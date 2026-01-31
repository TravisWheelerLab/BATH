/* Profile HMM construction from a multiple sequence alignment
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_mpi.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_msaweight.h"
#include "esl_msacluster.h"
#include "esl_stopwatch.h"
#include "esl_vectorops.h"

#ifdef HMMER_THREADS
#include <unistd.h>
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG	           *bg;
  P7_BUILDER       *bld;
} WORKER_INFO;

#ifdef HMMER_THREADS
typedef struct {
  int         nali;
  int         processed;
  ESL_MSA    *postmsa;
  ESL_MSA    *msa;
  ESL_SQ     *sq;
  P7_HMM     *hmm;
  int         force_single; /* FALSE by default,  TRUE if esl_opt_IsUsed(go, "--singlemx") ;  only matters for single sequences */
} WORK_ITEM;

typedef struct _pending_s {
  int         nali;
  ESL_MSA    *postmsa;
  ESL_MSA    *msa;
  ESL_SQ     *sq;
  P7_HMM     *hmm;
  struct _pending_s *next;
} PENDING_ITEM;
#endif /*HMMER_THREADS*/

struct cfg_s {

  char         *infile;         /* name of the input file we're building HMMs from  */
  int           fmt;        /* format code for inputfile */
  ESL_MSAFILE  *afp;            /* open infile if MSA */
  ESL_SQFILE   *sfp;            /* open infile is sequence */
  ESL_ALPHABET *abc;        /* digital alphabet */

  char         *hmmName;        /* hmm file name                           */
  char         *hmmfile;        /* file to write HMM to                    */
  FILE         *hmmfp;          /* HMM output file handle                  */

  int           nali;       /* which # alignment this is in file (only valid in serial mode)   */
  int           nnamed;     /* number of alignments that had their own names */
};

static char usage[]  = "[-options] <hmmfile_out> <msafile>";
static char banner[] = "profile HMM construction from multiple sequence alignments";

static void serial_loop  (WORKER_INFO *info, struct cfg_s *cfg, const ESL_GETOPTS *go);
#ifdef HMMER_THREADS
static void thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, struct cfg_s *cfg, const ESL_GETOPTS *go);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/

static int set_msa_name(struct cfg_s *cfg, char *errbuf, ESL_MSA *msa);
static void apply_fixed_gap_params(P7_HMM *hmm, double popen, double pextend);

/*
 * infile
 * hmmfile
 *
 */


/* p7_search_builder()
 * The usual version of hmmbuild, serial or threaded
 * For each MSA, build an HMM and save it.
 * 
 * A master can only return if it's successful. 
 * All errors are handled immediately and fatally with p7_Fail() or equiv.
 */
int
p7_search_builder(const ESL_GETOPTS *go, ESL_ALPHABET *abc, char *infile, char *hmmfile, int fmt)
{
  int              ncpus    = 0;
  int              infocnt  = 0;
  struct cfg_s     cfg;
  WORKER_INFO     *info     = NULL;
#ifdef HMMER_THREADS
  WORK_ITEM       *item     = NULL;
  ESL_THREADS     *threadObj= NULL;
  ESL_WORK_QUEUE  *queue    = NULL;
#endif
  double           popen;
  double           pextend;
  int              i;
  int              status;

  cfg.fmt     = fmt;
  cfg.infile  = infile;
  cfg.hmmfile = hmmfile;
  cfg.abc     = abc;
  cfg.hmmName = NULL;
  cfg.afp     = NULL;
  cfg.sfp     = NULL;
  cfg.hmmfp   = NULL; 

  cfg.nali       = 0;
  cfg.nnamed     = 0;

  if(fmt > 100) {
    status = esl_msafile_Open(&(cfg.abc), cfg.infile, NULL, cfg.fmt, NULL, &(cfg.afp));
    if (status != eslOK) esl_msafile_OpenFailure(cfg.afp, status);
  } else {
    status = esl_sqfile_Open(cfg.infile, cfg.fmt, NULL, &(cfg.sfp));
    if (status != eslOK) {
      if      (status == eslENOTFOUND) p7_Fail("Failed to open input sequence file %s for reading\n",          cfg.infile);
      else if (status == eslEFORMAT)   p7_Fail("Input sequence file %s is empty or misformatted\n",            cfg.infile);
      else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
      else if (status != eslOK)        p7_Fail("Unexpected error %d opening input sequence file %s\n", status, cfg.infile);
    } 
  }

  cfg.hmmfp = fopen(cfg.hmmfile, "w");
  if (cfg.hmmfp == NULL) p7_Fail("Failed to open HMM file %s for writing", cfg.hmmfile);

#ifdef HMMER_THREADS
  /* initialize thread data */
  ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
      queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus <= 0) ? 1 : ncpus;
  ESL_ALLOC(info, (ptrdiff_t) sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
  {
      info[i].bg = p7_bg_Create(cfg.abc);
      info[i].bld = p7_builder_Create(NULL, cfg.abc);

      if (info[i].bld == NULL)  p7_Fail("p7_builder_Create failed");

      //protein defaults
      popen   = esl_opt_IsUsed(go, "--popen")   ? esl_opt_GetReal(go, "--popen")   : 0.02;
      pextend = esl_opt_IsUsed(go, "--pextend") ? esl_opt_GetReal(go, "--pextend") : 0.4;

      /* Default matrix is stored in the --mx option, so it's always IsOn().
       * Check --mxfile first; then go to the --mx option and the default.
       */
      if ( esl_opt_IsUsed(go, "--singlemx") || cfg.sfp != NULL) {
        char  *mx      = esl_opt_GetString(go, "--mx");

        if (esl_opt_IsOn(go, "--mxfile")) status = p7_builder_SetScoreSystem (info[i].bld, esl_opt_GetString(go, "--mxfile"), NULL, popen, pextend, info[i].bg);
        else                              status = p7_builder_LoadScoreSystem(info[i].bld, mx,                                      popen, pextend, info[i].bg);
        if (status != eslOK) p7_Fail("Failed to set single query seq score system:\n%s\n", info[i].bld->errbuf);
      } else {
        if (esl_opt_IsUsed(go, "--popen") )  info[i].bld->popen   = popen;
        if (esl_opt_IsUsed(go, "--pextend")) info[i].bld->pextend = pextend;
      }

      /* special arguments for bathbuild */
      info[i].bld->fs         = (go != NULL && esl_opt_IsUsed (go, "--fs"))     ?  esl_opt_GetBoolean(go, "--fs")      : FALSE;
      info[i].bld->fsprob     = (go != NULL && esl_opt_IsUsed (go, "--fsprob")) ?  esl_opt_GetReal   (go, "--fsprob")  : 0.01;
      info[i].bld->ct         = (go != NULL && esl_opt_IsUsed (go, "--ct"))     ?  esl_opt_GetInteger(go, "--ct")      : 1;
      info[i].bld->w_len      = (go != NULL && esl_opt_IsOn (go, "--w_length")) ?  esl_opt_GetInteger(go, "--w_length"): -1;
      info[i].bld->w_beta     = (go != NULL && esl_opt_IsOn (go, "--w_beta"))   ?  esl_opt_GetReal   (go, "--w_beta")  : p7_DEFAULT_WINDOW_BETA;
      if ( info[i].bld->w_beta < 0 || info[i].bld->w_beta > 1  ) esl_fatal("Invalid window-length beta value\n");

#ifdef HMMER_THREADS
      info[i].queue = queue;
      if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
  }

#ifdef HMMER_THREADS
  for (i = 0; i < ncpus * 2; ++i)
  {
      ESL_ALLOC(item, sizeof(*item));

      item->nali      = 0;
      item->processed = FALSE;
      item->postmsa   = NULL;
      item->msa       = NULL;
      item->hmm       = NULL;

      status = esl_workqueue_Init(queue, item);
      if (status != eslOK) esl_fatal("Failed to add block to work queue");
  }
#endif

#ifdef HMMER_THREADS
  if (ncpus > 0)  thread_loop(threadObj, queue, &cfg, go);
  else            serial_loop(info, &cfg, go);
#else
  serial_loop(info, &cfg, go);
#endif

  for (i = 0; i < infocnt; ++i)
  {
      p7_bg_Destroy(info[i].bg);
      p7_builder_Destroy(info[i].bld);
  }

#ifdef HMMER_THREADS
  if (ncpus > 0)
  {
      esl_workqueue_Reset(queue);
      while (esl_workqueue_Remove(queue, (void **) &item) == eslOK)
      {
        free(item);
      }
      esl_workqueue_Destroy(queue);
      esl_threads_Destroy(threadObj);
  }
#endif

  free(info);

  if (cfg.afp)   esl_msafile_Close(cfg.afp);
  if (cfg.sfp)   esl_sqfile_Close(cfg.sfp);
  if (cfg.hmmfp) fclose(cfg.hmmfp);

  return eslOK;

 ERROR:
  return eslFAIL;
}


static void
serial_loop(WORKER_INFO *info, struct cfg_s *cfg, const ESL_GETOPTS *go)
{
  P7_BUILDER *bld         = NULL;
  ESL_MSA    *msa         = NULL;
  ESL_SQ     *sq          = NULL;
  P7_HMM     *hmm         = NULL;
  char        errmsg[eslERRBUFSIZE];
  int         status;


  if (cfg->fmt > 100) { 
    cfg->nali = 0;
    while ((status = esl_msafile_Read(cfg->afp, &msa)) != eslEOF)
      {
        if (status != eslOK) esl_msafile_ReadFailure(cfg->afp, status);
        cfg->nali++;  

        if ((status = set_msa_name(cfg, errmsg, msa)) != eslOK) p7_Fail("%s\n", errmsg); /* cfg->nnamed gets incremented in this call */


        /*         bg   new-HMM trarr gm   om  */
        if ( msa->nseq == 1 && esl_opt_IsUsed(go, "--singlemx")) {
          if ((status = esl_sq_FetchFromMSA(msa, 0, &sq)) != eslOK) p7_Fail("build failed: %s", bld->errbuf);
          if ((status = p7_SingleBuilder(info->bld, sq, info->bg, &hmm, NULL, NULL, NULL)) != eslOK) p7_Fail("build failed: %s", bld->errbuf);
          esl_sq_Destroy(sq);
          sq = NULL;
          hmm->eff_nseq = 1;
        } else {
          if ((status = p7_Builder(info->bld, msa, info->bg, &hmm, NULL, NULL, NULL, NULL )) != eslOK) p7_Fail("build failed: %s", bld->errbuf);

          //if not --singlemx, but the user set the popen/pextend flags, override the computed gap params now:
          if (info->bld->popen != -1 || info->bld->pextend != -1) {
            apply_fixed_gap_params(hmm, info->bld->popen, info->bld->pextend);
          }
        }
        if ((status = p7_hmmfile_WriteASCII(cfg->hmmfp, p7_BATH_3f, hmm)) != eslOK) p7_Fail("HMM save failed"); 
        p7_hmm_Destroy(hmm);
        esl_msa_Destroy(msa);
      }
  } else {
    
    esl_sqfile_SetDigital(cfg->sfp, cfg->abc);
    sq = esl_sq_CreateDigital(cfg->abc); 
    while ((status = esl_sqio_Read(cfg->sfp, sq)) != eslEOF)
      {  
         cfg->nali++;
         if (status != eslOK) p7_Fail("reading unaligned sequences from input file %s (%d)\n", cfg->infile, status); 
         if ((status = p7_SingleBuilder(info->bld, sq, info->bg, &hmm, NULL, NULL, NULL)) != eslOK) p7_Fail("build failed: %s", bld->errbuf);
     
         hmm->eff_nseq = 1;
         if ((status = p7_hmmfile_WriteASCII(cfg->hmmfp, p7_BATH_3f, hmm)) != eslOK) p7_Fail("HMM save failed"); 
         esl_sq_Reuse(sq);
         p7_hmm_Destroy(hmm);
      }
  }

  if(sq != NULL) esl_sq_Destroy(sq);
   
}

#ifdef HMMER_THREADS
static void
thread_loop(ESL_THREADS *obj, ESL_WORK_QUEUE *queue, struct cfg_s *cfg, const ESL_GETOPTS *go)
{
  int          status    = eslOK;
  int          sstatus   = eslOK;
  int          processed = 0;
  WORK_ITEM   *item;
  void        *newItem;

  int           next     = 1;
  PENDING_ITEM *top      = NULL;
  PENDING_ITEM *empty    = NULL;
  PENDING_ITEM *tmp      = NULL;

  char        errmsg[eslERRBUFSIZE];

   cfg->nali = 0;

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newItem);
  if (status != eslOK) esl_fatal("Work queue reader failed");
      
  /* Main loop: */
  item = (WORK_ITEM *) newItem;
  
  while (sstatus == eslOK) {
 
    if (cfg->fmt > 100) {
      item->sq = NULL;
      sstatus = esl_msafile_Read(cfg->afp, &item->msa);
      
      
      if (sstatus == eslOK) {
        item->nali = ++cfg->nali;
        if (set_msa_name(cfg, errmsg, item->msa) != eslOK) p7_Fail("%s\n", errmsg);
      }
      else if (sstatus == eslEOF && processed < cfg->nali) sstatus = eslOK;
      else if (sstatus != eslEOF) 
        esl_msafile_ReadFailure(cfg->afp, sstatus);

      if (sstatus == eslOK) {
        item->force_single = esl_opt_IsUsed(go, "--singlemx");
        status = esl_workqueue_ReaderUpdate(queue, item, &newItem);
        if (status != eslOK) esl_fatal("Work queue reader failed");
       
        /* process any results */
        item = (WORK_ITEM *) newItem;
        
        if (item->processed == TRUE) {
	  ++processed;

	  /* try to keep the input output order the same */
	  if (item->nali == next) {
        if ((sstatus = p7_hmmfile_WriteASCII(cfg->hmmfp, p7_BATH_3f, item->hmm)) != eslOK) p7_Fail("HMM save failed");

	    p7_hmm_Destroy(item->hmm);
	    esl_msa_Destroy(item->msa);
	    esl_msa_Destroy(item->postmsa);

	    ++next;

	    /* output any pending msa as long as the order
	     * remains the same as read in.
	     */
	    while (top != NULL && top->nali == next) {
          if ((sstatus = p7_hmmfile_WriteASCII(cfg->hmmfp, p7_BATH_3f, top->hmm)) != eslOK) p7_Fail("HMM save failed");

	      p7_hmm_Destroy(top->hmm);
	      esl_msa_Destroy(top->msa);
	      esl_msa_Destroy(top->postmsa);

	      tmp = top;
	      top = tmp->next;

	      tmp->next = empty;
	      empty     = tmp;
	    
	      ++next;
	    }
	  } else {
	    /* queue up the msa so the sequence order is the same in
	     * the .sto and .hmm
	     */
	    if (empty != NULL) {
	      tmp   = empty;
	      empty = tmp->next;
	    } else {
	      ESL_ALLOC(tmp, sizeof(PENDING_ITEM));
	    }

	    tmp->nali     = item->nali;
	    tmp->hmm      = item->hmm;
	    tmp->msa      = item->msa;
	    tmp->postmsa  = item->postmsa;

	    /* add the msa to the pending list */
	    if (top == NULL || tmp->nali < top->nali) {
	      tmp->next = top;
	      top       = tmp;
	    } else {
	      PENDING_ITEM *ptr = top;
	      while (ptr->next != NULL && tmp->nali > ptr->next->nali) {
	        ptr = ptr->next;
	      }
	      tmp->next = ptr->next;
	      ptr->next = tmp;
	    }
	  }

	  item->nali      = 0;
	  item->processed = FALSE;
	  item->hmm       = NULL;
	  item->msa       = NULL;
	  item->postmsa   = NULL;
        }
      }
    } else {  // input is unaligned sequences
      item->msa = NULL;
      esl_sqfile_SetDigital(cfg->sfp, cfg->abc);
      item->sq = esl_sq_CreateDigital(cfg->abc); 
      sstatus = esl_sqio_Read(cfg->sfp, item->sq);

      if (sstatus == eslOK) item->nali = ++cfg->nali;
      else if (sstatus == eslEOF) {
          if (processed < cfg->nali) sstatus = eslOK;
          if (item->sq != NULL) { 
            esl_sq_Destroy(item->sq);
            item->sq = NULL;
          }
      }
      else if (sstatus != eslEOF) 
        if (status != eslOK) p7_Fail("reading unaligned sequences from input file %s (%d)\n", cfg->infile, sstatus); 
	
       if (sstatus == eslOK) {
        item->force_single = esl_opt_IsUsed(go, "--singlemx");
        status = esl_workqueue_ReaderUpdate(queue, item, &newItem);
        if (status != eslOK) esl_fatal("Work queue reader failed");
        
        /* process any results */
        item = (WORK_ITEM *) newItem;
        if (item->processed == TRUE) {
	  ++processed;
       
	  /* try to keep the input output order the same */
	  if (item->nali == next) {
        if ((sstatus = p7_hmmfile_WriteASCII(cfg->hmmfp, p7_BATH_3f, item->hmm)) != eslOK) p7_Fail("HMM save failed");

	    if(item->hmm != NULL) p7_hmm_Destroy(item->hmm);
	    if(item->sq  != NULL) esl_sq_Destroy(item->sq);

	    ++next;

	    /* output any pending msa as long as the order
	     * remains the same as read in.
	     */
	    while (top != NULL && top->nali == next) {
          if ((sstatus = p7_hmmfile_WriteASCII(cfg->hmmfp, p7_BATH_3f, top->hmm)) != eslOK) p7_Fail("HMM save failed");

	      if(item->hmm != NULL) p7_hmm_Destroy(top->hmm);
	      if(item->sq  != NULL) esl_sq_Destroy(top->sq);

	      tmp = top;
	      top = tmp->next;

	      tmp->next = empty;
	      empty     = tmp;
	    
	      ++next;
	    }
	  } else {
	    /* queue up the msa so the sequence order is the same in
	     * the .sto and .hmm
	     */
	    if (empty != NULL) {
	      tmp   = empty;
	      empty = tmp->next;
	    } else {
	      ESL_ALLOC(tmp, sizeof(PENDING_ITEM));
	    }

	    tmp->nali     = item->nali;
	    tmp->hmm      = item->hmm;
	    tmp->sq       = item->sq;

	    /* add the msa to the pending list */
	    if (top == NULL || tmp->nali < top->nali) {
	      tmp->next = top;
	      top       = tmp;
	    } else {
	      PENDING_ITEM *ptr = top;
	      while (ptr->next != NULL && tmp->nali > ptr->next->nali) {
	        ptr = ptr->next;
	      }
	      tmp->next = ptr->next;
	      ptr->next = tmp;
	    }
	  }
    
	  item->nali      = 0;
	  item->processed = FALSE;
	  item->hmm       = NULL;
	  item->sq       = NULL;
	  item->postmsa   = NULL;
        }
      }
    }
  }
  if (top != NULL) esl_fatal("Top is not empty\n");

  while (empty != NULL) {
    tmp   = empty;
    empty = tmp->next;
    free(tmp);
  }

  status = esl_workqueue_ReaderUpdate(queue, item, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF)
    {
      /* wait for all the threads to complete */
      esl_threads_WaitForFinish(obj);
      esl_workqueue_Complete(queue);  
    }
  return;

 ERROR:
  p7_Fail("thread_loop failed: memory allocation problem");
}

static void 
pipeline_thread(void *arg)
{
  
  int           workeridx;
  int           status;

  WORK_ITEM    *item;
  void         *newItem;

  WORKER_INFO  *info;
  ESL_THREADS  *obj;
  ESL_SQ     *sq          = NULL;

  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newItem);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  item = (WORK_ITEM *) newItem;

  while (item->msa != NULL)
    {
      if ( item->msa->nseq == 1 && item->force_single) {
   
        status = esl_sq_FetchFromMSA(item->msa, 0, &sq);
        if (status != eslOK) p7_Fail("build failed: %s", info->bld->errbuf);

        status = p7_SingleBuilder(info->bld, sq, info->bg, &item->hmm, NULL, NULL, NULL);
        if (status != eslOK) p7_Fail("build failed: %s", info->bld->errbuf);

        esl_sq_Destroy(sq);
        sq = NULL;
        item->hmm->eff_nseq = 1;
      } else {
    
        status = p7_Builder(info->bld, item->msa, info->bg, &item->hmm, NULL, NULL, NULL, &item->postmsa);
        if (status != eslOK) p7_Fail("build failed: %s", info->bld->errbuf);
        //if not --singlemx, but the user set the popen/pextend flags, override the computed gap params now:
        if (info->bld->popen != -1 || info->bld->pextend != -1) {
          apply_fixed_gap_params(item->hmm, info->bld->popen, info->bld->pextend);
        }
      }
      item->processed = TRUE;

      status = esl_workqueue_WorkerUpdate(info->queue, item, &newItem);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      item = (WORK_ITEM *) newItem;
     
    }
    while (item->sq != NULL)
    {
        
        status = p7_SingleBuilder(info->bld, item->sq, info->bg, &item->hmm, NULL, NULL, NULL);
        if (status != eslOK) p7_Fail("build failed: %s", info->bld->errbuf);

        item->hmm->eff_nseq = 1;
        item->processed = TRUE;
      
      status = esl_workqueue_WorkerUpdate(info->queue, item, &newItem);
      if (status != eslOK) esl_fatal("Work queue worker failed");

      item = (WORK_ITEM *) newItem;
    }
   
  status = esl_workqueue_WorkerUpdate(info->queue, item, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
  
  return;
}
#endif   /* HMMER_THREADS */
 



/* set_msa_name() 
 * Make sure the alignment has a name; this name will
 * then be transferred to the model.
 * 
 * We can only do this for a single alignment in a file. For multi-MSA
 * files, each MSA is required to have a name already.
 *
 * Priority is:
 *      1. Use -n <name> if set, overriding any name the alignment might already have. 
 *      2. Use alignment's existing name, if non-NULL.
 *      3. Make a name, from alignment file name without path and without filename extension 
 *         (e.g. "/usr/foo/globins.slx" gets named "globins")
 * If none of these succeeds, return <eslEINVAL>.
 *         
 * If a multiple MSA database (e.g. Stockholm/Pfam), and we encounter
 * an MSA that doesn't already have a name, return <eslEINVAL> if nali > 1.
 * (We don't know we're in a multiple MSA database until we're on the second
 * alignment.)
 * 
 * Because we can't tell whether we've got more than one
 * alignment 'til we're on the second one, these fatal errors
 * only happen after the first HMM has already been built.
 * Oh well.
 */
static int
set_msa_name(struct cfg_s *cfg, char *errbuf, ESL_MSA *msa)
{
  char *name = NULL;
  int   status;

  if (cfg->nali == 1) /* first (only?) HMM in file: */
    {
      if  (cfg->hmmName != NULL)
	{
	  if ((status = esl_msa_SetName(msa, cfg->hmmName, -1)) != eslOK) return status;
	}
      else if (msa->name != NULL) 
	{
	  cfg->nnamed++;
	}
      else if (cfg->afp->bf->filename)
	{
	  if ((status = esl_FileTail(cfg->afp->bf->filename, TRUE, &name)) != eslOK) return status; /* TRUE=nosuffix */	  
	  if ((status = esl_msa_SetName(msa, name, -1))                    != eslOK) return status;
	  free(name);
	}
      else ESL_FAIL(eslEINVAL, errbuf, "Failed to set model name: msa has no name, no msa filename, and no -n");
    }
  else 
    {
      if (cfg->hmmName   != NULL) ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. You can't use -n with an alignment database.");
      else if (msa->name != NULL) cfg->nnamed++;
      else                        ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment in a multi MSA file; failed on #%d", cfg->nali+1);

      /* special kind of failure: the *first* alignment didn't have a name, and we used the filename to
       * construct one; now that we see a second alignment, we realize this was a boo-boo*/
      if (cfg->nnamed != cfg->nali)            ESL_FAIL(eslEINVAL, errbuf, "Oops. Wait. I need name annotation on each alignment in a multi MSA file; first MSA didn't have one");
    }
  return eslOK;
}


void
apply_fixed_gap_params(P7_HMM *hmm, double popen, double pextend){
  int k;
  for (k = 0; k <= hmm->M; k++)
  {
     if (popen != -1) {
        hmm->t[k][p7H_MM] = 1.0 - 2 * popen;
        hmm->t[k][p7H_MI] = popen;
        hmm->t[k][p7H_MD] = popen;
     }
     if (pextend != -1) {
        hmm->t[k][p7H_IM] = 1.0 - pextend;
        hmm->t[k][p7H_II] = pextend;
        hmm->t[k][p7H_DM] = 1.0 - pextend;
        hmm->t[k][p7H_DD] = pextend;
     }
  }

  /* Deal w/ special stuff at node M, overwriting a little of what we
   * just did.
   */
  if (popen != -1) {
    hmm->t[hmm->M][p7H_MM] = 1.0 - popen;
  }
  hmm->t[hmm->M][p7H_MD] = 0.;
  hmm->t[hmm->M][p7H_DM] = 1.0;
  hmm->t[hmm->M][p7H_DD] = 0.;

}

