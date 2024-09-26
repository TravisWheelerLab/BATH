/* bathsearch: search protein profile HMM(s) against a DNA sequence database. */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sq.h"
#include "esl_sqio.h"
#include "esl_stopwatch.h"
#include "esl_gencode.h"

#ifdef HMMER_THREADS
#include "esl_threads.h"
#include "esl_workqueue.h"
#endif /*HMMER_THREADS*/

#include "hmmer.h"

/* set the max residue count to 1/4 meg when reading a block */
#define BATH_MAX_RESIDUE_COUNT (1024 * 256)  /* 1/4 Mb */


typedef struct {
#ifdef HMMER_THREADS
  ESL_WORK_QUEUE   *queue;
#endif /*HMMER_THREADS*/
  P7_BG            *bg;	         /* null model                                                        */
  ESL_SQ           *ntsq;        /* DNA target sequence                                               */
  P7_PIPELINE      *pli;         /* work pipeline                                                     */
  P7_TOPHITS       *th;          /* top hit results                                                   */
  P7_OPROFILE      *om;          /* optimized query profile                                           */
  P7_PROFILE       *gm;		 /* non-optimized query profile                                       */
  P7_FS_PROFILE    *gm_fs;       /* non optimized frameshift query profile                            */
  P7_SCOREDATA     *scoredata;   /* used to create DNA windows from ORFs                              */
  ESL_GENCODE      *gcode;       /* used for translating ORFs                                         */
  ESL_GENCODE_WORKSTATE *wrk1;   /* used for intitial translation of taget DNA to ORFs                */ 
  ESL_GENCODE_WORKSTATE *wrk2;   /* used for secondary translation of DNA window for bias calcultaion */
} WORKER_INFO;


/* items used to keep track ot orignal taget lengths */
typedef struct {
  int    id;         /* internal sequence ID  */
  int    length;     /* length of sequence */
} ID_LENGTH;

typedef struct {
  ID_LENGTH  *id_lengths;
  int        count;
  int        size;
} ID_LENGTH_LIST;

static ID_LENGTH_LIST* init_id_length( int size );
static void            destroy_id_length( ID_LENGTH_LIST *list );
static int             add_id_length(ID_LENGTH_LIST *list, int id, int L);
static int             assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list);

#define REPOPTS     "-E,-T"//--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT"//--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL

static ESL_OPTIONS options[] = {
  /* name             type            default    env          range      toggles reqs  incomp          help                                                                           docgroup*/
  { "-h",             eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "show brief help on version and usage",                                        1 },
  /* Control of output */
  { "-o",             eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "direct output to file <f>, not stdout",                                       2 },
  { "--tblout",       eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "save parseable table of hits to file <f>",                                    2 },
  { "--fstblout",     eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "save table of frameshift locations to file <f>",                              2 },
  { "--hmmout",       eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "if input is alignment(s) or sequence(s) write produced hmms to file <f>",     2 },
  { "--splice",       eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,"--nofs", NULL,          "enable spliced alignments (requires SSI index 'esl-sftech --index <seqdb>')", 2 },
  { "--acc",          eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "prefer accessions over names in output",                                      2 },
  { "--noali",        eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "don't output alignments, so output is smaller",                               2 },
  { "--notrans",      eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "don't show the translated DNA sequence in  alignment",                        2 }, 
  { "--frameline",    eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "include frame of each codon in  alignment",                                   2 },
  { "--cigar",        eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,"--tblout", NULL,        "include alignment CIGAR string in table output (with --tblout)",              2 },
  { "--notextw",      eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL,"--textw",       "unlimit ASCII text output line width",                                        2 },
  { "--textw",        eslARG_INT,    "150",      NULL,       "n>=150",   NULL,   NULL,"--notextw",     "set max width of ASCII text output lines",                                    2 },
  /* Control of scoring system */
  { "--singlemx",     eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "use substitution score matrix w/ single-sequence MSA-format inputs",          3 },
  { "--popen",        eslARG_REAL,   "0.02",     NULL,       "0<=x<0.5", NULL,   NULL, NULL,           "gap open probability",                                                        3 },
  { "--pextend",      eslARG_REAL,   "0.4",      NULL,       "0<=x<1",   NULL,   NULL, NULL,           "gap extend probability",                                                      3 },
  { "--mx",           eslARG_STRING, "BLOSUM62", NULL,        NULL,      NULL,   NULL,"--mxfile",      "substitution score matrix choice (of some built-in matrices)",                3 },
  { "--mxfile",       eslARG_INFILE,  NULL,      NULL,        NULL,      NULL,   NULL,"--mx",          "read substitution score matrix from file <f>",                                3 },
  /* Control of reporting and inclusion thresholds */
  { "-E",             eslARG_REAL,   "10.0",     NULL,       "x>0",      NULL,   NULL, REPOPTS,        "report sequences <= this E-value threshold in output",                        4 },
  { "-T",             eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,   NULL, REPOPTS,        "report sequences >= this score threshold in output",                          4 },
  { "--incE",         eslARG_REAL,   "0.01",     NULL,       "x>0",      NULL,   NULL, INCOPTS,        "consider sequences <= this E-value threshold as significant",                 4 },
  { "--incT",         eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,   NULL, INCOPTS,        "consider sequences >= this score threshold as significant",                   4 },
  /* input formats */
  { "--qformat",      eslARG_STRING,  NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "assert query is in format <s> (can be seq or msa format)",                    5 },
  { "--qsingle_seqs", eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "force query to be read as individual sequences, even if in an msa format",    5 },
  { "--tformat",      eslARG_STRING,  NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "assert target <seqfile> is in format <s>: no autodetection",                  5 },

  /* Control of acceleration pipeline */
  { "--max",          eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"--F1,--F2,--F3","turn all heuristic filters off (less speed, more power)",                     7 },
  { "--F1",           eslARG_REAL,   "0.02",     NULL,        NULL,      NULL,   NULL,"--max",         "stage 1 (MSV) threshold: promote hits w/ P <= F1",                            7 },
  { "--F2",           eslARG_REAL,   "1e-3",     NULL,        NULL,      NULL,   NULL,"--max",         "stage 2 (Vit) threshold: promote hits w/ P <= F2",                            7 },
  { "--F3",           eslARG_REAL,   "1e-5",     NULL,        NULL,      NULL,   NULL,"--max",         "stage 3 (Fwd) threshold: promote hits w/ P <= F3",                            7 },
  { "--nobias",       eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL,"--max",         "turn off composition bias filter",                                            7 },
  { "--nonull2",      eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "turn off biased composition score corrections",                               7 },
  { "--fsonly",       eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"--nofs",        "send all potential hits to the frameshift aware pipeline",                    7 },
  { "--nofs",         eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"--fsonly",      "send all potential hits to the non-frameshift aware pipeline",                7 },
/* Other options */
  { "-Z",             eslARG_REAL,    FALSE,     NULL,       "x>=0",     NULL,   NULL, NULL,           "set database size (Megabases) to <x> for E-value calculations",               12 }, 
  { "--seed",         eslARG_INT,    "42",       NULL,       "n>=0",     NULL,   NULL, NULL,           "set RNG seed to <n> (if 0: one-time arbitrary seed)",                         12 },
  { "--w_beta",       eslARG_REAL,    NULL,      NULL,       "0>=x<=1",  NULL,   NULL, NULL,           "tail mass at which window length is determined",                              12 },
  { "--w_length",     eslARG_INT,     NULL,      NULL,       "x>=4",      NULL,   NULL, NULL,           "window length - essentially max expected hit length" ,                       12 },
  #ifdef HMMER_THREADS 
  { "--block_length", eslARG_INT,     NULL,      NULL,       "n>=50000", NULL,   NULL, NULL,           "length of blocks read from target database (threaded) ",                      12 },
  { "--cpu",          eslARG_INT,     p7_NCPU,  "HMMER_NCPU","n>=0",     NULL,   NULL, CPUOPTS,        "number of parallel CPU workers to use for multithreads",                      12 },
#endif
  /* Translation options */ 
  { "--ct",           eslARG_INT,    "1",        NULL,        NULL,      NULL,   NULL, NULL,           "use alt genetic code of NCBI translation table (see end of help)",            15 },
  { "-l",             eslARG_INT,    "20",       NULL,        NULL,      NULL,   NULL, NULL,           "minimum ORF length",                                                          15 },
  { "-m",             eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"-M",            "ORFs must initiate with AUG only",                                            15 },
  { "-M",             eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"-m",            "ORFs must start with allowed initiation codon",                               15 },
  { "--strand",       eslARG_STRING, "both",     NULL,        NULL,      NULL,   NULL, NULL,           "translate only forward strand 'plus' or reverse complement strand 'minus'",   15 },
  
  /* Restrict search to subset of database - hidden because these flags are
   *   (a) currently for internal use
   *   (b) probably going to change
   */
  { "--restrictdb_stkey", eslARG_STRING,"0",     NULL,        NULL,      NULL,   NULL, NULL,           "Search starts at the sequence with name <s> ",                                99 },
  { "--restrictdb_n",     eslARG_INT,   "-1",    NULL,        NULL,      NULL,   NULL, NULL,           "Search <j> target sequences (starting at --restrictdb_stkey)",                99 },
  { "--ssifile",          eslARG_STRING, NULL,   NULL,        NULL,      NULL,   NULL, NULL,           "restrictdb_x values require ssi file. Override default to <s>",               99 },

  /* stage-specific window length used for bias composition estimate,
   * hidden because they are confusing/expert options. May drag them out
   * into the daylight eventually
   */
  { "--B1",           eslARG_INT,    "110",      NULL,        NULL,      NULL,  NULL,"--max,--nobias", "window length for biased-composition modifier (SSV)",                         99 },
  { "--B2",           eslARG_INT,    "240",      NULL,        NULL,      NULL,  NULL,"--max,--nobias", "window length for biased-composition modifier (Vit)",                         99 },
  { "--B3",           eslARG_INT,    "1000",     NULL,        NULL,      NULL,  NULL,"--max,--nobias", "window length for biased-composition modifier (Fwd)",                         99 },
 
  /* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so that it doesn't print to help*/
  { "--domZ",         eslARG_REAL,    FALSE,     NULL,       "x>0",      NULL,  NULL, NULL,            "Not used",                                                                    99 },
  { "--domE",         eslARG_REAL,   "10.0",     NULL,       "x>0",      NULL,  NULL, DOMREPOPTS,      "Not used",                                                                    99 },
  { "--domT",         eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,  NULL, DOMREPOPTS,      "Not used",                                                                    99 },
  { "--incdomE",      eslARG_REAL,   "0.01",     NULL,       "x>0",      NULL,  NULL, INCDOMOPTS,      "Not used",                                                                    99 },
  { "--incdomT",      eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,  NULL, INCDOMOPTS,      "Not used",                                                                    99 },
  { "--crick",        eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,  NULL, NULL,            "only translate top strand",                                                   99 },
  { "--watson",       eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,  NULL, NULL,            "only translate bottom strand",                                                99 }, 
  { "--fs",           eslARG_REAL,   "0.01",     NULL,       "0<=x<=1",  NULL,  NULL, NULL,            "set the frameshift probabilty",                                               99 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

 /* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads).
 */
struct cfg_s {
  char            *dbfile;            /* target sequence database file                   */
  char            *queryfile;         /* query file (hmm, fasta, or some MSA)            */
  int              qfmt;  

  char             *firstseq_key;     /* name of the first sequence in the restricted db range */
  int              n_targetseq;       /* number of sequences in the restricted range */
};

static char usage[]  = "[options] <hmm, msa, or seq file> <seqdb>";
static char banner[] = "search protein profile(s) against DNA sequence database";

static int  serial_master(ESL_GETOPTS *go, struct cfg_s *cfg);
static int  serial_loop  (WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs);

#ifdef HMMER_THREADS
#define BLOCK_SIZE 1000

static int  thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseq);
static void pipeline_thread(void *arg);
#endif /*HMMER_THREADS*/


static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_hmmfile, char **ret_seqfile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK)  { if (printf("Failed to process environment: %s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK)  { if (printf("Failed to parse command line: %s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 100); /* 1= group; 2 = indentation; 100=textwidth*/

      if (puts("\nOptions directing output:")                                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 100); 

      if (puts("\nOptions controlling translation:")                      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 15, 2, 100); 

      if (puts("\nOptions controlling reporting and inclusion thresholds:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 100); 

      if (puts("\nOptions controlling acceleration heuristics:")             < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 100); 

      if (puts("\nOptions setting input formats:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 100);

      if (puts("\nOptions handling single sequence inputs:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 100);

      if (puts("\nOther expert options:")                                    < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 12, 2, 100); 
      
      if (puts("\nAvailable NCBI genetic code tables (for --ct <id>):")        < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_gencode_DumpAltCodeTable(stdout);

	  exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 2)     { if (puts("Incorrect number of command line arguments.")      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_hmmfile = esl_opt_GetArg(go, 1)) == NULL)  { if (puts("Failed to get <hmmfile> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_seqfile = esl_opt_GetArg(go, 2)) == NULL)  { if (puts("Failed to get <seqdb> argument on command line")   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* Validate any attempted use of stdin streams */
  if (strcmp(*ret_hmmfile, "-") == 0 && strcmp(*ret_seqfile, "-") == 0) 
    { if (puts("Either <hmmfile> or <seqdb> may be '-' (to read from stdin), but not both.") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  if (puts("\nwhere most common options are:")                                 < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80); /* 1= group; 2 = indentation; 80=textwidth*/
  if (printf("\nTo see more help on available options, do %s -h\n\n", argv[0]) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

static int
output_header(FILE *ofp, const ESL_GETOPTS *go, char *hmmfile, char *seqfile)
{
  
  if (                                                         fprintf(ofp, "# query HMM file:                                %s\n", hmmfile)                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (                                                         fprintf(ofp, "# target sequence database:                      %s\n", seqfile)                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (	                                                       fprintf(ofp, "# frameshift probability:                        %f\n", esl_opt_GetReal(go, "--fs"))                      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (                                                         fprintf(ofp, "# codon translation table:                       %d\n", esl_opt_GetInteger(go, "--ct"))                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")                              && fprintf(ofp, "# output directed to file:                       %s\n",      esl_opt_GetString(go, "-o"))                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")                        && fprintf(ofp, "# per-seq hits tabular output:                   %s\n",      esl_opt_GetString(go, "--tblout"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fstblout")                      && fprintf(ofp, "# frameshift tabular output:                     %s\n",      esl_opt_GetString(go, "--fstblout"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--hmmout")                        && fprintf(ofp, "# hmm output:                                    %s\n",      esl_opt_GetString(go, "--hmmout"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--splice")                        && fprintf(ofp, "# enable spliced alignments:                     yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--acc")                           && fprintf(ofp, "# prefer accessions over names:                  yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--noali")                         && fprintf(ofp, "# show alignments in output:                     no\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notextw")                       && fprintf(ofp, "# max ASCII text line length:                    unlimited\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--textw")                         && fprintf(ofp, "# max ASCII text line length:                    %d\n",      esl_opt_GetInteger(go, "--textw"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--notrans")                       && fprintf(ofp, "# show translated DNA sequence:                  no\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--singlemx")                      && fprintf(ofp, "# Use score matrix for 1-seq MSAs:               on\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--popen")                         && fprintf(ofp, "# gap open probability:                          %f\n",      esl_opt_GetReal  (go, "--popen"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pextend")                       && fprintf(ofp, "# gap extend probability:                        %f\n",      esl_opt_GetReal  (go, "--pextend"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--mx")                            && fprintf(ofp, "# subst score matrix (built-in):                 %s\n",      esl_opt_GetString(go, "--mx"))               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--mxfile")                        && fprintf(ofp, "# subst score matrix (file):                     %s\n",      esl_opt_GetString(go, "--mxfile"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-E")                              && fprintf(ofp, "# sequence reporting threshold:       E-value <= %g\n",      esl_opt_GetReal(go, "-E"))                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-T")                              && fprintf(ofp, "# sequence reporting threshold:         score >= %g\n",      esl_opt_GetReal(go, "-T"))                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incE")                          && fprintf(ofp, "# sequence inclusion threshold:       E-value <= %g\n",      esl_opt_GetReal(go, "--incE"))               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--incT")                          && fprintf(ofp, "# sequence inclusion threshold:         score >= %g\n",      esl_opt_GetReal(go, "--incT"))               < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--max")                           && fprintf(ofp, "# Max sensitivity mode:                          on [all heuristic filters off]\n")                       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F1")                            && fprintf(ofp, "# MSV filter P threshold:                     <= %g\n",      esl_opt_GetReal(go, "--F1"))                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F2")                            && fprintf(ofp, "# Vit filter P threshold:                     <= %g\n",      esl_opt_GetReal(go, "--F2"))                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--F3")                            && fprintf(ofp, "# Fwd filter P threshold:                     <= %g\n",      esl_opt_GetReal(go, "--F3"))                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")                        && fprintf(ofp, "# biased composition HMM filter:                 off\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nonull2")                       && fprintf(ofp, "# null2 bias corrections:                        off\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fsonly")                        && fprintf(ofp, "# Use only the frameshift aware pipeline\n")                                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--nofs")                          && fprintf(ofp, "# Use only the non-frameshift aware pipeline\n")                                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_stkey")              && fprintf(ofp, "# Restrict db to start at seq key:               %s\n",      esl_opt_GetString(go, "--restrictdb_stkey")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_n")                  && fprintf(ofp, "# Restrict db to # target seqs:                  %d\n",      esl_opt_GetInteger(go, "--restrictdb_n"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--ssifile")                       && fprintf(ofp, "# Override ssi file to:                          %s\n",      esl_opt_GetString(go, "--ssifile"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "-Z")                              && fprintf(ofp, "# database size is set to:                       %.1f Mb\n", esl_opt_GetReal(go, "-Z"))                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0               && fprintf(ofp, "# random number seed:                            one-time arbitrary\n")                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                                                  fprintf(ofp, "# random number seed set to:                     %d\n",      esl_opt_GetInteger(go, "--seed"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--qformat")                       && fprintf(ofp, "# query format asserted:                         %s\n",      esl_opt_GetString(go, "--qformat"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--qsingle_seqs")                  && fprintf(ofp, "# query contains individual seqs:                on\n")                                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--tformat")                       && fprintf(ofp, "# targ <seqfile> format asserted:                %s\n",      esl_opt_GetString(go, "--tformat"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_beta")                        && fprintf(ofp, "# window length beta value:                      %g\n",      esl_opt_GetReal(go, "--w_beta"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_length")                      && fprintf(ofp, "# window length :                                %d\n",      esl_opt_GetInteger(go, "--w_length"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
#ifdef HMMER_THREADS
  if (esl_opt_IsUsed(go, "--cpu")                           && fprintf(ofp, "# number of worker threads:                      %d\n",      esl_opt_GetInteger(go, "--cpu"))             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");  
#endif
  if (esl_opt_IsUsed(go, "-l")                              && fprintf(ofp, "# minimum ORF length:                            %d\n",      esl_opt_GetInteger(go, "-l"))                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-m")                              && fprintf(ofp, "# ORFs must initiate with AUG only:              yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-M")                              && fprintf(ofp, "# ORFs must start with allowed initiation codon: yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--strand")) {
     if     (!strcmp(esl_opt_GetString(go, "--strand"), "plus")   && fprintf(ofp, "# only translate the forward strand:             yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     else if(!strcmp(esl_opt_GetString(go, "--strand"), "minus")  && fprintf(ofp, "# only translate the reverse complement strand:  yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
     else if(!strcmp(esl_opt_GetString(go, "--strand"), "both")   && fprintf(ofp, "# translate both strands:                        yes\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (fprintf(ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n")                                                                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go       = NULL;	
  struct cfg_s     cfg;        
  int              status   = eslOK;

  impl_Init();                  /* processor specific initialization */
  p7_FLogsumInit();		/* we're going to use table-driven Logsum() approximations at times */

  /* Initialize what we can in the config structure (without knowing the alphabet yet) */
  cfg.queryfile    = NULL;
  cfg.dbfile       = NULL;
  cfg.qfmt         = eslSQFILE_UNKNOWN;
  cfg.firstseq_key = NULL;
  cfg.n_targetseq  = -1;
  process_commandline(argc, argv, &go, &cfg.queryfile, &cfg.dbfile);    

  if (esl_opt_IsOn(go, "--qformat")) { /* is this an msa or a single sequence file? */
    cfg.qfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--qformat")); // try single sequence format
    if (cfg.qfmt == eslSQFILE_UNKNOWN) {
      p7_Fail("%s is not a recognized input file format\n", esl_opt_GetString(go, "--qformat"));
    } else { /* disallow target-only formats */
      if (cfg.qfmt == eslSQFILE_NCBI    || cfg.qfmt == eslSQFILE_DAEMON ||
          cfg.qfmt == eslSQFILE_HMMPGMD || cfg.qfmt == eslSQFILE_FMINDEX )
        p7_Fail("%s is not a valid query format\n", esl_opt_GetString(go, "--qformat"));
    }
  }


  /* is the range restricted? */

#ifndef eslAUGMENT_SSI
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")  || esl_opt_IsUsed(go, "--ssifile")  )
    p7_Fail("Unable to use range-control options unless an SSI index file is available. See 'esl_sfetch --index'\n");
#else
  if (esl_opt_IsUsed(go, "--restrictdb_stkey") )
    if ((cfg.firstseq_key = esl_opt_GetString(go, "--restrictdb_stkey")) == NULL)  p7_Fail("Failure capturing --restrictdb_stkey\n");
  if (esl_opt_IsUsed(go, "--restrictdb_n") )
    cfg.n_targetseq = esl_opt_GetInteger(go, "--restrictdb_n");

  if ( cfg.n_targetseq != -1 && cfg.n_targetseq < 1 )
    p7_Fail("--restrictdb_n must be >= 1\n");

#endif

  status = serial_master(go, &cfg);

  esl_getopts_Destroy(go);

  return status;
}

/* create a set of ORFs for each DNA target sequence */
static int
do_sq_by_sequences(ESL_GENCODE *gcode, ESL_GENCODE_WORKSTATE *wrk, ESL_SQ *sq)
{
      esl_gencode_ProcessStart(gcode, wrk, sq);
      esl_gencode_ProcessPiece(gcode, wrk, sq);
      esl_gencode_ProcessEnd(wrk, sq);

  return eslOK;
}

/* query file open functions */
static int
bath_open_hmm_file(struct cfg_s *cfg,  P7_HMMFILE **hfp, char *errbuf, ESL_ALPHABET **abc, P7_HMM **hmm   ) {

  int status = p7_hmmfile_OpenE(cfg->queryfile, NULL, hfp, errbuf);

  if (status == eslENOTFOUND) {
    p7_Fail("File existence/permissions problem in trying to open query file %s.\n%s\n", cfg->queryfile, errbuf);
  } else if (status == eslOK) {
    //Successfully opened HMM file
    status = p7_hmmfile_Read(*hfp, abc, hmm);
    if (status != eslOK) p7_Fail("Error reading hmm from file %s (%d)\n", cfg->queryfile, status);
  }
    return status;
}

static int
bath_open_msa_file(struct cfg_s *cfg,  ESL_MSAFILE **qfp_msa, ESL_ALPHABET **abc, ESL_MSA **msa  ) {

  int status = esl_msafile_Open(abc, cfg->queryfile, NULL, cfg->qfmt, NULL, qfp_msa);

  if (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open query file %s.\n", cfg->queryfile);
  if (status == eslOK) {
    status = esl_msafile_Read(*qfp_msa, msa);
  }
  return status;
}

static int
bath_open_seq_file (struct cfg_s *cfg, ESL_SQFILE **qfp_sq, ESL_ALPHABET **abc, ESL_SQ **qsq, int used_qsingle_seqs) {

  int status = esl_sqfile_Open(cfg->queryfile, cfg->qfmt, NULL, qfp_sq);

  if (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open query file %s.\n", cfg->queryfile);
    if (status == eslOK) {
      if (*abc == NULL) {
        int q_type = eslUNKNOWN;
        status = esl_sqfile_GuessAlphabet(*qfp_sq, &q_type);
        if (  (*qfp_sq)->format == eslSQFILE_FASTA  /* we've guessed or been told it's a single sequence fasta file */
              && status == eslEFORMAT /* format error most likely to be due to presence of a gap character, so it's really an afa file */
              && used_qsingle_seqs  /* we were instructed to treat the input as single seqs, so override the fasta guess/instruction, and force single-sequence handling of afa file */
           ) {
           esl_sqfile_Close(*qfp_sq);
           status = esl_sqfile_Open(cfg->queryfile, eslMSAFILE_AFA, NULL, qfp_sq);
           if (status == eslOK && *abc == NULL)
             status = esl_sqfile_GuessAlphabet(*qfp_sq, &q_type);
        }
        if (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s):\n%s\n", (*qfp_sq)->filename, esl_sqfile_GetErrorBuf(*qfp_sq));
         if (q_type == eslUNKNOWN) p7_Fail("Unable to guess alphabet for the %s%s query file %s\n", (cfg->qfmt==eslUNKNOWN ? "" : esl_sqio_DecodeFormat(cfg->qfmt)), (cfg->qfmt==eslSQFILE_UNKNOWN ? "":"-formatted"), cfg->queryfile);
           *abc = esl_alphabet_Create(q_type);
      }
      if ((*abc)->type != eslAMINO) { 
        p7_Fail("Invalid alphabet type in the %s%squery file %s. Expect Amino Acid\n", (cfg->qfmt==eslUNKNOWN ? "" : esl_sqio_DecodeFormat(cfg->qfmt)), (cfg->qfmt==eslSQFILE_UNKNOWN ? "":"-formatted "), cfg->queryfile);
      }
      esl_sqfile_SetDigital(*qfp_sq, *abc);
      // read first sequence
      *qsq = esl_sq_CreateDigital(*abc);
      status = esl_sqio_Read(*qfp_sq, *qsq);
      if (status != eslOK) p7_Fail("reading sequence from file %s (%d): \n%s\n", cfg->queryfile, status, esl_sqfile_GetErrorBuf(*qfp_sq));
    }
    return status;
}

/* serial_master()
 * The serial version of bathsearch.
 * For each query HMM search the target database for hits.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().  We also use the
 * ESL_EXCEPTION and ERROR: mechanisms, but only because we know we're
 * using a fatal exception handler.
 */
static int
serial_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{

  int              i, d;
  
  /* output files */
  FILE            *ofp                      = stdout;            /* results output file (-o)                        */
  FILE            *tblfp                    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *fstblfp                  = NULL;              /* output stream for tabular per-ali (--fstblout)  */
  FILE            *hmmoutfp                 = NULL;              /* output stream for hmms (--hmmout),  only if input is an alignment file    */  
  char            *hmmfile                  = NULL;              /* file to write HMM to                            */
  int              force_single             = ( esl_opt_IsOn(go, "--singlemx") ? TRUE : FALSE );
  int              textw                    = 0;

 /* input files */
  P7_HMMFILE      *hfp                      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp                     = NULL;              /* open input sequence file                        */
  ESL_MSAFILE     *qfp_msa                  = NULL;              /* open query alifile                              */
  ESL_SQFILE      *qfp_sq                   = NULL;              /* open query seqfile                              */
  int              dbfmt                    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */

 /* query formats and HMM construction*/
  P7_HMM          *hmm                      = NULL;              /* one HMM query                                   */
  ESL_SQ          *qsq                      = NULL;              /* query sequence                                  */
  ESL_MSA         *msa                      = NULL;              /* query MSA */
  P7_BUILDER      *builder                  = NULL;
  int              msas_named               = 0;
  int              nquery                   = 0;

  /* alphabets and translation */
  float            indel_cost; 
  int              codon_table;
  ESL_ALPHABET    *abcAA                    = NULL;              /* AA  query  alphabet                                */
  ESL_ALPHABET    *abcDNA                   = NULL;              /* DNA target alphabet                              */
  ESL_GENCODE     *gcode                    = NULL;
 
 /* worker and worker items */ 
  WORKER_INFO     *info                     = NULL;
  P7_SCOREDATA    *scoredata                = NULL;              
  P7_FS_PROFILE   *gm_fs                    = NULL;
  P7_PROFILE      *gm                       = NULL;
  P7_OPROFILE     *om                       = NULL;       /* optimized query profile                  */

  /* post processing */
  int64_t          resCnt                   = 0;
  P7_TOPHITS      *tophits_accumulator      = NULL; /* to hold the top hits information from all 6 frame translations     */
  P7_PIPELINE     *pipelinehits_accumulator = NULL; /* to hold the pipeline hit information from all 6 frame translations */
  ID_LENGTH_LIST  *id_length_list           = NULL;
  ESL_STOPWATCH   *watch;
 

  /* multi threading */
  int              ncpus                    = 0; 
  int              infocnt                  = 0;
#ifdef HMMER_THREADS
  ESL_SQ_BLOCK    *block                    = NULL;
  ESL_THREADS     *threadObj                = NULL;
  ESL_WORK_QUEUE  *queue                    = NULL;
#endif

  /*error handeling */
  char             errbuf[eslERRBUFSIZE];
  int              status                   = eslOK;
  int              qhstatus                 = eslOK;
  int              sstatus                  = eslOK;
  int              ssistatus                = eslOK;  


  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");
 
  /* bathsearch accepts query files that are either hmm(s), msa(s), or sequence(s). The following 
   * code will follow the mandate of --qformat, and otherwise figure what the file type is. */

  /* (1) If we were told a specific query file type, just do what we're told */

  if (esl_sqio_IsAlignment(cfg->qfmt) /* msa file */ && !esl_opt_IsOn(go, "--qsingle_seqs") /* msa intent is not overridden */) {
      status = bath_open_msa_file(cfg, &qfp_msa, &abcAA, &msa);
      if (status != eslOK) p7_Fail("Error reading msa from the %s-formatted file %s (%d)\n", esl_sqio_DecodeFormat(cfg->qfmt), cfg->queryfile, status);
  } else if (cfg->qfmt != eslSQFILE_UNKNOWN /* sequence file */) {
      status = bath_open_seq_file(cfg, &qfp_sq, &abcAA, &qsq, esl_opt_IsOn(go, "--qsingle_seqs"));
      if (status != eslOK) p7_Fail("Error reading sequence from the %s-formatted file %s (%d)\n", esl_sqio_DecodeFormat(cfg->qfmt), cfg->queryfile, status);
  }

/* (2) Guessing query format.
 * First check if it's an HMM.  This fails easily if it's not,
 * and lets us either (a) give up right away if the input is piped (not rewindable),
 * or (b) continue guessing
 *
 * If it isn't an HMM, and it's a rewindable file, we'll check to see
 * if it's obviously an MSA file or obviously a sequence file
 * If not obvious, we'll force the user to tell us.
 * That looks like this:
 *      - Try to open it as an MSA file
 *         - if ok (i.e. it opens and passes the MSA check, including that
 *           all sequences are the same length)
 *            - if it's a FASTA format, it still might be a sequence file
 *              (note: a2m is FASTA-like, but explicitly a multiple sequence alignment)
 *                 - if the "MSA" is a single sequence, then rewind and call it
 *                   a sequence input.  Otherwise give "must specify" message
 *            - otherwise, it's an MSA;  proceed accordingly
 *         - if not ok (i.e. it's not an MSA file)
 *            - if it's anything, it must be a sequence file, proceed accordingly *
 */

 if ( cfg->qfmt == eslSQFILE_UNKNOWN ) {
      status = bath_open_hmm_file(cfg, &hfp, errbuf, &abcAA, &hmm);
      if (status != eslOK) { /* if it is eslOK, then it's an HMM, so we're done guessing */
          if (hfp!=NULL) { p7_hmmfile_Close(hfp); hfp=NULL;}
          if (strcmp(cfg->queryfile, "-") == 0 ) {
              /* we can't rewind a piped file, so we can't perform any more autodetection on the query format*/
              p7_Fail("Must specify query file format (--qformat) to read <query file> from stdin ('-')");
          } else {
              if (esl_opt_IsOn(go, "--qsingle_seqs")) { /* only try to open as a seq file*/
                  status = bath_open_seq_file(cfg, &qfp_sq, &abcAA, &qsq, esl_opt_IsOn(go, "--qsingle_seqs"));
                  if (status != eslOK) p7_Fail("Error reading query file %s (%d)\n", cfg->queryfile, status);
              } else { /* first try as an msa, then fall back to seq */
                  status = bath_open_msa_file(cfg, &qfp_msa, &abcAA, &msa);
                  if (status == eslOK) {
                      if (qfp_msa->format == eslMSAFILE_AFA) {
                          /* this could just be a sequence file with o single sequence (in which case, fall through
                           * to the "sequence" case), or with several same-sized sequences (in which case ask for guidance) */
                          if (msa->nseq > 1)
                              p7_Fail("Query file type could be either aligned or unaligned; please specify (--qformat [afa|fasta])");
                      } else {
                          /* if ok, and not fasta, then it's an MSA ... proceed */
                          cfg->qfmt = qfp_msa->format;
                      }
                  }
                  if (cfg->qfmt == eslSQFILE_UNKNOWN) { /* it's not an MSA, try seq */
                      if (qfp_msa) {
                          esl_msafile_Close(qfp_msa);
                          qfp_msa = NULL;
                          esl_msa_Destroy(msa);
                      }
                      status = bath_open_seq_file(cfg, &qfp_sq, &abcAA, &qsq, esl_opt_IsOn(go, "--qsingle_seqs"));
                      if (status != eslOK) p7_Fail("Error reading query file %s (%d)\n", cfg->queryfile, status);
                  }
              }
          }
      }  else {
          if (esl_opt_IsOn(go, "--qsingle_seqs"))
              p7_Fail("--qsingle_seqs flag is incompatible with an hmm-formatted query file\n");
      }
  }

  if (abcAA->type != eslAMINO)
     p7_Fail("Invalid alphabet type in query for bathsearch. Expect Amino Acid.\n"); 

  /* target format */
  if (esl_opt_IsOn(go, "--tformat")) {
    dbfmt = esl_sqio_EncodeFormat(esl_opt_GetString(go, "--tformat"));
    if (dbfmt == eslSQFILE_UNKNOWN) p7_Fail("%s is not a recognized sequence database file format\n", esl_opt_GetString(go, "--tformat"));
  }

  /* Open the target sequence database */
  status = esl_sqfile_Open(cfg->dbfile, dbfmt, p7_SEQDBENV, &dbfp);
  if      (status == eslENOTFOUND) p7_Fail("Failed to open sequence file %s for reading\n",          cfg->dbfile);
  else if (status == eslEFORMAT)   p7_Fail("Sequence file %s is empty or misformatted\n",            cfg->dbfile);
  else if (status == eslEINVAL)    p7_Fail("Can't autodetect format of a stdin or .gz seqfile");
  else if (status != eslOK)        p7_Fail("Unexpected error %d opening sequence file %s\n", status, cfg->dbfile);  

  /* if splicing is enabled check for SSI index*/
  if (esl_opt_IsUsed(go, "--splice")) {
    ssistatus = esl_sqfile_OpenSSI(dbfp, NULL);
    if (ssistatus != eslOK) p7_Fail("An SSI file is required for splicing. Create SSI using 'esl-sfetch --index %s' \n", cfg->dbfile);
  }


  if (esl_opt_IsUsed(go, "--restrictdb_stkey") || esl_opt_IsUsed(go, "--restrictdb_n")) {
    if (esl_opt_IsUsed(go, "--ssifile"))
      esl_sqfile_OpenSSI(dbfp, esl_opt_GetString(go, "--ssifile"));
    else
      esl_sqfile_OpenSSI(dbfp, NULL);
  }

  /* Open the results output files */
  if (esl_opt_IsOn(go, "-o"))          { if ((ofp      = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))    { if ((tblfp    = fopen(esl_opt_GetString(go, "--tblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--fstblout"))    { if ((fstblfp    = fopen(esl_opt_GetString(go, "--fstblout"),    "w")) == NULL)  esl_fatal("Failed to open tabular per-ali frameshift file %s for writing\n", esl_opt_GetString(go, "--fstblout")); }
  if (qfp_msa != NULL || qfp_sq != NULL) {
    if (esl_opt_IsOn(go, "--hmmout")) {
      hmmfile = esl_opt_GetString(go, "--hmmout");
      if ((hmmoutfp        = fopen(hmmfile,"w")) == NULL)        esl_fatal("Failed to open hmm output file %s for writing\n", hmmfile);
    }
  }


#ifdef HMMER_THREADS
  /* initialize thread data */
   ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
 
  if (ncpus > 0)
    {
      threadObj = esl_threads_Create(&pipeline_thread);
	  queue = esl_workqueue_Create(ncpus * 2);
    }
#endif

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, sizeof(*info) * infocnt);
  

   /*the query sequence will be DNA but will be translated to amino acids */
   abcDNA = esl_alphabet_Create(eslDNA); 

  /* Get translation variables */
  indel_cost = esl_opt_GetReal(go, "--fs");
  codon_table = esl_opt_GetInteger(go, "--ct");
 
  if (status == eslOK)
  {
    /* One-time initializations after alphabet <abc> becomes known */
    output_header(ofp, go, cfg->queryfile, cfg->dbfile);
    esl_sqfile_SetDigital(dbfp, abcDNA); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks

    for (i = 0; i < infocnt; ++i)
    {
      info[i].bg    = p7_bg_fs_Create(abcAA);
#ifdef HMMER_THREADS
      info[i].queue = queue;
#endif
    }

#ifdef HMMER_THREADS
    for (i = 0; i < ncpus * 2; ++i)
    {
      block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abcDNA);
      if (block == NULL)           esl_fatal("Failed to allocate sequence block");
      status = esl_workqueue_Init(queue, block);
      if (status != eslOK)          esl_fatal("Failed to add block to work queue");
    }
#endif
  }

  if (qfp_sq != NULL || qfp_msa  != NULL )  {  // need to convert query sequence / msa to HMM
    builder = p7_builder_Create(NULL, abcAA);
    if (builder == NULL)  p7_Fail("p7_builder_Create failed");
    builder->w_len      = (go != NULL && esl_opt_IsOn (go, "--w_length")) ?  esl_opt_GetInteger(go, "--w_length"): -1;
    builder->w_beta     = (go != NULL && esl_opt_IsOn (go, "--w_beta"))   ?  esl_opt_GetReal   (go, "--w_beta")    : p7_DEFAULT_WINDOW_BETA;
    if ( builder->w_beta < 0 || builder->w_beta > 1  ) esl_fatal("Invalid window-length beta value\n");  
    builder->fs         = (go != NULL)                                    ?  esl_opt_GetReal(go, "--fs"): 0.01;
  }

  if (qfp_sq != NULL || (qfp_msa != NULL && force_single )) {
    /* We'll use this scoring matrix whenever we have a single sequence (even in MSA format)
     * Default is stored in the --mx option, so it's always IsOn(). Check --mxfile first; 
     * then go to the --mx option and the default. */
    if (esl_opt_IsOn(go, "--mxfile")) status = p7_builder_SetScoreSystem (builder, esl_opt_GetString(go, "--mxfile"), NULL, esl_opt_GetReal(go, "--popen"), esl_opt_GetReal(go, "--pextend"), info->bg);
    else                              status = p7_builder_LoadScoreSystem(builder, esl_opt_GetString(go, "--mx"),           esl_opt_GetReal(go, "--popen"), esl_opt_GetReal(go, "--pextend"), info->bg);
    if (status != eslOK) p7_Fail("Failed to set single query seq score system:\n%s\n", builder->errbuf);
  }

  /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any aa   */
  gcode = esl_gencode_Create(abcDNA, abcAA);
  esl_gencode_Set(gcode, codon_table);  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m nor -M are set

  /* Outer loop: over each query HMM, alignment , or sequence in <query file>. */
  while (qhstatus == eslOK) 
  {
    gm_fs   = NULL;
    gm      = NULL;
    om      = NULL;       /* optimized query profile                  */

    if ( qfp_sq != NULL) {//  FASTA format, each query is a single sequence, they all have names
      //Turn sequence into an HMM     
      if ((qhstatus = p7_SingleBuilder(builder, qsq, info->bg, &hmm, NULL, NULL, NULL)) != eslOK) p7_Fail("build failed: %s", builder->errbuf);

    } else if ( qfp_msa != NULL ) {
      //deal with recently read MSA
      //if name isn't assigned, give it one (can only do this if there's a single unnamed alignment, so pick its filename)
      
      if (msa->name == NULL) {
        char *name = NULL;
        if (msas_named>0) p7_Fail("Name annotation is required for each alignment in a multi MSA file; failed on #%d", nquery+1);

        if (cfg->queryfile != NULL) {
          if ((status = esl_FileTail(cfg->queryfile, TRUE, &name)) != eslOK) return status; /* TRUE=nosuffix */
        } else { name = "Query"; }

        if ((status = esl_msa_SetName(msa, name, -1)) != eslOK) p7_Fail("Error assigning name to alignment");
        
        msas_named++;
        free(name);
      }

      //Turn sequence alignment into an HMM
      if (msa->nseq == 1 && force_single) {
        if (qsq!=NULL) esl_sq_Destroy(qsq);
        qsq = esl_sq_CreateDigitalFrom(msa->abc, (msa->sqname?msa->sqname[0]:"Query"), msa->ax[0], msa->alen, (msa->sqdesc?msa->sqdesc[0]:NULL), (msa->sqacc?msa->sqacc[0]:NULL), NULL);
        esl_abc_XDealign(qsq->abc, qsq->dsq,  qsq->dsq, &(qsq->n));
        if ((qhstatus = p7_SingleBuilder(builder, qsq, info->bg, &hmm, NULL, NULL, NULL)) != eslOK) p7_Fail("build failed: %s", builder->errbuf);
      } else {
        if ((qhstatus = p7_Builder(builder, msa, info->bg, &hmm, NULL, NULL, NULL, NULL)) != eslOK) p7_Fail("build failed: %s", builder->errbuf);
      }
    } else if(! (esl_opt_IsUsed(go, "--nofs"))) { //check that HMM is properly formated for bathsearch
      if( ! (hmm->evparam[p7_FTAUFS] != p7_EVPARAM_UNSET && hmm->fs && hmm->ct)) p7_Fail("HMM file %s not formated for bathsearch. Please run bathconvert.\n", cfg->queryfile);
     
      if( hmm->fs != indel_cost)  p7_Fail("Requested frameshift probability of %f does not match the frameshift probability in the HMM file %s. Please either run bathsearch with option '--fs %f' or run bathconvert with option '--fs %f'.\n", indel_cost, cfg->queryfile, hmm->fs, indel_cost);
      
      if( hmm->ct != esl_opt_GetInteger(go, "--ct"))  p7_Fail("Requested codon translation tabel ID %d does not match the codon translation tabel ID of the HMM file %s. Please either run bathsearch with option '--ct %d' or run bathconvert with option '--ct %d'.\n", codon_table, cfg->queryfile, hmm->ct, codon_table);
    } 

    if(hmm->max_length == -1)
      p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);

    if (hmmoutfp != NULL) {
      if(esl_opt_IsUsed(go, "--fs") || hmm->fs == 0.0)  hmm->fs = indel_cost;
      if(esl_opt_IsUsed(go, "--ct") || hmm->ct == 0)    hmm->ct = esl_opt_GetInteger(go, "--ct");	
      if ((status = p7_hmmfile_WriteASCII(hmmoutfp, p7_BATH_3f, hmm)) != eslOK) ESL_FAIL(status, errbuf, "HMM save failed");
    }

    nquery++;
    
    watch = esl_stopwatch_Create();
    esl_stopwatch_Start(watch);
	
    /* seqfile may need to be rewound (multiquery mode) */
    if (nquery > 1)
    {
      if (! esl_sqfile_IsRewindable(dbfp))
        esl_fatal("Target sequence file %s isn't rewindable; can't search it with multiple queries", cfg->dbfile);

      if (! esl_opt_IsUsed(go, "--restrictdb_stkey") )
        esl_sqfile_Position(dbfp, 0); //only re-set current position to 0 if we're not planning to set it in a moment
    }

    if ( cfg->firstseq_key != NULL ) { //it's tempting to want to do this once and capture the offset position for future passes, but ncbi files make this non-trivial, so this keeps it general
      sstatus = esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key);
      if (sstatus != eslOK)
        p7_Fail("Failure setting restrictdb_stkey to %d\n", cfg->firstseq_key);
    }

    if (fprintf(ofp, "Query:       %s  [M=%d]\n", hmm->name, hmm->M)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    if (hmm->acc)  { if (fprintf(ofp, "Accession:   %s\n", hmm->acc)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }
    if (hmm->desc) { if (fprintf(ofp, "Description: %s\n", hmm->desc) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

    /* Convert to an optimized model */
    gm_fs = p7_profile_fs_Create (hmm->M, abcAA);
    gm = p7_profile_Create (hmm->M, abcAA);
    om = p7_oprofile_Create(hmm->M, abcAA);
    p7_ProfileConfig(hmm, info->bg, gm, 100, p7_LOCAL); /* 100 is a dummy length for now; and MSVFilter requires local mode */
      
    p7_oprofile_Convert(gm, om);                                      /* convert <om> to <gm>*/
    p7_ProfileConfig_fs(hmm, info->bg, gcode, gm_fs, 100, p7_LOCAL);  /* build framshift aware codon HMM */
      
    /* Create processing pipeline and hit list accumulators */
    tophits_accumulator  = p7_tophits_Create(); 
    pipelinehits_accumulator = p7_pipeline_fs_Create(go, 100, 300, p7_SEARCH_SEQS);
    pipelinehits_accumulator->nmodels = 1;
    pipelinehits_accumulator->nnodes = hmm->M;

    scoredata = p7_hmm_ScoreDataCreate(om, NULL);

    for (i = 0; i < infocnt; ++i)
    {
      /* Create processing pipeline and hit list */
      info[i].gcode = gcode;
      info[i].wrk1 = esl_gencode_WorkstateCreate(go, gcode);
      info[i].wrk1->orf_block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abcAA);
      info[i].wrk2 = esl_gencode_WorkstateCreate(go, gcode);
      info[i].wrk2->orf_block = esl_sq_CreateDigitalBlock(BLOCK_SIZE, abcAA);
      info[i].th     = p7_tophits_Create();
      info[i].om     = p7_oprofile_Clone(om);
      info[i].gm     = p7_profile_Clone(gm);
      info[i].gm_fs  = p7_profile_fs_Clone(gm_fs);
      info[i].scoredata = p7_hmm_ScoreDataClone(scoredata, om->abc->Kp);
      info[i].pli = p7_pipeline_fs_Create(go, om->M, 300, p7_SEARCH_SEQS); /* L_hint = 300 is just a dummy for now */
      status = p7_pli_NewModel(info[i].pli, info[i].om, info[i].bg);
      if (status == eslEINVAL) p7_Fail(info->pli->errbuf);

      if      (strcmp(esl_opt_GetString(go, "--strand"), "both")  == 0) info[i].pli->strands = p7_STRAND_BOTH; 
      else if (strcmp(esl_opt_GetString(go, "--strand"), "plus")  == 0) info[i].pli->strands = p7_STRAND_TOPONLY;
      else if (strcmp(esl_opt_GetString(go, "--strand"), "minus") == 0) info[i].pli->strands = p7_STRAND_BOTTOMONLY;
      else     p7_Fail("Unrecognized argument for --strand ('%s'). Only 'both', 'plus', and 'minus' allowed.", esl_opt_GetString(go, "--strand"));

      if (  esl_opt_IsUsed(go, "--block_length") )
        info[i].pli->block_length = esl_opt_GetInteger(go, "--block_length");
      else
        info[i].pli->block_length = BATH_MAX_RESIDUE_COUNT;

#ifdef HMMER_THREADS
      if (ncpus > 0) esl_threads_AddThread(threadObj, &info[i]);
#endif
    }

    /* establish the id_lengths data structutre */
    id_length_list = init_id_length(1000);

#ifdef HMMER_THREADS
    if (ncpus > 0)  sstatus = thread_loop(info, id_length_list,threadObj, queue, dbfp, cfg->firstseq_key, cfg->n_targetseq);
    else
#endif
      sstatus = serial_loop(info, id_length_list, dbfp, cfg->firstseq_key, cfg->n_targetseq);

    switch(sstatus) {
      case eslEFORMAT:
        esl_fatal("Parse failed (sequence file %s):\n%s\n",
                   dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
        break;
      case eslEOF:
      case eslOK:
        /* do nothing */
         break;
      default:
        esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);
    }
      
    //need to re-compute e-values before merging (when list will be sorted)
    if (esl_opt_IsUsed(go, "-Z")) {
      resCnt = 1000000*esl_opt_GetReal(go, "-Z");
      if ( info[0].pli->strands == p7_STRAND_BOTH)
        resCnt *= 2;
    }
    else
    {
      for (i = 0; i < infocnt; ++i){
        resCnt += info[i].pli->nres;
      }
    }

    for (i = 0; i < infocnt; ++i)
      p7_tophits_ComputeBATHEvalues(info[i].th, resCnt, info[i].om->max_length);

    /* merge the results of the search results */
    for (i = 0; i < infocnt; ++i)
    {
      p7_tophits_Merge(tophits_accumulator, info[i].th);
      p7_pipeline_Merge(pipelinehits_accumulator, info[i].pli);
      p7_pipeline_fs_Destroy(info[i].pli);
      p7_tophits_Destroy(info[i].th);
      p7_oprofile_Destroy(info[i].om);
      p7_profile_Destroy(info[i].gm);
      p7_profile_fs_Destroy(info[i].gm_fs);
      p7_hmm_ScoreDataDestroy(info[i].scoredata);

      if(info[i].wrk1->orf_block != NULL)
      {
        esl_sq_DestroyBlock(info[i].wrk1->orf_block);
        info[i].wrk1->orf_block = NULL;
        esl_gencode_WorkstateDestroy(info[i].wrk1);
      }
      if(info[i].wrk2->orf_block != NULL)
      {
        esl_sq_DestroyBlock(info[i].wrk2->orf_block);
        info[i].wrk2->orf_block = NULL;
        esl_gencode_WorkstateDestroy(info[i].wrk2);
      }
    }

    /* Sort and remove duplicates */
    p7_tophits_SortBySeqidxAndAlipos(tophits_accumulator);
    assign_Lengths(tophits_accumulator, id_length_list);
    p7_tophits_RemoveDuplicates(tophits_accumulator, pipelinehits_accumulator->use_bit_cutoffs);

    /* Sort and remove hits bellow threshold */
    p7_tophits_SortBySortkey(tophits_accumulator);

    /* Set Z = 1 to prevent changing e-values. Correct Z 
     * was calcualted by p7_tophits_ComputeBathEvalues() */
    pipelinehits_accumulator->Z = 1;    
    p7_tophits_Threshold(tophits_accumulator, pipelinehits_accumulator);

    if (esl_opt_IsUsed(go, "--splice") && tophits_accumulator->N) {
      p7_splice_SpliceHits(tophits_accumulator, om, gm, gm_fs, go, gcode, dbfp, ofp);
      //SpliceHits(tophits_accumulator,dbfp,gm,om,gcode,go,ofp,textw);
    }
    
    /* Print the results.  */
    pipelinehits_accumulator->n_output = pipelinehits_accumulator->pos_output = 0; 
    for (i = 0; i < tophits_accumulator->N; i++) {
      if ( (tophits_accumulator->hit[i]->flags & p7_IS_REPORTED) || tophits_accumulator->hit[i]->flags & p7_IS_INCLUDED) {
        pipelinehits_accumulator->n_output++;
          
	  for(d = 0; d < tophits_accumulator->hit[i]->ndom; d++)
        pipelinehits_accumulator->pos_output += 1 + (tophits_accumulator->hit[i]->dcl[d].jali > tophits_accumulator->hit[i]->dcl[d].iali ? tophits_accumulator->hit[i]->dcl[d].jali - tophits_accumulator->hit[i]->dcl[d].iali : tophits_accumulator->hit[i]->dcl[d].iali - tophits_accumulator->hit[i]->dcl[d].jali) ;
      }
    }

    p7_tophits_Targets(ofp, tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    p7_tophits_Domains(ofp, tophits_accumulator, pipelinehits_accumulator, textw); if (fprintf(ofp, "\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

    if (tblfp)     p7_tophits_TabularTargets(tblfp,    hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));
    if (fstblfp)   p7_tophits_TabularFrameshifts(fstblfp,    hmm->name, hmm->acc, tophits_accumulator, pipelinehits_accumulator, (nquery == 1));

    esl_stopwatch_Stop(watch);
    p7_pli_Statistics(ofp, pipelinehits_accumulator, watch);
    if (fprintf(ofp, "//\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

    p7_pipeline_fs_Destroy(pipelinehits_accumulator);
    p7_tophits_Destroy(tophits_accumulator);
    p7_oprofile_Destroy(om);
    p7_profile_Destroy(gm);
    p7_profile_fs_Destroy(gm_fs);
    p7_hmm_Destroy(hmm);
    p7_hmm_ScoreDataDestroy(scoredata);
    destroy_id_length(id_length_list);
    if (qsq != NULL) esl_sq_Reuse(qsq);
      
    if (hfp != NULL) {
      qhstatus = p7_hmmfile_Read(hfp, &abcAA, &hmm);
    } else if (qfp_msa != NULL){
      esl_msa_Destroy(msa);
      qhstatus = esl_msafile_Read(qfp_msa, &msa);
    } else { // qfp_sq
      qhstatus = esl_sqio_Read(qfp_sq, qsq);
    } 
    
   if (qhstatus != eslOK && qhstatus != eslEOF) p7_Fail("reading from query file %s (%d)\n", cfg->queryfile, qhstatus);
  } /* end outer loop over query HMMs */

  if (hfp != NULL) {
    switch(qhstatus) {
      case eslEOD:       p7_Fail("read failed, HMM file %s may be truncated?", cfg->queryfile);      break;
      case eslEFORMAT:   p7_Fail("bad file format in HMM file %s",             cfg->queryfile);      break;
      case eslEINCOMPAT: p7_Fail("HMM file %s contains different alphabets",   cfg->queryfile);      break;
      case eslEOF:       /* do nothing. EOF is what we want. */                                    break;
      default:           p7_Fail("Unexpected error (%d) in reading HMMs from %s", qhstatus, cfg->queryfile);
    }
  } else if (qfp_msa != NULL){
    if (qhstatus != eslEOF ) esl_msafile_ReadFailure(qfp_msa, status);
  } else { // qfp_sq
    if      (qhstatus == eslEFORMAT) p7_Fail("Parse failed (sequence file %s):\n%s\n",
                qfp_sq->filename, esl_sqfile_GetErrorBuf(qfp_sq));
    else if (qhstatus != eslEOF)     p7_Fail("Unexpected error %d reading sequence file %s",
                qhstatus, qfp_sq->filename);
  }
  
  if (hmmoutfp != NULL)
    fclose(hmmoutfp);

  /* Terminate outputs... any last words? */
  if (tblfp)    p7_tophits_TabularTail(tblfp,    "bathsearch", p7_SEARCH_SEQS, cfg->queryfile, cfg->dbfile, go);
  if (fstblfp)  p7_tophits_TabularTail(fstblfp,  "bathsearch", p7_SEARCH_SEQS, cfg->queryfile, cfg->dbfile, go); 
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit */
  for (i = 0; i < infocnt; ++i)
    p7_bg_Destroy(info[i].bg);

#ifdef HMMER_THREADS
  if (ncpus > 0)
  {
    esl_workqueue_Reset(queue);
    while (esl_workqueue_Remove(queue, (void **) &block) == eslOK) { 
      esl_sq_DestroyBlock(block);
    }
    esl_workqueue_Destroy(queue);
    esl_threads_Destroy(threadObj);
  }
#endif
  
  free(info);

  if (hfp) p7_hmmfile_Close(hfp);
  if (qfp_msa) esl_msafile_Close(qfp_msa);
  if (qfp_sq)  esl_sqfile_Close(qfp_sq);
  
  if (builder) p7_builder_Destroy(builder);
  if (qsq)     esl_sq_Destroy(qsq);
  
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_gencode_Destroy(gcode);
  esl_stopwatch_Destroy(watch);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (fstblfp)         fclose(fstblfp);

  return eslOK;

ERROR:

  if (hfp) p7_hmmfile_Close(hfp);
  if (qfp_msa) esl_msafile_Close(qfp_msa);
  if (qfp_sq)  esl_sqfile_Close(qfp_sq);

  if (builder) p7_builder_Destroy(builder);
  if (qsq)     esl_sq_Destroy(qsq);

  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (fstblfp)       fclose(fstblfp);

  if (hmmfile != NULL) free (hmmfile);
  return eslFAIL;
}

static int
serial_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs)
{
  int  sstatus = eslOK;
  int seq_id = 0;

  ESL_ALPHABET *abcDNA = esl_alphabet_Create(eslDNA);
  ESL_SQ       *dbsq_dna    = esl_sq_CreateDigital(abcDNA);   /* (digital) nucleotide sequence, to be translated into ORFs  */
  sstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq_dna);

  while (sstatus == eslOK && (n_targetseqs==-1 || seq_id < n_targetseqs) ) 
  {
    dbsq_dna->idx = seq_id;
    if (dbsq_dna->n < 15) continue; /* do not process sequence of less than 5 codons */

    dbsq_dna->L = dbsq_dna->n; /* here, L is not the full length of the sequence in the db, just of the currently-active window;  required for esl_gencode machinations */
    
    if (info->pli->strands != p7_STRAND_BOTTOMONLY) 
    {
      info->pli->nres += dbsq_dna->n;
   
       /* translate DNA sequence to 3 frame ORFs */
      do_sq_by_sequences(info->gcode, info->wrk1, dbsq_dna);

      p7_Pipeline_BATH(info->pli, info->om, info->gm, info->gm_fs, info->scoredata, info->bg, info->th, info->pli->nseqs, dbsq_dna, info->wrk1->orf_block, info->wrk2, info->gcode, p7_NOCOMPLEMENT);
      p7_pipeline_fs_Reuse(info->pli); // prepare for next search

      esl_sq_ReuseBlock(info->wrk1->orf_block);    
    } 

    if (info->pli->strands != p7_STRAND_TOPONLY) 
    {   
      info->pli->nres += dbsq_dna->n;
  
      /* Reverse complement and translate DNA sequence to 3 frame ORFs */
      esl_sq_ReverseComplement(dbsq_dna);
      do_sq_by_sequences(info->gcode, info->wrk1, dbsq_dna);
	
      p7_Pipeline_BATH(info->pli, info->om, info->gm, info->gm_fs, info->scoredata, info->bg, info->th, info->pli->nseqs, dbsq_dna, info->wrk1->orf_block, info->wrk2, info->gcode, p7_COMPLEMENT); 
      p7_pipeline_fs_Reuse(info->pli); // prepare for next search
      
      esl_sq_ReuseBlock(info->wrk1->orf_block);
      
      /* Reverse sequence back to original */
      esl_sq_ReverseComplement(dbsq_dna);
    } 

    sstatus = esl_sqio_ReadWindow(dbfp, info->om->max_length, info->pli->block_length, dbsq_dna);
    
    if (sstatus == eslEOD) 
    { 
      /* no more left of this sequence ... move along to the next sequence. */
      add_id_length(id_length_list, dbsq_dna->idx, dbsq_dna->L);
      info->pli->nseqs++;
      esl_sq_Reuse(dbsq_dna);
      sstatus = esl_sqio_ReadWindow(dbfp, 0, info->pli->block_length, dbsq_dna);
      seq_id++;
    }
  }

  if(abcDNA) esl_alphabet_Destroy(abcDNA);
  if(dbsq_dna) esl_sq_Destroy(dbsq_dna);
  
  return sstatus;
}

#ifdef HMMER_THREADS
static int
thread_loop(WORKER_INFO *info, ID_LENGTH_LIST *id_length_list, ESL_THREADS *obj, ESL_WORK_QUEUE *queue, ESL_SQFILE *dbfp, char *firstseq_key, int n_targetseqs)
{
  int i;
  int           status   = eslOK;
  int           sstatus  = eslOK;
  int           eofCount = 0;
  int           seqid    = -1;
  int           abort    = FALSE; // in the case n_targetseqs != -1, a block may get abbreviated
  ESL_SQ_BLOCK *block;
  ESL_SQ       *tmpsq;
  void         *newBlock;
  
  tmpsq = esl_sq_CreateDigital(dbfp->abc);

  esl_workqueue_Reset(queue);
  esl_threads_WaitForStart(obj);

  status = esl_workqueue_ReaderUpdate(queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue reader failed");
  ((ESL_SQ_BLOCK *)newBlock)->complete = TRUE;

  /* Main loop: */
  while (sstatus == eslOK)
  {
    block = (ESL_SQ_BLOCK *) newBlock;     
	
    if (abort) {
      block->count = 0;
      sstatus = eslEOF;
    } else {
      sstatus = esl_sqio_ReadBlock(dbfp, block, info->pli->block_length, n_targetseqs, FALSE, TRUE);
    }

    block->first_seqidx = info->pli->nseqs;
    seqid = block->first_seqidx;

    for (i=0; i<block->count; i++) {
      block->list[i].idx = seqid;
      add_id_length(id_length_list, seqid, block->list[i].L); // NOLINT(cppcoreguidelines-narrowing-conversions)
      seqid++;

      if (       seqid == n_targetseqs // hit the sequence target
           && ( i<block->count-1 ||  block->complete ) // and either it's not the last sequence (so it's complete), or its complete
         ) 
      {
        abort = TRUE;
        block->count = i+1;
        break;
      }
    } 

    info->pli->nseqs += block->count  - ((abort || block->complete) ? 0 : 1);// if there's an incomplete sequence read into the block wait to count it until it's complete.	

    if (sstatus == eslEOF) {
      if (eofCount < esl_threads_GetWorkerCount(obj)) sstatus = eslOK;
      ++eofCount;
    } else if (!block->complete ) {
      /* The final sequence on the block was an incomplete window of the 
       * active sequence, so our next read will need a copy of it to 
       * correctly deal with overlapping regions. We capture a copy of the 
       * sequence here before sending it off to the pipeline to avoid odd 
       * race conditions that can occur otherwise. Copying the entire sequence 
       * isn't really necessary, and is a bit heavy-handed. Could accelerate 
       * if this proves to have any notable impact on speed. */
      esl_sq_Copy(block->list + (block->count - 1) , tmpsq);
    }

    if (sstatus == eslOK) {
      status = esl_workqueue_ReaderUpdate(queue, block, &newBlock);
      if (status != eslOK) esl_fatal("Work queue reader failed");

      /*newBlock needs all this information so the next ReadBlock call will know what to do */
      ((ESL_SQ_BLOCK *)newBlock)->complete = block->complete;
      if (!block->complete) {
        /* Push the captured copy of the previously-read sequence into the new block,
         * in preparation for ReadWindow  (double copy ... slower than necessary) */
        esl_sq_Copy(tmpsq, ((ESL_SQ_BLOCK *)newBlock)->list);

        if (  ((ESL_SQ_BLOCK *)newBlock)->list->n < info->om->max_length ) {
          /*no reason to search the final partial sequence on the block, as 
           * the next block will search this whole chunk */
          ((ESL_SQ_BLOCK *)newBlock)->list->C = ((ESL_SQ_BLOCK *)newBlock)->list->n;
          (((ESL_SQ_BLOCK *)newBlock)->count)--;
        } else {
          ((ESL_SQ_BLOCK *)newBlock)->list->C = info->om->max_length;
        }
      }
    }
  }

  status = esl_workqueue_ReaderUpdate(queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue reader failed");

  if (sstatus == eslEOF) {
    /* wait for all the threads to complete */
    esl_threads_WaitForFinish(obj);
    esl_workqueue_Complete(queue);
  }

  esl_sq_Destroy(tmpsq);

  return sstatus;
}

static void 
pipeline_thread(void *arg)
{
  int i;
  int status;
  int workeridx;
  WORKER_INFO   *info;
  ESL_THREADS   *obj;
  ESL_SQ_BLOCK  *block = NULL;
  void          *newBlock;

  impl_Init();
  obj = (ESL_THREADS *) arg;
  esl_threads_Started(obj, &workeridx);

  info = (WORKER_INFO *) esl_threads_GetData(obj, workeridx);

  status = esl_workqueue_WorkerUpdate(info->queue, NULL, &newBlock);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  /* loop until all blocks have been processed */
  block = (ESL_SQ_BLOCK *) newBlock;
 
  while (block->count > 0)
  {
    /* Main loop: */
    for (i = 0; i < block->count; ++i)
    {
      ESL_SQ *dnaSeq = block->list + i;
      dnaSeq->L = dnaSeq->n; /* here, L is not the full length of the sequence in the db, just of the currently-active window;  required for esl_gencode machinations */
     

      if (info->pli->strands != p7_STRAND_BOTTOMONLY) {

        info->pli->nres += dnaSeq->n;
        do_sq_by_sequences(info->gcode, info->wrk1, dnaSeq);
       
        p7_Pipeline_BATH(info->pli, info->om, info->gm, info->gm_fs, info->scoredata, info->bg, info->th, block->first_seqidx + i, dnaSeq, info->wrk1->orf_block, info->wrk2, info->gcode, p7_NOCOMPLEMENT);

        p7_pipeline_fs_Reuse(info->pli); // prepare for next search

        esl_sq_ReuseBlock(info->wrk1->orf_block);
      } 

      if (info->pli->strands != p7_STRAND_TOPONLY) {
        info->pli->nres += dnaSeq->n;
        esl_sq_ReverseComplement(dnaSeq);
        do_sq_by_sequences(info->gcode, info->wrk1, dnaSeq);
	
        p7_Pipeline_BATH(info->pli, info->om, info->gm, info->gm_fs, info->scoredata, info->bg, info->th, block->first_seqidx + i, dnaSeq, info->wrk1->orf_block, info->wrk2, info->gcode, p7_COMPLEMENT);

        p7_pipeline_fs_Reuse(info->pli); // prepare for next search

	esl_sq_ReuseBlock(info->wrk1->orf_block);
        esl_sq_ReverseComplement(dnaSeq);
      }
    }  
    status = esl_workqueue_WorkerUpdate(info->queue, block, &newBlock);
    if (status != eslOK) esl_fatal("Work queue worker failed");
   
    /* loop until all blocks have been processed */
    block = (ESL_SQ_BLOCK *) newBlock; 
  } 
  
  status = esl_workqueue_WorkerUpdate(info->queue, block, NULL);
  if (status != eslOK) esl_fatal("Work queue worker failed");

  esl_threads_Finished(obj, workeridx);
}
#endif   /* HMMER_THREADS */
 

static ID_LENGTH_LIST *
init_id_length( int size )
{
  int status;
  ID_LENGTH_LIST *list;

  ESL_ALLOC (list, sizeof(ID_LENGTH_LIST));
  list->count = 0;
  list->size  = size;
  list->id_lengths = NULL;

  ESL_ALLOC (list->id_lengths, size * sizeof(ID_LENGTH));

  return list;

ERROR:
  return NULL;
}

static void
destroy_id_length( ID_LENGTH_LIST *list )
{

  if (list != NULL) {
    if (list->id_lengths != NULL) free (list->id_lengths);
    free (list);
  }

}

static int
add_id_length(ID_LENGTH_LIST *list, int id, int L)
{
  int status;

  if (list->count > 0 && list->id_lengths[list->count-1].id == id) {
    /* the last time this gets updated, it'll have the sequence's actual length */
    list->id_lengths[list->count-1].length = L;
  } else {
    if (list->count == list->size) {
      list->size *= 10;
      ESL_REALLOC(list->id_lengths, list->size * sizeof(ID_LENGTH));
    }

    list->id_lengths[list->count].id     = id;
    list->id_lengths[list->count].length = L;

    list->count++;
  }
  return eslOK;

ERROR:
  return status;
}


static int
assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list) 
{

  int i;
  int j = 0;

  for (i=0; i<th->N; i++) {
    while (th->hit[i]->seqidx != id_length_list->id_lengths[j].id) { j++; } 
    th->hit[i]->dcl[0].ad->L = id_length_list->id_lengths[j].length;
  }

  return eslOK;
}


