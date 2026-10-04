/* bathsearch: search protein profile HMM(s) against a DNA sequence database.
 * The database is read and translated once, for all queries. */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

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
#endif /*HMMER_THREADS*/

#include "hmmer.h"
#include "p7_splice.h"

/* set the max residue count to 1/4 meg when reading a block */
#define BATH_MAX_RESIDUE_COUNT (1024 * 256)  /* 1/4 Mb */


typedef struct {
  P7_BG            *bg;	        /* null model                                                        */
  ESL_SQ           *ntsq;       /* DNA target sequence                                               */
  P7_PIPELINE      *pli;        /* work pipeline                                                     */
  P7_TOPHITS       *th;         /* top hit results                                                   */
  P7_OPROFILE      *om;         /* optimized query profile                                           */
  P7_PROFILE       *gm;         /* non optimized query profile                                       */
  P7_FS_PROFILE    *gm_fs5;     /* non optimized 5 codon length frameshift query profile             */
  P7_FS_OPROFILE   *om_fs5;     /* optimized 5 codon length frameshift query profile             */
  P7_FS_OPROFILE   *om_fs3;     /* optimized 3 codon length frameshift query profile             */
  P7_SCOREDATA     *scoredata;  /* used to create DNA windows from ORFs                              */
  ESL_GENCODE      *gcode;      /* used for translating ORFs                                         */
  ESL_GENCODE_WORKSTATE *wrk;   /* used for translation of taget DNA to ORFs                         */ 
  P7_HMM_WINDOWLIST     *hw;    /* exon seeds for splicing algorithms                                */

} WORKER_INFO;



static ID_LENGTH_LIST* init_id_length( int size );
static void            destroy_id_length( ID_LENGTH_LIST *list );
static int             add_id_length(ID_LENGTH_LIST *list, int id, int64_t L);
static int             assign_Lengths(P7_TOPHITS *th, ID_LENGTH_LIST *id_length_list);

#define REPOPTS     "-E,-T"//--cut_ga,--cut_nc,--cut_tc"
#define DOMREPOPTS  "--domE,--domT,--cut_ga,--cut_nc,--cut_tc"
#define INCOPTS     "--incE,--incT"//--cut_ga,--cut_nc,--cut_tc"
#define INCDOMOPTS  "--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"
#define THRESHOPTS  "-E,-T,--domE,--domT,--incE,--incT,,--incdomE,--incdomT,--cut_ga,--cut_nc,--cut_tc"

#define CPUOPTS     NULL
#define MPIOPTS     NULL

static ESL_OPTIONS options[] = {
  /* name             type            default    env          range      toggles reqs  incomp          help                                                                        docgroup*/
  { "-h",             eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "show brief help on version and usage",                                     1 },

  /* Algorithm options */
  { "--fs",           eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, "--splice",     "use frameshift alignment algorthims",                                      2 },
  { "--splice",       eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, "--fs",         "use spliced alignment algorithms ",                                        2 },

  /* Control of output */
  { "-o",             eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "direct output to file <f>, not stdout",                                    3 },
  { "--tblout",       eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "save parseable table of hits to file <f>",                                 3 },
  { "--exontblout",   eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,"--splice",NULL,         "save parseable table of exons to file <f>",                                3 },
  { "--fstblout",     eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL, "--fs", NULL,           "save table of frameshift locations to file <f>",                           3 },
  { "--hmmout",       eslARG_OUTFILE, NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "if input is alignment(s) or sequence(s) write produced hmms to file <f>",  3 },
  { "--acc",          eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "prefer accessions over names in output",                                   3 },
  { "--noali",        eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "don't output alignments, so output is smaller",                            3 },
  { "--notrans",      eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "don't show the translated DNA sequence in  alignment",                     3 }, 
  { "--frameline",    eslARG_NONE,    FALSE,     NULL,        NULL,      NULL, "--fs", NULL,           "include frame of each codon in  alignment",                                3 },
  { "--cigar",        eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,"--tblout", NULL,        "include alignment CIGAR string in table output (with --tblout)",           3 },
  { "--notextw",      eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL,"--textw",       "unlimit ASCII text output line width",                                     3 },
  { "--textw",        eslARG_INT,    "150",      NULL,       "n>=120",   NULL,   NULL,"--notextw",     "set max width of ASCII text output lines",                                 3 },

   /* Translation options */ 
  { "--ct",           eslARG_INT,    "1",        NULL,        NULL,      NULL,   NULL, NULL,           "use alt genetic code of NCBI translation table (see end of help)",         4 },
  { "-l",             eslARG_INT,    "20",       NULL,        NULL,      NULL,   NULL, NULL,           "minimum ORF length",                                                       4 },
  { "-m",             eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"-M",            "ORFs must initiate with AUG only",                                         4 },
  { "-M",             eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"-m",            "ORFs must start with allowed initiation codon",                            4 },
  { "--strand",       eslARG_STRING, "both",     NULL,        NULL,      NULL,   NULL, NULL,           "translate only forward strand 'plus' or reverse complement strand 'minus'",4 },

  /* Splicing options */
  { "--min_intron",   eslARG_INT,    "13",       NULL,   "13<=n<=50",    NULL,"--splice", NULL,         "minimum intron length",                                                    5 },
  { "--max_intron",   eslARG_INT,    "200000",   NULL, "10000<=n<=125000000",  NULL,"--splice", NULL,   "maximum intron length",                                                    5 },  

  /* Control of reporting and inclusion thresholds */
  { "-E",             eslARG_REAL,   "10.0",     NULL,       "x>0",      NULL,   NULL, REPOPTS,        "report sequences <= this E-value threshold in output",                     6 },
  { "-T",             eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,   NULL, REPOPTS,        "report sequences >= this score threshold in output",                       6 },
  { "--incE",         eslARG_REAL,   "0.01",     NULL,       "x>0",      NULL,   NULL, INCOPTS,        "consider sequences <= this E-value threshold as significant",              6 },
  { "--incT",         eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,   NULL, INCOPTS,        "consider sequences >= this score threshold as significant",                6 },

  /* Control of acceleration pipeline */
  { "--max",          eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL,"--F1,--F2,--F3,--F4","turn all heuristic filters off (less speed, more power)",             7 },
  { "--F1",           eslARG_REAL,   "0.02",     NULL,        NULL,      NULL,   NULL,"--max",         "stage 1 (MSV) threshold: promote hits w/ P <= F1",                         7 },
  { "--F2",           eslARG_REAL,   "1e-3",     NULL,        NULL,      NULL,   NULL,"--max",         "stage 2 (Vit) threshold: promote hits w/ P <= F2",                         7 },
  { "--F3",           eslARG_REAL,   "1e-5",     NULL,        NULL,      NULL,   NULL,"--max",         "stage 3 (Fwd) threshold: promote hits w/ P <= F3",                         7 },
  { "--F4",           eslARG_REAL,   "5e-4",     NULL,        NULL,      NULL,  "--fs","--max",         "stage 4 (FS-Fwd) threshold: promote hits w/ P <= F4",                     7 },
  { "--nobias",       eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL,"--max",         "turn off composition bias filter",                                         7 },
  { "--nonull2",      eslARG_NONE,    NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "turn off biased composition score corrections",                            7 },

  /* input formats */
  { "--qformat",      eslARG_STRING,  NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "assert query is in format <s> (can be seq or msa format)",                 8 },
  { "--tformat",      eslARG_STRING,  NULL,      NULL,        NULL,      NULL,   NULL, NULL,           "assert target <seqfile> is in format <s>: no autodetection",               8 },

  /* Control of scoring system */
  { "--singlemx",     eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, NULL,           "use substitution score matrix w/ single-sequence MSA-format inputs",       9 },
  { "--popen",        eslARG_REAL,   "0.02",     NULL,       "0<=x<0.5", NULL,   NULL, NULL,           "gap open probability",                                                     9 },
  { "--pextend",      eslARG_REAL,   "0.4",      NULL,       "0<=x<1",   NULL,   NULL, NULL,           "gap extend probability",                                                   9 },
  { "--mx",           eslARG_STRING, "BLOSUM62", NULL,        NULL,      NULL,   NULL,"--mxfile",      "substitution score matrix choice (of some built-in matrices)",             9 },
  { "--mxfile",       eslARG_INFILE,  NULL,      NULL,        NULL,      NULL,   NULL,"--mx",          "read substitution score matrix from file <f>",                             9 },

/* Other options */
  { "-Z",             eslARG_REAL,    FALSE,     NULL,       "x>=0",     NULL,   NULL, NULL,           "set database size (Megabases) to <x> for E-value calculations",           10 }, 
  { "--seed",         eslARG_INT,    "42",       NULL,       "n>=0",     NULL,   NULL, NULL,           "set RNG seed to <n> (if 0: one-time arbitrary seed)",                     10 },
  { "--w_beta",       eslARG_REAL,    NULL,      NULL,       "0>=x<=1",  NULL,   NULL, NULL,           "tail mass at which window length is determined",                          10 },
  { "--w_length",     eslARG_INT,     NULL,      NULL,       "x>=4",      NULL,   NULL, NULL,           "window length - essentially max expected hit length" ,                   10 },
  #ifdef HMMER_THREADS 
  { "--block_length", eslARG_INT,     NULL,      NULL,       "n>=50000", NULL,   NULL, NULL,           "length of blocks read from target database (threaded) ",                  10 },
  { "--cpu",          eslARG_INT,     p7_NCPU,  "HMMER_NCPU","n>=0",     NULL,   NULL, CPUOPTS,        "number of parallel CPU workers to use for multithreads",                  10 },
#endif
 
  /* Restrict search to subset of database - hidden because these flags are
   *   (a) currently for internal use
   *   (b) probably going to change
   */
  { "--restrictdb_stkey", eslARG_STRING,"0",     NULL,        NULL,      NULL,   NULL, NULL,           "Search starts at the sequence with name <s> ",                             99 },
  { "--restrictdb_n",     eslARG_INT,   "-1",    NULL,        NULL,      NULL,   NULL, NULL,           "Search <j> target sequences (starting at --restrictdb_stkey)",             99 },
  { "--ssifile",          eslARG_STRING, NULL,   NULL,        NULL,      NULL,   NULL, NULL,           "restrictdb_x values require ssi file. Override default to <s>",            99 },

  /* Not used, but retained because esl option-handling code errors if it isn't kept here.  Placed in group 99 so that it doesn't print to help*/
  { "--domZ",         eslARG_REAL,    FALSE,     NULL,       "x>0",      NULL,  NULL, NULL,            "Not used",                                                                 99 },
  { "--domE",         eslARG_REAL,   "10.0",     NULL,       "x>0",      NULL,  NULL, DOMREPOPTS,      "Not used",                                                                 99 },
  { "--domT",         eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,  NULL, DOMREPOPTS,      "Not used",                                                                 99 },
  { "--incdomE",      eslARG_REAL,   "0.01",     NULL,       "x>0",      NULL,  NULL, INCDOMOPTS,      "Not used",                                                                 99 },
  { "--incdomT",      eslARG_REAL,    FALSE,     NULL,        NULL,      NULL,  NULL, INCDOMOPTS,      "Not used",                                                                 99 },
  { "--crick",        eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,  NULL, NULL,            "only translate top strand",                                                99 },
  { "--watson",       eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,  NULL, NULL,            "only translate bottom strand",                                             99 }, 
  /* Hidden frameshift options - for debugging/testing */
  { "--fsonly",       eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,   NULL, "--splice",     "send all potential hits to the frameshift aware pipeline",                 99 },
  
  /* Hidden splicing options - for debugging/testing */
  { "--nodeinfo",     eslARG_NONE,    FALSE,     NULL,        NULL,      NULL,"--exontblout", NULL,    "additional info on node types for --exontblout",                           99 },
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

static int  search_master(ESL_GETOPTS *go, struct cfg_s *cfg);

#define BLOCK_SIZE 1000

static int  scan_search(ESL_GETOPTS *go, struct cfg_s *cfg, P7_HMMFILE *hfp, P7_HMM *hmm, ESL_ALPHABET **p_abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode,
                        int ncpus, ESL_SQFILE *dbfp, FILE *ofp, FILE *tblfp, FILE *exontblfp, FILE *fstblfp, int textw, ESL_STOPWATCH *watch);



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
      p7_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);
      if (puts("\nBasic options:")                                          < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 100); /* 1= group; 2 = indentation; 100=textwidth*/

      if (puts("\nAlgorithm options:")                                      < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 100);

      if (puts("\nOptions directing output:")                               < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 100); 

      if (puts("\nOptions controlling translation:")                        < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 100); 

      if (puts("\nOptions controlling splicing:")                           < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 100);

      if (puts("\nOptions controlling reporting and inclusion thresholds:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 100); 

      if (puts("\nOptions controlling acceleration heuristics:")            < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 7, 2, 100); 

      if (puts("\nOptions setting input formats:")                          < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 100);

      if (puts("\nOptions handling single sequence inputs:")                < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 100);

      if (puts("\nOther expert options:")                                   < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 100); 
      
      if (puts("\nAvailable NCBI genetic code tables (for --ct <id>):")     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
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
  if (                                                         fprintf(ofp, "# codon translation table:                       %d\n",      esl_opt_GetInteger(go, "--ct"))              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")                              && fprintf(ofp, "# output directed to file:                       %s\n",      esl_opt_GetString(go, "-o"))                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--tblout")                        && fprintf(ofp, "# per-seq hits tabular output:                   %s\n",      esl_opt_GetString(go, "--tblout"))           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--exontblout")                    && fprintf(ofp, "# per-seq exon tabular output:                   %s\n",      esl_opt_GetString(go, "--exontblout"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
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
  if (esl_opt_IsUsed(go, "--F4")                            && fprintf(ofp, "# ORF P threshold for FS FWD:                 <= %g\n",      esl_opt_GetReal(go, "--F4"))                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nobias")                        && fprintf(ofp, "# biased composition HMM filter:                 off\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--nonull2")                       && fprintf(ofp, "# null2 bias corrections:                        off\n")                                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fs")                            && fprintf(ofp, "# Use the frameshift aware algorithms\n")                                                                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fsonly")                        && fprintf(ofp, "# Use only the frameshift aware pipeline\n")                                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--restrictdb_stkey")              && fprintf(ofp, "# Restrict db to start at seq key:               %s\n",      esl_opt_GetString(go, "--restrictdb_stkey")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--restrictdb_n")                  && fprintf(ofp, "# Restrict db to # target seqs:                  %d\n",      esl_opt_GetInteger(go, "--restrictdb_n"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--ssifile")                       && fprintf(ofp, "# Override ssi file to:                          %s\n",      esl_opt_GetString(go, "--ssifile"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "-Z")                              && fprintf(ofp, "# database size is set to:                       %.1f Mb\n", esl_opt_GetReal(go, "-Z"))                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); 
  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0               && fprintf(ofp, "# random number seed:                            one-time arbitrary\n")                                   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if (                                                  fprintf(ofp, "# random number seed set to:                     %d\n",      esl_opt_GetInteger(go, "--seed"))            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--qformat")                       && fprintf(ofp, "# query format asserted:                         %s\n",      esl_opt_GetString(go, "--qformat"))          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
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

  status = search_master(go, &cfg);

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
bath_open_msa_file(struct cfg_s *cfg,  ESL_MSAFILE **qfp_msa, ESL_ALPHABET **abc) {

  
  ESL_MSA *msa = NULL;

  int status = esl_msafile_Open(abc, cfg->queryfile, NULL, cfg->qfmt, NULL, qfp_msa);
  if (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open query file %s.\n", cfg->queryfile);

  if (status == eslOK) {
    status = esl_msafile_Read(*qfp_msa, &msa);
  }

  if (status == eslOK && (*qfp_msa)->format == eslMSAFILE_AFA && cfg->qfmt != eslMSAFILE_AFA) {
    /* this could just be a sequence file with o single sequence (in which case, fall through
     * to the "sequence" case), or with several same-sized sequences (in which case ask for guidance) */
    if (msa != NULL && msa->nseq > 1)
      p7_Fail("Query file type could be either aligned or unaligned; please specify (--qformat [afa|fasta])");
  } else if (status == eslOK) {
    /* if ok, and not fasta, then it's an MSA ... proceed */
    cfg->qfmt = (*qfp_msa)->format;
  }
  else {
    esl_msafile_Close(*qfp_msa);
    *qfp_msa = NULL;
  }

  esl_msa_Destroy(msa);
  return status;
}

static int
bath_open_seq_file (struct cfg_s *cfg, ESL_SQFILE **qfp_sq, ESL_ALPHABET **abc) {
 
  ESL_SQ *qsq = NULL;
  int status = esl_sqfile_Open(cfg->queryfile, cfg->qfmt, NULL, qfp_sq);

  if (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open query file %s.\n", cfg->queryfile);
    if (status == eslOK) {
      if (*abc == NULL) {
        int q_type = eslUNKNOWN;
        status = esl_sqfile_GuessAlphabet(*qfp_sq, &q_type);
        if (status == eslEFORMAT) p7_Fail("Parse failed (sequence file %s):\n%s\n", (*qfp_sq)->filename, esl_sqfile_GetErrorBuf(*qfp_sq));
         if (q_type == eslUNKNOWN) p7_Fail("Unable to guess alphabet for the %s%s query file %s\n", (cfg->qfmt==eslUNKNOWN ? "" : esl_sqio_DecodeFormat(cfg->qfmt)), (cfg->qfmt==eslSQFILE_UNKNOWN ? "":"-formatted"), cfg->queryfile);
           *abc = esl_alphabet_Create(q_type);
      }
      if ((*abc)->type != eslAMINO) p7_Fail("Invalid alphabet type in the %s%squery file %s. Expect Amino Acid\n", (cfg->qfmt==eslUNKNOWN ? "" : esl_sqio_DecodeFormat(cfg->qfmt)), (cfg->qfmt==eslSQFILE_UNKNOWN ? "":"-formatted "), cfg->queryfile);

        esl_sqfile_SetDigital(*qfp_sq, *abc);
        // read first sequence
        qsq = esl_sq_CreateDigital(*abc);
        status = esl_sqio_Read(*qfp_sq, qsq);
        if (status != eslOK) p7_Fail("reading sequence from file %s (%d): \n%s\n", cfg->queryfile, status, esl_sqfile_GetErrorBuf(*qfp_sq));
    }
    if (qsq!=NULL) esl_sq_Destroy(qsq);
    return status;
}

/* search_master()
 * Open the query and target files, then search every query HMM
 * against the target database in one pass over it.
 * 
 * A master can only return if it's successful. All errors are handled
 * immediately and fatally with p7_Fail().  We also use the
 * ESL_EXCEPTION and ERROR: mechanisms, but only because we know we're
 * using a fatal exception handler.
 */
static int
search_master(ESL_GETOPTS *go, struct cfg_s *cfg)
{

  /* output files */
  FILE            *ofp                      = stdout;            /* results output file (-o)                        */
  FILE            *tblfp                    = NULL;              /* output stream for tabular per-seq (--tblout)    */
  FILE            *exontblfp                = NULL;               /* output stream for tabular per-seq (--exontblout)    */
  FILE            *fstblfp                  = NULL;              /* output stream for tabular per-ali (--fstblout)  */
  char            *hmmfile                  = NULL;              /* file to write HMM to                            */
  int              textw                    = 0;

 /* input files */
  P7_HMMFILE      *hfp                      = NULL;              /* open input HMM file                             */
  ESL_SQFILE      *dbfp                     = NULL;              /* open input sequence file                        */
  ESL_MSAFILE     *qfp_msa                  = NULL;              /* open query alifile                              */
  ESL_SQFILE      *qfp_sq                   = NULL;              /* open query seqfile                              */
  int              dbfmt                    = eslSQFILE_UNKNOWN; /* format code for sequence database file          */

 /* query formats and HMM construction*/
  P7_HMM          *hmm                      = NULL;              /* one HMM query                                   */

  /* alphabets and translation */
  int              codon_table;
  ESL_ALPHABET    *abcAA                    = NULL;              /* AA  query  alphabet                                */
  ESL_ALPHABET    *abcDNA                   = NULL;              /* DNA target alphabet                              */
  ESL_GENCODE     *gcode                    = NULL;
 
  ESL_STOPWATCH   *watch;

  /* multi threading */
  int              ncpus                    = 0; 

  /*error handeling */
  char             errbuf[eslERRBUFSIZE];
  int              status                   = eslOK;
  int              ssistatus                = eslOK;  

  if (esl_opt_GetBoolean(go, "--notextw")) textw = 0;
  else                                     textw = esl_opt_GetInteger(go, "--textw");
 
  /* bathsearch accepts query files that are either hmm(s), msa(s), or sequence(s). The following 
   * code will follow the mandate of --qformat, and otherwise figure what the file type is. */

  /* (1) If we were told a specific query file type, just do what we're told */

  if (esl_sqio_IsAlignment(cfg->qfmt) /* msa file */) {

    /* First check that the user has provided an output file for the converted HMMs */
    if (esl_opt_IsOn(go, "--hmmout")) {
      hmmfile = esl_opt_GetString(go, "--hmmout");
      
    }
    else {
	  ESL_ALLOC(hmmfile, sizeof(char) * 256);
	  snprintf(hmmfile, 256, "/tmp/hmmfile_%d.bhmm", getpid());
    }

    status = bath_open_msa_file(cfg, &qfp_msa, &abcAA);
    if (status != eslOK) p7_Fail("Error reading msa from the %s-formatted file %s (%d)\n", esl_sqio_DecodeFormat(cfg->qfmt), cfg->queryfile, status);
  } else if (cfg->qfmt != eslSQFILE_UNKNOWN /* sequence file */) {
      /* First check that the user has provided an output file for the converted HMMs */
      if (esl_opt_IsOn(go, "--hmmout")) {
        hmmfile = esl_opt_GetString(go, "--hmmout");
      }
      else {
	    ESL_ALLOC(hmmfile, sizeof(char) * 256);
        snprintf(hmmfile, 256, "/tmp/hmmfile_%d.bhmm", getpid());
      }
 
      status = bath_open_seq_file(cfg, &qfp_sq, &abcAA);
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
          
	        /* First check that the user has provided an output file for the converted HMMs */ 
        if (esl_opt_IsOn(go, "--hmmout")) {
          hmmfile = esl_opt_GetString(go, "--hmmout");
        }
        else {
	      ESL_ALLOC(hmmfile, sizeof(char) * 256);
          snprintf(hmmfile, 256, "/tmp/hmmfile_%d.bhmm", getpid());
        }           

        status = bath_open_msa_file(cfg, &qfp_msa, &abcAA);

        if (cfg->qfmt == eslSQFILE_UNKNOWN) { /* it's not an MSA, try seq */
          status = bath_open_seq_file(cfg, &qfp_sq, &abcAA);
          if (status != eslOK) p7_Fail("Error reading query file %s (%d)\n", cfg->queryfile, status);
        }
      }
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

  if (esl_opt_IsOn(go, "-o"))           { if ((ofp       = fopen(esl_opt_GetString(go, "-o"),           "w")) == NULL)  p7_Fail("Failed to open output file %s for writing\n",    esl_opt_GetString(go, "-o")); }
  if (esl_opt_IsOn(go, "--tblout"))     { if ((tblfp     = fopen(esl_opt_GetString(go, "--tblout"),     "w")) == NULL)  esl_fatal("Failed to open tabular per-seq output file %s for writing\n", esl_opt_GetString(go, "--tblout")); }
  if (esl_opt_IsOn(go, "--exontblout")) { if ((exontblfp = fopen(esl_opt_GetString(go, "--exontblout"), "w")) == NULL)  esl_fatal("Failed to open tabular per-seq exon output file %s for writing\n", esl_opt_GetString(go, "--exontblout")); }
  if (esl_opt_IsOn(go, "--fstblout"))   { if ((fstblfp   = fopen(esl_opt_GetString(go, "--fstblout"),   "w")) == NULL)  esl_fatal("Failed to open tabular per-ali frameshift file %s for writing\n", esl_opt_GetString(go, "--fstblout")); }
  

#ifdef HMMER_THREADS
  ncpus = ESL_MIN(esl_opt_GetInteger(go, "--cpu"), esl_threads_GetCPUCount());
#endif

   /*the query sequence will be DNA but will be translated to amino acids */
   abcDNA = esl_alphabet_Create(eslDNA); 

  /* Get translation variables */
  codon_table = esl_opt_GetInteger(go, "--ct");
 
  if (status == eslOK)
  {
    /* One-time initializations after alphabet <abc> becomes known */
    p7_banner(ofp, go->argv[0], banner);
    output_header(ofp, go, cfg->queryfile, cfg->dbfile);
    esl_sqfile_SetDigital(dbfp, abcDNA); //ReadBlock requires knowledge of the alphabet to decide how best to read blocks
  }

   /* Set up the genetic code. Default = NCBI 1, the standard code; allow ORFs to start at any aa   */
  gcode = esl_gencode_Create(abcDNA, abcAA);
  esl_gencode_Set(gcode, codon_table);  // default = 1, the standard genetic code

  if      (esl_opt_GetBoolean(go, "-m"))   esl_gencode_SetInitiatorOnlyAUG(gcode);
  else if (! esl_opt_GetBoolean(go, "-M")) esl_gencode_SetInitiatorAny(gcode);      // note this is the default, if neither -m nor -M are set
  
  /* If query is alignment or sequence build the HMMs */
  if ( qfp_msa != NULL || qfp_sq != NULL ) {

    if (qfp_msa) esl_msafile_Close(qfp_msa);
    if (qfp_sq)  esl_sqfile_Close(qfp_sq);

    p7_search_builder(go, abcAA, cfg->queryfile, hmmfile, cfg->qfmt);
    cfg->queryfile = hmmfile;
    
    status = bath_open_hmm_file(cfg, &hfp, errbuf, &abcAA, &hmm);   
    if (status != eslOK) p7_Fail("Error reading hmms from %s (%d)\n", cfg->queryfile, status);
  }

  watch = esl_stopwatch_Create();

  /* every query HMM, in one pass over the target */
  scan_search(go, cfg, hfp, hmm, &abcAA, abcDNA, gcode, ncpus, dbfp, ofp, tblfp, exontblfp, fstblfp, textw, watch);

    /* Terminate outputs... any last words? */
  if (tblfp)      p7_tophits_TabularTail(tblfp,      "bathsearch", p7_SEARCH_SEQS, cfg->queryfile, cfg->dbfile, go);
  if (fstblfp)    p7_tophits_TabularTail(fstblfp,    "bathsearch", p7_SEARCH_SEQS, cfg->queryfile, cfg->dbfile, go); 
  if (exontblfp)  p7_tophits_TabularTail(exontblfp,  "bathsearch", p7_SEARCH_SEQS, cfg->queryfile, cfg->dbfile, go);
  if (ofp)      { if (fprintf(ofp, "[ok]\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed"); }

  /* Cleanup - prepare for exit */
  if (hfp) p7_hmmfile_Close(hfp);
  esl_sqfile_Close(dbfp);
  esl_alphabet_Destroy(abcAA);
  esl_alphabet_Destroy(abcDNA);
  esl_gencode_Destroy(gcode);
  esl_stopwatch_Destroy(watch);
  
  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (exontblfp)     fclose(exontblfp);
  if (fstblfp)       fclose(fstblfp);

  if (!esl_opt_IsOn(go, "--hmmout") && hmmfile != NULL) free(hmmfile);

  return eslOK;

ERROR:

  if (hfp) p7_hmmfile_Close(hfp);
  if (ofp != stdout) fclose(ofp);
  if (tblfp)         fclose(tblfp);
  if (exontblfp)     fclose(exontblfp);
  if (fstblfp)       fclose(fstblfp);
  if (!esl_opt_IsOn(go, "--hmmout") && hmmfile != NULL) free(hmmfile);
  
  return eslFAIL;
}

/*****************************************************************
 * The search driver.
 * The reader fills a small ring of block slots. A free worker translates
 * a whole slot once (both strands) into shared, read-only ORFs and
 * reverse-complemented DNA. A unit of work is one query over one
 * translated slot; a query's state (profiles, pipeline, hit list) is a
 * lane, used by one worker at a time, so nothing is copied per thread.
 * A query has one lane unless workers would otherwise wait for it. A
 * slot is freed when every query has run on it.
 *
 * Without worker threads (--cpu 0) there is one slot and one lane per
 * query: each block is read, translated and searched with every query
 * in turn.
 *****************************************************************/

typedef struct {
  P7_HMM          *hmm;
  P7_BG           *bg;
  P7_PROFILE      *gm;
  P7_OPROFILE     *om;
  P7_FS_PROFILE   *gm_fs5;
  P7_FS_PROFILE   *gm_fs3;
  P7_FS_OPROFILE  *om_fs3;
  P7_FS_OPROFILE  *om_fs5;
  P7_SCOREDATA    *scoredata;
  const ESL_GENCODE *gcode;
} SCAN_QUERY;

enum { SLOT_FREE = 0, SLOT_READING, SLOT_READ, SLOT_TRANSLATING, SLOT_READY };

typedef struct {
  ESL_SQ_BLOCK   *block;        /* DNA windows from the reader                  */
  ESL_SQ        **rc;           /* [block->count] reverse-complemented windows  */
  int             rc_alloc;
  ESL_ORF_BLOCK  *orf[2];       /* ORFs of all windows: [0] top, [1] bottom     */
  int            *orf_start[2]; /* [count+1] first ORF of each window           */
  int             st_alloc;
  int             state;
  int             ndone;
  unsigned char  *qdone;        /* [nq] query has run on this slot              */
  int64_t         seq;
} SCAN_SLOT;

/* --fs: what a lane's pipeline hooks need to build its frameshift profiles on demand */
typedef struct {
  SCAN_QUERY     *q;
  P7_BG          *bg;           /* the lane's own                                  */
  P7_FS_OPROFILE *om_fs3;       /* the lane's own: allocated, written on first use */
  int             ready;
} SCAN_FSARG;

typedef struct {
#ifdef HMMER_THREADS
  pthread_mutex_t  lock;
  pthread_cond_t   cv;          /* workers wait here for something to run or translate */
  pthread_cond_t   cv_rd;       /* the reader waits here for a free slot               */
#endif
  int              nslots;
  SCAN_SLOT      *slot;
  int              nq;
  int              L;           /* lanes a query starts with: states that can run it at the same time */
  int              maxl;        /* most lanes a query can have; stride of q[], qbusy[], fsa[]         */
  int             *nl;          /* [nq] lanes each query has                            */
  WORKER_INFO     *q;           /* [nq*maxl] one state per query lane                   */
  unsigned char   *qbusy;       /* [nq*maxl]                                            */
  SCAN_FSARG     *fsa;         /* [nq*maxl]                                            */
  SCAN_QUERY      *Q;           /* [nq]                                                 */
  ESL_GETOPTS     *go;
  ESL_GENCODE     *gcode;
  ESL_ALPHABET    *abcDNA;
  ESL_ALPHABET    *abcAA;
  int              strands;
  int              block_length;
  int              use_fs;
  int              eof;
} SCAN;

typedef struct {
  SCAN                 *s;
  ESL_GENCODE_WORKSTATE *wrk;
  P7_PIPELINE           *ws;        /* scratch lent to the query being run   */
  P7_FS_PROFILE         *gm5;       /* --fs: 5-codon profiles for frameshift domain  */
  P7_FS_OPROFILE        *om5;       /*   definition, configured for gm5_q on demand  */
  void                  *gm5_q;
} SCAN_WORKER;

#ifdef HMMER_THREADS
static __thread SCAN_WORKER *scan_tl_wk = NULL;   /* the worker running on this thread */
#else
static SCAN_WORKER *scan_tl_wk = NULL;
#endif

/* A pipeline's DP matrices, domain definition and RNG are scratch, only live
 * during one p7_Pipeline_BATH() call (the RNG is reseeded per region), so
 * queries hold none and borrow their worker's for each unit of work. */
static void
scan_swap_scratch(P7_PIPELINE *a, P7_PIPELINE *b)
{
  ESL_SWAP(a->oxf,    b->oxf,    P7_OMX *);
  ESL_SWAP(a->oxb,    b->oxb,    P7_OMX *);
  ESL_SWAP(a->fwd,    b->fwd,    P7_OMX *);
  ESL_SWAP(a->bck,    b->bck,    P7_OMX *);
  ESL_SWAP(a->oxf_fs, b->oxf_fs, P7_OMX *);
  ESL_SWAP(a->oxb_fs, b->oxb_fs, P7_OMX *);
  ESL_SWAP(a->fwd_fs, b->fwd_fs, P7_OMX *);
  ESL_SWAP(a->bck_fs, b->bck_fs, P7_OMX *);
  ESL_SWAP(a->ov3,    b->ov3,    P7_OIVX *);
  ESL_SWAP(a->ov5,    b->ov5,    P7_OIVX *);
  ESL_SWAP(a->r,      b->r,      ESL_RANDOMNESS *);
  ESL_SWAP(a->ddef,   b->ddef,   P7_DOMAINDEF *);
}

/* --fs: build a lane's om_fs3 the first time a DNA window reaches
 * the frameshift stage. It was allocated up front but not written, so a lane
 * that never gets there never faults it in. */
static int
scan_fs_prepare(void *arg)
{
  SCAN_FSARG   *a = (SCAN_FSARG *) arg;
  P7_FS_PROFILE *gm_fs3;
  if (a->ready) return eslOK;
  gm_fs3 = p7_profile_fs_Create(a->q->hmm->M, a->q->hmm->abc, p7P_3CODONS);   /* only needed to make om_fs3 */
  p7_ProfileConfig_fs(a->q->hmm, a->bg, a->q->gcode, gm_fs3, 100, p7_LOCAL);
  p7_fs_oprofile_Convert(gm_fs3, a->om_fs3);
  p7_profile_fs_Destroy(gm_fs3);
  a->ready = TRUE;
  return eslOK;
}

/* --fs: frameshift domain definition is rare (a few windows per query
 * per genome), so its 5-codon profiles, most of a query's --fs memory, are
 * configured into the worker's own pair when needed instead of kept per lane.
 * Their state doesn't carry between windows: domain definition restores
 * gm_fs5's length and the pipeline configures om_fs5 for each window. */
static int
scan_fs5_get(void *arg, P7_FS_OPROFILE **ret_om, P7_FS_PROFILE **ret_gm)
{
  SCAN_FSARG  *a  = (SCAN_FSARG *) arg;
  SCAN_QUERY   *q  = a->q;
  SCAN_WORKER *wk = scan_tl_wk;

  if (wk->gm5_q != (void *) q) {   /* a pair of this query's size: one for the largest query, kept by every worker, was most of the --fs memory */
    p7_profile_fs_Destroy(wk->gm5);
    p7_fs_oprofile_Destroy(wk->om5);
    wk->gm5 = p7_profile_fs_Create(q->hmm->M, q->hmm->abc, p7P_5CODONS);
    wk->om5 = p7_fs_oprofile_Create(q->hmm->M, q->hmm->abc, p7P_5CODONS);
    p7_ProfileConfig_fs(q->hmm, a->bg, q->gcode, wk->gm5, 100, p7_LOCAL);
    p7_fs_oprofile_Convert(wk->gm5, wk->om5);
    wk->gm5_q = q;
  }
  *ret_om = wk->om5;
  *ret_gm = wk->gm5;
  return eslOK;
}

/* After a unit, a worker's four full DP matrices go back to their starting
 * size. They grow to the largest window a unit meets, and would otherwise
 * stay that large in every worker for the rest of the search. Replacing
 * them after every unit costs no measurable time. */
static void
scan_reset_scratch(P7_PIPELINE *ws)
{
  p7_omx_Destroy(ws->fwd);     ws->fwd    = p7_omx_Create(100, 100, 100);
  p7_omx_Destroy(ws->bck);     ws->bck    = p7_omx_Create(100, 100, 100);
  p7_omx_Destroy(ws->fwd_fs);  ws->fwd_fs = p7_omx_Create_dpf(100, 100, 100, p7G_NSCELLS_FS);
  p7_omx_Destroy(ws->bck_fs);  ws->bck_fs = p7_omx_Create_dpf(100, 100, 100, p7G_NSCELLS);
  if (ws->fwd == NULL || ws->bck == NULL || ws->fwd_fs == NULL || ws->bck_fs == NULL) p7_Fail("allocation failure");
}

static void
scan_free_scratch(P7_PIPELINE *pli)
{
  P7_PIPELINE *tmp = calloc(1, sizeof(P7_PIPELINE));
  scan_swap_scratch(pli, tmp);
  p7_pipeline_Destroy_BATH(tmp);
}

/* Build lane <l> of query <k>: the state one search of the query needs.
 * The master profiles stay unchanged for the splice step, so om is a
 * copy; lanes after the first copy everything the search changes. With --fs
 * the frameshift profiles come from the pipeline hooks, on first use. The
 * caller holds the lane, or runs before the workers start. */
static void
scan_lane_create(SCAN *s, int k, int l)
{
  SCAN_QUERY  *q  = s->Q + k;
  WORKER_INFO *qi = s->q + k * s->maxl + l;

  qi->bg        = p7_bg_Create(s->abcAA);
  qi->gcode     = s->gcode;
  qi->hw        = p7_hmmwindow_CreateList();
  qi->th        = p7_tophits_Create();
  qi->om        = p7_oprofile_Clone(q->om);
  qi->gm        = (l == 0) ? q->gm : p7_profile_Clone(q->gm);   /* read-only in the search */
  qi->gm_fs5    = NULL;
  qi->om_fs5    = NULL;
  qi->om_fs3    = NULL;
  if (s->use_fs) qi->om_fs3 = (l == 0) ? q->om_fs3 : p7_fs_oprofile_Create(q->hmm->M, s->abcAA, p7P_3CODONS);
  qi->scoredata = (l == 0) ? q->scoredata : p7_hmm_ScoreDataClone(q->scoredata, q->om->abc->Kp);
  qi->pli       = p7_pipeline_Create_BATH(s->go, 100, 100, p7_SEARCH_SEQS);   /* small: its DP memory is freed just below */
  if (p7_pli_NewModel(qi->pli, qi->om, qi->bg) == eslEINVAL) p7_Fail(qi->pli->errbuf);
  scan_free_scratch(qi->pli);
  if (s->use_fs) {   /* scan_fs_prepare() builds om_fs3, and scan_fs5_get() supplies the 5-codon pair */
    SCAN_FSARG *fa = s->fsa + k * s->maxl + l;
    fa->q = q; fa->bg = qi->bg; fa->om_fs3 = qi->om_fs3; fa->ready = FALSE;
    qi->pli->fs_prepare = scan_fs_prepare; qi->pli->fs5_get = scan_fs5_get; qi->pli->fs_arg = fa;
  }
  qi->pli->strands      = s->strands;
  qi->pli->block_length = s->block_length;
}

static void
scan_translate(SCAN *s, SCAN_WORKER *wk, SCAN_SLOT *sl)
{
  ESL_SQ_BLOCK *b = sl->block;
  int i;

  if (b->count + 1 > sl->st_alloc) {
    sl->st_alloc = b->count + 1;
    sl->orf_start[0] = realloc(sl->orf_start[0], sizeof(int) * sl->st_alloc);
    sl->orf_start[1] = realloc(sl->orf_start[1], sizeof(int) * sl->st_alloc);
  }
  if (b->count > sl->rc_alloc) {
    sl->rc = realloc(sl->rc, sizeof(ESL_SQ *) * b->count);
    for (i = sl->rc_alloc; i < b->count; i++) sl->rc[i] = esl_sq_CreateDigital(s->abcDNA);
    sl->rc_alloc = b->count;
  }
  esl_gencode_OrfBlockReuse(sl->orf[0]);
  esl_gencode_OrfBlockReuse(sl->orf[1]);

  for (i = 0; i < b->count; i++) {
    ESL_SQ *dna = b->list + i;
    dna->L = dna->n;
    sl->orf_start[0][i] = sl->orf[0]->count;
    sl->orf_start[1][i] = sl->orf[1]->count;
    if (s->strands != p7_STRAND_BOTTOMONLY) {
      wk->wrk->orf_block = sl->orf[0];
      do_sq_by_sequences(s->gcode, wk->wrk, dna);
    }
    if (s->strands != p7_STRAND_TOPONLY) {
      esl_sq_Reuse(sl->rc[i]);
      esl_sq_Copy(dna, sl->rc[i]);
      esl_sq_ReverseComplement(sl->rc[i]);
      wk->wrk->orf_block = sl->orf[1];
      do_sq_by_sequences(s->gcode, wk->wrk, sl->rc[i]);
    }
  }
  sl->orf_start[0][b->count] = sl->orf[0]->count;
  sl->orf_start[1][b->count] = sl->orf[1]->count;
  wk->wrk->orf_block = NULL;
}

static void
scan_run_query(SCAN *s, SCAN_WORKER *wk, SCAN_SLOT *sl, int li)
{
  WORKER_INFO  *qi = s->q + li;
  ESL_ORF_BLOCK v = { 0 };
  int           i, st, n, strand;

  for (i = 0; i < sl->block->count; i++)
    for (strand = 0; strand < 2; strand++) {
      ESL_SQ *dna;
      if (strand == 0 && s->strands == p7_STRAND_BOTTOMONLY) continue;
      if (strand == 1 && s->strands == p7_STRAND_TOPONLY)    continue;
      dna = (strand == 0) ? sl->block->list + i : sl->rc[i];
      st  = sl->orf_start[strand][i];
      n   = sl->orf_start[strand][i+1] - st;
      v.count = v.listSize = n;
      v.list = sl->orf[strand]->list + st;   /* shared, read-only: the pipeline writes nothing into ORFs */

      qi->pli->nres += dna->W;
      scan_swap_scratch(qi->pli, wk->ws);
      p7_Pipeline_BATH(qi->pli, qi->om, qi->gm, qi->om_fs3, qi->om_fs5, qi->gm_fs5, qi->scoredata, qi->bg, qi->th, sl->block->first_seqidx + i, dna, &v, s->gcode, qi->hw, strand ? p7_COMPLEMENT : p7_NOCOMPLEMENT);
      p7_pipeline_Reuse_BATH(qi->pli);
      scan_swap_scratch(qi->pli, wk->ws);
    }
  scan_reset_scratch(wk->ws);
}

#ifdef HMMER_THREADS
static void *
scan_worker(void *arg)
{
  SCAN_WORKER *wk = (SCAN_WORKER *) arg;
  SCAN        *s  = wk->s;
  scan_tl_wk = wk;
  int           j, k, l, best, pick_k, pick_li;

  impl_Init();
  pthread_mutex_lock(&s->lock);
  /* A worker that takes work wakes one more worker, which looks for work in its
   * turn. So new work wakes as many workers as it has use for, one after
   * another, and the idle ones don't all queue for the lock at every change. */
  while (1)
    {
      /* 1. the oldest translated slot with a query still to run on it and a free lane of that
       *    query. Results don't depend on the order a query sees slots in. */
      best = -1; pick_k = -1; pick_li = -1;
      for (j = 0; j < s->nslots; j++) {
        SCAN_SLOT *sl = s->slot + j;
        if (sl->state != SLOT_READY || sl->ndone == s->nq) continue;
        if (best >= 0 && sl->seq > s->slot[best].seq)    continue;
        for (k = 0; k < s->nq; k++) {
          if (sl->qdone[k]) continue;
          for (l = 0; l < s->nl[k]; l++) if (!s->qbusy[k * s->maxl + l]) break;
          if (l < s->nl[k]) { best = j; pick_k = k; pick_li = k * s->maxl + l; break; }
        }
      }
      if (best >= 0) {
        SCAN_SLOT *sl = s->slot + best;
        int li = pick_li;
        s->qbusy[li] = 1;
        sl->qdone[pick_k] = 1;
        pthread_cond_signal(&s->cv);
        pthread_mutex_unlock(&s->lock);
        scan_run_query(s, wk, sl, li);
        pthread_mutex_lock(&s->lock);
        s->qbusy[li] = 0;
        if (++sl->ndone == s->nq) { sl->state = SLOT_FREE; pthread_cond_signal(&s->cv_rd); }
        continue;
      }
      /* 2. a slot to translate */
      best = -1;
      for (j = 0; j < s->nslots; j++)
        if (s->slot[j].state == SLOT_READ && (best < 0 || s->slot[j].seq < s->slot[best].seq)) best = j;
      if (best >= 0) {
        SCAN_SLOT *sl = s->slot + best;
        sl->state = SLOT_TRANSLATING;
        pthread_cond_signal(&s->cv);
        pthread_mutex_unlock(&s->lock);
        scan_translate(s, wk, sl);
        pthread_mutex_lock(&s->lock);
        sl->state = SLOT_READY;
        continue;
      }
      /* 3. nothing to run or translate. A query that still has a slot to do is held back only
       *    by its busy lanes, and is what the other workers wait for: give it another lane. */
      best = -1; pick_k = -1;
      for (j = 0; j < s->nslots; j++) {
        SCAN_SLOT *sl = s->slot + j;
        if (sl->state != SLOT_READY || sl->ndone == s->nq) continue;
        if (best >= 0 && sl->seq > s->slot[best].seq)    continue;
        for (k = 0; k < s->nq; k++)
          if (!sl->qdone[k] && s->nl[k] < s->maxl) { best = j; pick_k = k; break; }
      }
      if (best >= 0) {
        SCAN_SLOT *sl = s->slot + best;
        int li;
        l  = s->nl[pick_k]++;
        li = pick_k * s->maxl + l;
        s->qbusy[li] = 1;
        sl->qdone[pick_k] = 1;
        pthread_cond_signal(&s->cv);
        pthread_mutex_unlock(&s->lock);
        scan_lane_create(s, pick_k, l);
        scan_run_query(s, wk, sl, li);
        pthread_mutex_lock(&s->lock);
        s->qbusy[li] = 0;
        if (++sl->ndone == s->nq) { sl->state = SLOT_FREE; pthread_cond_signal(&s->cv_rd); }
        continue;
      }
      /* 4. done when the reader is finished and every slot is free */
      if (s->eof) {
        for (j = 0; j < s->nslots; j++) if (s->slot[j].state != SLOT_FREE) break;
        if (j == s->nslots) { pthread_cond_broadcast(&s->cv); break; }
      }
      pthread_cond_wait(&s->cv, &s->lock);
    }
  pthread_mutex_unlock(&s->lock);
  return NULL;
}

#endif /*HMMER_THREADS*/

/* The reader calls this before it fills a slot. A slot gets its blocks the
 * first time it is used, so a small target doesn't pay for the whole ring. */
static void
scan_slot_ready(SCAN *s, SCAN_SLOT *sl)
{
  if (sl->block != NULL) return;
  sl->block  = esl_sq_CreateDigitalBlock(BLOCK_SIZE, s->abcDNA);
  sl->orf[0] = esl_gencode_OrfBlockCreate(BLOCK_SIZE);
  sl->orf[1] = esl_gencode_OrfBlockCreate(BLOCK_SIZE);
}

/* Read the next block of the target into <b>. A window the last block ended
 * in the middle of is carried over in <tmpsq>, with <C> residues of context.
 * Returns the read status; <*abort> is set once --restrictdb_n is reached. */
static int
scan_read_block(ESL_SQFILE *dbfp, ESL_SQ_BLOCK *b, ESL_SQ *tmpsq, int *prev_complete, ID_LENGTH_LIST *id_length_list,
                int block_length, int C, int n_targetseqs, int64_t *nseqs, int *abort)
{
  int64_t seqid;
  int     i, sstatus;

  b->complete = *prev_complete;
  if (! *prev_complete) {
    esl_sq_Copy(tmpsq, b->list);
    b->list->C = (b->list->n < C) ? b->list->n : C;
  }
  sstatus = esl_sqio_ReadBlock(dbfp, b, block_length, n_targetseqs, FALSE, TRUE);

  b->first_seqidx = *nseqs;
  seqid = *nseqs;
  for (i = 0; i < b->count; i++) {
    b->list[i].idx = seqid;
    add_id_length(id_length_list, seqid, b->list[i].L);
    seqid++;
    if (seqid == n_targetseqs && (i < b->count-1 || b->complete)) { *abort = TRUE; b->count = i+1; break; }
  }
  *nseqs += b->count - ((*abort || b->complete) ? 0 : 1);
  *prev_complete = b->complete;
  if (!b->complete && b->count > 0) esl_sq_Copy(b->list + (b->count - 1), tmpsq);
  return sstatus;
}

#ifdef HMMER_THREADS
/* One pass over the target for queries q[0..nq-1]; the reader runs here. */
static int
scan_pass(SCAN *s, int ncpus, ESL_SQFILE *dbfp, ID_LENGTH_LIST *id_length_list, int block_length, int C, int n_targetseqs, int64_t *ret_nseqs)
{
  pthread_t    *tid = malloc(sizeof(pthread_t) * ncpus);
  SCAN_WORKER *wk  = calloc(ncpus, sizeof(SCAN_WORKER));
  ESL_SQ       *tmpsq = esl_sq_CreateDigital(dbfp->abc);
  int           prev_complete = TRUE, abort = FALSE, sstatus = eslOK, j, t;
  int64_t       nseqs = 0, order = 0;

  for (j = 0; j < s->nslots; j++) { s->slot[j].state = SLOT_FREE; s->slot[j].ndone = 0; }
  s->eof = FALSE;
  for (t = 0; t < ncpus; t++) {
    wk[t].s   = s;
    wk[t].wrk = esl_gencode_WorkstateCreate(s->go, s->gcode);
    wk[t].ws  = p7_pipeline_Create_BATH(s->go, 100, 100, p7_SEARCH_SEQS);
    pthread_create(&tid[t], NULL, scan_worker, &wk[t]);
  }

  while (!abort)
    {
      SCAN_SLOT   *sl;
      pthread_mutex_lock(&s->lock);
      for (;;) {
        for (j = 0; j < s->nslots; j++) if (s->slot[j].state == SLOT_FREE) break;
        if (j < s->nslots) break;
        pthread_cond_wait(&s->cv_rd, &s->lock);
      }
      sl = s->slot + j;
      sl->state = SLOT_READING;
      pthread_mutex_unlock(&s->lock);

      scan_slot_ready(s, sl);
      sstatus = scan_read_block(dbfp, sl->block, tmpsq, &prev_complete, id_length_list, block_length, C, n_targetseqs, &nseqs, &abort);

      pthread_mutex_lock(&s->lock);
      if (sstatus != eslOK || sl->block->count == 0) {
        sl->state = SLOT_FREE;
        s->eof = TRUE;
        pthread_cond_broadcast(&s->cv);
        pthread_mutex_unlock(&s->lock);
        break;
      }
      memset(sl->qdone, 0, s->nq);
      sl->ndone = 0;
      sl->seq   = order++;
      sl->state = SLOT_READ;
      pthread_cond_signal(&s->cv);
      pthread_mutex_unlock(&s->lock);
    }
  pthread_mutex_lock(&s->lock);
  s->eof = TRUE;
  pthread_cond_broadcast(&s->cv);
  pthread_mutex_unlock(&s->lock);

  for (t = 0; t < ncpus; t++) {
    pthread_join(tid[t], NULL);
    esl_gencode_WorkstateDestroy(wk[t].wrk);
    p7_pipeline_Destroy_BATH(wk[t].ws);
    p7_profile_fs_Destroy(wk[t].gm5);
    p7_fs_oprofile_Destroy(wk[t].om5);
  }
  free(tid); free(wk);
  esl_sq_Destroy(tmpsq);
  *ret_nseqs = nseqs;
  return (sstatus == eslEOF || sstatus == eslOK) ? eslOK : sstatus;
}
#endif /*HMMER_THREADS*/

/* One pass over the target with no worker threads. */
static int
scan_pass_serial(SCAN *s, ESL_SQFILE *dbfp, ID_LENGTH_LIST *id_length_list, int block_length, int C, int n_targetseqs, int64_t *ret_nseqs)
{
  SCAN_WORKER  wk;
  SCAN_SLOT   *sl    = s->slot;
  ESL_SQ      *tmpsq = esl_sq_CreateDigital(dbfp->abc);
  int          prev_complete = TRUE, abort = FALSE, sstatus = eslOK, k;
  int64_t      nseqs = 0;

  memset(&wk, 0, sizeof(wk));
  wk.s   = s;
  wk.wrk = esl_gencode_WorkstateCreate(s->go, s->gcode);
  wk.ws  = p7_pipeline_Create_BATH(s->go, 100, 100, p7_SEARCH_SEQS);
  scan_tl_wk = &wk;
  scan_slot_ready(s, sl);

  while (!abort)
    {
      sstatus = scan_read_block(dbfp, sl->block, tmpsq, &prev_complete, id_length_list, block_length, C, n_targetseqs, &nseqs, &abort);
      if (sstatus != eslOK || sl->block->count == 0) break;
      scan_translate(s, &wk, sl);
      for (k = 0; k < s->nq; k++) scan_run_query(s, &wk, sl, k * s->maxl);
    }

  scan_tl_wk = NULL;
  esl_gencode_WorkstateDestroy(wk.wrk);
  p7_pipeline_Destroy_BATH(wk.ws);
  p7_profile_fs_Destroy(wk.gm5);
  p7_fs_oprofile_Destroy(wk.om5);
  esl_sq_Destroy(tmpsq);
  *ret_nseqs = nseqs;
  return (sstatus == eslEOF || sstatus == eslOK) ? eslOK : sstatus;
}

/* Search every query in one pass over the target */
static int
scan_search(ESL_GETOPTS *go, struct cfg_s *cfg, P7_HMMFILE *hfp, P7_HMM *hmm, ESL_ALPHABET **p_abcAA, ESL_ALPHABET *abcDNA, ESL_GENCODE *gcode,
             int ncpus, ESL_SQFILE *dbfp, FILE *ofp, FILE *tblfp, FILE *exontblfp, FILE *fstblfp, int textw, ESL_STOPWATCH *watch)
{
  ESL_ALPHABET      *abcAA   = *p_abcAA;
  int                use_fs  = (esl_opt_IsUsed(go, "--fs") || esl_opt_IsUsed(go, "--fsonly"));
  int                splice  = esl_opt_IsUsed(go, "--splice");
  int                codon_table = esl_opt_GetInteger(go, "--ct");
  int                batch;
  int                strands, block_length, qhstatus = eslOK, nout = 0;
  int                b0, nb, k, j, t, d, C, sstatus;
  int64_t            nseqs, resCnt;
  SCAN              s;
  SCAN_QUERY        *Q;
  ID_LENGTH_LIST    *id_length_list;
  P7_TOPHITS        *seed_hits;
  P7_HMM_WINDOWLIST *seed_accumulator;
  P7_FS_PROFILE     *gm_tr;
  P7_HMM           **H = NULL;
  int                nq = 0, qalloc = 0;

  if      (strcmp(esl_opt_GetString(go, "--strand"), "both")  == 0) strands = p7_STRAND_BOTH;
  else if (strcmp(esl_opt_GetString(go, "--strand"), "plus")  == 0) strands = p7_STRAND_TOPONLY;
  else                                                                strands = p7_STRAND_BOTTOMONLY;
  block_length = esl_opt_IsUsed(go, "--block_length") ? esl_opt_GetInteger(go, "--block_length") : BATH_MAX_RESIDUE_COUNT;

  /* read every query */
  while (qhstatus == eslOK) {
    if (use_fs) { //check that HMM is properly formated for bathsearch
      if( !(hmm->flags & p7H_STATS) )
        p7_Fail("HMM file %s has no E-value statistics, which bathsearch requires.\nRebuild with 'bathbuild --fs', or add them with 'bathconvert --fs new_file.bhmm %s'.\n", cfg->queryfile, cfg->queryfile);

      if( !(hmm->fsprob && hmm->ct)                      ||
          hmm->evparam[p7_FTAUFS3] == p7_EVPARAM_UNSET   ||
          hmm->evparam[p7_FTAUFS5] == p7_EVPARAM_UNSET )
        p7_Fail("HMM file %s has no frameshift statistics, which --fs requires.\nRebuild with 'bathbuild --fs', or add them with 'bathconvert --fs new_file.bhmm %s'.\n", cfg->queryfile, cfg->queryfile);

      /* frameshift E-values are computed from a specific codon table, so --fs/--fsonly
       * requires the HMM's table to match the one bathsearch is using */
      if( hmm->ct != codon_table)  p7_Fail("Requested codon translation tabel ID %d does not match the codon translation tabel ID of the HMM file %s. Please either run bathsearch with option '--ct %d' or run bathconvert with option '--ct %d'.\n", codon_table, cfg->queryfile, hmm->ct, codon_table);
    } else {
      if( !(hmm->flags & p7H_STATS) )
        p7_Fail("HMM file %s has no E-value statistics, which bathsearch requires.\nRebuild with 'bathbuild' (without --nostats), or add them with 'bathconvert --addstats new_file.bhmm %s'.\n", cfg->queryfile, cfg->queryfile);

      hmm->fs = FALSE;
      hmm->fsprob = 0.;
    }
    if (hmm->max_length == -1) p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);
    if (nq == qalloc) { qalloc = qalloc ? qalloc*2 : 64; H = realloc(H, sizeof(P7_HMM *) * qalloc); }
    H[nq++] = hmm;
    hmm = NULL;
    qhstatus = p7_hmmfile_Read(hfp, p_abcAA, &hmm);
    if (qhstatus != eslOK && qhstatus != eslEOF) p7_Fail("reading from query file %s (%d)\n", cfg->queryfile, qhstatus);
  }
  batch = nq;   /* all queries in one pass */
  s.L = ESL_MAX(1, (ncpus + batch - 1) / batch);   /* rounded up, so every worker has a lane */

#ifdef HMMER_THREADS
  pthread_mutex_init(&s.lock, NULL);
  pthread_cond_init(&s.cv, NULL);
  pthread_cond_init(&s.cv_rd, NULL);
#endif
  s.go = go; s.gcode = gcode; s.abcDNA = abcDNA; s.abcAA = abcAA; s.strands = strands;
  s.block_length = block_length; s.use_fs = use_fs;
  s.nslots = (ncpus > 0) ? ESL_MAX(4, 2*s.L + 2) : 1;
  s.maxl   = ESL_MAX(s.L, ESL_MIN(ncpus, s.nslots));   /* a query can run on one slot per lane */
  s.slot   = calloc(s.nslots, sizeof(SCAN_SLOT));
  for (j = 0; j < s.nslots; j++) s.slot[j].qdone = calloc(batch, 1);   /* its blocks come on first use: scan_slot_ready() */
  Q     = calloc(batch, sizeof(SCAN_QUERY));
  s.Q     = Q;
  s.q     = calloc(batch * s.maxl, sizeof(WORKER_INFO));
  s.qbusy = calloc(batch * s.maxl, 1);
  s.fsa   = calloc(batch * s.maxl, sizeof(SCAN_FSARG));
  s.nl    = calloc(batch, sizeof(int));

  for (b0 = 0; b0 < nq; b0 += batch)
    {
      nb = ESL_MIN(batch, nq - b0);
      s.nq = nb;
      C = 0;
      for (k = 0; k < nb; k++) {       /* one state per query */
        SCAN_QUERY  *q  = Q + k;
        P7_HMM      *h  = H[b0 + k];
        int          l;
        q->hmm = h;
        q->bg  = p7_bg_Create(abcAA);
        q->gcode    = gcode;
        q->gm_fs5 = NULL; q->gm_fs3 = NULL; q->om_fs3 = NULL; q->om_fs5 = NULL;
        q->gm = p7_profile_Create(h->M, abcAA);
        q->om = p7_oprofile_Create(h->M, abcAA);
        p7_ProfileConfig(h, q->bg, q->gm, 100, p7_LOCAL);
        p7_oprofile_Convert(q->gm, q->om);
        if (use_fs) q->om_fs3 = p7_fs_oprofile_Create(h->M, abcAA, p7P_3CODONS);   /* lane 0's; frameshift profiles only with --fs */
        q->scoredata = p7_hmm_ScoreDataCreate(q->om, NULL);
        p7_hmm_ScoreDataComputeRest(q->om, q->scoredata);
        C = ESL_MAX(C, q->om->max_length*3);

        s.nl[k] = s.L;
        for (l = 0; l < s.nl[k]; l++) scan_lane_create(&s, k, l);
      }

      if (b0 > 0 && esl_sqfile_Position(dbfp, 0) != eslOK) p7_Fail("can't rewind target file");
      if (cfg->firstseq_key != NULL && esl_sqfile_PositionByKey(dbfp, cfg->firstseq_key) != eslOK)
        p7_Fail("Failure setting restrictdb_stkey to %s\n", cfg->firstseq_key);
      esl_stopwatch_Start(watch);
      id_length_list = init_id_length(1000);
#ifdef HMMER_THREADS
      if (ncpus > 0) sstatus = scan_pass(&s, ncpus, dbfp, id_length_list, block_length, C, cfg->n_targetseq, &nseqs);
      else
#endif
                     sstatus = scan_pass_serial(&s, dbfp, id_length_list, block_length, C, cfg->n_targetseq, &nseqs);
      if (sstatus == eslEFORMAT) esl_fatal("Parse failed (sequence file %s):\n%s\n", dbfp->filename, esl_sqfile_GetErrorBuf(dbfp));
      else if (sstatus != eslOK) esl_fatal("Unexpected error %d reading sequence file %s", sstatus, dbfp->filename);

      for (k = 0; k < nb; k++)
        {
          SCAN_QUERY  *q  = Q + k;
          WORKER_INFO *ql = s.q + k * s.maxl;   /* this query's lanes */
          P7_HMM      *qh = q->hmm;
          P7_TOPHITS  *th;
          P7_PIPELINE *pl;
          int          l;

          if (fprintf(ofp, "Query:       %s  [M=%d]\n", qh->name, qh->M) < 0) p7_Fail("write failed");
          if (qh->acc)  fprintf(ofp, "Accession:   %s\n", qh->acc);
          if (qh->desc) fprintf(ofp, "Description: %s\n", qh->desc);

          ql[0].pli->nseqs = nseqs;   /* sequences are counted once, on the first lane */
          resCnt = 0;
          if (esl_opt_IsUsed(go, "-Z")) { resCnt = 1000000*esl_opt_GetReal(go, "-Z"); if (strands == p7_STRAND_BOTH) resCnt *= 2; }
          else for (l = 0; l < s.nl[k]; l++) resCnt += ql[l].pli->nres;
          for (l = 0; l < s.nl[k]; l++) p7_tophits_ComputeEvalues_BATH(ql[l].th, resCnt, ql[l].om->max_length*3);

          /* merge the lanes */
          th = p7_tophits_Create();
          pl = p7_pipeline_Create_BATH(go, 100, 300, p7_SEARCH_SEQS);
          pl->nmodels = 1; pl->nnodes = qh->M;
          seed_accumulator = splice ? p7_hmmwindow_CreateList() : NULL;
          for (l = 0; l < s.nl[k]; l++) {
            WORKER_INFO *qi = ql + l;
            p7_tophits_Merge(th, qi->th);
            p7_pipeline_Merge(pl, qi->pli);
            if (splice) p7_hmmwindow_Merge(seed_accumulator, qi->hw);
            p7_pipeline_Destroy_BATH(qi->pli);
            p7_tophits_Destroy(qi->th);
            p7_hmmwindow_DestroyList(qi->hw);
            p7_oprofile_Destroy(qi->om);
            if (qi->gm != q->gm) p7_profile_Destroy(qi->gm);
            p7_profile_fs_Destroy(qi->gm_fs5);
            if (l > 0) {
              p7_fs_oprofile_Destroy(qi->om_fs3);
              p7_fs_oprofile_Destroy(qi->om_fs5);
              p7_hmm_ScoreDataDestroy(qi->scoredata);
            }
            p7_bg_Destroy(qi->bg);
          }

          p7_tophits_SortBySeqidxAndAlipos(th);
          if (!splice) assign_Lengths(th, id_length_list);
          p7_tophits_RemoveDuplicates(th, pl->use_bit_cutoffs);
          p7_tophits_SortBySortkey(th);
          pl->Z = 1;
          p7_tophits_Threshold(th, pl);

          gm_tr = NULL;
          if (splice && th->N) {
            gm_tr = p7_profile_fs_Create(qh->M, abcAA, 1);
            p7_ProfileConfig_fs(qh, q->bg, gcode, gm_tr, 100, p7_UNILOCAL);
            p7_tophits_SortBySeqidxAndAlipos(th);
            p7_hmmwindow_RemoveDuplicates(seed_accumulator, th, pl->F3);
            seed_hits = p7_hmmwindow_GetSeedHits(seed_accumulator, th, qh, q->gm, dbfp, gcode, pl->F3, esl_opt_GetInteger(go, "--max_intron"));
            p7_splice_SpliceHits(th, seed_hits, q->om, q->gm, gm_tr, go, gcode, dbfp, id_length_list, resCnt);
            for (d = 0; d < seed_hits->N; d++) {
              p7_trace_fs_Destroy(seed_hits->unsrt[d].dcl->tr);
              free(seed_hits->unsrt[d].dcl->scores_per_pos);
              free(seed_hits->unsrt[d].dcl->k_per_pos);
            }
            p7_tophits_Destroy(seed_hits);
            assign_Lengths(th, id_length_list);
            p7_tophits_RemoveDuplicates(th, pl->use_bit_cutoffs);
            p7_tophits_SortBySortkey(th);
          }

          pl->n_output = pl->pos_output = 0;
          for (t = 0; t < th->N; t++)
            if ((th->hit[t]->flags & p7_IS_REPORTED) || th->hit[t]->flags & p7_IS_INCLUDED) {
              pl->n_output++;
              for (d = 0; d < th->hit[t]->ndom; d++)
                pl->pos_output += 1 + (th->hit[t]->dcl[d].jali > th->hit[t]->dcl[d].iali ? th->hit[t]->dcl[d].jali - th->hit[t]->dcl[d].iali : th->hit[t]->dcl[d].iali - th->hit[t]->dcl[d].jali);
            }
          p7_tophits_Targets(ofp, th, pl, textw); fprintf(ofp, "\n\n");
          p7_tophits_Domains(ofp, th, pl, textw); fprintf(ofp, "\n\n");
          if (tblfp)     p7_tophits_TabularTargets    (tblfp,     qh->name, qh->acc, th, pl, (nout == 0));
          if (exontblfp) p7_tophits_TabularExons      (exontblfp, qh->name, qh->acc, th, pl, (nout == 0), esl_opt_IsUsed(go, "--nodeinfo"));
          if (fstblfp)   p7_tophits_TabularFrameshifts(fstblfp,   qh->name, qh->acc, th, pl, (nout == 0));
          nout++;
          esl_stopwatch_Stop(watch);
          p7_pli_Statistics(ofp, pl, watch);
          fprintf(ofp, "//\n");

          p7_pipeline_Destroy_BATH(pl);
          p7_tophits_Destroy(th);
          p7_hmmwindow_DestroyList(seed_accumulator);
          p7_profile_fs_Destroy(gm_tr);
          p7_oprofile_Destroy(q->om);
          p7_profile_Destroy(q->gm);
          p7_profile_fs_Destroy(q->gm_fs5);
          p7_profile_fs_Destroy(q->gm_fs3);
          p7_fs_oprofile_Destroy(q->om_fs3);
          p7_fs_oprofile_Destroy(q->om_fs5);
          p7_hmm_ScoreDataDestroy(q->scoredata);
          p7_bg_Destroy(q->bg);
          p7_hmm_Destroy(qh);
        }
      destroy_id_length(id_length_list);
    }

  for (j = 0; j < s.nslots; j++) {
    if (s.slot[j].block != NULL) esl_sq_DestroyBlock(s.slot[j].block);
    esl_gencode_OrfBlockDestroy(s.slot[j].orf[0]);
    esl_gencode_OrfBlockDestroy(s.slot[j].orf[1]);
    for (k = 0; k < s.slot[j].rc_alloc; k++) esl_sq_Destroy(s.slot[j].rc[k]);
    free(s.slot[j].rc); free(s.slot[j].orf_start[0]); free(s.slot[j].orf_start[1]); free(s.slot[j].qdone);
  }
  free(s.slot); free(Q); free(s.q); free(s.qbusy); free(s.fsa); free(s.nl); free(H);
#ifdef HMMER_THREADS
  pthread_mutex_destroy(&s.lock);
  pthread_cond_destroy(&s.cv);
  pthread_cond_destroy(&s.cv_rd);
#endif
  return eslOK;
}

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
add_id_length(ID_LENGTH_LIST *list, int id, int64_t L)
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
    if(th->hit[i]->dcl[0].ad != NULL) th->hit[i]->dcl[0].ad->L = id_length_list->id_lengths[j].length; 
  }

  return eslOK;
}


