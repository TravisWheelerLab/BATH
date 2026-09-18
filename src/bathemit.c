/* bathemit: sample sequence(s) from a profile HMM.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_random.h"
#include "esl_sq.h"
#include "esl_sqio.h"

#include "hmmer.h"

/* Stringify p7P_FSPROB at compile time, so --fsprob's displayed default can
 * never drift out of sync with the constant it actually defaults to. */
#define bathemit_STR(x)  #x
#define bathemit_XSTR(x) bathemit_STR(x)

#define MODEOPTS "--local,--unilocal,--glocal,--uniglocal"
#define EMITOPTS "-a,-c,-C,-p"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs      incomp  help   docgroup*/
  { "-h",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,      NULL,  "show brief help on version and usage",             1 },
  { "-o",          eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,      NULL,  "send sequence output to file <f>, not stdout",     1 },
  { "-N",          eslARG_INT,      "1", NULL, "n>0",     NULL,      NULL,   "-c,-C",  "number of seqs to sample",                         1 },
  { "--seed",      eslARG_INT,      "0", NULL, "n>=0",    NULL,      NULL,      NULL,  "set RNG seed to <n>",                              1 },
  { "-a",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  EMITOPTS,  "emit alignment of sampled sequences",              2 },
  { "-c",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  EMITOPTS,  "emit simple majority-rule consensus sequence",     2 },
  { "-C",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  EMITOPTS,  "emit fancier consensus (needs --minl, --minu)",    2 },
  { "-p",          eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,  EMITOPTS,  "sample from profile, not core model",              2 },
  { "-L",          eslARG_INT,    "400", NULL, NULL,      NULL,      "-p",      NULL,  "with -p, set expected sequence length to <n>",     3 },
  { "--local",     eslARG_NONE,"default",NULL, NULL,   MODEOPTS,     "-p",      NULL,  "with -p, configure profile in multihit local mode",   3 },
  { "--unilocal",  eslARG_NONE,   FALSE, NULL, NULL,   MODEOPTS,     "-p",      NULL,  "with -p, configure profile in unihit local mode",     3 },
  { "--glocal",    eslARG_NONE,   FALSE, NULL, NULL,   MODEOPTS,     "-p",      NULL,  "with -p, configure profile in multihit glocal mode",  3 },
  { "--uniglocal", eslARG_NONE,   FALSE, NULL, NULL,   MODEOPTS,     "-p",      NULL,  "with -p, configure profile in unihit glocal mode",    3 },
  { "--minl",      eslARG_REAL,   "0.0", NULL, "0<=x<=1", NULL,      "-C",      NULL,  "show 'any' (X/N) below this consensus fraction",   4 },
  { "--minu",      eslARG_REAL,   "0.0", NULL, "0<=x<=1", NULL,      "-C",      NULL,  "show upper case at/above this consensus fraction", 4 },
  { "--dna",       eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,      NULL,  "emit DNA, back-translated via a codon table",      5 },
  { "--ct",        eslARG_INT,      "1", NULL, NULL,      NULL,   "--dna",      NULL,  "with --dna, use alt genetic code of NCBI transl table <n>", 5 },
  { "--fs",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,   "--dna",   "-c,-C",  "with --dna, inject frameshifts at rate --fsprob",  5 },
  { "--fsprob",    eslARG_REAL, bathemit_XSTR(p7P_FSPROB), NULL, "0.001<=x<=0.1", NULL, "--fs", NULL, "with --fs, per-codon frameshift rate", 5 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[]  = "[-options] <hmmfile>";
static char banner[] = "sample sequence(s) from a profile HMM";

/* One residue's rolled quasicodon, broken into the 6-slot frame used by
 * dna_align_fs(): [insert-before-1] pos1 [insert-1-2] pos2 [insert-2-3] pos3.
 * <has> is FALSE for a column where this sequence has no residue at all
 * (a plain amino-alignment gap); the three n_* counts are 0, 1, or 2. */
typedef struct {
  int     has;
  int     n_before; ESL_DSQ before[2];
  int     has1;     ESL_DSQ pos1;
  int     n_mid12;  ESL_DSQ mid12[2];
  int     has2;     ESL_DSQ pos2;
  int     n_mid23;  ESL_DSQ mid23[2];
  int     has3;     ESL_DSQ pos3;
} FS_SLOT;

static void cmdline_failure(char *argv0, char *format, ...);
static void cmdline_help(char *argv0, ESL_GETOPTS *go);

static void emit_sequences(ESL_GETOPTS *go, FILE *ofp, ESL_RANDOMNESS *r, P7_HMM *hmm,
                            P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA);
static void emit_one(FILE *ofp, ESL_RANDOMNESS *r, ESL_SQ *sq, const char *name,
                      P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA,
                      int do_fs, double fsprob);
static int  emit_fs_codon(ESL_RANDOMNESS *r, const ESL_DSQ *codon, ESL_DSQ *out, double fsprob);
static void emit_fs_slots(ESL_RANDOMNESS *r, const ESL_DSQ *codon, double fsprob, FS_SLOT *slot);
static ESL_MSA *dna_align(ESL_RANDOMNESS *r, const ESL_MSA *msa, const ESL_ALPHABET *abcAmino,
                           P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA);
static ESL_MSA *dna_align_fs(ESL_RANDOMNESS *r, const ESL_MSA *msa, const ESL_ALPHABET *abcAmino,
                              P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA, double fsprob);


int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go          = NULL;             /* command line processing                 */
  ESL_ALPHABET    *abc         = NULL;             /* sequence alphabet                       */
  ESL_ALPHABET    *abcDNA      = NULL;             /* DNA alphabet, only if --dna              */
  ESL_GENCODE     *gcode       = NULL;             /* genetic code, only if --dna              */
  P7_CODONTABLE   *codon_table = NULL;             /* amino->codon lookup, only if --dna       */
  ESL_RANDOMNESS  *r           = NULL;             /* source of randomness                    */
  char            *hmmfile     = NULL;             /* file to read HMM(s) from                */
  P7_HMMFILE      *hfp         = NULL;             /* open hmmfile                            */
  P7_HMM          *hmm         = NULL;             /* HMM to emit from                        */
  FILE            *ofp         = NULL;	            /* output stream                           */
  int              do_dna      = 0;
  int              nhmms       = 0;
  int              status;
  char             errbuf[eslERRBUFSIZE];

  go = esl_getopts_Create(options);
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) cmdline_failure(argv[0], "Failed to parse command line: %s\n", go->errbuf);
  if (esl_opt_VerifyConfig(go)               != eslOK) cmdline_failure(argv[0], "Error in configuration: %s\n",       go->errbuf);
  if (esl_opt_GetBoolean(go, "-h"))                    cmdline_help   (argv[0], go);
  if (esl_opt_ArgNumber(go) != 1)                      cmdline_failure(argv[0], "Incorrect number of command line arguments.\n");

  if ((hmmfile = esl_opt_GetArg(go, 1)) == NULL)       cmdline_failure(argv[0], "Failed to get <hmmfile> on cmdline: %s\n", go->errbuf);

  do_dna = esl_opt_GetBoolean(go, "--dna");

  if ( esl_opt_IsOn(go, "-o") ) {
    if ((ofp = fopen(esl_opt_GetString(go, "-o"), "w")) == NULL) esl_fatal("Failed to open output file %s", esl_opt_GetString(go, "-o"));
  } else ofp = stdout;

  r = esl_randomness_CreateFast(esl_opt_GetInteger(go, "--seed"));

  status = p7_hmmfile_OpenE(hmmfile, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile, errbuf);

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) != eslEOF)
    {
      if      (status == eslEFORMAT)    esl_fatal("Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
      else if (status == eslEINCOMPAT)  esl_fatal("HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(abc->type));
      else if (status != eslOK)         esl_fatal("Unexpected error in reading HMMs from %s\n",     hfp->fname);

      if (nhmms == 0 && do_dna) { 	/* first time initialization, now that alphabet known */
        if (abc->type != eslAMINO) esl_fatal("--dna requires a protein HMM; %s is %s, not amino\n", hfp->fname, esl_abc_DecodeType(abc->type));
        if ((abcDNA = esl_alphabet_Create(eslDNA))          == NULL) esl_fatal("failed to create DNA alphabet");
        if ((gcode  = esl_gencode_Create(abcDNA, abc))      == NULL) esl_fatal("failed to create genetic code");
        if (esl_gencode_Set(gcode, esl_opt_GetInteger(go, "--ct")) != eslOK) esl_fatal("failed to set genetic code translation table");
        if ((codon_table = p7_codontable_Create(gcode))     == NULL) esl_fatal("failed to create codon table");
      }
      nhmms++;

      emit_sequences(go, ofp, r, hmm, codon_table, abcDNA);

      p7_hmm_Destroy(hmm);
    }
  if (nhmms == 0) esl_fatal("Empty HMM file %s? No HMM data found.\n", hmmfile);

  if (esl_opt_IsOn(go, "-o")) { fclose(ofp); }
  if (codon_table) p7_codontable_Destroy(codon_table);
  if (gcode)       esl_gencode_Destroy(gcode);
  if (abcDNA)      esl_alphabet_Destroy(abcDNA);
  esl_randomness_Destroy(r);
  esl_alphabet_Destroy(abc);
  esl_getopts_Destroy(go);
  p7_hmmfile_Close(hfp);
  return eslOK;
}


static void
cmdline_failure(char *argv0, char *format, ...)
{
  va_list argp;
  printf("\nERROR: ");
  va_start(argp, format);
  vfprintf(stdout, format, argp);
  va_end(argp);
  esl_usage(stdout, argv0, usage);
  printf("\nTo see more help on available options, do %s -h\n\n", argv0);
  exit(1);
}

static void
cmdline_help(char *argv0, ESL_GETOPTS *go)
{
  p7_banner (stdout, argv0, banner);
  esl_usage (stdout, argv0, usage);
  puts("\nBasic options:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  puts("\nOptions controlling what to emit:");
  esl_opt_DisplayHelp(stdout, go, 2, 2, 100);
  puts("\nOptions controlling emission from profiles with -p:");
  esl_opt_DisplayHelp(stdout, go, 3, 2, 100);
  puts("\nOptions controlling fancy consensus emission with -C:");
  esl_opt_DisplayHelp(stdout, go, 4, 2, 100);
  puts("\nOptions controlling DNA output:");
  esl_opt_DisplayHelp(stdout, go, 5, 2, 100);
  puts("\nAvailable NCBI genetic code tables (for --ct <id>):");
  esl_gencode_DumpAltCodeTable(stdout);
  exit(0);
}

static void
emit_sequences(ESL_GETOPTS *go, FILE *ofp, ESL_RANDOMNESS *r, P7_HMM *hmm,
                P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA)
{
  ESL_SQ     *sq           = NULL;
  P7_TRACE   *tr           = NULL;
  P7_BG      *bg           = NULL;
  P7_PROFILE *gm           = NULL;
  int         do_fs        = esl_opt_GetBoolean(go, "--fs");
  int         do_profile   = esl_opt_GetBoolean(go, "-p");
  int         do_cons      = esl_opt_GetBoolean(go, "-c");
  int         do_fancycons = esl_opt_GetBoolean(go, "-C");
  int         do_align     = esl_opt_GetBoolean(go, "-a");
  double      fsprob       = esl_opt_GetReal(go, "--fsprob");
  int         N            = esl_opt_GetInteger(go, "-N");
  int         L            = esl_opt_GetInteger(go, "-L");
  int         mode         = p7_LOCAL;
  char        name[eslERRBUFSIZE];
  int         nseq;
  int         status;

  if      (esl_opt_GetBoolean(go, "--local"))     mode = p7_LOCAL;
  else if (esl_opt_GetBoolean(go, "--unilocal"))  mode = p7_UNILOCAL;
  else if (esl_opt_GetBoolean(go, "--glocal"))    mode = p7_GLOCAL;
  else if (esl_opt_GetBoolean(go, "--uniglocal")) mode = p7_UNIGLOCAL;

  if (p7_hmm_Validate(hmm, NULL, 0.0001) != eslOK) esl_fatal("whoops, HMM is bad!");

  if (do_cons || do_fancycons)
    {
      /* Consensus modes: exactly one deterministic output, no -N loop. */
      if (do_fancycons)
        {
          if ((sq = esl_sq_Create()) == NULL) esl_fatal("failed to allocate sequence");
          if (p7_emit_FancyConsensus(hmm, esl_opt_GetReal(go, "--minl"), esl_opt_GetReal(go, "--minu"), sq) != eslOK)
            esl_fatal("failed to create fancy consensus sequence");
          /* Only digitize (which discards -C's upper/lower case display) if we're
           * headed to --dna, where we need digital codes for codon lookup and case
           * is irrelevant. Plain protein output keeps the text-mode sq as-is. */
          if (codon_table && esl_sq_Digitize(hmm->abc, sq) != eslOK) esl_fatal("failed to digitize consensus sequence");
        }
      else
        {
          if ((sq = esl_sq_CreateDigital(hmm->abc)) == NULL) esl_fatal("failed to allocate sequence");
          if (p7_emit_SimpleConsensus(hmm, sq) != eslOK) esl_fatal("failed to create consensus sequence");
        }

      if (snprintf(name, sizeof(name), "%s-consensus", hmm->name) >= (int) sizeof(name))
        esl_fatal("model name too long\n");

      emit_one(ofp, r, sq, name, codon_table, abcDNA, do_fs, fsprob);
      esl_sq_Destroy(sq);
      return;
    }

  if (do_align)
    {
      /* Alignment mode: sample N sequences via the core model, each with
       * its own trace, then thread them through p7_tracealign_Seqs -- the
       * same alignment engine bathalign uses -- to assemble a genuine MSA.
       * Since each sequence's trace is known by construction (it's exactly
       * how the sequence was generated), no alignment inference is needed;
       * this just formats the already-known traces into a shared MSA. */
      ESL_MSA   *msa    = NULL;
      ESL_MSA   *dnamsa = NULL;
      ESL_SQ   **sqs    = NULL;
      P7_TRACE **trs    = NULL;
      int        i;

      ESL_ALLOC(sqs, sizeof(ESL_SQ *)   * N);
      ESL_ALLOC(trs, sizeof(P7_TRACE *) * N);
      for (i = 0; i < N; i++)
        {
          if ((sqs[i] = esl_sq_CreateDigital(hmm->abc)) == NULL) esl_fatal("failed to allocate seq");
          if ((trs[i] = p7_trace_Create())              == NULL) esl_fatal("failed to allocate trace");
        }

      for (i = 0; i < N; i++)
        {
          if (p7_CoreEmit(r, hmm, sqs[i], trs[i]) != eslOK) esl_fatal("Failed to emit sequence\n");
          if (snprintf(name, sizeof(name), "%s-sample%d", hmm->name, i+1) >= (int) sizeof(name))
            esl_fatal("model name too long\n");
          if (esl_sq_FormatName(sqs[i], "%s", name) != eslOK) esl_fatal("Failed to set sequence name\n");
        }

      if (p7_tracealign_Seqs(sqs, trs, N, hmm->M, p7_ALL_CONSENSUS_COLS, hmm, &msa) != eslOK)
        esl_fatal("Failed to build alignment\n");

      if (codon_table)
        {
          if (do_fs) dnamsa = dna_align_fs(r, msa, hmm->abc, codon_table, abcDNA, fsprob);
          else       dnamsa = dna_align   (r, msa, hmm->abc, codon_table, abcDNA);
          if (esl_msafile_Write(ofp, dnamsa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write alignment\n");
          esl_msa_Destroy(dnamsa);
        }
      else
        {
          if (esl_msafile_Write(ofp, msa, eslMSAFILE_STOCKHOLM) != eslOK) esl_fatal("Failed to write alignment\n");
        }

      for (i = 0; i < N; i++) { p7_trace_Destroy(trs[i]); esl_sq_Destroy(sqs[i]); }
      free(trs); free(sqs);
      esl_msa_Destroy(msa);
      return;
    }

  if ((sq = esl_sq_CreateDigital(hmm->abc)) == NULL)  esl_fatal("failed to allocate sequence");
  if ((tr = p7_trace_Create())              == NULL)  esl_fatal("failed to allocate trace");

  if (do_profile)
    {
      if ((bg = p7_bg_Create(hmm->abc))              == NULL)  esl_fatal("failed to create null model");
      if ((gm = p7_profile_Create(hmm->M, hmm->abc)) == NULL)  esl_fatal("failed to create profile");
      if (p7_ProfileConfig(hmm, bg, gm, L, mode)     != eslOK) esl_fatal("failed to configure profile");
      if (p7_bg_SetLength(bg, L)                     != eslOK) esl_fatal("failed to reconfig null model length");
      if (p7_profile_Validate(gm, NULL, 0.0001)      != eslOK) esl_fatal("whoops, profile is bad!");
    }

  for (nseq = 1; nseq <= N; nseq++)
    {
      if (do_profile) status = p7_ProfileEmit(r, hmm, gm, bg, sq, tr);
      else             status = p7_CoreEmit   (r, hmm, sq, tr);
      if (status)  esl_fatal("Failed to emit sequence\n");

      if (snprintf(name, sizeof(name), "%s-sample%d", hmm->name, nseq) >= (int) sizeof(name))
        esl_fatal("model name too long\n");

      emit_one(ofp, r, sq, name, codon_table, abcDNA, do_fs, fsprob);

      p7_trace_Reuse(tr);
      esl_sq_Reuse(sq);
    }

  esl_sq_Destroy(sq);
  p7_trace_Destroy(tr);
  if (bg) p7_bg_Destroy(bg);
  if (gm) p7_profile_Destroy(gm);
  return;

 ERROR:
  esl_fatal("allocation failure while assembling alignment\n");
}


/* emit_one()
 *
 * Given one already-emitted digital protein sequence <sq> (from any
 * mode: core/profile sampling or either consensus), write it out under
 * <name>, either as protein directly, or (if <codon_table> is non-NULL)
 * back-translated to DNA via the codon table, optionally with --fs
 * frameshift injection.
 *
 * A position holding the alphabet's "unknown" degenerate code (i.e. an
 * -c/-C consensus column with no confident call -- masked, or below
 * --minl) has no real codon to render or frameshift, so it's written
 * as a plain NNN regardless of --fs.
 */
static void
emit_one(FILE *ofp, ESL_RANDOMNESS *r, ESL_SQ *sq, const char *name,
          P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA,
          int do_fs, double fsprob)
{
  ESL_SQ  *dnasq   = NULL;
  ESL_DSQ *dna_dsq = NULL;
  ESL_DSQ  codon[3];
  int64_t  pos;
  int64_t  j;              /* next unfilled position in dna_dsq (1-offset) */
  int      nnuc;
  int      status;

  if (codon_table)
    {
      /* Back-translate the emitted amino acid sequence into DNA, one
       * uniformly-random synonymous codon per residue. With --fs, each
       * codon may instead be rendered as a 1,2,4, or 5-nucleotide
       * quasicodon (see emit_fs_codon()), so allow for the worst case
       * (every residue expanding to 5nt) when sizing the buffer. */
      ESL_ALLOC(dna_dsq, sizeof(ESL_DSQ) * (sq->n * (do_fs ? 5 : 3) + 2));
      dna_dsq[0] = eslDSQ_SENTINEL;
      j = 1;
      for (pos = 1; pos <= sq->n; pos++)
        {
          if (esl_abc_XIsUnknown(sq->abc, sq->dsq[pos]))
            {
              dna_dsq[j] = dna_dsq[j+1] = dna_dsq[j+2] = esl_abc_XGetUnknown(abcDNA);
              nnuc = 3;
            }
          else
            {
              status = p7_codontable_GetCodon(codon_table, r, sq->dsq[pos], codon);
              if (status != eslOK) esl_fatal("Failed to back-translate residue %" PRId64 " to a codon\n", pos);

              if (do_fs) nnuc = emit_fs_codon(r, codon, dna_dsq + j, fsprob);
              else
                {
                  dna_dsq[j] = codon[0]; dna_dsq[j+1] = codon[1]; dna_dsq[j+2] = codon[2];
                  nnuc = 3;
                }
            }
          j += nnuc;
        }
      dna_dsq[j] = eslDSQ_SENTINEL;

      if ((dnasq = esl_sq_CreateDigitalFrom(abcDNA, name, dna_dsq, j-1, NULL, NULL, NULL)) == NULL)
        esl_fatal("failed to create DNA sequence");

      status = esl_sqio_Write(ofp, dnasq, eslSQFILE_FASTA, FALSE);
      if (status != eslOK) esl_fatal("Failed to write sequence\n");

      esl_sq_Destroy(dnasq);
      free(dna_dsq);
    }
  else
    {
      status = esl_sq_FormatName(sq, "%s", name);
      if (status) esl_fatal("Failed to set sequence name\n");

      status = esl_sqio_Write(ofp, sq, eslSQFILE_FASTA, FALSE);
      if (status != eslOK) esl_fatal("Failed to write sequence\n");
    }

  return;

 ERROR:
  esl_fatal("allocation failure while back-translating to DNA\n");
}


/* emit_fs_codon()
 *
 * Given the 3 "true" nucleotides <codon> the codon table chose for one
 * residue, roll whether/how a frameshift renders it as a 1-, 2-, 3-, 4-,
 * or 5-nucleotide quasicodon, and write the result into <out> (must have
 * room for 5). Returns the number of nucleotides written.
 *
 * The permitted indel placements mirror enum p7p_rsc_indels in hmmer.h:
 * a deletion or insertion may fall at or before any internal codon
 * boundary, but never strictly after the codon's last position (that
 * placement is indistinguishable from, and is instead attributed to,
 * "before the next codon").
 *
 * <fsprob> sets the per-codon rate: a single-nucleotide shift (insertion or
 * deletion) occurs at rate <fsprob> in each direction, a double-nucleotide
 * shift at rate 0.5*<fsprob> in each direction, and the remainder is an
 * ordinary in-frame codon. Defaults to p7P_FSPROB (see --fsprob), the
 * fixed rate used elsewhere in BATH's own frameshift calibration.
 * --fsprob's declared range (0.001 to 0.1) keeps it well under the 1/3
 * bound past which the no-frameshift probability would go negative.
 */
static int
emit_fs_codon(ESL_RANDOMNESS *r, const ESL_DSQ *codon, ESL_DSQ *out, double fsprob)
{
  double p[5];
  int    len_idx;
  int    arr;

  p[0] = 0.5 * fsprob;      /* length 1: two deletions   (p7P_X__, p7P___X)            */
  p[1] = 1.0 * fsprob;      /* length 2: one deletion    (p7P_XX_, p7P_X_X, p7P__XX)   */
  p[2] = 1.0 - 3.0*fsprob;  /* length 3: standard codon  (p7P_XXX)                     */
  p[3] = 1.0 * fsprob;      /* length 4: one insertion   (p7P_xXXX, p7P_XxXX, p7P_XXxX)   */
  p[4] = 0.5 * fsprob;      /* length 5: two insertions  (p7P_xxXXX, p7P_XxxXX, p7P_XXxxX) */

  len_idx = esl_rnd_DChoose(r, p, 5);

  switch (len_idx)
    {
    case 0:  /* one nucleotide retained, two deleted */
      if (esl_rnd_Roll(r, 2) == 0) out[0] = codon[0];  /* p7P_X__ */
      else                         out[0] = codon[2];  /* p7P___X */
      return 1;

    case 1:  /* two nucleotides retained, one deleted */
      arr = esl_rnd_Roll(r, 3);
      if      (arr == 0) { out[0] = codon[0]; out[1] = codon[1]; }  /* p7P_XX_ */
      else if (arr == 1) { out[0] = codon[0]; out[1] = codon[2]; }  /* p7P_X_X */
      else                { out[0] = codon[1]; out[1] = codon[2]; } /* p7P__XX */
      return 2;

    case 2:  /* standard codon, no frameshift */
      out[0] = codon[0]; out[1] = codon[1]; out[2] = codon[2];      /* p7P_XXX */
      return 3;

    case 3:  /* one nucleotide inserted */
      arr = esl_rnd_Roll(r, 3);
      if (arr == 0)                                                /* p7P_xXXX */
        { out[0] = esl_rnd_Roll(r, p7P_MAXNUC); out[1] = codon[0]; out[2] = codon[1]; out[3] = codon[2]; }
      else if (arr == 1)                                           /* p7P_XxXX */
        { out[0] = codon[0]; out[1] = esl_rnd_Roll(r, p7P_MAXNUC); out[2] = codon[1]; out[3] = codon[2]; }
      else                                                         /* p7P_XXxX */
        { out[0] = codon[0]; out[1] = codon[1]; out[2] = esl_rnd_Roll(r, p7P_MAXNUC); out[3] = codon[2]; }
      return 4;

    default: /* len_idx == 4: two nucleotides inserted */
      arr = esl_rnd_Roll(r, 3);
      if (arr == 0)   /* p7P_xxXXX */
        {
          out[0] = esl_rnd_Roll(r, p7P_MAXNUC); out[1] = esl_rnd_Roll(r, p7P_MAXNUC);
          out[2] = codon[0]; out[3] = codon[1]; out[4] = codon[2];
        }
      else if (arr == 1)  /* p7P_XxxXX */
        {
          out[0] = codon[0];
          out[1] = esl_rnd_Roll(r, p7P_MAXNUC); out[2] = esl_rnd_Roll(r, p7P_MAXNUC);
          out[3] = codon[1]; out[4] = codon[2];
        }
      else  /* p7P_XXxxX */
        {
          out[0] = codon[0]; out[1] = codon[1];
          out[2] = esl_rnd_Roll(r, p7P_MAXNUC); out[3] = esl_rnd_Roll(r, p7P_MAXNUC);
          out[4] = codon[2];
        }
      return 5;
    }
}


/* dna_align()
 *
 * Given an amino acid alignment <msa> (text-mode, from p7_tracealign_Seqs
 * called *without* p7_DIGITIZE), build the DNA equivalent: every column
 * becomes a fixed-width 3nt group, a residue back-translated via
 * <codon_table> into a random synonymous codon, a gap copied through
 * verbatim (three copies of whatever character -- '.' or '-' -- the
 * engine used at that position).
 *
 * <msa> has to be text-mode, not digital: a digitized alignment only
 * carries a single gap code (Easel treats '.' and '-' as equivalent on
 * input), so the distinction this function preserves is already gone
 * once an alignment has been digitized.
 *
 * The #=GC RF line (if present) is expanded the same way as everything
 * else, each character repeated 3x to stay aligned under its column
 * group.
 *
 * Handles the no-frameshift case only, where every column group is
 * exactly 3nt wide; see dna_align_fs() for --fs, where a column group's
 * width varies with the widest quasicodon rendering used there.
 */
static ESL_MSA *
dna_align(ESL_RANDOMNESS *r, const ESL_MSA *msa, const ESL_ALPHABET *abcAmino,
          P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA)
{
  ESL_MSA *dnamsa  = NULL;
  ESL_DSQ  codon[3];
  ESL_DSQ  aa;
  int64_t  alen    = msa->alen;
  int64_t  dnalen  = alen * 3;   /* every column group is exactly 3nt wide; see above */
  int      idx;
  int64_t  apos, dpos;
  char     c;
  int      status;

  if ((dnamsa = esl_msa_Create(msa->nseq, dnalen)) == NULL) esl_fatal("failed to allocate DNA alignment");

  for (idx = 0; idx < msa->nseq; idx++)
    {
      dpos = 0;
      for (apos = 0; apos < alen; apos++)
        {
          c = msa->aseq[idx][apos];

          if (c == '.' || c == '-')
            {
              dnamsa->aseq[idx][dpos] = dnamsa->aseq[idx][dpos+1] = dnamsa->aseq[idx][dpos+2] = c;
            }
          else
            {
              aa = abcAmino->inmap[(int) c];
              status = p7_codontable_GetCodon(codon_table, r, aa, codon);
              if (status != eslOK) esl_fatal("Failed to back-translate alignment column %" PRId64 " to a codon\n", apos+1);
              dnamsa->aseq[idx][dpos]   = abcDNA->sym[codon[0]];
              dnamsa->aseq[idx][dpos+1] = abcDNA->sym[codon[1]];
              dnamsa->aseq[idx][dpos+2] = abcDNA->sym[codon[2]];
            }
          dpos += 3;
        }

      if (esl_msa_SetSeqName(dnamsa, idx, msa->sqname[idx], -1) != eslOK) esl_fatal("failed to set DNA alignment sequence name");
    }

  if (msa->rf)
    {
      ESL_ALLOC(dnamsa->rf, sizeof(char) * (dnalen+1));
      dpos = 0;
      for (apos = 0; apos < alen; apos++)
        {
          dnamsa->rf[dpos] = dnamsa->rf[dpos+1] = dnamsa->rf[dpos+2] = msa->rf[apos];
          dpos += 3;
        }
      dnamsa->rf[dnalen] = '\0';
    }

  return dnamsa;

 ERROR:
  esl_fatal("allocation failure while building DNA alignment\n");
  return NULL;
}


/* emit_fs_slots()
 *
 * Same roll as emit_fs_codon() (identical probabilities, identical
 * enum p7p_rsc_indels patterns -- see that function's header for the
 * derivation), but instead of flattening the result into a contiguous
 * nucleotide run, records which of the 6 slots (insert-before-1, pos1,
 * insert-1-2, pos2, insert-2-3, pos3) the roll actually used, for
 * dna_align_fs() to place into a shared column frame alongside other
 * sequences' rolls at the same alignment column.
 */
static void
emit_fs_slots(ESL_RANDOMNESS *r, const ESL_DSQ *codon, double fsprob, FS_SLOT *slot)
{
  double p[5];
  int    len_idx, arr;

  p[0] = 0.5 * fsprob;
  p[1] = 1.0 * fsprob;
  p[2] = 1.0 - 3.0*fsprob;
  p[3] = 1.0 * fsprob;
  p[4] = 0.5 * fsprob;

  slot->has      = TRUE;
  slot->n_before = slot->n_mid12 = slot->n_mid23 = 0;
  slot->has1     = slot->has2    = slot->has3    = FALSE;

  len_idx = esl_rnd_DChoose(r, p, 5);

  switch (len_idx)
    {
    case 0:  /* one nucleotide retained, two deleted */
      if (esl_rnd_Roll(r, 2) == 0) { slot->has1 = TRUE; slot->pos1 = codon[0]; }  /* p7P_X__ */
      else                         { slot->has3 = TRUE; slot->pos3 = codon[2]; }  /* p7P___X */
      break;

    case 1:  /* two nucleotides retained, one deleted */
      arr = esl_rnd_Roll(r, 3);
      if      (arr == 0) { slot->has1 = TRUE; slot->pos1 = codon[0]; slot->has2 = TRUE; slot->pos2 = codon[1]; }  /* XX_ */
      else if (arr == 1) { slot->has1 = TRUE; slot->pos1 = codon[0]; slot->has3 = TRUE; slot->pos3 = codon[2]; }  /* X_X */
      else                { slot->has2 = TRUE; slot->pos2 = codon[1]; slot->has3 = TRUE; slot->pos3 = codon[2]; } /* _XX */
      break;

    case 2:  /* standard codon, no frameshift */
      slot->has1 = TRUE; slot->pos1 = codon[0];
      slot->has2 = TRUE; slot->pos2 = codon[1];
      slot->has3 = TRUE; slot->pos3 = codon[2];
      break;

    case 3:  /* one nucleotide inserted */
      slot->has1 = TRUE; slot->pos1 = codon[0];
      slot->has2 = TRUE; slot->pos2 = codon[1];
      slot->has3 = TRUE; slot->pos3 = codon[2];
      arr = esl_rnd_Roll(r, 3);
      if (arr == 0)      { slot->n_before = 1; slot->before[0] = esl_rnd_Roll(r, p7P_MAXNUC); }  /* xXXX */
      else if (arr == 1) { slot->n_mid12  = 1; slot->mid12[0]  = esl_rnd_Roll(r, p7P_MAXNUC); }  /* XxXX */
      else                { slot->n_mid23  = 1; slot->mid23[0]  = esl_rnd_Roll(r, p7P_MAXNUC); }  /* XXxX */
      break;

    default: /* len_idx == 4: two nucleotides inserted */
      slot->has1 = TRUE; slot->pos1 = codon[0];
      slot->has2 = TRUE; slot->pos2 = codon[1];
      slot->has3 = TRUE; slot->pos3 = codon[2];
      arr = esl_rnd_Roll(r, 3);
      if (arr == 0)      { slot->n_before = 2; slot->before[0] = esl_rnd_Roll(r, p7P_MAXNUC); slot->before[1] = esl_rnd_Roll(r, p7P_MAXNUC); }  /* xxXXX */
      else if (arr == 1) { slot->n_mid12  = 2; slot->mid12[0]  = esl_rnd_Roll(r, p7P_MAXNUC); slot->mid12[1]  = esl_rnd_Roll(r, p7P_MAXNUC); }  /* XxxXX */
      else                { slot->n_mid23  = 2; slot->mid23[0]  = esl_rnd_Roll(r, p7P_MAXNUC); slot->mid23[1]  = esl_rnd_Roll(r, p7P_MAXNUC); }  /* XXxxX */
      break;
    }
}


/* dna_align_fs()
 *
 * The --fs sibling of dna_align(): builds the DNA alignment when
 * quasicodon renderings can vary in width per sequence. Every amino
 * alignment column becomes a "column group" of 6 slots --
 *   [insert-before-1] pos1 [insert-1-2] pos2 [insert-2-3] pos3
 * -- sized to the widest usage any sequence actually rolled at that
 * position (0, 1, or 2 for each insert sub-slot; the three pos slots
 * are always reserved, 1 wide each). A sequence with no residue at
 * that amino column is gapped across the whole group using whichever
 * character ('.' or '-') the amino alignment used there. A sequence
 * whose own roll used less of an insert sub-slot than the group's max
 * is right-justified within it (real nucleotides adjacent to the fixed
 * position slot that follows), gap-padded on the far side. A position
 * slot (pos1/pos2/pos3) the roll didn't fill (an in-codon deletion) is '-'.
 *
 * Three passes over the alignment: (1) roll every residue's slot
 * breakdown, (2) find each column's required group width from the
 * widest usage there, (3) render every row into the now-known frame.
 */
static ESL_MSA *
dna_align_fs(ESL_RANDOMNESS *r, const ESL_MSA *msa, const ESL_ALPHABET *abcAmino,
             P7_CODONTABLE *codon_table, const ESL_ALPHABET *abcDNA, double fsprob)
{
  ESL_MSA  *dnamsa  = NULL;
  FS_SLOT **grid    = NULL;   /* [idx][apos] */
  int      *mbefore = NULL;   /* [apos]: this column's insert-before-1 sub-slot width */
  int      *mmid12  = NULL;   /* [apos]: this column's insert-1-2 sub-slot width      */
  int      *mmid23  = NULL;   /* [apos]: this column's insert-2-3 sub-slot width      */
  int64_t  *gstart  = NULL;   /* [apos]: this column's group's starting DNA position  */
  int64_t   alen    = msa->alen;
  int64_t   dnalen;
  int       idx;
  int64_t   apos, dpos, k, pad;
  ESL_DSQ   codon[3];
  ESL_DSQ   aa;
  char      c;
  int       status;

  ESL_ALLOC(grid, sizeof(FS_SLOT *) * msa->nseq);
  for (idx = 0; idx < msa->nseq; idx++) { grid[idx] = NULL; }
  for (idx = 0; idx < msa->nseq; idx++) { ESL_ALLOC(grid[idx], sizeof(FS_SLOT) * alen); }
  ESL_ALLOC(mbefore, sizeof(int)     * alen);
  ESL_ALLOC(mmid12,  sizeof(int)     * alen);
  ESL_ALLOC(mmid23,  sizeof(int)     * alen);
  ESL_ALLOC(gstart,  sizeof(int64_t) * alen);

  /* Pass 1: roll every residue's quasicodon slot breakdown. */
  for (idx = 0; idx < msa->nseq; idx++)
    for (apos = 0; apos < alen; apos++)
      {
        c = msa->aseq[idx][apos];
        if (c == '.' || c == '-') { grid[idx][apos].has = FALSE; continue; }

        aa = abcAmino->inmap[(int) c];
        status = p7_codontable_GetCodon(codon_table, r, aa, codon);
        if (status != eslOK) esl_fatal("Failed to back-translate alignment column %" PRId64 " to a codon\n", apos+1);
        emit_fs_slots(r, codon, fsprob, &(grid[idx][apos]));
      }

  /* Pass 2: each column's group width = widest insert usage there, plus the 3 fixed position slots. */
  dnalen = 0;
  for (apos = 0; apos < alen; apos++)
    {
      mbefore[apos] = mmid12[apos] = mmid23[apos] = 0;
      for (idx = 0; idx < msa->nseq; idx++)
        if (grid[idx][apos].has)
          {
            if (grid[idx][apos].n_before > mbefore[apos]) mbefore[apos] = grid[idx][apos].n_before;
            if (grid[idx][apos].n_mid12  > mmid12[apos])  mmid12[apos]  = grid[idx][apos].n_mid12;
            if (grid[idx][apos].n_mid23  > mmid23[apos])  mmid23[apos]  = grid[idx][apos].n_mid23;
          }
      gstart[apos] = dnalen;
      dnalen += mbefore[apos] + 1 + mmid12[apos] + 1 + mmid23[apos] + 1;
    }

  if ((dnamsa = esl_msa_Create(msa->nseq, dnalen)) == NULL) esl_fatal("failed to allocate DNA alignment");

  /* Pass 3: render every row into the now-known column frame. */
  for (idx = 0; idx < msa->nseq; idx++)
    {
      for (apos = 0; apos < alen; apos++)
        {
          dpos = gstart[apos];

          if (! grid[idx][apos].has)
            {
              c = msa->aseq[idx][apos];   /* '.' or '-', copied across the whole group */
              for (k = 0; k < mbefore[apos]+1+mmid12[apos]+1+mmid23[apos]+1; k++) dnamsa->aseq[idx][dpos+k] = c;
              continue;
            }

          /* insert-before-1 sub-slot: right-justified against pos1, gap-padded on the left */
          pad = mbefore[apos] - grid[idx][apos].n_before;
          for (k = 0; k < pad; k++)                        dnamsa->aseq[idx][dpos+k]     = '-';
          for (k = 0; k < grid[idx][apos].n_before; k++)    dnamsa->aseq[idx][dpos+pad+k] = abcDNA->sym[grid[idx][apos].before[k]];
          dpos += mbefore[apos];

          dnamsa->aseq[idx][dpos] = grid[idx][apos].has1 ? abcDNA->sym[grid[idx][apos].pos1] : '-';
          dpos += 1;

          /* insert-1-2 sub-slot: right-justified against pos2 */
          pad = mmid12[apos] - grid[idx][apos].n_mid12;
          for (k = 0; k < pad; k++)                        dnamsa->aseq[idx][dpos+k]     = '-';
          for (k = 0; k < grid[idx][apos].n_mid12; k++)     dnamsa->aseq[idx][dpos+pad+k] = abcDNA->sym[grid[idx][apos].mid12[k]];
          dpos += mmid12[apos];

          dnamsa->aseq[idx][dpos] = grid[idx][apos].has2 ? abcDNA->sym[grid[idx][apos].pos2] : '-';
          dpos += 1;

          /* insert-2-3 sub-slot: right-justified against pos3 */
          pad = mmid23[apos] - grid[idx][apos].n_mid23;
          for (k = 0; k < pad; k++)                        dnamsa->aseq[idx][dpos+k]     = '-';
          for (k = 0; k < grid[idx][apos].n_mid23; k++)     dnamsa->aseq[idx][dpos+pad+k] = abcDNA->sym[grid[idx][apos].mid23[k]];
          dpos += mmid23[apos];

          dnamsa->aseq[idx][dpos] = grid[idx][apos].has3 ? abcDNA->sym[grid[idx][apos].pos3] : '-';
        }

      if (esl_msa_SetSeqName(dnamsa, idx, msa->sqname[idx], -1) != eslOK) esl_fatal("failed to set DNA alignment sequence name");
    }

  /* RF line: each original character repeated to fill its column's (now variable) group width. */
  if (msa->rf)
    {
      ESL_ALLOC(dnamsa->rf, sizeof(char) * (dnalen+1));
      for (apos = 0; apos < alen; apos++)
        {
          int width = mbefore[apos]+1+mmid12[apos]+1+mmid23[apos]+1;
          for (k = 0; k < width; k++) dnamsa->rf[gstart[apos]+k] = msa->rf[apos];
        }
      dnamsa->rf[dnalen] = '\0';
    }

  for (idx = 0; idx < msa->nseq; idx++) free(grid[idx]);
  free(grid); free(mbefore); free(mmid12); free(mmid23); free(gstart);

  return dnamsa;

 ERROR:
  esl_fatal("allocation failure while building frameshifted DNA alignment\n");
  return NULL;
}
