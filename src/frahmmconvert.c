/* frahmmconvert: converting HMMER formated profile HMM files to FraHMMER format.
 */
#include "p7_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_getopts.h"

#include "hmmer.h"

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,       NULL,    NULL, "show brief help on version and usage",                             1 },
  { "--fs",     eslARG_REAL,  "0.01",NULL, "0.001<=x<=0.05", NULL, NULL, NULL,  "set the frameshift probabilty",                 1 },
  { "-c",        eslARG_INT,      "1", NULL,   NULL,      NULL,        NULL,  NULL,  "use alt genetic code of NCBI transl table <n> ",        1 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[]  = "[-options] <hmmfile_out> <hmmfile_in>";
static char banner[] = "convert HMMER format pHMM to FraHMMER format";


int 
main(int argc, char **argv)
{
  ESL_GETOPTS   *go      = p7_CreateDefaultApp(options, 2, argc, argv, banner, usage);
  ESL_ALPHABET  *abc     = NULL;
  char          *hmmfile_out = esl_opt_GetArg(go, 1);
  char          *hmmfile_in  = esl_opt_GetArg(go, 2);
  P7_HMMFILE    *hfp     = NULL;
  P7_HMM        *hmm     = NULL;
  FILE          *ofp     = NULL;
  int            fmtcode = 7;	/* 7 = FraHMMER format */
  int            status;
  char           errbuf[eslERRBUFSIZE];
  float          fs;
  int            ct;
  P7_BG          *bg   = NULL;
  ESL_RANDOMNESS *r  = NULL;
  P7_FS_PROFILE  *gm_fs  = NULL;
  double          tau_fs;

  if (esl_opt_GetBoolean(go, "-h") == TRUE)
  {
     // p7_banner(stdout, argv[0], banner); This is HMMER banner - need to format FraHMMER version
     esl_usage(stdout, argv[0], usage);

     if (puts("\nBasic options:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);

     if (puts("\nAvailable NCBI genetic code tables (for -c <id>):")        < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_gencode_DumpAltCodeTable(stdout);
      exit(0);
  }
  
  status = p7_hmmfile_OpenE(hmmfile_in, NULL, &hfp, errbuf);
  if      (status == eslENOTFOUND) p7_Fail("File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile_in, errbuf);
  else if (status == eslEFORMAT)   p7_Fail("File format problem in trying to open HMM file %s.\n%s\n",                hmmfile_in, errbuf);
  else if (status != eslOK)        p7_Fail("Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile_in, errbuf);  

  ofp = fopen(hmmfile_out, "w");
  if (ofp == NULL) p7_Fail("Failed to open HMM file %s for writing", hmmfile_out);

  while ((status = p7_hmmfile_Read(hfp, &abc, &hmm)) == eslOK)
    {
      if(hmm->abc->type != eslAMINO) p7_Fail("Invalid alphabet type in the pHMM input file %s. Expect Amino Acid\n", hmmfile_in); 

      fs = esl_opt_GetReal(go, "--fs");
      ct = esl_opt_GetInteger(go, "-c");
 
      if(fs != hmm->fs || ct != hmm->ct)
        {
          hmm->fs = fs;
          hmm->ct = ct;
          r = esl_randomness_CreateFast(42);
          gm_fs = p7_profile_fs_Create (hmm->M, hmm->abc);
          bg = p7_bg_Create(hmm->abc);

          p7_fs_Tau(r, gm_fs, hmm, bg, 100, 200, hmm->fs, hmm->evparam[p7_FLAMBDA], 0.04, &tau_fs);
          hmm->evparam[p7_FTAUFS] = tau_fs;
        }
      if(hmm->max_length == -1)
	{
		p7_Builder_MaxLength(hmm, p7_DEFAULT_WINDOW_BETA);	
	} 

      p7_hmmfile_WriteASCII (ofp, fmtcode, hmm);

      p7_hmm_Destroy(hmm);
    }
  if      (status == eslEFORMAT)   p7_Fail("bad file format in HMM file %s",             hmmfile_in);
  else if (status == eslEINCOMPAT) p7_Fail("HMM file %s contains different alphabets",   hmmfile_in);
  else if (status != eslEOF)       p7_Fail("Unexpected error in reading HMMs from %s",   hmmfile_in);
  
  if(bg != NULL)    p7_bg_Destroy(bg);
  if(gm_fs != NULL) p7_profile_fs_Destroy(gm_fs);
  if(r != NULL)     esl_randomness_Destroy(r);
  if(ofp) fclose(ofp);
  if(hfp) p7_hmmfile_Close(hfp);
  if(abc) esl_alphabet_Destroy(abc);
  if(go) esl_getopts_Destroy(go);
  return 0;

 ERROR:
  if(bg != NULL)    p7_bg_Destroy(bg);
  if(gm_fs != NULL) p7_profile_fs_Destroy(gm_fs);
  if(r != NULL)     esl_randomness_Destroy(r);
  if(ofp) fclose(ofp);
  if(hfp) p7_hmmfile_Close(hfp);
  if(abc) esl_alphabet_Destroy(abc);
  if(go) esl_getopts_Destroy(go);
  exit(status);
}


