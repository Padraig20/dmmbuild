/* dmmbuild: building profile HMMs from multiple sequence alignments for DUMMER.
 *
 * Contents:
 *   1. Definitions of structures used in the program.
 *   2. Process command line arguments.
 *   3. Helper functions.
 *   4. Serial loop, build HMMs for each MSA.
 *   5. Master that builds HMMs for each MSA.
 *   6. Main program.
 */

#include "dummer.h"

/*****************************************************************
 * 1. Definitions of structures used in the program
 *****************************************************************/

/* struct cfg_s : "Global" application configuration shared by all threads/processes
 * 
 * This structure is passed to routines within main.c, as a means of semi-encapsulation
 * of shared data amongst different parallel processes (threads or MPI processes).
 */
struct cfg_s {
  FILE         *ofp;		/* output file (default is stdout) */

  char         *alifile;	/* name of the alignment file we're building HMMs from  */
  int           fmt;		/* format code for alifile */
  ESL_MSAFILE  *afp;            /* open alifile  */
  ESL_ALPHABET *abc;		/* digital alphabet */

  char         *hmmName;        /* hmm file name supplied from -n          */
  char         *hmmfile;        /* file to write HMM to                    */
  FILE         *hmmfp;          /* HMM output file handle                  */

  int           nali;		/* which # alignment this is in file (only valid in serial mode)   */
  int           nnamed;		/* number of alignments that had their own names */

  int           do_mpi;		/* TRUE if we're doing MPI parallelization */
  int           nproc;		/* how many MPI processes, total */
  int           my_rank;	/* who am I, in 0..nproc-1 */
  int           do_stall;	/* TRUE to stall the program until gdb attaches */
};

#define ALPHOPTS "--amino,--dna,--rna"                         /* Exclusive options for alphabet choice */
#define CONOPTS "--fast"                                /* Exclusive options for model construction                    */
#define EFFOPTS "--eent,--eentexp,--eclust,--eset,--enone"               /* Exclusive options for effective sequence number calculation */
#define WGTOPTS "--wgsc,--wblosum,--wpb,--wnone,--wgiven"      /* Exclusive options for relative weighting                    */

typedef struct {
  F4_BG	           *bg;
  F4_BUILDER       *bld;
} WORKER_INFO;

static char usage[]  = "[-options] <dmmfile_out> <msafile>";
static char banner[] = "profile HMM construction from multiple sequence alignments for DUMMER";

static ESL_OPTIONS options[] = {
  /* name           type      default  env  range     toggles      reqs   incomp  help   docgroup*/
  { "-h",        eslARG_NONE,   FALSE, NULL, NULL,      NULL,      NULL,    NULL, "show brief help on version and usage",                  1 },
  { "-n",        eslARG_STRING,  NULL, NULL, NULL,      NULL,      NULL,    NULL, "name the HMM <s>",                                      1 },
  { "-o",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "direct summary output to file <f>, not stdout",         1 },
  { "-O",        eslARG_OUTFILE,FALSE, NULL, NULL,      NULL,      NULL,    NULL, "resave annotated, possibly modified MSA to file <f>",   1 },
 
  /* Selecting the alphabet rather than autoguessing it */
  { "--amino",   eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is protein sequence data",              2 },
  { "--dna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is DNA sequence data",                  2 },
  { "--rna",     eslARG_NONE,   FALSE, NULL, NULL,   ALPHOPTS,    NULL,     NULL, "input alignment is RNA sequence data",                  2 },

  /* Alternate model construction strategies */
  { "--fast",    eslARG_NONE,"default",NULL, NULL,    CONOPTS,    NULL,     NULL, "assign cols w/ >= symfrac residues as consensus",       3 },
  //{ "--hand",    eslARG_NONE,   FALSE, NULL, NULL,    CONOPTS,    NULL,     NULL, "manual construction (requires reference annotation)",   3 },
  { "--symfrac", eslARG_REAL,   "0.5", NULL, "0<=x<=1", NULL,   "--fast",   NULL, "sets sym fraction controlling --fast construction",     3 },
  { "--fragthresh",eslARG_REAL, "0.5", NULL, "0<=x<=1", NULL,     NULL,     NULL, "if L <= x*alen, tag sequence as a fragment",            3 },

  /* Alternate relative sequence weighting strategies */
  /* --wme not implemented in HMMER3 yet */
  { "--wpb",     eslARG_NONE,"default",NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff position-based weights",                      4 },
  { "--wgsc",    eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Gerstein/Sonnhammer/Chothia tree weights",             4 },
  { "--wblosum", eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "Henikoff simple filter weights",                       4 },
  { "--wnone",   eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "don't do any relative weighting; set all to 1",        4 },
  { "--wgiven",  eslARG_NONE,   NULL,  NULL, NULL,    WGTOPTS,    NULL,      NULL, "use weights as given in MSA file",                     4 },
  { "--wid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--wblosum",   NULL, "for --wblosum: set identity cutoff",                   4 },

  /* Alternative effective sequence weighting strategies */
  { "--eent",    eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to achieve relative entropy target",  5 },
  { "--eclust",  eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "eff seq # is # of single linkage clusters",            5 },
  { "--enone",   eslARG_NONE,  FALSE,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "no effective seq # weighting: just use nseq",          5 },
  { "--eset",    eslARG_REAL,   NULL,  NULL, NULL,    EFFOPTS,    NULL,      NULL, "set eff seq # for all models to <x>",                  5 },
  { "--ere",     eslARG_REAL,   NULL,  NULL,"x>0",       NULL,    NULL,     NULL, "for --eent: set minimum rel entropy/position to <x>",  5 },
  { "--esigma",  eslARG_REAL, "45.0",  NULL,"x>0",       NULL,    NULL,     NULL, "for --eent: set sigma param to <x>",                   5 },
  { "--eid",     eslARG_REAL, "0.62",  NULL,"0<=x<=1",   NULL,"--eclust",    NULL, "for --eclust: set fractional identity cutoff to <x>",  5 },

  /* Alternative prior strategies */
  { "--pnone",   eslARG_NONE,  FALSE,  NULL, NULL,       NULL,  NULL,"--plaplace", "don't use any prior; parameters are frequencies",      9 },
  { "--plaplace",eslARG_NONE,  FALSE,  NULL, NULL,       NULL,  NULL,   "--pnone", "use a Laplace +1 prior",                               9 },

  /* Single sequence methods */
  //{ "--singlemx", eslARG_NONE,   FALSE, NULL,   NULL,   NULL,  NULL,           "",   "use substitution score matrix for single-sequence inputs",     10 },
  //{ "--mx",     eslARG_STRING, "BLOSUM62", NULL, NULL,   NULL, NULL,   "--mxfile",   "substitution score matrix (built-in matrices, with --singlemx)", 10 },
  //{ "--mxfile", eslARG_INFILE,     NULL, NULL,   NULL,   NULL, NULL,       "--mx",   "read substitution score matrix from file <f> (with --singlemx)", 10 },
  //{ "--popen",    eslARG_REAL,  NULL,  NULL,"0<=x<0.5",NULL, NULL,           "",   "force gap open prob. (w/ --singlemx, aa default 0.02, nt 0.031)",  10 },
  //{ "--pextend",  eslARG_REAL,  NULL,  NULL, "0<=x<1", NULL, NULL,           "",   "force gap extend prob. (w/ --singlemx, aa default 0.4, nt 0.75)",  10 },

  /* Control of E-value calibration */
  { "--EmL",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for MSV Gumbel mu fit",            6 },   
  { "--EmN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for MSV Gumbel mu fit",            6 },   
  { "--EvL",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Viterbi Gumbel mu fit",        6 },   
  { "--EvN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Viterbi Gumbel mu fit",        6 },   
  { "--EfL",     eslARG_INT,    "100", NULL,"n>0",       NULL,    NULL,      NULL, "length of sequences for Forward exp tail tau fit",     6 },   
  { "--EfN",     eslARG_INT,    "200", NULL,"n>0",       NULL,    NULL,      NULL, "number of sequences for Forward exp tail tau fit",     6 },   
  { "--Eft",     eslARG_REAL,  "0.04", NULL,"0<x<1",     NULL,    NULL,      NULL, "tail mass for Forward exponential tail tau fit",       6 },   

  /* Other options */
  { "--stall",   eslARG_NONE,       FALSE, NULL, NULL,    NULL,     NULL,    NULL, "arrest after start: for attaching debugger to process", 8 },
  { "--informat", eslARG_STRING,     NULL, NULL, NULL,    NULL,     NULL,    NULL, "assert input alifile is in format <s> (no autodetect)", 8 },
  { "--seed",     eslARG_INT,        "42", NULL, "n>=0",  NULL,     NULL,    NULL, "set RNG seed to <n> (if 0: one-time arbitrary seed)",   8 },
  { "--w_beta",   eslARG_REAL,       NULL, NULL, NULL,    NULL,     NULL,    NULL, "tail mass at which window length is determined",        8 },
  { "--w_length", eslARG_INT,        NULL, NULL, NULL,    NULL,     NULL,    NULL, "window length ",                                        8 },
  { "--maxinsertlen",  eslARG_INT,   NULL, NULL, "n>=5",  NULL,     NULL,    NULL, "pretend all inserts are length <= <n>",   8 },

  /*Expert-only option, hidden from view. Likely to be removed in the future.
    This is an experimental alternative method for weighting sequence counts.
    It produces non-uniform scaling across columns -- columns with higher
    counts are reduced to a greater extent than columns with low counts. It
    has proven useful in reducing alignment overextension with Dfam models -
    these are often based on MSAs reconstructed from heavily-fragmented
    sequences, and thus may have wildly different per-column coverage levels.
    Uniform scaling causes low-coverage regions to be scaled down to
    essentially no (weighted) observations, leading to low local relative
    entropy, which produces high overextension.

    However, we don't recommend using this method in general, as it
    may show odd behavior (or fail to even work) in the case of low
    overall observed counts (e.g. an alignemnt of 2 or 4 sequences)
    */
    { "--eentexp", eslARG_NONE,"default",NULL, NULL,    EFFOPTS,    NULL,      NULL, "adjust eff seq # to reach rel. ent. target using exp scaling",  99 },
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};

/*****************************************************************
 * 2. Process command line arguments.
 *****************************************************************/

static int
process_commandline(int argc, char **argv, ESL_GETOPTS **ret_go, char **ret_dmmfile, char **ret_alifile)
{
  ESL_GETOPTS *go = esl_getopts_Create(options);
  int          status;

  if (esl_opt_ProcessEnvironment(go)         != eslOK) { if (printf("Failed to process environment:\n%s\n", go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_ProcessCmdline(go, argc, argv) != eslOK) { if (printf("Failed to parse command line:\n%s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (esl_opt_VerifyConfig(go)               != eslOK) { if (printf("Failed to parse command line:\n%s\n",  go->errbuf) < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  /* help format: */
  if (esl_opt_GetBoolean(go, "-h") == TRUE) 
    {
      f4_banner(stdout, argv[0], banner);
      esl_usage(stdout, argv[0], usage);

      if (puts("\nBasic options:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 1, 2, 80);

      if (puts("\nOptions for selecting alphabet rather than guessing it:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 2, 2, 80);

      if (puts("\nAlternative model construction strategies:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 3, 2, 80);

      if (puts("\nAlternative relative sequence weighting strategies:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 4, 2, 80);

      if (puts("\nAlternative effective sequence weighting strategies:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 5, 2, 80);

      if (puts("\nAlternative prior strategies:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 9, 2, 80);

      if (puts("\nHandling single sequence inputs:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 10, 2, 80);

      if (puts("\nControl of E-value calibration:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 6, 2, 80);

      if (puts("\nOther options:") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed");
      esl_opt_DisplayHelp(stdout, go, 8, 2, 80);
      exit(0);
    }

  if (esl_opt_ArgNumber(go)                  != 2)    { if (puts("Incorrect number of command line arguments.")          < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_dmmfile = esl_opt_GetArg(go, 1)) == NULL) { if (puts("Failed to get <dmmfile_out> argument on command line") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if ((*ret_alifile = esl_opt_GetArg(go, 2)) == NULL) { if (puts("Failed to get <msafile> argument on command line")     < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }

  if (strcmp(*ret_dmmfile, "-") == 0) 
    { if (puts("Can't write <dmmfile_out> to stdout: don't use '-'")         < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  if (strcmp(*ret_alifile, "-") == 0 && ! esl_opt_IsOn(go, "--informat"))
    { if (puts("Must specify --informat to read <alifile> from stdin ('-')") < 0) ESL_XEXCEPTION_SYS(eslEWRITE, "write failed"); goto FAILURE; }
  
  *ret_go = go;
  return eslOK;
  
 FAILURE:  /* all errors handled here are user errors, so be polite.  */
  esl_usage(stdout, argv[0], usage);
  puts("\nwhere basic options are:");
  esl_opt_DisplayHelp(stdout, go, 1, 2, 80);
  printf("\nTo see more help on other available options, do:\n  %s -h\n\n", argv[0]);
  esl_getopts_Destroy(go);
  exit(1);  

 ERROR:
  if (go) esl_getopts_Destroy(go);
  exit(status);
}

/*****************************************************************
 * 3. Helper functions.
 *****************************************************************/

static int
output_header(const ESL_GETOPTS *go, const struct cfg_s *cfg)
{
  if (cfg->my_rank > 0)  return eslOK;

  f4_banner(cfg->ofp, go->argv[0], banner);
  
  if (fprintf(cfg->ofp, "# input alignment file:             %s\n", cfg->alifile) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (fprintf(cfg->ofp, "# output HMM file:                  %s\n", cfg->hmmfile) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "-n")           && fprintf(cfg->ofp, "# name (the single) HMM:            %s\n",        esl_opt_GetString(go, "-n"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-o")           && fprintf(cfg->ofp, "# output directed to file:          %s\n",        esl_opt_GetString(go, "-o"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "-O")           && fprintf(cfg->ofp, "# processed alignment resaved to:   %s\n",        esl_opt_GetString(go, "-O"))         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--amino")      && fprintf(cfg->ofp, "# input alignment is asserted as:   protein\n")                                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--dna")        && fprintf(cfg->ofp, "# input alignment is asserted as:   DNA\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--rna")        && fprintf(cfg->ofp, "# input alignment is asserted as:   RNA\n")                                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fast")       && fprintf(cfg->ofp, "# model architecture construction:  fast/heuristic\n")                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  //if (esl_opt_IsUsed(go, "--hand")       && fprintf(cfg->ofp, "# model architecture construction:  hand-specified by RF annotation\n")                < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--symfrac")    && fprintf(cfg->ofp, "# sym fraction for model structure: %.3f\n",      esl_opt_GetReal(go, "--symfrac"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--fragthresh") && fprintf(cfg->ofp, "# seq called frag if L <= x*alen:   %.3f\n",      esl_opt_GetReal(go, "--fragthresh")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--wpb")        && fprintf(cfg->ofp, "# relative weighting scheme:        Henikoff PB\n")                                    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--wgsc")       && fprintf(cfg->ofp, "# relative weighting scheme:        G/S/C\n")                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--wblosum")    && fprintf(cfg->ofp, "# relative weighting scheme:        BLOSUM filter\n")                                  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--wnone")      && fprintf(cfg->ofp, "# relative weighting scheme:        none\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--wid")        && fprintf(cfg->ofp, "# frac id cutoff for BLOSUM wgts:   %f\n",        esl_opt_GetReal(go, "--wid"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--eent")       && fprintf(cfg->ofp, "# effective seq number scheme:      entropy weighting\n")                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--eentexp")    && fprintf(cfg->ofp, "# effective seq number scheme:      entropy weighting using exponent-based scaling\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--eclust")     && fprintf(cfg->ofp, "# effective seq number scheme:      single linkage clusters\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--enone")      && fprintf(cfg->ofp, "# effective seq number scheme:      none\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--eset")       && fprintf(cfg->ofp, "# effective seq number:             set to %f\n", esl_opt_GetReal(go, "--eset"))       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--ere")        && fprintf(cfg->ofp, "# minimum rel entropy target:       %f bits\n",   esl_opt_GetReal(go, "--ere"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--esigma")     && fprintf(cfg->ofp, "# entropy target sigma parameter:   %f bits\n",   esl_opt_GetReal(go, "--esigma"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--eid")        && fprintf(cfg->ofp, "# frac id cutoff for --eclust:      %f\n",        esl_opt_GetReal(go, "--eid"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--pnone")      && fprintf(cfg->ofp, "# prior scheme:                     none\n")                                           < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--plaplace")   && fprintf(cfg->ofp, "# prior scheme:                     Laplace +1\n")                                     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--EmL")        && fprintf(cfg->ofp, "# seq length for MSV Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EmL"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--EmN")        && fprintf(cfg->ofp, "# seq number for MSV Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EmN"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--EvL")        && fprintf(cfg->ofp, "# seq length for Vit Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EvL"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--EvN")        && fprintf(cfg->ofp, "# seq number for Vit Gumbel mu fit: %d\n",        esl_opt_GetInteger(go, "--EvN"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--EfL")        && fprintf(cfg->ofp, "# seq length for Fwd exp tau fit:   %d\n",        esl_opt_GetInteger(go, "--EfL"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--EfN")        && fprintf(cfg->ofp, "# seq number for Fwd exp tau fit:   %d\n",        esl_opt_GetInteger(go, "--EfN"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--Eft")        && fprintf(cfg->ofp, "# tail mass for Fwd exp tau fit:    %f\n",        esl_opt_GetReal(go, "--Eft"))        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  //if (esl_opt_IsUsed(go, "--singlemx")   && fprintf(cfg->ofp, "# use score matrix for 1-seq MSAs:  on\n")                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  //if (esl_opt_IsUsed(go, "--popen")      && fprintf(cfg->ofp, "# gap open probability:             %f\n",         esl_opt_GetReal   (go, "--popen"))   < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  //if (esl_opt_IsUsed(go, "--pextend")    && fprintf(cfg->ofp, "# gap extend probability:           %f\n",         esl_opt_GetReal   (go, "--pextend")) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  //if (esl_opt_IsUsed(go, "--mx")         && fprintf(cfg->ofp, "# subst score matrix (built-in):    %s\n",         esl_opt_GetString (go, "--mx"))      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  //if (esl_opt_IsUsed(go, "--mxfile")     && fprintf(cfg->ofp, "# subst score matrix (file):        %s\n",         esl_opt_GetString (go, "--mxfile"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--maxinsertlen")  && fprintf(cfg->ofp, "# max insert length:                %d\n",         esl_opt_GetInteger (go, "--maxinsertlen"))  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");

  if (esl_opt_IsUsed(go, "--seed"))  {
    if (esl_opt_GetInteger(go, "--seed") == 0  && fprintf(cfg->ofp,"# random number seed:               one-time arbitrary\n")                        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
    else if                              (  fprintf(cfg->ofp,"# random number seed set to:        %d\n",         esl_opt_GetInteger(go, "--seed"))    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  }
  if (esl_opt_IsUsed(go, "--w_beta")     && fprintf(cfg->ofp, "# window length tail mass:          %g bits\n",   esl_opt_GetReal(go, "--w_beta"))     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  if (esl_opt_IsUsed(go, "--w_length")   && fprintf(cfg->ofp, "# window length :                   %d\n",        esl_opt_GetInteger(go, "--w_length"))< 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");


  if (fprintf(cfg->ofp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n\n") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "write failed");
  return eslOK;
}

static int
output_result(const struct cfg_s *cfg, char *errbuf, int msaidx, ESL_MSA *msa, F4_HMM *hmm, double entropy)
{
  int status;

  /* Special case: output the tabular results header. 
   * Arranged this way to keep the two fprintf()'s close together in the code,
   * so we can keep the data and labels properly sync'ed.
   */
  if (msa == NULL)
  {
    if (cfg->abc->type == eslAMINO) {
      if (fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s %8s %6s %s\n", " idx", "name",                 "nseq",  "alen",  "mlen",  "eff_nseq",  "re/pos",  "description")     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "output_result: write failed");
      if (fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s %8s %6s %s\n", "----", "--------------------", "-----", "-----", "-----", "--------",  "------",  "-----------") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "output_result: write failed");
    } else {
      if (fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s %5s %8s %6s %s\n", " idx", "name",                 "nseq",  "alen",  "mlen",  "W", "eff_nseq",  "re/pos",  "description")     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "output_result: write failed");
      if (fprintf(cfg->ofp, "#%4s %-20s %5s %5s %5s %5s %8s %6s %s\n", "----", "--------------------", "-----", "-----", "-----", "-----", "--------",  "------",  "-----------") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "output_result: write failed");
    }
    return eslOK;
  }
  if ((status = f4_hmmfile_WriteASCII(cfg->hmmfp, -1, hmm)) != eslOK) ESL_FAIL(status, errbuf, "HMM save failed");

	             /* #   name nseq alen M max_length eff_nseq re/pos description */
  if (cfg->abc->type == eslAMINO) {
    if (fprintf(cfg->ofp, "%-5d %-20s %5d %5" PRId64 " %5d %8.2f %6.3f %s\n",
          msaidx,
          (msa->name != NULL) ? msa->name : "",
          msa->nseq,
          msa->alen,
          hmm->M,
          hmm->eff_nseq,
          entropy,
          (msa->desc != NULL) ? msa->desc : "") < 0)
      ESL_EXCEPTION_SYS(eslEWRITE, "output_result: write failed");
  } else {
    if (fprintf(cfg->ofp, "%-5d %-20s %5d %5" PRId64 " %5d %5d %8.2f %6.3f %s\n",
          msaidx,
          (msa->name != NULL) ? msa->name : "",
          msa->nseq,
          msa->alen,
          hmm->M,
          hmm->max_length,
          hmm->eff_nseq,
          entropy,
          (msa->desc != NULL) ? msa->desc : "") < 0)
      ESL_EXCEPTION_SYS(eslEWRITE, "output_result: write failed");
  }

  return eslOK;
}


/* Function:  f4_MeanMatchRelativeEntropy()
 *
 * Purpose:   Calculate the mean relative entropy per match state emission
 *            distribution, in bits:
 *            
 *            \[
 *              \frac{1}{M} \sum_{k=1}^{M} \sum_x p_k(x) \log_2 \frac{p_k(x)}{f(x)}
 *            \]
 *       
 *            where $p_k(x)$ is emission probability for symbol $x$
 *            from match state $k$, and $f(x)$ is the null model's 
 *            background emission probability for $x$. 
 */
double
f4_MeanMatchRelativeEntropy(const F4_HMM *hmm, const F4_BG *bg)
{
  int    k;
  double KL = 0.;

  for (k = 1; k <= hmm->M; k++)
    KL += esl_vec_FRelEntropy(hmm->mat[k], bg->f, hmm->abc->K);
  KL /= (double) hmm->M;
  return KL;
}

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

  if (cfg->do_mpi == FALSE && cfg->nali == 1) /* first (only?) HMM in file: */
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

/*****************************************************************
 * 4. Serial loop, build HMMs for each MSA.
 *****************************************************************/

static void
serial_loop(WORKER_INFO *info, struct cfg_s *cfg, const ESL_GETOPTS *go)
{
  F4_BUILDER *bld         = NULL;
  ESL_MSA    *msa         = NULL;
  F4_HMM     *hmm         = NULL;
  char        errmsg[eslERRBUFSIZE];
  int         status;

  double      entropy;

  cfg->nali = 0;
  while ((status = esl_msafile_Read(cfg->afp, &msa)) != eslEOF)
    {
      if (status != eslOK) esl_msafile_ReadFailure(cfg->afp, status);
      cfg->nali++;  

      if ((status = set_msa_name(cfg, errmsg, msa)) != eslOK) f4_Fail("%s\n", errmsg); /* cfg->nnamed gets incremented in this call */

      /*         bg   new-HMM trarr gm   om  */
      if ((status = f4_Builder(info->bld, msa, info->bg, &hmm, NULL)) != eslOK) f4_Fail("build failed: %s", bld->errbuf);

      entropy = f4_MeanMatchRelativeEntropy(hmm, info->bg);
      if ((status = output_result(cfg, errmsg, cfg->nali, msa, hmm, entropy))         != eslOK) f4_Fail(errmsg);

      f4_hmm_Destroy(hmm);
      esl_msa_Destroy(msa);
    }
}

/*****************************************************************
 * 5. Master that builds HMMs for each MSA.
 *****************************************************************/

/* usual_master()
 * The usual version of dmmbuild (serial)
 * For each MSA, build an HMM and save it.
 * 
 * A master can only return if it's successful. 
 * All errors are handled immediately and fatally with f4_Fail() or equiv.
 */
static int
usual_master(const ESL_GETOPTS *go, struct cfg_s *cfg)
{
  int              ncpus    = 0;
  int              infocnt  = 0;
  WORKER_INFO     *info     = NULL;
  int              i;
  int              status;

  /* Open files, set alphabet.
   *   cfg->afp       - open alignment file for input
   *   cfg->abc       - alphabet expected or guessed in ali file
   *   cfg->hmmfp     - open HMM file for output
   *   cfp->ofp       - optional open output file, or stdout
   */
  if      (esl_opt_GetBoolean(go, "--amino"))   cfg->abc = esl_alphabet_Create(eslAMINO);
  else if (esl_opt_GetBoolean(go, "--dna"))     cfg->abc = esl_alphabet_Create(eslDNA);
  else if (esl_opt_GetBoolean(go, "--rna"))     cfg->abc = esl_alphabet_Create(eslRNA);
  else                                          cfg->abc = NULL;
  
  status = esl_msafile_Open(&(cfg->abc), cfg->alifile, NULL, cfg->fmt, NULL, &(cfg->afp));
  if (status != eslOK) esl_msafile_OpenFailure(cfg->afp, status);

  cfg->hmmfp = fopen(cfg->hmmfile, "w");
  if (cfg->hmmfp == NULL) f4_Fail("Failed to open HMM file %s for writing", cfg->hmmfile);

  if (esl_opt_IsUsed(go, "-o")) 
    {
      cfg->ofp = fopen(esl_opt_GetString(go, "-o"), "w");
      if (cfg->ofp == NULL) f4_Fail("Failed to open -o output file %s\n", esl_opt_GetString(go, "-o"));
    } 
  else cfg->ofp = stdout;


  /* Looks like the i/o is set up successfully...
   * Initial output to the user
   */
  output_header(go, cfg);                                  /* cheery output header                                */
  output_result(cfg, NULL, 0, NULL, NULL, 0.0);	   /* tabular results header (with no args, special-case) */

  infocnt = (ncpus == 0) ? 1 : ncpus;
  ESL_ALLOC(info, (ptrdiff_t) sizeof(*info) * infocnt);

  for (i = 0; i < infocnt; ++i)
  {
      info[i].bg = f4_bg_Create(cfg->abc);
      info[i].bld = f4_builder_Create(go, cfg->abc);

      if (info[i].bld == NULL)  f4_Fail("f4_builder_Create failed");

      //do this here instead of in f4_builder_Create(), because it's an hmmbuild-specific option

      if ( esl_opt_IsOn(go, "--maxinsertlen") )
        info[i].bld->max_insert_len    = esl_opt_GetInteger(go, "--maxinsertlen");

      /* special arguments for hmmbuild */
      info[i].bld->w_len      = (go != NULL && esl_opt_IsOn (go, "--w_length")) ?  esl_opt_GetInteger(go, "--w_length"): -1;
      info[i].bld->w_beta     = (go != NULL && esl_opt_IsOn (go, "--w_beta"))   ?  esl_opt_GetReal   (go, "--w_beta")    : f4_DEFAULT_WINDOW_BETA;
      if ( info[i].bld->w_beta < 0 || info[i].bld->w_beta > 1  ) esl_fatal("Invalid window-length beta value\n");
  }

  serial_loop(info, cfg, go);

  for (i = 0; i < infocnt; ++i)
  {
      f4_bg_Destroy(info[i].bg);
      f4_builder_Destroy(info[i].bld);
  }

  free(info);
  return eslOK;

 ERROR:
  return eslFAIL;
}

/*****************************************************************
 * 6. Main program.
 *****************************************************************/

int
main(int argc, char **argv)
{
  ESL_GETOPTS     *go = NULL;	/* command line processing                 */
  ESL_STOPWATCH   *w  = esl_stopwatch_Create();
  struct cfg_s     cfg;

  /* Set processor specific flags */
  //impl_Init();

  cfg.alifile     = NULL;
  cfg.hmmfile     = NULL;

  /* Parse the command line
   */
  process_commandline(argc, argv, &go, &cfg.hmmfile, &cfg.alifile);    

  /* Initialize what we can in the config structure (without knowing the alphabet yet).
   * Fields controlled by masters are set up in usual_master()
   */
  cfg.ofp         = NULL;	           
  cfg.fmt         = eslMSAFILE_UNKNOWN;    /* autodetect alignment format by default. */ 
  cfg.afp         = NULL;	           
  cfg.abc         = NULL;	           
  cfg.hmmfp       = NULL;	           

  cfg.nali       = 0;		           /* this counter is incremented in masters */
  cfg.nnamed     = 0;		           /* 0 or 1 if a single MSA; == nali if multiple MSAs */
  cfg.do_mpi     = FALSE;	           /* this gets reset below, if we init MPI */
  cfg.nproc      = 0;		           /* this gets reset below, if we init MPI */
  cfg.my_rank    = 0;		           /* this gets reset below, if we init MPI */
  cfg.do_stall   = esl_opt_GetBoolean(go, "--stall");
  cfg.hmmName    = esl_opt_GetString(go, "-n"); /* NULL by default */

  if (esl_opt_IsOn(go, "--informat")) {
    cfg.fmt = esl_msafile_EncodeFormat(esl_opt_GetString(go, "--informat"));
    if (cfg.fmt == eslMSAFILE_UNKNOWN) f4_Fail("%s is not a recognized input sequence file format\n", esl_opt_GetString(go, "--informat"));
  }

  /* Start timing. */
  esl_stopwatch_Start(w);

  usual_master(go, &cfg);
  esl_stopwatch_Stop(w);

  if (cfg.my_rank == 0) {
    fputc('\n', cfg.ofp);
    esl_stopwatch_Display(cfg.ofp, w, "# CPU time: ");
  }

  /* Clean up the shared cfg. 
   */
  if (cfg.my_rank == 0) {
    if (esl_opt_IsOn(go, "-o")) { fclose(cfg.ofp); }
    if (cfg.afp)   esl_msafile_Close(cfg.afp);
    if (cfg.abc)   esl_alphabet_Destroy(cfg.abc);
    if (cfg.hmmfp) fclose(cfg.hmmfp);
  }
  esl_getopts_Destroy(go);
  esl_stopwatch_Destroy(w);
  return 0;
}