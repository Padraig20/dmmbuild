/* Printing the result to an HMMER3 ASCII save file.
 * 
 * Contents:
 *    1. Helper functions for i/o stuff.
 *    2. Writing HMMER3 HMM files.
 */

#include "dummer.h"

/*****************************************************************
 * 1. Helper functions for i/o stuff
 *****************************************************************/

 /* These tags need to be in temporal order, so we can do tests
 * like "if (format >= p7_HMMFILE_3b) ..."
 */
enum p7_hmmfile_formats_e {
  f4_HMMFILE_20 = 0,
  f4_HMMFILE_3a = 1,
  f4_HMMFILE_3b = 2,
  f4_HMMFILE_3c = 3,
  f4_HMMFILE_3d = 4,
  f4_HMMFILE_3e = 5,
  f4_HMMFILE_3f = 6,
};

/* multiline()
 * 
 * Used to print the command log to ASCII save files.
 *
 * Given a record (like the comlog) that contains 
 * multiple lines, print it as multiple lines with
 * a given prefix. e.g.:
 *           
 * given:   "COM   ", "foo\nbar\nbaz"
 * print:   COM   1 foo
 *          COM   2 bar
 *          COM   3 baz
 *
 * If <s> is NULL, no-op. Otherwise <s> must be a <NUL>-terminated
 * string.  It does not matter if it ends in <\n> or not. <pfx>
 * must be a valid <NUL>-terminated string; it may be empty.
 *           
 * Args:     fp:   FILE to print to
 *           pfx:  prefix for each line
 *           s:    line to break up and print; tolerates a NULL
 *
 * Returns: <eslOK> on success.
 *
 * Throws:  <eslEWRITE> on write error.
 */
static int
multiline(FILE *fp, const char *pfx, char *s)
{
  char *sptr  = s;
  char *end   = NULL;
  size_t   n     = 0;
  int   nline = 1;

  do {
    end = strchr(sptr, '\n');

    if (end != NULL)                  /* if there's no \n left, end == NULL */
      {
  n = end - sptr;                       /* n chars exclusive of \n */
  if (fprintf(fp, "%s [%d] ", pfx, nline++) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fwrite(sptr, sizeof(char), n, fp)    != n) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");  /* using fwrite lets us write fixed # of chars   */
  if (fprintf(fp, "\n")                     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");  /* while writing \n w/ printf allows newline conversion */
  sptr += n + 1;                       /* +1 to get past \n */
      } 
    else 
      {
  if (fprintf(fp, "%s [%d] %s\n", pfx, nline++, sptr) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); /* last line */
      }
  } while (end != NULL  && *sptr != '\0');   /* *sptr == 0 if <s> terminates with a \n */
  return eslOK;
}

static int
printprob(FILE *fp, int fieldwidth, float p)
{
  if      (p == 0.0) { if (fprintf(fp, " %*s",   fieldwidth, "*")      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else if (p == 1.0) { if (fprintf(fp, " %*.5f", fieldwidth, 0.0)      < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else               { if (fprintf(fp, " %*.5f", fieldwidth, -logf(p)) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  return eslOK;
}

/*****************************************************************
 * 2. Writing HMMER3 HMM files.
 *****************************************************************/

/* Function:  f4_hmmfile_WriteASCII()
 * Synopsis:  Write a HMMER3 ASCII save file.
 *
 * Purpose:   Write a profile HMM <hmm> in an ASCII save file format to
 *            an open stream <fp>.
 *
 *            Legacy file formats in the 3.x release series are
 *            supported by specifying the <format> code. Pass <-1> to
 *            use the default current standard format; pass a valid
 *            code such as <f4_HMMFILE_3a> to select a specific
 *            format.
 *
 * Args:      fp     - open stream for writing
 *            format - -1 for default format, or a 3.x format code like <f4_HMMFILE_3a>
 *            hmm    - HMM to save
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <format> isn't a valid 3.0 format code.
 *            <eslEWRITE> on write error.
 */
int
f4_hmmfile_WriteASCII(FILE *fp, int format, F4_HMM *hmm)
{
  int k, x;
  int status;
  

  if (format == -1) format = f4_HMMFILE_3f;

  if      (format == f4_HMMFILE_3f)  { if (fprintf(fp, "HMMER3/f [%s | %s]\n",                             DUMMER_VERSION, DUMMER_DATE) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");}
  else if (format == f4_HMMFILE_3e)  { if (fprintf(fp, "HMMER3/e [%s | %s; reverse compatibility mode]\n", DUMMER_VERSION, DUMMER_DATE) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else if (format == f4_HMMFILE_3d)  { if (fprintf(fp, "HMMER3/d [%s | %s; reverse compatibility mode]\n", DUMMER_VERSION, DUMMER_DATE) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else if (format == f4_HMMFILE_3c)  { if (fprintf(fp, "HMMER3/c [%s | %s; reverse compatibility mode]\n", DUMMER_VERSION, DUMMER_DATE) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else if (format == f4_HMMFILE_3b)  { if (fprintf(fp, "HMMER3/b [%s | %s; reverse compatibility mode]\n", DUMMER_VERSION, DUMMER_DATE) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else if (format == f4_HMMFILE_3a)  { if (fprintf(fp, "HMMER3/a [%s | %s; reverse compatibility mode]\n", DUMMER_VERSION, DUMMER_DATE) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  else ESL_EXCEPTION(eslEINVAL, "invalid HMM file format code");
  
  if (fprintf(fp, "NAME  %s\n", hmm->name)                                                          < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (hmm->acc  && fprintf(fp, "ACC   %s\n", hmm->acc)                                              < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (hmm->desc && fprintf(fp, "DESC  %s\n", hmm->desc)                                             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fprintf(fp, "LENG  %d\n", hmm->M)                                                             < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (format >= f4_HMMFILE_3c && hmm->max_length > 0 && fprintf(fp, "MAXL  %d\n", hmm->max_length)  < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fprintf(fp, "ALPH  %s\n", esl_abc_DecodeType(hmm->abc->type))                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fprintf(fp, "RF    %s\n", (hmm->flags & f4H_RF)    ? "yes" : "no")                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (format >= f4_HMMFILE_3f && fprintf(fp, "MM    %s\n", (hmm->flags & f4H_MMASK) ? "yes" : "no") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (format >= f4_HMMFILE_3e && fprintf(fp, "CONS  %s\n", (hmm->flags & f4H_CONS)  ? "yes" : "no") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fprintf(fp, "CS    %s\n", (hmm->flags & f4H_CS)    ? "yes" : "no")                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fprintf(fp, "MAP   %s\n", (hmm->flags & f4H_MAP)   ? "yes" : "no")                            < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (hmm->ctime    != NULL)   { if (           fprintf  (fp, "DATE  %s\n", hmm->ctime)        < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  if (hmm->comlog   != NULL)   { if ( (status = multiline(fp, "COM  ",      hmm->comlog)) != eslOK) return status; }
  if (hmm->nseq     >  0)      { if (           fprintf  (fp, "NSEQ  %d\n", hmm->nseq)         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  if (hmm->eff_nseq >= 0)      { if (           fprintf  (fp, "EFFN  %f\n", hmm->eff_nseq)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  if (hmm->flags & f4H_CHKSUM) { if (           fprintf  (fp, "CKSUM %u\n", hmm->checksum)     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); } /* unsigned 32-bit */

  if (hmm->abc->type == eslRNA || hmm->abc->type == eslDNA ) {
    if ((hmm->flags & f4H_GA)  && fprintf(fp, "GA    %.2f\n", hmm->cutoff[f4_GA1]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    if ((hmm->flags & f4H_TC)  && fprintf(fp, "TC    %.2f\n", hmm->cutoff[f4_TC1]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    if ((hmm->flags & f4H_NC)  && fprintf(fp, "NC    %.2f\n", hmm->cutoff[f4_NC1]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  } else {
    if ((hmm->flags & f4H_GA)  && fprintf(fp, "GA    %.2f %.2f\n", hmm->cutoff[f4_GA1], hmm->cutoff[f4_GA2]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    if ((hmm->flags & f4H_TC)  && fprintf(fp, "TC    %.2f %.2f\n", hmm->cutoff[f4_TC1], hmm->cutoff[f4_TC2]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    if ((hmm->flags & f4H_NC)  && fprintf(fp, "NC    %.2f %.2f\n", hmm->cutoff[f4_NC1], hmm->cutoff[f4_NC2]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  }
  if (hmm->flags & f4H_STATS) {
    if (format == f4_HMMFILE_3a)  {        /* reverse compatibility */
      if (fprintf(fp, "STATS LOCAL     VLAMBDA %f\n", hmm->evparam[f4_MLAMBDA]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
      if (fprintf(fp, "STATS LOCAL         VMU %f\n", hmm->evparam[f4_MMU])     < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
      if (fprintf(fp, "STATS LOCAL        FTAU %f\n", hmm->evparam[f4_FTAU])    < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    } else {        /* default stats lines */
      if (fprintf(fp, "STATS LOCAL MSV      %8.4f %8.5f\n", hmm->evparam[f4_MMU],  hmm->evparam[f4_MLAMBDA]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
      if (fprintf(fp, "STATS LOCAL VITERBI  %8.4f %8.5f\n", hmm->evparam[f4_VMU],  hmm->evparam[f4_VLAMBDA]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
      if (fprintf(fp, "STATS LOCAL FORWARD  %8.4f %8.5f\n", hmm->evparam[f4_FTAU], hmm->evparam[f4_FLAMBDA]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    }
  }

  if (fprintf(fp, "DMM     ")                                         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  for (x = 0; x < hmm->abc->K; x++) 
    { if (fprintf(fp, "     %c   ", hmm->abc->sym[x])                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
  if (fputc('\n', fp)                                                 < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (fprintf(fp, "        %8s %8s %8s %8s %8s %8s %8s %8s %8s\n",
        "m->m", "m->i", "m->d", "i->m", "i->i", "i->d", "d->m", "d->d", "d->i") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  if (hmm->flags & f4H_COMPO) {
    if (fprintf(fp, "  COMPO ") < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    for (x = 0; x < hmm->abc->K; x++) 
      { if ( (status = printprob(fp, 8, hmm->compo[x])) != eslOK) return status; }
    if (fputc('\n', fp)         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  }

  /* node 0 is special: insert emissions, and B-> transitions */
  if (fputs("        ", fp) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  for (x = 0; x < hmm->abc->K;      x++) 
    { if ( (status = printprob(fp, 8, hmm->ins[0][x])) != eslOK) return status; }  
  if (fputc('\n', fp)       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");

  if (fputs("        ", fp) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  for (x = 0; x < f4H_NTRANSITIONS; x++) 
    { if ( (status = printprob(fp, 8, hmm->t[0][x])) != eslOK) return status; }    
  if (fputc('\n', fp)       < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  for (k = 1; k <= hmm->M; k++) {
    /* Line 1: k; match emissions; optional map, RF, MM, CS */
    if (fprintf(fp, " %6d ",  k) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    for (x = 0; x < hmm->abc->K; x++)
    { if ( (status = printprob(fp, 8, hmm->mat[k][x])) != eslOK) return status; }

    if (hmm->flags & f4H_MAP) { if (fprintf(fp, " %6d", hmm->map[k]) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
    else                      { if (fprintf(fp, " %6s", "-")         < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }

    if (format >= f4_HMMFILE_3e) {
      if (format >= f4_HMMFILE_3f && (hmm->flags & f4H_MMASK) && hmm->mm[k] == 'm' )
        x = tolower(hmm->abc->sym[hmm->abc->Kp-3]);
      else if (hmm->flags & f4H_CONS)
        x = hmm->consensus[k];
      else
        x = '-';
      if (fprintf(fp, " %c", x) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    }


    if (hmm->rf && hmm->rf[k] == ' ') ESL_EXCEPTION_SYS(eslEWRITE, "input alignment contains an RF line with spaces");
    if (fprintf(fp, " %c",   (hmm->flags & f4H_RF)    ? hmm->rf[k]        : '-') < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    if (format >= f4_HMMFILE_3f) { if (fprintf(fp, " %c",   (hmm->flags & f4H_MMASK) ? hmm->mm[k]       : '-') < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed"); }
    if (fprintf(fp, " %c\n", (hmm->flags & f4H_CS)    ? hmm->cs[k]        : '-') < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");

    /* Line 2:   insert emissions */
    if (fputs("        ", fp) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    for (x = 0; x < hmm->abc->K; x++)
    { if ( (status = printprob(fp, 8, hmm->ins[k][x])) != eslOK) return status; }
    /* Line 3:   transitions */
    if (fputs("\n        ", fp) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
    for (x = 0; x < f4H_NTRANSITIONS; x++)
    { if ( (status = printprob(fp, 8, hmm->t[k][x])) != eslOK) return status; }
    if (fputc('\n', fp) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  }
  if (fputs("//\n", fp) < 0) ESL_EXCEPTION_SYS(eslEWRITE, "hmm write failed");
  return eslOK;
}