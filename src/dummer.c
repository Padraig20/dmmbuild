#include <math.h>
#include <float.h>

#include "easel.h"
#include "esl_getopts.h"
#include "dummer.h"

/*****************************************************************
 * 1. Miscellaneous functions
 *****************************************************************/

/* Function:  f4_banner()
 * Synopsis:  print standard DUMMER application output header
 *
 * Returns:   (void)
 */
void
f4_banner(FILE *fp, const char *progname, char *banner)
{
  char *appname = NULL;

  if (esl_FileTail(progname, FALSE, &appname) != eslOK) esl_strdup(progname, -1, &appname);

  fprintf(fp, "# %s :: %s\n", appname, banner);
  fprintf(fp, "# DUMMER %s (%s); %s\n", DUMMER_VERSION, DUMMER_DATE, DUMMER_URL);
  fprintf(fp, "# %s\n", DUMMER_COPYRIGHT);
  fprintf(fp, "# %s\n", DUMMER_LICENSE);
  fprintf(fp, "# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");

  if (appname) free(appname);
  return;
}

/*****************************************************************
 * 2. Error functions
 *****************************************************************/

/* Function:  f4_Fail()
 * Synopsis:  Handle a user error (something that's the user's fault).
 */
void
f4_Fail(char *format, ...)
{
  va_list  argp;
                                /* format the error mesg */
  fprintf(stderr, "\nError: ");
  va_start(argp, format);
  vfprintf(stderr, format, argp);
  va_end(argp);
  fprintf(stderr, "\n");
  fflush(stderr);
  exit(1);
}

/* Function:  f4_AminoFrequencies()
 * Incept:    SRE, Fri Jan 12 13:46:41 2007 [Janelia]
 *
 * Purpose:   Fills a vector <f> with amino acid background frequencies,
 *            in [A..Y] alphabetic order, same order that Easel digital
 *            alphabet uses. Caller must provide <f> allocated for at
 *            least 20 floats.
 *            
 *            These were updated 4 Sept 2007, from Swiss-Prot 50.8,
 *            (Oct 2006), counting over 85956127 (86.0M) residues.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_AminoFrequencies(float *f)
{
  f[0] = 0.0787945;		/* A */
  f[1] = 0.0151600;		/* C */
  f[2] = 0.0535222;		/* D */
  f[3] = 0.0668298;		/* E */
  f[4] = 0.0397062;		/* F */
  f[5] = 0.0695071;		/* G */
  f[6] = 0.0229198;		/* H */
  f[7] = 0.0590092;		/* I */
  f[8] = 0.0594422;		/* K */
  f[9] = 0.0963728;		/* L */
  f[10]= 0.0237718;		/* M */
  f[11]= 0.0414386;		/* N */
  f[12]= 0.0482904;		/* P */
  f[13]= 0.0395639;		/* Q */
  f[14]= 0.0540978;		/* R */
  f[15]= 0.0683364;		/* S */
  f[16]= 0.0540687;		/* T */
  f[17]= 0.0673417;		/* V */
  f[18]= 0.0114135;		/* W */
  f[19]= 0.0304133;		/* Y */
  return eslOK;
}