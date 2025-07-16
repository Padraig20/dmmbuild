/* F4_BG: the null (background) model */

#include <string.h>
#include "dummer.h"

/* Function:  f4_bg_Create()
 * Synopsis:  Create a <F4_BG> null model object.
 *
 * Purpose:   Allocate a <F4_BG> object for digital alphabet <abc>,
 *            initializes it to appropriate default values, and
 *            returns a pointer to it.
 *            
 *            For protein models, default iid background frequencies
 *            are set (by <f4_AminoFrequencies()>) to average
 *            Swiss-Prot residue composition. For DNA, RNA and other
 *            alphabets, default frequencies are set to a uniform
 *            distribution.
 *            
 *            The model composition <bg->mcomp[]> is not initialized
 *            here; neither is the filter null model <bg->fhmm>.  To
 *            use the filter null model, caller will want to
 *            initialize these fields by calling
 *            <f4_bg_SetFilter()>.
 *
 * Throws:    <NULL> on allocation failure.
 */
F4_BG *
f4_bg_Create(const ESL_ALPHABET *abc)
{
  F4_BG *bg = NULL;
  int    status;

  ESL_ALLOC(bg, sizeof(F4_BG));
  bg->f     = NULL;
  bg->fhmm  = NULL;

  ESL_ALLOC(bg->f,     sizeof(float) * abc->K);
  if ((bg->fhmm = esl_hmm_Create(abc, 2)) == NULL) goto ERROR;

  if       (abc->type == eslAMINO)
    {
      if (f4_AminoFrequencies(bg->f) != eslOK) goto ERROR;
    }
  else
    esl_vec_FSet(bg->f, abc->K, 1. / (float) abc->K);

  bg->p1    = 350./351.;
  bg->omega = 1./256.;
  bg->abc   = abc;
  return bg;

 ERROR:
  f4_bg_Destroy(bg);
  return NULL;
}

/* Function:  f4_bg_Destroy()
 *
 * Purpose:   Frees a <F4_BG> object.
 *
 * Returns:   (void)
 */
void
f4_bg_Destroy(F4_BG *bg)
{
  if (bg != NULL) {
    if (bg->f     != NULL) free(bg->f);
    if (bg->fhmm  != NULL) esl_hmm_Destroy(bg->fhmm);
    free(bg);
  }
  return;
}