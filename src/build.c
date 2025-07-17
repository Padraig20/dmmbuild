/* Construction of new HMMs from multiple alignments.
 * 
 * f4_Fastmodelmaker() -- Krogh/Haussler heuristic
 * 
 * The maximum likelihood model construction algorithm that was in previous
 * HMMER versions has been deprecated, at least for the moment.
 * 
 * The meat of the model construction code is in matassign2hmm().
 * The two model construction strategies simply label which columns
 * are supposed to be match states, and then hand this info to
 * matassign2hmm().
 * 
 * Contents:
 *    1. Exported API: model construction routines.
 *    2. Private functions used in constructing models.
 */

#include "dummer.h"

static int do_modelmask( ESL_MSA *msa);
static int matassign2hmm(ESL_MSA *msa, int *matassign, F4_HMM **ret_hmm, F4_TRACE ***opt_tr);
static int annotate_model(F4_HMM *hmm, int *matassign, ESL_MSA *msa);

/*****************************************************************
 * 1. Exported API: model construction routines.
 *****************************************************************/

/* Function: f4_Fastmodelmaker()
 * 
 * Purpose:  Heuristic model construction.
 *           Construct an HMM from an alignment by a simple rule,
 *           based on the fractional occupancy of each columns w/
 *           residues vs gaps. Any column w/ a fractional
 *           occupancy of $\geq$ <symfrac> is assigned as a MATCH column;
 *           for instance, if thresh = 0.5, columns w/ $\geq$ 50\% 
 *           residues are assigned to match... roughly speaking.
 *           
 *           "Roughly speaking" because sequences may be weighted
 *           in the input <msa>, and because missing data symbols are
 *           ignored, in order to deal with sequence fragments.
 *
 *           The <msa> must be in digital mode. 
 *
 *           If the caller wants to designate any sequences as
 *           fragments, it does so by converting all N-terminal and
 *           C-terminal flanking gap symbols to missing data symbols.
 *
 *           Returns the HMM in counts form (ready for applying Dirichlet
 *           priors as the next step). Also returns fake traceback
 *           for each training sequence.
 *           
 *           Models must have at least one node, so if the <msa> defined 
 *           no consensus columns, a <eslENORESULT> error is returned.
 *           
 * Args:     msa       - multiple sequence alignment
 *           symfrac   - threshold for residue occupancy; >= assigns MATCH
 *           bld       - holds information on regions requiring masking, optionally NULL -> no masking
 *           ret_hmm   - RETURN: counts-form HMM
 *           opt_tr    - optRETURN: array of tracebacks for aseq's
 *           
 * Return:   <eslOK> on success. ret_hmm and opt_tr allocated here,
 *           and must be free'd by the caller (FreeTrace(tr[i]), free(tr),
 *           FreeHMM(hmm)).       
 *
 *           Returns <eslENORESULT> if no consensus columns were annotated;
 *           in this case, <ret_hmm> and <opt_tr> are returned NULL.
 *           
 * Throws:   <eslEMEM> on allocation failure; <eslEINVAL> if the 
 *           <msa> isn't in digital mode.
 */
int
f4_Fastmodelmaker(ESL_MSA *msa, float symfrac, F4_BUILDER *bld, F4_HMM **ret_hmm, F4_TRACE ***opt_tr)
{
  int      status;	     /* return status flag                  */
  int     *matassign = NULL; /* MAT state assignments if 1; 1..alen */
  int      idx;              /* counter over sequences              */
  int      apos;             /* counter for aligned columns         */
  float    r;		         /* weighted residue count              */
  float    totwgt;	     /* weighted residue+gap count          */

  if (! (msa->flags & eslMSA_DIGITAL)) ESL_XEXCEPTION(eslEINVAL, "need digital MSA");

  /* Allocations: matassign is 1..alen array of bit flags. */
  matassign = calloc(msa->alen + 1, sizeof(int));
  if (matassign == NULL) ESL_XEXCEPTION(eslEMEM, "allocation failed for matassign");

  /* Determine weighted sym freq in each column, set matassign[] accordingly. */
  for (apos = 1; apos <= msa->alen; apos++) 
    {  
      r = totwgt = 0.;
      for (idx = 0; idx < msa->nseq; idx++) 
      {
        if       (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) { r += msa->wgt[idx]; totwgt += msa->wgt[idx]; }
        else if  (esl_abc_XIsGap(msa->abc,     msa->ax[idx][apos])) {                     totwgt += msa->wgt[idx]; }
        else if  (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos])) continue;
      }
      if (r > 0. && r / totwgt >= symfrac) matassign[apos] = TRUE;
      else                                 matassign[apos] = FALSE;
    }


  /* Once we have matassign calculated, modelmakers behave
   * the same; matassign2hmm() does this stuff (traceback construction,
   * trace counting) and sets up ret_hmm and opt_tr.
   */
  if ((status = matassign2hmm(msa, matassign, ret_hmm, opt_tr)) != eslOK) {
    fprintf (stderr, "hmm construction error during trace counting\n");
    goto ERROR;
  }

  free(matassign);
  return eslOK;

 ERROR:
  if (matassign != NULL) free(matassign);
  return status;
}

/*-------------------- end, exported API -------------------------*/

/*****************************************************************
 * 2. Private functions used in constructing models.
 *****************************************************************/ 

/* Function: do_modelmask()
 *
 * Purpose:  If the given <msa> has a MM CS line, mask (turn to
 *           degenerate) residues in the msa positions associated
 *           with the marked position in the MM (marked with 'm')
 *
 * Return:   <eslOK> on success.
 *           <eslENORESULT> if error.
 */
static int
do_modelmask( ESL_MSA *msa)
{
  int i,j;

  if (msa->mm == NULL)  return eslOK;  //nothing to do

  for (i = 1; i <= msa->alen; i++) {
    for (j = 0; j < msa->nseq; j++) {
      if (msa->mm[i-1] == 'm') {
        if (msa->ax[j][i] != msa->abc->K && msa->ax[j][i] != msa->abc->Kp-1) // if not gap
          msa->ax[j][i] = msa->abc->Kp-3; //that's the degenerate "any character" (N for DNA, X for protein)
      }
    }
  }
  return eslOK;
}

/* Function: matassign2hmm()
 * 
 * Purpose:  Given an assignment of alignment columns to match vs.
 *           insert, finish the final part of the model construction 
 *           calculation that is constant between model construction
 *           algorithms.
 *           
 * Args:     msa       - multiple sequence alignment
 *           matassign - 1..alen bit flags for column assignments
 *           ret_hmm   - RETURN: counts-form HMM
 *           opt_tr    - optRETURN: array of tracebacks for aseq's
 *                         
 * Return:   <eslOK> on success.
 *           <eslENORESULT> if no consensus columns are identified.
 *
 *           ret_hmm and opt_tr alloc'ed here.
 */
static int
matassign2hmm(ESL_MSA *msa, int *matassign, F4_HMM **ret_hmm, F4_TRACE ***opt_tr)
{
  int        status;		/* return status                       */
  F4_HMM    *hmm = NULL;        /* RETURN: new hmm                     */
  F4_TRACE **tr  = NULL;        /* RETURN: 0..nseq-1 fake traces       */
  int      M;                   /* length of new model in match states */
  int      idx;                 /* counter over sequences              */
  int      apos;                /* counter for aligned columns         */
  char errbuf[eslERRBUFSIZE];

  /* apply the model mask in the 'GC MM' row */
  do_modelmask(msa);

  /* How many match states in the HMM? */
  for (M = 0, apos = 1; apos <= msa->alen; apos++) 
    if (matassign[apos]) M++;
  if (M == 0) { status = eslENORESULT; goto ERROR; }

  /* Make fake tracebacks for each seq */
  ESL_ALLOC(tr, sizeof(F4_TRACE *) * msa->nseq);
  if ((status = f4_trace_FauxFromMSA(msa, matassign, f4_MSA_COORDS, tr))        != eslOK) goto ERROR;

  for (idx = 0; idx < msa->nseq; idx++)
    {
      if ((status = f4_trace_Validate(tr[idx], msa->abc, msa->ax[idx], errbuf)) != eslOK) 
	ESL_XEXCEPTION(eslFAIL, "validation failed: %s", errbuf);
    }
  
  /* Build count model from tracebacks */
  if ((hmm    = f4_hmm_Create(M, msa->abc)) == NULL)  { status = eslEMEM; goto ERROR; }
  if ((status = f4_hmm_Zero(hmm))           != eslOK) goto ERROR;
  for (idx = 0; idx < msa->nseq; idx++) {
    if (tr[idx] == NULL) continue; /* skip rare examples of empty sequences */
    if ((status = f4_trace_Count(hmm, msa->ax[idx], msa->wgt[idx], tr[idx])) != eslOK) goto ERROR;
  }

  // now hmm has the initial profile (aka counts), which will now be refined via Baum-Welch

  double **letter_probs;    // [0..msa->alen-1][1..hmm->abc->K-1] letter probabilities
  double *background_probs; // [1..hmm->abc->K-1]                 background probabilities

  if ((status = letter_probs_build(hmm->M, msa->abc->K, &letter_probs, &background_probs)) != eslOK) goto ERROR;

  for (idx = 0; idx < msa->nseq; idx++) {
    if ((status = letter_probs_count(hmm->M, msa->alen, msa->ax[idx], msa->wgt[idx], msa->abc, letter_probs, background_probs, matassign)) != eslOK) goto ERROR;
  }

  if ((status = letter_probs_normalize(hmm->M, msa->abc->K, letter_probs, background_probs))             != eslOK) goto ERROR;

  if ((status = f4_trace_Estimate(hmm, msa, tr, letter_probs, background_probs)) != eslOK) goto ERROR;

  letter_probs_destroy(hmm->M, letter_probs, background_probs);

  hmm->nseq     = msa->nseq;
  hmm->eff_nseq = msa->nseq;

  /* Transfer annotation from the MSA to the new model */
  if ((status = annotate_model(hmm, matassign, msa)) != eslOK) goto ERROR;

  /* Reset #=RF line of alignment to reflect our assignment
   * of match, delete. matassign is valid from 1..alen and is off
   * by one from msa->rf.
   */
  if (msa->rf == NULL)  ESL_ALLOC(msa->rf, sizeof(char) * (msa->alen + 1));
  for (apos = 1; apos <= msa->alen; apos++)
    msa->rf[apos-1] = matassign[apos] ? 'x' : '.';
  msa->rf[msa->alen] = '\0';

  if (opt_tr  != NULL) *opt_tr  = tr; 
  else                  f4_trace_DestroyArray(tr, msa->nseq);
  *ret_hmm = hmm;
  return eslOK;

 ERROR:
  if (tr     != NULL) f4_trace_DestroyArray(tr, msa->nseq);
  if (hmm    != NULL) f4_hmm_Destroy(hmm);
  if (opt_tr != NULL) *opt_tr = NULL;
  if (letter_probs != NULL || background_probs != NULL) letter_probs_destroy(msa->alen, letter_probs, background_probs);
  *ret_hmm = NULL;
  return status;
}

/* Function: annotate_model()
 * 
 * Purpose:  Transfer rf, cs, and other optional annotation from the alignment
 *           to the new model.
 * 
 * Args:     hmm       - [M] new model to annotate
 *           matassign - which alignment columns are MAT; [1..alen]
 *           msa       - alignment, including annotation to transfer
 *           
 * Return:   <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error.
 */
static int
annotate_model(F4_HMM *hmm, int *matassign, ESL_MSA *msa)
{                      
  int   apos;			/* position in matassign, 1.alen  */
  int   k;			/* position in model, 1.M         */
  int   status;

  /* Reference coord annotation  */
  if (msa->rf != NULL) {
    ESL_ALLOC(hmm->rf, sizeof(char) * (hmm->M+2));
    hmm->rf[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++) 
      if (matassign[apos]) hmm->rf[k++] = msa->rf[apos-1]; /* watch off-by-one in msa's rf */
    hmm->rf[k] = '\0';
    hmm->flags |= f4H_RF;
  }

  /* Model mask annotation  */
  if (msa->mm != NULL) {
    ESL_ALLOC(hmm->mm, sizeof(char) * (hmm->M+2));
    hmm->mm[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->mm[k++] = ( msa->mm[apos-1] == '.' ? '-' : msa->mm[apos-1]) ;
    hmm->mm[k] = '\0';
    hmm->flags |= f4H_MMASK;
  }

  /* Consensus structure annotation */
  if (msa->ss_cons != NULL) {
    ESL_ALLOC(hmm->cs, sizeof(char) * (hmm->M+2));
    hmm->cs[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->cs[k++] = msa->ss_cons[apos-1];
    hmm->cs[k] = '\0';
    hmm->flags |= f4H_CS;
  }

  /* Surface accessibility annotation */
  if (msa->sa_cons != NULL) {
    ESL_ALLOC(hmm->ca, sizeof(char) * (hmm->M+2));
    hmm->ca[0] = ' ';
    for (apos = k = 1; apos <= msa->alen; apos++)
      if (matassign[apos]) hmm->ca[k++] = msa->sa_cons[apos-1];
    hmm->ca[k] = '\0';
    hmm->flags |= f4H_CA;
  }

  /* The alignment map (1..M in model, 1..alen in alignment) */
  ESL_ALLOC(hmm->map, sizeof(int) * (hmm->M+1));
  hmm->map[0] = 0;
  for (apos = k = 1; apos <= msa->alen; apos++)
    if (matassign[apos]) hmm->map[k++] = apos;
  hmm->flags |= f4H_MAP;

  return eslOK;

 ERROR:
  return status;
}