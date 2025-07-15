/* Model configuration: 
 * Converting a core model to a fully configured Fig4 search profile.
 * 
 * Contents:
 *     1. Routines in the exposed API.
 */

#include "dummer.h"

/*****************************************************************
 * 1. Routines in the exposed API.
 *****************************************************************/
 
/* Function:  f4_ProfileConfig()
 * Synopsis:  Configure a search profile.
 *
 * Purpose:   Given a model <hmm> with core probabilities, the null1
 *            model <bg>, a desired search <mode> (one of <f4_LOCAL>,
 *            <p4_GLOCAL>, <f4_UNILOCAL>, or <f4_UNIGLOCAL>), and an
 *            expected target sequence length <L>; configure the
 *            search model in <gm> with lod scores relative to the
 *            background frequencies in <bg>.
 *            
 * Returns:   <eslOK> on success; the profile <gm> now contains 
 *            scores and is ready for searching target sequences.
 *            
 * Throws:    <eslEMEM> on allocation error.
 */
int
f4_ProfileConfig(const F4_HMM *hmm, const F4_BG *bg, F4_PROFILE *gm, int L, int mode)
{
  int   k, x, z;	/* counters over states, residues, annotation */
  int   status;
  float *occ = NULL;
  float *tp, *rp;
  float  sc[f4_MAXCODE];
  float  Z;
 
  /* Contract checks */
  if (gm->abc->type != hmm->abc->type) ESL_XEXCEPTION(eslEINVAL, "HMM and profile alphabet don't match");
  if (hmm->M > gm->allocM)             ESL_XEXCEPTION(eslEINVAL, "profile too small to hold HMM");
  if (! (hmm->flags & f4H_CONS))       ESL_XEXCEPTION(eslEINVAL, "HMM must have a consensus to transfer to the profile");

  /* Copy some pointer references and other info across from HMM  */
  gm->M                = hmm->M;
  gm->max_length       = hmm->max_length;
  gm->mode             = mode;
  gm->roff             = -1;
  gm->eoff             = -1;
  gm->offs[f4_MOFFSET] = -1;
  gm->offs[f4_FOFFSET] = -1;
  gm->offs[f4_POFFSET] = -1;
  if (gm->name != NULL) free(gm->name);
  if (gm->acc  != NULL) free(gm->acc);
  if (gm->desc != NULL) free(gm->desc);
  if ((status = esl_strdup(hmm->name,   -1, &(gm->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->acc,    -1, &(gm->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->desc,   -1, &(gm->desc))) != eslOK) goto ERROR;
  if (hmm->flags & f4H_RF)    strcpy(gm->rf,        hmm->rf);
  if (hmm->flags & f4H_MMASK) strcpy(gm->mm,        hmm->mm);
  if (hmm->flags & f4H_CONS)  strcpy(gm->consensus, hmm->consensus); /* must be present, actually, so the flag test is just for symmetry w/ other optional HMM fields */
  if (hmm->flags & f4H_CS)    strcpy(gm->cs,        hmm->cs);
  for (z = 0; z < f4_NEVPARAM; z++) gm->evparam[z] = hmm->evparam[z];
  for (z = 0; z < f4_NCUTOFFS; z++) gm->cutoff[z]  = hmm->cutoff[z];
  for (z = 0; z < f4_MAXABET;  z++) gm->compo[z]   = hmm->compo[z];

  /* Entry scores. */
  if (f4_profile_IsLocal(gm))
    {
      /* Local mode entry:  occ[k] /( \sum_i occ[i] * (M-i+1))
       * (Reduces to uniform 2/(M(M+1)) for occupancies of 1.0)  */
      Z = 0.;
      ESL_ALLOC(occ, sizeof(float) * (hmm->M+1));

      if ((status = f4_hmm_CalculateOccupancy(hmm, occ, NULL)) != eslOK) goto ERROR;
      for (k = 1; k <= hmm->M; k++) 
	Z += occ[k] * (float) (hmm->M-k+1);
      for (k = 1; k <= hmm->M; k++) 
	f4P_TSC(gm, k-1, f4P_BM) = log(occ[k] / Z); /* note off-by-one: entry at Mk stored as [k-1][BM] */

      free(occ);
    }
  else	/* glocal modes: left wing retraction; must be in log space for precision */
    {
      Z = log(hmm->t[0][f4H_MD]);
      f4P_TSC(gm, 0, f4P_BM) = log(1.0 - hmm->t[0][f4H_MD]);
      for (k = 1; k < hmm->M; k++) 
	{
	   f4P_TSC(gm, k, f4P_BM) = Z + log(hmm->t[k][f4H_DM]);
	   Z += log(hmm->t[k][f4H_DD]);
	}
    }

  /* E state loop/move probabilities: nonzero for MOVE allows loops/multihits
   * N,C,J transitions are set later by length config 
   */
  if (f4_profile_IsMultihit(gm)) {
    gm->xsc[f4P_E][f4P_MOVE] = -eslCONST_LOG2;   
    gm->xsc[f4P_E][f4P_LOOP] = -eslCONST_LOG2;   
    gm->nj                   = 1.0f;
  } else {
    gm->xsc[f4P_E][f4P_MOVE] = 0.0f;   
    gm->xsc[f4P_E][f4P_LOOP] = -eslINFINITY;  
    gm->nj                   = 0.0f;
  }

  /* Transition scores. */
  for (k = 1; k < gm->M; k++) {
    tp = gm->tsc + k * f4P_NTRANS;
    tp[f4P_MM] = log(hmm->t[k][f4H_MM]);
    tp[f4P_MI] = log(hmm->t[k][f4H_MI]);
    tp[f4P_MD] = log(hmm->t[k][f4H_MD]);
    tp[f4P_IM] = log(hmm->t[k][f4H_IM]);
    tp[f4P_II] = log(hmm->t[k][f4H_II]);
    tp[f4P_DM] = log(hmm->t[k][f4H_DM]);
    tp[f4P_DD] = log(hmm->t[k][f4H_DD]);
  }
  
  /* Match emission scores. */
  sc[hmm->abc->K]     = -eslINFINITY; /* gap character */
  sc[hmm->abc->Kp-2]  = -eslINFINITY; /* nonresidue character */
  sc[hmm->abc->Kp-1]  = -eslINFINITY; /* missing data character */
  for (k = 1; k <= hmm->M; k++) {
    for (x = 0; x < hmm->abc->K; x++)
     sc[x] = log((double)hmm->mat[k][x] / bg->f[x]);

    esl_abc_FExpectScVec(hmm->abc, sc, bg->f); 

    for (x = 0; x < hmm->abc->Kp; x++) {
      rp = gm->rsc[x] + k * f4P_NR;
      rp[f4P_MSC] = sc[x];
    }
  }

  /* Insert emission scores */
  /* SRE, Fri Dec 5 08:41:08 2008: We currently hardwire insert scores
   * to 0, i.e. corresponding to the insertion emission probabilities
   * being equal to the background probabilities. Benchmarking shows
   * that setting inserts to informative emission distributions causes
   * more problems than it's worth: polar biased composition hits
   * driven by stretches of "insertion" occur, and are difficult to
   * correct for.
   */
  for (x = 0; x < gm->abc->Kp; x++)
    {
      for (k = 1; k < hmm->M; k++) f4P_ISC(gm, k, x) = 0.0f;
      f4P_ISC(gm, hmm->M, x) = -eslINFINITY;   /* init I_M to impossible.   */
    }
  for (k = 1; k <= hmm->M; k++) f4P_ISC(gm, k, gm->abc->K)    = -eslINFINITY; /* gap symbol */
  for (k = 1; k <= hmm->M; k++) f4P_ISC(gm, k, gm->abc->Kp-2) = -eslINFINITY; /* nonresidue symbol */
  for (k = 1; k <= hmm->M; k++) f4P_ISC(gm, k, gm->abc->Kp-1) = -eslINFINITY; /* missing data symbol */

  /* Remaining specials, [NCJ][MOVE | LOOP] are set by ReconfigLength()
   */
  gm->L = 0;			/* force ReconfigLength to reconfig */
  if ((status = f4_ReconfigLength(gm, L)) != eslOK) goto ERROR;
  return eslOK;

 ERROR:
  if (occ != NULL) free(occ);
  return status;
}


/* Function:  f4_ReconfigLength()
 * Synopsis:  Set the target sequence length of a model.
 *
 * Purpose:   Given a model already configured for scoring, in some
 *            particular algorithm mode; reset the expected length
 *            distribution of the profile for a new mean of <L>.
 *
 *            This doesn't affect the length distribution of the null
 *            model. That must also be reset, using <f4_bg_SetLength()>.
 *            
 *            We want this routine to run as fast as possible, because
 *            the caller needs to dynamically reconfigure the model
 *            for the length of each target sequence in a database
 *            search. The profile has precalculated <gm->nj>, 
 *            the number of times the J state is expected to be used,
 *            based on the E state loop transition in the current
 *            configuration.
 *
 * Returns:   <eslOK> on success; xsc[NCJ] scores are set here. These
 *            control the target length dependence of the model.
 */
int
f4_ReconfigLength(F4_PROFILE *gm, int L)
{
  float ploop, pmove;
  
  /* Configure N,J,C transitions so they bear L/(2+nj) of the total
   * unannotated sequence length L. 
   */
  pmove = (2.0f + gm->nj) / ((float) L + 2.0f + gm->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;
  gm->xsc[f4P_N][f4P_LOOP] =  gm->xsc[f4P_C][f4P_LOOP] = gm->xsc[f4P_J][f4P_LOOP] = log(ploop);
  gm->xsc[f4P_N][f4P_MOVE] =  gm->xsc[f4P_C][f4P_MOVE] = gm->xsc[f4P_J][f4P_MOVE] = log(pmove);
  gm->L = L;
  return eslOK;
}
