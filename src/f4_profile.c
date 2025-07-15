#include "dummer.h"

/* Function:  f4_profile_Create()
 * Synopsis:  Allocates a profile.
 *
 * Purpose:   Allocates for a profile of up to <M> nodes, for digital
 *            alphabet <abc>.
 *            
 *            Because this function might be in the critical path (in
 *            hmmscan, for example), we leave much of the model
 *            unintialized, including scores and length model
 *            probabilities. The <f4_ProfileConfig()> call is what
 *            sets these. 
 *            
 *            The alignment mode is set to <f4_NO_MODE>.  The
 *            reference pointer <gm->abc> is set to <abc>.
 *
 * Returns:   a pointer to the new profile.
 *
 * Throws:    <NULL> on allocation error.
 *
 * Xref:      STL11/125.
 */
F4_PROFILE *
f4_profile_Create(int allocM, const ESL_ALPHABET *abc)
{
  F4_PROFILE *gm = NULL;
  int         x;
  int         status;

  /* level 0 */
  ESL_ALLOC(gm, sizeof(F4_PROFILE));
  gm->tsc       = NULL;
  gm->rsc       = NULL;
  gm->rf        = NULL;
  gm->mm        = NULL;
  gm->cs        = NULL;
  gm->consensus = NULL;

  /* level 1 */
  ESL_ALLOC(gm->tsc,       sizeof(float)   * allocM * f4P_NTRANS); 
  ESL_ALLOC(gm->rsc,       sizeof(float *) * abc->Kp);
  ESL_ALLOC(gm->rf,        sizeof(char)    * (allocM+2)); /* yes, +2: each is (0)1..M, +trailing \0  */
  ESL_ALLOC(gm->mm,        sizeof(char)    * (allocM+2));
  ESL_ALLOC(gm->cs,        sizeof(char)    * (allocM+2));
  ESL_ALLOC(gm->consensus, sizeof(char)    * (allocM+2));
  gm->rsc[0] = NULL;
  
  /* level 2 */
  ESL_ALLOC(gm->rsc[0], sizeof(float) * abc->Kp * (allocM+1) * f4P_NR);
  for (x = 1; x < abc->Kp; x++) 
    gm->rsc[x] = gm->rsc[0] + x * (allocM+1) * f4P_NR;

  /* Initialize some edge pieces of memory that are never used,
   * and are only present for indexing convenience.
   */
  esl_vec_FSet(gm->tsc, f4P_NTRANS, -eslINFINITY);     /* node 0 nonexistent, has no transitions  */
  if (allocM > 1) {
    f4P_TSC(gm, 1, f4P_DM) = -eslINFINITY;             /* delete state D_1 is wing-retracted      */
    f4P_TSC(gm, 1, f4P_DD) = -eslINFINITY;
  }
  for (x = 0; x < abc->Kp; x++) {        
    f4P_MSC(gm, 0,      x) = -eslINFINITY;             /* no emissions from nonexistent M_0... */
    f4P_ISC(gm, 0,      x) = -eslINFINITY;             /* or I_0... */
    /* I_M is initialized in profile config, when we know actual M, not just allocated max M   */
  }
  x = esl_abc_XGetGap(abc);	                       /* no emission can emit/score gap characters */
  esl_vec_FSet(gm->rsc[x], (allocM+1)*f4P_NR, -eslINFINITY);
  x = esl_abc_XGetMissing(abc);	                      /* no emission can emit/score missing data characters */
  esl_vec_FSet(gm->rsc[x], (allocM+1)*f4P_NR, -eslINFINITY);

  /* Set remaining info  */
  gm->mode             = f4_NO_MODE;
  gm->L                = 0;
  gm->allocM           = allocM;
  gm->M                = 0;
  gm->max_length       = -1;
  gm->nj               = 0.0f;

  gm->roff             = -1;
  gm->eoff             = -1;
  gm->offs[f4_MOFFSET] = -1;
  gm->offs[f4_FOFFSET] = -1;
  gm->offs[f4_POFFSET] = -1;

  gm->name             = NULL;
  gm->acc              = NULL;
  gm->desc             = NULL;
  gm->rf[0]            = 0;     /* RF line is optional annotation; this flags that it's not set yet */
  gm->mm[0]            = 0;     /* likewise for MM annotation line */
  gm->cs[0]            = 0;     /* likewise for CS annotation line */
  gm->consensus[0]     = 0;
  
  for (x = 0; x < f4_NEVPARAM; x++) gm->evparam[x] = f4_EVPARAM_UNSET;
  for (x = 0; x < f4_NCUTOFFS; x++) gm->cutoff[x]  = f4_CUTOFF_UNSET;
  for (x = 0; x < f4_MAXABET;  x++) gm->compo[x]   = f4_COMPO_UNSET;

  gm->abc         = abc;
  return gm;

 ERROR:
  f4_profile_Destroy(gm);
  return NULL;
}

/* Function:  f4_profile_Destroy()
 * Synopsis:  Frees a profile.
 *
 * Purpose:   Frees a profile <gm>.
 *
 * Returns:   (void).
 */
void
f4_profile_Destroy(F4_PROFILE *gm)
{
  if (gm != NULL) {
    if (gm->rsc   != NULL && gm->rsc[0] != NULL) free(gm->rsc[0]);
    if (gm->tsc       != NULL) free(gm->tsc);
    if (gm->rsc       != NULL) free(gm->rsc);
    if (gm->name      != NULL) free(gm->name);
    if (gm->acc       != NULL) free(gm->acc);
    if (gm->desc      != NULL) free(gm->desc);
    if (gm->rf        != NULL) free(gm->rf);
    if (gm->mm        != NULL) free(gm->mm);
    if (gm->cs        != NULL) free(gm->cs);
    if (gm->consensus != NULL) free(gm->consensus);
    free(gm);
  }
  return;
}

/* Function:  f4_profile_IsLocal()
 * Synopsis:  Return TRUE if profile is in a local alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a local alignment mode.
 */
int
f4_profile_IsLocal(const F4_PROFILE *gm)
{
  if (gm->mode == f4_UNILOCAL || gm->mode == f4_LOCAL) return TRUE;
  return FALSE;
}

/* Function:  f4_profile_IsMultihit()
 * Synopsis:  Return TRUE if profile is in a multihit alignment mode.
 *
 * Purpose:   Return <TRUE> if profile is in a multihit alignment mode.
 */
int
f4_profile_IsMultihit(const F4_PROFILE *gm)
{
  if (gm->mode == f4_LOCAL || gm->mode == f4_GLOCAL) return TRUE;
  return FALSE;
}




/* Function:  f4_profile_GetT()
 *
 * Purpose:   Convenience function that looks up a transition score in
 *            profile <gm> for a transition from state type <st1> in
 *            node <k1> to state type <st2> in node <k2>. For unique
 *            state types that aren't in nodes (<f4T_S>, for example), the
 *            <k> value is ignored, though it would be customarily passed as 0.
 *            Return the transition score in <ret_tsc>.
 *            
 *            This function would almost always be called on profile
 *            traces, of course, but it's possible to call it
 *            on core traces (for example, if you were to try to 
 *            trace_Dump() during HMM construction, and you wanted
 *            to see detailed profile scores for that trace). Core traces
 *            can contain <f4T_X> "states" used solely to signal
 *            a sequence fragment, treated as missing data. Transitions
 *            involving <f4T_X> states are assigned zero score here.
 *            Other transitions that occur only in core traces
 *            (B->I0, B->D1, I_M->E) also silently get a zero score.
 *            This is safe, because we would only ever use this number
 *            for display, not as a log probability somewhere.
 *
 * Returns:   <eslOK> on success, and <*ret_tsc> contains the requested
 *            transition score.            
 * 
 * Throws:    <eslEINVAL> if a nonexistent transition is requested. Now
 *            <*ret_tsc> is set to $-\infty$.
 *            
 */
int
f4_profile_GetT(const F4_PROFILE *gm, char st1, int k1, char st2, int k2, float *ret_tsc)
{
  float tsc = 0.0f;
  int   status;

  /* Detect transitions that can only come from core traces;
   * return 0.0 as a special case (this is only done for displaying
   * "scores" in trace dumps, during debugging.)
   */
  if (st1 == f4T_X || st2 == f4T_X) return eslOK;
  if (st1 == f4T_B && st2 == f4T_I) return eslOK;
  if (st1 == f4T_B && st2 == f4T_D) return eslOK;
  if (st1 == f4T_I && st2 == f4T_E) return eslOK;

  /* Now we're sure this is a profile trace, as it should usually be. */
  switch (st1) {
  case f4T_S:  break;
  case f4T_T:  break;
  case f4T_N:
    switch (st2) {
    case f4T_B: tsc =  gm->xsc[f4P_N][f4P_MOVE]; break;
    case f4T_N: tsc =  gm->xsc[f4P_N][f4P_LOOP]; break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", f4_hmm_DecodeStatetype(st1), f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_B:
    switch (st2) {
    case f4T_M: tsc = f4P_TSC(gm, k2-1, f4P_BM); break; /* remember, B->Mk is stored in [k-1][f4P_BM] */
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", f4_hmm_DecodeStatetype(st1), f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_M:
    switch (st2) {
    case f4T_M: tsc = f4P_TSC(gm, k1, f4P_MM); break;
    case f4T_I: tsc = f4P_TSC(gm, k1, f4P_MI); break;
    case f4T_D: tsc = f4P_TSC(gm, k1, f4P_MD); break;
    case f4T_E: 
      if (k1 != gm->M && ! f4_profile_IsLocal(gm)) ESL_EXCEPTION(eslEINVAL, "local end transition (M%d of %d) in non-local model", k1, gm->M);
      tsc = 0.0f;		/* by def'n in H3 local alignment */
      break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", f4_hmm_DecodeStatetype(st1), k1, f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_D:
    switch (st2) {
    case f4T_M: tsc = f4P_TSC(gm, k1, f4P_DM); break;
    case f4T_D: tsc = f4P_TSC(gm, k1, f4P_DD); break;
    case f4T_E: 
      if (k1 != gm->M && ! f4_profile_IsLocal(gm)) ESL_EXCEPTION(eslEINVAL, "local end transition (D%d of %d) in non-local model", k1, gm->M);
      tsc = 0.0f;		/* by def'n in H3 local alignment */
      break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", f4_hmm_DecodeStatetype(st1), k1, f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_I:
    switch (st2) {
    case f4T_M: tsc = f4P_TSC(gm, k1, f4P_IM); break;
    case f4T_I: tsc = f4P_TSC(gm, k1, f4P_II); break;
    default:    ESL_XEXCEPTION(eslEINVAL, "bad transition %s_%d->%s", f4_hmm_DecodeStatetype(st1), k1, f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_E:
    switch (st2) {
    case f4T_C: tsc = gm->xsc[f4P_E][f4P_MOVE]; break;
    case f4T_J: tsc = gm->xsc[f4P_E][f4P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", f4_hmm_DecodeStatetype(st1), f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_J:
    switch (st2) {
    case f4T_B: tsc = gm->xsc[f4P_J][f4P_MOVE]; break;
    case f4T_J: tsc = gm->xsc[f4P_J][f4P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", f4_hmm_DecodeStatetype(st1), f4_hmm_DecodeStatetype(st2));
    }
    break;

  case f4T_C:
    switch (st2) {
    case f4T_T:  tsc = gm->xsc[f4P_C][f4P_MOVE]; break;
    case f4T_C:  tsc = gm->xsc[f4P_C][f4P_LOOP]; break;
    default:     ESL_XEXCEPTION(eslEINVAL, "bad transition %s->%s", f4_hmm_DecodeStatetype(st1), f4_hmm_DecodeStatetype(st2));
    }
    break;

  default: ESL_XEXCEPTION(eslEINVAL, "bad state type %d in traceback", st1);
  }

  *ret_tsc = tsc;
  return eslOK;

 ERROR:
  *ret_tsc = -eslINFINITY;
  return status;
}