/* F4_TRACE, the traceback structure.
 *
 * Contents:
 *   1. The F4_TRACE structure
 *   2. Debugging tools
 *   3. Creating traces by DP traceback
 *   4. Creating faux traces from existing MSAs
 *   5. Counting traces into new HMMs
 * 
 * Stylistic note: elements in a trace path are usually indexed by z.
 */

#include "dummer.h"


/*****************************************************************
 * 1. The F4_TRACE structure.
 *****************************************************************/

static F4_TRACE *trace_create_engine(int initial_nalloc, int initial_ndomalloc, int with_posteriors);

/* Function:  f4_trace_Create()
 * Synopsis:  Allocates a (growable, reusable) traceback.
 *
 * Purpose:   Allocates a traceback. 
 *  
 *            Tracebacks are growable. A reasonable initial internal
 *            allocation is made here, and routines that generate
 *            tracebacks will dynamically grow the trace as needed.
 *            
 *            Tracebacks are reusable. Usually a routine only
 *            allocates one, and reuses its memory over and over as
 *            new target sequences are aligned.
 *
 * Returns:   a pointer to the new <F4_TRACE> structure on success.
 *
 * Throws:    <NULL> on allocation error.
 */
F4_TRACE *
f4_trace_Create(void)
{
  int       initial_nalloc    = 256;
  int       initial_ndomalloc = 16;
  int       with_posteriors   = FALSE;
  return trace_create_engine(initial_nalloc, initial_ndomalloc, with_posteriors);
}

static F4_TRACE *
trace_create_engine(int initial_nalloc, int initial_ndomalloc, int with_posteriors)
{
  F4_TRACE *tr      = NULL;
  int       status;

  ESL_ALLOC(tr, sizeof(F4_TRACE));
  tr->st = NULL;
  tr->k  = NULL;
  tr->i  = NULL;
  tr->pp = NULL;
  tr->M  = 0;
  tr->L  = 0;
  tr->tfrom   = tr->tto   = NULL;
  tr->sqfrom  = tr->sqto  = NULL;
  tr->hmmfrom = tr->hmmto = NULL;

  /* The trace data itself */
  ESL_ALLOC(tr->st, sizeof(char) * initial_nalloc);
  ESL_ALLOC(tr->k,  sizeof(int)  * initial_nalloc);
  ESL_ALLOC(tr->i,  sizeof(int)  * initial_nalloc);
  if (with_posteriors)
    ESL_ALLOC(tr->pp, sizeof(float) * initial_nalloc);
  tr->N      = 0;
  tr->nalloc = initial_nalloc;

  /* The trace's index: table of domain start/stop coords */
  ESL_ALLOC(tr->tfrom,   sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->tto,     sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqfrom,  sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->sqto,    sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->hmmfrom, sizeof(int) * initial_ndomalloc);
  ESL_ALLOC(tr->hmmto,   sizeof(int) * initial_ndomalloc);
  tr->ndom      = 0;
  tr->ndomalloc = initial_ndomalloc;
  return tr;

 ERROR:
  if (tr != NULL) f4_trace_Destroy(tr);
  return NULL;
}

/* Function:  f4_trace_Grow()
 * Synopsis:  Grow the allocation for trace data.
 *
 * Purpose:   If <tr> can't fit another state, double its allocation for
 *            traceback data.
 *            
 *            This doesn't reallocate the domain index; see
 *            <f4_trace_GrowIndex()> or <f4_trace_GrowIndexTo()> for
 *            that.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the data in
 *            <tr> are unaffected.
 */
int
f4_trace_Grow(F4_TRACE *tr)
{
  void *tmp;
  int   status;
  
  if (tr->N < tr->nalloc) return eslOK;

  ESL_RALLOC(tr->st, tmp, sizeof(char) *2*tr->nalloc);
  ESL_RALLOC(tr->k,  tmp, sizeof(int)  *2*tr->nalloc);
  ESL_RALLOC(tr->i,  tmp, sizeof(int)  *2*tr->nalloc);
  if (tr->pp != NULL) ESL_RALLOC(tr->pp,  tmp, sizeof(float) *2*tr->nalloc);
  tr->nalloc *= 2;
  return eslOK;

 ERROR:
  return status;
}

/* Function:  f4_trace_Destroy()
 * Synopsis:  Frees a trace.
 *
 * Purpose:   Frees a trace structure <tr>.
 *
 * Returns:   (void)
 */
void 
f4_trace_Destroy(F4_TRACE *tr)
{
  if (tr == NULL) return;
  if (tr->st      != NULL) free(tr->st);
  if (tr->k       != NULL) free(tr->k);
  if (tr->i       != NULL) free(tr->i);
  if (tr->pp      != NULL) free(tr->pp);
  if (tr->tfrom   != NULL) free(tr->tfrom);
  if (tr->tto     != NULL) free(tr->tto);
  if (tr->sqfrom  != NULL) free(tr->sqfrom);
  if (tr->sqto    != NULL) free(tr->sqto);
  if (tr->hmmfrom != NULL) free(tr->hmmfrom);
  if (tr->hmmto   != NULL) free(tr->hmmto);
  free(tr);
  return;
}

/* Function:  f4_trace_DestroyArray()
 *
 * Purpose:   Frees an array of <N> trace structures, <tr[0..N-1]>.
 *
 * Returns:   (void)
 */
void 
f4_trace_DestroyArray(F4_TRACE **tr, int N)
{
  int idx;

  if (tr == NULL) return;
  for (idx = 0; idx < N; idx++)
    {
      if (tr[idx] == NULL) continue;
      f4_trace_Destroy(tr[idx]);
    }
  free(tr);
  return;
}

/*---------------------- end, F4_TRACE --------------------------*/

/*****************************************************************
 * 2. Debugging tools.
 *****************************************************************/

/* Function:  f4_trace_Validate()
 *
 * Purpose:   Validate the internal data in a trace structure <tr>
 *            representing an alignment of an HMM to a 
 *            digital sequence <sq>. The digital sequence may be either
 *            unaligned (usually) or aligned (in the case of "fake"
 *            tracebacks generated from an MSA during a
 *            model construction process). 
 *            
 *            We don't pass the HMM that the trace is associated with,
 *            because we might have constructed the trace during
 *            HMM construction when we don't have an HMM yet; but 
 *            we always have a digital sequence.
 *
 *            Intended for debugging, development, and testing
 *            purposes.
 *            
 * Args:      tr     - trace to validate
 *            abc    - alphabet corresponding to sequence <sq>
 *            sq     - digital sequence that <tr> is explaining
 *            errbuf - NULL, or an error message buffer allocated
 *                     for at least eslERRBUFSIZE chars.           
 *
 * Returns:   <eslOK> if trace appears fine.
 *            Returns <eslFAIL> if a problem is detected; if <errbuf> is
 *            provided (non-<NULL>), an informative message is formatted
 *            there to indicate the reason for the failure.
 */
int
f4_trace_Validate(const F4_TRACE *tr, const ESL_ALPHABET *abc, const ESL_DSQ *dsq, char *errbuf)
{
  int  z;			/* position in trace    */
  int  i;			/* position in sequence */
  int  k;			/* position in model */
  char prv;			/* type of the previous state */
  int  is_core;			/* TRUE if trace is a core trace, not profile */

  /* minimum trace length is a core's B->Mk->E. If we don't have at least that,
   * we're definitely in trouble
   */
  if (tr->N < 3)          ESL_FAIL(eslFAIL, errbuf, "trace is too short");
  if (tr->N > tr->nalloc) ESL_FAIL(eslFAIL, errbuf, "N of %d isn't sensible", tr->N);

  /* Determine if this is a core trace or a profile trace, so we can
   * construct validation tests appropriately.
   */
  if      (tr->st[0] == f4T_B) is_core = TRUE;
  else if (tr->st[0] == f4T_S) is_core = FALSE;
  else    ESL_FAIL(eslFAIL, errbuf, "first state neither S nor B");

  /* Verify "sentinels", the final states of the trace
   * (before we start looking backwards and forwards from each state in 
   * our main validation loop)
   */
  if (is_core  && tr->st[tr->N-1] != f4T_E) ESL_FAIL(eslFAIL, errbuf, "last state not E");
  if (!is_core && tr->st[tr->N-1] != f4T_T) ESL_FAIL(eslFAIL, errbuf, "last state not T");
  if (tr->k[0]        != 0)                 ESL_FAIL(eslFAIL, errbuf, "first state shouldn't have k set");
  if (tr->i[0]        != 0)                 ESL_FAIL(eslFAIL, errbuf, "first state shouldn't have i set");
  if (tr->k[tr->N-1]  != 0)                 ESL_FAIL(eslFAIL, errbuf, "last state shouldn't have k set");
  if (tr->i[tr->N-1]  != 0)                 ESL_FAIL(eslFAIL, errbuf, "last state shouldn't have i set");

  if (tr->pp != NULL && tr->pp[0]       != 0.0) ESL_FAIL(eslFAIL, errbuf, "first state doesn't emit; but post prob isn't 0");
  if (tr->pp != NULL && tr->pp[tr->N-1] != 0.0) ESL_FAIL(eslFAIL, errbuf, "last state doesn't emit; but post prob isn't 0");

  /* Main validation loop. */
  k = 0; 
  i = 1;
  for (z = 1; z < tr->N-1; z++)
    {
      for (; dsq[i] != eslDSQ_SENTINEL; i++) /* find next non-gap residue in dsq */
	if (esl_abc_XIsResidue(abc, dsq[i]) || esl_abc_XIsNonresidue(abc, dsq[i])) break; /* '*' included as emitted "residue"  */

      /* watch out for missing data states X: can only be one.
       * prv state might have to skip over one (but not more) missing data states
       */
      prv = (tr->st[z-1] == f4T_X)? tr->st[z-2] : tr->st[z-1];

      switch (tr->st[z]) {
      case f4T_S:
	ESL_FAIL(eslFAIL, errbuf, "S must be first state");
	break;
	
      case f4T_X:
	if (! is_core)       ESL_FAIL(eslFAIL, errbuf, "X state (missing data) only appears in core traces");
	if (prv != f4T_B && tr->st[z+1] != f4T_E)	/* only B->X and X->E are possible */
	  ESL_FAIL(eslFAIL, errbuf, "bad transition involving missing data (X state) not at start/end");
	break;

      case f4T_N:
	if (is_core)       ESL_FAIL(eslFAIL, errbuf, "core trace can't contain N");
	if (tr->k[z] != 0) ESL_FAIL(eslFAIL, errbuf, "no N should have k set");
	if (prv == f4T_S) { /* 1st N doesn't emit */
	  if (tr->i[z] != 0)                      ESL_FAIL(eslFAIL, errbuf, "first N shouldn't have i set");
	  if (tr->pp != NULL && tr->pp[z] != 0.0) ESL_FAIL(eslFAIL, errbuf, "first N can't have nonzero post prob");
	} else if (prv == f4T_N) { /* subsequent N's do */
	  if (tr->i[z] != i) ESL_FAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_FAIL(eslFAIL, errbuf, "bad transition to N; expected {S,N}->N");
	break;

      case f4T_B:
	if (tr->k[z] != 0)                      ESL_FAIL(eslFAIL, errbuf, "B shouldn't have k set");
	if (tr->i[z] != 0)                      ESL_FAIL(eslFAIL, errbuf, "B shouldn't have i set");
	if (tr->pp != NULL && tr->pp[z] != 0.0) ESL_FAIL(eslFAIL, errbuf, "B can't have nonzero post prob");
	if (prv != f4T_N && prv != f4T_J) 
	  ESL_FAIL(eslFAIL, errbuf, "bad transition to B; expected {N,J}->B");
	break;

      case f4T_M:
	if (prv == f4T_B) k = tr->k[z]; else k++; /* on a B->Mk entry, trust k; else verify */

	if (tr->k[z] != k) ESL_FAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != i) ESL_FAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	if (prv != f4T_B && prv != f4T_M && prv != f4T_D && prv != f4T_I)
	  ESL_FAIL(eslFAIL, errbuf, "bad transition to M; expected {B,M,D,I}->M");
	i++;
	break;

      case f4T_D:
	k++;
	if (tr->st[z-1] == f4T_X)  k = tr->k[z]; /* with fragments, a X->Ik case is possible */
	if (tr->k[z] != k)                      ESL_FAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != 0)                      ESL_FAIL(eslFAIL, errbuf, "D shouldn't have i set");
	if (tr->pp != NULL && tr->pp[z] != 0.0) ESL_FAIL(eslFAIL, errbuf, "D can't have nonzero post prob");
	if (is_core) {
	  if (prv != f4T_M && prv != f4T_D && prv != f4T_B && prv != f4T_I)
	    ESL_FAIL(eslFAIL, errbuf, "bad transition to D; expected {B,M,D,I}->D");
	} else {
	  if (prv != f4T_M && prv != f4T_D && prv != f4T_I)
	    ESL_FAIL(eslFAIL, errbuf, "bad transition to D; expected {M,D,I}->D");
	}
	break;
	
      case f4T_I:
	if (tr->st[z-1] == f4T_X)  k = tr->k[z]; /* with fragments, a X->Ik case is possible */
	if (tr->k[z] != k) ESL_FAIL(eslFAIL, errbuf, "expected k doesn't match trace's k");
	if (tr->i[z] != i) ESL_FAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	if (is_core) {
	  if (prv != f4T_B && prv != f4T_M && prv != f4T_I && prv != f4T_D)
	    ESL_FAIL(eslFAIL, errbuf, "bad transition to I; expected {B,M,D,I}->I");
	} else {
	  if (prv != f4T_M && prv != f4T_I && prv != f4T_D)
	    ESL_FAIL(eslFAIL, errbuf, "bad transition to I; expected {M,D,I}->I");
	}
	i++;
	break;

      case f4T_E:
	if (tr->k[z] != 0) ESL_FAIL(eslFAIL, errbuf, "E shouldn't have k set");
	if (tr->i[z] != 0) ESL_FAIL(eslFAIL, errbuf, "E shouldn't have i set");
	if (tr->pp != NULL && tr->pp[z] != 0.0) ESL_FAIL(eslFAIL, errbuf, "E can't have nonzero post prob");
	if (is_core) {
	  if (prv != f4T_M && prv != f4T_D && prv != f4T_I)
	    ESL_FAIL(eslFAIL, errbuf, "bad transition to E; expected {M,D,I}->E");
	} else {
	  if (prv != f4T_M && prv != f4T_D)
	    ESL_FAIL(eslFAIL, errbuf, "bad transition to E; expected {M,D}->E");
	}
	break;
	
      case f4T_J:
	if (tr->k[z] != 0) ESL_FAIL(eslFAIL, errbuf, "no J should have k set");
	if (prv == f4T_E) { /* 1st J doesn't emit */
	  if (tr->i[z] != 0)                      ESL_FAIL(eslFAIL, errbuf, "first J shouldn't have i set");
	  if (tr->pp != NULL && tr->pp[z] != 0.0) ESL_FAIL(eslFAIL, errbuf, "first J can't have nonzero post prob");
	} else if (prv == f4T_J) { /* subsequent J's do */
	  if (tr->i[z] != i) ESL_FAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_FAIL(eslFAIL, errbuf, "bad transition to J; expected {E,J}->J");
	break;

      case f4T_C:
	if (is_core)       ESL_FAIL(eslFAIL, errbuf, "core trace can't contain C");
	if (tr->k[z] != 0) ESL_FAIL(eslFAIL, errbuf, "no C should have k set");
	if (prv == f4T_E) { /* 1st C doesn't emit */
	  if (tr->i[z] != 0)                      ESL_FAIL(eslFAIL, errbuf, "first C shouldn't have i set");
	  if (tr->pp != NULL && tr->pp[z] != 0.0) ESL_FAIL(eslFAIL, errbuf, "first C can't have nonzero post prob");
	} else if (prv == f4T_C) { /* subsequent C's do */
	  if (tr->i[z] != i) ESL_FAIL(eslFAIL, errbuf, "expected i doesn't match trace's i");
	  i++;
	} else ESL_FAIL(eslFAIL, errbuf, "bad transition to C; expected {E,C}->C");
	break;
	
      case f4T_T:
	ESL_FAIL(eslFAIL, errbuf, "T must be last state");
	break;	
      }
    }

  /* Trace should have accounted for all residues in the dsq */
  for (; dsq[i] != eslDSQ_SENTINEL; i++) 
    if (esl_abc_XIsResidue(abc, dsq[i])) 
      ESL_FAIL(eslFAIL, errbuf, "trace didn't account for all residues in the sq");

  /* No k larger than M; no i-1 larger than L (i is sitting on dsq[n+1] sentinel right now) */
  if (k   > tr->M) ESL_FAIL(eslFAIL, errbuf, "M=%d, but k went to %d\n", tr->M, k);
  if (i-1 > tr->L) ESL_FAIL(eslFAIL, errbuf, "L=%d, but i went to %d\n", tr->L, i);

  return eslOK;
}

/*------------------ end, debugging tools -----------------------*/

/*****************************************************************
 * 3. Creating traces by DP traceback
 *****************************************************************/

/* Function:  f4_trace_Append()
 * Synopsis:  Add an element (state/residue) to a growing trace.
 *
 * Purpose:   Adds an element to a trace <tr> that is growing
 *            left-to-right. The element is defined by a state type
 *            <st> (such as <f4T_M>); a node index <k> (1..M for
 *            M,D,I main states; else 0); and a dsq position <i> (1..L
 *            for emitters, else 0).
 *            
 *            For CNJ states, which emit on transition, by convention
 *            we associate the emission with the downstream state; therefore
 *            the first state in any run of CNJ states has i=0. 
 *            
 *            Reallocates the trace (by doubling) if necessary.
 *            
 *            Caller can grow a trace right-to-left too, if it
 *            plans to call <f4_trace_Reverse()>. 
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on reallocation failure. The element is successfully
 *            added, but no more elements can be added before this trace is
 *            destroyed.
 *            
 *            <eslEINVAL> if you try to add an element to a trace whose
 *            reallocation has already failed.
 */
int
f4_trace_Append(F4_TRACE *tr, char st, int k, int i)
{
  int status;

  if ((status = f4_trace_Grow(tr)) != eslOK) return status;

  switch (st) {
    /* Emit-on-transition states: */
  case f4T_N: 
  case f4T_C: 
  case f4T_J: 
    tr->i[tr->N] = ( (tr->st[tr->N-1] == st) ? i : 0);
    tr->k[tr->N] = 0;
    break;
    /* Nonemitting states, outside main model: */
  case f4T_X:
  case f4T_S:
  case f4T_B:
  case f4T_E:
  case f4T_T: tr->i[tr->N] = 0; tr->k[tr->N] = 0; break;
    /* Nonemitting, but in main model (k valid) */
  case f4T_D: tr->i[tr->N] = 0; tr->k[tr->N] = k; break;
    /* Emitting states, with valid k position in model: */
  case f4T_M: 
  case f4T_I: tr->i[tr->N] = i; tr->k[tr->N] = k; break;
  default:    ESL_EXCEPTION(eslEINVAL, "no such state; can't append");
  }

  tr->st[tr->N] = st;
  tr->N++;
  return eslOK;
}

/*----------- end, creating traces by DP traceback ---------------*/


/*****************************************************************
 * 4. Creating faux traces from MSAs
 *****************************************************************/

/* Function:  f4_trace_FauxFromMSA()
 * Synopsis:  Create array of faux tracebacks from an existing MSA.
 *
 * Purpose:   Given an existing <msa> and an array <matassign> that
 *            flags the alignment columns that are assigned to consensus
 *            match states (matassign[1..alen] = 1|0); create an array
 *            of faux traces <tr[0..msa->nseq-1]>. <optflags> controls 
 *            optional behavior; it can be <f4_DEFAULT> or <f4_MSA_COORDS>,
 *            as explained below.
 *            
 *            The traces are core traces: they start/end with B/E,
 *            they may use I_0,I_M, and D_1 states. Any flanking
 *            insertions (outside the first/last consensus column) are
 *            assigned to I_0 and I_M.
 *            
 *            If the input alignment contains sequence fragments,
 *            caller should first convert leading/trailing gaps to
 *            missing data symbols. This hack causes entry/exit
 *            transitions to be encoded in the trace as B->X->{MDI}k
 *            and {MDI}k->X->E, rather than B->DDDD->Mk, Mk->DDDDD->E
 *            paths involving terminal deletions, and all functions
 *            that use traces, such as <f4_trace_Count()>, (should)
 *            ignore transitions involving <f4T_X> states.
 *            
 *            By default (<optflags = f4_DEFAULT>), the <i> coordinate
 *            in the faux tracebacks is <1..L>, relative to the
 *            unaligned raw sequences in <msa>, the way most H3 traces
 *            are supposed to be. In some cases (such as model
 *            construction from an MSA) it is convenient to reference
 *            residues in the MSA cooordinate system directly; setting
 *            <optflags = f4_MSA_COORDS> makes the traces come out
 *            with <i=1..alen> coords for residues.
 *            
 *            Important: an MSA may imply DI and ID transitions that
 *            are illegal in a core model. If the only purpose of the
 *            traces is to go straight back into alignment
 *            construction through a <f4_tracealign_*> function, this
 *            is ok, because the <f4_tracealign_*> routines can handle
 *            DI and ID transitions (enabling reconstruction of almost
 *            exactly the same input alignment, modulo unaligned
 *            insertions). This is what happens for <hmmalign
 *            --mapali>, for example. However, if the caller wants to
 *            use the traces for anything else, these illegal DI and
 *            ID transitions have to be removed first, and the caller
 *            should use <f4_trace_Doctor()> to do it.
 *
 * Args:      msa       - digital alignment
 *            matassign - flag for each alignment column, whether
 *                        it is consensus or not. matassign[1..alen] = 1|0; 
 *                        matassign[0] = 0
 *            optflags  - f4_DEFAULT | f4_MSA_COORDS 
 *            tr        - RETURN: caller provides 0..nseq-1 pointer 
 *                        array for holding returned traces.
 *
 * Returns:   <eslOK> on success, and tr[0..nseq-1] now point to newly
 *            created traces; caller is responsible for freeing these.
 *
 * Throws:    <eslEMEM> on allocation error.
 */
int
f4_trace_FauxFromMSA(ESL_MSA *msa, int *matassign, int optflags, F4_TRACE **tr)
{		      
  int  idx;			/* counter over seqs in MSA */
  int  k;                       /* position in HMM                 */
  int  apos;                    /* position in alignment columns 1..alen */
  int  rpos;			/* position in unaligned sequence residues 1..L */
  int  showpos;			/* coord to actually record: apos or rpos */
  int  status = eslOK;
 
  for (idx = 0; idx < msa->nseq; idx++) tr[idx] = NULL;
 
  for (idx = 0; idx < msa->nseq; idx++)
    {
      if ((tr[idx] = f4_trace_Create())                      == NULL) goto ERROR; 
      if ((status  = f4_trace_Append(tr[idx], f4T_B, 0, 0)) != eslOK) goto ERROR;

      for (k = 0, rpos = 1, apos = 1; apos <= msa->alen; apos++)
	{
	  showpos = (optflags & f4_MSA_COORDS) ? apos : rpos;

	  if (matassign[apos]) 
	    {			/* match or delete */
	      k++;
	      if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) 
		status = f4_trace_Append(tr[idx], f4T_M, k, showpos);
	      else if (esl_abc_XIsGap    (msa->abc, msa->ax[idx][apos])) 
		status = f4_trace_Append(tr[idx], f4T_D, k, 0);          
	      else if (esl_abc_XIsNonresidue(msa->abc, msa->ax[idx][apos]))
		status = f4_trace_Append(tr[idx], f4T_M, k, showpos); /* treat * as a residue! */
	      else if (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
		{
		  if (tr[idx]->st[tr[idx]->N-1] != f4T_X)
		    status = f4_trace_Append(tr[idx], f4T_X, k, 0); /* allow only one X in a row */
		}
	      else ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen");
	    }
	  else
	    { 			/* insert or nothing */
	      if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos]))
		status = f4_trace_Append(tr[idx], f4T_I, k, showpos);
	      else if (esl_abc_XIsNonresidue(msa->abc, msa->ax[idx][apos]))
		status = f4_trace_Append(tr[idx], f4T_I, k, showpos); /* treat * as a residue! */
	      else if (esl_abc_XIsMissing(msa->abc, msa->ax[idx][apos]))
		{ 
		  if (tr[idx]->st[tr[idx]->N-1] != f4T_X)
		    status = f4_trace_Append(tr[idx], f4T_X, k, 0);
		}
	      else if (! esl_abc_XIsGap(msa->abc, msa->ax[idx][apos]))
		ESL_XEXCEPTION(eslEINCONCEIVABLE, "can't happen");
	    }

	  if (esl_abc_XIsResidue(msa->abc, msa->ax[idx][apos])) rpos++; 
	  if (status != eslOK) goto ERROR;
	}
      if ((status = f4_trace_Append(tr[idx], f4T_E, 0, 0)) != eslOK) goto ERROR;
      /* k == M by construction; set tr->L = msa->alen since coords are w.r.t. ax */
      tr[idx]->M = k;
      tr[idx]->L = msa->alen;
    }
  return eslOK;


 ERROR:
  for (idx = 0; idx < msa->nseq; idx++) { f4_trace_Destroy(tr[idx]); tr[idx] = NULL; }
  return status; 
}

/*-------------- end, faux traces from MSAs ---------------------*/


/*****************************************************************
 * 5. Counting traces into new HMMs.
 *****************************************************************/

/* Function: f4_trace_Count()
 * 
 * Purpose:  Count a traceback into a count-based core HMM structure.
 *           (Usually as part of a model parameter re-estimation.)
 *           This is a special version of <p7_trace_Count()> but for f4-HMMs,
 *           which allows D->I and I->D transitions, and does not count
 *           based on transition occurrences, but on the actual parameters,
 *           alpha, beta, delta and epsilon that are used to later
 *           estimate the model parameters (transitions).
 * 
 *           The counts are usually only used as a starting point to be
 *           later refined by an EM algorithm, such as Baum-Welch.
 *           
 *           The traceback may either be a core traceback (as in model
 *           construction) or a profile traceback (as in model
 *           reestimation).
 *           
 *           If it is a profile traceback, we have to be careful how
 *           we translate an internal entry path from a score profile
 *           back to the core model. Sometimes a B->M_k transition is
 *           an internal entry from local alignment, and sometimes it
 *           is a wing-folded B->D_1..DDM_k alignment to the core
 *           model.
 *           
 *           This is one of the purposes of the special f4T_X
 *           'missing data' state in tracebacks. Local alignment entry
 *           is indicated by a B->X->{MDI}_k 'missing data' path, and
 *           direct B->M_k or M_k->E transitions in a traceback are
 *           interpreted as wing retraction in a glocal model.
 * 
 *           The <f4T_X> state is also used in core traces in model
 *           construction literally to mean missing data, in the
 *           treatment of sequence fragments.
 *
 * Args:     hmm   - counts-based HMM to count <tr> into
 *           tr    - alignment of seq to HMM
 *           dsq   - digitized sequence that traceback aligns to the HMM (1..L)
 *                   (or can be an ax, aligned digital seq)
 *           wt    - weight on this sequence
 *           
 * Return:   <eslOK> on success.
 *           Weighted count events are accumulated in hmm's mat[][], ins[][],
 *           tp[][] fields: the core probability model.
 *           
 * Throws:   <eslEINVAL> if something's corrupt in the trace; effect on hmm
 *           counts is undefined, because it may abort at any point in the trace.
 */
int
f4_trace_Count(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr)
{
  int z;			/* position index in trace */
  int i;			/* symbol position in seq */
  int st,st2;     		/* state type (cur, nxt)  */
  int k,k2,ktmp;		/* node index (cur, nxt)  */
  int z1 = 0;			/* left bound - may get set to an M position for a left fragment */
  int z2 = tr->N-1;		/* right bound, ditto for a right fragment. N-1 not N, because main loop accesses z,z+1 */
  
  /* If this is a core fragment trace (it has B->X and/or X->E) then
   * set z1 and/or z2 bound on first and/or last M state, so we don't
   * count incomplete flanking insertions. A fragment doesn't
   * necessarily have X's on both sides because of the way they get
   * set from ~'s in an input alignment.
   * 
   * A local alignment profile trace has B->X and X->E, and may have
   * >1 domain, but is guaranteed to be B->X->Mk, Mk->X->E, so
   * limiting trace counting to z1..z2 would have no effect... nonetheless,
   * we check, differentiating core vs. profile trace by the lead B vs S.
   * 
   * It's possible for a core trace to have no M's at all, just
   * B->(X)->III->(X)->E, as in bug #h82, so watch out for that; we don't
   * count anything in such a trace, even the II transitions, because
   * we don't get to see the complete length of the insertion (or the
   * IM transition), so we don't want to be estimating the I-state
   * geometric distribution from it.
   */
  if (tr->st[0] == f4T_B && tr->st[1] == f4T_X)
    for (z = 2; z < tr->N-1; z++)
      if (tr->st[z] == f4T_M) { z1 = z; break; }
  if (tr->st[tr->N-1] == f4T_E && tr->st[tr->N-2] == f4T_X)
    for (z = tr->N-3; z > 0; z--)
      if (tr->st[z] == f4T_M) { z2 = z; break; }
  
  for (z = z1; z < z2; z++) 
    {
      if (tr->st[z] == f4T_X) continue; /* skip missing data */

      /* pull some info into tmp vars for notational clarity later. */
      st  = tr->st[z]; 
      st2 = tr->st[z+1];
      k   = tr->k[z]; 
      k2  = tr->k[z+1];
      i   = tr->i[z];

      /* Emission counts. */
      if      (st == f4T_M) esl_abc_FCount(hmm->abc, hmm->mat[k], dsq[i], wt);
      else if (st == f4T_I) esl_abc_FCount(hmm->abc, hmm->ins[k], dsq[i], wt);

      /* Transition counts */
      if (st2 == f4T_X) continue; /* ignore transition to missing data */

      if (st == f4T_B) {
	if (st2 == f4T_M && k2 > 1)   /* wing-retracted B->DD->Mk path */
	  {
	    hmm->tp[0][f4H_DELTA] += wt; //MD
	    for (ktmp = 1; ktmp < k2-1; ktmp++) {
	      hmm->tp[ktmp][f4H_EPSILON] += wt; hmm->tp[ktmp][f4H_DELTA] += wt; hmm->tp[ktmp][f4H_EPSILONP] += wt; //DD
	      hmm->tp[ktmp][f4H_EPSILONP] += wt; hmm->tp[ktmp][f4H_GAMMA] += wt; //DM
      }
	  }
	else  {
	  switch (st2) {
	  case f4T_M: hmm->tp[0][f4H_GAMMA] += wt; break; //MM
	  case f4T_I: hmm->tp[0][f4H_ALPHA] += wt; break; //MI
	  case f4T_D: hmm->tp[0][f4H_DELTA] += wt; break; //MD
	  default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	  }
	}
      }
      else if (st == f4T_M) {
     	switch (st2) {
	case f4T_M: hmm->tp[0][f4H_GAMMA] += wt; break; //MM
	case f4T_I: hmm->tp[k][f4H_ALPHA] += wt; break; //MI
	case f4T_D: hmm->tp[k][f4H_DELTA] += wt; break; //MD
	case f4T_E: hmm->tp[0][f4H_GAMMA] += wt; break; /* k==M. A local alignment would've been Mk->X->E. */ //MM
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == f4T_I) {
	switch (st2) {
	case f4T_M: hmm->tp[k][f4H_BETAP] += wt; hmm->tp[k][f4H_GAMMA] += wt; break; //IM
	case f4T_I: hmm->tp[k][f4H_BETA] += wt; hmm->tp[k][f4H_ALPHA] += wt; hmm->tp[k][f4H_BETAP] += wt; break; //II
	case f4T_E: hmm->tp[k][f4H_GAMMA] += wt; hmm->tp[k][f4H_BETAP] += wt; break; /* k==M. */ //IM
  case f4T_D: hmm->tp[k][f4H_DELTA] += wt; hmm->tp[k][f4H_BETAP] += wt; break; /* Add this for f4-HMM */ //ID
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
      else if (st == f4T_D) {
	switch (st2) {
	case f4T_M: hmm->tp[k][f4H_EPSILONP] += wt; hmm->tp[k][f4H_GAMMA] += wt; break; //DM
	case f4T_D: hmm->tp[k][f4H_EPSILON] += wt; hmm->tp[k][f4H_DELTA] += wt; hmm->tp[k][f4H_EPSILONP] += wt; break; //DD
	case f4T_E: hmm->tp[k][f4H_GAMMA] += wt; hmm->tp[k][f4H_EPSILONP] += wt; break; /* k==M. A local alignment would've been Dk->X->E. */ //DM
  case f4T_I: hmm->tp[k][f4H_ALPHA] += wt; hmm->tp[k][f4H_EPSILONP] += wt; break; /* Add this for f4-HMM */ //DI
	default:     ESL_EXCEPTION(eslEINVAL, "bad transition in trace");
	}
      }
    } /* end loop over trace position */
  return eslOK;
}

/*--------------------- end, trace counting ---------------------*/