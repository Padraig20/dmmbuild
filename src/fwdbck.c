#include "dummer.h"

static int forward_engine (int do_full, const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, F4_OMX *fwd, float *opt_sc);

/* Function:  f4_ForwardParser()
 * Synopsis:  The Forward algorithm, linear memory parsing version.
 * Incept:    SRE, Fri Aug 15 19:05:26 2008 [Casa de Gatos]
 *
 * Purpose:   A linear $O(M+L)$ memory algorithm is
 *            used, keeping only the DP matrix values for the special
 *            (BENCJ) states. These are sufficient to do posterior
 *            decoding to identify high-probability regions where
 *            domains are.
 * 
 *            The caller must provide a suitably allocated "parsing"
 *            <ox> by calling <ox = f4_omx_Create(M, 0, L)> or
 *            <f4_omx_GrowTo(ox, M, 0, L)>.
 *            
 * Args:      dsq     - digital target sequence, 1..L
 *            L       - length of dsq in residues          
 *            om      - optimized profile
 *            ox      - RETURN: Forward DP matrix
 *            ret_sc  - RETURN: Forward score (in nats)          
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <ox> allocation is too small, or if the profile
 *            isn't in local alignment mode.
 *            <eslERANGE> if the score exceeds the limited range of
 *            a probability-space odds ratio.
 *            In either case, <*opt_sc> is undefined.
 */
int
f4_ForwardParser(const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, F4_OMX *ox, float *opt_sc)
{
  return forward_engine(FALSE, dsq, L, om, ox, opt_sc);
}

static int
forward_engine(int do_full, const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, F4_OMX *ox, float *opt_sc)
{
  register __m128 mpv, dpv, ipv;   /* previous row values                                       */
  register __m128 sv;		   /* temp storage of 1 curr row value in progress              */
  register __m128 dcv;		   /* delayed storage of D(i,q+1)                               */
  register __m128 xEv;		   /* E state: keeps max for Mk->E as we go                     */
  register __m128 xBv;		   /* B state: splatted vector of B[i-1] for B->Mk calculations */
  __m128   zerov;		   /* splatted 0.0's in a vector                                */
  float    xN, xE, xB, xC, xJ;	   /* special states' scores                                    */
  int i;			   /* counter over sequence positions 1..L                      */
  int q;			   /* counter over quads 0..nq-1                                */
  int j;			   /* counter over DD iterations (4 is full serialization)      */
  int Q       = f4O_NQF(om->M);	   /* segment length: # of vectors                              */
  __m128 *dpc = ox->dpf[0];        /* current row, for use in {MDI}MO(dpp,q) access macro       */
  __m128 *dpp;                     /* previous row, for use in {MDI}MO(dpp,q) access macro      */
  __m128 *rp;			   /* will point at om->rfv[x] for residue x[i]                 */
  __m128 *tp;			   /* will point into (and step thru) om->tfv                   */

  /* Initialization. */
  ox->M  = om->M;
  ox->L  = L;
  ox->has_own_scales = TRUE; 	/* all forward matrices control their own scalefactors */
  zerov  = _mm_setzero_ps();
  for (q = 0; q < Q; q++)
    MMO(dpc,q) = IMO(dpc,q) = DMO(dpc,q) = zerov;
  xE    = ox->xmx[f4X_E] = 0.;
  xN    = ox->xmx[f4X_N] = 1.;
  xJ    = ox->xmx[f4X_J] = 0.;
  xB    = ox->xmx[f4X_B] = om->xf[f4O_N][f4O_MOVE];
  xC    = ox->xmx[f4X_C] = 0.;

  ox->xmx[f4X_SCALE] = 1.0;
  ox->totscale       = 0.0;

  for (i = 1; i <= L; i++)
    {
      dpp   = dpc;                      
      dpc   = ox->dpf[do_full * i];     /* avoid conditional, use do_full as kronecker delta */
      rp    = om->rfv[dsq[i]];
      tp    = om->tfv;
      dcv   = _mm_setzero_ps();
      xEv   = _mm_setzero_ps();
      xBv   = _mm_set1_ps(xB);

      /* Right shifts by 4 bytes. 4,8,12,x becomes x,4,8,12.  Shift zeros on. */
      mpv   = esl_sse_rightshiftz_float(MMO(dpp,Q-1));
      dpv   = esl_sse_rightshiftz_float(DMO(dpp,Q-1));
      ipv   = esl_sse_rightshiftz_float(IMO(dpp,Q-1));
      
      for (q = 0; q < Q; q++)
	{
	  /* Calculate new MMO(i,q); don't store it yet, hold it in sv. */
	  sv   =                _mm_mul_ps(xBv, *tp);  tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(mpv, *tp)); tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
	  sv   = _mm_add_ps(sv, _mm_mul_ps(dpv, *tp)); tp++;
	  sv   = _mm_mul_ps(sv, *rp);                  rp++;
	  xEv  = _mm_add_ps(xEv, sv);
	  
	  /* Load {MDI}(i-1,q) into mpv, dpv, ipv;
	   * {MDI}MX(q) is then the current, not the prev row
	   */
	  mpv = MMO(dpp,q);
	  dpv = DMO(dpp,q);
	  ipv = IMO(dpp,q);

	  /* Do the delayed stores of {MD}(i,q) now that memory is usable */
	  MMO(dpc,q) = sv;
	  DMO(dpc,q) = dcv;

	  /* Calculate the next D(i,q+1) partially: M->D only;
           * delay storage, holding it in dcv
	   */
	  dcv   = _mm_mul_ps(sv, *tp); tp++;

	  /* Calculate and store I(i,q); assumes odds ratio for emission is 1.0 */
	  sv         =                _mm_mul_ps(mpv, *tp);  tp++;
	  IMO(dpc,q) = _mm_add_ps(sv, _mm_mul_ps(ipv, *tp)); tp++;
	}	  

      /* Now the DD paths. We would rather not serialize them but 
       * in an accurate Forward calculation, we have few options.
       */
      /* dcv has carried through from end of q loop above; store it 
       * in first pass, we add M->D and D->D path into DMX
       */
      /* We're almost certainly're obligated to do at least one complete 
       * DD path to be sure: 
       */
      dcv        = esl_sse_rightshiftz_float(dcv);
      DMO(dpc,0) = zerov;
      tp         = om->tfv + 7*Q;	/* set tp to start of the DD's */
      for (q = 0; q < Q; q++) 
	{
	  DMO(dpc,q) = _mm_add_ps(dcv, DMO(dpc,q));	
	  dcv        = _mm_mul_ps(DMO(dpc,q), *tp); tp++; /* extend DMO(q), so we include M->D and D->D paths */
	}

      /* now. on small models, it seems best (empirically) to just go
       * ahead and serialize. on large models, we can do a bit better,
       * by testing for when dcv (DD path) accrued to DMO(q) is below
       * machine epsilon for all q, in which case we know DMO(q) are all
       * at their final values. The tradeoff point is (empirically) somewhere around M=100,
       * at least on my desktop. We don't worry about the conditional here;
       * it's outside any inner loops.
       */
      if (om->M < 100)
	{			/* Fully serialized version */
	  for (j = 1; j < 4; j++)
	    {
	      dcv = esl_sse_rightshiftz_float(dcv);
	      tp  = om->tfv + 7*Q;	/* set tp to start of the DD's */
	      for (q = 0; q < Q; q++) 
		{ /* note, extend dcv, not DMO(q); only adding DD paths now */
		  DMO(dpc,q) = _mm_add_ps(dcv, DMO(dpc,q));	
		  dcv        = _mm_mul_ps(dcv, *tp);   tp++; 
		}	    
	    }
	} 
      else
	{			/* Slightly parallelized version, but which incurs some overhead */
	  for (j = 1; j < 4; j++)
	    {
	      register __m128 cv;	/* keeps track of whether any DD's change DMO(q) */

	      dcv = esl_sse_rightshiftz_float(dcv);
	      tp  = om->tfv + 7*Q;	/* set tp to start of the DD's */
	      cv  = zerov;
	      for (q = 0; q < Q; q++) 
		{ /* using cmpgt below tests if DD changed any DMO(q) *without* conditional branch */
		  sv         = _mm_add_ps(dcv, DMO(dpc,q));	
		  cv         = _mm_or_ps(cv, _mm_cmpgt_ps(sv, DMO(dpc,q))); 
		  DMO(dpc,q) = sv;	                                    /* store new DMO(q) */
		  dcv        = _mm_mul_ps(dcv, *tp);   tp++;            /* note, extend dcv, not DMO(q) */
		}	    
	      if (! _mm_movemask_ps(cv)) break; /* DD's didn't change any DMO(q)? Then done, break out. */
	    }
	}

      /* Add D's to xEv */
      for (q = 0; q < Q; q++) xEv = _mm_add_ps(DMO(dpc,q), xEv);

      /* Finally the "special" states, which start from Mk->E (->C, ->J->B) */
      /* The following incantation is a horizontal sum of xEv's elements  */
      /* These must follow DD calculations, because D's contribute to E in Forward
       * (as opposed to Viterbi)
       */
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(0, 3, 2, 1)));
      xEv = _mm_add_ps(xEv, _mm_shuffle_ps(xEv, xEv, _MM_SHUFFLE(1, 0, 3, 2)));
      _mm_store_ss(&xE, xEv);

      xN =  xN * om->xf[f4O_N][f4O_LOOP];
      xC = (xC * om->xf[f4O_C][f4O_LOOP]) +  (xE * om->xf[f4O_E][f4O_MOVE]);
      xJ = (xJ * om->xf[f4O_J][f4O_LOOP]) +  (xE * om->xf[f4O_E][f4O_LOOP]);
      xB = (xJ * om->xf[f4O_J][f4O_MOVE]) +  (xN * om->xf[f4O_N][f4O_MOVE]);
      /* and now xB will carry over into next i, and xC carries over after i=L */

      /* Sparse rescaling. xE above threshold? trigger a rescaling event.            */
      if (xE > 1.0e4)	/* that's a little less than e^10, ~10% of our dynamic range */
	{
	  xN  = xN / xE;
	  xC  = xC / xE;
	  xJ  = xJ / xE;
	  xB  = xB / xE;
	  xEv = _mm_set1_ps(1.0 / xE);
	  for (q = 0; q < Q; q++)
	    {
	      MMO(dpc,q) = _mm_mul_ps(MMO(dpc,q), xEv);
	      DMO(dpc,q) = _mm_mul_ps(DMO(dpc,q), xEv);
	      IMO(dpc,q) = _mm_mul_ps(IMO(dpc,q), xEv);
	    }
	  ox->xmx[i*f4X_NXCELLS+f4X_SCALE] = xE;
	  ox->totscale += log(xE);
	  xE = 1.0;		
	}
      else ox->xmx[i*f4X_NXCELLS+f4X_SCALE] = 1.0;

      /* Storage of the specials.  We could've stored these already
       * but using xE, etc. variables makes it easy to convert this
       * code to O(M) memory versions just by deleting storage steps.
       */
      ox->xmx[i*f4X_NXCELLS+f4X_E] = xE;
      ox->xmx[i*f4X_NXCELLS+f4X_N] = xN;
      ox->xmx[i*f4X_NXCELLS+f4X_J] = xJ;
      ox->xmx[i*f4X_NXCELLS+f4X_B] = xB;
      ox->xmx[i*f4X_NXCELLS+f4X_C] = xC;

    } /* end loop over sequence residues 1..L */

  /* finally C->T, and flip total score back to log space (nats) */
  /* On overflow, xC is inf or nan (nan arises because inf*0 = nan). */
  /* On an underflow (which shouldn't happen), we counterintuitively return infinity:
   * the effect of this is to force the caller to rescore us with full range.
   */
  if       (isnan(xC))        ESL_EXCEPTION(eslERANGE, "forward score is NaN");
  else if  (L>0 && xC == 0.0) ESL_EXCEPTION(eslERANGE, "forward score underflow (is 0.0)");     /* if L==0, xC *should* be 0.0; J5/118 */
  else if  (isinf(xC) == 1)   ESL_EXCEPTION(eslERANGE, "forward score overflow (is infinity)");

  if (opt_sc != NULL) *opt_sc = ox->totscale + log(xC * om->xf[f4O_C][f4O_MOVE]);
  return eslOK;
}