#include "dummer.h"

/* Function:  f4_oprofile_Create()
 * Synopsis:  Allocate an optimized profile structure.
 * Incept:    SRE, Sun Nov 25 12:03:19 2007 [Casa de Gatos]
 *
 * Purpose:   Allocate for profiles of up to <allocM> nodes for digital alphabet <abc>.
 *
 * Throws:    <NULL> on allocation error.
 */
F4_OPROFILE *
f4_oprofile_Create(int allocM, const ESL_ALPHABET *abc)
{
  int          status;
  F4_OPROFILE *om  = NULL;
  int          nqb = f4O_NQB(allocM); /* # of uchar vectors needed for query */
  int          nqw = f4O_NQW(allocM); /* # of sword vectors needed for query */
  int          nqf = f4O_NQF(allocM); /* # of float vectors needed for query */
  int          nqs = nqb + f4O_EXTRA_SB;
  int          x;

  /* level 0 */
  ESL_ALLOC(om, sizeof(F4_OPROFILE));
  om->rbv_mem = NULL;
  om->sbv_mem = NULL;
  om->rwv_mem = NULL;
  om->twv_mem = NULL;
  om->rfv_mem = NULL;
  om->tfv_mem = NULL;
  om->rbv     = NULL;
  om->sbv     = NULL;
  om->rwv     = NULL;
  om->twv     = NULL;
  om->rfv     = NULL;
  om->tfv     = NULL;
  om->clone   = 0;

  /* level 1 */
  ESL_ALLOC(om->rbv_mem, sizeof(__m128i) * nqb  * abc->Kp          +15); /* +15 is for manual 16-byte alignment */
  ESL_ALLOC(om->sbv_mem, sizeof(__m128i) * nqs  * abc->Kp          +15); 
  ESL_ALLOC(om->rwv_mem, sizeof(__m128i) * nqw  * abc->Kp          +15);                     
  ESL_ALLOC(om->twv_mem, sizeof(__m128i) * nqw  * f4O_NTRANS       +15);   
  ESL_ALLOC(om->rfv_mem, sizeof(__m128)  * nqf  * abc->Kp          +15);                     
  ESL_ALLOC(om->tfv_mem, sizeof(__m128)  * nqf  * f4O_NTRANS       +15);    

  ESL_ALLOC(om->rbv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->sbv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->rwv, sizeof(__m128i *) * abc->Kp); 
  ESL_ALLOC(om->rfv, sizeof(__m128  *) * abc->Kp); 

  /* align vector memory on 16-byte boundaries */
  om->rbv[0] = (__m128i *) (((unsigned long int) om->rbv_mem + 15) & (~0xf));
  om->sbv[0] = (__m128i *) (((unsigned long int) om->sbv_mem + 15) & (~0xf));
  om->rwv[0] = (__m128i *) (((unsigned long int) om->rwv_mem + 15) & (~0xf));
  om->twv    = (__m128i *) (((unsigned long int) om->twv_mem + 15) & (~0xf));
  om->rfv[0] = (__m128  *) (((unsigned long int) om->rfv_mem + 15) & (~0xf));
  om->tfv    = (__m128  *) (((unsigned long int) om->tfv_mem + 15) & (~0xf));

  /* set the rest of the row pointers for match emissions */
  for (x = 1; x < abc->Kp; x++) {
    om->rbv[x] = om->rbv[0] + (x * nqb);
    om->sbv[x] = om->sbv[0] + (x * nqs);
    om->rwv[x] = om->rwv[0] + (x * nqw);
    om->rfv[x] = om->rfv[0] + (x * nqf);
  }
  om->allocQ16  = nqb;
  om->allocQ8   = nqw;
  om->allocQ4   = nqf;

  /* Remaining initializations */
  om->tbm_b     = 0;
  om->tec_b     = 0;
  om->tjb_b     = 0;
  om->scale_b   = 0.0f;
  om->base_b    = 0;
  om->bias_b    = 0;

  om->scale_w      = 0.0f;
  om->base_w       = 0;
  om->ddbound_w    = 0;
  om->ncj_roundoff = 0.0f;	

  for (x = 0; x < f4_NOFFSETS; x++) om->offs[x]    = -1;
  for (x = 0; x < f4_NEVPARAM; x++) om->evparam[x] = f4_EVPARAM_UNSET;
  for (x = 0; x < f4_NCUTOFFS; x++) om->cutoff[x]  = f4_CUTOFF_UNSET;
  for (x = 0; x < f4_MAXABET;  x++) om->compo[x]   = f4_COMPO_UNSET;

  om->name      = NULL;
  om->acc       = NULL;
  om->desc      = NULL;

  /* in a F4_OPROFILE, we always allocate for the optional RF, CS annotation.  
   * we only rely on the leading \0 to signal that it's unused, but 
   * we initialize all this memory to zeros to shut valgrind up about 
   * fwrite'ing uninitialized memory in the io functions.
   */
  ESL_ALLOC(om->rf,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->mm,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->cs,          sizeof(char) * (allocM+2));
  ESL_ALLOC(om->consensus,   sizeof(char) * (allocM+2));
  memset(om->rf,       '\0', sizeof(char) * (allocM+2));
  memset(om->mm,       '\0', sizeof(char) * (allocM+2));
  memset(om->cs,       '\0', sizeof(char) * (allocM+2));
  memset(om->consensus,'\0', sizeof(char) * (allocM+2));

  om->abc        = abc;
  om->L          = 0;
  om->M          = 0;
  om->max_length = -1;
  om->allocM     = allocM;
  om->mode       = f4_NO_MODE;
  om->nj         = 0.0f;
  return om;

 ERROR:
  f4_oprofile_Destroy(om);
  return NULL;
}

/* Function:  f4_oprofile_Destroy()
 * Synopsis:  Frees an optimized profile structure.
 */
void
f4_oprofile_Destroy(F4_OPROFILE *om)
{
  if (om == NULL) return;

  if (om->clone == 0)
    {
      if (om->rbv_mem   != NULL) free(om->rbv_mem);
      if (om->sbv_mem   != NULL) free(om->sbv_mem);
      if (om->rwv_mem   != NULL) free(om->rwv_mem);
      if (om->twv_mem   != NULL) free(om->twv_mem);
      if (om->rfv_mem   != NULL) free(om->rfv_mem);
      if (om->tfv_mem   != NULL) free(om->tfv_mem);
      if (om->rbv       != NULL) free(om->rbv);
      if (om->sbv       != NULL) free(om->sbv);
      if (om->rwv       != NULL) free(om->rwv);
      if (om->rfv       != NULL) free(om->rfv);
      if (om->name      != NULL) free(om->name);
      if (om->acc       != NULL) free(om->acc);
      if (om->desc      != NULL) free(om->desc);
      if (om->rf        != NULL) free(om->rf);
      if (om->mm        != NULL) free(om->mm);
      if (om->cs        != NULL) free(om->cs);
      if (om->consensus != NULL) free(om->consensus);
    }

  free(om);
}

/* unbiased_byteify()
 * Convert original transition score to a rounded uchar cost
 * Transition scores for MSVFilter get this treatment.
 * e.g. a score of -2.1, with scale 3.0, becomes a cost of 6.
 * (A cost of +255 is our -infinity "prohibited event")
 */
static uint8_t 
unbiased_byteify(F4_OPROFILE *om, float sc)
{
  uint8_t b;

  sc  = -1.0f * roundf(om->scale_b * sc);       /* ugh. sc is now an integer cost represented in a float...    */
  b   = (sc > 255.) ? 255 : (uint8_t) sc;	/* and now we cast and saturate it to an unsigned char cost... */
  return b;
}

/* biased_byteify()
 * Converts original log-odds residue score to a rounded biased uchar cost.
 * Match emission scores for MSVFilter get this treatment.
 * e.g. a score of +3.2, with scale 3.0 and bias 12, becomes 2.
 *    3.2*3 = 9.6; rounded = 10; bias-10 = 2.
 * When used, we add the bias, then subtract this cost.
 * (A cost of +255 is our -infinity "prohibited event")
 */
static uint8_t
biased_byteify(F4_OPROFILE *om, float sc)
{
  uint8_t b;

  sc  = -1.0f * roundf(om->scale_b * sc);                          /* ugh. sc is now an integer cost represented in a float...           */
  b   = (sc > 255 - om->bias_b) ? 255 : (uint8_t) sc + om->bias_b; /* and now we cast, saturate, and bias it to an unsigned char cost... */
  return b;
}

/* wordify()
 * Converts log probability score to a rounded signed 16-bit integer cost.
 * Both emissions and transitions for ViterbiFilter get this treatment.
 * No bias term needed, because we use signed words. 
 *   e.g. a score of +3.2, with scale 500.0, becomes +1600.
 */
static int16_t 
wordify(F4_OPROFILE *om, float sc)
{
  sc  = roundf(om->scale_w * sc);
  if      (sc >=  32767.0) return  32767;
  else if (sc <= -32768.0) return -32768;
  else return (int16_t) sc;
}

/* sf_conversion():
 * Author: Bjarne Knudsen
 * 
 * Generates the SSVFilter() parts of the profile <om> scores
 * from the completed MSV score.  This includes calculating 
 * special versions of the match scores for using the the
 * ssv filter.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    (no abnormal error conditions)
 */
static int
sf_conversion(F4_OPROFILE *om)
{
  int     M   = om->M;		/* length of the query                                          */
  int     nq  = f4O_NQB(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  __m128i tmp;
  __m128i tmp2;

  /* We now want to fill out om->sbv with om->rbv - bias for use in the
   * SSV filter. The only challenge is that the om->rbv values are
   * unsigned and generally use the whole scale while the om->sbv
   * values are signed. To solve that problem we perform the following
   * calculation:
   *
   *   ((127 + bias) - rbv) ^ 127
   *
   * where the subtraction is unsigned saturated and the addition is
   * unsigned (it will not overflow, since bias is a small positive
   * number). The f(x) = x ^ 127 combined with a change from unsigned
   * to signed numbers have the same effect as f(x) = -x + 127. So if
   * we regard the above as signed instead of unsigned it is equal to:
   *
   *   -((127 + bias) - rbv) + 127 = rbv - bias
   *
   * which is what we want. The reason for this slightly complex idea
   * is that we wish the transformation to be fast, especially for
   * hmmscan where many models are loaded.
   */

  tmp = _mm_set1_epi8((int8_t) (om->bias_b + 127));
  tmp2  = _mm_set1_epi8(127);

  for (x = 0; x < om->abc->Kp; x++)
    {
      for (q = 0;  q < nq;            q++) om->sbv[x][q] = _mm_xor_si128(_mm_subs_epu8(tmp, om->rbv[x][q]), tmp2);
      for (q = nq; q < nq + f4O_EXTRA_SB; q++) om->sbv[x][q] = om->sbv[x][q % nq];
    }

  return eslOK;
}

/* mf_conversion(): 
 * 
 * This builds the MSVFilter() parts of the profile <om>, scores
 * in lspace uchars (16-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 */
static int
mf_conversion(const F4_PROFILE *gm, F4_OPROFILE *om)
{
  int     M   = gm->M;		/* length of the query                                          */
  int     nq  = f4O_NQB(M);     /* segment length; total # of striped vectors needed            */
  float   max = 0.0;		/* maximum residue score: used for unsigned emission score bias */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     z;			/* counter within elements of one SIMD minivector               */
  union { __m128i v; uint8_t i[16]; } tmp; /* used to align and load simd minivectors           */

  if (nq > om->allocQ16) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* First we determine the basis for the limited-precision MSVFilter scoring system. 
   * Default: 1/3 bit units, base offset 190:  range 0..255 => -190..65 => -63.3..21.7 bits
   * See J2/66, J4/138 for analysis.
   */
  for (x = 0; x < gm->abc->K; x++)  max = ESL_MAX(max, esl_vec_FMax(gm->rsc[x], (M+1)*2));
  om->scale_b = 3.0 / eslCONST_LOG2;                    /* scores in units of third-bits */
  om->base_b  = 190;
  om->bias_b  = unbiased_byteify(om, -1.0 * max);

  /* striped match costs: start at k=1.  */
  for (x = 0; x < gm->abc->Kp; x++)
  {
    for (q = 0, k = 1; q < nq; q++, k++)
    {
      for (z = 0; z < 16; z++) tmp.i[z] = ((k+ z*nq <= M) ? biased_byteify(om, f4P_MSC(gm, k+z*nq, x)) : 255);
      om->rbv[x][q]   = tmp.v;
    }
  }

  /* transition costs */
  om->tbm_b = unbiased_byteify(om, logf(2.0f / ((float) gm->M * (float) (gm->M+1)))); /* constant B->Mk penalty        */
  om->tec_b = unbiased_byteify(om, logf(0.5f));                                       /* constant multihit E->C = E->J */
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float) (gm->L+3))); /* this adopts the L setting of the parent profile */

  sf_conversion(om);

  return eslOK;
}


/* vf_conversion(): 
 * 
 * This builds the ViterbiFilter() parts of the profile <om>, scores
 * in lspace swords (8-way parallel), by rescaling, rounding, and
 * casting the scores in <gm>.
 * 
 * Returns <eslOK> on success;
 * throws <eslEINVAL> if <om> hasn't been allocated properly.
 */
static int
vf_conversion(const F4_PROFILE *gm, F4_OPROFILE *om)
{
  int     M   = gm->M;		/* length of the query                                          */
  int     nq  = f4O_NQW(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = f4O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     tg;			/* transition index in gm                                       */
  int     j;			/* counter in interleaved vector arrays in the profile          */
  int     ddtmp;		/* used in finding worst DD transition bound                    */
  int16_t  maxval;		/* used to prevent zero cost II                                 */
  int16_t  val;
  union { __m128i v; int16_t i[8]; } tmp; /* used to align and load simd minivectors            */

  if (nq > om->allocQ8) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* First set the basis for the limited-precision scoring system. 
   * Default: 1/500 bit units, base offset 12000:  range -32768..32767 => -44768..20767 => -89.54..41.53 bits
   * See J4/138 for analysis.
   */
  om->scale_w = 500.0 / eslCONST_LOG2;
  om->base_w  = 12000;

  /* striped match scores */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++)
      {
	for (z = 0; z < 8; z++) tmp.i[z] = ((k+ z*nq <= M) ? wordify(om, f4P_MSC(gm, k+z*nq, x)) : -32768);
	om->rwv[x][q]   = tmp.v;
      }

  /* Transition costs, all but the DD's. */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = f4O_BM; t <= f4O_II; t++) /* this loop of 7 transitions depends on the order in f4o_tsc_e */
	{
	  switch (t) {
	  case f4O_BM: tg = f4P_BM;  kb = k-1; maxval =  0; break; /* gm has tBMk stored off by one! start from k=0 not 1   */
	  case f4O_MM: tg = f4P_MM;  kb = k-1; maxval =  0; break; /* MM, DM, IM vectors are rotated by -1, start from k=0  */
	  case f4O_IM: tg = f4P_IM;  kb = k-1; maxval =  0; break;
	  case f4O_DM: tg = f4P_DM;  kb = k-1; maxval =  0; break;
	  case f4O_MD: tg = f4P_MD;  kb = k;   maxval =  0; break; /* the remaining ones are straight up  */
	  case f4O_MI: tg = f4P_MI;  kb = k;   maxval =  0; break; 
	  case f4O_II: tg = f4P_II;  kb = k;   maxval = -1; break; 
	  }

	  for (z = 0; z < 8; z++) {
	    val      = ((kb+ z*nq < M) ? wordify(om, f4P_TSC(gm, kb+ z*nq, tg)) : -32768);
	    tmp.i[z] = (val <= maxval) ? val : maxval; /* do not allow an II transition cost of 0, or hell may occur. */
	  }
	  om->twv[j++] = tmp.v;
	}
    }

  /* Finally the DD's, which are at the end of the optimized tsc vector; (j is already sitting there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 8; z++) tmp.i[z] = ((k+ z*nq < M) ? wordify(om, f4P_TSC(gm, k+ z*nq, f4P_DD)) : -32768);
      om->twv[j++] = tmp.v;
    }

  /* Specials. (Actually in same order in om and gm, but we copy in general form anyway.)  */
  /* VF CC,NN,JJ transitions hardcoded zero; -3.0 nat approximation used instead; this papers
   * over a length independence problem, where the approximation weirdly outperforms the
   * exact solution, probably indicating that the model's Pascal distribution is problematic,
   * and the "approximation" is in fact closer to the One True Model, the mythic H4 supermodel.
   * [xref J5/36] 
   */
  om->xw[f4O_E][f4O_LOOP] = wordify(om, gm->xsc[f4P_E][f4P_LOOP]);  
  om->xw[f4O_E][f4O_MOVE] = wordify(om, gm->xsc[f4P_E][f4P_MOVE]);
  om->xw[f4O_N][f4O_MOVE] = wordify(om, gm->xsc[f4P_N][f4P_MOVE]);
  om->xw[f4O_N][f4O_LOOP] = 0;                                        /* was wordify(om, gm->xsc[f4P_N][f4P_LOOP]); */
  om->xw[f4O_C][f4O_MOVE] = wordify(om, gm->xsc[f4P_C][f4P_MOVE]);
  om->xw[f4O_C][f4O_LOOP] = 0;                                        /* was wordify(om, gm->xsc[f4P_C][f4P_LOOP]); */
  om->xw[f4O_J][f4O_MOVE] = wordify(om, gm->xsc[f4P_J][f4P_MOVE]);
  om->xw[f4O_J][f4O_LOOP] = 0;                                        /* was wordify(om, gm->xsc[f4P_J][f4P_LOOP]); */

  om->ncj_roundoff = 0.0; /* goes along with NN=CC=JJ=0, -3.0 nat approximation */
                          /* otherwise, would be = om->scale_w * gm->xsc[f4P_N][f4P_LOOP] -  om->xw[f4O_N][f4O_LOOP];   */
			  /* see J4/150 for discussion of VF error suppression, superceded by the -3.0 nat approximation */

  /* Transition score bound for "lazy F" DD path evaluation (xref J2/52) */
  om->ddbound_w = -32768;	
  for (k = 2; k < M-1; k++) 
    {
      ddtmp         = (int) wordify(om, f4P_TSC(gm, k,   f4P_DD));
      ddtmp        += (int) wordify(om, f4P_TSC(gm, k+1, f4P_DM));
      ddtmp        -= (int) wordify(om, f4P_TSC(gm, k+1, f4P_BM));
      om->ddbound_w = ESL_MAX(om->ddbound_w, ddtmp);
    }

  return eslOK;
}


/* fb_conversion()
 * This builds the Forward/Backward part of the optimized profile <om>,
 * where we use odds ratios (not log-odds scores).
 */
static int
fb_conversion(const F4_PROFILE *gm, F4_OPROFILE *om)
{
  int     M   = gm->M;		/* length of the query                                          */
  int     nq  = f4O_NQF(M);     /* segment length; total # of striped vectors needed            */
  int     x;			/* counter over residues                                        */
  int     q;			/* q counts over total # of striped vectors, 0..nq-1            */
  int     k;			/* the usual counter over model nodes 1..M                      */
  int     kb;			/* possibly offset base k for loading om's TSC vectors          */
  int     z;			/* counter within elements of one SIMD minivector               */
  int     t;			/* counter over transitions 0..7 = f4O_{BM,MM,IM,DM,MD,MI,II,DD}*/
  int     tg;			/* transition index in gm                                       */
  int     j;			/* counter in interleaved vector arrays in the profile          */
  union { __m128 v; float x[4]; } tmp; /* used to align and load simd minivectors               */

  if (nq > om->allocQ4) ESL_EXCEPTION(eslEINVAL, "optimized profile is too small to hold conversion");

  /* striped match scores: start at k=1 */
  for (x = 0; x < gm->abc->Kp; x++)
    for (k = 1, q = 0; q < nq; q++, k++)
      {
	for (z = 0; z < 4; z++) tmp.x[z] = (k+ z*nq <= M) ? f4P_MSC(gm, k+z*nq, x) : -eslINFINITY;
	om->rfv[x][q] = esl_sse_expf(tmp.v);
      }


  /* Transition scores, all but the DD's. */
  for (j = 0, k = 1, q = 0; q < nq; q++, k++)
    {
      for (t = f4O_BM; t <= f4O_II; t++) /* this loop of 7 transitions depends on the order in the definition of f4o_tsc_e */
	{
	  switch (t) {
	  case f4O_BM: tg = f4P_BM;  kb = k-1; break; /* gm has tBMk stored off by one! start from k=0 not 1 */
	  case f4O_MM: tg = f4P_MM;  kb = k-1; break; /* MM, DM, IM quads are rotated by -1, start from k=0  */
	  case f4O_IM: tg = f4P_IM;  kb = k-1; break;
	  case f4O_DM: tg = f4P_DM;  kb = k-1; break;
	  case f4O_MD: tg = f4P_MD;  kb = k;   break; /* the remaining ones are straight up  */
	  case f4O_MI: tg = f4P_MI;  kb = k;   break; 
	  case f4O_II: tg = f4P_II;  kb = k;   break; 
	  }

	  for (z = 0; z < 4; z++) tmp.x[z] = (kb+z*nq < M) ? f4P_TSC(gm, kb+z*nq, tg) : -eslINFINITY;
	  om->tfv[j++] = esl_sse_expf(tmp.v);
	}
    }

  /* And finally the DD's, which are at the end of the optimized tfv vector; (j is already there) */
  for (k = 1, q = 0; q < nq; q++, k++)
    {
      for (z = 0; z < 4; z++) tmp.x[z] = (k+z*nq < M) ? f4P_TSC(gm, k+z*nq, f4P_DD) : -eslINFINITY;
      om->tfv[j++] = esl_sse_expf(tmp.v);
    }

  /* Specials. (These are actually in exactly the same order in om and
   *  gm, but we copy in general form anyway.)
   */
  om->xf[f4O_E][f4O_LOOP] = expf(gm->xsc[f4P_E][f4P_LOOP]);  
  om->xf[f4O_E][f4O_MOVE] = expf(gm->xsc[f4P_E][f4P_MOVE]);
  om->xf[f4O_N][f4O_LOOP] = expf(gm->xsc[f4P_N][f4P_LOOP]);
  om->xf[f4O_N][f4O_MOVE] = expf(gm->xsc[f4P_N][f4P_MOVE]);
  om->xf[f4O_C][f4O_LOOP] = expf(gm->xsc[f4P_C][f4P_LOOP]);
  om->xf[f4O_C][f4O_MOVE] = expf(gm->xsc[f4P_C][f4P_MOVE]);
  om->xf[f4O_J][f4O_LOOP] = expf(gm->xsc[f4P_J][f4P_LOOP]);
  om->xf[f4O_J][f4O_MOVE] = expf(gm->xsc[f4P_J][f4P_MOVE]);

  return eslOK;
}

/* Function:  f4_oprofile_Convert()
 * Synopsis:  Converts standard profile to an optimized one.
 * Incept:    SRE, Mon Nov 26 07:38:57 2007 [Janelia]
 *
 * Purpose:   Convert a standard profile <gm> to an optimized profile <om>,
 *            where <om> has already been allocated for a profile of at 
 *            least <gm->M> nodes and the same emission alphabet <gm->abc>.
 *
 * Args:      gm - profile to optimize
 *            om - allocated optimized profile for holding the result.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if <gm>, <om> aren't compatible. 
 *            <eslEMEM> on allocation failure.
 */
int
f4_oprofile_Convert(const F4_PROFILE *gm, F4_OPROFILE *om)
{
  int status, z;

  /* Set these first so they are available in the following calls */
  om->mode       = gm->mode;
  om->L          = gm->L;
  om->M          = gm->M;
  om->nj         = gm->nj;
  om->max_length = gm->max_length;

  if (gm->abc->type != om->abc->type)  ESL_EXCEPTION(eslEINVAL, "alphabets of the two profiles don't match");  
  if (gm->M         >  om->allocM)     ESL_EXCEPTION(eslEINVAL, "oprofile is too small");  

  if ((status =  mf_conversion(gm, om)) != eslOK) return status;   /* MSVFilter()'s information     */
  if ((status =  vf_conversion(gm, om)) != eslOK) return status;   /* ViterbiFilter()'s information */
  if ((status =  fb_conversion(gm, om)) != eslOK) return status;   /* ForwardFilter()'s information */

  if (om->name != NULL) free(om->name);
  if (om->acc  != NULL) free(om->acc);
  if (om->desc != NULL) free(om->desc);
  if ((status = esl_strdup(gm->name, -1, &(om->name))) != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->acc,  -1, &(om->acc)))  != eslOK) goto ERROR;
  if ((status = esl_strdup(gm->desc, -1, &(om->desc))) != eslOK) goto ERROR;
  strcpy(om->rf,        gm->rf);
  strcpy(om->mm,        gm->mm);
  strcpy(om->cs,        gm->cs);
  strcpy(om->consensus, gm->consensus);
  for (z = 0; z < f4_NEVPARAM; z++) om->evparam[z] = gm->evparam[z];
  for (z = 0; z < f4_NCUTOFFS; z++) om->cutoff[z]  = gm->cutoff[z];
  for (z = 0; z < f4_MAXABET;  z++) om->compo[z]   = gm->compo[z];

  return eslOK;

 ERROR:
  return status;
}

/* Function:  f4_oprofile_ReconfigMSVLength()
 * Synopsis:  Set the target sequence length of the MSVFilter part of the model.
 * Incept:    SRE, Tue Dec 16 13:39:17 2008 [Janelia]
 *
 * Purpose:   Given an  already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>, only for the part of the model that's used
 *            for the accelerated MSV filter.
 *            
 *            The acceleration pipeline uses this to defer reconfiguring the
 *            length distribution of the main model, mostly because hmmscan
 *            reads the model in two pieces, MSV part first, then the rest.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_oprofile_ReconfigMSVLength(F4_OPROFILE *om, int L)
{
  om->tjb_b = unbiased_byteify(om, logf(3.0f / (float) (L+3)));
  return eslOK;
}

/* Function:  f4_oprofile_ReconfigRestLength()
 * Synopsis:  Set the target sequence length of the main profile.
 * Incept:    SRE, Tue Dec 16 13:41:30 2008 [Janelia]
 *
 * Purpose:   Given an  already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>, for everything except the MSV filter part
 *            of the model.
 *            
 *            Calling <f4_oprofile_ReconfigMSVLength()> then
 *            <f4_oprofile_ReconfigRestLength()> is equivalent to
 *            just calling <f4_oprofile_ReconfigLength()>. The two
 *            part version is used in the acceleration pipeline.
 *
 * Returns:   <eslOK> on success.           
 */
int
f4_oprofile_ReconfigRestLength(F4_OPROFILE *om, int L)
{
  float pmove, ploop;
  
  pmove = (2.0f + om->nj) / ((float) L + 2.0f + om->nj); /* 2/(L+2) for sw; 3/(L+3) for fs */
  ploop = 1.0f - pmove;

  /* ForwardFilter() parameters: pspace floats */
  om->xf[f4O_N][f4O_LOOP] =  om->xf[f4O_C][f4O_LOOP] = om->xf[f4O_J][f4O_LOOP] = ploop;
  om->xf[f4O_N][f4O_MOVE] =  om->xf[f4O_C][f4O_MOVE] = om->xf[f4O_J][f4O_MOVE] = pmove;

  /* ViterbiFilter() parameters: lspace signed 16-bit ints */
  om->xw[f4O_N][f4O_MOVE] =  om->xw[f4O_C][f4O_MOVE] = om->xw[f4O_J][f4O_MOVE] = wordify(om, logf(pmove));
  /* om->xw[p7O_N][p7O_LOOP] =  om->xw[p7O_C][p7O_LOOP] = om->xw[p7O_J][p7O_LOOP] = wordify(om, logf(ploop)); */ /* 3nat approx in force: these stay 0 */
  /* om->ncj_roundoff        = (om->scale_w * logf(ploop)) - om->xw[p7O_N][p7O_LOOP];                         */ /* and this does too                  */

  om->L = L;
  return eslOK;
}

/* Function:  f4_oprofile_ReconfigLength()
 * Synopsis:  Set the target sequence length of a model.
 * Incept:    SRE, Thu Dec 20 09:56:40 2007 [Janelia]
 *
 * Purpose:   Given an already configured model <om>, quickly reset its
 *            expected length distribution for a new mean target sequence
 *            length of <L>. 
 *            
 *            This doesn't affect the length distribution of the null
 *            model. That must also be reset, using <f4_bg_SetLength()>.
 *            
 *            We want this routine to run as fast as possible, because
 *            this call is in the critical path: it must be called at
 *            each new target sequence in a database search.
 *
 * Returns:   <eslOK> on success. Costs/scores for N,C,J transitions are set
 *            here.
 */
int
f4_oprofile_ReconfigLength(F4_OPROFILE *om, int L)
{
  int status;
  if ((status = f4_oprofile_ReconfigMSVLength (om, L)) != eslOK) return status;
  if ((status = f4_oprofile_ReconfigRestLength(om, L)) != eslOK) return status;
  return eslOK;
}