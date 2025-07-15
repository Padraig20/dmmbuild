/* SSE implementation of an optimized profile structure.
 * 
 * Contents:
 *   1. The P7_OMX structure: a dynamic programming matrix
 * 
 * See also:
 *   p7_omx.ai - figure illustrating the layout of a P7_OMX.
 *
 * SRE, Sun Nov 25 11:26:48 2007 [Casa de Gatos]
 */

#include "dummer.h"

/* Function:  f4_omx_Destroy()
 * Synopsis:  Frees an optimized DP matrix.
 * Incept:    SRE, Tue Nov 27 09:11:42 2007 [Janelia]
 *
 * Purpose:   Frees optimized DP matrix <ox>.
 *
 * Returns:   (void)
 */
void
f4_omx_Destroy(F4_OMX *ox)
{
  if (ox == NULL) return;
  if (ox->x_mem   != NULL) free(ox->x_mem);
  if (ox->dp_mem  != NULL) free(ox->dp_mem);
  if (ox->dpf     != NULL) free(ox->dpf);
  if (ox->dpw     != NULL) free(ox->dpw);
  if (ox->dpb     != NULL) free(ox->dpb);
  free(ox);
  return;
}

/* Function:  f4_omx_Create()
 * Synopsis:  Create an optimized dynamic programming matrix.
 * Incept:    SRE, Tue Nov 27 08:48:20 2007 [Janelia]
 *
 * Purpose:   Allocates a reusable, resizeable <F4_OMX> for models up to
 *            size <allocM> and target sequences up to length
 *            <allocL/allocXL>, for use by any of the various optimized
 *            DP routines.
 *            
 *            To allocate the very memory-efficient one-row matrix
 *            used by *Filter() and *Score() functions that only
 *            calculate scores, <allocM=M>, <allocL=0>, and
 *            <allocXL=0>.
 *            
 *            To allocate the reasonably memory-efficient linear
 *            arrays used by *Parser() functions that only keep
 *            special (X) state scores, <allocM=M>, <allocL=0>,
 *            and <allocXL=L>.
 *            
 *            To allocate a complete matrix suitable for functions
 *            that need the whole DP matrix for traceback, sampling,
 *            posterior decoding, or reestimation, <allocM=M> and
 *            <allocL=allocXL=L>.
 *
 * Returns:   a pointer to the new <F4_OMX>.
 *
 * Throws:    <NULL> on allocation failure.
 */
F4_OMX *
f4_omx_Create(int allocM, int allocL, int allocXL)
{
  F4_OMX  *ox     = NULL;
  int      i;
  int      status;

  ESL_ALLOC(ox, sizeof(F4_OMX));
  ox->dp_mem = NULL;
  ox->dpb    = NULL;
  ox->dpw    = NULL;
  ox->dpf    = NULL;
  ox->xmx    = NULL;
  ox->x_mem  = NULL;

  /* DP matrix will be allocated for allocL+1 rows 0,1..L; allocQ4*f4X_NSCELLS columns */
  ox->allocR   = allocL+1;
  ox->validR   = ox->allocR;
  ox->allocQ4  = f4O_NQF(allocM);
  ox->allocQ8  = f4O_NQW(allocM);
  ox->allocQ16 = f4O_NQB(allocM);
  ox->ncells   = (int64_t) ox->allocR * (int64_t) ox->allocQ4 * 4;      /* # of DP cells allocated, where 1 cell contains MDI */

  ESL_ALLOC(ox->dp_mem, sizeof(__m128) * (int64_t) ox->allocR * (int64_t) ox->allocQ4 * f4X_NSCELLS + 15);  /* floats always dominate; +15 for alignment */
  ESL_ALLOC(ox->dpb,    sizeof(__m128i *) * ox->allocR);
  ESL_ALLOC(ox->dpw,    sizeof(__m128i *) * ox->allocR);
  ESL_ALLOC(ox->dpf,    sizeof(__m128  *) * ox->allocR);

  ox->dpb[0] = (__m128i *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
  ox->dpw[0] = (__m128i *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));
  ox->dpf[0] = (__m128  *) ( ( (unsigned long int) ((char *) ox->dp_mem + 15) & (~0xf)));

  for (i = 1; i <= allocL; i++) {
    ox->dpf[i] = ox->dpf[0] + (int64_t) i * (int64_t) ox->allocQ4  * f4X_NSCELLS;
    ox->dpw[i] = ox->dpw[0] + (int64_t) i * (int64_t) ox->allocQ8  * f4X_NSCELLS;
    ox->dpb[i] = ox->dpb[0] + (int64_t) i * (int64_t) ox->allocQ16;
  }

  ox->allocXR = allocXL+1;
  ESL_ALLOC(ox->x_mem,  sizeof(float) * ox->allocXR * f4X_NXCELLS + 15); 
  ox->xmx = (float *) ( ( (unsigned long int) ((char *) ox->x_mem  + 15) & (~0xf)));

  ox->M              = 0;
  ox->L              = 0;
  ox->totscale       = 0.0;
  ox->has_own_scales = TRUE;	/* most matrices are Forward, control their own scale factors */
  return ox;

 ERROR:
  f4_omx_Destroy(ox);
  return NULL;
}