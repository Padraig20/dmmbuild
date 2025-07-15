/* Creating profile HMMs from single sequences.
 * 
 * Contents:
 *   1. Routines in the exposed API.
 */

#include "dummer.h"

/*****************************************************************
 * 1. Routines in the exposed API.
 *****************************************************************/


/* Function:  f4_Seqmodel()
 * Synopsis:  Make a profile HMM from a single sequence.
 *
 * Purpose:   Make a profile HMM from a single sequence, for
 *            probabilistic Smith/Waterman alignment, DUMMER-style.
 *            
 *            The query is digital sequence <dsq> of length <M>
 *            residues in alphabet <abc>, named <name>. 
 *            
 *            The scoring system is given by <Q>, <f>, <popen>, and
 *            <pextend>. <Q> is a $K \times K$ matrix giving
 *            conditional residue probabilities $P(a \mid b)}$; these
 *            are typically obtained by reverse engineering a score
 *            matrix like BLOSUM62. <f> is a vector of $K$ background
 *            frequencies $p_a$. <popen> and <pextend> are the
 *            probabilities assigned to gap-open ($t_{MI}$ and
 *            $t_{MD}$) and gap-extend ($t_{II}$ and $t_{DD}$)
 *            transitions.
 *            
 * Args:      
 *
 * Returns:   <eslOK> on success, and a newly allocated HMM is returned
 *            in <ret_hmm>. 
 *
 * Throws:    <eslEMEM> on allocation error, and <*ret_hmm> is <NULL>.
 */
int
f4_Seqmodel(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
	    ESL_DMATRIX *Q, float *f, double popen, double pextend,
	    F4_HMM **ret_hmm)
{
  int     status;
  F4_HMM *hmm    = NULL;
  char   *logmsg = "[HMM created from a query sequence]";
  int     k;

  if ((hmm = f4_hmm_Create(M, abc)) == NULL) { status = eslEMEM; goto ERROR; }
  
  for (k = 0; k <= M; k++)
    {
      /* Use rows of P matrix as source of match emission vectors */
      if (k > 0) esl_vec_D2F(Q->mx[(int) dsq[k]], abc->K, hmm->mat[k]);

      /* Set inserts to background for now. This will be improved. */
      esl_vec_FCopy(f, abc->K, hmm->ins[k]);

      hmm->t[k][f4H_MM] = 1.0 - 2 * popen;
      hmm->t[k][f4H_MI] = popen;
      hmm->t[k][f4H_MD] = popen;
      hmm->t[k][f4H_IM] = 1.0 - pextend;
      hmm->t[k][f4H_II] = pextend;
      hmm->t[k][f4H_DM] = 1.0 - pextend;
      hmm->t[k][f4H_DD] = pextend;
    }

  /* Deal w/ special stuff at node M, overwriting a little of what we
   * just did. 
   */
  hmm->t[M][f4H_MM] = 1.0 - popen;
  hmm->t[M][f4H_MD] = 0.;
  hmm->t[M][f4H_DM] = 1.0;
  hmm->t[M][f4H_DD] = 0.;
  
  /* Add mandatory annotation
   */
  f4_hmm_SetName(hmm, name);
  f4_hmm_AppendComlog(hmm, 1, &logmsg);
  hmm->nseq     = 1;
  f4_hmm_SetCtime(hmm);
  hmm->checksum = 0;

  *ret_hmm = hmm;
  return eslOK;
  
 ERROR:
  if (hmm != NULL) f4_hmm_Destroy(hmm);
  *ret_hmm = NULL;
  return status;
}