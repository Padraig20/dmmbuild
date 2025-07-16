/* F4_FWDBCK, estimating parameter counts via Baum-Welch Expectation Maximization.
 *
 * Contents:
 *   1. The Baum-Welch structure.
 *   2. Letter probabilities and background probabilities.
 *   3. Forward and backward procedures.
 *   4. Parameter estimation.
 *   5. Estimating parameters from a traceback.
 */

#include "dummer.h"

/*****************************************************************
 * 1. The Baum-Welch structure.
 *****************************************************************/

 /* Function:  bw_build()
  * 
  * Purpose:  Allocate and initialize the data structures needed for
  *           the Baum-Welch algorithm.
  * 
  * Args:     M - length of the profile (number of states in the HMM)
  *           L - length of the sequence
  *           ret_W_bar, ret_Y_bar, ret_Z_bar - RETURN: pointers to the W_bar, Y_bar, Z_bar matrices
  *           ret_X, ret_Y, ret_Z - RETURN: pointers to the X, Y, Z matrices
  * 
  * Returns:  <eslOK> on success.
  *           <eslEMEM> on memory allocation failure.
  */
int
bw_build(int M, int L, 
       double ***ret_W_bar, double ***ret_Y_bar, double ***ret_Z_bar, 
       double ***ret_X, double ***ret_Y, double ***ret_Z)
{
  double **W_bar = malloc((M+2) * sizeof(double*));
  double **Y_bar = malloc((M+2) * sizeof(double*));
  double **Z_bar = malloc((M+2) * sizeof(double*));
  double **X     = malloc((M+2) * sizeof(double*));
  double **Y     = malloc((M+2) * sizeof(double*));
  double **Z     = malloc((M+2) * sizeof(double*));

  if (!W_bar || !Y_bar || !Z_bar || !X || !Y || !Z) {
    free(W_bar); free(Y_bar); free(Z_bar); free(X); free(Y); free(Z);
    return eslEMEM;
  }

  for (int k = 0; k <= M+1; k++) {
    W_bar[k] = calloc(L+2, sizeof(double));
    Y_bar[k] = calloc(L+2, sizeof(double));
    Z_bar[k] = calloc(L+2, sizeof(double));
    X[k]     = calloc(L+2, sizeof(double));
    Y[k]     = calloc(L+2, sizeof(double));
    Z[k]     = calloc(L+2, sizeof(double));

    if (!W_bar[k] || !Y_bar[k] || !Z_bar[k] || !X[k] || !Y[k] || !Z[k]) {
      for (int j = 0; j <= k; j++) {
        free(W_bar[j]); free(Y_bar[j]); free(Z_bar[j]);
        free(X[j]);     free(Y[j]);     free(Z[j]);
      }
      free(W_bar); free(Y_bar); free(Z_bar); free(X); free(Y); free(Z);
      return eslEMEM;
    }
  }

  if (ret_W_bar) *ret_W_bar = W_bar;
  if (ret_Y_bar) *ret_Y_bar = Y_bar;
  if (ret_Z_bar) *ret_Z_bar = Z_bar;
  if (ret_X)     *ret_X     = X;
  if (ret_Y)     *ret_Y     = Y;
  if (ret_Z)     *ret_Z     = Z;

  return eslOK;
}

/* Function: bw_destroy()
 * 
 * Purpose:  Free the data structures allocated by <bw_build()>.
 * 
 * Args:     M - length of the profile (number of states in the HMM)
 *           W_bar, Y_bar, Z_bar - matrices to free
 *           X, Y, Z - matrices to free
 * 
 * Returns:  (void)
 */
void
bw_destroy(int M, 
        double **W_bar, double **Y_bar, double **Z_bar, 
        double **X, double **Y, double **Z)
{
  if (W_bar) { for (int k = 0; k <= M+1; k++) free(W_bar[k]); free(W_bar); }
  if (Y_bar) { for (int k = 0; k <= M+1; k++) free(Y_bar[k]); free(Y_bar); }
  if (Z_bar) { for (int k = 0; k <= M+1; k++) free(Z_bar[k]); free(Z_bar); }
  if (X)     { for (int k = 0; k <= M+1; k++) free(X[k]);     free(X); }
  if (Y)     { for (int k = 0; k <= M+1; k++) free(Y[k]);     free(Y); }
  if (Z)     { for (int k = 0; k <= M+1; k++) free(Z[k]);     free(Z); }
}

/* Function: bw_zero()
 *
 * Purpose:  Initialize the W_bar, Y_bar, Z_bar, X, Y, Z matrices to zero.
 * 
 * Args:     M - length of the profile (number of states in the HMM)
 *           L - length of the sequence
 *           W_bar, Y_bar, Z_bar - matrices to zero
 *           X, Y, Z - matrices to zero
 * 
 * Returns:  <eslOK> on success.
 */
int
bw_zero(int M, int L, 
        double **W_bar, double **Y_bar, double **Z_bar, 
        double **X, double **Y, double **Z)
{
  for (int k = 0; k <= M+1; k++) {
    for (int l = 0; l <= L+1; l++) {
      W_bar[k][l] = 0.0;
      Y_bar[k][l] = 0.0;
      Z_bar[k][l] = 0.0;
      X[k][l]     = 0.0;
      Y[k][l]     = 0.0;
      Z[k][l]     = 0.0;
    }
  }
  return eslOK;
}

/*****************************************************************
 * 2. Letter probabilities and background probabilities.
 *****************************************************************/

/* Function: letter_probs_build()
 * 
 * Purpose:  Allocate and initialize the letter probabilities and background probabilities.
 *           The letter probabilities are a 2D array of size (profile_length+1) x (num_of_letters),
 *           where each row corresponds to a position in the profile and each column corresponds
 *           to a letter in the alphabet. The background probabilities are a 1D array of size (num_of_letters).
 * 
 * Args:     profile_length - length of the profile (number of states in the HMM)
 *           num_of_letters - number of letters in the alphabet
 *           ret_letter_probs - RETURN: pointer to the letter probabilities array
 *           ret_background_probs - RETURN: pointer to the background probabilities array
 * 
 * Returns:  <eslOK> on success.
 *           <eslEMEM> on memory allocation failure.
 */
int
letter_probs_build(int profile_length, int num_of_letters, double ***ret_letter_probs, double **ret_background_probs)
{
  double **letter_probs    = malloc((profile_length+1) * sizeof(double*));
  double *background_probs = calloc(num_of_letters, sizeof(double));


  if (!letter_probs || !background_probs) {
    free(letter_probs); free(background_probs);
    return eslEMEM;
  }

  for (int i = 0; i < profile_length+1; i++) {
    letter_probs[i] = calloc(num_of_letters, sizeof(double));
    if (!letter_probs[i]) {
      for (int j = 0; j < i; j++) free(letter_probs[j]);
      free(letter_probs); free(background_probs);
      return eslEMEM;
    }
  }

  if (ret_letter_probs) *ret_letter_probs = letter_probs;
  if (ret_background_probs) *ret_background_probs = background_probs;

  return eslOK;
}

/* Function: letter_probs_destroy()
 * 
 * Purpose:  Free the letter probabilities and background probabilities.
 * 
 * Args:     profile_length - length of the profile (number of states in the HMM)
 *           letter_probs - letter probabilities array to free
 *           background_probs - background probabilities array to free
 * 
 * Returns:  (void)
 */
void
letter_probs_destroy(int profile_length, double **letter_probs, double *background_probs)
{
  if (letter_probs) { for (int i = 0; i < profile_length+1; i++) free(letter_probs[i]); free(letter_probs); }
  if (background_probs) free(background_probs);
}

/* Function: letter_probs_normalize()
 * 
 * Purpose:  Calculate the letter probabilities and background probabilities by their counts.
 *           Each row of the letter probabilities is normalized to sum to 1,
 *           and the background probabilities are normalized to sum to 1.
 * 
 * Args:     profile_length - length of the profile (number of states in the HMM)
 *           num_of_letters - number of letters in the alphabet
 *           letter_probs - holds counts of letters at each position in the profile
 *           background_probs - holds counts of letters of all positions in the profile
 * 
 * Returns:  <eslOK> on success.
 */
int
letter_probs_normalize(int profile_length, int num_of_letters, double **letter_probs, double *background_probs)
{
  for (int i = 0; i < profile_length; i++) {
    double row_sum = 0.0;
    for (int k = 0; k < num_of_letters; k++) row_sum += letter_probs[i][k];
    if (row_sum > 0.0) {
      for (int k = 0; k < num_of_letters; k++) letter_probs[i][k] /= row_sum;
    }
  }
  double bg_sum = 0.0;
  for (int k = 0; k < num_of_letters; k++) bg_sum += background_probs[k];
  if (bg_sum > 0.0) {
    for (int k = 0; k < num_of_letters; k++) background_probs[k] /= bg_sum;
  }
  return eslOK;
}

/* Function: letter_probs_count()
 * 
 * Purpose:  Count the occurrences of letters in the sequence and update the letter probabilities
 *           and background probabilities accordingly.
 * 
 * Args:     profile_length - length of the profile (number of states in the HMM)
 *           dsq - sequence to count letters from
 *           wt - weight for the counts
 *           letter_probs - holds counts of letters at each position in the profile
 *           background_probs - holds counts of letters of all positions in the profile
 * 
 * Returns:  <eslOK> on success.
 */
int
letter_probs_count(int profile_length, ESL_DSQ *dsq, float wt, double **letter_probs, double *background_probs)
{
  for (int i = 0; i < profile_length; i++) {
    letter_probs[i][dsq[i]]  += wt;
    background_probs[dsq[i]] += wt;
  }
  return eslOK;
}

/*****************************************************************
 * 3. Forward and backward procedures.
 *****************************************************************/

/* Function: f4_fwd()
 *
 * Purpose:  Perform the forward algorithm for the f4-HMM.
 *           Corresponds to "Algorithm 3" in the paper.
 * 
 * Args:     hmm - the f4-HMM model
 *           dsq - the sequence to align
 *           wt - weight for the sequence
 *           tr - the traceback structure (not used in this function)
 *           X, Y, Z - matrices for the forward algorithm (output stored here)
 *           M - length of the profile (number of states)
 *           N - length of the sequence
 *           letter_probs - letter probabilities for the sequence
 *           background_probs - background probabilities for the alphabet
 *           w_sum - pointer to accumulate the sum of weights
 * 
 * Returns:  <eslOK> on success.
 *           <eslEMEM> on memory allocation failure.
 * 
 * Note:     Note that the X, Y and Z matrices in Algorithm 3 also include i,j = -1
 *           and i = m, j = n. To account for that, index 0 in the array corresponds
 *           to -1 in the algorithm, and index m/n corresponds to m+1/n+1. I.e. the
 *           array was "shifted" towards the right.
 */
int
f4_fwd(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr, double **X, double **Y, double **Z, int M, int N, double **letter_probs, double *background_probs, double *w_sum)
{
  double S;                                                              // score
  double a_prime, b_prime, d_prime, e_prime;                             // parameters for BW
  double alpha_prob, beta_prob, delta_prob, epsilon_prob, epsilon_prob1; // counts to probabilities
  double w;

  for (int i = 1; i <= M+1; i++) {
    for (int j = 1; j <= N+1; j++) {

      alpha_prob   = hmm->tp[i-1][f4H_ALPHA]   / (hmm->tp[i-1][f4H_ALPHA]   + hmm->tp[i-1][f4H_DELTA] + hmm->tp[i-1][f4H_GAMMA]);
      beta_prob    = hmm->tp[i-1][f4H_BETA]    / (hmm->tp[i-1][f4H_BETA]    + hmm->tp[i-1][f4H_BETAP]);
      delta_prob   = hmm->tp[i-1][f4H_DELTA]   / (hmm->tp[i-1][f4H_ALPHA]   + hmm->tp[i-1][f4H_DELTA] + hmm->tp[i-1][f4H_GAMMA]);
      epsilon_prob = hmm->tp[i-1][f4H_EPSILON] / (hmm->tp[i-1][f4H_EPSILON] + hmm->tp[i-1][f4H_EPSILONP]);

      if (isnan(alpha_prob))   alpha_prob   = 0.0;
      if (isnan(beta_prob))    beta_prob    = 0.0;
      if (isnan(delta_prob))   delta_prob   = 0.0;
      if (isnan(epsilon_prob)) epsilon_prob = 0.0;

      a_prime = alpha_prob * (1.0 - beta_prob);
      b_prime = beta_prob;

      if (i == M+1) { // can use arbitrary values for S_m, d_m, e_m
        S       = 0.0;
        d_prime = 0.0;
        e_prime = 0.0;
      } else {
        epsilon_prob1 = hmm->tp[i][f4H_EPSILONP] / (hmm->tp[i][f4H_EPSILON] + hmm->tp[i][f4H_EPSILONP]);
        if (isnan(epsilon_prob1)) epsilon_prob1 = 0.0;

        S = (1 - alpha_prob - delta_prob) * letter_probs[i-1][dsq[j-1]] / background_probs[dsq[j-1]];

        d_prime = delta_prob   * (1 - epsilon_prob1);
        e_prime = epsilon_prob * (1 - epsilon_prob1) / (1 - epsilon_prob);
        if (isnan(e_prime)) e_prime = 0.0; // Avoid NaN issues
      }

      w = X[i-1][j-1] + Y[i-1][j] + Z[i][j-1] + 1;
      *w_sum += w;

      /* Sanity check: ensure no NaN values in the parameters */
      /*if (isnan(w) || isnan(S) || isnan(a_prime) || isnan(b_prime) || isnan(d_prime) || isnan(e_prime)) {
        fprintf(stderr, "Error: NaN detected at i=%d, j=%d (w=%f, S=%f, a_prime=%f, b_prime=%f, d_prime=%f, e_prime=%f)\n", 
            i, j, w, S, a_prime, b_prime, d_prime, e_prime);
        abort();
      }
      */
      
      X[i][j] = S * w;
      Y[i][j] = d_prime * w + e_prime * Y[i-1][j];
      Z[i][j] = a_prime * w + b_prime * Z[i][j-1];
    }
  }

  return eslOK;
}

/* Function: f4_bwd()
 *
 * Purpose:  Perform the backward algorithm for the f4-HMM.
 *           Corresponds to "Algorithm S5" in the paper.
 * 
 * Args:     hmm - the f4-HMM model
 *           dsq - the sequence to align
 *           wt - weight for the sequence
 *           tr - the traceback structure (not used in this function)
 *           W_bar, Y_bar, Z_bar - matrices for the backward algorithm (output stored here)
 *           M - length of the profile (number of states)
 *           N - length of the sequence
 *           letter_probs - letter probabilities for the sequence
 *           background_probs - background probabilities for the alphabet
 * 
 * Note:     Note that the X, Y and Z matrices in Algorithm S5 also include i,j = -1
 *           and i = m, j = n. To account for that, index 0 in the array corresponds
 *           to -1 in the algorithm, and index m/n corresponds to m+1/n+1. I.e. the
 *           array was "shifted" towards the right.
 */
int
f4_bwd(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr, double **W_bar, double **Y_bar, double **Z_bar, int M, int N, double **letter_probs, double *background_probs)
{
  double S;                                                              // score
  double a_prime, b_prime, d_prime, e_prime;                             // parameters for BW
  double alpha_prob, beta_prob, delta_prob, epsilon_prob, epsilon_prob1; // counts to probabilities
  double x;

  for (int i = M; i >= 0; i--) {
    for (int j = N; j >= 0; j--) {

      // calculate probabilities based on transition counts
      alpha_prob   = hmm->tp[i][f4H_ALPHA]   / (hmm->tp[i][f4H_ALPHA]   + hmm->tp[i][f4H_DELTA] + hmm->tp[i][f4H_GAMMA]);
      beta_prob    = hmm->tp[i][f4H_BETA]    / (hmm->tp[i][f4H_BETA]    + hmm->tp[i][f4H_BETAP]);
      delta_prob   = hmm->tp[i][f4H_DELTA]   / (hmm->tp[i][f4H_ALPHA]   + hmm->tp[i][f4H_DELTA] + hmm->tp[i][f4H_GAMMA]);
      epsilon_prob = hmm->tp[i][f4H_EPSILON] / (hmm->tp[i][f4H_EPSILON] + hmm->tp[i][f4H_EPSILONP]);

      if (isnan(alpha_prob))   alpha_prob   = 0.0;
      if (isnan(beta_prob))    beta_prob    = 0.0;
      if (isnan(delta_prob))   delta_prob   = 0.0;
      if (isnan(epsilon_prob)) epsilon_prob = 0.0;

      a_prime = alpha_prob * (1.0 - beta_prob);
      b_prime = beta_prob;

      if (i == M+1) { // can use arbitrary values for S_m, d_m, e_m
        S       = 0.0;
        d_prime = 0.0;
        e_prime = 0.0;
      } else {
        epsilon_prob1 = hmm->tp[i][f4H_EPSILONP] / (hmm->tp[i][f4H_EPSILON] + hmm->tp[i][f4H_EPSILONP]);
        if (isnan(epsilon_prob1)) epsilon_prob1 = 0.0;

        S = (1 - alpha_prob - delta_prob) * letter_probs[i][dsq[j]] / background_probs[dsq[j]];

        d_prime = delta_prob   * (1 - epsilon_prob1);
        e_prime = epsilon_prob * (1 - epsilon_prob1) / (1 - epsilon_prob);
        if (isnan(e_prime)) e_prime = 0.0; // avoid NaN issues
      }

      x = S * W_bar[i+1][j+1];

      /* Sanity check: ensure no NaN values in the parameters */
      /*
      if (isnan(x) || isnan(S) || isnan(a_prime) || isnan(b_prime) || isnan(d_prime) || isnan(e_prime)) {
        printf("Error: NaN detected at i=%d, j=%d (w=%f, S=%f, a_prime=%f, b_prime=%f, d_prime=%f, e_prime=%f)\n", 
            i, j, x, S, a_prime, b_prime, d_prime, e_prime);
        abort();
      }
      */
      
      W_bar[i][j] = x + d_prime * Y_bar[i+1][j] + a_prime * Z_bar[i][j+1] + 1;
      Y_bar[i][j] = W_bar[i][j] + e_prime * Y_bar[i+1][j];
      Z_bar[i][j] = W_bar[i][j] + b_prime * Z_bar[i][j+1];
    }
  }

  return eslOK;
}

/*****************************************************************
 * 4. Parameter estimation.
 *****************************************************************/

/* Function: f4_calculate_parameters()
 *
 * Purpose:  Estimate the parameters of the f4-HMM using the counts from the forward and backward procedures.
 *           Essentially, we calculate the expected counts of the parameters based on the forward and backward matrices.
 * 
 * Args:     hmm - the f4-HMM model (output will be stored in hmm->tp)
 *           dsq - the sequence to align
 *           wt - weight for the sequence
 *           tr - the traceback structure (not used in this function)
 *           W_bar, Y_bar, Z_bar - matrices for the forward algorithm
 *           X, Y, Z - matrices for the backward algorithm
 *           v - the sum of weights from the forward algorithm
 *           letter_probs - letter probabilities for the sequence
 *           termination_condition - pointer to an integer that will be updated with the termination condition
 * 
 * Returns:  <eslOK> on success.
 * 
 * Note:     The termination condition is updated to indicate whether the parameters have
 *           converged under a specific threshold.
 */
int
f4_calculate_parameters(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr, 
  double **W_bar, double **Y_bar, double **Z_bar,
  double **X, double **Y, double **Z,
  double v, double **letter_probs,
  int *termination_condition)
{
  int M = hmm->M;
  int N = tr->L;

  double alpha, beta, delta, epsilon;
  double gamma, betap, epsilonp;
  double epsilonp_i1, epsilon_i1;

  // estimate counts for every position
  for (int i = 1; i <= M+1; i++) {
    // expected count of (1 - beta)
    betap = 0.0;
    for (int j = 1; j <= N+1; j++) {
      betap += Z[i][j-1] * W_bar[i][j];
    }
    betap /= v;

    // expected count of beta
    beta = 0.0;
    for (int j = 1; j <= N+1; j++) {
      beta += Z[i][j-1] * Z_bar[i][j];
    }
    beta /= v;
    beta -= betap;

    // expected count of (1 - epsilon)
    epsilonp = 0.0;
    for (int j = 1; j <= N+1; j++) {
      epsilonp += Y[i-1][j] * W_bar[i][j];
    }
    epsilonp /= v;

    // expected count of epsilon
    epsilon = 0.0;
    for (int j = 1; j <= N+1; j++) {
      epsilon += Y[i-1][j] * Y_bar[i][j];
    }
    epsilon /= v;
    epsilon -= epsilonp;

    // expected count of alpha
    alpha = betap;

    // expected count of delta
    delta = 0.0;
    // requires epsilonp_i+1 and epsilon_i+1
    epsilonp_i1 = 0.0;
    epsilon_i1  = 0.0;
    if (i < M+1) { // epsilon_m may be arbitrary, therefore 0 is fine
      // expected count of epsilonp_i+1
      for (int j = 1; j <= N+1; j++) {
        epsilonp_i1 += Y[i][j] * W_bar[i+1][j]; // i++
      }
      epsilonp_i1 /= v;
      // expected count of epsilon_i+1
      for (int j = 1; j <= N+1; j++) {
        epsilon_i1 += Y[i][j] * Y_bar[i+1][j]; // i++
      }
      epsilon_i1 /= v;
      epsilon_i1 -= epsilonp_i1;
    }
    // now we can calculate delta
    delta = epsilonp_i1 + epsilon_i1 - epsilon;

    // expected count of gamma
    gamma = 0.0;
    // requires expected counts of letter_probs
    for (int letter = 0; letter < hmm->abc->K; letter++) {
      double letter_sum = 0.0;
      for (int j = 1; j <= N; j++) {
        if (dsq[j-1] == letter && i < M+1) {
          letter_sum += X[i][j] * W_bar[i+1][j+1];
        }
      }
      letter_sum /= v;
      // now we can calculate gamma
      gamma += letter_sum;
    }

    double threshold = 1e-5;
    *termination_condition &= 
      fabs(alpha    * wt - hmm->tp[i-1][f4H_ALPHA])    < threshold &&
      fabs(beta     * wt - hmm->tp[i-1][f4H_BETA])     < threshold &&
      fabs(delta    * wt - hmm->tp[i-1][f4H_DELTA])    < threshold &&
      fabs(epsilon  * wt - hmm->tp[i-1][f4H_EPSILON])  < threshold &&
      fabs(gamma    * wt - hmm->tp[i-1][f4H_GAMMA])    < threshold &&
      fabs(betap    * wt - hmm->tp[i-1][f4H_BETAP])    < threshold &&
      fabs(epsilonp * wt - hmm->tp[i-1][f4H_EPSILONP]) < threshold;

    // update the HMM parameters
    hmm->tp[i-1][f4H_ALPHA]    = alpha    * wt;
    hmm->tp[i-1][f4H_BETA]     = beta     * wt;
    hmm->tp[i-1][f4H_DELTA]    = delta    * wt;
    hmm->tp[i-1][f4H_EPSILON]  = epsilon  * wt;
    hmm->tp[i-1][f4H_GAMMA]    = gamma    * wt;
    hmm->tp[i-1][f4H_BETAP]    = betap    * wt;
    hmm->tp[i-1][f4H_EPSILONP] = epsilonp * wt;
  }
  
  return eslOK;
}

/*****************************************************************
 * 5. Estimating parameters from a traceback.
 *****************************************************************/

/* Function: f4_trace_Estimate()
 * 
 * Purpose:  Instead of directly counting a traceback into a count-based
 *           core HMM structure, this function estimates the parameters
 *           of the core HMM from a traceback via the Baum-Welch
 *           algorithm.
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
 *
 * Notes:    We stop after 50 iterations or when the termination condition is met.
 *           Possibly, therefore, the parameters may not be fully converged.
 */
int
f4_trace_Estimate(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr, double **letter_probs, double *background_probs)
{
  int M = hmm->M; // profile size
  int N = tr->L;  // sequence length

  double **W_bar, **Y_bar, **Z_bar;  // backward
  double **X, **Y, **Z;              // forward
  double v;                          // variable to accumulate the sum of weights

  bw_build(M, N, &W_bar, &Y_bar, &Z_bar, &X, &Y, &Z);

  int termination_condition = 1;
  int num_iterations        = 0;

  do {

    v = 0.0;

    if (bw_zero(M, N, W_bar, Y_bar, Z_bar, X, Y, Z) != eslOK) {
      bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);
      return eslEMEM;
    }

    /* Forward pass, calculate X, Y, Z */
    if (f4_fwd(hmm, dsq, wt, tr, X, Y, Z, M, N, letter_probs, background_probs, &v) != eslOK) {
      bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);
      return eslEINVAL;
    }

    /* Backward pass, calculate W_bar, Y_bar, Z_bar */
    if (f4_bwd(hmm, dsq, wt, tr, W_bar, Y_bar, Z_bar, M, N, letter_probs, background_probs) != eslOK) {
      bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);
      return eslEINVAL;
    }

    /* Calculate and update parameters in hmm, determine termination condition */
    if (f4_calculate_parameters(hmm, dsq, wt, tr, W_bar, Y_bar, Z_bar, X, Y, Z, v, letter_probs, &termination_condition) != eslOK) {
      bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);
      return eslEINVAL;
    }

    num_iterations++;

  } while (!termination_condition && num_iterations < 50);

  bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);
  
  return eslOK;
}