/* F4_FWDBCK, estimating parameter counts via Baum-Welch Expectation Maximization.
 *
 * Contents:
 *   0. Helper Functions.
 *   1. The Baum-Welch structure.
 *   2. Letter probabilities and background probabilities.
 *   3. Forward and backward procedures.
 *   4. Parameter estimation.
 *   5. Estimating parameters from a traceback.
 */

#include "dummer.h"

/*****************************************************************
 * 0. Helper Functions.
 *****************************************************************/

/* Function: is_gap_or_missing()
 * 
 * Purpose:  Check if a letter is a gap or missing data in the alphabet.
 * 
 * Args:     letter - the letter to check
 *           abc - pointer to the alphabet structure
 * 
 * Returns:  1 (true) if the letter is a gap or missing data, 0 otherwise.
 * 
 * Notes:    We assume that the alphabet of the MSA is correct. If not,
 *           that's on the user :( Sorry about the memory leaks, but fix
 *           the MSA!
 */
int
is_gap_or_missing(ESL_DSQ letter, ESL_ALPHABET *abc) {
  if (letter >= abc->K && letter < abc->Kp) {
    return 1; // gap or missing data
  } else if (letter < 0 || letter >= abc->Kp) {
    fprintf(stderr, "Whoa. That does not look good. This letter should not even be in the alphabet: %d\n", letter);
    fprintf(stderr, "The alphabet size is %d, and the total size of the alphabet (including gaps, missing data, ...) is %d.\n", abc->K, abc->Kp);
    fprintf(stderr, "Hell breaks loose. This script will abort now. Sorry about that. Fix your MSA and come back!\n");
    abort(); // this is a serious error, we can not continue
  }
  return 0; // valid symbol
}

/* Function: get_sequence_without_gaps()
 * 
 * Purpose:  Extract the sequence from the MSA without gaps or missing data.
 * 
 * Args:     msa - pointer to the MSA structure
 *           ret_dsq - return pointer to the array of sequences
 *                     without gaps [0..nseq-1][1..length]
 *           ret_length - return pointer to the array of lengths of the
 *                        sequences without gaps [0..nseq-1]
 * 
 * Returns:  <eslOK> on success.
 *           <eslEMEM> on memory allocation failure.
 * 
 * Notes:    Both arrays will be allocated and must be freed by the *caller*.
 */
int
get_sequence_without_gaps(ESL_MSA *msa, ESL_DSQ ***ret_dsq, int **ret_length) {
  int i, j, k;
  int length;
  ESL_DSQ *dsq;

  *ret_dsq    = malloc(msa->nseq * sizeof(ESL_DSQ*));
  *ret_length = calloc(msa->nseq, sizeof(int));

  for (i = 0; i < msa->nseq; i++) {
    length = 0;

    for (j = 1; j <= msa->alen; j++)
      if (!is_gap_or_missing(msa->ax[i][j], msa->abc))
        length++;
    
    (*ret_length)[i] = length;
  }

  for (i = 0; i < msa->nseq; i++) {
    k = 0;
    dsq = malloc(sizeof(ESL_DSQ) * ((*ret_length)[i] + 1)); // +1 for null terminator
    if (dsq == NULL) goto ERROR;

    for (j = 1; j <= msa->alen; j++)
      if (!is_gap_or_missing(msa->ax[i][j], msa->abc))
        dsq[k++] = msa->ax[i][j];

    dsq[k] = eslDSQ_SENTINEL; // null terminator
    (*ret_dsq)[i] = dsq;
  }
  
  return eslOK;

 ERROR:
  if (dsq != NULL) free(dsq);
  return eslEMEM;   
}

/* Function: free_sequence_without_gaps()
 *
 * Purpose:  Free the memory allocated for the sequences without gaps
 *           as well as the lengths.
 * 
 * Args:     dsq - pointer to the array of sequences without gaps
 *           length - pointer to the array of lengths of the sequences
 *           nseq - number of sequences
 * 
 * Returns:  (void)
 */
void
free_sequence_without_gaps(ESL_DSQ **dsq, int *length, int nseq) {
  int i;
  for (i = 0; i < nseq; i++) {
    if (dsq[i] != NULL) free(dsq[i]);
  }
  free(dsq);
  free(length);
}

/*****************************************************************
 * 1. The Baum-Welch structure.
 *****************************************************************/

 /* Function: bw_build()
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
 *           N - length of the sequence
 *           W_bar, Y_bar, Z_bar - matrices to zero
 *           X, Y, Z - matrices to zero
 * 
 * Returns:  (void)
 */
void
bw_zero(int M, int N, 
        double **W_bar, double **Y_bar, double **Z_bar, 
        double **X, double **Y, double **Z)
{
  for (int k = 0; k <= M+1; k++) {
    for (int l = 0; l <= N+1; l++) {
      W_bar[k][l] = 0.0;
      Y_bar[k][l] = 0.0;
      Z_bar[k][l] = 0.0;
      X[k][l]     = 0.0;
      Y[k][l]     = 0.0;
      Z[k][l]     = 0.0;
    }
  }
}

/* Function: param_counts_save_from_to()
 *
 * Purpose:  Save the parameter counts into another HMM structure.
 * 
 * Args:     from - the HMM structure from which to save the parameter counts
 *           to   - the HMM structure where to save the parameter counts
 * 
 * Returns:  (void)
 */
void
param_counts_save_from_to(F4_HMM *from, F4_HMM *to)
{
  int k;

  for (k = 0; k <= to->M; k++) {
    to->tp[k][f4H_ALPHA]    = from->tp[k][f4H_ALPHA];
    to->tp[k][f4H_BETA]     = from->tp[k][f4H_BETA];
    to->tp[k][f4H_DELTA]    = from->tp[k][f4H_DELTA];
    to->tp[k][f4H_EPSILON]  = from->tp[k][f4H_EPSILON];
    to->tp[k][f4H_GAMMA]    = from->tp[k][f4H_GAMMA];
    to->tp[k][f4H_BETAP]    = from->tp[k][f4H_BETAP];
    to->tp[k][f4H_EPSILONP] = from->tp[k][f4H_EPSILONP];
  }
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
 *           ret_letter_probs - RETURN: pointer to the letter probabilities array [1..profile_length][0..K-1]
 *           ret_background_probs - RETURN: pointer to the background probabilities array [0..K-1]
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

  for (int i = 1; i <= profile_length; i++) {
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
  if (letter_probs) { for (int i = 1; i <= profile_length; i++) free(letter_probs[i]); free(letter_probs); }
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
  for (int i = 1; i <= profile_length; i++) {
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
  } else {
    fprintf(stderr, "Warning: Background probabilities sum to zero. This may indicate an issue with the input data.\n");
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
 *           letter_probs - holds counts of letters at each position in the profile [1..profile_length][0..K-1]
 *           background_probs - holds counts of letters of all positions in the profile [0..K-1]
 * 
 * Returns:  <eslOK> on success.
 *           <eslFAIL> on invalid symbol in the sequence OR profile has not been extracted successfully.
 */
int
letter_probs_count(int profile_length, int seq_length, ESL_DSQ *dsq, float wt, ESL_ALPHABET *abc,
                   double **letter_probs, double *background_probs, int *matassign)
{
  int idx = 0;

  for (int i = 1; i <= seq_length; i++) {

    if (matassign[i]) {
      ++idx;
    } else {
      continue; // skip positions not assigned to the profile
    }

    if (is_gap_or_missing(dsq[i], abc)) {
      continue; // this is just a gap or missing symbol, skip it
    }

    letter_probs[idx][dsq[i]]  += wt;
    background_probs[dsq[i]]   += wt;
  }

  if (idx != profile_length) {
    fprintf(stderr, "Error: Profile length (%d) does not match the number of assigned positions (%d).\n", profile_length, idx);
    return eslFAIL;
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
 *           <eslFAIL> on invalid symbol in the sequence.
 * 
 * Note:     Note that the X, Y and Z matrices in Algorithm 3 also include i,j = -1
 *           and i = m, j = n. To account for that, index 0 in the array corresponds
 *           to -1 in the algorithm, and index m/n corresponds to m+1/n+1. I.e. the
 *           array was "shifted" towards the right.
 */
int
f4_fwd(F4_HMM *hmm, ESL_DSQ *dsq, F4_TRACE *tr, ESL_ALPHABET *abc,
       double **X, double **Y, double **Z, int M, int N,
       double **letter_probs, double *background_probs, double *w_sum)
{
  double S;                                                              // score
  double a_prime, b_prime, d_prime, e_prime;                             // parameters for BW
  double alpha_prob, beta_prob, delta_prob, epsilon_prob, epsilon_prob1; // counts to probabilities
  double w;
  int letter;

  for (int i = 1; i <= M+1; i++) {

    alpha_prob    = hmm->tp[i-1][f4H_ALPHA];
    beta_prob     = hmm->tp[i-1][f4H_BETA];
    delta_prob    = hmm->tp[i-1][f4H_DELTA];
    epsilon_prob  = hmm->tp[i-1][f4H_EPSILON];

    epsilon_prob1 = (i == M+1) ? 0.0 : hmm->tp[i][f4H_EPSILON];

    a_prime = alpha_prob * (1.0 - beta_prob);
    b_prime = beta_prob;

    // can use arbitrary values for S_m, d_m, e_m
    d_prime = (i == M+1) ? 0.0 : delta_prob   * (1 - epsilon_prob1);
    e_prime = (i == M+1) ? 0.0 : epsilon_prob * (1 - epsilon_prob1) / (1 - epsilon_prob);
    if (isnan(e_prime)) e_prime = 0.0; // Avoid NaN issues

    for (int j = 1; j <= N+1; j++) {

      if (j == N+1) letter = 0; // can use arbitrary values for theta_n
      else letter = dsq[j-1];

      S = (i == M+1) ? 0.0 : (1 - alpha_prob - delta_prob) * letter_probs[i][letter] / background_probs[letter];

      w = X[i-1][j-1] + Y[i-1][j] + Z[i][j-1] + 1.0; // score is 1.0
      *w_sum += w;
      
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
 *           tr - the traceback structure (not used in this function)
 *           W_bar, Y_bar, Z_bar - matrices for the backward algorithm (output stored here)
 *           M - length of the profile (number of states)
 *           N - length of the sequence
 *           letter_probs - letter probabilities for the sequence
 *           background_probs - background probabilities for the alphabet
 * 
 * Returns:  <eslOK> on success.
 *           <eslEMEM> on memory allocation failure.
 *           <eslFAIL> on invalid symbol in the sequence.
 * 
 * Note:     Note that the X, Y and Z matrices in Algorithm S5 also include i,j = -1
 *           and i = m, j = n. To account for that, index 0 in the array corresponds
 *           to -1 in the algorithm, and index m/n corresponds to m+1/n+1. I.e. the
 *           array was "shifted" towards the right.
 */
int
f4_bwd(F4_HMM *hmm, ESL_DSQ *dsq, F4_TRACE *tr, ESL_ALPHABET *abc,
       double **W_bar, double **Y_bar, double **Z_bar, int M, int N,
       double **letter_probs, double *background_probs)
{
  double S;                                                              // score
  double a_prime, b_prime, d_prime, e_prime;                             // parameters for BW
  double alpha_prob, beta_prob, delta_prob, epsilon_prob, epsilon_prob1; // counts to probabilities
  double x;
  int letter;

  for (int i = M; i >= 0; i--) {

    alpha_prob   = hmm->tp[i][f4H_ALPHA];
    beta_prob    = hmm->tp[i][f4H_BETA];
    delta_prob   = hmm->tp[i][f4H_DELTA];
    epsilon_prob = hmm->tp[i][f4H_EPSILON];

    epsilon_prob1 = (i == M) ? 0.0 : hmm->tp[i+1][f4H_EPSILON];

    a_prime = alpha_prob * (1.0 - beta_prob);
    b_prime = beta_prob;

    // can use arbitrary values for S_m, d_m, e_m
    d_prime = (i == M) ? 0.0 : delta_prob   * (1 - epsilon_prob1);
    e_prime = (i == M) ? 0.0 : epsilon_prob * (1 - epsilon_prob1) / (1 - epsilon_prob);
    if (isnan(e_prime)) e_prime = 0.0; // avoid NaN issues

    for (int j = N; j >= 0; j--) {

      if (j == N) letter = 0; // can use arbitrary values for theta_n
      else letter = dsq[j];

      S = (i == M) ? 0.0 : (1 - alpha_prob - delta_prob) * letter_probs[i+1][letter] / background_probs[letter];

      x = S * W_bar[i+1][j+1];
      
      W_bar[i][j] = x + d_prime * Y_bar[i+1][j] + a_prime * Z_bar[i][j+1] + 1.0; // score is 1.0
      Y_bar[i][j] = W_bar[i][j] + e_prime * Y_bar[i+1][j];
      Z_bar[i][j] = W_bar[i][j] + b_prime * Z_bar[i][j+1];
    }
  }

  return eslOK;
}

int determine_termination_condition(F4_HMM *old, F4_HMM *new)
{

  for (int i = 0; i <= old->M; i++) {
    if (
      fabs(old->tp[i][f4H_ALPHA]    - new->tp[i][f4H_ALPHA])    >= f4_BW_CONVERGE ||
      fabs(old->tp[i][f4H_BETA]     - new->tp[i][f4H_BETA])     >= f4_BW_CONVERGE ||
      fabs(old->tp[i][f4H_DELTA]    - new->tp[i][f4H_DELTA])    >= f4_BW_CONVERGE ||
      fabs(old->tp[i][f4H_EPSILON]  - new->tp[i][f4H_EPSILON])  >= f4_BW_CONVERGE ||
      fabs(old->tp[i][f4H_GAMMA]    - new->tp[i][f4H_GAMMA])    >= f4_BW_CONVERGE ||
      fabs(old->tp[i][f4H_BETAP]    - new->tp[i][f4H_BETAP])    >= f4_BW_CONVERGE ||
      fabs(old->tp[i][f4H_EPSILONP] - new->tp[i][f4H_EPSILONP]) >= f4_BW_CONVERGE
    ) {
      return 0;  // not converged
    }
  }

  return 1;  // converged
}

/*****************************************************************
 * 4. Parameter estimation.
 *****************************************************************/

/* Function: f4_calculate_parameters()
 *
 * Purpose:  Estimate the parameters of the f4-HMM using the counts from the forward and backward procedures.
 *           Essentially, we calculate the expected counts of the parameters based on the forward and backward matrices.
 * 
 * Args:     hmm - the f4-HMM model (parameter transition probabilities are read from here)
 *           N - length of the sequence
 *           wt - weight for the sequence
 *           tr - the traceback structure (not used in this function)
 *           W_bar, Y_bar, Z_bar - matrices for the forward algorithm
 *           X, Y, Z - matrices for the backward algorithm
 *           v - the sum of weights from the forward algorithm
 *           letter_probs - letter probabilities for the sequence
 *           param_counts - matrix to accumulate the estimated parameters (output will be stored here)
 * 
 * Returns:  <eslOK> on success.
 */
int
f4_calculate_parameters(F4_HMM *hmm, int N, float wt, 
  double **W_bar, double **Y_bar, double **Z_bar,
  double **X, double **Y, double **Z,
  double v, double **letter_probs, F4_HMM *param_counts)
{
  int M = hmm->M;
  int i, j;

  double alpha, beta, delta, epsilon;
  double gamma, betap, epsilonp;
  double epsilonp_i1, epsilon_i1;

  // estimate counts for every position
  for (i = 1; i <= M+1; i++) {

    betap    = 0.0;  // expected count of (1 - beta)
    beta     = 0.0;
    epsilonp = 0.0;  // expected count of (1 - epsilon)
    epsilon  = 0.0;

    for (j = 1; j <= N+1; j++) {
      betap    += Z[i][j-1] * W_bar[i-1][j-1];
      beta     += Z[i][j-1] * Z_bar[i-1][j-1];
      epsilonp += Y[i-1][j] * W_bar[i-1][j-1];
      epsilon  += Y[i-1][j] * Y_bar[i-1][j-1];
    }

    betap /= v;
    beta  /= v;
    beta  -= betap;
    if (beta < -DBL_EPSILON) printf("Warning: Negative beta count at position %d.\n", i);

    epsilonp /= v;
    epsilon  /= v;
    epsilon  -= epsilonp;
    if (epsilon < -DBL_EPSILON) printf("Warning: Negative epsilon count at position %d.\n", i);

    // expected count of alpha
    alpha = betap;

    // expected count of gamma
    // gamma_m is arbitrary
    gamma = 0.0;

    // expected count of delta
    delta = 0.0;
    // requires epsilonp_i+1 and epsilon_i+1
    if (i < M+1) { // epsilon_m may be arbitrary, therefore 0 is fine
      epsilonp_i1 = 0.0;
      epsilon_i1  = 0.0;
      for (int j = 1; j <= N+1; j++) {
        epsilonp_i1 += Y[i][j] * W_bar[i][j-1]; // i++
        epsilon_i1  += Y[i][j] * Y_bar[i][j-1]; // i++
      }
      epsilonp_i1 /= v;
      epsilon_i1  /= v;
      epsilon_i1  -= epsilonp_i1;
      
      // now we can calculate delta
      delta = epsilonp_i1 + epsilon_i1 - epsilon;
      if (delta < -DBL_EPSILON) printf("Warning: Negative delta count at position %d. %g + %g - %g = %g\n", i,
				       epsilonp_i1, epsilon_i1, epsilon, delta);
      for (j = 1; j <= N; j++) {
        gamma += X[i][j] * W_bar[i][j];
      }
      gamma /= v;
    }

    // update the HMM parameters
    param_counts->tp[i-1][f4H_ALPHA]    += alpha    * wt;
    param_counts->tp[i-1][f4H_BETA]     += beta     * wt;
    param_counts->tp[i-1][f4H_DELTA]    += delta    * wt;
    param_counts->tp[i-1][f4H_EPSILON]  += epsilon  * wt;
    param_counts->tp[i-1][f4H_GAMMA]    += gamma    * wt;
    param_counts->tp[i-1][f4H_BETAP]    += betap    * wt;
    param_counts->tp[i-1][f4H_EPSILONP] += epsilonp * wt;
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
 *           msa   - Multiple Sequence Alignment (MSA) structure
 *           tr    - array of all alignments of seq to HMM
 *           pri   - prior probabilities for the HMM
 *           letter_probs - letter probabilities for the sequence
 *           background_probs - background probabilities for the alphabet
 *           
 * Return:   <eslOK> on success.
 *           Weighted count events are accumulated in hmm's mat[][], ins[][],
 *           tp[][] fields: the core probability model.
 *           
 * Throws:   <eslEINVAL> if something's corrupt in the trace; effect on hmm
 *           counts is undefined, because it may abort at any point in the trace.
 *           <eslEMEM> if memory allocation fails.
 *
 * Notes:    We stop after a certain number of iterations or when the termination condition is met.
 *           Possibly, therefore, the parameters may not be fully converged.
 */
int
f4_trace_Estimate(F4_HMM *hmm, ESL_MSA *msa, F4_TRACE **tr, const F4_PRIOR *pri, double **letter_probs, double *background_probs)
{
  int M = hmm->M;     // profile size
  int N = tr[0]->L;   // sequence length (all are the same anyways)
  float wt;           // weight for the sequences
  int idx;            // index for sequences in the MSA

  double **W_bar = NULL, **Y_bar = NULL, **Z_bar = NULL; // backward
  double **X = NULL, **Y = NULL, **Z = NULL;             // forward
  F4_HMM *param_counts_new = NULL;                       // current parameter counts for the HMM
  F4_HMM *param_counts_old = NULL;                       // previous parameter counts for the HMM
  double v;                                              // variable to accumulate the sum of weights

  ESL_DSQ *dsq;     // sequence without gaps
  ESL_DSQ **dsqs;   // array of sequences without gaps
  int seq_length;   // length of sequence without gaps
  int *seq_lengths; // lengths of sequences without gaps

  int status = eslEINVAL; // status code for error handling

  if ((status = get_sequence_without_gaps(msa, &dsqs, &seq_lengths)) != eslOK)
    goto ERROR;

  if ((status = bw_build(M, N, &W_bar, &Y_bar, &Z_bar, &X, &Y, &Z)) != eslOK)
    goto ERROR;

  if ((param_counts_new = f4_hmm_Create(M, msa->abc)) == NULL) {
    status = eslEMEM;
    goto ERROR;
  }

  if ((param_counts_old = f4_hmm_Create(M, msa->abc)) == NULL) {
    status = eslEMEM;
    goto ERROR;
  }

  param_counts_save_from_to(hmm, param_counts_old);

  if ((status = estimate_parameters(param_counts_old, pri)) != eslOK) 
    goto ERROR;

  /* Needed to determine termination. */
  int termination_condition;
  int num_iterations = 0;

  do {

    /* We aggregate the expected counts over all sequences. */
    f4_hmm_Zero(param_counts_new);

    for (idx = 0; idx < msa->nseq; idx++) {

      bw_zero(M, N, W_bar, Y_bar, Z_bar, X, Y, Z);

      wt = msa->wgt[idx];

      dsq = dsqs[idx]; // get the sequence without gaps
      seq_length = seq_lengths[idx];

      /* Forward pass, calculate X, Y, Z and the aggregated v (i.e. sum of w-values). */
      v = 0.0;
      if ((status = f4_fwd(param_counts_old, dsq, tr[idx], msa->abc, X, Y, Z, M, seq_length, letter_probs, background_probs, &v)) != eslOK)
        goto ERROR;
      
      /* Backward pass, calculate W_bar, Y_bar, Z_bar. */
      if ((status = f4_bwd(param_counts_old, dsq, tr[idx], msa->abc, W_bar, Y_bar, Z_bar, M, seq_length, letter_probs, background_probs)) != eslOK)
        goto ERROR;

      /* Calculate and update parameters in new parameter counts. */
      if ((status = f4_calculate_parameters(param_counts_old, seq_length, wt, W_bar, Y_bar, Z_bar, X, Y, Z, v, letter_probs, param_counts_new)) != eslOK)
        goto ERROR;
    }

    /* Save the new estimated parameter counts into the HMM. */
    param_counts_save_from_to(param_counts_new, hmm);

    /* Turn the parameter counts into probabilities via priors. */
    if ((status = estimate_parameters(param_counts_new, pri)) != eslOK) 
      goto ERROR;

    /* Determine termination condition. Then overwrite the old probabilities with the new ones. */
    termination_condition = determine_termination_condition(param_counts_old, param_counts_new);
    param_counts_save_from_to(param_counts_new, param_counts_old);
    num_iterations++;

    printf("\rIteration %d / %d: %s", num_iterations, f4_BW_MAXITER, termination_condition ? "    converged" : "not converged");
    fflush(stdout);

  } while (!termination_condition && num_iterations < f4_BW_MAXITER);

  printf("\nError tolerance: %g, iterations: %d, %s\n", f4_BW_CONVERGE, num_iterations, termination_condition ? "Baum-Welch converged" : "Baum-Welch did not converge");

  bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);
  f4_hmm_Destroy(param_counts_new);
  f4_hmm_Destroy(param_counts_old);
  free_sequence_without_gaps(dsqs, seq_lengths, msa->nseq);
  
  return eslOK;

  ERROR:
  if (W_bar || Y_bar || Z_bar || X || Y || Z)
    bw_destroy(M, W_bar, Y_bar, Z_bar, X, Y, Z);

  if (param_counts_new)
    f4_hmm_Destroy(param_counts_new);

  if (param_counts_old)
    f4_hmm_Destroy(param_counts_old);

  if (dsqs || seq_lengths)
    free_sequence_without_gaps(dsqs, seq_lengths, msa->nseq);

  return status;
}