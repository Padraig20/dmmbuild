/* The Fig-4 core HMM data structure.
 * 
 * Contents:
 *   1. The F4_HMM object: allocation, initialization, destruction.
 *   2. Convenience routines for setting fields in an HMM.
 *   3. Renormalization and rescaling counts in core HMMs.
 *   4. Debugging and development code.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "dummer.h"

/*****************************************************************
 * 1. The F4_HMM object: allocation, initialization, destruction.
 *****************************************************************/

/* Function:  f4_hmm_Create()
 * Synopsis:  Allocate a new <F4_HMM>.
 *
 * Purpose:   Allocate a <F4_HMM> of <M> nodes, for symbol
 *            alphabet <abc>, and return a pointer to it.
 *            
 *            The HMM only keeps a copy of the <abc> alphabet
 *            pointer. The caller is responsible for providing the
 *            alphabet, keeping it around while the HMM is in use,
 *            and (eventually) free'ing the alphabet when it's
 *            not needed any more. (Basically, just a step removed
 *            from keeping the alphabet as a global.)
 *
 * Throws:    <NULL> on allocation failure.
 */
F4_HMM *
f4_hmm_Create(int M, const ESL_ALPHABET *abc) 
{
  F4_HMM *hmm = NULL;

  if ((hmm = f4_hmm_CreateShell()) == NULL) return NULL;
  f4_hmm_CreateBody(hmm, M, abc);
  return hmm;
}  


/* Function:  f4_hmm_CreateShell()
 * Synopsis:  Allocate the ``shell'' of a <F4_HMM>.
 *
 * Purpose:   Allocate the shell of a <F4_HMM>: everything that
 *            doesn't depend on knowing the number of nodes M. 
 *            
 *            HMM input (<hmmio.c>) uses two-step shell/body
 *            allocation because it has to read for a ways from the
 *            HMM file before it reads the model size M or the
 *            alphabet type.
 *
 * Returns:   a pointer to the new <F4_HMM> on success.
 *
 * Throws:    <NULL> on allocation failure.
 */
F4_HMM *
f4_hmm_CreateShell(void) 
{
  F4_HMM *hmm = NULL;
  int     z;
  int     status;

  ESL_ALLOC(hmm, sizeof(F4_HMM));
  hmm->M          = 0;
  hmm->t          = NULL;
  hmm->tp         = NULL;
  hmm->mat        = NULL;
  hmm->ins        = NULL;

  hmm->name       = NULL;
  hmm->acc        = NULL;
  hmm->desc       = NULL;
  hmm->rf         = NULL;
  hmm->mm         = NULL;
  hmm->consensus  = NULL;
  hmm->cs         = NULL;
  hmm->ca         = NULL;
  hmm->comlog     = NULL; 
  hmm->nseq       = -1;
  hmm->eff_nseq   = -1.0;
  hmm->max_length = -1;
  hmm->ctime      = NULL;
  hmm->map        = NULL;
  hmm->checksum   = 0;

  for (z = 0; z < f4_NCUTOFFS; z++) hmm->cutoff[z]  = f4_CUTOFF_UNSET;
  for (z = 0; z < f4_NEVPARAM; z++) hmm->evparam[z] = f4_EVPARAM_UNSET;
  for (z = 0; z < f4_MAXABET;  z++) hmm->compo[z]   = f4_COMPO_UNSET;

  hmm->offset   = 0;
  hmm->flags    = 0;
  hmm->abc      = NULL;
  return hmm;

 ERROR:
  return NULL;
}  

/* Function:  f4_hmm_CreateBody()
 * Synopsis:  Allocate the "body" of a <F4_HMM>.
 *
 * Purpose:   Given an allocated shell <hmm>, and a now-known number
 *            of nodes <M> and alphabet <abc>, allocate
 *            the remainder of it for that many nodes.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure; in this case, the HMM
 *            is likely corrupted, and the caller should destroy it.
 */
int
f4_hmm_CreateBody(F4_HMM *hmm, int M, const ESL_ALPHABET *abc) 
{
  int k;
  int status;

  hmm->abc = abc;
  hmm->M   = M;

  /* level 1 */
  ESL_ALLOC(hmm->t,    (M+1) * sizeof(float *));
  ESL_ALLOC(hmm->tp,   (M+1) * sizeof(float *));
  ESL_ALLOC(hmm->mat,  (M+1) * sizeof(float *));
  ESL_ALLOC(hmm->ins,  (M+1) * sizeof(float *));
  hmm->t[0]   = NULL;
  hmm->mat[0] = NULL;
  hmm->ins[0] = NULL;

  /* level 2 */
  ESL_ALLOC(hmm->t[0],   (f4H_NTRANSITIONS*(M+1)) * sizeof(float));
  ESL_ALLOC(hmm->tp[0],  (f4H_NPARAMS*(M+1))      * sizeof(float));
  ESL_ALLOC(hmm->mat[0], (abc->K*(M+1))           * sizeof(float));
  ESL_ALLOC(hmm->ins[0], (abc->K*(M+1))           * sizeof(float));
  for (k = 1; k <= M; k++) {
    hmm->mat[k] = hmm->mat[0] + k * hmm->abc->K;
    hmm->ins[k] = hmm->ins[0] + k * hmm->abc->K;
    hmm->t[k]   = hmm->t[0]   + k * f4H_NTRANSITIONS;
    hmm->tp[k]  = hmm->tp[0]  + k * f4H_NPARAMS;
  }

  /* Enforce conventions on unused but allocated distributions, so
   * Compare() tests succeed unless memory was corrupted.
   */
  if ((status = f4_hmm_Zero(hmm)) != eslOK) goto ERROR;
  hmm->mat[0][0]    = 1.0;
  hmm->t[0][f4H_DM] = 1.0; 

  /* Optional allocation, status flag dependent */
  if (hmm->flags & f4H_RF)    ESL_ALLOC(hmm->rf,         (M+2) * sizeof(char));
  if (hmm->flags & f4H_MMASK) ESL_ALLOC(hmm->mm,         (M+2) * sizeof(char));
  if (hmm->flags & f4H_CONS)  ESL_ALLOC(hmm->consensus,  (M+2) * sizeof(char));
  if (hmm->flags & f4H_CS)    ESL_ALLOC(hmm->cs,         (M+2) * sizeof(char));
  if (hmm->flags & f4H_CA)    ESL_ALLOC(hmm->ca,         (M+2) * sizeof(char));
  if (hmm->flags & f4H_MAP)   ESL_ALLOC(hmm->map,        (M+1) * sizeof(int));
  
  return eslOK;

 ERROR:
  return status;
}  


/* Function:  f4_hmm_Destroy()
 * Synopsis:  Free a <F4_HMM>.
 *
 * Purpose:   Frees both the shell and body of an <hmm>.
 *            Works even if the <hmm> is damaged (incompletely allocated)
 *            or even <NULL>.
 *
 * Note:      Remember, leave reference pointers like abc, gm, and
 *            bg alone. These are under the application's control not ours.
 *
 * Returns:   (void).
 */
void
f4_hmm_Destroy(F4_HMM *hmm)
{
  if (hmm == NULL) return;

  if (hmm->mat) {  if (hmm->mat[0]) free(hmm->mat[0]); free(hmm->mat); }
  if (hmm->ins) {  if (hmm->ins[0]) free(hmm->ins[0]); free(hmm->ins); }
  if (hmm->t)   {  if (hmm->t[0])   free(hmm->t[0]);   free(hmm->t);   }
  if (hmm->tp)  {  if (hmm->tp[0])  free(hmm->tp[0]);  free(hmm->tp);  }

  if (hmm->name)      free(hmm->name);
  if (hmm->acc)       free(hmm->acc);
  if (hmm->desc)      free(hmm->desc);
  if (hmm->rf)        free(hmm->rf);
  if (hmm->mm)        free(hmm->mm);
  if (hmm->consensus) free(hmm->consensus);
  if (hmm->cs)        free(hmm->cs);
  if (hmm->ca)        free(hmm->ca);
  if (hmm->comlog)    free(hmm->comlog);
  if (hmm->ctime)     free(hmm->ctime);
  if (hmm->map)       free(hmm->map);

  free(hmm);
  return;
}

/* Function:  f4_hmm_CopyParameters()
 * Synopsis:  Copy parameters from one HMM to another.
 *
 * Purpose:   Copy parameters of <src> to <dest>. The HMM <dest> must
 *            be allocated by the caller for the same 
 *            alphabet and M as <src>. 
 *            
 *            Both core and search model parameters are copied.
 *
 *            No annotation is copied.  This is because several
 *            annotation fields are variable-length strings that
 *            require individual allocations.  The
 *            <f4_hmm_CopyParameters()> function is for cases where we
 *            have to repeatedly reset the parameters of a model - for
 *            example, in entropy weighting.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_hmm_CopyParameters(const F4_HMM *src, F4_HMM *dest)
{
  int k;
  for (k = 0; k <= src->M; k++) {
    esl_vec_FCopy(src->t[k],   f4H_NTRANSITIONS, dest->t[k]);
    esl_vec_FCopy(src->tp[k],  f4H_NPARAMS,      dest->tp[k]);
    esl_vec_FCopy(src->mat[k], src->abc->K,      dest->mat[k]);
    esl_vec_FCopy(src->ins[k], src->abc->K,      dest->ins[k]);
  }
  return eslOK;
}

/* Function:  f4_hmm_Clone()
 * Synopsis:  Make an exact duplicate of an HMM.
 *
 * Purpose:   Duplicates an hmm.
 * 
 *            Note: does not duplicate the objects the HMM refers to,
 *            if any (profile, null model, or alphabet); only copies
 *            the reference pointers.
 * 
 * Returns:   a pointer to the duplicate.
 *
 * Throws:    <NULL> on allocation failure.
 */
F4_HMM *
f4_hmm_Clone(const F4_HMM *hmm)
{
  int     status;
  F4_HMM *new = NULL;
  int     z;

  if ((new = f4_hmm_Create(hmm->M, hmm->abc)) == NULL) goto ERROR;
  f4_hmm_CopyParameters(hmm, new);
  
  if ((status = esl_strdup(hmm->name,   -1, &(new->name)))   != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->acc,    -1, &(new->acc)))    != eslOK) goto ERROR;
  if ((status = esl_strdup(hmm->desc,   -1, &(new->desc)))   != eslOK) goto ERROR;

  if ((hmm->flags & f4H_RF)    && (status = esl_strdup(hmm->rf,        -1, &(new->rf)))        != eslOK) goto ERROR;
  if ((hmm->flags & f4H_MMASK) && (status = esl_strdup(hmm->mm,        -1, &(new->mm)))        != eslOK) goto ERROR;
  if ((hmm->flags & f4H_CONS)  && (status = esl_strdup(hmm->consensus, -1, &(new->consensus))) != eslOK) goto ERROR;
  if ((hmm->flags & f4H_CS)    && (status = esl_strdup(hmm->cs,        -1, &(new->cs)))        != eslOK) goto ERROR;
  if ((hmm->flags & f4H_CA)    && (status = esl_strdup(hmm->ca,        -1, &(new->ca)))        != eslOK) goto ERROR;
  if ((hmm->comlog != NULL)    && (status = esl_strdup(hmm->comlog,    -1, &(new->comlog)))    != eslOK) goto ERROR;
  if ((hmm->ctime  != NULL)    && (status = esl_strdup(hmm->ctime,     -1, &(new->ctime)))     != eslOK) goto ERROR;
  if (hmm->flags & f4H_MAP) {
    ESL_ALLOC(new->map, sizeof(int) * (hmm->M+1));
    esl_vec_ICopy(hmm->map, hmm->M+1, new->map);
  }
  new->nseq       = hmm->nseq;
  new->eff_nseq   = hmm->eff_nseq;
  new->max_length = hmm->max_length;
  new->checksum   = hmm->checksum;

  for (z = 0; z < f4_NEVPARAM; z++) new->evparam[z] = hmm->evparam[z];
  for (z = 0; z < f4_NCUTOFFS; z++) new->cutoff[z]  = hmm->cutoff[z];
  for (z = 0; z < f4_MAXABET;  z++) new->compo[z]   = hmm->compo[z];

  new->offset   = hmm->offset;
  new->flags    = hmm->flags;
  new->abc      = hmm->abc;
  return new;

 ERROR:
  if (new != NULL) f4_hmm_Destroy(new);
  return NULL;
}


/* Function:  f4_hmm_Zero()
 * Synopsis:  Set all parameters to zero (including model composition).
 *
 * Purpose:   Zeroes all counts/probabilities fields in core model,
 *            including emissions, transitions, and model
 *            composition.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_hmm_Zero(F4_HMM *hmm)
{
  int k;

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FSet(hmm->t[k],   f4H_NTRANSITIONS, 0.);  
    esl_vec_FSet(hmm->tp[k],  f4H_NPARAMS,      0.);
    esl_vec_FSet(hmm->mat[k], hmm->abc->K, 0.);  
    esl_vec_FSet(hmm->ins[k], hmm->abc->K, 0.);  
  }
  esl_vec_FSet(hmm->compo, f4_MAXABET, 0.);
  return eslOK;
}

/*****************************************************************
 * 2. Convenience routines for setting fields in an HMM.
 *****************************************************************/ 

/* Function: f4_hmm_SetName()
 * Synopsis: Set or change the name of an HMM.
 * 
 * Purpose:  Set or change the name of a Plan7 HMM to <name>.
 *           Any trailing whitespace (including newline) is chopped off.     
 *      
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
f4_hmm_SetName(F4_HMM *hmm, char *name)
{
  int   status;
  void *tmp;
  int   n;

  if (name == NULL) {
    if (hmm->name != NULL) free(hmm->name); 
    hmm->name = NULL;
  } else {
    n = strlen(name);
    ESL_RALLOC(hmm->name, tmp, sizeof(char)*(n+1));
    strcpy(hmm->name, name);
    if ((status = esl_strchop(hmm->name, n)) != eslOK) goto ERROR;
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: f4_hmm_SetAccession()
 * Synopsis: Set or change the accession of an HMM.
 * 
 * Purpose:  Set or change the accession number of a fig-4 HMM to <acc>,
 *           and raise the <F4_ACC> flag. Trailing whitespace (including newline) 
 *           is chopped.  
 *           
 *           If <acc> is <NULL>, unset the HMM's accession (if any) and drop 
 *           the <F4_ACC> flag.
 *
 * Returns:  <eslOK> on success.
 *
 * Throws:   <eslEMEM> on allocation error, and original name (if any) 
 *           remains.
 */
int
f4_hmm_SetAccession(F4_HMM *hmm, char *acc)
{
  int   status;
  void *tmp;
  int   n;

  if (acc == NULL) {
    if (hmm->acc != NULL) free(hmm->acc); 
    hmm->acc = NULL;
    hmm->flags &= ~f4H_ACC;	/* legacy */
  } else {
    n = strlen(acc);
    ESL_RALLOC(hmm->acc, tmp, sizeof(char)*(n+1));
    strcpy(hmm->acc, acc);
    if ((status = esl_strchop(hmm->acc, n)) != eslOK) goto ERROR;
    hmm->flags |= f4H_ACC;	/* legacy */
  }
  return eslOK;

 ERROR:
  return status;
}

/* Function: f4_hmm_SetDescription()
 * Synopsis: Set or change the description line of an HMM.
 * 
 * Purpose:  Set or change the description line of a fig4 HMM. 
 *           Trailing whitespace (including newline) is chopped.
 */
int
f4_hmm_SetDescription(F4_HMM *hmm, char *desc)
{
  int   status;
  void *tmp;
  int   n;

  if (desc == NULL) 
    {
      if (hmm->desc != NULL) free(hmm->desc); 
      hmm->desc   = NULL;
      hmm->flags &= ~f4H_DESC;	/* legacy */
    }
  else
    {
      n = strlen(desc);
      ESL_RALLOC(hmm->desc, tmp, sizeof(char)*(n+1));
      strcpy(hmm->desc, desc);
      if ((status = esl_strchop(hmm->desc, n)) != eslOK) goto ERROR;
      hmm->flags |= f4H_DESC;	/* legacy */
    }
  return eslOK;

 ERROR:
  return status;
}

/* Function: f4_hmm_AppendComlog()
 * Synopsis: Concatenate and append command line to the command line log.
 * 
 * Purpose:  Concatenate command line options and append as a new line in the
 *           command line log. Command line log is multiline, with each line
 *           ending in newline char, except for last line.
 *           
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure.          
 */
int
f4_hmm_AppendComlog(F4_HMM *hmm, int argc, char **argv)
{
  int   status;
  void *tmp;
  int   n;
  int   i;

  /* figure out length of added command line, and (re)allocate comlog */
  n = argc-1;	/* account for 1 space per arg, except last one */
  for (i = 0; i < argc; i++)
    n += strlen(argv[i]);

  if (hmm->comlog != NULL) {
    n += strlen(hmm->comlog) + 1; /* +1 for the \n we're going to add to the old comlog */
    ESL_RALLOC(hmm->comlog, tmp, sizeof(char)* (n+1));
    strcat(hmm->comlog, "\n");
  } else {
    ESL_ALLOC(hmm->comlog, sizeof(char)* (n+1));
    *(hmm->comlog) = '\0'; /* need this to make strcat work */
  }

  for (i = 0; i < argc-1; i++)
    {
      strcat(hmm->comlog, argv[i]);
      strcat(hmm->comlog, " ");
    }
  strcat(hmm->comlog, argv[argc-1]);
  return eslOK;

 ERROR:
  return status;
}

/* Function: f4_hmm_SetCtime()
 * Synopsis: Timestamp an HMM.
 * 
 * Purpose:  Set the <ctime> field in a new HMM to the current time.
 * 
 * Returns:  <eslOK> on success.
 * 
 * Throws:   <eslEMEM> on allocation failure. 
 *           <eslESYS> if system calls fail to obtain (or format) the time.
 *           
 * Notes:    This function calls <ctime_r()>, supposedly a part of the
 *           ISO/IEC 9945-1:1996 (POSIX.1) standard, but not ANSI
 *           C99, so we have potential portability problems here.
 *
 *           A known one: <ctime_r()> is by default a three-argument
 *           call on Solaris 10 systems. Our autoconf script sets
 *           -D_POSIX_PTHREAD_SEMANTICS on Solaris systems to fix this
 *           issue, requesting Solaris to use a compliant version of 
 *           ctime_r().
 *
 *           We might want to use strftime() instead; that's what 
 *           POSIX 2008 recommends; but we'd still need localtime_r() or
 *           its equivalent, and that has its own portability issues.
 *           
 *           Note to porters: it really doesn't matter what this
 *           timestamp is. HMMER doesn't look at it, it's for human
 *           notetaking. If you have to, set it to an empty string.
 *
 * TODO:     Oi. Time is complicated. Easel should give us an
 *           easy and portable call to generate time stamps like this;
 *           an esl_time module, perhaps?
 */
int
f4_hmm_SetCtime(F4_HMM *hmm)
{
  char    *s = NULL;
  time_t   date;
  int      status;

  ESL_ALLOC(s, 32);
  if ((date = time(NULL)) == -1)               { status = eslESYS; goto ERROR; }
  if (ctime_r(&date, s) == NULL)               { status = eslESYS; goto ERROR; }
  if ((status = esl_strchop(s, -1)) != eslOK)  {                   goto ERROR; }
  
  if (hmm->ctime != NULL) free(hmm->ctime);
  hmm->ctime = s;
  return eslOK;

 ERROR:
  if (s) free(s);
  return status;
}

/* Function:  f4_hmm_SetComposition()
 * Synopsis:  Calculate and set model composition, <hmm->compo[]>
 *
 * Purpose:   Calculates the mean residue composition emitted by
 *            model <hmm>, and set <hmm->compo[]> to it.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEMEM> on allocation failure, in which case 
 *            values in <hmm->compo[]> are unchanged.
 *            
 * Note:      In principle, you should be able to normalize
 *            hmm->compo[] by dividing thru by the sum of
 *            mocc[] and iocc[] vectors, and that's what the
 *            3.0 release version did. This allowed f4_hmm_Validate()
 *            to check compo[] for summation to 1.0, as a smoke
 *            check for bugs here. The problem with that 
 *            is numerical roundoff error accumulation, when
 *            hmm->M is large [bug #h84]. To fix #h84, we
 *            simply renormalize compo[], rather than the fancier
 *            previous version. This avoids error accumulation,
 *            but it also guarantees that compo[] will trivially
 *            pass the hmm_Validation() step; it's not really
 *            validating the SetComposition() calculation at all.                                 
 *            (For description of #h84, error analysis, and the fix,
 *            xref J7/7; SRE, Tue Nov  2 14:32:29 2010)
 */
int
f4_hmm_SetComposition(F4_HMM *hmm)
{
  float *mocc = NULL;
  float *iocc = NULL;
  int    k;
  int    status;

  ESL_ALLOC(mocc, sizeof(float) * (hmm->M+1));
  ESL_ALLOC(iocc, sizeof(float) * (hmm->M+1));

  f4_hmm_CalculateOccupancy(hmm, mocc, iocc);
  esl_vec_FSet(hmm->compo, hmm->abc->K, 0.0);
  esl_vec_FAddScaled(hmm->compo, hmm->ins[0], iocc[0], hmm->abc->K);
  for (k = 1; k <= hmm->M; k++)
    {
      esl_vec_FAddScaled(hmm->compo, hmm->mat[k], mocc[k], hmm->abc->K);
      esl_vec_FAddScaled(hmm->compo, hmm->ins[k], iocc[k], hmm->abc->K);
    }

  esl_vec_FNorm(hmm->compo, hmm->abc->K);
  hmm->flags  |= f4H_COMPO;

  free(mocc);
  free(iocc);
  return eslOK;

 ERROR:
  if (mocc != NULL) free(mocc);
  if (iocc != NULL) free(iocc);
  return status;
}

/* Function:  f4_hmm_SetConsensus()
 * Synopsis:  Set the consensus residue line of the HMM.
 *
 * Purpose:   Sets the consensus annotation line of the model <hmm>.
 *            
 *            Behavior differs, depending on whether this is a
 *            single-sequence model (i.e. phmmer) or a standard
 *            model of a multiple sequence alignment. If <sq> is
 *            non-<NULL> this is a single-sequence model and <sq> is
 *            the digital sequence it was built from. If <sq> is <NULL>
 *            this is a standard multiple-sequence model.
 *            
 *            In a standard model, the most likely (highest emission
 *            probability) residue is the consensus at each position.
 *            In a single-sequence model, the consensus is the
 *            sequence itself.
 *            
 *            In both cases, if the emission probability is $\geq$
 *            certain threshold, the residue is upper cased. The
 *            threshold is arbitrarily set to 0.9 for nucleic acid
 *            alphabets (<eslDNA>, <eslRNA>) and 0.5 for amino acid
 *            alphabets (<eslAMINO>) and all other alphabets.
 *            
 *            The special handling of single-sequence models avoids
 *            a counterintuitive case where the most likely residue is
 *            not the original residue. For example, under the
 *            BLOSUM62 matrix, given an observed M, the most likely
 *            aligned residue is an L, not an M. (Because L is so much
 *            more likely a priori than M.)
 *
 * Args:      hmm    - model with valid probability parameters mat[1..M][x]
 *            sq     - NULL if a standard model;
 *                     or the query sequence for a single-sequence model.
 *           
 * Returns:   <eslOK> on success. The <f4H_CONS> flag on the <hmm> is raised
 *            if it wasn't already. The <hmm->consensus> line is set.
 *
 * Throws:    <eslEMEM> on allocation error. The <f4H_CONS> is dropped, even
 *            if it was up to begin with, and the <hmm->consensus> is <NULL>,
 *            even if we had one to begin with.
 */
int
f4_hmm_SetConsensus(F4_HMM *hmm, ESL_SQ *sq)
{
  int   k, x;
  float mthresh;
  int   status;
  
  /* allocation, if needed */
  if (! hmm->consensus) ESL_ALLOC(hmm->consensus, sizeof(char) * (hmm->M+2));

  /* set our arbitrary threshold for upper/lower casing */
  if      (hmm->abc->type == eslAMINO) mthresh = 0.5;
  else if (hmm->abc->type == eslDNA)   mthresh = 0.9;
  else if (hmm->abc->type == eslRNA)   mthresh = 0.9;
  else                                 mthresh = 0.5;

  hmm->consensus[0] = ' ';
  for (k = 1; k <= hmm->M; k++) 
    {
      x = (sq ?  sq->dsq[k] : esl_vec_FArgMax(hmm->mat[k], hmm->abc->K));
      hmm->consensus[k] = ((hmm->mat[k][x] >= mthresh) ? toupper(hmm->abc->sym[x]) : tolower(hmm->abc->sym[x]));
    }
  hmm->consensus[hmm->M+1] = '\0';
  hmm->flags  |= f4H_CONS;	
  return eslOK;

 ERROR:
  if (hmm->consensus) free(hmm->consensus);
  hmm->consensus = NULL;
  hmm->flags    &= (~f4H_CONS);	
  return status;
}
/*---------------- end, internal-setting routines ---------------*/

/*****************************************************************
 * 3. Renormalization and rescaling counts in core HMMs.
 *****************************************************************/ 

/* Function:  f4_hmm_Scale()
 * Synopsis:  In a model containing counts, rescale counts by a factor.
 *
 * Purpose:   Given a counts-based model <hmm>, scale core
 *            by a multiplicative factor of <scale>, where <scale> is
 *            often <eff_nseq/nseq> for absolute sequence weighting.
 *            Only affects core probability model emissions and 
 *            transitions (<t>, <ins>, and <mat>).
 *
 * Args:      hmm        - counts based HMM.
 *            scale      - scaling factor (e.g. eff_nseq/nseq); 1.0=no scaling.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_hmm_Scale(F4_HMM *hmm, double scale)
{
  int k;

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FScale(hmm->t[k],   f4H_NTRANSITIONS, scale); 
    esl_vec_FScale(hmm->tp[k],  f4H_NPARAMS,      scale); 
    esl_vec_FScale(hmm->mat[k], hmm->abc->K,      scale);  
    esl_vec_FScale(hmm->ins[k], hmm->abc->K,      scale);  
  }
  return eslOK;
}


/* Function:  f4_hmm_ScaleExponential()
 * Synopsis:  In a model containing counts, rescale counts by an exponential factor.
 *
 * Purpose:   Given a counts-based model <hmm>, scale core by an
 *            exponential factor <exp>. This should be thought of as
 *            an alternative to f4_hmm_Scale(). Let C_i be the total
 *            observed count in column i, and F be the scale. In
 *            f4_hmm_Scale, the updated total observed count would be
 *            C_i = C_i * F  (i.e. the scaling factor is uniform across
 *            all columns). In this function, C_i = C_i ^ F. The result
 *            is a non-uniform scaling across columns -- columns with
 *            higher C_i will be reduced to a greater extent than will
 *            columns with low counts.
 *
 *            Consider the case where one column has 30 observations and a
 *            bunch of others have 300. This can happen when heavily-
 *            fragmented sequences are used to reconstruct a family MSA, as
 *            in Dfam models ... but isn't likely to have been seen in Pfam
 *            alignments. Though the column with 30 observations isn't nearly
 *            as complete as the one with 300, it still has enough that we
 *            shouldn't be willing to discount the observations entirely -
 *            something that might happen if uniform entropy weighting needs
 *            to push the average observations down 10-fold in order to achieve
 *            the desired avg relative entropy.
 *            e.g.
 *
 *             C_i    F  ->  C_i
 *               3   .8      2.4
 *              30   .8       15
 *             300   .8       96
 *               3   .7      2.2
 *              30   .7       11
 *             300   .7       54
 *               3   .6      1.9
 *              30   .6      7.6
 *             300   .6       30
 *
 *            Note: the observed counts will never drop below 1 in this case.
 *
 *            After computing the per-column total scale for column i, that
 *            scale is applied to the core probability model emissions and
 *            transitions (<t>, <ins>, and <mat>) for position i.
 *
 * Args:      hmm     - counts based HMM.
 *            exp     - exponential factor; 1.0=no scaling.
 *            ret_scaleavg - returns the mean of the per-column scale factors corresponding
 *                           to the factor exp.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_hmm_ScaleExponential(F4_HMM *hmm, double exp)
{

  int k;
  for (k = 1; k <= hmm->M; k++) {

    float count = esl_vec_FSum(hmm->mat[k], hmm->abc->K);
    float new_count = pow(count, exp);
    double scale = count>0 ? new_count / count : 1.0;  /* if no counts in the column (strange, but possible), just use default freqs*/

    esl_vec_FScale(hmm->t[k],   f4H_NTRANSITIONS, scale);
    esl_vec_FScale(hmm->tp[k],  f4H_NPARAMS,      scale);
    esl_vec_FScale(hmm->mat[k], hmm->abc->K,      scale);
    esl_vec_FScale(hmm->ins[k], hmm->abc->K,      scale);
  }
  return eslOK;
}

/* Function: f4_hmm_Renormalize()
 * Synopsis: Renormalize all parameter vectors (emissions/transitions).
 * 
 * Purpose:  Take a core HMM in counts form, and renormalize
 *           all probability vectors in the core probability model. Enforces
 *           fig4 restrictions on nonexistent transitions.
 *
 *           Leaves other flags (stats and profile) alone, so caller
 *           needs to be wary. Renormalizing a probability model that
 *           has stats and profile scores wouldn't usually invalidate
 *           those data; and if we're renormalizing a counts model, we
 *           shouldn't have stats or profile scores yet anyway.
 *           
 * Args:     hmm - the model to renormalize.
 *                 
 * Return:   <eslOK> on success.
 */                          
int
f4_hmm_Renormalize(F4_HMM *hmm)
{
  int   k;			/* counter for model position */

  for (k = 0; k <= hmm->M; k++) {
    esl_vec_FNorm(hmm->mat[k], hmm->abc->K);
    esl_vec_FNorm(hmm->ins[k], hmm->abc->K);
    esl_vec_FNorm(F4H_TMAT(hmm, k), f4H_NTMAT);	/* TMX */
    esl_vec_FNorm(F4H_TDEL(hmm, k), f4H_NTDEL);	/* TIX */
    esl_vec_FNorm(F4H_TINS(hmm, k), f4H_NTINS);	/* TDX */
  }
  /* If t[M][TD*] distribution was all zeros, we just made TDD nonzero. Oops.
   * Re-enforce t's on that final delete state. */
  hmm->t[hmm->M][f4H_DM] = 1.0;
  hmm->t[hmm->M][f4H_DD] = 0.0;

  /* Rare: if t[M][TM*] distribution was all zeros (all final transitions
   * were D_M -> E) then we just made nonexistent M_M->D_M+1 transition nonzero.
   * Fix that too.
   */
  if (hmm->t[hmm->M][f4H_MD] > 0.) {
    hmm->t[hmm->M][f4H_MD] = 0.;
    hmm->t[hmm->M][f4H_MM] = 0.5;
    hmm->t[hmm->M][f4H_MI] = 0.5;
  }

  return eslOK;
}

/*****************************************************************
 * 4. Other routines in the API.
 *****************************************************************/

/* Function:  f4_hmm_CalculateOccupancy()
 * Synopsis:  Calculate match occupancy and insert expected use count vectors.
 *
 * Purpose:   Calculate a vector <mocc[1..M]> containing probability
 *            that each match state is used in a sampled glocal path through
 *            the model. Caller provides allocated space (<M+1> floats)
 *            for <mocc>.
 *            
 *            Caller may optionally provide an array <iocc[0..M]> as
 *            well, which (if provided) will be set to contain the
 *            expected number of times that a sampled path would contain
 *            each insert state.
 *
 * Returns:   <eslOK> on success.
 */
int
f4_hmm_CalculateOccupancy(const F4_HMM *hmm, float *mocc, float *iocc)
{
  int k;

  mocc[0] = 0.;			                     /* no M_0 state */
  mocc[1] = hmm->t[0][f4H_MI] + hmm->t[0][f4H_MM];   /* initialize w/ 1 - B->D_1 */
  for (k = 2; k <= hmm->M; k++)
      mocc[k] = mocc[k-1] * (hmm->t[k-1][f4H_MM] + hmm->t[k-1][f4H_MI]) +
        (1.0-mocc[k-1]) * hmm->t[k-1][f4H_DM];
  if (iocc != NULL) {
    iocc[0] = hmm->t[0][f4H_MI] / hmm->t[0][f4H_IM];
    for (k = 1; k <= hmm->M; k++)
      iocc[k] = mocc[k] * hmm->t[k][f4H_MI] / hmm->t[k][f4H_IM];
  }

  return eslOK;
}

/*---------------- end of the rest of the API -------------------*/
