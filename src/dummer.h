#include <stdio.h>		
#include <stddef.h>             // ptrdiff_t 
#include <string.h>
#include <math.h>
#include <float.h>
#include <ctype.h>

#include "easel.h"
#include "esl_alphabet.h"	/* ESL_DSQ, ESL_ALPHABET */
#include "esl_dmatrix.h"	/* ESL_DMATRIX           */
#include "esl_dirichlet.h"  /* ESL_DIRICHLET         */
#include "esl_getopts.h"	/* ESL_GETOPTS           */
#include "esl_gumbel.h"  /* ESL_GUMBEL            */
#include "esl_histogram.h"      /* ESL_HISTOGRAM         */
#include "esl_hmm.h"	        /* ESL_HMM               */
#include "esl_keyhash.h"        /* ESL_KEYHASH           */
#include "esl_mixdchlet.h"	/* ESL_MIXDCHLET         */
#include "esl_msa.h"		/* ESL_MSA               */
#include "esl_msafile.h"    /* ESL_MSAFILE           */
#include "esl_msacluster.h" /* ESL_MSACLUSTER        */
#include "esl_random.h"		/* ESL_RANDOMNESS        */
#include "esl_rand64.h" /* ESL_RAND64 */
#include "esl_rootfinder.h"     /* ESL_ROOTFINDER        */
#include "esl_sq.h"		/* ESL_SQ                */
#include "esl_sse.h"    /* ESL_SSE               */
#include "esl_scorematrix.h"    /* ESL_SCOREMATRIX       */
#include "esl_stopwatch.h"      /* ESL_STOPWATCH         */
#include "esl_msaweight.h"      /* ESL_MSAWEIGHT         */
#include "esl_vectorops.h"      /* ESL_VECTOROPS         */

#include <xmmintrin.h>    /* SSE  */
#include <emmintrin.h>    /* SSE2 */

/*****************************************************************
 * 0. Configuration & Constants
 *****************************************************************/

#define DUMMER_VERSION "0.1"
#define  DUMMER_DATE "Jul 2025"
#define  DUMMER_COPYRIGHT "Copyright (C) 2025 The University of Tokyo."
#define  DUMMER_LICENSE "Freely distributed under the MIT open source license."
#define  DUMMER_URL "https://gitlab.com/mcfrith/seq-position-probs"

/* Search modes. */
#define f4_NO_MODE   0
#define f4_LOCAL     1		/* multihit local:  "fs" mode   */
#define f4_GLOCAL    2		/* multihit glocal: "ls" mode   */
#define f4_UNILOCAL  3		/* unihit local: "sw" mode      */
#define f4_UNIGLOCAL 4		/* unihit glocal: "s" mode      */

#define f4_IsLocal(mode)  (mode == f4_LOCAL || mode == f4_UNILOCAL)
#define f4_IsMulti(mode)  (mode == f4_LOCAL || mode == f4_GLOCAL)

/* In calculating Q, the number of vectors we need in a row, we have
 * to make sure there's at least 2, or a striped implementation fails.
 */
#define f4O_NQB(M)   ( ESL_MAX(2, ((((M)-1) / 16) + 1)))   /* 16 uchars  */
#define f4O_NQW(M)   ( ESL_MAX(2, ((((M)-1) / 8)  + 1)))   /*  8 words   */
#define f4O_NQF(M)   ( ESL_MAX(2, ((((M)-1) / 4)  + 1)))   /*  4 floats  */

#define f4O_EXTRA_SB 17 /* see ssvfilter.c (in HMMER repo) for explanation */

/* Relative entropy target defaults:
 * For proteins, hmmbuild's effective sequence number calculation
 * aims to achieve a certain relative entropy per match emission.
 * (= average score per match emission).
 * These are empirically tuned constants,
 */
#define f4_ETARGET_AMINO  0.59 /* bits,  from the work of Steve Johnson. */
#define f4_ETARGET_DNA    0.62 /* bits,  from the work of Travis Wheeler and Robert Hubley. */
#define f4_ETARGET_OTHER  1.0  /* bits */ /* if you define your own alphabet, set this */

enum f4_evparams_e {    f4_MMU  = 0, f4_MLAMBDA = 1,     f4_VMU = 2,  f4_VLAMBDA = 3, f4_FTAU = 4, f4_FLAMBDA = 5 };
enum f4_cutoffs_e  {     f4_GA1 = 0,     f4_GA2 = 1,     f4_TC1 = 2,      f4_TC2 = 3,  f4_NC1 = 4,     f4_NC2 = 5 };
enum f4_offsets_e  { f4_MOFFSET = 0, f4_FOFFSET = 1, f4_POFFSET = 2 };

/* Option flags when creating multiple alignments with f4_tracealign_*() */
#define f4_DEFAULT             0
#define f4_DIGITIZE            (1<<0)
#define f4_ALL_CONSENSUS_COLS  (1<<1)
#define f4_TRIM                (1<<2)

/*****************************************************************
 * 1. Miscellaneous functions
 *****************************************************************/

extern void f4_banner(FILE *fp, const char *progname, char *banner);
extern void f4_Die(char *format, ...);
extern void f4_Fail(char *format, ...);
extern char f4_alidisplay_EncodePostProb(float p);
extern int f4_AminoFrequencies(float *f);

/*****************************************************************
 * 2. F4_HMM: a core fig-4 HMM model.
 *****************************************************************/

/* Bit flags used in <hmm->flags>: optional annotation in an HMM
 * 
 * Flags marked with ! may not be changed nor used for other meanings,
 * because they're codes used by HMMER2 (and earlier) that must be
 * preserved for reverse compatibility with old HMMER files.
 * 
 * Why use flags? (So I don't ask this question of myself again:)
 *   1. The way we allocate an HMM, we need to know if we're allocating
 *      M-width annotation fields (RF, CS, CA, MAP) before we read the
 *      annotation from the file.
 *   2. Historically, H2 used flags, so we still need to read H2 flags
 *      for backwards compatibility; so we may as well keep using them.
 */
#define f4H_HASBITS (1<<0)    /* obsolete (was: model has log-odds scores)       !*/
#define f4H_DESC    (1<<1)    /* description exists (legacy; xref SRE:J5/114)    !*/
#define f4H_RF      (1<<2)    /* #RF annotation available                        !*/
#define f4H_CS      (1<<3)    /* #CS annotation available                        !*/
#define f4H_XRAY    (1<<4)    /* obsolete (was: structural data available)       !*/
#define f4H_HASPROB (1<<5)    /* obsolete (was: model in probability form)       !*/
#define f4H_HASDNA  (1<<6)    /* obsolete (was: protein HMM->DNA seq params set) !*/
#define f4H_STATS   (1<<7)    /* model has E-value statistics calibrated         !*/
#define f4H_MAP     (1<<8)    /* alignment map is available                      !*/
#define f4H_ACC     (1<<9)    /* accession is available (legacy; xref SRE:J5/114)!*/
#define f4H_GA      (1<<10)   /* gathering thresholds available                  !*/
#define f4H_TC      (1<<11)   /* trusted cutoffs available                       !*/
#define f4H_NC      (1<<12)   /* noise cutoffs available                         !*/
#define f4H_CA      (1<<13)   /* surface accessibilities available               !*/
#define f4H_COMPO   (1<<14)   /* model-specific residue composition available     */
#define f4H_CHKSUM  (1<<15)   /* model has an alignment checksum                  */
#define f4H_CONS    (1<<16)   /* consensus residue line available                 */
#define f4H_MMASK   (1<<17)   /* #MM annotation available                        !*/

/* Indices of Fig4 main model state transitions, hmm->t[k][] */
enum f4h_transitions_e {
  f4H_MM = 0,
  f4H_MI = 1,
  f4H_MD = 2,
  f4H_IM = 3,
  f4H_II = 4,
  f4H_ID = 5, /* Add this for f4-HMM */
  f4H_DM = 6,
  f4H_DD = 7 ,
  f4H_DI = 8, /* Add this for f4-HMM */
};
#define f4H_NTRANSITIONS 9 /* Change this for f4-HMM */

/* Parameters used for estimating transition probabilities if fig-4 HMM */
enum f4h_params_e {
  f4H_ALPHA     = 0,
  f4H_DELTA     = 1,
  f4H_GAMMA     = 2,
  f4H_BETA      = 3,
  f4H_BETAP     = 4,
  f4H_EPSILON   = 5,
  f4H_EPSILONP  = 6
};
#define f4H_NPARAMS 7

/* How the hmm->t[k] vector is interpreted as separate probability vectors. */
#define F4H_TMAT(hmm, k) ((hmm)->t[k])
#define F4H_TINS(hmm, k) ((hmm)->t[k]+3)
#define F4H_TDEL(hmm, k) ((hmm)->t[k]+6)
#define f4H_NTMAT 3
#define f4H_NTDEL 3 /* Change this for f4-HMM */
#define f4H_NTINS 3 /* Change this for f4-HMM */

#define f4_NEVPARAM 6	/* number of statistical parameters stored in models                      */
#define f4_NCUTOFFS 6	/* number of Pfam score cutoffs stored in models                          */
#define f4_NOFFSETS 3	/* number of disk offsets stored in models for hmmscan's fast model input */

/* Option flags when creating faux traces with f4_trace_FauxFromMSA() */
#define f4_MSA_COORDS	       (1<<0) /* default: i = unaligned seq residue coords     */

/* The symbol alphabet is handled by ESL_ALPHABET objects, which
 * dynamically allocate; but sometimes HMMER uses statically-allocated
 * space, and it's useful to know a reasonable maximum for
 * symbol alphabet size.
 */
#define f4_MAXABET    20      /* maximum size of alphabet (4 or 20)              */
#define f4_MAXCODE    29      /* maximum degenerate alphabet size (18 or 29)     */

#define f4_EVPARAM_UNSET -99999.0f  /* if evparam[0] is unset, then all unset                         */
#define f4_CUTOFF_UNSET  -99999.0f  /* if cutoff[XX1] is unset, then cutoff[XX2] unset, XX={GA,TC,NC} */
#define f4_COMPO_UNSET   -1.0f      /* if compo[0] is unset, then all unset                           */

/* Some notes:
 *   0. The model might be either in counts or probability form.
 *   1. t[0] is special: t[0][TMM,TMI,TMD] are the begin->M_1,I_0,D_1 entry probabilities,
 *      t[0][TIM,TII] are the I_0 transitions, and delete state 0 doesn't
 *      exist. Therefore D[0] transitions and mat[0] emissions are unused.
 *      To simplify some normalization code, we adopt a convention that these are set
 *      to valid probability distributions: 1.0 for t[0][TDM] and mat[0][0],
 *      and 0 for the rest.
 *   2. t[M] is also special: TMD and TDD are 0 because there is no next delete state;
 *      TDM is therefore 1.0 by definition. TMM and TDM are interpreted as the
 *      M->E and D->E end transitions. t[M][TDM] must be 1.0, therefore.
 */
typedef struct f4_hmm_s {
  /*::cexcerpt::fig4_core::begin::*/
  int     M;                    /* length of the model (# nodes)                           */
  float **t;                    /* transition prob's. t[(0),1..M][0..f4H_NTRANSITIONS-1]   */
  float **mat;                  /* match emissions.  mat[1..M][0..K-1]                     */ 
  float **ins;                  /* insert emissions. ins[1..M][0..K-1]                     */
  float **tp;                   /* parameter transitions for f4-hmm. tp[0..M][0..3]        */
  /*::cexcerpt::fig4_core::end::*/

  /* Annotation. Everything but <name> is optional. Flags are set when
   * optional values are set. All the char *'s are proper nul-terminated
   * strings, not just arrays. (hmm->map is an int array).
   */
  char    *name;                 /* name of the model                     (mandatory)      */ /* String, \0-terminated   */
  char    *acc;	                 /* accession number of model (Pfam)      (optional: NULL) */ /* String, \0-terminated   */
  char    *desc;                 /* brief (1-line) description of model   (optional: NULL) */ /* String, \0-terminated   */
  char    *rf;                   /* reference line from alignment 1..M    (f4H_RF)         */ /* String; 0=' ', M+1='\0' */
  char    *mm;                   /* model mask line from alignment 1..M   (f4H_MM)         */ /* String; 0=' ', M+1='\0' */
  char    *consensus;		 /* consensus residue line        1..M    (f4H_CONS)       */ /* String; 0=' ', M+1='\0' */
  char    *cs;                   /* consensus structure line      1..M    (f4H_CS)         */ /* String; 0=' ', M+1='\0' */
  char    *ca;	                 /* consensus accessibility line  1..M    (f4H_CA)         */ /* String; 0=' ', M+1='\0' */

  char    *comlog;               /* command line(s) that built model      (optional: NULL) */ /* String, \0-terminated   */
  int      nseq;	         /* number of training sequences          (optional: -1)   */
  float    eff_nseq;             /* effective number of seqs (<= nseq)    (optional: -1)   */
  int	   max_length;           /* upper bound length, all but 1e-7 prob (optional: -1)   */
  char    *ctime;	         /* creation date                         (optional: NULL) */
  int     *map;	                 /* map of alignment cols onto model 1..M (f4H_MAP)        */ /* Array; map[0]=0 */
  uint32_t checksum;             /* checksum of training sequences        (f4H_CHKSUM)     */
  float    evparam[f4_NEVPARAM]; /* E-value params                        (f4H_STATS)      */
  float    cutoff[f4_NCUTOFFS];  /* Pfam score cutoffs                    (f4H_{GA,TC,NC}) */
  float    compo[f4_MAXABET];    /* model bg residue comp                 (f4H_COMPO)      */

  off_t    offset;               /* HMM record offset on disk                              */
  const ESL_ALPHABET *abc;       /* ptr to alphabet info (hmm->abc->K is alphabet size)    */
  int      flags;                /* status flags                                           */
} F4_HMM;

/*****************************************************************
 * 3. F4_BG: a null (background) model.
 *****************************************************************/

/* This really contains three different things: 
 *     
 *   - the "null1" model, a one-state HMM consisting of background
 *     frequencies <f> and a parameter <p1> for a target-length
 *     dependent geometric;
 *     
 *   - the "bias filter" <fhmm> a two-state HMM composed from null1's
 *     background <f> and the model's mean composition <compo>. This
 *     model is constructed dynamically, every time a new profile is 
 *     considered;
 *     
 *   - a single term <omega> that's needed by the "null2" model to set
 *     a balance between the null1 and null2 scoring terms.  The null2
 *     model is otherwise defined by construction, in f4_domaindef.c.
 *
 */
typedef struct f4_bg_s {
  float   *f;		/* null1 background residue frequencies [0..K-1]: set at initialization    */
  float    p1;		/* null1's transition prob: f4_bg_SetLength() sets this from target seq L  */

  ESL_HMM *fhmm;	/* bias filter: f4_bg_SetFilter() sets this, from model's mean composition */

  float    omega;	/* the "prior" on null2/null3: set at initialization (one omega for both null types)  */

  const ESL_ALPHABET *abc;	/* reference to alphabet in use: set at initialization             */
} F4_BG;

/*****************************************************************
 * 2. F4_PRIOR: mixture Dirichlet prior for profile HMMs
 *****************************************************************/

typedef struct f4_prior_s {
  ESL_MIXDCHLET *tm;		/*  match transitions */
  ESL_MIXDCHLET *ti;		/* insert transitions */
  ESL_MIXDCHLET *td;		/* delete transitions */
  ESL_MIXDCHLET *em;		/*  match emissions   */
  ESL_MIXDCHLET *ei;		/* insert emissions   */
  ESL_MIXDCHLET *pradg; /* parameter transitions - alpha, delta, gamma */
  ESL_MIXDCHLET *prb;   /* parameter transitions - beta, beta' */
  ESL_MIXDCHLET *pre;   /* parameter transitions - epsilon, epsilon' */
} F4_PRIOR;

/*****************************************************************
 * 3. F4_BUILDER: pipeline for new HMM construction
 *****************************************************************/

#define f4_DEFAULT_WINDOW_BETA  1e-7

enum f4_archchoice_e { f4_ARCH_FAST = 0 };
enum f4_wgtchoice_e  { f4_WGT_NONE  = 0, f4_WGT_GIVEN = 1, f4_WGT_GSC    = 2, f4_WGT_PB       = 3, f4_WGT_BLOSUM = 4 };
enum f4_effnchoice_e { f4_EFFN_NONE = 0, f4_EFFN_SET  = 1, f4_EFFN_CLUST = 2, f4_EFFN_ENTROPY = 3, f4_EFFN_ENTROPY_EXP = 4 };

typedef struct f4_builder_s {
  /* Model architecture                                                                            */
  enum f4_archchoice_e arch_strategy;    /* choice of model architecture determination algorithm   */
  float                symfrac;	         /* residue occ thresh for fast architecture determination */
  float                fragthresh;	 /* if L <= fragthresh*alen, seq is called a fragment      */

  /* Relative sequence weights                                                                     */
  enum f4_wgtchoice_e  wgt_strategy;     /* choice of relative sequence weighting algorithm        */
  double               wid;		 /* %id threshold for BLOSUM relative weighting            */

  /* Effective sequence number                                                                     */
  enum f4_effnchoice_e effn_strategy;    /* choice of effective seq # determination algorithm      */
  double               re_target;	 /* rel entropy target for effn eweighting, if set; or -1.0*/
  double               esigma;		 /* min total rel ent parameter for effn entropy weights   */
  double               eid;		 /* %id threshold for effn clustering                      */
  double               eset;		 /* effective sequence number, if --eset; or -1.0          */

  /* Run-to-run variation due to random number generation                                          */
  ESL_RANDOMNESS      *r;	         /* RNG for E-value calibration simulations                */
  int                  do_reseeding;	 /* TRUE to reseed, making results reproducible            */

  /* E-value parameter calibration                                                                 */
  int                  EmL;            	 /* length of sequences generated for MSV fitting          */
  int                  EmN;	         /* # of sequences generated for MSV fitting               */
  int                  EvL;            	 /* length of sequences generated for Viterbi fitting      */
  int                  EvN;	         /* # of sequences generated for Viterbi fitting           */
  int                  EfL;	         /* length of sequences generated for Forward fitting      */
  int                  EfN;	         /* # of sequences generated for Forward fitting           */
  double               Eft;	         /* tail mass used for Forward fitting                     */

  /* Choice of prior                                                                               */
  F4_PRIOR            *prior;	         /* choice of prior when parameterizing from counts        */
  int                  max_insert_len;

  /* Optional: information used for parameterizing single sequence queries                         */
  ESL_SCOREMATRIX     *S;		 /* residue score matrix                                   */
  ESL_DMATRIX         *Q;	         /* Q->mx[a][b] = P(b|a) residue probabilities             */
  double               popen;         	 /* gap open probability                                   */
  double               pextend;          /* gap extend probability                                 */

  double               w_beta;    /*beta value used to compute W (window length)   */
  int                  w_len;     /*W (window length)  explicitly set */

  const ESL_ALPHABET  *abc;		 /* COPY of alphabet                                       */
  char errbuf[eslERRBUFSIZE];            /* informative message on model construction failure      */
} F4_BUILDER;


/*****************************************************************
 * 4. F4_TRACE:  a traceback (alignment of seq to profile).
 *****************************************************************/

/* Traceback structure for alignment of a model to a sequence.
 *
 * A traceback only makes sense in a triplet (tr, gm, dsq), for a
 * given profile or HMM (with nodes 1..M) and a given digital sequence
 * (with positions 1..L).
 * 
 * A traceback may be relative to a profile (usually) or to a core
 * model (as a special case in model construction; see build.c). You
 * can tell the difference by looking at the first statetype,
 * tr->st[0]; if it's a f4T_S, it's for a profile, and if it's f4T_B,
 * it's for a core model.
 * 
 * A "profile" trace uniquely has S,N,C,T,J states and their
 * transitions; it also can have B->Mk and Mk->E internal entry/exit
 * transitions for local alignments. It may not contain X states.
 *
 * A "core" trace may contain I0, IM, and D1 states and their
 * transitions. A "core" trace can also have B->X->{MDI}k and
 * {MDI}k->X->E transitions as a special hack in a build procedure, to
 * deal with the case of a local alignment fragment implied by an
 * input alignment, which is "impossible" for a core model.
 * X "states" only appear in core traces, and only at these
 * entry/exit places; some code depends on this.
 *   
 * A profile's N,C,J states emit on transition, not on state, so a
 * path of N emits 0 residues, NN emits 1 residue, NNN emits 2
 * residues, and so on. By convention, the trace always associates an
 * emission-on-transition with the trailing (destination) state, so
 * the first N, C, or J is stored in a trace as a nonemitter (i=0).
 *
 * A i coords in a traceback are usually 1..L with respect to an
 * unaligned digital target sequence, but in the special case of
 * traces faked from existing MSAs (as in hmmbuild), the coords may
 * be 1..alen relative to an MSA's columns.
 */

/* State types */
enum f4t_statetype_e {
  f4T_BOGUS =  0,
  f4T_M     =  1,
  f4T_D     =  2,
  f4T_I     =  3,
  f4T_S     =  4,
  f4T_N     =  5,
  f4T_B     =  6, 
  f4T_E     =  7,
  f4T_C     =  8, 
  f4T_T     =  9, 
  f4T_J     = 10,
  f4T_X     = 11, 	/* missing data: used esp. for local entry/exits */
};
#define f4T_NSTATETYPES 12

typedef struct f4_trace_s {
  int    N;		/* length of traceback                       */
  int    nalloc;        /* allocated length of traceback             */
  char  *st;		/* state type code                   [0..N-1]*/
  int   *k;		/* node index; 1..M if M,D,I; else 0 [0..N-1]*/
  int   *i;		/* pos emitted in dsq, 1..L; else 0  [0..N-1]*/
  float *pp;		/* posterior prob of x_i; else 0     [0..N-1]*/
  int    M;		/* model length M (maximum k)                */
  int    L;		/* sequence length L (maximum i)             */

  /* The following section is data generated by "indexing" a trace's domains */
  int   ndom;		/* number of domains in trace (= # of B or E states) */
  int  *tfrom,   *tto;	/* locations of B/E states in trace (0..tr->N-1)     */
  int  *sqfrom,  *sqto;	/* first/last M-emitted residue on sequence (1..L)   */
  int  *hmmfrom, *hmmto;/* first/last M state on model (1..M)                */
  int   ndomalloc;	/* current allocated size of these stacks            */

} F4_TRACE;

/*****************************************************************
 * 5. F4_PROFILE: a scoring profile, and its implicit model.
 *****************************************************************/

/* Indices for special state types in the length model, gm->xsc[x][]
 */
enum f4p_xstates_e { 
  f4P_E = 0,
  f4P_N = 1,
  f4P_J = 2,
  f4P_C = 3
};
#define f4P_NXSTATES 4

/* Indices for transitions from the length modeling scores gm->xsc[][x]
 */
enum f4p_xtransitions_e {
  f4P_LOOP = 0,
  f4P_MOVE = 1
};
#define f4P_NXTRANS 2

/* Indices for transition scores gm->tsc[k][] */
/* order is optimized for dynamic programming */
enum f4p_tsc_e {
  f4P_MM = 0, 
  f4P_IM = 1, 
  f4P_DM = 2, 
  f4P_BM = 3, 
  f4P_MD = 4, 
  f4P_DD = 5, 
  f4P_MI = 6, 
  f4P_II = 7,
  f4P_ID = 8, /* Add this for f4-HMM */
  f4P_DI = 9, /* Add this for f4-HMM */
};
#define f4P_NTRANS 8

/* Indices for residue emission score vectors
 */
enum f4p_rsc_e {
  f4P_MSC = 0, 
  f4P_ISC = 1
};
#define f4P_NR 2

/* Accessing transition, emission scores */
/* _BM is specially stored off-by-one: [k-1][f4P_BM] is score for entering at Mk */
#define f4P_TSC(gm, k, s) ((gm)->tsc[(k) * f4P_NTRANS + (s)])
#define f4P_MSC(gm, k, x) ((gm)->rsc[x][(k) * f4P_NR + f4P_MSC])
#define f4P_ISC(gm, k, x) ((gm)->rsc[x][(k) * f4P_NR + f4P_ISC])


typedef struct f4_profile_s {
  float  *tsc;          /* transitions  [0.1..M-1][0..f4P_NTRANS-1], hand-indexed  */
  float **rsc;          /* emissions [0..Kp-1][0.1..M][f4P_NR], hand-indexed       */
  float   xsc[f4P_NXSTATES][f4P_NXTRANS]; /* special transitions [NECJ][LOOP,MOVE] */

  int     mode;        	/* configured algorithm mode (e.g. f4_LOCAL)               */ 
  int     L;		/* current configured target seq length                    */
  int     allocM;	/* max # of nodes allocated in this structure              */
  int     M;		/* number of nodes in the model                            */
  int     max_length;	/* calculated upper bound on emitted seq length            */
  float   nj;		/* expected # of uses of J; precalculated from loop config */

  /* Info, most of which is a copy from parent HMM:                                       */
  char  *name;			/* unique name of model                                   */
  char  *acc;			/* unique accession of model, or NULL                     */
  char  *desc;                  /* brief (1-line) description of model, or NULL           */
  char  *rf;                    /* reference line from alignment 1..M; *rf=0 means unused */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line      1..M, *cs=0 means unused */
  char  *consensus;		/* consensus residues to display in alignments, 1..M      */
  float  evparam[f4_NEVPARAM]; 	/* parameters for determining E-values, or UNSET          */
  float  cutoff[f4_NCUTOFFS]; 	/* per-seq/per-domain bit score cutoffs, or UNSET         */
  float  compo[f4_MAXABET];	/* per-model HMM filter composition, or UNSET             */

  /* Disk offset information for hmmpfam's fast model retrieval                           */
  off_t  offs[f4_NOFFSETS];     /* f4_{MFP}OFFSET, or -1                                  */

  off_t  roff;                  /* record offset (start of record); -1 if none            */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown           */

  const ESL_ALPHABET *abc;	/* copy of pointer to appropriate alphabet                */
} F4_PROFILE;

/*****************************************************************
 * 6. SSE implementation of MSV, Viterbi, and Forward routines.
 * TODO this will need to be changed to include the DI and ID transitions.
 *****************************************************************/

#define f4O_NXSTATES  4    /* special states stored: ENJC                       */
#define f4O_NXTRANS   2         /* special states all have 2 transitions: move, loop */
#define f4O_NTRANS    8    /* 7 core transitions + BMk entry                    */
enum f4o_xstates_e      { f4O_E    = 0, f4O_N    = 1,  f4O_J  = 2,  f4O_C  = 3 };
enum f4o_xtransitions_e { f4O_MOVE = 0, f4O_LOOP = 1 };
enum f4o_tsc_e          { f4O_BM   = 0, f4O_MM   = 1,  f4O_IM = 2,  f4O_DM = 3, f4O_MD   = 4, f4O_MI   = 5,  f4O_II = 6,  f4O_DD = 7 };

typedef struct f4_oprofile_s {
  /* MSVFilter uses scaled, biased uchars: 16x unsigned byte vectors                 */
  __m128i **rbv;         /* match scores [x][q]: rm, rm[0] are allocated      */
  __m128i **sbv;         /* match scores for ssvfilter                        */
  uint8_t   tbm_b;    /* constant B->Mk cost:    scaled log 2/M(M+1)       */
  uint8_t   tec_b;    /* constant E->C  cost:    scaled log 0.5            */
  uint8_t   tjb_b;    /* constant NCJ move cost: scaled log 3/(L+3)        */
  float     scale_b;    /* typically 3 / log2: scores scale to 1/3 bits      */
  uint8_t   base_b;            /* typically +190: offset of uchar scores            */
  uint8_t   bias_b;    /* positive bias to emission scores, make them >=0   */

  /* ViterbiFilter uses scaled swords: 8x signed 16-bit integer vectors              */
  __m128i **rwv;    /* [x][q]: rw, rw[0] are allocated  [Kp][Q8]         */
  __m128i  *twv;    /* transition score blocks          [8*Q8]           */
  int16_t   xw[f4O_NXSTATES][f4O_NXTRANS]; /* NECJ state transition costs            */
  float     scale_w;            /* score units: typically 500 / log(2), 1/500 bits   */
  int16_t   base_w;             /* offset of sword scores: typically +12000          */
  int16_t   ddbound_w;    /* threshold precalculated for lazy DD evaluation    */
  float     ncj_roundoff;  /* missing precision on NN,CC,JJ after rounding      */

  /* Forward, Backward use IEEE754 single-precision floats: 4x vectors               */
  __m128 **rfv;         /* [x][q]:  rf, rf[0] are allocated [Kp][Q4]         */
  __m128  *tfv;          /* transition probability blocks    [8*Q4]           */
  float    xf[f4O_NXSTATES][f4O_NXTRANS]; /* NECJ transition costs                   */

  /* Our actual vector mallocs, before we align the memory                           */
  __m128i  *rbv_mem;
  __m128i  *sbv_mem;
  __m128i  *rwv_mem;
  __m128i  *twv_mem;
  __m128   *tfv_mem;
  __m128   *rfv_mem;
  
  /* Disk offset information for hmmpfam's fast model retrieval                      */
  off_t  offs[f4_NOFFSETS];     /* f4_{MFP}OFFSET, or -1                             */

  /* Disk offset bookkeeping for h3f:                                                */
  off_t  roff;                  /* record offset (start of record); -1 if none       */
  off_t  eoff;                  /* offset to last byte of record; -1 if unknown      */

  /* Information, annotation copied from parent profile:                             */
  char  *name;      /* unique name of model                              */
  char  *acc;      /* unique accession of model, or NULL                */
  char  *desc;                  /* brief (1-line) description of model, or NULL      */
  char  *rf;                    /* reference line           1..M; *ref=0: unused     */
  char  *mm;                    /* modelmask line           1..M; *ref=0: unused     */
  char  *cs;                    /* consensus structure line 1..M, *cs=0: unused      */
  char  *consensus;    /* consensus residues for ali display, 1..M          */
  float  evparam[f4_NEVPARAM];   /* parameters for determining E-values, or UNSET     */
  float  cutoff[f4_NCUTOFFS];   /* per-seq/per-dom bit cutoffs, or UNSET             */
  float  compo[f4_MAXABET];  /* per-model HMM filter composition, or UNSET        */
  const ESL_ALPHABET *abc;  /* copy of ptr to alphabet information               */

  /* Information about current configuration, size, allocation                       */
  int    L;      /* current configured target seq length              */
  int    M;      /* model length                                      */
  int    max_length;    /* upper bound on emitted sequence length            */
  int    allocM;    /* maximum model length currently allocated for      */
  int    allocQ4;    /* f4_NQF(allocM): alloc size for tf, rf             */
  int    allocQ8;    /* f4_NQW(allocM): alloc size for tw, rw             */
  int    allocQ16;    /* f4_NQB(allocM): alloc size for rb                 */
  int    mode;      /* currently must be f4_LOCAL                        */
  float  nj;      /* expected # of J's: 0 or 1, uni vs. multihit       */

  int    clone;                 /* this optimized profile structure is just a copy   */
                                /* of another profile structre.  all pointers of     */
                                /* this structure should not be freed.               */
} F4_OPROFILE;

/*****************************************************************
 * 7. F4_GMX: a "generic" dynamic programming matrix
 *****************************************************************/

enum f4g_scells_e {
  f4G_M = 0,
  f4G_I = 1,
  f4G_D = 2,
};
#define f4G_NSCELLS 3

enum f4g_xcells_e {
  f4G_E  = 0,
  f4G_N  = 1,
  f4G_J  = 2,
  f4G_B  = 3,
  f4G_C  = 4
};
#define f4G_NXCELLS 5

typedef struct f4_gmx_s {
  int  M;		/* actual model dimension (model 1..M)    */
  int  L;		/* actual sequence dimension (seq 1..L)   */
  
  int      allocR;      /* current allocated # of rows : L+1 <= validR <= allocR                */
  int      validR;	/* # of rows actually pointing at DP memory                             */
  int      allocW;	/* current set row width :  M+1 <= allocW                               */
  uint64_t ncells;	/* total # of allocated cells in 2D matrix : ncells >= (validR)(allocW) */

  float **dp;           /* logically [0.1..L][0.1..M][0..f4G_NSCELLS-1]; indexed [i][k*f4G_NSCELLS+s] */
  float  *xmx;          /* logically [0.1..L][0..f4G_NXCELLS-1]; indexed [i*f4G_NXCELLS+s]            */

  float  *dp_mem;
} F4_GMX;

/* Macros below implement indexing idioms for generic DP routines.
 * They require the following setup, for profile <gm> and matrix <gx>:
 *   float const *tsc = gm->tsc;
 *   float      **dp  = gx->dp;
 *   float       *xmx = gx->xmx;
 * and for each row i (target residue x_i in digital seq <dsq>):
 *   float const *rsc = gm->rsc[dsq[i]];
 */
#define MMX(i,k) (dp[(i)][(k) * f4G_NSCELLS + f4G_M])
#define IMX(i,k) (dp[(i)][(k) * f4G_NSCELLS + f4G_I])
#define DMX(i,k) (dp[(i)][(k) * f4G_NSCELLS + f4G_D])
#define XMX(i,s) (xmx[(i) * f4G_NXCELLS + (s)])

#define TSC(s,k) (tsc[(k) * f4P_NTRANS + (s)])
#define MSC(k)   (rsc[(k) * f4P_NR     + f4P_MSC])
#define ISC(k)   (rsc[(k) * f4P_NR     + f4P_ISC])

/* Flags that control F4_GMX debugging dumps */
#define f4_HIDE_SPECIALS (1<<0)
#define f4_SHOW_LOG      (1<<1)

/*****************************************************************
 * 8. F4_OMX: a one-row dynamic programming matrix
 *****************************************************************/

enum f4x_scells_e { f4X_M = 0, f4X_D = 1, f4X_I = 2 };
#define f4X_NSCELLS 3

/* Besides ENJBC states, we may also store a rescaling factor on each row  */
enum f4x_xcells_e { f4X_E = 0, f4X_N = 1, f4X_J = 2, f4X_B = 3, f4X_C = 4, f4X_SCALE = 5 }; 
#define f4X_NXCELLS 6

/* 
 * 
 * dpf[][] 
 *    to access M(i,k) for i=0,1..L; k=1..M:  dpf[i][(k-1)/4 + f4X_M].element[(k-1)%4]
 * 
 * xmx[] arrays for individual special states:
 *    xmx[ENJBC] = [0 1 2 3][4 5 6 7]..[L-2 L-1 L x]     XRQ >= (L/4)+1
 *    to access B[i] for example, for i=0..L:   xmx[B][i/4].x[i%4]  (quad i/4; element i%4).
 */  
typedef struct f4_omx_s {
  int       M;      /* current actual model dimension                              */
  int       L;      /* current actual sequence dimension                           */

  /* The main dynamic programming matrix for M,D,I states                                      */
  __m128  **dpf;    /* striped DP matrix for [0,1..L][0..Q-1][MDI], float vectors  */
  __m128i **dpw;    /* striped DP matrix for [0,1..L][0..Q-1][MDI], sword vectors  */
  __m128i **dpb;    /* striped DP matrix for [0,1..L][0..Q-1] uchar vectors        */
  void     *dp_mem;    /* DP memory shared by <dpb>, <dpw>, <dpf>                     */
  int       allocR;    /* current allocated # rows in dp{uf}. allocR >= validR >= L+1 */
  int       validR;    /* current # of rows actually pointing at DP memory            */
  int       allocQ4;    /* current set row width in <dpf> quads:   allocQ4*4 >= M      */
  int       allocQ8;    /* current set row width in <dpw> octets:  allocQ8*8 >= M      */
  int       allocQ16;   /* current set row width in <dpb> 16-mers: allocQ16*16 >= M    */
  int64_t   ncells;     /* current allocation size of <dp_mem>, in accessible cells    */

  /* The X states (for full,parser; or NULL, for scorer)                                       */
  float    *xmx;          /* logically [0.1..L][ENJBCS]; indexed [i*f4X_NXCELLS+s]       */
  void     *x_mem;    /* X memory before 16-byte alignment                           */
  int       allocXR;    /* # of rows allocated in each xmx[] array; allocXR >= L+1     */
  float     totscale;    /* log of the product of all scale factors (0.0 if unscaled)   */
  int       has_own_scales;  /* TRUE to use own scale factors; FALSE if scales provided     */

  /* Parsers,scorers only hold a row at a time, so to get them to dump full matrix, it
   * must be done during a DP calculation, after each row is calculated 
   */
  int     debugging;    /* TRUE if we're in debugging mode                             */
  FILE   *dfp;      /* output stream for diagnostics                               */
} F4_OMX;

/* ?MXo(q) access macros work for either uchar or float, so long as you
 * init your "dp" to point to the appropriate array.
 */
#define MMXo(q)   (dp[(q) * f4X_NSCELLS + f4X_M])
#define DMXo(q)   (dp[(q) * f4X_NSCELLS + f4X_D])
#define IMXo(q)   (dp[(q) * f4X_NSCELLS + f4X_I])
#define XMXo(i,s) (xmx[(i) * f4X_NXCELLS + s])

/* and this version works with a ptr to the approp DP row. */
#define MMO(dp,q) ((dp)[(q) * f4X_NSCELLS + f4X_M])
#define DMO(dp,q) ((dp)[(q) * f4X_NSCELLS + f4X_D])
#define IMO(dp,q) ((dp)[(q) * f4X_NSCELLS + f4X_I])

/*****************************************************************
 * 9. Routines in DUMMER's exposed API.
 *****************************************************************/

/* dmmbuild.c */
extern double f4_MeanMatchRelativeEntropy(const F4_HMM *hmm, const F4_BG *bg);

/* f4_bg.c */
extern F4_BG *f4_bg_Create(const ESL_ALPHABET *abc);
extern F4_BG *f4_bg_CreateUniform(const ESL_ALPHABET *abc);
extern F4_BG *f4_bg_Clone(const F4_BG *bg);
extern int    f4_bg_Dump(FILE *ofp, const F4_BG *bg);
extern void   f4_bg_Destroy(F4_BG *bg);
extern int    f4_bg_SetLength(F4_BG *bg, int L);
extern int    f4_bg_NullOne(const F4_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

extern int    f4_bg_Read(char *bgfile, F4_BG *bg, char *errbuf);
extern int    f4_bg_Write(FILE *fp, F4_BG *bg);

extern int    f4_bg_SetFilter  (F4_BG *bg, int M, const float *compo);
extern int    f4_bg_FilterScore(F4_BG *bg, const ESL_DSQ *dsq, int L, float *ret_sc);

/* f4_builder.c */
extern F4_BUILDER *f4_builder_Create(const ESL_GETOPTS *go, const ESL_ALPHABET *abc);
extern int         f4_builder_LoadScoreSystem(F4_BUILDER *bld, const char *matrix,                  double popen, double pextend, F4_BG *bg);
extern int         f4_builder_SetScoreSystem (F4_BUILDER *bld, const char *mxfile, const char *env, double popen, double pextend, F4_BG *bg);
extern void        f4_builder_Destroy(F4_BUILDER *bld);

extern int f4_Builder      (F4_BUILDER *bld, ESL_MSA *msa, F4_BG *bg, F4_HMM **opt_hmm, F4_TRACE ***opt_trarr, F4_PROFILE **opt_gm, F4_OPROFILE **opt_om, ESL_MSA **opt_postmsa);
extern int f4_SingleBuilder(F4_BUILDER *bld, ESL_SQ *sq,   F4_BG *bg, F4_HMM **opt_hmm, F4_TRACE  **opt_tr,    F4_PROFILE **opt_gm, F4_OPROFILE **opt_om); 
extern int f4_Builder_MaxLength      (F4_HMM *hmm, double emit_thresh);

/* f4_hmmfile.c */
extern int f4_hmmfile_WriteASCII(FILE *fp, int format, F4_HMM *hmm);

/* f4_hmm.c */
/*      1. The F4_HMM object: allocation, initialization, destruction. */
extern F4_HMM *f4_hmm_Create(int M, const ESL_ALPHABET *abc);
extern F4_HMM *f4_hmm_CreateShell(void);
extern int     f4_hmm_CreateBody(F4_HMM *hmm, int M, const ESL_ALPHABET *abc);
extern void    f4_hmm_Destroy(F4_HMM *hmm);
extern int     f4_hmm_CopyParameters(const F4_HMM *src, F4_HMM *dest);
extern F4_HMM *f4_hmm_Clone(const F4_HMM *hmm);
extern int     f4_hmm_Zero(F4_HMM *hmm);
extern char    f4_hmm_EncodeStatetype(char *typestring);
extern char   *f4_hmm_DecodeStatetype(char st);
/*      2. Convenience routines for setting fields in an HMM. */
extern int     f4_hmm_SetName       (F4_HMM *hmm, char *name);
extern int     f4_hmm_SetAccession  (F4_HMM *hmm, char *acc);
extern int     f4_hmm_SetDescription(F4_HMM *hmm, char *desc);
extern int     f4_hmm_AppendComlog  (F4_HMM *hmm, int argc, char **argv);
extern int     f4_hmm_SetCtime      (F4_HMM *hmm);
extern int     f4_hmm_SetComposition(F4_HMM *hmm);
extern int     f4_hmm_SetConsensus  (F4_HMM *hmm, ESL_SQ *sq);
/*      3. Renormalization and rescaling counts in core HMMs. */
extern int     f4_hmm_Scale      (F4_HMM *hmm, double scale);
extern int     f4_hmm_ScaleExponential(F4_HMM *hmm, double exp);
extern int     f4_hmm_Renormalize(F4_HMM *hmm);
/*      4. Other routines in the API */
extern int     f4_hmm_CalculateOccupancy(const F4_HMM *hmm, float *mocc, float *iocc);

/* f4_prior.c */
extern F4_PRIOR  *f4_prior_CreateAmino(void);
extern F4_PRIOR  *f4_prior_CreateNucleic(void);
extern F4_PRIOR  *f4_prior_CreateLaplace(const ESL_ALPHABET *abc);
extern void       f4_prior_Destroy(F4_PRIOR *pri);

extern int        p7_ParameterEstimation(F4_HMM *hmm, const F4_PRIOR *pri); //TODO will be removed in future
extern int        f4_ParameterEstimation(F4_HMM *hmm, const F4_PRIOR *pri);

/* f4_trace.c */
extern F4_TRACE *f4_trace_Create(void);
extern F4_TRACE *f4_trace_CreateWithPP(void);
extern int  f4_trace_Reuse(F4_TRACE *tr);
extern int  f4_trace_Grow(F4_TRACE *tr);
extern int  f4_trace_GrowIndex(F4_TRACE *tr);
extern int  f4_trace_GrowTo(F4_TRACE *tr, int N);
extern int  f4_trace_GrowIndexTo(F4_TRACE *tr, int ndom);
extern void f4_trace_Destroy(F4_TRACE *tr);
extern void f4_trace_DestroyArray(F4_TRACE **tr, int N);

extern int  f4_trace_GetDomainCount   (const F4_TRACE *tr, int *ret_ndom);
extern int  f4_trace_GetStateUseCounts(const F4_TRACE *tr, int *counts);
extern int  f4_trace_GetDomainCoords  (const F4_TRACE *tr, int which, int *ret_i1, int *ret_i2,
				       int *ret_k1, int *ret_k2);

extern int   f4_trace_Validate(const F4_TRACE *tr, const ESL_ALPHABET *abc, const ESL_DSQ *dsq, char *errbuf);
extern int   f4_trace_Dump(FILE *fp, const F4_TRACE *tr, const F4_PROFILE *gm, const ESL_DSQ *dsq);
extern int   f4_trace_Compare(F4_TRACE *tr1, F4_TRACE *tr2, float pptol);
extern int   f4_trace_Score(F4_TRACE *tr, ESL_DSQ *dsq, F4_PROFILE *gm, float *ret_sc);
extern int   f4_trace_SetPP(F4_TRACE *tr, const F4_GMX *pp);
extern float f4_trace_GetExpectedAccuracy(const F4_TRACE *tr);

extern int  f4_trace_Append(F4_TRACE *tr, char st, int k, int i);
extern int  f4_trace_AppendWithPP(F4_TRACE *tr, char st, int k, int i, float pp);
extern int  f4_trace_Reverse(F4_TRACE *tr);
extern int  f4_trace_Index(F4_TRACE *tr);

extern int  f4_trace_FauxFromMSA(ESL_MSA *msa, int *matassign, int optflags, F4_TRACE **tr);

extern int  p7_trace_Count(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr);
extern int  f4_trace_Count(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr);

extern int  f4_trace_Estimate(F4_HMM *hmm, ESL_DSQ *dsq, float wt, F4_TRACE *tr, double **letter_probs, double *background_probs);
extern int  letter_probs_build(int profile_length, int num_of_letters, double ***ret_letter_probs, double **ret_background_probs);
extern int  letter_probs_count(int sequence_length, ESL_DSQ *dsq, float wt, double **letter_probs, double *background_probs);
extern int letter_probs_normalize(int profile_length, int num_of_letters, double **letter_probs, double *background_probs);
extern void letter_probs_destroy(int sequence_length, double **letter_probs, double *background_probs);

/* seqmodel.c */
extern int f4_Seqmodel(const ESL_ALPHABET *abc, ESL_DSQ *dsq, int M, char *name,
		       ESL_DMATRIX *P, float *f, double popen, double pextend,
		       F4_HMM **ret_hmm);

/* eweight.c */
extern int f4_EntropyWeight(const F4_HMM *hmm, const F4_BG *bg, const F4_PRIOR *pri, double infotarget, double *ret_Neff);
extern int f4_EntropyWeight_exp(const F4_HMM *hmm, const F4_BG *bg, const F4_PRIOR *pri, double etarget, double *ret_exp);

/* evalues.c */
extern int f4_Calibrate(F4_HMM *hmm, F4_BUILDER *cfg_b, ESL_RANDOMNESS **byp_rng, F4_BG **byp_bg, F4_PROFILE **byp_gm, F4_OPROFILE **byp_om);
extern int f4_Lambda(F4_HMM *hmm, F4_BG *bg, double *ret_lambda);
extern int f4_MSVMu     (ESL_RANDOMNESS *r, F4_OPROFILE *om, F4_BG *bg, int L, int N, double lambda,               double *ret_mmu);
extern int f4_ViterbiMu (ESL_RANDOMNESS *r, F4_OPROFILE *om, F4_BG *bg, int L, int N, double lambda,               double *ret_vmu);
extern int f4_Tau       (ESL_RANDOMNESS *r, F4_OPROFILE *om, F4_BG *bg, int L, int N, double lambda, double tailp, double *ret_tau);

/* tracealign.c */
extern int f4_tracealign_MSA(const ESL_MSA *premsa, F4_TRACE **tr, int M, int optflags, ESL_MSA **ret_postmsa);

/* build.c */
extern int f4_Fastmodelmaker(ESL_MSA *msa, float symfrac, F4_BUILDER *bld, F4_HMM **ret_hmm, F4_TRACE ***ret_tr);

/* modelconfig.c */
extern int f4_ProfileConfig(const F4_HMM *hmm, const F4_BG *bg, F4_PROFILE *gm, int L, int mode);
extern int f4_ReconfigLength(F4_PROFILE *gm, int L);

/* f4_oprofile.c */
extern void f4_oprofile_Destroy(F4_OPROFILE *om);
extern F4_OPROFILE * f4_oprofile_Create(int allocM, const ESL_ALPHABET *abc);
extern int f4_oprofile_Convert(const F4_PROFILE *gm, F4_OPROFILE *om);

/* f4_profile.c */
extern void f4_profile_Destroy(F4_PROFILE *gm);
extern F4_PROFILE * f4_profile_Create(int allocM, const ESL_ALPHABET *abc);
extern int f4_oprofile_ReconfigLength(F4_OPROFILE *om, int L);
extern int f4_profile_GetT(const F4_PROFILE *gm, char st1, int k1, char st2, int k2, float *ret_tsc);
extern int f4_profile_IsLocal(const F4_PROFILE *gm);
extern int f4_profile_IsMultihit(const F4_PROFILE *gm);

/* vitfilter.c */
extern int f4_ViterbiFilter(const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, F4_OMX *ox, float *ret_sc);

/* msvfilter.c */
extern int f4_MSVFilter(const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, F4_OMX *ox, float *ret_sc);

/* ssvfilter.c */
extern int f4_SSVFilter(const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, float *ret_sc);

/* f4_omx.c */
extern F4_OMX * f4_omx_Create(int allocM, int allocL, int allocXL);
extern void f4_omx_Destroy(F4_OMX *ox);

/* fwdbck.c */
extern int f4_ForwardParser(const ESL_DSQ *dsq, int L, const F4_OPROFILE *om, F4_OMX *ox, float *opt_sc);