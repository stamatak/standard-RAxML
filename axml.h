/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *
 *  and
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 *
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses
 *  with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include <assert.h>
#include <stdint.h>


#ifdef _USE_PTHREADS

#include <pthread.h>

#endif


#define MAX_TIP_EV     0.999999999 /* max tip vector value, sum of EVs needs to be smaller than 1.0, otherwise the numerics break down */
#define smoothings     32          /* maximum smoothing passes through tree */
#define iterations     10          /* maximum iterations of iterations per insert */
#define newzpercycle   1           /* iterations of makenewz per tree traversal */
#define nmlngth        256         /* number of characters in species name */
#define deltaz         0.00001     /* test of net branch length change in update */
#define defaultz       0.9         /* value of z assigned as starting point */
#define unlikely       -1.0E300    /* low likelihood for initialization */


#define SUMMARIZE_LENGTH -3
#define SUMMARIZE_LH     -2
#define NO_BRANCHES      -1

#define MASK_LENGTH 32
#define VECTOR_LENGTH (NUM_BRANCHES / MASK_LENGTH)
#define GET_BITVECTOR_LENGTH(x) ((x % MASK_LENGTH) ? (x / MASK_LENGTH + 1) : (x / MASK_LENGTH))

#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define zmax (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */

#define twotothe256  \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0
                                                     /*  2**256 (exactly)  */

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

#define badRear         -1

#define NUM_BRANCHES   128

#define TRUE             1
#define FALSE            0



#define LIKELIHOOD_EPSILON 0.0000001

#define AA_SCALE 10.0
#define AA_SCALE_PLUS_EPSILON 10.001

/* ALPHA_MIN is critical -> numerical instability, eg for 4 discrete rate cats                    */
/* and alpha = 0.01 the lowest rate r_0 is                                                        */
/* 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487             */
/* which leads to numerical problems Table for alpha settings below:                              */
/*                                                                                                */
/* 0.010000 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487    */
/* 0.010000 yielded nasty numerical bugs in at least one case !                                   */
/* 0.020000 0.00000000000000000000000000000044136090435925743185910935350715027016962154188875    */
/* 0.030000 0.00000000000000000000476844846859006690412039180149775802624789852441798419292220    */
/* 0.040000 0.00000000000000049522423236954066431210260930029681736928018820007024736185030633    */
/* 0.050000 0.00000000000050625351310359203371872643495343928538368616365517027588794007897377    */
/* 0.060000 0.00000000005134625283884191118711474021861409372524676086868566926568746566772461    */
/* 0.070000 0.00000000139080650074206434685544624965062437960128249869740102440118789672851562    */
/* 0.080000 0.00000001650681201563587066858709818343436959153791576682124286890029907226562500    */
/* 0.090000 0.00000011301977332931251259273962858978301859735893231118097901344299316406250000    */
/* 0.100000 0.00000052651925834844387815526344648331402709118265192955732345581054687500000000    */


#define ALPHA_MIN    0.02
#define ALPHA_MAX    1000.0

#define RATE_MIN     0.0000001
#define RATE_MAX     1000000.0

#define INVAR_MIN    0.0001
#define INVAR_MAX    0.9999

#define TT_MIN       0.0000001
#define TT_MAX       1000000.0

#define FREQ_MIN     0.001

/* 
   previous values between 0.001 and 0.000001

   TO AVOID NUMERICAL PROBLEMS WHEN FREQ == 0 IN PARTITIONED MODELS, ESPECIALLY WITH AA 
   previous value of FREQ_MIN was: 0.000001, but this seemed to cause problems with some 
   of the 7-state secondary structure models with some rather exotic small toy test datasets,
   on the other hand 0.001 caused problems with some of the 16-state secondary structure models

   For some reason the frequency settings seem to be repeatedly causing numerical problems
   
*/

#define ITMAX 100



#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))


#define LOG(x)  log(x)
#define EXP(x)  exp(x)






#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

#define programName        "RAxML"
#define programVersion     "7.3.0"
#define programDate        "June 2011"


#define  TREE_EVALUATION            0
#define  BIG_RAPID_MODE             1
#define  CALC_BIPARTITIONS          3
#define  SPLIT_MULTI_GENE           4
#define  CHECK_ALIGNMENT            5
#define  PER_SITE_LL                6
#define  PARSIMONY_ADDITION         7
#define  CLASSIFY_ML                9
#define  DISTANCE_MODE              11
#define  GENERATE_BS                12
#define  COMPUTE_ELW                13
#define  BOOTSTOP_ONLY              14
#define  COMPUTE_LHS                17
#define  COMPUTE_BIPARTITION_CORRELATION 18
#define  COMPUTE_RF_DISTANCE        20
#define  MORPH_CALIBRATOR           21
#define  CONSENSUS_ONLY             22
#define  MESH_TREE_SEARCH           23
#define  FAST_SEARCH                24        
#define  EPA_SITE_SPECIFIC_BIAS     27
#define  SH_LIKE_SUPPORTS           28
#define  CLASSIFY_MP                29
#define  ANCESTRAL_STATES           30

#define M_GTRCAT         1
#define M_GTRGAMMA       2
#define M_BINCAT         3
#define M_BINGAMMA       4
#define M_PROTCAT        5
#define M_PROTGAMMA      6
#define M_32CAT          7
#define M_32GAMMA        8
#define M_64CAT          9
#define M_64GAMMA        10

#define DAYHOFF      0
#define DCMUT        1
#define JTT          2
#define MTREV        3
#define WAG          4
#define RTREV        5
#define CPREV        6
#define VT           7
#define BLOSUM62     8
#define MTMAM        9
#define LG           10
#define MTART        11
#define MTZOA        12
#define PMB          13
#define HIVB         14
#define HIVW         15
#define JTTDCMUT     16
#define FLU          17 
#define PROT_FILE    18
#define GTR_UNLINKED 19
#define GTR          20  /* GTR always needs to be the last one */

#define NUM_PROT_MODELS 21


/* bipartition stuff */

#define BIPARTITIONS_ALL       0
#define GET_BIPARTITIONS_BEST  1
#define DRAW_BIPARTITIONS_BEST 2
#define BIPARTITIONS_BOOTSTOP  3
#define BIPARTITIONS_RF  4



/* bootstopping stuff */

#define BOOTSTOP_PERMUTATIONS 100
#define START_BSTOP_TEST      10

#define FC_THRESHOLD          99
#define FC_SPACING            50
#define FC_LOWER              0.99
#define FC_INIT               20

#define FREQUENCY_STOP 0
#define MR_STOP        1
#define MRE_STOP       2
#define MRE_IGN_STOP   3

#define MR_CONSENSUS 0
#define MRE_CONSENSUS 1
#define STRICT_CONSENSUS 2



/* bootstopping stuff end */


#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2

#define MIN_MODEL        -1
#define BINARY_DATA      0
#define DNA_DATA         1
#define AA_DATA          2
#define SECONDARY_DATA   3
#define SECONDARY_DATA_6 4
#define SECONDARY_DATA_7 5
#define GENERIC_32       6
#define GENERIC_64       7
#define MAX_MODEL        8

#define SEC_6_A 0
#define SEC_6_B 1
#define SEC_6_C 2
#define SEC_6_D 3
#define SEC_6_E 4

#define SEC_7_A 5
#define SEC_7_B 6
#define SEC_7_C 7
#define SEC_7_D 8
#define SEC_7_E 9
#define SEC_7_F 10

#define SEC_16   11
#define SEC_16_A 12
#define SEC_16_B 13
#define SEC_16_C 14
#define SEC_16_D 15
#define SEC_16_E 16
#define SEC_16_F 17
#define SEC_16_I 18
#define SEC_16_J 19
#define SEC_16_K 20

#define ORDERED_MULTI_STATE 0
#define MK_MULTI_STATE      1
#define GTR_MULTI_STATE     2





#define CAT         0
#define GAMMA       1
#define GAMMA_I     2




typedef  int boolean;


typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

struct ent
{
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int amountTips;
  int *supportVector;
  unsigned int bipNumber;
  unsigned int bipNumber2;
  unsigned int supportFromTreeset[2]; 
  struct ent *next;
};

typedef struct ent entry;

typedef unsigned int hashNumberType;

typedef unsigned int parsimonyNumber;

/*typedef uint_fast32_t parsimonyNumber;*/

#define PCF 32

/*
  typedef uint64_t parsimonyNumber;

  #define PCF 16


typedef unsigned char parsimonyNumber;

#define PCF 2
*/

typedef struct
{
  hashNumberType tableSize;
  entry **table;
  hashNumberType entryCount;
}
  hashtable;


struct stringEnt
{
  int nodeNumber;
  char *word;
  struct stringEnt *next;
};

typedef struct stringEnt stringEntry;
 
typedef struct
{
  hashNumberType tableSize;
  stringEntry **table;
}
  stringHashtable;




typedef struct ratec
{
  double accumulatedSiteLikelihood;
  double rate;
}
  rateCategorize;


typedef struct
{
  int tipCase;
  int pNumber;
  int qNumber;
  int rNumber;
  double qz[NUM_BRANCHES];
  double rz[NUM_BRANCHES];
} traversalInfo;

typedef struct
{
  traversalInfo *ti;
  int count;
} traversalData;


struct noderec;

typedef struct epBrData
{
  int    *countThem;
  int    *executeThem;
  unsigned int *parsimonyScore;
  double *branches;
  double *distalBranches; 
  double *likelihoods;
  double originalBranchLength;
  char branchLabel[64];
  int leftNodeNumber;
  int rightNodeNumber;
  int *leftScaling;
  int *rightScaling;
  double branchLengths[NUM_BRANCHES];
  double *left;
  double *right;
  int branchNumber;
  int jointLabel;
} epaBranchData;

typedef struct
{
  epaBranchData *epa;

  unsigned int *vector; 
  int support;   
  struct noderec *oP;
  struct noderec *oQ;
} branchInfo;








typedef struct
{
  boolean valid;
  int partitions;
  int *partitionList;
}
  linkageData;

typedef struct
{
  int entries;
  linkageData* ld;
}
  linkageList;


typedef  struct noderec
{
  unsigned int    isPresent[NUM_BRANCHES / MASK_LENGTH];
  struct noderec  *backs[NUM_BRANCHES];
  char            xs[NUM_BRANCHES];
  branchInfo      *bInf;
  double           z[NUM_BRANCHES];
  struct noderec  *next;
  struct noderec  *back;
  hashNumberType   hash;
  int              support;
  int              number;
  char             x;
}
  node, *nodeptr;

typedef struct
  {
    double lh;
    double pendantBranch;
    double distalBranch;    
    int number;
  }
  info;

typedef struct bInf {
  double likelihood;
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;




typedef  struct
{
  int              numsp;
  int              sites;
  unsigned char             **y;
  unsigned char             *y0;
  unsigned char             *yBUF;
  int              *wgt;
} rawdata;

typedef  struct {
  int             *alias;       /* site representing a pattern */
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory;
  int              endsite;     /* # of sequence patterns */
  double          *patrat;      /* rates per pattern */
  double          *patratStored; 
} cruncheddata;




typedef struct {
  int     states;
  int     maxTipStates;
  size_t    lower;
  size_t     upper;
  size_t     width;
  int     dataType;
  int     protModels;
  boolean usePredefinedProtFreqs;
  int     mxtips;
  int     numberOfCategories;
  int             **expVector;
  double          **xVector;
  size_t             *xSpaceVector;
  size_t             *expSpaceVector;
 
  unsigned char            **yVector;
 
  char   *partitionName;
  char   proteinSubstitutionFileName[2048];
  double externalAAMatrix[420];
  double *sumBuffer;
 
  double *gammaRates;

  double *EIGN;
  double *EV;



  double *EI;



  

  double *left;
  double *right;

  


  double *frequencies;
  double *tipVector;
 
  double *substRates;
  double *perSiteLL;
  
  double *perSiteRates;

  double *wr;
  double *wr2;

  

  unsigned int    *globalScaler;
  double          *globalScalerDouble;
  int    *perSiteAAModel;
  int    *wgt;
  int    *invariant;
  int    *rateCategory;
  int    *symmetryVector;
  int    *frequencyGrouping;
  boolean nonGTR;
  double alpha;
  double propInvariant;

  int gapVectorLength;
  unsigned int *gapVector;
  double *gapColumn;

  size_t initialGapVectorSize;

  size_t parsimonyLength;
  parsimonyNumber *parsVect; 

} pInfo;



typedef struct 
{
  int left;
  int right;
  double likelihood;
} lhEntry;


typedef struct 
{
  int count;
  int size;
  lhEntry *entries;
} lhList;


typedef struct List_{
  void *value; 			
  struct List_ *next; 
} List;


typedef struct {
  double EIGN[19];             
  double EV[400];                
  double EI[380];
  double substRates[190];        
  double frequencies[20];      
  double tipVector[460];
  double fracchange[1];
  double left[1600];
  double right[1600];
} siteAAModels;

typedef  struct  {
  siteAAModels siteProtModel[2 * (NUM_PROT_MODELS - 3)];

  boolean estimatePerSiteAA;

  boolean useGappedImplementation;
  boolean saveMemory;
  
  int    *resample;

  int numberOfBranches;
  int    numberOfTipsForInsertion;
  int    *inserts;
  int    branchCounter;
  
  int *ti; 
  
  int numberOfTrees; 

  stringHashtable  *nameHash;

  pInfo            *partitionData;
  pInfo            *initialPartitionData;
  pInfo            *extendedPartitionData;

  int              *dataVector;
  int              *initialDataVector;
  int              *extendedDataVector;

  int              *patternPosition;
  int              *columnPosition;

  char             *secondaryStructureInput;

  boolean          *executeModel;

  double           *perPartitionLH;
  double           *storedPerPartitionLH;

  traversalData td[NUM_BRANCHES];

  unsigned int *parsimonyScore;

  int              maxCategories;

  double           *wr;
  double           *wr2;
  double           *sumBuffer;
  double           *perSiteLL;  
  double           coreLZ[NUM_BRANCHES];
  int              modelNumber;
  int              multiBranch;
  int              numBranches;
  int              bootStopCriterion;
  int              consensusType;
  double           wcThreshold;

  int              multiGene;

  nodeptr          startVector[NUM_BRANCHES];

  nodeptr          removeNodes[NUM_BRANCHES];
  nodeptr          leftNodes[NUM_BRANCHES];
  nodeptr          rightNodes[NUM_BRANCHES]; 
  nodeptr          storedBacks[NUM_BRANCHES];
  lhList          *likelihoodList[NUM_BRANCHES];
 
  double           zLeft[NUM_BRANCHES];
  double           zRight[NUM_BRANCHES];
  double           zDown[NUM_BRANCHES];

  boolean          useFastScaling;
 
  branchInfo	   *bInf;

  int              multiStateModel;


  size_t innerNodes;

  boolean curvatOK[NUM_BRANCHES];
  /* the stuff below is shared among DNA and AA, span does
     not change depending on datatype */

  double           *invariants;
  double           *fracchanges;

  /* model stuff end */

  unsigned char             **yVector;
  int              secondaryStructureModel;
  int              discreteRateCategories;
  int              originalCrunchedLength;
  int              fullSites;
  int              *originalModel;
  int              *originalDataVector;
  int              *originalWeights;
  int              *secondaryStructurePairs;


  double            *partitionContributions;
  double            fracchange;
  double            lhCutoff;
  double            lhAVG;
  unsigned long     lhDEC;
  unsigned long     itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;
  int               rateHetModel;

  double           startLH;
  double           endLH;
  double           likelihood;
  double          *likelihoods;
  int             *invariant;
  node           **nodep;
  node            *start;
  int              mxtips;
  int              mxtipsVector[NUM_BRANCHES];
  int              *model;

  int              *constraintVector;
  int              numberOfSecondaryColumns;
  boolean          searchConvergenceCriterion;
  int              ntips;
  int              nextnode;
  int              NumberOfModels;
  int              parsimonyLength;
  
  int              checkPointCounter;
  int              treeID;
  int              numberOfOutgroups;
  int             *outgroupNums;
  char           **outgroups;
  boolean          useEpaHeuristics;
  double           fastEPAthreshold;
  boolean          bigCutoff;
  boolean          partitionSmoothed[NUM_BRANCHES];
  boolean          partitionConverged[NUM_BRANCHES];
  boolean          rooted;
  boolean          grouped;
  boolean          constrained;
  boolean          doCutoff;
  boolean          catOnly;
  rawdata         *rdta;
  cruncheddata    *cdta;

  char **nameList;
  char *tree_string;
  int treeStringLength;
  unsigned int bestParsimony;
  double bestOfNode;
  nodeptr removeNode;
  nodeptr insertNode;

  double zqr[NUM_BRANCHES];
  double currentZQR[NUM_BRANCHES];

  double currentLZR[NUM_BRANCHES];
  double currentLZQ[NUM_BRANCHES];
  double currentLZS[NUM_BRANCHES];
  double currentLZI[NUM_BRANCHES];
  double lzs[NUM_BRANCHES];
  double lzq[NUM_BRANCHES];
  double lzr[NUM_BRANCHES];
  double lzi[NUM_BRANCHES];

 
  int mr_thresh;

  boolean wasRooted;
  nodeptr leftRootNode;
  nodeptr rightRootNode;
  int rootLabel;
  
 

#ifdef _USE_PTHREADS

  double *ancestralStates;

  hashtable *h;
 
    
   

  double *temporaryVector;  
  int    *temporaryScaling;
  double *temporarySumBuffer;
  size_t contiguousVectorLength;
  size_t contiguousScalingLength;  
 
 

  int *contiguousRateCategory;
  int *contiguousWgt;
  int *contiguousInvariant;  

  unsigned char **contiguousTips;

  double *contiguousWR;
  double *contiguousWR2;
  
  
  unsigned char *y_ptr;
  
  double *wrPtr;
  double *wr2Ptr;

  

  double *perSiteLLPtr;
  int    *wgtPtr;
  int    *invariantPtr;
  int    *rateCategoryPtr;

  int threadID;
  double lower_spacing;
  double upper_spacing;
  double *lhs;

  /* stuff for parallel MRE */

  /* added by aberer */  
  entry **sbi;
  entry **sbw;
  int *len;   

  /* mre */
  int sectionEnd;
  int bipStatusLen;
  entry **entriesOfSection;
  int recommendedAmountJobs;
 

  /* used for printBip */
  boolean *hasAncestor; 
  List **listOfDirectChildren;
  unsigned int bitVectorLength;                 /* we now need this also in sequential mode */
  entry **consensusBips;
  int consensusBipLen;          /* also used in mre */  
  int maxBips;
  int *bipStatus;


  /* for parallel hash insert */
  pthread_mutex_t** mutexesForHashing; 


 

#endif


} tree;


/***************************************************************/

typedef struct
{
  double z[NUM_BRANCHES];
  nodeptr p, q;
  int cp, cq;
}
  connectRELL, *connptrRELL;

typedef  struct
{
  connectRELL     *connect; 
  int             start;
  double          likelihood;
}
  topolRELL;


typedef  struct
{
  int max;
  topolRELL **t;
}
  topolRELL_LIST;


/**************************************************************/



typedef struct conntyp {
    double           z[NUM_BRANCHES];           /* branch length */
    node            *p, *q;       /* parent and child sectors */
    void            *valptr;      /* pointer to value of subtree */
    int              descend;     /* pointer to first connect of child */
    int              sibling;     /* next connect from same parent */
    } connect, *connptr;

typedef  struct {
    double           likelihood;
  int              initialTreeNumber;
    connect         *links;       /* pointer to first connect (start) */
    node            *start;
    int              nextlink;    /* index of next available connect */
                                  /* tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;
    int              scrNum;      /* position in sorted list of scores */
    int              tplNum;      /* position in sorted list of trees */

    } topol;

typedef struct {
    double           best;        /* highest score saved */
    double           worst;       /* lowest score saved */
    topol           *start;       /* starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /* maximum topologies to save */
    int              nvalid;      /* number of topologies saved */
    int              ninit;       /* number of topologies initialized */
    int              numtrees;    /* number of alternatives tested */
    boolean          improved;
    } bestlist;

typedef  struct {
  int              categories;
  int              model;
  int              bestTrav;
  int              max_rearrange;
  int              stepwidth;
  int              initial;
  boolean          initialSet;
  int              mode;
  long             boot;
  long             rapidBoot;
  boolean          bootstrapBranchLengths;
  boolean          restart;
  boolean          useWeightFile;
  boolean          useMultipleModel;
  boolean          constraint;
  boolean          grouping;
  boolean          randomStartingTree;
  boolean          useInvariant;
  int            protEmpiricalFreqs;
  int            proteinMatrix;
  int            checkpoints;
  int            startingTreeOnly;
  int            multipleRuns;
  long           parsimonySeed;
  boolean        perGeneBranchLengths;
  boolean        likelihoodTest;
  boolean        outgroup;
  boolean        permuteTreeoptimize;
  boolean        allInOne;
  boolean        generateBS;
  boolean        bootStopping;
  boolean        useExcludeFile;
  boolean        userProteinModel;
  boolean        computeELW;
  boolean        computeDistance;
  boolean        compressPatterns;
  boolean        useSecondaryStructure; 
  double         likelihoodEpsilon;
  double         gapyness;
  int            similarityFilterMode;
  double        externalAAMatrix[420];
  boolean       useFloat;
  boolean       readTaxaOnly;
  int           meshSearch;  
  boolean       veryFast;
  boolean       useBinaryModelFile;
  boolean       leaveDropMode;
  int           slidingWindowSize;
} analdef;


typedef struct 
{
  int leftLength;
  int rightLength;
  int eignLength;
  int evLength;
  int eiLength;
  int substRatesLength;
  int frequenciesLength;
  int tipVectorLength;
  int symmetryVectorLength;
  int frequencyGroupingLength;

  boolean nonGTR;

  int undetermined;

  const char *inverseMeaning;

  int states;

  boolean smoothFrequencies;

  const unsigned  int *bitVector;

} partitionLengths;

/****************************** FUNCTIONS ****************************************************/

extern void computePlacementBias(tree *tr, analdef *adef);

extern int lookupWord(char *s, stringHashtable *h);

extern void getDataTypeString(tree *tr, int model, char typeOfData[1024]);

extern unsigned int genericBitCount(unsigned int* bitVector, unsigned int bitVectorLength);
extern int countTips(nodeptr p, int numsp);
extern entry *initEntry(void);
extern void computeRogueTaxa(tree *tr, char* treeSetFileName, analdef *adef);
extern unsigned int precomputed16_bitcount(unsigned int n);


extern double evaluateGenericMulti (tree *tr, nodeptr p, int model);
extern void setupPointerMesh(tree *tr);
extern void determineFullTraversalMulti(nodeptr p, tree *tr);
extern void computeTraversalInfoMulti(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model);
extern double evaluateIterativeMulti(tree *, boolean writeVector);
extern void newviewIterativeMulti (tree *tr);
extern void getxsnode (nodeptr p, int model);
extern void findNext(nodeptr p, tree *tr, nodeptr *result);




extern partitionLengths * getPartitionLengths(pInfo *p);
extern boolean getSmoothFreqs(int dataType);
extern const unsigned int *getBitVector(int dataType);
extern int getUndetermined(int dataType);
extern int getStates(int dataType);
extern char getInverseMeaning(int dataType, unsigned char state);
extern void printModelParams(tree *tr, analdef *adef);
extern double gettime ( void );
extern int gettimeSrand ( void );
extern double randum ( long *seed );

extern void getxnode ( nodeptr p );
extern void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
extern void hookupDefault ( nodeptr p, nodeptr q, int numBranches);
extern boolean whitechar ( int ch );
extern void errorExit ( int e );
extern void printResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBootstrapResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBipartitionResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printLog ( tree *tr, analdef *adef, boolean finalPrint );
extern void printStartingTree ( tree *tr, analdef *adef, boolean finalPrint );
extern void writeInfoFile ( analdef *adef, tree *tr, double t );
extern int main ( int argc, char *argv[] );
extern void calcBipartitions ( tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName );
extern void initReversibleGTR (tree *tr, int model);
extern double LnGamma ( double alpha );
extern double IncompleteGamma ( double x, double alpha, double ln_gamma_alpha );
extern double PointNormal ( double prob );
extern double PointChi2 ( double prob, double v );
extern void makeGammaCats (double alpha, double *gammaRates, int K);
extern void initModel ( tree *tr, rawdata *rdta, cruncheddata *cdta, analdef *adef );
extern void doAllInOne ( tree *tr, analdef *adef );

extern void classifyML(tree *tr, analdef *adef);
extern void classifyMP(tree *tr, analdef *adef);
extern void markTips(nodeptr p, int *perm, int maxTips);
extern char *Tree2StringClassify(char *treestr, tree *tr, int *inserts, 
				 boolean  originalTree, boolean jointLabels, boolean likelihood);


extern void doBootstrap ( tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta );
extern void doInference ( tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta );
extern void resetBranches ( tree *tr );
extern void modOpt ( tree *tr, analdef *adef , boolean resetModel, double likelihoodEpsilon, boolean testGappedImplementation);


extern void parsePartitions ( analdef *adef, rawdata *rdta, tree *tr);
extern void computeBOOTRAPID (tree *tr, analdef *adef, long *radiusSeed);
extern void optimizeRAPID ( tree *tr, analdef *adef );
extern void thoroughOptimization ( tree *tr, analdef *adef, topolRELL_LIST *rl, int index );
extern int treeOptimizeThorough ( tree *tr, int mintrav, int maxtrav);

extern int checker ( tree *tr, nodeptr p );
extern int randomInt ( int n );
extern void makePermutation ( int *perm, int n, analdef *adef );
extern boolean tipHomogeneityChecker ( tree *tr, nodeptr p, int grouping );
extern void makeRandomTree ( tree *tr, analdef *adef );
extern void nodeRectifier ( tree *tr );
extern void makeParsimonyTreeThorough(tree *tr, analdef *adef);
extern void makeParsimonyTree ( tree *tr, analdef *adef );
extern void makeParsimonyTreeFast(tree *tr, analdef *adef, boolean full);
extern void makeParsimonyTreeIncomplete ( tree *tr, analdef *adef );
extern void makeParsimonyInsertions(tree *tr, nodeptr startNodeQ, nodeptr startNodeR);



extern FILE *myfopen(const char *path, const char *mode);


extern boolean initrav ( tree *tr, nodeptr p );
extern void initravPartition ( tree *tr, nodeptr p, int model );
extern boolean update ( tree *tr, nodeptr p );
extern boolean smooth ( tree *tr, nodeptr p );
extern boolean smoothTree ( tree *tr, int maxtimes );
extern boolean localSmooth ( tree *tr, nodeptr p, int maxtimes );
extern boolean localSmoothMulti(tree *tr, nodeptr p, int maxtimes, int model);
extern void initInfoList ( int n );
extern void freeInfoList ( void );
extern void insertInfoList ( nodeptr node, double likelihood );
extern boolean smoothRegion ( tree *tr, nodeptr p, int region );
extern boolean regionalSmooth ( tree *tr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( tree *tr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( tree *tr, nodeptr p );
extern boolean insertBIG ( tree *tr, nodeptr p, nodeptr q, int numBranches);
extern boolean insertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern boolean testInsertBIG ( tree *tr, nodeptr p, nodeptr q );
extern void addTraverseBIG ( tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( tree *tr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern double treeOptimizeRapid ( tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt);
extern boolean testInsertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( tree *tr );
extern int determineRearrangementSetting ( tree *tr, analdef *adef, bestlist *bestT, bestlist *bt );
extern void computeBIGRAPID ( tree *tr, analdef *adef, boolean estimateModel);
extern boolean treeEvaluate ( tree *tr, double smoothFactor );
extern boolean treeEvaluatePartition ( tree *tr, double smoothFactor, int model );

extern void meshTreeSearch(tree *tr, analdef *adef, int thorough);

extern void initTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void freeTL ( topolRELL_LIST *rl);
extern void restoreTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, tree *tr, int index );

extern int  saveBestTree (bestlist *bt, tree *tr);
extern int  recallBestTree (bestlist *bt, int rank, tree *tr);
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern boolean freeBestTree ( bestlist *bt );


extern char *Tree2String ( char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, 
			   boolean rellTree, boolean finalPrint, analdef *adef, int perGene, boolean branchLabelSupport, boolean printSHSupport);
extern void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission);



extern int treeReadLen (FILE *fp, tree *tr, boolean readBranches, boolean readNodeLabels, boolean topologyOnly, analdef *adef, boolean completeTree);
extern boolean treeReadLenMULT ( FILE *fp, tree *tr, analdef *adef );

extern void getStartingTree ( tree *tr, analdef *adef );
extern double treeLength(tree *tr, int model);

extern void computeBootStopOnly(tree *tr, char *bootStrapFileName, analdef *adef);
extern boolean bootStop(tree *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength);
extern void computeConsensusOnly(tree *tr, char* treeSetFileName, analdef *adef);
extern double evaluatePartialGeneric (tree *, int i, double ki, int _model);
extern double evaluateGeneric (tree *tr, nodeptr p);
extern void newviewGeneric (tree *tr, nodeptr p);
extern void newviewGenericMulti (tree *tr, nodeptr p, int model);
extern void newviewGenericMasked (tree *tr, nodeptr p);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask);
extern void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (tree *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (tree *tr, nodeptr p, int model);
extern double evaluateGenericVector (tree *tr, nodeptr p);
extern void categorizeGeneric (tree *tr, nodeptr p);
extern double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern boolean isTip(int number, int maxTips);
extern void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches);



extern void   newviewIterative(tree *);

extern double evaluateIterative(tree *, boolean writeVector);

extern void *malloc_aligned( size_t size, size_t align);
extern double FABS(double x);




extern void makenewzIterative(tree *);
extern void execCore(tree *, volatile double *dlnLdlz, volatile double *d2lnLdlz2);



extern void determineFullTraversal(nodeptr p, tree *tr);
/*extern void optRateCat(tree *, int i, double lower_spacing, double upper_spacing, double *lhs);*/

extern unsigned int evaluateParsimonyIterative(tree *);
extern void newviewParsimonyIterative(tree *);

extern unsigned int evaluateParsimonyIterativeFast(tree *);
extern void newviewParsimonyIterativeFast(tree *);


extern void initravParsimonyNormal(tree *tr, nodeptr p);

extern double evaluateGenericInitrav (tree *tr, nodeptr p);
extern double evaluateGenericInitravPartition(tree *tr, nodeptr p, int model);
extern void evaluateGenericVectorIterative(tree *, int startIndex, int endIndex);
extern void categorizeIterative(tree *, int startIndex, int endIndex);
extern void onlyInitrav(tree *tr, nodeptr p);
extern void onlyInitravPartition(tree *tr, nodeptr p, int model);

extern void fixModelIndices(tree *tr, int endsite, boolean fixRates);
extern void calculateModelOffsets(tree *tr);
extern void gammaToCat(tree *tr);
extern void catToGamma(tree *tr, analdef *adef);
extern void handleExcludeFile(tree *tr, analdef *adef, rawdata *rdta);
extern void printBaseFrequencies(tree *tr);
extern nodeptr findAnyTip(nodeptr p, int numsp);

extern void parseProteinModel(double *externalAAMatrix, char *fileName);
extern int filexists(char *filename);
extern void computeFullTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches);

extern void computeNextReplicate(tree *tr, long *seed, int *originalRateCategories, int *originalInvariant, boolean isRapid, boolean fixRates);
/*extern void computeNextReplicate(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);*/

extern void reductionCleanup(tree *tr, int *originalRateCategories, int *originalInvariant);
extern void parseSecondaryStructure(tree *tr, analdef *adef, int sites);
extern void printPartitions(tree *tr);
extern void calcDiagptable(double z, int data, int numberOfCategories, double *rptr, double *EIGN, double *diagptable);
extern void compareBips(tree *tr, char *bootStrapFileName, analdef *adef);
extern void computeRF(tree *tr, char *bootStrapFileName, analdef *adef);


extern  unsigned int **initBitVector(tree *tr, unsigned int *vectorLength);
extern hashtable *copyHashTable(hashtable *src, unsigned int vectorLength);
extern hashtable *initHashTable(unsigned int n);
extern void cleanupHashTable(hashtable *h, int state);
extern double convergenceCriterion(hashtable *h, int mxtips);
extern void freeBitVectors(unsigned int **v, int n);
extern void freeHashTable(hashtable *h);
extern stringHashtable *initStringHashTable(hashNumberType n);
extern void addword(char *s, stringHashtable *h, int nodeNumber);


extern void printBothOpen(const char* format, ... );
extern void printBothOpenMPI(const char* format, ... );
extern void initRateMatrix(tree *tr);

extern void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf,
				    int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF);

extern int getIncrement(tree *tr, int model);

extern void fastSearch(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta);
extern void shSupports(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta);

extern FILE *getNumberOfTrees(tree *tr, char *fileName, analdef *adef);

extern void writeBinaryModel(tree *tr);
extern void readBinaryModel(tree *tr);
extern void treeEvaluateRandom (tree *tr, double smoothFactor);
extern void treeEvaluateProgressive(tree *tr);

extern void testGapped(tree *tr);

extern boolean issubset(unsigned int* bipA, unsigned int* bipB, unsigned int vectorLen);
extern boolean compatible(entry* e1, entry* e2, unsigned int bvlen);



extern int *permutationSH(tree *tr, int nBootstrap, long _randomSeed);

extern void updatePerSiteRates(tree *tr, boolean scaleRates);

extern void newviewIterativeAncestral(tree *tr);
extern void newviewGenericAncestral(tree *tr, nodeptr p, boolean atRoot);
extern void computeAncestralStates(tree *tr, double referenceLikelihood, analdef *adef);
extern void makeP_Flex(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, const int numStates);

#ifdef _WAYNE_MPI

extern boolean computeBootStopMPI(tree *tr, char *bootStrapFileName, analdef *adef, double *pearsonAverage);

#endif


#ifdef _USE_PTHREADS

extern size_t getContiguousVectorLength(tree *tr);

extern void makenewzClassify(tree *tr, int maxiter, double *result, double *z0, double *x1_start, double *x2_start,
			     unsigned char *tipX1,  unsigned char *tipX2, int tipCase, boolean *partitionConverged);

extern void    newviewClassify(tree *tr, branchInfo *bInf, double *z);

extern void addTraverseRobIterative(tree *tr, int branchNumber);
extern void insertionsParsimonyIterative(tree *tr, int branchNumber);

extern void newviewClassifySpecial(tree *tr, double *x1_start, double *x2_start, double *x3_start, int *ex1, int *ex2, int *ex3,
				   unsigned char *tipX1,  unsigned char *tipX2, int tipCase, double *pz, double *qz);
extern double evalCL(tree *tr, double *x2, int *ex2, unsigned char *tip, double *pz);

extern void testInsertThoroughIterative(tree *tr, int branchNumber);


/* parallel MRE stuff */

#define MRE_POSSIBLE_CANDIDATE  1
#define MRE_EXCLUDED  2
#define MRE_ADDED  3    
#define SECTION_CONSTANT 1
#define MRE_MIN_AMOUNT_JOBS_PER_THREAD 5


/* work tags for parallel regions */

#define THREAD_NEWVIEW                0
#define THREAD_EVALUATE               1
#define THREAD_MAKENEWZ               2
#define THREAD_MAKENEWZ_FIRST         3
#define THREAD_RATE_CATS              4
#define THREAD_EVALUATE_VECTOR        7
#define THREAD_ALLOC_LIKELIHOOD       8
#define THREAD_COPY_RATE_CATS         9
#define THREAD_COPY_INVAR             10
#define THREAD_COPY_INIT_MODEL        11
#define THREAD_FIX_MODEL_INDICES      12
#define THREAD_INIT_PARTITION         13
#define THREAD_OPT_INVAR              14
#define THREAD_OPT_ALPHA              15
#define THREAD_OPT_RATE               16
#define THREAD_RESET_MODEL            17
#define THREAD_COPY_ALPHA             18
#define THREAD_COPY_RATES             19
#define THREAD_CAT_TO_GAMMA           20
#define THREAD_GAMMA_TO_CAT           21
#define THREAD_NEWVIEW_MASKED         22
#define THREAD_COPY_PARAMS            26
#define THREAD_INIT_EPA                     28
#define THREAD_GATHER_LIKELIHOOD            29
#define THREAD_INSERT_CLASSIFY              30
#define THREAD_INSERT_CLASSIFY_THOROUGH     31
#define THREAD_GATHER_PARSIMONY             32
#define THREAD_PARSIMONY_INSERTIONS         34
#define THREAD_PREPARE_EPA_PARSIMONY        35
#define THREAD_CLEANUP_EPA_PARSIMONY        36
#define THREAD_USE_GAPPED                   38
#define THREAD_PREPARE_BIPS_FOR_PRINT       39
#define THREAD_MRE_COMPUTE                  40
#define THREAD_NEWVIEW_ANCESTRAL            41
#define THREAD_GATHER_ANCESTRAL             42

/*

parallel tree parsing abandoned ...

  #define THREAD_FILL_HASH_FOR_CONSENSUS      41

parallel dropset comp. currently doesn't scale well 

  #define THREAD_FIND_BEST_DROPSET            42
  #define THREAD_CALC_DROPSETS                43

*/

/*
  Pthreads-based MP computations don't really scale for the 
  SSE3-based optimized version 

  #define THREAD_FAST_EVALUATE_PARSIMONY        XX
  #define THREAD_FAST_NEWVIEW_PARSIMONY         XX
  #define THREAD_INIT_FAST_PARSIMONY            XX
*/

typedef struct
{
  tree *tr;
  int threadNumber;
}
  threadData;

void threadMakeVector(tree *tr, int tid);
void threadComputeAverage(tree *tr, int tid);
void threadComputePearson(tree *tr, int tid);
extern void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid);
extern void masterBarrier(int jobType, tree *tr);

#endif



