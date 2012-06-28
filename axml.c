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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models".
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifdef WIN32
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>

#ifdef  _WAYNE_MPI
#include <mpi.h>
#endif



#ifdef _USE_PTHREADS
#include <pthread.h>

#endif

#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
#include <xmmintrin.h>
/*
  special bug fix, enforces denormalized numbers to be flushed to zero,
  without this program is a tiny bit faster though.
  #include <emmintrin.h> 
  #define MM_DAZ_MASK    0x0040
  #define MM_DAZ_ON    0x0040
  #define MM_DAZ_OFF    0x0000
*/
#endif

#include "axml.h"
#include "globalVariables.h"


#define _PORTABLE_PTHREADS


/***************** UTILITY FUNCTIONS **************************/


double FABS(double x)
{
  /*  if(x < -1.0E-10)
      assert(0);*/
  
  /* if(x < 0.0)
     printf("%1.40f\n", x); */

  return fabs(x);
}

void *malloc_aligned(size_t size) 
{
  void *ptr = (void *)NULL;
  int res;
  

#if defined (__APPLE__)
  /* 
     presumably malloc on MACs always returns 
     a 16-byte aligned pointer
  */

  ptr = malloc(size);
  
  if(ptr == (void*)NULL) 
   assert(0);

#else
  res = posix_memalign( &ptr, BYTE_ALIGNMENT, size );

  if(res != 0) 
    assert(0);
#endif 

  /*
    to ensure that the allocated pages are mapped 
    correctly on the distributed shared memory system:

    for(i=0; i<N; i++) 
    // or i+=PAGE_SIZE
    huge[i] = 0.0; // mapping takes place here!
    

   */

   
  return ptr;
}




FILE *getNumberOfTrees(tree *tr, char *fileName, analdef *adef)
{
  FILE 
    *f = myfopen(fileName, "r");

  int 
    trees = 0,
    ch;

  while((ch = fgetc(f)) != EOF)
    if(ch == ';')
      trees++;

  assert(trees > 0);

  tr->numberOfTrees = trees;

  if(!adef->allInOne)   
    printBothOpen("\n\nFound %d trees in File %s\n\n", trees, fileName);


  rewind(f);

  return f;
}

static void printBoth(FILE *f, const char* format, ... )
{
  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);
}

void printBothOpen(const char* format, ... )
{
  FILE *f = myfopen(infoFileName, "ab");

  va_list args;
  va_start(args, format);
  vfprintf(f, format, args );
  va_end(args);

  va_start(args, format);
  vprintf(format, args );
  va_end(args);

  fclose(f);
}

void printBothOpenMPI(const char* format, ... )
{
#ifdef _WAYNE_MPI
  if(processID == 0)
#endif
    {
      FILE *f = myfopen(infoFileName, "ab");

      va_list args;
      va_start(args, format);
      vfprintf(f, format, args );
      va_end(args);
      
      va_start(args, format);
      vprintf(format, args );
      va_end(args);
      
      fclose(f);
    }
}


boolean getSmoothFreqs(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].smoothFrequencies;
}

const unsigned int *getBitVector(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].bitVector;
}


int getStates(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].states;
}

int getUndetermined(int dataType)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return pLengths[dataType].undetermined;
}



char getInverseMeaning(int dataType, unsigned char state)
{
  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  return  pLengths[dataType].inverseMeaning[state];
}

partitionLengths *getPartitionLengths(pInfo *p)
{
  int 
    dataType  = p->dataType,
    states    = p->states,
    tipLength = p->maxTipStates;

  assert(states != -1 && tipLength != -1);

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  pLength.leftLength = pLength.rightLength = states * states;
  pLength.eignLength = states -1;
  pLength.evLength   = states * states;
  pLength.eiLength   = states * states - states;
  pLength.substRatesLength = (states * states - states) / 2;
  pLength.frequenciesLength = states;
  pLength.tipVectorLength   = tipLength * states;
  pLength.symmetryVectorLength = (states * states - states) / 2;
  pLength.frequencyGroupingLength = states;
  pLength.nonGTR = FALSE;

  return (&pLengths[dataType]); 
}



static boolean isCat(analdef *adef)
{
  if(adef->model == M_PROTCAT || adef->model == M_GTRCAT || adef->model == M_BINCAT || adef->model == M_32CAT || adef->model == M_64CAT)
    return TRUE;
  else
    return FALSE;
}

static boolean isGamma(analdef *adef)
{
  if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA || adef->model == M_BINGAMMA || 
     adef->model == M_32GAMMA || adef->model == M_64GAMMA)
    return TRUE;
  else
    return FALSE;
}


static int stateAnalyzer(tree *tr, int model, int maxStates)
{
  boolean
    counter[256],
    previous,
    inputError = FALSE;
  
  int
    lower = tr->partitionData[model].lower,
    upper = tr->partitionData[model].upper,
    i,
    j,
    states = 0;

  for(i = 0; i < 256; i++)
    counter[i] = FALSE;

  for(i = 0; i < tr->rdta->numsp; i++)
    {
      unsigned char *yptr =  &(tr->rdta->y0[((size_t)i) * ((size_t)tr->originalCrunchedLength)]);

      for(j = lower; j < upper; j++)
	if(yptr[j] != getUndetermined(GENERIC_32))
	  counter[yptr[j]] = TRUE;		
    
    }
  
  for(i = 0; i < maxStates; i++)
    {      
      if(counter[i])
	states++;
    }
  

  previous = counter[0];
  
  for(i = 1; i < 256; i++)
    {
      if(previous == FALSE && counter[i] == TRUE)
	{
	  inputError = TRUE;
	  break;
	}     		
      else
	{
	  if(previous == TRUE && counter[i] ==  FALSE)
	    previous = FALSE;
	}
    }
  
  if(inputError)
    {
      printf("Multi State Error, characters must be used in the order they are available, i.e.\n");
      printf("0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V\n");
      printf("You are using the following characters: \n");
      for(i = 0; i < 256; i++)
	if(counter[i])
	  printf("%c ", inverseMeaningGeneric32[i]);
      printf("\n");
      exit(-1);
    }

  return states;
}




static void setRateHetAndDataIncrement(tree *tr, analdef *adef)
{
  int model;

  if(isCat(adef))
    tr->rateHetModel = CAT;
  else
    {
      if(adef->useInvariant)
	tr->rateHetModel = GAMMA_I;
      else
	tr->rateHetModel = GAMMA;
    }

  switch(tr->rateHetModel)
    {
    case GAMMA:
    case GAMMA_I:
      tr->discreteRateCategories = 4;      
      break;
    case CAT:
      if((adef->boot && !adef->bootstrapBranchLengths) || (adef->mode == CLASSIFY_ML) || (tr->catOnly))	
	tr->discreteRateCategories = 1; 	
      else
	tr->discreteRateCategories = 4;
      break;
    default:
      assert(0);
    }

  if(adef->bootstrapBranchLengths)
    assert(tr->discreteRateCategories == 4);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	states = -1,
	maxTipStates = getUndetermined(tr->partitionData[model].dataType) + 1;
      
      switch(tr->partitionData[model].dataType)
	{
	case BINARY_DATA:
	case DNA_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	  states = getStates(tr->partitionData[model].dataType);	 
	  break;	
	case GENERIC_32:
	case GENERIC_64:
	  states = stateAnalyzer(tr, model, getStates(tr->partitionData[model].dataType));	 	 	 	  
	  break;
	default:
	  assert(0);
	}

      tr->partitionData[model].states       = states;
      tr->partitionData[model].maxTipStates = maxTipStates;
    }
}


double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
{
  long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
  double res;

  mult0 = 1549;
  seed0 = *seed & 4095;
  sum  = mult0 * seed0;
  newseed0 = sum & 4095;
  sum >>= 12;
  seed1 = (*seed >> 12) & 4095;
  mult1 =  406;
  sum += mult0 * seed1 + mult1 * seed0;
  newseed1 = sum & 4095;
  sum >>= 12;
  seed2 = (*seed >> 24) & 255;
  sum += mult0 * seed2 + mult1 * seed1;
  newseed2 = sum & 255;

  *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
  res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));

  return res;
}

int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"rb");

  if(fp)
    {
      res = 1;
      fclose(fp);
    }
  else
    res = 0;

  return res;
}


FILE *myfopen(const char *path, const char *mode)
{
  FILE *fp = fopen(path, mode);

  if(strcmp(mode,"r") == 0 || strcmp(mode,"rb") == 0)
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s you want to open for reading does not exist, exiting ...\n", path);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }
  else
    {
      if(fp)
	return fp;
      else
	{
	  if(processID == 0)
	    printf("The file %s RAxML wants to open for writing or appending can not be opened [mode: %s], exiting ...\n",
		   path, mode);
	  errorExit(-1);
	  return (FILE *)NULL;
	}
    }


}





/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/


boolean isTip(int number, int maxTips)
{
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}





void getxsnode (nodeptr p, int model)  
{  
  assert(p->xs[model] || p->next->xs[model] || p->next->next->xs[model]);
  assert(p->xs[model] + p->next->xs[model] + p->next->next->xs[model] == 1);
  
  assert(p == p->next->next->next);

  p->xs[model] = 1;
  
  if(p->next->xs[model])
    {      
      p->next->xs[model] = 0;
      return;
    }
  else
    {
      p->next->next->xs[model] = 0;
      return;
    }  

  assert(0);
}



void getxnode (nodeptr p)
{
  nodeptr  s;

  if ((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }

  assert(p->x);
}





void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;
}


/***********************reading and initializing input ******************/

static void getnums (rawdata *rdta)
{
  if (fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2)
    {
      if(processID == 0)
	printf("ERROR: Problem reading number of species and sites\n");
      errorExit(-1);
    }

  if (rdta->numsp < 4)
    {
      if(processID == 0)
	printf("TOO FEW SPECIES\n");
      errorExit(-1);
    }

  if (rdta->sites < 1)
    {
      if(processID == 0)
	printf("TOO FEW SITES\n");
      errorExit(-1);
    }

  return;
}





boolean whitechar (int ch)
{
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r');
}


static void uppercase (int *chptr)
{
  int  ch;

  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
}




static void getyspace (rawdata *rdta)
{
  size_t size = 4 * ((size_t)(rdta->sites / 4 + 1));
  int    i;
  unsigned char *y0;

  rdta->y = (unsigned char **) malloc((rdta->numsp + 1) * sizeof(unsigned char *));
  assert(rdta->y);   

  y0 = (unsigned char *) malloc(((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));
  assert(y0);   

  rdta->y0 = y0;

  for (i = 0; i <= rdta->numsp; i++)
    {
      rdta->y[i] = y0;
      y0 += size;
    }

  return;
}


static unsigned int KISS32(void)
{
  static unsigned int 
    x = 123456789, 
    y = 362436069,
    z = 21288629,
    w = 14921776,
    c = 0;

  unsigned int t;

  x += 545925293;
  y ^= (y<<13); 
  y ^= (y>>17); 
  y ^= (y<<5);
  t = z + w + c; 
  z = w; 
  c = (t>>31); 
  w = t & 2147483647;

  return (x+y+w);
}

static boolean setupTree (tree *tr, analdef *adef)
{
  nodeptr  p0, p, q;
  int
    i,
    j,
    k,
    tips,
    inter; 
  
  tr->brLenScaler = 1.0;
  tr->storedBrLens = (double*)NULL;

  if(!adef->readTaxaOnly)
    {
      tr->bigCutoff = FALSE;

      tr->patternPosition = (int*)NULL;
      tr->columnPosition = (int*)NULL;

      tr->maxCategories = MAX(4, adef->categories);

      tr->partitionContributions = (double *)malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->partitionContributions[i] = -1.0;

      tr->perPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);
      tr->storedPerPartitionLH = (double *)malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  tr->perPartitionLH[i] = 0.0;
	  tr->storedPerPartitionLH[i] = 0.0;
	}

      if(adef->grouping)
	tr->grouped = TRUE;
      else
	tr->grouped = FALSE;

      if(adef->constraint)
	tr->constrained = TRUE;
      else
	tr->constrained = FALSE;

      tr->treeID = 0;
    }

  tips  = tr->mxtips;
  inter = tr->mxtips - 1;

  if(!adef->readTaxaOnly)
    {
      tr->yVector      = (unsigned char **)  malloc((tr->mxtips + 1) * sizeof(unsigned char *));

      tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double));
      tr->likelihoods  = (double *)malloc(adef->multipleRuns * sizeof(double));
    }

  tr->numberOfTrees = -1;

 

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char)); 

  /*TODO, must that be so long ?*/

  if(!adef->readTaxaOnly)
    {
      
      if(tr->multiGene)
	{
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      tr->td[i].count = 0;
	      tr->td[i].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
	    }
	}
      else
	{
	  tr->td[0].count = 0;
	  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
	}

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->fracchanges[i] = -1.0;
      tr->fracchange = -1.0;

      tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));

      tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));
    }

  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }

  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++)
    {
      p = p0++;

      p->hash   =  KISS32(); /* hast table stuff */
      p->x      =  0;
      p->number =  i;
      p->next   =  p;
      p->back   = (node *)NULL;
      p->bInf   = (branchInfo *)NULL;

      
      for(k = 0; k < NUM_BRANCHES; k++)
	{
	  p->xs[k]    = 0;
	  p->backs[k] = (nodeptr)NULL;
	}

      for(k = 0; k < VECTOR_LENGTH; k++)
	p->isPresent[k] = 0;

      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++)
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++)
	{	 
	  p = p0++;
	  if(j == 1)
	    p->x = 1;
	  else
	    p->x =  0;
	  p->number = i;
	  p->next   = q;
	  p->bInf   = (branchInfo *)NULL;
	  p->back   = (node *) NULL;
	  p->hash   = 0;

	  if(j == 1)
	    for(k = 0; k < NUM_BRANCHES; k++)
	      {
		p->xs[k]    = 1;
		p->backs[k] = (nodeptr)NULL;
	      }
	  else
	    for(k = 0; k < NUM_BRANCHES; k++)
	      {
		p->xs[k]    = 0;
		p->backs[k] = (nodeptr)NULL;
	      }

	  for(k = 0; k < VECTOR_LENGTH; k++)
	    p->isPresent[k] = 0;


	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->startVector[i]  = (node *) NULL;

  tr->ntips       = 0;
  tr->nextnode    = 0;

  if(!adef->readTaxaOnly)
    {
      for(i = 0; i < tr->numBranches; i++)
	tr->partitionSmoothed[i] = FALSE;
    }

  return TRUE;
}


static void checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
    {
      boolean valid;

      switch(buffer[i])
	{
	case '\0':
	case '\t':
	case '\n':
	case '\r':
	case ' ':
	case ':':
	case ',':
	case '(':
	case ')':
	case ';':
	case '[':
	case ']':
	case '\'':
	  valid = FALSE;
	  break;
	default:
	  valid = TRUE;
	}

      if(!valid)
	{
	  printf("ERROR: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n", buffer, i, buffer[i]);
	  printf("Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\", \"\'\" \n");
	  printf("Exiting\n");
	  exit(-1);
	}

    }
  assert(buffer[len - 1] == '\0');
}

static boolean getdata(analdef *adef, rawdata *rdta, tree *tr)
{
  int   
    i, 
    j, 
    basesread, 
    basesnew, 
    ch, my_i, meaning,
    len,
    meaningAA[256], 
    meaningDNA[256], 
    meaningBINARY[256],
    meaningGeneric32[256],
    meaningGeneric64[256];
  
  boolean  
    allread, 
    firstpass;
  
  char 
    buffer[nmlngth + 2];
  
  unsigned char
    genericChars32[32] = {'0', '1', '2', '3', '4', '5', '6', '7', 
			  '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
			  'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
			  'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'};  
  unsigned long 
    total = 0,
    gaps  = 0;

  for (i = 0; i < 256; i++)
    {      
      meaningAA[i]          = -1;
      meaningDNA[i]         = -1;
      meaningBINARY[i]      = -1;
      meaningGeneric32[i]   = -1;
      meaningGeneric64[i]   = -1;
    }

  /* generic 32 data */

  for(i = 0; i < 32; i++)
    meaningGeneric32[genericChars32[i]] = i;
  meaningGeneric32['-'] = getUndetermined(GENERIC_32);
  meaningGeneric32['?'] = getUndetermined(GENERIC_32);

  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20; /* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21; /*21 glutamine glutamic 5 and 6*/

  meaningAA['X'] = 
    meaningAA['?'] = 
    meaningAA['*'] = 
    meaningAA['-'] = 
    getUndetermined(AA_DATA);

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;  
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9; 
  meaningDNA['Y'] = 10;

  meaningDNA['N'] = 
    meaningDNA['O'] = 
    meaningDNA['X'] = 
    meaningDNA['-'] = 
    meaningDNA['?'] = 
    getUndetermined(DNA_DATA);

  /* BINARY DATA */

  meaningBINARY['0'] = 1;
  meaningBINARY['1'] = 2;
  
  meaningBINARY['-'] = 
    meaningBINARY['?'] = 
    getUndetermined(BINARY_DATA);


  /*******************************************************************/

  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread)
    {
      for (i = 1; i <= tr->mxtips; i++)
	{
	  if (firstpass)
	    {
	      ch = getc(INFILE);
	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);

	      my_i = 0;

	      do
		{
		  buffer[my_i] = ch;
		  ch = getc(INFILE);
		  my_i++;
		  if(my_i >= nmlngth)
		    {
		      if(processID == 0)
			{
			  printf("Taxon Name to long at taxon %d, adapt constant nmlngth in\n", i);
			  printf("axml.h, current setting %d\n", nmlngth);
			}
		      errorExit(-1);
		    }
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');

	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);
	      
	      ungetc(ch, INFILE);

	      buffer[my_i] = '\0';
	      len = strlen(buffer) + 1;
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * len);
	      strcpy(tr->nameList[i], buffer);
	    }

	  j = basesread;

	  while ((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
	    {
	      uppercase(& ch);

	      assert(tr->dataVector[j + 1] != -1);

	      switch(tr->dataVector[j + 1])
		{
		case BINARY_DATA:
		  meaning = meaningBINARY[ch];
		  break;
		case DNA_DATA:
		case SECONDARY_DATA:
		case SECONDARY_DATA_6:
		case SECONDARY_DATA_7:
		  /*
		     still dealing with DNA/RNA here, hence just act if as they where DNA characters
		     corresponding column merging for sec struct models will take place later
		  */
		  meaning = meaningDNA[ch];
		  break;
		case AA_DATA:
		  meaning = meaningAA[ch];
		  break;
		case GENERIC_32:
		  meaning = meaningGeneric32[ch];
		  break;
		case GENERIC_64:
		  meaning = meaningGeneric64[ch];
		  break;
		default:
		  assert(0);
		}

	      if (meaning != -1)
		{
		  j++;
		  rdta->y[i][j] = ch;		 
		}
	      else
		{
		  if(!whitechar(ch))
		    {
		      printf("ERROR: Bad base (%c) at site %d of sequence %d\n",
			     ch, j + 1, i);
		      return FALSE;
		    }
		}
	    }

	  if (ch == EOF)
	    {
	      printf("ERROR: End-of-file at site %d of sequence %d\n", j + 1, i);
	      return  FALSE;
	    }

	  if (! firstpass && (j == basesread))
	    i--;
	  else
	    {
	      if (i == 1)
		basesnew = j;
	      else
		if (j != basesnew)
		  {
		    printf("ERROR: Sequences out of alignment\n");
		    printf("%d (instead of %d) residues read in sequence %d %s\n",
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= rdta->sites);
    }

  for(j = 1; j <= tr->mxtips; j++)
    for(i = 1; i <= rdta->sites; i++)
      {
	assert(tr->dataVector[i] != -1);

	switch(tr->dataVector[i])
	  {
	  case BINARY_DATA:
	    meaning = meaningBINARY[rdta->y[j][i]];
	    if(meaning == getUndetermined(BINARY_DATA))
	      gaps++;
	    break;

	  case SECONDARY_DATA:
	  case SECONDARY_DATA_6:
	  case SECONDARY_DATA_7:
	    assert(tr->secondaryStructurePairs[i - 1] != -1);
	    assert(i - 1 == tr->secondaryStructurePairs[tr->secondaryStructurePairs[i - 1]]);
	    /*
	       don't worry too much about undetermined column count here for sec-struct, just count
	       DNA/RNA gaps here and worry about the rest later-on, falling through to DNA again :-)
	    */
	  case DNA_DATA:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == getUndetermined(DNA_DATA))
	      gaps++;
	    break;

	  case AA_DATA:
	    meaning = meaningAA[rdta->y[j][i]];
	    if(meaning == getUndetermined(AA_DATA))
	      gaps++;
	    break;

	  case GENERIC_32:
	    meaning = meaningGeneric32[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_32))
	      gaps++;
	    break;

	  case GENERIC_64:
	    meaning = meaningGeneric64[rdta->y[j][i]];
	    if(meaning == getUndetermined(GENERIC_64))
	      gaps++;
	    break;
	  default:
	    assert(0);
	  }

	total++;
	rdta->y[j][i] = meaning;
      }

  adef->gapyness = (double)gaps / (double)total;

  


  return  TRUE;
}



static void inputweights (rawdata *rdta)
{
  int i, w, fres;
  FILE *weightFile;
  int *wv = (int *)malloc(sizeof(int) *  rdta->sites);

  weightFile = myfopen(weightFileName, "rb");

  i = 0;

  while((fres = fscanf(weightFile,"%d", &w)) != EOF)
    {
      if(!fres)
	{
	  if(processID == 0)
	    printf("error reading weight file probably encountered a non-integer weight value\n");
	  errorExit(-1);
	}
      wv[i] = w;
      i++;
    }

  if(i != rdta->sites)
    {
      if(processID == 0)
	printf("number %d of weights not equal to number %d of alignment columns\n", i, rdta->sites);
      errorExit(-1);
    }

  for(i = 1; i <= rdta->sites; i++)
    rdta->wgt[i] = wv[i - 1];

  fclose(weightFile);
  free(wv);
}



static void getinput(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int i;

  if(!adef->readTaxaOnly)
    {
      INFILE = myfopen(seq_file, "rb");
  
      getnums(rdta);
    }

  tr->mxtips            = rdta->numsp;
  
  if(!adef->readTaxaOnly)
    {
      rdta->wgt             = (int *)    malloc((rdta->sites + 1) * sizeof(int));
      cdta->alias           = (int *)    malloc((rdta->sites + 1) * sizeof(int));
      cdta->aliaswgt        = (int *)    malloc((rdta->sites + 1) * sizeof(int));
      cdta->rateCategory    = (int *)    malloc((rdta->sites + 1) * sizeof(int));
      tr->model             = (int *)    calloc((rdta->sites + 1), sizeof(int));
      tr->initialDataVector  = (int *)    malloc((rdta->sites + 1) * sizeof(int));
      tr->extendedDataVector = (int *)    malloc((rdta->sites + 1) * sizeof(int));     
      cdta->patrat          = (double *) malloc((rdta->sites + 1) * sizeof(double));
      cdta->patratStored    = (double *) malloc((rdta->sites + 1) * sizeof(double));      
      tr->wr                = (double *) malloc((rdta->sites + 1) * sizeof(double)); 
      tr->wr2               = (double *) malloc((rdta->sites + 1) * sizeof(double)); 


      if(!adef->useWeightFile)
	{
	  for (i = 1; i <= rdta->sites; i++)
	    rdta->wgt[i] = 1;
	}
      else
	{
	  assert(!adef->useSecondaryStructure);
	  inputweights(rdta);
	}
    }

  tr->multiBranch = 0;
  tr->numBranches = 1;

  if(!adef->readTaxaOnly)
    {
      if(adef->useMultipleModel)
	{
	  int ref;
	  
	  parsePartitions(adef, rdta, tr);	  
	  
	  for(i = 1; i <= rdta->sites; i++)
	    {
	      ref = tr->model[i];
	      tr->initialDataVector[i] = tr->initialPartitionData[ref].dataType;
	    }
	}
      else
	{
	  int 
	    dataType = -1;
	  
	  tr->initialPartitionData  = (pInfo*)malloc(sizeof(pInfo));
	  tr->initialPartitionData[0].partitionName = (char*)malloc(128 * sizeof(char));
	  strcpy(tr->initialPartitionData[0].partitionName, "No Name Provided");
	  
	  tr->initialPartitionData[0].protModels = adef->proteinMatrix;
	  if(adef->protEmpiricalFreqs)
	    tr->initialPartitionData[0].usePredefinedProtFreqs  = FALSE;
	  else
	    tr->initialPartitionData[0].usePredefinedProtFreqs  = TRUE;
	  
	  

	  tr->NumberOfModels = 1;
	  
	  if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
	    dataType = AA_DATA;
	  if(adef->model == M_GTRCAT || adef->model == M_GTRGAMMA)
	    dataType = DNA_DATA;
	  if(adef->model == M_BINCAT || adef->model == M_BINGAMMA)
	    dataType = BINARY_DATA;
	  if(adef->model == M_32CAT || adef->model == M_32GAMMA)
	    dataType = GENERIC_32;
	  if(adef->model == M_64CAT || adef->model == M_64GAMMA)
	    dataType = GENERIC_64;
	     
	     

	  assert(dataType == BINARY_DATA || dataType == DNA_DATA || dataType == AA_DATA || 
		 dataType == GENERIC_32  || dataType == GENERIC_64);

	  tr->initialPartitionData[0].dataType = dataType;
	  
	  if(dataType == AA_DATA && adef->userProteinModel)
	    {
	      tr->initialPartitionData[0].protModels = PROT_FILE;
	      tr->initialPartitionData[0].usePredefinedProtFreqs  = TRUE;
	      strcpy(tr->initialPartitionData[0].proteinSubstitutionFileName, proteinModelFileName);
	    }
	  
	  for(i = 0; i <= rdta->sites; i++)
	    {
	      tr->initialDataVector[i] = dataType;
	      tr->model[i]      = 0;
	    }
	}

      if(adef->useSecondaryStructure)
	{
	  memcpy(tr->extendedDataVector, tr->initialDataVector, (rdta->sites + 1) * sizeof(int));
	  
	  tr->extendedPartitionData =(pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);
	  
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      tr->extendedPartitionData[i].partitionName = (char*)malloc((strlen(tr->initialPartitionData[i].partitionName) + 1) * sizeof(char));
	      strcpy(tr->extendedPartitionData[i].partitionName, tr->initialPartitionData[i].partitionName);
	      strcpy(tr->extendedPartitionData[i].proteinSubstitutionFileName, tr->initialPartitionData[i].proteinSubstitutionFileName);
	      tr->extendedPartitionData[i].dataType   = tr->initialPartitionData[i].dataType;	      
	      tr->extendedPartitionData[i].protModels = tr->initialPartitionData[i].protModels;
	      tr->extendedPartitionData[i].usePredefinedProtFreqs  = tr->initialPartitionData[i].usePredefinedProtFreqs;
	    }
	  
	  parseSecondaryStructure(tr, adef, rdta->sites);
	  
	  tr->dataVector    = tr->extendedDataVector;
	  tr->partitionData = tr->extendedPartitionData;
	}
      else
	{
	  tr->dataVector    = tr->initialDataVector;
	  tr->partitionData = tr->initialPartitionData;
	}
     
      

      for(i = 0; i < tr->NumberOfModels; i++)
	if(tr->partitionData[i].dataType == AA_DATA && tr->partitionData[i].protModels == PROT_FILE)
	  parseProteinModel(tr->partitionData[i].externalAAMatrix, tr->partitionData[i].proteinSubstitutionFileName);
      
      

      tr->executeModel   = (boolean *)malloc(sizeof(boolean) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->executeModel[i] = TRUE;

      getyspace(rdta);
    } 

  setupTree(tr, adef);


  if(!adef->readTaxaOnly)
    {
      if(!getdata(adef, rdta, tr))
	{
	  printf("Problem reading alignment file \n");
	  errorExit(1);
	}
      tr->nameHash = initStringHashTable(10 * tr->mxtips);
      for(i = 1; i <= tr->mxtips; i++)
	addword(tr->nameList[i], tr->nameHash, i);

      fclose(INFILE);
    }
}



static unsigned char buildStates(int secModel, unsigned char v1, unsigned char v2)
{
  unsigned char new = 0;

  switch(secModel)
    {
    case SECONDARY_DATA:
      new = v1;
      new = new << 4;
      new = new | v2;
      break;
    case SECONDARY_DATA_6:
      {
	int
	  meaningDNA[256],
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};

	const unsigned char
	  finalBinaryStates[6] = {1, 2, 4, 8, 16, 32};

	unsigned char
	  intermediateBinaryStates[6];

	int length = 6;

	for(i = 0; i < 256; i++)
	  meaningDNA[i] = -1;

	meaningDNA['A'] =  1;
	meaningDNA['B'] = 14;
	meaningDNA['C'] =  2;
	meaningDNA['D'] = 13;
	meaningDNA['G'] =  4;
	meaningDNA['H'] = 11;
	meaningDNA['K'] = 12;
	meaningDNA['M'] =  3;
	meaningDNA['N'] = 15;
	meaningDNA['O'] = 15;
	meaningDNA['R'] =  5;
	meaningDNA['S'] =  6;
	meaningDNA['T'] =  8;
	meaningDNA['U'] =  8;
	meaningDNA['V'] =  7;
	meaningDNA['W'] =  9;
	meaningDNA['X'] = 15;
	meaningDNA['Y'] = 10;
	meaningDNA['-'] = 15;
	meaningDNA['?'] = 15;

	for(i = 0; i < length; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];

	    new = n1;
	    new = new << 4;
	    new = new | n2;

	    intermediateBinaryStates[i] = new;
	  }

	new = v1;
	new = new << 4;
	new = new | v2;

	for(i = 0; i < length; i++)
	  {
	    if(new == intermediateBinaryStates[i])
	      break;
	  }
	if(i < length)
	  new = finalBinaryStates[i];
	else
	  {
	    new = 0;
	    for(i = 0; i < length; i++)
	      {
		if(v1 & meaningDNA[allowedStates[i][0]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    new |= finalBinaryStates[i];
		  }
		if(v2 & meaningDNA[allowedStates[i][1]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    new |= finalBinaryStates[i];
		  }
	      }
	  }	
      }
      break;
    case SECONDARY_DATA_7:
      {
	int
	  meaningDNA[256],
	  i;

	const unsigned char
	  allowedStates[6][2] = {{'A','T'}, {'C', 'G'}, {'G', 'C'}, {'G','T'}, {'T', 'A'}, {'T', 'G'}};

	const unsigned char
	  finalBinaryStates[7] = {1, 2, 4, 8, 16, 32, 64};

	unsigned char
	  intermediateBinaryStates[7];

	for(i = 0; i < 256; i++)
	  meaningDNA[i] = -1;

	meaningDNA['A'] =  1;
	meaningDNA['B'] = 14;
	meaningDNA['C'] =  2;
	meaningDNA['D'] = 13;
	meaningDNA['G'] =  4;
	meaningDNA['H'] = 11;
	meaningDNA['K'] = 12;
	meaningDNA['M'] =  3;
	meaningDNA['N'] = 15;
	meaningDNA['O'] = 15;
	meaningDNA['R'] =  5;
	meaningDNA['S'] =  6;
	meaningDNA['T'] =  8;
	meaningDNA['U'] =  8;
	meaningDNA['V'] =  7;
	meaningDNA['W'] =  9;
	meaningDNA['X'] = 15;
	meaningDNA['Y'] = 10;
	meaningDNA['-'] = 15;
	meaningDNA['?'] = 15;
	

	for(i = 0; i < 6; i++)
	  {
	    unsigned char n1 = meaningDNA[allowedStates[i][0]];
	    unsigned char n2 = meaningDNA[allowedStates[i][1]];

	    new = n1;
	    new = new << 4;
	    new = new | n2;

	    intermediateBinaryStates[i] = new;
	  }

	new = v1;
	new = new << 4;
	new = new | v2;

	for(i = 0; i < 6; i++)
	  {
	    /* exact match */
	    if(new == intermediateBinaryStates[i])
	      break;
	  }
	if(i < 6)
	  new = finalBinaryStates[i];
	else
	  {
	    /* distinguish between exact mismatches and partial mismatches */

	    for(i = 0; i < 6; i++)
	      if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		break;
	    if(i < 6)
	      {
		/* printf("partial mismatch\n"); */

		new = 0;
		for(i = 0; i < 6; i++)
		  {
		    if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		      {
			/*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
			new |= finalBinaryStates[i];
		      }
		    else
		      new |=  finalBinaryStates[6];
		  }
	      }
	    else
	      new = finalBinaryStates[6];
	  }	
      }
      break;
    default:
      assert(0);
    }

  return new;

}



static void adaptRdataToSecondary(tree *tr, rawdata *rdta)
{
  int *alias = (int*)calloc(rdta->sites, sizeof(int));
  int i, j, realPosition;  

  for(i = 0; i < rdta->sites; i++)
    alias[i] = -1;

  for(i = 0, realPosition = 0; i < rdta->sites; i++)
    {
      int partner = tr->secondaryStructurePairs[i];
      if(partner != -1)
	{
	  assert(tr->dataVector[i+1] == SECONDARY_DATA || tr->dataVector[i+1] == SECONDARY_DATA_6 || tr->dataVector[i+1] == SECONDARY_DATA_7);

	  if(i < partner)
	    {
	      for(j = 1; j <= rdta->numsp; j++)
		{
		  unsigned char v1 = rdta->y[j][i+1];
		  unsigned char v2 = rdta->y[j][partner+1];

		  assert(i+1 < partner+1);

		  rdta->y[j][i+1] = buildStates(tr->dataVector[i+1], v1, v2);
		}
	      alias[realPosition] = i;
	      realPosition++;
	    }
	}
      else
	{
	  alias[realPosition] = i;
	  realPosition++;
	}
    }

  assert(rdta->sites - realPosition == tr->numberOfSecondaryColumns / 2);

  rdta->sites = realPosition;

  for(i = 0; i < rdta->sites; i++)
    {
      assert(alias[i] != -1);
      tr->model[i+1]    = tr->model[alias[i]+1];
      tr->dataVector[i+1] = tr->dataVector[alias[i]+1];
      rdta->wgt[i+1] =  rdta->wgt[alias[i]+1];

      for(j = 1; j <= rdta->numsp; j++)
	rdta->y[j][i+1] = rdta->y[j][alias[i]+1];
    }

  free(alias);
}

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  gap, i, j, jj, jg, k, n, nsp;
  int  
    *index, 
    *category = (int*)NULL;

  boolean  flip, tied;
  unsigned char  **data;

  if(adef->useSecondaryStructure)
    {
      assert(tr->NumberOfModels > 1 && adef->useMultipleModel);

      adaptRdataToSecondary(tr, rdta);
    }

  if(adef->useMultipleModel)    
    category      = tr->model;
  

  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;


  if(adef->compressPatterns)
    {
      for (gap = n / 2; gap > 0; gap /= 2)
	{
	  for (i = gap + 1; i <= n; i++)
	    {
	      j = i - gap;

	      do
		{
		  jj = index[j];
		  jg = index[j+gap];
		  if(adef->useMultipleModel)
		    {		     		      
		      assert(category[jj] != -1 &&
			     category[jg] != -1);
		     
		      flip = (category[jj] > category[jg]);
		      tied = (category[jj] == category[jg]);		     
		    }
		  else
		    {
		      flip = 0;
		      tied = 1;
		    }

		  for (k = 1; (k <= nsp) && tied; k++)
		    {
		      flip = (data[k][jj] >  data[k][jg]);
		      tied = (data[k][jj] == data[k][jg]);
		    }

		  if (flip)
		    {
		      index[j]     = jg;
		      index[j+gap] = jj;
		      j -= gap;
		    }
		}
	      while (flip && (j > 0));
	    }
	}
    }
}


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  i, sitei, j, sitej, k;
  boolean  tied;
  int 
    *aliasModel = (int*)NULL,
    *aliasSuperModel = (int*)NULL;

  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)malloc(sizeof(int) * (rdta->sites + 1));
      aliasModel      = (int*)malloc(sizeof(int) * (rdta->sites + 1));
    } 

  i = 0;
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;

  if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES)
    {
      int i;

      tr->patternPosition = (int*)malloc(sizeof(int) * rdta->sites);
      tr->columnPosition  = (int*)malloc(sizeof(int) * rdta->sites);

      for(i = 0; i < rdta->sites; i++)
	{
	  tr->patternPosition[i] = -1;
	  tr->columnPosition[i]  = -1;
	}
    }

  

  i = 0;
  for (j = 1; j <= rdta->sites; j++)
    {
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];
      if(!adef->compressPatterns)
	tied = 0;
      else
	{
	  if(adef->useMultipleModel)
	    {	     
	      tied = (tr->model[sitei] == tr->model[sitej]);
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);
	    }
	  else
	    tied = 1;
	}

      for (k = 1; tied && (k <= rdta->numsp); k++)
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);

      if (tied)
	{
	  if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d from column %d also at site %d\n", i, sitei, sitej);*/
	    }


	  cdta->aliaswgt[i] += rdta->wgt[sitej];
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else
	{
	  if (cdta->aliaswgt[i] > 0) i++;

	  if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES)
	    {
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	      /*printf("Pattern %d is from cloumn %d\n", i, sitej);*/
	    }

	  cdta->aliaswgt[i] = rdta->wgt[sitej];
	  cdta->alias[i] = sitej;
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
    }

  cdta->endsite = i;
  if (cdta->aliaswgt[i] > 0) cdta->endsite++;

  if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES)
    {
      for(i = 0; i < rdta->sites; i++)
	{
	  int p  = tr->patternPosition[i];
	  int c  = tr->columnPosition[i];

	  assert(p >= 0 && p < cdta->endsite);
	  assert(c >= 1 && c <= rdta->sites);
	}
    }


  if(adef->useMultipleModel)
    {
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->model[i]      = aliasModel[i];
	  tr->dataVector[i] = aliasSuperModel[i];
	}
    }

  if(adef->useMultipleModel)
    {
      free(aliasModel);
      free(aliasSuperModel);
    }     
}


static boolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int  i;

  for (i = 1; i <= rdta->sites; i++)
    cdta->alias[i] = i;

  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef);

  return TRUE;
}




static boolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  i, j, model, fullSites = 0, modelCounter;

  unsigned char
    *y    = (unsigned char *)malloc(((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char)),
    *yBUF = (unsigned char *)malloc( ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));

  for (i = 1; i <= rdta->numsp; i++)
    for (j = 0; j < cdta->endsite; j++)
      y[(((size_t)(i - 1)) * ((size_t)cdta->endsite)) + j] = rdta->y[i][cdta->alias[j]];

  free(rdta->y0);
  free(rdta->y);

  rdta->y0 = y;
  memcpy(yBUF, y, ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));
  rdta->yBUF = yBUF;

  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->useMultipleModel)
    {
      tr->partitionData[0].lower = 0;

      model        = tr->model[0];
      modelCounter = 0;
     
      i            = 1;

      while(i <  cdta->endsite)
	{	  
	  if(tr->model[i] != model)
	    {	     
	      tr->partitionData[modelCounter].upper     = i;
	      tr->partitionData[modelCounter + 1].lower = i;

	      model = tr->model[i];	     
	      modelCounter++;
	    }
	  i++;
	}

      tr->partitionData[tr->NumberOfModels - 1].upper = cdta->endsite;      
    
      for(i = 0; i < tr->NumberOfModels; i++)		  
	tr->partitionData[i].width      = tr->partitionData[i].upper -  tr->partitionData[i].lower;
	 
      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;
      i            = 1;
	
      while(i < cdta->endsite)
	{	 
	  if(tr->model[i] != model)
	    {
	      model = tr->model[i];
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}      
    }
  else
    {
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = cdta->endsite;
      tr->partitionData[0].width =  tr->partitionData[0].upper -  tr->partitionData[0].lower;
    }

  tr->rdta       = rdta;
  tr->cdta       = cdta;

  tr->invariant          = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalDataVector = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalModel      = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalWeights    = (int *)malloc(cdta->endsite * sizeof(int));

  memcpy(tr->originalModel, tr->model,            cdta->endsite * sizeof(int));
  memcpy(tr->originalDataVector, tr->dataVector,  cdta->endsite * sizeof(int));
  memcpy(tr->originalWeights, tr->cdta->aliaswgt, cdta->endsite * sizeof(int));

  tr->originalCrunchedLength = tr->cdta->endsite;
  for(i = 0; i < tr->cdta->endsite; i++)
    fullSites += tr->cdta->aliaswgt[i];

  tr->fullSites = fullSites;

  for(i = 0; i < rdta->numsp; i++)
    tr->yVector[i + 1] = &(rdta->y0[((size_t)tr->originalCrunchedLength) * ((size_t)i)]);

  return TRUE;
}


static void makeMissingData(tree *tr)
{
  if(tr->multiGene)
    {
      int 
	model, 
	i, 
	j;
      
      double
	totalWidth = 0.0,
	missingWidth = 0.0;
      
      unsigned char 
	undetermined;       

#ifdef _USE_PTHREADS
      assert(0);
#endif

      assert(tr->NumberOfModels > 1 && tr->numBranches > 1);

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 
	    countMissing = 0,
	    width = tr->partitionData[model].upper - tr->partitionData[model].lower;	

	  tr->mxtipsVector[model] = 0;
	  
	  undetermined = getUndetermined(tr->partitionData[model].dataType);	  

	  for(i = 1; i <= tr->mxtips; i++)
	    {
	      unsigned char *tip = tr->partitionData[model].yVector[i];	      
	      

	      assert(width > 0);

	      for(j = 0; j < width; j++)
		if(tip[j] != undetermined)
		  break;

	      if(j == width)				 
		countMissing++;		
	      else
		{
		  tr->nodep[i]->isPresent[model / MASK_LENGTH] |= mask32[model % MASK_LENGTH];
		  if(!tr->startVector[model])
		    {
		      tr->startVector[model] =  tr->nodep[i];
		      /*printf("placing VR into terminal branch %d\n", i);*/
		    }
		}
	    }

	  tr->mxtipsVector[model] = tr->mxtips - countMissing;
	  assert( tr->mxtipsVector[model] + countMissing == tr->mxtips);

	  printBothOpen("Partition %d has %d missing taxa and %d present taxa\n\n", model, countMissing, tr->mxtipsVector[model]);
	  assert(countMissing < tr->mxtips);

	  totalWidth   += (double)(tr->mxtips) * (double)(width);
	  missingWidth += (double)(countMissing) * (double)(width);
	}

      printBothOpen("Percentage of gene-sampling induced gappyness in this alignment: %2.2f%s\n", 100 * (missingWidth / totalWidth), "%");

    }
}






static int sequenceSimilarity(unsigned char *tipJ, unsigned char *tipK, int n)
{
  int i;

  for(i = 0; i < n; i++)
    if(*tipJ++ != *tipK++)
      return 0;

  return 1;
}

static void checkSequences(tree *tr, rawdata *rdta, analdef *adef)
{
  int n = tr->mxtips + 1;
  int i, j;
  int *omissionList     = (int *)calloc(n, sizeof(int));
  int *undeterminedList = (int *)calloc((rdta->sites + 1), sizeof(int));
  int *modelList        = (int *)malloc((rdta->sites + 1)* sizeof(int));
  int count = 0;
  int countNameDuplicates = 0;
  int countUndeterminedColumns = 0;
  int countOnlyGaps = 0;
  int modelCounter = 1;
  unsigned char *tipI, *tipJ;

  for(i = 1; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	if(strcmp(tr->nameList[i], tr->nameList[j]) == 0)
	  {
	    countNameDuplicates++;
	    if(processID == 0)
	      printBothOpen("Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
	  }
    }

  if(countNameDuplicates > 0)
    {
      if(processID == 0)
	printBothOpen("ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
      errorExit(-1);
    }

  if(adef->checkForUndeterminedSequences)
    {
      for(i = 1; i < n; i++)
	{
	  j = 1;
	  
	  while(j <= rdta->sites)
	    {	  
	      if(rdta->y[i][j] != getUndetermined(tr->dataVector[j]))
		break;	  	  
	      
	      j++;
	    }
	  
	  if(j == (rdta->sites + 1))
	    {
	      if(processID == 0)
		printBothOpen("ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",
			      tr->nameList[i]);
	      
	      countOnlyGaps++;
	    }
	}

      if(countOnlyGaps > 0)
	{
	  if(processID == 0)
	    printBothOpen("ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  
	  errorExit(-1);
	}
    }

  for(i = 0; i <= rdta->sites; i++)
    modelList[i] = -1;

  for(i = 1; i <= rdta->sites; i++)
    {
      j = 1;

      while(j < n)
	{
	  if(rdta->y[j][i] != getUndetermined(tr->dataVector[i]))
	    break;

	  
	  j++;
	}

      if(j == n)
	{
	  undeterminedList[i] = 1;

	  if(processID == 0)
	    printBothOpen("IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);

	  countUndeterminedColumns++;
	}
      else
	{
	  if(adef->useMultipleModel)
	    {
	      modelList[modelCounter] = tr->model[i];
	      modelCounter++;
	    }
	}
    }


  for(i = 1; i < n; i++)
    {
      if(omissionList[i] == 0)
	{
	  tipI = &(rdta->y[i][1]);

	  for(j = i + 1; j < n; j++)
	    {
	      if(omissionList[j] == 0)
		{
		  tipJ = &(rdta->y[j][1]);
		  if(sequenceSimilarity(tipI, tipJ, rdta->sites))
		    {
		      if(processID == 0)
			printBothOpen("\n\nIMPORTANT WARNING: Sequences %s and %s are exactly identical\n", tr->nameList[i], tr->nameList[j]);

		      omissionList[j] = 1;
		      count++;
		    }
		}
	    }
	}
    }

  if(count > 0 || countUndeterminedColumns > 0)
    {
      char noDupFile[2048];
      char noDupModels[2048];
      char noDupSecondary[2048];

      if(count > 0 &&processID == 0)
	{
	  printBothOpen("\nIMPORTANT WARNING\n");

	  printBothOpen("Found %d %s that %s exactly identical to other sequences in the alignment.\n", count, (count == 1)?"sequence":"sequences", (count == 1)?"is":"are");

	  printBothOpen("Normally they should be excluded from the analysis.\n\n");
	}

      if(countUndeterminedColumns > 0 && processID == 0)
	{
	  printBothOpen("\nIMPORTANT WARNING\n");

	  printBothOpen("Found %d %s that %s only undetermined values which will be treated as missing data.\n",
			countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");

	  printBothOpen("Normally these columns should be excluded from the analysis.\n\n");
	}

      strcpy(noDupFile, seq_file);
      strcat(noDupFile, ".reduced");

      strcpy(noDupModels, modelFileName);
      strcat(noDupModels, ".reduced");

      strcpy(noDupSecondary, secondaryStructureFileName);
      strcat(noDupSecondary, ".reduced");

      if(processID == 0)
	{
	  if(adef->useSecondaryStructure)
	    {
	      if(countUndeterminedColumns && !filexists(noDupSecondary))
		{
		  FILE *newFile = myfopen(noDupSecondary, "wb");
		  int count;

		  printBothOpen("\nJust in case you might need it, a secondary structure file with \n");
		  printBothOpen("structure assignments for undetermined columns removed is printed to file %s\n",noDupSecondary);

		  for(i = 1, count = 0; i <= rdta->sites; i++)
		    {
		      if(undeterminedList[i] == 0)
			fprintf(newFile, "%c", tr->secondaryStructureInput[i - 1]);
		      else
			count++;
		    }

		  assert(count == countUndeterminedColumns);

		  fprintf(newFile,"\n");

		  fclose(newFile);
		}
	      else
		{
		  if(countUndeterminedColumns)
		    {
		      printBothOpen("\nA secondary structure file with model assignments for undetermined\n");
		      printBothOpen("columns removed has already been printed to  file %s\n",noDupSecondary);
		    }
		}
	    }


	  if(adef->useMultipleModel && !filexists(noDupModels) && countUndeterminedColumns)
	    {
	      FILE *newFile = myfopen(noDupModels, "wb");

	      printBothOpen("\nJust in case you might need it, a mixed model file with \n");
	      printBothOpen("model assignments for undetermined columns removed is printed to file %s\n",noDupModels);

	      for(i = 0; i < tr->NumberOfModels; i++)
		{
		  boolean modelStillExists = FALSE;

		  for(j = 1; (j <= rdta->sites) && (!modelStillExists); j++)
		    {
		      if(modelList[j] == i)
			modelStillExists = TRUE;
		    }

		  if(modelStillExists)
		    {
		      int k = 1;
		      int lower, upper;
		      int parts = 0;


		      switch(tr->partitionData[i].dataType)
			{
			case AA_DATA:
			  {
			    char 
			      AAmodel[1024];

			    if(tr->partitionData[i].protModels != PROT_FILE)
			      {
				strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
				if(tr->partitionData[i].usePredefinedProtFreqs == FALSE)
				  strcat(AAmodel, "F");
				
				fprintf(newFile, "%s, ", AAmodel);
			      }
			    else
			      fprintf(newFile, "[%s], ", tr->partitionData[i].proteinSubstitutionFileName);
			  }
			  break;
			case DNA_DATA:
			  fprintf(newFile, "DNA, ");
			  break;
			case BINARY_DATA:
			  fprintf(newFile, "BIN, ");
			  break;
			case GENERIC_32:
			  fprintf(newFile, "MULTI, ");
			  break;
			case GENERIC_64:
			  fprintf(newFile, "CODON, ");
			  break;
			default:
			  assert(0);
			}

		      fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);

		      while(k <= rdta->sites)
			{
			  if(modelList[k] == i)
			    {
			      lower = k;
			      while((modelList[k + 1] == i) && (k <= rdta->sites))
				k++;
			      upper = k;

			      if(lower == upper)
				{
				  if(parts == 0)
				    fprintf(newFile, "%d", lower);
				  else
				    fprintf(newFile, ",%d", lower);
				}
			      else
				{
				  if(parts == 0)
				    fprintf(newFile, "%d-%d", lower, upper);
				  else
				    fprintf(newFile, ",%d-%d", lower, upper);
				}
			      parts++;
			    }
			  k++;
			}
		      fprintf(newFile, "\n");
		    }
		}
	      fclose(newFile);
	    }
	  else
	    {
	      if(adef->useMultipleModel)
		{
		  printBothOpen("\nA mixed model file with model assignments for undetermined\n");
		  printBothOpen("columns removed has already been printed to  file %s\n",noDupModels);
		}
	    }


	  if(!filexists(noDupFile))
	    {
	      FILE *newFile;

	      printBothOpen("Just in case you might need it, an alignment file with \n");
	      if(count && !countUndeterminedColumns)
		printBothOpen("sequence duplicates removed is printed to file %s\n", noDupFile);
	      if(!count && countUndeterminedColumns)
		printBothOpen("undetermined columns removed is printed to file %s\n", noDupFile);
	      if(count && countUndeterminedColumns)
		printBothOpen("sequence duplicates and undetermined columns removed is printed to file %s\n", noDupFile);

	      newFile = myfopen(noDupFile, "wb");

	      fprintf(newFile, "%d %d\n", tr->mxtips - count, rdta->sites - countUndeterminedColumns);

	      for(i = 1; i < n; i++)
		{
		  if(!omissionList[i])
		    {
		      fprintf(newFile, "%s ", tr->nameList[i]);
		      tipI =  &(rdta->y[i][1]);

		      for(j = 0; j < rdta->sites; j++)
			{
			  if(undeterminedList[j + 1] == 0)			    
			    fprintf(newFile, "%c", getInverseMeaning(tr->dataVector[j + 1], tipI[j]));			      			     			 
			}

		      fprintf(newFile, "\n");
		    }
		}

	      fclose(newFile);
	    }
	  else
	    {
	      if(count && !countUndeterminedColumns)
		printBothOpen("An alignment file with sequence duplicates removed has already\n");
	      if(!count && countUndeterminedColumns)
		printBothOpen("An alignment file with undetermined columns removed has already\n");
	      if(count && countUndeterminedColumns)
		printBothOpen("An alignment file with undetermined columns and sequence duplicates removed has already\n");

	      printBothOpen("been printed to file %s\n",  noDupFile);
	    }
	}
    }

  free(undeterminedList);
  free(omissionList);
  free(modelList);
}







static void generateBS(tree *tr, analdef *adef)
{
  int i, j, k, w;
  int count;
  char outName[1024], buf[16];
  FILE *of;

  assert(adef->boot != 0);

  for(i = 0; i < adef->multipleRuns; i++)
    {
      computeNextReplicate(tr, &adef->boot, (int*)NULL, (int*)NULL, FALSE, FALSE);

      count = 0;
      for(j = 0; j < tr->cdta->endsite; j++)
	count += tr->cdta->aliaswgt[j];

      assert(count == tr->rdta->sites);

      strcpy(outName, workdir);
      strcat(outName, seq_file);
      strcat(outName, ".BS");
      sprintf(buf, "%d", i);
      strcat(outName, buf);
      printf("Printing replicate %d to %s\n", i, outName);

      of = myfopen(outName, "wb");

      fprintf(of, "%d %d\n", tr->mxtips, count);

      for(j = 1; j <= tr->mxtips; j++)
	{
	  unsigned char *tip   =  tr->yVector[tr->nodep[j]->number];
	  fprintf(of, "%s ", tr->nameList[j]);

	  for(k = 0; k < tr->cdta->endsite; k++)
	    {
	      for(w = 0; w < tr->cdta->aliaswgt[k]; w++)
		fprintf(of, "%c", getInverseMeaning(tr->dataVector[k], tip[k]));	      
	    }

	  fprintf(of, "\n");
	}
      fclose(of);
    }
}





static void splitMultiGene(tree *tr, rawdata *rdta)
{
  int i, l;
  int n = rdta->sites + 1;
  int *modelFilter = (int *)malloc(sizeof(int) * n);
  int length, k;
  unsigned char *tip;
  FILE *outf;
  char outFileName[2048];
  char buf[16];

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      strcpy(outFileName, seq_file);
      sprintf(buf, "%d", i);
      strcat(outFileName, ".GENE.");
      strcat(outFileName, buf);
      outf = myfopen(outFileName, "wb");
      length = 0;
      for(k = 1; k < n; k++)
	{
	  if(tr->model[k] == i)
	    {
	      modelFilter[k] = 1;
	      length++;
	    }
	  else
	    modelFilter[k] = -1;
	}

      fprintf(outf, "%d %d\n", rdta->numsp, length);

      for(l = 1; l <= rdta->numsp; l++)
	{
	  fprintf(outf, "%s ", tr->nameList[l]);

	  tip = &(rdta->y[l][0]);

	  for(k = 1; k < n; k++)
	    {
	      if(modelFilter[k] == 1)		
		fprintf(outf, "%c", getInverseMeaning(tr->dataVector[k], tip[k]));		 	     
	    }
	  fprintf(outf, "\n");

	}

      fclose(outf);

      printf("Wrote individual gene/partition alignment to file %s\n", outFileName);
    }

  free(modelFilter);
  printf("Wrote all %d individual gene/partition alignments\n", tr->NumberOfModels);
  printf("Exiting normally\n");
}


static int countTaxaInTopology(void)
{
  FILE *f = myfopen(tree_file, "rb");   

  int
    c,   
    taxaCount = 0;

  while((c = fgetc(f)) != EOF)
    {
      if(c == '(' || c == ',')
	{
	  c = fgetc(f);
	  if(c ==  '(' || c == ',')
	    ungetc(c, f);
	  else
	    {	      	      	  	      
	      do
		{		
		  c = fgetc(f);
		}
	      while(c != ':' && c != ')' && c != ',');	    

	      taxaCount++;	     	     
	    
	      ungetc(c, f);
	    }
	}
    }
 
  printBothOpen("Found a total of %d taxa in tree file %s\n", taxaCount, tree_file);

  fclose(f);

  return taxaCount;
}







static void allocPartitions(tree *tr)
{
  int
    i,
    maxCategories = tr->maxCategories;

 


  for(i = 0; i < tr->NumberOfModels; i++)
    {
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[i]));

      size_t
	k,
	width = tr->partitionData[i].width;      

      tr->partitionData[i].perSiteAAModel = (int *)malloc(sizeof(int) * width);
      for(k = 0; k < width; k++)
	tr->partitionData[i].perSiteAAModel[k] = WAG;
      

      tr->partitionData[i].wr = (double *)malloc(sizeof(double) * width);
      tr->partitionData[i].wr2 = (double *)malloc(sizeof(double) * width);

     

      if(tr->useFastScaling)
	{
	  tr->partitionData[i].globalScaler    = (unsigned int *)calloc(2 * tr->mxtips, sizeof(unsigned int));  	  
	}

      tr->partitionData[i].left              = (double *)malloc_aligned(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[i].right             = (double *)malloc_aligned(pl->rightLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[i].EIGN              = (double*)malloc(pl->eignLength * sizeof(double));
      tr->partitionData[i].EV                = (double*)malloc_aligned(pl->evLength * sizeof(double));
      tr->partitionData[i].EI                = (double*)malloc(pl->eiLength * sizeof(double));
      tr->partitionData[i].substRates        = (double *)malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[i].frequencies       = (double*)malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[i].tipVector         = (double *)malloc_aligned(pl->tipVectorLength * sizeof(double));
      tr->partitionData[i].symmetryVector    = (int *)malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[i].frequencyGrouping = (int *)malloc(pl->frequencyGroupingLength  * sizeof(int));
      tr->partitionData[i].perSiteRates      = (double *)malloc(sizeof(double) * tr->maxCategories);
      tr->partitionData[i].unscaled_perSiteRates = (double *)malloc(sizeof(double) * tr->maxCategories);
      
      
      tr->partitionData[i].nonGTR = FALSE;
       
      

      tr->partitionData[i].gammaRates = (double*)malloc(sizeof(double) * 4);
      tr->partitionData[i].yVector = (unsigned char **)malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

           
      tr->partitionData[i].xVector = (double **)malloc(sizeof(double*) * tr->innerNodes);     
      tr->partitionData[i].xSpaceVector = (size_t *)calloc(tr->innerNodes, sizeof(size_t));	
           
      tr->partitionData[i].expVector      = (int **)malloc(sizeof(int*) * tr->innerNodes);
      tr->partitionData[i].expSpaceVector = (size_t *)calloc(tr->innerNodes, sizeof(size_t));

      tr->partitionData[i].mxtips  = tr->mxtips;

     


#ifndef _USE_PTHREADS    
      {
	int j;

	for(j = 1; j <= tr->mxtips; j++)
	  tr->partitionData[i].yVector[j] = &(tr->yVector[j][tr->partitionData[i].lower]);
      }
#endif

    }
}

#ifndef _USE_PTHREADS





static void allocNodex (tree *tr)
{
  size_t
    i,   
    model,
    offset,
    memoryRequirements = 0;

  allocPartitions(tr);

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;

      int 
	undetermined, 
	j;

      memoryRequirements += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * width;              
	
      tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
      
      tr->partitionData[model].gapVector = (unsigned int*)calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));


      tr->partitionData[model].initialGapVectorSize = tr->partitionData[model].gapVectorLength * 2 * tr->mxtips * sizeof(int);
	
      /* always multiply by 4 due to frequent switching between CAT and GAMMA in standard RAxML */
      
      tr->partitionData[model].gapColumn = (double *)malloc_aligned(((size_t)tr->innerNodes) *
								    ((size_t)4) * 
								    ((size_t)(tr->partitionData[model].states)) *
								    sizeof(double));		  		
	
      undetermined = getUndetermined(tr->partitionData[model].dataType);

      for(j = 1; j <= tr->mxtips; j++)
	for(i = 0; i < width; i++)
	  if(tr->partitionData[model].yVector[j][i] == undetermined)
	    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];      
    }

  tr->perSiteLL       = (double *)malloc((size_t)tr->cdta->endsite * sizeof(double));
  assert(tr->perSiteLL != NULL);

  tr->sumBuffer  = (double *)malloc_aligned(memoryRequirements * sizeof(double));
  assert(tr->sumBuffer != NULL);
 
  offset = 0;

  /* C-OPT for initial testing tr->NumberOfModels will be 1 */

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t 
	lower = tr->partitionData[model].lower,
	width = tr->partitionData[model].upper - lower;

      /* TODO all of this must be reset/adapted when fixModelIndices is called ! */

      
      tr->partitionData[model].sumBuffer       = &tr->sumBuffer[offset];
     

      tr->partitionData[model].perSiteLL    = &tr->perSiteLL[lower];        


      tr->partitionData[model].wgt          = &tr->cdta->aliaswgt[lower];
      tr->partitionData[model].invariant    = &tr->invariant[lower];
      tr->partitionData[model].rateCategory = &tr->cdta->rateCategory[lower];

      offset += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * width;      
    }

  for(i = 0; i < tr->innerNodes; i++)
    {     
      for(model = 0; model < (size_t)tr->NumberOfModels; model++)
	{	 	  
	  tr->partitionData[model].expVector[i] = (int*)NULL;
	  tr->partitionData[model].xVector[i]   = (double*)NULL;		  			      		  		      		      		  	    	    	  	 
	}
    }
}

#endif


static void initAdef(analdef *adef)
{  
  adef->useSecondaryStructure  = FALSE;
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;
  adef->max_rearrange          = 21;
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
  adef->mode                   = BIG_RAPID_MODE;
  adef->categories             = 25;
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE;
  adef->randomStartingTree     = FALSE;
  adef->parsimonySeed          = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;
  adef->outgroup               = FALSE;
  adef->useInvariant           = FALSE;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE;
  adef->allInOne               = FALSE;
  adef->likelihoodTest         = FALSE;
  adef->perGeneBranchLengths   = FALSE;
  adef->generateBS             = FALSE;
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
  adef->compressPatterns       = TRUE; 
  adef->readTaxaOnly           = FALSE;
  adef->meshSearch             = 0;
  adef->useBinaryModelFile     = FALSE;
  adef->leaveDropMode          = FALSE;
  adef->slidingWindowSize      = 100;
  adef->checkForUndeterminedSequences = TRUE;
}




static int modelExists(char *model, analdef *adef)
{
  int i;
  char thisModel[1024];

  /********** BINARY ********************/

   if(strcmp(model, "BINGAMMAI\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }

  if(strcmp(model, "BINGAMMA\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "BINCAT\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "BINCATI\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = TRUE;
      return 1;
    }

  /*********** 32 state ****************************/

  if(strcmp(model, "MULTIGAMMAI\0") == 0)
    {
      adef->model = M_32GAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }

  if(strcmp(model, "MULTIGAMMA\0") == 0)
    {
      adef->model = M_32GAMMA;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "MULTICAT\0") == 0)
    {
      adef->model = M_32CAT;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "MULTICATI\0") == 0)
    {
      adef->model = M_32CAT;
      adef->useInvariant = TRUE;
      return 1;
    }

  /*********** 64 state ****************************/

  if(strcmp(model, "CODONGAMMAI\0") == 0)
    {
      adef->model = M_64GAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }

  if(strcmp(model, "CODONGAMMA\0") == 0)
    {
      adef->model = M_64GAMMA;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "CODONCAT\0") == 0)
    {
      adef->model = M_64CAT;
      adef->useInvariant = FALSE;
      return 1;
    }

  if(strcmp(model, "CODONCATI\0") == 0)
    {
      adef->model = M_64CAT;
      adef->useInvariant = TRUE;
      return 1;
    }


  /*********** DNA **********************/

  if(strcmp(model, "GTRGAMMAI\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }

  if(strcmp(model, "GTRGAMMA\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = FALSE;
      return 1;
    }

  

  if(strcmp(model, "GTRCAT\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useInvariant = FALSE;
      return 1;
    }

   

  if(strcmp(model, "GTRCATI\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useInvariant = TRUE;
      return 1;
    }




  /*************** AA GTR ********************/

  /* TODO empirical FREQS */

  if(strcmp(model, "PROTCATGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      return 1;
    }

  if(strcmp(model, "PROTCATIGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = TRUE;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAGTR\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAIGTR\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = TRUE;
      adef->protEmpiricalFreqs = 1;
      return 1;
    }
  
  /*************** AA GTR_UNLINKED ********************/

  if(strcmp(model, "PROTCATGTR_UNLINKED\0") == 0)    
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");

      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      return 1;
    }

  if(strcmp(model, "PROTCATIGTR_UNLINKED\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = TRUE;
      adef->protEmpiricalFreqs = 1;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAGTR_UNLINKED\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAIGTR_UNLINKED\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = TRUE;
      return 1;
    }
  
  /****************** AA ************************/

  for(i = 0; i < NUM_PROT_MODELS - 2; i++)
    {
      /* check CAT */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  return 1;
	}

      /* check CATF */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  return 1;
	}


      /* check CATI */

      strcpy(thisModel, "PROTCATI");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->useInvariant = TRUE;
	  return 1;
	}

      /* check CATIF */

      strcpy(thisModel, "PROTCATI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}


      /****************check GAMMA ************************/

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->useInvariant = FALSE;
	  return 1;
	}

     


      /*check GAMMAI*/

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->useInvariant = TRUE;
	  return 1;
	}


      /* check GAMMAmodelF */

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = FALSE;
	  return 1;
	}

     
      /* check GAMMAImodelF */

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}

    }

  /*********************************************************************************/



  return 0;
}



static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }

static void checkOutgroups(tree *tr, analdef *adef)
{
  if(adef->outgroup)
    {
      boolean found;
      int i, j;

      for(j = 0; j < tr->numberOfOutgroups; j++)
	{
	  found = FALSE;
	  for(i = 1; (i <= tr->mxtips) && !found; i++)
	    {
	      if(strcmp(tr->nameList[i], tr->outgroups[j]) == 0)
		{
		  tr->outgroupNums[j] = i;
		  found = TRUE;
		}
	    }
	  if(!found)
	    {
	      printf("Error, the outgroup name \"%s\" you specified can not be found in the alignment, exiting ....\n", tr->outgroups[j]);
	      errorExit(-1);
	    }
	}
    }

}

static void parseOutgroups(char *outgr, tree *tr)
{
  int count = 1, i, k;
  char name[nmlngth];

  i = 0;
  while(outgr[i] != '\0')
    {
      if(outgr[i] == ',')
	count++;
      i++;
    }

  tr->numberOfOutgroups = count;

  tr->outgroups = (char **)malloc(sizeof(char *) * count);

  for(i = 0; i < tr->numberOfOutgroups; i++)
    tr->outgroups[i] = (char *)malloc(sizeof(char) * nmlngth);

  tr->outgroupNums = (int *)malloc(sizeof(int) * count);

  i = 0;
  k = 0;
  count = 0;
  while(outgr[i] != '\0')
    {
      if(outgr[i] == ',')
	{
	  name[k] = '\0';
	  strcpy(tr->outgroups[count], name);
	  count++;
	  k = 0;
	}
      else
	{
	  name[k] = outgr[i];
	  k++;
	}
      i++;
    }

  name[k] = '\0';
  strcpy(tr->outgroups[count], name);

  /*for(i = 0; i < tr->numberOfOutgroups; i++)
    printf("%d %s \n", i, tr->outgroups[i]);*/


  /*printf("%s \n", name);*/
}


/*********************************** OUTGROUP STUFF END *********************************************************/


static void printVersionInfo(void)
{
  printf("\n\nThis is %s version %s released by Alexandros Stamatakis in %s.\n\n",  programName, programVersion, programDate);
  printf("With greatly appreciated code contributions by:\n");
  printf("Andre Aberer (HITS)\n");     
  printf("Simon Berger (HITS)\n"); 
  printf("Nick Pattengale (Sandia)\n"); 
  printf("Wayne Pfeiffer (SDSC)\n");
  printf("Akifumi S. Tanabe (Univ. Tsukuba)\n\n");
}

static void printMinusFUsage(void)
{
  printf("\n");
  printf("              \"-f a\": rapid Bootstrap analysis and search for best-scoring ML tree in one program run\n");  

  printf("              \"-f A\": compute marginal ancestral states on a ROOTED reference tree provided with \"t\"\n");

  printf("              \"-f b\": draw bipartition information on a tree provided with \"-t\" based on multiple trees\n");
  printf("                      (e.g., from a bootstrap) in a file specifed by \"-z\"\n");

  printf("              \"-f B\": optimize br-len scaler and other model parameters (GTR, alpha, etc.) on a tree provided with \"-t\".\n");
  printf("                      The tree needs to contain branch lengths. The branch lengths will not be optimized, just scaled by a single common value.\n");


  printf("              \"-f c\": check if the alignment can be properly read by RAxML\n");

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("              \"-f e\": optimize model+branch lengths for given input tree under GAMMA/GAMMAI only\n");

  printf("              \"-f E\": execute very fast experimental tree search, at present only for testing\n");

  printf("              \"-f F\": execute fast experimental tree search, at present only for testing\n");

  printf("              \"-f g\": compute per site log Likelihoods for one ore more trees passed via\n");
  printf("                      \"-z\" and write them to a file that can be read by CONSEL\n");
 
  printf("              \"-f h\": compute log likelihood test (SH-test) between best tree passed via \"-t\"\n");
  printf("                      and a bunch of other trees passed via \"-z\" \n");  

  printf("              \"-f i\": EXPERIMENTAL do not use for real tree inferences: conducts a single cycle of fast lazy SPR moves\n");
  printf("                      on a given input tree, to be used in combination with -C and -M \n");
  
  printf("              \"-f I\": EXPERIMENTAL do not use for real tree inferences: conducts a single cycle of thorough lazy SPR moves\n");
  printf("                      on a given input tree, to be used in combination with -C and -M \n");

  printf("              \"-f j\": generate a bunch of bootstrapped alignment files from an original alignemnt file.\n");
  printf("                      You need to specify a seed with \"-b\" and the number of replicates with \"-#\" \n"); 

  printf("              \"-f J\": Compute SH-like support values on a given tree passed via \"-t\".\n"); 

  printf("              \"-f m\": compare bipartitions between two bunches of trees passed via \"-t\" and \"-z\" \n");
  printf("                      respectively. This will return the Pearson correlation between all bipartitions found\n");
  printf("                      in the two tree files. A file called RAxML_bipartitionFrequencies.outpuFileName\n");
  printf("                      will be printed that contains the pair-wise bipartition frequencies of the two sets\n");

  printf("              \"-f n\": compute the log likelihood score of all trees contained in a tree file provided by\n");
  printf("                      \"-z\" under GAMMA or GAMMA+P-Invar\n");

  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");

  printf("              \"-f p\": perform pure stepwise MP addition of new sequences to an incomplete starting tree and exit\n");

  printf("              \"-f q\": fast quartet calculator\n");

  printf("              \"-f r\": compute pairwise Robinson-Foulds (RF) distances between all pairs of trees in a tree file passed via \"-z\" \n");
  printf("                      if the trees have node labales represented as integer support values the program will also compute two flavors of\n");
  printf("                      the weighted Robinson-Foulds (WRF) distance\n");

  printf("              \"-f s\": split up a multi-gene partitioned alignment into the respective subalignments \n");

  printf("              \"-f S\": compute site-specific placement bias using a leave one out test inspired by the evolutionary placement algorithm\n");

  printf("              \"-f t\": do randomized tree searches on one fixed starting tree\n");

  printf("              \"-f T\": do final thorough optimization of ML tree from rapid bootstrap search in stand-alone mode\n");

  printf("              \"-f u\": execute morphological weight calibration using maximum likelihood, this will return a weight vector.\n");
  printf("                      you need to provide a morphological alignment and a reference tree via \"-t\" \n");    

  printf("              \"-f v\": classify a bunch of environmental sequences into a reference tree using thorough read insertions\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");

  printf("              \"-f w\": compute ELW test on a bunch of trees passed via \"-z\" \n");

  printf("              \"-f x\": compute pair-wise ML distances, ML model parameters will be estimated on an MP \n");
  printf("                      starting tree or a user-defined tree passed via \"-t\", only allowed for GAMMA-based\n");
  printf("                      models of rate heterogeneity\n");

  printf("              \"-f y\": classify a bunch of environmental sequences into a reference tree using parsimony\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");
  
  printf("\n");
  printf("              DEFAULT for \"-f\": new rapid hill climbing\n");

  printf("\n");
}


static void printREADME(void)
{
  printVersionInfo();
  printf("\n");
  printf("Please also consult the RAxML-manual\n");
  printf("\nTo report bugs send an email to stamatak@cs.tum.edu\n");
  printf("Please send me all input files, the exact invocation, details of the HW and operating system,\n");
  printf("as well as all error messages printed to screen.\n\n\n");

  printf("raxmlHPC[-SSE3|-PTHREADS|-PTHREADS-SSE3|-HYBRID|-HYBRID-SSE3]\n");
  printf("      -s sequenceFileName -n outputFileName -m substitutionModel\n");
  printf("      [-a weightFileName] [-A secondaryStructureSubstModel]\n");
  printf("      [-b bootstrapRandomNumberSeed] [-B wcCriterionThreshold]\n");
  printf("      [-c numberOfCategories] [-C] [-d] [-D]\n");
  printf("      [-e likelihoodEpsilon] [-E excludeFileName]\n");
  printf("      [-f a|A|b|B|c|d|e|E|F|g|h|i|I|j|J|m|n|o|p|q|r|s|S|t|T|u|v|w|x|y] [-F]\n");
  printf("      [-g groupingFileName] [-G placementThreshold] [-h]\n");
  printf("      [-i initialRearrangementSetting] [-I autoFC|autoMR|autoMRE|autoMRE_IGN]\n");
  printf("      [-j] [-J MR|MR_DROP|MRE|STRICT|STRICT_DROP] [-k] [-K] [-M]\n");
  printf("      [-o outGroupName1[,outGroupName2[,...]]][-O]\n");
  printf("      [-p parsimonyRandomSeed] [-P proteinModel]\n");
  printf("      [-q multipleModelFileName] [-r binaryConstraintTree]\n");
  printf("      [-R binaryModelParamFile] [-S secondaryStructureFile] [-t userStartingTree]\n");
  printf("      [-T numberOfThreads] [-u] [-U] [-v] [-V] [-w outputDirectory] [-W slidingWindowSize]\n");
  printf("      [-x rapidBootstrapRandomNumberSeed] [-X] [-y]\n");
  printf("      [-z multipleTreesFile] [-#|-N numberOfRuns|autoFC|autoMR|autoMRE|autoMRE_IGN]\n");
  printf("\n");
  printf("      -a      Specify a column weight file name to assign individual weights to each column of \n");
  printf("              the alignment. Those weights must be integers separated by any type and number \n");
  printf("              of whitespaces whithin a separate file, see file \"example_weights\" for an example.\n");
  printf("\n");
  printf("      -A      Specify one of the secondary structure substitution models implemented in RAxML.\n");
  printf("              The same nomenclature as in the PHASE manual is used, available models: \n");
  printf("              S6A, S6B, S6C, S6D, S6E, S7A, S7B, S7C, S7D, S7E, S7F, S16, S16A, S16B\n");
  printf("\n");
  printf("              DEFAULT: 16-state GTR model (S16)\n");
  printf("\n");
  printf("      -b      Specify an integer number (random seed) and turn on bootstrapping\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -B      specify a floating point number between 0.0 and 1.0 that will be used as cutoff threshold \n");
  printf("              for the MR-based bootstopping criteria. The recommended setting is 0.03.\n");
  printf("\n");
  printf("              DEFAULT: 0.03 (recommended empirically determined setting)\n");
  printf("\n");
  printf("      -c      Specify number of distinct rate catgories for RAxML when modelOfEvolution\n");
  printf("              is set to GTRCAT or GTRMIX\n");
  printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
  printf("              categories to accelerate computations. \n");
  printf("\n");
  printf("              DEFAULT: 25\n");
  printf("\n");
  printf("      -C      Conduct model parameter optimization on gappy, partitioned multi-gene alignments with per-partition\n");
  printf("              branch length estimates (-M enabled) using the fast method with pointer meshes described in:\n");
  printf("              Stamatakis and Ott: \"Efficient computation of the phylogenetic likelihood function on multi-gene alignments and multi-core processors\"\n");
  printf("              WARNING: We can not conduct useful tree searches using this method yet! Does not work with Pthreads version.\n");
  printf("\n");
  printf("      -d      start ML optimization from random starting tree \n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -D      ML search convergence criterion. This will break off ML searches if the relative \n");
  printf("              Robinson-Foulds distance between the trees obtained from two consecutive lazy SPR cycles\n");
  printf("              is smaller or equal to 1%s. Usage recommended for very large datasets in terms of taxa.\n", "%");
  printf("              On trees with more than 500 taxa this will yield execution time improvements of approximately 50%s\n",  "%");
  printf("              While yielding only slightly worse trees.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -e      set model optimization precision in log likelihood units for final\n");
  printf("              optimization of tree topology under MIX/MIXI or GAMMA/GAMMAI\n");
  printf("\n");
  printf("              DEFAULT: 0.1   for models not using proportion of invariant sites estimate\n");
  printf("                       0.001 for models using proportion of invariant sites estimate\n");
  printf("\n");
  printf("      -E      specify an exclude file name, that contains a specification of alignment positions you wish to exclude.\n");
  printf("              Format is similar to Nexus, the file shall contain entries like \"100-200 300-400\", to exclude a\n");
  printf("              single column write, e.g., \"100-100\", if you use a mixed model, an appropriatly adapted model file\n");
  printf("              will be written.\n");
  printf("\n");
  printf("      -f      select algorithm:\n");

  printMinusFUsage();

  printf("\n");
  printf("      -F      enable ML tree searches under CAT model for very large trees without switching to \n");
  printf("              GAMMA in the end (saves memory).\n");
  printf("              This option can also be used with the GAMMA models in order to avoid the thorough optimization \n");
  printf("              of the best-scoring ML tree in the end.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -g      specify the file name of a multifurcating constraint tree\n");
  printf("              this tree does not need to be comprehensive, i.e. must not contain all taxa\n");
  printf("\n");
  printf("      -G      enable the ML-based evolutionary placement algorithm heuristics\n");
  printf("              by specifiyng a threshold value (fraction of insertion branches to be evaluated\n");
  printf("              using slow insertions under ML).\n");
  printf("\n");
  printf("      -h      Display this help message.\n");
  printf("\n");
  printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
  printf("              changes phase\n");
  printf("\n");
  printf("      -I      a posteriori bootstopping analysis. Use:\n");
  printf("             \"-I autoFC\" for the frequency-based criterion\n");
  printf("             \"-I autoMR\" for the majority-rule consensus tree criterion\n");
  printf("             \"-I autoMRE\" for the extended majority-rule consensus tree criterion\n");
  printf("             \"-I autoMRE_IGN\" for metrics similar to MRE, but include bipartitions under the threshold whether they are compatible\n");
  printf("                              or not. This emulates MRE but is faster to compute.\n");
  printf("              You also need to pass a tree file containg several bootstrap replicates via \"-z\" \n"); 
  printf("\n");
  printf("      -j      Specifies that intermediate tree files shall be written to file during the standard ML and BS tree searches.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -J      Compute majority rule consensus tree with \"-J MR\" or extended majority rule consensus tree with \"-J MRE\"\n");
  printf("              or strict consensus tree with \"-J STRICT\".\n");
  printf("              Options \"-J STRICT_DROP\" and \"-J MR_DROP\" will execute an algorithm that identifies dropsets which contain\n");
  printf("              rogue taxa as proposed by Pattengale et al. in the paper \"Uncovering hidden phylogenetic consensus\".\n");
  printf("              You will also need to provide a tree file containing several UNROOTED trees via \"-z\"\n");
  printf("\n");
  printf("      -k      Specifies that bootstrapped trees should be printed with branch lengths.\n");
  printf("              The bootstraps will run a bit longer, because model parameters will be optimized\n");
  printf("              at the end of each run under GAMMA or GAMMA+P-Invar respectively.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");  
  printf("      -K      Specify one of the multi-state substitution models (max 32 states) implemented in RAxML.\n");
  printf("              Available models are: ORDERED, MK, GTR\n");
  printf("\n");
  printf("              DEFAULT: GTR model \n");
  printf("\n");
  printf("      -m      Model of Binary (Morphological), Nucleotide, Multi-State, or Amino Acid Substitution: \n");
  printf("\n");
  printf("              BINARY:\n\n");
  printf("                \"-m BINCAT\"         : Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                      automatically under BINGAMMA, depending on the tree search option\n");
  printf("                \"-m BINCATI\"        : Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                      automatically under BINGAMMAI, depending on the tree search option \n");
  printf("                \"-m BINGAMMA\"       : GAMMA model of rate \n");
  printf("                                      heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m BINGAMMAI\"      : Same as BINGAMMA, but with estimate of proportion of invariable sites\n");
  printf("\n");
  printf("              NUCLEOTIDES:\n\n");
  printf("                \"-m GTRCAT\"         : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                      under GTRGAMMA, depending on the tree search option\n");  
  printf("                \"-m GTRCATI\"        : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                      under GTRGAMMAI, depending on the tree search option\n");
  printf("                \"-m GTRGAMMA\"       : GTR + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                      heterogeneity (alpha parameter will be estimated)\n");  
  printf("                \"-m GTRGAMMAI\"      : Same as GTRGAMMA, but with estimate of proportion of invariable sites \n");
  printf("\n");
  printf("              MULTI-STATE:\n\n");
  printf("                \"-m MULTICAT\"         : Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                      automatically under MULTIGAMMA, depending on the tree search option\n");
  printf("                \"-m MULTICATI\"        : Optimization of site-specific\n");
  printf("                                      evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                      rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                      automatically under MULTIGAMMAI, depending on the tree search option \n");
  printf("                \"-m MULTIGAMMA\"       : GAMMA model of rate \n");
  printf("                                      heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m MULTIGAMMAI\"      : Same as MULTIGAMMA, but with estimate of proportion of invariable sites\n");
  printf("\n");
  printf("                You can use up to 32 distinct character states to encode multi-state regions, they must be used in the following order:\n");
  printf("                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V\n");
  printf("                i.e., if you have 6 distinct character states you would use 0, 1, 2, 3, 4, 5 to encode these.\n");
  printf("                The substitution model for the multi-state regions can be selected via the \"-K\" option\n");
  printf("\n");
  printf("              AMINO ACIDS:\n\n");
  printf("                \"-m PROTCATmatrixName[F]\"         : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                    evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                    rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                    automatically under PROTGAMMAmatrixName[f], depending on the tree search option\n");  
  printf("                \"-m PROTCATImatrixName[F]\"        : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                    evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                    rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                    automatically under PROTGAMMAImatrixName[f], depending on the tree search option\n");
  printf("                \"-m PROTGAMMAmatrixName[F]\"       : specified AA matrix + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                                    heterogeneity (alpha parameter will be estimated)\n");  
  printf("                \"-m PROTGAMMAImatrixName[F]\"      : Same as PROTGAMMAmatrixName[F], but with estimate of proportion of invariable sites \n");
  printf("\n");
  printf("                Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, DUMMY, DUMMY2, GTR_UNLINKED, GTR\n");
  printf("                With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
  printf("                Please note that for mixed models you can in addition specify the per-gene AA model in\n");
  printf("                the mixed model file (see manual for details). Also note that if you estimate AA GTR parameters on a partitioned\n");
  printf("                dataset, they will be linked (estimated jointly) across all partitions to avoid over-parametrization\n");
  printf("\n");
  printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");
  printf("              Branch lengths for individual partitions will be printed to separate files\n");
  printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
  printf("\n"),
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -n      Specifies the name of the output file.\n");
  printf("\n");
  printf("      -o      Specify the name of a single outgrpoup or a comma-separated list of outgroups, eg \"-o Rat\" \n");
  printf("              or \"-o Rat,Mouse\", in case that multiple outgroups are not monophyletic the first name \n");
  printf("              in the list will be selected as outgroup, don't leave spaces between taxon names!\n"); 
  printf("\n");
  printf("      -O      Disable check for completely undetermined sequence in alignment.\n");
  printf("              The program will not exit with an error message when \"-O\" is specified.\n");
  printf("\n");
  printf("              DEFAULT: check enabled\n");
  printf("\n");
  printf("      -p      Specify a random number seed for the parsimony inferences. This allows you to reproduce your results\n");
  printf("              and will help me debug the program.\n");
  printf("\n");
  printf("      -P      Specify the file name of a user-defined AA (Protein) substitution model. This file must contain\n");
  printf("              420 entries, the first 400 being the AA substitution rates (this must be a symmetric matrix) and the\n");
  printf("              last 20 are the empirical base frequencies\n");
  printf("\n");
  printf("      -q      Specify the file name which contains the assignment of models to alignment\n");
  printf("              partitions for multiple models of substitution. For the syntax of this file\n");
  printf("              please consult the manual.\n"); 
  printf("\n");
  printf("      -r      Specify the file name of a binary constraint tree.\n");
  printf("              this tree does not need to be comprehensive, i.e. must not contain all taxa\n");
  printf("\n");
  printf("      -R      Specify the file name of a binary model parameter file that has previously been generated\n");
  printf("              with RAxML using the -f e tree evaluation option. The file name should be: \n");
  printf("              RAxML_binaryModelParameters.runID\n");
  printf("\n");
  printf("      -s      Specify the name of the alignment data file in PHYLIP format\n");
  printf("\n");
  printf("      -S      Specify the name of a secondary structure file. The file can contain \".\" for \n");
  printf("              alignment columns that do not form part of a stem and characters \"()<>[]{}\" to define \n");
  printf("              stem regions and pseudoknots\n");
  printf("\n");
  printf("      -t      Specify a user starting tree file name in Newick format\n");
  printf("\n");
  printf("      -T      PTHREADS VERSION ONLY! Specify the number of threads you want to run.\n");
  printf("              Make sure to set \"-T\" to at most the number of CPUs you have on your machine,\n");
  printf("              otherwise, there will be a huge performance decrease!\n");
  printf("\n");
  printf("      -u      use the median for the discrete approximation of the GAMMA model of rate heterogeneity\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -U      Try to save memory by using SEV-based implementation for gap columns on large gappy alignments\n");
  printf("              The technique is described here: http://www.biomedcentral.com/1471-2105/12/470\n");
  printf("              This will only work for DNA and/or PROTEIN data and only with the SSE3 or AVX-vextorized version of the code.\n");
  printf("\n");
  printf("      -v      Display version information\n");
  printf("\n");
  printf("      -V      Disable rate heterogeneity among sites model and use one without rate heterogeneity instead.\n");
  printf("              Only works if you specify the CAT model of rate heterogeneity.\n");
  printf("\n");
  printf("              DEFAULT: use rate heterogeneity\n");
  printf("\n");
  printf("      -w      FULL (!) path to the directory into which RAxML shall write its output files\n");
  printf("\n");
  printf("              DEFAULT: current directory\n");
  printf("\n");
  printf("      -W      Sliding window size for leave-one-out site-specific placement bias algorithm\n");
  printf("              only effective when used in combination with \"-f S\" \n");
  printf("\n");
  printf("              DEFAULT: 100 sites\n");
  printf("\n");
  printf("      -x      Specify an integer number (random seed) and turn on rapid bootstrapping\n");
  printf("              CAUTION: unlike in version 7.0.4 RAxML will conduct rapid BS replicates under \n");
  printf("              the model of rate heterogeneity you specified via \"-m\" and not by default under CAT\n");
  printf("\n");
  printf("      -X      EXPERIMENTAL OPTION: This option will do a per-site estimate of protein substitution models\n");
  printf("              by looping over all given, fixed models LG, WAG, JTT, etc and using their respective base frequencies to independently\n");
  printf("              assign a prot subst. model to each site via ML optimization\n");
  printf("              At present this option only works with the GTR+GAMMA model, unpartitioned datasets, and in the sequential\n");
  printf("              version only.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -y      If you want to only compute a parsimony starting tree with RAxML specify \"-y\",\n");
  printf("              the program will exit after computation of the starting tree\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n"); 
  printf("      -z      Specify the file name of a file containing multiple trees e.g. from a bootstrap\n");
  printf("              that shall be used to draw bipartition values onto a tree provided with \"-t\",\n");
  printf("              It can also be used to compute per site log likelihoods in combination with \"-f g\"\n");
  printf("              and to read a bunch of trees for a couple of other options (\"-f h\", \"-f m\", \"-f n\").\n");
  printf("\n");
  printf("      -#|-N   Specify the number of alternative runs on distinct starting trees\n");
  printf("              In combination with the \"-b\" option, this will invoke a multiple boostrap analysis\n");
  printf("              Note that \"-N\" has been added as an alternative since \"-#\" sometimes caused problems\n");
  printf("              with certain MPI job submission systems, since \"-#\" is often used to start comments.\n");
  printf("              If you want to use the bootstopping criteria specify \"-# autoMR\" or \"-# autoMRE\" or \"-# autoMRE_IGN\"\n");
  printf("              for the majority-rule tree based criteria (see -I option) or \"-# autoFC\" for the frequency-based criterion.\n");
  printf("              Bootstopping will only work in combination with \"-x\" or \"-b\"\n");
  printf("\n");
  printf("              DEFAULT: 1 single analysis\n");
  printf("\n\n\n\n");

}




static void analyzeRunId(char id[128])
{
  int i = 0;

  while(id[i] != '\0')
    {    
      if(i >= 128)
	{
	  printf("Error: run id after \"-n\" is too long, it has %d characters please use a shorter one\n", i);
	  assert(0);
	}
      
      if(id[i] == '/')
	{
	  printf("Error character %c not allowed in run ID\n", id[i]);
	  assert(0);
	}


      i++;
    }

  if(i == 0)
    {
      printf("Error: please provide a string for the run id after \"-n\" \n");
      assert(0);
    }

}

static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  boolean
    bad_opt    =FALSE,
    resultDirSet = FALSE;

  char
    resultDir[1024] = "",
    aut[256],         
    *optarg,
    model[2048] = "",
    secondaryModel[2048] = "",
    multiStateModel[2048] = "",
    modelChar;

  double 
    likelihoodEpsilon,    
    wcThreshold,
    fastEPAthreshold;
  
  int  
    optind = 1,        
    c,
    nameSet = 0,
    alignmentSet = 0,
    multipleRuns = 0,
    constraintSet = 0,
    treeSet = 0,
    groupSet = 0,
    modelSet = 0,
    treesSet  = 0;

  boolean
    bSeedSet = FALSE,
    xSeedSet = FALSE,
    multipleRunsSet = FALSE;

  run_id[0] = 0;
  workdir[0] = 0;
  seq_file[0] = 0;
  tree_file[0] = 0;
  model[0] = 0;
  weightFileName[0] = 0;
  modelFileName[0] = 0;

  /*********** tr inits **************/

#ifdef _USE_PTHREADS
  NumberOfThreads = 0;
#endif
  

  tr->useFastScaling = TRUE;

 
  tr->bootStopCriterion = -1;
  tr->wcThreshold = 0.03;
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->catOnly = FALSE;
  tr->multiGene = 0;
  tr->useEpaHeuristics = FALSE;
  tr->fastEPAthreshold = -1.0;
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->useGappedImplementation = FALSE;
  tr->saveMemory = FALSE;
  tr->estimatePerSiteAA = FALSE;
  tr->useGammaMedian = FALSE;
  tr->noRateHet = FALSE;
  
  /********* tr inits end*************/


  while(!bad_opt &&
	((c = mygetopt(argc,argv,"R:T:E:N:B:L:P:S:A:G:H:I:J:K:W:l:x:z:g:r:e:a:b:c:f:i:m:t:w:s:n:o:q:#:p:vudyjhkMDFCQUXOV", &optind, &optarg))!=-1))
    {
    switch(c)
      {
      case 'V':
	tr->noRateHet = TRUE;
	break;
      case 'u':
	tr->useGammaMedian = TRUE;
	break;
      case 'O':
	adef->checkForUndeterminedSequences = FALSE;
	break;
      case 'X':
	tr->estimatePerSiteAA = TRUE;
	tr->useFastScaling    = FALSE;
	break;
      case 'W':
	sscanf(optarg,"%d", &(adef->slidingWindowSize));
	if(adef->slidingWindowSize <= 0)
	  {
	    printf("You can't use a sliding window size smaller than 1, you specified %d\n", adef->slidingWindowSize);
	    exit(-1);
	  }
	if(adef->slidingWindowSize <= 10)	  
	  {
	    printf("You specified a small sliding window size of %d sites\n", adef->slidingWindowSize);	
	    printf("Are you sure you want to do this?\n");
	  }
	if(adef->slidingWindowSize >= 500)	  
	  {
	    printf("You specified a large sliding window size of %d sites\n", adef->slidingWindowSize);	
	    printf("Are you sure you want to do this?\n");
	  }
	break;
      case 'U':
	tr->saveMemory = TRUE;
#if (!defined(__SIM_SSE3) && !defined(__AVX))	
	printf("\nmemory saving option -U does only work with the AVX and SSE3 vectorized versions of the code\n");
	printf("please remove this option and execute the program again\n");
	printf("exiting ....\n\n");
	errorExit(0);
#endif
	break;
      case 'R':
	adef->useBinaryModelFile = TRUE;
	strcpy(binaryModelParamsInputFileName, optarg);
	break;          
      case 'K':
	{
	  const char *modelList[3] = { "ORDERED", "MK", "GTR"};
	  const int states[3] = {ORDERED_MULTI_STATE, MK_MULTI_STATE, GTR_MULTI_STATE};
	  int i;

	  sscanf(optarg, "%s", multiStateModel);

	  for(i = 0; i < 3; i++)
	    if(strcmp(multiStateModel, modelList[i]) == 0)
	      break;

	  if(i < 3)
	    tr->multiStateModel = states[i];
	  else
	    {
	      printf("The multi-state model %s you want to use does not exist, exiting .... \n", multiStateModel);
	      errorExit(0);
	    }
	  

	}
	break;
      case 'A':
	{
	  const char *modelList[21] = { "S6A", "S6B", "S6C", "S6D", "S6E", "S7A", "S7B", "S7C", "S7D", "S7E", "S7F", "S16", "S16A", "S16B", "S16C",
				      "S16D", "S16E", "S16F", "S16I", "S16J", "S16K"};
	  int i;

	  sscanf(optarg, "%s", secondaryModel);

	  for(i = 0; i < 21; i++)
	    if(strcmp(secondaryModel, modelList[i]) == 0)
	      break;

	  if(i < 21)
	    tr->secondaryStructureModel = i;
	  else
	    {
	      printf("The secondary structure model %s you want to use does not exist, exiting .... \n", secondaryModel);
	      errorExit(0);
	    }
	}
	break;
      case 'B':
	sscanf(optarg,"%lf", &wcThreshold);
	tr->wcThreshold = wcThreshold;
	if(wcThreshold <= 0.0 || wcThreshold >= 1.0)
	  {
	    printf("\nBootstrap threshold must be set to values between 0.0 and 1.0, you just set it to %f\n", wcThreshold);
	    exit(-1);
	  }
	if(wcThreshold < 0.01 || wcThreshold > 0.05)
	  {
	    printf("\n\nWARNING, reasonable settings for Bootstopping threshold with MR-based criteria range between 0.01 and 0.05.\n");
	    printf("You are just setting it to %f, the most reasonable empirically determined setting is 0.03 \n\n", wcThreshold);
	  }
	break;
      case 'C':
	tr->multiGene = 1;
	break;
      case 'D':
	tr->searchConvergenceCriterion = TRUE;
	break;
      case 'E':
	strcpy(excludeFileName, optarg);
	adef->useExcludeFile = TRUE;
	break;
      case 'F':
	tr->catOnly = TRUE;
	break;
      case 'G':
	tr->useEpaHeuristics = TRUE;
	
	sscanf(optarg,"%lf", &fastEPAthreshold);
	tr->fastEPAthreshold = fastEPAthreshold;
	
	if(fastEPAthreshold <= 0.0 || fastEPAthreshold >= 1.0)
	  {
	    printf("\nHeuristic EPA threshold must be set to values between 0.0 and 1.0, you just set it to %f\n", fastEPAthreshold);
	    exit(-1);
	  }
	if(fastEPAthreshold < 0.015625 || fastEPAthreshold > 0.5)
	  {
	    printf("\n\nWARNING, reasonable settings for heuristic EPA threshold range between 0.015625 (1/64) and 0.5 (1/2).\n");
	    printf("You are just setting it to %f\n\n", fastEPAthreshold);
	  }	
#ifdef _USE_PTHREADS
	tr->useFastScaling = FALSE;
#endif	
	break;	
     
      case 'I':     
	adef->readTaxaOnly = TRUE;
	adef->mode = BOOTSTOP_ONLY;
	if((sscanf(optarg,"%s", aut) > 0) && ((strcmp(aut, "autoFC") == 0) || (strcmp(aut, "autoMR") == 0) || 
					      (strcmp(aut, "autoMRE") == 0) || (strcmp(aut, "autoMRE_IGN") == 0)))
	  {
	    if((strcmp(aut, "autoFC") == 0))	   
	      tr->bootStopCriterion = FREQUENCY_STOP;
	    if((strcmp(aut, "autoMR") == 0))		  	    
	      tr->bootStopCriterion = MR_STOP;	   
	    if((strcmp(aut, "autoMRE") == 0))	   
	      tr->bootStopCriterion = MRE_STOP;
	    if((strcmp(aut, "autoMRE_IGN") == 0))
	      tr->bootStopCriterion = MRE_IGN_STOP;
	  }
	else
	  {
	    if(processID == 0)	      
	      printf("Use -I a posteriori bootstop option either as \"-I autoFC\" or \"-I autoMR\" or \"-I autoMRE\" or \"-I autoMRE_IGN\"\n");	       	      
	    errorExit(0);
	  }
	break;     
      case 'J':	
	adef->readTaxaOnly = TRUE;
	adef->mode = CONSENSUS_ONLY;
	
	if((sscanf(optarg,"%s", aut) > 0) && ((strcmp(aut, "MR") == 0) || (strcmp(aut, "MRE") == 0) || (strcmp(aut, "STRICT") == 0) || 
					      (strcmp(aut, "STRICT_DROP") == 0) || (strcmp(aut, "MR_DROP") == 0)))
	  {
	    if((strcmp(aut, "MR") == 0))	   
	      tr->consensusType = MR_CONSENSUS;
	    if((strcmp(aut, "MR_DROP") == 0))	   
	      {
		tr->consensusType = MR_CONSENSUS;
		adef->leaveDropMode = TRUE;
	      }

	    if((strcmp(aut, "MRE") == 0))		  	    
	      tr->consensusType = MRE_CONSENSUS;

	    
	    if((strcmp(aut, "STRICT") == 0))		  	    
	      tr->consensusType = STRICT_CONSENSUS;	
	    if((strcmp(aut, "STRICT_DROP") == 0))		
	      {
		tr->consensusType = STRICT_CONSENSUS; 
		adef->leaveDropMode = TRUE;
	      }
	  }
	else
	  {
	    if(processID == 0)	      
	      printf("Use -J consensus tree option either as \"-J MR\" or \"-J MRE\" or \"-J STRICT\" or \"-J MR_DROP\"  or \"-J STRICT_DROP\"\n");	       	      
	    errorExit(0);
	  }	
	      
	break;
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;
      case 'P':
	strcpy(proteinModelFileName, optarg);
	adef->userProteinModel = TRUE;
	/*parseProteinModel(adef->externalAAMatrix, proteinModelFileName);*/
	break;
      case 'S':
	adef->useSecondaryStructure = TRUE;
	strcpy(secondaryStructureFileName, optarg);
	break;
      case 'T':
#ifdef _USE_PTHREADS
	sscanf(optarg,"%d", &NumberOfThreads);
#else
	if(processID == 0)
	  {
	    printf("Option -T does not have any effect with the sequential or parallel MPI version.\n");
	    printf("It is used to specify the number of threads for the Pthreads-based parallelization\n");
	  }
#endif
	break;                  
      case 'o':
	{
	  char *outgroups;
	  outgroups = (char*)malloc(sizeof(char) * (strlen(optarg) + 1));
	  strcpy(outgroups, optarg);
	  parseOutgroups(outgroups, tr);
	  free(outgroups);
	  adef->outgroup = TRUE;
	}
	break;
      case 'k':
	adef->bootstrapBranchLengths = TRUE;
	break;
      case 'z':
	strcpy(bootStrapFile, optarg);
	treesSet = 1;
	break;
      case 'd':
	adef->randomStartingTree = TRUE;
	break;
      case 'g':
	strcpy(tree_file, optarg);
	adef->grouping = TRUE;
	adef->restart  = TRUE;
	groupSet = 1;
	break;
      case 'r':
	strcpy(tree_file, optarg);
	adef->restart = TRUE;
	adef->constraint = TRUE;
	constraintSet = 1;
	break;
      case 'e':
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;
	break;
      case 'q':
	strcpy(modelFileName,optarg);
	adef->useMultipleModel = TRUE;
        break;
      case 'p':
	sscanf(optarg,"%ld", &(adef->parsimonySeed));	
	if(adef->parsimonySeed <= 0)
	  {
	    printf("Parsimony seed specified via -p must be greater than zero\n");
	    errorExit(-1);
	  }
	break;
      case 'N':
      case '#':
	if(sscanf(optarg,"%d", &multipleRuns) > 0)
	  {
	    adef->multipleRuns = multipleRuns;
	  }
	else
	  {
	    if((sscanf(optarg,"%s", aut) > 0) && ((strcmp(aut, "autoFC") == 0) || (strcmp(aut, "autoMR") == 0) || 
						  (strcmp(aut, "autoMRE") == 0) || (strcmp(aut, "autoMRE_IGN") == 0)))
						  
	      {
		adef->bootStopping = TRUE;
		adef->multipleRuns = 1000;

		if((strcmp(aut, "autoFC") == 0))	   
		  tr->bootStopCriterion = FREQUENCY_STOP;
		if((strcmp(aut, "autoMR") == 0))		  	    
		  tr->bootStopCriterion = MR_STOP;	   
		if((strcmp(aut, "autoMRE") == 0))	   
		  tr->bootStopCriterion = MRE_STOP;
		if((strcmp(aut, "autoMRE_IGN") == 0))
		  tr->bootStopCriterion = MRE_IGN_STOP;
	      }
	    else
	      {
		if(processID == 0)
		  {
		    printf("Use -# or -N option either with an integer, e.g., -# 100 or with -# autoFC or -# autoMR or -# autoMRE or -# autoMRE_IGN\n");
		    printf("or -N 100 or  -N autoFC or -N autoMR or -N autoMRE or -N autoMRE_IGN respectively, note that auto will not work for the\n");
		    printf("MPI-based parallel version\n");
		  }
		errorExit(0);
	      }
	  }
	multipleRunsSet = TRUE;
	break;
      case 'v':
	printVersionInfo();
	errorExit(0);
      case 'y':
	adef->startingTreeOnly = 1;
	break;     
      case 'h':
	printREADME();
	errorExit(0);
      case 'j':
	adef->checkpoints = 1;
	break;
      case 'a':
	strcpy(weightFileName,optarg);
	adef->useWeightFile = TRUE;
        break;
      case 'b':
	sscanf(optarg,"%ld", &adef->boot);
	if(adef->boot <= 0)
	  {
	    printf("Bootstrap seed specified via -b must be greater than zero\n");
	    errorExit(-1);
	  }
	bSeedSet = TRUE;
	break;
      case 'x':
	sscanf(optarg,"%ld", &adef->rapidBoot);
	if(adef->rapidBoot <= 0)
	  {
	    printf("Bootstrap seed specified via -x must be greater than zero\n");
	    errorExit(-1);
	  }
	xSeedSet = TRUE;
	break;
      case 'c':
	sscanf(optarg, "%d", &adef->categories);
	break;     
      case 'f':
	sscanf(optarg, "%c", &modelChar);
	switch(modelChar)
	  {
	  case 'A':
	    adef->mode = ANCESTRAL_STATES; 
	    /*adef->compressPatterns  = FALSE;*/
	    break;
	  case 'a':
	    adef->allInOne = TRUE;
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;
	  case 'b':
	    adef->readTaxaOnly = TRUE;
	    adef->mode = CALC_BIPARTITIONS;
	    break;
	  case 'B':
	    adef->mode = OPTIMIZE_BR_LEN_SCALER;
	    break;
	  case 'c':
	    adef->mode = CHECK_ALIGNMENT;
	    break;
	  case 'd':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;
	  case 'e':
	    adef->mode = TREE_EVALUATION;
	    break;
	  case 'F':
	    adef->mode = FAST_SEARCH;
	    adef->veryFast = FALSE;
	    break;
	  case 'E':
	    adef->mode = FAST_SEARCH;
	    adef->veryFast = TRUE;
	    break;
	  case 'g':
	    tr->useFastScaling = FALSE;
	    adef->mode = PER_SITE_LL;
	    break;
	  case 'h':
	    adef->mode = TREE_EVALUATION;
	    adef->likelihoodTest = TRUE;
	    tr->useFastScaling = FALSE;
	    break;
	  case 'i':
	    adef->mode = MESH_TREE_SEARCH;
	    adef->meshSearch = 0;
	    break;
	  case 'I':
	    adef->mode = MESH_TREE_SEARCH;
	    adef->meshSearch = 1;
	    break;
	  case 'j':
	    adef->mode = GENERATE_BS;
	    adef->generateBS = TRUE;
	    break;
	  case 'J':
	    adef->mode = SH_LIKE_SUPPORTS; 
	    tr->useFastScaling = FALSE;
	    break;
	  case 'm': 
	    adef->readTaxaOnly = TRUE;	    
	    adef->mode = COMPUTE_BIPARTITION_CORRELATION;
	    break;
	  case 'n':
	    adef->mode = COMPUTE_LHS;
	    break;
	  case 'o':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = FALSE;
	    break;
	  case 'q':
	    adef->mode = QUARTET_CALCULATION;
	    break;
	  case 'p':
	    adef->mode =  PARSIMONY_ADDITION;
	    break;	 
	  case 'r':
	    adef->readTaxaOnly = TRUE;
	    adef->mode = COMPUTE_RF_DISTANCE;
	    break;	  
	  case 's':
	    adef->mode = SPLIT_MULTI_GENE;
	    break;
	  case 'S':
	    adef->mode = EPA_SITE_SPECIFIC_BIAS;
	    tr->useFastScaling = FALSE;
	    adef->compressPatterns  = FALSE;
	    break;
	  case 't':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    adef->permuteTreeoptimize = TRUE;
	    break;
	  case 'T':
	    adef->mode = THOROUGH_OPTIMIZATION;
	    break;
	  case 'u':
	    adef->mode = MORPH_CALIBRATOR;
	    tr->useFastScaling = FALSE;
	    adef->compressPatterns  = FALSE;	    
	    break;	  
	  case 'v':	    
	    adef->mode = CLASSIFY_ML;	   
	   	   
#ifdef _PAVLOS
	    adef->compressPatterns  = FALSE; 
#endif
#ifdef _USE_PTHREADS
	    tr->useFastScaling = FALSE;
#endif
	    break;
	  case 'y':
	    adef->mode = CLASSIFY_MP;
	    break;
	  case 'w':	    
	    adef->mode = COMPUTE_ELW;
	    adef->computeELW = TRUE;
	    break;
	  case 'x':
	    adef->mode = DISTANCE_MODE;
	    adef->computeDistance = TRUE;
	    break;	  	  	  	  	  	  	     
	  default:
	    {
	      if(processID == 0)
		{
		  printf("Error select one of the following algorithms via -f :\n");
		  printMinusFUsage();
		}
	      errorExit(-1);
	    }
	  }
	break;
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;
      case 'n':
        strcpy(run_id,optarg);
	analyzeRunId(run_id);
	nameSet = 1;
        break;
      case 'w':
        strcpy(resultDir, optarg);
	resultDirSet = TRUE;
        break;
      case 't':
	strcpy(tree_file, optarg);
	adef->restart = TRUE;
	treeSet = 1;
	break;
      case 's':
	strcpy(seq_file, optarg);
	alignmentSet = 1;
	break;
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, adef) == 0)
	  {
	    if(processID == 0)
	      {
		printf("Model %s does not exist\n\n", model);
                printf("For BINARY data use: BINCAT                or BINGAMMA                or\n");
		printf("                     BINCATI               or BINGAMMAI                 \n");
		printf("For DNA data use:    GTRCAT                or GTRGAMMA                or\n");
		printf("                     GTRCATI               or GTRGAMMAI                 \n");
		printf("For AA data use:     PROTCATmatrixName[F]  or PROTGAMMAmatrixName[F]  or\n");
		printf("                     PROTCATImatrixName[F] or PROTGAMMAImatrixName[F]   \n");
		printf("The AA substitution matrix can be one of the following: \n");
		printf("DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, LG, MTART, MTZOA, PMB, HIVB, HIVW, JTTDCMUT, FLU, GTR\n\n");
		printf("With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
		printf("Please note that for mixed models you can in addition specify the per-gene model in\n");
		printf("the mixed model file (see manual for details)\n");
	      }
	    errorExit(-1);
	  }
	else
	  modelSet = 1;
	break;
      default:
	errorExit(-1);
    }
  }

 

#ifdef _USE_PTHREADS
  if(NumberOfThreads < 2)
    {
      printf("\nThe number of threads is currently set to %d\n", NumberOfThreads);
      printf("Specify the number of threads to run via -T numberOfThreads\n");
      printf("NumberOfThreads must be set to an integer value greater than 1\n\n");
      errorExit(-1);
    }
#endif

 

  if(bSeedSet && xSeedSet)
    {
      printf("Error, you can't seed random seeds by using -x and -b at the same time\n");
      printf("use either -x or -b, exiting ......\n");
      errorExit(-1);
    }

  if(bSeedSet || xSeedSet)
    {
      if(!multipleRunsSet)
	{
	  printf("Error, you have specified a random number seed via -x or -b for some sort of bootstrapping,\n");
	  printf("but you have not specified a number of replicates via -N or -#, exiting ....\n");
	  errorExit(-1);
	}
      
      if(adef->multipleRuns == 1)
	{
	  printf("WARNING, you have specified a random number seed via -x or -b for some sort of bootstrapping,\n");
	  printf("but you have specified a number of replicates via -N or -# euqal to one\n");
	  printf("Are you really sure that this is what you want to do?\n");
	}

     
    }


  if(adef->computeELW)
    {
      if(processID == 0)
	{
	  if(adef->boot == 0)
	    {
	      printf("Error, you must specify a bootstrap seed via \"-b\" to compute ELW statistics\n");
	      errorExit(-1);
	    }

	  if(adef->multipleRuns < 2)
	    {
	      printf("Error, you must specify the number of BS replicates via \"-#\" or \"-N\" to compute ELW statistics\n");
	      printf("it should be larger than 1, recommended setting is 100\n");
	      errorExit(-1);
	    }

	  if(!treesSet)
	    {
	      printf("Error, you must specify an input file containing several candidate trees\n");
	      printf("via \"-z\" to compute ELW statistics.\n");
	      errorExit(-1);
	    }

	  if(!isGamma(adef))
	    {
	      printf("Error ELW test can only be conducted undetr GAMMA or GAMMA+P-Invar models\n");
	      errorExit(-1);
	    }
	}
    }


  if(isGamma(adef) && tr->noRateHet)
    {
      printf("\n\nError: using a  model without any rate heterogeneity (enabled via \"-V\") only works if you specify a CAT model\n");
      printf("via the \"-m\" switch, exiting ....\n\n");
      errorExit(-1);
    }
  
  if(((!adef->boot) && (!adef->rapidBoot)) && adef->bootStopping)
    {
      if(processID == 0)
	{
	  printf("Can't use automatic bootstopping without actually doing a Bootstrap\n");
	  printf("Specify either -x randomNumberSeed (rapid) or -b randomNumberSeed (standard)\n");
	  errorExit(-1);
	}
    }

  if(adef->boot && adef->rapidBoot)
    {
      if(processID == 0)
	{
	  printf("Can't use standard and rapid BOOTSTRAP simultaneously\n");
	  errorExit(-1);
	}
    }

  if(adef->rapidBoot)
    {
      if(processID == 0 && (adef->restart || treesSet) && !(groupSet || constraintSet))
	{
	  printf("Error, starting tree(s) will be ignored by rapid Bootstrapping\n");
	  errorExit(-1);
	}
    }

  if(adef->allInOne && (adef->rapidBoot == 0))
    {
      if(processID == 0)
	{
	  printf("Error, to carry out an ML search after a rapid BS inference you must specify a random number seed with -x\n");
	  errorExit(-1);
	}
    }


  

  if(adef->mode == PER_SITE_LL)
    {
      if(!isGamma(adef))
	{
	  if(processID == 0)
	    printf("\n ERROR: Computation of per-site log LHs is only allowed under GAMMA model of rate heterogeneity!\n");
	  errorExit(-1);
	}

      if(!treesSet)
	{
	  if(processID == 0)
	    printf("\n ERROR: For Computation of per-site log LHs you need to specify several input trees with \"-z\"\n");
	  errorExit(-1);
	}
    }



  if(adef->mode == SPLIT_MULTI_GENE && (!adef->useMultipleModel))
    {
      if(processID == 0)
	{
	  printf("\n  Error, you are trying to split a multi-gene alignment into individual genes with the \"-f s\" option\n");
	  printf("Without specifying a multiple model file with \"-q modelFileName\" \n");
	}
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS && !treesSet)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition computation mode you must specify a file containing multiple trees with the \"-z\" option\n");
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS && !adef->restart)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition computation mode you must specify a tree on which bipartition information will be drawn with the \"-t\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    {
      if(processID == 0)
	printf("\n Error, you must specify a model of substitution with the \"-m\" option\n");
      errorExit(-1);
    }

  if(adef->computeDistance)
    {
      if(isCat(adef))
	{
	  if(processID == 0)
	    printf("\n Error pairwise distance computation only allowed for GAMMA-based models of rate heterogeneity\n");
	  errorExit(-1);
	}

      if(adef->restart)
	{
	  if(adef->randomStartingTree)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for random starting trees\n");
	      errorExit(-1);
	    }

	  if(adef->constraint)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for binary backbone  constraint tree\n");
	      errorExit(-1);
	    }

	  if(adef->grouping)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for constraint tree\n");
	      errorExit(-1);
	    }

	}

      if(adef->boot || adef->rapidBoot)
	{
	  if(processID == 0)
	    printf("\n Bootstrapping not implemented for pairwise distance computation\n");
	  errorExit(-1);
	}
    }








  if(!adef->restart && adef->mode == PARSIMONY_ADDITION)
    {
       if(processID == 0)
	 {
	   printf("\n You need to specify an incomplete binary input tree with \"-t\" to execute \n");
	   printf(" RAxML MP stepwise addition with \"-f p\"\n");
	 }
      errorExit(-1);
    }



  if(adef->restart && adef->randomStartingTree)
    {
      if(processID == 0)
	{
	  if(adef->constraint)
	    {
	      printf("\n Error you specified a binary constraint tree with -r AND the computation\n");
	      printf("of a random starting tree with -d for the same run\n");
	    }
	  else
	    {
	      if(adef->grouping)
		{
		  printf("\n Error you specified a multifurcating constraint tree with -g AND the computation\n");
		  printf("of a random starting tree with -d for the same run\n");
		}
	      else
		{
		  printf("\n Error you specified a starting tree with -t AND the computation\n");
		  printf("of a random starting tree with -d for the same run\n");
		}
	    }
	}
      errorExit(-1);
    }

  if(adef->outgroup && adef->mode == ANCESTRAL_STATES)
    {
      if(processID == 0)
	{
	  printf("\n Specifying an outgroup for ancestral state reconstruction is not allowed\n");
	  printf(" You already need to specify a rooted input tree for computing ancestral states anyway.\n\n");
	}
      errorExit(-1);
    }

  if(!treeSet && adef->mode == ANCESTRAL_STATES)
    {
      if(processID == 0)
	printf("\n Error you need to specify a ROOTED binary reference tree for ancestral state computations\n");
      errorExit(-1);
    }
  
  if(treeSet && constraintSet)
    {
      if(processID == 0)
	printf("\n Error you specified a binary constraint tree AND a starting tree for the same run\n");
      errorExit(-1);
    }


  if(treeSet && groupSet)
    {
      if(processID == 0)
	printf("\n Error you specified a multifurcating constraint tree AND a starting tree for the same run\n");
      errorExit(-1);
    }


  if(groupSet && constraintSet)
    {
      if(processID == 0)
	printf("\n Error you specified a bifurcating constraint tree AND a multifurcating constraint tree for the same run\n");
      errorExit(-1);
    }

  if(adef->restart && adef->startingTreeOnly)
    {
      if(processID == 0)
	{
	  printf("\n Error conflicting options: you want to compute only a parsimony starting tree with -y\n");
	  printf(" while you actually specified a starting tree with -t %s\n", tree_file);
	}
      errorExit(-1);
    }

  if((adef->mode == TREE_EVALUATION || adef->mode == OPTIMIZE_BR_LEN_SCALER) && (!adef->restart))
    {
      if(processID == 0)
	printf("\n Error: please specify a treefile for the tree you want to evaluate with -t\n");
      errorExit(-1);
    }

#ifdef _WAYNE_MPI

  if(adef->mode == SPLIT_MULTI_GENE)
    {
      if(processID == 0)
	printf("Multi gene alignment splitting (-f s) not implemented for the MPI-Version\n");
      errorExit(-1);
    }

  if(adef->mode == TREE_EVALUATION)
    {
      if(processID == 0)
	printf("Tree Evaluation mode (-f e) not implemented for the MPI-Version\n");
      errorExit(-1);
    }
  
  if(adef->mode == OPTIMIZE_BR_LEN_SCALER)
    {
      if(processID == 0)
	printf("Branch length scaler optimization mode (-f B) not implemented for the MPI-Version\n");
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS)
    {
      if(processID == 0)
	 printf("Computation of bipartitions (-f b) not implemented for the MPI-Version\n");
      errorExit(-1);
    }

  if(adef->multipleRuns == 1)
    {
      if(processID == 0)
	{
	  printf("Error: you are running the parallel MPI program but only want to compute one tree\n");
	  printf("For the MPI version you must specify a number of trees greater than 1 with the -# or -N option\n");
	}
      errorExit(-1);
    }

#endif

   if((adef->mode == TREE_EVALUATION || adef->mode == OPTIMIZE_BR_LEN_SCALER) && (isCat(adef)))
     {
       if(processID == 0)
	 {
	   printf("\n Warning: tree evaluation with CAT model of rate heterogeneity\n");
	   printf("Only compare likelihood values for identical rate category assignments\n");	  
	   printf("CAT-based Branch lengths are on average shorter by factor 0.5 than GAMMA-based branch lengths\n");
	   printf("... but highly correlated with GAMMA branch lengths\n");
	 }      
     }

  if(!nameSet)
    {
      if(processID == 0)
	printf("\n Error: please specify a name for this run with -n\n");
      errorExit(-1);
    }

  if(! alignmentSet && !adef->readTaxaOnly)
    {
      if(processID == 0)
	printf("\n Error: please specify an alignment for this run with -s\n");
      errorExit(-1);
    }

  
  {
#ifdef WIN32
    const 
      char *separator = "\\";
#else
    const 
      char *separator = "/";
#endif

    if(resultDirSet)
      {
	char 
	  dir[1024] = "";
	
#ifndef WIN32
	if(resultDir[0] != separator[0])
	  strcat(dir, separator);
#endif
	
	strcat(dir, resultDir);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	strcpy(workdir, dir);
      }
    else
      {
	char 
	  dir[1024] = "",
	  *result = getcwd(dir, sizeof(dir));
	
	assert(result != (char*)NULL);
	
	if(dir[strlen(dir) - 1] != separator[0]) 
	  strcat(dir, separator);
	
	strcpy(workdir, dir);		
      }
  }

  return;
}




void errorExit(int e)
{

#ifdef _WAYNE_MPI
  MPI_Finalize();
#endif

  exit(e);

}



static void makeFileNames(void)
{
  int infoFileExists = 0;


  strcpy(permFileName,         workdir);
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);
  strcpy(checkpointFileName,   workdir);
  strcpy(infoFileName,         workdir);
  strcpy(randomFileName,       workdir);
  strcpy(bootstrapFileName,    workdir);
  strcpy(bipartitionsFileName, workdir);
  strcpy(bipartitionsFileNameBranchLabels, workdir);
  strcpy(ratesFileName,        workdir);
  strcpy(lengthFileName,       workdir);
  strcpy(lengthFileNameModel,  workdir);
  strcpy(perSiteLLsFileName,  workdir);
  strcpy(binaryModelParamsOutputFileName,  workdir);

  strcat(permFileName,         "RAxML_parsimonyTree.");
  strcat(resultFileName,       "RAxML_result.");
  strcat(logFileName,          "RAxML_log.");
  strcat(checkpointFileName,   "RAxML_checkpoint.");
  strcat(infoFileName,         "RAxML_info.");
  strcat(randomFileName,       "RAxML_randomTree.");
  strcat(bootstrapFileName,    "RAxML_bootstrap.");
  strcat(bipartitionsFileName, "RAxML_bipartitions.");
  strcat(bipartitionsFileNameBranchLabels, "RAxML_bipartitionsBranchLabels.");
  strcat(ratesFileName,        "RAxML_perSiteRates.");
  strcat(lengthFileName,       "RAxML_treeLength.");
  strcat(lengthFileNameModel,  "RAxML_treeLengthModel.");
  strcat(perSiteLLsFileName,   "RAxML_perSiteLLs.");
  strcat(binaryModelParamsOutputFileName,   "RAxML_binaryModelParameters.");

  strcat(permFileName,         run_id);
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);
  strcat(checkpointFileName,   run_id);
  strcat(infoFileName,         run_id);
  strcat(randomFileName,       run_id);
  strcat(bootstrapFileName,    run_id);
  strcat(bipartitionsFileName, run_id);
  strcat(bipartitionsFileNameBranchLabels, run_id);  
  strcat(ratesFileName,        run_id);
  strcat(lengthFileName,       run_id);
  strcat(lengthFileNameModel,  run_id);
  strcat(perSiteLLsFileName,   run_id);
  strcat(binaryModelParamsOutputFileName, run_id);

#ifdef _WAYNE_MPI  
  {
    char buf[64];
    
    strcpy(bootstrapFileNamePID, bootstrapFileName);
    strcat(bootstrapFileNamePID, ".PID.");
    sprintf(buf, "%d", processID);
    strcat(bootstrapFileNamePID, buf);
  }
#endif

  if(processID == 0)
    {
      infoFileExists = filexists(infoFileName);

      if(infoFileExists)
	{
	  printf("RAxML output files with the run ID <%s> already exist \n", run_id);
	  printf("in directory %s ...... exiting\n", workdir);

	  exit(-1);
	}
    }
}




 




/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


void printBaseFrequencies(tree *tr)
{
  if(processID == 0)
    {
      int 
	model;

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int i;

	  printBothOpen("Partition: %d with name: %s\n", model, tr->partitionData[model].partitionName);
	  printBothOpen("Base frequencies: ");
	  
	  for(i = 0; i < tr->partitionData[model].states; i++)
	    printBothOpen("%1.3f ", tr->partitionData[model].frequencies[i]);
	  
	  printBothOpen("\n\n");
	}	      
    }
}

static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{
  if(processID == 0)
    {
      int i, model;
      FILE *infoFile = myfopen(infoFileName, "ab");
      char modelType[128];

      if(!adef->readTaxaOnly)
	{
	  if(adef->useInvariant)
	    strcpy(modelType, "GAMMA+P-Invar");
	  else
	    strcpy(modelType, "GAMMA");
	}
     
      printBoth(infoFile, "\n\nThis is %s version %s released by Alexandros Stamatakis in %s.\n\n",  programName, programVersion, programDate);
      printBoth(infoFile, "With greatly appreciated code contributions by:\n");
      printBoth(infoFile, "Andre Aberer (HITS)\n");     
      printBoth(infoFile, "Simon Berger (HITS)\n");     
      printBoth(infoFile, "Nick Pattengale (Sandia)\n"); 
      printBoth(infoFile, "Wayne Pfeiffer (SDSC)\n");
      printBoth(infoFile, "Akifumi S. Tanabe (Univ. Tsukuba)\n\n");
      
      if(!adef->readTaxaOnly)
	{
	  if(!adef->compressPatterns)
	    printBoth(infoFile, "\nAlignment has %d columns\n\n",  tr->cdta->endsite);
	  else
	    printBoth(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->cdta->endsite);
	  
	  if(adef->useInvariant)
	    printBoth(infoFile, "Found %d invariant alignment patterns that correspond to %d columns \n", tr->numberOfInvariableColumns, tr->weightOfInvariableColumns);

	  printBoth(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %3.2f%s\n", 100.0 * adef->gapyness, "%");
	}

      switch(adef->mode)
	{	
	case DISTANCE_MODE:
	  printBoth(infoFile, "\nRAxML Computation of pairwise distances\n\n");
	  break;
	case TREE_EVALUATION :
	  printBoth(infoFile, "\nRAxML Model Optimization up to an accuracy of %f log likelihood units\n\n", adef->likelihoodEpsilon);
	  break;
	case  BIG_RAPID_MODE:
	  if(adef->rapidBoot)
	    {
	      if(adef->allInOne)
		printBoth(infoFile, "\nRAxML rapid bootstrapping and subsequent ML search\n\n");
	      else
		printBoth(infoFile,  "\nRAxML rapid bootstrapping algorithm\n\n");
	    }
	  else
	    printBoth(infoFile, "\nRAxML rapid hill-climbing mode\n\n");
	  break;
	case CALC_BIPARTITIONS:
	  printBoth(infoFile, "\nRAxML Bipartition Computation: Drawing support values from trees in file %s onto tree in file %s\n\n",
		    bootStrapFile, tree_file);
	  break;
	case PER_SITE_LL:
	  printBoth(infoFile, "\nRAxML computation of per-site log likelihoods\n");
	  break;
	case PARSIMONY_ADDITION:
	  printBoth(infoFile, "\nRAxML stepwise MP addition to incomplete starting tree\n\n");
	  break;
	case CLASSIFY_ML:
	  printBoth(infoFile, "\nRAxML likelihood-based placement algorithm\n\n");
	  break;
	case CLASSIFY_MP:
	  printBoth(infoFile, "\nRAxML parsimony-based placement algorithm\n\n");
	  break;
	case GENERATE_BS:
	  printBoth(infoFile, "\nRAxML BS replicate generation\n\n");
	  break;
	case COMPUTE_ELW:
	  printBoth(infoFile, "\nRAxML ELW test\n\n");
	  break;
	case BOOTSTOP_ONLY:
	  printBoth(infoFile, "\nRAxML a posteriori Bootstrap convergence assessment\n\n");
	  break;
	case CONSENSUS_ONLY:
	  if(adef->leaveDropMode)
	    printBoth(infoFile, "\nRAxML rogue taxa computation by Andre Aberer (HITS)\n\n");
	  else
	    printBoth(infoFile, "\nRAxML consensus tree computation\n\n");
	  break;
	case COMPUTE_LHS:
	  printBoth(infoFile, "\nRAxML computation of likelihoods for a set of trees\n\n");
	  break;
	case COMPUTE_BIPARTITION_CORRELATION:
	  printBoth(infoFile, "\nRAxML computation of bipartition support correlation on two sets of trees\n\n");
	  break;
	case COMPUTE_RF_DISTANCE:
	  printBoth(infoFile, "\nRAxML computation of RF distances for all pairs of trees in a set of trees\n\n");
	  break;
	case MORPH_CALIBRATOR:
	  printBoth(infoFile, "\nRAxML morphological calibrator using Maximum Likelihood\n\n");
	  break;		  
	case MESH_TREE_SEARCH:
	  printBoth(infoFile, "\nRAxML experimental mesh tree search\n\n");
	  break;
	case FAST_SEARCH:
	  printBoth(infoFile, "\nRAxML experimental very fast tree search\n\n");
	  break;
	case SH_LIKE_SUPPORTS:
	  printBoth(infoFile, "\nRAxML computation of SH-like support values on a given tree\n\n");
	  break;
	case EPA_SITE_SPECIFIC_BIAS:
	  printBoth(infoFile, "\nRAxML exprimental site-specfific phylogenetic placement bias analysis algorithm\n\n");
	  break;

	case ANCESTRAL_STATES:
	  printBoth(infoFile, "\nRAxML marginal ancestral state computation\n\n");
	  break;
	case  QUARTET_CALCULATION:
	  printBoth(infoFile, "\nRAxML quartet computation\n\n");
	  break;
	case THOROUGH_OPTIMIZATION:
	  printBoth(infoFile, "\nRAxML thorough tree optimization\n\n");
	  break;
	case OPTIMIZE_BR_LEN_SCALER :
	  printBoth(infoFile, "\nRAxML Branch length scaler and other model parameter optimization up to an accuracy of %f log likelihood units\n\n", adef->likelihoodEpsilon);
	  break;
	default:
	  assert(0);
	}

    
      if(!adef->readTaxaOnly)
	{
	  if(adef->perGeneBranchLengths)
	    printBoth(infoFile, "Using %d distinct models/data partitions with individual per partition branch length optimization\n\n\n", tr->NumberOfModels);
	  else
	    printBoth(infoFile, "Using %d distinct models/data partitions with joint branch length optimization\n\n\n", tr->NumberOfModels);
	}    

      if(adef->mode == BIG_RAPID_MODE)
	{
	  if(adef->rapidBoot)
	    {
	      if(adef->allInOne)
		printBoth(infoFile, "\nExecuting %d rapid bootstrap inferences and thereafter a thorough ML search \n\n", adef->multipleRuns);
	      else
		printBoth(infoFile, "\nExecuting %d rapid bootstrap inferences\n\n", adef->multipleRuns);
	    }
	  else
	    {
	      if(adef->boot)
		printBoth(infoFile, "Executing %d non-parametric bootstrap inferences\n\n", adef->multipleRuns);
	      else
		{
		  char treeType[1024];

		  if(adef->restart)
		    strcpy(treeType, "user-specifed");
		  else
		    {
		      if(adef->randomStartingTree)
			strcpy(treeType, "distinct complete random");
		      else
			strcpy(treeType, "distinct randomized MP");
		    }

		  printBoth(infoFile, "Executing %d inferences on the original alignment using %d %s trees\n\n",
			    adef->multipleRuns, adef->multipleRuns, treeType);
		}
	    }
	}


      if(!adef->readTaxaOnly)
	{	  
	  printBoth(infoFile, "All free model parameters will be estimated by RAxML\n");

	  
	  if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
	    printBoth(infoFile, "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n\n", modelType);
	  else
	    {
	      printBoth(infoFile, "ML estimate of %d per site rate categories\n\n", adef->categories);
	      if(adef->mode != CLASSIFY_ML && adef->mode != CLASSIFY_MP)
		printBoth(infoFile, "Likelihood of final tree will be evaluated and optimized under %s\n\n", modelType);
	    }
	  
	  if(adef->mode != CLASSIFY_ML && adef->mode != CLASSIFY_MP)
	    printBoth(infoFile, "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n",
		      modelType, adef->likelihoodEpsilon);
	  

	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      printBoth(infoFile, "Partition: %d\n", model);
	      printBoth(infoFile, "Alignment Patterns: %d\n", tr->partitionData[model].upper - tr->partitionData[model].lower);
	      printBoth(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);
	      
	      switch(tr->partitionData[model].dataType)
		{
		case DNA_DATA:
		  printBoth(infoFile, "DataType: DNA\n");		  
		  printBoth(infoFile, "Substitution Matrix: GTR\n");
		  break;
		case AA_DATA:
		  assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
		  printBoth(infoFile, "DataType: AA\n");
		  if(tr->partitionData[model].protModels != PROT_FILE)
		    {		     		     
		      printBoth(infoFile, "Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);		      
		      printBoth(infoFile, "Using %s base frequencies\n", (tr->partitionData[model].usePredefinedProtFreqs == TRUE)?"fixed":"empirical");		      
		    }
		  else
		    {
		       printBoth(infoFile, "Substitution Matrix File name: %s\n", tr->partitionData[model].proteinSubstitutionFileName);
		       printBoth(infoFile, "Using base frequencies as provided in the model file\n");
		    }
		  break;
		case BINARY_DATA:
		  printBoth(infoFile, "DataType: BINARY/MORPHOLOGICAL\n");		  
		  printBoth(infoFile, "Substitution Matrix: Uncorrected\n");
		  break;
		case SECONDARY_DATA:
		  printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");		  
		  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		  break;
		case SECONDARY_DATA_6:
		  printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");		  
		  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		  break;
		case SECONDARY_DATA_7:
		  printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");		 
		  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		  break;
		case GENERIC_32:
		  printBoth(infoFile, "DataType: Multi-State with %d distinct states in use (maximum 32)\n",tr->partitionData[model].states);		  
		  switch(tr->multiStateModel)
		    {
		    case ORDERED_MULTI_STATE:
		      printBoth(infoFile, "Substitution Matrix: Ordered Likelihood\n");
		      break;
		    case MK_MULTI_STATE:
		      printBoth(infoFile, "Substitution Matrix: MK model\n");
		      break;
		    case GTR_MULTI_STATE:
		      printBoth(infoFile, "Substitution Matrix: GTR\n");
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case GENERIC_64:
		  printBoth(infoFile, "DataType: Codon\n");		  
		  break;		
		default:
		  assert(0);
		}
	      printBoth(infoFile, "\n\n\n");
	    }
	}

      printBoth(infoFile, "\n");

      printBoth(infoFile, "RAxML was called as follows:\n\n");
      for(i = 0; i < argc; i++)
	printBoth(infoFile,"%s ", argv[i]);
      printBoth(infoFile,"\n\n\n");

      fclose(infoFile);
    }
}

void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "", treeID[64] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
    {    
    case MESH_TREE_SEARCH:    
    case MORPH_CALIBRATOR:
      break;
    case TREE_EVALUATION:


      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE);

      logFile = myfopen(temporaryFileName, "wb");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

      if(adef->perGeneBranchLengths)
	printTreePerGene(tr, adef, temporaryFileName, "wb");


      break;
    case BIG_RAPID_MODE:
      if(!adef->boot)
	{
	  if(adef->multipleRuns > 1)
	    {
	      sprintf(treeID, "%d", tr->treeID);
	      strcat(temporaryFileName, ".RUN.");
	      strcat(temporaryFileName, treeID);
	    }


	  if(finalPrint)
	    {
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		case GAMMA_I:
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef,
			      SUMMARIZE_LH, FALSE, FALSE);

		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);

		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "wb");
		  break;
		case CAT:
		  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
			      NO_BRANCHES, FALSE, FALSE);

		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);

		  break;
		default:
		  assert(0);
		}
	    }
	  else
	    {
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef,
			  NO_BRANCHES, FALSE, FALSE);
	      logFile = myfopen(temporaryFileName, "wb");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);
	    }
	}
      break;
    default:
      printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
      exit(-1);
      break;
    }
}

void printBootstrapResult(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE 
    *logFile;
#ifdef _WAYNE_MPI
  char 
    *fileName = bootstrapFileNamePID;
#else
  char 
    *fileName = bootstrapFileName;
#endif

  if(adef->mode == BIG_RAPID_MODE && (adef->boot || adef->rapidBoot))
    {
      if(adef->bootstrapBranchLengths)
	{
	  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE);

	  logFile = myfopen(fileName, "ab");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	  
	  if(adef->perGeneBranchLengths)
	    printTreePerGene(tr, adef, fileName, "ab");
	}
      else
	{
	  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE);
	  
	  logFile = myfopen(fileName, "ab");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	}
    }
  else
    {
      printf("FATAL ERROR in  printBootstrapResult\n");
      exit(-1);
    }
}



void printBipartitionResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0 || adef->allInOne)
    {
      FILE *logFile;

      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, TRUE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE);
      logFile = myfopen(bipartitionsFileName, "ab");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, TRUE, FALSE);

      logFile = myfopen(bipartitionsFileNameBranchLabels, "ab");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);

    }
}



void printLog(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "", checkPoints[1024] = "", treeID[64] = "";
  double lh, t;

  lh = tr->likelihood;
  t = gettime() - masterTime;

  strcpy(temporaryFileName, logFileName);
  strcpy(checkPoints,       checkpointFileName);

  switch(adef->mode)
    {
    case TREE_EVALUATION:
      logFile = myfopen(temporaryFileName, "ab");

      printf("%f %f\n", t, lh);
      fprintf(logFile, "%f %f\n", t, lh);

      fclose(logFile);
      break;
    case BIG_RAPID_MODE:
      if(adef->boot || adef->rapidBoot)
	{
	  /* testing only printf("%f %f\n", t, lh);*/
	  /* NOTHING PRINTED so far */
	}
      else
	{
	  if(adef->multipleRuns > 1)
	    {
	      sprintf(treeID, "%d", tr->treeID);
	      strcat(temporaryFileName, ".RUN.");
	      strcat(temporaryFileName, treeID);

	      strcat(checkPoints, ".RUN.");
	      strcat(checkPoints, treeID);
	    }


	  if(!adef->checkpoints)
	    {
	      logFile = myfopen(temporaryFileName, "ab");

	      fprintf(logFile, "%f %f\n", t, lh);

	      fclose(logFile);
	    }
	  else
	    {
	      logFile = myfopen(temporaryFileName, "ab");

	      fprintf(logFile, "%f %f %d\n", t, lh, tr->checkPointCounter);

	      fclose(logFile);

	      strcat(checkPoints, ".");

	      sprintf(treeID, "%d", tr->checkPointCounter);
	      strcat(checkPoints, treeID);

	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE);

	      logFile = myfopen(checkPoints, "ab");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);

	      tr->checkPointCounter++;
	    }
	}
      break;   
    case MORPH_CALIBRATOR:
      break;
    default:
      assert(0);
    }
}



void printStartingTree(tree *tr, analdef *adef, boolean finalPrint)
{
  if(adef->boot)
    {
      /* not printing starting trees for bootstrap */
    }
  else
    {
      FILE *treeFile;
      char temporaryFileName[1024] = "", treeID[64] = "";

      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE);

      if(adef->randomStartingTree)
	strcpy(temporaryFileName, randomFileName);
      else
	strcpy(temporaryFileName, permFileName);

      if(adef->multipleRuns > 1)
	{
	  sprintf(treeID, "%d", tr->treeID);
	  strcat(temporaryFileName, ".RUN.");
	  strcat(temporaryFileName, treeID);
	}

      treeFile = myfopen(temporaryFileName, "ab");
      fprintf(treeFile, "%s", tr->tree_string);
      fclose(treeFile);
    }
}

void writeInfoFile(analdef *adef, tree *tr, double t)
{

    {      
      switch(adef->mode)
	{
	case MESH_TREE_SEARCH:
	  break;
	case TREE_EVALUATION:
	  break;
	case BIG_RAPID_MODE:
	  if(adef->boot || adef->rapidBoot)
	    {
	      if(!adef->initialSet)	
		printBothOpen("Bootstrap[%d]: Time %f seconds, bootstrap likelihood %f, best rearrangement setting %d\n", tr->treeID, t, tr->likelihood,  adef->bestTrav);		
	      else	
		printBothOpen("Bootstrap[%d]: Time %f seconds, bootstrap likelihood %f\n", tr->treeID, t, tr->likelihood);		
	    }
	  else
	    {
	      int model;
	      char modelType[128];

	      switch(tr->rateHetModel)
		{
		case GAMMA_I:
		  strcpy(modelType, "GAMMA+P-Invar");
		  break;
		case GAMMA:
		  strcpy(modelType, "GAMMA");
		  break;
		case CAT:
		  strcpy(modelType, "CAT");
		  break;
		default:
		  assert(0);
		}

	      if(!adef->initialSet)		
		printBothOpen("Inference[%d]: Time %f %s-based likelihood %f, best rearrangement setting %d\n",
			      tr->treeID, t, modelType, tr->likelihood,  adef->bestTrav);		 
	      else		
		printBothOpen("Inference[%d]: Time %f %s-based likelihood %f\n",
			      tr->treeID, t, modelType, tr->likelihood);		 

	      {
		FILE *infoFile = myfopen(infoFileName, "ab");

		for(model = 0; model < tr->NumberOfModels; model++)
		  {
		    fprintf(infoFile, "alpha[%d]: %f ", model, tr->partitionData[model].alpha);
		    if(adef->useInvariant)
		      fprintf(infoFile, "invar[%d]: %f ", model, tr->partitionData[model].propInvariant);

		    if(tr->partitionData[model].dataType == DNA_DATA)
		      {
			int 
			  k,
			  states = tr->partitionData[model].states,
			  rates = ((states * states - states) / 2);
			
			fprintf(infoFile, "rates[%d] ac ag at cg ct gt: ", model);
			for(k = 0; k < rates; k++)
			  fprintf(infoFile, "%f ", tr->partitionData[model].substRates[k]);
		      }		 

		  }

		fprintf(infoFile, "\n");
		fclose(infoFile);
	      }
	    }
	  break;
	default:
	  assert(0);
	}      
    }
}

static void printFreqs(int n, double *f, char **names)
{
  int k;

  for(k = 0; k < n; k++)
    printBothOpen("freq pi(%s): %f\n", names[k], f[k]);
}

static void printRatesDNA_BIN(int n, double *r, char **names)
{
  int i, j, c;

  for(i = 0, c = 0; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	{
	  if(i == n - 2 && j == n - 1)
	    printBothOpen("rate %s <-> %s: %f\n", names[i], names[j], 1.0);
	  else
	    printBothOpen("rate %s <-> %s: %f\n", names[i], names[j], r[c]);
	  c++;
	}
    }
}

static void printRatesRest(int n, double *r, char **names)
{
  int i, j, c;

  for(i = 0, c = 0; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	{
	  printBothOpen("rate %s <-> %s: %f\n", names[i], names[j], r[c]);
	  c++;
	}
    }
}


void getDataTypeString(tree *tr, int model, char typeOfData[1024])
{
  switch(tr->partitionData[model].dataType)
    {
    case AA_DATA:
      strcpy(typeOfData,"AA");
      break;
    case DNA_DATA:
      strcpy(typeOfData,"DNA");
      break;
    case BINARY_DATA:
      strcpy(typeOfData,"BINARY/MORPHOLOGICAL");
      break;
    case SECONDARY_DATA:
      strcpy(typeOfData,"SECONDARY 16 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_6:
      strcpy(typeOfData,"SECONDARY 6 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case SECONDARY_DATA_7:
      strcpy(typeOfData,"SECONDARY 7 STATE MODEL USING ");
      strcat(typeOfData, secondaryModelList[tr->secondaryStructureModel]);
      break;
    case GENERIC_32:
      strcpy(typeOfData,"Multi-State");
      break;
    case GENERIC_64:
      strcpy(typeOfData,"Codon"); 
      break;
    default:
      assert(0);
    }
}



void printModelParams(tree *tr, analdef *adef)
{
  int
    model;

  double
    *f = (double*)NULL,
    *r = (double*)NULL;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      double tl;
      char typeOfData[1024];

      getDataTypeString(tr, model, typeOfData);      

      printBothOpen("Model Parameters of Partition %d, Name: %s, Type of Data: %s\n",
		    model, tr->partitionData[model].partitionName, typeOfData);
      printBothOpen("alpha: %f\n", tr->partitionData[model].alpha);

      if(adef->useInvariant)
	printBothOpen("invar: %f\n", tr->partitionData[model].propInvariant);

      if(adef->perGeneBranchLengths)
	tl = treeLength(tr, model);
      else
	tl = treeLength(tr, 0);

      printBothOpen("Tree-Length: %f\n", tl);

      f = tr->partitionData[model].frequencies;
      r = tr->partitionData[model].substRates;

      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  {
	    char *freqNames[20] = {"A", "R", "N ","D", "C", "Q", "E", "G",
				   "H", "I", "L", "K", "M", "F", "P", "S",
				   "T", "W", "Y", "V"};

	    printRatesRest(20, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(20, f, freqNames);
	  }
	  break;
	case GENERIC_32:
	  {
	    char *freqNames[32] = {"0", "1", "2", "3", "4", "5", "6", "7", 
				   "8", "9", "A", "B", "C", "D", "E", "F",
				   "G", "H", "I", "J", "K", "L", "M", "N",
				   "O", "P", "Q", "R", "S", "T", "U", "V"}; 

	    printRatesRest(32, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(32, f, freqNames);
	  }
	  break;
	case GENERIC_64:
	  assert(0);
	  break;
	case DNA_DATA:
	  {
	    char *freqNames[4] = {"A", "C", "G", "T"};

	    printRatesDNA_BIN(4, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(4, f, freqNames);
	  }
	  break;
	case SECONDARY_DATA_6:
	   {
	    char *freqNames[6] = {"AU", "CG", "GC", "GU", "UA", "UG"};

	    printRatesRest(6, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(6, f, freqNames);
	  }
	  break;
	case SECONDARY_DATA_7:
	  {
	    char *freqNames[7] = {"AU", "CG", "GC", "GU", "UA", "UG", "REST"};

	    printRatesRest(7, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(7, f, freqNames);
	  }
	  break;
	case SECONDARY_DATA:
	  {
	    char *freqNames[16] = {"AA", "AC", "AG", "AU", "CA", "CC", "CG", "CU",
				   "GA", "GC", "GG", "GU", "UA", "UC", "UG", "UU"};

	    printRatesRest(16, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(16, f, freqNames);
	  }
	  break;
	case BINARY_DATA:
	  {
	    char *freqNames[2] = {"0", "1"};

	    printRatesDNA_BIN(2, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(2, f, freqNames);
	  }
	  break;
	default:
	  assert(0);
	}

      printBothOpen("\n");
    }
}

static void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(processID == 0)
    {
      double t;

      t = gettime() - masterTime;

      switch(adef->mode)
	{
	case MESH_TREE_SEARCH:
	  break;
	case TREE_EVALUATION :
	  printBothOpen("\n\nOverall Time for Tree Evaluation %f\n", t);
	  printBothOpen("Final GAMMA  likelihood: %f\n", tr->likelihood);

	  {
	    boolean
	      linkedProteinGTR = FALSE;
	    
	    int
	      model,
	      params = 0,
	      paramsBrLen = 0;

	    for(model = 0; model < tr->NumberOfModels; model++)
	      {
		switch(tr->partitionData[model].dataType)
		  {
		  case AA_DATA:	 
		    if(tr->partitionData[model].protModels == GTR_UNLINKED)
		      params += 189;
		    
		    if(tr->partitionData[model].protModels == GTR)
		      linkedProteinGTR = TRUE;
		    
		    if(!tr->partitionData[model].usePredefinedProtFreqs)
		      params += 19;
		    break;
		  case GENERIC_32:
		    {
		      int 
			states = tr->partitionData[model].states;

		      /* frequencies */
		      
		      params += (states - 1);

		      switch(tr->multiStateModel)
			{
			case ORDERED_MULTI_STATE:			 
			  break;
			case MK_MULTI_STATE:
			  params += (states - 1);
			  break;
			case GTR_MULTI_STATE:
			  params += ((((states * states) - states) / 2) - 1);
			  break;
			default:
			  assert(0);
			}
		    break;
		  case GENERIC_64:
		    assert(0);
		    break;
		  case DNA_DATA:
		    params += 5 + 3;
		    break;
		  case SECONDARY_DATA_6:	  	      
		  case SECONDARY_DATA_7:	 		 
		  case SECONDARY_DATA: 
		    {
		      int 
			states = tr->partitionData[model].states;
		      
		      switch(tr->secondaryStructureModel)
			{
			case SEC_6_A:
			  params += ((((states * states) - states) / 2) - 1); /*rates*/
			  params += (states - 1); /* frequencies */
			  break;
			case SEC_6_B:
			  params += 1; /*rates */
			  params += 5; /* frequencies */	     
			  break;
			case SEC_6_C:
			  params += 1; /*rates */
			  params += 2; /* frequencies */	     
			  break;
			case SEC_6_D:
			  params += 1; /*rates */
			  params += 1; /* frequencies */	     
			  break;
			case SEC_6_E:
			  params += 1; /*rates */
			  params += 5; /* frequencies */	      
			  break;
			case SEC_7_A:
			  params += ((((states * states) - states) / 2) - 1); /*rates*/
			  params += (states - 1); /* frequencies */		
			  break;
			case SEC_7_B:
			  params += 20; /*rates */
			  params += 3; /* frequencies */		
			  break;
			case SEC_7_C:
			  params += 9; /*rates */
			  params += 6; /* frequencies */	     
			  break;
			case SEC_7_D:	
			  params += 3; /*rates */
			  params += 6; /* frequencies */	     
			  break;	      	   
			case SEC_7_E:
			  params += 1; /*rates */
			  params += 6; /* frequencies */	     
			  break;				  
			case SEC_7_F:
			  params += 3; /*rates */
			  params += 3; /* frequencies */	     
			  break;				     
			case SEC_16:
			  params += ((((states * states) - states) / 2) - 1); /*rates*/
			  params += (states - 1); /* frequencies */
			  break;
			case SEC_16_A:	
			  params += 4; /*rates */
			  params += 15; /* frequencies */	      
			  break;
			case SEC_16_B:
			  params += 0; /*rates */
			  params += 15; /* frequencies */	      
			  break;	     	    
			case SEC_16_C:	      
			case SEC_16_D:
			case SEC_16_E:
			case SEC_16_F:
			case SEC_16_I:
			case SEC_16_J:
			case SEC_16_K:
			  assert(0);
			default:
			  assert(0);
			}	 
		    }
		    break;
		    case BINARY_DATA:
		      params += 1;
		      break;
		    default:
		      assert(0);
		    }
		  }
		
		if(adef->useInvariant)
		  params += 2;
		else /* GAMMA */
		  params += 1;
	      }
	    
	    if(linkedProteinGTR)
	      params += 189;

	    if(tr->multiBranch)
	      paramsBrLen = params + tr->NumberOfModels * (2 * tr->mxtips - 3);
	    else
	      paramsBrLen = params + 2 * tr->mxtips - 3;

	    printBothOpen("\n");

	   
	    printBothOpen("Number of free parameters for AIC-TEST(BR-LEN): %d\n",    paramsBrLen);
	    printBothOpen("Number of free parameters for AIC-TEST(NO-BR-LEN): %d\n", params);
	    
	    
	    printBothOpen("\n\n");
	    
	    printModelParams(tr, adef);
	    
	    printBothOpen("Final tree written to:                 %s\n", resultFileName);
	    printBothOpen("Execution Log File written to:         %s\n", logFileName);
	 
	  }
	  break;
	case  BIG_RAPID_MODE:
	  if(adef->boot)
	    {
	      printBothOpen("\n\nOverall Time for %d Bootstraps %f\n", adef->multipleRuns, t);
	      printBothOpen("\n\nAverage Time per Bootstrap %f\n", (double)(t/((double)adef->multipleRuns)));
	      printBothOpen("All %d bootstrapped trees written to: %s\n", adef->multipleRuns, bootstrapFileName);
	    }
	  else
	    {
	      if(adef->multipleRuns > 1)
		{
		  double avgLH = 0;
		  double bestLH = unlikely;
		  int i, bestI  = 0;

		  for(i = 0; i < adef->multipleRuns; i++)
		    {
		      avgLH   += tr->likelihoods[i];
		      if(tr->likelihoods[i] > bestLH)
			{
			  bestLH = tr->likelihoods[i];
			  bestI  = i;
			}
		    }
		  avgLH /= ((double)adef->multipleRuns);

		  printBothOpen("\n\nOverall Time for %d Inferences %f\n", adef->multipleRuns, t);
		  printBothOpen("Average Time per Inference %f\n", (double)(t/((double)adef->multipleRuns)));
		  printBothOpen("Average Likelihood   : %f\n", avgLH);
		  printBothOpen("\n");
		  printBothOpen("Best Likelihood in run number %d: likelihood %f\n\n", bestI, bestLH);

		  if(adef->checkpoints)
		    printBothOpen("Checkpoints written to:                 %s.RUN.%d.* to %d.*\n", checkpointFileName, 0, adef->multipleRuns - 1);
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			printBothOpen("Random starting trees written to:       %s.RUN.%d to %d\n", randomFileName, 0, adef->multipleRuns - 1);
		      else
			printBothOpen("Parsimony starting trees written to:    %s.RUN.%d to %d\n", permFileName, 0, adef->multipleRuns - 1);
		    }
		  printBothOpen("Final trees written to:                 %s.RUN.%d to %d\n", resultFileName,  0, adef->multipleRuns - 1);
		  printBothOpen("Execution Log Files written to:         %s.RUN.%d to %d\n", logFileName, 0, adef->multipleRuns - 1);
		  printBothOpen("Execution information file written to:  %s\n", infoFileName);
		}
	      else
		{
		  printBothOpen("\n\nOverall Time for 1 Inference %f\n", t);
		  printBothOpen("Likelihood   : %f\n", tr->likelihood);
		  printBothOpen("\n\n");

		  if(adef->checkpoints)
		  printBothOpen("Checkpoints written to:                %s.*\n", checkpointFileName);
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			printBothOpen("Random starting tree written to:       %s\n", randomFileName);
		      else
			printBothOpen("Parsimony starting tree written to:    %s\n", permFileName);
		    }
		  printBothOpen("Final tree written to:                 %s\n", resultFileName);
		  printBothOpen("Execution Log File written to:         %s\n", logFileName);
		  printBothOpen("Execution information file written to: %s\n",infoFileName);
		}
	    }

	  break;
	case CALC_BIPARTITIONS:
	  printBothOpen("\n\nTime for Computation of Bipartitions %f\n", t);
	  printBothOpen("Tree with bipartitions written to file:  %s\n", bipartitionsFileName);
	  printBothOpen("Tree with bipartitions as branch labels written to file:  %s\n", bipartitionsFileNameBranchLabels);	  
	  printBothOpen("Execution information file written to :  %s\n",infoFileName);
	  break;
	case PER_SITE_LL:
	  printBothOpen("\n\nTime for Optimization of per-site log likelihoods %f\n", t);
	  printBothOpen("Per-site Log Likelihoods written to File %s in Tree-Puzzle format\n",  perSiteLLsFileName);
	  printBothOpen("Execution information file written to :  %s\n",infoFileName);

	  break;
	case PARSIMONY_ADDITION:
	  printBothOpen("\n\nTime for MP stepwise addition %f\n", t);
	  printBothOpen("Execution information file written to :  %s\n",infoFileName);
	  printBothOpen("Complete parsimony tree written to:      %s\n", permFileName);
	  break;
	case ANCESTRAL_STATES:
	  printBothOpen("\n\nTime for marginal ancestral state computation: %f\n\n", t);
	  break;
	case QUARTET_CALCULATION:
	  printBothOpen("\n\nOverall Time for quartet computation: %f\n\n", t);
	  break;
	case THOROUGH_OPTIMIZATION:
	  printBothOpen("\n\nTime for thorough tree optimization: %f\n\n", t);
	  break;
	case OPTIMIZE_BR_LEN_SCALER:
	  printBothOpen("\n\nTime for branch length scaler and remaining model parameters optimization: %f\n\n", t);
	  break;
	default:
	  assert(0);
	}
    }

}


/************************************************************************************/


#ifdef _USE_PTHREADS






static void computeFraction(tree *localTree, int tid, int n)
{
  int
    i,
    model;

  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      int width = 0;

      for(i = localTree->partitionData[model].lower; i < localTree->partitionData[model].upper; i++)
	if(i % n == tid)
	      width++;

      localTree->partitionData[model].width = width;
    }
}



static void threadFixModelIndices(tree *tr, tree *localTree, int tid, int n)
{
  size_t
    model,
    j,
    i,
    globalCounter = 0,
    localCounter  = 0,
    offset,
    countOffset,
    myLength = 0,
    memoryRequirements = 0;

  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      localTree->partitionData[model].lower      = tr->partitionData[model].lower;
      localTree->partitionData[model].upper      = tr->partitionData[model].upper;
    }

  computeFraction(localTree, tid, n);

  for(model = 0, offset = 0, countOffset = 0; model < (size_t)localTree->NumberOfModels; model++)
    {           
      localTree->partitionData[model].sumBuffer    = &localTree->sumBuffer[offset];      
      localTree->partitionData[model].perSiteLL    = &localTree->perSiteLLPtr[countOffset];          
      localTree->partitionData[model].wgt          = &localTree->wgtPtr[countOffset];
      localTree->partitionData[model].invariant    = &localTree->invariantPtr[countOffset];
      localTree->partitionData[model].rateCategory = &localTree->rateCategoryPtr[countOffset];     

      countOffset += localTree->partitionData[model].width;

      offset += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * (size_t)(localTree->partitionData[model].width);      
    }

  myLength           = countOffset;
  memoryRequirements = offset;


  /* figure in data */   

  for(i = 0; i < (size_t)localTree->mxtips; i++)
    {
      for(model = 0, offset = 0, countOffset = 0; model < (size_t)localTree->NumberOfModels; model++)
	{
	  localTree->partitionData[model].yVector[i+1]   = &localTree->y_ptr[i * myLength + countOffset];
	  countOffset +=  localTree->partitionData[model].width;
	}
      assert(countOffset == myLength);
    }

 

  for(model = 0, globalCounter = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      for(localCounter = 0, i = (size_t)localTree->partitionData[model].lower;  i < (size_t)localTree->partitionData[model].upper; i++)
	{
	  if(i % (size_t)n == (size_t)tid)
	    {
	      localTree->partitionData[model].wgt[localCounter]          = tr->cdta->aliaswgt[globalCounter];	      
	      localTree->partitionData[model].invariant[localCounter]    = tr->invariant[globalCounter];
	      localTree->partitionData[model].rateCategory[localCounter] = tr->cdta->rateCategory[globalCounter];	      

	      for(j = 1; j <= (size_t)localTree->mxtips; j++)
	       localTree->partitionData[model].yVector[j][localCounter] = tr->yVector[j][globalCounter]; 	     

	      localCounter++;
	    }
	  globalCounter++;
	}
    }
  
  for(model = 0; model < (size_t)localTree->NumberOfModels; model++)
    {
      int        
	undetermined = getUndetermined(localTree->partitionData[model].dataType);
      
      size_t
	width =  localTree->partitionData[model].width;
      
      localTree->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
      
      memset(localTree->partitionData[model].gapVector, 0, localTree->partitionData[model].initialGapVectorSize);

      for(j = 1; j <= (size_t)(localTree->mxtips); j++)
	for(i = 0; i < width; i++)
	  if(localTree->partitionData[model].yVector[j][i] == undetermined)
	    localTree->partitionData[model].gapVector[localTree->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];
    }
}


static void initPartition(tree *tr, tree *localTree, int tid)
{
  int model;

  localTree->threadID = tid; 

  if(tid > 0)
    {
      int totalLength = 0;
      
      localTree->useGammaMedian          = tr->useGammaMedian;
      localTree->saveMemory              = tr->saveMemory;
      localTree->useGappedImplementation = tr->useGappedImplementation;
      localTree->innerNodes              = tr->innerNodes;
      localTree->useFastScaling          = tr->useFastScaling;
      localTree->maxCategories           = tr->maxCategories;
     
      localTree->originalCrunchedLength  = tr->originalCrunchedLength;
      localTree->NumberOfModels          = tr->NumberOfModels;
      localTree->mxtips                  = tr->mxtips;
      localTree->multiBranch             = tr->multiBranch;
      localTree->multiGene               = tr->multiGene;
      assert(localTree->multiGene == 0);
      localTree->numBranches             = tr->numBranches;
      localTree->lhs                     = (double*)malloc(sizeof(double)   * localTree->originalCrunchedLength);
      localTree->executeModel            = (boolean*)malloc(sizeof(boolean) * localTree->NumberOfModels);
      localTree->perPartitionLH          = (double*)malloc(sizeof(double)   * localTree->NumberOfModels);
      localTree->storedPerPartitionLH    = (double*)malloc(sizeof(double)   * localTree->NumberOfModels);

      localTree->fracchanges = (double*)malloc(sizeof(double)   * localTree->NumberOfModels);
      localTree->partitionContributions = (double*)malloc(sizeof(double)   * localTree->NumberOfModels);

      localTree->partitionData = (pInfo*)malloc(sizeof(pInfo) * localTree->NumberOfModels);

      /* extend for multi-branch */
      localTree->td[0].count = 0;
      localTree->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * localTree->mxtips);

      localTree->cdta               = (cruncheddata*)malloc(sizeof(cruncheddata));
      localTree->cdta->patrat       = (double*)malloc(sizeof(double) * localTree->originalCrunchedLength);
      localTree->cdta->patratStored = (double*)malloc(sizeof(double) * localTree->originalCrunchedLength);      

      localTree->discreteRateCategories = tr->discreteRateCategories;     

      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  localTree->partitionData[model].numberOfCategories    = tr->partitionData[model].numberOfCategories;
	  localTree->partitionData[model].states     = tr->partitionData[model].states;
	  localTree->partitionData[model].maxTipStates    = tr->partitionData[model].maxTipStates;
	  localTree->partitionData[model].dataType   = tr->partitionData[model].dataType;
	  localTree->partitionData[model].protModels = tr->partitionData[model].protModels;
	  localTree->partitionData[model].usePredefinedProtFreqs  = tr->partitionData[model].usePredefinedProtFreqs;
	  localTree->partitionData[model].mxtips     = tr->partitionData[model].mxtips;
	  localTree->partitionData[model].lower      = tr->partitionData[model].lower;
	  localTree->partitionData[model].upper      = tr->partitionData[model].upper;
	  localTree->executeModel[model]             = TRUE;
	  localTree->perPartitionLH[model]           = 0.0;
	  localTree->storedPerPartitionLH[model]     = 0.0;
	  totalLength += (localTree->partitionData[model].upper -  localTree->partitionData[model].lower);		  
	}

      assert(totalLength == localTree->originalCrunchedLength);
    }

  for(model = 0; model < localTree->NumberOfModels; model++)
    localTree->partitionData[model].width        = 0;
}




static void allocNodex(tree *tr, int tid, int n)
{
  size_t   
    model,
    memoryRequirements = 0,
    myLength = 0;

  computeFraction(tr, tid, n);

  allocPartitions(tr);

  for(model = 0; model < (size_t)tr->NumberOfModels; model++)
    {
      size_t 
	width = tr->partitionData[model].width;

      int 
	i;

      myLength += width;

      memoryRequirements += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * width;
     
      tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
      
      tr->partitionData[model].gapVector = (unsigned int*)calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));
      
      tr->partitionData[model].initialGapVectorSize = tr->partitionData[model].gapVectorLength * 2 * tr->mxtips * sizeof(int);
      
      /* always multiply by 4 due to frequent switching between CAT and GAMMA in standard RAxML */

      tr->partitionData[model].gapColumn = (double *)malloc_aligned(
								    ((size_t)(tr->innerNodes)) *
								    ((size_t)(4)) * 
								    ((size_t)(tr->partitionData[model].states)) *
								    sizeof(double));		             
      for(i = 0; i < tr->innerNodes; i++)
	{
	  tr->partitionData[model].xVector[i]   = (double*)NULL;     
	  tr->partitionData[model].expVector[i]   = (int*)NULL;
	}
    }

  if(tid == 0)
    {
      tr->perSiteLL       = (double *)malloc((size_t)tr->cdta->endsite * sizeof(double));
      assert(tr->perSiteLL != NULL);
    }
  
  tr->sumBuffer  = (double *)malloc_aligned(memoryRequirements * sizeof(double));
  assert(tr->sumBuffer != NULL);
   
  tr->y_ptr = (unsigned char *)malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);  

  tr->perSiteLLPtr     = (double*) malloc(myLength * sizeof(double));
  assert(tr->perSiteLLPtr != NULL);

  tr->wgtPtr           = (int*)    malloc(myLength * sizeof(int));
  assert(tr->wgtPtr != NULL);  

  tr->invariantPtr     = (int*)    malloc(myLength * sizeof(int));
  assert(tr->invariantPtr != NULL);

  tr->rateCategoryPtr  = (int*)    malloc(myLength * sizeof(int));
  assert(tr->rateCategoryPtr != NULL);
}






inline static void sendTraversalInfo(tree *localTree, tree *tr)
{
  /* the one below is a hack we are re-assigning the local pointer to the global one
     the memcpy version below is just for testing and preparing the
     fine-grained MPI BlueGene version */

  if(1)
    {     
      localTree->td[0] = tr->td[0];
    }
  else
    {
      localTree->td[0].count = tr->td[0].count;
      memcpy(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));
    }
}


static void collectDouble(double *dst, double *src, tree *tr, int n, int tid)
{
  int model, i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	{
	  if(i % n == tid)
	    dst[i] = src[i];
	}
    }
}


static void broadcastPerSiteRates(tree *tr, tree *localTree)
{
  int
    i = 0,
    model = 0;

  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      localTree->partitionData[model].numberOfCategories = tr->partitionData[model].numberOfCategories;

      for(i = 0; i < localTree->partitionData[model].numberOfCategories; i++)
	{
	  localTree->partitionData[model].perSiteRates[i] = tr->partitionData[model].perSiteRates[i];
	  localTree->partitionData[model].unscaled_perSiteRates[i] = tr->partitionData[model].unscaled_perSiteRates[i];
	}
    }

}

static void execFunction(tree *tr, tree *localTree, int tid, int n)
{
  double volatile result;
  int
    i,
    currentJob,
    parsimonyResult,
    model,
    localCounter,
    globalCounter;

  currentJob = threadJob >> 16;

  switch(currentJob)
    {     
    case THREAD_INIT_PARTITION:
      initPartition(tr, localTree, tid);
      break;
    case THREAD_ALLOC_LIKELIHOOD:
      allocNodex(localTree, tid, n);
      threadFixModelIndices(tr, localTree, tid, n);
      break;
    case THREAD_FIX_MODEL_INDICES:
      threadFixModelIndices(tr, localTree, tid, n);
      break;    
    case THREAD_EVALUATE:
      sendTraversalInfo(localTree, tr);
      result = evaluateIterative(localTree, FALSE);

      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_NEWVIEW_MASKED:
      sendTraversalInfo(localTree, tr);
      memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);
      newviewIterative(localTree);
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_NEWVIEW:
      sendTraversalInfo(localTree, tr);
      newviewIterative(localTree);
      break;
    case THREAD_MAKENEWZ_FIRST:
      {
	volatile double
	  dlnLdlz[NUM_BRANCHES],
	  d2lnLdlz2[NUM_BRANCHES];

	sendTraversalInfo(localTree, tr);
	if(tid > 0)
	  {
	    memcpy(localTree->coreLZ,   tr->coreLZ,   sizeof(double) *  localTree->numBranches);
	    memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);
	  }

	makenewzIterative(localTree);	
	execCore(localTree, dlnLdlz, d2lnLdlz2);

	if(!tr->multiBranch)
	  {
	    reductionBuffer[tid]    = dlnLdlz[0];
	    reductionBufferTwo[tid] = d2lnLdlz2[0];
	  }
	else
	  {
	    for(i = 0; i < localTree->NumberOfModels; i++)
	      {
		reductionBuffer[tid * localTree->NumberOfModels + i]    = dlnLdlz[i];
		reductionBufferTwo[tid * localTree->NumberOfModels + i] = d2lnLdlz2[i];
	      }
	  }

	if(tid > 0)
	  {
	    for(model = 0; model < localTree->NumberOfModels; model++)
	      localTree->executeModel[model] = TRUE;
	  }
      }
      break;
    case THREAD_MAKENEWZ:
      {
	volatile double
	  dlnLdlz[NUM_BRANCHES],
	  d2lnLdlz2[NUM_BRANCHES];

	memcpy(localTree->coreLZ,   tr->coreLZ,   sizeof(double) *  localTree->numBranches);
	memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);
	
	execCore(localTree, dlnLdlz, d2lnLdlz2);

	if(!tr->multiBranch)
	  {
	    reductionBuffer[tid]    = dlnLdlz[0];
	    reductionBufferTwo[tid] = d2lnLdlz2[0];
	  }
	else
	  {
	    for(i = 0; i < localTree->NumberOfModels; i++)
	      {
		reductionBuffer[tid * localTree->NumberOfModels + i]    = dlnLdlz[i];
		reductionBufferTwo[tid * localTree->NumberOfModels + i] = d2lnLdlz2[i];
	      }
	  }
	if(tid > 0)
	  {
	    for(model = 0; model < localTree->NumberOfModels; model++)
	      localTree->executeModel[model] = TRUE;
	  }
      }
      break;
    case THREAD_COPY_RATES:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {	      
	      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
	      
	      memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        pl->eignLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          pl->evLength * sizeof(double));		  
	      memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          pl->eiLength * sizeof(double));
	      memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   pl->tipVectorLength * sizeof(double));
	      
	     	
	    }
	}
      break;
    case THREAD_OPT_RATE:
      if(tid > 0)
	{
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));

	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
	      
	      memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        pl->eignLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          pl->evLength * sizeof(double));		  
	      memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          pl->eiLength * sizeof(double));
	      memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   pl->tipVectorLength * sizeof(double));
	      
	      	     
	    }
	}

      result = evaluateIterative(localTree, FALSE);


      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;


      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_COPY_INVAR:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	}
      break;
    case THREAD_OPT_INVAR:
      if(tid > 0)
	{
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	}

      result = evaluateIterative(localTree, FALSE);

      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_COPY_ALPHA:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	    }
	}
      break;
    case THREAD_OPT_ALPHA:
      if(tid > 0)
	{
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	}

      result = evaluateIterative(localTree, FALSE);


      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid *  localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    case THREAD_RESET_MODEL:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));

	      memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        pl->eignLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          pl->evLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          pl->eiLength * sizeof(double));
	      memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  pl->substRatesLength * sizeof(double));
	      memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, pl->frequenciesLength * sizeof(double));
	      memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   pl->tipVectorLength * sizeof(double));
	      
	       	      	      

	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	      localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	    }
	}
      break;     
    case THREAD_COPY_INIT_MODEL:
      if(tid > 0)
	{
	  localTree->rateHetModel       = tr->rateHetModel;

	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));

	      memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        pl->eignLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          pl->evLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          pl->eiLength * sizeof(double));
	      memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  pl->substRatesLength * sizeof(double));
	      memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, pl->frequenciesLength * sizeof(double));
	      memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   pl->tipVectorLength * sizeof(double));
	      
	          

	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	      localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	      localTree->partitionData[model].lower      = tr->partitionData[model].lower;
	      localTree->partitionData[model].upper      = tr->partitionData[model].upper;
	      
	      localTree->partitionData[model].numberOfCategories      = tr->partitionData[model].numberOfCategories;
	    }

	  memcpy(localTree->cdta->patrat,        tr->cdta->patrat,      localTree->originalCrunchedLength * sizeof(double));
	  memcpy(localTree->cdta->patratStored, tr->cdta->patratStored, localTree->originalCrunchedLength * sizeof(double));	  
	}     

       for(model = 0; model < localTree->NumberOfModels; model++)
	 {
	   int localIndex;
	   for(i = localTree->partitionData[model].lower, localIndex = 0; i <  localTree->partitionData[model].upper; i++)
	     if(i % n == tid)
	       {
		 localTree->partitionData[model].wgt[localIndex]          = tr->cdta->aliaswgt[i];				 		
		 localTree->partitionData[model].invariant[localIndex]    = tr->invariant[i];		

		 localIndex++;
	       }	  
	 }
      break;    
    case THREAD_RATE_CATS:
      sendTraversalInfo(localTree, tr);
      if(tid > 0)
	{
	  localTree->lower_spacing = tr->lower_spacing;
	  localTree->upper_spacing = tr->upper_spacing;
	}

      optRateCatPthreads(localTree, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);

      if(tid > 0)
	{
	  collectDouble(tr->cdta->patrat,       localTree->cdta->patrat,         localTree, n, tid);
	  collectDouble(tr->cdta->patratStored, localTree->cdta->patratStored,   localTree, n, tid);
	  collectDouble(tr->lhs,                localTree->lhs,                  localTree, n, tid);
	}
      break;
    case THREAD_COPY_RATE_CATS:
      if(tid > 0)
	{	
	  memcpy(localTree->cdta->patrat,       tr->cdta->patrat,         localTree->originalCrunchedLength * sizeof(double));
	  memcpy(localTree->cdta->patratStored, tr->cdta->patratStored,   localTree->originalCrunchedLength * sizeof(double));
	  broadcastPerSiteRates(tr, localTree);
	}

      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  localTree->partitionData[model].numberOfCategories = tr->partitionData[model].numberOfCategories;

	  for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++)
	    {
	      if(i % n == tid)
		{		 
		  localTree->partitionData[model].rateCategory[localCounter] = tr->cdta->rateCategory[i];
		  localTree->partitionData[model].wr[localCounter]             = tr->wr[i];
		  localTree->partitionData[model].wr2[localCounter]            = tr->wr2[i];

		
		 
		  localCounter++;
		}
	    }
	}
      break;
    case THREAD_CAT_TO_GAMMA:
      if(tid > 0)
	localTree->rateHetModel = tr->rateHetModel;
      break;
    case THREAD_GAMMA_TO_CAT:
      if(tid > 0)
	localTree->rateHetModel = tr->rateHetModel;
      break;
    case THREAD_EVALUATE_VECTOR:
      sendTraversalInfo(localTree, tr);
      result = evaluateIterative(localTree, TRUE);            
     
      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid * localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}

      for(model = 0, globalCounter = 0; model < localTree->NumberOfModels; model++)
	{
	  for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++)
	    {
	      if(i % n == tid)
		{
		  tr->perSiteLL[globalCounter] =  localTree->partitionData[model].perSiteLL[localCounter];
		  localCounter++;
		}
	      globalCounter++;
	    }
	}
      break;    
    case THREAD_COPY_PARAMS:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[model]));
	      
	      memcpy(localTree->partitionData[model].EIGN,        tr->partitionData[model].EIGN,        pl->eignLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EV,          tr->partitionData[model].EV,          pl->evLength * sizeof(double));
	      memcpy(localTree->partitionData[model].EI,          tr->partitionData[model].EI,          pl->eiLength * sizeof(double));
	      memcpy(localTree->partitionData[model].substRates,  tr->partitionData[model].substRates,  pl->substRatesLength * sizeof(double));
	      memcpy(localTree->partitionData[model].frequencies, tr->partitionData[model].frequencies, pl->frequenciesLength * sizeof(double));
	      memcpy(localTree->partitionData[model].tipVector,   tr->partitionData[model].tipVector,   pl->tipVectorLength * sizeof(double));
	      
	     	     
	    }
	}
      break;
    case THREAD_INIT_EPA:     
      if(tid > 0)
	{
	  localTree->leftRootNode             = tr->leftRootNode;
	  localTree->rightRootNode            = tr->rightRootNode;
	  localTree->wasRooted                = tr->wasRooted;
	  localTree->bInf                     = tr->bInf;	 
	  localTree->numberOfBranches         = tr->numberOfBranches;
	  localTree->contiguousVectorLength   = tr->contiguousVectorLength;
	  localTree->contiguousScalingLength  = tr->contiguousScalingLength;
	  localTree->inserts                  = tr->inserts;
	  localTree->numberOfTipsForInsertion = tr->numberOfTipsForInsertion;	
	  localTree->fracchange = tr->fracchange;
	  memcpy(localTree->partitionContributions, tr->partitionContributions, sizeof(double) * localTree->NumberOfModels);
	  memcpy(localTree->fracchanges, tr->fracchanges, sizeof(double) * localTree->NumberOfModels);	 
	}                                                

      localTree->temporarySumBuffer = (double *)malloc_aligned(sizeof(double) * localTree->contiguousVectorLength);
      localTree->temporaryVector  = (double *)malloc_aligned(sizeof(double) * localTree->contiguousVectorLength);      

      localTree->temporaryScaling = (int *)malloc(sizeof(int) * localTree->contiguousScalingLength);
                 
      
      localTree->contiguousWgt          = (int*)malloc(sizeof(int) * localTree->contiguousScalingLength);
      localTree->contiguousInvariant    = (int*)malloc(sizeof(int) * localTree->contiguousScalingLength);	  
      
     
      memcpy(localTree->contiguousWgt         , tr->cdta->aliaswgt,     sizeof(int) * localTree->contiguousScalingLength);
      memcpy(localTree->contiguousInvariant   , tr->invariant,          sizeof(int) * localTree->contiguousScalingLength);
      
      if(tid > 0)
	broadcastPerSiteRates(tr, localTree);

      localTree->contiguousWR           = (double*)malloc(sizeof(double) * localTree->contiguousScalingLength);
      localTree->contiguousWR2          = (double*)malloc(sizeof(double) * localTree->contiguousScalingLength);
      localTree->contiguousRateCategory = (int*)malloc(sizeof(int) * localTree->contiguousScalingLength);
      
      memcpy(localTree->contiguousWR, tr->wr, sizeof(double) * localTree->contiguousScalingLength);
      memcpy(localTree->contiguousWR2, tr->wr2, sizeof(double) * localTree->contiguousScalingLength);
      memcpy(localTree->contiguousRateCategory, tr->cdta->rateCategory, sizeof(int) * localTree->contiguousScalingLength);           
     
      localTree->contiguousTips = tr->yVector;	  	
	 
      break;         
    case THREAD_GATHER_LIKELIHOOD:
      {	
	int 
	  branchCounter = tr->branchCounter;

	double
	  *leftContigousVector = localTree->bInf[branchCounter].epa->left,
	  *rightContigousVector = localTree->bInf[branchCounter].epa->right;
      
	int
	  *leftContigousScalingVector = localTree->bInf[branchCounter].epa->leftScaling,
	  *rightContigousScalingVector = localTree->bInf[branchCounter].epa->rightScaling,	
	  globalColumnCount = 0,
	  globalCount       = 0,
	  rightNumber = localTree->bInf[branchCounter].epa->rightNodeNumber,
	  leftNumber  = localTree->bInf[branchCounter].epa->leftNodeNumber;	


	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    size_t
	      blockRequirements;

	    double
	      *leftStridedVector  =  (double *)NULL,
	      *rightStridedVector =  (double *)NULL;

	    int
	      *leftStridedScalingVector  =  (int *)NULL,
	      *rightStridedScalingVector =  (int *)NULL,
	     
	      localColumnCount = 0,
	      localCount = 0;	   

	    if(!isTip(leftNumber, localTree->mxtips))
	      {
		leftStridedVector        = localTree->partitionData[model].xVector[leftNumber - localTree->mxtips - 1];
		leftStridedScalingVector = localTree->partitionData[model].expVector[leftNumber - localTree->mxtips - 1];
	      }	   

	    if(!isTip(rightNumber, localTree->mxtips))
	      {
		rightStridedVector        = localTree->partitionData[model].xVector[rightNumber - localTree->mxtips - 1];
		rightStridedScalingVector = localTree->partitionData[model].expVector[rightNumber - localTree->mxtips - 1];
	      }	    

	    assert(!(isTip(leftNumber, localTree->mxtips) && isTip(rightNumber, localTree->mxtips)));	   

	    blockRequirements = (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states);	   

	    for(globalColumnCount = localTree->partitionData[model].lower; globalColumnCount < localTree->partitionData[model].upper; globalColumnCount++)
	      {	
		if(globalColumnCount % n == tid)
		  {		    
		    if(leftStridedVector)
		      {
			memcpy(&leftContigousVector[globalCount], &leftStridedVector[localCount], sizeof(double) * blockRequirements);
			leftContigousScalingVector[globalColumnCount] = leftStridedScalingVector[localColumnCount];
		      }
		    
		    if(rightStridedVector)
		      {
			memcpy(&rightContigousVector[globalCount], &rightStridedVector[localCount], sizeof(double) * blockRequirements);
			rightContigousScalingVector[globalColumnCount] = rightStridedScalingVector[localColumnCount];
		      }
		   
		    localColumnCount++;
		    localCount += blockRequirements;
		  }

	

		globalCount += blockRequirements;
	      }	    

	    assert(localColumnCount == localTree->partitionData[model].width);
	    assert(localCount == (localTree->partitionData[model].width * (int)blockRequirements));

	  }
      }
      break;                 
    case THREAD_INSERT_CLASSIFY:          
    case THREAD_INSERT_CLASSIFY_THOROUGH:       
      { 
	int
	  branchNumber;
	
	boolean 
	  done = FALSE;	       

	while(!done)
	  {	      	      	    		    
	    pthread_mutex_lock(&mutex);
	      
	    if(NumberOfJobs == 0)
	      done = TRUE;
	    else
	      {		  
		branchNumber = localTree->numberOfBranches - NumberOfJobs;		 
		NumberOfJobs--;		 
	      }
	      	   
	    pthread_mutex_unlock(&mutex);
	      
	    if(!done)
	      {		 		 		 		 	 		      
		switch(currentJob)
		  {
		  case THREAD_INSERT_CLASSIFY:
		    addTraverseRobIterative(localTree, branchNumber);
		    break;		  
		  case  THREAD_INSERT_CLASSIFY_THOROUGH:		    
		    testInsertThoroughIterative(localTree, branchNumber);		   
		    break;   		 		 
		  default:
		    assert(0);
		  }
		  		  		
	      }	    
	    }
      }
      break;
    case  THREAD_USE_GAPPED:
      localTree->useGappedImplementation = tr->useGappedImplementation;
      break;     
    case THREAD_PREPARE_BIPS_FOR_PRINT:
      {       
	int 
	  i, 
	  j;	
	
	boolean 
	  done = FALSE;

	while(!done)
	  {
	    pthread_mutex_lock(&mutex);
	    
	    if(NumberOfJobs == 0)
	      done = TRUE;
	    else
	      {
		i = tr->consensusBipLen - NumberOfJobs;
		NumberOfJobs--;
	      }
	    
	    pthread_mutex_unlock(&mutex);

	    if( ! done)	      
	      {
		for(j = i+1; j < tr->consensusBipLen; j++)
		  {
		    if(tr->consensusBips[i]->amountTips < tr->consensusBips[j]->amountTips &&
		       issubset(tr->consensusBips[i]->bitVector, tr->consensusBips[j]->bitVector, tr->bitVectorLength))
		      { 
			/* i is child of j */		    
			List *elem = (List*) malloc(sizeof(List));
			
			elem->value = calloc(1, sizeof(int));
			
			*(int*)elem->value = i;
			
			pthread_mutex_lock(tr->mutexesForHashing[j]); /* LOCKED */
			
			tr->hasAncestor[i] = TRUE;

			elem->next = tr->listOfDirectChildren[j];
			tr->listOfDirectChildren[j] = elem;
			
			pthread_mutex_unlock(tr->mutexesForHashing[j]); /* UNLOCKED */
			
			break;	/* each node has only 1 parent -> nothing more to do */
		      }
		  }
	      }
	  }       
      }
      break;
    case THREAD_MRE_COMPUTE: 	
       {	
	if(tid > 0)
	  {		
	    /* worker threads */	  	  
	    boolean done = FALSE; 	  	  
	    int localEntryCount = (int) tr->h->entryCount; /* problem? */
	    while(!done )
	      {
		int acquiredJobs =  0; 
		int jobId = -1;
		
		/* get new job */
		
		pthread_mutex_lock(&mutex) ; /* START LOCK */
		
		if( NumberOfJobs == 0 )
		  { 
		    /* finish  */
		    done = TRUE;
		  }
		else 
		  if( localEntryCount - NumberOfJobs + tr->recommendedAmountJobs  < tr->sectionEnd)
		    { 
		      /* try to acquire the recommended amount of jobs */
		      jobId = localEntryCount - NumberOfJobs;
		      acquiredJobs = tr->recommendedAmountJobs; 
		      NumberOfJobs -= acquiredJobs;
		    }
		  else 
		    if( localEntryCount - NumberOfJobs  < (signed int)tr->sectionEnd)
		      { 
			/* at least get one job */
			jobId = tr->h->entryCount - NumberOfJobs;
			acquiredJobs = 1; 
			NumberOfJobs--;
		      } 
		
		pthread_mutex_unlock(&mutex); /* END LOCK */
		
		if(*(tr->len) >= tr->maxBips)
		  break;

		/* check all  */
		while(acquiredJobs > 0)
		  {
		    boolean 
		      compatflag = TRUE;
		    
		    entry 
		      *currentEntry = tr->sbw[jobId];
		    
		    int k; 

		    if(!(tr->mr_thresh < currentEntry->supportFromTreeset[0])) 
		      {
			for(k = *(tr->len); k > 0; k--)
			  {
			    if(! compatible(tr->sbi[k-1], currentEntry, tr->bitVectorLength))
			      {
				compatflag = FALSE;
				break;
			      }
			  }	
		      }
		    if(compatflag) 
		      tr->bipStatus[jobId - tr->sectionEnd + tr->bipStatusLen] = MRE_POSSIBLE_CANDIDATE; /* ready to check */ 
		    else
		      tr->bipStatus[jobId - tr->sectionEnd + tr->bipStatusLen] = MRE_EXCLUDED; /* can be omitted */
		    
		    acquiredJobs--;
		    jobId++;
		  }
	      }
	  }      
	else			
	  /* master thread */
	  {	    
	    /* check in a looping manner, if bipartitions could be added */
	    
	    int 
	      highestToCheck,
	      tmpCounter = 0;	    
	    
	    double 
	      density = 0.0;
	    
	    while(TRUE)	      
	      {		
		/* get highest bip to check  */
		highestToCheck = 0;
		while(highestToCheck < tr->bipStatusLen)
		  { 
		    /* waits busily as long as there is nothing to do */
		    /* printf("%d is highest to check\n", highestToCheck); */
		    if( ! tr->bipStatus[highestToCheck] )
		      highestToCheck = 0; 
		    else 
		      if(tr->bipStatus[highestToCheck] == MRE_POSSIBLE_CANDIDATE)
			break;
		      else
			highestToCheck++;
		  } 
		
		/* try to finish */
		if( tmpCounter >= tr->maxBips || 
		    (highestToCheck == tr->bipStatusLen	/* end of buffer that is examined */
		     && tr->sectionEnd == tr->h->entryCount /* the end of the buffer is also the hashtable  */
		     && tr->bipStatus[highestToCheck-1] > MRE_POSSIBLE_CANDIDATE))
		  { 
		    /* the last entry in buffer was already processed */
		    *(tr->len) = tmpCounter; /* for the workers to finish */
		    break;		   /* master says goodbye */
		  }

		/* reset section (resp. the buffer to be checked) */ 
		else 
		  if( highestToCheck == tr->bipStatusLen)
		    { 
		      int 
			newSectionEnd, 
			min, 
			max; 
		      
		      *(tr->len) = tmpCounter; /* reset counter for workers  */
		      tr->entriesOfSection = &(tr->sbw[tr->sectionEnd ]); 
		      
		      /* find new section end: tries to find a new window
			 size (and resp. sectionEnd) s.t. the expected
			 amount of work for master and workers is the same.
		      */		  
		      density /= tr->bipStatusLen; 

		      /* I am not entirely sure, if this makes the code really incredible faster...  */
		      max = 5 * (NumberOfThreads-1); 
		      min = 1; 
		      tr->recommendedAmountJobs = (int)(max + (min - max) * density); /* recommend an amount of jobs to be calculate per thread between min and max  */

		      if(density)
			{
			  int 
			    tmp = MAX((2 * tmpCounter * SECTION_CONSTANT / (NumberOfThreads * density)), /* the above discussed formula */
				      NumberOfThreads * MRE_MIN_AMOUNT_JOBS_PER_THREAD ); /* we need at least a bit work  */
			  newSectionEnd  = MIN(tr->sectionEnd + tmp, tr->h->entryCount);
			}
		      else
			newSectionEnd = tr->h->entryCount; 
		  
		      density = 0.0;

		      tr->bipStatusLen = newSectionEnd - tr->sectionEnd;
		      free(tr->bipStatus);
		      /* printf("%d\n" ,tr->bipStatusLen); */
		      tr->bipStatus = (int*)calloc(tr->bipStatusLen, sizeof(int));
		      tr->sectionEnd = newSectionEnd; 		  		  
		      continue;
		    }
		  		
		assert( tr->bipStatus[highestToCheck] == MRE_POSSIBLE_CANDIDATE);

		for(i = highestToCheck; i > 0; i--) /* checking new bip */
		  {	
		    assert(tr->bipStatus[i-1] == MRE_ADDED || tr->bipStatus[i-1] == MRE_EXCLUDED);
		    
		    if(tr->bipStatus[i-1] == MRE_ADDED 
		       && ! compatible(tr->entriesOfSection[i-1], tr->entriesOfSection[highestToCheck], tr->bitVectorLength))
		      {			
			tr->bipStatus[highestToCheck] = MRE_EXCLUDED;
			break;
		      }
		  } 
		
		if(i == 0)	/* accepting */
		  {
		    tr->bipStatus[highestToCheck] = MRE_ADDED;
		    tr->sbi[tmpCounter] = tr->entriesOfSection[highestToCheck];
		    tmpCounter++;
		    density++;
		  }
	      }
	  }     
       }      
      break;         
    case  THREAD_NEWVIEW_ANCESTRAL: 
      sendTraversalInfo(localTree, tr);
      newviewIterativeAncestral(localTree);
      break;
    case THREAD_GATHER_ANCESTRAL:      
      {
	double
	  *contigousVector = tr->ancestralStates;
      
	int	 
	  globalColumnCount = 0,
	  globalCount       = 0;	

	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    size_t
	      blockRequirements = (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states);
	    
	    int	     	     
	      localColumnCount = 0,
	      localCount = 0;	   

	    double 
	      *stridedVector = localTree->partitionData[model].sumBuffer;	    	   	    

	    for(globalColumnCount = localTree->partitionData[model].lower; globalColumnCount < localTree->partitionData[model].upper; globalColumnCount++)
	      {	
		if(globalColumnCount % n == tid)
		  {
		    memcpy(&contigousVector[globalCount], &stridedVector[localCount], sizeof(double) * blockRequirements);		
		    		   		   
		    localColumnCount++;
		    localCount += blockRequirements;
		  }	

		globalCount += blockRequirements;
	      }	    

	    assert(localColumnCount == localTree->partitionData[model].width);
	    assert(localCount == (localTree->partitionData[model].width * (int)blockRequirements));
	  }	
      }      
      break;
    case THREAD_OPT_SCALER:
      if(tid > 0)	
	memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));	 
	
      result = evaluateIterative(localTree, FALSE);

      if(localTree->NumberOfModels > 1)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    reductionBuffer[tid *  localTree->NumberOfModels + model] = localTree->perPartitionLH[model];
	}
      else
	reductionBuffer[tid] = result;

      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    localTree->executeModel[model] = TRUE;
	}
      break;
    default:
      printf("Job %d\n", currentJob);
      assert(0);
    }
}




void masterBarrier(int jobType, tree *tr)
{
  const int 
    n = NumberOfThreads;
  
  int 
    i, 
    sum;

  jobCycle = !jobCycle;
  threadJob = (jobType << 16) + jobCycle;

  execFunction(tr, tr, 0, n);

 
  do
    {
      for(i = 1, sum = 1; i < n; i++)
	sum += barrierBuffer[i];
    }
  while(sum < n);  

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
}

#ifndef _PORTABLE_PTHREADS

static void pinToCore(int tid)
{
  cpu_set_t cpuset;
         
  CPU_ZERO(&cpuset);    
  CPU_SET(tid, &cpuset);

  if(pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset) != 0)
    {
      printBothOpen("\n\nThere was a problem finding a physical core for thread number %d to run on.\n", tid);
      printBothOpen("Probably this happend because you are trying to run more threads than you have cores available,\n");
      printBothOpen("which is a thing you should never ever do again, good bye .... \n\n");
      assert(0);
    }
}

#endif

static void *likelihoodThread(void *tData)
{
  threadData *td = (threadData*)tData;
  tree
    *tr = td->tr,
    *localTree = (tree *)malloc(sizeof(tree));
  int
    myCycle = 0;

  const int 
    n = NumberOfThreads,
    tid = td->threadNumber;

#ifndef _PORTABLE_PTHREADS
  pinToCore(tid);
#endif
 
  printf("\nThis is RAxML Worker Pthread Number: %d\n", tid);

  while(1)
    {
      while (myCycle == threadJob);
      myCycle = threadJob;

      execFunction(tr, localTree, tid, n);

   
      barrierBuffer[tid] = 1;     
    }

  return (void*)NULL;
}

static void startPthreads(tree *tr)
{
  pthread_t *threads;
  pthread_attr_t attr;
  int rc, t;
  threadData *tData;

  jobCycle        = 0;
  threadJob       = 0;

  printf("\nThis is the RAxML Master Pthread\n");

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);

  pthread_mutex_init(&mutex , (pthread_mutexattr_t *)NULL);

  threads    = (pthread_t *)malloc(NumberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc(NumberOfThreads * sizeof(threadData));
  reductionBuffer          = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferTwo       = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferThree     = (volatile double *)malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferParsimony = (volatile int *)malloc(sizeof(volatile int) *  NumberOfThreads);

  
  barrierBuffer            = (volatile char *)malloc(sizeof(volatile char) *  NumberOfThreads);
  
  for(t = 0; t < NumberOfThreads; t++)
    barrierBuffer[t] = 0;

 
  branchInfos              = (volatile branchInfo **)malloc(sizeof(volatile branchInfo *) * NumberOfThreads);
 
  for(t = 1; t < NumberOfThreads; t++)
    {
      tData[t].tr  = tr;
      tData[t].threadNumber = t;
      rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t]));
      if(rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
    }
}



#endif


/*************************************************************************************************************************************************************/

static int elwCompare(const void *p1, const void *p2)
{
  elw *rc1 = (elw *)p1;
  elw *rc2 = (elw *)p2;

  double i = rc1->weight;
  double j = rc2->weight;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}

static int elwCompareLikelihood(const void *p1, const void *p2)
{
  elw *rc1 = (elw *)p1;
  elw *rc2 = (elw *)p2;

  double i = rc1->lh;
  double j = rc2->lh;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}

static void computeLHTest(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int
    i;
  
  double 
    bestLH, 
    currentLH, 
    weightSum = 0.0;
  
  FILE 
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);
  
  double 
    *bestVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);

  for(i = 0; i < tr->cdta->endsite; i++)
    weightSum += (double)(tr->cdta->aliaswgt[i]);

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);
  printBothOpen("Model optimization, best Tree: %f\n", tr->likelihood);
  bestLH = tr->likelihood;  

  evaluateGenericVector(tr, tr->start);
  memcpy(bestVector, tr->perSiteLL, tr->cdta->endsite * sizeof(double)); 

  for(i = 0; i < tr->numberOfTrees; i++)
    {
      int 	
	j;
	 
      double 
	temp, 
	wtemp, 
	sum = 0.0, 
	sum2 = 0.0, 
	sd;           

      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE);
     
      treeEvaluate(tr, 2);
      tr->start = tr->nodep[1];      

      currentLH = tr->likelihood;
      
      if(currentLH > bestLH)	
	printBothOpen("Better tree found %d at %f\n", i, currentLH);	 

      evaluateGenericVector(tr, tr->start);         

      sum = 0.0;
      sum2 = 0.0;

      for (j = 0; j < tr->cdta->endsite; j++)
	{
	  temp  = bestVector[j] - tr->perSiteLL[j];
	  wtemp = tr->cdta->aliaswgt[j] * temp;
	  sum  += wtemp;
	  sum2 += wtemp * temp;
	}

      sd = sqrt( weightSum * (sum2 - sum*sum / weightSum) / (weightSum - 1) );
      /* this is for a 5% p level */
	 
      printBothOpen("Tree: %d Likelihood: %f D(LH): %f SD: %f Significantly Worse: %s (5%s), %s (2%s), %s (1%s)\n", 
		    i, currentLH, currentLH - bestLH, sd, 
		    (sum > 1.95996 * sd) ? "Yes" : " No", "%",
		    (sum > 2.326 * sd) ? "Yes" : " No", "%",
		    (sum > 2.57583 * sd) ? "Yes" : " No", "%");	       
    }

 
  free(bestVector);
  fclose(treeFile);
  exit(0);
}

static void computePerSiteLLs(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int   
    i;
  
  FILE 
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef),
    *tlf = myfopen(perSiteLLsFileName, "wb");

  double
    *unsortedSites = (double*)malloc(sizeof(double) * tr->rdta->sites);

  

  fprintf(tlf, "  %d  %d\n", tr->numberOfTrees, tr->rdta->sites);

  for(i = 0; i < tr->numberOfTrees; i++)
    {      
      int 	
	k, 
	j;
           
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE);
      assert(tr->ntips == tr->mxtips);
      
      if(i == 0) 
	{
	  if(adef->useBinaryModelFile)
	    {
	      readBinaryModel(tr);
	      evaluateGenericInitrav(tr, tr->start);
	      treeEvaluate(tr, 2);
	    }
	  else
	    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);	
	}
      else
	treeEvaluate(tr, 2);     

      tr->start = tr->nodep[1];     
      
      evaluateGenericVector(tr, tr->start);                
      
      printBothOpen("Tree %d: %f\n", i, tr->likelihood);

      fprintf(tlf, "tr%d\t", i + 1);     

      for(j = 0; j < tr->cdta->endsite; j++)
	{
	  for(k = 0; k < tr->rdta->sites; k++)	    
	    if(j == tr->patternPosition[k])		
	      unsortedSites[tr->columnPosition[k] - 1] = tr->perSiteLL[j];		 	          
	}      

      for(j = 0; j < tr->rdta->sites; j++)	  	
	fprintf(tlf, "%f ", unsortedSites[j]);	   	             	                   

      fprintf(tlf, "\n");
    }

  fclose(treeFile);

  free(unsortedSites); 
  fclose(tlf);   
}


static double cumulativeTreeLength(tree *tr, analdef *adef)
{
  double tl = 0.0;

  if(adef->perGeneBranchLengths)
    {
      int 
	accWgt = 0,
	model;
      
      double 
	accLength = 0.0;

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int
	    wgt = 0,
	    i,
	    lower,
	    upper;

	  double 
	    tlm;

	  tlm = treeLength(tr, model);

	  lower = tr->partitionData[model].lower;
	  upper = tr->partitionData[model].upper;

	  for(i = lower; i < upper; i++)
	    wgt += tr->cdta->aliaswgt[i];

	  accLength += ((double)wgt) * tlm;
	  accWgt += wgt;	  
	}

      tl = accLength / ((double)accWgt);

    }
  else
    tl = treeLength(tr, 0);


  return tl;
}

static void computeAllLHs(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int    
    i;

  double
    bestLH = unlikely;
  
  bestlist 
    *bestT;
  
  FILE 
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef),
    *result = myfopen(resultFileName, "wb");

  elw 
    *list;
  
  INFILE = getNumberOfTrees(tr, bootStrapFileName, adef);
 
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);   
 
  list = (elw *)malloc(sizeof(elw) * tr->numberOfTrees); 
   
  for(i = 0; i < tr->numberOfTrees; i++)
    {            
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE);
      resetBranches(tr); 
      
      if(i == 0)
	{
	  testGapped(tr);

	  if(adef->useBinaryModelFile)
	    {
	      readBinaryModel(tr);
	      evaluateGenericInitrav(tr, tr->start);
	      treeEvaluate(tr, 2);
	    }
	  else
	    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);
	  
	  printBothOpen("Model optimization on first Tree: %f\n", tr->likelihood);	  	 
	}
      else
      	{	       
	  evaluateGenericInitrav(tr, tr->start);
           
	  /*
	    treeEvaluateProgressive(tr);
	    treeEvaluateRandom(tr, 2);      
	  */
	  
	  treeEvaluate(tr, 2);
	}            

      list[i].tree = i;
      list[i].lh   = tr->likelihood;

      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE,
		  TRUE, adef, SUMMARIZE_LH, FALSE, FALSE);

      fprintf(result, "%s", tr->tree_string);

      saveBestTree(bestT, tr);

      if(tr->likelihood > bestLH)
	bestLH   = tr->likelihood;
      
      printBothOpen("Tree %d Likelihood %f Tree-Length %f\n", i, tr->likelihood, cumulativeTreeLength(tr, adef));    
    }

  qsort(list, tr->numberOfTrees, sizeof(elw), elwCompareLikelihood);

  printBothOpen("\n");
  for(i = 0; i < tr->numberOfTrees; i++)
    printBothOpen("%d %f\n", list[i].tree, list[i].lh);

  printBothOpen("\n");
 
  /*
    recallBestTree(bestT, 1, tr);
    evaluateGeneric(tr, tr->start);
    printf("Model optimization, %f <-> %f\n", bestLH, tr->likelihood);
    fprintf(infoFile, "Model optimization, %f <-> %f\n", bestLH, tr->likelihood);
    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
    treeEvaluate(tr, 2);
    printf("Model optimization, %f <-> %f\n", bestLH, tr->likelihood);
    fprintf(infoFile, "Model optimization, %f <-> %f\n", bestLH, tr->likelihood);
  */

  printBothOpen("\nAll evaluated trees with branch lengths written to File: %s\n", resultFileName);
  printBothOpen("\nTotal execution time: %f\n", gettime() - masterTime);

  
  fclose(result);
  exit(0);
}




static void computeELW(tree *tr, analdef *adef, char *bootStrapFileName)
{
  FILE 
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  int
    position = 0,   
    bestIndex = -1,  
    i,
    k,
    *originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int)),
    *originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int)),
    *countBest;

  long 
    startSeed;

  double
    best = unlikely,
    **lhs,
    **lhweights,
    sum = 0.0;

  elw
    *bootweights,
    **rankTest;

  initModel(tr, tr->rdta, tr->cdta, adef);      

  if(tr->numberOfTrees < 2)
    {
      printBothOpen("Error, there is only one tree in file %s which you want to use to conduct an ELW test\n", bootStrapFileName);

      exit(-1);
    }  

  bootweights = (elw *)malloc(sizeof(elw) * tr->numberOfTrees);

  rankTest = (elw **)malloc(sizeof(elw *) * adef->multipleRuns);

  for(k = 0; k < adef->multipleRuns; k++)
    rankTest[k] = (elw *)malloc(sizeof(elw) * tr->numberOfTrees);

  lhs = (double **)malloc(sizeof(double *) * tr->numberOfTrees);

  for(k = 0; k < tr->numberOfTrees; k++)
    lhs[k] = (double *)calloc(adef->multipleRuns, sizeof(double));


  lhweights = (double **)malloc(sizeof(double *) * tr->numberOfTrees);

  for(k = 0; k < tr->numberOfTrees; k++)
    lhweights[k] = (double *)calloc(adef->multipleRuns, sizeof(double));

  countBest = (int*)calloc(adef->multipleRuns, sizeof(int));

  /* read in the first tree and optimize ML params on it */  
  
  treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE);   
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);
  rewind(treeFile);

  printBothOpen("Model optimization, first Tree: %f\n", tr->likelihood);

  memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
  memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

  assert(adef->boot > 0);

  /* TODO this is ugly, should be passed as param to computenextreplicate() */

  startSeed = adef->boot;


  /*
     now read the trees one by one, do a couple of BS replicates and re-compute their likelihood
     for every replicate
  */

  /* loop over all trees */

  for(i = 0; i < tr->numberOfTrees; i++)
    {     
      position = 0;

      /* read in new tree */
     
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE); 
      treeEvaluate(tr, 2.0);
      printBothOpen("Original tree %d likelihood %f\n", i, tr->likelihood);

      if(tr->likelihood > best)
	{
	  best      = tr->likelihood;
	  bestIndex = i;
	}
      /* reset branches to default values */

      resetBranches(tr);

      /* reset BS random seed, we want to use the same replicates for every tree */

      adef->rapidBoot = startSeed;

      for(k = 0; k < adef->multipleRuns; k++)
	{
	  /* compute the next BS replicate, i.e., re-sample alignment columns */

	  computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE, TRUE);

	  /* if this is the first replicate for this tree do a slightly more thorough br-len opt */
	  /* we don't re-estimate ML model params (except branches) for every replicate to make things a bit faster */

	  if(k == 0)
	    treeEvaluate(tr, 2.0);
	  else
	    treeEvaluate(tr, 0.5);	  

	  /* store the likelihood of replicate k for tree i */
	  lhs[i][k] = tr->likelihood;

	  rankTest[k][i].lh   = tr->likelihood;
	  rankTest[k][i].tree = i;
	}

      /* restore the original alignment to start BS procedure for the next tree */

      reductionCleanup(tr, originalRateCategories, originalInvariant);
    }

  assert(bestIndex >= 0 && best != unlikely);

  printBothOpen("Best-Scoring tree is tree %d with score %f\n", bestIndex, best);


  /* now loop over all replicates */

  for(k = 0; k < adef->multipleRuns; k++)
    {
      /* find best score for this replicate */

      for(i = 0, best = unlikely; i < tr->numberOfTrees; i++)
	if(lhs[i][k] > best)
	  best = lhs[i][k];

      /* compute exponential weights w.r.t. the best likelihood for replicate k */

      for(i = 0; i < tr->numberOfTrees; i++)
	lhweights[i][k] = exp(lhs[i][k] - best);

      /* sum over all exponential weights */

      for(i = 0, sum = 0.0; i < tr->numberOfTrees; i++)
	sum += lhweights[i][k];

      /* and normalize by the sum */

      for(i = 0; i < tr->numberOfTrees; i++)
	lhweights[i][k] = lhweights[i][k] / sum;

    }

  /* now loop over all trees */

  for(i = 0; i < tr->numberOfTrees; i++)
    {

      /* loop to sum over all replicate weights for tree i  */

      for(k = 0, sum = 0.0; k < adef->multipleRuns; k++)
	sum += lhweights[i][k];

      /* set the weight and the index of the respective tree */

      bootweights[i].weight = sum / ((double)adef->multipleRuns);
      bootweights[i].tree   = i;
    }

  /* now just sort the tree collection by weights */

  qsort(bootweights, tr->numberOfTrees, sizeof(elw), elwCompare);

  printBothOpen("Tree\t Posterior Probability \t Cumulative posterior probability\n");

  /* loop over the sorted array of trees and print out statistics */

  for(i = 0, sum = 0.0; i < tr->numberOfTrees; i++)
    {
      sum += bootweights[i].weight;

      printBothOpen("%d\t\t %f \t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum);
    }


  /*
    if(0)  
    {
    // now compute the super-duper rank test
    
    printBothOpen("\n\nNow also computing the super-duper rank test, though I still don't\n");
    printBothOpen("understand what it actually means. What this thing does is to initially determine\n");
    printBothOpen("the best-scoring ML tree on the original alignment and then the scores of the input\n");
    printBothOpen("trees on the number of specified Bootstrap replicates. Then it sorts the scores of the trees\n");
    printBothOpen("for every bootstrap replicate and determines the rank of the best-scoring tree on every BS\n");
    printBothOpen("replicate. It then prints out how many positions in the sorted lists of thz BS replicates \n");
    printBothOpen("must be included in order for the best scoring tree to appear 95 and 99 times respectively.\n");
    printBothOpen("This gives some intuition about how variable the score order of the trees will be under\n");
    printBothOpen("slight alterations of the data.\n\n");
    
    // sort all BS replicates accodring to likelihood scores
    
    for(i = 0; i < adef->multipleRuns; i++)
    qsort(rankTest[i], tr->numberOfTrees, sizeof(elw), elwCompareLikelihood);
    
    
    // search for our best-scoring tree in every sorted array of likelihood scores
    
    for(i = 0; i < adef->multipleRuns; i++)
    {
    for(k = 0; k < tr->numberOfTrees; k++)
    {
    if(rankTest[i][k].tree == bestIndex)
    countBest[k]++;
    }
    }
    
    for(k = 0; k < tr->numberOfTrees; k++)
    {
    if(k > 0)
    countBest[k] += countBest[k - 1];
    
    printBothOpen("Number of Occurences of best-scoring tree for %d BS replicates up to position %d in sorted list: %d\n",
    adef->multipleRuns, k, countBest[k]);
    
    if(cutOff95 == -1 && countBest[k] <= (int)((double)adef->multipleRuns * 0.95 + 0.5))
    cutOff95 = k;
    
    if(cutOff99 == -1 && countBest[k] <= (int)((double)adef->multipleRuns * 0.99 + 0.5))
    cutOff99 = k;
    }
    
    assert(countBest[k-1] == adef->multipleRuns);
    assert(cutOff95 >= 0 && cutOff99 >= 0);
    
    printBothOpen("\n95%s cutoff reached after including %d out of %d sorted likelihood columns\n", "%", countBest[cutOff95], adef->multipleRuns);
    
    printBothOpen("99%s cutoff reached after including %d out of %d sorted likelihood columns\n\n", "%", countBest[cutOff99], adef->multipleRuns);
    }
  */
    
  printBothOpen("\nTotal execution time: %f\n\n", gettime() - masterTime);
  
  free(originalRateCategories);
  free(originalInvariant);
  fclose(treeFile);  
 
  exit(0);
}



static void computeDistances(tree *tr, analdef *adef)
{
  int i, j, modelCounter;
  double z0[NUM_BRANCHES];
  double result[NUM_BRANCHES];
  double t;
  char distanceFileName[1024];

  FILE
    *out;

  strcpy(distanceFileName,         workdir);
  strcat(distanceFileName,         "RAxML_distances.");
  strcat(distanceFileName,         run_id);

  out = myfopen(distanceFileName, "wb"); 
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);

  printBothOpen("\nLog Likelihood Score after parameter optimization: %f\n\n", tr->likelihood);
  printBothOpen("\nComputing pairwise ML-distances ...\n");

  for(modelCounter = 0; modelCounter < tr->NumberOfModels; modelCounter++)
    z0[modelCounter] = defaultz;

  t = gettime();

  for(i = 1; i <= tr->mxtips; i++)
    for(j = i + 1; j <= tr->mxtips; j++)
      {
	double z, x;

	makenewzGenericDistance(tr, 10, z0, result, i, j);

	if(tr->multiBranch)
	  {
	    int k;

	    for(k = 0, x = 0.0; k < tr->numBranches; k++)
	      {
		assert(tr->partitionContributions[k] != -1.0);
		assert(tr->fracchanges[k] != -1.0);
		z = result[k];
		if (z < zmin)
		  z = zmin;
		x += (-log(z) * tr->fracchanges[k]) * tr->partitionContributions[k];
	      }
	  }
	else
	  {
	    z = result[0];
	    if (z < zmin)
	      z = zmin;
	    x = -log(z) * tr->fracchange;
	  }

	/*printf("%s-%s \t %f\n", tr->nameList[i], tr->nameList[j], x);*/
	fprintf(out, "%s %s \t %f\n", tr->nameList[i], tr->nameList[j], x);
      }

  fclose(out);

  t = gettime() - t;

  printBothOpen("\nTime for pair-wise ML distance computation of %d distances: %f seconds\n",
		 (tr->mxtips * tr->mxtips - tr->mxtips) / 2, t);
  printBothOpen("\nDistances written to file: %s\n", distanceFileName);



  exit(0);
}



static void morphologicalCalibration(tree *tr, analdef *adef)
{
  int 
    replicates = adef->multipleRuns,
    i,     
    *significanceCounter = (int*)malloc(sizeof(int) * tr->cdta->endsite); 

  double 
    *reference  = (double*)malloc(sizeof(double) *  tr->cdta->endsite);

  char    
    integerFileName[1024] = "";

  FILE 
    *integerFile;

  if(replicates == 1)
    {
      printBothOpen("You did not specify the number of random trees to be generated by \"-#\" !\n");
      printBothOpen("Automatically setting it to 100.\n");
      replicates = 100;
    }      

  printBothOpen("Likelihood on Reference tree: %f\n\n", tr->likelihood);  

  evaluateGenericVector(tr, tr->start);

  for(i = 0; i < tr->cdta->endsite; i++)    
    significanceCounter[i] = 0;             

  memcpy(reference, tr->perSiteLL, tr->cdta->endsite * sizeof(double));

  for(i = 0; i < replicates; i++)
    {    
      int k;
      
      printBothOpen("Testing Random Tree [%d]\n", i);
      makeRandomTree(tr, adef);
      evaluateGenericInitrav(tr, tr->start);
      treeEvaluate(tr, 2);
      
      /*
	don't really need modOpt here
	modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
      */
      
      evaluateGenericVector(tr, tr->start);
            
      
      for(k = 0; k < tr->cdta->endsite; k++)	
	if(tr->perSiteLL[k] <= reference[k])
	  significanceCounter[k] = significanceCounter[k] + 1;	        
    }
   
  strcpy(integerFileName,         workdir);
  strcat(integerFileName,         "RAxML_weights.");
  strcat(integerFileName,         run_id);

  integerFile = myfopen(integerFileName, "wb");  

  for(i = 0; i < tr->cdta->endsite; i++)   
    fprintf(integerFile, "%d ", significanceCounter[i]);
    
  fclose(integerFile);
 
  printBothOpen("RAxML calibrated integer weight file written to: %s\n", integerFileName);

  exit(0);
}









static void extractTaxaFromTopology(tree *tr, rawdata *rdta, cruncheddata *cdta)
{
  FILE *f = myfopen(bootStrapFile, "rb");

  char 
    **nameList,
    buffer[nmlngth + 2]; 

  int
    i = 0,
    c,
    taxaSize = 1024,
    taxaCount = 0;
   
  nameList = (char**)malloc(sizeof(char*) * taxaSize);  

  while((c = fgetc(f)) != ';')
    {
      if(c == '(' || c == ',')
	{
	  c = fgetc(f);
	  if(c ==  '(' || c == ',')
	    ungetc(c, f);
	  else
	    {	      
	      i = 0;	      	     
	      
	      do
		{
		  buffer[i++] = c;
		  c = fgetc(f);
		}
	      while(c != ':' && c != ')' && c != ',');
	      buffer[i] = '\0';	    

	      for(i = 0; i < taxaCount; i++)
		{
		  if(strcmp(buffer, nameList[i]) == 0)
		    {
		      printf("A taxon labelled by %s appears twice in the first tree of tree collection %s, exiting ...\n", buffer, bootStrapFile);
		      exit(-1);
		    }
		}	     
	     
	      if(taxaCount == taxaSize)
		{		  
		  taxaSize *= 2;
		  nameList = (char **)realloc(nameList, sizeof(char*) * taxaSize);		 
		}
	      
	      nameList[taxaCount] = (char*)malloc(sizeof(char) * (strlen(buffer) + 1));
	      strcpy(nameList[taxaCount], buffer);
	     
	      taxaCount++;
			    
	      ungetc(c, f);
	    }
	}   
    }
  
  printf("Found a total of %d taxa in first tree of tree collection %s\n", taxaCount, bootStrapFile);
  printf("Expecting all remaining trees in collection to have the same taxon set\n");

  rdta->numsp = taxaCount;

  tr->nameList = (char **)malloc(sizeof(char *) * (taxaCount + 1));  
  for(i = 1; i <= taxaCount; i++)
    tr->nameList[i] = nameList[i - 1];
  
  free(nameList);

  tr->rdta       = rdta;
  tr->cdta       = cdta;

  if (rdta->numsp < 4)
    {    
      printf("TOO FEW SPECIES, tree contains only %d species\n", rdta->numsp);
      assert(0);
    }

  tr->nameHash = initStringHashTable(10 * taxaCount);
  for(i = 1; i <= taxaCount; i++)
    addword(tr->nameList[i], tr->nameHash, i);

  fclose(f);
}


static void myfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, stream);

  assert(bytes_written = nmemb);
}


void writeBinaryModel(tree *tr)
{
  int   
    model; 
  
  FILE 
    *f = myfopen(binaryModelParamsOutputFileName, "w"); 

  /* cdta */   

  myfwrite(tr->cdta->rateCategory, sizeof(int), tr->rdta->sites + 1, f);
  myfwrite(tr->cdta->patrat, sizeof(double), tr->rdta->sites + 1, f);
  myfwrite(tr->cdta->patratStored, sizeof(double), tr->rdta->sites + 1, f);
  
  
  /* pInfo */
   
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	dataType = tr->partitionData[model].dataType;
      
      myfwrite(tr->partitionData[model].gammaRates, sizeof(double), 4, f);
      myfwrite(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
      myfwrite(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
      myfwrite(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);  

      myfwrite(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
      myfwrite(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);  
      myfwrite(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);     
      
      myfwrite(&(tr->partitionData[model].alpha), sizeof(double),  1, f);
      myfwrite(&(tr->partitionData[model].propInvariant), sizeof(double), 1, f);
    }

  printBothOpen("\nModel parameters (binary file format) written to: %s\n", binaryModelParamsOutputFileName);
  
  fclose(f);
}

static void myfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t
    bytes_read;
  
  bytes_read = fread(ptr, size, nmemb, stream);

  assert(bytes_read == nmemb);
}

void readBinaryModel(tree *tr)
{
  int  
    model;

  FILE 
    *f;


  printBothOpen("\nRAxML is reading a binary model file and not optimizing model params\n");

  f = fopen(binaryModelParamsInputFileName, "r");   

  /* cdta */   

  myfread(tr->cdta->rateCategory, sizeof(int),    (size_t)(tr->rdta->sites + 1), f);
  myfread(tr->cdta->patrat,       sizeof(double), (size_t)(tr->rdta->sites + 1), f);
  myfread(tr->cdta->patratStored, sizeof(double), (size_t)(tr->rdta->sites + 1), f);
  
  
  /* pInfo */
   
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	dataType = tr->partitionData[model].dataType;
      
      myfread(tr->partitionData[model].gammaRates, sizeof(double), 4, f);
      
      myfread(tr->partitionData[model].EIGN, sizeof(double), (size_t)(pLengths[dataType].eignLength), f);
      myfread(tr->partitionData[model].EV,   sizeof(double), (size_t)(pLengths[dataType].evLength), f);
      myfread(tr->partitionData[model].EI,   sizeof(double), (size_t)(pLengths[dataType].eiLength), f);  

      myfread(tr->partitionData[model].frequencies, sizeof(double),  (size_t)(pLengths[dataType].frequenciesLength), f);
      myfread(tr->partitionData[model].tipVector,   sizeof(double),  (size_t)(pLengths[dataType].tipVectorLength),   f);  
      myfread(tr->partitionData[model].substRates,  sizeof(double),  (size_t)(pLengths[dataType].substRatesLength),  f);     
      
      myfread(&(tr->partitionData[model].alpha),         sizeof(double), 1, f);
      myfread(&(tr->partitionData[model].propInvariant), sizeof(double), 1, f);
    }

#ifdef _USE_PTHREADS
  /* TODO need to add stuff if we want this to work for CAT models as well */
  masterBarrier(THREAD_RESET_MODEL, tr);
#endif  

  fclose(f);
}
  

void testGapped(tree *tr)
{
  if((!tr->saveMemory) && (!tr->estimatePerSiteAA) && tr->rateHetModel == GAMMA && tr->useGappedImplementation == FALSE)
    {
      int 
	i;
      
      double 
	gappedTime,
	ungappedTime;
      
      printBothOpen("Testing which likelihood implementation to use\n");
      
      tr->useGappedImplementation = FALSE;
      ungappedTime = gettime();
      for(i = 0; i < 8; i++)
	evaluateGenericInitrav(tr, tr->start);
      ungappedTime = gettime() - ungappedTime;
      
      tr->useGappedImplementation = TRUE;
      
#ifdef _USE_PTHREADS
      masterBarrier(THREAD_USE_GAPPED, tr);
#endif
      
      gappedTime = gettime();
      for(i = 0; i < 8; i++)
	evaluateGenericInitrav(tr, tr->start);
      gappedTime = gettime() - gappedTime;
      
      
      printBothOpen("Standard Implementation full tree traversal time: %f\n", ungappedTime);
      
      printBothOpen("Subtree Equality Vectors for gap columns full tree traversal time: %f\n", gappedTime);
      
      if((0.8 * ungappedTime) <= gappedTime) /* purely empirical */
	{
	  tr->useGappedImplementation = FALSE;
	  printBothOpen("... using standard implementation\n\n");
	}
      else
	{
	  tr->useGappedImplementation = TRUE;
	  printBothOpen("... using SEV-based implementation\n\n");
	}
      
#ifdef _USE_PTHREADS
      masterBarrier(THREAD_USE_GAPPED, tr);
#endif
    }
}

static int iterated_bitcount(unsigned int n)
{
    int 
      count=0;    
    
    while(n)
      {
        count += n & 0x1u ;    
        n >>= 1 ;
      }
    
    return count;
}

static char bits_in_16bits [0x1u << 16];

static void compute_bits_in_16bits(void)
{
    unsigned int i;    
    
    for (i = 0; i < (0x1u<<16); i++)
        bits_in_16bits[i] = iterated_bitcount(i);
    
    return ;
}

unsigned int precomputed16_bitcount (unsigned int n)
{
  /* works only for 32-bit int*/
    
    return bits_in_16bits [n         & 0xffffu]
        +  bits_in_16bits [(n >> 16) & 0xffffu] ;
}

/* function to compute the likelihood on quartets */


static double quartetLikelihood(tree *tr, nodeptr p1, nodeptr p2, nodeptr p3, nodeptr p4, nodeptr q1, nodeptr q2)
{
  /* 
     build a quartet tree, where q1 and q2 are the inner nodes and p1, p2, p3, p4
     are the tips of the quartet where the sequence data is located.

     initially set all branch lengths to the default value.
  */

  /* 
     for the tree and node data structure used, please see one of the last chapter's of Joe 
     Felsensteins book. 
  */

  hookupDefault(q1, q2, tr->numBranches);
  
  hookupDefault(q1->next,       p1, tr->numBranches);
  hookupDefault(q1->next->next, p2, tr->numBranches);
  
  hookupDefault(q2->next,       p3, tr->numBranches);
  hookupDefault(q2->next->next, p4, tr->numBranches);
  
  /* now compute the likelihood vectors at the two inner nodes of the tree,
     here the virtual root is located between the two inner nodes q1 and q2.
  */

  newviewGeneric(tr, q1);
  newviewGeneric(tr, q2);
  
  /* call a function that is also used for NNIs that iteratively optimizes all 
     5 branch lengths in the tree.

     Note that 16 is an important tuning parameter, this integer value determines 
     how many times we visit all branches until we give up further optimizing the branch length 
     configuration.
  */

  nniSmooth(tr, q1, 16);

  /* now compute the log likelihood of the tree for the virtual root located between inner nodes q1 and q2 */
  
  /* debugging code 
     {
    double l;
  */
  
  evaluateGeneric(tr, q1->back->next->next);
  
  /* debugging code 
     
     l = tr->likelihood;

     newviewGeneric(tr, q1);
     newviewGeneric(tr, q2);
     evaluateGeneric(tr, q1);
     
   
     assert(ABS(l - tr->likelihood) < 0.00001);
     }
  */

  return (tr->likelihood);
}

static void computeAllThreeQuartets(tree *tr, nodeptr q1, nodeptr q2, int t1, int t2, int t3, int t4, FILE *f)
{
  /* set the tip nodes to different sequences 
     with the tip indices t1, t2, t3, t4 */
	       
  nodeptr 
    p1 = tr->nodep[t1],
    p2 = tr->nodep[t2],
    p3 = tr->nodep[t3], 
    p4 = tr->nodep[t4];
  
  double 
    l;
  
  /* first quartet */	    
  
  /* compute the likelihood of tree ((p1, p2), (p3, p4)) */
  
  l = quartetLikelihood(tr, p1, p2, p3, p4, q1, q2);
  
  fprintf(f, "%d %d | %d %d: %f\n", p1->number, p2->number, p3->number, p4->number, l);
  
  /* second quartet */	    
  
  /* compute the likelihood of tree ((p1, p3), (p2, p4)) */
  
  l = quartetLikelihood(tr, p1, p3, p2, p4, q1, q2);
  
  fprintf(f, "%d %d | %d %d: %f\n", p1->number, p3->number, p2->number, p4->number, l);
  
  /* third quartet */	    
  
  /* compute the likelihood of tree ((p1, p4), (p2, p3)) */
  
  l = quartetLikelihood(tr, p1, p4, p2, p3, q1, q2);
  
  fprintf(f, "%d %d | %d %d: %f\n", p1->number, p4->number, p2->number, p3->number, l);	    	   
}

static void computeQuartets(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  /* some indices for generating quartets in an arbitrary way */

  int
    i,
    t1, 
    t2, 
    t3, 
    t4;

  double
    fraction,
    t;

  unsigned long int
    randomQuartets = (unsigned long int)(adef->multipleRuns),
    quartetCounter = 0,
    numberOfQuartets = ((unsigned long int)tr->mxtips * ((unsigned long int)tr->mxtips - 1) * ((unsigned long int)tr->mxtips - 2) * ((unsigned long int)tr->mxtips - 3)) / 24;

  /* use two inner nodes for building quartet trees */

  nodeptr 	
    q1 = tr->nodep[tr->mxtips + 1],
    q2 = tr->nodep[tr->mxtips + 2];


  char 
    quartetFileName[1024];

  FILE 
    *f;

  strcpy(quartetFileName,         workdir);
  strcat(quartetFileName,         "RAxML_quartets.");
  strcat(quartetFileName,         run_id);
  
  f = myfopen(quartetFileName, "w");

  /* initialize model parameters */

  initModel(tr, rdta, cdta, adef);
      
  /* get a starting tree: either reads in a tree or computes a randomized stepwise addition parsimony tree */

  getStartingTree(tr, adef);
  
  /* optimize model parameters on that comprehensive tree that can subsequently be used for qyartet building */

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);

  printBothOpen("Time for parsing input tree or building parsimony tree and optimizing model parameters: %f\n\n", gettime() - masterTime); 

  if(randomQuartets > numberOfQuartets)
    randomQuartets = 1;

  if(randomQuartets == 1)
    printBothOpen("There are %u quartet sets for which RAxML will evaluate all %u quartet trees\n", numberOfQuartets, numberOfQuartets * 3);
  else
    {
      /* cast from unsigned long int to double may be dangeruous for very large integer values */

      fraction = (double)randomQuartets / (double)numberOfQuartets;

      printBothOpen("There are %u quartet sets for which RAxML will randomly sub-sambple %u sets (%f\%), i.e., compute %u quartet trees\n", numberOfQuartets, randomQuartets, 100 * fraction, randomQuartets * 3);
    }

  fprintf(f, "Taxon names and indices:\n\n");

  for(i = 1; i <= tr->mxtips; i++)
    {
      fprintf(f, "%s %d\n", tr->nameList[i], i);
      assert(tr->nodep[i]->number == i);
    }

  fprintf(f, "\n\n");

  t = gettime();
  
  /* do a loop to generate some quartets to test.
     note that tip nodes/sequences in RAxML are indexed from 1,...,n
     and not from 0,...,n-1 as one might expect 
     
     tr->mxtips is the maximum number of tips in the alignment/tree
  */

  if(randomQuartets == 1)
    {
      for(t1 = 1; t1 <= tr->mxtips; t1++)
	for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
	  for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
	    for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
	      {
		computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f);
		quartetCounter++;
	      }
      
      assert(quartetCounter == numberOfQuartets);
    }
  else
    {
      for(t1 = 1; t1 <= tr->mxtips; t1++)
	for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
	  for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
	    for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
	      {
		double
		  r = randum(&adef->parsimonySeed);

		if(r < fraction)
		  {
		    computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f);
		    quartetCounter++;
		  }

		if(quartetCounter == randomQuartets)
		  goto DONE;
	      }
      
    DONE:
      assert(quartetCounter == randomQuartets);
    }

  t = gettime() - t;

  printBothOpen("\nPure quartet computation time: %f secs\n", t);
  
  printBothOpen("\nAll quartets and corresponding likelihoods written to file %s\n", quartetFileName);

  fclose(f);
}

static void thoroughTreeOptimization(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  char 
    bestTreeFileName[1024]; 

  FILE 
    *f;
  
  initModel(tr, rdta, cdta, adef);
      
  getStartingTree(tr, adef);  

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);

  Thorough = 1;
  tr->doCutoff = FALSE;  
	 
  printBothOpen("\nStart likelihood: %f\n\n", tr->likelihood);

  treeOptimizeThorough(tr, 1, 10);
  evaluateGenericInitrav(tr, tr->start);
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);

  printBothOpen("End likelihood: %f\n\n", tr->likelihood);

  printModelParams(tr, adef);    
  
  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_bestTree.");
  strcat(bestTreeFileName,         run_id);

  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE);
  f = myfopen(bestTreeFileName, "wb");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);

  printBothOpen("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
}

int main (int argc, char *argv[])
{
  rawdata      *rdta;
  cruncheddata *cdta;
  tree         *tr;
  analdef      *adef;
  int
    i,
    countGTR = 0,
    countOtherModel = 0;

#if (defined(_USE_PTHREADS) && !defined(_PORTABLE_PTHREADS))  
  pinToCore(0);
#endif 
 
#ifdef _WAYNE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  printf("\nThis is RAxML MPI Process Number: %d\n", processID);
#else
  processID = 0;
#endif

  masterTime = gettime();

  globalArgc = argc;
  globalArgv = (char **)malloc(sizeof(char *) * argc);
  for(i = 0; i < argc; i++)
    globalArgv[i] = argv[i];



#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))

  /* 
     David Defour's command  
     _mm_setcsr( _mm_getcsr() | (_MM_FLUSH_ZERO_ON | MM_DAZ_ON));  
  */

   _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);

#endif 

  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  cdta = (cruncheddata *)malloc(sizeof(cruncheddata));
  tr   = (tree *)malloc(sizeof(tree));

  /* initialize lookup table for fast bit counter */

  compute_bits_in_16bits();

  initAdef(adef);
  get_args(argc,argv, adef, tr); 
  
  if(adef->readTaxaOnly)  
    extractTaxaFromTopology(tr, rdta, cdta);   
 
  getinput(adef, rdta, cdta, tr);

  checkOutgroups(tr, adef);
  makeFileNames();

#ifdef _WAYNE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif

  if(adef->useInvariant && adef->likelihoodEpsilon > 0.001)
    {
      printBothOpen("\nYou are using a proportion of Invariable sites estimate, although I don't\n");
      printBothOpen("like it. The likelihood epsilon \"-f e\" will be automatically lowered to 0.001\n");
      printBothOpen("to avoid unfavorable effects caused by simultaneous optimization of alpha and P-Invar\n");

      adef->likelihoodEpsilon = 0.001;
    }


  /*
     switch back to model without secondary structure for all this
     checking stuff
  */

  if(adef->useSecondaryStructure)
    {
      tr->dataVector    = tr->initialDataVector;
      tr->partitionData = tr->initialPartitionData;
      tr->NumberOfModels--;
    }

  if(adef->useExcludeFile)
    {
      handleExcludeFile(tr, adef, rdta);
      exit(0);
    }

 
  if(!adef->readTaxaOnly && adef->mode != FAST_SEARCH && adef->mode != SH_LIKE_SUPPORTS)
    checkSequences(tr, rdta, adef);
  

  if(adef->mode == SPLIT_MULTI_GENE)
    {
      splitMultiGene(tr, rdta);
      exit(0);
    }

  if(adef->mode == CHECK_ALIGNMENT)
    {
      printf("Alignment format can be read by RAxML \n");
      exit(0);
    }

  /*
     switch back to model with secondary structure for all this
     checking stuff
  */

  if(adef->useSecondaryStructure && !adef->readTaxaOnly)
    {
      tr->dataVector    = tr->extendedDataVector;
      tr->partitionData = tr->extendedPartitionData;
      tr->NumberOfModels++;
      /* might as well free the initial structures here */

    }
  
  if(!adef->readTaxaOnly)
    {
      int 
	countNonSev = 0;

      makeweights(adef, rdta, cdta, tr);
      makevalues(rdta, cdta, tr, adef);      

      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  if(!(tr->partitionData[i].dataType == AA_DATA || tr->partitionData[i].dataType == DNA_DATA))
	    countNonSev++;

	  if(tr->partitionData[i].dataType == AA_DATA)
	    {
	      if(tr->partitionData[i].protModels == GTR || tr->partitionData[i].protModels == GTR_UNLINKED)
		countGTR++;
	      else
		countOtherModel++;
	    }
	}

      if(tr->saveMemory && countNonSev > 0)
	{
	  printf("\nError, you want to use the SEV-based memory saving technique for large gappy datasets with missing data.\n");
	  printf("However, this is only implelemented for DNA and protein data partitions, one of your partitions is neither DNA\n");
	  printf("nor protein data ... exiting to prevent bad things from happening ;-) \n\n");

	  errorExit(-1);
	}


      if(countGTR > 0 && countOtherModel > 0)
	{
	  printf("Error, it is only allowed to conduct partitioned AA analyses\n");
	  printf("with a GTR model of AA substitution, if not all AA partitions are assigned\n");
	  printf("the GTR or GTR_UNLINKED model.\n\n");
	  
	  printf("The following partitions do not use GTR:\n");
	  
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      if(tr->partitionData[i].dataType == AA_DATA && (tr->partitionData[i].protModels != GTR || tr->partitionData[i].protModels != GTR_UNLINKED))
		printf("Partition %s\n", tr->partitionData[i].partitionName);
	    }
	  printf("exiting ...\n");
	  errorExit(-1);
	}

      if(countGTR > 0 && tr->NumberOfModels > 1)
	{
	  FILE *info = myfopen(infoFileName, "ab");

	  printBoth(info, "You are using the GTR model of AA substitution!\n");
	  printBoth(info, "GTR parameters for AA substiution will automatically be estimated\n");
	  printBoth(info, "either jointly (GTR params will be linked) or independently (when using GTR_UNLINKED) across all partitions.\n");
	  printBoth(info, "WARNING: you may be over-parametrizing the model!\n\n\n");

	  fclose(info);
	}
    }

  if(adef->mode == CLASSIFY_ML || adef->mode == CLASSIFY_MP)              
    tr->innerNodes = (size_t)(countTaxaInTopology() - 1);   
  else
    tr->innerNodes = tr->mxtips;

  
  setRateHetAndDataIncrement(tr, adef);

#ifdef _USE_PTHREADS
  startPthreads(tr);
  masterBarrier(THREAD_INIT_PARTITION, tr);
  if(!adef->readTaxaOnly)  
    masterBarrier(THREAD_ALLOC_LIKELIHOOD, tr);
#else
  if(!adef->readTaxaOnly)  
    allocNodex(tr);    
#endif

  makeMissingData(tr);

  printModelAndProgramInfo(tr, adef, argc, argv);

  switch(adef->mode)
    {  
    case CLASSIFY_MP:
      getStartingTree(tr, adef);
      assert(0);
      break;
    case CLASSIFY_ML:
      initModel(tr, rdta, cdta, adef);
      getStartingTree(tr, adef);
      exit(0);
      break;
    case GENERATE_BS:
      generateBS(tr, adef);
      exit(0);
      break;
    case COMPUTE_ELW:
      computeELW(tr, adef, bootStrapFile);
      exit(0);
      break;
    case COMPUTE_LHS:
      initModel(tr, rdta, cdta, adef);
      computeAllLHs(tr, adef, bootStrapFile);
      exit(0);
      break;
    case COMPUTE_BIPARTITION_CORRELATION:
      compareBips(tr, bootStrapFile, adef);
      exit(0);
      break;
    case COMPUTE_RF_DISTANCE:
      computeRF(tr, bootStrapFile, adef);
      exit(0);
      break;
    case BOOTSTOP_ONLY:
      computeBootStopOnly(tr, bootStrapFile, adef);
      exit(0);
      break;
    case CONSENSUS_ONLY:      
      if(adef->leaveDropMode)
	computeRogueTaxa(tr, bootStrapFile, adef);
      else
	computeConsensusOnly(tr, bootStrapFile, adef);
      exit(0);
      break;
    case DISTANCE_MODE:
      initModel(tr, rdta, cdta, adef);
      getStartingTree(tr, adef);
      computeDistances(tr, adef);
      break;
    case  PARSIMONY_ADDITION:
      initModel(tr, rdta, cdta, adef);
      getStartingTree(tr, adef);
      printStartingTree(tr, adef, TRUE);
      break;
    case PER_SITE_LL:
      initModel(tr, rdta, cdta, adef);
      computePerSiteLLs(tr, adef, bootStrapFile);
      break;
    case TREE_EVALUATION:
      initModel(tr, rdta, cdta, adef);
      
      getStartingTree(tr, adef);      
      
      if(adef->likelihoodTest)
	computeLHTest(tr, adef, bootStrapFile);
      else
	{	 	 
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);	  
	  writeBinaryModel(tr);
	  printLog(tr, adef, TRUE);
	  printResult(tr, adef, TRUE);
	}
  
      break;
    case ANCESTRAL_STATES:
      initModel(tr, rdta, cdta, adef);
      
      getStartingTree(tr, adef);
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);
       
      evaluateGenericInitrav(tr, tr->start);                                       	                  
      
      computeAncestralStates(tr, tr->likelihood, adef);
      break;
    case  QUARTET_CALCULATION:                                             	                        
      computeQuartets(tr, adef, rdta, cdta);
      break;
    case THOROUGH_OPTIMIZATION:
      thoroughTreeOptimization(tr, adef, rdta, cdta);
      break;
    case CALC_BIPARTITIONS:      
      calcBipartitions(tr, adef, tree_file, bootStrapFile);
      break;
    case BIG_RAPID_MODE:
      if(adef->boot)
	doBootstrap(tr, adef, rdta, cdta);
      else
	{
	  if(adef->rapidBoot)
	    {
	      initModel(tr, rdta, cdta, adef);
	      doAllInOne(tr, adef);
	    }
	  else	    	    
	    doInference(tr, adef, rdta, cdta);	     	
	}
      break;
    case MORPH_CALIBRATOR:
      initModel(tr, rdta, cdta, adef);
      getStartingTree(tr, adef);
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);
      morphologicalCalibration(tr, adef);
      break;    
    case MESH_TREE_SEARCH:
      initModel(tr, rdta, cdta, adef); 
      getStartingTree(tr, adef); 
      meshTreeSearch(tr, adef, adef->meshSearch);
      /* TODO */
      break;
    case FAST_SEARCH:
      fastSearch(tr, adef, rdta, cdta);
      exit(0);
    case SH_LIKE_SUPPORTS:
      shSupports(tr, adef, rdta, cdta);
      break;    
    case EPA_SITE_SPECIFIC_BIAS:
      initModel(tr, rdta, cdta, adef);
      getStartingTree(tr, adef);
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);
      computePlacementBias(tr, adef);
      break;
    case OPTIMIZE_BR_LEN_SCALER:
      initModel(tr, rdta, cdta, adef);
      
      getStartingTree(tr, adef);      
            	 	 
      modOpt(tr, adef, FALSE, adef->likelihoodEpsilon, FALSE);	  
      
      printBothOpen("Likelihood: %f\n", tr->likelihood);

      break;
    default:
      assert(0);
    }

  finalizeInfoFile(tr, adef);

#ifdef _WAYNE_MPI
  MPI_Finalize();
#endif

  return 0;
}


