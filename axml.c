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
#define WIN32_LEAN_AND_MEAN // skips unwanted headers like socket etc.
#include <windows.h>
#include <direct.h>
#endif

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#endif

#include <time.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdarg.h>
#include <limits.h>
#include <inttypes.h>
#include <getopt.h>
//#include <stdbool.h>


#if (defined(_WAYNE_MPI) || defined (_QUARTET_MPI))
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

static void checkStdoutFlush(void)
{
  /* If stdout is redirected, other processes monitoring RAxML's output
     (e.g., via tail, or a pipe) do not receive any standard output until
     stdio gets around to flushing the file, which may be a long time.
     To provide more continuous feeding of RAxML output to these processes,
     we force a flush of the stdout stream once per second.
     (Dave Swofford 16july2016)
  */
  
  static clock_t 
    lastFlush;
  
  clock_t
    now = clock();
  
  if(now - lastFlush > CLOCKS_PER_SEC)
    {
      fflush(stdout);
      lastFlush = now;
    }
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
  checkStdoutFlush();
}

void printBothOpen(const char* format, ... )
{
#ifdef _QUARTET_MPI
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
      checkStdoutFlush();
      
      fclose(f);
    }     
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
      checkStdoutFlush();
      
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

unsigned char getUndetermined(int dataType)
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
  //pLength.optimizeBaseFrequencies = FALSE;

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
#ifdef WIN32 // WINDOWS build
	FILETIME tm;
	ULONGLONG t;
#if defined(NTDDI_WIN8) && NTDDI_VERSION >= NTDDI_WIN8 // >= WIN8
	GetSystemTimePreciseAsFileTime( &tm );
#else // < WIN8
	GetSystemTimeAsFileTime( &tm );
#endif
	t = ((ULONGLONG)tm.dwHighDateTime << 32) | (ULONGLONG)tm.dwLowDateTime;
	return (double)t / 10000000.0;
#else // Unixoid build
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}



double randum (int64_t  *seed)
{
  int64_t  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
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

#ifdef _BASTIEN
  for(i = 0; i < numBranches; i++)
    p->secondDerivativeValid[i] = q->secondDerivativeValid[i] = FALSE;
#endif
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;

  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;

#ifdef _BASTIEN
  for(i = 0; i < numBranches; i++)
    p->secondDerivativeValid[i] = q->secondDerivativeValid[i] = FALSE;
#endif
}


/***********************reading and initializing input ******************/

static void rax_getline_insptr_valid(char **lineptr, size_t *n, size_t ins_ptr ) 
{
  const size_t 
    n_inc = 1024;

  if(ins_ptr >= *n) 
    {
      assert( *n <= (SSIZE_MAX - n_inc));

      *n += n_inc;
          
      *lineptr = (char*)rax_realloc((void*)(*lineptr), *n * sizeof(char), FALSE);
    
      assert(*lineptr != 0);
  }
}

ssize_t rax_getline(char **lineptr, size_t *n, FILE *h) 
{
  size_t 
    ins_ptr = 0;
  
  /* this implementation does not conform to the standard regarding error checking (i.e., asserts on errors ) */
  
  assert(h != (FILE*)NULL);

  if(*lineptr == (char *)NULL)     
    *n = 0;    

  while(1) 
    {
      int 
	c = fgetc(h);

      /* handle EOF: if no character has been read on the current line throw an error. 
	 Otherwise treat as end-of-line. Don't know if this is correct, 
	 as I don't have the POSIX standard and the linux manpage is unclear. */
      
    if(c == EOF) 
      {
	if(ins_ptr == 0) 	  
	  return -1;	  
	else 
	  c = '\n';
	  //break;	  
      }

    if(c == '\r') 
      {
	//this is the original GNU implementation
	/* windows line-end: must be followed by a '\n'. Don't tolerate anything else. */
	//c = fgetc(h);
	//assert(c == '\n');

	//fixed to essentialy replace windows line endings by '\n'
	c = '\n';
      }

    /* insert character (including '\n') into buffer */
    rax_getline_insptr_valid(lineptr, n, ins_ptr);
    (*lineptr)[ins_ptr] = c;
    ++ins_ptr;

    if(c == '\n')      
      break;    
  }

  /* null-terminate */
  rax_getline_insptr_valid( lineptr, n, ins_ptr );
  (*lineptr)[ins_ptr] = 0;

  return ((ssize_t)ins_ptr);
}


static void getnums (rawdata *rdta, analdef *adef)
{
  if(fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2)
    {      
      char 
	*line = NULL;

      size_t 
	len = 0;

      ssize_t 
	read;     

      int
	sequenceLength = 0,       
	sequences = 0,
	taxa = 0,
	sites =0;
      
      if(processID == 0)
	{
	  printf("\nRAxML can't, parse the alignment file as phylip file \n");
	  printf("it will now try to parse it as FASTA file\n\n");
	}

      while((read = rax_getline(&line, &len, INFILE)) != -1) 
	{
	  ssize_t
	    i = 0;
	  	  	  
	  while((i < read - 1) && (line[i] == ' ' || line[i] == '\t'))
	    i++;
	
	  if(line[i] == '>')
	    {	      
	      if(taxa == 1)   		
		sequenceLength = sites;
	       
	      if(taxa > 0)
		{
		  if(sites == 0 && processID == 0)
		    {
		      printf("Fasta parsing error, RAxML was expecting sequence data before: %s\n", line);
		      errorExit(-1);
		    }
		  assert(sites > 0);
		  sequences++;		 
		}
	      
	      if(taxa > 0)
		{
		  if(sequenceLength != sites && processID == 0)
		    {
		      printf("Fasta parsing error, RAxML expects an alignment.\n");
		      printf("the sequence before taxon %s: seems to have a different length\n", line);
		      errorExit(-1);
		    }
		  assert(sequenceLength == sites);	     
		}
	      
	      taxa++;
	     
	      sites = 0;
	    }
	  else
	    {	     
	      while(i < read - 1)
		{
		  if(!(line[i] == ' ' || line[i] == '\t'))
		    {		    
		      sites++;
		    }		 
		  i++;
		}
	      //printf("sites %d %s\n", sites, line);
	    }	  
	}

      if(sites > 0)
	sequences++;
      if(taxa != sequences && processID == 0)
	{
	  printf("Fasta parsing error, the number of taxa %d and sequences %d are not equal!\n", taxa, sequences);
	  errorExit(-1);
	}
      assert(taxa == sequences);
      
      if(sequenceLength != sites && processID == 0)
	{
	  printf("Fasta parsing error, RAxML expects an alignment.\n");
	  printf("the last sequence in the alignment seems to have a different length\n");
	  errorExit(-1);
	}
      
      assert(sites == sequenceLength);

      if(line)
	rax_free(line);

      rewind(INFILE);

      adef->alignmentFileType = FASTA;

      rdta->numsp = taxa;
      rdta->sites = sites;
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

  rdta->y = (unsigned char **) rax_malloc((rdta->numsp + 1) * sizeof(unsigned char *));
  assert(rdta->y);   

  y0 = (unsigned char *) rax_malloc(((size_t)(rdta->numsp + 1)) * size * sizeof(unsigned char));
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
    tips,
    inter; 
  
  
  
  tr->storedBrLens = (double*)NULL;

  if(!adef->readTaxaOnly)
    {
      tr->bigCutoff = FALSE;

      tr->patternPosition = (int*)NULL;
      tr->columnPosition = (int*)NULL;

      tr->maxCategories = MAX(4, adef->categories);

      tr->partitionContributions = (double *)rax_malloc(sizeof(double) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->partitionContributions[i] = -1.0;

      tr->perPartitionLH = (double *)rax_malloc(sizeof(double) * tr->NumberOfModels);
      tr->storedPerPartitionLH = (double *)rax_malloc(sizeof(double) * tr->NumberOfModels);

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
      tr->yVector      = (unsigned char **)  rax_malloc((tr->mxtips + 1) * sizeof(unsigned char *));     

      tr->likelihoods  = (double *)rax_malloc(adef->multipleRuns * sizeof(double));
    }

  tr->numberOfTrees = -1;

  tr->treeStringLength = 
    2 * (size_t)tr->mxtips + //parentheses
    2 * (size_t)tr->mxtips * 64 + //branche lengths with : and . and branch labels and
    (size_t)tr->mxtips + //commas
    1 + //closing semicolon 
    (size_t)tr->mxtips * nmlngth; //taxon names

  //tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;

  //printf("tips %d Tree String Length %d old length %d\n", tr->mxtips, tr->treeStringLength,tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2 );

  tr->tree_string  = (char*)rax_calloc(tr->treeStringLength, sizeof(char)); 

  if(!adef->readTaxaOnly)
    {
           
      tr->td[0].count = 0;
      tr->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * tr->mxtips);	

     
      tr->constraintVector = (int *)rax_malloc((2 * tr->mxtips) * sizeof(int));

      tr->nameList = (char **)rax_malloc(sizeof(char *) * (tips + 1));
    }

  if (!(p0 = (nodeptr) rax_malloc((tips + 3*inter) * sizeof(node))))
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (!(tr->nodep = (nodeptr *) rax_malloc((2*tr->mxtips) * sizeof(nodeptr))))
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

	 

	  


	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;

  

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

  //printf("checking taxon name\n");

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

static void printParsingErrorContext(FILE *f)
{
  const int64_t
    contextWidth = 20;

  int64_t
    i,
    currentPos = ftell(f),
    contextPos = MAX(currentPos - contextWidth, 0);
  
  fseek(f, MAX(currentPos - contextWidth, 0), SEEK_SET);
  
  printf("Printing error context:\n\n");
  
  for(i = contextPos; i < currentPos + contextWidth; i++)
    {
      int 
	ch = getc(f);
      if(ch != EOF)
	printf("%c", ch);
      else
	break;
    }
  
  printf("\n\n");
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
  uint64_t
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
      for(i = 1; i <= tr->mxtips; i++)
	{
	  if(firstpass)
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
			  printf("Taxon Name too long at taxon %d, adapt constant nmlngth in\n", i);
			  printf("axml.h, current setting %d\n", nmlngth);
			}
		      errorExit(-1);
		    }
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');	      

	      buffer[my_i] = '\0';
	      len = strlen(buffer) + 1;
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)rax_malloc(sizeof(char) * len);
	      strcpy(tr->nameList[i], buffer);	      

	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r')
		ch = getc(INFILE);
	      	      
	      ungetc(ch, INFILE);
	    }
	  
	  j = basesread;

	  while((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r'))
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

		      printParsingErrorContext(INFILE);		      
		      			  
		      return FALSE;
		    }
		}
	    }

	  if (ch == EOF)
	    {
	      printf("ERROR: End-of-file at site %d of sequence %d\n", j + 1, i);
	      
	      printParsingErrorContext(INFILE);

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

		     printParsingErrorContext(INFILE);

		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') 
	    ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= rdta->sites);
    }

  if(ch != EOF)
    {
      int
	garbageCount = 0;  

      do
	{
	  if(!whitechar(ch))
	    garbageCount++;

	  if(garbageCount > 0)
	    {
	      if(garbageCount == 1)
		printf("\nOups, there is garbage at the end of your file:\n\n");
	      if(garbageCount < 1000)
		printf("%c", ch);	     
	      if(garbageCount == 1000)
		printf("\n .... and so on\n");
	    }
	}
      while((ch = getc(INFILE)) != EOF);
      
      if(garbageCount > 0)
	{
	  printf("\n\nRAxML correctly finished parsing your PHYLIP file,\n");
	  printf("but there seems to be garbage at the end of the file, i.e., non-whitespace characters!\n");
	  printf("RAxML will exit now\n\n");
	  errorExit(-1);
	}
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

static void parseFasta(analdef *adef, rawdata *rdta, tree *tr)
{
  int 
    index,
    meaning,
    meaningAA[256], 
    meaningDNA[256], 
    meaningBINARY[256],
    meaningGeneric32[256],
    meaningGeneric64[256];
    
  char 
    buffer[nmlngth + 2];
  
  unsigned char
    genericChars32[32] = {'0', '1', '2', '3', '4', '5', '6', '7', 
			  '8', '9', 'A', 'B', 'C', 'D', 'E', 'F',
			  'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N',
			  'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V'};  
  uint64_t
    total = 0,
    gaps  = 0;

  for(index = 0; index < 256; index++)
    {      
      meaningAA[index]          = -1;
      meaningDNA[index]         = -1;
      meaningBINARY[index]      = -1;
      meaningGeneric32[index]   = -1;
      meaningGeneric64[index]   = -1;
    }

  /* generic 32 data */

  for(index = 0; index < 32; index++)
    meaningGeneric32[genericChars32[index]] = index;
  
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

  {
    char 
      *line = NULL;

    size_t 
      len = 0;
    
    ssize_t 
      read;     
    
    int
      sequenceLength = 0,       
      sequences = 0,
      taxa = 0,
      sites = 0;
    
         
    while((read = rax_getline(&line, &len, INFILE)) != -1) 
	{
	  ssize_t
	    i = 0;
	  	  	  	  
	  while((i < read - 1) && (line[i] == ' ' || line[i] == '\t'))
	    i++;
	
	  if(line[i] == '>')
	    {
	      int
		nameCount = 0,
		nameLength;

	      
	      
	      if(taxa == 1)   		
		sequenceLength = sites;
	       
	      if(taxa > 0)
		{
		  assert(sites > 0);
		  sequences++;		 
		}
	      
	      if(taxa > 0)
		assert(sequenceLength == sites);	     
	      
	      taxa++;
	     
	      i++;
	      
	      while((i < read - 1) && (line[i] == ' ' || line[i] == '\t'))
		i++;

	      while((i < read - 1) && !(line[i] == ' ' || line[i] == '\t'))
		{		  
		  buffer[nameCount] = line[i];
		  nameCount++;
		  i++;
		}

	      if(nameCount >= nmlngth)
		{
		  if(processID == 0)
		    {
		      printf("Taxon Name too long at taxon %d, adapt constant nmlngth in\n", taxa);
		      printf("axml.h, current setting %d\n", nmlngth);
		    }
		  errorExit(-1);
		}

	      buffer[nameCount] = '\0';	      
	      nameLength = strlen(buffer) + 1;
	      checkTaxonName(buffer, nameLength);
	      tr->nameList[taxa] = (char *)rax_malloc(sizeof(char) * nameLength);
	      strcpy(tr->nameList[taxa], buffer);

	      sites = 0;
	    }
	  else
	    {	     
	      while(i < read - 1)
		{
		  if(!(line[i] == ' ' || line[i] == '\t'))
		    {	
		      int 
			ch = line[i];
		      
		      uppercase(&ch);

		      assert(tr->dataVector[sites + 1] != -1);
		      
		      switch(tr->dataVector[sites + 1])
			{
			case BINARY_DATA:
			  meaning = meaningBINARY[ch];
			  break;
			case DNA_DATA:
			case SECONDARY_DATA:
			case SECONDARY_DATA_6:
			case SECONDARY_DATA_7:			 
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
			rdta->y[taxa][sites + 1] = ch;		 			
		      else
			{
			  if(processID == 0)
			    {
			      printf("ERROR: Bad base (%c) at site %d of sequence %d\n",
				     ch, sites + 1, taxa);
			    }
			  errorExit(-1);
			}		    	       	   

		      sites++;
		    }
		  i++;
		}
	    }	  
	}

      if(sites > 0)
	sequences++;

      /* the assertions below should never fail, the have already been checked in getNums */

      assert(taxa == sequences);    
      assert(sites == sequenceLength);

      if(line)
	rax_free(line);
  }


  {
    int 
      i,
      j;

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
  }
    
  adef->gapyness = (double)gaps / (double)total;
    
  return;
}



static void inputweights (rawdata *rdta)
{
  int i, w, fres;
  FILE *weightFile;
  int *wv = (int *)rax_malloc(sizeof(int) *  rdta->sites);

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
  rax_free(wv);
}



static void getinput(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  int i;

  if(!adef->readTaxaOnly)
    {
      INFILE = myfopen(seq_file, "rb");
  
      getnums(rdta, adef);
    }

  tr->mxtips            = rdta->numsp;
  
  if(!adef->readTaxaOnly)
    {
      rdta->wgt             = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int)); 
      cdta->alias           = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
      cdta->aliaswgt        = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
      cdta->rateCategory    = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
      tr->model             = (int *)    rax_calloc((rdta->sites + 1), sizeof(int));
      tr->initialDataVector  = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));
      tr->extendedDataVector = (int *)    rax_malloc((rdta->sites + 1) * sizeof(int));     
      cdta->patrat          = (double *) rax_malloc((rdta->sites + 1) * sizeof(double));
      cdta->patratStored    = (double *) rax_malloc((rdta->sites + 1) * sizeof(double));      
     


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
  else
    tr->NumberOfModels = 0;

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
	  
	  tr->initialPartitionData  = (pInfo*)rax_malloc(sizeof(pInfo));
	  tr->initialPartitionData[0].partitionName = (char*)rax_malloc(128 * sizeof(char));
	  strcpy(tr->initialPartitionData[0].partitionName, "No Name Provided");
	  
	  tr->initialPartitionData[0].protModels = adef->proteinMatrix;
	  if(adef->protEmpiricalFreqs)
	    tr->initialPartitionData[0].usePredefinedProtFreqs  = FALSE;
	  else
	    tr->initialPartitionData[0].usePredefinedProtFreqs  = TRUE;
	  
	  if(adef->optimizeBaseFrequencies)
	    {	     
	      tr->initialPartitionData[0].optimizeBaseFrequencies  = TRUE;
	      tr->initialPartitionData[0].usePredefinedProtFreqs  = FALSE;
	    }
	  else
	    tr->initialPartitionData[0].optimizeBaseFrequencies  = FALSE;

	  if(adef->ascertainmentBias)
	    tr->initialPartitionData[0].ascBias = TRUE;
	  else
	    tr->initialPartitionData[0].ascBias = FALSE;


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
	      tr->initialPartitionData[0].optimizeBaseFrequencies = FALSE;
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
	  
	  tr->extendedPartitionData =(pInfo*)rax_malloc(sizeof(pInfo) * tr->NumberOfModels);
	  
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      tr->extendedPartitionData[i].partitionName = (char*)rax_malloc((strlen(tr->initialPartitionData[i].partitionName) + 1) * sizeof(char));
	      strcpy(tr->extendedPartitionData[i].partitionName, tr->initialPartitionData[i].partitionName);
	      strcpy(tr->extendedPartitionData[i].proteinSubstitutionFileName, tr->initialPartitionData[i].proteinSubstitutionFileName);
	      strcpy(tr->extendedPartitionData[i].ascFileName, tr->initialPartitionData[i].ascFileName);
	      tr->extendedPartitionData[i].dataType   = tr->initialPartitionData[i].dataType;	      
	      tr->extendedPartitionData[i].protModels = tr->initialPartitionData[i].protModels;
	      tr->extendedPartitionData[i].usePredefinedProtFreqs  = tr->initialPartitionData[i].usePredefinedProtFreqs;
	      tr->extendedPartitionData[i].optimizeBaseFrequencies  = tr->initialPartitionData[i].optimizeBaseFrequencies;
	      tr->extendedPartitionData[i].ascBias  = tr->initialPartitionData[i].ascBias;
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
      
      

      tr->executeModel   = (boolean *)rax_malloc(sizeof(boolean) * tr->NumberOfModels);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->executeModel[i] = TRUE;

      getyspace(rdta);
    } 

  setupTree(tr, adef);


  if(!adef->readTaxaOnly)
    {
      switch(adef->alignmentFileType)
	{
	case PHYLIP:
	  if(!getdata(adef, rdta, tr))
	    {
	      printf("Problem reading alignment file \n");
	      errorExit(1);
	    }
	  break;
	case FASTA:
	  parseFasta(adef, rdta, tr);
	  break;
	default:
	  assert(0);
	}
      
      tr->nameHash = initStringHashTable(10 * tr->mxtips);
      for(i = 1; i <= tr->mxtips; i++)
	{	 
	  addword(tr->nameList[i], tr->nameHash, i);
	}

      fclose(INFILE);
    }
}



static unsigned char buildStates(int secModel, unsigned char v1, unsigned char v2)
{
  unsigned char newChar = 0;

  switch(secModel)
    {
    case SECONDARY_DATA:
      newChar = v1;
      newChar = newChar << 4;
      newChar = newChar | v2;
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

	    newChar = n1;
	    newChar = newChar << 4;
	    newChar = newChar | n2;

	    intermediateBinaryStates[i] = newChar;
	  }

	newChar = v1;
	newChar = newChar << 4;
	newChar = newChar | v2;

	for(i = 0; i < length; i++)
	  {
	    if(newChar == intermediateBinaryStates[i])
	      break;
	  }
	if(i < length)
	  newChar = finalBinaryStates[i];
	else
	  {
	    newChar = 0;
	    for(i = 0; i < length; i++)
	      {
		if(v1 & meaningDNA[allowedStates[i][0]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    newChar |= finalBinaryStates[i];
		  }
		if(v2 & meaningDNA[allowedStates[i][1]])
		  {
		    /*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
		    newChar |= finalBinaryStates[i];
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

	    newChar = n1;
	    newChar = newChar << 4;
	    newChar = newChar | n2;

	    intermediateBinaryStates[i] = newChar;
	  }

	newChar = v1;
	newChar = newChar << 4;
	newChar = newChar | v2;

	for(i = 0; i < 6; i++)
	  {
	    /* exact match */
	    if(newChar == intermediateBinaryStates[i])
	      break;
	  }
	if(i < 6)
	  newChar = finalBinaryStates[i];
	else
	  {
	    /* distinguish between exact mismatches and partial mismatches */

	    for(i = 0; i < 6; i++)
	      if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		break;
	    if(i < 6)
	      {
		/* printf("partial mismatch\n"); */

		newChar = 0;
		for(i = 0; i < 6; i++)
		  {
		    if((v1 & meaningDNA[allowedStates[i][0]]) && (v2 & meaningDNA[allowedStates[i][1]]))
		      {
			/*printf("Adding %c%c\n", allowedStates[i][0], allowedStates[i][1]);*/
			newChar |= finalBinaryStates[i];
		      }
		    else
		      newChar |=  finalBinaryStates[6];
		  }
	      }
	    else
	      newChar = finalBinaryStates[6];
	  }	
      }
      break;
    default:
      assert(0);
    }

  return newChar;

}



static void adaptRdataToSecondary(tree *tr, rawdata *rdta)
{
  int *alias = (int*)rax_calloc(rdta->sites, sizeof(int));
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

  rax_free(alias);
}

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  
    gap, 
    i, 
    j, 
    jj, 
    jg, 
    k, 
    n, 
    nsp,  
    *index, 
    *category = (int*)NULL;

  boolean  
    flip, 
    tied;
  
  unsigned char  
    **data;

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


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef, int countAscBias)
{
  boolean  
    tied;
  
  int   
    i, 
    sitei, 
    j, 
    sitej, 
    k,
    *aliasModel = (int*)NULL,
    *aliasSuperModel = (int*)NULL,
    undeterminedSites = 0;

 

  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)rax_malloc(sizeof(int) * (rdta->sites + 1));
      aliasModel      = (int*)rax_malloc(sizeof(int) * (rdta->sites + 1));
    } 

  i = 0;
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;

  //if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES)
  {
    int i;
    
    tr->patternPosition = (int*)rax_malloc(sizeof(int) * rdta->sites);
    tr->columnPosition  = (int*)rax_malloc(sizeof(int) * rdta->sites);
    
    for(i = 0; i < rdta->sites; i++)
      {
	tr->patternPosition[i] = -1;
	tr->columnPosition[i]  = -1;
      }
  }

  i = 0;

  for (j = 1; j <= rdta->sites; j++)
    {
      int 
	allGap = TRUE;

      unsigned char 
	undetermined;
      
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];

      undetermined = getUndetermined(tr->dataVector[sitej]);

      for(k = 1; k <= rdta->numsp; k++)
	{	 
	  if(rdta->y[k][sitej] != undetermined)
	    {
	      allGap = FALSE;
	      break;
	    }
	}

      if(allGap)      
	undeterminedSites++;
      
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

      assert(!(tied && allGap));
      
      if(tied && !allGap)
	{	  
	  tr->patternPosition[j - 1] = i;
	  tr->columnPosition[j - 1] = sitej;
	  
	  cdta->aliaswgt[i] += rdta->wgt[sitej];

	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else
	{
	  if(!allGap)
	    {
	      if(cdta->aliaswgt[i] > 0) 
		i++;
	      	      
	      tr->patternPosition[j - 1] = i;
	      tr->columnPosition[j - 1] = sitej;
	
	      cdta->aliaswgt[i] = rdta->wgt[sitej];
	      cdta->alias[i] = sitej;
	      
	      if(adef->useMultipleModel)
		{
		  aliasModel[i]      = tr->model[sitej];
		  aliasSuperModel[i] = tr->dataVector[sitej];
		}
	    }	
	}
    }

  cdta->endsite = i;

  if (cdta->aliaswgt[i] > 0) 
    cdta->endsite++;

  if(adef->mode == PER_SITE_LL || adef->mode == ANCESTRAL_STATES || (countAscBias > 0))
    {
      if(undeterminedSites > 0)
	{
	  printBothOpen("You are trying to infer per site likelihoods or ancestral states or\n");
	  printBothOpen("do calculations with an ascertainment bias correction\n");
	  printBothOpen("on an alignment containing %d sites consisting only of undetermined\n", undeterminedSites);
	  printBothOpen("characters. Please remove them first and then re-run RAxML!\n");

	  errorExit(-1);
	}

      for(i = 0; i < rdta->sites; i++)
	{
	  int 
	    p  = tr->patternPosition[i],
	    c  = tr->columnPosition[i];
	  	  
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
      rax_free(aliasModel);
      rax_free(aliasSuperModel);
    }     

  if(undeterminedSites > 0)    
    printBothOpen("\nAlignment has %d completely undetermined sites that will be automatically removed from the input data\n\n", undeterminedSites);

  //exit(-1);
}


static boolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr, int countAscBias)
{
  int  i;

  for (i = 1; i <= rdta->sites; i++)
    cdta->alias[i] = i;

  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef, countAscBias);

  return TRUE;
}




static boolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)
{
  int  i, j, model, fullSites = 0, modelCounter;

  unsigned char
    *y    = (unsigned char *)rax_malloc(((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char)),
    *yBUF = (unsigned char *)rax_malloc( ((size_t)rdta->numsp) * ((size_t)cdta->endsite) * sizeof(unsigned char));

  for (i = 1; i <= rdta->numsp; i++)
    for (j = 0; j < cdta->endsite; j++)
      y[(((size_t)(i - 1)) * ((size_t)cdta->endsite)) + j] = rdta->y[i][cdta->alias[j]];

  rax_free(rdta->y0);
  rax_free(rdta->y);

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

      if(modelCounter <  tr->NumberOfModels - 1)
	{
	  printf("\nYou specified %d partitions, but after parsing and pre-processing ExaML only found %d partitions\n", tr->NumberOfModels, modelCounter + 1);
	  printf("Presumably one or more partitions vanished because they consisted entirely of undetermined characters.\n");
	  printf("Please fix your data!\n\n");
	  exit(-1);
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

  tr->invariant          = (int *)rax_malloc(cdta->endsite * sizeof(int));
  tr->originalDataVector = (int *)rax_malloc(cdta->endsite * sizeof(int));
  tr->originalModel      = (int *)rax_malloc(cdta->endsite * sizeof(int));
  tr->originalWeights    = (int *)rax_malloc(cdta->endsite * sizeof(int));

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
  int *omissionList     = (int *)rax_calloc(n, sizeof(int));
  int *undeterminedList = (int *)rax_calloc((rdta->sites + 1), sizeof(int));
  int *modelList        = (int *)rax_malloc((rdta->sites + 1)* sizeof(int));
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

	  if(processID == 0 && !adef->silent)
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
		      if(processID == 0 && !adef->silent)
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

      if(count > 0 && processID == 0 && !adef->silent)
	{
	  printBothOpen("\nIMPORTANT WARNING\n");

	  printBothOpen("Found %d %s that %s exactly identical to other sequences in the alignment.\n", count, (count == 1)?"sequence":"sequences", (count == 1)?"is":"are");

	  printBothOpen("Normally they should be excluded from the analysis.\n\n");
	}

      if(countUndeterminedColumns > 0 && processID == 0 && !adef->silent)
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
				if(tr->partitionData[i].ascBias)
				  {
				    strcpy(AAmodel, "ASC_");
				    strcat(AAmodel, protModels[tr->partitionData[i].protModels]);
				  }
				else
				  strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
				if(tr->partitionData[i].usePredefinedProtFreqs == FALSE)
				  strcat(AAmodel, "F");

				if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				  strcat(AAmodel, "X");

				assert(!(tr->partitionData[i].optimizeBaseFrequencies && tr->partitionData[i].usePredefinedProtFreqs));
				
				fprintf(newFile, "%s, ", AAmodel);
			      }
			    else
			      fprintf(newFile, "[%s], ", tr->partitionData[i].proteinSubstitutionFileName);
			  }
			  break;
			case DNA_DATA:	
			  if(tr->partitionData[i].ascBias)
			    {
			      if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				fprintf(newFile, "ASC_DNAX, ");
			      else
				fprintf(newFile, "ASC_DNA, ");
			    }
			  else
			    {
			      if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				fprintf(newFile, "DNAX, ");
			      else
				fprintf(newFile, "DNA, ");
			    }
			  break;
			case BINARY_DATA:	
			   if(tr->partitionData[i].ascBias)
			     {
			       if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				 fprintf(newFile, "ASC_BINX, ");
			       else
				 fprintf(newFile, "ASC_BIN, ");
			     }
			   else
			     {
			       if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				 fprintf(newFile, "BINX, ");
			       else
				 fprintf(newFile, "BIN, ");
			     }
			  break;
			case GENERIC_32:
			  if(tr->partitionData[i].ascBias)
			    {
			      if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				fprintf(newFile, "ASC_MULTIX, ");
			      else
				fprintf(newFile, "ASC_MULTI, ");
			    }
			  else
			    {
			       if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				fprintf(newFile, "MULTIX, ");
			       else
				fprintf(newFile, "MULTI, ");
			    }
			  break;
			case GENERIC_64:
			  if(tr->partitionData[i].ascBias)
			    {
			      if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				fprintf(newFile, "ASC_CODONX, ");
			      else
				fprintf(newFile, "ASC_CODON, ");
			    }
			  else
			    {
			       if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
				fprintf(newFile, "CODONX, ");
			      else
				fprintf(newFile, "CODON, ");
			    }
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


	  if(adef->printIdenticalSequences == TRUE)
	    {
	      count = 0;
	    }

	  if(!filexists(noDupFile))
	    {
	      FILE 
		*newFile;

	      if(adef->silent && (count || countUndeterminedColumns))
		printBothOpen("\nIMPORTANT WARNING: Alignment validation warnings have been suppressed. Found %d duplicate %s and %d undetermined %s\n\n", 
			      count, count > 1 ? "sequences" : "sequence", countUndeterminedColumns, countUndeterminedColumns > 1 ? "columns" : "column");
	      
	      //if(adef->printIdenticalSequences)
	      //	count = 0;
 	      
	      if(count > 0 || countUndeterminedColumns > 0)
		{
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
		      if(!omissionList[i] || count == 0)
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

  rax_free(undeterminedList);
  rax_free(omissionList);
  rax_free(modelList);
}





static void printPartitionFile(tree *tr, analdef *adef, char* newPartitionFile)
{
  if(adef->useMultipleModel && !filexists(newPartitionFile))
    {
      FILE 
	*newFile = myfopen(newPartitionFile, "wb");
     
      int 	
	i,
	l = 1,
	partitions = 0;

      printBothOpen("\n\nA partitioned model file with model assignments for bootstrap alignments \n");
      printBothOpen("is printed to file %s\n",newPartitionFile);
      printBothOpen("IMPORTANT: You MUST use this new model file and NOT the original one when running RAxML and ExaML on these bootstrapped alignments!\n\n");


      for(i = 1; i < tr->cdta->endsite; i++)
	assert(tr->model[i] >= tr->model[i-1]);	
	       
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  int 	   
	    lower, 
	    upper;

	    switch(tr->partitionData[i].dataType)
	      {
	      case AA_DATA:
		{
		  char
		    AAmodel[1024];

		  if(tr->partitionData[i].protModels != PROT_FILE)
		    {
		      if(tr->partitionData[i].ascBias)
			{
			  strcpy(AAmodel, "ASC_");
			  strcat(AAmodel, protModels[tr->partitionData[i].protModels]);
			}
		      else
			strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
		      if(tr->partitionData[i].usePredefinedProtFreqs == FALSE)
			strcat(AAmodel, "F");

		      if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
			strcat(AAmodel, "X");

		      assert(!(tr->partitionData[i].optimizeBaseFrequencies && tr->partitionData[i].usePredefinedProtFreqs));

		      fprintf(newFile, "%s, ", AAmodel);
		    }
		  else
		    fprintf(newFile, "[%s], ", tr->partitionData[i].proteinSubstitutionFileName);
		}
		break;
	      case DNA_DATA:
		if(tr->partitionData[i].ascBias)
		  {
		    if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		      fprintf(newFile, "ASC_DNAX, ");
		    else
		      fprintf(newFile, "ASC_DNA, ");
		  }
		else
		  {
		    if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		      fprintf(newFile, "DNAX, ");
		    else
		      fprintf(newFile, "DNA, ");
		  }
		break;
	      case BINARY_DATA:
		 if(tr->partitionData[i].ascBias)
		   {
		     if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		       fprintf(newFile, "ASC_BINX, ");
		     else
		       fprintf(newFile, "ASC_BIN, ");
		   }
		 else
		   {
		     if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		       fprintf(newFile, "BINX, ");
		     else
		       fprintf(newFile, "BIN, ");
		   }
		break;
	      case GENERIC_32:
		if(tr->partitionData[i].ascBias)
		  {
		    if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		      fprintf(newFile, "ASC_MULTIX, ");
		    else
		      fprintf(newFile, "ASC_MULTI, ");
		  }
		else
		  {
		     if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		      fprintf(newFile, "MULTIX, ");
		     else
		      fprintf(newFile, "MULTI, ");
		  }
		break;
	      case GENERIC_64:
		if(tr->partitionData[i].ascBias)
		  {
		    if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		      fprintf(newFile, "ASC_CODONX, ");
		    else
		      fprintf(newFile, "ASC_CODON, ");
		  }
		else
		  {
		     if(tr->partitionData[i].optimizeBaseFrequencies == TRUE)
		      fprintf(newFile, "CODONX, ");
		    else
		      fprintf(newFile, "CODON, ");
		  }
		break;
	      default:
		assert(0);
	      }

	    fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);

	    int 
	      k = 0;
	    
	    while(k < tr->cdta->endsite)
	      {
		if(tr->model[k] == i)
		  {
		    lower = l;
		    
		    do
		      {
			l += tr->cdta->aliaswgt[k];
		      }
		    while((++k < tr->cdta->endsite) && (tr->model[k] == i) );
		    
		    upper = l-1;

		    if(lower == upper)		      		       
		      fprintf(newFile, "%d", lower);		 
		    else
		      {
			assert(lower < upper);			
			fprintf(newFile, "%d-%d", lower, upper);		  
		      }		   
		    partitions++;
		  }
		else
		  k++;
	      }
            //printf("k: %d, cdta: %d\n", k, tr->cdta->endsite);
            assert(k == tr->cdta->endsite);	    
	    fprintf(newFile, "\n");
	}
      
      assert(partitions == tr->NumberOfModels);
      //printf("l:%d, rdta: %d\n", l, tr->rdta->sites);
      assert(l == tr->rdta->sites + 1);
      //assert(parts == tr->NumberOfModels);
      fclose(newFile);      
    }
  else
    {
      if(adef->useMultipleModel)
	{
	  printBothOpen("\nA partitioned model file with model assignments for bootstrap alignments\n");
	  printBothOpen("has already been printed to  file %s\n",newPartitionFile);
	}
    }
}

static void generateBS(tree *tr, analdef *adef)
{
  int 
    i, 
    j, 
    k, 
    w;
  
  char outName[1024], partName[1024], buf[16];
  FILE *of;

  assert(adef->boot != 0);

  {
    int 
      i,
      w = 0;

    for(i = 0; i < tr->cdta->endsite; i++)
      w += tr->cdta->aliaswgt[i];        

    if(w < tr->rdta->sites)
      {
	printBothOpen("Error in BS replicate generation. Apparently your input alignment contains %d completely undetermined sites.\n", tr->rdta->sites - w);
	printBothOpen("RAxML cowardly refuses to generate BS replicate MSAs on original MSAs containing entirely undetermined sites.\n\n");
	errorExit(-1);
      }
  }
  

  for(i = 0; i < adef->multipleRuns; i++)
    {
      int 
	count = 0;

      computeNextReplicate(tr, &adef->boot, (int*)NULL, (int*)NULL, FALSE, FALSE, adef);

      count = 0;
      for(j = 0; j < tr->cdta->endsite; j++)
	count += tr->cdta->aliaswgt[j];

      assert(count == tr->fullSites);

      /* generate model file name */
      strcpy(partName, workdir);
      strcat(partName, modelFileName);
      strcat(partName, ".BS");
      sprintf(buf, "%d", i);
      strcat(partName, buf);
      
      printPartitionFile(tr, adef, partName);
      /*******/
       

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
  int *modelFilter = (int *)rax_malloc(sizeof(int) * n);
  int length, k;
  unsigned char *tip;
  FILE *outf;
  char outFileName[2048];
  
  /* char buf[16]; */

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      strcpy(outFileName, seq_file);

      /*sprintf(buf, "%d", i);*/
      /*strcat(outFileName, ".GENE.");*/
      
      strcat(outFileName, ".");
      strcat(outFileName, tr->partitionData[i].partitionName);
      strcat(outFileName, ".phy");
      
      /*strcat(outFileName, buf);*/
      
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

  rax_free(modelFilter);
  printf("Wrote all %d individual gene/partition alignments\n", tr->NumberOfModels);
  printf("Exiting normally\n");
}


static int countTaxaInTopology(void)
{
  FILE 
    *f = myfopen(tree_file, "rb");   

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
      const partitionLengths 
	*pl = getPartitionLengths(&(tr->partitionData[i]));
            
      if(tr->useFastScaling)	
	tr->partitionData[i].globalScaler    = (unsigned int *)rax_calloc(2 * tr->mxtips, sizeof(unsigned int));  	         

      
      tr->partitionData[i].left              = (double *)rax_malloc(pl->leftLength * (maxCategories + 1) * sizeof(double));
      tr->partitionData[i].right             = (double *)rax_malloc(pl->rightLength * (maxCategories + 1) * sizeof(double));      
      tr->partitionData[i].EIGN              = (double*)rax_malloc(pl->eignLength * sizeof(double));
      tr->partitionData[i].EV                = (double*)rax_malloc(pl->evLength * sizeof(double));
      tr->partitionData[i].EI                = (double*)rax_malloc(pl->eiLength * sizeof(double));
      tr->partitionData[i].substRates        = (double *)rax_malloc(pl->substRatesLength * sizeof(double));
      tr->partitionData[i].frequencies       = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[i].freqExponents     = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
      tr->partitionData[i].tipVector         = (double *)rax_malloc(pl->tipVectorLength * sizeof(double));

      tr->partitionData[i].invariableFrequencies  = (double *)rax_malloc(pl->states * sizeof(double));      
      

#ifdef _HET
      
      tr->partitionData[i].EIGN_TIP          = (double*)rax_malloc(pl->eignLength * sizeof(double));
      tr->partitionData[i].EV_TIP            = (double*)rax_malloc(pl->evLength * sizeof(double));
      tr->partitionData[i].EI_TIP            = (double*)rax_malloc(pl->eiLength * sizeof(double));
      tr->partitionData[i].substRates_TIP    = (double *)rax_malloc(pl->substRatesLength * sizeof(double));      
      tr->partitionData[i].tipVector_TIP     = (double *)rax_malloc(pl->tipVectorLength * sizeof(double));

#endif
      


      if(tr->partitionData[i].protModels == LG4 || tr->partitionData[i].protModels == LG4X)      
	{	  	  
	  int 
	    k;
	  
	  for(k = 0; k < 4; k++)
	    {	    
	      tr->partitionData[i].EIGN_LG4[k]              = (double*)rax_malloc(pl->eignLength * sizeof(double));
	      tr->partitionData[i].rawEIGN_LG4[k]              = (double*)rax_malloc(pl->eignLength * sizeof(double));	      
	      tr->partitionData[i].EV_LG4[k]                = (double*)rax_malloc(pl->evLength * sizeof(double));
	      tr->partitionData[i].EI_LG4[k]                = (double*)rax_malloc(pl->eiLength * sizeof(double));
	      tr->partitionData[i].substRates_LG4[k]        = (double *)rax_malloc(pl->substRatesLength * sizeof(double));
	      tr->partitionData[i].frequencies_LG4[k]       = (double*)rax_malloc(pl->frequenciesLength * sizeof(double));
	      tr->partitionData[i].tipVector_LG4[k]         = (double *)rax_malloc(pl->tipVectorLength * sizeof(double));
	    }
	}

      

      tr->partitionData[i].symmetryVector    = (int *)rax_malloc(pl->symmetryVectorLength  * sizeof(int));
      tr->partitionData[i].frequencyGrouping = (int *)rax_malloc(pl->frequencyGroupingLength  * sizeof(int));
      tr->partitionData[i].perSiteRates      = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
      tr->partitionData[i].unscaled_perSiteRates = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
      
      
      tr->partitionData[i].nonGTR = FALSE;     
      

      tr->partitionData[i].gammaRates = (double*)rax_malloc(sizeof(double) * 4);
      tr->partitionData[i].yVector = (unsigned char **)rax_malloc(sizeof(unsigned char*) * (tr->mxtips + 1));

           
      tr->partitionData[i].xVector = (double **)rax_malloc(sizeof(double*) * tr->innerNodes);     
      tr->partitionData[i].xSpaceVector = (size_t *)rax_calloc(tr->innerNodes, sizeof(size_t));	
           
      tr->partitionData[i].expVector      = (int **)rax_malloc(sizeof(int*) * tr->innerNodes);
      tr->partitionData[i].expSpaceVector = (size_t *)rax_calloc(tr->innerNodes, sizeof(size_t));

      tr->partitionData[i].mxtips  = tr->mxtips;

     
      //andre-opt
      tr->partitionData[i].presenceMap = (unsigned int *)rax_calloc((size_t)tr->mxtips + 1 , sizeof(unsigned int));

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
	
      //asc
          
      if(tr->partitionData[model].ascBias)
	{	 
	  tr->partitionData[model].ascOffset = 4 * tr->partitionData[model].states * tr->partitionData[model].states;

	  tr->partitionData[model].ascVector = (double *)rax_malloc(((size_t)tr->innerNodes) *
								    ((size_t)tr->partitionData[model].ascOffset) * 
								    sizeof(double));
	  	 
	  tr->partitionData[model].ascExpVector = (int *)rax_calloc(((size_t)tr->innerNodes) * ((size_t)tr->partitionData[model].states),
								 sizeof(int));
	 
	  
	  tr->partitionData[model].ascSumBuffer = (double *)rax_malloc(((size_t)tr->partitionData[model].ascOffset) *
								       sizeof(double));
	  
	}
      
      //asc

      tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
      
      tr->partitionData[model].gapVector = (unsigned int*)rax_calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));


      tr->partitionData[model].initialGapVectorSize = tr->partitionData[model].gapVectorLength * 2 * tr->mxtips * sizeof(int);
	
      /* always multiply by 4 due to frequent switching between CAT and GAMMA in standard RAxML */
      
      tr->partitionData[model].gapColumn = (double *)rax_malloc(((size_t)tr->innerNodes) *
								    ((size_t)4) * 
								    ((size_t)(tr->partitionData[model].states)) *
								    sizeof(double));		  		
	
      

      undetermined = getUndetermined(tr->partitionData[model].dataType);

      for(j = 1; j <= tr->mxtips; j++)
	for(i = 0; i < width; i++)
	  if(tr->partitionData[model].yVector[j][i] == undetermined)
	    tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];      
    }

  tr->perSiteLL       = (double *)rax_malloc((size_t)tr->cdta->endsite * sizeof(double));
  assert(tr->perSiteLL != NULL);

  tr->sumBuffer  = (double *)rax_malloc(memoryRequirements * sizeof(double));
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
  adef->constraintSeed         = 0;
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
  adef->useBinaryModelFile     = FALSE;
  adef->leaveDropMode          = FALSE;
  adef->slidingWindowSize      = 100;
  adef->checkForUndeterminedSequences = TRUE;
  adef->useQuartetGrouping = FALSE;
  adef->alignmentFileType = PHYLIP;
  adef->calculateIC = FALSE;
  adef->verboseIC = FALSE;
  adef->stepwiseAdditionOnly = FALSE;
  adef->optimizeBaseFrequencies = FALSE;
  adef->ascertainmentBias = FALSE;
  adef->rellBootstrap = FALSE;
  adef->mesquite = FALSE;
  adef->silent = FALSE;
  adef->noSequenceCheck = FALSE;
  adef->useBFGS = TRUE;
  adef->setThreadAffinity = FALSE;
  adef->bootstopPermutations = 100;
  adef->fcThreshold = 99;
  adef->sampleQuartetsWithoutReplacement = FALSE;
  adef->printIdenticalSequences = FALSE; 
}

static int modelExists(char *model, analdef *adef)
{
  int 
    i;
  
  char 
    thisModel[1024];

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

  if(strcmp(model, "ASC_BINGAMMA\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
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

   if(strcmp(model, "ASC_BINCAT\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  


  if(strcmp(model, "BINGAMMAX\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_BINGAMMAX\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "BINCATX\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_BINCATX\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

   if(strcmp(model, "BINGAMMAIX\0") == 0)
    {
      adef->model = M_BINGAMMA;
      adef->useInvariant = TRUE;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

   if(strcmp(model, "BINCATIX\0") == 0)
    {
      adef->model = M_BINCAT;
      adef->useInvariant = TRUE;
      adef->optimizeBaseFrequencies = TRUE;
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
  
  if(strcmp(model, "ASC_MULTIGAMMA\0") == 0)
    {
      adef->model = M_32GAMMA;
      adef->useInvariant = FALSE; 
      adef->ascertainmentBias = TRUE;
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

   if(strcmp(model, "ASC_MULTICAT\0") == 0)
    {
      adef->model = M_32CAT;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  


  if(strcmp(model, "MULTIGAMMAX\0") == 0)
    {
      adef->model = M_32GAMMA;
      adef->useInvariant = FALSE; 
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_MULTIGAMMAX\0") == 0)
    {
      adef->model = M_32GAMMA;
      adef->useInvariant = FALSE; 
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "MULTICATX\0") == 0)
    {
      adef->model = M_32CAT;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_MULTICATX\0") == 0)
    {
      adef->model = M_32CAT;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }
  
   if(strcmp(model, "MULTIGAMMAIX\0") == 0)
    {
      adef->model = M_32GAMMA;
      adef->useInvariant = TRUE;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

    if(strcmp(model, "MULTICATIX\0") == 0)
      {
	adef->model = M_32CAT;
	adef->useInvariant = TRUE; 
	adef->optimizeBaseFrequencies = TRUE;
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

  if(strcmp(model, "ASC_CODONGAMMA\0") == 0)
    {
      adef->model = M_64GAMMA;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
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

  if(strcmp(model, "ASC_CODONCAT\0") == 0)
    {
      adef->model = M_64CAT;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

 


  if(strcmp(model, "CODONGAMMAX\0") == 0)
    {
      adef->model = M_64GAMMA;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;       
      return 1;
    }

  if(strcmp(model, "ASC_CODONGAMMAX\0") == 0)
    {
      adef->model = M_64GAMMA;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "CODONCATX\0") == 0)
    {
      adef->model = M_64CAT;
      adef->useInvariant = FALSE; 
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_CODONCATX\0") == 0)
    {
      adef->model = M_64CAT;
      adef->useInvariant = FALSE; 
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  
   if(strcmp(model, "CODONGAMMAIX\0") == 0)
    {
      adef->model = M_64GAMMA;
      adef->useInvariant = TRUE; 
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

   if(strcmp(model, "CODONCATIX\0") == 0)
    {
      adef->model = M_64CAT;
      adef->useInvariant = TRUE;
      adef->optimizeBaseFrequencies = TRUE;
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

  if(strcmp(model, "ASC_GTRGAMMA\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
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

   if(strcmp(model, "ASC_GTRCAT\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useInvariant = FALSE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }
   
 


  if(strcmp(model, "GTRGAMMAX\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;

      //printf("opt base freqs!\n");
      return 1;
    }
  

  if(strcmp(model, "ASC_GTRGAMMAX\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = FALSE;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      //printf("opt base freqs!\n");
      return 1;
    }

  if(strcmp(model, "GTRCATX\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useInvariant = FALSE; 
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_GTRCATX\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useInvariant = FALSE; 
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

   if(strcmp(model, "GTRGAMMAIX\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = TRUE;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

   if(strcmp(model, "GTRCATI\0") == 0)
     {
       adef->model = M_GTRCAT;
       adef->useInvariant = TRUE;
       adef->optimizeBaseFrequencies = TRUE;
       return 1;
     }

   if(strcmp(model, "GTRCATIX\0") == 0)
     {
       adef->model = M_GTRCAT;
       adef->useInvariant = TRUE;
       adef->optimizeBaseFrequencies = TRUE;
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

  if(strcmp(model, "ASC_PROTGAMMAGTR\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      adef->ascertainmentBias = TRUE;
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
  
   if(strcmp(model, "PROTCATGTRX\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

   if(strcmp(model, "ASC_PROTCATGTRX\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAGTRX\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0;
       adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }


  if(strcmp(model, "ASC_PROTGAMMAGTRX\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0;
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "PROTGAMMAIGTRX\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      adef->useInvariant = TRUE; 
      adef->protEmpiricalFreqs = 0;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  if(strcmp(model, "PROTCATIGTRX\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useInvariant = TRUE;
      adef->protEmpiricalFreqs = 0;
      adef->optimizeBaseFrequencies = TRUE;
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

  if(strcmp(model, "PROTCATGTR_UNLINKED_X\0") == 0)    
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");

      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0; 
      adef->optimizeBaseFrequencies = TRUE;
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

  if(strcmp(model, "ASC_PROTCATGTR_UNLINKED\0") == 0)    
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");

      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_PROTCATGTR_UNLINKED_X\0") == 0)    
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");

      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0; 
      adef->optimizeBaseFrequencies = TRUE;
      adef->ascertainmentBias = TRUE;
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

  if(strcmp(model, "PROTGAMMAGTR_UNLINKED_X\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0;
       adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

  
   if(strcmp(model, "ASC_PROTGAMMAGTR_UNLINKED\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 1;
      adef->ascertainmentBias = TRUE;
      return 1;
    }

  if(strcmp(model, "ASC_PROTGAMMAGTR_UNLINKED_X\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = FALSE;
      adef->protEmpiricalFreqs = 0;
       adef->optimizeBaseFrequencies = TRUE;
       adef->ascertainmentBias = TRUE;
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
  
  if(strcmp(model, "PROTGAMMAIGTR_UNLINKED_X\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = TRUE;
      adef->protEmpiricalFreqs = 0;
      adef->optimizeBaseFrequencies = TRUE;
      return 1;
    }

   if(strcmp(model, "PROTCATIGTR_UNLINKED_X\0") == 0)
    {
      printf("Advisory: GTR_UNLINKED only has an effect if specified in the partition file\n");
      
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR_UNLINKED;
      adef->useInvariant = TRUE;
      adef->protEmpiricalFreqs = 0;
       adef->optimizeBaseFrequencies = TRUE;
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

      /* check CATX */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "X");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 0;
	  adef->optimizeBaseFrequencies = TRUE;
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

      /* check CATIX */

      strcpy(thisModel, "PROTCATI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "X");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE; 
	  adef->optimizeBaseFrequencies = TRUE;
	  return 1;
	}

      /**************check CAT ASC ***********************************/


       /* check CAT */

      strcpy(thisModel, "ASC_PROTCAT");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->ascertainmentBias = TRUE;
	  return 1;
	}

      /* check CATF */

      strcpy(thisModel, "ASC_PROTCAT");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->ascertainmentBias = TRUE;
	  return 1;
	}

      /* check CATX */

      strcpy(thisModel, "ASC_PROTCAT");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "X");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 0;
	  adef->optimizeBaseFrequencies = TRUE;
	  adef->ascertainmentBias = TRUE;
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

      strcpy(thisModel, "ASC_PROTGAMMA");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->useInvariant = FALSE;
	  adef->ascertainmentBias = TRUE;
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

      strcpy(thisModel, "ASC_PROTGAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = FALSE;
	  adef->ascertainmentBias = TRUE;
	  return 1;
	}

     
       /* check GAMMAmodelX */

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "X");

      if(strcmp(model, thisModel) == 0)
	{	 
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 0;
	  adef->useInvariant = FALSE;
	  adef->optimizeBaseFrequencies = TRUE;	  
	  return 1;
	}

      strcpy(thisModel, "ASC_PROTGAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "X");

      if(strcmp(model, thisModel) == 0)
	{	 
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 0;
	  adef->useInvariant = FALSE;
	  adef->optimizeBaseFrequencies = TRUE;
	  adef->ascertainmentBias = TRUE;
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

      /* check GAMMAImodelX */

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "X");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 0;
	  adef->useInvariant = TRUE;
	  adef->optimizeBaseFrequencies = TRUE;
	  return 1;
	}

    }

  /*********************************************************************************/



  return 0;
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

  tr->outgroups = (char **)rax_malloc(sizeof(char *) * count);

  for(i = 0; i < tr->numberOfOutgroups; i++)
    tr->outgroups[i] = (char *)rax_malloc(sizeof(char) * nmlngth);

  tr->outgroupNums = (int *)rax_malloc(sizeof(int) * count);

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


static void printVersionInfo(boolean terminal, FILE *infoFile)
{
  char 
    text[12][1024];

  int 
    i;

  sprintf(text[0], "\n\nThis is %s version %s released by Alexandros Stamatakis on %s.\n\n",  programName, programVersion, programDate);
  sprintf(text[1], "With greatly appreciated code contributions by:\n");
  sprintf(text[2], "Andre Aberer      (HITS)\n");     
  sprintf(text[3], "Simon Berger      (HITS)\n"); 
  sprintf(text[4], "Alexey Kozlov     (HITS)\n"); 
  sprintf(text[5], "Kassian Kobert    (HITS)\n"); 
  sprintf(text[6], "David Dao         (KIT and HITS)\n");
  sprintf(text[7], "Sarah Lutteropp   (KIT and HITS)\n");
  sprintf(text[8], "Nick Pattengale   (Sandia)\n"); 
  sprintf(text[9], "Wayne Pfeiffer    (SDSC)\n");
  sprintf(text[10], "Akifumi S. Tanabe (NRIFS)\n");  
  sprintf(text[11], "Charlie Taylor    (UF)\n\n");
  

  for(i = 0; i < 12; i++)
    {
      if(terminal)    
	printf("%s", text[i]);
      else     
	printBoth(infoFile, text[i]);
    }
  
}

static void printMinusFUsage(void)
{
  printf("\n");
  printf("              \"-f a\": rapid Bootstrap analysis and search for best-scoring ML tree in one program run\n");  

  printf("              \"-f A\": compute marginal ancestral states on a ROOTED reference tree provided with \"-t\"\n");

  printf("              \"-f b\": draw bipartition information on a tree provided with \"-t\" based on multiple trees\n");
  printf("                      (e.g., from a bootstrap) in a file specified by \"-z\"\n");

  printf("              \"-f B\": optimize br-len scaler and other model parameters (GTR, alpha, etc.) on a tree provided with \"-t\".\n");
  printf("                      The tree needs to contain branch lengths. The branch lengths will not be optimized, just scaled by a single common value.\n");


  printf("              \"-f c\": check if the alignment can be properly read by RAxML\n");

  printf("              \"-f C\": ancestral sequence test for Jiajie, users will also need to provide a list of taxon names via -Y separated by whitespaces\n");

  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("                      DEFAULT: ON\n");

  printf("              \"-f D\": rapid hill-climbing with RELL bootstraps\n");

  printf("              \"-f e\": optimize model+branch lengths for given input tree under GAMMA/GAMMAI only\n");

  

  printf("              \"-f E\": execute very fast experimental tree search, at present only for testing\n");

  printf("              \"-f F\": execute fast experimental tree search, at present only for testing\n");

  printf("              \"-f g\": compute per site log Likelihoods for one or more trees passed via\n");
  printf("                      \"-z\" and write them to a file that can be read by CONSEL\n");
  printf("                      The model parameters will be estimated on the first tree only!\n");
  
  printf("              \"-f G\": compute per site log Likelihoods for one or more trees passed via\n");
  printf("                      \"-z\" and write them to a file that can be read by CONSEL.\n");
  printf("                      The model parameters will be re-estimated for each tree\n");
 
  printf("              \"-f h\": compute log likelihood test (SH-test) between best tree passed via \"-t\"\n");
  printf("                      and a bunch of other trees passed via \"-z\" \n");  
  printf("                      The model parameters will be estimated on the first tree only!\n");

  printf("              \"-f H\": compute log likelihood test (SH-test) between best tree passed via \"-t\"\n");
  printf("                      and a bunch of other trees passed via \"-z\" \n");  
  printf("                      The model parameters will be re-estimated for each tree\n");

  printf("              \"-f i\": calculate IC and TC scores (Salichos and Rokas 2013) on a tree provided with \"-t\" based on multiple trees\n");
  printf("                      (e.g., from a bootstrap) in a file specified by \"-z\"\n");

  printf("              \"-f I\": a simple tree rooting algorithm for unrooted trees.\n");
  printf("                      It roots the tree by rooting it at the branch that best balances the subtree lengths\n");
  printf("                      (sum over branches in the subtrees) of the left and right subtree.\n");
  printf("                      A branch with an optimal balance does not always exist!\n");
  printf("                      You need to specify the tree you want to root via \"-t\".\n"); 

  printf("              \"-f j\": generate a bunch of bootstrapped alignment files from an original alignemnt file.\n");
  printf("                      You need to specify a seed with \"-b\" and the number of replicates with \"-#\" \n");   

  printf("              \"-f J\": Compute SH-like support values on a given tree passed via \"-t\".\n"); 

  printf("              \"-f k\": Fix long branch lengths in partitioned data sets with missing data using the\n");
  printf("                      branch length stealing algorithm.\n");  
  printf("                      This option only works in conjunction with \"-t\", \"-M\", and \"-q\".\n");
  printf("                      It will print out a tree with shorter branch lengths, but having the same likelihood score.\n");

  printf("              \"-f m\": compare bipartitions between two bunches of trees passed via \"-t\" and \"-z\" \n");
  printf("                      respectively. This will return the Pearson correlation between all bipartitions found\n");
  printf("                      in the two tree files. A file called RAxML_bipartitionFrequencies.outpuFileName\n");
  printf("                      will be printed that contains the pair-wise bipartition frequencies of the two sets\n");

  printf("              \"-f n\": compute the log likelihood score of all trees contained in a tree file provided by\n");
  printf("                      \"-z\" under GAMMA or GAMMA+P-Invar\n");
  printf("                      The model parameters will be estimated on the first tree only!\n");

  printf("              \"-f N\": compute the log likelihood score of all trees contained in a tree file provided by\n");
  printf("                      \"-z\" under GAMMA or GAMMA+P-Invar\n");
  printf("                      The model parameters will be re-estimated for each tree\n");


  printf("              \"-f o\": old and slower rapid hill-climbing without heuristic cutoff\n");

  printf("              \"-f p\": perform pure stepwise MP addition of new sequences to an incomplete starting tree and exit\n");

  printf("              \"-f P\": perform a phylogenetic placement of sub trees specified in a file passed via \"-z\" into a given reference tree\n");
  printf("                      in which these subtrees are contained that is passed via \"-t\" using the evolutionary placement algorithm.\n");

  printf("              \"-f q\": fast quartet calculator\n");

  printf("              \"-f r\": compute pairwise Robinson-Foulds (RF) distances between all pairs of trees in a tree file passed via \"-z\" \n");
  printf("                      if the trees have node labales represented as integer support values the program will also compute two flavors of\n");
  printf("                      the weighted Robinson-Foulds (WRF) distance\n");

  printf("              \"-f R\": compute all pairwise Robinson-Foulds (RF) distances between a large reference tree  passed via \"-t\" \n");
  printf("                      and many smaller trees (that must have a subset of the taxa of the large tree) passed via \"-z\".\n");
  printf("                      This option is intended for checking the plausibility of very large phylogenies that can not be inspected\n");
  printf("                      visually any more.\n");

  printf("              \"-f s\": split up a multi-gene partitioned alignment into the respective subalignments \n");

  printf("              \"-f S\": compute site-specific placement bias using a leave one out test inspired by the evolutionary placement algorithm\n");

  printf("              \"-f t\": do randomized tree searches on one fixed starting tree\n");

  printf("              \"-f T\": do final thorough optimization of ML tree from rapid bootstrap search in stand-alone mode\n");

  printf("              \"-f u\": execute morphological weight calibration using maximum likelihood, this will return a weight vector.\n");
  printf("                      you need to provide a morphological alignment and a reference tree via \"-t\" \n");    

  printf("              \"-f v\": classify a bunch of environmental sequences into a reference tree using thorough read insertions\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");

  printf("              \"-f V\": classify a bunch of environmental sequences into a reference tree using thorough read insertions\n");
  printf("                      you will need to start RAxML with a non-comprehensive reference tree and an alignment containing all sequences (reference + query)\n");
  printf("                      WARNING: this is a test implementation for more efficient handling of multi-gene/whole-genome datasets!\n");

  printf("              \"-f w\": compute ELW test on a bunch of trees passed via \"-z\" \n");
  printf("                      The model parameters will be estimated on the first tree only!\n");

  printf("              \"-f W\": compute ELW test on a bunch of trees passed via \"-z\" \n");
  printf("                      The model parameters will be re-estimated for each tree\n");

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
  printVersionInfo(TRUE, (FILE*)NULL);
  printf("\n");
  printf("Please also consult the RAxML-manual\n");
  printf("\nPlease report bugs via the RAxML google group!\n");
  printf("Please send us all input files, the exact invocation, details of the HW and operating system,\n");
  printf("as well as all error messages printed to screen.\n\n\n");

  printf("raxmlHPC[-SSE3|-AVX|-PTHREADS|-PTHREADS-SSE3|-PTHREADS-AVX|-HYBRID|-HYBRID-SSE3|HYBRID-AVX]\n");
  printf("      -s sequenceFileName -n outputFileName -m substitutionModel\n");
  printf("      [-a weightFileName] [-A secondaryStructureSubstModel]\n");
  printf("      [-b bootstrapRandomNumberSeed] [-B wcCriterionThreshold]\n");
  printf("      [-c numberOfCategories] [-C] [-d] [-D]\n");
  printf("      [-e likelihoodEpsilon] [-E excludeFileName]\n");
  printf("      [-f a|A|b|B|c|C|d|D|e|E|F|g|G|h|H|i|I|j|J|k|m|n|N|o|p|P|q|r|R|s|S|t|T|u|v|V|w|W|x|y] [-F]\n");
  printf("      [-g groupingFileName] [-G placementThreshold] [-h] [-H]\n");
  printf("      [-i initialRearrangementSetting] [-I autoFC|autoMR|autoMRE|autoMRE_IGN]\n");
  printf("      [-j] [-J MR|MR_DROP|MRE|STRICT|STRICT_DROP|T_<PERCENT>] [-k] [-K] \n");
  printf("      [-L MR|MRE|T_<PERCENT>] [-M]\n");
  printf("      [-o outGroupName1[,outGroupName2[,...]]][-O]\n");
  printf("      [-p parsimonyRandomSeed] [-P proteinModel]\n");
  printf("      [-q multipleModelFileName] [-r binaryConstraintTree]\n");
  printf("      [-R binaryModelParamFile] [-S secondaryStructureFile] [-t userStartingTree]\n");
  printf("      [-T numberOfThreads] [-u] [-U] [-v] [-V] [-w outputDirectory] [-W slidingWindowSize]\n");
  printf("      [-x rapidBootstrapRandomNumberSeed] [-X] [-y] [-Y quartetGroupingFileName|ancestralSequenceCandidatesFileName]\n");
  printf("      [-z multipleTreesFile] [-#|-N numberOfRuns|autoFC|autoMR|autoMRE|autoMRE_IGN]\n");
#ifdef _WAYNE_MPI
  printf("      [--silent][--no-seq-check][--no-bfgs]\n");
#else
  printf("      [--mesquite][--silent][--no-seq-check][--no-bfgs]\n");
#endif
#ifdef _NICK
  printf("      [--asc-corr=stamatakis|felsenstein|lewis|goldman1|goldman2|goldman3]\n");
#else
  printf("      [--asc-corr=stamatakis|felsenstein|lewis]\n");
#endif
  printf("      [--flag-check][--auto-prot=ml|bic|aic|aicc]\n");
  printf("      [--epa-keep-placements=number][--epa-accumulated-threshold=threshold]\n");
  printf("      [--epa-prob-threshold=threshold]\n");
  printf("      [--JC69][--K80][--HKY85]\n");
#if (defined(_WAYNE_MPI) && defined(_USE_PTHREADS))
  printf("      [--set-thread-affinity]\n");
#endif
  printf("      [--bootstop-perms=number]\n");
  printf("\n");
  printf("      [--quartets-without-replacement]\n");
  printf("\n");
  printf("      [---without-replacement]\n");
  printf("\n");
  printf("      [--print-identical-sequences]\n");
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
  printf("      -c      Specify number of distinct rate catgories for RAxML when model of rate heterogeneity\n");
  printf("              is set to CAT\n");
  printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
  printf("              categories to accelerate computations. \n");
  printf("\n");
  printf("              DEFAULT: 25\n");
  printf("\n"); 
  printf("      -C      Enable verbose output for the \"-L\" and \"-f i\" options. This will produce more, as well as more verbose output files\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
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
  printf("              optimization of tree topology\n");
  printf("\n");
  printf("              DEFAULT: 0.1   for models not using proportion of invariant sites estimate\n");
  printf("                       0.001 for models using proportion of invariant sites estimate\n");
  printf("\n");
  printf("      -E      specify an exclude file name, that contains a specification of alignment positions you wish to exclude.\n");
  printf("              Format is similar to Nexus, the file shall contain entries like \"100-200 300-400\", to exclude a\n");
  printf("              single column write, e.g., \"100-100\", if you use a mixed model, an appropriately adapted model file\n");
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
  printf("      -H      Disable pattern compression.\n");
  printf("\n");
  printf("              DEFAULT: ON\n");
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
  printf("              or strict consensus tree with \"-J STRICT\". For a custom consensus threshold >= 50%%, specify T_<NUM>, where 100 >= NUM >= 50.\n");
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
  printf("      -L     Compute consensus trees labelled by IC supports and the overall TC value as proposed in Salichos and Rokas 2013.\n");
  printf("             Compute a majority rule consensus tree with \"-L MR\" or an extended majority rule consensus tree with \"-L MRE\".\n");
  printf("             For a custom consensus threshold >= 50%%, specify \"-L T_<NUM>\", where 100 >= NUM >= 50.\n");  
  printf("             You will of course also need to provide a tree file containing several UNROOTED trees via \"-z\"!\n");
  printf("\n");
  printf("      -m      Model of Binary (Morphological), Nucleotide, Multi-State, or Amino Acid Substitution: \n");
  printf("\n");
  printf("              BINARY:\n\n");
  printf("                \"-m BINCAT[X]\"       : Optimization of site-specific\n");
  printf("                                       evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                       rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                       automatically under BINGAMMA, depending on the tree search option.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m BINCATI[X]\"      : Optimization of site-specific\n");
  printf("                                       evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                       rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                       automatically under BINGAMMAI, depending on the tree search option.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m ASC_BINCAT[X]\"   : Optimization of site-specific\n");
  printf("                                       evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                       rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                       automatically under BINGAMMA, depending on the tree search option.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                                       The ASC prefix willl correct the likelihood for ascertainment bias.\n");  
  printf("                \"-m BINGAMMA[X]\"     : GAMMA model of rate heterogeneity (alpha parameter will be estimated).\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m ASC_BINGAMMA[X]\" : GAMMA model of rate heterogeneity (alpha parameter will be estimated).\n");
  printf("                                       The ASC prefix willl correct the likelihood for ascertainment bias.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m BINGAMMAI[X]\"    : Same as BINGAMMA, but with estimate of proportion of invariable sites.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("\n");
  printf("              NUCLEOTIDES:\n\n");
  printf("                \"-m GTRCAT[X]\"       : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                       evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                       rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                       under GTRGAMMA, depending on the tree search option.\n");  
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m GTRCATI[X]\"      : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                       evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                       rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                       under GTRGAMMAI, depending on the tree search option.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m ASC_GTRCAT[X]\"   : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                       evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                       rate categories for greater computational efficiency.  Final tree might be evaluated\n");
  printf("                                       under GTRGAMMA, depending on the tree search option.\n"); 
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                                       The ASC prefix willl correct the likelihood for ascertainment bias.\n");  
  printf("                \"-m GTRGAMMA[X]\"     : GTR + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                       heterogeneity (alpha parameter will be estimated).\n");  
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m ASC_GTRGAMMA[X]\" : GTR + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                       heterogeneity (alpha parameter will be estimated).\n");  
  printf("                                       The ASC prefix willl correct the likelihood for ascertainment bias.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m GTRGAMMAI[X]\"    : Same as GTRGAMMA, but with estimate of proportion of invariable sites.\n");
  printf("                                       With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("\n");
  printf("              MULTI-STATE:\n\n");
  printf("                \"-m MULTICAT[X]\"       : Optimization of site-specific\n");
  printf("                                         evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                         rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                         automatically under MULTIGAMMA, depending on the tree search option.\n"); 
  printf("                                         With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m MULTICATI[X]\"      : Optimization of site-specific\n");
  printf("                                         evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                         rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                         automatically under MULTIGAMMAI, depending on the tree search option.\n");
  printf("                                         With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n"); 
  printf("                \"-m ASC_MULTICAT[X]\"   : Optimization of site-specific\n");
  printf("                                         evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                         rate categories for greater computational efficiency. Final tree might be evaluated\n");
  printf("                                         automatically under MULTIGAMMA, depending on the tree search option.\n"); 
  printf("                                         With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                                         The ASC prefix willl correct the likelihood for ascertainment bias.\n");  
  printf("                \"-m MULTIGAMMA[X]\"     : GAMMA model of rate heterogeneity (alpha parameter will be estimated).\n");
  printf("                                         With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n"); 
  printf("                \"-m ASC_MULTIGAMMA[X]\" : GAMMA model of rate heterogeneity (alpha parameter will be estimated).\n");
  printf("                                         The ASC prefix willl correct the likelihood for ascertainment bias.\n");
  printf("                                         With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n"); 
  printf("                \"-m MULTIGAMMAI[X]\"    : Same as MULTIGAMMA, but with estimate of proportion of invariable sites.\n");
  printf("                                         With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("\n");
  printf("                You can use up to 32 distinct character states to encode multi-state regions, they must be used in the following order:\n");
  printf("                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, A, B, C, D, E, F, G, H, I, J, K, L, M, N, O, P, Q, R, S, T, U, V\n");
  printf("                i.e., if you have 6 distinct character states you would use 0, 1, 2, 3, 4, 5 to encode these.\n");
  printf("                The substitution model for the multi-state regions can be selected via the \"-K\" option\n");
  printf("\n");
  printf("              AMINO ACIDS:\n\n");
  printf("                \"-m PROTCATmatrixName[F|X]\"       : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                    evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                    rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                    automatically under PROTGAMMAmatrixName[F|X], depending on the tree search option.\n");  
  printf("                                                    With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m PROTCATImatrixName[F|X]\"      : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                    evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                    rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                    automatically under PROTGAMMAImatrixName[F|X], depending on the tree search option.\n");
  printf("                                                    With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m ASC_PROTCATmatrixName[F|X]\"   : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                    evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                    rate categories for greater computational efficiency.   Final tree might be evaluated\n");
  printf("                                                    automatically under PROTGAMMAmatrixName[F|X], depending on the tree search option.\n");  
  printf("                                                    With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n"); 
  printf("                                                    The ASC prefix willl correct the likelihood for ascertainment bias.\n");
  printf("                \"-m PROTGAMMAmatrixName[F|X]\"     : specified AA matrix + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                                    heterogeneity (alpha parameter will be estimated).\n"); 
  printf("                                                    With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("                \"-m ASC_PROTGAMMAmatrixName[F|X]\" : specified AA matrix + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                                    heterogeneity (alpha parameter will be estimated).\n"); 
  printf("                                                    The ASC prefix willl correct the likelihood for ascertainment bias.\n");
  printf("                                                    With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n"); 
  printf("                \"-m PROTGAMMAImatrixName[F|X]\"    : Same as PROTGAMMAmatrixName[F|X], but with estimate of proportion of invariable sites.\n");
  printf("                                                    With the optional \"X\" appendix you can specify a ML estimate of base frequencies.\n");  
  printf("\n");
  printf("                Available AA substitution models:\n");
  printf("                ");
  
  {
    int 
      i;
    
    for(i = 0; i < NUM_PROT_MODELS - 1; i++)
      {
	if(i > 0 && (i % 8 == 0))
	  {
	    printf("\n");
	    printf("                ");
	  }
	printf("%s, ", protModels[i]);
      }
    
    printf("%s\n", protModels[i]);
  }
  
  printf("                With the optional \"F\" appendix you can specify if you want to use empirical base frequencies.\n");
  printf("                AUTOF and AUTOX are not supported any more, if you specify AUTO it will test prot subst. models with and without empirical\n");
  printf("                base frequencies now!\n");
  printf("                Please note that for partitioned models you can in addition specify the per-gene AA model in\n");
  printf("                the partition file (see manual for details). Also note that if you estimate AA GTR parameters on a partitioned\n");
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
  printf("      -o      Specify the name of a single outgroup or a comma-separated list of outgroups, eg \"-o Rat\" \n");
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
  printf("      -X      Same as the \"-y\" option below, however the parsimony search is more superficial.\n");
  printf("              RAxML will only do a randomized stepwise addition order parsimony tree reconstruction\n");
  printf("              without performing any additional SPRs.\n");
  printf("              This may be helpful for very broad whole-genome datasets, since this can generate topologically\n");
  printf("              more different starting trees.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");   
  printf("      -y      If you want to only compute a parsimony starting tree with RAxML specify \"-y\",\n");
  printf("              the program will exit after computation of the starting tree\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n"); 
  printf("      -Y      Pass a quartet grouping file name defining four groups from which to draw quartets\n");
  printf("              The file input format must contain 4 groups in the following form:\n");
  printf("              (Chicken, Human, Loach), (Cow, Carp), (Mouse, Rat, Seal), (Whale, Frog);\n");
  printf("              Only works in combination with -f q !\n");
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
  printf("\n");  
#ifndef _WAYNE_MPI
  printf("      --mesquite Print output files that can be parsed by Mesquite.\n");
  printf("\n");
  printf("              DEFAULT: Off\n");
  printf("\n");
#endif
  printf("      --silent Disables printout of warnings related to identical sequences and entirely undetermined sites in the alignment\n");
  printf("\n");
  printf("              DEFAULT: Off\n");
  printf("\n");
  printf("      --no-seq-check Disables checking the input MSA for identical sequences and entirely undetermined sites.\n");
  printf("                     Enabling this option may save time, in particular for large phylogenomic alignments.\n");
  printf("                     Before using this, make sure to check the alignment using the \"-f c\" option!\n");
  printf("\n");
  printf("              DEFAULT: Off \n");
  printf("\n");
  printf("      --no-bfgs Disables automatic usage of BFGS method to optimize GTR rates on unpartitioned DNA datasets\n");
  printf("\n");
  printf("              DEFAULT: BFGS on\n");
  printf("\n");
#ifdef _NICK
  printf("      --asc-corr Allows to specify the type of ascertainment bias correction you wish to use. There are %d\n", NUM_ASC_CORRECTIONS);
#else
  printf("      --asc-corr Allows to specify the type of ascertainment bias correction you wish to use. There are %d\n", 3);
#endif
  printf("                 types available:\n");
  printf("                 --asc-corr=lewis: the standard correction by Paul Lewis\n");
  printf("                 --asc-corr=felsenstein: a correction introduced by Joe Felsenstein that allows to explicitely specify\n");
  printf("                                         the number of invariable sites (if known) one wants to correct for.\n");
  printf("                 --asc-corr=stamatakis: a correction introduced by myself that allows to explicitely specify\n");
  printf("                                        the number of invariable sites for each character (if known) one wants to correct for.\n");
#ifdef _NICK
  printf("                 --asc-corr=goldman1: 1st correction proposed by Nick Goldman\n");
  printf("                 --asc-corr=goldman2: 2nd correction proposed by Nick Goldman\n");
  printf("                 --asc-corr=goldman3: 3rd correction proposed by Nick Goldman\n");
  printf("                 For further details about the Goldman corrections please refer to the manual!\n");
#endif
  printf("\n");
  printf("      --flag-check When using this option, RAxML will only check if all command line flags specifed are available and then exit\n");
  printf("                   with a message listing all invalid command line flags or with a message stating that all flags are valid.\n");
  printf("\n");
  printf("      --auto-prot=ml|bic|aic|aicc When using automatic protein model selection you can chose the criterion for selecting these models.\n");
  printf("                  RAxML will test all available prot subst. models except for LG4M, LG4X and GTR-based models, with and without empirical base frequencies.\n");
  printf("                  You can chose between ML score based selection and the BIC, AIC, and AICc criteria.\n");
  printf("\n");
  printf("                  DEFAULT: ml\n");
  printf("\n");
  printf("      --epa-keep-placements=number specify the number of potential placements you want to keep for each read in the EPA algorithm.\n");
  printf("                  Note that, the actual values printed will also depend on the settings for --epa-prob-threshold=threshold !\n");
  printf("\n");
  printf("                  DEFAULT: 7\n");
  printf("\n");  
  printf("      --epa-prob-threshold=threshold specify a percent threshold for including potential placements of a read depending on the \n");
  printf("                  maximum placement weight for this read. If you set this value to 0.01 placements that have a placement weight of 1 per cent of\n");
  printf("                  the maximum placement will still be printed to file if the setting of --epa-keep-placements allows for it\n");
  printf("\n");
  printf("                  DEFAULT: 0.01\n");
  printf("\n");
  printf("      --epa-accumulated-threshold=threshold specify an accumulated likelihood weight threshold for which different placements of read are printed\n");
  printf("                  to file. Placements for a read will be printed until the sum of their placement weights has reached the threshold value.\n");
  printf("                  Note that, this option can neither be used in combination with --epa-prob-threshold nor with --epa-keep-placements!\n");
  printf("\n");
  printf("      --JC69 specify that all DNA partitions will evolve under the Jukes-Cantor model, this overrides all other model specifications for DNA partitions.\n");
  printf("\n");
  printf("                  DEFAULT: Off\n");
  printf("\n");
  printf("      --K80 specify that all DNA partitions will evolve under the K80 model, this overrides all other model specifications for DNA partitions.\n");
  printf("\n");
  printf("                  DEFAULT: Off\n");
  printf("\n");
   printf("      --HKY85 specify that all DNA partitions will evolve under the HKY85 model, this overrides all other model specifications for DNA partitions.\n");
  printf("\n");
  printf("                  DEFAULT: Off\n");
  printf("\n");  
#if (defined(_WAYNE_MPI) && defined(_USE_PTHREADS))
  printf("\n");
  printf("      --set-thread-affinity specify that thread-to-core affinity shall be set by RAxML for the hybrid MPI-PThreads version\n");
  printf("\n");
  printf("                  DEFAULT: Off\n");  
  printf("\n");
#endif
  printf("\n");
  printf("      --bootstop-perms=number specify the number of permutations to be conducted for the bootstopping/bootstrap convergence test.\n");
  printf("                  The allowed minimum number is 100!\n");
  printf("\n");
  printf("                  DEFAULT: 100\n"); 
  printf("\n");
  printf("      --quartets-without-replacement specify that quartets are randomly subsampled, but without replacement.\n");
   printf("\n");
  printf("                  DEFAULT: random sampling with replacements\n"); 
  printf("\n");
  printf("      --print-identical-sequences specify that RAxML shall automatically generate a .reduced alignment with all\n");
  printf("                  undetermined columns removed, but without removing exactly identical sequences\n");
  printf("\n");
  printf("                  DEFAULT: identical sequences will also be removed in the .reduced file\n"); 
  printf("\n");
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
    disablePatternCompression = FALSE,   
    resultDirSet = FALSE, 
    epaSet = FALSE;

  char
    resultDir[1024] = "",
    aut[256],           
    model[2048] = "",
    secondaryModel[2048] = "",
    multiStateModel[2048] = "",
    modelChar;

  double 
    likelihoodEpsilon,    
    wcThreshold,
    fastEPAthreshold;
  
  int
    fOptionCount = 0,
    invalidOptions = 0,
    i,
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
    multipleRunsSet = FALSE,
    yFileSet = FALSE,
    flagCheck = FALSE;

  FILE 
    *flagCheckFile;

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
  

  tr->doSubtreeEPA = FALSE;
  tr->useFastScaling = TRUE; 
  tr->bootStopCriterion = -1;
  tr->wcThreshold = 0.03;
  tr->doCutoff = TRUE;
  tr->secondaryStructureModel = SEC_16; /* default setting */
  tr->searchConvergenceCriterion = FALSE;
  tr->catOnly = FALSE;
  tr->useEpaHeuristics = FALSE;
  tr->fastEPAthreshold = -1.0;
  tr->multiStateModel  = GTR_MULTI_STATE;
  tr->saveMemory = FALSE;
  tr->useGammaMedian = FALSE;
  tr->noRateHet = FALSE;
  tr->perPartitionEPA = FALSE;
  tr->useBrLenScaler = FALSE;
  tr->ascertainmentCorrectionType = NOT_DEFINED;
  tr->autoProteinSelectionType = AUTO_ML;
  
  //EPA related stuff 

  tr->numberOfEPAEntries = 7;
  tr->accumulatedEPACutoff = 0.95;
  tr->useAccumulatedEPACutoff = FALSE;
  tr->probThresholdEPA = 0.01;

  //JC and K80 and HKY85

  tr->useK80 = FALSE;
  tr->useJC69 = FALSE;
  tr->useHKY85 = FALSE;

#ifdef _BASTIEN
  tr->doBastienStuff = FALSE;
#endif

  /********* tr inits end*************/

  for(i = 1; i < argc; i++)
    if(strcmp(argv[i], "--flag-check") == 0)
      {
	flagCheck = TRUE; 
	flagCheckFile = myfopen("RAxML_flagCheck", "wb");
	opterr = 0;
	break;
      }


  static int flag;

  while(1)
    {      
      static struct 
	option long_options[18] =
	{	 
	  {"mesquite",                  no_argument,       &flag, 1},
	  {"silent",                    no_argument,       &flag, 1},
	  {"no-seq-check",              no_argument,       &flag, 1},
	  {"no-bfgs",                   no_argument,       &flag, 1},
	  {"asc-corr",                  required_argument, &flag, 1},
	  {"flag-check",                no_argument,       &flag, 1},
	  {"auto-prot",                 required_argument, &flag, 1},
	  {"epa-keep-placements",       required_argument, &flag, 1},
	  {"epa-accumulated-threshold", required_argument, &flag, 1},
	  {"epa-prob-threshold",        required_argument, &flag, 1}, 
	  {"JC69",                      no_argument,       &flag, 1},
	  {"K80",                       no_argument,       &flag, 1},
	  {"HKY85",                     no_argument,       &flag, 1},	 	 	 
	  {"set-thread-affinity",       no_argument,       &flag, 1},
	  {"bootstop-perms",            required_argument, &flag, 1},
	  {"quartets-without-replacement", no_argument,    &flag, 1},
	  {"print-identical-sequences", no_argument,       &flag, 1},
	  {0, 0, 0, 0}
	};
      
      int 
	option_index;
      
      flag = 0;
      
      c = getopt_long(argc,argv, "R:T:E:N:B:L:P:S:Y:A:G:I:J:K:W:l:x:z:g:r:e:a:b:c:f:i:m:t:w:s:n:o:q:#:p:vudyjhHkMDFQUOVCX", long_options, &option_index/*&optind, &optarg*/);
         
      if(c == -1)
	break;          
      
      if(flag > 0)
	{
	  switch(option_index)
	    {
	    case 0:	    
	      adef->mesquite = TRUE;
#ifdef _WAYNE_MPI
	     if(processID == 0)
	       printf("\nMesquite output option not supported for hybrid version!\n");
	     errorExit(-1);
#endif
	      //printf("%s mes\n", long_options[option_index].name);
	      break;
	    case 1:
	      adef->silent = TRUE;
	      //printf("%s sil\n", long_options[option_index].name);
	      break;
	    case 2:
	      adef->noSequenceCheck = TRUE;
	      break;
	    case 3:
	      adef->useBFGS = FALSE;
	      break;
	    case 4:
	      {		
		char 
		  *ascModels[NUM_ASC_CORRECTIONS] = {"lewis", "felsenstein", "stamatakis", "goldman1", "goldman2", "goldman3"};

		int 
		  k;

		for(k = 0; k < NUM_ASC_CORRECTIONS; k++)		  
		  if(strcmp(optarg, ascModels[k]) == 0)
		    break;

		if(k == NUM_ASC_CORRECTIONS)
		  {
		    printf("\nError, unknown ascertainment correction type, you can specify one of the following corrections:\n\n");
		    for(k = 0; k < NUM_ASC_CORRECTIONS; k++)
		      printf("--asc-corr=%s\n", ascModels[k]);
		    printf("\n");
		    errorExit(-1);
		  }
		else
		  {
		    switch(k)
		      {
		      case 0:
			tr->ascertainmentCorrectionType = LEWIS_CORRECTION;
			break;
		      case 1:
			tr->ascertainmentCorrectionType = FELSENSTEIN_CORRECTION;
			break;
		      case 2:
			tr->ascertainmentCorrectionType = STAMATAKIS_CORRECTION;
			break;
		      case 3:
			tr->ascertainmentCorrectionType = GOLDMAN_CORRECTION_1;
			break;
		      case 4:
			tr->ascertainmentCorrectionType = GOLDMAN_CORRECTION_2;
			break;
		      case 5:
			tr->ascertainmentCorrectionType = GOLDMAN_CORRECTION_3;
			break;
		      default:
			assert(0);
		      }
		  }
	      }
	      break;
	    case 5:
	      break;
	    case 6:
	      {
		char 
		  *autoModels[4] = {"ml", "bic", "aic", "aicc"};

		int 
		  k;

		for(k = 0; k < 4; k++)		  
		  if(strcmp(optarg, autoModels[k]) == 0)
		    break;

		if(k == 4)
		  {
		    printf("\nError, unknown protein model selection type, you can specify one of the following selection criteria:\n\n");
		    for(k = 0; k < 4; k++)
		      printf("--auto-prot=%s\n", autoModels[k]);
		    printf("\n");
		    errorExit(-1);
		  }
		else
		  {
		    switch(k)
		      {
		      case 0:
			tr->autoProteinSelectionType = AUTO_ML;
			break;
		      case 1:
			tr->autoProteinSelectionType = AUTO_BIC;
			break;
		      case 2:
			tr->autoProteinSelectionType = AUTO_AIC;
			break;
		      case 3:
			tr->autoProteinSelectionType = AUTO_AICC;
			break;
		      default:
			assert(0);
		      }
		  }
	      }
	      break;
	    case 7:	      	 
	      if(sscanf(optarg,"%u", &(tr->numberOfEPAEntries)) != 1)
		{
		  printf("\nError parsing number of EPA placements to print, RAxML expects a positive integer value\n\n");
		  errorExit(-1);
		}
	     	     
	      if(tr->numberOfEPAEntries == 0)
		{
		  printf("\nError parsing number of EPA placements to print, RAxML expects an integer value larger than 0\n\n");
		  errorExit(-1);
		}

	      epaSet = TRUE;
	      break;
	    case 8:
	      if(sscanf(optarg, "%lf", &(tr->accumulatedEPACutoff)) != 1)
		{
		  printf("\nError parsing accumulated EPA placement weight cutoff, RAxML expects a floating point value > 0.0 and <= 1.0\n\n");
		  errorExit(-1);
		}

	      if(tr->accumulatedEPACutoff <= 0.0 || tr->accumulatedEPACutoff > 1.0)
		{
		  printf("\nError parsing accumulated EPA placement weight cutoff, RAxML expects a floating point value > 0.0 and <= 1.0\n\n");
		  errorExit(-1);
		}
	     
	      tr->useAccumulatedEPACutoff = TRUE;
	      break;
	    case 9:
	      if(sscanf(optarg, "%lf", &(tr->probThresholdEPA)) != 1)
		{
		  printf("\nError parsing EPA relative probability cutoff, RAxML expects a floating point value >= 0.0 and <= 1.0\n\n");
		  errorExit(-1);
		}
	      if(tr->probThresholdEPA < 0.0 || tr->probThresholdEPA > 1.0)
		{
		  printf("\nError parsing accumulated epa placement weight cutoff, RAxML expects a floating point value >= 0.0 and <= 1.0\n\n");
		  errorExit(-1);
		}	    

	      epaSet = TRUE;
	      break;
	    case 10:
	      tr->useJC69 = TRUE;
	      break;
	    case 11:
	      tr->useK80 = TRUE;
	      break;
	    case 12:
	       tr->useHKY85 = TRUE;
	      break;
	      break;	   
	    case 13:
#if (defined(_WAYNE_MPI) && defined(_USE_PTHREADS))
	      adef->setThreadAffinity = TRUE;
#else
	      printf("Warning: flag --set-thread-affinity has no effect if you don't use the hybrid MPI-PThreads version\n");
#endif
	      break;
	    case 14:
	      {
		int 
		  perms = -1;
		
		if(sscanf(optarg,"%d", &perms) != 1 || (perms < 100))
		  {
		    printf("\nError parsing number of bootstop permutations to execute, RAxML expects a positive integer value larger or equal to 100\n\n");
		    errorExit(-1);
		  }
		
		adef->bootstopPermutations = perms;
		adef->fcThreshold = perms - round((double)perms / 100.0);		
	      }
	      break;
	    case 15:
	      adef->sampleQuartetsWithoutReplacement = TRUE;
	      break;
	    case 16:
	      adef->printIdenticalSequences = TRUE;
	      break;
	    default:
	      if(flagCheck)
		{
		  printf("Option %s not supported\n", argv[optind - 1]);
		  fprintf(flagCheckFile, "Option %s not supported\n", argv[optind - 1]);
		  invalidOptions++;
		}
	      else
		assert(0);
	    }

	  if((tr->useK80 + tr->useJC69 + tr->useHKY85) > 1)
	    {
	      printf("\nYou can't use \"--JC69\" and \"--K80\" and \"--HKY85\" options simultaneously!\n");
	      printf("You can either use  \"--JC69\" or \"--K80\" or \"--HKY85\"!\n");
	      printf("exiting .....\n\n");
	      errorExit(-1);
	    }

	}
      else
	switch(c)
	  {           
	  case 'Y':
	    adef->useQuartetGrouping = TRUE;
	    yFileSet = TRUE;
	    strcpy(quartetGroupingFileName, optarg);
	    break;
	  case 'V':
	    tr->noRateHet = TRUE;
	    break;
	  case 'u':
	    tr->useGammaMedian = TRUE;
	    break;
	  case 'O':
	    adef->checkForUndeterminedSequences = FALSE;
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
	    adef->calculateIC = FALSE;
	    
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
		if( (sscanf( optarg, "%s", aut) > 0)  && optarg[0] == 'T' && optarg[1] == '_')
		  {
		    tr->consensusType = USER_DEFINED;
		    sscanf(optarg + 2,"%d", &tr->consensusUserThreshold);
		    
		    if(tr->consensusUserThreshold < 50 || tr->consensusUserThreshold > 100)
		      {
			printf("Please specify a custom threshold c, with 50 <= c <= 100\n" );
			errorExit(0); 
		      }
		  }
		else
		  {
		    if(processID == 0)	      
		      printf("Use -J consensus tree option either as \"-J MR\" or \"-J MRE\" or \"-J STRICT\" or \"-J MR_DROP\"  or \"-J STRICT_DROP\" or T_<NUM>, where NUM >= 50\n");
		    errorExit(0);
		  }	
	      }	     
	    break;
	  case 'C':
	    adef->verboseIC = TRUE;
	    break;
	  case 'L':
	    adef->readTaxaOnly = TRUE;
	    adef->mode = CONSENSUS_ONLY;
	    adef->leaveDropMode = FALSE;
	    adef->calculateIC = TRUE;
	    
	    if((sscanf(optarg,"%s", aut) > 0) && ((strcmp(aut, "MR") == 0) || (strcmp(aut, "MRE") == 0)))
	      {
		if((strcmp(aut, "MR") == 0))	   
		  tr->consensusType = MR_CONSENSUS;	    
		
		if((strcmp(aut, "MRE") == 0))		  	    
		  tr->consensusType = MRE_CONSENSUS;	    	   		    
	      }
	    else
	      {
		if((sscanf( optarg, "%s", aut) > 0)  && optarg[0] == 'T' && optarg[1] == '_')
		  {
		    tr->consensusType = USER_DEFINED;
		    sscanf(optarg + 2,"%d", &tr->consensusUserThreshold);
		    
		    if(tr->consensusUserThreshold < 50 || tr->consensusUserThreshold > 100)
		      {
			printf("Please specify a custom threshold c, with 50 <= c <= 100\n" );
			errorExit(0); 
		      }
		  }
		else
		  {
		    if(processID == 0)	      
		      printf("Use -L consensus tree option including IC/TC score computation either as \"-L MR\" or \"-L MRE\" or \"-L T_<NUM>\", where NUM >= 50\n");
		    errorExit(0);
		  }	
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
	      char 
		*outgroups = (char*)rax_malloc(sizeof(char) * (strlen(optarg) + 1));
	      strcpy(outgroups, optarg);
	      parseOutgroups(outgroups, tr);
	      rax_free(outgroups);
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
	    sscanf(optarg,"%" PRId64, &(adef->parsimonySeed));	
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
	    printVersionInfo(TRUE, (FILE*)NULL);
	    errorExit(0);
	  case 'y':
	    adef->stepwiseAdditionOnly = FALSE;
	    adef->startingTreeOnly = 1;
	    break;
	  case 'X':
	    adef->stepwiseAdditionOnly = TRUE;
	    adef->startingTreeOnly = 1;
	    break;     
	  case 'h':
	    printREADME();
	    errorExit(0);
	  case 'H':
	    disablePatternCompression = TRUE;
	    break;
	  case 'j':
	    adef->checkpoints = 1;
	    break;
	  case 'a':
	    strcpy(weightFileName,optarg);
	    adef->useWeightFile = TRUE;
	    break;
	  case 'b':
	    sscanf(optarg,"%" PRId64, &adef->boot);
	    if(adef->boot <= 0)
	      {
		printf("Bootstrap seed specified via -b must be greater than zero\n");
		errorExit(-1);
	      }
	    bSeedSet = TRUE;
	    break;
	  case 'x':
	    sscanf(optarg,"%" PRId64, &adef->rapidBoot);
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
	    fOptionCount++;
	    if(fOptionCount > 1) 
	      {
		printf("\nError: only one of the various \"-f \" options can be used per RAxML run!\n");
		printf("They are mutually exclusive! exiting ...\n\n");
		errorExit(-1);
	      }
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
		adef->perGeneBranchLengths = TRUE;
		tr->useBrLenScaler = TRUE;
		break;
	      case 'c':
		adef->mode = CHECK_ALIGNMENT;
		break;
	      case 'C':
		adef->mode = ANCESTRAL_SEQUENCE_TEST;
		tr->useFastScaling = FALSE;
		break;
	      case 'd':
		adef->mode = BIG_RAPID_MODE;
		tr->doCutoff = TRUE;
		break;
	      case 'D':
		adef->mode = BIG_RAPID_MODE;
		tr->doCutoff = TRUE;	
		tr->useFastScaling = FALSE;
		adef->rellBootstrap = TRUE;
		break;
	      case 'e':
		adef->mode = TREE_EVALUATION;
		break; 
	      case 'E':
		adef->mode = FAST_SEARCH;
		adef->veryFast = TRUE;
		break;
	      case 'F':
		adef->mode = FAST_SEARCH;
		adef->veryFast = FALSE;
		break;	 
	      case 'g':
		tr->useFastScaling = FALSE;
		tr->optimizeAllTrees = FALSE;    
		adef->mode = PER_SITE_LL;
		break;
	      case 'G':
		tr->useFastScaling = FALSE;
		tr->optimizeAllTrees = TRUE;
		adef->mode = PER_SITE_LL;
		break;
	      case 'h':
		tr->optimizeAllTrees = FALSE;
		adef->mode = TREE_EVALUATION;
		adef->likelihoodTest = TRUE;
		tr->useFastScaling = FALSE;
		break;	 
	      case 'H': 
		tr->optimizeAllTrees = TRUE;
		adef->mode = TREE_EVALUATION;
		adef->likelihoodTest = TRUE;
		tr->useFastScaling = FALSE;
		break;
	      case 'i':	    
		adef->readTaxaOnly = TRUE;
		adef->mode = CALC_BIPARTITIONS_IC;
		break;
	      case 'I':
		adef->mode = ROOT_TREE;
		adef->readTaxaOnly = TRUE;
		break;
	      case 'j':
		adef->mode = GENERATE_BS;
		adef->generateBS = TRUE;
		break;
	      case 'J':
		adef->mode = SH_LIKE_SUPPORTS; 
		tr->useFastScaling = FALSE;
		break;
	      case 'k':
		adef->mode = STEAL_BRANCH_LENGTHS;
		tr->useFastScaling = FALSE;	    
		adef->compressPatterns  = FALSE;
		break;
	      case 'm': 
		adef->readTaxaOnly = TRUE;	    
		adef->mode = COMPUTE_BIPARTITION_CORRELATION;
		break;
	      case 'n': 
		tr->optimizeAllTrees = FALSE;
		adef->mode = COMPUTE_LHS;
		break;
	      case 'N':
		tr->optimizeAllTrees = TRUE;
		adef->mode = COMPUTE_LHS;
		break;	    
	      case 'o':
		adef->mode = BIG_RAPID_MODE;
		tr->doCutoff = FALSE;
		break;
	      case 'p':
		adef->mode =  PARSIMONY_ADDITION;
		break;
	      case 'P':
		adef->mode =  SUBTREE_EPA;
		tr->doSubtreeEPA = TRUE;
		break;
	      case 'q':
		adef->mode = QUARTET_CALCULATION;
		break;	  	 
	      case 'r':
		adef->readTaxaOnly = TRUE;
		adef->mode = COMPUTE_RF_DISTANCE;
		break;	  
	      case 'R':
		adef->readTaxaOnly = TRUE;
		adef->mode = PLAUSIBILITY_CHECKER;
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
		tr->perPartitionEPA = FALSE;
		
#ifdef _USE_PTHREADS
		tr->useFastScaling = FALSE;
#endif
		break;
		
	      case 'V':
		adef->mode = CLASSIFY_ML;	   	   	   	    
		tr->perPartitionEPA = TRUE;
		
#ifdef _USE_PTHREADS
		tr->useFastScaling = FALSE;
#endif	    
		break;	  
	      case 'w':	    
		adef->mode = COMPUTE_ELW;
		adef->computeELW = TRUE;
		tr->optimizeAllTrees = FALSE;
		break;
	      case 'W':
		adef->mode = COMPUTE_ELW;
		adef->computeELW = TRUE;
		tr->optimizeAllTrees = TRUE;
		break;
	      case 'x':
		adef->mode = DISTANCE_MODE;
		adef->computeDistance = TRUE;
		break;
	      case 'y':
		adef->mode = CLASSIFY_MP;
		break;
	      default:
		if(flagCheck)
		  {
		    printf("Option -f %s not supported\n", optarg);
		    fprintf(flagCheckFile, "Option -f %s not supported\n", optarg);
		    invalidOptions++;
		  }
		else
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
	    {
	      int 
		result;
	      strcpy(model,optarg);
	      
	      result = modelExists(model, adef);

	      if(adef->proteinMatrix == AUTO && !(adef->optimizeBaseFrequencies) && adef->protEmpiricalFreqs)
		{
		  printf("\nError: Option AUTOF has been deprecated, exiting\n\n");
		   errorExit(-1);
		}
	       if(adef->proteinMatrix == AUTO && adef->optimizeBaseFrequencies)
		{
		  printf("\nError: Option AUTOX has been deprecated, exiting\n\n");
		   errorExit(-1);
		}
	      
	      

	      if(result == 0)
		{		  		  
		  if(processID == 0)
		    {		    
		      printf("Model %s does not exist\n\n", model);
		      printf("For BINARY data use:  BINCAT[X]             or BINGAMMA[X]             or\n");
		      printf("                      BINCATI[X]            or BINGAMMAI[X]            or\n");
		      printf("                      ASC_BINGAMMA[X]       or ASC_BINCAT[X]\n");
		      printf("For DNA data use:     GTRCAT[X]             or GTRGAMMA[X]             or\n");
		      printf("                      GTRCATI[X]            or GTRGAMMAI[X]            or\n");
		      printf("                      ASC_GTRGAMMA[X]       or ASC_GTRCAT[X]\n");
		      printf("For Multi-state data: MULTICAT[X]          or MULTIGAMMA[X]           or\n");
		      printf("                      MULTICATI[X]         or MULTIGAMMAI[X]          or\n");		
		      printf("                      ASC_MULTIGAMMA[X]    or ASC_MULTICAT[X]\n");
		      printf("For AA data use:      PROTCATmatrixName[F|X]  or PROTGAMMAmatrixName[F|X]  or\n");
		      printf("                      PROTCATImatrixName[F|X] or PROTGAMMAImatrixName[F|X] or\n");
		      printf("                      ASC_PROTGAMMAmatrixName[X] or  ASC_PROTCATmatrixName[X]\n");
		      printf("The AA substitution matrix can be one of the following: \n");
		      
		      {
			int 
			  i;
			
			for(i = 0; i < NUM_PROT_MODELS - 1; i++)
			  {
			    if(i % 8 == 0)	  
			      printf("\n");
			    printf("%s, ", protModels[i]);
			  }
			
			printf("%s\n\n", protModels[i]);
		      }
		      
		      printf("With the optional \"F\" appendix you can specify if you want to use empirical base frequencies.\n");
		      printf("With the optional \"X\" appendix you can specify that you want to do a ML estimate of base frequencies.\n");
		      printf("Please note that for mixed models you can in addition specify the per-gene model in\n");
		      printf("the mixed model file (see manual for details).\n");
		    }
		  errorExit(-1);
		}
	      else
		modelSet = 1;
	    }
	    break;
	  default:
	    if(flagCheck)
	      {
		printf("Option %s not supported\n", argv[optind-1]);
		fprintf(flagCheckFile, "Option %s not supported\n", argv[optind-1]);
		invalidOptions++;		
	      }
	    else
	      errorExit(-1);
	  }
    }
  
  if(flagCheck)
    {
      if(invalidOptions == 0)
	{
	  printf("All options supported\n");
	  fprintf(flagCheckFile, "All options supported\n");
	  fclose(flagCheckFile);
	  exit(0);
	}
      else
	{
	  fclose(flagCheckFile);
	  exit(-1);
	}
     
    }
  
  if(tr->useAccumulatedEPACutoff && epaSet)
    {
      printf("\nError: the EPA flag: --epa-accumulated-threshold can neither be used in combination with --epa-keep-placements nor --epa-prob-threshold\n\n");
      errorExit(-1);
    }
    
  if(disablePatternCompression)
    adef->compressPatterns = FALSE;

  if(adef->mode == SUBTREE_EPA)
    {
      if(!treeSet)
	{
	  printf("Error: for subtree placements you need to specify a tree file via \"-t\"\n");
	  errorExit(-1);
	}
      if(!treesSet)
	{
	  printf("Error: for subtree placements you need to specify a subtree specification file via \"-z\"\n");
	  errorExit(-1);
	}
    }
  
  if(adef->mode == STEAL_BRANCH_LENGTHS)
    {
      if(!treeSet)
	{
	  printf("Error: for branch length stealing you need to specify a tree file via \"-t\"\n");
	  errorExit(-1);
	}
      if(!adef->useMultipleModel)
	{
	  printf("Error: for branch length stealing  you need to specify a partition file via \"-q\"\n");
	  errorExit(-1);
	}
      if(!adef->perGeneBranchLengths)
	{
	  printf("Error: for branch length stealing you need to specify a per partition branch length estimate via \"-M\"\n");
	  errorExit(-1);
	}      
    }

#ifdef _USE_PTHREADS
  if(NumberOfThreads < 2)
    {
      printf("\nWARNING: The number of threads is currently set to %d\n", NumberOfThreads);
      printf("You can specify the number of threads to run via -T numberOfThreads\n");
      printf("NumberOfThreads must be set to an integer value greater than 1\n\n");
      printf("RAxML, will now set the number of threads automatically to 2 !\n\n");
      NumberOfThreads = 2;
      //errorExit(-1);
    }
#endif

#ifdef _QUARTET_MPI
  if(adef->mode != QUARTET_CALCULATION)
    {
      if(processID == 0)
	{
	  printf("you are using the dedicated RAxML MPI version for parallel quartet computations\n");
	  printf("However you are not using the quartet option \"-f q\", raxml will exit now ...\n");
	}
      
      errorExit(-1);
    }

  if(!adef->useBinaryModelFile)
    {
       if(processID == 0)
	{
	  printf("you are using the dedicated RAxML MPI version for parallel quartet computations\n");
	  printf("However you must provide a binary model file via \"-R\" when using the MPI version, raxml will exit now ...\n");
	}
      
      errorExit(-1);
    }
#endif

  if(adef->mode ==  ANCESTRAL_SEQUENCE_TEST && !yFileSet)
    {
      if(!yFileSet)
	{
	  printf("Error, for using the ancestral sequence test you have to provide a ancestral taxon name\n");
	  printf("candidate file via \"-Y\" \n");
	  errorExit(-1);
	}

      if(!treeSet)
	{
	  printf("Error, for using the ancestral sequence test you have to provide a tree file\n");
	  printf("via \"-t\" \n");
	  errorExit(-1);
	}
    }

  if(tr->catOnly && adef->rapidBoot)
    {
      printf("Error, you can not use \"-F\" in conjunction with the rapid bootstrapping option!\n");
      printf("it will only work with standard ML tree searches\n");
      errorExit(-1);
    }

  if(tr->catOnly && adef->boot)
    {
      printf("Error, you can not use \"-F\" in conjunction with the standard bootstrapping option!\n");
      printf("it will only work with standard ML tree searches\n");
      errorExit(-1);
    }
     

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


 
  if(adef->rellBootstrap && adef->parsimonySeed == 0)
    {
      if(processID == 0)
	{
	  printf("Error, you must specify a random number seed via \"-p\" in conjunction with the RELL bootstrap\n");
	  errorExit(-1);
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

  if(adef->grouping)
    {
      if(processID == 0)
	{
	  if(adef->parsimonySeed == 0)
	    {
	      printf("\nERROR: you must specify a random number seed via \"-p\" when using multi-furcating constraint trees\n");
	      errorExit(-1);
	    }
	}
      
      adef->constraintSeed = adef->parsimonySeed;
      assert(adef->constraintSeed > 0);
    }

  if(tr->searchConvergenceCriterion && ((adef->rapidBoot > 0) || (adef->allInOne)))
    {
      if(processID == 0)
	{
	  printf("\nError: the tree search convergence criterion \"-D\" has no effect in conjunction with the: \n");
	  printf("\"-x\" or \"-f a\" options.\n\n");
	  errorExit(-1);
	}
    }

  if(adef->mode == SH_LIKE_SUPPORTS)
    {
      if(processID == 0)
	{
	  if(adef->parsimonySeed == 0)
	    {
	      printf("\nERROR: you must specify a random number seed via \"-p\" to calculate SH-like supports\n");
	      errorExit(-1);
	    }
	}
         
      assert(adef->parsimonySeed > 0);
    }

  if(adef->mode == BOOTSTOP_ONLY ||  adef->bootStopping == TRUE)
    {
      if(processID == 0)
	{
	  if(adef->parsimonySeed == 0)
	    {
	      printf("\nERROR: you must specify a random number seed via \"-p\" for bootstopping\n");
	      errorExit(-1);
	    }
	}
         
      assert(adef->parsimonySeed > 0);
    }

   if(adef->randomStartingTree)
    {
      if(processID == 0)
	{
	  if(adef->parsimonySeed == 0)
	    {
	      printf("\nERROR: you must specify a random number seed via \"-p\" when using random starting trees\n");
	      errorExit(-1);
	    }
	}
         
      assert(adef->parsimonySeed > 0);
    }


  if(adef->mode == FAST_SEARCH && (adef->grouping || adef->constraint))
    {
      if(processID == 0)
	printf("\n ERROR: Fast ML search algorithms -f F and -f E can not take as input constraint trees specified via -g or -r, since they will be ignored\n");
      errorExit(-1);
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

  if(adef->mode == ROOT_TREE && !treeSet)
    {
      if(processID == 0)
	printf("\n  Error, for the tree rooting algorithm you need to specify a file containing the tree you want to root via \"-t\"\n");
      errorExit(-1);
    }

  if((adef->mode == CALC_BIPARTITIONS || adef->mode == CALC_BIPARTITIONS_IC) && !treesSet)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition and IC computation mode you must specify a file containing multiple trees with the \"-z\" option\n");
      errorExit(-1);
    }

  if((adef->mode == CALC_BIPARTITIONS || adef->mode == CALC_BIPARTITIONS_IC) && !adef->restart)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition and IC computation mode you must specify a tree on which bipartition information will be drawn with the \"-t\" option\n");
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

  if(adef->mode == CALC_BIPARTITIONS_IC)
    {
      if(processID == 0)
	 printf("Computation of IC and TC scores (-f i) not implemented for the MPI-Version\n");
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
	
	printf("Warning, you specified a working directory via \"-w\"\n");
	printf("Keep in mind that RAxML only accepts absolute path names, not relative ones!\n");

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

#if (defined(_WAYNE_MPI) || defined (_QUARTET_MPI))
  MPI_Finalize();
#endif

  exit(e);

}



static void makeFileNames(void)
{
  int infoFileExists = 0;

  strcpy(verboseSplitsFileName, workdir);
  strcpy(permFileName,         workdir);
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);
  strcpy(checkpointFileName,   workdir);
  strcpy(infoFileName,         workdir);
  strcpy(randomFileName,       workdir);
  strcpy(bootstrapFileName,    workdir);
  strcpy(bipartitionsFileName, workdir);
  strcpy(bipartitionsFileNameBranchLabels, workdir);
  strcpy(icFileNameBranchLabels, workdir);
  strcpy(icFileNameBranchLabelsUniform, workdir);
  strcpy(icFileNameBranchLabelsStochastic, workdir);
  strcpy(ratesFileName,        workdir);
  strcpy(lengthFileName,       workdir);
  strcpy(lengthFileNameModel,  workdir);
  strcpy(perSiteLLsFileName,  workdir);
  strcpy(binaryModelParamsOutputFileName,  workdir);
  strcpy(rellBootstrapFileName, workdir);
  strcpy(mesquiteModel, workdir);
  strcpy(mesquiteTrees, workdir);
  strcpy(mesquiteMLTrees, workdir);
  strcpy(mesquiteMLLikes, workdir);
  
  strcat(verboseSplitsFileName,             "RAxML_verboseSplits.");
  strcat(permFileName,                      "RAxML_parsimonyTree.");
  strcat(resultFileName,                    "RAxML_result.");
  strcat(logFileName,                       "RAxML_log.");
  strcat(checkpointFileName,                "RAxML_checkpoint.");
  strcat(infoFileName,                      "RAxML_info.");
  strcat(randomFileName,                    "RAxML_randomTree.");
  strcat(bootstrapFileName,                 "RAxML_bootstrap.");
  strcat(bipartitionsFileName,              "RAxML_bipartitions.");
  strcat(bipartitionsFileNameBranchLabels,  "RAxML_bipartitionsBranchLabels.");
  strcat(icFileNameBranchLabels,            "RAxML_IC_Score_BranchLabels.");
  strcat(icFileNameBranchLabelsStochastic,  "RAxML_Corrected_Probabilistic_IC_Score_BranchLabels.");
  strcat(icFileNameBranchLabelsUniform,     "RAxML_Corrected_Lossless_IC_Score_BranchLabels.");
  strcat(ratesFileName,                     "RAxML_perSiteRates.");
  strcat(lengthFileName,                    "RAxML_treeLength.");
  strcat(lengthFileNameModel,               "RAxML_treeLengthModel.");
  strcat(perSiteLLsFileName,                "RAxML_perSiteLLs.");
  strcat(binaryModelParamsOutputFileName,   "RAxML_binaryModelParameters.");
  strcat(rellBootstrapFileName,             "RAxML_rellBootstrap.");
  strcat(mesquiteModel,                     "RAxML_mesquiteModel.");
  strcat(mesquiteTrees,                     "RAxML_mesquiteTrees.");
  strcat(mesquiteMLTrees,                   "RAxML_mesquite_ML_Trees.");
  strcat(mesquiteMLLikes,                   "RAxML_mesquite_ML_Likes.");

  strcat(verboseSplitsFileName,            run_id);
  strcat(permFileName,                     run_id);
  strcat(resultFileName,                   run_id);
  strcat(logFileName,                      run_id);
  strcat(checkpointFileName,               run_id);
  strcat(infoFileName,                     run_id);
  strcat(randomFileName,                   run_id);
  strcat(bootstrapFileName,                run_id);
  strcat(bipartitionsFileName,             run_id);
  strcat(bipartitionsFileNameBranchLabels, run_id);  
  strcat(icFileNameBranchLabels,           run_id); 
  strcat(icFileNameBranchLabelsUniform,    run_id);
  strcat(icFileNameBranchLabelsStochastic, run_id);
  strcat(ratesFileName,                    run_id);
  strcat(lengthFileName,                   run_id);
  strcat(lengthFileNameModel,              run_id);
  strcat(perSiteLLsFileName,               run_id);
  strcat(binaryModelParamsOutputFileName,  run_id);
  strcat(rellBootstrapFileName,            run_id);
  strcat(mesquiteModel,                    run_id);
  strcat(mesquiteTrees,                    run_id);
  strcat(mesquiteMLTrees,                  run_id);
  strcat(mesquiteMLLikes,                  run_id);

#ifdef _WAYNE_MPI  
  {
    char 
      buf[64];
    
    sprintf(buf, "%d", processID);

    strcpy(bootstrapFileNamePID, bootstrapFileName);
    strcat(bootstrapFileNamePID, ".PID.");   
    strcat(bootstrapFileNamePID, buf);

    strcpy(rellBootstrapFileNamePID, rellBootstrapFileName);
    strcat(rellBootstrapFileNamePID, ".PID.");   
    strcat(rellBootstrapFileNamePID, buf);

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

	  if(tr->partitionData[model].optimizeBaseFrequencies)
	    printBothOpen("Initial base frequencies, prior to ML estimate: ");
	  else
	    printBothOpen("Base frequencies: ");
	  
	  if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
	    {
	      int
		k;
	      
	      printBothOpen("\n");
	      
	      for(k = 0; k < 4; k++)
		{
		  printBothOpen("LG4 %d: ", k);
		  for(i = 0; i < tr->partitionData[model].states; i++)
		    printBothOpen("%1.3f ", tr->partitionData[model].frequencies_LG4[k][i]);
		  printBothOpen("\n");
		}
	    }
	  else
	    {
	      for(i = 0; i < tr->partitionData[model].states; i++)
		printBothOpen("%1.3f ", tr->partitionData[model].frequencies[i]);
	    }
	  
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
     
      printVersionInfo(FALSE, infoFile);

      
      
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
	case STEAL_BRANCH_LENGTHS:
	  printBoth(infoFile, "\nRAxML branch length stealing\n\n");
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
	case CALC_BIPARTITIONS_IC:
	  printBoth(infoFile, "\nRAxML IC and TC score Computation: Computing IC and TC scores induced by trees in file %s w.r.t. tree in file %s\n\n",
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
	case FAST_SEARCH:
	  printBoth(infoFile, "\nRAxML experimental very fast tree search\n\n");
	  break;
	case SH_LIKE_SUPPORTS:
	  printBoth(infoFile, "\nRAxML computation of SH-like support values on a given tree\n\n");
	  break;
	case EPA_SITE_SPECIFIC_BIAS:
	  printBoth(infoFile, "\nRAxML experimental site-specfific phylogenetic placement bias analysis algorithm\n\n");
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
	case ANCESTRAL_SEQUENCE_TEST:
	  printBoth(infoFile, "\nRAxML ancestral sequence test for Jiajie\n\n");
	  break;
	case PLAUSIBILITY_CHECKER:
	  printBoth(infoFile, "\nRAxML large-tree plausibility-checker\n\n");
	  break;
	case ROOT_TREE:
	  printBoth(infoFile, "\nRAxML tree rooting algorithm\n\n");
	  break;
	case SUBTREE_EPA:
	  printBoth(infoFile, "\nRAxML Evolutionary Placement Algorithm for (taxonomic) subtrees\n\n");
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
		    strcpy(treeType, "user-specified");
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
	    printBoth(infoFile, "%s model of rate heterogeneity, ML estimate of alpha-parameter\n\n", modelType);
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
		  {
		    char 
		      *matrices[4] = {"GTR", "JC69", "K80", "HKY85"};
		    
		    int 
		      index = -1;

		    printBoth(infoFile, "DataType: DNA\n");
		    
		    if(tr->useJC69)
		      index = 1;
		    else
		      {
			if(tr->useK80)
			  index = 2;
			else
			  {
			    if(tr->useHKY85)
			      index = 3;
			    else
			      index = 0;
			  }
		      }
		    
		    printBoth(infoFile, "Substitution Matrix: %s\n", matrices[index]);
		    
		    if(tr->partitionData[model].optimizeBaseFrequencies)
		      printBoth(infoFile, "Base frequencies: ML estimate\n");
		  }
		  break;
		case AA_DATA:
		  assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < NUM_PROT_MODELS);
		  printBoth(infoFile, "DataType: AA\n");
		  if(tr->partitionData[model].protModels != PROT_FILE)
		    {		     		     
		      printBoth(infoFile, "Substitution Matrix: %s\n", protModels[tr->partitionData[model].protModels]);
		      if(!tr->partitionData[model].optimizeBaseFrequencies)
			printBoth(infoFile, "Using %s base frequencies\n", (tr->partitionData[model].usePredefinedProtFreqs == TRUE)?"fixed":"empirical");		      
		      else
			printBoth(infoFile, "Using ML estimate of base frequencies\n");
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
		  if(tr->partitionData[model].optimizeBaseFrequencies)
		    printBoth(infoFile, "Base frequencies: ML estimate\n");
		  break;
		case SECONDARY_DATA:
		  printBoth(infoFile, "DataType: SECONDARY STRUCTURE\n");		  
		  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]); 
		  if(tr->partitionData[model].optimizeBaseFrequencies)
		    printBoth(infoFile, "Base frequencies: ML estimate\n");
		  break;
		case SECONDARY_DATA_6:
		  printBoth(infoFile, "DataType: SECONDARY STRUCTURE 6 STATE\n");		  
		  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]);
		  if(tr->partitionData[model].optimizeBaseFrequencies)
		    printBoth(infoFile, "Base frequencies: ML estimate\n");
		  break;
		case SECONDARY_DATA_7:
		  printBoth(infoFile, "DataType: SECONDARY STRUCTURE 7 STATE\n");		 
		  printBoth(infoFile, "Substitution Matrix: %s\n", secondaryModelList[tr->secondaryStructureModel]); 
		  if(tr->partitionData[model].optimizeBaseFrequencies)
		    printBoth(infoFile, "Base frequencies: ML estimate\n");
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
		  if(tr->partitionData[model].optimizeBaseFrequencies)
		     printBoth(infoFile, "Base frequencies: ML estimate\n");
		  break;
		case GENERIC_64:
		  printBoth(infoFile, "DataType: Codon\n"); 
		  if(tr->partitionData[model].optimizeBaseFrequencies)
		    printBoth(infoFile, "Base frequencies: ML estimate\n");		  
		  break;		
		default:
		  assert(0);
		}

	      if(tr->partitionData[model].ascBias)
		printBoth(infoFile, "Correcting likelihood for ascertainment bias\n");

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
    case MORPH_CALIBRATOR:
      break;
    case TREE_EVALUATION:


      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);

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
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);

		  logFile = myfopen(temporaryFileName, "wb");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);

		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "wb");
		  
		  break;
		case CAT:
		  if(adef->mesquite)
		    Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
		  else
		    Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE, FALSE, FALSE);

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
	      if(adef->mesquite)
		Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
	      else
		Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE, FALSE, FALSE);
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
	  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);	  

	  logFile = myfopen(fileName, "ab");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	  
	  if(adef->perGeneBranchLengths)
	    printTreePerGene(tr, adef, fileName, "ab");
	}
      else
	{
	  if(adef->mesquite)
	    Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
	  else
	    Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE, FALSE, FALSE);
	  
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



void printBipartitionResult(tree *tr, analdef *adef, boolean finalPrint, boolean printIC, char *fileName)
{
  if(processID == 0 || adef->allInOne)
    {
      FILE 
	*logFile;
     
      if(!printIC)
	{
	  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, TRUE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE, printIC, FALSE);
            
	  logFile = myfopen(bipartitionsFileName, "ab");
      
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);
	}
      
     
      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, TRUE, FALSE, printIC, FALSE);
      
      
      logFile = myfopen(fileName, "ab");     
      
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


	  if(!adef->checkpoints && !adef->mesquite)
	    {
	      logFile = myfopen(temporaryFileName, "ab");

	      fprintf(logFile, "%f %f\n", t, lh);

	      fclose(logFile);
	    }
	  else
	    {
	      if(adef->mesquite)
		{ 
		  char 
		    temporaryFileName2[1024] = "";
		  
		  logFile = myfopen(temporaryFileName, "ab");

		  fprintf(logFile, "%f %f\n", t, lh);
		  
		  fclose(logFile);

		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
		  
		  //now print raxml result 
		  
		  strcpy(temporaryFileName2, resultFileName);
		  
		  if(adef->multipleRuns > 1)
		    {
		      char
			treeID[64] = "";
		      sprintf(treeID, "%d", tr->treeID);
		      strcat(temporaryFileName2, ".RUN.");
		      strcat(temporaryFileName2, treeID);
		    }

		  //printf("Mesquite printing intemediate tree %s to file %s\n", tr->tree_string, temporaryFileName2);

		  logFile = myfopen(temporaryFileName2, "wb");
		  
		  fprintf(logFile, "%s", tr->tree_string);
		  
		  fclose(logFile);
		}
	      else
		{
		  assert(adef->checkpoints);

		  logFile = myfopen(temporaryFileName, "ab");
		  
		  fprintf(logFile, "%f %f %d\n", t, lh, tr->checkPointCounter);
		  
		  fclose(logFile);
		  
		  strcat(checkPoints, ".");
		  
		  sprintf(treeID, "%d", tr->checkPointCounter);
		  strcat(checkPoints, treeID);
		  
		  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE, FALSE, FALSE);
		  
		  logFile = myfopen(checkPoints, "ab");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		  
		  tr->checkPointCounter++;
		}
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

      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES, FALSE, FALSE, FALSE, FALSE);

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
		FILE 
		  *infoFile = myfopen(infoFileName, "ab");		

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

		    if(tr->partitionData[model].optimizeBaseFrequencies)
		      {
			int
			  k,
			  states = tr->partitionData[model].states;
			
			fprintf(infoFile, "ML estimate base freqs[%d]: ", model);
			
			for(k = 0; k < states; k++)
			  fprintf(infoFile, "%f ", tr->partitionData[model].frequencies[k]);
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

      if(tr->useBrLenScaler)
	printBothOpen("Branch length scaler: %f\n", tr->partitionData[model].brLenScaler);

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
	    char *freqNames[20] = {"A", "R", "N","D", "C", "Q", "E", "G",
				   "H", "I", "L", "K", "M", "F", "P", "S",
				   "T", "W", "Y", "V"};

	    if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
	      {
		int 
		  k;
		
		printBothOpen("\nLG4X rates: ");
		for(k = 0; k < 4; k++)
		  printBothOpen("%f ", tr->partitionData[model].gammaRates[k]);		 

		printBothOpen("\n\nLG4X weights: ");
		for(k = 0; k < 4; k++)
		  printBothOpen("%f ", tr->partitionData[model].weights[k]);	

		printBothOpen("\n\n");

		for(k = 0; k < 4; k++)
		  {
		    printBothOpen("LGM %d\n", k);
		    printRatesRest(20, tr->partitionData[model].substRates_LG4[k], freqNames);
		    printBothOpen("\n");
		    printFreqs(20, tr->partitionData[model].frequencies_LG4[k], freqNames);
		  }
	      }
	    else
	      {
		printRatesRest(20, r, freqNames);
		printBothOpen("\n");
		printFreqs(20, f, freqNames);
	      }
	  }
	  break;
	case GENERIC_32:
	  {
	    char *freqNames[32] = {"0", "1", "2", "3", "4", "5", "6", "7", 
				   "8", "9", "A", "B", "C", "D", "E", "F",
				   "G", "H", "I", "J", "K", "L", "M", "N",
				   "O", "P", "Q", "R", "S", "T", "U", "V"}; 

	    printRatesRest(tr->partitionData[model].states, r, freqNames);
	    printBothOpen("\n");
	    printFreqs(tr->partitionData[model].states, f, freqNames);
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
	case TREE_EVALUATION :	
	case OPTIMIZE_BR_LEN_SCALER:
	  
	  if(adef->mode == OPTIMIZE_BR_LEN_SCALER)
	    printBothOpen("\n\nOverall Time for Tree Evaluation with branch length scalers: %f\n", t);
	  else	    
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
		    
		    if(!tr->partitionData[model].usePredefinedProtFreqs || tr->partitionData[model].optimizeBaseFrequencies)
		      params += 19;
		    if(tr->partitionData[model].protModels == LG4X)
		      params += 6;
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
		      
		      assert(!tr->partitionData[model].optimizeBaseFrequencies);

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
		
		if(adef->useInvariant)
		  params += 2;
		else /* GAMMA */
		  params += 1;
	      }
	    
	    if(linkedProteinGTR)
	      params += 189;

	    if(adef->mode == TREE_EVALUATION)
	      {
		if(tr->multiBranch)
		  paramsBrLen = params + tr->NumberOfModels * (2 * tr->mxtips - 3);
		else
		  paramsBrLen = params + 2 * tr->mxtips - 3;
	      }
	    else
	      {			
		paramsBrLen = params + tr->NumberOfModels;		
	      }
	    
	    printBothOpen("\n");

	   
	    printBothOpen("Number of free parameters for AIC-TEST(BR-LEN): %d\n",    paramsBrLen);
	    printBothOpen("Number of free parameters for AIC-TEST(NO-BR-LEN): %d\n", params);
	    
	    
	    printBothOpen("\n\n");
	    
	    printModelParams(tr, adef);
	    
	    if(adef->mode == TREE_EVALUATION)
	      {
		printBothOpen("Final tree written to:                 %s\n", resultFileName);
		printBothOpen("Execution Log File written to:         %s\n", logFileName);
	      }
	 
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
		  double 
		    avgLH = 0.0,
		    bestLH = unlikely;
		  
		  int 
		    i, 
		    bestI  = 0;

		  for(i = 0; i < adef->multipleRuns; i++)
		    {
		      avgLH += tr->likelihoods[i];
		      
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
		  printBothOpen("Execution information file written to: %s\n", infoFileName);
		}

	      if(adef->mesquite)
		{
		  printBothOpen("Mesquite tree file written to:   %s\n", mesquiteTrees);
		  printBothOpen("Mesquite model file written to:   %s\n\n", mesquiteModel);
		}
	    }
	  
	  break;
	case CALC_BIPARTITIONS:
	  printBothOpen("\n\nTime for Computation of Bipartitions %f\n", t);
	  printBothOpen("Tree with bipartitions written to file:  %s\n", bipartitionsFileName);
	  printBothOpen("Tree with bipartitions as branch labels written to file:  %s\n", bipartitionsFileNameBranchLabels);	  
	  printBothOpen("Execution information file written to :  %s\n",infoFileName);
	  break;
	case CALC_BIPARTITIONS_IC:
	  printBothOpen("\n\nTime for Computation of TC and IC scores %f\n\n", t);
	  
	  if(tr->corrected_IC_Score)
	    {
	      printBothOpen("Tree with corrected (for partial gene trees) stochastic IC scores as branch labels written to file:  %s\n\n", icFileNameBranchLabelsStochastic);
	      printBothOpen("Tree with corrected (for partial gene trees) uniform IC scores as branch labels written to file:  %s\n\n", icFileNameBranchLabelsUniform);
	      
	    }
	  else
	    printBothOpen("Tree with IC scores as branch labels written to file:  %s\n", icFileNameBranchLabels);	  
	  
	  printBothOpen("Execution information file written to :  %s\n\n",infoFileName);
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
	case ROOT_TREE:
	  printBothOpen("\n\nTime for tree rooting: %f\n\n", t);
	  break;
	default:
	  assert(0);
	}
    }

}


/************************************************************************************/
static void setupPresenceMask(tree *tr)
{
  int 
    model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      int 
	j;

      for(j = 1; j <= tr->mxtips; j++)
	{
	  unsigned int 
	    presenceMask = 0,	  
	    i;

	  for(i = 0; i < tr->partitionData[model].width; i++)	      	   	 	             	  		  	    	 		
	    presenceMask = presenceMask | mask32[tr->partitionData[model].yVector[j][i]];       
	  
	  tr->partitionData[model].presenceMap[j] = presenceMask;
	  /*
	    #ifdef _USE_PTHREADS
	    printf("Thread %d Taxon %d has a total of %d states present\n", tr->threadID, j, BIT_COUNT(presenceMask));
	    #else
	    printf("Taxon %d has a total of %d states present\n", j, BIT_COUNT(presenceMask));
	    #endif
	  */
	}
    }  
}

/**************************************************************************************************/

#ifdef _USE_PTHREADS






static void computeFraction(tree *localTree, int tid, int n)
{
  int  
    model;

  size_t
    i;

  for(model = 0; model < localTree->NumberOfModels; model++)
    {
      int width = 0;

      for(i = localTree->partitionData[model].lower; i < localTree->partitionData[model].upper; i++)
	if(i % (size_t)n == (size_t)tid)
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
    myLength = 0;

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
      localTree->innerNodes              = tr->innerNodes;
      localTree->useFastScaling          = tr->useFastScaling;
      localTree->perPartitionEPA         = tr->perPartitionEPA;
      localTree->maxCategories           = tr->maxCategories;
     
      localTree->originalCrunchedLength  = tr->originalCrunchedLength;
      localTree->NumberOfModels          = tr->NumberOfModels;
      localTree->mxtips                  = tr->mxtips;
      localTree->multiBranch             = tr->multiBranch;
         
      localTree->nameList                = tr->nameList;
      localTree->numBranches             = tr->numBranches;
      localTree->lhs                     = tr->lhs;
      localTree->executeModel            = (boolean*)rax_malloc(sizeof(boolean) * localTree->NumberOfModels);
      localTree->perPartitionLH          = (double*)rax_malloc(sizeof(double)   * localTree->NumberOfModels);
      localTree->storedPerPartitionLH    = (double*)rax_malloc(sizeof(double)   * localTree->NumberOfModels);

     

      localTree->partitionContributions = (double*)rax_malloc(sizeof(double)   * localTree->NumberOfModels);

      localTree->partitionData = (pInfo*)rax_malloc(sizeof(pInfo) * localTree->NumberOfModels);

      /* not required any more */
      //localTree->td[0].count = 0;
      //localTree->td[0].ti    = (traversalInfo *)rax_malloc(sizeof(traversalInfo) * localTree->mxtips);

      localTree->cdta               = (cruncheddata*)rax_malloc(sizeof(cruncheddata));
      localTree->cdta->patrat       = tr->cdta->patrat;
      localTree->cdta->patratStored = tr->cdta->patratStored;

      localTree->discreteRateCategories = tr->discreteRateCategories;     
      localTree->ascertainmentCorrectionType = tr->ascertainmentCorrectionType;

      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  localTree->partitionData[model].numberOfCategories    = tr->partitionData[model].numberOfCategories;
	  localTree->partitionData[model].states     = tr->partitionData[model].states;
	  localTree->partitionData[model].maxTipStates    = tr->partitionData[model].maxTipStates;
	  localTree->partitionData[model].dataType   = tr->partitionData[model].dataType;
	  localTree->partitionData[model].protModels = tr->partitionData[model].protModels;
	  localTree->partitionData[model].usePredefinedProtFreqs  = tr->partitionData[model].usePredefinedProtFreqs;
	  localTree->partitionData[model].optimizeBaseFrequencies  = tr->partitionData[model].optimizeBaseFrequencies;
	  localTree->partitionData[model].ascBias  = tr->partitionData[model].ascBias;
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
	width = tr->partitionData[model].width,
	i;

      myLength += width;

      memoryRequirements += (size_t)(tr->discreteRateCategories) * (size_t)(tr->partitionData[model].states) * width;
     
      tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;
      
      tr->partitionData[model].gapVector = (unsigned int*)rax_calloc(tr->partitionData[model].gapVectorLength * 2 * tr->mxtips, sizeof(unsigned int));
      
      tr->partitionData[model].initialGapVectorSize = tr->partitionData[model].gapVectorLength * 2 * tr->mxtips * sizeof(int);
      
      /* always multiply by 4 due to frequent switching between CAT and GAMMA in standard RAxML */

      tr->partitionData[model].gapColumn = (double *)rax_malloc(
								    ((size_t)(tr->innerNodes)) *
								    ((size_t)(4)) * 
								    ((size_t)(tr->partitionData[model].states)) *
								    sizeof(double));

       //asc
          
      if(tr->partitionData[model].ascBias && tid == 0)
	{	 
	  tr->partitionData[model].ascOffset = 4 * tr->partitionData[model].states * tr->partitionData[model].states;
	  
	  tr->partitionData[model].ascVector = (double *)rax_malloc(((size_t)tr->innerNodes) *
								    ((size_t)tr->partitionData[model].ascOffset) * 
								    sizeof(double));

	  tr->partitionData[model].ascExpVector = (int *)rax_calloc(((size_t)tr->innerNodes) * ((size_t)tr->partitionData[model].states),
								 sizeof(int));	  
	  
	  tr->partitionData[model].ascSumBuffer = (double *)rax_malloc(((size_t)tr->partitionData[model].ascOffset) *
								       sizeof(double));	  
	}
      
      //asc
		             
      for(i = 0; i < tr->innerNodes; i++)
	{
	  tr->partitionData[model].xVector[i]   = (double*)NULL;     
	  tr->partitionData[model].expVector[i]   = (int*)NULL;
	}
    }

  if(tid == 0)
    {
      tr->perSiteLL       = (double *)rax_malloc((size_t)tr->cdta->endsite * sizeof(double));
      assert(tr->perSiteLL != NULL);
    }
  
  tr->sumBuffer  = (double *)rax_malloc(memoryRequirements * sizeof(double));
  assert(tr->sumBuffer != NULL);
   
  tr->y_ptr = (unsigned char *)rax_malloc(myLength * (size_t)(tr->mxtips) * sizeof(unsigned char));
  assert(tr->y_ptr != NULL);  

  tr->perSiteLLPtr     = (double*) rax_malloc(myLength * sizeof(double));
  assert(tr->perSiteLLPtr != NULL);

  tr->wgtPtr           = (int*)    rax_malloc(myLength * sizeof(int));
  assert(tr->wgtPtr != NULL);  

  tr->invariantPtr     = (int*)    rax_malloc(myLength * sizeof(int));
  assert(tr->invariantPtr != NULL);

  tr->rateCategoryPtr  = (int*)    rax_malloc(myLength * sizeof(int));
  assert(tr->rateCategoryPtr != NULL);
}






inline static void sendTraversalInfo(tree *localTree, tree *tr)
{   
  localTree->td[0] = tr->td[0];  
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

static void copyLG4(tree *localTree, tree *tr, int model, const partitionLengths *pl)
{
  if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
    {
      int 
	k;
          
      for(k = 0; k < 4; k++)
	{
	   memcpy(localTree->partitionData[model].EIGN_LG4[k],        tr->partitionData[model].EIGN_LG4[k],        pl->eignLength * sizeof(double));
	   memcpy(localTree->partitionData[model].rawEIGN_LG4[k],        tr->partitionData[model].rawEIGN_LG4[k],        pl->eignLength * sizeof(double));
	   memcpy(localTree->partitionData[model].EV_LG4[k],          tr->partitionData[model].EV_LG4[k],          pl->evLength * sizeof(double));
	   memcpy(localTree->partitionData[model].EI_LG4[k],          tr->partitionData[model].EI_LG4[k],          pl->eiLength * sizeof(double));
	   memcpy(localTree->partitionData[model].substRates_LG4[k],  tr->partitionData[model].substRates_LG4[k],  pl->substRatesLength * sizeof(double));
	   memcpy(localTree->partitionData[model].frequencies_LG4[k], tr->partitionData[model].frequencies_LG4[k], pl->frequenciesLength * sizeof(double));
	   memcpy(localTree->partitionData[model].tipVector_LG4[k],   tr->partitionData[model].tipVector_LG4[k],   pl->tipVectorLength * sizeof(double));
	}
    }
}

static void execFunction(tree *tr, tree *localTree, int tid, int n)
{
  double volatile result;

  size_t
    i;

  int 
    currentJob,
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
      setupPresenceMask(localTree);
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
	    for(i = 0; i < (size_t)localTree->NumberOfModels; i++)
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

	if(tid > 0)
	  {
	    memcpy(localTree->coreLZ,   tr->coreLZ,   sizeof(double) *  localTree->numBranches);
	    memcpy(localTree->executeModel, tr->executeModel, sizeof(boolean) * localTree->NumberOfModels);
	  }
	
	execCore(localTree, dlnLdlz, d2lnLdlz2);

	if(!tr->multiBranch)
	  {
	    reductionBuffer[tid]    = dlnLdlz[0];
	    reductionBufferTwo[tid] = d2lnLdlz2[0];
	  }
	else
	  {
	    for(i = 0; i < (size_t)localTree->NumberOfModels; i++)
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
	      
	      copyLG4(localTree, tr, model, pl);
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
	      
	      copyLG4(localTree, tr, model, pl);	     
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
	    {
	      localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	    }
	}
      break;
    case THREAD_OPT_INVAR:
      if(tid > 0)
	{
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
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
	      
	      copyLG4(localTree, tr, model, pl);

	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	      localTree->partitionData[model].brLenScaler = tr->partitionData[model].brLenScaler;
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
	      
	      copyLG4(localTree, tr, model, pl);

	      memcpy(localTree->partitionData[model].weights, tr->partitionData[model].weights, sizeof(double) * 4);
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	      localTree->partitionData[model].alpha = tr->partitionData[model].alpha;
	      localTree->partitionData[model].brLenScaler = tr->partitionData[model].brLenScaler;
	      localTree->partitionData[model].propInvariant = tr->partitionData[model].propInvariant;
	      localTree->partitionData[model].lower      = tr->partitionData[model].lower;
	      localTree->partitionData[model].upper      = tr->partitionData[model].upper;
	      	      
	      localTree->partitionData[model].numberOfCategories      = tr->partitionData[model].numberOfCategories;
	      
	      localTree->partitionData[model].invariableWeight = tr->partitionData[model].invariableWeight;
	      memcpy(localTree->partitionData[model].invariableFrequencies,   tr->partitionData[model].invariableFrequencies,   pl->states * sizeof(double));
	    }	   
	}     

       for(model = 0; model < localTree->NumberOfModels; model++)
	 {
	   int 
	     localIndex;
	   
	   for(i = localTree->partitionData[model].lower, localIndex = 0; i <  localTree->partitionData[model].upper; i++)
	     if(i % (size_t)n == (size_t)tid)
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
	  localTree->lhs = tr->lhs;
	  localTree->lower_spacing = tr->lower_spacing;
	  localTree->upper_spacing = tr->upper_spacing;
	}
     
      optRateCatPthreads(localTree, localTree->lower_spacing, localTree->upper_spacing, localTree->lhs, n, tid);

      
      break;
    case THREAD_COPY_RATE_CATS:
      if(tid > 0)		
	broadcastPerSiteRates(tr, localTree);       

      for(model = 0; model < localTree->NumberOfModels; model++)
	{
	  localTree->partitionData[model].numberOfCategories = tr->partitionData[model].numberOfCategories;

	  for(localCounter = 0, i = localTree->partitionData[model].lower;  i < localTree->partitionData[model].upper; i++)
	    {
	      if(i % (size_t)n == (size_t)tid)
		{		 
		  localTree->partitionData[model].rateCategory[localCounter] = tr->cdta->rateCategory[i];		  				 
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
	      if(i % (size_t)n == (size_t)tid)
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
	      
	      copyLG4(localTree, tr, model, pl);
	     	     
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
	 
	  
	  memcpy(localTree->partitionContributions, tr->partitionContributions, sizeof(double) * localTree->NumberOfModels);
	  
	  


	  if(localTree->perPartitionEPA)
	    {
	      localTree->readPartition = (int *)rax_malloc(sizeof(int) * (size_t)localTree->numberOfTipsForInsertion);
	      memcpy(localTree->readPartition, tr->readPartition, sizeof(int) * (size_t)localTree->numberOfTipsForInsertion);
	    }

	}                                                

      localTree->temporarySumBuffer = (double *)rax_malloc(sizeof(double) * localTree->contiguousVectorLength);
      localTree->temporaryVector  = (double *)rax_malloc(sizeof(double) * localTree->contiguousVectorLength);      

      localTree->temporaryScaling = (int *)rax_malloc(sizeof(int) * localTree->contiguousScalingLength);
                 
      
      localTree->contiguousWgt          = (int*)rax_malloc(sizeof(int) * localTree->contiguousScalingLength);
      localTree->contiguousInvariant    = (int*)rax_malloc(sizeof(int) * localTree->contiguousScalingLength);	  
      
     
      memcpy(localTree->contiguousWgt         , tr->cdta->aliaswgt,     sizeof(int) * localTree->contiguousScalingLength);
      memcpy(localTree->contiguousInvariant   , tr->invariant,          sizeof(int) * localTree->contiguousScalingLength);
      
      if(tid > 0)
	broadcastPerSiteRates(tr, localTree);

     
      localTree->contiguousRateCategory = (int*)rax_malloc(sizeof(int) * localTree->contiguousScalingLength);
      
     
      memcpy(localTree->contiguousRateCategory, tr->cdta->rateCategory, sizeof(int) * localTree->contiguousScalingLength);           
     
      localTree->contiguousTips = tr->yVector;	  	
	 

      //NUMA optimizations, I adapted the original suggestion by Alexey Kozlov
      //initially all vectors in the branch data structure are allocated in function 
      //allocBranchX() in callssify.c 
      //Hence I can call the NUMA first touches here because this parallel region is invoked 
      //immediately after the call to allocBranchX(tr); in classify.c 
      //doing an actual memset might be an overkill, normally it would be sufficient 
      //to set/write a single byte per page size
      {
	size_t 
	  vectorChunkSize = tr->contiguousVectorLength / n,
	  scalingChunkSize = tr->contiguousScalingLength / n,
	  vectorOffset = vectorChunkSize * tid,
	  scalingOffset = scalingChunkSize * tid;
	
	int 
	  branchNumber;
	
	if(tid == n-1)
	  {
	    vectorChunkSize = tr->contiguousVectorLength - vectorOffset;
	    scalingChunkSize = tr->contiguousScalingLength - scalingOffset;
	    
	    assert(vectorChunkSize > 0 && scalingChunkSize > 0);
	  }
	
	for(branchNumber = 0; branchNumber < tr->numberOfBranches; branchNumber++)
	  {
	    double
	      *leftContigousVector = localTree->bInf[branchNumber].epa->left,
	      *rightContigousVector = localTree->bInf[branchNumber].epa->right;
	    
	    int
	      *leftContigousScalingVector = localTree->bInf[branchNumber].epa->leftScaling,
	      *rightContigousScalingVector = localTree->bInf[branchNumber].epa->rightScaling;
	    
	    
	    memset(&leftContigousVector[vectorOffset], 0, sizeof(double) * vectorChunkSize);
	    memset(&leftContigousScalingVector[scalingOffset], 0, sizeof(int) * scalingChunkSize);
	    
	    memset(&rightContigousVector[vectorOffset], 0, sizeof(double) * vectorChunkSize);
	    memset(&rightContigousScalingVector[scalingOffset], 0, sizeof(int) * scalingChunkSize);	  
	  }
      }      
      break;         
    case THREAD_FREE_VECTORS:
      {	
	size_t 
	  i;
      
	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    for(i = 0; i < localTree->innerNodes; i++)
	      {	    		      
		rax_free(localTree->partitionData[model].xVector[i]);
		rax_free(localTree->partitionData[model].expVector[i]);				
	      }
	  }
      }
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
	  rightNumber = localTree->bInf[branchCounter].epa->rightNodeNumber,
	  leftNumber  = localTree->bInf[branchCounter].epa->leftNodeNumber;	

	size_t
	  globalColumnCount = 0,
	  globalCount       = 0;

	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    size_t
	      blockRequirements;

	    double
	      *leftStridedVector  =  (double *)NULL,
	      *rightStridedVector =  (double *)NULL;

	    int
	      *leftStridedScalingVector  =  (int *)NULL,
	      *rightStridedScalingVector =  (int *)NULL;

	    size_t
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
		if(globalColumnCount % (size_t)n == (size_t)tid)
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
    case THREAD_PREPARE_BIPS_FOR_PRINT:
      {       
	int 
	  i = 0, 
	  j = 0;	
	
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
		entry 
		  *bipA = tr->consensusBips[i] ; 

		unsigned int 
		  firstIndex = 0;
		
		while(firstIndex < tr->bitVectorLength && bipA->bitVector[firstIndex] == 0 )
		  firstIndex++;		


		for(j = i + 1; j < tr->consensusBipLen; j++)
		  {
		    entry
		      *bipB = tr->consensusBips[j]; 
		    
		    if(bipA->amountTips < bipB->amountTips &&
		       issubset(bipA->bitVector, bipB->bitVector, tr->bitVectorLength, firstIndex))
		      { 
			/* i is child of j */		    
			IdList 
			  *elem = (IdList*) rax_calloc(1,sizeof(IdList));
			elem->value = i;
			
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

		    if(!((unsigned int)tr->mr_thresh < currentEntry->supportFromTreeset[0])) 
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
		     && (unsigned int)tr->sectionEnd == tr->h->entryCount /* the end of the buffer is also the hashtable  */
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
			  newSectionEnd  = MIN(tr->sectionEnd + tmp, (int)(tr->h->entryCount));
			}
		      else
			newSectionEnd = tr->h->entryCount; 
		  
		      density = 0.0;

		      tr->bipStatusLen = newSectionEnd - tr->sectionEnd;
		      rax_free(tr->bipStatus);
		      /* printf("%d\n" ,tr->bipStatusLen); */
		      tr->bipStatus = (int*)rax_calloc(tr->bipStatusLen, sizeof(int));
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
      
	size_t	 
	  globalColumnCount = 0,
	  globalCount       = 0;	

	for(model = 0; model < localTree->NumberOfModels; model++)
	  {
	    size_t
	      rateHet,	    	 
	      blockRequirements;
	    
	   
	    size_t
	      localColumnCount = 0,
	      localCount = 0;	   

	    double 
	      *stridedVector = localTree->partitionData[model].sumBuffer;
	    
	    if(tr->rateHetModel == CAT)
	      rateHet = 1;
	    else
	      rateHet = 4;
	    
	    blockRequirements = (size_t)(rateHet) * (size_t)(tr->partitionData[model].states);

	    for(globalColumnCount = localTree->partitionData[model].lower; globalColumnCount < localTree->partitionData[model].upper; globalColumnCount++)
	      {	
		if(globalColumnCount % (size_t)n == (size_t)tid)
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
	{
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));	 	 
	
	  for(model = 0; model < localTree->NumberOfModels; model++)	    	      
	    localTree->partitionData[model].brLenScaler = tr->partitionData[model].brLenScaler;
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
    case THREAD_COPY_LG4X_RATES:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      memcpy(localTree->partitionData[model].weights,    tr->partitionData[model].weights,    sizeof(double) * 4);
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	    }
	}
      break;
    case THREAD_OPT_LG4X_RATES:
      if(tid > 0)
	{
	  memcpy(localTree->executeModel, tr->executeModel, localTree->NumberOfModels * sizeof(boolean));

	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      memcpy(localTree->partitionData[model].weights,    tr->partitionData[model].weights,    sizeof(double) * 4);
	      memcpy(localTree->partitionData[model].gammaRates, tr->partitionData[model].gammaRates, sizeof(double) * 4);
	    }
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
    case THREAD_COPY_LG4X_EIGN:
      if(tid > 0)
	{
	  for(model = 0; model < localTree->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
		{
		  memcpy(localTree->partitionData[model].EIGN_LG4[0],    tr->partitionData[model].EIGN_LG4[0],    sizeof(double) * 19);
		  memcpy(localTree->partitionData[model].EIGN_LG4[1],    tr->partitionData[model].EIGN_LG4[1],    sizeof(double) * 19);
		  memcpy(localTree->partitionData[model].EIGN_LG4[2],    tr->partitionData[model].EIGN_LG4[2],    sizeof(double) * 19);
		  memcpy(localTree->partitionData[model].EIGN_LG4[3],    tr->partitionData[model].EIGN_LG4[3],    sizeof(double) * 19);
		}
	    }
	}
      break;
    case THREAD_SETUP_PRESENCE_MAP:
      setupPresenceMask(localTree);
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
    *localTree = (tree *)rax_malloc(sizeof(tree));
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

#if (defined(_WAYNE_MPI) && defined(_USE_PTHREADS))

static int setPthreadAffinity(void)
/*
 * Note: By setting the cpuset for the threads to include
 * all available cpus, we allow the kernel to decide on
 * which cpu/core each thread should execute.  Without this,
 * in the hybrid build of the code, the threads inherited
 * the OpenMPI-imposed affinity of the parent MPI process
 * forcing them to time-slice on a single core.  This
 * negated any potential performance improvement from
 * multiple threads on multi-core systems.
 *
 * If you set the processor affinity of the threads by
 * undef'ing _PORTABLE_PTHREADS, the tid is used to
 * create the cpuset.  This won't work in batch
 * environments with shared systems because the
 * processors corresponding to the tid may not be
 * available to the user's job on a given system.
 */
{
  int         
    maxCores, 
    cpuNum, 
    rc;
  
  int64_t        
    systemCores;
  
  cpu_set_t  
    *cpuSetPtr;
  
  size_t      
    setSize;
  /*
   * Get the number of cores available on the system. If
   * it is more than the glibc cpuset interface can handle,
   * set maxCores to the current max (CPU_SETSIZE).
   */
  
  systemCores = sysconf(_SC_NPROCESSORS_ONLN);

  if(systemCores > CPU_SETSIZE)
    maxCores = CPU_SETSIZE;  
  else
    maxCores = (int)systemCores;
 
  /* TODO */

  //printf("Number of cores available is %d\n", maxCores);

  /*
   * Allocate and configure the cpuset with all available
   * processors i.e. cores.
   */
  
  setSize   = CPU_ALLOC_SIZE(maxCores);
  cpuSetPtr = CPU_ALLOC(maxCores);
  
  if(cpuSetPtr == NULL) 
    {
      perror("CPU_ALLOC");
      rc = -1;
      return(rc);
    }
  CPU_ZERO_S(setSize, cpuSetPtr);
  
  for(cpuNum=0; cpuNum<maxCores; cpuNum++)
    CPU_SET_S(cpuNum, setSize, cpuSetPtr);
 
  /*
   * Set the processor affinity of the calling process using the
   * configured cpu set (i.e. bitmask).
   */
  
  printf("Modifying processor affinity for RAxML Master Pthread.\n");
  rc = sched_setaffinity(0, sizeof(cpu_set_t), cpuSetPtr);
  
  if(rc != 0)    
    perror("sched_setaffinity"); 

  CPU_FREE(cpuSetPtr);
  
  return rc;
}

#endif

static void startPthreads(tree *tr, analdef *adef)
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

  threads    = (pthread_t *)rax_malloc(NumberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)rax_malloc(NumberOfThreads * sizeof(threadData));
  
  
  reductionBuffer          = (volatile double *)rax_malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferTwo       = (volatile double *)rax_malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferThree     = (volatile double *)rax_malloc(sizeof(volatile double) *  NumberOfThreads * tr->NumberOfModels);
  reductionBufferParsimony = (volatile int *)rax_malloc(sizeof(volatile int) *  NumberOfThreads);

  
  barrierBuffer            = (volatile char *)rax_malloc(sizeof(volatile char) *  NumberOfThreads);
  
#if (defined(_WAYNE_MPI) && defined(_USE_PTHREADS))
  /*
   * Set the processor affinity of this (the master) thread.  Subsequent
   * threads spawned via pthread_create() will inherit the affinity of
   * their parent thread.
   */
  
  if(adef->setThreadAffinity)
    {
      if(setPthreadAffinity() != 0) 
	{
	  printf("Warning: Error setting thread processor affinity.\n");
	  printf("         This may adversely impact multi-core performance.\n");
	}
    }
#endif

  for(t = 0; t < NumberOfThreads; t++)
    barrierBuffer[t] = 0;

 
  branchInfos              = (volatile branchInfo **)rax_malloc(sizeof(volatile branchInfo *) * NumberOfThreads);
 
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
    *bestVector = (double*)rax_malloc(sizeof(double) * tr->cdta->endsite);

  for(i = 0; i < tr->cdta->endsite; i++)
    weightSum += (double)(tr->cdta->aliaswgt[i]);

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
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

      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE, FALSE);
     

      if(tr->optimizeAllTrees)
	{
	  treeEvaluate(tr, 1);
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);
	}
      else
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

 
  rax_free(bestVector);
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
    *unsortedSites = (double*)rax_malloc(sizeof(double) * tr->rdta->sites);

  

  fprintf(tlf, "  %d  %d\n", tr->numberOfTrees, tr->rdta->sites);

  for(i = 0; i < tr->numberOfTrees; i++)
    {      
      int 	
	k, 
	j;
           
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE, FALSE);
      assert(tr->ntips == tr->mxtips);
      
      if(i == 0) 
	{
	  if(adef->useBinaryModelFile)
	    {
	      readBinaryModel(tr, adef);
	      evaluateGenericInitrav(tr, tr->start);
	      treeEvaluate(tr, 2);
	    }
	  else
	    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);	
	}
      else
	{
	  if(tr->optimizeAllTrees)
	    {
	      treeEvaluate(tr, 1);
	      evaluateGenericInitrav(tr, tr->start);
	      modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);
	    }
	  else
	    treeEvaluate(tr, 2);     
	}

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

  rax_free(unsortedSites); 
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
 
  bestT = (bestlist *) rax_malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);   
 
  list = (elw *)rax_malloc(sizeof(elw) * tr->numberOfTrees); 
   
  for(i = 0; i < tr->numberOfTrees; i++)
    {            
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE, FALSE);
      resetBranches(tr); 
      
      if(i == 0)
	{	

	  if(adef->useBinaryModelFile)
	    {
	      readBinaryModel(tr, adef);
	      evaluateGenericInitrav(tr, tr->start);
	      treeEvaluate(tr, 2);
	    }
	  else
	    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
	  
	  printBothOpen("Model optimization on first Tree: %f\n", tr->likelihood);	  	 
	}
      else
      	{	       
	  evaluateGenericInitrav(tr, tr->start);
           
	  /*
	    treeEvaluateProgressive(tr);
	    treeEvaluateRandom(tr, 2);      
	  */
	  if(tr->optimizeAllTrees)
	    {
	      treeEvaluate(tr, 1);
	      evaluateGenericInitrav(tr, tr->start);
	      modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);
	    }
	  else
	    treeEvaluate(tr, 2);
	}            

      list[i].tree = i;
      list[i].lh   = tr->likelihood;

      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);

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
    bestIndex = -1,  
    i,
    k,
    *originalRateCategories = (int*)rax_malloc(tr->cdta->endsite * sizeof(int)),
    *originalInvariant      = (int*)rax_malloc(tr->cdta->endsite * sizeof(int));

  int64_t 
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

  bootweights = (elw *)rax_malloc(sizeof(elw) * tr->numberOfTrees);

  rankTest = (elw **)rax_malloc(sizeof(elw *) * adef->multipleRuns);

  for(k = 0; k < adef->multipleRuns; k++)
    rankTest[k] = (elw *)rax_malloc(sizeof(elw) * tr->numberOfTrees);

  lhs = (double **)rax_malloc(sizeof(double *) * tr->numberOfTrees);

  for(k = 0; k < tr->numberOfTrees; k++)
    lhs[k] = (double *)rax_calloc(adef->multipleRuns, sizeof(double));


  lhweights = (double **)rax_malloc(sizeof(double *) * tr->numberOfTrees);

  for(k = 0; k < tr->numberOfTrees; k++)
    lhweights[k] = (double *)rax_calloc(adef->multipleRuns, sizeof(double));

  /* read in the first tree and optimize ML params on it */  
  
  treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE, FALSE);   
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
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

      /* read in new tree */
     
      treeReadLen(treeFile, tr, FALSE, FALSE, FALSE, adef, TRUE, FALSE);  
      
      if(tr->optimizeAllTrees)
	{
	  treeEvaluate(tr, 1);
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);
	}
      else
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

	  computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE, TRUE, adef);

	  evaluateGenericInitrav(tr, tr->start);

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
  
  rax_free(originalRateCategories);
  rax_free(originalInvariant);
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
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);

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

#ifdef _BASTIEN
	assert(0);
#endif 

	if(tr->multiBranch)
	  {
	    int k;

	    for(k = 0, x = 0.0; k < tr->numBranches; k++)
	      {
		assert(tr->partitionContributions[k] != -1.0);

		z = result[k];

		if (z < zmin)
		  z = zmin;
		x += -log(z) * tr->partitionContributions[k];
	      }
	  }
	else
	  {
	    z = result[0];
	    if (z < zmin)
	      z = zmin;
	    x = -log(z);
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
    *significanceCounter = (int*)rax_malloc(sizeof(int) * tr->cdta->endsite); 

  double 
    *reference  = (double*)rax_malloc(sizeof(double) *  tr->cdta->endsite);

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




static int sortLex(const void *a, const void *b)
{
  int 
    i = 0; 
  
  char 
    *aPtr = *(char**)a,
    *bPtr = *(char**)b; 
  
  while((aPtr[i] != '\0') && (bPtr[i] != '\0') && (aPtr[i] == bPtr[i]))
    i++;

  if((aPtr[i] == '\0') || (bPtr[i] == '\0'))
    return (bPtr[i] == '\0');
  
  return (aPtr[i] > bPtr[i]);   
}


static void extractTaxaFromTopology(tree *tr, rawdata *rdta, cruncheddata *cdta, char fileName[1024])
{
  FILE 
    *f = myfopen(fileName, "rb");

  char 
    **nameList,
    buffer[nmlngth + 2]; 

  int
    i = 0,
    c,
    taxaSize = 1024,
    taxaCount = 0;
   
  nameList = (char**)rax_malloc(sizeof(char*) * taxaSize);  

  while((c = fgetc(f)) != ';')
    {
      
      if(c == '(' || c == ',')
	{
	  c = fgetc(f);
	  if(c ==  '(' || c == ',')
	    {
	      ungetc(c, f);
	    }
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

	      if(taxaCount == taxaSize)
		{		  
		  taxaSize *= 2;
		  nameList = (char **)rax_realloc(nameList, sizeof(char*) * taxaSize, FALSE);		 
		}
	      
	      nameList[taxaCount] = (char*)rax_malloc(sizeof(char) * (strlen(buffer) + 1));
	      strcpy(nameList[taxaCount], buffer);
	     
	      taxaCount++;
			    
	      ungetc(c, f);
	    }
	}   
    }


  /* BEGIN ensuring no taxon occurs twice */
  {
    char 
      **taxList = (char **)rax_malloc(sizeof(char *) * (size_t)taxaCount); 
    
    for(i = 0; i < taxaCount; ++i)
      taxList[i] = nameList[i]; 
  
    qsort(taxList, taxaCount, sizeof(char**), sortLex); 
    
    for(i = 1; i < taxaCount; ++i)
      {	
	if(strcmp(taxList[i], taxList[i-1]) == 0)
	  {
	    printf("\n\nA taxon labelled by %s appears twice in the first tree of tree collection %s, exiting ...\n\n", taxList[i], bootStrapFile);
	    exit(-1);
	  }
      }

    rax_free(taxList);
  }
  /* END */

  
  printf("Found a total of %d taxa in first tree of tree collection %s\n", taxaCount, bootStrapFile);
  printf("Expecting all remaining trees in collection to have the same taxon set\n");

  rdta->numsp = taxaCount;

  tr->nameList = (char **)rax_malloc(sizeof(char *) * (taxaCount + 1));  
  for(i = 1; i <= taxaCount; i++)
    tr->nameList[i] = nameList[i - 1];
  
  rax_free(nameList);

  tr->rdta       = rdta;
  tr->cdta       = cdta;

  if (rdta->numsp < 4)
    {    
      printf("TOO FEW SPECIES, tree contains only %d species\n", rdta->numsp);
      assert(0);
    }

  tr->nameHash = initStringHashTable(10 * taxaCount);
  for(i = 1; i <= taxaCount; i++)   
    {
      printf("add [%s]\n", tr->nameList[i]);
      addword(tr->nameList[i], tr->nameHash, i);
    }

  fclose(f);
}


static void myfwrite(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  size_t  
    bytes_written = fwrite(ptr, size, nmemb, stream);

  assert(bytes_written = nmemb);
}


static void writeLG4(tree *tr, int model, int dataType, FILE *f, partitionLengths pLengths[MAX_MODEL])
{
  if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
    {
      int 
	k;

      for(k = 0; k < 4; k++)
	{
	  myfwrite(tr->partitionData[model].EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	  myfwrite(tr->partitionData[model].rawEIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	  myfwrite(tr->partitionData[model].EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	  myfwrite(tr->partitionData[model].EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	  myfwrite(tr->partitionData[model].frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	  myfwrite(tr->partitionData[model].tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	  myfwrite(tr->partitionData[model].substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	}
    }
}


void writeBinaryModel(tree *tr, analdef *adef)
{
  int   
    model,
    progVers = programVersionInt; 
  
  FILE 
    *f = myfopen(binaryModelParamsOutputFileName, "w"); 

  /* number of tips in the tree that generates the model data */

  myfwrite(&(tr->ntips), sizeof(int), 1, f);

  /* pattern compression */

  myfwrite(&adef->compressPatterns, sizeof(boolean), 1, f);

  /* rate heterogeneity model */

  myfwrite(&tr->rateHetModel, sizeof(int), 1, f);

  /* prog version */

  myfwrite(&progVers, sizeof(int), 1, f);

  /* cdta */   

  myfwrite(tr->cdta->rateCategory, sizeof(int), tr->rdta->sites + 1, f);
  myfwrite(tr->cdta->patrat, sizeof(double), tr->rdta->sites + 1, f);
  myfwrite(tr->cdta->patratStored, sizeof(double), tr->rdta->sites + 1, f);

  /* partition contributions */

  myfwrite(tr->partitionContributions, sizeof(double), tr->NumberOfModels, f);

  
    
  /* pInfo */
   
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	dataType = tr->partitionData[model].dataType;

      myfwrite(tr->partitionData[model].weightExponents, sizeof(double), 4, f);
      myfwrite(tr->partitionData[model].weights,         sizeof(double), 4, f);

      myfwrite(tr->partitionData[model].gammaRates, sizeof(double), 4, f);

      myfwrite(tr->partitionData[model].EIGN, sizeof(double), pLengths[dataType].eignLength, f);
      myfwrite(tr->partitionData[model].EV, sizeof(double),  pLengths[dataType].evLength, f);
      myfwrite(tr->partitionData[model].EI, sizeof(double),  pLengths[dataType].eiLength, f);  

      myfwrite(tr->partitionData[model].frequencies, sizeof(double),  pLengths[dataType].frequenciesLength, f);
      myfwrite(tr->partitionData[model].freqExponents, sizeof(double),  pLengths[dataType].frequenciesLength, f);

      myfwrite(tr->partitionData[model].tipVector, sizeof(double),  pLengths[dataType].tipVectorLength, f);  
      myfwrite(tr->partitionData[model].substRates, sizeof(double),  pLengths[dataType].substRatesLength, f);           
      myfwrite(&(tr->partitionData[model].alpha), sizeof(double),  1, f);
      myfwrite(&(tr->partitionData[model].propInvariant), sizeof(double), 1, f);

      myfwrite(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);

      myfwrite(&(tr->partitionData[model].protModels), sizeof(int), 1, f);
      myfwrite(&(tr->partitionData[model].autoProtModels), sizeof(int), 1, f);

      myfwrite(tr->partitionData[model].perSiteRates,          sizeof(double), tr->partitionData[model].numberOfCategories, f);
      myfwrite(tr->partitionData[model].unscaled_perSiteRates, sizeof(double), tr->partitionData[model].numberOfCategories, f);   

      writeLG4(tr, model, dataType, f, pLengths);
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


static void readLG4(tree *tr, int model, int dataType, FILE *f, partitionLengths pLengths[MAX_MODEL])
{
  if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
    {
      int 
	k;

      for(k = 0; k < 4; k++)
	{
	  myfread(tr->partitionData[model].EIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	  myfread(tr->partitionData[model].rawEIGN_LG4[k], sizeof(double), pLengths[dataType].eignLength, f);
	  myfread(tr->partitionData[model].EV_LG4[k], sizeof(double),  pLengths[dataType].evLength, f);
	  myfread(tr->partitionData[model].EI_LG4[k], sizeof(double),  pLengths[dataType].eiLength, f);    
	  myfread(tr->partitionData[model].frequencies_LG4[k], sizeof(double),  pLengths[dataType].frequenciesLength, f);
	  myfread(tr->partitionData[model].tipVector_LG4[k], sizeof(double),  pLengths[dataType].tipVectorLength, f);  
	  myfread(tr->partitionData[model].substRates_LG4[k], sizeof(double),  pLengths[dataType].substRatesLength, f);    
	}
    }
}

void readBinaryModel(tree *tr, analdef *adef)
{
  boolean 
    compressPatterns;

  int 
    progVers,
    rateHetModel,
    model;

  FILE 
    *f;


  printBothOpen("\nRAxML is reading a binary model file and not optimizing model params\n");

  f = fopen(binaryModelParamsInputFileName, "r");   

  /* number of tips in the tree that generated the model data */

  myfread(&(tr->binaryFile_ntips), sizeof(int), 1, f);

  /* pattern compression */

  myfread(&compressPatterns, sizeof(boolean), 1, f);
  if(tr->rateHetModel == CAT && adef->compressPatterns && adef->mode == CLASSIFY_ML)
    {     
      printf("\n\nError: You need to disable site pattern compression by specifying the \"-H\" command line option\n");
      printf("when generating and reading binary model checkpoints for the EPA under CAT!\n\n");
      errorExit(-1);
    }

  if(compressPatterns != adef->compressPatterns)
    {
      printf("\n\nError you may need to disable pattern compression via the \"-H\" command line option!\n");
      printf("Or, when using CAT, disable it in the call you use to generate the binary model file.\n\n");
      errorExit(-1);
    }

  /* rate heterogeneity model */

  
  myfread(&rateHetModel, sizeof(int), 1, f);

  if(rateHetModel != tr->rateHetModel)
    {
      char 
	*modelNames[3] = {"CAT", "GAMMA", "GAMMAI"};
    
      printf("\n\nError: Rate heterogeneity models between binary model file that uses %s and the current command line that uses %s don't match \n\n\n", 
	     modelNames[rateHetModel], modelNames[tr->rateHetModel]);     
      errorExit(-1);
    } 

  /* program version  */

  myfread(&progVers, sizeof(int), 1, f);

   if(progVers != programVersionInt)
    {       
      printf("Error: Program versions between binary model file: %d and the current RAxML executable: %d don't match \n\n\n", 
	     progVers, programVersionInt);     
      errorExit(-1);
    } 

  /* cdta */   

  myfread(tr->cdta->rateCategory, sizeof(int),    (size_t)(tr->rdta->sites + 1), f);
  myfread(tr->cdta->patrat,       sizeof(double), (size_t)(tr->rdta->sites + 1), f);
  myfread(tr->cdta->patratStored, sizeof(double), (size_t)(tr->rdta->sites + 1), f);

  /* partition contributions */

  myfread(tr->partitionContributions, sizeof(double), tr->NumberOfModels, f);
  
 
  
  /* pInfo */
   
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	dataType = tr->partitionData[model].dataType;
      
      myfread(tr->partitionData[model].weightExponents, sizeof(double), 4, f);
      myfread(tr->partitionData[model].weights,         sizeof(double), 4, f);

      myfread(tr->partitionData[model].gammaRates, sizeof(double), 4, f);
      
      myfread(tr->partitionData[model].EIGN, sizeof(double), (size_t)(pLengths[dataType].eignLength), f);
      myfread(tr->partitionData[model].EV,   sizeof(double), (size_t)(pLengths[dataType].evLength), f);
      myfread(tr->partitionData[model].EI,   sizeof(double), (size_t)(pLengths[dataType].eiLength), f);  

      myfread(tr->partitionData[model].frequencies, sizeof(double),  (size_t)(pLengths[dataType].frequenciesLength), f);
      myfread(tr->partitionData[model].freqExponents, sizeof(double),  (size_t)(pLengths[dataType].frequenciesLength), f);     

      myfread(tr->partitionData[model].tipVector,   sizeof(double),  (size_t)(pLengths[dataType].tipVectorLength),   f);  
      myfread(tr->partitionData[model].substRates,  sizeof(double),  (size_t)(pLengths[dataType].substRatesLength),  f);     
      
      myfread(&(tr->partitionData[model].alpha),         sizeof(double), 1, f);
      myfread(&(tr->partitionData[model].propInvariant), sizeof(double), 1, f);

      myfread(&(tr->partitionData[model].numberOfCategories), sizeof(int), 1, f);
      
      myfread(&(tr->partitionData[model].protModels), sizeof(int), 1, f);
      myfread(&(tr->partitionData[model].autoProtModels), sizeof(int), 1, f);

      myfread(tr->partitionData[model].perSiteRates,          sizeof(double), tr->partitionData[model].numberOfCategories, f);
      myfread(tr->partitionData[model].unscaled_perSiteRates, sizeof(double), tr->partitionData[model].numberOfCategories, f);      

      readLG4(tr, model, dataType, f, pLengths);
    }

#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_INIT_MODEL, tr); 
  //masterBarrier(THREAD_RESET_MODEL, tr);
#endif  

  if(tr->rateHetModel == CAT)
    {
#ifdef _USE_PTHREADS
      masterBarrier(THREAD_COPY_RATE_CATS, tr);
#else
      {
	size_t 
	  i;
	int
	  model;	  

	for(model = 0; model < tr->NumberOfModels; model++)
	  {	  
	    int 
	      localCounter = 0;

	    for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++, localCounter++)	      
	      tr->partitionData[model].rateCategory[localCounter] = tr->cdta->rateCategory[i];	      
	  }      
      }
#endif
    }

  fclose(f);
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

    assert(sizeof(unsigned int) == 4);

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

/* functions to compute likelihoods on quartets */


/*** functions by Sarah for drawing quartets without replacement ***/

/*
Given the following nested for-loops:
for b = a+1 to n-2
  for c = b+1 to n-1
    for d = c+1 to n
How many iterations do we have for a given a?
*/
static uint64_t f2(int n, int a) 
{
  long double nDouble = n;
  long double aDouble = a;
  long double res = (nDouble - aDouble) * (nDouble - 1 - aDouble) * (nDouble - 2 - aDouble) / 6;
  return round(res);
};

/*
Given the following nested for-loops:
for c = b+1 to n-1
  for d = c+1 to n
How many iterations do we have for a given b?
*/
static uint64_t f3(int n, int b) 
{
  long double nDouble = n;
  long double bDouble = b;
  long double res = (nDouble - bDouble) * (nDouble - 1 - bDouble) / 2;
  return round(res);
};

/*
Given the following for-loop:
for d = c+1 to n
How many iterations do we have for a given c?
*/
static uint64_t f4(int n, int c) 
{
  return (n-c);
};

static void preprocessQuartetPrefix(int numberOfTaxa, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4)
{
  int 
    i,
    n = numberOfTaxa;
  
  /*
  Given the following nested for-loops:
  it = 0;
  for a = 1 to n-3
    for b = a+1 to n-2
      for c = b+1 to n-1
        for d = c+1 to n
          it++;
  prefixSumF2[i]: first value of it that belongs to a = i+1 
  prefixSumF3[i]: first value of it that belongs to b = i+2
  prefixSumF4[i]: first value of it that belongs to c = i+3
  */
  prefixSumF2[0] = 1;
  prefixSumF3[0] = 1;
  prefixSumF4[0] = 1;
  
  for (i = 1; i < n - 3; ++i) 
    {
      prefixSumF2[i] = prefixSumF2[i - 1] + f2(n, i);
      prefixSumF3[i] = prefixSumF3[i - 1] + f3(n, i+1);
      prefixSumF4[i] = prefixSumF4[i - 1] + f4(n, i+2);
  }
}

/*
Binary search in sorted array of size n-2. Returns the index of the greatest value in array that is <= z.
*/
static unsigned int binarySearch(uint64_t* array, uint64_t z, int n)
{
  unsigned int 
    first = 0,
    last = n-3,
    middle = (first + last) / 2, 
    lastSmallerOrEqual = 0;
  
  while(first <= last)
    {
      if(array[middle] < z)
	{
	  first = middle + 1;
	  lastSmallerOrEqual = middle;
	}
      else 
	{
	  if (array[middle] > z)	  
	    last = middle-1;	 
	  else 
	    { 
	      // array[middle] == z
	      lastSmallerOrEqual = middle;
	      break;
	    }
	}
      
      middle = (first + last)/2;
    }

  return lastSmallerOrEqual;
}

/**
Map an integer value z to a quartet (t1,t2,t3,t4).

@param numberOfTaxa The number of taxa in the tree.
@param z A value encoding a quartet (t1,t2,t3,t4). 
*/
static void mapNumberToQuartet(int numberOfTaxa, uint64_t z, int *t1, int *t2, int *t3, int *t4, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4)
{
  /*
  Given the following nested for-loops:
  z = 0;
  for t1 = 1 to numberOfTaxa-3
    for t2 = t1+1 to numberOfTaxa-2
      for t3 = t2+1 to numberOfTaxa-1
        for t4 = t3+1 to numberOfTaxa
          z++;
  Find the quartet (t1,t2,t3,t4) that belongs to the given value of z.
  */
  
  uint64_t    
    wantedT1 = z;

  // find the first value of z that belongs to t1
  *t1 = binarySearch(prefixSumF2, z, numberOfTaxa) + 1;

  uint64_t 
    foundT1 = prefixSumF2[*t1 - 1];
  
  if(wantedT1 == foundT1) 
    {
      *t2 = *t1+1;
      *t3 = *t1+2;
      *t4 = *t1+3;
      return;
    }
  
  uint64_t 
    wantedT2 = (prefixSumF3[*t1 - 1]) + (wantedT1 - foundT1);
  
  // find the first value of z that belongs to t2
  *t2 = binarySearch(prefixSumF3, wantedT2, numberOfTaxa) + 2;

  uint64_t 
    foundT2 = prefixSumF3[*t2 - 2];
  
  if(wantedT2 == foundT2) 
    {
      *t3 = *t2 + 1;
      *t4 = *t2 + 2;
      return;
    }
  
  uint64_t 
    wantedT3 = (prefixSumF4[*t2 - 2]) + (wantedT2 - foundT2);
  
  // find the first value of z that belongs to t3
  *t3 = binarySearch(prefixSumF4, wantedT3, numberOfTaxa) + 3;

  uint64_t 
    foundT3 = prefixSumF4[*t3 - 3];
  
  if (wantedT3 == foundT3) 
    {
      *t4 = *t3 + 1;
      return;
    }

  // find the value of z that belongs to t4
  *t4 = wantedT3 - foundT3 + *t3 + 1;
}


/* a parser error function */

static void parseError(int c)
{
  printf("Quartet grouping parser expecting symbol: %c\n", c);
  assert(0);
}

/* parser for the taxon grouping format, one has to specify 4 groups in a newick-like 
   format from which quartets (a substantially smaller number compared to ungrouped quartets) 
   will be drawn */

static void groupingParser(char *quartetGroupFileName, int *groups[4], int groupSize[4], tree *tr)
{
  FILE 
    *f = myfopen(quartetGroupFileName, "r");
  
  int 
    taxonCounter = 0,
    n,
    state = 0,
    groupCounter = 0,
    ch,
    i;

  printf("%s\n", quartetGroupFileName);

  for(i = 0; i < 4; i++)
    {
      groups[i] = (int*)rax_malloc(sizeof(int) * (tr->mxtips + 1));
      groupSize[i] = 0;
    }
  
  while((ch = getc(f)) != EOF)
    {
      if(!whitechar(ch))
	{
	  switch(state)
	    {
	    case 0:
	      if(ch != '(')
		parseError('(');
	      state = 1;
	      break;
	    case 1:
	      ungetc(ch, f);
	      n = treeFindTipName(f, tr, FALSE);  
	      if(n <= 0 || n > tr->mxtips)		
		printf("parsing error, raxml is expecting to read a taxon name, found \"%c\" instead\n", ch);		
	      assert(n > 0 && n <= tr->mxtips);	     
	      taxonCounter++;
	      groups[groupCounter][groupSize[groupCounter]] = n;
	      groupSize[groupCounter] = groupSize[groupCounter] + 1;	    
	      state = 2;
	      break;
	    case 2:
	      if(ch == ',')
		state = 1;
	      else
		{
		  if(ch == ')')
		    {
		      groupCounter++;
		      state = 3;
		    }
		  else
		    parseError('?');
		}
	      break;
	    case 3:
	      if(groupCounter == 4)
		{
		  if(ch == ';')
		    state = 4;
		  else
		    parseError(';');
		}
	      else
		{
		  if(ch != ',')
		    parseError(',');
		  state = 0;
		}
	      break; 
	    case 4:
	      printf("Error: extra char after ; %c\n", ch);
	      assert(0);
	    default:
	      assert(0);
	    }
	}
    }

  assert(state == 4);
  assert(groupCounter == 4); 
  assert(taxonCounter == tr->mxtips);

  printBothOpen("Successfully parsed quartet groups\n\n");

  /* print out the taxa that have been assigned to the 4 groups */

  for(i = 0; i < 4; i++)
    {
      int 
	j;
      
      printBothOpen("group %d has %d members\n", i, groupSize[i]);

      for(j = 0; j < groupSize[i]; j++)
	printBothOpen("%s\n", tr->nameList[groups[i][j]]);

      printBothOpen("\n");
    }

  fclose(f);
}


static double quartetLikelihood(tree *tr, nodeptr p1, nodeptr p2, nodeptr p3, nodeptr p4, nodeptr q1, nodeptr q2, analdef *adef, boolean firstQuartet)
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

  
#ifdef __BLACKRIM 
  if(firstQuartet)
    {
      tr->start = q1->next->back;
  
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
    }
#endif
  
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

#ifdef _QUARTET_MPI

typedef struct 
{
  int a1;
  int b1;
  int c1; 
  int d1;
  
  int a2;
  int b2;
  int c2; 
  int d2;

  int a3;
  int b3;
  int c3; 
  int d3;

  double l1;
  double l2;
  double l3;
} quartetResult;

#define QUARTET_MESSAGE_SIZE sizeof(quartetResult)  
#define QUARTET_MESSAGE 0
#define I_AM_DONE       1

static void startQuartetMaster(tree *tr, FILE *f)
{
  quartetResult 
    *qr = (quartetResult *)rax_malloc(sizeof(quartetResult));
  
  MPI_Status 
    status,
    recvStatus;

  int 
    dummy,
    workersDone = 0;
  
  assert(processID == 0);

  while(1)
    {
      MPI_Probe(MPI_ANY_SOURCE,  MPI_ANY_TAG, MPI_COMM_WORLD, &status);

      switch(status.MPI_TAG)
	{
	case QUARTET_MESSAGE:
	  MPI_Recv((void *)(qr), QUARTET_MESSAGE_SIZE, MPI_BYTE, status.MPI_SOURCE, QUARTET_MESSAGE, MPI_COMM_WORLD, &recvStatus);
	  fprintf(f, "%d %d | %d %d: %f\n", qr->a1, qr->b1, qr->c1, qr->d1, qr->l1);
	  fprintf(f, "%d %d | %d %d: %f\n", qr->a2, qr->b2, qr->c2, qr->d2, qr->l2);
	  fprintf(f, "%d %d | %d %d: %f\n", qr->a3, qr->b3, qr->c3, qr->d3, qr->l3);
	  break;
	case I_AM_DONE:
	  MPI_Recv(&dummy, 1, MPI_INT, status.MPI_SOURCE, I_AM_DONE, MPI_COMM_WORLD, &recvStatus);
	  workersDone++;
	  if(workersDone == processes -1)
	    goto END_IT;
	  break;
	default:
	  assert(0);
	}
    }
    
 END_IT:
  rax_free(qr);
  return;      
}

#endif



static void computeAllThreeQuartets(tree *tr, nodeptr q1, nodeptr q2, int t1, int t2, int t3, int t4, FILE *f, analdef *adef)
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

#ifdef _QUARTET_MPI
  quartetResult 
    *qr = (quartetResult *)rax_malloc(sizeof(quartetResult));
#endif
  
  /* first quartet */	    
  
  /* compute the likelihood of tree ((p1, p2), (p3, p4)) */
  
  l = quartetLikelihood(tr, p1, p2, p3, p4, q1, q2, adef, TRUE);
 
#ifndef _QUARTET_MPI
  fprintf(f, "%d %d | %d %d: %f\n", p1->number, p2->number, p3->number, p4->number, l);
#else
  qr->a1 = p1->number;
  qr->b1 = p2->number;
  qr->c1 = p3->number;
  qr->d1 = p4->number;
  qr->l1 = l;
#endif
  /* second quartet */	    
  
  /* compute the likelihood of tree ((p1, p3), (p2, p4)) */
  
  l = quartetLikelihood(tr, p1, p3, p2, p4, q1, q2, adef, FALSE);

#ifndef _QUARTET_MPI  
  fprintf(f, "%d %d | %d %d: %f\n", p1->number, p3->number, p2->number, p4->number, l);
#else
  qr->a2 = p1->number;
  qr->b2 = p3->number;
  qr->c2 = p2->number;
  qr->d2 = p4->number;
  qr->l2 = l;
#endif
  /* third quartet */	    
  
  /* compute the likelihood of tree ((p1, p4), (p2, p3)) */
  
  l = quartetLikelihood(tr, p1, p4, p2, p3, q1, q2, adef, FALSE);
  
#ifndef _QUARTET_MPI
  fprintf(f, "%d %d | %d %d: %f\n", p1->number, p4->number, p2->number, p3->number, l);	    	   
#else
  qr->a3 = p1->number;
  qr->b3 = p4->number;
  qr->c3 = p2->number;
  qr->d3 = p3->number;
  qr->l3 = l;

  MPI_Send((void *)qr, QUARTET_MESSAGE_SIZE, MPI_BYTE, 0, QUARTET_MESSAGE, MPI_COMM_WORLD);

  assert(processID > 0);
  rax_free(qr);
#endif
}

/* the three quartet options: all quartets, randomly sub-sample a certain number n of quartets, 
   subsample all quartets from 4 pre-defined groups of quartets */

#define ALL_QUARTETS 0
#define RANDOM_QUARTETS 1
#define GROUPED_QUARTETS 2

/**
Sample random quartets in ascending order using the methodA algorithm from J. S. Vitter, "An efficient algorithm for sequential random sampling". The runtime of this algorithm is O(numberOfQuartets).

@param tr The tree.
@param numberOfTaxa The number of taxa in the tree.
@param seed
@param numberOfQuartets The total number of different quartets that exist for numberOfTaxa taxa.
@param randomQuartets The number of quartets to sample.
@param q1
@param q2
@param prefixSumF2
@param prefixSumF3
@param prefixSumF4
@param f
@param adef
@param actVal The value of the last drawn random number representing a quartet.
*/
static void sampleQuartetsWithoutReplacementA(tree *tr, int numberOfTaxa, int64_t seed, uint64_t numberOfQuartets, uint64_t randomQuartets, nodeptr q1, nodeptr q2, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4, FILE *f, analdef *adef, uint64_t actVal)
{
  int64_t 
    myseed = seed;

  uint64_t    
    sampleSize = randomQuartets,
    quartetCounter = 0,
    top = numberOfQuartets - sampleSize,
    s;
  
  int 
    t1,
    t2,
    t3,
    t4;

  double 
    NReal = (double)numberOfQuartets, 
    v, 
    quot; 
  
  while(sampleSize >= 2)
    {
      v = randum(&myseed);
      s = 0;
      quot = top / NReal;
    
      while (quot > v)
	{
	  s++; 
	  top--; 
	  NReal--;
	  quot = (quot * top) / NReal;
	}
    // Skip over the next s records and select the following one for the sample
      actVal += s+1;
      mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
      computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
      quartetCounter++;
      
      NReal--;
      sampleSize--;
    }
  
  // Special case sampleSize == 1
  s = trunc(round(NReal) * randum(&myseed));
  // Skip over the next s records and select the following one for the sample
  actVal += s+1;
  
  mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
  #ifdef _QUARTET_MPI
				  //MPI version very simple and naive way to determine which processor 
				  //is going to do the likelihood calculations for this quartet
				  if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
  computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
  quartetCounter++;

  assert(quartetCounter == randomQuartets);
}

/**
Sample random quartets in ascending order using the methodD algorithm from J. S. Vitter, "An efficient algorithm for sequential random sampling". The runtime of this algorithm is O(randomQuartets). The main idea of the algorithm is to decide ho many quartets to skip instead of testing for each quartet whether to take it or not.

@param tr The tree.
@param numberOfTaxa The number of taxa in the tree.
@param seed
@param numberOfQuartets The total number of different quartets that exist for numberOfTaxa taxa.
@param randomQuartets The number of quartets to sample.
@param q1
@param q2
@param prefixSumF2
@param prefixSumF3
@param prefixSumF4
@param f
@param adef
@param actVal The value of the last drawn random number representing a quartet.
*/
static void sampleQuartetsWithoutReplacementD(tree *tr, int numberOfTaxa, int64_t seed, uint64_t numberOfQuartets, uint64_t randomQuartets, nodeptr q1, nodeptr q2, uint64_t *prefixSumF2, uint64_t *prefixSumF3, uint64_t *prefixSumF4, FILE *f, analdef *adef, uint64_t actVal)
{
  int64_t       
    myseed = seed;
  
  uint64_t
    sampleSize = randomQuartets,
    quartetCounter = 0,
    s,
    qu1,
    threshold,
    t,
    limit;
    
  int 
    t1,
    t2,
    t3,
    t4;
    
  double
    negalphainv = -1.0/13,
    nreal = sampleSize,
    ninv = 1.0 / nreal,
    Nreal = numberOfQuartets,
    vprime = exp(log(randum(&myseed)) * ninv),
    qu1real,
    nmin1inv,
    x,
    u, 
    negSreal,
    y1,
    y2,
    top,
    bottom;
    
  qu1 = -sampleSize + 1 + numberOfQuartets;
  qu1real = -nreal + 1.0 + Nreal;
  threshold = -negalphainv * sampleSize;

  while((sampleSize > 1) && (threshold < numberOfQuartets))
  {
    nmin1inv = 1.0 / (-1.0 + nreal);
    while(TRUE)
      {
	while (TRUE)
	  // step D2: Generate U and X
	  {
	    x = Nreal * (-vprime + 1.0);
	    s = trunc(x);
	    if (s < qu1) break;
	    vprime = exp(log(randum(&myseed)) * ninv);
	  }
	u = randum(&myseed);
	negSreal = (double) s * (-1);
	// step D3: Accept?
	y1 = exp(log(u * Nreal / qu1real) * nmin1inv);
	vprime = y1 * (-x / Nreal + 1.0) * (qu1real / (negSreal + qu1real));
	if (vprime <= 1.0) break; // Accept! test (2.8) is true
	// step D4: Accept?
	y2 = 1.0;
	top = -1.0 + Nreal;
	if(-1 + sampleSize > s)
	  {
	    bottom = -nreal + Nreal;
	    limit = -s + numberOfQuartets;
	  }
	else
	  {
	    bottom = -1.0 + negSreal + Nreal;
	    limit = qu1;
	  }
	for (t = -1 + numberOfQuartets; t >= limit; t--)
	  {
	    y2 = (y2 * top)/bottom;
	    top = -1.0 + top;
	    bottom = -1.0 + bottom;
	  }
	
	if(Nreal / (-x + Nreal) >= y1 * exp(log(y2) * nmin1inv))
	  {
	    // Accept!
	    vprime = exp(log(randum(&myseed)) * nmin1inv);
	    break;
	  }
	vprime = exp(log(randum(&myseed)) * ninv);
      }
    // Step D5: Select the (s+1)st record
    // Skip over the next s records and select the following one for the sample
    actVal += s+1;
    mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
    computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
    quartetCounter++;
    
    numberOfQuartets = -s + (-1 + numberOfQuartets);
    Nreal = negSreal + (-1.0 + Nreal);
    sampleSize--;
    nreal = nreal - 1.0;
    ninv = nmin1inv;
    qu1 = qu1 - s;
    qu1real += negSreal;
    threshold += negalphainv;
  }
  if (sampleSize > 1)
    {
      // Use Method A to finish the sampling
      assert(quartetCounter == randomQuartets - sampleSize);
      sampleQuartetsWithoutReplacementA(tr, numberOfTaxa, seed, numberOfQuartets, sampleSize, q1, q2, prefixSumF2, prefixSumF3, prefixSumF4, f, adef, actVal);
    }
  else // Special case sampleSize == 1
    {
      s = trunc(numberOfQuartets * vprime);
      // Skip over the next s records and select the following one for the sample
      actVal += s+1;
      mapNumberToQuartet(numberOfTaxa, actVal, &t1, &t2, &t3, &t4, prefixSumF2, prefixSumF3, prefixSumF4);
      #ifdef _QUARTET_MPI
				  //MPI version very simple and naive way to determine which processor 
				  //is going to do the likelihood calculations for this quartet
				  if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
      computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
      quartetCounter++;
      assert(quartetCounter == randomQuartets);
    }
}

static void computeQuartets(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  /* some indices for generating quartets in an arbitrary way */

  int
    flavor = ALL_QUARTETS, //type of quartet calculation 
    i, 
    t1, 
    t2, 
    t3, 
    t4, 
    *groups[4],
    groupSize[4];    

  double
    fraction = 0.0, //fraction of random quartets to compute
    t;

  uint64_t
    randomQuartets = (uint64_t)(adef->multipleRuns), //number of random quartets to compute 
    quartetCounter = 0, 
    //total number of possible quartets, note that we count the following ((A,B),(C,D)), ((A,C),(B,D)), ((A,D),(B,C)) as one quartet here 
    numberOfQuartets = ((uint64_t)tr->mxtips * ((uint64_t)tr->mxtips - 1) * ((uint64_t)tr->mxtips - 2) * ((uint64_t)tr->mxtips - 3)) / 24; 
  
  /* use two inner tree nodes for building quartet trees */

  nodeptr 	
    q1 = tr->nodep[tr->mxtips + 1],
    q2 = tr->nodep[tr->mxtips + 2];

  char 
    quartetFileName[1024];

  FILE 
    *f;
       
  /***********************************/
 
  


  /* build output file name */
    
  strcpy(quartetFileName,         workdir);
  strcat(quartetFileName,         "RAxML_quartets.");
  strcat(quartetFileName,         run_id);
  
  /* open output file */

#ifdef _QUARTET_MPI
  if(processID == 0)
#endif
    f = myfopen(quartetFileName, "w");

  /* initialize model parameters */

  initModel(tr, rdta, cdta, adef);

  

  if(!adef->useBinaryModelFile)
    {
#ifdef _QUARTET_MPI
      //the parallel version requires a pre-computed model parameter file as input!
      assert(0);
#endif

      /* get a starting tree on which we optimize the likelihood model parameters: either reads in a tree or computes a randomized stepwise addition parsimony tree */

      getStartingTree(tr, adef);
   
      /* optimize model parameters on that comprehensive tree that can subsequently be used for evaluation of quartet likelihoods */

#ifndef __BLACKRIM //if BLACKRIM is defined, the model parameters will be optimized for each quartet individually
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
#endif

      printBothOpen("Time for parsing input tree or building parsimony tree and optimizing model parameters: %f\n\n", gettime() - masterTime); 

#ifndef __BLACKRIM
      printBothOpen("Tree likelihood: %f\n\n", tr->likelihood);
#endif
    }
  else
    {
      //if a binary model parameter file has been specified, we just read the model parameters from there 
      readBinaryModel(tr, adef);

      printBothOpen("Time for reading model parameters: %f\n\n", gettime() - masterTime); 
    }


  /* figure out which flavor of quartets we want to compute */

  if(adef->useQuartetGrouping)
    {
      //quartet grouping evaluates all possible quartets from four disjoint 
      //sets of user-specified taxon names 

      flavor = GROUPED_QUARTETS;
      
      //parse the four disjoint sets of taxon names specified by the user from file      
      groupingParser(quartetGroupingFileName, groups, groupSize, tr);

#ifdef __BLACKRIM     
      //special implementation where we only sub-sample randomly from the quartets 
      //defined by the four user-specified taxon sets 
      numberOfQuartets =  (uint64_t)groupSize[0] * (uint64_t)groupSize[1] * (uint64_t)groupSize[2] * (uint64_t)groupSize[3];

      if(randomQuartets > numberOfQuartets)
	randomQuartets = 1;

      fraction = (double)randomQuartets / (double)numberOfQuartets;     

      //printf("%d %d %f\n", numberOfQuartets, randomQuartets, fraction);
#endif
    }
  else
    {
      //if the user specified more random quartets to sample than there actually 
      //exist for the number of taxa, then fix this.
      if(randomQuartets > numberOfQuartets)
	randomQuartets = 1;
  
      if(randomQuartets == 1)   
	//change flavor if randomQuartets > possibleQuartets
	flavor = ALL_QUARTETS;
      else
	{      
	  //compute the fraction of random quartets to sample 
	  //there may be an issue here with the unit64_t <-> double cast
	  fraction = (double)randomQuartets / (double)numberOfQuartets;      
	  flavor = RANDOM_QUARTETS;
	}
    }

  /* print some output on what we are doing*/

  switch(flavor)
    {
    case ALL_QUARTETS:
      printBothOpen("There are %" PRIu64 " quartet sets for which RAxML will evaluate all %" PRIu64 " quartet trees\n", numberOfQuartets, numberOfQuartets * 3);
      break;
    case RANDOM_QUARTETS:
      printBothOpen("There are %" PRIu64 " quartet sets for which RAxML will randomly sub-sambple %" PRIu64 " sets (%f per cent), i.e., compute %" PRIu64 " quartet trees\n", 
		    numberOfQuartets, randomQuartets, 100 * fraction, randomQuartets * 3);
      break;
    case GROUPED_QUARTETS:  
#ifdef __BLACKRIM          
      printBothOpen("There are 4 quartet groups from which RAxML will evaluate the three alternatives for %u out of all possible %u quartet trees\n", 
		    (unsigned int)randomQuartets, 
		    (unsigned int)groupSize[0] * (unsigned int)groupSize[1] * (unsigned int)groupSize[2] * (unsigned int)groupSize[3]);
#else
      printBothOpen("There are 4 quartet groups from which RAxML will evaluate all %u quartet trees\n", 
		    (unsigned int)groupSize[0] * (unsigned int)groupSize[1] * (unsigned int)groupSize[2] * (unsigned int)groupSize[3] * 3);
#endif
      break;
    default:
      assert(0);
    }

  /* print taxon name to taxon number correspondance table to output file */
#ifdef _QUARTET_MPI
  if(processID == 0)   
#endif
    {
      fprintf(f, "Taxon names and indices:\n\n");

      for(i = 1; i <= tr->mxtips; i++)
	{
	  fprintf(f, "%s %d\n", tr->nameList[i], i);
	  assert(tr->nodep[i]->number == i);
	}
      
      fprintf(f, "\n\n");
    }


  t = gettime();
  
  /* do a loop to generate some quartets to test.
     note that tip nodes/sequences in RAxML are indexed from 1,...,n
     and not from 0,...,n-1 as one might expect 
     
     tr->mxtips is the maximum number of tips in the alignment/tree
  */


  //now do the respective quartet evaluations by switching over the three distinct flavors 

#ifdef _QUARTET_MPI
  if(processID > 0)   
#endif
    {
      switch(flavor)
	{
	case ALL_QUARTETS:
	  {
	    assert(randomQuartets == 1);
	    
	    /* compute all possible quartets */
	    
	    for(t1 = 1; t1 <= tr->mxtips; t1++)
	      for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
		for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
		  for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
		    {
#ifdef _QUARTET_MPI
		      if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
			computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
		      quartetCounter++;
		    }
	    
	    assert(quartetCounter == numberOfQuartets);
	  }
	  break;
	case RANDOM_QUARTETS:
	  {	 
	    //code contributed by Sarah for drawing quartets without replacement :-)
        // Sample random quartets without replacement in O(randomQuartets * log(tr->mxtips)) time and O(tr->mxtips) space.
        // This is achieved by drawing random numbers in ascending order and using prefix sums to map a number to a
        // quartet (t1,t2,t3,t4) using the lexicographical ordering of the quartets. For each quartet, it is required
        // that 1 <= t1 < t2 < t3 < t4 <= tr->mxtips.
	    
	    if(adef->sampleQuartetsWithoutReplacement)
	      {
		uint64_t
		  *prefixSumF2 = (uint64_t*)rax_malloc(sizeof(uint64_t) * (size_t)(tr->mxtips - 2)),
		  *prefixSumF3 = (uint64_t*)rax_malloc(sizeof(uint64_t) * (size_t)(tr->mxtips - 2)),
		  *prefixSumF4 = (uint64_t*)rax_malloc(sizeof(uint64_t) * (size_t)(tr->mxtips - 2));

		preprocessQuartetPrefix(tr->mxtips, prefixSumF2, prefixSumF3, prefixSumF4);

		if (randomQuartets >= numberOfQuartets/13) // decide for each quartet whether to take it or not
		  sampleQuartetsWithoutReplacementA(tr, tr->mxtips, adef->parsimonySeed, numberOfQuartets, randomQuartets, q1, q2, prefixSumF2, prefixSumF3, prefixSumF4, f, adef, 0);
		else // decide how many quartets to skip before taking the next one
		  sampleQuartetsWithoutReplacementD(tr, tr->mxtips, adef->parsimonySeed, numberOfQuartets, randomQuartets, q1, q2, prefixSumF2, prefixSumF3, prefixSumF4, f, adef, 0);

		rax_free(prefixSumF2);
		rax_free(prefixSumF3);
		rax_free(prefixSumF4);
	      }
	    else
	      {
		//endless loop ta make sure we randomly sub-sample exactly as many quartets as the user specified

		//This is not very elegant, but it works, note however, that especially when the number of 
		//random quartets to be sampled is large, that is, close to the total number of quartets 
		//some quartets may be sampled twice by pure chance. To randomly sample unique quartets 
		//using hashes or bitmaps to store which quartets have already been sampled is not memory efficient.
		//Insetad, we need to use a random number generator that can generate a unique series of random numbers 
		//and then have a function f() that maps those random numbers to the corresponding index quartet (t1, t2, t3, t4),
		//see above 
		
		do
		  {	      
		    //loop over all quartets 
		    for(t1 = 1; t1 <= tr->mxtips; t1++)
		      for(t2 = t1 + 1; t2 <= tr->mxtips; t2++)
			for(t3 = t2 + 1; t3 <= tr->mxtips; t3++)
			  for(t4 = t3 + 1; t4 <= tr->mxtips; t4++)
			    {
			      //chose a random number
			      double
				r = randum(&adef->parsimonySeed);
			      
			      //if the random number is smaller than the fraction of quartets to subsample
			      //evaluate the likelihood of the current quartet
			      if(r < fraction)
				{
#ifdef _QUARTET_MPI
				  //MPI version very simple and naive way to determine which processor 
				  //is going to do the likelihood calculations for this quartet
				  if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
				    //function that computes the likelihood for all three possible unrooted trees 
				    //defined by the given quartet of taxa 
				    computeAllThreeQuartets(tr, q1, q2, t1, t2, t3, t4, f, adef);
				  //increment quartet counter that counts how many quartets we have evaluated
				  quartetCounter++;
				}
			      
			      //exit endless loop if we have randomly sub-sampled as many quartets as the user specified
			      if(quartetCounter == randomQuartets)
				goto DONE;
			    }
		  }
		while(1);
		
	      DONE:
		assert(quartetCounter == randomQuartets);	  
	      }
	  }
	  break;
	case GROUPED_QUARTETS:
	  {
	    /* compute all quartets that can be built out of the four pre-defined groups */
	    
	    for(t1 = 0; t1 < groupSize[0]; t1++)
	      for(t2 = 0; t2 < groupSize[1]; t2++)
		for(t3 = 0; t3 < groupSize[2]; t3++)
		  for(t4 = 0; t4 < groupSize[3]; t4++)
		    {
		      int
			i1 = groups[0][t1],
			i2 = groups[1][t2],
			i3 = groups[2][t3],
			i4 = groups[3][t4];
		      
#ifdef __BLACKRIM
		      double
			r = randum(&adef->parsimonySeed);
		      
		      if(r < fraction)
			{
#endif
			  
#ifdef _QUARTET_MPI
			  if((quartetCounter % (uint64_t)(processes - 1)) == (uint64_t)(processID - 1))
#endif
			    computeAllThreeQuartets(tr, q1, q2, i1, i2, i3, i4, f, adef);
			  quartetCounter++;
#ifdef __BLACKRIM
			}
		      if(quartetCounter == randomQuartets)
			goto DONE_GROUPED;
#endif
		    }
#ifdef __BLACKRIM   
	  DONE_GROUPED:
	    printBothOpen("\nComputed %" PRIu64 " random quartets for grouping\n", quartetCounter);
	    assert(quartetCounter == randomQuartets);
#else
	    printBothOpen("\nComputed all %" PRIu64 " possible grouped quartets\n", quartetCounter);
#endif	    
	    
	  }
	  break;
	default:
	  assert(0);
	}
    }
#ifdef _QUARTET_MPI
  if(processID == 0)
    startQuartetMaster(tr, f);
  else
    {
      int 
	dummy;
      
      MPI_Send(&dummy, 1, MPI_INT, 0, I_AM_DONE, MPI_COMM_WORLD);
    }
#endif
  
 

  t = gettime() - t;

  printBothOpen("\nPure quartet computation time: %f secs\n", t);
  
  printBothOpen("\nAll quartets and corresponding likelihoods written to file %s\n", quartetFileName);

#ifdef _QUARTET_MPI
  if(processID == 0)
#endif
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

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);

  Thorough = 1;
  tr->doCutoff = FALSE;  
	 
  printBothOpen("\nStart likelihood: %f\n\n", tr->likelihood);

  treeOptimizeThorough(tr, 1, 10);
  evaluateGenericInitrav(tr, tr->start);
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);

  printBothOpen("End likelihood: %f\n\n", tr->likelihood);

  printModelParams(tr, adef);    
  
  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_bestTree.");
  strcat(bestTreeFileName,         run_id);

  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);
  f = myfopen(bestTreeFileName, "wb");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);

  printBothOpen("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
}

static void evaluateSD(tree *tr, double bestLH, double *bestVector, double weightSum, int configuration, int i, FILE *f)
{
  double
    sum = 0.0,
    sum2 = 0.0,
    sd,
    currentLH;

  int 
    k;

  evaluateGenericInitrav(tr, tr->start);	  
  evaluateGenericVector(tr, tr->start); 	 	  
  
  currentLH = tr->likelihood;

  printBothOpen("Configuration %d Likelihood: %f\n", configuration, tr->likelihood);       

  fprintf(f, "tr%d\t", configuration);

  if(currentLH > bestLH)	
    printBothOpen("WARNING tree with ancestral sequence taxon %s has a better likelihood %f > %f than the reference tree!\n", tr->nameList[i], currentLH, bestLH);

  for (k = 0; k < tr->cdta->endsite; k++)
    {
      int 
	w;
      
      double 
	temp = bestVector[k] - tr->perSiteLL[k],
	wtemp = tr->cdta->aliaswgt[k] * temp;
     
      for(w = 0; w < tr->cdta->aliaswgt[k]; w++)
	fprintf(f, "%f ",  tr->perSiteLL[k]);
 
      sum  += wtemp;
      sum2 += wtemp * temp;
    }

  fprintf(f, "\n");

  sd = sqrt( weightSum * (sum2 - sum * sum / weightSum) / (weightSum - 1) );
	 
  printBothOpen("Ancestral Taxon: %s Likelihood: %f D(LH): %f SD: %f \nSignificantly Worse: %s (5%s), %s (2%s), %s (1%s)\n", 
		tr->nameList[i], currentLH, currentLH - bestLH, sd, 
		(sum > 1.95996 * sd) ? "Yes" : " No", "%",
		(sum > 2.326 * sd) ? "Yes" : " No", "%",
		(sum > 2.57583 * sd) ? "Yes" : " No", "%");

  printBothOpen("\n");
}

static void ancestralSequenceTest(tree *tr)
{
  FILE 
    *f = myfopen(quartetGroupingFileName, "r");

  int
    ch,
    i,
    *candidateAncestorList = (int *)rax_calloc((tr->mxtips + 1), sizeof(int)),
    numberOfCandidateAncestors = 0;  

  double 
    bestLH, 
    weightSum = 0.0,
    *bestVector = (double*)rax_malloc(sizeof(double) * tr->cdta->endsite);

  assert(tr->useFastScaling == FALSE);

  for(i = 0; i < tr->cdta->endsite; i++)
    weightSum += (double)(tr->cdta->aliaswgt[i]);
  
  evaluateGenericInitrav(tr, tr->start);
  evaluateGenericVector(tr, tr->start);
  
  bestLH = tr->likelihood;
  
  memcpy(bestVector, tr->perSiteLL, tr->cdta->endsite * sizeof(double));

  printBothOpen("Likelihood of reference tree: %f\n\n\n", tr->likelihood);  

  while((ch = getc(f)) != EOF)
    {
      if(!whitechar(ch))
	{	  
	  int 
	    n;
	  
	  ungetc(ch, f);
	  
	  n = treeFindTipName(f, tr, FALSE);  
	  
	  if(n <= 0 || n > tr->mxtips)		
	    printf("parsing error, raxml is expecting to read a taxon name that is contained in the reference tree you passed!\n");		
	  
	  assert(n > 0 && n <= tr->mxtips);	
	  
	  candidateAncestorList[n] = 1;
	  numberOfCandidateAncestors++;	     
	}
    }

  fclose(f);

  for(i = 1; i <= tr->mxtips; i++)
    {
      if(candidateAncestorList[i])
	{
	  nodeptr 	    
	    p = tr->nodep[i],
	    q = p->back,
	    l = q->next,
	    r = q->next->next;

	  int
	    k;	 

	  double 	  
	    attachmentBranch[NUM_BRANCHES],
	    leftBranch[NUM_BRANCHES],
	    rightBranch[NUM_BRANCHES];

	  FILE 
	    *f;

	  char 
	    fileName[1024];

	  strcpy(fileName, workdir);
	  strcat(fileName, "RAxML_ancestralTest.");
	  strcat(fileName, tr->nameList[i]);
	  strcat(fileName, ".");
	  strcat(fileName, run_id);

	  f  = myfopen(fileName, "w");	
 
	  fprintf(f, "  3  %d\n", tr->rdta->sites);

	  assert(strcmp(tr->nameList[i], tr->nameList[p->number]) == 0);

	  printBothOpen("Checking if %s is a candidate ancestor\n\n", tr->nameList[i]);
	  printBothOpen("Per site log likelihoods for the three configurations will be written to file %s\n\n", fileName);
	  
	  memcpy(attachmentBranch, p->z, sizeof(double) * NUM_BRANCHES);
	  memcpy(leftBranch, l->z, sizeof(double) * NUM_BRANCHES);
	  memcpy(rightBranch, r->z, sizeof(double) * NUM_BRANCHES);


	  //configuration 1

	  for(k = 0; k < NUM_BRANCHES; k++)
	    p->z[k] = q->z[k] = zmax;
	    	   
	  evaluateSD(tr, bestLH, bestVector, weightSum, 1, i, f);	

	  memcpy(p->z, attachmentBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(p->back->z, attachmentBranch, sizeof(double) * NUM_BRANCHES);

	  evaluateGenericInitrav(tr, tr->start);
	  assert(tr->likelihood == bestLH);

	  //configuration 2

	  for(k = 0; k < NUM_BRANCHES; k++)
	    {
	      p->z[k] = q->z[k] = zmax;
	      l->z[k] = l->back->z[k] = zmax;
	    }

	  evaluateSD(tr, bestLH, bestVector, weightSum, 2, i, f);	

	  memcpy(p->z, attachmentBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(p->back->z, attachmentBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(l->z, leftBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(l->back->z, leftBranch, sizeof(double) * NUM_BRANCHES);

	  evaluateGenericInitrav(tr, tr->start);
	  assert(tr->likelihood == bestLH);

	  //configuration 3

	  for(k = 0; k < NUM_BRANCHES; k++)
	    {
	      p->z[k] = q->z[k] = zmax;
	      r->z[k] = r->back->z[k] = zmax;
	    }

	  evaluateSD(tr, bestLH, bestVector, weightSum, 3, i, f);	

	  memcpy(p->z, attachmentBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(p->back->z, attachmentBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(r->z, rightBranch, sizeof(double) * NUM_BRANCHES);
	  memcpy(r->back->z, rightBranch, sizeof(double) * NUM_BRANCHES);
	  
	  evaluateGenericInitrav(tr, tr->start);
	  assert(tr->likelihood == bestLH);
	  
	  printBothOpen("\n\n");
	  fclose(f);
	}
    }

  printBothOpen("good-bye\n\n");

  rax_free(candidateAncestorList);
  rax_free(bestVector);
  exit(0);  
}

static double distancesInitial(nodeptr p, double *distances, tree *tr, boolean fullTraversal)
{
  if(isTip(p->number, tr->mxtips))
    return p->z[0];
  else
    {
      double 
	acc = 0.0;
      
      nodeptr 
	q;                
     
      if(fullTraversal || !p->x)
	{
	  q = p->next;      

	  while(q != p)
	    {	 
	      acc += distancesInitial(q->back, distances, tr, fullTraversal);
	      q = q->next;
	    }

	  distances[p->number] = acc;
	  p->x = 1;
	  p->next->x = 0;
	  p->next->next->x = 0;
	}
      else
	acc = distances[p->number];

      return acc + p->z[0];
    }
}



static void distancesNewview(nodeptr p, double *distances, tree *tr, nodeptr *rootBranch, double *minimum)
{ 
  nodeptr 
    q;                
      
  double 
    left = 0.0,
    right = 0.0;
  
  if(isTip(p->number, tr->mxtips))
    {          
      q = p->back;  

      if(!isTip(q->number, tr->mxtips))
	{
	  if(!q->x)
	    distancesInitial(q, distances, tr, FALSE);		  
	  left = distances[q->number];	  
	}                 

      if(left <= p->z[0])
	{	  
	  //the balanced root is in this branch
	  *rootBranch = p;
	  *minimum = 0.0;
	}
      else
	{	  	  	    
	  double 
	    diff = left - p->z[0];	    
	  
	  if(diff < *minimum)
	    { 	     
	      *minimum = diff;
	      *rootBranch = p;
	    }
	}	
    }
  else
    {          
      q = p->back;  

      if(!isTip(q->number, tr->mxtips))
	{
	  if(!q->x)
	    distancesInitial(q, distances, tr, FALSE);	
	  
	  left = distances[q->number];	  
	}
      else
	left = 0.0;
      
      if(!isTip(p->number, tr->mxtips))
	{
	  if(!p->x)
	    distancesInitial(p, distances, tr, FALSE);	
	  
	  right = distances[p->number];	  
	}
      else
	right = 0.0;
                 
      if(ABS(left - right) <= p->z[0])
	{	 
	  *rootBranch = p;
	  *minimum = 0.0;
	}
      else
	{
	  double
	    diff;

	  if(left > right)	   
	    diff = left - (right + p->z[0]);	    	    
	  else	    
	    diff = right - (left + p->z[0]);	  

	  if(*minimum > diff)
	    {	      
	      *minimum = diff;
	      *rootBranch = p;
	    }
	}

      q = p->next;

      while(q != p)
	{
	  distancesNewview(q->back, distances, tr, rootBranch, minimum);
	  q = q->next;
	}      
    }
}

static void printTreeRec(FILE *f, nodeptr p, tree *tr, boolean rootDescendant, boolean printBranchLabels)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(rootDescendant)
	fprintf(f, "%s", tr->nameList[p->number]);
      else
	fprintf(f, "%s:%f", tr->nameList[p->number], p->z[0]);
    }
  else
    {
      fprintf(f, "(");
      printTreeRec(f, p->next->back, tr, FALSE, printBranchLabels);
      fprintf(f, ",");
      printTreeRec(f, p->next->next->back, tr, FALSE, printBranchLabels);
      
      if(rootDescendant)
	fprintf(f, ")");
      else
	{
	  if(printBranchLabels && !isTip(p->number, tr->mxtips) && !isTip(p->back->number, tr->mxtips))
	    {
	      assert(p->support == p->back->support);
	      fprintf(f, "):%f[%d]", p->z[0], p->support);
	    }
	  else
	    fprintf(f, "):%f", p->z[0]);
	}
    }
}

static void printTree(nodeptr p, tree *tr, double *distances, FILE *f, boolean printBranchLabels)
{
  double
    leftRoot,
    rightRoot,  
    thisBranch = p->z[0],   
    left = 0.0,
    right = 0.0;

  nodeptr 
    q = p->back;

  if(!isTip(p->number, tr->mxtips))
    {
      if(!p->x)
	distancesInitial(p, distances, tr, FALSE);	
	  
      left = distances[p->number];	  
    }
  else
    left = 0.0;

   if(!isTip(q->number, tr->mxtips))
    {
      if(!q->x)
	distancesInitial(q, distances, tr, FALSE);	
	  
      right = distances[q->number];	  
    }
   else
     left = 0.0;

   //printf("left %f right %f thisBranch %f\n", left, right, thisBranch);

   if(ABS(left - right) <= thisBranch)
     {      
       if(left < right)
	 {
	   leftRoot = (right + thisBranch - left) / 2.0;
	   rightRoot = thisBranch - leftRoot;
	 }
       else
	 {
	   rightRoot = (left + thisBranch - right) / 2.0;
	   leftRoot = thisBranch - rightRoot;
	 }	                            
     }
   else
     {
       if(left < right)
	 {
	   leftRoot  = thisBranch;
	   rightRoot = 0.0;
	 }
       else
	 {
	   leftRoot  = 0.0;
	   rightRoot = thisBranch;	   
	 }
     }

   //descend into right subtree and print it

   fprintf(f, "(");
   printTreeRec(f, p, tr, TRUE, printBranchLabels);

   //finished right subtree, print attachment branch of right subtree
   //noew descent into left subtree
   
   if(printBranchLabels && !isTip(p->number, tr->mxtips) && !isTip(q->number, tr->mxtips))
     {
       assert(p->support == q->support);
       fprintf(f, ":%f[%d], ", leftRoot, p->support);
     }
   else
     fprintf(f, ":%f, ", leftRoot);
   printTreeRec(f, q, tr, TRUE, printBranchLabels);
   
   //finished left subtree, now print its branch to the root node 
   //and we are done 

   if(printBranchLabels && !isTip(p->number, tr->mxtips) && !isTip(q->number, tr->mxtips))
     {
       assert(p->support == q->support);
       fprintf(f, ":%f[%d]);", rightRoot, q->support);
     }
   else
     fprintf(f, ":%f);", rightRoot);
}

static void rootTree(tree *tr, analdef *adef)
{
  int 
    i;

  double
    checkDistances,
    minimum,
    *distances = (double *)rax_malloc(sizeof(double) * 2 * tr->mxtips);
  
  char 
    rootedTreeFile[1024];

  FILE 
    *f = myfopen(tree_file, "r");
   
  nodeptr 
    rootBranch;

  boolean 
    printBranchLabels = FALSE;

  for(i = 0; i < 2 * tr->mxtips; i++)
    distances[i] = 0.0;
  
  strcpy(rootedTreeFile,         workdir);
  strcat(rootedTreeFile,         "RAxML_rootedTree.");
  strcat(rootedTreeFile,         run_id);

  treeReadLen(f, tr, TRUE, FALSE, TRUE, adef, TRUE, TRUE);

  if(tr->branchLabelCounter > 0)
    {
      assert(tr->branchLabelCounter == (tr->ntips - 3));
      printBranchLabels = TRUE;
      printBothOpen("\nYour input tree contains branch labels, these will also be printed in the output tree ...\n\n");
    }

  fclose(f);

  minimum = checkDistances = distancesInitial(tr->start->back, distances, tr, TRUE);

  //printf("Tree Lenght: %f\n", checkDistances); 
  
  f = myfopen(rootedTreeFile, "w");

  distancesNewview(tr->start->back, distances, tr, &rootBranch, &minimum);

  printTree(rootBranch, tr, distances, f, printBranchLabels);
  fprintf(f, "\n");
  
  fclose(f);

  printBothOpen("RAxML-rooted tree using subtree length-balance printed to file:\n%s\n",  rootedTreeFile);

  rax_free(distances);
}

static boolean partitionHasInvariantSites(tree *tr, int model)
{
  unsigned int   
    patternCount = 0,
    i;

  int
    j,
    *patternIndices = (int *)rax_malloc(sizeof(int) * (tr->partitionData[model].upper - tr->partitionData[model].lower));

  const unsigned int 
    *bitVector = getBitVector(tr->partitionData[model].dataType),
    undetermined = getUndetermined(tr->partitionData[model].dataType);

  for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)	      
    {		 
      unsigned int 
	encoding = undetermined;            	  	

      for(j = 1; j <= tr->mxtips; j++)		
	encoding = encoding & bitVector[tr->yVector[j][i]];       

      if(encoding > 0)    	
	patternIndices[patternCount++] = i;
    }

  if(patternCount > 0)
    {     	
      printBothOpen("Partition %d with name \"%s\" is to be analyzed using ascertainment bias correction, but it has %u invariable site(s)!\n", 
		    model, tr->partitionData[model].partitionName, patternCount);
      printBothOpen("This is is not allowed! RAxML will print the offending site(s) and then exit.\n\n");

      for(i = 0; i < patternCount; i++)	
	{
	  int 
	    k;
	  
	  printBothOpen("Pattern: ");
	  for(j = 1; j <= tr->mxtips; j++)
	    printBothOpen("%c", getInverseMeaning(tr->partitionData[model].dataType, tr->yVector[j][patternIndices[i]]));
	  printBothOpen("\n");

	  printBothOpen("Pattern occurs at the following sites of the input alignment: \n");
	    
	  for(k = 0; k < tr->rdta->sites; k++)	    
	    if(patternIndices[i] == tr->patternPosition[k])
	      printBothOpen("Site %d \n", tr->columnPosition[k]); 
	  
	  printBothOpen("\n");
	}	 	
	     	   
      rax_free(patternIndices);

      return TRUE;	        
    }
  else
    {
      rax_free(patternIndices);

      return FALSE;
    }
}

static void checkAscBias(tree *tr)
{
  int 
    model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      if(tr->partitionData[model].ascBias)
	{
	  boolean
	    hasInvariantSites = partitionHasInvariantSites(tr, model);
      
	  if(hasInvariantSites)
	    {
	      printBothOpen("\n\nFor partition %s you specified that the likelihood score shall be corrected for invariant sites\n", tr->partitionData[model].partitionName);
	      printBothOpen("via an ascertainment bias correction. However, some sites in this partition are already invariant.\n");
	      printBothOpen("This is not allowed, please remove all invariant sites and try again, exiting ... \n\n\n");
	      errorExit(-1);
	    }	  
	}
    }

  if(tr->rateHetModel == CAT)
    for(model = 0; model < tr->NumberOfModels; model++)
      {
	if(tr->partitionData[model].ascBias && !tr->noRateHet)
	  {
	    printBothOpen("\nWARNING: you specified that you want to use an ascertainment bias correction for partition %d\n", model);
	    printBothOpen("for the CAT model of rate heterogeneity. Are you sure that you don't want to use a model without any rate \n");
	    printBothOpen("heterogeneity modeling via the \"-V\" command line switch?\n");
	  }
      }
}


static void readAscFiles(tree *tr)
 {     
   if(tr->ascertainmentCorrectionType == STAMATAKIS_CORRECTION || tr->ascertainmentCorrectionType == FELSENSTEIN_CORRECTION || tr->ascertainmentCorrectionType == GOLDMAN_CORRECTION_3)
     {
       int 
	 model;
       
       for(model = 0; model < tr->NumberOfModels; model++)
	 {
	   if(tr->partitionData[model].ascBias)
	     {
	       if(filexists(tr->partitionData[model].ascFileName))
		 {
		   FILE 
		     *f = myfopen(tr->partitionData[model].ascFileName, "rb");
		   
		   switch(tr->ascertainmentCorrectionType)
		     {
		     case STAMATAKIS_CORRECTION:
		       {
			 boolean
			   errorDetected = FALSE;
			 
			 unsigned int
			   *weights = (unsigned int*)rax_malloc((size_t)tr->partitionData[model].states * sizeof(unsigned int));
			 
			 int
			   i;
			 
			 for(i = 0; i < tr->partitionData[model].states; i++)
			   if(fscanf(f, "%u", &weights[i]) != 1)
			     {
			       errorDetected = TRUE;
			       break;
			     }
			 
			 if(errorDetected)
			   {
			     printf("\nProblem reading number of invariable sites in ascertainment correction file %s\n", tr->partitionData[model].ascFileName);
			     printf("for stamatakis ascertainment bias correction for partition %d with name %s.\n\n", model, tr->partitionData[model].partitionName);
			     errorExit(-1);
			   }
			 
			 for(i = 0; i < tr->partitionData[model].states; i++)
			   tr->partitionData[model].invariableFrequencies[i] = (double)weights[i];
			 
			 rax_free(weights);
		       }
		       break;
		     case FELSENSTEIN_CORRECTION:
		     case GOLDMAN_CORRECTION_3:
		       {
			 unsigned int 
			   length;
			 
			 if(fscanf(f, "%u", &length) != 1)
			   {
			     printf("\nProblem reading number of invariable sites in ascertainment correction file %s\n", tr->partitionData[model].ascFileName);
			     printf("for felsenstein ascertainment bias correction for partition %d with name %s.\n\n", model, tr->partitionData[model].partitionName);
			     errorExit(-1);
			   }
			 
			 tr->partitionData[0].invariableWeight = (double)length;
		       }
		       break;
		     default:
		       assert(0);
		     }
		   
		   fclose(f);
		 }
	       else
		 {
		   printf("You specified that you want to use a stamatakis or felsenstein ascertainment bias correction for partition %d with name %s.\n", 
			  model, tr->partitionData[model].partitionName);
		   printf("but did not specify a correction file for this partition in the partition file!\n\n");
		   errorExit(-1);
		 }
	     }	      			      
	 }
     }
 }
      
     

/******* branch length stealing code -f k option **************************************************/


/* function to determine if tip with tip index number for model/partition model 
   only consists of missing data */

static boolean tipHasData(tree *tr, int model, int number)
{
  unsigned char 	
    undetermined = getUndetermined(tr->partitionData[model].dataType),
    *tip = tr->partitionData[model].yVector[number];

  size_t
    i;

  boolean
    data = FALSE;
  
  for(i = 0; i < tr->partitionData[model].width && (!data); i++)
    if(tip[i] != undetermined)
      data = TRUE; 

  return data;
}

/* recursive function to determine if the subtree rooted at p has taxa with sequence data for partition model 
   in at least one of the tips */

static boolean hasData(tree *tr, nodeptr p, int model)
{
  if(isTip(p->number, tr->mxtips))         
    return tipHasData(tr, model, p->number);   
  else   
    return (hasData(tr, p->next->back, model) || hasData(tr, p->next->next->back, model));   
}
  

/* recursive function to adapt branch lengths of subtrees of partitions with missing data using 
   branch length info from those partitions that have data ! */

static void adaptBranchLengths(tree *tr, nodeptr p, int *count)
{  
  int 
    *missingData = (int *)rax_calloc(tr->NumberOfModels, sizeof(int)),
    wgtsum = 0,
    model,
    partitionsWithData = 0,
    partitionsWithoutData = 0;

  double
    branchLength = 0.0;
    
  /* first we check if the branch defined by p and p->back 
     has data on both sides (at least one tip with data
     for all partitions.

     When this is not the case we need to "steal" an approximate branch length from 
     one of the other partitions.
  */

  //increment the number of branches we have visited
  *count = *count + 1;

  //compute the number of sites per partition for which we have data on both sides of the branch

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(hasData(tr, p, model) && hasData(tr, p->back, model))
	{
	  wgtsum += (tr->partitionData[model].upper - tr->partitionData[model].lower);
	  partitionsWithData++;
	  missingData[model] = 0;
	}  
      else
	{
	  missingData[model] = 1;
	  partitionsWithoutData++;
	}
    }

  //make sure that there is at least one partition from which we can steal the branch length 
  assert(partitionsWithData > 0);
  

  //now compute the branch length average over all partitions 
  //that have data on both sides of this branch 
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(missingData[model] == 0 &&  partitionsWithoutData > 0)
	{
	  double 
	    factor = (double)(tr->partitionData[model].upper - tr->partitionData[model].lower) / (double)wgtsum,
	    z = p->z[model];	      
	  
	  //printf("factor %f\n", factor);

	  if(z < zmin)
	    z = zmin;
	  
	  z = -log(z);
	  
	  branchLength += z * factor;

	  //printf("br-len: %f %f\n", p->z[model], branchLength);
	}
    }     
  
  //and assign this average branch length to the partitions that don't have data present across this branch

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(missingData[model] == 1)
	{
	  double
	    targetBranch = exp(-(branchLength));

	  //printf("adapted one branch in part %d %1.40f -> %1.40f\n", model, p->z[model], targetBranch);

	  if(targetBranch < zmin)
	    targetBranch = zmin;

	  if(targetBranch > zmax)
	    targetBranch = zmax;

	  p->z[model] = targetBranch;
	  p->back->z[model] = targetBranch;	  
	}  
    }

  rax_free(missingData);
        
  //now handle the remaining branches recursively
  
  if(!isTip(p->number, tr->mxtips))
    {
      adaptBranchLengths(tr, p->next->back, count);
      adaptBranchLengths(tr, p->next->next->back, count);      
    }

}

static void stealBranchLengths(tree *tr, analdef *adef)
{
  char 
    fileName[1024];

  FILE
    *f;
  
  int
    model,
    count = 0;

  double
    wgtsum = 0.0,
    treeLengthBefore = 0.0,
    treeLengthAfter = 0.0;
  
  
  //-M must be set in the command line and we need to use a partition file

  assert(tr->numBranches == tr->NumberOfModels && tr->NumberOfModels > 0);

  //initially optimize the model parameters of the given tree 

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
  printBothOpen("After model optimization on the tree: %f with %d taxa\n", tr->likelihood, tr->mxtips);

  assert(!isTip(tr->start->back->number, tr->mxtips));            

  //get the initial likelihood   

  evaluateGenericInitrav(tr, tr->start);

  for(model = 0; model < tr->NumberOfModels; model++)   
    wgtsum += (double)(tr->partitionData[model].upper - tr->partitionData[model].lower);
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      double 
	factor = (double)(tr->partitionData[model].upper - tr->partitionData[model].lower) / wgtsum;           

      assert(factor == tr->partitionContributions[model]);
      
      //printf("%d %f factor %f\n", model, treeLength(tr, model), factor);

      treeLengthBefore += treeLength(tr, model) * factor;
    }

  printBothOpen("Likelihood before br-len stealing %f\n\n", tr->likelihood);
  printBothOpen("Tree length before br-len stealing %f\n\n", treeLengthBefore);
  
  //Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);

  printf("%s \n", tr->tree_string);
  
  //now adapt branch lengths
  adaptBranchLengths(tr, tr->start->back, &count);
  
  //make sure that we visited all branches
  assert(count == 2 * tr->mxtips - 3);


  //re-calculate likelihood, note that the branch length adaptation doesn't alter 
  //the likelihood (or shouldn't if this is correctly implemented) since we only changed
  //branch lengths for those parts of the tree that have missing data
  evaluateGenericInitrav(tr, tr->start);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      double 
	factor = (double)(tr->partitionData[model].upper - tr->partitionData[model].lower) / wgtsum;           

      treeLengthAfter += treeLength(tr, model) * factor;
    }
  
  printBothOpen("Likelihood after br-len stealing %f\n\n", tr->likelihood);
  printBothOpen("Tree length after br-len stealing %f\n\n", treeLengthAfter);
  
  strcpy(fileName, workdir);
  strcat(fileName,             "RAxML_stolenBranchLengths.");
  strcat(fileName,            run_id);

  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE, FALSE, FALSE);

  f = myfopen(fileName, "wb");

  fprintf(f, "%s", tr->tree_string);
  
  fclose(f);

  printBothOpen("Tree with stolen branch lengths written to file %s, giasou file mou.\n\n", fileName);

  exit(0);
}






/*******************************************************/


int main (int argc, char *argv[])
{
#if (defined(_WAYNE_MPI) || defined (_QUARTET_MPI))
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);
  MPI_Comm_size(MPI_COMM_WORLD, &processes);
  printf("\nThis is RAxML MPI Process Number: %d\n", processID);
#else
  processID = 0;
#endif
  {
    rawdata      *rdta;
    cruncheddata *cdta;
    tree         *tr;
    analdef      *adef;
    int
      i,
      countGTR = 0,
      countOtherModel = 0,
      countAscBias = 0;
    
#if (defined(_USE_PTHREADS) && !defined(_PORTABLE_PTHREADS))  
    pinToCore(0);
#endif 
    
    

    masterTime = gettime();

    globalArgc = argc;
    globalArgv = (char **)rax_malloc(sizeof(char *) * argc);
    for(i = 0; i < argc; i++)
      globalArgv[i] = argv[i];
    
    
    
#if ! (defined(__ppc) || defined(__powerpc__) || defined(PPC))
    
    /* 
       David Defour's command  
       _mm_setcsr( _mm_getcsr() | (_MM_FLUSH_ZERO_ON | MM_DAZ_ON));  
    */
    
    _mm_setcsr( _mm_getcsr() | _MM_FLUSH_ZERO_ON);
    
#endif 
    
    adef = (analdef *)rax_malloc(sizeof(analdef));
    rdta = (rawdata *)rax_malloc(sizeof(rawdata));
    cdta = (cruncheddata *)rax_malloc(sizeof(cruncheddata));
    tr   = (tree *)rax_malloc(sizeof(tree));
    
    /* initialize lookup table for fast bit counter */
    
    compute_bits_in_16bits();
    
    initAdef(adef);
    get_args(argc,argv, adef, tr); 
    
    
    if(adef->readTaxaOnly)  
      {
	if(adef->mode == PLAUSIBILITY_CHECKER || adef->mode == ROOT_TREE || adef->mode == CALC_BIPARTITIONS_IC)
	  extractTaxaFromTopology(tr, rdta, cdta, tree_file);   
	else
	  extractTaxaFromTopology(tr, rdta, cdta, bootStrapFile);
      }
    
    getinput(adef, rdta, cdta, tr);
    
    

    checkOutgroups(tr, adef);
    makeFileNames();
    
#if (defined(_WAYNE_MPI) || defined (_QUARTET_MPI))
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
    
    
    if(!adef->noSequenceCheck && !adef->readTaxaOnly && adef->mode != FAST_SEARCH && adef->mode != SH_LIKE_SUPPORTS)
      checkSequences(tr, rdta, adef);
    else
      printBothOpen("\n\nWARNING: RAxML is not checking sequences for duplicate seqs and sites with missing data!\n\n");
    
    if(tr->NumberOfModels == 1 && tr->partitionData[0].dataType == DNA_DATA && adef->useBFGS)
      printBothOpen("\n\nUsing BFGS method to optimize GTR rate parameters, to disable this specify \"--no-bfgs\" \n\n");

    
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
	/* might as well rax_free the initial structures here */
	
      }
  
    if(!adef->readTaxaOnly)
      {
	int        
	  countNonSev = 0,	 
	  countLG4 =0;
	
	assert(countAscBias == 0);

	for(i = 0; i < tr->NumberOfModels; i++)	  
	  if(tr->partitionData[i].ascBias)
	    countAscBias++;

	makeweights(adef, rdta, cdta, tr, countAscBias);
	makevalues(rdta, cdta, tr, adef);      
	
	for(i = 0; i < tr->NumberOfModels; i++)
	  {	    
	    if(!(tr->partitionData[i].dataType == AA_DATA || tr->partitionData[i].dataType == DNA_DATA))
	      countNonSev++;	    
	       
	    
	    if(tr->partitionData[i].protModels == LG4 || tr->partitionData[i].protModels == LG4X)
	      {
		countLG4++;

		if(tr->partitionData[i].optimizeBaseFrequencies)
		  {
		    FILE 
		      *info = myfopen(infoFileName, "ab");
	    
		    printBoth(info, "\n\nYou are using the %s model of AA substitution for partition %d.\n", (tr->partitionData[i].protModels == LG4)?"LG4M":"LG4X", i);
		    printBoth(info, "and specified that you want to use a joint ML estimate of the base frequencies shared across all 4\n");
		    printBoth(info, "substitution matrices of the LG4 model.\n");
		    printBoth(info, "WARNING: This does not correspond to the idea behind the LG4 model!\n");
		    printBoth(info, "         which has different sets of base ferquencies for each substitution matrix.\n");
		    printBoth(info, "         The ML estimate of base freqs with %s was simply implemented for convenience in RAxML.\n", (tr->partitionData[i].protModels == LG4)?"LG4M":"LG4X");
		    printBoth(info, "Olivier Gascuel asked me to tell you that he doesn't like it to be used like this!\n");
		    printBoth(info, "Since he is a nice guy you should use the original model %s instead.\n\n\n", (tr->partitionData[i].protModels == LG4)?"LG4M":"LG4X");
	    
		    fclose(info);		    
		  }
		else
		  if(!tr->partitionData[i].usePredefinedProtFreqs)
		    {
		      FILE 
			*info = myfopen(infoFileName, "ab");
	    
		      printBoth(info, "\n\nYou are using the %s model of AA substitution for partition %d.\n", (tr->partitionData[i].protModels == LG4)?"LG4M":"LG4X", i);
		      printBoth(info, "and specified that you want to use joint empirical base frequencies (drawn from the alignment) shared across all 4\n");
		      printBoth(info, "substitution matrices of the LG4 model.\n");
		      printBoth(info, "WARNING: This does not correspond to the idea behind the LG4 model!\n");
		      printBoth(info, "         which has different sets of base ferquencies for each substitution matrix.\n");
		      printBoth(info, "         The option to use empirical base freqs with %s was simply implemented for convenience in RAxML.\n", (tr->partitionData[i].protModels == LG4)?"LG4M":"LG4X");
		      printBoth(info, "Olivier Gascuel asked me to tell you that he doesn't like it to be used like this!\n");
		      printBoth(info, "Since he is a nice guy you'd want to use the original model %s instead.\n\n\n", (tr->partitionData[i].protModels == LG4)?"LG4M":"LG4X");
	    
		      fclose(info);
		    }
	      }
	    
	    if(tr->partitionData[i].dataType == AA_DATA)
	      {
		if(tr->partitionData[i].protModels == GTR || tr->partitionData[i].protModels == GTR_UNLINKED)
		  countGTR++;
		else
		  countOtherModel++;
	      }
	  }
	
	if(countLG4 > 0)
	  {
	    if(tr->saveMemory)
	      {
		printf("Error: the LG4 substitution model does not work in combination with the \"-U\" memory saving flag!\n\n");	  
		errorExit(-1);
	      }
	    
	    if(adef->useInvariant)
	      {
		printf("Error: the LG4 substitution model does not work for proportion of invariavble sites estimates!\n\n");	  
		errorExit(-1);
	      }
	    
	    if(isCat(adef))
	      {
		printf("Error: the LG4 substitution model does not work with the CAT model of rate heterogeneity!\n\n");	  
		errorExit(-1);	    
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

    if(tr->rateHetModel == GAMMA_I && tr->saveMemory)
      {
	 printf("\nError: Memory saving option \"-U\" not implemented for models with proportion\n");
	 printf("of variable site estimates, reomve \"-U\" from your command line and re-run, exiting now.\n\n");
	 errorExit(-1);
      }

    if(countAscBias > 0 && !adef->readTaxaOnly)
      {
	if(tr->ascertainmentCorrectionType == NOT_DEFINED)
	  {
	    printf("\nError, one or more of your partitions use an ascertainment bias correction,\n");
	    printf("but you have not specified which correction type you want to use via \"--asc-corr\" \n\n");
	    errorExit(-1);
	  }

	checkAscBias(tr);				      					
      }

#ifdef _USE_PTHREADS
    startPthreads(tr, adef);
    masterBarrier(THREAD_INIT_PARTITION, tr);
    if(!adef->readTaxaOnly)  
      masterBarrier(THREAD_ALLOC_LIKELIHOOD, tr);
#else
    if(!adef->readTaxaOnly)  
      allocNodex(tr);
#endif
    
    
    readAscFiles(tr);

    if(!adef->readTaxaOnly) 
      {
#ifdef _USE_PTHREADS
	masterBarrier(THREAD_SETUP_PRESENCE_MAP, tr);
#else     
	setupPresenceMask(tr);
#endif
      }
    
    printModelAndProgramInfo(tr, adef, argc, argv);
      
#ifdef _BASTIEN
    if(adef->mode != TREE_EVALUATION)
      assert(0);
#endif


      switch(adef->mode)
	{  
	case SUBTREE_EPA:
	  if(adef->useBinaryModelFile)      
	    readBinaryModel(tr, adef);	       
	  else
	    initModel(tr, rdta, cdta, adef);
	  
	  getStartingTree(tr, adef);

	  subtreeEPA(tr, adef);
	  assert(0);
	  break;
	case CLASSIFY_MP:
	  getStartingTree(tr, adef);
	  assert(0);
	  break;
	case CLASSIFY_ML:
	  if(adef->useBinaryModelFile)      
	    readBinaryModel(tr, adef);	       
	  else
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
	    computeConsensusOnly(tr, bootStrapFile, adef, adef->calculateIC);
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
	case STEAL_BRANCH_LENGTHS:
	  initModel(tr, rdta, cdta, adef);      
	  getStartingTree(tr, adef); 
	  stealBranchLengths(tr, adef);      
	  break;
	case TREE_EVALUATION:
	  initModel(tr, rdta, cdta, adef);
	  
	  getStartingTree(tr, adef);      
	  
	  if(adef->likelihoodTest)
	    computeLHTest(tr, adef, bootStrapFile);
	  else
	    { 
	      if(adef->useBinaryModelFile)	 	 
		{
		  readBinaryModel(tr, adef);
		  evaluateGenericInitrav(tr, tr->start);	      
		  treeEvaluate(tr, 2);
		}
	      else
		{	      
		  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);	  
#ifdef _BASTIEN
		  printf("likelihood of current tree: %f\n", tr->likelihood);

		  tr->doBastienStuff = TRUE;
		  printf("do a full traversal\n");
		  evaluateGenericInitrav(tr, tr->start);
		  exit(0);
#endif
		  writeBinaryModel(tr, adef);
		}
	      
	      printLog(tr, adef, TRUE);
	      printResult(tr, adef, TRUE);
	    }
	  
	  break;
	case ANCESTRAL_STATES:
	  initModel(tr, rdta, cdta, adef);
	  
	  getStartingTree(tr, adef);
	  
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
	  
	  evaluateGenericInitrav(tr, tr->start);                                       	                  
	  
	  computeAncestralStates(tr, tr->likelihood);
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
	case CALC_BIPARTITIONS_IC:
	  calcBipartitions_IC_Global(tr, adef, tree_file, bootStrapFile);
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
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
	  morphologicalCalibration(tr, adef);
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
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon);
	  computePlacementBias(tr, adef);
	  break;
	case OPTIMIZE_BR_LEN_SCALER:
	  initModel(tr, rdta, cdta, adef);
	  
	  getStartingTree(tr, adef);      
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);	  
	  
	  printBothOpen("Likelihood: %f\n", tr->likelihood);
	  
	  break;
	case ANCESTRAL_SEQUENCE_TEST:
	  initModel(tr, rdta, cdta, adef);
	  
	  getStartingTree(tr, adef);  
	  
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, adef->likelihoodEpsilon);	
	  
	  ancestralSequenceTest(tr);
	  break;
	case PLAUSIBILITY_CHECKER:
	  plausibilityChecker(tr, adef);
	  exit(0);
	  break;
	case ROOT_TREE:
	  rootTree(tr, adef);    
	  break;
	default:
	  assert(0);
	}
      
      finalizeInfoFile(tr, adef);

#if (defined(_WAYNE_MPI) || defined (_QUARTET_MPI))
      MPI_Finalize();
#endif
  }

  return 0;
}


