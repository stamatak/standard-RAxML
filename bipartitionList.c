/*  RAxML-HPC, a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright March 2006 by Alexandros Stamatakis
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
 *
 *  When publishing work that is based on the results from RAxML please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#ifndef WIN32  
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include "axml.h"
#include "rmq.h" //include range minimum queries for fast plausibility checker

#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>

#endif

#ifdef _USE_PTHREADS
#include <pthread.h>
#endif

#ifdef _WAYNE_MPI
#include <mpi.h>
extern int processID;
extern int processes;
#endif


extern FILE *INFILE;
extern char run_id[128];
extern char workdir[1024];
extern char bootStrapFile[1024];
extern char tree_file[1024];
extern char infoFileName[1024];
extern char resultFileName[1024];
extern char verboseSplitsFileName[1024];
extern char bipartitionsFileNameBranchLabels[1024];
extern char icFileNameBranchLabelsStochastic[1024];
extern char icFileNameBranchLabelsUniform[1024];
extern char icFileNameBranchLabels[1024];

extern double masterTime;

extern const unsigned int mask32[32];

extern volatile branchInfo      **branchInfos;
extern volatile int NumberOfThreads;
extern volatile int NumberOfJobs;

static void checkTreeNumber(int numberOfTrees, char *fileName)
{
  if(numberOfTrees <= 1)
    {
      printf("RAxML is expecting to read more than one tree in file %s for this operation on a set of trees!\n", fileName);
      printf("The program will exit now\n");      
      exit(-1);
    }
}

static void mre(hashtable *h, boolean icp, entry*** sbi, int* len, int which, int n, unsigned int vectorLength, boolean sortp, tree *tr, boolean bootStopping);


entry *initEntry(void)
{
  entry *e = (entry*)rax_malloc(sizeof(entry));

  e->bitVector     = (unsigned int*)NULL;
  e->treeVector    = (unsigned int*)NULL;
  e->supportVector = (int*)NULL;
  e->bipNumber  = 0;
  e->bipNumber2 = 0;
  e->supportFromTreeset[0] = 0;
  e->supportFromTreeset[1] = 0;
  e->next       = (entry*)NULL;

  //Partial IC/TC modification
  e->taxonMask     = (unsigned int*)NULL;
  e->adjustedSupport =  0;
  e->tempSupport     =  0;
  e->tempSupportFrom = -1;
  e->coveredNumber   =  0;
  e->covered         = FALSE;
  //Partial IC/TC end

  return e;
} 




hashtable *initHashTable(hashNumberType n)
{
  /* 
     init with primes 
    
     static const hashNumberType initTable[] = {53, 97, 193, 389, 769, 1543, 3079, 6151, 12289, 24593, 49157, 98317,
     196613, 393241, 786433, 1572869, 3145739, 6291469, 12582917, 25165843,
     50331653, 100663319, 201326611, 402653189, 805306457, 1610612741};
  */

  /* init with powers of two */

  static const  hashNumberType initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384,
					      32768, 65536, 131072, 262144, 524288, 1048576, 2097152,
					      4194304, 8388608, 16777216, 33554432, 67108864, 134217728,
					      268435456, 536870912, 1073741824, 2147483648U};
  
  hashtable *h = (hashtable*)rax_malloc(sizeof(hashtable));
  
  hashNumberType
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (hashNumberType)-1;    

  assert(n <= maxSize);

  i = 0;

  while(initTable[i] < n && i < primeTableLength)
    i++;

  assert(i < primeTableLength);

  tableSize = initTable[i];

  /* printf("Hash table init with size %u\n", tableSize); */

  h->table = (entry**)rax_calloc(tableSize, sizeof(entry*));
  h->tableSize = tableSize;  
  h->entryCount = 0;  

  return h;
}




void freeHashTable(hashtable *h)
{
  hashNumberType
    i,
    entryCount = 0;
   

  for(i = 0; i < h->tableSize; i++)
    {
      if(h->table[i] != NULL)
	{
	  entry *e = h->table[i];
	  entry *previous;	 

	  do
	    {
	      previous = e;
	      e = e->next;

	      if(previous->bitVector)
		rax_free(previous->bitVector);

	      //Partial IC/TC modif 

	      if(previous->taxonMask)
		rax_free(previous->taxonMask);
	      
	      //Partial IC/TC modif end 

	      if(previous->treeVector)
		rax_free(previous->treeVector);

	      if(previous->supportVector)
		rax_free(previous->supportVector);
	      
	      rax_free(previous);	      
	      entryCount++;
	    }
	  while(e != NULL);	  
	}

    }

  assert(entryCount == h->entryCount);
 
  rax_free(h->table);
}



void cleanupHashTable(hashtable *h, int state)
{
  hashNumberType
    k,
    entryCount = 0,
    removeCount = 0;
 
  assert(state == 1 || state == 0);

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];
	  entry *start     = (entry*)NULL;
	  entry *lastValid = (entry*)NULL;
	  	  
	  do
	    {	   	 	      	
	      if(state == 0)
		{
		  e->treeVector[0] = e->treeVector[0] & 2;	
		  assert(!(e->treeVector[0] & 1));
		}
	      else
		{
		  e->treeVector[0] = e->treeVector[0] & 1;
		  assert(!(e->treeVector[0] & 2));
		}
	      
	      if(e->treeVector[0] != 0)
		{
		  if(!start)
		    start = e;
		  lastValid = e;
		  e = e->next;
		}	  
	      else
		{
		  entry *remove = e;
		  e = e->next;
		  
		  removeCount++;

		  if(lastValid)		    		    
		    lastValid->next = remove->next;

		  if(remove->bitVector)
		    rax_free(remove->bitVector);
		  
		  if(remove->treeVector)
		    rax_free(remove->treeVector);
		  
		  if(remove->supportVector)
		    rax_free(remove->supportVector);

		  //Partial IC/TC

		  if(remove->taxonMask)
		    rax_free(remove->taxonMask);

		  //Partial IC/TC end

		  rax_free(remove);		 
		}
	      
	      entryCount++;	     	     
	    }
	  while(e != NULL);	 

	  if(!start)
	    {
	      assert(!lastValid);
	      h->table[k] = NULL;
	    }
	  else
	    {
	      h->table[k] = start;
	    }	 	 
	}    
    }

  assert(entryCount ==  h->entryCount);  

  h->entryCount -= removeCount;
}

unsigned int **initBitVector(tree *tr, unsigned int *vectorLength)
{
  unsigned int **bitVectors = (unsigned int **)rax_malloc(sizeof(unsigned int*) * 2 * tr->mxtips);
  int i;

  if(tr->mxtips % MASK_LENGTH == 0)
    *vectorLength = tr->mxtips / MASK_LENGTH;
  else
    *vectorLength = 1 + (tr->mxtips / MASK_LENGTH); 
  
  for(i = 1; i <= tr->mxtips; i++)
    {
      bitVectors[i] = (unsigned int *)rax_calloc(*vectorLength, sizeof(unsigned int));
      bitVectors[i][(i - 1) / MASK_LENGTH] |= mask32[(i - 1) % MASK_LENGTH];
    }
  
  for(i = tr->mxtips + 1; i < 2 * tr->mxtips; i++) 
    bitVectors[i] = (unsigned int *)rax_malloc(sizeof(unsigned int) * *vectorLength);

  return bitVectors;
}

void freeBitVectors(unsigned int **v, int n)
{
  int i;

  for(i = 1; i < n; i++)
    rax_free(v[i]);
}







/* compute bit-vectors representing bipartitions/splits for a multi-furcating tree */

static void newviewBipartitionsMultifurcating(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength)
{
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q,
      firstDescendant;
    
    unsigned int       
      *vector = bitVectors[p->number];
   
    unsigned 
      int i;           
    
    int
      x_set = 0,
      number;
    
    /* Set the directional x token to the correct element in the 
       cyclically linked list representing an inner node */

    q = p->next;
    
    if(p->x)
      x_set++;

    p->x = 1;   
    number = p->number;

    while(q != p)
      {
	if(q->x) 
	  x_set++;
	q->x = 0;
	assert(q->number == number);
	q = q->next;
      }
   
    assert(x_set == 1);

    /* get the first connecting branch of the node to initialize the 
       bipartition hash value and the bit vector representing this split.
     */

    firstDescendant = p->next->back;      
    
    /* if this is not a tip, we first need to recursively 
       compute the hash values and bit-vectors for the subtree rooted at 
       firstDescendant */

    if(!isTip(firstDescendant->number, numsp))
      {	
	if(!firstDescendant->x)
	  newviewBipartitionsMultifurcating(bitVectors, firstDescendant, numsp, vectorLength);
      }
	
    /* initialize the bit vector of the current split by the bit vector of the first descandant */
    
    for(i = 0; i < vectorLength; i++)
      vector[i] = bitVectors[firstDescendant->number][i];
	
    /* initialize the hash key of the current split by the hash key of the first descandant */

    p->hash = firstDescendant->hash;
    
    /* handle all other descendants of this inner node potentially representing a multi-furcation */

    q = p->next->next;

    while(q != p)
      {
	/* update the has key by xoring with the current hash with the hash of this descendant */
	
	p->hash = p->hash ^ q->back->hash;
	
	if(!isTip(q->back->number, numsp))
	  {	   
	    if(!q->back->x)
	      newviewBipartitionsMultifurcating(bitVectors, q->back, numsp, vectorLength);
	  }
	    	
	/* update the bit-vector representing the current split by applying a bitwise with the bit vector of the descendant */
	
	for(i = 0; i < vectorLength; i++)
	  vector[i] = bitVectors[q->back->number][i] | vector[i];	    	 
	
	q = q->next;
      }     
  }     
}


static void insertHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int bipNumber, hashNumberType position)
{
  entry *e = initEntry();

  e->bipNumber = bipNumber; 
  /*e->bitVector = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int)); */

  e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
  memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));
 
  memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
  
  if(h->table[position] != NULL)
    {
      e->next = h->table[position];
      h->table[position] = e;           
    }
  else
    h->table[position] = e;

  h->entryCount =  h->entryCount + 1;
}



static int countHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position)
{ 
  if(h->table[position] == NULL)         
    return -1;
  {
    entry *e = h->table[position];     

    do
      {	 
	unsigned int i;

	for(i = 0; i < vectorLength; i++)
	  if(bitVector[i] != e->bitVector[i])
	    goto NEXT;
	   
	return (e->bipNumber);	 
      NEXT:
	e = e->next;
      }
    while(e != (entry*)NULL); 
     
    return -1;   
  }

}

static void insertHashAll(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber,  hashNumberType position)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      if(treeNumber == 0)
		e->bipNumber = 	e->bipNumber  + 1;
	      else
		e->bipNumber2 = e->bipNumber2 + 1;
	      return;
	    }
	  
	  e = e->next;	 
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
  
      /*e->bitVector  = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int)); */
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));


      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

      if(treeNumber == 0)	
	e->bipNumber  = 1;       	
      else		 
	e->bipNumber2 = 1;
	
      e->next = h->table[position];
      h->table[position] = e;              
    }
  else
    {
      entry *e = initEntry(); 
  
      /*e->bitVector  = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int)); */

      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

      if(treeNumber == 0)	
	e->bipNumber  = 1;	  	
      else    
	e->bipNumber2 = 1;	

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}



static void insertHashBootstop(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber, int treeVectorLength, hashNumberType position)
{    
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 

      e->bipNumber = h->entryCount;
       
      /*e->bitVector  = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int));*/
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));


      e->treeVector = (unsigned int*)rax_calloc(treeVectorLength, sizeof(unsigned int));
      
      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 

      e->bipNumber = h->entryCount;

      /*e->bitVector  = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int));*/

      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      e->treeVector = (unsigned int*)rax_calloc(treeVectorLength, sizeof(unsigned int));

      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}

static void insertHashRF(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, int treeNumber, int treeVectorLength, hashNumberType position, int support, 
			 boolean computeWRF, unsigned int bLink)
{     
  if(h->table[position] != NULL)
    {
      entry *e = h->table[position];     

      do
	{	 
	  unsigned int i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {
	      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
	      if(computeWRF)
		{
		  e->supportVector[treeNumber] = support;
		 
		  assert(0 <= treeNumber && treeNumber < treeVectorLength * MASK_LENGTH);
		}
	      return;
	    }
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
       
      /*e->bitVector  = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int));*/
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));


      e->treeVector = (unsigned int*)rax_calloc(treeVectorLength, sizeof(unsigned int));
      if(computeWRF)
	e->supportVector = (int*)rax_calloc(treeVectorLength * MASK_LENGTH, sizeof(int));

      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      if(computeWRF)
	{
	  e->supportVector[treeNumber] = support;
	 
	  assert(0 <= treeNumber && treeNumber < treeVectorLength * MASK_LENGTH);
	}

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->bLink=bLink;
      
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry *e = initEntry(); 
       
      /*e->bitVector  = (unsigned int*)rax_calloc(vectorLength, sizeof(unsigned int)); */

      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));

      e->treeVector = (unsigned int*)rax_calloc(treeVectorLength, sizeof(unsigned int));
      if(computeWRF)	
	e->supportVector = (int*)rax_calloc(treeVectorLength * MASK_LENGTH, sizeof(int));


      e->treeVector[treeNumber / MASK_LENGTH] |= mask32[treeNumber % MASK_LENGTH];
      if(computeWRF)
	{
	  e->supportVector[treeNumber] = support;
	 
	  assert(0 <= treeNumber && treeNumber < treeVectorLength * MASK_LENGTH);
	}

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      
      e->bLink=bLink;//bLink: Branch link, needed for Partial TC/IC correction to be able to link Bipartitions to bInfos.
      
      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}



void bitVectorInitravSpecial(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf, 
			     int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF)
{
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;          

      do 
	{
	  bitVectorInitravSpecial(bitVectors, q->back, numsp, vectorLength, h, treeNumber, function, bInf, countBranches, treeVectorLength, traverseOnly, computeWRF);
	  q = q->next;
	}
      while(q != p);
           
      //newviewBipartitions(bitVectors, p, numsp, vectorLength);
      newviewBipartitionsMultifurcating(bitVectors, p, numsp, vectorLength);

      assert(p->x);

      if(traverseOnly)
	{
	  if(!(isTip(p->back->number, numsp)))
	    *countBranches =  *countBranches + 1;
	  return;
	}

      if(!(isTip(p->back->number, numsp)))
	{
	  unsigned int *toInsert  = bitVectors[p->number];
	  hashNumberType position = p->hash % h->tableSize;
	 
	  assert(!(toInsert[0] & 1));	 

	  switch(function)
	    {
	    case BIPARTITIONS_ALL:	      
	      insertHashAll(toInsert, h, vectorLength, treeNumber, position);
	      *countBranches =  *countBranches + 1;	
	      break;
	    case GET_BIPARTITIONS_BEST:	   	     
	      insertHash(toInsert, h, vectorLength, *countBranches, position);	     
	      
	      p->bInf            = &bInf[*countBranches];
	      p->back->bInf      = &bInf[*countBranches];        
	      p->bInf->support   = 0;	  	 
	      p->bInf->oP = p;
	      p->bInf->oQ = p->back;
	      
	      *countBranches =  *countBranches + 1;		
	      break;
	    case DRAW_BIPARTITIONS_BEST:	     
	      {
		int found = countHash(toInsert, h, vectorLength, position);
		if(found >= 0)
		  bInf[found].support =  bInf[found].support + 1;
		*countBranches =  *countBranches + 1;
	      }	      
	      break;
	    case BIPARTITIONS_BOOTSTOP:	      
	      insertHashBootstop(toInsert, h, vectorLength, treeNumber, treeVectorLength, position);
	      *countBranches =  *countBranches + 1;
	      break;
	    case BIPARTITIONS_RF:
	      if(computeWRF)
		assert(p->support == p->back->support);
	      insertHashRF(toInsert, h, vectorLength, treeNumber, treeVectorLength, position, p->support, computeWRF, 0);
	      *countBranches =  *countBranches + 1;
	      break;  	
	    case BIPARTITIONS_PARTIAL_TC:
	       if(computeWRF)
		assert(p->support == p->back->support);
	      insertHashRF(toInsert, h, vectorLength, treeNumber, treeVectorLength, position, p->support, computeWRF, *countBranches );

	      p->bInf            = &bInf[*countBranches];
	      p->back->bInf      = &bInf[*countBranches];        
	      p->bInf->support   = 0;	  	 
	      p->bInf->oP = p;
	      p->bInf->oQ = p->back;

	      *countBranches =  *countBranches + 1;
	      break;  
	    default:
	      assert(0);
	    }	  	  
	}
      
    }
}

static void linkBipartitions(nodeptr p, tree *tr, branchInfo *bInf, int *counter, int numberOfTrees)
{
  if(isTip(p->number, tr->mxtips))    
    {
      assert(p->bInf == (branchInfo*) NULL && p->back->bInf == (branchInfo*) NULL);      
      return;
    }
  else
    {
      nodeptr q;          
      
      q = p->next;

      while(q != p)
	{
	  linkBipartitions(q->back, tr, bInf, counter, numberOfTrees);	
	  q = q->next;
	}
     
      if(!(isTip(p->back->number, tr->mxtips)))
	{
	  double support;

	  p->bInf       = &bInf[*counter];
	  p->back->bInf = &bInf[*counter]; 

	  support = ((double)(p->bInf->support)) / ((double) (numberOfTrees));
	  p->bInf->support = (int)(0.5 + support * 100.0);	 	       	  

	  assert(p->bInf->oP == p);
	  assert(p->bInf->oQ == p->back);
	  
	  *counter = *counter + 1;
	}


      return;
    }
}


static int readSingleTree(tree *tr, char *fileName, analdef *adef, boolean readBranches, boolean readNodeLabels, boolean completeTree)
{ 
  FILE 
    *f = myfopen(fileName, "r");

  int 
    numberOfTaxa,
    ch,
    trees = 0;

  while((ch = fgetc(f)) != EOF)
    if(ch == ';')
      trees++;
    
  assert(trees == 1);

  printBothOpen("\n\nFound 1 tree in File %s\n\n", fileName);

  rewind(f);

  treeReadLen(f, tr, readBranches, readNodeLabels, TRUE, adef, completeTree, FALSE);
  
  numberOfTaxa = tr->ntips;

  fclose(f);

  return numberOfTaxa;
}

/*************** function for drawing bipartitions on a bifurcating tree ***********/

void calcBipartitions(tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName)
{  
  branchInfo 
    *bInf;
  unsigned int vLength;

  int 
    numberOfTaxa = 0,
    branchCounter = 0,
    counter = 0,  
    numberOfTrees = 0,
    i;

  unsigned int 
    **bitVectors = initBitVector(tr, &vLength);
  
  hashtable *h = 
    initHashTable(tr->mxtips * 10);

  FILE 
    *treeFile; 
  
  numberOfTaxa = readSingleTree(tr, bestTreeFileName, adef, FALSE, FALSE, TRUE);    
  
  bInf = (branchInfo*)rax_malloc(sizeof(branchInfo) * (tr->mxtips - 3));

  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, GET_BIPARTITIONS_BEST, bInf, &branchCounter, 0, FALSE, FALSE);   

  if(numberOfTaxa != tr->mxtips)
    {
      printBothOpen("The number of taxa in the reference tree file \"%s\" is %d and\n",  bestTreeFileName, numberOfTaxa);
      printBothOpen("is not equal to the number of taxa in the bootstrap tree file \"%s\" which is %d.\n", bootStrapFileName, tr->mxtips);
      printBothOpen("RAxML will exit now with an error ....\n\n");
    }
 
  assert((int)h->entryCount == (tr->mxtips - 3));  
  assert(branchCounter == (tr->mxtips - 3));
  
  treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  numberOfTrees = tr->numberOfTrees;

  checkTreeNumber(numberOfTrees, bootStrapFileName);

  for(i = 0; i < numberOfTrees; i++)
    {                
      int 
	bCount = 0;
      
      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE);
      assert(tr->ntips == tr->mxtips);
      
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, DRAW_BIPARTITIONS_BEST, bInf, &bCount, 0, FALSE, FALSE);      
      assert(bCount == tr->mxtips - 3);      
    }     
  
  fclose(treeFile);
   
  readSingleTree(tr, bestTreeFileName, adef, TRUE, FALSE, TRUE); 
   
  linkBipartitions(tr->nodep[1]->back, tr, bInf, &counter, numberOfTrees);

  assert(counter == branchCounter);

  printBipartitionResult(tr, adef, TRUE, FALSE, bipartitionsFileNameBranchLabels);    

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h); 

  rax_free(bInf);
}


/****** functions for IC computation ***********************************************/

static void insertHash_IC(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position)
{    
  if(h->table[position] != NULL)
    {
      entry 
	*e = h->table[position];     

      do
	{	 
	  unsigned int 
	    i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)
	    {	     
	      e->bipNumber = 	e->bipNumber  + 1;	     
	      return;
	    }
	  
	  e = e->next;	 
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
        
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
   
      e->bipNumber  = 1;      
	
      e->next = h->table[position];
      h->table[position] = e;              
    }
  else
    {
      entry 
	*e = initEntry(); 
  
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));
      memset(e->bitVector, 0, vectorLength * sizeof(unsigned int));
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);

     
      e->bipNumber  = 1;	
      
      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}



static unsigned int findHash_IC(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position)
{
  if(h->table[position] == NULL)         
    return 0;
  {
    entry *e = h->table[position];     

    do
      {	 
	unsigned int 
	  i;

	for(i = 0; i < vectorLength; i++)
	  if(bitVector[i] != e->bitVector[i])
	    goto NEXT;
	   
	return e->bipNumber;	 
      NEXT:
	e = e->next;
      }
    while(e != (entry*)NULL); 
     
    return 0;   
  }
}

static boolean compatibleIC(unsigned int *A, unsigned int *C, unsigned int bvlen)
{
  unsigned int 
    i;
  
  for(i = 0; i < bvlen; i++)        
    if(A[i] & C[i])
      break;
          
  if(i == bvlen)
    return TRUE;
  
  for(i = 0; i < bvlen; i++)
    if(A[i] & ~C[i])
      break;
   
  if(i == bvlen)  
    return TRUE;  
  
  /*
    not required -> ask Andre
  for(i = 0; i < bvlen; i++)
    if(~A[i] & ~C[i])
      break;
   
  if(i == bvlen)  
  return TRUE; */

  
  for(i = 0; i < bvlen; i++)
    if(~A[i] & C[i])
      break;
  
  if(i == bvlen)
    return TRUE;  
  else
    return FALSE;
}

static int sortByBipNumber(const void *a, const void *b)
{       
  int        
    ca = ((*((entry **)a))->bipNumber),
    cb = ((*((entry **)b))->bipNumber);
  
  if (ca == cb) 
    return 0;
  
  return ((ca < cb)?1:-1);
}

static int sortByTreeSet(const void *a, const void *b)
{       
  int        
    ca = ((*((entry **)a))->supportFromTreeset)[0],
    cb = ((*((entry **)b))->supportFromTreeset)[0];
  
  if (ca == cb) 
    return 0;
  
  return ((ca < cb)?1:-1);
}

#define BIP_FILTER

static unsigned int countIncompatibleBipartitions(unsigned int *toInsert, hashtable *h,  hashNumberType vectorLength, unsigned int *maxima, unsigned int *maxCounter, boolean bipNumber, 
						  unsigned int numberOfTrees, unsigned int **maximaBitVectors)
{
  unsigned int   
    threshold = numberOfTrees / 20,
    max = 0,
    entryVectorSize = h->entryCount,
    entryVectorElements = 0,
    k,
    uncompatible = 0;

  entry
    **entryVector = (entry**)rax_malloc(sizeof(entry*) * entryVectorSize);

#ifdef BIP_FILTER
     boolean
       *disregard = (boolean*)rax_malloc(sizeof(boolean) * entryVectorSize);
#endif
    

  for(k = 0; k < entryVectorSize; k++)
    {
      entryVector[k] = (entry*)NULL;
#ifdef BIP_FILTER   
      disregard[k] = FALSE;
#endif
    }

  for(k = 0; k < h->tableSize; k++)	     
    {    
      if(h->table[k] != NULL)
	{
	  entry 
	    *e = h->table[k];	   	  

	  do
	    {	    		  
	      if(!compatibleIC(toInsert, e->bitVector, vectorLength))
		{
		  unsigned int
		    support;

		  if(bipNumber)
		    support = e->bipNumber;
		  else
		    support = e->supportFromTreeset[0];		 		  

		  if(support > max)
		    max = support;		  		
		  
		  uncompatible += support;

		  entryVector[entryVectorElements] = e;
		  entryVectorElements++;
		  assert(entryVectorElements < entryVectorSize);
		}
		  
	      e = e->next;
	    }
	  while(e != NULL);
	}
    }

  if(entryVectorElements == 0)    
    uncompatible = 0;    
  else
    {
      //printf("threshold %d\n", threshold);

      int 
	candidates = 0;

      if(bipNumber)
	{
	  qsort(entryVector, entryVectorElements, sizeof(entry *), sortByBipNumber);
	  assert(max == entryVector[0]->bipNumber);	 
	}
      else
	{
	  qsort(entryVector, entryVectorElements, sizeof(entry *), sortByTreeSet);
	  assert(max == entryVector[0]->supportFromTreeset[0]);	  
	}
      
      for(k = 0; k < entryVectorElements; k++)
	{      
	  unsigned int 
	    j,
	    support;
	  
	  boolean
	    uncompat = TRUE;
	  
	  entry
	    *referenceEntry = entryVector[k];            
	  
	  if(bipNumber)
	    support = entryVector[k]->bipNumber;
	  else
	    support = entryVector[k]->supportFromTreeset[0];
	  
	  if(k > 0)
	    {
	      if(support > threshold)
		{
		  candidates++;
		  
		  //printf("%d\n", support);
#ifdef BIP_FILTER
		  for(j = 0; j < k; j++)
		    {

		      if(!disregard[j])
			{
			  entry 
			    *checkEntry = entryVector[j];
			  
			  if(compatibleIC(checkEntry->bitVector, referenceEntry->bitVector, vectorLength)) 
			    {
			      uncompat = FALSE;
			      break;
			    }
			}
		    }
#endif
		}
	      else
		uncompat = FALSE;
	    }
	  
	  if(uncompat)
	    {	     
#ifdef BIP_FILTER
	      disregard[k] = FALSE;
#endif
	      maximaBitVectors[*maxCounter] = entryVector[k]->bitVector;	  	
	      maxima[*maxCounter] = support;		 	  
	      *maxCounter = *maxCounter + 1;	  
	    }
#ifdef BIP_FILTER
	  else
	    disregard[k] = TRUE;
#endif
	}

      if(0)
	{
	  char 
	    n[1024];
	  
	  strcpy(n, "trace.");
	  strcat(n, run_id);
	  
	  FILE 
	    *f = fopen(n, "a");
	
	  fprintf(f, "conf bips %d candidates %d\n\n", *maxCounter, candidates);
	
	  fclose(f);
	}
    }
    
  rax_free(entryVector);
#ifdef BIP_FILTER
  rax_free(disregard);
#endif
    
  return uncompatible;
}

static double computeIC_Value(unsigned int supportedBips, unsigned int *maxima, unsigned int numberOfTrees, unsigned int maxCounter, boolean computeIC_All, boolean warnNegativeIC)
{
  unsigned int 	
    loopLength,
    i,
    totalBipsAll = supportedBips;

  double 
    ic,
    n = 1 + (double)maxCounter,
    supportFreq;
  
  boolean
    negativeIC = FALSE;    
  
  if(computeIC_All)
    {
      loopLength = maxCounter;
      n = 1 + (double)maxCounter;
    }
  else
    {
      loopLength = 1;
      n = 2.0;
    }

  // should never enter this function when the bip is supported by 100%

  assert(supportedBips < numberOfTrees);

  // must be larger than 0 in this case

  //assert(maxCounter > 0);

  // figure out if the competing bipartition is higher support 
  // can happen for MRE and when drawing IC values on best ML tree 

  if(maxima[0] > supportedBips)
    {
      negativeIC = TRUE;
      
      if(warnNegativeIC)
	 {
	   printBothOpen("\nMax conflicting bipartition frequency: %d is larger than frequency of the included bipartition: %d\n", maxima[0], supportedBips);
	   printBothOpen("This is interesting, but not unexpected when computing extended Majority Rule consensus trees.\n");
	   printBothOpen("Please send an email with the input files and command line\n");
	   printBothOpen("to Alexandros.Stamatakis@gmail.com.\n");
	   printBothOpen("Thank you :-)\n\n");   
	 }      
    }

  for(i = 0; i < loopLength; i++)
    totalBipsAll += maxima[i];  

  //neither support for the actual bipartition, nor for the conflicting ones
  //I am not sure that this will happen, but anyway
  if(totalBipsAll == 0)
    return 0.0;
  
  supportFreq = (double)supportedBips / (double)totalBipsAll;
  
  if(supportedBips == 0)
    ic = log(n);
  else    
    ic = log(n) + supportFreq * log(supportFreq);

  for(i = 0; i < loopLength; i++)
    {
      assert(maxima[i] > 0);

      supportFreq =  (double)maxima[i] / (double)totalBipsAll;
      
      if(maxima[i] != 0)
	ic += supportFreq * log(supportFreq);
    }
  
  ic /= log(n);
  
  if(negativeIC)
    return (-ic);
  else
    return ic;
}

static void printSplit(FILE *f, FILE *v, unsigned int *bitVector, tree *tr, double support, double ic, unsigned int frequency)
{
  int   
    i,
    countLeftTaxa = 0,
    countRightTaxa = 0,
    taxa = 0,
    totalTaxa = 0;

  for(i = 0; i < tr->mxtips; i++)        
    if((bitVector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]))    	
      countLeftTaxa++;

  countRightTaxa = tr->mxtips - countLeftTaxa;

  fprintf(f, "((");

  for(i = 0; i < tr->mxtips; i++)    
    {
      if((bitVector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]))    	
	{
	  fprintf(v, "*");	
	  fprintf(f, "%s", tr->nameList[i+1]);	  	  
	  taxa++;
	  totalTaxa++;
	  if(taxa < countLeftTaxa)
	    fprintf(f, ", ");
	}   
      else
	fprintf(v, "-");

      if((i + 1) % 5 == 0)
	fprintf(v, " ");
    } 

  fprintf(v, "\t%u/%f/%f\n", frequency, support * 100.0, ic);

  fprintf(f, "),(");

  taxa = 0;

  for(i = 0; i < tr->mxtips; i++)    
    {
      if(!(bitVector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]))    	
	{
	  fprintf(f, "%s", tr->nameList[i+1]);	  
	  taxa++;
	  totalTaxa++;
	  if(taxa < countRightTaxa)
	    fprintf(f, ", ");
	}   
    }

  assert(totalTaxa == tr->mxtips);

  fprintf(f, "));\n");
}

static void printFullySupportedSplit(tree *tr, unsigned int *bitVector, unsigned int numberOfTrees)
{
  FILE    
    *v = myfopen(verboseSplitsFileName, "a");

  int 
    i;

  fprintf(v, "partition: \n");

  for(i = 0; i < tr->mxtips; i++)    
    {
      if((bitVector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]))    		
	fprintf(v, "*");	  	
      else
	fprintf(v, "-");

      if((i+1) % 5 == 0)
	fprintf(v, " ");
    } 

  fprintf(v, "\t%u/%f/%f\n\n\n", numberOfTrees, 100.0, 1.0);

  fclose(v);
}

static void printVerboseTaxonNames(tree *tr)
{
  FILE 
    *f = myfopen(verboseSplitsFileName, "w");

  int 
    i;

  fprintf(f, "\n");

  for(i = 1; i <= tr->mxtips; i++)
    fprintf(f, "%d. %s \n", i, tr->nameList[i]);

  fprintf(f, "\n");

  fclose(f);
  
}

static void printVerboseIC(tree *tr, unsigned int supportedBips, unsigned int *toInsert, unsigned int maxCounter, unsigned int *maxima, 
			   unsigned int **maximaBitVectors, unsigned int numberOfTrees, int counter, double ic)
{
  unsigned int 
    i;

  double 
    support = (double)supportedBips / (double)numberOfTrees;

  FILE 
    *f,
    *v = myfopen(verboseSplitsFileName, "a");

  char 
    fileName[1024],
    id[64];

  sprintf(id, "%d", counter);
  strcpy(fileName, workdir);
  strcat(fileName, "RAxML_verboseIC.");
  strcat(fileName, run_id);
  strcat(fileName, ".");
  strcat(fileName, id);

  f = myfopen(fileName, "w");

  printBothOpen("Support for split number %d in tree: %f\n", counter, support);

  fprintf(v, "partition: \n");

  printSplit(f, v, toInsert, tr, support, ic, supportedBips);  

  for(i = 0; i < maxCounter; i++)
    {
      support = (double)maxima[i] / (double)numberOfTrees;

      printBothOpen("Support for conflicting split number %u: %f\n", i, support);

      printSplit(f, v, maximaBitVectors[i], tr, support, ic, maxima[i]);        
    }
 
  printBothOpen("All Newick-formatted splits for this bipartition have been written to file %s\n", fileName);
  printBothOpen("\n\n");

  fprintf(v, "\n\n");

  fclose(f);
  fclose(v);
}




static void bitVectorInitravIC(tree *tr, unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h, int treeNumber, int function, branchInfo *bInf, 
			       int *countBranches, int treeVectorLength, boolean traverseOnly, boolean computeWRF, double *tc, double *tcAll, boolean verboseIC)
{
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;          

      do 
	{
	  bitVectorInitravIC(tr, bitVectors, q->back, numsp, vectorLength, h, treeNumber, function, bInf, countBranches, treeVectorLength, traverseOnly, computeWRF, tc, tcAll, verboseIC);
	  q = q->next;
	}
      while(q != p);
           
      newviewBipartitionsMultifurcating(bitVectors, p, numsp, vectorLength);
      //newviewBipartitions(bitVectors, p, numsp, vectorLength);
      
      assert(p->x);

      if(traverseOnly)
	{
	  if(!(isTip(p->back->number, numsp)))
	    *countBranches =  *countBranches + 1;
	  return;
	}

      if(!(isTip(p->back->number, numsp)))
	{
	  unsigned int *toInsert  = bitVectors[p->number];
	  hashNumberType position = p->hash % h->tableSize;
	 
	  assert(!(toInsert[0] & 1));	 

	  switch(function)
	    {	  
	    case GATHER_BIPARTITIONS_IC:
	      insertHash_IC(toInsert, h, vectorLength, position);
	      *countBranches =  *countBranches + 1;
	      break;
	    case FIND_BIPARTITIONS_IC:
	      {
		unsigned int
		  maxCounter = 0,
		  *maxima = (unsigned int *)rax_calloc(h->entryCount, sizeof(unsigned int)),
		  **maximaBitVectors = (unsigned int **)rax_calloc(h->entryCount, sizeof(unsigned int *)),
		  numberOfTrees = (unsigned int)treeVectorLength,		  
		  supportedBips;
		
		double 		  
		  ic,
		  icAll;		

		//obtain the support for the bipartitions in the tree 
		supportedBips = findHash_IC(toInsert, h, vectorLength, position);

		//printf("supported: %d trees %d\n", supportedBips, numberOfTrees);

		if(supportedBips == numberOfTrees)
		  {
		    ic = 1.0;
		    icAll = 1.0;

		    if(verboseIC)
		      printFullySupportedSplit(tr, toInsert, numberOfTrees);
		  }
		else
		  {
		    unsigned int 
		      incompatibleBipartitions = 0;

		    //find all incompatible bipartitions in the hash table and also 
		    //get the conflicting bipartition with maximum support
		    incompatibleBipartitions = countIncompatibleBipartitions(toInsert, h, vectorLength, maxima, &maxCounter, TRUE, numberOfTrees, maximaBitVectors);

		    if(incompatibleBipartitions == 0) 
		      {
			ic = 1.0;
			icAll = 1.0;

			printBothOpen("WARNING, returning an IC score of 1.0, while only %d out of %d trees support the current bipartition\n", supportedBips, numberOfTrees);
			printBothOpen("The IC is still 1.0, but some input trees do not contain information about this bipartition!\n\n");

			if(verboseIC)
			  printFullySupportedSplit(tr, toInsert, numberOfTrees);
		      }
		    else
		      {
			//this number must be smaller or equal to the total number of trees		   
			
			assert(supportedBips + maxima[0] <= numberOfTrees);	     
			
			//now compute the IC score for this bipartition 
			
			ic    = computeIC_Value(supportedBips, maxima, numberOfTrees, maxCounter, FALSE, FALSE);
			icAll = computeIC_Value(supportedBips, maxima, numberOfTrees, maxCounter, TRUE, FALSE);	
			
			if(verboseIC)		      
			  printVerboseIC(tr, supportedBips, toInsert, maxCounter, maxima, maximaBitVectors, numberOfTrees, *countBranches, ic);
		      }
		  }

		//printf("%d %d %d %d IC %f\n", supportedBips, unsupportedBips, max, treeVectorLength, ic);
		
		p->bInf            = &bInf[*countBranches];
		p->back->bInf      = &bInf[*countBranches];        	     
		p->bInf->oP = p;
		p->bInf->oQ = p->back;

		p->bInf->ic    = ic;
		p->bInf->icAll = icAll;

		//increment tc
		*tc    += ic;
		*tcAll += icAll;
			
		rax_free(maxima);
		rax_free(maximaBitVectors);
	      }
	      *countBranches =  *countBranches + 1;
	      break;
	    default:
	      assert(0);
	    }	  	  
	}
      
    }
}


static boolean allTreesEqualSize(tree *tr, int numberOfTrees, FILE *treeFile)
{
  int 
    i; 

  char    
    buffer[nmlngth + 2]; 

  boolean
    allSameSize = TRUE;

  for(i = 0; i < numberOfTrees; i++)
    {
      int 
	c,
	taxaCount = 0;

      while((c = fgetc(treeFile)) != ';')
	{
	  if(c == '(' || c == ',')
	    {
	      c = fgetc(treeFile);
	      if(c ==  '(' || c == ',')
		ungetc(c, treeFile);
	      else
		{
		  int 
		    lookup = -1,
		    j = 0;	      	     
		  
		  do
		    {
		      buffer[j++] = c;
		      c = fgetc(treeFile);
		    }
		  while(c != ':' && c != ')' && c != ',');
		  buffer[j] = '\0';	    			      		  
		  
		  lookup = lookupWord(buffer, tr->nameHash);

		  if(lookup <= 0)
		    {
		      printf("ERROR: Cannot find tree species: %s\n", buffer);
		      assert(0);
		    }

		  taxaCount++;
		  
		  ungetc(c, treeFile);
		}
	    }   
	}

      printf("Tree %d: found %d taxa, reference tree has %d taxa\n", i, taxaCount, tr->mxtips);

     
      assert(taxaCount <= tr->mxtips);

      allSameSize = allSameSize && (tr->mxtips == taxaCount);
    }

  rewind(treeFile);

  //TODO 
  return allSameSize;
}




static void calcBipartitions_IC(tree *tr, analdef *adef, char *bestTreeFileName, FILE *treeFile, int numberOfTrees)
{  
  branchInfo 
    *bInf;
  
  unsigned int 
    vLength,
    **bitVectors = initBitVector(tr, &vLength);
  
  int 
    bCount = 0,
    //numberOfTaxa = 0,
    //numberOfTrees = 0,
    i;
     
  hashtable *h = 
    initHashTable(tr->mxtips * 10);

  //FILE 
  //  *treeFile; 
  
  double
    rtc = 0.0,
    rtcAll = 0.0,
    tc = 0.0,
    tcAll = 0.0;

  //multifurcations 
  tree
    *inputTree = (tree *)rax_malloc(sizeof(tree));

  //numberOfTaxa = readSingleTree(tr, bestTreeFileName, adef, FALSE, FALSE, TRUE);    
  
  bInf = (branchInfo*)rax_malloc(sizeof(branchInfo) * (tr->mxtips - 3));

  if(adef->verboseIC)
    printVerboseTaxonNames(tr);

  /*if(numberOfTaxa != tr->mxtips)
    {
      printBothOpen("The number of taxa in the reference tree file \"%s\" is %d and\n",  bestTreeFileName, numberOfTaxa);
      printBothOpen("is not equal to the number of taxa in the bootstrap tree file \"%s\" which is %d.\n", bootStrapFileName, tr->mxtips);
      printBothOpen("RAxML will exit now with an error ....\n\n");
    }  
  */
  
  //treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  //numberOfTrees = tr->numberOfTrees;
 
  //multifurcations  
  allocateMultifurcations(tr, inputTree);

  for(i = 0; i < numberOfTrees; i++)
    {     
      int 
	numberOfSplits;                
      
      bCount = 0;      
      
      //printf("tree %d\n", i);

      inputTree->start = (nodeptr)NULL;

      numberOfSplits = readMultifurcatingTree(treeFile, inputTree, adef, FALSE);

      //treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE);
      assert(inputTree->ntips == inputTree->mxtips);
      assert(tr->mxtips == inputTree->mxtips);

      bitVectorInitravIC(inputTree, bitVectors, inputTree->nodep[1]->back, inputTree->mxtips, vLength, h, 0, GATHER_BIPARTITIONS_IC, (branchInfo*)NULL, &bCount, 0, FALSE, FALSE, &tc, &tcAll, FALSE);      
      
      assert(bCount == numberOfSplits);      
    }     
  
  fclose(treeFile);
   
  readSingleTree(tr, bestTreeFileName, adef, TRUE, FALSE, TRUE); 
   
  bCount = 0;
  tc = 0.0;

  bitVectorInitravIC(tr, bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, FIND_BIPARTITIONS_IC, bInf, &bCount, numberOfTrees, FALSE, FALSE, &tc, &tcAll, adef->verboseIC); 

  rtc = tc / (double)(tr->mxtips - 3);

  assert(bCount == tr->mxtips - 3); 
  assert(tc <= (double)(tr->mxtips - 3));
  
  printBothOpen("Tree certainty for this tree: %f\n", tc);
  printBothOpen("Relative tree certainty for this tree: %f\n\n", rtc);

  rtcAll = tcAll / (double)(tr->mxtips - 3);

  printBothOpen("Tree certainty including all conflicting bipartitions (TCA) for this tree: %f\n", tcAll);
  printBothOpen("Relative tree certainty including all conflicting bipartitions (TCA) for this tree: %f\n\n", rtcAll);

  if(adef->verboseIC)
    printBothOpen("Verbose PHYLIP-style formatted bipartition information written to file: %s\n\n",  verboseSplitsFileName);

  printBipartitionResult(tr, adef, TRUE, TRUE, icFileNameBranchLabels);    

  //multifurcations
  freeMultifurcations(inputTree);
  rax_free(inputTree); 

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h); 

  rax_free(bInf);
}

static void calcBipartitions_IC_PartialGeneTrees(tree *tr, analdef *adef, char *bootStrapFileName, FILE *treeFile, int numberOfTrees, int numberOfTaxa);

void calcBipartitions_IC_Global(tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName)
{
  int    
    numberOfTrees,
    numberOfTaxa = readSingleTree(tr, bestTreeFileName, adef, TRUE, FALSE, TRUE);

  FILE 
    *treeFile; 
  

  //must be identical because tr->mxtips was extracted from the reference tree

  assert(tr->mxtips == numberOfTaxa);

  treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  numberOfTrees = tr->numberOfTrees;
  checkTreeNumber(numberOfTrees, bootStrapFileName);

  if(allTreesEqualSize(tr, numberOfTrees, treeFile))
    {
      tr->corrected_IC_Score = FALSE;
      printBothOpen("All trees in tree set have the same size, executing standard IC/TC algorithm ... \n");
      calcBipartitions_IC(tr, adef, bestTreeFileName, treeFile, numberOfTrees);
    }
  else
    {
      tr->corrected_IC_Score = TRUE;
      printBothOpen("The trees in tree set don't have the same size, executing IC/TC algorithm with correction ... \n");
      
      /* 
	 TODO: why has the random number seed check been removed?

	 printBothOpen("This RAxML option has not been activated yet .... exiting now ...\n");
	 
	 exit(-1);
	 
	 if(adef->parsimonySeed <= 0)
	 {
	 printBothOpen("Error: You need to specify a random number seed with \"-p\" for calculating TC/IC scores\n");
	 printBothOpen("on partial gene trees, RAxML will exit now ....\n\n");
	 
	 assert(adef->parsimonySeed > 0);
	 }
      */
      
      calcBipartitions_IC_PartialGeneTrees(tr, adef, bootStrapFileName, treeFile, numberOfTrees, numberOfTaxa);
    }
}


/*******************/

static double testFreq(double *vect1, double *vect2, int n);



void compareBips(tree *tr, char *bootStrapFileName, analdef *adef)
{
  int 
    numberOfTreesAll = 0, 
    numberOfTreesStop = 0,
    i; 
  unsigned int k, entryCount;
  double *vect1, *vect2, p, avg1 = 0.0, avg2 = 0.0, scaleAll, scaleStop;
  int 
    bipAll = 0,
    bipStop = 0;
  char bipFileName[1024];
  FILE 
    *outf,
    *treeFile;
  
  unsigned 
    int vLength;
  
  unsigned int **bitVectors = initBitVector(tr, &vLength);
  hashtable *h = initHashTable(tr->mxtips * 100);    
  uint64_t c1 = 0;
  uint64_t c2 = 0;   


  /*********************************************************************************************************/
  
  treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);  
  numberOfTreesAll = tr->numberOfTrees;
              
  checkTreeNumber(numberOfTreesAll, bootStrapFileName);

  for(i = 0; i < numberOfTreesAll; i++)
    { 
      int 
	bCounter = 0;
      
      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE);                
      assert(tr->mxtips == tr->ntips); 

      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, BIPARTITIONS_ALL, (branchInfo*)NULL, &bCounter, 0, FALSE, FALSE);
      assert(bCounter == tr->mxtips - 3);      
    }
	  
  fclose(treeFile);


  /*********************************************************************************************************/   

  treeFile = getNumberOfTrees(tr, tree_file, adef);
  
  numberOfTreesStop = tr->numberOfTrees;   
  checkTreeNumber(numberOfTreesStop, tree_file);

  for(i = 0; i < numberOfTreesStop; i++)
    {              
      int 
	bCounter = 0;


      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE);      
      assert(tr->mxtips == tr->ntips);
      
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 1, BIPARTITIONS_ALL, (branchInfo*)NULL, &bCounter, 0, FALSE, FALSE);
      assert(bCounter == tr->mxtips - 3);     
    }
	  
  fclose(treeFile);  

  /***************************************************************************************************/
      
  vect1 = (double *)rax_malloc(sizeof(double) * h->entryCount);
  vect2 = (double *)rax_malloc(sizeof(double) * h->entryCount);

  strcpy(bipFileName,         workdir);  
  strcat(bipFileName,         "RAxML_bipartitionFrequencies.");
  strcat(bipFileName,         run_id);

  outf = myfopen(bipFileName, "wb");


  scaleAll  = 1.0 / (double)numberOfTreesAll;
  scaleStop = 1.0 / (double)numberOfTreesStop;

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {
      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];

	  do
	    {
	      c1 += e->bipNumber;
	      c2 += e->bipNumber2;
	      vect1[entryCount] = ((double)e->bipNumber) * scaleAll;
	      if(vect1[entryCount] > 0)
		bipAll++;
	      vect2[entryCount] = ((double)e->bipNumber2) * scaleStop;
	      if(vect2[entryCount] > 0)
		bipStop++;
	      fprintf(outf, "%f %f\n", vect1[entryCount], vect2[entryCount]);
	      entryCount++;
	      e = e->next;
	    }
	  while(e != NULL);
	}     
    }
  
  printBothOpen("%" PRIu64 "%" PRIu64 "\n", c1, c2);

  assert(entryCount == h->entryCount);

  fclose(outf);

  p = testFreq(vect1, vect2, h->entryCount);

  for(k = 0; k < h->entryCount; k++)
    {
      avg1 += vect1[k];
      avg2 += vect2[k];
    }

  avg1 /= ((double)h->entryCount);
  avg2 /= ((double)h->entryCount);
  
  
  printBothOpen("Average [%s]: %1.40f [%s]: %1.40f\n", bootStrapFileName, avg1, tree_file, avg2);
  printBothOpen("Pearson: %f Bipartitions in [%s]: %d Bipartitions in [%s]: %d Total Bipartitions: %d\n", p, bootStrapFileName, bipAll, tree_file, bipStop, h->entryCount);

  printBothOpen("\nFile containing pair-wise bipartition frequencies written to %s\n\n", bipFileName);
  
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h);

  rax_free(vect1);
  rax_free(vect2);

  exit(0);
}

/*************************************************************************************************************/

void computeRF(tree *tr, char *bootStrapFileName, analdef *adef)
{
  int     
    numberOfUniqueTrees = 0,
    treeVectorLength, 
    numberOfTrees = 0, 
    i,
    j, 
    *rfMat,
    *wrfMat,
    *wrf2Mat;

  unsigned int
    vLength; 

  unsigned int 
    k, 
    entryCount,    
    **bitVectors = initBitVector(tr, &vLength);
  
  char rfFileName[1024];

  boolean 
    computeWRF = FALSE,
    *unique;

  double 
    maxRF, 
    avgRF,
    avgWRF,
    avgWRF2;

  FILE 
    *outf,
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  hashtable *h = (hashtable *)NULL;   
  
  numberOfTrees = tr->numberOfTrees;

  checkTreeNumber(numberOfTrees, bootStrapFileName);

  unique = (boolean *)rax_malloc(sizeof(boolean) * numberOfTrees);

  for(i = 0; i < numberOfTrees; i++)
    unique[i] = TRUE;

  h = initHashTable(tr->mxtips * 2 * numberOfTrees); 

  if(numberOfTrees % MASK_LENGTH == 0)
    treeVectorLength = numberOfTrees / MASK_LENGTH;
  else
    treeVectorLength = 1 + (numberOfTrees / MASK_LENGTH);

 
  rfMat = (int*)rax_calloc(numberOfTrees * numberOfTrees, sizeof(int));
  wrfMat = (int*)rax_calloc(numberOfTrees * numberOfTrees, sizeof(int));
  wrf2Mat = (int*)rax_calloc(numberOfTrees * numberOfTrees, sizeof(int));

  for(i = 0; i < numberOfTrees; i++)
    { 
      int 
	bCounter = 0,      
	lcount   = 0;
      
      lcount = treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE); 

      
      
      assert(tr->mxtips == tr->ntips); 
      
      if(i == 0)
	{
	  assert(lcount == 0 || lcount == tr->mxtips - 3);
	  if(lcount == 0)
	    computeWRF = FALSE;
	  else
	    computeWRF = TRUE;
	}
     
      if(computeWRF)
	assert(lcount == tr->mxtips - 3);     	

      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, i, BIPARTITIONS_RF, (branchInfo *)NULL,
			      &bCounter, treeVectorLength, FALSE, computeWRF);
     
      assert(bCounter == tr->mxtips - 3);          
    }
	  
  fclose(treeFile);  

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {    
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];

	  do
	    {
	      unsigned int *vector = e->treeVector;
	      	    	      
	      if(!computeWRF)
		{	
		  i = 0;
		  while(i < numberOfTrees)
		    {
		      if(vector[i / MASK_LENGTH] == 0)
			i += MASK_LENGTH;
		      else
			{		     
			  if((vector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]) > 0)
			    {			
			      int *r = &rfMat[i * numberOfTrees];

			      for(j = 0; j < numberOfTrees; j++)
				if((vector[j / MASK_LENGTH] & mask32[j % MASK_LENGTH]) == 0)
				  r[j]++;			     			    
			    }
			  i++;
			}
		    }	 		  		  		  
		}
	      else
		{
		  int *supportVector = e->supportVector;

		  i = 0;

		  while(i < numberOfTrees)
		    {
		      if(vector[i / MASK_LENGTH] == 0)
			i += MASK_LENGTH;
		      else
			{		     
			  if((vector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH]) > 0)
			    {			
			      int 
				*r = &rfMat[i * numberOfTrees],
				*w   = &wrfMat[i * numberOfTrees],
				*w2  = &wrf2Mat[i * numberOfTrees],
				support = supportVector[i];
			      
			      

			      for(j = 0; j < numberOfTrees; j++)
				{
				  if((vector[j / MASK_LENGTH] & mask32[j % MASK_LENGTH]) == 0)
				    {
				      r[j]++;			     			    
				      w[j] += ABS(support - supportVector[j]);
				      w2[j] += ABS(support - supportVector[j]);
				    }
				  else
				    {
				      w2[j] += ABS(support - supportVector[j]);
				    }

				}
			    }
			  i++;
			}
		    }			  
		}      
	      
	      entryCount++;
	      e = e->next;
	    }
	  while(e != NULL);
	}

     
    }
  assert(entryCount == h->entryCount);  


  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  outf = myfopen(rfFileName, "wb");
  
  maxRF  = ((double)(2 * (tr->mxtips - 3)));
  avgRF  = 0.0;
  avgWRF = 0.0;
  avgWRF2 = 0.0;

  if(!computeWRF)
    {
      for(i = 0; i < numberOfTrees; i++)    
	for(j = i + 1; j < numberOfTrees; j++)
	  rfMat[i * numberOfTrees + j] += rfMat[j * numberOfTrees + i];
    }
  else
    {
      for(i = 0; i < numberOfTrees; i++)    
	for(j = i + 1; j < numberOfTrees; j++)
	  {
	    rfMat[i * numberOfTrees + j] += rfMat[j * numberOfTrees + i];
	    wrfMat[i * numberOfTrees + j] += wrfMat[j * numberOfTrees + i];
	    wrf2Mat[i * numberOfTrees + j] += wrf2Mat[j * numberOfTrees + i];
	  }
    }

  for(i = 0; i < numberOfTrees; i++)
    for(j = i + 1; j < numberOfTrees; j++)
      {
	int    rf = rfMat[i * numberOfTrees + j];
	double rrf = (double)rf / maxRF;

	if(rf == 0 && unique[i])
	  unique[j] = FALSE;
	
	if(computeWRF)
	  {
	    double wrf = wrfMat[i * numberOfTrees + j] / 100.0;
	    double rwrf = wrf / maxRF;
	    double wrf2 = wrf2Mat[i * numberOfTrees + j] / 100.0;
	    double rwrf2 = wrf2 / maxRF;
	
	    fprintf(outf, "%d %d: %d %f, %f %f, %f %f\n", i, j, rf, rrf, wrf, rwrf, wrf2, rwrf2);
	    avgWRF  += rwrf;
	    avgWRF2 += rwrf2; 
	  }
	else
	  fprintf(outf, "%d %d: %d %f\n", i, j, rf, rrf);
	
	avgRF += rrf;
      }
  
  fclose(outf);

  for(i = 0; i < numberOfTrees; i++)
    if(unique[i])
      numberOfUniqueTrees++;
  
  printBothOpen("\n\nAverage relative RF in this set: %f\n", avgRF / ((double)((numberOfTrees * numberOfTrees - numberOfTrees) / 2)));
  printBothOpen("\n\nNumber of unique trees in this tree set: %d\n", numberOfUniqueTrees);
  
  if(computeWRF)
    {
      printBothOpen("\n\nAverage relative WRF in this set: %f\n", avgWRF / ((double)((numberOfTrees * numberOfTrees - numberOfTrees) / 2)));
      printBothOpen("\n\nAverage relative WRF2 in this set: %f\n", avgWRF2 / ((double)((numberOfTrees * numberOfTrees - numberOfTrees) / 2)));
      printBothOpen("\nFile containing all %d pair-wise RF, WRF and WRF2 distances written to file %s\n\n", (numberOfTrees * numberOfTrees - numberOfTrees) / 2, rfFileName);
    }    
  else
    printBothOpen("\nFile containing all %d pair-wise RF distances written to file %s\n\n", (numberOfTrees * numberOfTrees - numberOfTrees) / 2, rfFileName);

  rax_free(unique);
  rax_free(rfMat);
  rax_free(wrfMat);
  rax_free(wrf2Mat);
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h);  

  exit(0);
}

/********************plausibility checker **********************************/



/* we can not use the default hash numbers generated e.g., in the RF code, based on the tree shape.
   we need to compute a hash on the large/long vector that has as many bits as the big tree has taxa */

static hashNumberType oat_hash(unsigned char *p, int len)
{
  unsigned int 
    h = 0;
  int 
    i;
  
  for(i = 0; i < len; i++) 
    {
      h += p[i];
      h += ( h << 10 );
      h ^= ( h >> 6 );
    }
  
  h += ( h << 3 );
  h ^= ( h >> 11 );
  h += ( h << 15 );
  
  return h;
}

/* function that re-hashes bipartitions from the large tree into the new hash table */

static void insertHashPlausibility(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position)
{     
  if(h->table[position] != NULL)
    {
      entry 
	*e = h->table[position];     

      do
	{	 
	  unsigned int 
	    i;
	  
	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)	 	    	     
	    return;	   	    
	  
	  e = e->next;
	}
      while(e != (entry*)NULL); 

      e = initEntry(); 
            
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));                
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);
     
      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      entry 
	*e = initEntry(); 
             
      e->bitVector = (unsigned int*)rax_malloc(vectorLength * sizeof(unsigned int));      
      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * vectorLength);     

      h->table[position] = e;
    }

  h->entryCount =  h->entryCount + 1;
}

/* this function is called while we parse the small trees and extract the bipartitions, it will just look 
   if the bipartition stored in bitVector is already in the hastable, if it is in there this means that 
   this bipartition is also present in the big tree */

static int findHash(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position)
{ 
  if(h->table[position] == NULL)         
    return 0;
  {
    entry *e = h->table[position];     

    do
      {	 
	unsigned int i;

	for(i = 0; i < vectorLength; i++)
	  if(bitVector[i] != e->bitVector[i])
	    goto NEXT;
	   
	return 1;	 
      NEXT:
	e = e->next;
      }
    while(e != (entry*)NULL); 
     
    return 0;   
  }
}




/************************************  new fast code by David Dao *****************************************************************************/

/* Calculates the size of every subtree and extract its bipartition by seperating the subtree from the small tree 
 This Algorithm works on the induced bifurcating subtree only and needs therefore no multifurcation adaption 
 It then counts, how many bipartition is shared with the reference small tree*/
static int rec_findBipartitions(unsigned int ** bitvectors, int* seq, int arraysize, int* translate, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction)
{
  int 
    i, 
    j,
    o, 
    taxon,
    found = 0,
    before;
    
  unsigned int 
    k, 
    *toInsert = (unsigned int*)NULL,
    *bipartition = (unsigned int*)NULL;
    
  hashNumberType 
    position;
  
  int 
    numberOfSplits = 0, /* stop after extracting n-3 splits */ 
    firstTaxon = taxonToReduction[first - 1] + 1,   
    *V = (int *)rax_malloc((arraysize) * sizeof(int)),
    *bipartitions = (int *)rax_malloc((ntips - 3) * sizeof(int));

  for(i = arraysize - 1; i >= 0; i--) 
    {
      V[i] = 1;
    
      /* Extract Bipartiton from inner node! */
      if(!isTip(translate[seq[i]],numsp) && (numberOfSplits < (ntips - 3)))
	{
	  /* we can be sure of bifurcating trees, therefore j < deg(i) = 2 */
	  for(j = 0; j < 2; j++) 	    
	    V[i] = V[i] + V[i + V[i]];	   

	  /* find out the memory efficient index for this taxon which lies between 0.. (ntips - 1) instead of 1..mxtips  
	     We have to subtract 1 because array starts from 0.. */  
	  int 
	    index = taxonToReduction[translate[seq[i]] - 1];

	  /* save bipartition */
	  bipartitions[numberOfSplits] = index;

	  toInsert = bitvectors[index];

	  //printf("We are working with %i \n", index);
	  /* Calculate Bipartition */
	  for(j = 1; j < V[i]; j++) 
	    {   
	      /* look if Preorderlabel at index i + j is a tip! */
	      if(isTip(translate[seq[i + j]],numsp))
		{
		  /* translate taxon number to a number between 0 and smalltreesize */      
		  taxon = taxonToReduction[translate[seq[i + j]] - 1] + 1;

		  /* set bit to one */
		  toInsert[(taxon-1)  / MASK_LENGTH] |= mask32[(taxon-1) % MASK_LENGTH];
		}
	      else
		{      
		  before = taxonToReduction[translate[seq[i+j]] - 1];

		  bipartition = bitvectors[before];

		  for(k = 0; k < vLength; k++) 
		    toInsert[k] |= bipartition[k];

		  /* jump to the next relevant taxon */ 
		  j = j + V[i+j] - 1;
		}        
	    }
	  /* count number of splits and stop at n-3 */
	  numberOfSplits += 1;
	}
    }
  
  for(i=0;i < (ntips - 3); i++) 
    {
      toInsert = bitvectors[bipartitions[i]];

      /* if bitvector contains first taxon, use its complement */
      if(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) 
	{                 
	  /* Padding the last bits! */
	  if(ntips % MASK_LENGTH != 0) 
	    {            
	      for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--)         
		toInsert[vLength - 1] |= mask32[o-1];       
	    }
	  
	  for(k=0;k < vLength; k++)         
	    toInsert[k] = ~toInsert[k];         
	}
      
      assert(!(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]));
      
      /* compute hash */
      position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vLength);
      
      position = position % hash->tableSize;
      
      /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
	 to the taxa contained in the small tree
      */
      
      found = found + findHash(toInsert, hash, vLength, position);              
    }

  rax_free(V);
  rax_free(bipartitions);
  
  return found;
}


/* method adapted for multifurcating trees, changes are: 
we now need information about the degree of an inner node, because it is not 2 anymore 
we also can have less than n-3 splits and therefore need a new parameter telling us the actual split number */
static void rec_extractBipartitionsMulti(unsigned int** bitvectors, int* seq, int arraysize, int numsp, unsigned int vLength, int ntips, int first, hashtable* hash, int* taxonToReduction, int* taxonHasDegree, int maxSplits)
{
  int 
    i,
    j,
    o,
    taxon,
    numberOfSplits = 0,
    firstTaxon = taxonToReduction[first - 1] + 1,
    *V = (int *)rax_malloc((arraysize) * sizeof(int)),
    *bipartitions = (int*)rax_malloc((ntips - 3) * sizeof(int));

  hashNumberType 
    position;
  

  unsigned int 
    k,
    *toInsert = (unsigned int*)NULL,
    *bipartition = (unsigned int*)NULL;


  for(i = arraysize - 1; i >= 0; i--) 
    {
      V[i] = 1;
      /* instead of n-3 we stop after maxSplits */
      if(!isTip(seq[i],numsp) && (numberOfSplits < maxSplits))
	{ 
	  /* instead of j < 2, for multifurcating trees, they can have an abitrarily degree. 
	     We save this information in an array called taxonHasDegree and look it up quickly*/
	  for(j = 0; j < taxonHasDegree[seq[i] - 1] ; j++)      
	    V[i] = V[i] + V[i + V[i]];              
	  
	  int 
	    index = taxonToReduction[seq[i] - 1];

	  bipartitions[numberOfSplits] = index;

	  toInsert = bitvectors[index];

	  /* Extract Bipartition */
	  for(j = 1; j < V[i]; j++) 
	    {
	      if(isTip(seq[i + j],numsp))
		{
		  taxon = taxonToReduction[seq[i + j] - 1] + 1;
		  
		  toInsert[(taxon-1) / MASK_LENGTH] |= mask32[(taxon-1) % MASK_LENGTH];
		}
	      else
		{		  
		  int 
		    before = taxonToReduction[seq[i+j] - 1];
		  
		  bipartition = bitvectors[before];
		  
		  for(k = 0; k < vLength; k++) 
		    toInsert[k] |= bipartition[k];
		  
		  /* jump to the next subtree */
		  j = j + V[i+j] - 1;
		}        
	    }
	  
	  numberOfSplits += 1;          
	}
    }
  /* we now iterate to maxSplits instead of n-3 */
  for(i=0;i < maxSplits; i++) 
    {
      toInsert = bitvectors[bipartitions[i]];
      
      if(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]) 
	{          
	  /* Padding the last bits! */
	  if(ntips % MASK_LENGTH != 0) 
	    {            
	      for(o = MASK_LENGTH; o > (ntips % MASK_LENGTH); o--)         
		toInsert[vLength - 1] |= mask32[o-1];       
	    }

	  for(k=0;k < vLength; k++)         
	    toInsert[k] = ~toInsert[k];         
	}
        
      assert(!(toInsert[(firstTaxon-1) / MASK_LENGTH] & mask32[(firstTaxon-1) % MASK_LENGTH]));

      assert(vLength > 0);

      for(k=0; k < vLength; k++)    
	/* compute hash */
	position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vLength);
    
      position = position % hash->tableSize;
    
      /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
	 to the taxa contained in the small tree
      */
      insertHashPlausibility(toInsert, hash, vLength, position);     
    }

  rax_free(V);
  rax_free(bipartitions);  
}




/*Preordertraversal of the big tree using bitVectorInitrav as reference and taking start->back node, 
number of tips and start->number as parameter and delivers a TaxonToPreOrderLabel and LabelToTaxon Array*/
static void preOrderTraversal(nodeptr p, int numsp, int start, int* array, int* backarray, int* pos)
{
  if(isTip(p->number, numsp))
    {
      array[p->number - 1] = *pos; 
      
      backarray[*pos] = p->number;
      
      *pos = *pos + 1;
      
      return;
    }
  else
    {
      nodeptr q = p->next;
      
      array[p->number - 1] = *pos; 
      
      backarray[*pos] = p->number;
      
      *pos = *pos + 1;
      
      /* get start element */
      if(p->back->number == start) 	
	preOrderTraversal(p->back, numsp, start, array, backarray, pos);       

      do
        {
          preOrderTraversal(q->back, numsp, start, array, backarray, pos);

          q = q->next;
        }
      while(q != p); 
    }
}



/*extract all smalltree taxa and store a list of all Taxon*/
static void rec_extractTaxa(int* smallTreeTaxa, int* taxonToReduction, nodeptr p, int numsp, int* pos, int* pos2)
{ 
  if(isTip(p->number, numsp))
    {
      smallTreeTaxa[*pos] = p->number; 
      taxonToReduction[p->number - 1] = *pos;
      *pos = *pos + 1;
      
      return;
    }
  else
    {    
      nodeptr 
	q = p->next;
      
      taxonToReduction[p->number - 1] = *pos2;      
      *pos2 = *pos2 + 1;
      
      do
	{
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, q->back, numsp, pos, pos2);
	  q = q->next;
	}
      while(q != p);            
  }
}

/* traverses the reference small tree and additionally extracting for every node its degree. It is stored in int* deg */
static void rec_preOrderTraversalMulti(nodeptr p, int numsp, int start, int* backarray, int* deg, int* pos)
{
  int degree = 0;

  if(isTip(p->number, numsp))
    { 
      backarray[*pos] = p->number;
      
      *pos = *pos + 1;
      
      deg[p->number - 1] = degree;
      return;
    }
  else
    {
      nodeptr 
	q = p->next;

      backarray[*pos] = p->number;

      *pos = *pos + 1;

      if(p->back->number == start) 
	{
	  rec_preOrderTraversalMulti(p->back, numsp, start, backarray, deg, pos);
	  degree += 1;
	} 

      do
        {
          rec_preOrderTraversalMulti(q->back, numsp, start, backarray, deg, pos);

          q = q->next;

          degree += 1;
         }
      while(q != p); 

      deg[p->number - 1] = degree;
    }
}



/* special function inits bitvectors to store bipartitions for dynamic programming*/
static unsigned int **rec_initBitVector(tree *tr, unsigned int vectorLength)
{
  unsigned int 
    **bitVectors = (unsigned int **)rax_malloc(sizeof(unsigned int*) * 2 * tr->ntips - 1);
  
  int 
    i;
  
  for(i = 0; i < (2 * tr->ntips - 1); i++)    
    bitVectors[i] = (unsigned int *)rax_calloc(vectorLength, sizeof(unsigned int));
        
  return bitVectors;
}


static void rec_freeBitVector(tree *tr, unsigned int **bitVectors)
{  
  int 
    i;
  
  for(i = 0; i < (2 * tr->ntips - 1); i++)    
    rax_free(bitVectors[i]);
}

/*euler traversal for binary and rooted trees*/                                                                                                                                                              
static void eulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex)
{
  array[*pos] = reference[p->number - 1];

  if (isTip(p->number, numsp)) 
    {
      if (taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;
  
  if(!isTip(p->number, numsp))
    {
      eulerTour(p->next->back, numsp, array, reference, pos, taxonToEulerIndex);
      
      array[*pos] = reference[p->number - 1]; 
      
      *pos = *pos + 1;
      
      eulerTour(p->next->next->back, numsp, array, reference, pos, taxonToEulerIndex);
      
      array[*pos] = reference[p->number - 1];
      
      *pos = *pos + 1;
    }
}

/*For unrooted Trees there is a special case for the arbitrary root which has degree 3 */
static void unrootedEulerTour(nodeptr p, int numsp, int* array, int* reference, int* pos, int* taxonToEulerIndex)
{
  array[*pos] = reference[p->number - 1];
  
  if (isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;
  
  if(!isTip(p->number, numsp))    
    eulerTour(p->back, numsp,array,reference, pos, taxonToEulerIndex);   

  array[*pos] = reference[p->number - 1];

  if(isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;

  if(!isTip(p->number, numsp))    
    eulerTour(p->next->back, numsp,array,reference, pos, taxonToEulerIndex);    

  array[*pos] = reference[p->number - 1];

  if(isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;

  if(!isTip(p->number, numsp))    
    eulerTour(p->next->next->back, numsp,array,reference, pos, taxonToEulerIndex);

  array[*pos] = reference[p->number - 1];

  if(isTip(p->number, numsp)) 
    {
      if(taxonToEulerIndex[p->number - 1] == -1) 
	taxonToEulerIndex[p->number - 1] = *pos;
    }

  *pos = *pos + 1;
}

//function for built-in quicksort

static int sortIntegers(const void *a, const void *b)
{
  int 
    ia = *(int *)(a),
    ib = *(int *)(b);

  if(ia == ib)
    return 0;

  if(ib > ia)
    return -1;
  else
    return 1;
}



/**********************************************************************************/


/************************************* efficient multifurcating algorithm **************************************/

/* Edits compared to the bifurcating algorithm: 
Array storing the degree of every node, called taxonHasDeg,
rec_preOrderTraversalMulti now extracts information about degree of inner nodes,
rec_extractBipartitionsMulti needs two additional parameters: array storing the degree of the nodes and max number of splits */

void plausibilityChecker(tree *tr, analdef *adef)
{
  FILE 
    *treeFile,
    *rfFile;
  
  tree 
    *smallTree = (tree *)rax_malloc(sizeof(tree));

  char 
    rfFileName[1024];

  int
    numberOfTreesAnalyzed = 0,
    i;

  double 
    avgRF = 0.0,
    sumEffectivetime = 0.0;

  /* set up an output file name */

  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  rfFile = myfopen(rfFileName, "wb");  

  assert(adef->mode ==  PLAUSIBILITY_CHECKER);

  /* open the big reference tree file and parse it */

  treeFile = myfopen(tree_file, "r");

  printBothOpen("Parsing reference tree %s\n", tree_file);

  treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);

  assert(tr->mxtips == tr->ntips);
  
  /*************************************************************************************/
  /* Preprocessing Step */

  double 
    preprocesstime = gettime();
  
  /* taxonToLabel[2*tr->mxtips - 2]; 
  Array storing all 2n-2 labels from the preordertraversal: (Taxonnumber - 1) -> (Preorderlabel) */
  int 
    *taxonToLabel  = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int)),

    /* taxonHasDeg[2*tr->mxtips - 2] 
    Array used to store the degree of every taxon, is needed to extract Bipartitions from multifurcating trees 
    (Taxonnumber - 1) -> (degree of node(Taxonnumber)) */

    *taxonHasDeg = (int *)rax_calloc((2*tr->mxtips - 2),sizeof(int)),

    /* taxonToReduction[2*tr->mxtips - 2]; 
  Array used for reducing bitvector and speeding up extraction: 
  (Taxonnumber - 1) -> (0..1 (increment count of taxa appearing in small tree))
  (Taxonnumber - 1) -> (0..1 (increment count of inner nodes appearing in small tree)) */

    *taxonToReduction = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
    
  int 
    newcount = 0; //counter used for correct traversals

  /* labelToTaxon[2*tr->mxtips - 2];
  is used to translate between Perorderlabel and p->number: (Preorderlabel) -> (Taxonnumber) */
  int 
    *labelToTaxon = (int *)rax_malloc((2*tr->mxtips - 2) * sizeof(int));
  
  /* Preorder-Traversal of the large tree */
  preOrderTraversal(tr->start->back,tr->mxtips, tr->start->number, taxonToLabel, labelToTaxon, &newcount);

  newcount = 0; //counter set to 0 to be now used for Eulertraversal

  /* eulerIndexToLabel[4*tr->mxtips - 5]; 
  Array storing all 4n-5 PreOrderlabels created during eulertour: (Eulerindex) -> (Preorderlabel) */
  int* 
    eulerIndexToLabel = (int *)rax_malloc((4*tr->mxtips - 5) * sizeof(int));

  /* taxonToEulerIndex[tr->mxtips]; 
  Stores all indices of the first appearance of a taxa in the eulerTour: (Taxonnumber - 1) -> (Index of the Eulertour where Taxonnumber first appears) 
  is used for efficient computation of the Lowest Common Ancestor during Reconstruction Step
  */
  int*
    taxonToEulerIndex  = (int *)rax_malloc((tr->mxtips) * sizeof(int));

  /* Init taxonToEulerIndex and taxonToReduction */
  int 
    ix;

  for(ix = 0; ix < tr->mxtips; ++ix)    
    taxonToEulerIndex[ix] = -1;    
  
  for(ix = 0; ix < (2*tr->mxtips - 2); ++ix)    
    taxonToReduction[ix] = -1;    


  /* Eulertraversal of the large tree*/
  unrootedEulerTour(tr->start->back,tr->mxtips, eulerIndexToLabel, taxonToLabel, &newcount, taxonToEulerIndex);

  /* Creating RMQ Datastructure for efficient retrieval of LCAs, using Johannes Fischers Library rewritten in C
  Following Files: rmq.h,rmqs.c,rmqs.h are included in Makefile.RMQ.gcc */
  RMQ_succinct(eulerIndexToLabel,4*tr->mxtips - 5);

  double 
    preprocessendtime = gettime() - preprocesstime;

  /* Proprocessing Step End */
  /*************************************************************************************/

  printBothOpen("The reference tree has %d tips\n", tr->ntips);

  fclose(treeFile);
  
  /* now see how many small trees we have */

  treeFile = getNumberOfTrees(tr, bootStrapFile, adef);

  checkTreeNumber(tr->numberOfTrees, bootStrapFile);

  /* allocate a data structure for parsing the potentially mult-furcating tree */

  allocateMultifurcations(tr, smallTree);

  /* loop over all small trees */

  for(i = 0; i < tr->numberOfTrees;  i++)
    {      
      int
	numberOfSplits = readMultifurcatingTree(treeFile, smallTree, adef, TRUE);
      
      if(numberOfSplits > 0)
	{
	  int
	    firstTaxon;           

	  double
	    rec_rf,
	    maxRF;

	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Small tree %d has %d tips and %d bipartitions\n", i, smallTree->ntips, numberOfSplits);    
	  
	  /* compute the maximum RF distance for computing the relative RF distance later-on */
	  
	  /* note that here we need to pay attention, since the RF distance is not normalized 
	     by 2 * (n-3) but we need to account for the fact that the multifurcating small tree 
	     will potentially contain less bipartitions. 
	     Hence the normalization factor is obtained as n-3 + numberOfSplits, where n-3 is the number 
	     of bipartitions of the pruned down large reference tree for which we know that it is 
	     bifurcating/strictly binary */
	  
	  maxRF = (double)(2 * numberOfSplits);
	  
	  /* now get the index of the first taxon of the small tree.
	     we will use this to unambiguously store the bipartitions 
	  */
	  
	  firstTaxon = smallTree->start->number;
	  
	  /***********************************************************************************/
	  /* Reconstruction Step */
	  
	  double 
	    time_start = gettime();
	  
	  /* Init hashtable to store Bipartitions of the induced subtree */
	  /* 
	     using smallTree->ntips instead of smallTree->mxtips yields faster code 
	     e.g. 120 versus 128 seconds for 20,000 small trees on my laptop 
	   */
	  hashtable
	    *s_hash = initHashTable(smallTree->ntips * 4);
	  
	  /* smallTreeTaxa[smallTree->ntips]; 
	     Stores all taxa numbers from smallTree into an array called smallTreeTaxa: (Index) -> (Taxonnumber)  */
	  int* 
	    smallTreeTaxa = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
	  
	  /* counter is set to 0 for correctly extracting taxa of the small tree */
	  newcount = 0; 
	  
	  int 
	    newcount2 = 0;
	  
	  /* seq2[2*smallTree->ntips - 2]; 
	     stores PreorderSequence of the reference smalltree: (Preorderindex) -> (Taxonnumber) */
	  int* 
	    seq2 = (int *)rax_malloc((2*smallTree->ntips - 2) * sizeof(int));
	  /* used to store the vectorLength of the bitvector */
	  unsigned int 
	    vectorLength;
	  
	  /* extract all taxa of the smalltree and store it into an array, 
	     also store all counts of taxa and nontaxa in taxonToReduction */
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, smallTree->start, smallTree->mxtips, &newcount, &newcount2);
	  
	  rec_extractTaxa(smallTreeTaxa, taxonToReduction, smallTree->start->back, smallTree->mxtips, &newcount, &newcount2);
	  
	  /* counter is set to 0 to correctly preorder traverse the small tree */
	  newcount = 0;
	  
	  /* Preordertraversal of the small tree and save its sequence into seq2 for later extracting the bipartitions, it
	     also stores information about the degree of every node */
	  
	  rec_preOrderTraversalMulti(smallTree->start->back,smallTree->mxtips, smallTree->start->number, seq2, taxonHasDeg, &newcount);
	  
	  /* calculate the bitvector length */
	  if(smallTree->ntips % MASK_LENGTH == 0)
	    vectorLength = smallTree->ntips / MASK_LENGTH;
	  else
	    vectorLength = 1 + (smallTree->ntips / MASK_LENGTH); 
	  
	  unsigned int 
	    **bitVectors = rec_initBitVector(smallTree, vectorLength);
	  
	  /* store all non trivial bitvectors using an subtree approach for the induced subtree and 
	     store it into a hashtable, this method was changed for multifurcation */
	  rec_extractBipartitionsMulti(bitVectors, seq2, newcount,tr->mxtips, vectorLength, smallTree->ntips, 
				       firstTaxon, s_hash, taxonToReduction, taxonHasDeg, numberOfSplits);
	  
	  /* counter is set to 0 to be used for correctly storing all EulerIndices */
	  newcount = 0; 
	  
	  /* smallTreeTaxonToEulerIndex[smallTree->ntips]; 
	     Saves all first Euler indices for all Taxons appearing in small Tree: 
	     (Index) -> (Index of the Eulertour where the taxonnumber of the small tree first appears) */
	  int* 
	    smallTreeTaxonToEulerIndex = (int *)rax_malloc((smallTree->ntips) * sizeof(int));
	  
	  /* seq[(smallTree->ntips*2) - 1] 
	     Stores the Preordersequence of the induced small tree */
	  int* 
	    seq = (int *)rax_malloc((2*smallTree->ntips - 1) * sizeof(int));
	  
	  
	  /* iterate through all small tree taxa */
	  for(ix = 0; ix < smallTree->ntips; ix++) 
	    {        
	      int 
		taxanumber = smallTreeTaxa[ix];
	      
	      /* To create smallTreeTaxonToEulerIndex we filter taxonToEulerIndex for taxa in the small tree*/
	      smallTreeTaxonToEulerIndex[newcount] = taxonToEulerIndex[taxanumber-1]; 
	      
	      /* Saves all Preorderlabel of the smalltree taxa in seq*/
	      seq[newcount] = taxonToLabel[taxanumber-1];
	      
	      newcount++;
	    }
	  
	  /* sort the euler indices to correctly calculate LCA */
	  //quicksort(smallTreeTaxonToEulerIndex,0,newcount - 1);             
	  
	  qsort(smallTreeTaxonToEulerIndex, newcount, sizeof(int), sortIntegers);
	  
	  //printf("newcount2 %i \n", newcount2);      
	  /* Iterate through all small tree taxa */
	  for(ix = 1; ix < newcount; ix++)
	    {  
	      /* query LCAs using RMQ Datastructure */
	      seq[newcount - 1 + ix] =  eulerIndexToLabel[query(smallTreeTaxonToEulerIndex[ix - 1],smallTreeTaxonToEulerIndex[ix])]; 	 
	      
	      /* Used for dynamic programming. We save an index for every inner node:
		 For example the reference tree has 3 inner nodes which we saves them as 0,1,2.
		 Now we calculate for example 5 LCA to construct the induced subtree, which are also inner nodes. 
		 Therefore we mark them as 3,4,5,6,7  */
	      
	      taxonToReduction[labelToTaxon[seq[newcount - 1 + ix]] - 1] = newcount2;
	      
	      newcount2 += 1;
	    }
	  
	  /* sort to construct the Preordersequence of the induced subtree */
	  //quicksort(seq,0,(2*smallTree->ntips - 2));
	  
	  qsort(seq, (2 * smallTree->ntips - 2) + 1, sizeof(int), sortIntegers);
	  
	  /* calculates all bipartitions of the reference small tree and count how many bipartition it shares with the induced small tree */
	  int 
	    rec_bips = rec_findBipartitions(bitVectors, seq,(2*smallTree->ntips - 1), labelToTaxon, tr->mxtips, vectorLength, smallTree->ntips, firstTaxon, s_hash, taxonToReduction);
	  
	  /* Reconstruction Step End */
	  /***********************************************************************************/
	  
	  double 
	    effectivetime = gettime() - time_start;
	  
	  /*
	    if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Reconstruction time: %.10f secs\n\n", effectivetime);
	  */
	  
	  /* compute the relative RF */
	  
	  rec_rf = (double)(2 * (numberOfSplits - rec_bips)) / maxRF;
	  
	  assert(numberOfSplits >= rec_bips);	  	 

	  avgRF += rec_rf;
	  sumEffectivetime += effectivetime;
	  
	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Relative RF tree %d: %f\n\n", i, rec_rf);
	  
	  fprintf(rfFile, "%d %f\n", i, rec_rf);
	  
	  /* free masks and hast table for this iteration */
	  rec_freeBitVector(smallTree, bitVectors);
	  rax_free(bitVectors);
	  
	  freeHashTable(s_hash);
	  rax_free(s_hash);
	  
	  rax_free(smallTreeTaxa);
	  rax_free(seq);
	  rax_free(seq2);
	  rax_free(smallTreeTaxonToEulerIndex);

	  numberOfTreesAnalyzed++;
	}
    }
  
  printBothOpen("Number of small trees skipped: %d\n\n", tr->numberOfTrees - numberOfTreesAnalyzed);
  
  printBothOpen("Average RF distance %f\n\n", avgRF / (double)numberOfTreesAnalyzed);
  
  printBothOpen("Large Tree: %i, Number of SmallTrees analyzed: %i \n\n", tr->mxtips, numberOfTreesAnalyzed); 
  
  printBothOpen("Total execution time: %f secs\n\n", gettime() - masterTime);
   
  printBothOpen("File containing all %d pair-wise RF distances written to file %s\n\n", numberOfTreesAnalyzed, rfFileName);

  printBothOpen("execution stats:\n\n");
  printBothOpen("Accumulated time Effective algorithm: %.5f sec \n", sumEffectivetime);
  printBothOpen("Average time for effective: %.10f sec \n",sumEffectivetime / (double)numberOfTreesAnalyzed);
  printBothOpen("Preprocessingtime: %0.5f sec \n\n", preprocessendtime);
 

  fclose(treeFile);
  fclose(rfFile);    
  
  /* free the data structure used for parsing the potentially multi-furcating tree */

  freeMultifurcations(smallTree);
  rax_free(smallTree);

  rax_free(taxonToLabel);
  rax_free(taxonToEulerIndex);
  rax_free(labelToTaxon);
  rax_free(eulerIndexToLabel);
  rax_free(taxonToReduction);
  rax_free(taxonHasDeg);
}



/********************************************************/

double convergenceCriterion(hashtable *h, int mxtips)
{
  int      
    rf = 0; 

  unsigned int 
    k = 0, 
    entryCount = 0;
  
  double    
    rrf;  

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {      
      if(h->table[k] != NULL)
	{
	  entry *e = h->table[k];

	  do
	    {
	      unsigned int *vector = e->treeVector;	     
	      if(((vector[0] & 1) > 0) + ((vector[0] & 2) > 0) == 1)
		rf++;	     
	      
	      entryCount++;
	      e = e->next;
	    }
	  while(e != NULL);
	}     
    }

  assert(entryCount == h->entryCount);  
      
  rrf = (double)rf/((double)(2 * (mxtips - 3)));  
 
  return rrf;
}




/*************************************************************************************************************/

static void permute(unsigned int *perm, unsigned int n, int64_t *seed)
{
  unsigned int  i, j, k;
 
  for (i = 0; i < n; i++) 
    {
      k =  (int)((double)(n - i) * randum(seed));
      j        = perm[i];    
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
      /*assert(i + k < n);*/
    }
}





static double testFreq(double *vect1, double *vect2, int n)
{
  int 
    i;
  
  boolean 
    allEqual = TRUE;

  double
    avg1 = 0.0, 
    avg2 = 0.0,
    sum_xy = 0.0, 
    sum_x  = 0.0, 
    sum_y  = 0.0,
    corr   = 0.0;
 
  for(i = 0; i < n; i++)
    {	     
      allEqual = allEqual && (vect1[i] == vect2[i]);

      avg1 += vect1[i];
      avg2 += vect2[i];
    }
      
  avg1 /= ((double)n);
  avg2 /= ((double)n); 
      
  for(i = 0; i < n; i++)
    {
      sum_xy += ((vect1[i] - avg1) * (vect2[i] - avg2));
      sum_x  += ((vect1[i] - avg1) * (vect1[i] - avg1));
      sum_y  += ((vect2[i] - avg2) * (vect2[i] - avg2));	 
    }       

  if(allEqual)
    return 1.0;

  if(sum_x == 0.0 || sum_y == 0.0) 
    return 0.0;

  corr = sum_xy / (sqrt(sum_x) * sqrt(sum_y));
   
  /*
    #ifndef WIN32
    if(isnan(corr))
    {
    printf("Numerical Error pearson correlation is not a number\n");
    assert(0);
    }
    #endif
  */

  return corr;
}

static double frequencyCriterion(int numberOfTrees, hashtable *h, int *countBetter, int bootstopPermutations, analdef *adef)
{
  int 
    k, 
    l;
    
  int64_t 
    seed = adef->parsimonySeed;

  double     
    result, 
    avg = 0, 
    *vect1, 
    *vect2; 

  unsigned int
    *perm =  (unsigned int *)rax_malloc(sizeof(unsigned int) * numberOfTrees),
    j;

  assert(*countBetter == 0);
  assert(seed > 0);
	  
#ifdef _WAYNE_MPI
  seed = seed + 10000 * processID;
#endif
	  
  for(j = 0; j < (unsigned int)numberOfTrees; j++)
    perm[j] = j;
	  
  for(k = 0; k < bootstopPermutations; k++)
    {   	      		      	      
      unsigned int entryCount = 0;

      permute(perm, numberOfTrees, &seed);
      
     

      vect1 = (double *)rax_calloc(h->entryCount, sizeof(double));
      vect2 = (double *)rax_calloc(h->entryCount, sizeof(double));	     

        
      
      for(j = 0; j < h->tableSize; j++)
	{		
	  if(h->table[j] != NULL)
	    {		
	      entry *e = h->table[j];
	      
	      do
		{
		  unsigned int *set = e->treeVector;       
		  
		  for(l = 0; l < numberOfTrees; l++)
		    {			     
		      if((set[l / MASK_LENGTH] != 0) && (set[l / MASK_LENGTH] & mask32[l % MASK_LENGTH]))
			{
			  if(perm[l] % 2 == 0)
			    vect1[entryCount] = vect1[entryCount] + 1.0;
			  else			
			    vect2[entryCount] = vect2[entryCount] + 1.0;
			}
		    }
		  entryCount++;
		  e = e->next;
		}
	      while(e != NULL);
	    }			     
	}		    	
      
      
      
      
      assert(entryCount == h->entryCount);
      
      

      result = testFreq(vect1, vect2, entryCount);
	  
          
      
      if(result >= FC_LOWER)
	*countBetter = *countBetter + 1;
	      
      avg += result;
	      
      rax_free(vect1);		  
      rax_free(vect2);		 
    }

  rax_free(perm);
	  
  avg /= bootstopPermutations;
	
      

  return avg;
}




static double wcCriterion(int numberOfTrees, hashtable *h, int *countBetter, double *wrf_thresh_avg, double *wrf_avg, tree *tr, unsigned int vectorLength, int bootstopPermutations, analdef *adef)
{
  int 
    k, 
    l,   
    wrf,
    mr_thresh = ((double)numberOfTrees/4.0);
   
  unsigned int 
    *perm =  (unsigned int *)rax_malloc(sizeof(unsigned int) * numberOfTrees),
    j;

  int64_t 
    seed = adef->parsimonySeed;  
  
  double 
    wrf_thresh = 0.0,
    pct_avg = 0.0;

#ifdef _WAYNE_MPI
  seed = seed + 10000 * processID;
#endif

  assert(seed > 0);

  assert(*countBetter == 0 && *wrf_thresh_avg == 0.0 && *wrf_avg == 0.0);
	   	  
  for(j = 0; j < (unsigned int)numberOfTrees; j++)
    perm[j] = j;
	  
  for(k = 0; k < bootstopPermutations; k++)
    {   	      		           
      int mcnt1 = 0;			  
      int mcnt2 = 0;
      unsigned int entryCount = 0;
      double halfOfConsideredBips = 0.0;

      entry ** sortedByKeyA = (entry **)NULL;
      entry ** sortedByKeyB = (entry **)NULL;
      int lenA, lenB;
      boolean ignoreCompatibilityP;

      int iA, iB;
      wrf = 0;
      	      
      permute(perm, numberOfTrees, &seed);      
    
      for(j = 0; j < h->tableSize; j++)
	{		
	  if(h->table[j] != NULL)
	    {
	      entry *e = h->table[j];

	        do
		  {
		    int cnt1 = 0;
		    int cnt2 = 0;

		    unsigned int *set = e->treeVector;

		    for(l = 0; l < numberOfTrees; l++)
		      {			     
			if((set[l / MASK_LENGTH] != 0) && (set[l / MASK_LENGTH] & mask32[l % MASK_LENGTH]))
			  {			    
			    if(perm[l] % 2 == 0)
			      cnt1++;
			    else			
			      cnt2++;
			  }			     
		      }
		    
		    switch(tr->bootStopCriterion)
		      {
		      case MR_STOP:
			if(cnt1 <= mr_thresh)			      
			  cnt1 = 0;
		       
			if(cnt2 <= mr_thresh)	    
			  cnt2 = 0;

			if(cnt1 > 0)			      
			  mcnt1++;
			
			if(cnt2 > 0)			      
			  mcnt2++;

			wrf += ((cnt1 > cnt2) ? cnt1 - cnt2 : cnt2 - cnt1);
			break;
		      case MRE_STOP:
		      case MRE_IGN_STOP:
			e->supportFromTreeset[0] = cnt1;
			e->supportFromTreeset[1] = cnt2;
			break;
		      default:
			assert(0);
		      }

		    entryCount++;
		    e = e->next;
		  }
		while(e != NULL);
	    }	  	  	  	
	}	
      
      assert(entryCount == h->entryCount);

      if((tr->bootStopCriterion == MRE_STOP) || (tr->bootStopCriterion == MRE_IGN_STOP))
	{
	
	  
	  if (tr->bootStopCriterion == MRE_IGN_STOP)
	    ignoreCompatibilityP = TRUE;
	  else
	    ignoreCompatibilityP = FALSE;
	    
	 
	 
	  mre(h, ignoreCompatibilityP, &sortedByKeyA, &lenA, 0, tr->mxtips, vectorLength, TRUE, tr, TRUE);
	  mre(h, ignoreCompatibilityP, &sortedByKeyB, &lenB, 1, tr->mxtips, vectorLength, TRUE, tr, TRUE);
	  
	   
	  mcnt1 = lenA;
	  mcnt2 = lenB;
	  
	  iA = iB = 0;
	  
	  while(iA < mcnt1 || iB < mcnt2)
	    {	    
	      if( iB == mcnt2 || (iA < mcnt1 && sortedByKeyA[iA] < sortedByKeyB[iB]) )
		{
		  wrf += sortedByKeyA[iA]->supportFromTreeset[0];
		  iA++;
		}
	      else
		{
		  if( iA == mcnt1 || (iB < mcnt2 && sortedByKeyB[iB] < sortedByKeyA[iA]) )
		    {
		      wrf += sortedByKeyB[iB]->supportFromTreeset[1];
		      iB++;
		    }
		  else
		    {
		      int cnt1, cnt2;
		      
		      assert (sortedByKeyA[iA] == sortedByKeyB[iB]);
		      
		      cnt1 = sortedByKeyA[iA]->supportFromTreeset[0];
		      cnt2 = sortedByKeyB[iB]->supportFromTreeset[1];
		      
		      wrf += ((cnt1 > cnt2) ? cnt1 - cnt2 : cnt2 - cnt1);
		      
		      iA++; 
		      iB++;
		    }
		}
	    }

	  rax_free(sortedByKeyA);
	  rax_free(sortedByKeyB);
	  
	  assert (iA == mcnt1);
	  assert (iB == mcnt2);	  
	}

      halfOfConsideredBips = ( ((((double)numberOfTrees/2.0) * (double)mcnt1)) + ((((double)numberOfTrees/2.0) * (double)mcnt2)) );

      /* 
	 wrf_thresh is the 'custom' threshold computed for this pair
	 of majority rules trees (i.e. one of the BS_PERMS splits),
	 and simply takes into account the resolution of the two trees
      */

      wrf_thresh = (tr->wcThreshold) * halfOfConsideredBips;    
      
      /*
	we count this random split as 'succeeding' when
	 the wrf between maj rules trees is exceeded
	 by its custom threshold
      */

      if((double)wrf <= wrf_thresh)			        
	*countBetter = *countBetter + 1;

      /* 
	 here we accumulate outcomes and thresholds, because
	 we're not going to stop until the avg dist is less
	 than the avg threshold
      */

      pct_avg += (double)wrf / halfOfConsideredBips  * 100.0;
      *wrf_avg += (double)wrf;
      *wrf_thresh_avg += wrf_thresh;
    }
 
  rax_free(perm);

  pct_avg /= (double)bootstopPermutations; 
  *wrf_avg /= (double)bootstopPermutations; 
  *wrf_thresh_avg /= (double)bootstopPermutations;   

  /*printf("%d \t\t %f \t\t %d \t\t\t\t %f\n", numberOfTrees, *wrf_avg, *countBetter, *wrf_thresh_avg);	  	      */

  return pct_avg; 
}	  






void computeBootStopOnly(tree *tr, char *bootStrapFileName, analdef *adef)
{
  int numberOfTrees = 0, i;
  boolean stop = FALSE;
  double avg;
  int checkEvery;
  int treesAdded = 0;
  hashtable *h = initHashTable(tr->mxtips * FC_INIT * 10); 
  unsigned int 
    treeVectorLength, 
    vectorLength;
  unsigned int **bitVectors = initBitVector(tr, &vectorLength);   
 

  FILE 
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  checkTreeNumber(tr->numberOfTrees, bootStrapFileName);

  assert((FC_SPACING % 2 == 0) && (adef->fcThreshold < adef->bootstopPermutations));
   
  numberOfTrees = tr->numberOfTrees;
 
  
  printBothOpen("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  
  assert(sizeof(unsigned char) == 1);
  
  if(numberOfTrees % MASK_LENGTH == 0)
    treeVectorLength = numberOfTrees / MASK_LENGTH;
  else
    treeVectorLength = 1 + (numberOfTrees / MASK_LENGTH);  
 
  checkEvery = FC_SPACING;
        
  switch(tr->bootStopCriterion)
    {
    case FREQUENCY_STOP:
      printBothOpen("# Trees \t Average Pearson Coefficient \t # Permutations: pearson >= %f\n", 
		    FC_LOWER);
      break;
    case MR_STOP:
    case MRE_STOP:
    case MRE_IGN_STOP:
      printBothOpen("# Trees \t Avg WRF in %s \t # Perms: wrf <= %1.2f %s\n","%", 100.0 * tr->wcThreshold, "%");
      break;
    default:    
      assert(0);
    }
  
  for(i = 1; i <= numberOfTrees && !stop; i++)
    {                  
      int 
	bCount = 0;           
          
      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE); 
      assert(tr->mxtips == tr->ntips);
      
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			      &bCount, treeVectorLength, FALSE, FALSE);
      
      assert(bCount == tr->mxtips - 3);
                 
      treesAdded++;	
            
      if((i > START_BSTOP_TEST) && (i % checkEvery == 0))
	{ 
	  int countBetter = 0;
	  	 
	  switch(tr->bootStopCriterion)
	    {
	    case FREQUENCY_STOP:
	      avg = frequencyCriterion(i, h, &countBetter, adef->bootstopPermutations, adef);	  	  	  
	      printBothOpen("%d \t\t\t %f \t\t\t\t %d\n", i, avg, countBetter);
	  	  
	      stop = (countBetter >= adef->fcThreshold && avg >= FC_LOWER);	  	 
	      break;
	    case MR_STOP:
	    case MRE_STOP:
	    case MRE_IGN_STOP:
	      {
		double 
		  wrf_thresh_avg = 0.0,
		  wrf_avg = 0.0;
		avg = wcCriterion(i, h, &countBetter, &wrf_thresh_avg, &wrf_avg, tr, vectorLength, adef->bootstopPermutations, adef);
		printBothOpen("%d \t\t %1.2f \t\t\t %d\n", i, avg, countBetter);	       
		
		stop = (countBetter >= adef->fcThreshold && wrf_avg <= wrf_thresh_avg);
	      }
	      break;
	    default:
	      assert(0);
	    }	 
	}	 	   
      
    }
  
 

  if(stop)              
    printBothOpen("Converged after %d replicates\n", treesAdded);           
  else    
    printBothOpen("Bootstopping test did not converge after %d trees\n", treesAdded);     

  fclose(treeFile); 
  
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h);

  

 
  exit(0);
}

#ifdef _WAYNE_MPI

boolean computeBootStopMPI(tree *tr, char *bootStrapFileName, analdef *adef, double *pearsonAverage)
{
  boolean
    stop = FALSE;

  int 
    bootStopPermutations = 0,
    numberOfTrees = 0, 
    i,
    countBetter = 0;
  
  unsigned int
    treeVectorLength, 
    vectorLength;

  double 
    avg;
  
  hashtable 
    *h = initHashTable(tr->mxtips * FC_INIT * 10); 
 
  unsigned 
    int **bitVectors = initBitVector(tr, &vectorLength);   

  
  FILE 
    *treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);
  
  numberOfTrees = tr->numberOfTrees;

  checkTreeNumber(numberOfTrees, bootStrapFileName);

  if(numberOfTrees % 2 != 0)
    numberOfTrees--;
   
  /*printf("\n\nProcess %d Found %d trees in File %s\n\n", processID, numberOfTrees, bootStrapFileName);*/
  
  assert(sizeof(unsigned char) == 1);
  

  if(adef->bootstopPermutations % processes == 0)
    bootStopPermutations = adef->bootstopPermutations / processes;
  else
    bootStopPermutations = 1 + (adef->bootstopPermutations / processes);
  
  /*printf("Perms %d\n",  bootStopPermutations);*/

  if(numberOfTrees % MASK_LENGTH == 0)
    treeVectorLength = numberOfTrees / MASK_LENGTH;
  else
    treeVectorLength = 1 + (numberOfTrees / MASK_LENGTH);    
  
  for(i = 1; i <= numberOfTrees; i++)
    {                  
      int 
	bCount = 0;          
     
      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE); 
      assert(tr->mxtips == tr->ntips);
      
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			      &bCount, treeVectorLength, FALSE, FALSE);    
      assert(bCount == tr->mxtips - 3);
    }                                     
	     
  switch(tr->bootStopCriterion)
    {
    case FREQUENCY_STOP:
      {
	double 
	  allOut[2],
	  allIn[2];
	
	avg = frequencyCriterion(numberOfTrees, h, &countBetter, bootStopPermutations, adef);	  	  	  
	
	/*printf("%d \t\t\t %f \t\t\t\t %d\n", numberOfTrees, avg, countBetter);*/
	
	allOut[0] = (double)countBetter;
	allOut[1] = avg;

	MPI_Allreduce(allOut, allIn, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	/*printf("%d %f %f\n", processID, allIn[0], allIn[1]);*/

	stop = (((int)allIn[0]) >= adef->fcThreshold && (allIn[1] / ((double)processes)) >= FC_LOWER);

	*pearsonAverage = (allIn[1] / ((double)processes));
      }
      break;
    case MR_STOP:
    case MRE_STOP:
    case MRE_IGN_STOP:
      {
	double 
	  allOut[4],
	  allIn[4];
	
	double 
	  wrf_thresh_avg = 0.0,
	  wrf_avg = 0.0;
	
	avg = wcCriterion(numberOfTrees, h, &countBetter, &wrf_thresh_avg, &wrf_avg, tr, vectorLength, bootStopPermutations, adef);
	
	/*printf("%d %1.2f  %d %f %f\n", numberOfTrees, avg, countBetter, wrf_thresh_avg, wrf_avg);*/

	allOut[0] = (double)countBetter;
	allOut[1] = wrf_thresh_avg;
	allOut[2] = wrf_avg;
	allOut[3] = avg;

	MPI_Allreduce(allOut, allIn, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	
	/*printf("%d %f %f %f\n", processID, allIn[0], allIn[1], allIn[2]);*/
	
	stop = (((int)allIn[0]) >= adef->fcThreshold && (allIn[2] / ((double)processes)) <= (allIn[1] / ((double)processes)));

	*pearsonAverage = (allIn[3] / ((double)processes));
      }
      break;
    default:
      assert(0);
    }
        
  fclose(treeFile); 
  
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h);

  return stop;
}

#endif

boolean bootStop(tree *tr, hashtable *h, int numberOfTrees, double *pearsonAverage, unsigned int **bitVectors, int treeVectorLength, unsigned int vectorLength, analdef *adef)
{
  int 
    n = numberOfTrees + 1,
    bCount = 0;

  assert((FC_SPACING % 2 == 0) && (adef->fcThreshold < adef->bootstopPermutations));
  assert(tr->mxtips == tr->rdta->numsp);

  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, numberOfTrees, BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			  &bCount, treeVectorLength, FALSE, FALSE);
  assert(bCount == tr->mxtips - 3); 

  if((n > START_BSTOP_TEST) && (n % FC_SPACING == 0))
    {     
      int countBetter = 0;

      switch(tr->bootStopCriterion)
	{
	case FREQUENCY_STOP:
	  *pearsonAverage = frequencyCriterion(n, h, &countBetter, adef->bootstopPermutations, adef);	  	        	       

	  if(countBetter >= adef->fcThreshold && *pearsonAverage >= FC_LOWER)
	    return TRUE;
	  else
	    return FALSE;
	  break;
	case MR_STOP:
	case MRE_STOP:
	case MRE_IGN_STOP:
	 {	   
	   double 
	     wrf_thresh_avg = 0.0,
	     wrf_avg = 0.0;
	   
	   *pearsonAverage = wcCriterion(n, h, &countBetter, &wrf_thresh_avg, &wrf_avg, tr, vectorLength, adef->bootstopPermutations, adef);
	  
	   if(countBetter >= adef->fcThreshold && wrf_avg <= wrf_thresh_avg)
	     return TRUE;
	   else
	     return FALSE;
	 } 
	default:
	  assert(0);
	}
    }
  else
    return FALSE;
}




/* consensus stuff */

boolean compatible(entry* e1, entry* e2, unsigned int bvlen)
{
  unsigned int i;

  unsigned int 
    *A = e1->bitVector,
    *C = e2->bitVector;
  
  for(i = 0; i < bvlen; i++)        
    if(A[i] & C[i])
      break;
          
  if(i == bvlen)
    return TRUE;
  
  for(i = 0; i < bvlen; i++)
    if(A[i] & ~C[i])
      break;
   
  if(i == bvlen)  
    return TRUE;  
  
  for(i = 0; i < bvlen; i++)
    if(~A[i] & C[i])
      break;
  
  if(i == bvlen)
    return TRUE;  
  else
    return FALSE;
}



static int sortByWeight(const void *a, const void *b, int which)
{
  /* recall, we want to sort descending, instead of ascending */
     
  int 
    ca,
    cb;
    
  ca = ((*((entry **)a))->supportFromTreeset)[which];
  cb = ((*((entry **)b))->supportFromTreeset)[which];
  
  if (ca == cb) 
    return 0;
  
  return ((ca<cb)?1:-1);
}

static int sortByIndex(const void *a, const void *b)
{
  if ( (*((entry **)a)) == (*((entry **)b)) ) return 0;
  return (( (*((entry **)a)) < (*((entry **)b)) )?-1:1);
}

static int _sortByWeight0(const void *a, const void *b)
{
  return sortByWeight(a,b,0);
}

static int _sortByWeight1(const void *a, const void *b)
{
  return sortByWeight(a,b,1);
}

boolean issubset(unsigned int* bipA, unsigned int* bipB, unsigned int vectorLen, unsigned int firstIndex)
{
  unsigned int 
    i;
  
  for(i = firstIndex; i < vectorLen; i++)
    if((bipA[i] & bipB[i]) != bipA[i])    
      return FALSE;
        
  return TRUE; 
}


static void mre(hashtable *h, boolean icp, entry*** sbi, int* len, int which, int n, unsigned int vectorLength, boolean sortp, tree *tr, boolean bootStopping)
{
  entry 
    **sbw;
  
  unsigned int   
    i = 0,
    j = 0;
  
  sbw = (entry **) rax_calloc(h->entryCount, sizeof(entry *));
   
  for(i = 0; i < h->tableSize; i++) /* copy hashtable h to list sbw */
    {		
      if(h->table[i] != NULL)
	{
	  entry 
	    *e = h->table[i];
	  
	  do
	    {
	      sbw[j] = e;
	      j++;
	      e = e->next;
	    }
	  
	  while(e != NULL);
	}
    }

  assert(h->entryCount == j);

  if(which == 0)  		/* sort the sbw list */
    qsort(sbw, h->entryCount, sizeof(entry *), _sortByWeight0);
  else
    qsort(sbw, h->entryCount, sizeof(entry *), _sortByWeight1);

  *sbi = (entry **)rax_calloc(n - 3, sizeof(entry *));

  *len = 0;

  if(icp == FALSE)
    {
      
#ifdef _USE_PTHREADS
      /*
	We only deploy the parallel version of MRE when not using it 
	in conjunction with bootstopping for the time being.
	When bootstopping it is probably easier and more efficient to
	parallelize over the permutations 
      */

      if(!bootStopping)
	{
	  //printf("Parallel region \n" );

	  tr->h = h;
	  NumberOfJobs = tr->h->entryCount;
	  tr->sectionEnd = MIN(NumberOfJobs, NumberOfThreads * MRE_MIN_AMOUNT_JOBS_PER_THREAD); //NumberOfThreads * MRE_MIN_AMOUNT_JOBS_PER_THREAD;
	  tr->len = len;
	  tr->sbi = (*sbi);
	  tr->maxBips = n - 3;	
	  tr->recommendedAmountJobs = 1;
	  tr->bitVectorLength = vectorLength;
	  tr->sbw = sbw;
	  tr->entriesOfSection = tr->sbw;
	  tr->bipStatus = (int*)rax_calloc(tr->sectionEnd, sizeof(int));
	  tr->bipStatusLen = tr->sectionEnd;
	  masterBarrier(THREAD_MRE_COMPUTE, tr);
	}
      else
#endif
	{
	  for(i = 0; i < h->entryCount && (*len) < n-3; i++)
	    {			
	      boolean 
		compatflag = TRUE;
	      
	      entry 
		*currentEntry = sbw[i];	  	     
	      
	      assert(*len < n-3);
	      
	      if(currentEntry->supportFromTreeset[which] <= ((unsigned int)tr->mr_thresh))
		{
		  int k;
		  
		  for(k = (*len); k > 0; k--)
		    {
		      if( ! compatible((*sbi)[k-1], currentEntry, vectorLength))
			{
			  compatflag = FALSE;	    
			  break;
			}
		    }
		}
	      
	      if(compatflag)
		{	      
		  (*sbi)[*len] = sbw[i];
		  (*len)++;
		}         
	    }
	}
    }
  else 
    {      
      for(i = 0; i < (unsigned int)(n-3); i++)
	{
	  (*sbi)[i] = sbw[i];
	  (*len)++;
	}
    }


  rax_free(sbw);

  if (sortp == TRUE)
    qsort(*sbi, (*len), sizeof(entry *), sortByIndex);

  return;
}


/* if we encounter the first bits that are set, then we can determine,
   whether bip a is a subset of bip b.  We already know, that A has
   more bits set than B and that both bips are compatible to each
   other. Thus, if A & B is true (and A contains bits), then A MUST be
   a proper subset of B (given the setting). */

/* check different versions of this ! */




static int sortByAmountTips(const void *a, const void *b)
{				
  entry 
    *A = (*(entry **)a),
    *B = (*(entry **)b);
  
  if((unsigned int)A->amountTips == (unsigned int)B->amountTips)
    return 0; 

  return (((unsigned int)A->amountTips < (unsigned int)B->amountTips) ?  -1 : 1); 
}


/******* IC function *******************/


static void calculateIC(tree *tr, hashtable *h, unsigned int *bitVector, unsigned int vectorLength, int trees, unsigned int supportedBips, double *ic, double *icAll, boolean verboseIC, int counter)
{
  unsigned int
    maxCounter = 0,
    *maxima = (unsigned int *)rax_calloc(h->entryCount, sizeof(unsigned int)),
    **maximaBitVectors = (unsigned int **)rax_calloc(h->entryCount, sizeof(unsigned int *)),
    numberOfTrees = (unsigned int)trees;

  *ic = 0.0,
  *icAll = 0.0;

  //if the support is 100% we don't need to consider any conflicting bipartitions and can save some time

  if(supportedBips == numberOfTrees)
    {
      *ic = 1.0;
      *icAll = 1.0;
      
      if(verboseIC)
	printFullySupportedSplit(tr, bitVector, numberOfTrees);
    }
  else
    {
      unsigned int 
	incompatibleBipartitions;

      //search conflicting bipartitions 

      incompatibleBipartitions = countIncompatibleBipartitions(bitVector, h, vectorLength, maxima, &maxCounter, FALSE, numberOfTrees, maximaBitVectors);           

      if(incompatibleBipartitions == 0)
	{
	   *ic = 1.0;
	   *icAll = 1.0;
      
	   printBothOpen("WARNING, returning an IC score of 1.0, while only %d out of %d trees support the current bipartition\n", supportedBips, numberOfTrees);
	   printBothOpen("The IC is still 1.0, but some input trees do not contain information about this bipartition!\n\n");

	   if(verboseIC)
	     printFullySupportedSplit(tr, bitVector, numberOfTrees);
	}
      else
	{
	  //make sure that the sum of raw supports is not higher than the number of trees 
	  
	  assert(supportedBips + maxima[0] <= numberOfTrees);
	  
	  *ic    = computeIC_Value(supportedBips, maxima, numberOfTrees, maxCounter, FALSE, FALSE);
	  *icAll = computeIC_Value(supportedBips, maxima, numberOfTrees, maxCounter, TRUE, FALSE);  
	  
	  if(verboseIC)		      
	    printVerboseIC(tr, supportedBips, bitVector, maxCounter, maxima, maximaBitVectors, numberOfTrees, counter, *ic);
	}
    }

  //printf("IC %f %f IC-all %f %fmaxima: %u\n", ic, _ic, icAll, _icAll, maxCounter);

  rax_free(maxima);
  rax_free(maximaBitVectors);
}

/******* IC function end ***************/

static void printBipsRecursive(tree *tr, FILE *outf, int consensusBipLen, entry **consensusBips, int numberOfTrees, 
			       int currentBipIdx, IdList **listOfDirectChildren, int bitVectorLength, int numTips, 
			       char **nameList, entry *currentBip, boolean *printed, boolean topLevel, unsigned int *printCounter, hashtable 
			       *h, boolean computeIC, double *tc, double *tcAll, boolean verboseIC)
{
  IdList 
    *idx; 
  
  int 
    i;
  
  unsigned int 
    *currentBitVector = (unsigned int*)rax_calloc(bitVectorLength,  sizeof(unsigned int));  

  /* open bip */
  if(*printed)
    fprintf(outf, ",");
  *printed = FALSE;
    
  if(!topLevel)
    fprintf(outf, "(");   

  /* determine tips that are not in sub bipartitions */
  for(i = 0; i < bitVectorLength; i++)
    {  	  
      idx = listOfDirectChildren[currentBipIdx]; 
      currentBitVector[i] = currentBip->bitVector[i];
      
      while(idx)
	{
	  currentBitVector[i] = currentBitVector[i] & ~ consensusBips[idx->value]->bitVector[i]; 
	  idx = idx->next;
	}
    }

  /* print out those tips that are direct leafs of the current bip */
  for(i = 0; i < numTips; i++)
    {
      if(currentBitVector[i / MASK_LENGTH] & mask32[i % MASK_LENGTH])
	{
	  if(*printed){fprintf(outf, ",");};
	  fprintf(outf, "%s", nameList[i+1]);    
	  *printed = TRUE;
	}
    }

  /* process all sub bips */    
  idx = listOfDirectChildren[currentBipIdx]; 
  while(idx)
    {
    
      if(*printed)
	{
	  fprintf(outf, ",");
	  *printed = FALSE;
	} 
      
      printBipsRecursive(tr, outf, consensusBipLen, consensusBips, numberOfTrees, 
			 idx->value, listOfDirectChildren, bitVectorLength, numTips, nameList, 
			 consensusBips[idx->value], printed, FALSE, printCounter, h, computeIC, tc, tcAll, verboseIC);
      *printed  = TRUE;
      idx = idx->next; 
    }

  /* close the bipartition */
  if(currentBipIdx != consensusBipLen)
    {
      if(computeIC)
	{
	  double 
	    ic,
	    icAll;
	  
	  calculateIC(tr, h, currentBip->bitVector, bitVectorLength, numberOfTrees, currentBip->supportFromTreeset[0], &ic, &icAll, verboseIC, *printCounter);

	  *tc    += ic;
	  *tcAll += icAll;

	  fprintf(outf,"):1.0[%1.2f,%1.2f]", ic, icAll);
	}
      else
	{
	  double 
	    support = ((double)(currentBip->supportFromTreeset[0])) / ((double) (numberOfTrees));
	  
	  int 
	    branchLabel = (int)(0.5 + support * 100.0);
	  
	  fprintf(outf,"):1.0[%d]", branchLabel);
	}
      
      *printCounter = *printCounter + 1;
    }
  else
    fprintf(outf, ");\n");
  
  rax_free(currentBitVector); 
}




static void printSortedBips(entry **consensusBips, const int consensusBipLen, const int numTips, const unsigned int vectorLen, 
			    const int numberOfTrees, FILE *outf, char **nameList , tree *tr, unsigned int *printCounter, hashtable *h, boolean computeIC, boolean verboseIC)
{
  int 
    i;

  double
    tc = 0.0,
    tcAll = 0.0;

  IdList 
    **listOfDirectChildren = (IdList**) rax_calloc(consensusBipLen + 1, sizeof(IdList*)); /* reserve one more: the last one is the bip with all species */
  
  boolean 
    *hasAncestor = (boolean*) rax_calloc(consensusBipLen, sizeof(boolean)),
    *printed = (boolean*)rax_calloc(1, sizeof(boolean));
  
  entry 
    *topBip; 

  /* sort the consensusBips by the amount of tips they contain */
  
  for( i = 0; i < consensusBipLen; i++)
    consensusBips[i]->amountTips = genericBitCount(consensusBips[i]->bitVector, vectorLen);  

  qsort(consensusBips, consensusBipLen, sizeof(entry *), &sortByAmountTips);

  /* create an artificial entry for the top */
  topBip = (entry *)rax_malloc(sizeof(entry));
  topBip->bitVector = rax_calloc(sizeof(unsigned int), vectorLen);  
  
  for(i = 1; i < numTips ; i++)
    topBip->bitVector[i / MASK_LENGTH] |= mask32[i % MASK_LENGTH];  

 

  /* find the parent of each bip (in the tree they represent) and construct some kind of hashtable this way */
#ifdef _USE_PTHREADS

  //printf("Parallel region 2\n");

  NumberOfJobs = consensusBipLen;
  tr->consensusBipLen = consensusBipLen; 
  tr->consensusBips = consensusBips;
  tr->mxtips = numTips; /* don't need this ? */
  tr->hasAncestor = hasAncestor; 
  tr->listOfDirectChildren = listOfDirectChildren;
  tr->bitVectorLength = vectorLen;   
  tr->mutexesForHashing = (pthread_mutex_t**) rax_malloc(consensusBipLen * sizeof(pthread_mutex_t*));  
  
  for(i = 0; i < consensusBipLen; i++)
    {
      tr->mutexesForHashing[i] = (pthread_mutex_t*) rax_malloc(sizeof(pthread_mutex_t));
      pthread_mutex_init(tr->mutexesForHashing[i], (pthread_mutexattr_t *)NULL);
    }
  
  masterBarrier(THREAD_PREPARE_BIPS_FOR_PRINT, tr);

  /* cleanup */
  for(i = 0; i < consensusBipLen; i++)
    rax_free(tr->mutexesForHashing[i]);
  rax_free(tr->mutexesForHashing);

  /* restore the old variables - necessary? */
  
  hasAncestor = tr->hasAncestor;
  listOfDirectChildren = tr->listOfDirectChildren; 
#else 
  {
    int 
      j; 

    for(i = 0; i < consensusBipLen; i++)
      {
	entry 
	  *bipA = consensusBips[i]; 
	
	/* find first index  */
	unsigned int 
	  firstIndex = 0;
	
	while(firstIndex < vectorLen && bipA->bitVector[firstIndex] == 0 )
	  firstIndex++;
	
	for(j = i + 1; j < consensusBipLen; j++)
	  {		    
	    if((unsigned int)consensusBips[i]->amountTips < (unsigned int)consensusBips[j]->amountTips
	       && issubset(consensusBips[i]->bitVector, consensusBips[j]->bitVector, vectorLen, firstIndex))
	      { 	      
		IdList
		  *elem = (IdList*) rax_calloc(1,sizeof(IdList)); 
		elem->value = i; 
		elem->next = listOfDirectChildren[j];
		listOfDirectChildren[j] = elem;
		hasAncestor[i] = TRUE;
		break;
	      }
	  }
      }    
  }
#endif

  /****************************************************************/
  /* print the bips during a DFS search on the ancestor hashtable */
  /****************************************************************/

  /* insert these toplevel bips into the last field of the array */
  for(i = 0; i < consensusBipLen; i++)
    if( ! hasAncestor[i])
      {
	IdList
	  *elem  = rax_calloc(1,sizeof(IdList)); 
	elem->value = i; 
	elem->next = listOfDirectChildren[consensusBipLen];  
	listOfDirectChildren[consensusBipLen] = elem;       
      }
  
  /* start dfs search at the top level */
  printBipsRecursive(tr, outf, 
		     consensusBipLen, consensusBips, 
		     numberOfTrees, consensusBipLen,
		     listOfDirectChildren, vectorLen, 
		     numTips, nameList, 
		     topBip, printed, TRUE, printCounter, h, computeIC, &tc, &tcAll, verboseIC);

  if(computeIC)
    {
      double
	rtcAll = tcAll / (double)(tr->mxtips - 3),
	rtc    = tc    / (double)(tr->mxtips - 3);
      
      printBothOpen("Tree certainty for this tree: %f\n", tc);
      
      /* Leonida: for consensus trees I also calculate the relative tree certainty by dividing 
	 the tc by the total number of bipartitions of a tree with n (tr->mxtips) taxa 
	 to penalize the consensi for potentially being unresolved */
	 
      printBothOpen("Relative tree certainty for this tree: %f\n\n", rtc);

      printBothOpen("Tree certainty including all conflicting bipartitions (TCA) for this tree: %f\n", tcAll);
      printBothOpen("Relative tree certainty including all conflicting bipartitions (TCA) for this tree: %f\n\n", rtcAll);
    }

  rax_free(topBip->bitVector);
  rax_free(topBip);
  rax_free(printed);
  rax_free(hasAncestor);

  for( i = 0; i < consensusBipLen + 1; ++i)
    {
      IdList
	*iter = listOfDirectChildren[i]; 
      
      while(iter != NULL)
	{
	  IdList *nxt = iter->next; 
	  rax_free(iter); 
	  iter = nxt; 
	}
    }
  
  rax_free(listOfDirectChildren);
}




void computeConsensusOnly(tree *tr, char *treeSetFileName, analdef *adef, boolean computeIC)
{        
  hashtable 
    *h = initHashTable(tr->mxtips * FC_INIT * 10);

  hashNumberType
    entries = 0;

  unsigned int  
    printCounter  = 0,
    numberOfTrees = 0, 
    i, 
    j,     
    treeVectorLength, 
    vectorLength;

  int
    consensusBipsLen = 0;  

  unsigned int    
    **bitVectors = initBitVector(tr, &vectorLength);

  entry 
    **consensusBips;

  char 
    someChar[1024],
    consensusFileName[1024];   
  
  FILE 
    *outf,
    *treeFile = getNumberOfTrees(tr, treeSetFileName, adef);

  numberOfTrees = tr->numberOfTrees; 

  checkTreeNumber(numberOfTrees, treeSetFileName);

  tr->mr_thresh = ((double)numberOfTrees / 2.0);   

  assert(sizeof(unsigned char) == 1);
 
  treeVectorLength = GET_BITVECTOR_LENGTH(numberOfTrees);

  /* read the trees and process the bipartitions */ 

  //modified for multifurcations !

  {
    tree
      *inputTree = (tree *)rax_malloc(sizeof(tree));

    allocateMultifurcations(tr, inputTree);

    for(i = 1; i <= numberOfTrees; i++)
      {                  
	int 
	  numberOfSplits,
	  bCount = 0;
	
	numberOfSplits = readMultifurcatingTree(treeFile, inputTree, adef, FALSE);
	//treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE);               
	
	assert(inputTree->mxtips == inputTree->ntips);
	
	bitVectorInitravSpecial(bitVectors, inputTree->nodep[1]->back, inputTree->mxtips, vectorLength, h, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
				&bCount, treeVectorLength, FALSE, FALSE);
	
	assert(bCount == numberOfSplits);                     
      }

    freeMultifurcations(inputTree);
    rax_free(inputTree);     
  }
  
  fclose(treeFile);
  
  if(tr->consensusType == MR_CONSENSUS || tr->consensusType == STRICT_CONSENSUS || tr->consensusType == USER_DEFINED)
    {
      consensusBips = (entry **)rax_calloc(tr->mxtips - 3, sizeof(entry *));
      consensusBipsLen = 0;   
    }
  
  for(j = 0; j < (unsigned int)h->tableSize; j++) /* determine support of the bips */
    {		
      if(h->table[j] != NULL)
	{
	  entry *e = h->table[j];
	  
	  do
	    {	
	      unsigned int 
		cnt = genericBitCount(e->treeVector, treeVectorLength);

		if((tr->consensusType == MR_CONSENSUS     && cnt > (unsigned int)tr->mr_thresh) || 
		 (tr->consensusType == STRICT_CONSENSUS && cnt == numberOfTrees) ||
		 (tr->consensusType ==  USER_DEFINED    && cnt > (numberOfTrees * tr->consensusUserThreshold) / 100))
		{
		  consensusBips[consensusBipsLen] = e;
		  consensusBipsLen++;
		}

	      e->supportFromTreeset[0] = cnt;
	      e = e->next;
	      entries++;
	    }
	  while(e != NULL);		
	}	  	        
    }	
  
  assert(h->entryCount == entries);
  
  if(tr->consensusType == MR_CONSENSUS || tr->consensusType == STRICT_CONSENSUS || tr->consensusType == USER_DEFINED)
    assert(consensusBipsLen <= (tr->mxtips - 3));
   
  if(tr->consensusType == MRE_CONSENSUS)   
    mre(h, FALSE, &consensusBips, &consensusBipsLen, 0, tr->mxtips, vectorLength, FALSE , tr, FALSE);  

  /* printf("Bips NEW %d\n", consensusBipsLen); */

  strcpy(consensusFileName,         workdir);  
  
  switch(tr->consensusType)
    {
    case MR_CONSENSUS:
      if(computeIC)
	strcat(consensusFileName,         "RAxML_MajorityRuleConsensusTree_IC.");
      else
	strcat(consensusFileName,         "RAxML_MajorityRuleConsensusTree.");
      break;
    case MRE_CONSENSUS:
      if(computeIC)
	strcat(consensusFileName,         "RAxML_MajorityRuleExtendedConsensusTree_IC.");
      else
	strcat(consensusFileName,         "RAxML_MajorityRuleExtendedConsensusTree.");
      break;
    case STRICT_CONSENSUS:
      assert(!computeIC);
      strcat(consensusFileName,         "RAxML_StrictConsensusTree.");
      break;
    case USER_DEFINED :        
      if(computeIC)
	sprintf(someChar,         "RAxML_Threshold-%d-ConsensusTree_IC.", tr->consensusUserThreshold);
      else
	sprintf(someChar,         "RAxML_Threshold-%d-ConsensusTree.", tr->consensusUserThreshold);
      strcat(consensusFileName, someChar);       
      break; 
    default:
      assert(0);
    }
  
  strcat(consensusFileName,         run_id);

  outf = myfopen(consensusFileName, "wb");

  fprintf(outf, "(%s,", tr->nameList[1]);
  
  if(computeIC)
    {
      if(adef->verboseIC)
	printVerboseTaxonNames(tr);
      printSortedBips(consensusBips, consensusBipsLen, tr->mxtips, vectorLength, numberOfTrees, outf, tr->nameList, tr, &printCounter, h, computeIC, adef->verboseIC);
    }
  else
    printSortedBips(consensusBips, consensusBipsLen, tr->mxtips, vectorLength, numberOfTrees, outf, tr->nameList, tr, &printCounter, h, computeIC, FALSE);

  assert(printCounter ==  (unsigned int)consensusBipsLen);

  /* ????? fprintf(outf, ");\n"); */
  
  fclose(outf);
  
  if(adef->verboseIC && computeIC)
    printBothOpen("Verbose PHYLIP-style formatted bipartition information written to file: %s\n\n",  verboseSplitsFileName);
  
  switch(tr->consensusType)
    {
    case MR_CONSENSUS:
      if(computeIC)	
	printBothOpen("RAxML Majority Rule consensus tree with IC values written to file: %s\n\n", consensusFileName);	         
      else
	printBothOpen("RAxML Majority Rule consensus tree written to file: %s\n", consensusFileName);
      break;
    case MRE_CONSENSUS:
      if(computeIC)
	printBothOpen("RAxML extended Majority Rule consensus tree with IC values written to file: %s\n", consensusFileName);
      else
	printBothOpen("RAxML extended Majority Rule consensus tree written to file: %s\n", consensusFileName);
      break;
    case STRICT_CONSENSUS:
      printBothOpen("RAxML strict consensus tree written to file: %s\n", consensusFileName);
      break;
    case USER_DEFINED: 
      if(computeIC)
	printBothOpen("RAxML consensus tree with threshold %d with IC values written to file: %s\n", tr->consensusUserThreshold,  consensusFileName);
      else
	printBothOpen("RAxML consensus tree with threshold %d written to file: %s\n", tr->consensusUserThreshold,  consensusFileName);
      break;
    default:
      assert(0);
    }
  
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h);
  rax_free(consensusBips);

  exit(0);   
}


/***************************** start of TC/IC stuff for partial trees **********************************/

//function to obtain a bit mask of the present partial tree 
//the bit mask has at least as many bits (may also have padding bits!) 
//as the reference tree has taxa, taxa that are present are set to 1
//taxa that are absent are set to 0

static void setupMask(unsigned int *smallTreeMask, nodeptr p, int numsp)
{
  //if this is a tip set the bit to 1
  //note that the taxon numbering startes at 1 in RAxML, but we use 
  //the bit-mask as 0-based array 
  
  if(isTip(p->number, numsp))
    smallTreeMask[(p->number - 1) / MASK_LENGTH] |= mask32[(p->number - 1) % MASK_LENGTH];
  else
    {    
      nodeptr 
	q = p->next;

      /* I had to change this function to account for mult-furcating trees.
	 In this case an inner node can have more than 3 cyclically linked 
	 elements, because there might be more than 3 outgoing branches 
	 from an inner node */

      while(q != p)
	{
	  setupMask(smallTreeMask, q->back, numsp);
	  q = q->next;
	}      
    }
}

//dynamic array data structure to store the tree masks from the various trees 
//TODO move the struct below to axml.h at some point
typedef struct 
{
  //number of slots used
  size_t entries;
  
  //number of slots available
  size_t length;
  
  //array with pointers to bit masks 
  unsigned int **bitMasks;
} taxonMaskArray;

//function to initialize the dynamic array for storing bit masks

static taxonMaskArray* initTaxonMaskArray(unsigned int vLength)
{
  //initialize the structure itself 
  taxonMaskArray 
    *b = rax_malloc(sizeof(taxonMaskArray));
    
  //make sure that the number of 32-bit integers that we need 
  //for storing the bits has been set
  assert(vLength != 0);

  //set number of available slots to 100
  b->length = 100;

  //set number of used slots to 0
  b->entries = 0;

  //allocate the array to hold the bit masks 
  b->bitMasks = (unsigned int **)rax_malloc(sizeof(unsigned int *) * (size_t)b->length);

  return b;
}

//function to add a taxon bit mask to the dynamic array 

static void addTaxonMask(unsigned int vLength, unsigned int *mask, taxonMaskArray *b)
{
  //is the bit mask already stored?
  boolean 
    found = FALSE;
  
  unsigned int 
    i, 
    j;

  //loop over entries/used slots to see if the bit mask is already stored.

  for(i = 0; i < b->entries; i++)
    {
      //check mask to add and stored mask i for identity 
      for(j = 0; j < vLength; j++)
	if(mask[j] != b->bitMasks[i][j])	 
	  break;

      //if all 32-bit integers are identical, the bit mask is already stored in our data structure
      if(j == vLength)
	{
	  found = TRUE;
	  break;
	}
    }
  
  //if we have found the bit mask 
  //we can just exit
  if(found)
    return;
  else
    {
      //need to add the bit mask 

      //if there are still free slots available we can just add the new bit mask 
      if(b->entries < b->length)
	{
	  //allocate space for the bit mask 
	  b->bitMasks[b->entries] = (unsigned int *)rax_malloc(sizeof(unsigned int) * (size_t)vLength);

	  //copy the bit mask 
	  memcpy(b->bitMasks[b->entries], mask, sizeof(unsigned int) * (size_t)vLength);

	  //increment the number of used slots 
	  b->entries = b->entries + 1;
	}
      else
	{
	  //not enough space, increment the number of slots 
	  b->length = b->length * 2;

	  //reallocate the space for pointers to the bit masks 
	  b->bitMasks = (unsigned int **)rax_realloc(b->bitMasks, sizeof(unsigned int *) * (size_t)b->length, FALSE);

	  //now assign space for the new bit mask 
	  b->bitMasks[b->entries] = (unsigned int *)rax_malloc(sizeof(unsigned int) * (size_t)vLength);

	  //copy the bit mask
	  memcpy(b->bitMasks[b->entries], mask, sizeof(unsigned int) * (size_t)vLength);

	  //increment the number of used slots 
	  b->entries = b->entries + 1;
	}
    }
}

//function to free bit mask array 

static void freeTaxonMask(taxonMaskArray *b)
{
  unsigned int 
    i;
  
  //free all bit masks stored at used slots 
  for(i = 0; i <  b->entries; i++)
    rax_free(b->bitMasks[i]);
    
  //free bit mask array 
  rax_free(b->bitMasks);
    
  return;
}

//TODO add struct below to axml.h at some point

typedef struct 
{
  //number of slots used
  size_t entries;
  
  //number of slots available
  size_t length;
  
  //array with pointers to taxon masks 
  unsigned int **taxonMasks;
  
  //array with pointers to bit masks 
  unsigned int **bitMasks;
  
  //array with pointers to frequencies 
  double *frequencies;
  
} referenceMaskArray;


//insert bipartition into hash table

static void insertHashPartial(unsigned int *bitVector, hashtable *h, unsigned int vectorLength, hashNumberType position, unsigned int *mask)
{     
  //there is already an entry in the has table at this position,
  //so check

  if(h->table[position] != NULL)
    {
      entry 
	*e = h->table[position];     

      do
	{	 
	  unsigned int 
	    i;

	  boolean 
	    maskIdentical = FALSE,
	    bitVectorIdentical = FALSE;
	  
	  //check for identity of bipartitions 

	  for(i = 0; i < vectorLength; i++)
	    if(bitVector[i] != e->bitVector[i])
	      break;
	  
	  if(i == vectorLength)	 	    	     
	    bitVectorIdentical = TRUE;	   	    

	  //check for identity of taxon masks

	  for(i = 0; i < vectorLength; i++)
	    if(mask[i] != e->taxonMask[i])
	      break;
	  
	  if(i == vectorLength)	 	    	     
	    maskIdentical = TRUE;

	  //if the bipartition vector and the taxon mask are identical
	  //increment the frequency of this (taxonMask,bipartition) pair and exit 

	  if(maskIdentical && bitVectorIdentical)
	    {
	      e->supportFromTreeset[0] = e->supportFromTreeset[0] + 1;
	      return;
	    }
	  
	  //otherwise keep searching 

	  e = e->next;
	}
      while(e != (entry*)NULL); 

      //nothing found, so we need to create a new entry by chaining 

      e = initEntry(); 
            
      //allocate space and copy the taxon mask and the bipartition vector

      e->bitVector = (unsigned int*)rax_malloc((size_t)vectorLength * sizeof(unsigned int)); 
      e->taxonMask = (unsigned int*)rax_malloc((size_t)vectorLength * sizeof(unsigned int));

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * (size_t)vectorLength);
      memcpy(e->taxonMask, mask,      sizeof(unsigned int) * (size_t)vectorLength);
      
      //set frequency of this  (taxonMask,bipartition) pair to 1

      e->supportFromTreeset[0] = 1;

      //update the linked list

      e->next = h->table[position];
      h->table[position] = e;          
    }
  else
    {
      //no entry at the current has index, so create one 

      entry 
	*e = initEntry(); 
             
      //allocate space and copy the taxon mask and the bipartition bit vector 

      e->bitVector = (unsigned int*)rax_malloc((size_t)vectorLength * sizeof(unsigned int));      
      e->taxonMask = (unsigned int*)rax_malloc((size_t)vectorLength * sizeof(unsigned int)); 

      memcpy(e->bitVector, bitVector, sizeof(unsigned int) * (size_t)vectorLength);  
      memcpy(e->taxonMask, mask,      sizeof(unsigned int) * (size_t)vectorLength);
      
      //set frequency of occurence of this pair to 1

      e->supportFromTreeset[0] = 1;

      //set the hash table index to pint to this entry

      h->table[position] = e;
    }

  //increment hash table entry counter 

  h->entryCount =  h->entryCount + 1;
}

//function to traverse the partial trees, extract their non-trivial bipartitions and store them, together with the 
//taxon bit mask in the hash table 

static void bitVectorInitravPartial(tree *tr, unsigned int **bitVectors, unsigned int *taxonMask, nodeptr p,  unsigned int vectorLength, hashtable *h, int *countBranches, int numsp)
{
  //if this is a tip it's a trivial bipartition do nothing 
  if(isTip(p->number, numsp))
    return;
  else
    {
      nodeptr q = p->next;          

      //handle all other nodes recursively 
      do 
	{
	  bitVectorInitravPartial(tr, 
				  bitVectors, 
				  taxonMask, 
				  q->back, vectorLength, h, countBranches, numsp);
	  q = q->next;
	}
      while(q != p);
           

      //obtain the bipartition of the present branch in  very similar way we 
      //compute likelihoods, this is the dedicated function for multifurcating trees 
      newviewBipartitionsMultifurcating(bitVectors, p, numsp, vectorLength);
      
      //make sure that the node value (the bipartition for branch p <-> p->back is there 
      
      assert(p->x);
      
      //if the back pointer is not a tip insert the bipartition into the hash table 

      if(!(isTip(p->back->number, numsp)))
	{
	  unsigned int 
	    *toInsert  = bitVectors[p->number],
	    *correctedVector = (unsigned int*)rax_malloc(sizeof(unsigned int) * (size_t)vectorLength);

	  hashNumberType 
	    position;

	  //need to make a copy of the bit-vector since we may have to reverse it 

	  memcpy(correctedVector, toInsert, sizeof(unsigned int) * (size_t)vectorLength);	 

	  //in RAxML to store the bit vectors representing bipartitions in a canonical form
	  //we always have the restriction that the very first bit need to be set to 0
	  //if this is not the case we need to complement the bit vector 

	  if(correctedVector[0] & 1)
	    {
	      unsigned int
		j,
		count1 = 0,
		count2 = 0;
	      
	      //count the bits in non-reverted vector 
	      for(j = 0; j < vectorLength; j++)
		count1 += BIT_COUNT(correctedVector[j]);  

	      //revert the vector by computing the bit-wise complement,
	      //note that we also need to apply the taxon mask again 
	      //to make sure that the complemented vector only contains 
	      //taxa that are in the tree 
	      for(j = 0; j < vectorLength; j++)		
		correctedVector[j] = (~toInsert[j]) & taxonMask [j];		

	      //count the number of taxa in thi sbipartition
	      for(j = 0; j < vectorLength; j++)
		count2 += BIT_COUNT(correctedVector[j]);

	      //printf("reverting %d\n", tr->start->number);

	      //make sure that everything is correct
	      assert(count1 + count2 == (unsigned int)tr->ntips);
	    }
	 
	 
	  //check that the first biot in the bit vector representing the bipartition is 
	  //not set 
	 
	  assert(!(correctedVector[0] & 1));	 

	  //compute hash index 

	  position = oat_hash((unsigned char *)correctedVector, sizeof(unsigned int) * vectorLength);
	  position = position % h->tableSize;
	  
	  //insert into hash table 

	  insertHashPartial(correctedVector, h, vectorLength, position, taxonMask);
	  
	  //increment the non-trivial bipartition counter 
	  *countBranches =  *countBranches + 1;	 

	  //free intermediate storage 
	  rax_free(correctedVector);
	}
      
    }
}







//simply invert the bits in a bitvector (length is actual length in bits and may be larger than MASK_LENGTH)
static void invertBits(unsigned int * v, int length)
{
  int 
    i,
    j = ceil(((double)length)/((double)MASK_LENGTH));
  
  for(i = 0; i < j; i++)
    v[i] = ~v[i];
}

//simply set the bit at position i to 0 (i may be larger than MASK_LENGTH)
static void unsetBit(unsigned int *v, int i)
{  
  int
    j = floor(((double)i)/((double)MASK_LENGTH));
  
  v[j] &= ~mask32[i % MASK_LENGTH];  
}

//return value of bit at position i (i may be larger than MASK_LENGTH)
static int readBit(unsigned int * v, int i)
{  
  int
    j = floor(((double)i)/((double)MASK_LENGTH));
  
  if(!(v[j] & mask32[i % MASK_LENGTH]))
    return 0;
  
  return 1;
}

static unsigned int bitCountLong(unsigned int *bitVector, int length)
{
 unsigned int 
   j,
   count   = 0,
   vLength = ceil(((double)length)/((double)MASK_LENGTH));
 
 for(j = 0; j < vLength; j++)
   count += BIT_COUNT(bitVector[j]);
 
 return count;
}


static int sortEntryByMask(const void *a, const void *b)
{       

  entry 
    *ca =(*(entry **) a),
    *cb =(*(entry **) b);
  
  if(ca->covered == TRUE && cb->covered == TRUE) 
    return 0;
  
  if(ca->covered == TRUE)
    return 1;
  
  if(cb->covered == TRUE)
    return -1;
  
  if(ca->amountTips ==  cb->amountTips) 
    return 0;
  
  return (((unsigned int)ca->amountTips < (unsigned int)cb->amountTips) ?  1 : -1); 
}

static int sortEntryByAdjustedSupport(const void *a, const void *b)
{       
  entry 
    *ca =(*(entry **) a),
    *cb =(*(entry **) b);
  
  if(ca->covered == TRUE && cb->covered == TRUE) 
    return 0;
  
  if(ca->covered == TRUE)
    return 1;
  
  if(cb->covered == TRUE)
    return -1;
  
  if(ca->adjustedSupport ==  cb->adjustedSupport) 
    return 0;
  
  return ((ca->adjustedSupport < cb->adjustedSupport) ?  1 : -1); 
}

static int sortEntryByTempSupport(const void *a, const void *b)
{       
  entry 
    *ca =(*(entry **) a),
    *cb =(*(entry **) b);
  
  if(ca->covered == TRUE && cb->covered == TRUE) 
    return 0;
  
  if(ca->covered == TRUE)
    return 1;
  
  if(cb->covered == TRUE)
    return -1;
  
  if(ca->tempSupport ==  cb->tempSupport) 
    return 0;
  
  return ((ca->tempSupport < cb->tempSupport) ?  1 : -1); 
}

static double adjustScore(int big, int small, int support, int mxtips)
{
  unsigned int
	    limit,
	    pow;
  
  unsigned long 
	    result;
  
   assert(big >= small);
   
   limit = 2 * (BIT_COUNT(UINT_MAX));
   
   pow = mxtips - small;
   
   if(pow >= limit)
     result = ULONG_MAX; //TODO implement actual scaling
   else //calculate 2^pow
     result = 1 << pow; 
   
   return ((double)support / (double)result); 
}

//checks whether a full and one partial bipartition are compatible. ie X1 cup Y1 != {}, or X1 cup Y2 != {}, or X2 cup Y1 != {}, or X2 cup Y2 != {}. lenght is actual length in bits.

static boolean areCompatibleBipMask(unsigned int* A,unsigned int* C, unsigned int * mask, int length)
{
  unsigned int 
    i,
    vL = ceil(((double)length)/ ((double)MASK_LENGTH));
  
  for(i = 0; i < vL; i++)
    {        
      if((A[i] & C[i]) & mask[i])
	break;
    }
  
  if(i == vL)
    return TRUE;
  
  for(i = 0; i < vL; i++)
    {
      if((A[i] & ~C[i]) & mask[i])
	break;
    }
  
  if(i == vL)  
    return TRUE;  
  
  for(i = 0; i < vL; i++)
    {
      if((~A[i] & ~C[i]) & mask[i])
	break;
    }
  
  if(i == vL)  
    return TRUE;  
  
  for(i = 0; i < vL; i++)
    {
      if((~A[i] & C[i]) & mask[i])
	break;
    }
  
  if(i == vL)
    return TRUE;  
  else
    return FALSE;
}

//checks whether two partial bipartitions are compatible. ie X1 cup Y1 != {}, or X1 cup Y2 != {}, or X2 cup Y1 != {}, or X2 cup Y2 != {}
static boolean areCompatibleBipMaskMask(unsigned int* A,unsigned int* C, unsigned int * mask1, unsigned int * mask2, int length)
{
  unsigned int 
    i,
    vL = ceil(((double)length)/ ((double)MASK_LENGTH));
  
  for(i = 0; i < vL; i++)
    {        
      if((A[i] & C[i]) & (mask1[i] & mask2[i]))
	break;
    }
  
  if(i == vL)
    return TRUE;
  
  for(i = 0; i < vL; i++)
    {
      if((A[i] & ~C[i]) & (mask1[i] & mask2[i]))
	break;
    }
  
  if(i == vL)  
    return TRUE;  
  
  for(i = 0; i < vL; i++)
    {
      if((~A[i] & ~C[i]) & (mask1[i] & mask2[i]))
	break;
    }
  
  if(i == vL)  
    return TRUE;  
  
  for(i = 0; i < vL; i++)
    {
      if((~A[i] & C[i]) & (mask1[i] & mask2[i]))
	break;
    }
  
  if(i == vL)
    return TRUE;  
  else
    return FALSE;
  
}


static double computeIC_ValueDouble(double refFrequency, double confFrequency, boolean warnNegativeIC)
{  
   double 
      sum = refFrequency + confFrequency,
      relativeRef,
      relativeConf,
      ic;
  
   boolean
     negativeIC = FALSE;    

    
    //neither support for the actual bipartition, nor for the conflicting ones
   if(sum == 0)
     return 0.0;
    
   relativeRef  =  refFrequency  / sum;
   relativeConf =  confFrequency / sum;
   
   // figure out if the competing bipartition is higher support 
   // can happen for MRE and when drawing IC values on best ML tree 
   
   if(confFrequency > refFrequency)
     {
       negativeIC = TRUE;
       
       if(warnNegativeIC)
	 {
	   printBothOpen("\nMax conflicting bipartition frequency: %f is larger than frequency of the included bipartition: %f\n", confFrequency, refFrequency);
	   printBothOpen("This is interesting, but not unexpected when computing extended Majority Rule consensus trees.\n");
	   printBothOpen("Please send an email with the input files and command line\n");
	   printBothOpen("to Alexandros.Stamatakis@gmail.com.\n");
	   printBothOpen("Thank you :-)\n\n");   
	 }      
     }
    
   ic = log(2);
   
   if(refFrequency > 0)
     ic += relativeRef * log(relativeRef);
   
   if(confFrequency > 0)
     ic += relativeConf * log(relativeConf);
      
   //go from log(x) to log_2(x)
   ic /= log(2);
  
   if(negativeIC)
     return (-ic);
   else
     return ic;
}


static double computeICA_ValueDouble(double refFrequency, entry ** conflicts, unsigned int numberOfConflicts, boolean warnNegativeIC)
{  
  //if no conflicting bipartition is present, set score to 1
  if(numberOfConflicts == 0)    
    return 1.0;
  
  {
    unsigned int 	
      n   = 2,
      i;

    double 
      sum    = refFrequency + conflicts[0]->adjustedSupport,
      relSum = refFrequency + conflicts[0]->adjustedSupport,
      relativeRef,
      ica;
  
    boolean
      negativeIC = FALSE;    
   
      
    // figure out if the competing bipartition is higher support 
    // can happen for MRE and when drawing IC values on best ML tree 
    
    if(conflicts[0]->adjustedSupport > refFrequency)
      {
	negativeIC = TRUE;
	
	if(warnNegativeIC)
	  {
	    printBothOpen("\nMax conflicting bipartition frequency: %f is larger than frequency of the included bipartition: %f\n", conflicts[0]->adjustedSupport, refFrequency);
	    printBothOpen("This is interesting, but not unexpected when computing extended Majority Rule consensus trees.\n"); 
	  }      
      }
      
    for(i = 1; i < numberOfConflicts; i++)      
      sum += conflicts[i]->adjustedSupport;       
    
    
    for(i = 1; i < numberOfConflicts; i++)
      {
	//Ensure a threshold of 5% for the bipartitions to be considered in the ica calculation.
	//Include at least one conflicting bipartition if any exist though.
	//(Note, the way we calculate the threshold here is slightly differen from
	//what we do for full gene trees.)
	if(conflicts[i]->adjustedSupport / sum < 0.05)
	  break;
	
	relSum += conflicts[i]->adjustedSupport;
	n++;	
      }
          
    //neither support for the actual bipartition, nor for the conflicting ones
    if(relSum == 0)
      return 0.0;
    
    relativeRef = refFrequency / relSum;
        
    ica = log(n);
    
    if(relativeRef > 0)
      ica += relativeRef * log(relativeRef);

      
    for(i = 0; i < n-1; i++)
      {  
	if(conflicts[i]->adjustedSupport > 0)
          ica += (conflicts[i]->adjustedSupport / relSum) * log(conflicts[i]->adjustedSupport / relSum);      
      }
      
    //go from log(x) to log_n(x)
    ica /= log(n);
  
    if(negativeIC)
      return (-ica);
    else
      return ica;
  }
}

//static double computeICAFromTemp_ValueDouble(double refFrequency, entry ** conflicts, unsigned int numberOfConflicts, boolean warnNegativeIC)
static double computeICAFromTemp_ValueDouble(double refFrequency, entry ** conflicts, unsigned int numberOfConflicts, boolean warnNegativeIC, unsigned int numTrees)
{  
  //if no conflicting bipartition is present, set score to 1
  if(numberOfConflicts == 0)    
    return 1.0;
  
  {
    unsigned int 	
      n   = 2,
      i;

    double 
      sum,
      //sum    = refFrequency + conflicts[0]->tempSupport,
      relSum = refFrequency + conflicts[0]->tempSupport,
      relativeRef,
      ica;
  
    boolean
      negativeIC = FALSE;    
   
      
    // figure out if the competing bipartition is higher support 
    // can happen for MRE and when drawing IC values on best ML tree 
    
    //if(conflicts[0]->adjustedSupport > refFrequency)
    if(conflicts[0]->tempSupport > refFrequency)
      {
	negativeIC = TRUE;
	
	if(warnNegativeIC)
	  {
	   // printBothOpen("\nMax conflicting bipartition frequency: %f is larger than frequency of the included bipartition: %f\n", conflicts[0]->adjustedSupport, refFrequency);
	    printBothOpen("\nMax conflicting bipartition frequency: %f is larger than frequency of the included bipartition: %f\n", conflicts[0]->tempSupport, refFrequency);
	    printBothOpen("This is interesting, but not unexpected when computing extended Majority Rule consensus trees.\n");
	  }      
      }
    
    /*
    for(i = 1; i < numberOfConflicts; i++)
    {
    
     sum += conflicts[i]->tempSupport;  
    }
    */
        
    sum=numTrees;
    
    for(i = 1; i < numberOfConflicts; i++)
      {
	//Ensure a threshold of 5% for the bipartitions to be considered in the ica calculation.
	//Include at least one conflicting bipartition if any exist though.
	//(Note, the way we calculate the threshold here is slightly differen from
	//what we do for full gene trees.)
	if(conflicts[i]->tempSupport / sum < 0.05)
	  break;
	
	relSum += conflicts[i]->tempSupport;
	n++;      
      }
          
    //neither support for the actual bipartition, nor for the conflicting ones
    if(relSum == 0)
      return 0.0;
    
    relativeRef = refFrequency / relSum;
        
    ica = log(n);
    
    if(relativeRef > 0)
      ica += relativeRef * log(relativeRef);
      
    for(i = 0; i < n-1; i++)
      {  
	// if(conflicts[i]->adjustedSupport > 0)
	if(conflicts[i]->tempSupport > 0)
	  //     ica += (conflicts[i]->adjustedSupport / relSum) * log(conflicts[i]->adjustedSupport / relSum);      
	  ica += (conflicts[i]->tempSupport / relSum) * log(conflicts[i]->tempSupport / relSum);      
      }
    
    
    //go from log(x) to log_n(x)
    ica /= log(n);
  
    if(negativeIC)
      return (-ica);
    else
      return ica;
  }
}


//Given some reference bipartition refBip, find all conflicting bipartitions. That is, iteratively add bipartitions that are not compatible with anything already added.


static unsigned int findConflictingBipartitions(entry ** conflictingBips, unsigned int * refBip , entry ** candidates, unsigned int candidateLength, unsigned int mxtips)
{
  unsigned int
    conflicts = 0,
    k,
    j;
  
  for(j = 0; j < candidateLength; j++)
    {	
      //index corresponds to decreasing ordering with respect to potentialFrequencies
      //unsigned int 
      //  id = sortedList[j].id;
      
      //Easy access
      unsigned int 
	*confBip = candidates[j]->bitVector,
	*confMask=candidates[j]->taxonMask;
      
      //To be added, the new bipartition must not be compatible with the reference bipartition..
      if(!areCompatibleBipMask(refBip, confBip, confMask, mxtips))
	{
	  //..nor with any of the already added bipartitions.
	  for(k = 0; k < conflicts; k++)		
	    if(areCompatibleBipMaskMask(conflictingBips[k]->bitVector, confBip, confMask, conflictingBips[k]->taxonMask, mxtips))
	      break;		  
	  
	  //If it was not compatible, add it (i.e. the index) to the list.
	  if(k == conflicts)
	    {
	      conflictingBips[conflicts] = candidates[j];
	      conflicts++;
	    }	      	      		
	  
	}
      
      if(conflicts >= 20)
	break;
    }

  //return length of conflictingBips to be considered.
  return conflicts;
}

/*
 static unsigned int findConflictingBipartitions(unsigned int * conflictingBips, unsigned int * refBip , entry ** candidates, unsigned int candidateLength, unsigned int mxtips)
{
  unsigned int
            conflicts = 0,
            k,
            j;
  
      for(j = 0; j < candidateLength; j++)
	{	
	  //index corresponds to decreasing ordering with respect to potentialFrequencies
	  //unsigned int 
	  //  id = sortedList[j].id;
	  
	  //Easy access
	  unsigned int 
	    *confBip = candidates[j]->bitVector,
	    *confMask=candidates[j]->taxonMask;
		
	  //To be added, the new bipartition must not be compatible with the reference bipartition..
	  if(!areCompatibleBipMask(refBip, confBip, confMask, mxtips))
	    {
	      //..nor with any of the already added bipartitions.
	      for(k = 0; k < conflicts; k++)		
		if(areCompatibleBipMaskMask(candidates[conflictingBips[k]]->bitVector, confBip, confMask, candidates[conflictingBips[k]]->taxonMask, mxtips))
		  break;		  

	      //If it was not compatible, add it (i.e. the index) to the list.
	      if(k == conflicts)
		{
		  conflictingBips[conflicts] = j;
		  conflicts++;
		}	      	      		
		
	    }
	    
	    if(conflicts >= 20)
	      break;
	}
	//return length of conflictingBips to be considered.
	return conflicts;
}
 */

//masks a bipartition and rehashes it to find subpartitions if they exist. length is number of actual bits.
static entry * getSubPart(unsigned int * refBip, unsigned int * mask, hashtable* genesHash, int length)
{  
  unsigned int 
    vLength = ceil((double)length/(double)MASK_LENGTH),
    l, 	      
    *bipCopy = (unsigned int *)rax_malloc(sizeof(unsigned int) * (size_t)vLength); //allocate space for a copy of the reference bipartition because we need to mask it 
    
  hashNumberType 
    position; //hash table index
  
  entry 
    *e; //hash table element
  
  //copy the reference bipartition 
  memcpy(bipCopy, refBip, sizeof(unsigned int) * (size_t)vLength); 
  
  //mask the reference bipartition 
  for(l = 0; l < vLength; l++)
    bipCopy[l] = bipCopy[l] & mask[l];
  
  //check that we have a canonic representation of the masked ref. bip.
  assert(!(bipCopy[0] & 1));
    
  //compute hash table index for the masked reference bipartition 
  position = oat_hash((unsigned char *)bipCopy, sizeof(unsigned int) * vLength);
  position = position % genesHash->tableSize;
    
  e = genesHash->table[position];
    
  if(e != (entry*)NULL)
    {
      //if there is one or more elements here, we need to search the linked list 
      //for a matching bipartition/mask pair 
      //we don't necessarily need to find one !
      
      do 
	{	     
	  boolean 
	    maskIdentical = FALSE,
	    bitVectorIdentical = FALSE;
	
	  //make sure the bip. stored in the has table 
	  //is also in canonic form
	  assert(!(e->bitVector[0] & 1));
	
	  //check for identity of bipartitions 				
	  for(l = 0; l < vLength; l++)
	    if(bipCopy[l] != e->bitVector[l])
	      break;
	  
	  if(l == vLength)	 	    	     
	    bitVectorIdentical = TRUE;	   	    
	  
	  //check for identity of taxon masks		    
	  for(l = 0; l < vLength; l++)
	    if(mask[l] != e->taxonMask[l])
	      break;
	  
	  if(l == vLength)	 	    	     
	    maskIdentical = TRUE;
	  
	  //if the bipartition vector and the taxon mask are identical
	  //just print out the frequency 
	  
	  if(maskIdentical && bitVectorIdentical)	
	    {
	      rax_free(bipCopy); 
	      return e;      		  	    
	    }
	  e = e->next;
	}
      while(e != (entry*)NULL); 
    }
        
  if((e == (entry*)NULL) && (readBit(mask, 0) == 0))
    {
      //consider the case where the first bit is unset and thus subset relations may be flipped.
      
      invertBits(bipCopy, length);
      unsetBit(bipCopy, 0);//only allowed since "readBit(taxonMasks->bitMasks[k], 0)==0"
      
      //mask the reference bipartition 
      for(l = 0; l < vLength; l++)
	bipCopy[l] = bipCopy[l] & mask[l];
      
      //again check that we have a canonic representation of the masked ref. bip.
      assert(!(bipCopy[0] & 1));
      
      //compute hash table index for the masked reference bipartition 
      position = oat_hash((unsigned char *)bipCopy, sizeof(unsigned int) * vLength);
      position = position % genesHash->tableSize;
      
      e = genesHash->table[position];
            
      if(e != (entry*)NULL)
	{
	  //if there is one or more elements here, we need to search the linked list 
	  //for a matching bipartition/mask pair 
	  //we don't necessarily need to find one !
	  
	  do 
	    {	     
	      boolean 
		maskIdentical = FALSE,
		bitVectorIdentical = FALSE;
	      
	      //make sure the bip. stored in the has table 
	      //is also in canonic form
	      assert(!(e->bitVector[0] & 1));
	      
	      //check for identity of bipartitions 				
	      for(l = 0; l < vLength; l++)
		if(bipCopy[l] != e->bitVector[l])
		  break;
	      
	      if(l == vLength)	 	    	     
		bitVectorIdentical = TRUE;	   	    
	      
	      //check for identity of taxon masks		    
	      for(l = 0; l < vLength; l++)
		if(mask[l] != e->taxonMask[l])
		  break;
	      
	      if(l == vLength)	 	    	     
		maskIdentical = TRUE;
	      	      
	      if(maskIdentical && bitVectorIdentical)
		{
		  rax_free(bipCopy); 
		  return e;
		}
			    			  
	      e = e->next;
	    }
	  while(e != (entry*)NULL); 
	}	
    }
    
  rax_free(bipCopy); 
  
  return (entry*)NULL;
}


static double computeTempFreq(int id, tree * tr, unsigned int * referenceBipartition, entry ** conflictingBips, unsigned int conflicts, hashtable * genesHash, 
			      taxonMaskArray * taxonMasks, unsigned int vLength)
{
  double
    refSupport = 0;

  unsigned int 
    k,
    i,
    *refBip;
  
  entry 
      *f,
      *e;
	           
  //Find sub partitions of the reference partition to count coveredNumber.   
  
  refBip = referenceBipartition;

  for(k = 0; k < (unsigned int) taxonMasks->entries; k++)
    {	
      //mask and rehash bipartition
      f = getSubPart(refBip, taxonMasks->bitMasks[k], genesHash, tr->mxtips);
      if(f != (entry*)NULL)
	{
	  assert(f->tempSupportFrom < id);
	  f->tempSupportFrom = id;
	  f->coveredNumber   = 1;
	}
    }
  
  //Do the same for conflicting bipartition.
  for(i = 0; i < conflicts; i++)
    {
      e = conflictingBips[i];
      
      //e cannot be covered as it would not be a conflicting bipartition otherwise.
      assert(!e->covered);
      
      refBip = e->bitVector;
      
      for(k = 0; k < (unsigned int) taxonMasks->entries; k++)
	{	
	  //check size of mask (This step is not strictly needed, as it is also included in the subset test below)
	  if(e->amountTips < bitCountLong(taxonMasks->bitMasks[k], tr->mxtips) )
	    continue;
	  
	  //check if mask is a subset of mask of e. If not, continue
	  if(!issubset(taxonMasks->bitMasks[k], e->taxonMask, vLength, 0))
	    continue;
	  
	  //mask and rehash bipartition
	  f = getSubPart(refBip, taxonMasks->bitMasks[k], genesHash, tr->mxtips);
	  if(f != (entry*)NULL)
	    {
	      if(f->tempSupportFrom != id)
		{
		  f->tempSupportFrom = id;
		  f->coveredNumber   = 1;
		} else {
		f->coveredNumber++;
	      }
	      
	    }
	} 
    }
  
	  
  // Now calculate temporary adjusted frequency (by distributing the frequency of sub bipartitions to those passed to this function.) 
  //NOTE: We would not strictly need to iterate through the same lists again if we store a map of which bipartition is covered by which bipartition.

  refBip = referenceBipartition;
  for(k = 0; k < (unsigned int) taxonMasks->entries; k++)
    {	
      //mask and rehash bipartition
      f = getSubPart(refBip, taxonMasks->bitMasks[k], genesHash, tr->mxtips);
      if(f != (entry*)NULL)
	{
	  assert(f->tempSupportFrom == id && f->coveredNumber > 0);
	  refSupport += f->supportFromTreeset[0] / f->coveredNumber;        
	}
    }
  
	  
  //	printf("supports: %f ",refSupport);
	  
  for(i = 0; i < conflicts; i++)
    {
      e = conflictingBips[i];
      e->tempSupport = 0;
      
      refBip = e->bitVector;
      
      for(k = 0; k < (unsigned int) taxonMasks->entries; k++)
	{	
	  //check size of mask (This step is not strictly needed, as it is also included in the subset test below)
	  if(e->amountTips < bitCountLong(taxonMasks->bitMasks[k], tr->mxtips) )
	    continue;
	  
	  //check if mask is a subset of mask of e. If not, continue
	  if(!issubset(taxonMasks->bitMasks[k], e->taxonMask, vLength, 0))
	    continue;
	  
	  //mask and rehash bipartition
	  f = getSubPart(refBip, taxonMasks->bitMasks[k], genesHash, tr->mxtips);
	  if(f != (entry*)NULL)
	    {
	      assert(f->tempSupportFrom == id && f->coveredNumber > 0);
	      
	      e->tempSupport += f->supportFromTreeset[0] / f->coveredNumber;
	    }
	} 
      
      //  printf("%f ", e->tempSupport);
    }
  
  //printf("\n\n");
  return refSupport;
}

static void calcBipartitions_IC_PartialGeneTrees(tree *tr, analdef *adef, char *bootStrapFileName, FILE *treeFile, int numberOfTrees, int numberOfTaxa)
{  
  branchInfo 
    *bInf;
  
  unsigned int
    k = 0,
    entryCount = 0,
    vLength = 0,   
    **bitVectors, 
    //array for easy access to the bipartitions of the reference tree
    **referenceBipartitions,      
    //Numerical links between bInfos and referenceBipartitions
    *referenceLinks = (unsigned int *)rax_malloc(sizeof(unsigned int ) * (size_t)(tr->mxtips - 3));
  
  int 
    bCount = 0,
    i;
     
  hashtable 
    //hash table containing the bips of the reference tree
    *h = initHashTable(tr->mxtips * 10),
    //hash table containing bips of the partial gene trees 
    *genesHash = initHashTable(tr->mxtips * 10);

  //data structure for parsing the partial gene trees
  tree   
    *inputTree = (tree *)rax_malloc(sizeof(tree));

  //data structure to store all taxon masks induced by the partial gene trees
  taxonMaskArray
    *taxonMasks;

  referenceBipartitions = (unsigned int **)rax_malloc(sizeof(unsigned int *) * (tr->mxtips - 3));

  //make sure that everything has been parsed correctly
  
  assert(numberOfTaxa == tr->mxtips);  
  
  //initialize data structure for extracting bipartitions from all trees 

  bitVectors = initBitVector(tr, &vLength);
    
  //initialize array for storing taxon masks 

  taxonMasks = initTaxonMaskArray(vLength);

  //extract the bipartitions of the reference tree and store them in hash table h

  bInf = (branchInfo*)rax_malloc(sizeof(branchInfo) * (tr->mxtips - 3)); 

  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, BIPARTITIONS_PARTIAL_TC, bInf,
			  &bCount, 1, FALSE, FALSE);

  //make sure this is a fully resolved bifurcating tree 

  assert(bCount == (tr->mxtips - 3));
  
  //loop over the hash table to set up the array pointing to the bips of the reference tree

  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
    {    
      if(h->table[k] != NULL)
	{
	  entry 
	    *e = h->table[k];
		  
	  /* we resolve collisions by chaining, hence the loop here */
		  
	  do
	    {
	      unsigned int 
		*bitVector = e->bitVector; 

	      //set the array to point to the bipartition 

	      referenceBipartitions[entryCount] = bitVector;
	      referenceLinks[entryCount]=e->bLink;

	      entryCount++;
	      	      			 		     		      
	      e = e->next;
	    }
	  while(e != NULL);
	}
    }
	  
  /* 
     make sure that we got all non-trivial bipartitions of the original tree 
     and stored them in referenceBipartitions 
   */
	  
  assert(entryCount == (unsigned int)(tr->mxtips - 3));      
  
  //get the number of trees in the partial gene tree set 

  //  treeFile = getNumberOfTrees(tr, bootStrapFileName, adef);

  //numberOfTrees = tr->numberOfTrees;

  //check that we have enough trees to do something meaningful, i.e., more than 1 tree 

  checkTreeNumber(numberOfTrees, bootStrapFileName);

  //initialize data structure required for parsing potentially multifurcating partial gene trees
  allocateMultifurcations(tr, inputTree);

  printBothOpen("Looping over partial gene trees\n");
  
  //loop over partial gene trees
  for(i = 0; i < numberOfTrees; i++)
    {     
      int 
	numberOfSplits;

      printBothOpen("\r %d / %d (%.2f %%)",i + 1,numberOfTrees, 100.0 * (i + 1.0) / ( numberOfTrees ));
      
      inputTree->start = (nodeptr)NULL;
      
      //parse partial gene tree

      numberOfSplits = readMultifurcatingTree(treeFile, inputTree, adef, FALSE);  //myabe set to FALSE for complete cleanup ???
      
      //if number of non-trivial bipartitions in this tree is greater than 0, analyze it further
      //otherwise it's not interesting 

      if(numberOfSplits > 0)
	{
	  int
	    taxa = 0,
	    splits = 0;
	  
	  unsigned int 
	    j,
	    //allocate tree mask 
	    *smallTreeMask = (unsigned int *)rax_calloc(vLength, sizeof(unsigned int));            		  	        

	  //  printf("Tree %d Splits %d Taxa %d\n", i, numberOfSplits, inputTree->ntips);

	  //first get the bit mask of taxa in this tree 
	  
	  setupMask(smallTreeMask, inputTree->start,       inputTree->mxtips);
	  setupMask(smallTreeMask, inputTree->start->back, inputTree->mxtips);
	  
	  //do a check that the number of set bits in this bit mask is actually 
	  //identical to the number of taxa in the tree 

	  for(j = 0; j < vLength; j++)
	    taxa += BIT_COUNT(smallTreeMask[j]);
	  assert(taxa == inputTree->ntips);

	  //add taxon mask to the dynamic array of bitmasks 

	  addTaxonMask(vLength, smallTreeMask, taxonMasks);

	  //now get the bipartitions of this tree and hash them into the hash table geneHash

	  bitVectorInitravPartial(inputTree, bitVectors, smallTreeMask, inputTree->start->back,  vLength, genesHash, &splits, inputTree->mxtips);

	  //make sure that we hashed all non-trivial bipartitions 
	  
	  assert(splits == numberOfSplits);
	  
	  //free the tree mask
	  rax_free(smallTreeMask);
	}           
    }     
  
  printBothOpen("\n");
  //we are done with initialization
  // 1. genesHash contains all non-trivial bipartitions of the partial gene trees, they have been stored as taxon mask and 
  //    bipartition vector pairs and their frequencies 
  // 2. the dynamic array taxonMasks contains all unique taxon masks found while analyzing the partial gene trees
  // 3. the array  referenceBipartitions conatins pointers to all bipartitions of the reference tree 

  {
     double 
       * referenceFreq = (double *)rax_malloc(sizeof(double) * (size_t)(tr->mxtips - 3)),
       * ic            = (double *)rax_malloc(sizeof(double) * (size_t)(tr->mxtips - 3)),
       * ica           = (double *)rax_malloc(sizeof(double) * (size_t)(tr->mxtips - 3)),
       * icNoLoss      = (double *)rax_malloc(sizeof(double) * (size_t)(tr->mxtips - 3)),
       * icaNoLoss     = (double *)rax_malloc(sizeof(double) * (size_t)(tr->mxtips - 3)),
       tc        = 0.0,
       tca       = 0.0,
       tcNoLoss  = 0.0,
       tcaNoLoss = 0.0;
     
     unsigned int
       * refBip,
       coveredBips      = 0,
       newlyCoveredBips = 0,
       skipped          = 0, //variable only needed for printout below
       entryCountPartial;
     
     entry
       *  e,
       *  f,
       ** partialList = (entry **)rax_malloc(sizeof(entry *) * (genesHash->entryCount));
     
     //Set up list of partial bipartitions (entries in hashtable) NOTE should be done when adding entries to the hashtable
     printBothOpen("Setting up list of partial bipartitions from hashtable.\n");
     
     for(i = 0, entryCountPartial = 0; (unsigned int)i < genesHash->tableSize; i++)	     
       {    
	 printBothOpen("\r %d / %d (%.2f %%)",i,genesHash->tableSize - 1, 100.0 * i / (genesHash->tableSize -1));
	 
	 if(genesHash->table[i] != NULL)
	   {
	     e = genesHash->table[i];
	     
	     /* we resolve collisions by chaining, hence the loop here */
	     
	     do
	       {
		 e->amountTips = bitCountLong(e->taxonMask, tr->mxtips);
		 partialList[entryCountPartial] = e;
		 entryCountPartial++;
		 
		 e = e->next;
	       }
	     while(e != NULL);
	   }
       }
     
     printBothOpen("\n");
      
     printBothOpen("Adjusting score for reference bipartitions\n"); 

     //Adjust the score for reference bipartitions
     for(i = 0; i < tr->mxtips - 3; i++) 
       referenceFreq[i] = 0.0;
     
     // loop through all reference partitions first
     for(i = 0; i < tr->mxtips - 3; i++)
       {	     
	 printBothOpen("\r %d / %d (%.2f %%) ",i  + 1,tr->mxtips - 3 , 100.0 * (i +1.0) / (tr->mxtips -3));
	 refBip = referenceBipartitions[i];
	 
	 for(k = 0; k < (unsigned int) taxonMasks->entries; k++)
	   {	
	     //mask and rehash bipartition
	     f = getSubPart(refBip, taxonMasks->bitMasks[k], genesHash, tr->mxtips);
	     if(f != (entry*)NULL)
	       {
		 if(!f->covered)
		   coveredBips++;
		 
		 f->covered=TRUE;
		 referenceFreq[i]+=adjustScore(tr->mxtips, f->amountTips, f->supportFromTreeset[0], tr->mxtips);
	       }
	   }
	 
       }
     
     printBothOpen("\n");
	
     printBothOpen("Sorting bipartitions.\n");

     qsort(partialList, entryCountPartial, sizeof(entry *), sortEntryByMask);
  
     printBothOpen("Adjusting support\n");
     
     for(i = 0; (unsigned int)i + coveredBips < entryCountPartial; i++)
       {
	 printBothOpen("\r %d / %d (%.2f %%)",i - skipped ,entryCountPartial - coveredBips - newlyCoveredBips -1, 100.0 *( i -skipped) / (entryCountPartial - coveredBips - newlyCoveredBips -1));
	 e = partialList[i];
	  
	 //If there is a strict super-bipartition, the adjusted support is not needed.
	 if(e->covered)
	   {
	     skipped++;
	     continue;
	   }
	 
	 refBip = e->bitVector;
	 
	 for(k = 0; k < (unsigned int) taxonMasks->entries; k++)
	   {	
	     //check size of mask (This step is not strictly needed, as it is also included in the subset test below)
	     if(e->amountTips < bitCountLong(taxonMasks->bitMasks[k], tr->mxtips) )
	       continue;
	     
	     //check if mask is a subset of mask of e. If not, continue
	     if(!issubset(taxonMasks->bitMasks[k], e->taxonMask, vLength, 0))
	       continue;
	     
	     //mask and rehash bipartition
	     f = getSubPart(refBip, taxonMasks->bitMasks[k], genesHash, tr->mxtips);
	     if(f != (entry*)NULL)
	       {
		 if(!f->covered && e->amountTips > f->amountTips)
		   {
		     //Mark f as having a super-bipartition
		     newlyCoveredBips++; 
		     f->covered=TRUE;
		   }
		 //adjust support
		 e->adjustedSupport += adjustScore(e->amountTips, f->amountTips, f->supportFromTreeset[0], tr->mxtips);
	       }
	   } 
	 
	 e->tempSupport = e->adjustedSupport;
       }
     
     printBothOpen("\n");
     printBothOpen("Sorting by adjusted Support\n");

     qsort(partialList, entryCountPartial - coveredBips, sizeof(entry *), sortEntryByAdjustedSupport);
  
#ifdef _DEBUG_PARTIAL_IC
     int il;
     //for(il=0 ; il < entryCountPartial - coveredBips +1; il++)
     for(il = 0 ; il < 30; il++)
       {
	 printf("%3d: %2d, %2d, %f",il, partialList[il]->amountTips, partialList[il]->supportFromTreeset[0], partialList[il]->adjustedSupport);
	 if(!partialList[il]->covered)
	   printf(" FALSE");
	 else
	   printf(" TRUE");
	 printf("\n");
       }
#endif
     
     coveredBips += newlyCoveredBips;

     unsigned int
       //number of conflicting bipartitions
       conflicts;
     //list of entries of conflicting bipartitions. 20 entries are sufficient, as we imply a 5% (1/20) threhold.
     //unsigned int * conflictingBips = (unsigned int *)rax_malloc(sizeof(unsigned int) * 20);
     entry
       **  conflictingBips = (entry **)rax_malloc(sizeof(entry *) * 20);
     
     
     // loop through all reference partitions to obtain ic/ica scores
     for(i = 0; i < tr->mxtips - 3; i++)
       {
	 // for each reference bip, find conflicting bipartitions and store index in "conflictingBips".	      
	 conflicts = 0;
	 refBip = referenceBipartitions[i];
	 	 
	 conflicts = findConflictingBipartitions(conflictingBips,refBip , partialList, (entryCountPartial - coveredBips), tr->mxtips);
	 assert(conflicts < 21);                             
      
	 //calc ic,ica
	 ic[i]        = 1;
	 ica[i]       = 1;
	 icNoLoss[i]  = 1;
	 icaNoLoss[i] = 1;
	 
	 double 
	   tempReferenceFreq;
      
	 if(conflicts > 0)
	   {	     
	     ic[i]  = computeIC_ValueDouble(referenceFreq[i], conflictingBips[0]->adjustedSupport , FALSE);      
	     ica[i] = computeICA_ValueDouble(referenceFreq[i], conflictingBips, conflicts, FALSE);
	     
	     //assert they have the same sign
	     assert(ic[i]*ica[i] >= 0);
	     
	     //computeTempFreq will set e->tempSupport to the correct values for all e in conflicting bips & reference bip.
	     //Before this, tempSupport was equal to adjustedSupport.
	     tempReferenceFreq = computeTempFreq(i, tr, refBip, conflictingBips, conflicts, genesHash, taxonMasks, vLength);
	     
	     qsort(conflictingBips, conflicts, sizeof(entry *), sortEntryByTempSupport);
	     
	     icNoLoss[i]  = computeIC_ValueDouble(tempReferenceFreq, conflictingBips[0]->tempSupport , FALSE);
	     icaNoLoss[i] = computeICAFromTemp_ValueDouble(tempReferenceFreq, conflictingBips, conflicts, FALSE, tr->numberOfTrees);	           
	   }
      
#ifdef _DEBUG_PARTIAL_IC      
	 if(conflicts > 0){
	   unsigned int j;
	   printf("----------------------------------------------------------------\n");
	   printf("ref (%d) adjusted support : %f, ref temp Support %f \n",i,referenceFreq[i], tempReferenceFreq);
	   printf("conflicts %u \t:\t",conflicts);
	   for(j=0; j < conflicts;j++)
	     {printf(" %4u ", conflictingBips[j]->amountTips);
	       //printf("\t%u ", conflictingBips[j]);
	     }
	   printf("\n");
	   
	   printf("support f tree \t:\t");
	   for(j=0; j < conflicts;j++)
	     {printf(" %4u ", conflictingBips[j]->supportFromTreeset[0]);
	       //printf("\t%u ", conflictingBips[j]);
	     }
	   printf("\n");
	   
	   printf("adjusted     \t:\t");
	   for(j=0; j < conflicts;j++)
	     {printf(" %4.2f ", conflictingBips[j]->adjustedSupport);
	       //printf("\t%u ", conflictingBips[j]);
	     }
	   printf("\n");
	   
	   printf("tempadjusted \t:\t");
	   for(j=0; j < conflicts;j++)
	     {printf(" %4.2f ", conflictingBips[j]->tempSupport);
	       //printf("\t%u ", conflictingBips[j]);
	     }
	   printf("\n");     
	 }
#endif    
	 //Add to the tc/tcAll
	 tc  += ic[i];		
	 tca += ica[i];	
	 
	 tcNoLoss  += icNoLoss[i];
	 tcaNoLoss += icaNoLoss[i];	 
       }
    
     rax_free(conflictingBips);
       
     //set ic/icAll for actual branches (where sizes of sub bipartitions matter) 
     for(i = 0; i < (tr->mxtips - 3); i++)
       {
	 bInf[referenceLinks[i]].ic = ic[i];
	 bInf[referenceLinks[i]].icAll = ica[i];
       }
     //print to tree file   
     printBipartitionResult(tr, adef, TRUE, TRUE, icFileNameBranchLabelsStochastic); 
     
     assert(tc <= (double)(tr->mxtips - 3));
     
     printBothOpen("Tree certainty under stochastic bipartition adjustment for this tree: %f\n", tc);
     printBothOpen("Relative tree certainty under stochastic bipartition adjustment for this tree: %f\n\n", tc / (double)(tr->mxtips - 3));
     
     printBothOpen("Tree certainty including all conflicting bipartitions (TCA) under stochastic bipartition adjustment for this tree: %f\n", tca);
     printBothOpen("Relative tree certainty including all conflicting bipartitions (TCA) under stochastic bipartition adjustment for this tree: %f\n\n", tca / (double)(tr->mxtips - 3));
      
     //set ic/icAll for actual branches (completely distributed to conflicting bipartitions) 
     for(i = 0; i < (tr->mxtips - 3); i++)
       {
	 bInf[referenceLinks[i]].ic = icNoLoss[i];
	 bInf[referenceLinks[i]].icAll = icaNoLoss[i];
       }
     //print to tree file   
     printBipartitionResult(tr, adef, TRUE, TRUE, icFileNameBranchLabelsUniform); 
     
     assert(tc <= (double)(tr->mxtips - 3));
     
     printBothOpen("Tree certainty under uniform bipartition adjustment for this tree: %f\n", tcNoLoss);
     printBothOpen("Relative tree certainty under uniform bipartition adjustment for this tree: %f\n\n", tcNoLoss / (double)(tr->mxtips - 3));
     
     printBothOpen("Tree certainty including all conflicting bipartitions (TCA) under uniform bipartition adjustment for this tree: %f\n", tcaNoLoss);
     printBothOpen("Relative tree certainty including all conflicting bipartitions (TCA) under uniform bipartition adjustment for this tree: %f\n\n", tcaNoLoss / (double)(tr->mxtips - 3));
     
     
     
#ifdef _DEBUG_PARTIAL_IC
     printf("probability ic\t");
     
     for(i = 0 ; i < (unsigned int)(tr->mxtips - 3); i++)    
       printf("& %4d ", (int)(ic[i]*100)); 
     
     printf("\n");
     
     printf("probability ica\t");
     
     for(i = 0; i < (unsigned int)(tr->mxtips - 3); i++)    
       printf("& %4d ", (int)(ica[i]*100)); 
     
     printf("\n");
     printf("\n");
     
     printf("Lossless ic  \t");
     
     for(i = 0; i < (unsigned int)(tr->mxtips - 3); i++)    
       printf("& %4d ", (int)(icNoLoss[i]*100)); 
     
     printf("\n");
     printf("Lossless ica \t");
     
     for(i = 0; i < (unsigned int)(tr->mxtips - 3); i++)     
       printf("& %4d ", (int)(icaNoLoss[i]*100));      
     
     
     printf("\n");
     printf("\n");
#endif    
     
    /* consensus tree building from modified bip freqs, do not erase, to be continued 
     *	{
     *	  entry 
     **consensusBips = (entry **)rax_calloc(tr->mxtips - 3, sizeof(entry *));
     
     int
     consensusBipsLen = 0;  
     
     //TODO how do we set the threshold  
     
     //threshold definition correct? I use tree counts, Kassian uses real frequencies or doesn't he?
     
     //why is mr_tresh an integer and not a double????
     
     tr->mr_thresh = ((double)numberOfTrees / 2.0);          
   
     printf("trees %d threshold %f\n", numberOfTrees, tr->mr_thresh);
     
     printf(" entries in potential bips: %d\n", fullBips->entries);
     
     for(i = 0; i < fullBips->entries; i++)
     {
     int 
     cnt = fullBips->potentialFrequencies[i]; //is potential frequencies the correct variable?
     
     if((tr->consensusType == MR_CONSENSUS     && cnt > (unsigned int)tr->mr_thresh) || 
     (tr->consensusType == STRICT_CONSENSUS && cnt == numberOfTrees) ||
     (tr->consensusType ==  USER_DEFINED    && cnt > (numberOfTrees * tr->consensusUserThreshold) / 100))
     {
     entry 
     e = initEntry();
     
     //has vector length been set?
     
     e->bitVector = fullBips->bitMasks[i];
     //make sure that we have a canonic representation of the bit vector!
     assert(!(e->bitVector[0] & 1));
     //TODO check bitmask correct -> correct part of bits set to 1!
     
     
     consensusBips[consensusBipsLen] = e;
     consensusBipsLen++;
     }
     }
     
     
     
     //
     }
     
     //TODO pass this mask to consensus tree building algorithm 
     // type is defined as referenceMaskArray 
     
     */
    
  
     rax_free(icNoLoss);
     rax_free(icaNoLoss);
     rax_free(referenceFreq);
     rax_free(ic);
     rax_free(ica);            
     rax_free(partialList);
  }
  
  fclose(treeFile);
  
  //free data structures 
  
  freeMultifurcations(inputTree);
  rax_free(inputTree); 
  
  rax_free(referenceBipartitions);
  
  freeTaxonMask(taxonMasks);
  rax_free(taxonMasks);
  
  freeHashTable(genesHash);
  rax_free(genesHash); 
  
  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h); 
  
  rax_free(referenceLinks);
  rax_free(bInf);
}






