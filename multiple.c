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


#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"


#ifdef _WAYNE_MPI
#include <mpi.h>
extern int processID;
extern int processes;
#endif

extern int  optimizeRatesInvocations;
extern int  optimizeRateCategoryInvocations;
extern int  optimizeAlphaInvocations;
extern int  optimizeInvarInvocations;
extern int  checkPointCounter;
extern int  Thorough;
extern int  partCount;
extern char tree_file[1024];
extern const unsigned int mask32[32];
extern double masterTime;

extern FILE   *INFILE, *permutationFile, *logFile, *infoFile;

extern char seq_file[1024];
extern char permFileName[1024], resultFileName[1024], 
  logFileName[1024], checkpointFileName[1024], infoFileName[1024], run_id[128], workdir[1024], bootStrapFile[1024], bootstrapFileName[1024], 
  bipartitionsFileName[1024],bipartitionsFileNameBranchLabels[1024]; 




void catToGamma(tree *tr, analdef *adef)
{
  assert(tr->rateHetModel == CAT);  
  
  if(adef->useInvariant)
    tr->rateHetModel = GAMMA_I;
  else
    tr->rateHetModel = GAMMA;

  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_CAT_TO_GAMMA, tr); 
#endif
}

void gammaToCat(tree *tr)
{

  assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
   
  tr->rateHetModel = CAT;

 
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_GAMMA_TO_CAT, tr); 
#endif  
}


static void singleBootstrap(tree *tr, int i, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  tr->treeID = i;
  tr->checkPointCounter = 0;
     
  computeNextReplicate(tr, &adef->boot, (int*)NULL, (int*)NULL, FALSE, FALSE);
  
  initModel(tr, rdta, cdta, adef);                                   
         
  getStartingTree(tr, adef);

  if(i == 0)
    testGapped(tr);

  computeBIGRAPID(tr, adef, TRUE);

  if(adef->bootstrapBranchLengths)
    {
      switch(tr->rateHetModel)
	{
	case GAMMA:	  
	case GAMMA_I:      
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);		      	    	    
	  break;
	case CAT:    	      	 		  
	  tr->likelihood = unlikely;	       
	  catToGamma(tr, adef);	  
	  initModel(tr, rdta, cdta, adef);
	  if(i == 0)
	    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);	 
	  else
	    modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);
	  gammaToCat(tr); 		  	      	        	
	  break;
	default:
	  assert(0);
	}
    }  


   printBootstrapResult(tr, adef, TRUE);    
}





/***************************** EXPERIMENTAL FUNCTIONS ********************************************************/




static int compareTopolRell(const void *p1, const void *p2)
{
  topolRELL **rc1 = (topolRELL **)p1;
  topolRELL **rc2 = (topolRELL **)p2;
 
  double i = (*rc1)->likelihood;
  double j = (*rc2)->likelihood;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


void fixModelIndices(tree *tr, int endsite, boolean fixRates)
{
  int model, i;

  assert(tr->NumberOfModels > 0);   

  tr->partitionData[0].lower = 0;
     
  model = tr->model[0];
  i = 1;

  while(i < endsite)
    {
      if(tr->model[i] != model)
	{	      
	  tr->partitionData[model].upper = i;
	  tr->partitionData[model + 1].lower = i;
	  model = tr->model[i];
	}
      i++;
    }       
  
  tr->partitionData[tr->NumberOfModels - 1].upper = endsite;
  
  for(model = 0; model < tr->NumberOfModels; model++)    
    tr->partitionData[model].width = tr->partitionData[model].upper -  tr->partitionData[model].lower;
 
  
#ifndef _USE_PTHREADS
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	j,
	lower =  tr->partitionData[model].lower;
      
      /* SOS what about sumBuffer? */
      /* tr->partitionData[model].sumBuffer    = &tr->sumBuffer[offset]; */
      tr->partitionData[model].perSiteLL    = &tr->perSiteLL[lower];     
      tr->partitionData[model].wgt          = &tr->cdta->aliaswgt[lower];
      

      tr->partitionData[model].invariant    = &tr->invariant[lower];
      tr->partitionData[model].rateCategory = &tr->cdta->rateCategory[lower];
      

      for(j = 1; j <= tr->mxtips; j++)
	tr->partitionData[model].yVector[j] = &(tr->yVector[j][tr->partitionData[model].lower]);

      
      {
	int 
	  width =  tr->partitionData[model].width,
	  undetermined = getUndetermined(tr->partitionData[model].dataType),
	  j;		         		     
	
	tr->partitionData[model].gapVectorLength = ((int)width / 32) + 1;

	memset(tr->partitionData[model].gapVector, 0, tr->partitionData[model].initialGapVectorSize);

	for(j = 1; j <= tr->mxtips; j++)
	  for(i = 0; i < width; i++)
	    if(tr->partitionData[model].yVector[j][i] == undetermined)
	      tr->partitionData[model].gapVector[tr->partitionData[model].gapVectorLength * j + i / 32] |= mask32[i % 32];
      }

    }
#else
  masterBarrier(THREAD_FIX_MODEL_INDICES, tr);
#endif

  
  if(fixRates)
    updatePerSiteRates(tr, TRUE);
}

void reductionCleanup(tree *tr, int *originalRateCategories, int *originalInvariant)
{
  tr->cdta->endsite = tr->originalCrunchedLength;

  memcpy(tr->cdta->aliaswgt, tr->originalWeights, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->model, tr->originalModel, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->dataVector, tr->originalDataVector,  sizeof(int) * tr->cdta->endsite);

  memcpy(tr->cdta->rateCategory, originalRateCategories, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->invariant,          originalInvariant,      sizeof(int) * tr->cdta->endsite);                      
  
  updatePerSiteRates(tr, TRUE);

  memcpy(tr->rdta->y0, tr->rdta->yBUF, ((size_t)tr->rdta->numsp) * ((size_t)tr->cdta->endsite) * sizeof(char));  
      
  tr->cdta->endsite = tr->originalCrunchedLength;
  fixModelIndices(tr, tr->originalCrunchedLength, TRUE);      
}






void computeNextReplicate(tree *tr, long *randomSeed, int *originalRateCategories, int *originalInvariant, boolean isRapid, boolean fixRates)
{ 
  int pos, nonzero, j, model, w;   
  int *weightBuffer, endsite;                
  int *weights, i, l;  


  for(j = 0; j < tr->originalCrunchedLength; j++)
    tr->cdta->aliaswgt[j] = 0;

      	  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      nonzero = 0;        	 
      
      for (j = 0; j < tr->originalCrunchedLength; j++)  
	{
	  if(tr->originalModel[j] == model)
	    nonzero += tr->originalWeights[j];
	}				          
      
      weightBuffer = (int *)calloc(nonzero, sizeof(int));	 
     
      for (j = 0; j < nonzero; j++)
	weightBuffer[(int) (nonzero*randum(randomSeed))]++;                  
      
      pos = 0;	      
      
      for(j = 0; j < tr->originalCrunchedLength; j++) 
	{
	  if(model == tr->originalModel[j])
	    {
	      for(w = 0; w < tr->originalWeights[j]; w++)	  	  	 
		{
		  tr->cdta->aliaswgt[j] += weightBuffer[pos];
		  pos++;		      
		}				   
	    }
	}  

      free(weightBuffer);	  
      
    }       

  endsite = 0;
  
  for (j = 0; j < tr->originalCrunchedLength; j++)         
    if(tr->cdta->aliaswgt[j] > 0)
      endsite++;    
  
  weights = tr->cdta->aliaswgt;

  for(i = 0; i < tr->rdta->numsp; i++)
    {     
      unsigned char *yPos    = &(tr->rdta->y0[((size_t)tr->originalCrunchedLength) * ((size_t)i)]);
      unsigned char *origSeq = &(tr->rdta->yBUF[((size_t)tr->originalCrunchedLength) * ((size_t)i)]);
      int l, j;
      
      for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)      
	if(tr->cdta->aliaswgt[j] > 0)	  	    
	  yPos[l++] = origSeq[j];	                   
    }

  for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)
    {     
      if(weights[j])	
	{
	  tr->cdta->aliaswgt[l]     = tr->cdta->aliaswgt[j];
	  tr->dataVector[l]         = tr->originalDataVector[j];
	  tr->model[l]              = tr->originalModel[j];

	  if(isRapid)
	    {	      
	      tr->cdta->rateCategory[l] = originalRateCategories[j];
	      tr->invariant[l]          = originalInvariant[j];
	    }
	  l++;
	}
    }

  tr->cdta->endsite = endsite;
  fixModelIndices(tr, endsite, fixRates);
}


  




static pInfo *allocParams(tree *tr)
{
  int i;
  pInfo *partBuffer = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      const partitionLengths *pl = getPartitionLengths(&(tr->partitionData[i]));

      partBuffer[i].EIGN = (double*)malloc(pl->eignLength * sizeof(double));
      partBuffer[i].EV   = (double*)malloc(pl->evLength * sizeof(double));
      partBuffer[i].EI   = (double*)malloc(pl->eiLength * sizeof(double));	  
      partBuffer[i].substRates = (double *)malloc(pl->substRatesLength * sizeof(double));	  
      partBuffer[i].frequencies =  (double*)malloc(pl->frequenciesLength * sizeof(double));	  
      partBuffer[i].tipVector   = (double *)malloc(pl->tipVectorLength * sizeof(double));
      
     
    }

  return partBuffer;      
}

static void freeParams(int numberOfModels, pInfo *partBuffer)
{
  int i;

  for(i = 0; i < numberOfModels; i++)
    {
      free(partBuffer[i].EIGN); 
      free(partBuffer[i].EV);   
      free(partBuffer[i].EI);   
      free(partBuffer[i].substRates);
      free(partBuffer[i].frequencies); 
      free(partBuffer[i].tipVector);  
    }
      
}

static void copyParams(int numberOfModels, pInfo *dst, pInfo *src, tree *tr)
{
  int i;

  assert(src != dst);

  for(i = 0; i < numberOfModels; i++)
    {
      const partitionLengths *pl = getPartitionLengths(&src[i]);
      
      dst[i].dataType = src[i].dataType;

       memcpy(dst[i].EIGN,        src[i].EIGN,        pl->eignLength * sizeof(double));
       memcpy(dst[i].EV,          src[i].EV,          pl->evLength * sizeof(double));
       memcpy(dst[i].EI,          src[i].EI,          pl->eiLength * sizeof(double));	  
       memcpy(dst[i].substRates,  src[i].substRates,  pl->substRatesLength * sizeof(double));	  
       memcpy(dst[i].frequencies, src[i].frequencies, pl->frequenciesLength * sizeof(double));	  
       memcpy(dst[i].tipVector,   src[i].tipVector,   pl->tipVectorLength * sizeof(double));
       
       
    }
  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_PARAMS, tr);
#endif    
}


#ifdef _WAYNE_MPI

static void copyFiles(FILE *source, FILE *destination)
{
  size_t 
    c;
  
  char 
    buf[1024];
  
  assert(processID == 0);

  while((c = fread(buf, sizeof(char), 1024, source)) > 0)
    fwrite(buf, sizeof(char), c, destination);     	
}

static void concatenateBSFiles(int processes)
{
  if(processID == 0)
    {
      int i;
      
      FILE 
	*destination = myfopen(bootstrapFileName, "w"),
	*source = (FILE*)NULL;
      
      char 
	sourceName[1024];          
      
      strcpy(sourceName, bootstrapFileName);
      strcat(sourceName, ".PID.");
      
      for(i = 0; i < processes; i++)
	{
	  char 
	    buf[64],
	    temporary[1024];
	  
	  sprintf(buf, "%d", i);
	  strcpy(temporary, sourceName);
	  strcat(temporary, buf);
	  
	  source = myfopen(temporary, "r");
	  
	  copyFiles(source, destination);
	  
	  fclose(source);
	}
      
      fclose(destination);
    }
}

static void removeBSFiles(int processes)
{
  if(processID == 0)
    {
      int 
	i;
      
      char 
	sourceName[1024];            
      
      strcpy(sourceName, bootstrapFileName);
      strcat(sourceName, ".PID.");
      
      for(i = 0; i < processes; i++)
	{
	  char 
	    buf[64],
	    temporary[1024];
	  
	  sprintf(buf, "%d", i);
	  strcpy(temporary, sourceName);
	  strcat(temporary, buf);
	  
	  remove(temporary);
	}
    }
}

#endif


void doAllInOne(tree *tr, analdef *adef)
{
  int i, n, sites, bestIndex, bootstrapsPerformed;

#ifdef _WAYNE_MPI
  int 
    bootStopTests = 1,
    j,
    bootStrapsPerProcess = 0;
#endif 

  double loopTime; 
  int      *originalRateCategories;
  int      *originalInvariant;
#ifdef _WAYNE_MPI
  int      slowSearches, fastEvery;
#else
  int      slowSearches, fastEvery = 5;
#endif
  int treeVectorLength = -1;
  topolRELL_LIST *rl;  
  double bestLH, mlTime, overallTime;  
  long radiusSeed = adef->rapidBoot;
  FILE *f;
  char bestTreeFileName[1024];  
  hashtable *h = (hashtable*)NULL;
  unsigned int **bitVectors = (unsigned int**)NULL;
  boolean bootStopIt = FALSE;
  double pearsonAverage = 0.0;
  pInfo *catParams         = allocParams(tr);
  pInfo *gammaParams = allocParams(tr);
  unsigned int vLength;

  n = adef->multipleRuns; 

#ifdef _WAYNE_MPI
  if(n % processes != 0)
    n = processes * ((n / processes) + 1);
#endif

  if(adef->bootStopping)
    {    
      h = initHashTable(tr->mxtips * 100);

      treeVectorLength = adef->multipleRuns;
      
      bitVectors = initBitVector(tr, &vLength);          
    }

  rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
  initTL(rl, tr, n);
     
  originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));      
  originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));

  sites = tr->cdta->endsite;             

  initModel(tr, tr->rdta, tr->cdta, adef);

  if(adef->grouping)
    printBothOpen("\n\nThe topologies of all Bootstrap and ML trees will adhere to the constraint tree specified in %s\n", tree_file);
  if(adef->constraint)
    printBothOpen("\n\nThe topologies of all Bootstrap and ML trees will adhere to the bifurcating backbone constraint tree specified in %s\n", tree_file);
 

#ifdef _WAYNE_MPI
  long parsimonySeed0 = adef->parsimonySeed;
  long replicateSeed0 = adef->rapidBoot;
  n = n / processes;
#endif
 
  for(i = 0; i < n && !bootStopIt; i++)
    {  
#ifdef _WAYNE_MPI
      j = i + n * processID;
      tr->treeID = j;
#else              
      tr->treeID = i;
#endif

      tr->checkPointCounter = 0;
        
      loopTime = gettime();  

#ifdef _WAYNE_MPI
      if(i == 0)
        {
          if(parsimonySeed0 != 0)
            adef->parsimonySeed = parsimonySeed0 + 10000 * processID;
          adef->rapidBoot = replicateSeed0 + 10000 * processID;
          radiusSeed = adef->rapidBoot;
        }
#endif          
     
      if(i % 10 == 0)
	{
	  if(i > 0)	    	    
	    reductionCleanup(tr, originalRateCategories, originalInvariant);	    	  

	  if(adef->grouping || adef->constraint)
	    {
	      FILE *f = myfopen(tree_file, "rb");	

	      assert(adef->restart);
	      partCount = 0;
	      if (! treeReadLenMULT(f, tr, adef))
		exit(-1);
	     
	      fclose(f);
	    }
	  else
	    makeParsimonyTree(tr, adef);
	  
	  tr->likelihood = unlikely;
	  if(i == 0)
	    {
	      double t;
	          
	      onlyInitrav(tr, tr->start);
	      treeEvaluate(tr, 1);	     	
	     	      
	      t = gettime();    	      

	      modOpt(tr, adef, FALSE, 5.0, TRUE);	    
#ifdef _WAYNE_MPI
	      printBothOpen("\nTime for BS model parameter optimization on Process %d: %f seconds\n", processID, gettime() - t);	     
#else
	      printBothOpen("\nTime for BS model parameter optimization %f\n", gettime() - t);
#endif
	      
	      memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
	      memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

	      if(adef->bootstrapBranchLengths)
		{
		  if(tr->rateHetModel == CAT)
		    {
		      copyParams(tr->NumberOfModels, catParams, tr->partitionData, tr);		      
		      assert(tr->cdta->endsite == tr->originalCrunchedLength);		 
		      catToGamma(tr, adef);		      
		      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);
		      copyParams(tr->NumberOfModels, gammaParams, tr->partitionData, tr);		      
		      gammaToCat(tr);
		      copyParams(tr->NumberOfModels, tr->partitionData, catParams, tr);		      
		    }
		  else
		    {		  
		      assert(tr->cdta->endsite == tr->originalCrunchedLength);		 		     		     		      		     
		    }
		}
	    }	  	  
	}

      computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE, TRUE); 
      resetBranches(tr);

     

      evaluateGenericInitrav(tr, tr->start);
    
      treeEvaluate(tr, 1);    	             
     
      computeBOOTRAPID(tr, adef, &radiusSeed);  
#ifdef _WAYNE_MPI
      saveTL(rl, tr, j);
#else                      	  
      saveTL(rl, tr, i);
#endif

      if(adef->bootstrapBranchLengths)
	{
	  double lh = tr->likelihood;
	  int    endsite;
	 
	  if(tr->rateHetModel == CAT)
	    {
	      copyParams(tr->NumberOfModels, tr->partitionData, gammaParams, tr);	      
	      /*endsite = tr->cdta->endsite;
		tr->cdta->endsite = tr->originalCrunchedLength;*/
	      catToGamma(tr, adef);
	      /*tr->cdta->endsite = endsite;*/
	      
	      resetBranches(tr);
	      onlyInitrav(tr, tr->start);
	      treeEvaluate(tr, 2.0);
	  
	      /*endsite = tr->cdta->endsite;
		tr->cdta->endsite = tr->originalCrunchedLength;*/
	      gammaToCat(tr);
	      /*tr->cdta->endsite = endsite;	 	    */
	
	      copyParams(tr->NumberOfModels, tr->partitionData, catParams, tr);	      
	      tr->likelihood = lh;
	    }
	  else
	    {	     
	      treeEvaluate(tr, 2.0);
	      tr->likelihood = lh;
	    }
	}
      
      printBootstrapResult(tr, adef, TRUE); 

      loopTime = gettime() - loopTime; 
      writeInfoFile(adef, tr, loopTime); 
     
      if(adef->bootStopping)
#ifdef _WAYNE_MPI
	{
	  int 
	    nn = (i + 1) * processes;

	  if((nn > START_BSTOP_TEST) && 
	     (i * processes < FC_SPACING * bootStopTests) &&
	     ((i + 1) * processes >= FC_SPACING * bootStopTests)
	     )	     
	    {
	      MPI_Barrier(MPI_COMM_WORLD);
	                    
	      concatenateBSFiles(processes);                
	      
              MPI_Barrier(MPI_COMM_WORLD);	      
	      
	      bootStopIt = computeBootStopMPI(tr, bootstrapFileName, adef, &pearsonAverage);
	      bootStopTests++;
	    }
	}	
#else	
	bootStopIt = bootStop(tr, h, i, &pearsonAverage, bitVectors, treeVectorLength, vLength);
#endif


    }  
 
#ifdef _WAYNE_MPI      
  MPI_Barrier(MPI_COMM_WORLD);
  
  bootstrapsPerformed = i * processes; 
  bootStrapsPerProcess = i;   
      
  concatenateBSFiles(processes);
  removeBSFiles(processes);  
  
  MPI_Barrier(MPI_COMM_WORLD); 
#else
  bootstrapsPerformed = i;
#endif

  freeParams(tr->NumberOfModels, catParams);
  free(catParams);

  freeParams(tr->NumberOfModels, gammaParams);
  free(gammaParams);

  if(adef->bootStopping)
    {
      freeBitVectors(bitVectors, 2 * tr->mxtips);
      free(bitVectors);
      freeHashTable(h);
      free(h);      
    }

 
  {      
    double t;

    printBothOpenMPI("\n\n");
    
    if(adef->bootStopping)
      {
	if(bootStopIt)
	  {
	    switch(tr->bootStopCriterion)
	      {
	      case FREQUENCY_STOP:
		printBothOpenMPI("Stopped Rapid BS search after %d replicates with FC Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);	      
		break;
	      case MR_STOP:
		printBothOpenMPI("Stopped Rapid BS search after %d replicates with MR-based Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);	     
		break;
	      case MRE_STOP:
		printBothOpenMPI("Stopped Rapid BS search after %d replicates with MRE-based Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);	     
		break;
	      case MRE_IGN_STOP:
		printBothOpenMPI("Stopped Rapid BS search after %d replicates with MRE_IGN-based Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);	     
		break;
	      default:
		assert(0);
	      }
	  }
	else
	  { 
	    switch(tr->bootStopCriterion)	     
	      {
	      case FREQUENCY_STOP:
		printBothOpenMPI("Rapid BS search did not converge after %d replicates with FC Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);
		break;
	      case MR_STOP:
		printBothOpenMPI("Rapid BS search did not converge after %d replicates with MR-based Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
		break;
	      case MRE_STOP:
		printBothOpenMPI("Rapid BS search did not converge after %d replicates with MRE-based Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
		break;
	      case MRE_IGN_STOP:
		printBothOpenMPI("Rapid BS search did not converge after %d replicates with MR_IGN-based Bootstopping criterion\n", bootstrapsPerformed);
		printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
		break;
	      default:
		assert(0);
	      }
	  }
      }
    

    t = gettime() - masterTime;

    printBothOpenMPI("Overall Time for %d Rapid Bootstraps %f seconds\n", bootstrapsPerformed, t);     
    printBothOpenMPI("Average Time per Rapid Bootstrap %f seconds\n", (double)(t/((double)bootstrapsPerformed)));  
        
    if(!adef->allInOne)     
      {
	printBothOpenMPI("All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);

#ifdef _WAYNE_MPI      	 
	MPI_Finalize();
#endif
	exit(0);
      }
  }
 
  
  /* ML-search */ 

  mlTime = gettime();
  double t = mlTime;
  
  printBothOpenMPI("\nStarting ML Search ...\n\n"); 

  /***CLEAN UP reduction stuff */  

  reductionCleanup(tr, originalRateCategories, originalInvariant);  

  /****/     	   
  
#ifdef _WAYNE_MPI 
  restoreTL(rl, tr, n * processID); 
#else
  restoreTL(rl, tr, 0);
#endif

  resetBranches(tr);

  

  evaluateGenericInitrav(tr, tr->start);   

  

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);  

#ifdef _WAYNE_MPI
  
  if(bootstrapsPerformed <= 100)
    fastEvery = 5;
  else
    fastEvery = bootstrapsPerformed / 20;

  for(i = 0; i < bootstrapsPerformed; i++)
    rl->t[i]->likelihood = unlikely;

  for(i = 0; i < bootStrapsPerProcess; i++)
    {            
      j = i + n * processID;
    
      if(i % fastEvery == 0)
	{	 
	  restoreTL(rl, tr, j); 	 	    	   	
	  
	  resetBranches(tr);	 

	  evaluateGenericInitrav(tr, tr->start);
	  	  
	  treeEvaluate(tr, 1); 		 
	  	  
	  optimizeRAPID(tr, adef);	  			         	  
	  
	  saveTL(rl, tr, j);  
	}    
    }     
#else
  for(i = 0; i < bootstrapsPerformed; i++)
    {            
      rl->t[i]->likelihood = unlikely;
    
      if(i % fastEvery == 0)
	{
	 
	  
	  restoreTL(rl, tr, i); 	 	    	   	
	  
	  resetBranches(tr);	 

	  evaluateGenericInitrav(tr, tr->start);
	  	  
	  treeEvaluate(tr, 1); 		 
	  	  
	  optimizeRAPID(tr, adef);	  			         	  
	  
	 

	  saveTL(rl, tr, i); 	 
	}    
    }     
#endif
 
  printBothOpenMPI("Fast ML optimization finished\n\n"); 
  t = gettime() - t;
  
#ifdef _WAYNE_MPI
  printBothOpen("Fast ML search on Process %d: Time %f seconds\n\n", processID, t);
  j = n * processID;

  qsort(&(rl->t[j]), n, sizeof(topolRELL*), compareTopolRell);

  restoreTL(rl, tr, j);
#else
  printBothOpen("Fast ML search Time: %f seconds\n\n", t);
  qsort(&(rl->t[0]), bootstrapsPerformed, sizeof(topolRELL*), compareTopolRell);
       
  restoreTL(rl, tr, 0);
#endif
  t = gettime();
  
  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);

  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);     
  
  slowSearches = bootstrapsPerformed / 5;
  if(bootstrapsPerformed % 5 != 0)
    slowSearches++;

  slowSearches  = MIN(slowSearches, 10); 

#ifdef _WAYNE_MPI
   if(processes > 1)
    {
      if(slowSearches % processes == 0)
        slowSearches = slowSearches / processes;
      else
        slowSearches = (slowSearches / processes) + 1;
    }
   
   for(i = 0; i < slowSearches; i++)
    {           
      j = i + n * processID;
      restoreTL(rl, tr, j);     
      rl->t[j]->likelihood = unlikely;  
      
      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 1.0);   
      
      thoroughOptimization(tr, adef, rl, j); 
   }   
#else
  for(i = 0; i < slowSearches; i++)
    {           
      restoreTL(rl, tr, i);     
      rl->t[i]->likelihood = unlikely;  
      
      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 1.0);   
      
      thoroughOptimization(tr, adef, rl, i); 	 

   }
#endif
  
  

  /*************************************************************************************************************/  
  
  if(tr->rateHetModel == CAT) 
    {      
      catToGamma(tr, adef);    
      modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE); 
    }

  bestIndex = -1;
  bestLH = unlikely;
    
#ifdef _WAYNE_MPI
  for(i = 0; i < slowSearches; i++)
    { 
      j = i + n * processID;
      restoreTL(rl, tr, j);
      resetBranches(tr);

      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 2);
      
      printBothOpen("Slow ML Search %d Likelihood: %f\n", j, tr->likelihood);
      
      if(tr->likelihood > bestLH)
	{
	  bestLH = tr->likelihood;
	  bestIndex = j;
	}
    }
  /*printf("processID = %d, bestIndex = %d; bestLH = %f\n", processID, bestIndex, bestLH);*/
#else
  for(i = 0; i < slowSearches; i++)
    { 
      restoreTL(rl, tr, i);
      resetBranches(tr);

      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 2);
      
      printBothOpen("Slow ML Search %d Likelihood: %f\n", i, tr->likelihood);
      
      if(tr->likelihood > bestLH)
	{
	  bestLH = tr->likelihood;
	  bestIndex = i;
	}
    }
#endif
  
  printBothOpenMPI("Slow ML optimization finished\n\n");

  t = gettime() - t;

#ifdef _WAYNE_MPI
  printBothOpen("Slow ML search on Process %d: Time %f seconds\n", processID, t);
#else
  printBothOpen("Slow ML search Time: %f seconds\n", t);
#endif
  
  t = gettime();
  
  restoreTL(rl, tr, bestIndex);
  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);
 
  treeEvaluate(tr, 2); 
         
  Thorough = 1;
  tr->doCutoff = FALSE;  
	 
  treeOptimizeThorough(tr, 1, 10);
  evaluateGenericInitrav(tr, tr->start);
  
  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, FALSE);
  t = gettime() - t;

#ifdef _WAYNE_MPI
  printBothOpen("Thorough ML search on Process %d: Time %f seconds\n", processID, t);
#else
  printBothOpen("Thorough ML search Time: %f seconds\n", t);
#endif

#ifdef _WAYNE_MPI
  bestLH = tr->likelihood;

  printf("\nprocessID = %d, bestLH = %f\n", processID,  bestLH);

  if(processes > 1)
    {
      double *buffer;
      int bestProcess;

      buffer = (double *)malloc(sizeof(double) * processes);
      for(i = 0; i < processes; i++)
        buffer[i] = unlikely;
      buffer[processID] = bestLH;
      for(i = 0; i < processes; i++)
        MPI_Bcast(&buffer[i], 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
      bestLH = buffer[0];
      bestProcess = 0;
      for(i = 1; i < processes; i++)
        if(buffer[i] > bestLH)
          {
             bestLH = buffer[i];
             bestProcess = i;
          }
      free(buffer);

      if(processID != bestProcess)
        {
          MPI_Finalize();
          exit(0);
        }
    }
#endif

  printBothOpen("\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);   
  printBothOpen("\nModel Information:\n\n");
  
  printModelParams(tr, adef);    
  
  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_bestTree.");
  strcat(bestTreeFileName,         run_id);
   
  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE);
  f = myfopen(bestTreeFileName, "wb");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);

  if(adef->perGeneBranchLengths)
    printTreePerGene(tr, adef, bestTreeFileName, "w");

  
  overallTime = gettime() - masterTime;
  mlTime    = gettime() - mlTime;

  printBothOpen("\nML search took %f secs or %f hours\n", mlTime, mlTime / 3600.0); 
  printBothOpen("\nCombined Bootstrap and ML search took %f secs or %f hours\n", overallTime, overallTime / 3600.0);   
  printBothOpen("\nDrawing Bootstrap Support Values on best-scoring ML tree ...\n\n");
      
  
  freeTL(rl);   
  free(rl);       
  
  calcBipartitions(tr, adef, bestTreeFileName, bootstrapFileName);    
  

  overallTime = gettime() - masterTime;

  printBothOpen("Program execution info written to %s\n", infoFileName);
  printBothOpen("All %d bootstrapped trees written to: %s\n\n", bootstrapsPerformed, bootstrapFileName);
  printBothOpen("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    printBothOpen("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
		  tr->NumberOfModels - 1);    
  printBothOpen("Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
  printBothOpen("Best-scoring ML tree with support values as branch labels written to: %s\n\n", bipartitionsFileNameBranchLabels);
  printBothOpen("Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);

#ifdef _WAYNE_MPI
  MPI_Finalize();
#endif      

  exit(0); 
}


/*******************************************EXPERIMENTAL FUNCTIONS END *****************************************************/





void doBootstrap(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int 
    bootstrapsPerformed,
    i, 
    n, 
    treeVectorLength = -1;
  unsigned int
    vLength = 0;

#ifdef _WAYNE_MPI 
  int 
    j,
    bootStopTests = 1;
#endif

  double loopTime, pearsonAverage;
  hashtable *h = (hashtable*)NULL;
  unsigned int **bitVectors = (unsigned int **)NULL;
  boolean bootStopIt = FALSE; 
  

  n = adef->multipleRuns; 

#ifdef _WAYNE_MPI
  if(n % processes != 0)
    n = processes * ((n / processes) + 1);
  adef->multipleRuns = n;
#endif

  if(adef->bootStopping)
    {    
      h = initHashTable(tr->mxtips * 100);
      bitVectors = initBitVector(tr, &vLength);    

      treeVectorLength = adef->multipleRuns;        
    }           

#ifdef _WAYNE_MPI
  long parsimonySeed0 = adef->parsimonySeed;
  long replicateSeed0 = adef->rapidBoot;
  n = n / processes;
#endif

  for(i = 0; i < n && !bootStopIt; i++)
    {    
      loopTime = gettime();
                    
#ifdef _WAYNE_MPI
      if(i == 0)
        {
          if(parsimonySeed0 != 0)
            adef->parsimonySeed = parsimonySeed0 + 10000 * processID;
          adef->rapidBoot = replicateSeed0 + 10000 * processID;
        }
      j  = i + n*processID;
      singleBootstrap(tr, j, adef, rdta, cdta);
#else      
      singleBootstrap(tr, i, adef, rdta, cdta);              
#endif
      loopTime = gettime() - loopTime;     
      writeInfoFile(adef, tr, loopTime);  

      if(adef->bootStopping)
#ifdef _WAYNE_MPI
	{
	  int 
	    nn = (i + 1) * processes;

	  if((nn > START_BSTOP_TEST) && 
	     (i * processes < FC_SPACING * bootStopTests) &&
	     ((i + 1) * processes >= FC_SPACING * bootStopTests)
	     )	     
	    {
	      MPI_Barrier(MPI_COMM_WORLD);
	                    
	      concatenateBSFiles(processes);               
	      
              MPI_Barrier(MPI_COMM_WORLD);
	     
	      bootStopIt = computeBootStopMPI(tr, bootstrapFileName, adef, &pearsonAverage);
	      bootStopTests++;
	    }
	}
#else	      	
      bootStopIt = bootStop(tr, h, i, &pearsonAverage, bitVectors, treeVectorLength, vLength);
#endif
    }      

#ifdef _WAYNE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  
  bootstrapsPerformed = i * processes;
  
  if(processID == 0)
    {      
      if(!adef->bootStopping)
	concatenateBSFiles(processes);

      removeBSFiles(processes);      
    }
  
  MPI_Barrier(MPI_COMM_WORLD); 
#else
  bootstrapsPerformed = i;
#endif
  
  adef->multipleRuns = bootstrapsPerformed;
  
  if(adef->bootStopping)
    {
      freeBitVectors(bitVectors, 2 * tr->mxtips);
      free(bitVectors);
      freeHashTable(h);
      free(h);
      
       
      if(bootStopIt)
	{
	  switch(tr->bootStopCriterion)
	    {
	    case FREQUENCY_STOP:
	      printBothOpenMPI("Stopped Standard BS search after %d replicates with FC Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);	      
	      break;
	    case MR_STOP:
	      printBothOpenMPI("Stopped Standard BS search after %d replicates with MR-based Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);	     
	      break;
	    case MRE_STOP:
	      printBothOpenMPI("Stopped Standard BS search after %d replicates with MRE-based Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);	     
	      break;
	    case MRE_IGN_STOP:
	      printBothOpenMPI("Stopped Standard BS search after %d replicates with MRE_IGN-based Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);	     
	      break;
	    default:
	      assert(0);
	    }
	}
      else
	{
	  switch(tr->bootStopCriterion)
	    {
	    case FREQUENCY_STOP:
	      printBothOpenMPI("Standard BS search did not converge after %d replicates with FC Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("Pearson Average of %d random splits: %f\n",BOOTSTOP_PERMUTATIONS , pearsonAverage);
	      break;
	    case MR_STOP:
	      printBothOpenMPI("Standard BS search did not converge after %d replicates with MR-based Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    case MRE_STOP:
	      printBothOpenMPI("Standard BS search did not converge after %d replicates with MRE-based Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    case MRE_IGN_STOP:
	      printBothOpenMPI("Standard BS search did not converge after %d replicates with MR_IGN-based Bootstopping criterion\n", bootstrapsPerformed);
	      printBothOpenMPI("WRF Average of %d random splits: %f\n", BOOTSTOP_PERMUTATIONS, pearsonAverage);
	      break;
	    default:
	      assert(0);
	    }
	}     
    }
}

void doInference(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n;

#ifdef _WAYNE_MPI
  int j;
#endif

  double loopTime;
  topolRELL_LIST *rl = (topolRELL_LIST *)NULL; 
  int 
    best = -1,
    newBest = -1;
  double 
    bestLH = unlikely; 
  FILE *f;
  char bestTreeFileName[1024]; 
  double overallTime;

  n = adef->multipleRuns;
     
#ifdef _WAYNE_MPI
  if(n % processes != 0)
    n = processes * ((n / processes) + 1);
#endif 

  if(!tr->catOnly)
    {
      rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
      initTL(rl, tr, n);
    }

#ifdef _WAYNE_MPI
  long parsimonySeed0 = adef->parsimonySeed;
  n = n / processes;
#endif

  for(i = 0; i < n; i++)
    { 
#ifdef _WAYNE_MPI 
      if(i == 0)
        { 
          if(parsimonySeed0 != 0) 
            adef->parsimonySeed = parsimonySeed0 + 10000 * processID;
        }
      j = i + n * processID;
      tr->treeID = j;
#else    
      tr->treeID = i;
#endif

      tr->checkPointCounter = 0;
         
      loopTime = gettime();
                                             
      initModel(tr, rdta, cdta, adef); 

      if(i == 0)
	printBaseFrequencies(tr);
     
      getStartingTree(tr, adef); 

      if(i == 0)
	testGapped(tr);
                       
      computeBIGRAPID(tr, adef, TRUE);  

#ifdef _WAYNE_MPI
      if(tr->likelihood > bestLH)
	{
	  best = j;
	  bestLH = tr->likelihood;
	}

      if(!tr->catOnly)
	saveTL(rl, tr, j);
#else
      if(tr->likelihood > bestLH)
	{
	  best = i;
	  bestLH = tr->likelihood;
	}

      if(!tr->catOnly)
	saveTL(rl, tr, i);
#endif

      loopTime = gettime() - loopTime; 
      writeInfoFile(adef, tr, loopTime);
     
    }     
 
  assert(best >= 0);

#ifdef _WAYNE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  n = n * processes;
#endif

  if(tr->catOnly)
    {
      printBothOpenMPI("\n\nNOT conducting any final model optimizations on all %d trees under CAT-based model ....\n", n);
      printBothOpenMPI("\nREMEMBER that CAT-based likelihood scores are meaningless!\n\n", n);        
#ifdef _WAYNE_MPI
      if(processID != 0)
        {
          MPI_Finalize();
          exit(0);
        }
#endif
    }
  else
    {
      printBothOpenMPI("\n\nConducting final model optimizations on all %d trees under GAMMA-based models ....\n\n", n);
 
#ifdef _WAYNE_MPI
      n = n / processes;
#endif

      if(tr->rateHetModel == GAMMA ||  tr->rateHetModel == GAMMA_I)
	{
	  restoreTL(rl, tr, best);
	  onlyInitrav(tr, tr->start);
	  if(!adef->useBinaryModelFile)
	    modOpt(tr, adef, FALSE, adef->likelihoodEpsilon, FALSE); 
	  else
	    {
	      readBinaryModel(tr);
	      evaluateGenericInitrav(tr, tr->start);
	      treeEvaluate(tr, 2);
	    }
	  bestLH = tr->likelihood;
	  tr->likelihoods[best] = tr->likelihood;
	  saveTL(rl, tr, best);
	  tr->treeID = best; 
	  printResult(tr, adef, TRUE);
	  newBest = best;      
	  
	  for(i = 0; i < n; i++)
	    {
#ifdef _WAYNE_MPI
	      j = i + n * processID;
	      if(j != best)
		{
		  restoreTL(rl, tr, j);
		  onlyInitrav(tr, tr->start);
		  treeEvaluate(tr, 1);
		  tr->likelihoods[j] = tr->likelihood;
		  
		  if(tr->likelihood > bestLH)
		    {
		      newBest = j;
		      bestLH = tr->likelihood;		  
		      saveTL(rl, tr, j);
		    }
		  tr->treeID = j;
		  printResult(tr, adef, TRUE);
		}
	      if(n == 1 && processes == 1)
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);	   
	      else	    
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", j, tr->likelihoods[j], resultFileName, j);
#else	  
	      if(i != best)
		{
		  restoreTL(rl, tr, i);
		  onlyInitrav(tr, tr->start);
		  treeEvaluate(tr, 1);
		  tr->likelihoods[i] = tr->likelihood;
		  
		  if(tr->likelihood > bestLH)
		    {
		      newBest = i;
		      bestLH = tr->likelihood;		  
		      saveTL(rl, tr, i);
		    }
		  tr->treeID = i;
		  printResult(tr, adef, TRUE);
		}

	      
	      if(n == 1)
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);	   
	      else	    
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", i, tr->likelihoods[i], resultFileName, i);
#endif	    	 
	    }    
	}
      else
	{     
	  catToGamma(tr, adef);
	  
#ifdef _WAYNE_MPI
	  for(i = 0; i < n; i++)
            {
              j = i + n*processID;
	      rl->t[j]->likelihood = unlikely;
            }  
#else
	  for(i = 0; i < n; i++)
	    rl->t[i]->likelihood = unlikely;
#endif
	  
	  initModel(tr, rdta, cdta, adef);
	  
	  restoreTL(rl, tr, best);      
	  
	  resetBranches(tr);
	  onlyInitrav(tr, tr->start);
	  modOpt(tr, adef, TRUE, adef->likelihoodEpsilon, TRUE);      
	  tr->likelihoods[best] = tr->likelihood;
	  bestLH = tr->likelihood;     
	  saveTL(rl, tr, best);
	  tr->treeID = best;
	  printResult(tr, adef, TRUE);
	  newBest = best;
	  
	  for(i = 0; i < n; i++)
	    {
#ifdef _WAYNE_MPI
	      j = i + n*processID;
	      if(j != best)
		{
		  restoreTL(rl, tr, j);	    
		  resetBranches(tr);
		  onlyInitrav(tr, tr->start);
		  treeEvaluate(tr, 2);
		  tr->likelihoods[j] = tr->likelihood;
		  
		  if(tr->likelihood > bestLH)
		    { 
		      newBest = j;
		      bestLH = tr->likelihood;		
		      saveTL(rl, tr, j);	  
		    }
		  tr->treeID = j;
		  printResult(tr, adef, TRUE);
		} 
	      
	      if(n == 1 && processes == 1)	    
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);
	      else
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", j, tr->likelihoods[j], resultFileName, j);
#else
	      if(i != best)
		{
		  restoreTL(rl, tr, i);	    
		  resetBranches(tr);
		  onlyInitrav(tr, tr->start);
		  treeEvaluate(tr, 2);
		  tr->likelihoods[i] = tr->likelihood;
		  
		  if(tr->likelihood > bestLH)
		    { 
		      newBest = i;
		      bestLH = tr->likelihood;		
		      saveTL(rl, tr, i);	  
		    }
		  tr->treeID = i;
		  printResult(tr, adef, TRUE);
		} 
	      
	      if(n == 1)	    
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s\n", i, tr->likelihoods[i], resultFileName);
	      else
		printBothOpen("Inference[%d] final GAMMA-based Likelihood: %f tree written to file %s.RUN.%d\n", i, tr->likelihoods[i], resultFileName, i);	   	  
#endif
	    }
	}     
    
      assert(newBest >= 0);

#ifdef _WAYNE_MPI
      if(processes > 1)
	{
	  double *buffer;
	  int bestProcess;
	  
	  buffer = (double *)malloc(sizeof(double) * processes);
	  for(i = 0; i < processes; i++)
	    buffer[i] = unlikely;
	  buffer[processID] = bestLH;
	  for(i = 0; i < processes; i++)
	    MPI_Bcast(&buffer[i], 1, MPI_DOUBLE, i, MPI_COMM_WORLD);
	  bestLH = buffer[0];
	  bestProcess = 0;
	  for(i = 1; i < processes; i++)
	    if(buffer[i] > bestLH)
	      {
		bestLH = buffer[i];
		bestProcess = i;
	      }
	  free(buffer);
	  
	  if(processID != bestProcess)
	    {
	      MPI_Finalize();
	      exit(0);
	    }
	}
#endif

      restoreTL(rl, tr, newBest);
  
      onlyInitrav(tr, tr->start);
      printBothOpen("\n\nStarting final GAMMA-based thorough Optimization on tree %d likelihood %f .... \n\n", newBest, tr->likelihoods[newBest]);

      Thorough = 1;
      tr->doCutoff = FALSE; 
      treeOptimizeThorough(tr, 1, 10); 
      evaluateGenericInitrav(tr, tr->start);
  
      printBothOpen("Final GAMMA-based Score of best tree %f\n\n", tr->likelihood); 
    

      strcpy(bestTreeFileName, workdir); 
      strcat(bestTreeFileName, "RAxML_bestTree.");
      strcat(bestTreeFileName,         run_id);
      
     
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH, FALSE, FALSE);
      
      f = myfopen(bestTreeFileName, "wb");
      fprintf(f, "%s", tr->tree_string);
      fclose(f);

      if(adef->perGeneBranchLengths)
	printTreePerGene(tr, adef, bestTreeFileName, "w");
    }
  
  overallTime = gettime() - masterTime;

  printBothOpen("Program execution info written to %s\n", infoFileName);
  
  if(!tr->catOnly)
    {
      printBothOpen("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);

      if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
	printBothOpen("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
		      tr->NumberOfModels - 1);  
    }
   
  printBothOpen("Overall execution time: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);    

  if(!tr->catOnly)
    {
      freeTL(rl);   
      free(rl); 
    }
  
#ifdef _WAYNE_MPI
  MPI_Finalize();
#endif
  exit(0);
}

