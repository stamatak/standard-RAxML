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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
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


static const double MNBRAK_GOLD =    1.618034;
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =      1.e-5;
static const double BRENT_CGOLD =   0.3819660;

extern int optimizeRatesInvocations;
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];
extern char *protModels[NUM_PROT_MODELS];


#ifdef _USE_PTHREADS
extern volatile int             NumberOfThreads;
extern volatile double          *reductionBuffer;
#endif

/* TODO remove at some point */
#define _DEBUG_MODEL_OPTIMIZATION 

#define ALPHA_F    0
#define INVAR_F    1
#define RATE_F     2
#define SCALER_F   3
#define LXRATE_F   4
#define LXWEIGHT_F 5
#define FREQ_F     6


static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, tree *tr, linkageList *ll, double lim_inf, double lim_sup);

static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double lim_inf, double lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, linkageList *ll, double modelEpsilon);

static void optParamGeneric(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double lim_inf, double lim_sup, int whichParameterType);

static void updateWeights(tree *tr, int model, int rate, double value);

/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/




static void setRateModel(tree *tr, int model, double rate, int position)
{
  int
    states   = tr->partitionData[model].states,
    numRates = (states * states - states) / 2;

  if(tr->partitionData[model].dataType == DNA_DATA)
    assert(position >= 0 && position < (numRates - 1));
  else
    assert(position >= 0 && position < numRates);

  assert(tr->partitionData[model].dataType != BINARY_DATA); 

  if(!(tr->partitionData[model].dataType == SECONDARY_DATA ||
       tr->partitionData[model].dataType == SECONDARY_DATA_6 ||
       tr->partitionData[model].dataType == SECONDARY_DATA_7))
    assert(rate >= RATE_MIN && rate <= RATE_MAX);

  if(tr->partitionData[model].nonGTR)
    {    
      int 
	i, 
	k = tr->partitionData[model].symmetryVector[position];

      assert(tr->partitionData[model].dataType == SECONDARY_DATA ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_6 ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_7);

      if(k == -1)
	tr->partitionData[model].substRates[position] = 0.0;
      else
	{
	  if(k == tr->partitionData[model].symmetryVector[numRates - 1])
	    {
	      for(i = 0; i < numRates - 1; i++)
		if(tr->partitionData[model].symmetryVector[i] == k)
		  tr->partitionData[model].substRates[position] = 1.0;
	    }
	  else
	    {
	      for(i = 0; i < numRates - 1; i++)
		{
		  if(tr->partitionData[model].symmetryVector[i] == k)
		    tr->partitionData[model].substRates[i] = rate; 
		}	      	     
	    }
	}
    }
  else
    tr->partitionData[model].substRates[position] = rate;
}





static linkageList* initLinkageList(int *linkList, tree *tr)
{
  int 
    k,
    partitions,
    numberOfModels = 0,
    i,
    pos;
  
  linkageList* 
    ll = (linkageList*)rax_malloc(sizeof(linkageList));
      
  for(i = 0; i < tr->NumberOfModels; i++)
    {
      if(linkList[i] > numberOfModels)
	numberOfModels = linkList[i];
    }

  numberOfModels++;
  
  ll->entries = numberOfModels;
  ll->ld      = (linkageData*)rax_malloc(sizeof(linkageData) * numberOfModels);


  for(i = 0; i < numberOfModels; i++)
    {
      ll->ld[i].valid = TRUE;
      partitions = 0;

      for(k = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  partitions++;	    

      ll->ld[i].partitions = partitions;
      ll->ld[i].partitionList = (int*)rax_malloc(sizeof(int) * partitions);
      
      for(k = 0, pos = 0; k < tr->NumberOfModels; k++)	
	if(linkList[k] == i)
	  ll->ld[i].partitionList[pos++] = k;
    }

  return ll;
}


static linkageList* initLinkageListGTR(tree *tr)
{
  int
    i,
    *links = (int*)rax_malloc(sizeof(int) * tr->NumberOfModels),
    firstAA = tr->NumberOfModels + 2,
    countGTR = 0,
    countUnlinkedGTR = 0,
    countOtherModel = 0;
  
  linkageList* 
    ll;

  for(i = 0; i < tr->NumberOfModels; i++)
    {     
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  if(tr->partitionData[i].protModels == GTR)
	    {
	      if(i < firstAA)
		firstAA = i;
	      countGTR++;
	    }
	  else
	    {
	      if(tr->partitionData[i].protModels == GTR_UNLINKED)
		countUnlinkedGTR++;
	      else
		countOtherModel++;
	    }
	}
    }
  
  /* 
     TODO need to think what we actually want here ! 
     Shall we mix everything: linked, unlinked WAG etc?
  */

  assert((countGTR > 0 && countOtherModel == 0) || (countGTR == 0 && countOtherModel > 0) ||  (countGTR == 0 && countOtherModel == 0));

  if(countGTR == 0)
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	links[i] = i;
    }
  else
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  switch(tr->partitionData[i].dataType)
	    {	   
	    case DNA_DATA:
	    case BINARY_DATA:
	    case GENERIC_32:
	    case GENERIC_64:
	    case SECONDARY_DATA:
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7: 
	      links[i] = i;
	      break;
	    case AA_DATA:	  
	      links[i] = firstAA;
	      break;
	    default:
	      assert(0);
	    }
	}
    }
  

  ll = initLinkageList(links, tr);

  rax_free(links);
  
  return ll;
}

static void changeModelParameters(int index, int rateNumber, double value, int whichParameterType, tree *tr)
{
  switch(whichParameterType)
    {
    case RATE_F:
      setRateModel(tr, index, value, rateNumber);  
      initReversibleGTR(tr, index);		 
      break;
    case ALPHA_F:
      tr->partitionData[index].alpha = value;
      makeGammaCats(tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4, tr->useGammaMedian);
      break;
    case INVAR_F:
      tr->partitionData[index].propInvariant = value;
      break;
    case SCALER_F:     
      tr->partitionData[index].brLenScaler = value;
      scaleBranches(tr, FALSE);			        
      break;
    case LXRATE_F:
      tr->partitionData[index].gammaRates[rateNumber] = value;
      break;
    case LXWEIGHT_F:
      updateWeights(tr, index, rateNumber, value);
      break;	
    case FREQ_F:
      {
	int
	  states = tr->partitionData[index].states,
	  j;

	double 
	  w = 0.0;

	tr->partitionData[index].freqExponents[rateNumber] = value;

	for(j = 0; j < states; j++)
	  w += exp(tr->partitionData[index].freqExponents[j]);

	for(j = 0; j < states; j++)	    	    
	  tr->partitionData[index].frequencies[j] = exp(tr->partitionData[index].freqExponents[j]) / w;
	
	/*
	  for(j = 0; j < states; j++)
	  printf("%f ", tr->partitionData[index].frequencies[j]);
	  printf("\n");
	*/

	initReversibleGTR(tr, index);
      }
      break;
    default:
      assert(0);
    }
}


static void freeLinkageList( linkageList* ll)
{
  int i;    

  for(i = 0; i < ll->entries; i++)    
    rax_free(ll->ld[i].partitionList);         

  rax_free(ll->ld);
  rax_free(ll);   
}


static void updateWeights(tree *tr, int model, int rate, double value)
{
  int 
    j;

  double 
    w = 0.0;

  tr->partitionData[model].weightExponents[rate] = value;

  for(j = 0; j < 4; j++)
    w += exp(tr->partitionData[model].weightExponents[j]);

  for(j = 0; j < 4; j++)	    	    
    tr->partitionData[model].weights[j] = exp(tr->partitionData[model].weightExponents[j]) / w;
}


static void optimizeWeights(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels)
{
  int 
    i;
  
  double 
    initialLH = 0.0,
    finalLH   = 0.0;

  evaluateGeneric(tr, tr->start);
 
  initialLH = tr->likelihood;
  //printf("W: %f %f [%f] ->", tr->perPartitionLH[0], tr->perPartitionLH[1], initialLH);

  for(i = 0; i < 4; i++)   
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, i, -1000000.0, 200.0, LXWEIGHT_F);
    //optLG4X_Weights(tr, ll, numberOfModels, i, modelEpsilon);
  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_LG4X_RATES, tr);
#endif

  evaluateGenericInitrav(tr, tr->start); 

  finalLH = tr->likelihood;

  if(finalLH < initialLH)
    printf("Final: %f initial: %f\n", finalLH, initialLH);
  assert(finalLH >= initialLH);

  //printf("%f %f [%f]\n",  tr->perPartitionLH[0], tr->perPartitionLH[1], finalLH);
}
static void evaluateChange(tree *tr, int rateNumber, double *value, double *result, boolean* converged, int whichFunction, int numberOfModels, linkageList *ll, double modelEpsilon)
{ 
  int 
    i, 
    k, 
    pos;

  boolean
    atLeastOnePartition = FALSE;
   
  for(i = 0, pos = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  if(converged[pos])
	    {
	      //if parameter optimizations for this specific model have converged 
	      //set executeModel to FALSE 

	      for(k = 0; k < ll->ld[i].partitions; k++)
		tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;
	    }
	  else
	    {
	      atLeastOnePartition = TRUE;

	      for(k = 0; k < ll->ld[i].partitions; k++)
		{
		  int 
		    index = ll->ld[i].partitionList[k];

		  changeModelParameters(index, rateNumber, value[pos], whichFunction, tr);		    		  
		}
	    }
	  pos++;
	}
      else
	{
	  // if this partition is not being optimized anyway (e.g., we may be optimizing GTR rates for all DNA partitions,
	  // but there are also a couple of Protein partitions with fixed models like WAG, JTT, etc.) set executeModel to FALSE
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    tr->executeModel[ll->ld[i].partitionList[k]] = FALSE;	     
	}      
    }

  assert(pos == numberOfModels);

  //some error checks for individual model parameters
  //and individual pre-processing 

  switch(whichFunction)
    {
    case SCALER_F:      
      assert(ll->entries == tr->NumberOfModels);
      assert(ll->entries == tr->numBranches);
      scaleBranches(tr, FALSE);
      break;
    case RATE_F: 
      assert(rateNumber != -1);  
      if(tr->useBrLenScaler)
	determineFullTraversal(tr->start, tr);
      break;
    case ALPHA_F:	     
      break;
    case INVAR_F:     
      break;
    case LXRATE_F:
      assert(rateNumber != -1);
    case LXWEIGHT_F:
      assert(rateNumber != -1);
      break;
    case FREQ_F:
      break;
    default:
      assert(0);
    }	
  

  

#ifdef _USE_PTHREADS
  switch(whichFunction)
    {
    case RATE_F:
      masterBarrier(THREAD_OPT_RATE, tr);
      break;
    case ALPHA_F:	
      masterBarrier(THREAD_OPT_ALPHA, tr);
      break;
    case INVAR_F:
      masterBarrier(THREAD_OPT_INVAR, tr);
      break;
    case SCALER_F:
      determineFullTraversal(tr->start, tr);	
      masterBarrier(THREAD_OPT_SCALER, tr);
      break;
    case LXRATE_F:
      masterBarrier(THREAD_OPT_LG4X_RATES, tr);
      break;
    case LXWEIGHT_F:
      masterBarrier(THREAD_OPT_LG4X_RATES, tr);
      break;
    case FREQ_F:
      masterBarrier(THREAD_OPT_RATE, tr);
      break;
    default:
      assert(0);
    }	
  {
    volatile double 
      result,
      partitionResult;
		    
    int 
      j;
	
    result = 0.0;
    
    for(j = 0; j < tr->NumberOfModels; j++)
      {
	for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
	  partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
	result +=  partitionResult;
	tr->perPartitionLH[j] = partitionResult;
      }
  }
#else
  switch(whichFunction)
    {
    case RATE_F:
    case ALPHA_F:
    case SCALER_F:
    case LXRATE_F: 
    case FREQ_F:
      evaluateGenericInitrav(tr, tr->start);
      //printf("%f \n", tr->likelihood);
      break;
    case LXWEIGHT_F:
    case INVAR_F:
      evaluateGeneric(tr, tr->start);  
      break;
    default:
      assert(0);
    }
#endif

 

  //nested optimization for LX4 model, now optimize the weights!

  if(whichFunction == LXRATE_F && atLeastOnePartition)
    {
      boolean 
	*buffer = (boolean*)rax_malloc(tr->NumberOfModels * sizeof(boolean));
	    
      memcpy(buffer, tr->executeModel, sizeof(boolean) * tr->NumberOfModels);
	    
      for(i = 0; i < tr->NumberOfModels; i++)
	tr->executeModel[i] = FALSE;
	    
      for(i = 0, pos = 0; i < ll->entries; i++)	
	{  
	  int 
	    index = ll->ld[i].partitionList[0];	    	      
	    
	  if(ll->ld[i].valid)		  	    	      		   	    
	    tr->executeModel[index] = TRUE;	    
	}

      optimizeWeights(tr, modelEpsilon, ll, numberOfModels);      
	    
      memcpy(tr->executeModel, buffer, sizeof(boolean) * tr->NumberOfModels);
	    
      rax_free(buffer);
    }
    
  
  for(i = 0, pos = 0; i < ll->entries; i++)	
    {
      if(ll->ld[i].valid)
	{
	  result[pos] = 0.0;
	  
	  for(k = 0; k < ll->ld[i].partitions; k++)
	    {
	      int 
		index = ll->ld[i].partitionList[k];

	      assert(tr->perPartitionLH[index] <= 0.0);
	      
	      result[pos] -= tr->perPartitionLH[index];
	      
	    }
	  pos++;
	}

      //set execute model for ALL partitions to true again 
      //for consistency 

      for(k = 0; k < ll->ld[i].partitions; k++)
	{
	  int 
	    index = ll->ld[i].partitionList[k];	  
	  tr->executeModel[index] = TRUE;
	}	  
    }
  
  assert(pos == numberOfModels);   
}




static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, tree *tr, linkageList *ll, double lim_inf, double lim_sup)
{
  int iter, i;
  double 
    *a     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *b     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *d     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *etemp = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fu    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fv    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fw    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fx    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *p     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *q     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *r     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *tol1  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *tol2  = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *u     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *v     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *w     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *x     = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *xm    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *e     = (double *)rax_malloc(sizeof(double) * numberOfModels);
  boolean *converged = (boolean *)rax_malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;
  
  for(i = 0; i < numberOfModels; i++)    
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      e[i] = 0.0;
      d[i] = 0.0;
    }

  for(i = 0; i < numberOfModels; i++)
    {
      a[i]=((ax[i] < cx[i]) ? ax[i] : cx[i]);
      b[i]=((ax[i] > cx[i]) ? ax[i] : cx[i]);
      x[i] = w[i] = v[i] = bx[i];
      fw[i] = fv[i] = fx[i] = fb[i];
    }

  for(i = 0; i < numberOfModels; i++)
    {      
      assert(a[i] >= lim_inf && a[i] <= lim_sup);
      assert(b[i] >= lim_inf && b[i] <= lim_sup);
      assert(x[i] >= lim_inf && x[i] <= lim_sup);
      assert(v[i] >= lim_inf && v[i] <= lim_sup);
      assert(w[i] >= lim_inf && w[i] <= lim_sup);
    }
  
  

  for(iter = 1; iter <= ITMAX; iter++)
    {
      allConverged = TRUE;

      for(i = 0; i < numberOfModels && allConverged; i++)
	allConverged = allConverged && converged[i];

      if(allConverged)
	{
	  rax_free(converged);
	  rax_free(a);
	  rax_free(b);
	  rax_free(d);
	  rax_free(etemp);
	  rax_free(fu);
	  rax_free(fv);
	  rax_free(fw);
	  rax_free(fx);
	  rax_free(p);
	  rax_free(q);
	  rax_free(r);
	  rax_free(tol1);
	  rax_free(tol2);
	  rax_free(u);
	  rax_free(v);
	  rax_free(w);
	  rax_free(x);
	  rax_free(xm);
	  rax_free(e);
	  return;
	}     

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {	     	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
  
	      xm[i] = 0.5 * (a[i] + b[i]);
	      tol2[i] = 2.0 * (tol1[i] = tol * fabs(x[i]) + BRENT_ZEPS);
	  
	      if(fabs(x[i] - xm[i]) <= (tol2[i] - 0.5 * (b[i] - a[i])))
		{		 
		  result[i] =  -fx[i];
		  xmin[i]   = x[i];
		  converged[i] = TRUE;		  
		}
	      else
		{
		  if(fabs(e[i]) > tol1[i])
		    {		     
		      r[i] = (x[i] - w[i]) * (fx[i] - fv[i]);
		      q[i] = (x[i] - v[i]) * (fx[i] - fw[i]);
		      p[i] = (x[i] - v[i]) * q[i] - (x[i] - w[i]) * r[i];
		      q[i] = 2.0 * (q[i] - r[i]);
		      if(q[i] > 0.0)
			p[i] = -p[i];
		      q[i] = fabs(q[i]);
		      etemp[i] = e[i];
		      e[i] = d[i];
		      if((fabs(p[i]) >= fabs(0.5 * q[i] * etemp[i])) || (p[i] <= q[i] * (a[i]-x[i])) || (p[i] >= q[i] * (b[i] - x[i])))
			d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i] : b[i] - x[i]));
		      else
			{
			  d[i] = p[i] / q[i];
			  u[i] = x[i] + d[i];
			  if( u[i] - a[i] < tol2[i] || b[i] - u[i] < tol2[i])
			    d[i] = SIGN(tol1[i], xm[i] - x[i]);
			}
		    }
		  else
		    {		     
		      d[i] = BRENT_CGOLD * (e[i] = (x[i] >= xm[i] ? a[i] - x[i]: b[i] - x[i]));
		    }
		  u[i] = ((fabs(d[i]) >= tol1[i]) ? (x[i] + d[i]): (x[i] +SIGN(tol1[i], d[i])));
		}

	      if(!converged[i])
		assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
                 
      evaluateChange(tr, rateNumber, u, fu, converged, whichFunction, numberOfModels, ll, tol);

      for(i = 0; i < numberOfModels; i++)
	{
	  if(!converged[i])
	    {
	      if(fu[i] <= fx[i])
		{
		  if(u[i] >= x[i])
		    a[i] = x[i];
		  else
		    b[i] = x[i];
		  
		  SHFT(v[i],w[i],x[i],u[i]);
		  SHFT(fv[i],fw[i],fx[i],fu[i]);
		}
	      else
		{
		  if(u[i] < x[i])
		    a[i] = u[i];
		  else
		    b[i] = u[i];
		  
		  if(fu[i] <= fw[i] || w[i] == x[i])
		    {
		      v[i] = w[i];
		      w[i] = u[i];
		      fv[i] = fw[i];
		      fw[i] = fu[i];
		    }
		  else
		    {
		      if(fu[i] <= fv[i] || v[i] == x[i] || v[i] == w[i])
			{
			  v[i] = u[i];
			  fv[i] = fu[i];
			}
		    }	    
		}
	      
	      assert(a[i] >= lim_inf && a[i] <= lim_sup);
	      assert(b[i] >= lim_inf && b[i] <= lim_sup);
	      assert(x[i] >= lim_inf && x[i] <= lim_sup);
	      assert(v[i] >= lim_inf && v[i] <= lim_sup);
	      assert(w[i] >= lim_inf && w[i] <= lim_sup);
	      assert(u[i] >= lim_inf && u[i] <= lim_sup);
	    }
	}
    }

  rax_free(converged);
  rax_free(a);
  rax_free(b);
  rax_free(d);
  rax_free(etemp);
  rax_free(fu);
  rax_free(fv);
  rax_free(fw);
  rax_free(fx);
  rax_free(p);
  rax_free(q);
  rax_free(r);
  rax_free(tol1);
  rax_free(tol2);
  rax_free(u);
  rax_free(v);
  rax_free(w);
  rax_free(x);
  rax_free(xm);
  rax_free(e);

  printf("\n. Too many iterations in BRENT !");
  assert(0);
}



static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double lim_inf, double lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, linkageList *ll, double modelEpsilon)
{
  double 
    *ulim = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *u    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *r    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *q    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *fu   = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *dum  = (double *)rax_malloc(sizeof(double) * numberOfModels), 
    *temp = (double *)rax_malloc(sizeof(double) * numberOfModels);
  
  int 
    i,
    *state    = (int *)rax_malloc(sizeof(int) * numberOfModels),
    *endState = (int *)rax_malloc(sizeof(int) * numberOfModels);

  boolean *converged = (boolean *)rax_malloc(sizeof(boolean) * numberOfModels);
  boolean allConverged;

  for(i = 0; i < numberOfModels; i++)
    converged[i] = FALSE;

  for(i = 0; i < numberOfModels; i++)
    {
      state[i] = 0;
      endState[i] = 0;

      u[i] = 0.0;

      param[i] = ax[i];

      if(param[i] > lim_sup) 	
	param[i] = ax[i] = lim_sup;
      
      if(param[i] < lim_inf) 
	param[i] = ax[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
   
  
  evaluateChange(tr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll, modelEpsilon);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup) 
	param[i] = bx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = bx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
  evaluateChange(tr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll, modelEpsilon);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
	{	  
	  SHFT(dum[i],ax[i],bx[i],dum[i]);
	  SHFT(dum[i],fa[i],fb[i],dum[i]);
	}
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup) 
	param[i] = cx[i] = lim_sup;
      if(param[i] < lim_inf) 
	param[i] = cx[i] = lim_inf;

      assert(param[i] >= lim_inf && param[i] <= lim_sup);
    }
  
 
  evaluateChange(tr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll, modelEpsilon);

   while(1) 
     {       
       allConverged = TRUE;

       for(i = 0; i < numberOfModels && allConverged; i++)
	 allConverged = allConverged && converged[i];

       if(allConverged)
	 {
	   for(i = 0; i < numberOfModels; i++)
	     {	       
	       if(ax[i] > lim_sup) 
		 ax[i] = lim_sup;
	       if(ax[i] < lim_inf) 
		 ax[i] = lim_inf;

	       if(bx[i] > lim_sup) 
		 bx[i] = lim_sup;
	       if(bx[i] < lim_inf) 
		 bx[i] = lim_inf;
	       
	       if(cx[i] > lim_sup) 
		 cx[i] = lim_sup;
	       if(cx[i] < lim_inf) 
		 cx[i] = lim_inf;
	     }

	   rax_free(converged);
	   rax_free(ulim);
	   rax_free(u);
	   rax_free(r);
	   rax_free(q);
	   rax_free(fu);
	   rax_free(dum); 
	   rax_free(temp);
	   rax_free(state);   
	   rax_free(endState);
	   return 0;
	   
	 }

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {
	       switch(state[i])
		 {
		 case 0:
		   endState[i] = 0;
		   if(!(fb[i] > fc[i]))		         
		     converged[i] = TRUE;		       		     
		   else
		     {
		   
		       if(ax[i] > lim_sup) 
			 ax[i] = lim_sup;
		       if(ax[i] < lim_inf) 
			 ax[i] = lim_inf;
		       if(bx[i] > lim_sup) 
			 bx[i] = lim_sup;
		       if(bx[i] < lim_inf) 
			 bx[i] = lim_inf;
		       if(cx[i] > lim_sup) 
			 cx[i] = lim_sup;
		       if(cx[i] < lim_inf) 
			 cx[i] = lim_inf;
		       
		       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
		       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
		       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
			 (2.0*SIGN(MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
		       
		       ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
		       
		       if(u[i] > lim_sup) 
			 u[i] = lim_sup;
		       if(u[i] < lim_inf) 
			 u[i] = lim_inf;
		       if(ulim[i] > lim_sup) 
			 ulim[i] = lim_sup;
		       if(ulim[i] < lim_inf) 
			 ulim[i] = lim_inf;
		       
		       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
			 {
			   param[i] = u[i];
			   if(param[i] > lim_sup) 			     
			     param[i] = u[i] = lim_sup;
			   if(param[i] < lim_inf)
			     param[i] = u[i] = lim_inf;
			   endState[i] = 1;
			 }
		       else 
			 {
			   if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
			     {
			       param[i] = u[i];
			       if(param[i] > lim_sup) 
				 param[i] = u[i] = lim_sup;
			       if(param[i] < lim_inf) 
				 param[i] = u[i] = lim_inf;
			       endState[i] = 2;
			     }		  	       
			   else
			     {
			       if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
				 {
				   u[i] = ulim[i];
				   param[i] = u[i];	
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = ulim[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = ulim[i] = lim_inf;
				   endState[i] = 0;
				 }		  		
			       else 
				 {		  
				   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
				   param[i] = u[i];
				   endState[i] = 0;
				   if(param[i] > lim_sup) 
				     param[i] = u[i] = lim_sup;
				   if(param[i] < lim_inf) 
				     param[i] = u[i] = lim_inf;
				 }
			     }	  
			 }
		     }
		   break;
		 case 1:
		   endState[i] = 0;
		   break;
		 case 2:
		   endState[i] = 3;
		   break;
		 default:
		   assert(0);
		 }
	       assert(param[i] >= lim_inf && param[i] <= lim_sup);
	     }
	 }
             
       evaluateChange(tr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll, modelEpsilon);

       for(i = 0; i < numberOfModels; i++)
	 {
	   if(!converged[i])
	     {	       
	       switch(endState[i])
		 {
		 case 0:
		   fu[i] = temp[i];
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 case 1:
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {
		       ax[i]=(bx[i]);
		       bx[i]=u[i];
		       fa[i]=(fb[i]);
		       fb[i]=fu[i]; 
		       converged[i] = TRUE;		      
		     } 
		   else 
		     {
		       if (fu[i] > fb[i]) 
			 {
			   assert(u[i] >= lim_inf && u[i] <= lim_sup);
			   cx[i]=u[i];
			   fc[i]=fu[i];
			   converged[i] = TRUE;			  
			 }
		       else
			 {		   
			   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
			   param[i] = u[i];
			   if(param[i] > lim_sup) {param[i] = u[i] = lim_sup;}
			   if(param[i] < lim_inf) {param[i] = u[i] = lim_inf;}	  
			   state[i] = 1;		 
			 }		  
		     }
		   break;
		 case 2: 
		   fu[i] = temp[i];
		   if (fu[i] < fc[i]) 
		     {		     
		       SHFT(bx[i],cx[i],u[i], cx[i]+MNBRAK_GOLD*(cx[i]-bx[i]));
		       state[i] = 2;
		     }	   
		   else
		     {
		       state[i] = 0;
		       SHFT(ax[i],bx[i],cx[i],u[i]);
		       SHFT(fa[i],fb[i],fc[i],fu[i]);
		     }
		   break;	   
		 case 3:		  
		   SHFT(fb[i],fc[i],fu[i], temp[i]);
		   SHFT(ax[i],bx[i],cx[i],u[i]);
		   SHFT(fa[i],fb[i],fc[i],fu[i]);
		   state[i] = 0;
		   break;
		 default:
		   assert(0);
		 }
	     }
	 }
    }
   

   assert(0);
   rax_free(converged);
   rax_free(ulim);
   rax_free(u);
   rax_free(r);
   rax_free(q);
   rax_free(fu);
   rax_free(dum); 
   rax_free(temp);
   rax_free(state);   
   rax_free(endState);

  

   return(0);
}









static void optInvar(tree *tr, double modelEpsilon, linkageList *ll)
{
  optParamGeneric(tr, modelEpsilon, ll, tr->NumberOfModels, -1, INVAR_MIN, INVAR_MAX, INVAR_F);
}








/*******************************************************************************************************************/
/*******************generic optimization ******************************************************************************************/






static void optParamGeneric(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double lim_inf, double lim_sup, int whichParameterType)
{
  int
    l,
    k, 
    j, 
    pos;
    
  double 
    *startRates   = (double *)rax_malloc(sizeof(double) * numberOfModels * 4),
    *startWeights = (double *)rax_malloc(sizeof(double) * numberOfModels * 4),
    *startExponents = (double *)rax_malloc(sizeof(double) * numberOfModels * 4),
    *startValues = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *startLH    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *endLH      = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_a         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_b         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_c         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fa        = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fb        = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_fc        = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *_param     = (double *)rax_malloc(sizeof(double) * numberOfModels),    
    *_x         = (double *)rax_malloc(sizeof(double) * numberOfModels); 

  
  if(whichParameterType == LXWEIGHT_F)
    evaluateGeneric(tr, tr->start);
  else
    evaluateGenericInitrav(tr, tr->start);
  
#ifdef  _DEBUG_MODEL_OPTIMIZATION
  double
    initialLH = tr->likelihood;
#endif

  /* 
     at this point here every worker has the traversal data it needs for the 
     search 
  */

  for(l = 0, pos = 0; l < ll->entries; l++)
    {
      if(ll->ld[l].valid)
	{
	  endLH[pos] = unlikely;
	  startLH[pos] = 0.0;

	  for(j = 0; j < ll->ld[l].partitions; j++)
	    {
	      int 
		index = ll->ld[l].partitionList[j];
	      
	      startLH[pos] += tr->perPartitionLH[index];
	      
	      switch(whichParameterType)
		{
		case ALPHA_F:
		  startValues[pos] = tr->partitionData[index].alpha;
		  break;
		case RATE_F:
		  startValues[pos] = tr->partitionData[index].substRates[rateNumber];      
		  break;
		case INVAR_F:
		  startValues[pos] = tr->partitionData[index].propInvariant;
		  break;
		case SCALER_F:
		  startValues[pos] = tr->partitionData[index].brLenScaler;
		  break;
		case LXRATE_F:
		  assert(rateNumber >= 0 && rateNumber < 4);
		  startValues[pos] = tr->partitionData[index].gammaRates[rateNumber];
		  memcpy(&startRates[pos * 4],   tr->partitionData[index].gammaRates, 4 * sizeof(double)); 
		  memcpy(&startExponents[pos * 4], tr->partitionData[index].weightExponents, 4 * sizeof(double));
		  memcpy(&startWeights[pos * 4], tr->partitionData[index].weights,    4 * sizeof(double));
		  break;
		case LXWEIGHT_F:
		  assert(rateNumber >= 0 && rateNumber < 4);
		  startValues[pos] = tr->partitionData[index].weightExponents[rateNumber];
		  memcpy(&startRates[pos * 4],   tr->partitionData[index].gammaRates, 4 * sizeof(double)); 
		  memcpy(&startExponents[pos * 4], tr->partitionData[index].weightExponents, 4 * sizeof(double));
		  memcpy(&startWeights[pos * 4], tr->partitionData[index].weights,    4 * sizeof(double));
		  break;
		case FREQ_F:
		  startValues[pos] = tr->partitionData[index].freqExponents[rateNumber];
		  break;
		default:
		  assert(0);
		}
		  
	    }
	  pos++;
	}
    }  

  assert(pos == numberOfModels);
   
  for(k = 0, pos = 0; k < ll->entries; k++)
    {
      if(ll->ld[k].valid)
	{	 	 	  
	  _a[pos] = startValues[pos] + 0.1;
	  _b[pos] = startValues[pos] - 0.1;
	      
	  if(_a[pos] < lim_inf) 
	    _a[pos] = lim_inf;
	  
	  if(_a[pos] > lim_sup) 
	    _a[pos] = lim_sup;
	      
	  if(_b[pos] < lim_inf) 
	    _b[pos] = lim_inf;
	  
	  if(_b[pos] > lim_sup) 
	    _b[pos] = lim_sup;    

	  pos++;
	}
    }                    	     

  assert(pos == numberOfModels);

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, rateNumber, whichParameterType, tr, ll, modelEpsilon);
      
  for(k = 0; k < numberOfModels; k++)
    {
      assert(_a[k] >= lim_inf && _a[k] <= lim_sup);
      assert(_b[k] >= lim_inf && _b[k] <= lim_sup);	  
      assert(_c[k] >= lim_inf && _c[k] <= lim_sup);	    
    }      

  brentGeneric(_a, _b, _c, _fb, modelEpsilon, _x, endLH, numberOfModels, whichParameterType, rateNumber, tr,  ll, lim_inf, lim_sup);
		      
  for(k = 0, pos = 0; k < ll->entries; k++)
    {
      if(ll->ld[k].valid)
	{ 
	  if(startLH[pos] > endLH[pos])
	    {
	      //printf("revert\n");
	      //if the initial likelihood was better than the likelihodo after optimization, we set the values back 
	      //to their original values 

	      for(j = 0; j < ll->ld[k].partitions; j++)
		{
		  int 
		    index = ll->ld[k].partitionList[j];
		  
		  if(whichParameterType == LXWEIGHT_F || whichParameterType == LXRATE_F)
		    {
		      memcpy(tr->partitionData[index].weights,         &startWeights[pos * 4], sizeof(double) * 4);
		      memcpy(tr->partitionData[index].gammaRates,      &startRates[pos * 4], sizeof(double) * 4);
		      memcpy(tr->partitionData[index].weightExponents, &startExponents[pos * 4], 4 * sizeof(double));
		    }

		  changeModelParameters(index, rateNumber, startValues[pos], whichParameterType, tr);		 
		}
	    }
	  else
	    {
	      //	      printf("accept\n");

	      //otherwise we set the value to the optimized value 
	      //this used to be a bug in standard RAxML, before I fixed it 
	      //I was not using _x[pos] as value that needs to be set 

	      for(j = 0; j < ll->ld[k].partitions; j++)
		{
		  int 
		    index = ll->ld[k].partitionList[j];

		  changeModelParameters(index, rateNumber, _x[pos], whichParameterType, tr);		  		  
		}
	    }
	  pos++;
	}
    }

#ifdef _USE_PTHREADS
  switch(whichParameterType)
    {
    case FREQ_F:
    case RATE_F:      
      masterBarrier(THREAD_COPY_RATES, tr);
      break;
    case ALPHA_F:     
      masterBarrier(THREAD_COPY_ALPHA, tr);
      break;
    case INVAR_F:
      masterBarrier(THREAD_COPY_INVAR, tr);	 
      break;
    case SCALER_F: 
      //nothing to do here      
      break;
    case LXRATE_F:
    case LXWEIGHT_F:
      masterBarrier(THREAD_COPY_LG4X_RATES, tr);
      break;
    default:
      assert(0);
    }
#endif  

  //some individual post-processing 

  switch(whichParameterType)
    {
    case RATE_F:              
    case ALPHA_F:       
    case INVAR_F:     
      break;
    case SCALER_F:       
      evaluateGenericInitrav(tr, tr->start);
      break;
    case LXRATE_F:
    case LXWEIGHT_F: 
    case FREQ_F:  
      break;
    default:
      assert(0);
    }
  
 
  assert(pos == numberOfModels);

  rax_free(startLH);
  rax_free(endLH);
  rax_free(_a);
  rax_free(_b);
  rax_free(_c);
  rax_free(_fa);
  rax_free(_fb);
  rax_free(_fc);
  rax_free(_param);
  rax_free(_x);  
  rax_free(startValues);
  rax_free(startRates);
  rax_free(startWeights);
  rax_free(startExponents);

#ifdef _DEBUG_MODEL_OPTIMIZATION
  evaluateGenericInitrav(tr, tr->start);

  if(tr->likelihood < initialLH)
    printf("%f %f\n", tr->likelihood, initialLH);
  assert(tr->likelihood >= initialLH);
#endif

}


static void optFreqs(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{ 
  int 
    rateNumber;

  double
    freqMin = -1000000.0,
    freqMax = 200.0;
  
  for(rateNumber = 0; rateNumber < states; rateNumber++)   
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, freqMin, freqMax, FREQ_F);   
}

static void optBaseFreqs(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    states,
    dnaPartitions = 0,
    aaPartitions  = 0,
    binPartitions = 0,
    multPartitions = 0;

  /* first do DNA */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:	
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;	 
	  if(tr->partitionData[ll->ld[i].partitionList[0]].optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = TRUE;
	      dnaPartitions++;  	    
	    }
	  else
	     ll->ld[i].valid = FALSE;
	  break;       
	case BINARY_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  if(dnaPartitions > 0)
    optFreqs(tr, modelEpsilon, ll, dnaPartitions, states);
  
  /* then AA */

  
  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case AA_DATA:	  
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states; 	      
	  if(tr->partitionData[ll->ld[i].partitionList[0]].optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = TRUE;
	      aaPartitions++;		
	    }
	  else
	    ll->ld[i].valid = FALSE; 
	  break;
	case DNA_DATA:
	case BINARY_DATA:      
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:	    
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}	 
    }

  if(aaPartitions > 0)      
    optFreqs(tr, modelEpsilon, ll, aaPartitions, states);

  /* then binary */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case BINARY_DATA:	  
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states; 	      
	  if(tr->partitionData[ll->ld[i].partitionList[0]].optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = TRUE;
	      binPartitions++;		
	    }
	  else
	    ll->ld[i].valid = FALSE; 
	  break;
	case DNA_DATA:	  
	case AA_DATA:      
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:	    
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}	 
    }

  if(binPartitions > 0)      
    optFreqs(tr, modelEpsilon, ll, binPartitions, states);

  /* then multi */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case GENERIC_32:	  
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states; 	      
	  if(tr->partitionData[ll->ld[i].partitionList[0]].optimizeBaseFrequencies)
	    {
	      ll->ld[i].valid = TRUE;
	      multPartitions++;		
	    }
	  else
	    ll->ld[i].valid = FALSE; 
	  break;
	case DNA_DATA:
	case AA_DATA:
	case BINARY_DATA:      
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:      
	case GENERIC_64:	    
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}	 
    }

  if(multPartitions > 0)      
    optFreqs(tr, modelEpsilon, ll, multPartitions, states);

  /* done */


  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}



static void optRates(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1;

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, RATE_MIN, RATE_MAX, RATE_F);
}



static boolean AAisGTR(tree *tr)
{
  int 
    i, 
    count = 0;

  for(i = 0; i < tr->NumberOfModels; i++)   
    {
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  count++;
	  if(tr->partitionData[i].protModels != GTR)
	    return FALSE;
	}
    }

  if(count == 0)
    return FALSE;

  return TRUE;
}

static boolean AAisUnlinkedGTR(tree *tr)
{
  int 
    i, 
    count = 0;

  for(i = 0; i < tr->NumberOfModels; i++)   
    {
      if(tr->partitionData[i].dataType == AA_DATA)
	{
	  count++;
	  if(tr->partitionData[i].protModels != GTR_UNLINKED)
	    return FALSE;
	}
    }

  if(count == 0)
    return FALSE;

  return TRUE;
}

static void optRatesGeneric(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    dnaPartitions = 0,
    aaPartitionsLinked  = 0,
    aaPartitionsUnlinked = 0,
    secondaryPartitions = 0,
    secondaryModel = -1,
    states = -1;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* first do DNA */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:	
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  ll->ld[i].valid = TRUE;
	  dnaPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }   

  if(dnaPartitions > 0)
    optRates(tr, modelEpsilon, ll, dnaPartitions, states);
  

  /* then SECONDARY */

   for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	  /* we only have one type of secondary data models in one analysis */
	case SECONDARY_DATA_6:
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA_6;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA_7: 
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA_7;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case SECONDARY_DATA:
	  states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	  secondaryModel = SECONDARY_DATA;
	  ll->ld[i].valid = TRUE;
	  secondaryPartitions++;  
	  break;
	case BINARY_DATA:
	case AA_DATA:	
	case DNA_DATA:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = FALSE;
	  break;
	default:
	  assert(0);
	}      
    }

  
   
   if(secondaryPartitions > 0)
     {
       assert(secondaryPartitions == 1);

       switch(secondaryModel)
	 {
	 case SECONDARY_DATA:
	   optRates(tr, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_6:
	   optRates(tr, modelEpsilon, ll, secondaryPartitions, states);
	   break;
	 case SECONDARY_DATA_7:
	   optRates(tr, modelEpsilon, ll, secondaryPartitions, states);
	   break; 
	 default:
	   assert(0);
	 }
     }

  /* then AA */

  if(AAisGTR(tr))
    {
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case AA_DATA:
	      states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	      ll->ld[i].valid = TRUE;
	      aaPartitionsLinked++;
	      break;
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	      ll->ld[i].valid = FALSE;
	      break;
	    default:
	      assert(0);
	    }	 
	}

      assert(aaPartitionsLinked == 1);     
      
      optRates(tr, modelEpsilon, ll, aaPartitionsLinked, states);
    }
  
   if(AAisUnlinkedGTR(tr))
    {
      aaPartitionsUnlinked = 0;

      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case AA_DATA:
	      states = tr->partitionData[ll->ld[i].partitionList[0]].states;
	      ll->ld[i].valid = TRUE;
	      aaPartitionsUnlinked++;
	      break;
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	      ll->ld[i].valid = FALSE;
	      break;
	    default:
	      assert(0);
	    }	 
	}

      assert(aaPartitionsUnlinked >= 1);     
      
      optRates(tr, modelEpsilon, ll, aaPartitionsUnlinked, states);
    }
  /* then multi-state */

  /* 
     now here we have to be careful, because every multi-state partition can actually 
     have a distinct number of states, so we will go to every multi-state partition separately,
     parallel efficiency for this will suck, but what the hell .....
  */

  if(tr->multiStateModel == GTR_MULTI_STATE)
    {     
      for(i = 0; i < ll->entries; i++)
	{
	  switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	    {
	    case GENERIC_32:
	      {
		int k;
		
		states = tr->partitionData[ll->ld[i].partitionList[0]].states;			      

		ll->ld[i].valid = TRUE;
		
		for(k = 0; k < ll->entries; k++)
		  if(k != i)
		    ll->ld[k].valid = FALSE;
		
		optRates(tr, modelEpsilon, ll, 1, states);
	      }
	      break;
	    case AA_DATA:	    
	    case DNA_DATA:	    
	    case BINARY_DATA:
	    case SECONDARY_DATA:	
	    case SECONDARY_DATA_6:
	    case SECONDARY_DATA_7:
	    case GENERIC_64:
	      break;
	    default:
	      assert(0);
	    }	 
	}           
    }

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}



/*** WARNING: the re-scaling of the branch lengths will only work if 
     it is done after the rate optimization that modifies fracchange 
***/


static void optLG4X(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels)
{
  int 
    i;

  double
    lg4xScaler,
    *lg4xScalers = (double *)rax_calloc(tr->NumberOfModels, sizeof(double)),
    *modelWeights = (double *)rax_calloc(tr->NumberOfModels, sizeof(double)),
    wgtsum = 0.0;

  for(i = 0; i < 4; i++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, i, LG4X_RATE_MIN, LG4X_RATE_MAX, LXRATE_F);
    
  for(i = 0; i < tr->NumberOfModels; i++)
    lg4xScalers[i] = 1.0;

  for(i = 0; i < ll->entries; i++)
    {
      if(ll->ld[i].valid)
	{
	  int
	    j,
	    index = ll->ld[i].partitionList[0];
	  
	  double
	    averageRate = 0.0;
	  
	  assert(ll->ld[i].partitions == 1);
	  
	  for(j = 0; j < 4; j++)
	    averageRate += tr->partitionData[index].gammaRates[j];	  
	  
	  averageRate /= 4.0;
	  
	  lg4xScalers[index] = averageRate;
	}
    }

  if(tr->NumberOfModels > 1)
    {
      for(i = 0; i < tr->NumberOfModels; i++)
	tr->fracchanges[i] = tr->rawFracchanges[i] * (1.0 / lg4xScalers[i]);
    }

  for(i = 0; i < tr->cdta->endsite; i++)
    {
      modelWeights[tr->model[i]]  += (double)tr->cdta->aliaswgt[i];
      wgtsum                      += (double)tr->cdta->aliaswgt[i];
    }

  lg4xScaler = 0.0;

  for(i = 0; i < tr->NumberOfModels; i++)
    {
      double 
	fraction = modelWeights[i] / wgtsum; 
      
      lg4xScaler += (fraction * lg4xScalers[i]); 
    }

  tr->fracchange = tr->rawFracchange * (1.0 / lg4xScaler);

  rax_free(lg4xScalers);
  rax_free(modelWeights);
}

static void optAlphasGeneric(tree *tr, double modelEpsilon, linkageList *ll)
{
  int 
    i,
    non_LG4X_Partitions = 0,
    LG4X_Partitions  = 0;

  /* assumes homogeneous super-partitions, that either contain DNA or AA partitions !*/
  /* does not check whether AA are all linked */

  /* first do non-LG4X partitions */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:			  	
	case BINARY_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  ll->ld[i].valid = TRUE;
	  non_LG4X_Partitions++;
	  break;
	case AA_DATA:	  
	  if(tr->partitionData[ll->ld[i].partitionList[0]].protModels == LG4X)
	    {
	      LG4X_Partitions++;	      
	      ll->ld[i].valid = FALSE;
	    }
	  else
	    {
	      ll->ld[i].valid = TRUE;
	      non_LG4X_Partitions++;
	    }
	  break;
	default:
	  assert(0);
	}      
    }   

 

  if(non_LG4X_Partitions > 0)
    optParamGeneric(tr, modelEpsilon, ll, non_LG4X_Partitions, -1, ALPHA_MIN, ALPHA_MAX, ALPHA_F);
  
 

  /* then LG4x partitions */

  for(i = 0; i < ll->entries; i++)
    {
      switch(tr->partitionData[ll->ld[i].partitionList[0]].dataType)
	{
	case DNA_DATA:			  	
	case BINARY_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case GENERIC_32:
	case GENERIC_64:
	  {
	    //check that alpha values do not get too large and warn users that 
	    //they should maybe use a model that is not correcting for rate heterogeneity 
	    
	    int 
	      k;

	    double const 
	      alphaMaybeNoRateHet = 10.0;
	    	  	  	  
	    for(k = 0; k < ll->ld[i].partitions; k++)
	      {
		int index = ll->ld[i].partitionList[k];

		if(tr->partitionData[index].alpha >= alphaMaybeNoRateHet)
		  {
		    printBothOpen("\nWARNING the alpha parameter with a value of %f estimated by RAxML for partition number %d with the name \"%s\"\n",  
				  tr->partitionData[index].alpha, 
				  index, 
				  tr->partitionData[index].partitionName);
		    printBothOpen("is larger than %f. You should do a model test and confirm that you actually need to incorporate a model of rate heterogeneity!\n",  alphaMaybeNoRateHet);
		    printBothOpen("You can run inferences with a plain substitution model (without rate heterogeneity) by specifyng the CAT model and the \"-V\" option!\n\n");		    
		  }
	      }

	    ll->ld[i].valid = FALSE;
	  }
	  break;
	case AA_DATA:	  
	  if(tr->partitionData[ll->ld[i].partitionList[0]].protModels == LG4X)	      
	    ll->ld[i].valid = TRUE;	   
	  else	    
	    ll->ld[i].valid = FALSE;	   	    
	  break;
	default:
	  assert(0);
	}      
    }   
  
  

  if(LG4X_Partitions > 0)
    optLG4X(tr, modelEpsilon, ll, LG4X_Partitions);


  

  for(i = 0; i < ll->entries; i++)
    ll->ld[i].valid = TRUE;
}




/*********************FUNCTIONS FOR APPROXIMATE MODEL OPTIMIZATION ***************************************/






static int catCompare(const void *p1, const void *p2)
{
 rateCategorize *rc1 = (rateCategorize *)p1;
 rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}



static void categorizePartition(tree *tr, rateCategorize *rc, int model, int lower, int upper)
{
  int
    zeroCounter,
    i, 
    k;
  
  double 
    diff, 
    min;

  for (i = lower, zeroCounter = 0; i < upper; i++, zeroCounter++) 
      {
	double
	  temp = tr->cdta->patrat[i];

	int
	  found = 0;
	
	for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
	  {
	    if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	      {
		found = 1;
		tr->cdta->rateCategory[i] = k;				
		break;
	      }
	  }
	
	if(!found)
	  {
	    min = fabs(temp - rc[0].rate);
	    tr->cdta->rateCategory[i] = 0;

	    for(k = 1; k < tr->partitionData[model].numberOfCategories; k++)
	      {
		diff = fabs(temp - rc[k].rate);

		if(diff < min)
		  {
		    min = diff;
		    tr->cdta->rateCategory[i] = k;
		  }
	      }
	  }
      }

  for(k = 0; k < tr->partitionData[model].numberOfCategories; k++)
    tr->partitionData[model].unscaled_perSiteRates[k] = rc[k].rate; 
}


#ifdef _USE_PTHREADS

void optRateCatPthreads(tree *tr, double lower_spacing, double upper_spacing, double *lhs, int n, int tid)
{
  int 
    model;
     
  size_t
    localIndex,
    i;

  for(model = 0; model < tr->NumberOfModels; model++)
    {               
      for(i = tr->partitionData[model].lower, localIndex = 0;  i < tr->partitionData[model].upper; i++)
	{
	  if(i % ((size_t)n) == ((size_t)tid))
	    {
	      
	      double initialRate, initialLikelihood, 
		leftLH, rightLH, leftRate, rightRate, v;
	      const double epsilon = 0.00001;
	      int k;	      
	      
	     

	      tr->cdta->patrat[i] = tr->cdta->patratStored[i];     
	      initialRate = tr->cdta->patrat[i];
	      
	      

	      initialLikelihood = evaluatePartialGeneric(tr, localIndex, initialRate, model); /* i is real i ??? */
	     
	      
	      leftLH = rightLH = initialLikelihood;
	      leftRate = rightRate = initialRate;
	      
	      k = 1;
	      
	      while((initialRate - k * lower_spacing > 0.0001) && 
		    ((v = evaluatePartialGeneric(tr, localIndex, initialRate - k * lower_spacing, model)) 
		     > leftLH) && 
		    (fabs(leftLH - v) > epsilon))  
		{	  
#ifndef WIN32
		  if(isnan(v))
		    assert(0);
#endif
		  
		  leftLH = v;
		  leftRate = initialRate - k * lower_spacing;
		  k++;	  
		}      
	      
	      k = 1;
	      
	      while(((v = evaluatePartialGeneric(tr, localIndex, initialRate + k * upper_spacing, model)) > rightLH) &&
		    (fabs(rightLH - v) > epsilon))    	
		{
#ifndef WIN32
		  if(isnan(v))
		    assert(0);
#endif     
		  rightLH = v;
		  rightRate = initialRate + k * upper_spacing;	 
		  k++;
		}           
	      
	      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
		{
		  if(rightLH > leftLH)	    
		    {	     
		      tr->cdta->patrat[i] = rightRate;
		      lhs[i] = rightLH;
		    }
		  else
		    {	      
		      tr->cdta->patrat[i] = leftRate;
		      lhs[i] = leftLH;
		    }
		}
	      else
		lhs[i] = initialLikelihood;
	      
	     

	      tr->cdta->patratStored[i] = tr->cdta->patrat[i];
	      localIndex++;
	    }
	}
      assert(localIndex == tr->partitionData[model].width);    
    }
}



#else


static void optRateCatModel(tree *tr, int model, double lower_spacing, double upper_spacing, double *lhs)
{
  int lower = tr->partitionData[model].lower;
  int upper = tr->partitionData[model].upper;
  int i;
  for(i = lower; i < upper; i++)
    {
      double initialRate, initialLikelihood, 
	leftLH, rightLH, leftRate, rightRate, v;
      const double epsilon = 0.00001;
      int k;
      
      tr->cdta->patrat[i] = tr->cdta->patratStored[i];     
      initialRate = tr->cdta->patrat[i];
      
      initialLikelihood = evaluatePartialGeneric(tr, i, initialRate, model); 
      
      
      leftLH = rightLH = initialLikelihood;
      leftRate = rightRate = initialRate;
      
      k = 1;
      
      while((initialRate - k * lower_spacing > 0.0001) && 
	    ((v = evaluatePartialGeneric(tr, i, initialRate - k * lower_spacing, model)) 
	     > leftLH) && 
	    (fabs(leftLH - v) > epsilon))  
	{	  
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif
	  
	  leftLH = v;
	  leftRate = initialRate - k * lower_spacing;
	  k++;	  
	}      
      
      k = 1;
      
      while(((v = evaluatePartialGeneric(tr, i, initialRate + k * upper_spacing, model)) > rightLH) &&
	    (fabs(rightLH - v) > epsilon))    	
	{
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif     
	  rightLH = v;
	  rightRate = initialRate + k * upper_spacing;	 
	  k++;
	}           
  
      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
	{
	  if(rightLH > leftLH)	    
	    {	     
	      tr->cdta->patrat[i] = rightRate;
	      lhs[i] = rightLH;
	    }
	  else
	    {	      
	      tr->cdta->patrat[i] = leftRate;
	      lhs[i] = leftLH;
	    }
	}
      else
	lhs[i] = initialLikelihood;
      
      tr->cdta->patratStored[i] = tr->cdta->patrat[i];
    }

}


#endif



/* 
   set scaleRates to FALSE everywhere such that 
   per-site rates are not scaled to obtain an overall mean rate 
   of 1.0
*/

void updatePerSiteRates(tree *tr, boolean scaleRates)
{
  int 
    i,
    model;

  if(tr->multiBranch)
    {            
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 	       
	    lower = tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper;
	  
	  if(scaleRates)
	    {
	      double 
		scaler = 0.0,       
		accRat = 0.0; 

	      int 
		accWgt     = 0;
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].unscaled_perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	   
	  
	      accRat /= ((double)accWgt);
	  
	      scaler = 1.0 / ((double)accRat);
	  	  
	      for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		tr->partitionData[model].perSiteRates[i] = scaler * tr->partitionData[model].unscaled_perSiteRates[i];	    

	      accRat = 0.0;	 
	      
	      for(i = lower; i < upper; i++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}	         

	      accRat /= ((double)accWgt);	  

	      if(!(ABS(1.0 - accRat) < 1.0E-5))	   
		printBothOpen("An assertion will fail: accumulated rate categories: %1.40f\n", accRat);

	      assert(ABS(1.0 - accRat) < 1.0E-5);
	    }                

	          	  
	  

	}
    }
  else
    {
      int
	accWgt = 0;

      double 
	scaler = 0.0,       
	accRat = 0.0; 

      if(scaleRates)
	{
	  for(model = 0, accRat = 0.0, accWgt = 0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].unscaled_perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);
		  
		  accWgt += w;
		  
		  accRat += (w * rate);
		}	      
	    }
	  
	 

	  accRat /= ((double)accWgt);
	  
	  scaler = 1.0 / ((double)accRat);
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      for(i = 0; i < tr->partitionData[model].numberOfCategories; i++)
		tr->partitionData[model].perSiteRates[i] = scaler * tr->partitionData[model].unscaled_perSiteRates[i];
	    }

	  for(model = 0, accRat = 0.0; model < tr->NumberOfModels; model++)
	    {
	      int 
		localCount = 0,
		lower = tr->partitionData[model].lower,
		upper = tr->partitionData[model].upper;
	      
	      for(i = lower, localCount = 0; i < upper; i++, localCount++)
		{
		  int 
		    w = tr->cdta->aliaswgt[i];
		  
		  double
		    rate = tr->partitionData[model].perSiteRates[tr->cdta->rateCategory[i]];
		  
		  assert(0 <= tr->cdta->rateCategory[i] && tr->cdta->rateCategory[i] < tr->maxCategories);	      
		  
		  accRat += (w * rate);
		}
	    }           

	  accRat /= ((double)accWgt);	  

	  if(!(ABS(1.0 - accRat) < 1.0E-5))	   
	    printBothOpen("An assertion will fail: accumulated rate categories: %1.40f\n", accRat);
	    
	  assert(ABS(1.0 - accRat) < 1.0E-5);
	}
         
       

    }
  
      
#ifdef _USE_PTHREADS
      masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif               
}

static void optimizeRateCategories(tree *tr, int _maxCategories)
{
  assert(_maxCategories > 0);

  if(_maxCategories > 1)
    {
      double  
	temp,  
	lower_spacing, 
	upper_spacing,
	initialLH = tr->likelihood,	
	*ratStored = (double *)rax_malloc(sizeof(double) * tr->cdta->endsite),
	*lhs =       (double *)rax_malloc(sizeof(double) * tr->cdta->endsite),
	**oldCategorizedRates = (double **)rax_malloc(sizeof(double *) * tr->NumberOfModels),
	**oldUnscaledCategorizedRates = (double **)rax_malloc(sizeof(double *) * tr->NumberOfModels);

      int  
	i,
	k,
	maxCategories = _maxCategories,
	*oldCategory =  (int *)rax_malloc(sizeof(int) * tr->cdta->endsite),
	model,
	*oldNumbers = (int *)rax_malloc(sizeof(int) * tr->NumberOfModels);
  
      assert(isTip(tr->start->number, tr->rdta->numsp));         
      
      determineFullTraversal(tr->start, tr);

      if(optimizeRateCategoryInvocations == 1)
	{
	  lower_spacing = 0.5 / ((double)optimizeRateCategoryInvocations);
	  upper_spacing = 1.0 / ((double)optimizeRateCategoryInvocations);
	}
      else
	{
	  lower_spacing = 0.05 / ((double)optimizeRateCategoryInvocations);
	  upper_spacing = 0.1 / ((double)optimizeRateCategoryInvocations);
	}
      
      if(lower_spacing < 0.001)
	lower_spacing = 0.001;
      
      if(upper_spacing < 0.001)
	upper_spacing = 0.001;
      
      optimizeRateCategoryInvocations++;

      memcpy(oldCategory, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);	     
      memcpy(ratStored,   tr->cdta->patratStored, sizeof(double) * tr->cdta->endsite);

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldNumbers[model]          = tr->partitionData[model].numberOfCategories;

	  oldCategorizedRates[model] = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
	  oldUnscaledCategorizedRates[model] = (double *)rax_malloc(sizeof(double) * tr->maxCategories);
	  
	  memcpy(oldCategorizedRates[model], tr->partitionData[model].perSiteRates, tr->maxCategories * sizeof(double));	  	 	  
	  memcpy(oldUnscaledCategorizedRates[model], tr->partitionData[model].unscaled_perSiteRates, tr->maxCategories * sizeof(double));
	}      
      
#ifdef _USE_PTHREADS
      tr->lhs = lhs;
      tr->lower_spacing = lower_spacing;
      tr->upper_spacing = upper_spacing;
      masterBarrier(THREAD_RATE_CATS, tr);
#else
      for(model = 0; model < tr->NumberOfModels; model++)      
	optRateCatModel(tr, model, lower_spacing, upper_spacing, lhs);
#endif     

      for(model = 0; model < tr->NumberOfModels; model++)
	{     
	  int 
	    where = 1,
	    found = 0,
	    width = tr->partitionData[model].upper -  tr->partitionData[model].lower,
	    upper = tr->partitionData[model].upper,
	    lower = tr->partitionData[model].lower;
	    
	  rateCategorize 
	    *rc = (rateCategorize *)rax_malloc(sizeof(rateCategorize) * width);		 
	
	  for (i = 0; i < width; i++)
	    {
	      rc[i].accumulatedSiteLikelihood = 0.0;
	      rc[i].rate = 0.0;
	    }  
	
	  rc[0].accumulatedSiteLikelihood = lhs[lower];
	  rc[0].rate = tr->cdta->patrat[lower];
	
	  tr->cdta->rateCategory[lower] = 0;
	
	  for (i = lower + 1; i < upper; i++) 
	    {
	      temp = tr->cdta->patrat[i];
	      found = 0;
	    
	      for(k = 0; k < where; k++)
		{
		  if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
		    {
		      found = 1;						
		      rc[k].accumulatedSiteLikelihood += lhs[i];	
		      break;
		    }
		}
	    
	      if(!found)
		{	    
		  rc[where].rate = temp;	    
		  rc[where].accumulatedSiteLikelihood += lhs[i];	    
		  where++;
		}
	    }
	
	  qsort(rc, where, sizeof(rateCategorize), catCompare);
	
	  if(where < maxCategories)
	    {
	      tr->partitionData[model].numberOfCategories = where;
	      categorizePartition(tr, rc, model, lower, upper);
	    }
	  else
	    {
	      tr->partitionData[model].numberOfCategories = maxCategories;	
	      categorizePartition(tr, rc, model, lower, upper);
	    }
	
	  rax_free(rc);
	}
        	
      updatePerSiteRates(tr, TRUE);	

      evaluateGenericInitrav(tr, tr->start);
      
      if(tr->likelihood < initialLH)
	{	 		  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      tr->partitionData[model].numberOfCategories = oldNumbers[model];
	      memcpy(tr->partitionData[model].perSiteRates, oldCategorizedRates[model], tr->maxCategories * sizeof(double));
	      memcpy(tr->partitionData[model].unscaled_perSiteRates, oldUnscaledCategorizedRates[model], tr->maxCategories * sizeof(double));
	    }	      
	  
	  memcpy(tr->cdta->patratStored, ratStored, sizeof(double) * tr->cdta->endsite);
	  memcpy(tr->cdta->rateCategory, oldCategory, sizeof(int) * tr->cdta->endsite);	     
	  
	  updatePerSiteRates(tr, FALSE);
	  
	  evaluateGenericInitrav(tr, tr->start);

	  /* printf("REVERT: %1.40f %1.40f\n", initialLH, tr->likelihood); */

	  assert(initialLH == tr->likelihood);
	}
          
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  rax_free(oldCategorizedRates[model]);
	  rax_free(oldUnscaledCategorizedRates[model]);
	}
                   
      rax_free(oldCategorizedRates);
      rax_free(oldUnscaledCategorizedRates);
      rax_free(oldCategory);
      rax_free(ratStored);       
      rax_free(lhs); 
      rax_free(oldNumbers);
    }
}
  






/*****************************************************************************************************/

void resetBranches(tree *tr)
{
  nodeptr  p, q;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
      for(i = 0; i < tr->numBranches; i++)
	p->z[i] = defaultz;
	
      q = p->next;
      while(q != p)
	{	
	  for(i = 0; i < tr->numBranches; i++)
	    q->z[i] = defaultz;	    
	  q = q->next;
	}
      p++;
    }
}

static double fixZ(double z)
{
  if(z > zmax)
    return zmax;
  
  if(z < zmin)
    return zmin;
  
  return z;
}

static double getFracChange(tree *tr, int model)
{
  if(tr->NumberOfModels == 1)
    return tr->fracchange;
  else
    return tr->fracchanges[model];
}
void scaleBranches(tree *tr, boolean fromFile)
{
  nodeptr  
    p;
  
  int  
    model,
    i,
    nodes, 
    count = 0;

  double 
    z;
  
  if(!tr->storedBrLens)
    tr->storedBrLens = (double *)rax_malloc(sizeof(double) * (2 * tr->mxtips - 3) * 2);

  assert(tr->numBranches == tr->NumberOfModels);
  
  nodes = tr->mxtips  +  tr->mxtips - 2;
  p = tr->nodep[1];

  for(i = 1; i <= nodes; i++)
    {      
      p = tr->nodep[i];
      
      if(fromFile)
	{
	  tr->storedBrLens[count] = p->z[0];
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      z = exp(-p->z[0] / getFracChange(tr, model));

	      z = fixZ(z);	     

	      p->z[model] = z;
	    }
	}
      else
	{	
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      z = tr->partitionData[model].brLenScaler * tr->storedBrLens[count];
	     
	      z = exp(-z / getFracChange(tr, model));
	      
	      z = fixZ(z);

	      p->z[model] = z;
	    }
	}
      count++;
	
	
      if(i > tr->mxtips)
	{	
	  if(fromFile)
	    {
	      tr->storedBrLens[count] = p->next->z[0];
	      
	      for(model = 0; model < tr->NumberOfModels; model++)
		{
		  z = exp(-p->next->z[0] / getFracChange(tr, model));
		  
		  z = fixZ(z);

		  p->next->z[model] = z;

		}
	    }
	  else
	    {	      
	      for(model = 0; model < tr->NumberOfModels; model++)
		{
		  z = tr->partitionData[model].brLenScaler * tr->storedBrLens[count];
		  z = exp(-z / getFracChange(tr, model));		  		 

		  z = fixZ(z);

		  p->next->z[model] = z;
		}
	    }
	  count++;
	  
	  if(fromFile)
	    {
	      tr->storedBrLens[count] = p->next->next->z[0];
	      
	      for(model = 0; model < tr->NumberOfModels; model++)
		{
		  z = exp(-p->next->next->z[0] / getFracChange(tr, model));
		  
		  z = fixZ(z);		  
		  
		  p->next->next->z[model] = z;
		}
	    }
	  else	  
	    {	     
	       for(model = 0; model < tr->NumberOfModels; model++)
		{
		  z = tr->partitionData[model].brLenScaler * tr->storedBrLens[count];
		  
		  z = exp(-z / getFracChange(tr, model));		 

		  z = fixZ(z);

		  p->next->next->z[model] = z;
		}
	    }
	  count++;
	}	  
    }
  
  assert(count == (2 * tr->mxtips - 3) * 2);
}



static void printAAmatrix(tree *tr, double epsilon)
{
  

  if(AAisGTR(tr) || AAisUnlinkedGTR(tr))
    {
      int 
	model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionData[model].dataType == AA_DATA) 
	    {
	      char 
		gtrFileName[1024];
		/*epsilonStr[1024];*/
	      
	      FILE 
		*gtrFile;
	      
	      double 
		*rates = tr->partitionData[model].substRates,
		*f     = tr->partitionData[model].frequencies,
		q[20][20];
	      
	      int    
		r = 0,
		i, 
		j;

	      assert(tr->partitionData[model].protModels == GTR || tr->partitionData[model].protModels == GTR_UNLINKED);

	      /*sprintf(epsilonStr, "%f", epsilon);*/

	      strcpy(gtrFileName, workdir);
	      strcat(gtrFileName, "RAxML_proteinGTRmodel.");
	      strcat(gtrFileName, run_id);
	      
	      /*
		strcat(gtrFileName, "_");
		strcat(gtrFileName, epsilonStr);
	      */
	      
	      strcat(gtrFileName, "_Partition_");
	      strcat(gtrFileName, tr->partitionData[model].partitionName);

	      gtrFile = myfopen(gtrFileName, "wb");

	      for(i = 0; i < 20; i++)
		for(j = 0; j < 20; j++)
		  q[i][j] = 0.0;

	      for(i = 0; i < 19; i++)
		for(j = i + 1; j < 20; j++)
		  q[i][j] = rates[r++];

	      for(i = 0; i < 20; i++)
		for(j = 0; j <= i; j++)
		  {
		    if(i == j)
		      q[i][j] = 0.0;
		    else
		      q[i][j] = q[j][i];
		  }
	   
	      for(i = 0; i < 20; i++)
		{
		  for(j = 0; j < 20; j++)		
		    fprintf(gtrFile, "%1.80f ", q[i][j]);
		
		  fprintf(gtrFile, "\n");
		}
	      for(i = 0; i < 20; i++)
		fprintf(gtrFile, "%1.80f ", f[i]);
	      fprintf(gtrFile, "\n");

	      fclose(gtrFile);
	      
	      if(tr->partitionData[model].protModels == GTR)
		printBothOpen("\nPrinted linked AA GTR matrix that achieved an overall improvement of %f log likelihood units for partition %s to file %s\n\n", epsilon, tr->partitionData[model].partitionName, gtrFileName);	      	    
	      else
		printBothOpen("\nPrinted unlinked AA GTR matrix that achieved an overall improvement of %f log likelihood units for partition %s to file %s\n\n", epsilon, tr->partitionData[model].partitionName, gtrFileName);
	    }

	}	  
    }
}



static void optScaler(tree *tr, double modelEpsilon, linkageList *ll)
{  
  optParamGeneric(tr, modelEpsilon, ll, ll->entries, -1, 0.01, 100.0, SCALER_F);
}

static void autoProtein(tree *tr)
{
  int 
    countAutos = 0,   
    model;  
  
  for(model = 0; model < tr->NumberOfModels; model++)	      
    if(tr->partitionData[model].protModels == AUTO)
      countAutos++;

  if(countAutos > 0)
    {
      int 
	i,
	numProteinModels = AUTO,
	*bestIndex = (int*)rax_malloc(sizeof(int) * tr->NumberOfModels),
	*oldIndex  = (int*)rax_malloc(sizeof(int) * tr->NumberOfModels);

      double
	startLH,
	*bestScores = (double*)rax_malloc(sizeof(double) * tr->NumberOfModels);    

      topolRELL_LIST 
	*rl = (topolRELL_LIST *)rax_malloc(sizeof(topolRELL_LIST));

      initTL(rl, tr, 1);
      saveTL(rl, tr, 0);

      evaluateGenericInitrav(tr, tr->start); 

      startLH = tr->likelihood;

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldIndex[model] = tr->partitionData[model].autoProtModels;
	  bestIndex[model] = -1;
	  bestScores[model] = unlikely;
	}
      
      for(i = 0; i < numProteinModels; i++)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {	   
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = i;
		  initReversibleGTR(tr, model);  
		}
	    }

#ifdef _USE_PTHREADS	
	  masterBarrier(THREAD_COPY_RATES, tr);	   
#endif
      
	  resetBranches(tr);
	  evaluateGenericInitrav(tr, tr->start);  
	  treeEvaluate(tr, 0.5);     
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{		  
		  if(tr->perPartitionLH[model] > bestScores[model])
		    {
		      bestScores[model] = tr->perPartitionLH[model];
		      bestIndex[model] = i;		      
		    }
		}	      
	    }       
	}           
      
      printBothOpen("Automatic protein model assignment algorithm:\n\n");

      for(model = 0; model < tr->NumberOfModels; model++)
	{	   
	  if(tr->partitionData[model].protModels == AUTO)
	    {
	      tr->partitionData[model].autoProtModels = bestIndex[model];
	      initReversibleGTR(tr, model);  
	      printBothOpen("\tPartition: %d best-scoring AA model: %s likelihood %f\n", model, protModels[tr->partitionData[model].autoProtModels], bestScores[model]);
	    }	 
	}

      printBothOpen("\n\n");

#ifdef _USE_PTHREADS	
      masterBarrier(THREAD_COPY_RATES, tr);	   
#endif
          
      resetBranches(tr);
      evaluateGenericInitrav(tr, tr->start); 
      treeEvaluate(tr, 2.0);    
      
      if(tr->likelihood < startLH)
	{	
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = oldIndex[model];
		  initReversibleGTR(tr, model);
		}
	    }
	  
#ifdef _USE_PTHREADS	
	  masterBarrier(THREAD_COPY_RATES, tr);	   
#endif 
	  restoreTL(rl, tr, 0);	
	  evaluateGenericInitrav(tr, tr->start);              
	}
      
      assert(tr->likelihood >= startLH);
      
      freeTL(rl);   
      rax_free(rl); 
      
      rax_free(oldIndex);
      rax_free(bestIndex);
      rax_free(bestScores);
    }
}


//#define _DEBUG_MOD_OPT

void modOpt(tree *tr, analdef *adef, boolean resetModel, double likelihoodEpsilon)
{ 
  int i, model, catOpt = 0; 

  boolean
    changedRoot = FALSE;

  double 
    inputLikelihood,
    currentLikelihood,
    modelEpsilon = 0.0001;
  linkageList 
    *alphaList,
    *invarList,
    *rateList,
    *scalerList,
    *freqList; 
  /*
    int linkedAlpha[4] = {0, 0, 0, 0};   
    int linkedInvar[4] = {0, 0, 0, 0}; 
    int linkedRates[4] = {0, 0, 0, 0};
  */  
  int 
    *unlinked = (int *)rax_malloc(sizeof(int) * tr->NumberOfModels),
    *linked =  (int *)rax_malloc(sizeof(int) * tr->NumberOfModels);
  
  assert(!adef->useBinaryModelFile);
 
  modelEpsilon = 0.0001;
  
  
  for(i = 0; i < tr->NumberOfModels; i++)
    {
      unlinked[i] = i;
      linked[i] = 0;
    }
  
  alphaList = initLinkageList(unlinked, tr);
  invarList = initLinkageList(unlinked, tr);
  rateList  = initLinkageListGTR(tr);
  scalerList = initLinkageList(unlinked, tr);
  freqList = initLinkageList(unlinked, tr);
    
  if(!(adef->mode == CLASSIFY_ML))
    {
      if(tr->start != tr->nodep[1])
	changedRoot = TRUE;      
      tr->start = tr->nodep[1];
    }
  
  if(resetModel)
    {
      

      initRateMatrix(tr);
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{     	  
	  if(adef->useInvariant)
	    {
	      int lower, upper;
	      int count = 0;
	      int total = 0;
	      
	      lower = tr->partitionData[model].lower;
	      upper = tr->partitionData[model].upper;
	      	      
	      for(i = lower; i < upper; i++)
		{
		  if(tr->invariant[i] < 4) 		
		    count += tr->cdta->aliaswgt[i];		  		
		  total += tr->cdta->aliaswgt[i];
		}
	      
	      tr->partitionData[model].propInvariant = ((double)count)/((double) total);
	    }   
	  
	  tr->partitionData[model].alpha = 1.0;     
	  
	  initReversibleGTR(tr, model);      
	  
	  makeGammaCats(tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4, tr->useGammaMedian); 
	}
#ifdef _USE_PTHREADS     
      masterBarrier(THREAD_RESET_MODEL ,tr);    
#endif
      
      resetBranches(tr);
      
      evaluateGenericInitrav(tr, tr->start); 
            
      treeEvaluate(tr, 0.25);        

      evaluateGenericInitrav(tr, tr->start);
    }  
  
  inputLikelihood = tr->likelihood;

  evaluateGenericInitrav(tr, tr->start); 

  if(!changedRoot)
    assert(inputLikelihood == tr->likelihood);
  
  /* no need for individual models here, just an init on params equal for all partitions*/
  
  do
    {           
      currentLikelihood = tr->likelihood;      

#ifdef _DEBUG_MOD_OPT
      printf("start: %1.40f\n", currentLikelihood);
#endif

      optRatesGeneric(tr, modelEpsilon, rateList);         

      evaluateGenericInitrav(tr, tr->start); 
      
#ifdef _DEBUG_MOD_OPT
      printf("after rates %1.40f\n", tr->likelihood);
#endif

      autoProtein(tr);
      
      if(adef->mode != OPTIMIZE_BR_LEN_SCALER)
	treeEvaluate(tr, 0.0625);                     	            
      else 	
	optScaler(tr, modelEpsilon, scalerList);     

#ifdef _DEBUG_MOD_OPT
      evaluateGenericInitrav(tr, tr->start); 
      printf("after br-len 1 %1.40f\n", tr->likelihood);
#endif

      evaluateGenericInitrav(tr, tr->start);

      optBaseFreqs(tr, modelEpsilon, freqList);
      
      evaluateGenericInitrav(tr, tr->start);
      
      treeEvaluate(tr, 0.0625);

#ifdef _DEBUG_MOD_OPT
      evaluateGenericInitrav(tr, tr->start); 
      printf("after optBaseFreqs 1 %f\n", tr->likelihood);
#endif 

      switch(tr->rateHetModel)
	{	  
	case GAMMA_I:
	  optAlphasGeneric(tr, modelEpsilon, alphaList);

#ifdef _DEBUG_MOD_OPT
	  evaluateGenericInitrav(tr, tr->start); 
	  printf("after alphas %1.40f\n", tr->likelihood);
#endif

	  optInvar(tr, modelEpsilon, invarList);

#ifdef _DEBUG_MOD_OPT
	  evaluateGenericInitrav(tr, tr->start); 
	  printf("after invar %1.40f\n", tr->likelihood);
#endif

	  if(adef->mode != OPTIMIZE_BR_LEN_SCALER)		      	    	   	 
	    treeEvaluate(tr, 0.1);  
	  else 	  
	    optScaler(tr, modelEpsilon, scalerList);	
	  
#ifdef _DEBUG_MOD_OPT
	  evaluateGenericInitrav(tr, tr->start); 
	  printf("after br-len 2 %1.40f\n", tr->likelihood);
#endif

	  break;
	case GAMMA:      
	  optAlphasGeneric(tr, modelEpsilon, alphaList); 

#ifdef _DEBUG_MOD_OPT
	  evaluateGenericInitrav(tr, tr->start); 
	  printf("after alphas %1.40f\n", tr->likelihood);
#endif

	 
	  onlyInitrav(tr, tr->start); 
	  
	  evaluateGeneric(tr, tr->start); 
	  
	  if(adef->mode != OPTIMIZE_BR_LEN_SCALER)	 	 
	    treeEvaluate(tr, 0.1);
	  else 	   
	    optScaler(tr, modelEpsilon, scalerList);

#ifdef _DEBUG_MOD_OPT
	  evaluateGenericInitrav(tr, tr->start); 
	  printf("after br-len 3 %1.40f\n", tr->likelihood);
#endif

	 
	  break;	  
	case CAT:
	  if(!tr->noRateHet)
	    {
	      if(catOpt < 3)
		{	      	
		  evaluateGenericInitrav(tr, tr->start);
		  optimizeRateCategories(tr, adef->categories);	      	     	      	      
		  catOpt++;
		}
	    }
	 

#ifdef _DEBUG_MOD_OPT
	  evaluateGenericInitrav(tr, tr->start); 
	  printf("after cat-opt %f\n", tr->likelihood);
#endif	  

	  break;	  
	default:
	  assert(0);
	}       

      if(tr->likelihood < currentLikelihood)
	{
	 
	  if(fabs(tr->likelihood - currentLikelihood) > MIN(0.0000001, likelihoodEpsilon))
	    {
	      printf("%1.40f %1.40f\n", tr->likelihood, currentLikelihood);
	      assert(0);
	    }
	}
      
      printAAmatrix(tr, fabs(currentLikelihood - tr->likelihood));    
    }
  while(fabs(currentLikelihood - tr->likelihood) > likelihoodEpsilon);  
  
  rax_free(unlinked);
  rax_free(linked);
  freeLinkageList(alphaList);
  freeLinkageList(rateList);
  freeLinkageList(invarList);  
  freeLinkageList(scalerList);
  freeLinkageList(freqList);


}




/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/



static double branchLength(int model, double *z, tree *tr)
{
  double x;
  
  x = z[model];
  assert(x > 0);
  if (x < zmin) 
    x = zmin;  
  
 
  assert(x <= zmax);
  
  if(!tr->multiBranch)             
    x = -log(x) * tr->fracchange;       
  else
    x = -log(x) * tr->fracchanges[model];

  return x;

}


static double treeLengthRec(nodeptr p, tree *tr, int model)
{  
  double 
    x = branchLength(model, p->z, tr);

  if(isTip(p->number, tr->rdta->numsp))  
    return x;    
  else
    {
      double acc = 0;
      nodeptr q;                
     
      q = p->next;      

      while(q != p)
	{
	  acc += treeLengthRec(q->back, tr, model);
	  q = q->next;
	}

      return acc + x;
    }
}

double treeLength(tree *tr, int model)
{ 
  /* printf("SCALER: %f\n", tr->partitionData[model].brLenScaler); */

  return treeLengthRec(tr->start->back, tr, model);
}





