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
#ifdef _HET
#define RATE_F_HET 7
#endif

static boolean optimizeRatesBFGS(tree *tr);
static void setRateModel(tree *tr, int model, double rate, int position);

static void brentGeneric(double *ax, double *bx, double *cx, double *fb, double tol, double *xmin, double *result, int numberOfModels, 
			 int whichFunction, int rateNumber, tree *tr, linkageList *ll, double *lim_inf, double *lim_sup);

static int brakGeneric(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		       double *fc, double *lim_inf, double *lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, linkageList *ll);

static void optParamGeneric(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double lim_inf, double lim_sup, int whichParameterType);

static void updateWeights(tree *tr, int model, int rate, double value);

/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/





#ifdef _HET
static void setRateModel(tree *tr, int model, double rate, int position, boolean isHet)
#else
static void setRateModel(tree *tr, int model, double rate, int position)
#endif
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

  if(tr->partitionData[model].nonGTR || (tr->partitionData[model].dataType == DNA_DATA && (tr->useK80 || tr->useHKY85)))
    {    
      int 
	i, 
	k = tr->partitionData[model].symmetryVector[position];

      assert(tr->partitionData[model].dataType == SECONDARY_DATA ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_6 ||
	     tr->partitionData[model].dataType == SECONDARY_DATA_7 ||
	     tr->partitionData[model].dataType == DNA_DATA);

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
    {
#ifdef _HET
      if(isHet)		  
	tr->partitionData[model].substRates_TIP[position] = rate;	
      else       
	tr->partitionData[model].substRates[position] = rate;	
#else
      tr->partitionData[model].substRates[position] = rate;
#endif
    }
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
#ifdef _HET
    case RATE_F:
      setRateModel(tr, index, value, rateNumber, FALSE);  
      initReversibleGTR(tr, index);		 
      break;
    case RATE_F_HET:
      setRateModel(tr, index, value, rateNumber, TRUE);  
      initReversibleGTR(tr, index);		 
      break; 
#else
    case RATE_F:
      setRateModel(tr, index, value, rateNumber);  
      initReversibleGTR(tr, index);		 
      break;
#endif
    case ALPHA_F:
      tr->partitionData[index].alpha = value;
      makeGammaCats(tr->rateHetModel, tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4, tr->useGammaMedian, tr->partitionData[index].propInvariant);
      break;
    case INVAR_F:
      tr->partitionData[index].propInvariant = value;
      makeGammaCats(tr->rateHetModel, tr->partitionData[index].alpha, tr->partitionData[index].gammaRates, 4, tr->useGammaMedian, tr->partitionData[index].propInvariant);
      break;
    case SCALER_F:     
      tr->partitionData[index].brLenScaler = value;
      scaleBranches(tr, FALSE);			        
      break;
    case LXRATE_F:
      tr->partitionData[index].gammaRates[rateNumber] = value;         
      scaleLG4X_EIGN(tr, index);    
      break;
    case LXWEIGHT_F:
      updateWeights(tr, index, rateNumber, value);      
      scaleLG4X_EIGN(tr, index);
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

void scaleLG4X_EIGN(tree *tr, int model)
{
  double 
    acc = 0.0;

  int 
    i, 
    l;
          
  for(i = 0; i < 4; i++)	     
    acc += tr->partitionData[model].weights[i] *  tr->partitionData[model].gammaRates[i];

  acc = 1.0 / acc;

  /*
    printf("update %f %f %f %f %f\n", acc, tr->partitionData[model].gammaRates[0], tr->partitionData[model].gammaRates[1], tr->partitionData[model].gammaRates[2], 
    tr->partitionData[model].gammaRates[3]);

    printf("weigths: %f %f %f %f\n", tr->partitionData[model].weights[0], tr->partitionData[model].weights[1], tr->partitionData[model].weights[2], 
    tr->partitionData[model].weights[3]);
  */

  for(i = 0; i < 4; i++)
    for(l = 0; l < 19; l++)
	tr->partitionData[model].EIGN_LG4[i][l] = tr->partitionData[model].rawEIGN_LG4[i][l] * acc;
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_LG4X_EIGN, tr);
#endif

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

  evaluateGenericInitrav(tr, tr->start);
 
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
static void evaluateChange(tree *tr, int rateNumber, double *value, double *result, boolean* converged, int whichFunction, int numberOfModels, linkageList *ll)
{ 
  int 
    i, 
    k, 
    pos;  
   
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
#ifdef _HET
    case  RATE_F_HET:
#endif
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
#ifdef _HET
    case RATE_F_HET:
      assert(0); 
      // not implemented for Pthreads!
      //needs an own barrier
      break;
#endif
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
#ifdef _HET
    case RATE_F_HET:
#endif
    case ALPHA_F:
    case SCALER_F:
    case LXRATE_F: 
    case FREQ_F: 
    case LXWEIGHT_F:
    case INVAR_F:
      evaluateGenericInitrav(tr, tr->start);           
      break;
    default:
      assert(0);
    }
#endif

 

  
  
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
			 int whichFunction, int rateNumber, tree *tr, linkageList *ll, double *lim_inf, double *lim_sup)
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
      assert(a[i] >= lim_inf[i] && a[i] <= lim_sup[i]);
      assert(b[i] >= lim_inf[i] && b[i] <= lim_sup[i]);
      assert(x[i] >= lim_inf[i] && x[i] <= lim_sup[i]);
      assert(v[i] >= lim_inf[i] && v[i] <= lim_sup[i]);
      assert(w[i] >= lim_inf[i] && w[i] <= lim_sup[i]);
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
	      assert(a[i] >= lim_inf[i] && a[i] <= lim_sup[i]);
	      assert(b[i] >= lim_inf[i] && b[i] <= lim_sup[i]);
	      assert(x[i] >= lim_inf[i] && x[i] <= lim_sup[i]);
	      assert(v[i] >= lim_inf[i] && v[i] <= lim_sup[i]);
	      assert(w[i] >= lim_inf[i] && w[i] <= lim_sup[i]);
  
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
		assert(u[i] >= lim_inf[i] && u[i] <= lim_sup[i]);
	    }
	}
                 
      evaluateChange(tr, rateNumber, u, fu, converged, whichFunction, numberOfModels, ll);

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
	      
	      assert(a[i] >= lim_inf[i] && a[i] <= lim_sup[i]);
	      assert(b[i] >= lim_inf[i] && b[i] <= lim_sup[i]);
	      assert(x[i] >= lim_inf[i] && x[i] <= lim_sup[i]);
	      assert(v[i] >= lim_inf[i] && v[i] <= lim_sup[i]);
	      assert(w[i] >= lim_inf[i] && w[i] <= lim_sup[i]);
	      assert(u[i] >= lim_inf[i] && u[i] <= lim_sup[i]);
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
		       double *fc, double *lim_inf, double *lim_sup, 
		       int numberOfModels, int rateNumber, int whichFunction, tree *tr, linkageList *ll)
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

      if(param[i] > lim_sup[i]) 	
	param[i] = ax[i] = lim_sup[i];
      
      if(param[i] < lim_inf[i]) 
	param[i] = ax[i] = lim_inf[i];

      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
    }
   
  
  evaluateChange(tr, rateNumber, param, fa, converged, whichFunction, numberOfModels, ll);


  for(i = 0; i < numberOfModels; i++)
    {
      param[i] = bx[i];
      if(param[i] > lim_sup[i]) 
	param[i] = bx[i] = lim_sup[i];
      if(param[i] < lim_inf[i]) 
	param[i] = bx[i] = lim_inf[i];

      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
    }
  
  evaluateChange(tr, rateNumber, param, fb, converged, whichFunction, numberOfModels, ll);

  for(i = 0; i < numberOfModels; i++)  
    {
      if (fb[i] > fa[i]) 
	{	  
	  SHFT(dum[i],ax[i],bx[i],dum[i]);
	  SHFT(dum[i],fa[i],fb[i],dum[i]);
	}
      
      cx[i] = bx[i] + MNBRAK_GOLD * (bx[i] - ax[i]);
      
      param[i] = cx[i];
      
      if(param[i] > lim_sup[i]) 
	param[i] = cx[i] = lim_sup[i];
      if(param[i] < lim_inf[i]) 
	param[i] = cx[i] = lim_inf[i];

      assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
    }
  
 
  evaluateChange(tr, rateNumber, param, fc, converged, whichFunction, numberOfModels,  ll);

   while(1) 
     {       
       allConverged = TRUE;

       for(i = 0; i < numberOfModels && allConverged; i++)
	 allConverged = allConverged && converged[i];

       if(allConverged)
	 {
	   for(i = 0; i < numberOfModels; i++)
	     {	       
	       if(ax[i] > lim_sup[i]) 
		 ax[i] = lim_sup[i];
	       if(ax[i] < lim_inf[i]) 
		 ax[i] = lim_inf[i];

	       if(bx[i] > lim_sup[i]) 
		 bx[i] = lim_sup[i];
	       if(bx[i] < lim_inf[i]) 
		 bx[i] = lim_inf[i];
	       
	       if(cx[i] > lim_sup[i]) 
		 cx[i] = lim_sup[i];
	       if(cx[i] < lim_inf[i]) 
		 cx[i] = lim_inf[i];
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
		   
		       if(ax[i] > lim_sup[i]) 
			 ax[i] = lim_sup[i];
		       if(ax[i] < lim_inf[i]) 
			 ax[i] = lim_inf[i];
		       if(bx[i] > lim_sup[i]) 
			 bx[i] = lim_sup[i];
		       if(bx[i] < lim_inf[i]) 
			 bx[i] = lim_inf[i];
		       if(cx[i] > lim_sup[i]) 
			 cx[i] = lim_sup[i];
		       if(cx[i] < lim_inf[i]) 
			 cx[i] = lim_inf[i];
		       
		       r[i]=(bx[i]-ax[i])*(fb[i]-fc[i]);
		       q[i]=(bx[i]-cx[i])*(fb[i]-fa[i]);
		       u[i]=(bx[i])-((bx[i]-cx[i])*q[i]-(bx[i]-ax[i])*r[i])/
			 (2.0*SIGN(MAX(fabs(q[i]-r[i]),MNBRAK_TINY),q[i]-r[i]));
		       
		       ulim[i]=(bx[i])+MNBRAK_GLIMIT*(cx[i]-bx[i]);
		       
		       if(u[i] > lim_sup[i]) 
			 u[i] = lim_sup[i];
		       if(u[i] < lim_inf[i]) 
			 u[i] = lim_inf[i];
		       if(ulim[i] > lim_sup[i]) 
			 ulim[i] = lim_sup[i];
		       if(ulim[i] < lim_inf[i]) 
			 ulim[i] = lim_inf[i];
		       
		       if ((bx[i]-u[i])*(u[i]-cx[i]) > 0.0)
			 {
			   param[i] = u[i];
			   if(param[i] > lim_sup[i]) 			     
			     param[i] = u[i] = lim_sup[i];
			   if(param[i] < lim_inf[i])
			     param[i] = u[i] = lim_inf[i];
			   endState[i] = 1;
			 }
		       else 
			 {
			   if ((cx[i]-u[i])*(u[i]-ulim[i]) > 0.0) 
			     {
			       param[i] = u[i];
			       if(param[i] > lim_sup[i]) 
				 param[i] = u[i] = lim_sup[i];
			       if(param[i] < lim_inf[i]) 
				 param[i] = u[i] = lim_inf[i];
			       endState[i] = 2;
			     }		  	       
			   else
			     {
			       if ((u[i]-ulim[i])*(ulim[i]-cx[i]) >= 0.0) 
				 {
				   u[i] = ulim[i];
				   param[i] = u[i];	
				   if(param[i] > lim_sup[i]) 
				     param[i] = u[i] = ulim[i] = lim_sup[i];
				   if(param[i] < lim_inf[i]) 
				     param[i] = u[i] = ulim[i] = lim_inf[i];
				   endState[i] = 0;
				 }		  		
			       else 
				 {		  
				   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
				   param[i] = u[i];
				   endState[i] = 0;
				   if(param[i] > lim_sup[i]) 
				     param[i] = u[i] = lim_sup[i];
				   if(param[i] < lim_inf[i]) 
				     param[i] = u[i] = lim_inf[i];
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
	       assert(param[i] >= lim_inf[i] && param[i] <= lim_sup[i]);
	     }
	 }
             
       evaluateChange(tr, rateNumber, param, temp, converged, whichFunction, numberOfModels, ll);

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
			   assert(u[i] >= lim_inf[i] && u[i] <= lim_sup[i]);
			   cx[i]=u[i];
			   fc[i]=fu[i];
			   converged[i] = TRUE;			  
			 }
		       else
			 {		   
			   u[i]=(cx[i])+MNBRAK_GOLD*(cx[i]-bx[i]);
			   param[i] = u[i];
			   if(param[i] > lim_sup[i]) {param[i] = u[i] = lim_sup[i];}
			   if(param[i] < lim_inf[i]) {param[i] = u[i] = lim_inf[i];}	  
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

/* the two functions below have been added to guarantee that the minimum ML-estimated base frequency and 
   the maximum estimated base frequency remain within the predifined limits of FREQ_MIN [minimum]
   and 1.0 - FREQ_MIN * (states - 1) [maximum].

   This is required because the frequencies are optimized using the xponential function, that is:

   f[i] = e^w[i] / (e^w[0] + ... + e^w[n-1])

   such that the sum of the f[i] is 1.0.

   Thus, the quantity that is actually being optimized are the w[i] and not the frequencies.

   Since the ML estimate of frequencies is done one after the other and independently for each partition,
   the maximum and minimum allowed values for w[i] need to be recalculated at each iteration of the numerical optimization routine 
   depending on the (fixed) values of the w[j] where j = 0...n-1 and j != i

   The equations for the minimum are obtained by transforming the follwing inequality:

   find a value of w[i] such that:

   FREQ_MIN < e^w[i] / e^w[i] + c, where c := sum over all w[j] for j != i which are constant

   after some transformations this becomes:

   (FREQ_MIN * c) / (1 - FREQ_MIN) < e^w[i]

   taking the logarithm this is:

   log(FREQ_MIN) + log(c) - log(1 - FREQ_MIN) < w[i]

   For the maximum allowed frequency which is: 

   max := 1.0 - FREQ_MIN * (states - 1)

   the derivation is analoguous.

*/

static double minFreq(int index, int whichFreq, tree *tr, double absoluteMin)
{
  double 
    min = 0.0,
    *w = tr->partitionData[index].freqExponents,
    c = 0.0;

  int
    states = tr->partitionData[index].states,
    i;

  for(i = 0; i < states; i++)
    if(i != whichFreq)
      c += exp(w[i]);

  min = log(FREQ_MIN) + log(c) - log (1.0 - FREQ_MIN);

  if(0)
    {
      double
	check = exp(min) / (exp(min) + c);
      
      printf("check %f\n", check);    

      printf("min: %f \n", min);
    }
  
  return MAX(min, absoluteMin);
}

static double maxFreq(int index, int whichFreq, tree *tr, double absoluteMax)
{
  double 
    max = 0.0,
    *w = tr->partitionData[index].freqExponents,
    c = 0.0;

  int
    states = tr->partitionData[index].states,
    i;

  for(i = 0; i < states; i++)
    if(i != whichFreq)
      c += exp(w[i]);

  max = log(1.0 - ((double)(states - 1) * FREQ_MIN)) + log(c) - log ((double)(states - 1) * FREQ_MIN);

  if(0)
    {
      double
	check = exp(max) / (exp(max) + c);
      
      printf("check max %f\n", check);    
      
      printf("max: %f \n", max);
    }
  
  return MIN(max, absoluteMax);
}

//#define _DEBUG_MODEL_OPTIMIZATION

static void optParamGeneric(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int rateNumber, double _lim_inf, double _lim_sup, int whichParameterType)
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
    *_x         = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *lim_inf    = (double *)rax_malloc(sizeof(double) * numberOfModels),
    *lim_sup    = (double *)rax_malloc(sizeof(double) * numberOfModels);

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
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].alpha;
		  break;
		case RATE_F:
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].substRates[rateNumber];      
		  break;
#ifdef _HET
		case RATE_F_HET:
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].substRates_TIP[rateNumber];
		  break;
#endif
		case INVAR_F:
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].propInvariant;
		  break;
		case SCALER_F:
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  startValues[pos] = tr->partitionData[index].brLenScaler;
		  break;
		case LXRATE_F:		 
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  assert(rateNumber >= 0 && rateNumber < 4);
		  startValues[pos] = tr->partitionData[index].gammaRates[rateNumber];
		  memcpy(&startRates[pos * 4],   tr->partitionData[index].gammaRates, 4 * sizeof(double)); 
		  memcpy(&startExponents[pos * 4], tr->partitionData[index].weightExponents, 4 * sizeof(double));
		  memcpy(&startWeights[pos * 4], tr->partitionData[index].weights,    4 * sizeof(double));
		  break;
		case LXWEIGHT_F: 		  
		  lim_inf[pos] = _lim_inf;
		  lim_sup[pos] = _lim_sup;
		  assert(rateNumber >= 0 && rateNumber < 4);
		  startValues[pos] = tr->partitionData[index].weightExponents[rateNumber];		  
		  break;
		case FREQ_F:
		  lim_inf[pos] = minFreq(index, rateNumber, tr, _lim_inf);
		  lim_sup[pos] = maxFreq(index, rateNumber, tr, _lim_sup);
		  //lim_inf[pos] = _lim_inf;
		  //lim_sup[pos] = _lim_sup;
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
	      
	  if(_a[pos] < lim_inf[pos]) 
	    _a[pos] = lim_inf[pos];
	  
	  if(_a[pos] > lim_sup[pos]) 
	    _a[pos] = lim_sup[pos];
	      
	  if(_b[pos] < lim_inf[pos]) 
	    _b[pos] = lim_inf[pos];
	  
	  if(_b[pos] > lim_sup[pos]) 
	    _b[pos] = lim_sup[pos];    

	  pos++;
	}
    }                    	     

  assert(pos == numberOfModels);

  brakGeneric(_param, _a, _b, _c, _fa, _fb, _fc, lim_inf, lim_sup, numberOfModels, rateNumber, whichParameterType, tr, ll);
      
  for(k = 0; k < numberOfModels; k++)
    {
      assert(_a[k] >= lim_inf[k] && _a[k] <= lim_sup[k]);
      assert(_b[k] >= lim_inf[k] && _b[k] <= lim_sup[k]);	  
      assert(_c[k] >= lim_inf[k] && _c[k] <= lim_sup[k]);	    
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
		  		 
		  changeModelParameters(index, rateNumber, startValues[pos], whichParameterType, tr);		 
		}
	    }
	  else
	    {
	      //printf("accept\n");

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
#ifdef _HET
    case RATE_F_HET:
      assert(0);
      break;
#endif

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
#ifdef _HET
    case RATE_F_HET:
#endif
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
  rax_free(lim_inf);
  rax_free(lim_sup);

#ifdef _DEBUG_MODEL_OPTIMIZATION
  evaluateGenericInitrav(tr, tr->start);

  if(tr->likelihood < initialLH)
    printf("%1.40f %1.40f\n", tr->likelihood, initialLH);
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


#ifdef _HET
static void optRates(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1;

  evaluateGenericInitrav(tr, tr->start);
  //printf("\nL1 %f\n", tr->likelihood);

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, RATE_MIN, RATE_MAX, RATE_F);
  
  //treeEvaluate(tr, 0.0625);

  evaluateGenericInitrav(tr, tr->start);
  //printf("L2 %f\n\n", tr->likelihood);
  
 
  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, RATE_MIN, RATE_MAX, RATE_F_HET);

  evaluateGenericInitrav(tr, tr->start);
  //printf("L3 %f\n\n", tr->likelihood);
}
#else
static void optRates(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels, int states)
{
  int
    rateNumber,
    numberOfRates = ((states * states - states) / 2) - 1;

  for(rateNumber = 0; rateNumber < numberOfRates; rateNumber++)
    optParamGeneric(tr, modelEpsilon, ll, numberOfModels, rateNumber, RATE_MIN, RATE_MAX, RATE_F);
}
#endif



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
	  if(tr->useJC69)
	    ll->ld[i].valid = FALSE;
	  else
	    {
	      ll->ld[i].valid = TRUE;
	      dnaPartitions++;  
	    }
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






static void optLG4X(tree *tr, double modelEpsilon, linkageList *ll, int numberOfModels)
{
  int 
    i;

  for(i = 0; i < 4; i++)
    {
      optParamGeneric(tr, modelEpsilon, ll, numberOfModels, i, LG4X_RATE_MIN, LG4X_RATE_MAX, LXRATE_F);
      optimizeWeights(tr, modelEpsilon, ll, numberOfModels);   
    }
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
	      z = exp(-p->z[0]);

	      z = fixZ(z);	     

	      p->z[model] = z;
	    }
	}
      else
	{	
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      z = tr->partitionData[model].brLenScaler * tr->storedBrLens[count];
	     
	      z = exp(-z);
	      
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
		  z = exp(-p->next->z[0]);
		  
		  z = fixZ(z);

		  p->next->z[model] = z;

		}
	    }
	  else
	    {	      
	      for(model = 0; model < tr->NumberOfModels; model++)
		{
		  z = tr->partitionData[model].brLenScaler * tr->storedBrLens[count];
		  z = exp(-z);	  		 

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
		  z = exp(-p->next->next->z[0]);
		  
		  z = fixZ(z);		  
		  
		  p->next->next->z[model] = z;
		}
	    }
	  else	  
	    {	     
	       for(model = 0; model < tr->NumberOfModels; model++)
		{
		  z = tr->partitionData[model].brLenScaler * tr->storedBrLens[count];
		  
		  z = exp(-z);	 

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

static void optModel(tree *tr, int numProteinModels, int *bestIndex, double *bestScores, boolean empiricalFreqs)
{
  int
    i,
    model;
    
  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      bestIndex[model] = -1;
      bestScores[model] = unlikely;
    }
      
  for(i = 0; i < numProteinModels; i++)
    {
      for(model = 0; model < tr->NumberOfModels; model++)
	{	   
	  if(tr->partitionData[model].protModels == AUTO)
	    { 
	      if(empiricalFreqs)
		tr->partitionData[model].usePredefinedProtFreqs = FALSE;
	      else
		tr->partitionData[model].usePredefinedProtFreqs = TRUE;
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

      //printf("Subst Model %d Freqs: %s like %f %f\n", i, (empiricalFreqs == TRUE)?"empirical":"fixed", tr->likelihood, tr->perPartitionLH[0]);
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->partitionData[model].protModels == AUTO)
	    {	
	      /*int k;
	      for(k = 0; k < 20; k++)
		printf("%f ", tr->partitionData[model].frequencies[k]);
		printf("\n");*/

 
	      if(tr->perPartitionLH[model] > bestScores[model])
		{
		  bestScores[model] = tr->perPartitionLH[model];
		  bestIndex[model] = i;		      
		}
	    }	      
	}       
    }             
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
	numProteinModels = AUTO,
	*bestIndex = (int*)rax_malloc(sizeof(int) * tr->NumberOfModels),
	*oldIndex  = (int*)rax_malloc(sizeof(int) * tr->NumberOfModels),
	*bestIndexEmpFreqs = (int*)rax_malloc(sizeof(int) * tr->NumberOfModels);

      boolean
	*oldFreqs =  (boolean*)rax_malloc(sizeof(boolean) * tr->NumberOfModels);

      double
	startLH,
	*bestScores         = (double*)rax_malloc(sizeof(double) * tr->NumberOfModels),
	*bestScoresEmpFreqs = (double*)rax_malloc(sizeof(double) * tr->NumberOfModels);

      topolRELL_LIST 
	*rl = (topolRELL_LIST *)rax_malloc(sizeof(topolRELL_LIST));

      char
	*autoModels[4] = {"ML", "BIC", "AIC", "AICc"};

      initTL(rl, tr, 1);
      saveTL(rl, tr, 0);

      evaluateGenericInitrav(tr, tr->start); 

      startLH = tr->likelihood;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  oldIndex[model] = tr->partitionData[model].autoProtModels;
	  oldFreqs[model] = tr->partitionData[model].usePredefinedProtFreqs;
	}

      
      
      optModel(tr, numProteinModels, bestIndex, bestScores, FALSE);
      
      optModel(tr, numProteinModels, bestIndexEmpFreqs, bestScoresEmpFreqs, TRUE);      
     
      printBothOpen("Automatic protein model assignment algorithm using %s criterion:\n\n", autoModels[tr->autoProteinSelectionType]);

      for(model = 0; model < tr->NumberOfModels; model++)
	{	   
	  if(tr->partitionData[model].protModels == AUTO)
	    {	  
	      size_t
		k;

	      int 	       
		bestIndexFixed = bestIndex[model],
		bestIndexEmp = bestIndexEmpFreqs[model];
	      
	      double
		bestLhFixed = bestScores[model],
		bestLhEmp   = bestScoresEmpFreqs[model],
		samples = 0.0,		
		freeParamsFixed = 0.0,
		freeParamsEmp = 0.0;
	      
	      //actually get the number of sites, not the number of patterns!
	      for(k = tr->partitionData[model].lower; k < tr->partitionData[model].upper; k++)
		samples += (double)tr->cdta->aliaswgt[k];

	      //printf("sample %f\n", samples);

	      //include branches in parameter count
	      //catch case where tree contains less taxa than tr->mxtips e.g. in the EPA 
	      assert(tr->ntips <= tr->mxtips);
	      freeParamsFixed = freeParamsEmp = (2 * tr->ntips - 3);
	      freeParamsEmp += 19.0;

	      //printf("%f %f\n", bestLhFixed, bestLhEmp);
	      
	      /* AIC: 2 * (k - lnL)
		 AICc: AIC + (2 * k * (k + 1))/(n - k - 1)
		 BIC: -2 * lnL + k * ln(n) */
	      
	      switch(tr->rateHetModel)
		{
		case CAT:
		  freeParamsFixed += (double)tr->partitionData[model].numberOfCategories;
		  freeParamsEmp += (double)tr->partitionData[model].numberOfCategories;
		  break;
		case GAMMA: 
		  freeParamsFixed += 1.0;
		  freeParamsEmp += 1.0;
		  break;
		case GAMMA_I:
		  freeParamsFixed += 2.0;
		  freeParamsEmp += 2.0;
		  break;
		default:
		  assert(0);
		}
		    

	      switch(tr->autoProteinSelectionType)
		{
		case AUTO_ML:	
		  if(bestLhFixed > bestLhEmp)
		    {
		      tr->partitionData[model].autoProtModels = bestIndexFixed;
		      tr->partitionData[model].usePredefinedProtFreqs = TRUE;
		    }
		  else
		    {
		      tr->partitionData[model].autoProtModels = bestIndexEmp;
		      tr->partitionData[model].usePredefinedProtFreqs = FALSE;
		    }
		  break;
		case AUTO_BIC:
		  { 
		    //BIC: -2 * lnL + k * ln(n)
		    double
		      bicFixed = -2.0 * bestLhFixed + freeParamsFixed * log(samples),
		      bicEmp   = -2.0 * bestLhEmp   + freeParamsEmp   * log(samples);

		    if(bicFixed < bicEmp)
		      {
			tr->partitionData[model].autoProtModels = bestIndexFixed;
			tr->partitionData[model].usePredefinedProtFreqs = TRUE;
		      }
		    else
		      {
			tr->partitionData[model].autoProtModels = bestIndexEmp;
			tr->partitionData[model].usePredefinedProtFreqs = FALSE;
		      }		   
		  }
		  break;
		case AUTO_AIC:
		  {
		    //AIC: 2 * (k - lnL)
		    double
		      aicFixed = 2.0 * (freeParamsFixed - bestLhFixed),
		      aicEmp   = 2.0 * (freeParamsEmp   - bestLhEmp);
		    
		    if(aicFixed < aicEmp)
		      {
			tr->partitionData[model].autoProtModels = bestIndexFixed;
			tr->partitionData[model].usePredefinedProtFreqs = TRUE;
		      }
		    else
		      {
			tr->partitionData[model].autoProtModels = bestIndexEmp;
			tr->partitionData[model].usePredefinedProtFreqs = FALSE;
		      }	
		  }
		  break;
		case AUTO_AICC:
		  { 
		    //AICc: AIC + (2 * k * (k + 1))/(n - k - 1)
		    double
		      aiccFixed,
		      aiccEmp;

		    /* 
		     * Even though samples and freeParamsFixed are fp variables, they are actually integers.
		     * That's why we are comparing with a 0.5 threshold.
		     */
		    
		    if(fabs(samples - freeParamsFixed - 1.0) < 0.5) 		      
		      aiccFixed = 0.0;
		    else 
		      aiccFixed = (2.0 * (freeParamsFixed - bestLhFixed)) + ((2.0 * freeParamsFixed * (freeParamsFixed + 1.0)) / (samples - freeParamsFixed - 1.0));

		    if(fabs(samples - freeParamsEmp - 1.0) < 0.5)
		      aiccEmp = 0.0;
		    else 
		      aiccEmp   = (2.0 * (freeParamsEmp   - bestLhEmp))   + ((2.0 * freeParamsEmp   * (freeParamsEmp + 1.0))   / (samples - freeParamsEmp   - 1.0));

		    if(aiccFixed < aiccEmp)
		      {
			tr->partitionData[model].autoProtModels = bestIndexFixed;
			tr->partitionData[model].usePredefinedProtFreqs = TRUE;
		      }
		    else
		      {
			tr->partitionData[model].autoProtModels = bestIndexEmp;
			tr->partitionData[model].usePredefinedProtFreqs = FALSE;
		      }	
		  }
		  break;
		default:
		  assert(0);
		}

	      initReversibleGTR(tr, model);  
	      printBothOpen("\tPartition: %d best-scoring AA model: %s likelihood %f with %s base frequencies\n", 
			    model, protModels[tr->partitionData[model].autoProtModels], 
			    (tr->partitionData[model].usePredefinedProtFreqs == TRUE)?bestLhFixed:bestLhEmp, (tr->partitionData[model].usePredefinedProtFreqs == TRUE)?"fixed":"empirical");
		  
	    }	 
	}

      printBothOpen("\n\n");

#ifdef _USE_PTHREADS	
      masterBarrier(THREAD_COPY_RATES, tr);	   
#endif
          
      resetBranches(tr);
      evaluateGenericInitrav(tr, tr->start); 
      treeEvaluate(tr, 2.0);    

      //printf("exit %f\n", tr->likelihood);
      
      if(tr->likelihood < startLH)
	{	
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->partitionData[model].protModels == AUTO)
		{
		  tr->partitionData[model].autoProtModels = oldIndex[model];
		  tr->partitionData[model].usePredefinedProtFreqs = oldFreqs[model];
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
      rax_free(bestIndexEmpFreqs);
      rax_free(bestScoresEmpFreqs);
      rax_free(oldFreqs);
    }
}


//#define _DEBUG_MOD_OPT



static void checkTolerance(double l1, double l2)
{
  if(l1 < l2)
    {   
      double 
	
	tolerance = fabs(MAX(l1, l2) * 0.000000000001);      

      if(fabs(l1 - l2) > MIN(0.1, tolerance))
	{
	  printf("Likelihood problem in model optimization l1: %1.40f l2: %1.40f tolerance: %1.40f\n", l1, l2, tolerance);
	  assert(0);	
	}
    }
}


//#define  _DEBUG_MOD_OPT

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
    
#ifndef __BLACKRIM
  if(!(adef->mode == CLASSIFY_ML))
    {
      if(tr->start != tr->nodep[1])
	changedRoot = TRUE;      
      tr->start = tr->nodep[1];
    }
#endif
  
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
	  
	  makeGammaCats(tr->rateHetModel, tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4, tr->useGammaMedian, tr->partitionData[model].propInvariant); 
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

      if(tr->NumberOfModels == 1 && tr->partitionData[0].dataType == DNA_DATA && adef->useBFGS && !(tr->useJC69 || tr->useK80 || tr->useHKY85))
	{	  	 
	  if(optimizeRatesBFGS(tr) == FALSE)
	    {
	      adef->useBFGS = FALSE;
	      optRatesGeneric(tr, modelEpsilon, rateList);
	    }
	}
      else
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
      if(adef->mode != OPTIMIZE_BR_LEN_SCALER)
	printf("after br-len 1 %f\n", tr->likelihood);
      else
	printf("after opt-scaler 1 %f\n", tr->likelihood);     
#endif

      evaluateGenericInitrav(tr, tr->start);

      optBaseFreqs(tr, modelEpsilon, freqList);
      
      evaluateGenericInitrav(tr, tr->start);           

#ifdef _DEBUG_MOD_OPT
      evaluateGenericInitrav(tr, tr->start); 
      printf("after optBaseFreqs 1 %f\n", tr->likelihood);
#endif     

      if(adef->mode != OPTIMIZE_BR_LEN_SCALER)
	treeEvaluate(tr, 0.0625);                     	            
      else 	
	optScaler(tr, modelEpsilon, scalerList); 

#ifdef _DEBUG_MOD_OPT
      evaluateGenericInitrav(tr, tr->start); 
      if(adef->mode != OPTIMIZE_BR_LEN_SCALER)
	printf("after br-len 2 %f\n", tr->likelihood);
      else
	printf("after opt-scaler 2%f\n", tr->likelihood);
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
	  if(adef->mode != OPTIMIZE_BR_LEN_SCALER)
	    printf("after br-len 3 %f\n", tr->likelihood);
	  else
	    printf("after opt-scaler 3%f\n", tr->likelihood);	  
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
      
      //printf("%f \n", tr->likelihood);

      checkTolerance(tr->likelihood, currentLikelihood);
                  
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
    x = -log(x);       
  else
    x = -log(x);

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
      double 
	acc = 0.0;
      
      nodeptr 
	q = p->next;      

      while(q != p)
	{
	  acc += treeLengthRec(q->back, tr, model);
	  q = q->next;
	}

      return (acc + x);
    }
}

double treeLength(tree *tr, int model)
{ 
  /* printf("SCALER: %f\n", tr->partitionData[model].brLenScaler); */

  return treeLengthRec(tr->start->back, tr, model);
}

/********************** bfgs optimization ********************************/
/* most of the code below taken from IQ-Tree http://www.cibiv.at/software/iqtree/
   
   Bui Quang Minh, Minh Anh Thi Nguyen, and Arndt von Haeseler (2013) 
   Ultrafast approximation for phylogenetic bootstrap. 
   Mol. Biol. Evol., 30:1188-1195. (free reprint, DOI: 10.1093/molbev/mst024)
   
*/

static double targetFunk(double *x, int n, tree *tr)
{
  int
    model,
    i = 1;  
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	k;
      
      for(k = 0; k < 5; k++, i++)       
	setRateModel(tr, model, x[i], k);         
      
      initReversibleGTR(tr, model); 
    }

  assert(i == n + 1);
  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_RATES, tr);	
#endif
  
     
  //TODO when optimizing some quantities we actually need to 
  //only evaluate at the root without re-traversing the entire tree 
  
  evaluateGenericInitrav(tr, tr->start);    
  
  //minh confirm that we are actually really trying to minimize, hence 
  //reverting the sign of the lnl is correct below?
  ////MINH: Yes correct, this is minization, returning negative logl does the job 

  return (-1.0 * tr->likelihood);
}

#define ALF 1.0e-4
#define TOLX 1.0e-7
static double maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ?\
        (maxarg1) : (maxarg2))

static void fixBound(double *x, double *lower, double *upper, int n) 
{
  int 
    i;

  for (i = 1; i <= n; i++) 
    {
      if(x[i] < lower[i])
	x[i] = lower[i];
      else 
	if(x[i] > upper[i])
	  x[i] = upper[i];
    }
}

//minh please confirm that *f just is a single value and not potentially an array?
////MINH: yes, *f points to 1 single value

static void lnsrch(int n, double *xold, double fold, double *g, double *p, double *x,
                   double *f, double stpmax, int *check, double *lower, double *upper, tree *tr) 
{
  int i;
  
  double 
    a,alam,alam2=0,alamin,b,disc,f2=0,fold2=0,rhs1,rhs2,slope,sum,temp,
    test,tmplam;

  boolean 
    first_time = TRUE;
  
  *check=0;   

  for (sum=0.0,i=1;i<=n;i++) 
    sum += p[i]*p[i];
  
  sum=sqrt(sum);

  if(sum > stpmax)
    for(i=1;i<=n;i++) 
      p[i] *= stpmax/sum;
  
  for(slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  
  test=0.0;
  
  for (i=1;i<=n;i++) 
    {
      temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
      if(temp > test) 
	test=temp;
    }
  
  alamin=TOLX/test;
  alam=1.0;
  
  /*
    int rep = 0;
    do {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    if (!checkRange(x))
    alam *= 0.5;
    else
    break;
    rep++;
    } while (rep < 10);
  */
 
  
  for (;;) 
    {      
      for(i = 1;i <= n; i++) 		 
	x[i] = xold[i] + alam * p[i]; 	         
      
      fixBound(x, lower, upper, n);

      //minh is the commented code below needed?
      ////MINH: this is indeed not needed
      //checkRange(x);
      
      
      *f=targetFunk(x, n, tr);
      
      if(alam < alamin) 
	{
	  for (i=1;i<=n;i++) x[i]=xold[i];
	  *check=1;
	  return;
	} 
      else 
	if (*f <= fold+ALF*alam*slope) 
	  return;
	else 
	  {
	    if (first_time)
	      tmplam = -slope/(2.0*(*f-fold-slope));
	    else 
	      {
		rhs1 = *f-fold-alam*slope;
		rhs2=f2-fold2-alam2*slope;
		a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
		b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
		if (a == 0.0) 
		  tmplam = -slope/(2.0*b);
		else 
		  {
		    disc=b*b-3.0*a*slope;
		    if (disc<0.0) //nrerror("Roundoff problem in lnsrch.");
		      tmplam = 0.5 * alam;
		    else 
		      if (b <= 0.0) 
			tmplam=(-b+sqrt(disc))/(3.0*a);
		      else 
			tmplam = -slope/(b+sqrt(disc));
		  }
		if (tmplam>0.5*alam)
		  tmplam=0.5*alam;
	      }
	  }
      
      alam2=alam;
      
      f2 = *f;
      
      fold2=fold;
      
      alam = FMAX(tmplam,0.1*alam);
      
      first_time = FALSE;
    }
}
#undef ALF
#undef TOLX


const int MAX_ITER = 3;

static void dfpmin(double *p, int n, double *lower, double *upper, double gtol, int *iter, double *fret, tree *tr,  boolean *bfgsConverged);
static double derivativeFunk(double *x, double *dfx, int n, tree *tr);

//minh is guess over-written by this function?
////MINH: yes, guess will be overwritten by the optimized values. IMPORTANT: all arrays guess,lower,upper,bound_check are indexed from 1 to ndim (not starting from 0!). 
// I don't know why, ask Numerical Recipes guys ;-)
//minh what is the exact meaining of gtol, does it refer to x or f(x)?
////MINH: gtol is the tolerance for the first derivative f'(x). It stops when |f'(x)| < gtol, because optimization means to find the root of f'(x)

static double minimizeMultiDimen(double *guess, int ndim, double *lower, double *upper, boolean *bound_check, double gtol, tree *tr, boolean *bfgsConverged) 
{
  int 
    i, 
    iter,
    count = 0;
    
  double 
    fret, 
    minf = -1.0 * unlikely, 
    //10000000.0, 
    //minh is this some maximum value for (-log likelihood) of the tree does the number need to be large?
    ////MINH: yes it needs to be very large initially, as minf stores the function value at the minimum point
    *minx = (double*)rax_malloc(sizeof(double) * (ndim + 1));
  
  boolean 
    restart;

  int64_t
    seed = 12345;
  
  do 
    {
      dfpmin(guess, ndim, lower, upper, gtol, &iter, &fret, tr, bfgsConverged);
      
      if (fret < minf) 
	{
	  minf = fret;
	  for(i = 1; i <= ndim; i++)
	    minx[i] = guess[i];
	}
      
      count++;
      // restart the search if at the boundary
      // it's likely to end at a local optimum at the boundary
      restart = FALSE;
				
      for (i = 1; i <= ndim; i++)
	if (bound_check[i])
	  if (fabs(guess[i]-lower[i]) < 1e-4 || fabs(guess[i]-upper[i]) < 1e-4) 
	    {
	      restart = TRUE;
	      break;
	    }
		
      if (!restart)
	break;
      
      if (count == MAX_ITER)
	break;
      
      do 
	{
	  for (i = 1; i <= ndim; i++) 
	    {
	      //minh our randum() function yields values between [0,1), confirm that this is correct?
	      ////MINH: correct!
	      guess[i] = randum(&seed) * (upper[i] - lower[i])/3 + lower[i];
	    }
	} 
      while (FALSE);
      
      printf("Restart estimation at the boundary... \n");
      
    } 
  while (count < MAX_ITER);
  
  if(count > 1) 
    {
      for (i = 1; i <= ndim; i++)
	guess[i] = minx[i];
      fret = minf;
    }
				
  rax_free(minx);
	
  return fret;
}


#define ITMAX_BFGS 500
static double sqrarg;
#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0



static void freeAll(double *xi, double *pnew, double **hessin, double *hdg, double *g, double *dg, int n)
{
  int 
    i;

  rax_free(xi);
  rax_free(pnew);
  rax_free(hdg);
  rax_free(g);
  rax_free(dg);

  for(i = 0; i <= n; i++)
    rax_free(hessin[i]);

  rax_free(hessin);
}



static void dfpmin(double *p, int n, double *lower, double *upper, double gtol, int *iter, double *fret, tree *tr, boolean *bfgsConverged) 
{
  int 
    check,
    i,
    its,
    j;
  
  double 
    den,
    fac,
    fad,
    fae,
    fp,
    stpmax,
    sum=0.0,
    sumdg,
    sumxi,
    temp,
    test,
    *dg = (double*)rax_malloc(sizeof(double) * (n + 1)),
    *g  = (double*)rax_malloc(sizeof(double) * (n + 1)),
    *hdg = (double*)rax_malloc(sizeof(double) * (n + 1)),
    **hessin = (double**)rax_malloc(sizeof(double *) * (n + 1)),
    *pnew = (double*)rax_malloc(sizeof(double) * (n + 1)),
    *xi = (double*)rax_malloc(sizeof(double) * (n + 1));

  for(i = 0; i <= n; i++)
    hessin[i] = (double*)rax_malloc(sizeof(double) * (n + 1));

  fp = derivativeFunk(p, g, n, tr);
  
  for (i=1;i<=n;i++) 
    {
      for (j=1;j<=n;j++) 
	hessin[i][j]=0.0;
      hessin[i][i]=1.0;
      xi[i] = -g[i];
      sum += p[i]*p[i];
    }
  
 

  stpmax = STPMX * FMAX(sqrt(sum),(double)n);
  
  for(its = 1; its <= ITMAX_BFGS; its++) 
    {     
      *iter = its;
     
      lnsrch(n, p, fp, g, xi, pnew, fret, stpmax, &check, lower, upper, tr);
        
      fp = *fret;
      
      for (i=1;i<=n;i++) 
	{
	  xi[i]=pnew[i]-p[i];
	  p[i]=pnew[i];
	}
     
      test=0.0;
      
      for (i=1;i<=n;i++) 
	{
	  temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
	  if (temp > test) 
	    test=temp;
	}
      if (test < TOLX) 
	{
	  freeAll(xi, pnew, hessin, hdg, g, dg, n);
	  return;
	}
		
      for (i=1;i<=n;i++) 
	dg[i]=g[i];
      
      derivativeFunk(p, g, n, tr);            

      test=0.0;
      
      den = FMAX(fabs(*fret),1.0);
      //corrected according to this post here:
      //http://www.nr.com/forum/showthread.php?t=1327
      
      //den=FMAX(*fret,1.0);
      
      for (i=1;i<=n;i++) 
	{
	  temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
	  if (temp > test) 
	    test=temp;
	}
      
      if (test < gtol) 
	{
	  freeAll(xi, pnew, hessin, hdg, g, dg, n);
	  return;
	}
      
      for (i=1;i<=n;i++) 
	dg[i]=g[i]-dg[i];
      
      for (i=1;i<=n;i++) 
	{
	  hdg[i]=0.0;
	  for (j=1;j<=n;j++) 
	    hdg[i] += hessin[i][j]*dg[j];
	}
     
      fac=fae=sumdg=sumxi=0.0;      
      
      for (i=1;i<=n;i++) 
	{
	  fac += dg[i]*xi[i];
	  fae += dg[i]*hdg[i];
	  sumdg += SQR(dg[i]);
	  sumxi += SQR(xi[i]);
	}      

      if(fac*fac > EPS*sumdg*sumxi)
	{
	  fac=1.0/fac;
	  fad=1.0/fae;
	  for (i=1;i<=n;i++) 
	    dg[i] = fac * xi[i] - fad * hdg[i];
	  
	  for (i=1;i<=n;i++) 
	    {
	      for (j=1;j<=n;j++) 
		{
		  hessin[i][j] += fac*xi[i]*xi[j]
		    -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
		}
	    }
	}
      
      for (i=1;i<=n;i++) 
	{
	  xi[i]=0.0;
	  for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
	}           
    }

  //printf("too many iterations in dfpmin\n");
  //assert(0);
  
  printBothOpen("\n\n BFGS required too many iterations, reverting to Brent-based optimizer\n\n");
  *bfgsConverged = FALSE;
  
  freeAll(xi, pnew, hessin, hdg, g, dg, n);
}

#undef ITMAX_BFGS
#undef SQR
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
#undef FMAX
#define ERROR_X 1.0e-4

/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
static double derivativeFunk(double *x, double *dfx, int n, tree *tr) 
{ 
  double  
    h, 
    temp,
    fx = targetFunk(x, n, tr);

  int 
    dim;   

  for(dim = 1; dim <= n; dim++)
    {
      temp = x[dim];
      
      h = ERROR_X * fabs(temp);
      if (h == 0.0) 
	h = ERROR_X;
      x[dim] = temp + h;      
      h = x[dim] - temp;
      double t = targetFunk(x, n, tr);
      dfx[dim] = (t - fx) / h;
           	           
      x[dim] = temp;
    }
  
  return fx;
}






static boolean optimizeRatesBFGS(tree *tr)
{
  int 
    model,
    i = 0,
    nGTR = 5 * tr->NumberOfModels;



  double  
    startLH,
    endLH,
    *guessGTR = (double*)rax_malloc(sizeof(double) * (nGTR + 1)),
    *lowerGTR = (double*)rax_malloc(sizeof(double) * (nGTR + 1)),
    *upperGTR = (double*)rax_malloc(sizeof(double) * (nGTR + 1)),
    *rateBuffer = (double*)rax_malloc(sizeof(double) * 6);

  boolean
    bfgsConverged = TRUE,
    *bound_check_GTR = (boolean*)rax_malloc(sizeof(boolean) * (nGTR + 1));
  
  assert(tr->NumberOfModels == 1);
  assert(tr->partitionData[0].dataType == DNA_DATA);
  
  memcpy(rateBuffer, tr->partitionData[0].substRates, sizeof(double) * 6);

  evaluateGenericInitrav(tr, tr->start);
  startLH = tr->likelihood;

  for(model = 0, i = 1; model < tr->NumberOfModels; model++)
    {
      int 
	k;
      
      for(k = 0; k < 5; k++, i++)
	{
	  guessGTR[i] = tr->partitionData[model].substRates[k];	
	  bound_check_GTR[i] = FALSE;
	  //added this max here to prevent num problems when the boundary value is set lower than 
	  //the error margin!
	  lowerGTR[i] = MAX(RATE_MIN + RATE_MIN * ERROR_X, ERROR_X);
	  upperGTR[i] = RATE_MAX - RATE_MAX * ERROR_X;
	  //printf("%1.20f %1.20f %1.20f %1.20f\n", lowerGTR[i], upperGTR[i], RATE_MIN, RATE_MAX);
	}
    }

  assert(i == nGTR + 1);
 
  minimizeMultiDimen(guessGTR, nGTR, lowerGTR, upperGTR, bound_check_GTR, 0.0001, tr,  &bfgsConverged);  
  
  memcpy(tr->partitionData[0].substRates, &guessGTR[1], sizeof(double) * 5);   
  
  initReversibleGTR(tr, 0);
  
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_RATES, tr);	

  
#endif

  evaluateGenericInitrav(tr, tr->start);
  endLH = tr->likelihood; 
  
  if(endLH < startLH)
    {
      printBothOpen("Reverting BFGS ... \n");
      memcpy(tr->partitionData[0].substRates, rateBuffer, sizeof(double) * 6);
      initReversibleGTR(tr, 0);
#ifdef _USE_PTHREADS
      masterBarrier(THREAD_COPY_RATES, tr);	
#endif
      evaluateGenericInitrav(tr, tr->start);
      assert(startLH == tr->likelihood);
    }

  rax_free(guessGTR);
  rax_free(lowerGTR);
  rax_free(upperGTR);
  rax_free(bound_check_GTR);
  rax_free(rateBuffer);

  return bfgsConverged;
}
