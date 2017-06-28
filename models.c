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

#include "axml.h"

extern int optimizeRatesInvocations;
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeTTRatioInvocations;
extern int optimizeInvarInvocations;

extern const unsigned int bitVectorSecondary[256];
extern const unsigned int bitVector32[33];
extern const unsigned int bitVectorAA[23];
extern const unsigned int bitVectorIdentity[256];
extern const unsigned int mask32[32];
extern const partitionLengths pLengths[MAX_MODEL];

#ifdef _USE_PTHREADS
extern volatile int NumberOfThreads;
#endif


extern char *protModels[NUM_PROT_MODELS];



static void smoothFreqs(const int n, double *pfreqs, double *dst, pInfo *partitionData)
{
  //double
  //  acc = 0.0;

  int 
    countScale = 0, 
    l,
    loopCounter = 0;  
  

  
  for(l = 0; l < n; l++)
    if(pfreqs[l] < FREQ_MIN)
      countScale++;
  

  /*for(l = 0; l < n; l++)
    if(pfreqs[l] == 0.0)
    countScale++;*/

  if(countScale > 0)
    {	     
      while(countScale > 0)
	{
	  double 
	    correction = 0.0,
	    factor = 1.0;
	  
	  for(l = 0; l < n; l++)
	    {
	      if(pfreqs[l] == 0.0)		  
		correction += FREQ_MIN;		   		  
	      else
		if(pfreqs[l] < FREQ_MIN)		    
		  {
		    correction += (FREQ_MIN - pfreqs[l]);
		    factor -= (FREQ_MIN - pfreqs[l]);
		  }
	    }		      	    	    	  	  

	  countScale = 0;
	  
	  for(l = 0; l < n; l++)
	    {		    
	      if(pfreqs[l] >= FREQ_MIN)		      
		pfreqs[l] = pfreqs[l] - (pfreqs[l] * correction * factor);	
	      else
		pfreqs[l] = FREQ_MIN;
	      
	      if(pfreqs[l] < FREQ_MIN)
		countScale++;
	    }
	  assert(loopCounter < 100);
	  loopCounter++;
	}		    
    }

  for(l = 0; l < n; l++)
    {
      //acc += pfreqs[l];
      dst[l] = pfreqs[l];
    }

  //printf("sum %f\n", acc);

  if(partitionData->nonGTR)
    {
      int k;

      assert(partitionData->dataType == SECONDARY_DATA_7 || partitionData->dataType == SECONDARY_DATA_6 || partitionData->dataType == SECONDARY_DATA);
       
      for(l = 0; l < n; l++)
	{
	  int count = 1;	
	  
	  for(k = 0; k < n; k++)
	    {
	      if(k != l && partitionData->frequencyGrouping[l] == partitionData->frequencyGrouping[k])
		{
		  count++;
		  dst[l] += pfreqs[k];
		}
	    }
	  dst[l] /= ((double)count);
	}            
     }  
}
	    

static void genericBaseFrequencies(tree *tr, const int numFreqs, rawdata *rdta, cruncheddata *cdta, int lower, int upper, int model, boolean smoothFrequencies,
				   const unsigned int *bitMask)
{
  double 
    wj, 
    acc,
    pfreqs[64], 
    sumf[64],   
    temp[64];
 
  int 
    statesPresent[64],
    countStatesPresent = 0,
    i, 
    j, 
    k, 
    l;

  unsigned char  *yptr;  
	  

  for(i = 0; i < numFreqs; i++)
    statesPresent[i] = 0;

  for(i = 0; i < rdta->numsp; i++) 
    {
      yptr = &(rdta->y0[((size_t)i) * (tr->originalCrunchedLength)]);
      
      for(j = lower; j < upper; j++) 
	{
	  unsigned int
	    state = 0,
	    stateCount = 0,
	    code = bitMask[yptr[j]];

	  if(precomputed16_bitcount(code) == 1)
	    {
	      for(k = 0; k < numFreqs; k++)
		if(code & mask32[k])
		  {
		    state = k;
		    stateCount++;
		  }

	      assert(stateCount == 1);

	      statesPresent[state] = 1;
	    }		    		    	  	
	}
    }
	      
  for(i = 0, countStatesPresent = 0; i < numFreqs; i++)
    if(statesPresent[i] == 1)
      countStatesPresent++;


  if(tr->partitionData[model].optimizeBaseFrequencies)
    {   
      for(l = 0; l < numFreqs; l++)	    
	tr->partitionData[model].frequencies[l] = 1.0 / ((double)numFreqs);
    }
  else
    {
      
      for(l = 0; l < numFreqs; l++)	    
	pfreqs[l] = 1.0 / ((double)numFreqs);
      
      for (k = 1; k <= 8; k++) 
	{	     	   	    	      			
	  for(l = 0; l < numFreqs; l++)
	    sumf[l] = 0.0;
	  
	  for (i = 0; i < rdta->numsp; i++) 
	    {		 
	      yptr =  &(rdta->y0[((size_t)i) * ((size_t)tr->originalCrunchedLength)]);
	      
	      for(j = lower; j < upper; j++) 
		{
		  unsigned int code = bitMask[yptr[j]];
		  assert(code >= 1);
		  
		  for(l = 0; l < numFreqs; l++)
		    {
		      if((code >> l) & 1)
			temp[l] = pfreqs[l];
		      else
			temp[l] = 0.0;
		    }		      	      
		  
		  for(l = 0, acc = 0.0; l < numFreqs; l++)
		    {
		      if(temp[l] != 0.0)
			acc += temp[l];
		    }
		  
		  wj = ((double)cdta->aliaswgt[j]) / acc;
		  
		  for(l = 0; l < numFreqs; l++)
		    {
		      if(temp[l] != 0.0)		    
			sumf[l] += wj * temp[l];			     				   			     		   
		    }
		}
	    }	    	      
	  
	  for(l = 0, acc = 0.0; l < numFreqs; l++)
	    {
	      if(sumf[l] != 0.0)
		acc += sumf[l];
	    }
	  
	  for(l = 0; l < numFreqs; l++)
	    pfreqs[l] = sumf[l] / acc;	     
	}

      if(countStatesPresent < numFreqs)
	{
	  printf("Partition %s number %d has a problem, the number of expected states is %d the number of states that are present is %d.\n", 
		 tr->partitionData[model].partitionName, model, numFreqs, countStatesPresent);
	  printf("Please go and fix your data!\n\n");
	}
      
      if(smoothFrequencies)         
	{
	  //printf("smoothing!\n");
	  smoothFreqs(numFreqs, pfreqs,  tr->partitionData[model].frequencies, &(tr->partitionData[model]));	   
	}
      else    
	{
	  boolean 
	    zeroFreq = FALSE;
	  
	  char 
	    typeOfData[1024];
	  
	  getDataTypeString(tr, model, typeOfData);  
	  
	  for(l = 0; l < numFreqs; l++)
	    {
	      if(pfreqs[l] == 0.0)
		{
		  printBothOpen("Empirical base frequency for state number %d is equal to zero in %s data partition %s\n", l, typeOfData, tr->partitionData[model].partitionName);
		  printBothOpen("Since this is probably not what you want to do, RAxML will soon exit.\n\n");
		  zeroFreq = TRUE;
		}
	    }
	  
	  if(zeroFreq)
	    exit(-1);
	  
	  for(l = 0; l < numFreqs; l++)
	    {
	      assert(pfreqs[l] > 0.0);
	      tr->partitionData[model].frequencies[l] = pfreqs[l];
	    }     
	}  
    }
}








static void baseFrequenciesGTR(rawdata *rdta, cruncheddata *cdta, tree *tr)
{  
  int 
    model,
    lower,
    upper,
    states;

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      lower = tr->partitionData[model].lower;
      upper = tr->partitionData[model].upper;	  	 
      states = tr->partitionData[model].states;
	
      switch(tr->partitionData[model].dataType)
	{
	case GENERIC_32:
	  switch(tr->multiStateModel)
	    {
	    case ORDERED_MULTI_STATE:
	    case MK_MULTI_STATE:	      
	      {	       
		int i;
		double 
		  freq = 1.0 / (double)states,
		  acc = 0.0;

		for(i = 0; i < states; i++)
		  {
		    acc += freq;
		    tr->partitionData[model].frequencies[i] = freq;
		    /*printf("%f \n", freq);*/
		  }
		/*printf("Frequency Deviation: %1.60f\n", acc);*/
	      }
	      break;
	     case GTR_MULTI_STATE:
	      genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, TRUE,
				     bitVector32);
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case GENERIC_64:	 
	  assert(0);
	  break;
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	case SECONDARY_DATA:
	case AA_DATA:	
	case BINARY_DATA:
	  {
	    //printf("Smooth freqs: %d\n", getSmoothFreqs(tr->partitionData[model].dataType));
	    genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, 
				   getSmoothFreqs(tr->partitionData[model].dataType),
				   getBitVector(tr->partitionData[model].dataType));	  	 
	  }
	  break;
	case DNA_DATA:
	  if(tr->useK80 || tr->useJC69)
	    {	   
	      assert(!tr->partitionData[model].optimizeBaseFrequencies);

	      tr->partitionData[model].frequencies[0] = 0.25;
	      tr->partitionData[model].frequencies[1] = 0.25;
	      tr->partitionData[model].frequencies[2] = 0.25;
	      tr->partitionData[model].frequencies[3] = 0.25;
	    }
	  else
	    genericBaseFrequencies(tr, states, rdta, cdta, lower, upper, model, 
				   getSmoothFreqs(tr->partitionData[model].dataType),
				   getBitVector(tr->partitionData[model].dataType));
	  break;
	default:
	  assert(0);     
	}      
    }
  
  return;
}

static void putWAG(double *ext_initialRates)
{ 
  double
    scaler,
    q[20][20],
    daa[400];

  int 
    i,
    j,
    r;

  daa[ 1*20+ 0] =  55.15710; daa[ 2*20+ 0] =  50.98480; daa[ 2*20+ 1] =  63.53460; 
  daa[ 3*20+ 0] =  73.89980; daa[ 3*20+ 1] =  14.73040; daa[ 3*20+ 2] = 542.94200; 
  daa[ 4*20+ 0] = 102.70400; daa[ 4*20+ 1] =  52.81910; daa[ 4*20+ 2] =  26.52560; 
  daa[ 4*20+ 3] =   3.02949; daa[ 5*20+ 0] =  90.85980; daa[ 5*20+ 1] = 303.55000; 
  daa[ 5*20+ 2] = 154.36400; daa[ 5*20+ 3] =  61.67830; daa[ 5*20+ 4] =   9.88179; 
  daa[ 6*20+ 0] = 158.28500; daa[ 6*20+ 1] =  43.91570; daa[ 6*20+ 2] =  94.71980; 
  daa[ 6*20+ 3] = 617.41600; daa[ 6*20+ 4] =   2.13520; daa[ 6*20+ 5] = 546.94700; 
  daa[ 7*20+ 0] = 141.67200; daa[ 7*20+ 1] =  58.46650; daa[ 7*20+ 2] = 112.55600; 
  daa[ 7*20+ 3] =  86.55840; daa[ 7*20+ 4] =  30.66740; daa[ 7*20+ 5] =  33.00520; 
  daa[ 7*20+ 6] =  56.77170; daa[ 8*20+ 0] =  31.69540; daa[ 8*20+ 1] = 213.71500; 
  daa[ 8*20+ 2] = 395.62900; daa[ 8*20+ 3] =  93.06760; daa[ 8*20+ 4] =  24.89720; 
  daa[ 8*20+ 5] = 429.41100; daa[ 8*20+ 6] =  57.00250; daa[ 8*20+ 7] =  24.94100; 
  daa[ 9*20+ 0] =  19.33350; daa[ 9*20+ 1] =  18.69790; daa[ 9*20+ 2] =  55.42360; 
  daa[ 9*20+ 3] =   3.94370; daa[ 9*20+ 4] =  17.01350; daa[ 9*20+ 5] =  11.39170; 
  daa[ 9*20+ 6] =  12.73950; daa[ 9*20+ 7] =   3.04501; daa[ 9*20+ 8] =  13.81900; 
  daa[10*20+ 0] =  39.79150; daa[10*20+ 1] =  49.76710; daa[10*20+ 2] =  13.15280; 
  daa[10*20+ 3] =   8.48047; daa[10*20+ 4] =  38.42870; daa[10*20+ 5] =  86.94890; 
  daa[10*20+ 6] =  15.42630; daa[10*20+ 7] =   6.13037; daa[10*20+ 8] =  49.94620; 
  daa[10*20+ 9] = 317.09700; daa[11*20+ 0] =  90.62650; daa[11*20+ 1] = 535.14200; 
  daa[11*20+ 2] = 301.20100; daa[11*20+ 3] =  47.98550; daa[11*20+ 4] =   7.40339; 
  daa[11*20+ 5] = 389.49000; daa[11*20+ 6] = 258.44300; daa[11*20+ 7] =  37.35580; 
  daa[11*20+ 8] =  89.04320; daa[11*20+ 9] =  32.38320; daa[11*20+10] =  25.75550; 
  daa[12*20+ 0] =  89.34960; daa[12*20+ 1] =  68.31620; daa[12*20+ 2] =  19.82210; 
  daa[12*20+ 3] =  10.37540; daa[12*20+ 4] =  39.04820; daa[12*20+ 5] = 154.52600; 
  daa[12*20+ 6] =  31.51240; daa[12*20+ 7] =  17.41000; daa[12*20+ 8] =  40.41410; 
  daa[12*20+ 9] = 425.74600; daa[12*20+10] = 485.40200; daa[12*20+11] =  93.42760; 
  daa[13*20+ 0] =  21.04940; daa[13*20+ 1] =  10.27110; daa[13*20+ 2] =   9.61621; 
  daa[13*20+ 3] =   4.67304; daa[13*20+ 4] =  39.80200; daa[13*20+ 5] =   9.99208; 
  daa[13*20+ 6] =   8.11339; daa[13*20+ 7] =   4.99310; daa[13*20+ 8] =  67.93710; 
  daa[13*20+ 9] = 105.94700; daa[13*20+10] = 211.51700; daa[13*20+11] =   8.88360; 
  daa[13*20+12] = 119.06300; daa[14*20+ 0] = 143.85500; daa[14*20+ 1] =  67.94890; 
  daa[14*20+ 2] =  19.50810; daa[14*20+ 3] =  42.39840; daa[14*20+ 4] =  10.94040; 
  daa[14*20+ 5] =  93.33720; daa[14*20+ 6] =  68.23550; daa[14*20+ 7] =  24.35700; 
  daa[14*20+ 8] =  69.61980; daa[14*20+ 9] =   9.99288; daa[14*20+10] =  41.58440; 
  daa[14*20+11] =  55.68960; daa[14*20+12] =  17.13290; daa[14*20+13] =  16.14440; 
  daa[15*20+ 0] = 337.07900; daa[15*20+ 1] = 122.41900; daa[15*20+ 2] = 397.42300; 
  daa[15*20+ 3] = 107.17600; daa[15*20+ 4] = 140.76600; daa[15*20+ 5] = 102.88700; 
  daa[15*20+ 6] =  70.49390; daa[15*20+ 7] = 134.18200; daa[15*20+ 8] =  74.01690; 
  daa[15*20+ 9] =  31.94400; daa[15*20+10] =  34.47390; daa[15*20+11] =  96.71300; 
  daa[15*20+12] =  49.39050; daa[15*20+13] =  54.59310; daa[15*20+14] = 161.32800; 
  daa[16*20+ 0] = 212.11100; daa[16*20+ 1] =  55.44130; daa[16*20+ 2] = 203.00600; 
  daa[16*20+ 3] =  37.48660; daa[16*20+ 4] =  51.29840; daa[16*20+ 5] =  85.79280; 
  daa[16*20+ 6] =  82.27650; daa[16*20+ 7] =  22.58330; daa[16*20+ 8] =  47.33070; 
  daa[16*20+ 9] = 145.81600; daa[16*20+10] =  32.66220; daa[16*20+11] = 138.69800; 
  daa[16*20+12] = 151.61200; daa[16*20+13] =  17.19030; daa[16*20+14] =  79.53840; 
  daa[16*20+15] = 437.80200; daa[17*20+ 0] =  11.31330; daa[17*20+ 1] = 116.39200; 
  daa[17*20+ 2] =   7.19167; daa[17*20+ 3] =  12.97670; daa[17*20+ 4] =  71.70700; 
  daa[17*20+ 5] =  21.57370; daa[17*20+ 6] =  15.65570; daa[17*20+ 7] =  33.69830; 
  daa[17*20+ 8] =  26.25690; daa[17*20+ 9] =  21.24830; daa[17*20+10] =  66.53090; 
  daa[17*20+11] =  13.75050; daa[17*20+12] =  51.57060; daa[17*20+13] = 152.96400; 
  daa[17*20+14] =  13.94050; daa[17*20+15] =  52.37420; daa[17*20+16] =  11.08640; 
  daa[18*20+ 0] =  24.07350; daa[18*20+ 1] =  38.15330; daa[18*20+ 2] = 108.60000; 
  daa[18*20+ 3] =  32.57110; daa[18*20+ 4] =  54.38330; daa[18*20+ 5] =  22.77100; 
  daa[18*20+ 6] =  19.63030; daa[18*20+ 7] =  10.36040; daa[18*20+ 8] = 387.34400; 
  daa[18*20+ 9] =  42.01700; daa[18*20+10] =  39.86180; daa[18*20+11] =  13.32640; 
  daa[18*20+12] =  42.84370; daa[18*20+13] = 645.42800; daa[18*20+14] =  21.60460; 
  daa[18*20+15] =  78.69930; daa[18*20+16] =  29.11480; daa[18*20+17] = 248.53900; 
  daa[19*20+ 0] = 200.60100; daa[19*20+ 1] =  25.18490; daa[19*20+ 2] =  19.62460; 
  daa[19*20+ 3] =  15.23350; daa[19*20+ 4] = 100.21400; daa[19*20+ 5] =  30.12810; 
  daa[19*20+ 6] =  58.87310; daa[19*20+ 7] =  18.72470; daa[19*20+ 8] =  11.83580; 
  daa[19*20+ 9] = 782.13000; daa[19*20+10] = 180.03400; daa[19*20+11] =  30.54340; 
  daa[19*20+12] = 205.84500; daa[19*20+13] =  64.98920; daa[19*20+14] =  31.48870; 
  daa[19*20+15] =  23.27390; daa[19*20+16] = 138.82300; daa[19*20+17] =  36.53690; 
  daa[19*20+18] =  31.47300; 

  for(i = 0; i < 20; i++)
    for(j = 0; j < 20; j++)
      q[i][j] = 0.0;

  for (i=0; i<20; i++)  
    for (j=0; j<i; j++)               
      daa[j*20+i] = daa[i*20+j];

  for(i = 0; i < 19; i++)
    for(j = i + 1; j < 20; j++)      
      q[i][j] = daa[i * 20 + j];

  
  /*
    for (i=0; i<20; i++) 
    {
      for (j=0; j<20; j++)
	printf("%1.2f ", q[i][j]);
      printf("\n");
    }
    printf("\n");

    printf("%f\n", q[18][19]);
  */

  scaler = 1.0 / q[18][19];

  

  for(i = 0; i < 19; i++)
    for(j = i + 1; j < 20; j++)      
      q[i][j] *= scaler;

  for(i = 0, r = 0; i < 19; i++)          
    for(j = i + 1; j < 20; j++)      
      ext_initialRates[r++] = q[i][j];           
      
  /*
    for (i=0; i<20; i++) 
    {
      for (j=0; j<20; j++)
	printf("%1.2f ", q[i][j]);
      printf("\n");
    }
    printf("\n");
  */

}

static void makeAASubstMat(double *daa, double *f, double *rates, double *freqs)
{
  int 
    i, j, r = 0;

  for(i = 1; i < 20; i++)
    for(j = 0; j < i; j++)
      {
	daa[i * 20 + j] = rates[r];
	r++;
      }
  
  assert(r == 190);
  
  for(i = 0; i < 20; i++)
    f[i] = freqs[i];
}

static void initProtMat(double f[20], int proteinMatrix, double *ext_initialRates, int model, tree *tr, int lg4_index)
{ 
  double q[20][20];
  double daa[400], max, temp;
  int i, j, r;
  double *initialRates = ext_initialRates;
  double scaler;

  {
    switch(proteinMatrix)
      {
	case PROT_FILE:
	  memcpy(daa, tr->partitionData[model].externalAAMatrix,         400 * sizeof(double));
	  memcpy(f,   &(tr->partitionData[model].externalAAMatrix[400]), 20  * sizeof(double));
	  break;
	case DAYHOFF:
	  {	
	    daa[ 1*20+ 0] =   27.00; daa[ 2*20+ 0] =   98.00; daa[ 2*20+ 1] =   32.00; daa[ 3*20+ 0] =  120.00;
	    daa[ 3*20+ 1] =    0.00; daa[ 3*20+ 2] =  905.00; daa[ 4*20+ 0] =   36.00; daa[ 4*20+ 1] =   23.00;
	    daa[ 4*20+ 2] =    0.00; daa[ 4*20+ 3] =    0.00; daa[ 5*20+ 0] =   89.00; daa[ 5*20+ 1] =  246.00;
	    daa[ 5*20+ 2] =  103.00; daa[ 5*20+ 3] =  134.00; daa[ 5*20+ 4] =    0.00; daa[ 6*20+ 0] =  198.00;
	    daa[ 6*20+ 1] =    1.00; daa[ 6*20+ 2] =  148.00; daa[ 6*20+ 3] = 1153.00; daa[ 6*20+ 4] =    0.00;
	    daa[ 6*20+ 5] =  716.00; daa[ 7*20+ 0] =  240.00; daa[ 7*20+ 1] =    9.00; daa[ 7*20+ 2] =  139.00;
	    daa[ 7*20+ 3] =  125.00; daa[ 7*20+ 4] =   11.00; daa[ 7*20+ 5] =   28.00; daa[ 7*20+ 6] =   81.00;
	    daa[ 8*20+ 0] =   23.00; daa[ 8*20+ 1] =  240.00; daa[ 8*20+ 2] =  535.00; daa[ 8*20+ 3] =   86.00;
	    daa[ 8*20+ 4] =   28.00; daa[ 8*20+ 5] =  606.00; daa[ 8*20+ 6] =   43.00; daa[ 8*20+ 7] =   10.00;
	    daa[ 9*20+ 0] =   65.00; daa[ 9*20+ 1] =   64.00; daa[ 9*20+ 2] =   77.00; daa[ 9*20+ 3] =   24.00;
	    daa[ 9*20+ 4] =   44.00; daa[ 9*20+ 5] =   18.00; daa[ 9*20+ 6] =   61.00; daa[ 9*20+ 7] =    0.00;
	    daa[ 9*20+ 8] =    7.00; daa[10*20+ 0] =   41.00; daa[10*20+ 1] =   15.00; daa[10*20+ 2] =   34.00;
	    daa[10*20+ 3] =    0.00; daa[10*20+ 4] =    0.00; daa[10*20+ 5] =   73.00; daa[10*20+ 6] =   11.00;
	    daa[10*20+ 7] =    7.00; daa[10*20+ 8] =   44.00; daa[10*20+ 9] =  257.00; daa[11*20+ 0] =   26.00;
	    daa[11*20+ 1] =  464.00; daa[11*20+ 2] =  318.00; daa[11*20+ 3] =   71.00; daa[11*20+ 4] =    0.00;
	    daa[11*20+ 5] =  153.00; daa[11*20+ 6] =   83.00; daa[11*20+ 7] =   27.00; daa[11*20+ 8] =   26.00;
	    daa[11*20+ 9] =   46.00; daa[11*20+10] =   18.00; daa[12*20+ 0] =   72.00; daa[12*20+ 1] =   90.00;
	    daa[12*20+ 2] =    1.00; daa[12*20+ 3] =    0.00; daa[12*20+ 4] =    0.00; daa[12*20+ 5] =  114.00;
	    daa[12*20+ 6] =   30.00; daa[12*20+ 7] =   17.00; daa[12*20+ 8] =    0.00; daa[12*20+ 9] =  336.00;
	    daa[12*20+10] =  527.00; daa[12*20+11] =  243.00; daa[13*20+ 0] =   18.00; daa[13*20+ 1] =   14.00;
	    daa[13*20+ 2] =   14.00; daa[13*20+ 3] =    0.00; daa[13*20+ 4] =    0.00; daa[13*20+ 5] =    0.00;
	    daa[13*20+ 6] =    0.00; daa[13*20+ 7] =   15.00; daa[13*20+ 8] =   48.00; daa[13*20+ 9] =  196.00;
	    daa[13*20+10] =  157.00; daa[13*20+11] =    0.00; daa[13*20+12] =   92.00; daa[14*20+ 0] =  250.00;
	    daa[14*20+ 1] =  103.00; daa[14*20+ 2] =   42.00; daa[14*20+ 3] =   13.00; daa[14*20+ 4] =   19.00;
	    daa[14*20+ 5] =  153.00; daa[14*20+ 6] =   51.00; daa[14*20+ 7] =   34.00; daa[14*20+ 8] =   94.00;
	    daa[14*20+ 9] =   12.00; daa[14*20+10] =   32.00; daa[14*20+11] =   33.00; daa[14*20+12] =   17.00;
	    daa[14*20+13] =   11.00; daa[15*20+ 0] =  409.00; daa[15*20+ 1] =  154.00; daa[15*20+ 2] =  495.00;
	    daa[15*20+ 3] =   95.00; daa[15*20+ 4] =  161.00; daa[15*20+ 5] =   56.00; daa[15*20+ 6] =   79.00;
	    daa[15*20+ 7] =  234.00; daa[15*20+ 8] =   35.00; daa[15*20+ 9] =   24.00; daa[15*20+10] =   17.00;
	    daa[15*20+11] =   96.00; daa[15*20+12] =   62.00; daa[15*20+13] =   46.00; daa[15*20+14] =  245.00;
	    daa[16*20+ 0] =  371.00; daa[16*20+ 1] =   26.00; daa[16*20+ 2] =  229.00; daa[16*20+ 3] =   66.00;
	    daa[16*20+ 4] =   16.00; daa[16*20+ 5] =   53.00; daa[16*20+ 6] =   34.00; daa[16*20+ 7] =   30.00;
	    daa[16*20+ 8] =   22.00; daa[16*20+ 9] =  192.00; daa[16*20+10] =   33.00; daa[16*20+11] =  136.00;
	    daa[16*20+12] =  104.00; daa[16*20+13] =   13.00; daa[16*20+14] =   78.00; daa[16*20+15] =  550.00;
	    daa[17*20+ 0] =    0.00; daa[17*20+ 1] =  201.00; daa[17*20+ 2] =   23.00; daa[17*20+ 3] =    0.00;
	    daa[17*20+ 4] =    0.00; daa[17*20+ 5] =    0.00; daa[17*20+ 6] =    0.00; daa[17*20+ 7] =    0.00;
	    daa[17*20+ 8] =   27.00; daa[17*20+ 9] =    0.00; daa[17*20+10] =   46.00; daa[17*20+11] =    0.00;
	    daa[17*20+12] =    0.00; daa[17*20+13] =   76.00; daa[17*20+14] =    0.00; daa[17*20+15] =   75.00;
	    daa[17*20+16] =    0.00; daa[18*20+ 0] =   24.00; daa[18*20+ 1] =    8.00; daa[18*20+ 2] =   95.00;
	    daa[18*20+ 3] =    0.00; daa[18*20+ 4] =   96.00; daa[18*20+ 5] =    0.00; daa[18*20+ 6] =   22.00;
	    daa[18*20+ 7] =    0.00; daa[18*20+ 8] =  127.00; daa[18*20+ 9] =   37.00; daa[18*20+10] =   28.00;
	    daa[18*20+11] =   13.00; daa[18*20+12] =    0.00; daa[18*20+13] =  698.00; daa[18*20+14] =    0.00;
	    daa[18*20+15] =   34.00; daa[18*20+16] =   42.00; daa[18*20+17] =   61.00; daa[19*20+ 0] =  208.00;
	    daa[19*20+ 1] =   24.00; daa[19*20+ 2] =   15.00; daa[19*20+ 3] =   18.00; daa[19*20+ 4] =   49.00;
	    daa[19*20+ 5] =   35.00; daa[19*20+ 6] =   37.00; daa[19*20+ 7] =   54.00; daa[19*20+ 8] =   44.00;
	    daa[19*20+ 9] =  889.00; daa[19*20+10] =  175.00; daa[19*20+11] =   10.00; daa[19*20+12] =  258.00;
	    daa[19*20+13] =   12.00; daa[19*20+14] =   48.00; daa[19*20+15] =   30.00; daa[19*20+16] =  157.00;
	    daa[19*20+17] =    0.00; daa[19*20+18] =   28.00;	    	    


	    /*
	      f[ 0] = 0.087000; f[ 1] = 0.041000; f[ 2] = 0.040000; f[ 3] = 0.047000;
	      f[ 4] = 0.034000; f[ 5] = 0.038000; f[ 6] = 0.050000; f[ 7] = 0.089000;
	      f[ 8] = 0.034000; f[ 9] = 0.037000; f[10] = 0.085000; f[11] = 0.080000;
	      f[12] = 0.014000; f[13] = 0.040000; f[14] = 0.051000; f[15] = 0.070000;
	      f[16] = 0.058000; f[17] = 0.011000; f[18] = 0.030000; f[19] = 0.064000;
	    */
	    f[ 0] = 0.087127; f[ 1] = 0.040904; f[ 2] = 0.040432; f[ 3] = 0.046872;
	    f[ 4] = 0.033474; f[ 5] = 0.038255; f[ 6] = 0.049530; f[ 7] = 0.088612;
	    f[ 8] = 0.033618; f[ 9] = 0.036886; f[10] = 0.085357; f[11] = 0.080482;
	    f[12] = 0.014753; f[13] = 0.039772; f[14] = 0.050680; f[15] = 0.069577;
	    f[16] = 0.058542; f[17] = 0.010494; f[18] = 0.029916; f[19] = 0.064717;
	  }
	  break;
	case DCMUT:
	  {	
	    daa[ 1*20+ 0] =   26.78280; daa[ 2*20+ 0] =   98.44740; daa[ 2*20+ 1] =   32.70590; daa[ 3*20+ 0] =  119.98050; 
	    daa[ 3*20+ 1] =    0.00000; daa[ 3*20+ 2] =  893.15150; daa[ 4*20+ 0] =   36.00160; daa[ 4*20+ 1] =   23.23740; 
	    daa[ 4*20+ 2] =    0.00000; daa[ 4*20+ 3] =    0.00000; daa[ 5*20+ 0] =   88.77530; daa[ 5*20+ 1] =  243.99390; 
	    daa[ 5*20+ 2] =  102.85090; daa[ 5*20+ 3] =  134.85510; daa[ 5*20+ 4] =    0.00000; daa[ 6*20+ 0] =  196.11670; 
	    daa[ 6*20+ 1] =    0.00000; daa[ 6*20+ 2] =  149.34090; daa[ 6*20+ 3] = 1138.86590; daa[ 6*20+ 4] =    0.00000; 
	    daa[ 6*20+ 5] =  708.60220; daa[ 7*20+ 0] =  238.61110; daa[ 7*20+ 1] =    8.77910; daa[ 7*20+ 2] =  138.53520; 
	    daa[ 7*20+ 3] =  124.09810; daa[ 7*20+ 4] =   10.72780; daa[ 7*20+ 5] =   28.15810; daa[ 7*20+ 6] =   81.19070; 
	    daa[ 8*20+ 0] =   22.81160; daa[ 8*20+ 1] =  238.31480; daa[ 8*20+ 2] =  529.00240; daa[ 8*20+ 3] =   86.82410; 
	    daa[ 8*20+ 4] =   28.27290; daa[ 8*20+ 5] =  601.16130; daa[ 8*20+ 6] =   43.94690; daa[ 8*20+ 7] =   10.68020; 
	    daa[ 9*20+ 0] =   65.34160; daa[ 9*20+ 1] =   63.26290; daa[ 9*20+ 2] =   76.80240; daa[ 9*20+ 3] =   23.92480; 
	    daa[ 9*20+ 4] =   43.80740; daa[ 9*20+ 5] =   18.03930; daa[ 9*20+ 6] =   60.95260; daa[ 9*20+ 7] =    0.00000; 
	    daa[ 9*20+ 8] =    7.69810; daa[10*20+ 0] =   40.64310; daa[10*20+ 1] =   15.49240; daa[10*20+ 2] =   34.11130; 
	    daa[10*20+ 3] =    0.00000; daa[10*20+ 4] =    0.00000; daa[10*20+ 5] =   73.07720; daa[10*20+ 6] =   11.28800; 
	    daa[10*20+ 7] =    7.15140; daa[10*20+ 8] =   44.35040; daa[10*20+ 9] =  255.66850; daa[11*20+ 0] =   25.86350; 
	    daa[11*20+ 1] =  461.01240; daa[11*20+ 2] =  314.83710; daa[11*20+ 3] =   71.69130; daa[11*20+ 4] =    0.00000; 
	    daa[11*20+ 5] =  151.90780; daa[11*20+ 6] =   83.00780; daa[11*20+ 7] =   26.76830; daa[11*20+ 8] =   27.04750; 
	    daa[11*20+ 9] =   46.08570; daa[11*20+10] =   18.06290; daa[12*20+ 0] =   71.78400; daa[12*20+ 1] =   89.63210; 
	    daa[12*20+ 2] =    0.00000; daa[12*20+ 3] =    0.00000; daa[12*20+ 4] =    0.00000; daa[12*20+ 5] =  112.74990; 
	    daa[12*20+ 6] =   30.48030; daa[12*20+ 7] =   17.03720; daa[12*20+ 8] =    0.00000; daa[12*20+ 9] =  333.27320; 
	    daa[12*20+10] =  523.01150; daa[12*20+11] =  241.17390; daa[13*20+ 0] =   18.36410; daa[13*20+ 1] =   13.69060; 
	    daa[13*20+ 2] =   13.85030; daa[13*20+ 3] =    0.00000; daa[13*20+ 4] =    0.00000; daa[13*20+ 5] =    0.00000; 
	    daa[13*20+ 6] =    0.00000; daa[13*20+ 7] =   15.34780; daa[13*20+ 8] =   47.59270; daa[13*20+ 9] =  195.19510; 
	    daa[13*20+10] =  156.51600; daa[13*20+11] =    0.00000; daa[13*20+12] =   92.18600; daa[14*20+ 0] =  248.59200; 
	    daa[14*20+ 1] =  102.83130; daa[14*20+ 2] =   41.92440; daa[14*20+ 3] =   13.39400; daa[14*20+ 4] =   18.75500; 
	    daa[14*20+ 5] =  152.61880; daa[14*20+ 6] =   50.70030; daa[14*20+ 7] =   34.71530; daa[14*20+ 8] =   93.37090; 
	    daa[14*20+ 9] =   11.91520; daa[14*20+10] =   31.62580; daa[14*20+11] =   33.54190; daa[14*20+12] =   17.02050; 
	    daa[14*20+13] =   11.05060; daa[15*20+ 0] =  405.18700; daa[15*20+ 1] =  153.15900; daa[15*20+ 2] =  488.58920; 
	    daa[15*20+ 3] =   95.60970; daa[15*20+ 4] =  159.83560; daa[15*20+ 5] =   56.18280; daa[15*20+ 6] =   79.39990; 
	    daa[15*20+ 7] =  232.22430; daa[15*20+ 8] =   35.36430; daa[15*20+ 9] =   24.79550; daa[15*20+10] =   17.14320; 
	    daa[15*20+11] =   95.45570; daa[15*20+12] =   61.99510; daa[15*20+13] =   45.99010; daa[15*20+14] =  242.72020; 
	    daa[16*20+ 0] =  368.03650; daa[16*20+ 1] =   26.57450; daa[16*20+ 2] =  227.16970; daa[16*20+ 3] =   66.09300; 
	    daa[16*20+ 4] =   16.23660; daa[16*20+ 5] =   52.56510; daa[16*20+ 6] =   34.01560; daa[16*20+ 7] =   30.66620; 
	    daa[16*20+ 8] =   22.63330; daa[16*20+ 9] =  190.07390; daa[16*20+10] =   33.10900; daa[16*20+11] =  135.05990; 
	    daa[16*20+12] =  103.15340; daa[16*20+13] =   13.66550; daa[16*20+14] =   78.28570; daa[16*20+15] =  543.66740; 
	    daa[17*20+ 0] =    0.00000; daa[17*20+ 1] =  200.13750; daa[17*20+ 2] =   22.49680; daa[17*20+ 3] =    0.00000; 
	    daa[17*20+ 4] =    0.00000; daa[17*20+ 5] =    0.00000; daa[17*20+ 6] =    0.00000; daa[17*20+ 7] =    0.00000; 
	    daa[17*20+ 8] =   27.05640; daa[17*20+ 9] =    0.00000; daa[17*20+10] =   46.17760; daa[17*20+11] =    0.00000; 
	    daa[17*20+12] =    0.00000; daa[17*20+13] =   76.23540; daa[17*20+14] =    0.00000; daa[17*20+15] =   74.08190; 
	    daa[17*20+16] =    0.00000; daa[18*20+ 0] =   24.41390; daa[18*20+ 1] =    7.80120; daa[18*20+ 2] =   94.69400; 
	    daa[18*20+ 3] =    0.00000; daa[18*20+ 4] =   95.31640; daa[18*20+ 5] =    0.00000; daa[18*20+ 6] =   21.47170; 
	    daa[18*20+ 7] =    0.00000; daa[18*20+ 8] =  126.54000; daa[18*20+ 9] =   37.48340; daa[18*20+10] =   28.65720; 
	    daa[18*20+11] =   13.21420; daa[18*20+12] =    0.00000; daa[18*20+13] =  695.26290; daa[18*20+14] =    0.00000; 
	    daa[18*20+15] =   33.62890; daa[18*20+16] =   41.78390; daa[18*20+17] =   60.80700; daa[19*20+ 0] =  205.95640; 
	    daa[19*20+ 1] =   24.03680; daa[19*20+ 2] =   15.80670; daa[19*20+ 3] =   17.83160; daa[19*20+ 4] =   48.46780; 
	    daa[19*20+ 5] =   34.69830; daa[19*20+ 6] =   36.72500; daa[19*20+ 7] =   53.81650; daa[19*20+ 8] =   43.87150; 
	    daa[19*20+ 9] =  881.00380; daa[19*20+10] =  174.51560; daa[19*20+11] =   10.38500; daa[19*20+12] =  256.59550; 
	    daa[19*20+13] =   12.36060; daa[19*20+14] =   48.50260; daa[19*20+15] =   30.38360; daa[19*20+16] =  156.19970; 
	    daa[19*20+17] =    0.00000; daa[19*20+18] =   27.93790;   	    	   

	    /*f[ 0] = 0.08700; f[ 1] = 0.04100; f[ 2] = 0.04000; f[ 3] = 0.04700;
	    f[ 4] = 0.03300; f[ 5] = 0.03800; f[ 6] = 0.04900; f[ 7] = 0.08900;
	    f[ 8] = 0.03400; f[ 9] = 0.03700; f[10] = 0.08500; f[11] = 0.08000;
	    f[12] = 0.01500; f[13] = 0.04000; f[14] = 0.05200; f[15] = 0.06900;
	    f[16] = 0.05900; f[17] = 0.01000; f[18] = 0.03000; f[19] = 0.06500;*/
	    
	    f[ 0] = 0.087127; f[ 1] = 0.040904; f[ 2] = 0.040432; f[ 3] = 0.046872;
	    f[ 4] = 0.033474; f[ 5] = 0.038255; f[ 6] = 0.049530; f[ 7] = 0.088612;
	    f[ 8] = 0.033619; f[ 9] = 0.036886; f[10] = 0.085357; f[11] = 0.080481;
	    f[12] = 0.014753; f[13] = 0.039772; f[14] = 0.050680; f[15] = 0.069577;
	    f[16] = 0.058542; f[17] = 0.010494; f[18] = 0.029916; f[19] = 0.064717;

	  }
	  break;
	case JTT:
	  {
	    daa[ 1*20+ 0] =   58.00; daa[ 2*20+ 0] =   54.00; daa[ 2*20+ 1] =   45.00; daa[ 3*20+ 0] =   81.00;
	    daa[ 3*20+ 1] =   16.00; daa[ 3*20+ 2] =  528.00; daa[ 4*20+ 0] =   56.00; daa[ 4*20+ 1] =  113.00;
	    daa[ 4*20+ 2] =   34.00; daa[ 4*20+ 3] =   10.00; daa[ 5*20+ 0] =   57.00; daa[ 5*20+ 1] =  310.00;
	    daa[ 5*20+ 2] =   86.00; daa[ 5*20+ 3] =   49.00; daa[ 5*20+ 4] =    9.00; daa[ 6*20+ 0] =  105.00;
	    daa[ 6*20+ 1] =   29.00; daa[ 6*20+ 2] =   58.00; daa[ 6*20+ 3] =  767.00; daa[ 6*20+ 4] =    5.00;
	    daa[ 6*20+ 5] =  323.00; daa[ 7*20+ 0] =  179.00; daa[ 7*20+ 1] =  137.00; daa[ 7*20+ 2] =   81.00;
	    daa[ 7*20+ 3] =  130.00; daa[ 7*20+ 4] =   59.00; daa[ 7*20+ 5] =   26.00; daa[ 7*20+ 6] =  119.00;
	    daa[ 8*20+ 0] =   27.00; daa[ 8*20+ 1] =  328.00; daa[ 8*20+ 2] =  391.00; daa[ 8*20+ 3] =  112.00;
	    daa[ 8*20+ 4] =   69.00; daa[ 8*20+ 5] =  597.00; daa[ 8*20+ 6] =   26.00; daa[ 8*20+ 7] =   23.00;
	    daa[ 9*20+ 0] =   36.00; daa[ 9*20+ 1] =   22.00; daa[ 9*20+ 2] =   47.00; daa[ 9*20+ 3] =   11.00;
	    daa[ 9*20+ 4] =   17.00; daa[ 9*20+ 5] =    9.00; daa[ 9*20+ 6] =   12.00; daa[ 9*20+ 7] =    6.00;
	    daa[ 9*20+ 8] =   16.00; daa[10*20+ 0] =   30.00; daa[10*20+ 1] =   38.00; daa[10*20+ 2] =   12.00;
	    daa[10*20+ 3] =    7.00; daa[10*20+ 4] =   23.00; daa[10*20+ 5] =   72.00; daa[10*20+ 6] =    9.00;
	    daa[10*20+ 7] =    6.00; daa[10*20+ 8] =   56.00; daa[10*20+ 9] =  229.00; daa[11*20+ 0] =   35.00;
	    daa[11*20+ 1] =  646.00; daa[11*20+ 2] =  263.00; daa[11*20+ 3] =   26.00; daa[11*20+ 4] =    7.00;
	    daa[11*20+ 5] =  292.00; daa[11*20+ 6] =  181.00; daa[11*20+ 7] =   27.00; daa[11*20+ 8] =   45.00;
	    daa[11*20+ 9] =   21.00; daa[11*20+10] =   14.00; daa[12*20+ 0] =   54.00; daa[12*20+ 1] =   44.00;
	    daa[12*20+ 2] =   30.00; daa[12*20+ 3] =   15.00; daa[12*20+ 4] =   31.00; daa[12*20+ 5] =   43.00;
	    daa[12*20+ 6] =   18.00; daa[12*20+ 7] =   14.00; daa[12*20+ 8] =   33.00; daa[12*20+ 9] =  479.00;
	    daa[12*20+10] =  388.00; daa[12*20+11] =   65.00; daa[13*20+ 0] =   15.00; daa[13*20+ 1] =    5.00;
	    daa[13*20+ 2] =   10.00; daa[13*20+ 3] =    4.00; daa[13*20+ 4] =   78.00; daa[13*20+ 5] =    4.00;
	    daa[13*20+ 6] =    5.00; daa[13*20+ 7] =    5.00; daa[13*20+ 8] =   40.00; daa[13*20+ 9] =   89.00;
	    daa[13*20+10] =  248.00; daa[13*20+11] =    4.00; daa[13*20+12] =   43.00; daa[14*20+ 0] =  194.00;
	    daa[14*20+ 1] =   74.00; daa[14*20+ 2] =   15.00; daa[14*20+ 3] =   15.00; daa[14*20+ 4] =   14.00;
	    daa[14*20+ 5] =  164.00; daa[14*20+ 6] =   18.00; daa[14*20+ 7] =   24.00; daa[14*20+ 8] =  115.00;
	    daa[14*20+ 9] =   10.00; daa[14*20+10] =  102.00; daa[14*20+11] =   21.00; daa[14*20+12] =   16.00;
	    daa[14*20+13] =   17.00; daa[15*20+ 0] =  378.00; daa[15*20+ 1] =  101.00; daa[15*20+ 2] =  503.00;
	    daa[15*20+ 3] =   59.00; daa[15*20+ 4] =  223.00; daa[15*20+ 5] =   53.00; daa[15*20+ 6] =   30.00;
	    daa[15*20+ 7] =  201.00; daa[15*20+ 8] =   73.00; daa[15*20+ 9] =   40.00; daa[15*20+10] =   59.00;
	    daa[15*20+11] =   47.00; daa[15*20+12] =   29.00; daa[15*20+13] =   92.00; daa[15*20+14] =  285.00;
	    daa[16*20+ 0] =  475.00; daa[16*20+ 1] =   64.00; daa[16*20+ 2] =  232.00; daa[16*20+ 3] =   38.00;
	    daa[16*20+ 4] =   42.00; daa[16*20+ 5] =   51.00; daa[16*20+ 6] =   32.00; daa[16*20+ 7] =   33.00;
	    daa[16*20+ 8] =   46.00; daa[16*20+ 9] =  245.00; daa[16*20+10] =   25.00; daa[16*20+11] =  103.00;
	    daa[16*20+12] =  226.00; daa[16*20+13] =   12.00; daa[16*20+14] =  118.00; daa[16*20+15] =  477.00;
	    daa[17*20+ 0] =    9.00; daa[17*20+ 1] =  126.00; daa[17*20+ 2] =    8.00; daa[17*20+ 3] =    4.00;
	    daa[17*20+ 4] =  115.00; daa[17*20+ 5] =   18.00; daa[17*20+ 6] =   10.00; daa[17*20+ 7] =   55.00;
	    daa[17*20+ 8] =    8.00; daa[17*20+ 9] =    9.00; daa[17*20+10] =   52.00; daa[17*20+11] =   10.00;
	    daa[17*20+12] =   24.00; daa[17*20+13] =   53.00; daa[17*20+14] =    6.00; daa[17*20+15] =   35.00;
	    daa[17*20+16] =   12.00; daa[18*20+ 0] =   11.00; daa[18*20+ 1] =   20.00; daa[18*20+ 2] =   70.00;
	    daa[18*20+ 3] =   46.00; daa[18*20+ 4] =  209.00; daa[18*20+ 5] =   24.00; daa[18*20+ 6] =    7.00;
	    daa[18*20+ 7] =    8.00; daa[18*20+ 8] =  573.00; daa[18*20+ 9] =   32.00; daa[18*20+10] =   24.00;
	    daa[18*20+11] =    8.00; daa[18*20+12] =   18.00; daa[18*20+13] =  536.00; daa[18*20+14] =   10.00;
	    daa[18*20+15] =   63.00; daa[18*20+16] =   21.00; daa[18*20+17] =   71.00; daa[19*20+ 0] =  298.00;
	    daa[19*20+ 1] =   17.00; daa[19*20+ 2] =   16.00; daa[19*20+ 3] =   31.00; daa[19*20+ 4] =   62.00;
	    daa[19*20+ 5] =   20.00; daa[19*20+ 6] =   45.00; daa[19*20+ 7] =   47.00; daa[19*20+ 8] =   11.00;
	    daa[19*20+ 9] =  961.00; daa[19*20+10] =  180.00; daa[19*20+11] =   14.00; daa[19*20+12] =  323.00;
	    daa[19*20+13] =   62.00; daa[19*20+14] =   23.00; daa[19*20+15] =   38.00; daa[19*20+16] =  112.00;
	    daa[19*20+17] =   25.00; daa[19*20+18] =   16.00;
	    	    
	    /*f[ 0] = 0.07700; f[ 1] = 0.05200; f[ 2] = 0.04200; f[ 3] = 0.05100;
	    f[ 4] = 0.02000; f[ 5] = 0.04100; f[ 6] = 0.06200; f[ 7] = 0.07300;
	    f[ 8] = 0.02300; f[ 9] = 0.05400; f[10] = 0.09200; f[11] = 0.05900;
	    f[12] = 0.02400; f[13] = 0.04000; f[14] = 0.05100; f[15] = 0.06900;
	    f[16] = 0.05800; f[17] = 0.01400; f[18] = 0.03200; f[19] = 0.06600;*/

	    f[ 0] = 0.076748; f[ 1] = 0.051691; f[ 2] = 0.042645; f[ 3] = 0.051544;
	    f[ 4] = 0.019803; f[ 5] = 0.040752; f[ 6] = 0.061830; f[ 7] = 0.073152;
	    f[ 8] = 0.022944; f[ 9] = 0.053761; f[10] = 0.091904; f[11] = 0.058676;
	    f[12] = 0.023826; f[13] = 0.040126; f[14] = 0.050901; f[15] = 0.068765;
	    f[16] = 0.058565; f[17] = 0.014261; f[18] = 0.032102; f[19] = 0.066004;
	  }
	  break;
	case  MTREV:
	  {
	    daa[ 1*20+ 0] =   23.18; daa[ 2*20+ 0] =   26.95; daa[ 2*20+ 1] =   13.24; daa[ 3*20+ 0] =   17.67;
	    daa[ 3*20+ 1] =    1.90; daa[ 3*20+ 2] =  794.38; daa[ 4*20+ 0] =   59.93; daa[ 4*20+ 1] =  103.33;
	    daa[ 4*20+ 2] =   58.94; daa[ 4*20+ 3] =    1.90; daa[ 5*20+ 0] =    1.90; daa[ 5*20+ 1] =  220.99;
	    daa[ 5*20+ 2] =  173.56; daa[ 5*20+ 3] =   55.28; daa[ 5*20+ 4] =   75.24; daa[ 6*20+ 0] =    9.77;
	    daa[ 6*20+ 1] =    1.90; daa[ 6*20+ 2] =   63.05; daa[ 6*20+ 3] =  583.55; daa[ 6*20+ 4] =    1.90;
	    daa[ 6*20+ 5] =  313.56; daa[ 7*20+ 0] =  120.71; daa[ 7*20+ 1] =   23.03; daa[ 7*20+ 2] =   53.30;
	    daa[ 7*20+ 3] =   56.77; daa[ 7*20+ 4] =   30.71; daa[ 7*20+ 5] =    6.75; daa[ 7*20+ 6] =   28.28;
	    daa[ 8*20+ 0] =   13.90; daa[ 8*20+ 1] =  165.23; daa[ 8*20+ 2] =  496.13; daa[ 8*20+ 3] =  113.99;
	    daa[ 8*20+ 4] =  141.49; daa[ 8*20+ 5] =  582.40; daa[ 8*20+ 6] =   49.12; daa[ 8*20+ 7] =    1.90;
	    daa[ 9*20+ 0] =   96.49; daa[ 9*20+ 1] =    1.90; daa[ 9*20+ 2] =   27.10; daa[ 9*20+ 3] =    4.34;
	    daa[ 9*20+ 4] =   62.73; daa[ 9*20+ 5] =    8.34; daa[ 9*20+ 6] =    3.31; daa[ 9*20+ 7] =    5.98;
	    daa[ 9*20+ 8] =   12.26; daa[10*20+ 0] =   25.46; daa[10*20+ 1] =   15.58; daa[10*20+ 2] =   15.16;
	    daa[10*20+ 3] =    1.90; daa[10*20+ 4] =   25.65; daa[10*20+ 5] =   39.70; daa[10*20+ 6] =    1.90;
	    daa[10*20+ 7] =    2.41; daa[10*20+ 8] =   11.49; daa[10*20+ 9] =  329.09; daa[11*20+ 0] =    8.36;
	    daa[11*20+ 1] =  141.40; daa[11*20+ 2] =  608.70; daa[11*20+ 3] =    2.31; daa[11*20+ 4] =    1.90;
	    daa[11*20+ 5] =  465.58; daa[11*20+ 6] =  313.86; daa[11*20+ 7] =   22.73; daa[11*20+ 8] =  127.67;
	    daa[11*20+ 9] =   19.57; daa[11*20+10] =   14.88; daa[12*20+ 0] =  141.88; daa[12*20+ 1] =    1.90;
	    daa[12*20+ 2] =   65.41; daa[12*20+ 3] =    1.90; daa[12*20+ 4] =    6.18; daa[12*20+ 5] =   47.37;
	    daa[12*20+ 6] =    1.90; daa[12*20+ 7] =    1.90; daa[12*20+ 8] =   11.97; daa[12*20+ 9] =  517.98;
	    daa[12*20+10] =  537.53; daa[12*20+11] =   91.37; daa[13*20+ 0] =    6.37; daa[13*20+ 1] =    4.69;
	    daa[13*20+ 2] =   15.20; daa[13*20+ 3] =    4.98; daa[13*20+ 4] =   70.80; daa[13*20+ 5] =   19.11;
	    daa[13*20+ 6] =    2.67; daa[13*20+ 7] =    1.90; daa[13*20+ 8] =   48.16; daa[13*20+ 9] =   84.67;
	    daa[13*20+10] =  216.06; daa[13*20+11] =    6.44; daa[13*20+12] =   90.82; daa[14*20+ 0] =   54.31;
	    daa[14*20+ 1] =   23.64; daa[14*20+ 2] =   73.31; daa[14*20+ 3] =   13.43; daa[14*20+ 4] =   31.26;
	    daa[14*20+ 5] =  137.29; daa[14*20+ 6] =   12.83; daa[14*20+ 7] =    1.90; daa[14*20+ 8] =   60.97;
	    daa[14*20+ 9] =   20.63; daa[14*20+10] =   40.10; daa[14*20+11] =   50.10; daa[14*20+12] =   18.84;
	    daa[14*20+13] =   17.31; daa[15*20+ 0] =  387.86; daa[15*20+ 1] =    6.04; daa[15*20+ 2] =  494.39;
	    daa[15*20+ 3] =   69.02; daa[15*20+ 4] =  277.05; daa[15*20+ 5] =   54.11; daa[15*20+ 6] =   54.71;
	    daa[15*20+ 7] =  125.93; daa[15*20+ 8] =   77.46; daa[15*20+ 9] =   47.70; daa[15*20+10] =   73.61;
	    daa[15*20+11] =  105.79; daa[15*20+12] =  111.16; daa[15*20+13] =   64.29; daa[15*20+14] =  169.90;
	    daa[16*20+ 0] =  480.72; daa[16*20+ 1] =    2.08; daa[16*20+ 2] =  238.46; daa[16*20+ 3] =   28.01;
	    daa[16*20+ 4] =  179.97; daa[16*20+ 5] =   94.93; daa[16*20+ 6] =   14.82; daa[16*20+ 7] =   11.17;
	    daa[16*20+ 8] =   44.78; daa[16*20+ 9] =  368.43; daa[16*20+10] =  126.40; daa[16*20+11] =  136.33;
	    daa[16*20+12] =  528.17; daa[16*20+13] =   33.85; daa[16*20+14] =  128.22; daa[16*20+15] =  597.21;
	    daa[17*20+ 0] =    1.90; daa[17*20+ 1] =   21.95; daa[17*20+ 2] =   10.68; daa[17*20+ 3] =   19.86;
	    daa[17*20+ 4] =   33.60; daa[17*20+ 5] =    1.90; daa[17*20+ 6] =    1.90; daa[17*20+ 7] =   10.92;
	    daa[17*20+ 8] =    7.08; daa[17*20+ 9] =    1.90; daa[17*20+10] =   32.44; daa[17*20+11] =   24.00;
	    daa[17*20+12] =   21.71; daa[17*20+13] =    7.84; daa[17*20+14] =    4.21; daa[17*20+15] =   38.58;
	    daa[17*20+16] =    9.99; daa[18*20+ 0] =    6.48; daa[18*20+ 1] =    1.90; daa[18*20+ 2] =  191.36;
	    daa[18*20+ 3] =   21.21; daa[18*20+ 4] =  254.77; daa[18*20+ 5] =   38.82; daa[18*20+ 6] =   13.12;
	    daa[18*20+ 7] =    3.21; daa[18*20+ 8] =  670.14; daa[18*20+ 9] =   25.01; daa[18*20+10] =   44.15;
	    daa[18*20+11] =   51.17; daa[18*20+12] =   39.96; daa[18*20+13] =  465.58; daa[18*20+14] =   16.21;
	    daa[18*20+15] =   64.92; daa[18*20+16] =   38.73; daa[18*20+17] =   26.25; daa[19*20+ 0] =  195.06;
	    daa[19*20+ 1] =    7.64; daa[19*20+ 2] =    1.90; daa[19*20+ 3] =    1.90; daa[19*20+ 4] =    1.90;
	    daa[19*20+ 5] =   19.00; daa[19*20+ 6] =   21.14; daa[19*20+ 7] =    2.53; daa[19*20+ 8] =    1.90;
	    daa[19*20+ 9] = 1222.94; daa[19*20+10] =   91.67; daa[19*20+11] =    1.90; daa[19*20+12] =  387.54;
	    daa[19*20+13] =    6.35; daa[19*20+14] =    8.23; daa[19*20+15] =    1.90; daa[19*20+16] =  204.54;
	    daa[19*20+17] =    5.37; daa[19*20+18] =    1.90;
	    
	    
	    f[ 0] = 0.072000; f[ 1] = 0.019000; f[ 2] = 0.039000; f[ 3] = 0.019000;
	    f[ 4] = 0.006000; f[ 5] = 0.025000; f[ 6] = 0.024000; f[ 7] = 0.056000;
	    f[ 8] = 0.028000; f[ 9] = 0.088000; f[10] = 0.169000; f[11] = 0.023000;
	    f[12] = 0.054000; f[13] = 0.061000; f[14] = 0.054000; f[15] = 0.072000;
	    f[16] = 0.086000; f[17] = 0.029000; f[18] = 0.033000; f[19] = 0.043000;
	  }
	  break;
	case WAG:
	  {
	    daa[ 1*20+ 0] =  55.15710; daa[ 2*20+ 0] =  50.98480; daa[ 2*20+ 1] =  63.53460; 
	    daa[ 3*20+ 0] =  73.89980; daa[ 3*20+ 1] =  14.73040; daa[ 3*20+ 2] = 542.94200; 
	    daa[ 4*20+ 0] = 102.70400; daa[ 4*20+ 1] =  52.81910; daa[ 4*20+ 2] =  26.52560; 
	    daa[ 4*20+ 3] =   3.02949; daa[ 5*20+ 0] =  90.85980; daa[ 5*20+ 1] = 303.55000; 
	    daa[ 5*20+ 2] = 154.36400; daa[ 5*20+ 3] =  61.67830; daa[ 5*20+ 4] =   9.88179; 
	    daa[ 6*20+ 0] = 158.28500; daa[ 6*20+ 1] =  43.91570; daa[ 6*20+ 2] =  94.71980; 
	    daa[ 6*20+ 3] = 617.41600; daa[ 6*20+ 4] =   2.13520; daa[ 6*20+ 5] = 546.94700; 
	    daa[ 7*20+ 0] = 141.67200; daa[ 7*20+ 1] =  58.46650; daa[ 7*20+ 2] = 112.55600; 
	    daa[ 7*20+ 3] =  86.55840; daa[ 7*20+ 4] =  30.66740; daa[ 7*20+ 5] =  33.00520; 
	    daa[ 7*20+ 6] =  56.77170; daa[ 8*20+ 0] =  31.69540; daa[ 8*20+ 1] = 213.71500; 
	    daa[ 8*20+ 2] = 395.62900; daa[ 8*20+ 3] =  93.06760; daa[ 8*20+ 4] =  24.89720; 
	    daa[ 8*20+ 5] = 429.41100; daa[ 8*20+ 6] =  57.00250; daa[ 8*20+ 7] =  24.94100; 
	    daa[ 9*20+ 0] =  19.33350; daa[ 9*20+ 1] =  18.69790; daa[ 9*20+ 2] =  55.42360; 
	    daa[ 9*20+ 3] =   3.94370; daa[ 9*20+ 4] =  17.01350; daa[ 9*20+ 5] =  11.39170; 
	    daa[ 9*20+ 6] =  12.73950; daa[ 9*20+ 7] =   3.04501; daa[ 9*20+ 8] =  13.81900; 
	    daa[10*20+ 0] =  39.79150; daa[10*20+ 1] =  49.76710; daa[10*20+ 2] =  13.15280; 
	    daa[10*20+ 3] =   8.48047; daa[10*20+ 4] =  38.42870; daa[10*20+ 5] =  86.94890; 
	    daa[10*20+ 6] =  15.42630; daa[10*20+ 7] =   6.13037; daa[10*20+ 8] =  49.94620; 
	    daa[10*20+ 9] = 317.09700; daa[11*20+ 0] =  90.62650; daa[11*20+ 1] = 535.14200; 
	    daa[11*20+ 2] = 301.20100; daa[11*20+ 3] =  47.98550; daa[11*20+ 4] =   7.40339; 
	    daa[11*20+ 5] = 389.49000; daa[11*20+ 6] = 258.44300; daa[11*20+ 7] =  37.35580; 
	    daa[11*20+ 8] =  89.04320; daa[11*20+ 9] =  32.38320; daa[11*20+10] =  25.75550; 
	    daa[12*20+ 0] =  89.34960; daa[12*20+ 1] =  68.31620; daa[12*20+ 2] =  19.82210; 
	    daa[12*20+ 3] =  10.37540; daa[12*20+ 4] =  39.04820; daa[12*20+ 5] = 154.52600; 
	    daa[12*20+ 6] =  31.51240; daa[12*20+ 7] =  17.41000; daa[12*20+ 8] =  40.41410; 
	    daa[12*20+ 9] = 425.74600; daa[12*20+10] = 485.40200; daa[12*20+11] =  93.42760; 
	    daa[13*20+ 0] =  21.04940; daa[13*20+ 1] =  10.27110; daa[13*20+ 2] =   9.61621; 
	    daa[13*20+ 3] =   4.67304; daa[13*20+ 4] =  39.80200; daa[13*20+ 5] =   9.99208; 
	    daa[13*20+ 6] =   8.11339; daa[13*20+ 7] =   4.99310; daa[13*20+ 8] =  67.93710; 
	    daa[13*20+ 9] = 105.94700; daa[13*20+10] = 211.51700; daa[13*20+11] =   8.88360; 
	    daa[13*20+12] = 119.06300; daa[14*20+ 0] = 143.85500; daa[14*20+ 1] =  67.94890; 
	    daa[14*20+ 2] =  19.50810; daa[14*20+ 3] =  42.39840; daa[14*20+ 4] =  10.94040; 
	    daa[14*20+ 5] =  93.33720; daa[14*20+ 6] =  68.23550; daa[14*20+ 7] =  24.35700; 
	    daa[14*20+ 8] =  69.61980; daa[14*20+ 9] =   9.99288; daa[14*20+10] =  41.58440; 
	    daa[14*20+11] =  55.68960; daa[14*20+12] =  17.13290; daa[14*20+13] =  16.14440; 
	    daa[15*20+ 0] = 337.07900; daa[15*20+ 1] = 122.41900; daa[15*20+ 2] = 397.42300; 
	    daa[15*20+ 3] = 107.17600; daa[15*20+ 4] = 140.76600; daa[15*20+ 5] = 102.88700; 
	    daa[15*20+ 6] =  70.49390; daa[15*20+ 7] = 134.18200; daa[15*20+ 8] =  74.01690; 
	    daa[15*20+ 9] =  31.94400; daa[15*20+10] =  34.47390; daa[15*20+11] =  96.71300; 
	    daa[15*20+12] =  49.39050; daa[15*20+13] =  54.59310; daa[15*20+14] = 161.32800; 
	    daa[16*20+ 0] = 212.11100; daa[16*20+ 1] =  55.44130; daa[16*20+ 2] = 203.00600; 
	    daa[16*20+ 3] =  37.48660; daa[16*20+ 4] =  51.29840; daa[16*20+ 5] =  85.79280; 
	    daa[16*20+ 6] =  82.27650; daa[16*20+ 7] =  22.58330; daa[16*20+ 8] =  47.33070; 
	    daa[16*20+ 9] = 145.81600; daa[16*20+10] =  32.66220; daa[16*20+11] = 138.69800; 
	    daa[16*20+12] = 151.61200; daa[16*20+13] =  17.19030; daa[16*20+14] =  79.53840; 
	    daa[16*20+15] = 437.80200; daa[17*20+ 0] =  11.31330; daa[17*20+ 1] = 116.39200; 
	    daa[17*20+ 2] =   7.19167; daa[17*20+ 3] =  12.97670; daa[17*20+ 4] =  71.70700; 
	    daa[17*20+ 5] =  21.57370; daa[17*20+ 6] =  15.65570; daa[17*20+ 7] =  33.69830; 
	    daa[17*20+ 8] =  26.25690; daa[17*20+ 9] =  21.24830; daa[17*20+10] =  66.53090; 
	    daa[17*20+11] =  13.75050; daa[17*20+12] =  51.57060; daa[17*20+13] = 152.96400; 
	    daa[17*20+14] =  13.94050; daa[17*20+15] =  52.37420; daa[17*20+16] =  11.08640; 
	    daa[18*20+ 0] =  24.07350; daa[18*20+ 1] =  38.15330; daa[18*20+ 2] = 108.60000; 
	    daa[18*20+ 3] =  32.57110; daa[18*20+ 4] =  54.38330; daa[18*20+ 5] =  22.77100; 
	    daa[18*20+ 6] =  19.63030; daa[18*20+ 7] =  10.36040; daa[18*20+ 8] = 387.34400; 
	    daa[18*20+ 9] =  42.01700; daa[18*20+10] =  39.86180; daa[18*20+11] =  13.32640; 
	    daa[18*20+12] =  42.84370; daa[18*20+13] = 645.42800; daa[18*20+14] =  21.60460; 
	    daa[18*20+15] =  78.69930; daa[18*20+16] =  29.11480; daa[18*20+17] = 248.53900; 
	    daa[19*20+ 0] = 200.60100; daa[19*20+ 1] =  25.18490; daa[19*20+ 2] =  19.62460; 
	    daa[19*20+ 3] =  15.23350; daa[19*20+ 4] = 100.21400; daa[19*20+ 5] =  30.12810; 
	    daa[19*20+ 6] =  58.87310; daa[19*20+ 7] =  18.72470; daa[19*20+ 8] =  11.83580; 
	    daa[19*20+ 9] = 782.13000; daa[19*20+10] = 180.03400; daa[19*20+11] =  30.54340; 
	    daa[19*20+12] = 205.84500; daa[19*20+13] =  64.98920; daa[19*20+14] =  31.48870; 
	    daa[19*20+15] =  23.27390; daa[19*20+16] = 138.82300; daa[19*20+17] =  36.53690; 
	    daa[19*20+18] =  31.47300; 
	    	   
	    /* f[0]  = 0.08700; f[1]  = 0.04400; f[2]  = 0.03900; f[3]  = 0.05700;
	    f[4]  = 0.01900; f[5]  = 0.03700; f[6]  = 0.05800; f[7]  = 0.08300;
	    f[8]  = 0.02400; f[9]  = 0.04900; f[10] = 0.08600; f[11] = 0.06200;
	    f[12] = 0.02000; f[13] = 0.03800; f[14] = 0.04600; f[15] = 0.07000;
	    f[16] = 0.06100; f[17] = 0.01400; f[18] = 0.03500; f[19] = 0.07100;   
	    */

	    f[0] = 0.0866279; f[1] =  0.043972; f[2] =  0.0390894; f[3] =  0.0570451;
	    f[4] =  0.0193078; f[5] =  0.0367281; f[6] =  0.0580589; f[7] =  0.0832518;
	    f[8] =  0.0244313; f[9] =  0.048466; f[10] =  0.086209; f[11] = 0.0620286;
	    f[12] = 0.0195027; f[13] =  0.0384319; f[14] =  0.0457631; f[15] = 0.0695179;
	    f[16] =  0.0610127; f[17] =  0.0143859; f[18] =  0.0352742; f[19] =  0.0708957;

	  }
	  break;
	case RTREV:
	  {
	    daa[1*20+0]= 34;         daa[2*20+0]= 51;         daa[2*20+1]= 35;         daa[3*20+0]= 10;         
	    daa[3*20+1]= 30;         daa[3*20+2]= 384;        daa[4*20+0]= 439;        daa[4*20+1]= 92;         
	    daa[4*20+2]= 128;        daa[4*20+3]= 1;          daa[5*20+0]= 32;         daa[5*20+1]= 221;        
	    daa[5*20+2]= 236;        daa[5*20+3]= 78;         daa[5*20+4]= 70;         daa[6*20+0]= 81;         
	    daa[6*20+1]= 10;         daa[6*20+2]= 79;         daa[6*20+3]= 542;        daa[6*20+4]= 1;          
	    daa[6*20+5]= 372;        daa[7*20+0]= 135;        daa[7*20+1]= 41;         daa[7*20+2]= 94;         
	    daa[7*20+3]= 61;         daa[7*20+4]= 48;         daa[7*20+5]= 18;         daa[7*20+6]= 70;         
	    daa[8*20+0]= 30;         daa[8*20+1]= 90;         daa[8*20+2]= 320;        daa[8*20+3]= 91;         
	    daa[8*20+4]= 124;        daa[8*20+5]= 387;        daa[8*20+6]= 34;         daa[8*20+7]= 68;         
	    daa[9*20+0]= 1;          daa[9*20+1]= 24;         daa[9*20+2]= 35;         daa[9*20+3]= 1;          
	    daa[9*20+4]= 104;        daa[9*20+5]= 33;         daa[9*20+6]= 1;          daa[9*20+7]= 1;          
	    daa[9*20+8]= 34;         daa[10*20+0]= 45;        daa[10*20+1]= 18;        daa[10*20+2]= 15;        
	    daa[10*20+3]= 5;         daa[10*20+4]= 110;       daa[10*20+5]= 54;        daa[10*20+6]= 21;        
	    daa[10*20+7]= 3;         daa[10*20+8]= 51;        daa[10*20+9]= 385;       daa[11*20+0]= 38;        
	    daa[11*20+1]= 593;       daa[11*20+2]= 123;       daa[11*20+3]= 20;        daa[11*20+4]= 16;        
	    daa[11*20+5]= 309;       daa[11*20+6]= 141;       daa[11*20+7]= 30;        daa[11*20+8]= 76;        
	    daa[11*20+9]= 34;        daa[11*20+10]= 23;       daa[12*20+0]= 235;       daa[12*20+1]= 57;        
	    daa[12*20+2]= 1;         daa[12*20+3]= 1;         daa[12*20+4]= 156;       daa[12*20+5]= 158;       
	    daa[12*20+6]= 1;         daa[12*20+7]= 37;        daa[12*20+8]= 116;       daa[12*20+9]= 375;       
	    daa[12*20+10]= 581;      daa[12*20+11]= 134;      daa[13*20+0]= 1;         daa[13*20+1]= 7;         
	    daa[13*20+2]= 49;        daa[13*20+3]= 1;         daa[13*20+4]= 70;        daa[13*20+5]= 1;         
	    daa[13*20+6]= 1;         daa[13*20+7]= 7;         daa[13*20+8]= 141;       daa[13*20+9]= 64;        
	    daa[13*20+10]= 179;      daa[13*20+11]= 14;       daa[13*20+12]= 247;      daa[14*20+0]= 97;        
	    daa[14*20+1]= 24;        daa[14*20+2]= 33;        daa[14*20+3]= 55;        daa[14*20+4]= 1;         
	    daa[14*20+5]= 68;        daa[14*20+6]= 52;        daa[14*20+7]= 17;        daa[14*20+8]= 44;        
	    daa[14*20+9]= 10;        daa[14*20+10]= 22;       daa[14*20+11]= 43;       daa[14*20+12]= 1;        
	    daa[14*20+13]= 11;       daa[15*20+0]= 460;       daa[15*20+1]= 102;       daa[15*20+2]= 294;       
	    daa[15*20+3]= 136;       daa[15*20+4]= 75;        daa[15*20+5]= 225;       daa[15*20+6]= 95;        
	    daa[15*20+7]= 152;       daa[15*20+8]= 183;       daa[15*20+9]= 4;         daa[15*20+10]= 24;       
	    daa[15*20+11]= 77;       daa[15*20+12]= 1;        daa[15*20+13]= 20;       daa[15*20+14]= 134;      
	    daa[16*20+0]= 258;       daa[16*20+1]= 64;        daa[16*20+2]= 148;       daa[16*20+3]= 55;        
	    daa[16*20+4]= 117;       daa[16*20+5]= 146;       daa[16*20+6]= 82;        daa[16*20+7]= 7;         
	    daa[16*20+8]= 49;        daa[16*20+9]= 72;        daa[16*20+10]= 25;       daa[16*20+11]= 110;      
	    daa[16*20+12]= 131;      daa[16*20+13]= 69;       daa[16*20+14]= 62;       daa[16*20+15]= 671;      
	    daa[17*20+0]= 5;         daa[17*20+1]= 13;        daa[17*20+2]= 16;        daa[17*20+3]= 1;         
	    daa[17*20+4]= 55;        daa[17*20+5]= 10;        daa[17*20+6]= 17;        daa[17*20+7]= 23;        
	    daa[17*20+8]= 48;        daa[17*20+9]= 39;        daa[17*20+10]= 47;       daa[17*20+11]= 6;        
	    daa[17*20+12]= 111;      daa[17*20+13]= 182;      daa[17*20+14]= 9;        daa[17*20+15]= 14;       
	    daa[17*20+16]= 1;        daa[18*20+0]= 55;        daa[18*20+1]= 47;        daa[18*20+2]= 28;        
	    daa[18*20+3]= 1;         daa[18*20+4]= 131;       daa[18*20+5]= 45;        daa[18*20+6]= 1;         
	    daa[18*20+7]= 21;        daa[18*20+8]= 307;       daa[18*20+9]= 26;        daa[18*20+10]= 64;       
	    daa[18*20+11]= 1;        daa[18*20+12]= 74;       daa[18*20+13]= 1017;     daa[18*20+14]= 14;       
	    daa[18*20+15]= 31;       daa[18*20+16]= 34;       daa[18*20+17]= 176;      daa[19*20+0]= 197;       
	    daa[19*20+1]= 29;        daa[19*20+2]= 21;        daa[19*20+3]= 6;         daa[19*20+4]= 295;       
	    daa[19*20+5]= 36;        daa[19*20+6]= 35;        daa[19*20+7]= 3;         daa[19*20+8]= 1;         
	    daa[19*20+9]= 1048;      daa[19*20+10]= 112;      daa[19*20+11]= 19;       daa[19*20+12]= 236;      
	    daa[19*20+13]= 92;       daa[19*20+14]= 25;       daa[19*20+15]= 39;       daa[19*20+16]= 196;      
	    daa[19*20+17]= 26;       daa[19*20+18]= 59;       
	    
	    f[0]= 0.0646;           f[1]= 0.0453;           f[2]= 0.0376;           f[3]= 0.0422;           
	    f[4]= 0.0114;           f[5]= 0.0606;           f[6]= 0.0607;           f[7]= 0.0639;           
	    f[8]= 0.0273;           f[9]= 0.0679;           f[10]= 0.1018;          f[11]= 0.0751;          
	    f[12]= 0.015;           f[13]= 0.0287;          f[14]= 0.0681;          f[15]= 0.0488;          
	    f[16]= 0.0622;          f[17]= 0.0251;          f[18]= 0.0318;          f[19]= 0.0619;	    	    
	  }
	  break;
	case CPREV:
	  {
	    daa[1*20+0]= 105;        daa[2*20+0]= 227;        daa[2*20+1]= 357;        daa[3*20+0]= 175;        
	    daa[3*20+1]= 43;         daa[3*20+2]= 4435;       daa[4*20+0]= 669;        daa[4*20+1]= 823;        
	    daa[4*20+2]= 538;        daa[4*20+3]= 10;         daa[5*20+0]= 157;        daa[5*20+1]= 1745;       
	    daa[5*20+2]= 768;        daa[5*20+3]= 400;        daa[5*20+4]= 10;         daa[6*20+0]= 499;        
	    daa[6*20+1]= 152;        daa[6*20+2]= 1055;       daa[6*20+3]= 3691;       daa[6*20+4]= 10;         
	    daa[6*20+5]= 3122;       daa[7*20+0]= 665;        daa[7*20+1]= 243;        daa[7*20+2]= 653;        
	    daa[7*20+3]= 431;        daa[7*20+4]= 303;        daa[7*20+5]= 133;        daa[7*20+6]= 379;        
	    daa[8*20+0]= 66;         daa[8*20+1]= 715;        daa[8*20+2]= 1405;       daa[8*20+3]= 331;        
	    daa[8*20+4]= 441;        daa[8*20+5]= 1269;       daa[8*20+6]= 162;        daa[8*20+7]= 19;         
	    daa[9*20+0]= 145;        daa[9*20+1]= 136;        daa[9*20+2]= 168;        daa[9*20+3]= 10;         
	    daa[9*20+4]= 280;        daa[9*20+5]= 92;         daa[9*20+6]= 148;        daa[9*20+7]= 40;         
	    daa[9*20+8]= 29;         daa[10*20+0]= 197;       daa[10*20+1]= 203;       daa[10*20+2]= 113;       
	    daa[10*20+3]= 10;        daa[10*20+4]= 396;       daa[10*20+5]= 286;       daa[10*20+6]= 82;        
	    daa[10*20+7]= 20;        daa[10*20+8]= 66;        daa[10*20+9]= 1745;      daa[11*20+0]= 236;       
	    daa[11*20+1]= 4482;      daa[11*20+2]= 2430;      daa[11*20+3]= 412;       daa[11*20+4]= 48;        
	    daa[11*20+5]= 3313;      daa[11*20+6]= 2629;      daa[11*20+7]= 263;       daa[11*20+8]= 305;       
	    daa[11*20+9]= 345;       daa[11*20+10]= 218;      daa[12*20+0]= 185;       daa[12*20+1]= 125;       
	    daa[12*20+2]= 61;        daa[12*20+3]= 47;        daa[12*20+4]= 159;       daa[12*20+5]= 202;       
	    daa[12*20+6]= 113;       daa[12*20+7]= 21;        daa[12*20+8]= 10;        daa[12*20+9]= 1772;      
	    daa[12*20+10]= 1351;     daa[12*20+11]= 193;      daa[13*20+0]= 68;        daa[13*20+1]= 53;        
	    daa[13*20+2]= 97;        daa[13*20+3]= 22;        daa[13*20+4]= 726;       daa[13*20+5]= 10;        
	    daa[13*20+6]= 145;       daa[13*20+7]= 25;        daa[13*20+8]= 127;       daa[13*20+9]= 454;       
	    daa[13*20+10]= 1268;     daa[13*20+11]= 72;       daa[13*20+12]= 327;      daa[14*20+0]= 490;       
	    daa[14*20+1]= 87;        daa[14*20+2]= 173;       daa[14*20+3]= 170;       daa[14*20+4]= 285;       
	    daa[14*20+5]= 323;       daa[14*20+6]= 185;       daa[14*20+7]= 28;        daa[14*20+8]= 152;       
	    daa[14*20+9]= 117;       daa[14*20+10]= 219;      daa[14*20+11]= 302;      daa[14*20+12]= 100;      
	    daa[14*20+13]= 43;       daa[15*20+0]= 2440;      daa[15*20+1]= 385;       daa[15*20+2]= 2085;      
	    daa[15*20+3]= 590;       daa[15*20+4]= 2331;      daa[15*20+5]= 396;       daa[15*20+6]= 568;       
	    daa[15*20+7]= 691;       daa[15*20+8]= 303;       daa[15*20+9]= 216;       daa[15*20+10]= 516;      
	    daa[15*20+11]= 868;      daa[15*20+12]= 93;       daa[15*20+13]= 487;      daa[15*20+14]= 1202;     
	    daa[16*20+0]= 1340;      daa[16*20+1]= 314;       daa[16*20+2]= 1393;      daa[16*20+3]= 266;       
	    daa[16*20+4]= 576;       daa[16*20+5]= 241;       daa[16*20+6]= 369;       daa[16*20+7]= 92;        
	    daa[16*20+8]= 32;        daa[16*20+9]= 1040;      daa[16*20+10]= 156;      daa[16*20+11]= 918;      
	    daa[16*20+12]= 645;      daa[16*20+13]= 148;      daa[16*20+14]= 260;      daa[16*20+15]= 2151;     
	    daa[17*20+0]= 14;        daa[17*20+1]= 230;       daa[17*20+2]= 40;        daa[17*20+3]= 18;        
	    daa[17*20+4]= 435;       daa[17*20+5]= 53;        daa[17*20+6]= 63;        daa[17*20+7]= 82;        
	    daa[17*20+8]= 69;        daa[17*20+9]= 42;        daa[17*20+10]= 159;      daa[17*20+11]= 10;       
	    daa[17*20+12]= 86;       daa[17*20+13]= 468;      daa[17*20+14]= 49;       daa[17*20+15]= 73;       
	    daa[17*20+16]= 29;       daa[18*20+0]= 56;        daa[18*20+1]= 323;       daa[18*20+2]= 754;       
	    daa[18*20+3]= 281;       daa[18*20+4]= 1466;      daa[18*20+5]= 391;       daa[18*20+6]= 142;       
	    daa[18*20+7]= 10;        daa[18*20+8]= 1971;      daa[18*20+9]= 89;        daa[18*20+10]= 189;      
	    daa[18*20+11]= 247;      daa[18*20+12]= 215;      daa[18*20+13]= 2370;     daa[18*20+14]= 97;       
	    daa[18*20+15]= 522;      daa[18*20+16]= 71;       daa[18*20+17]= 346;      daa[19*20+0]= 968;       
	    daa[19*20+1]= 92;        daa[19*20+2]= 83;        daa[19*20+3]= 75;        daa[19*20+4]= 592;       
	    daa[19*20+5]= 54;        daa[19*20+6]= 200;       daa[19*20+7]= 91;        daa[19*20+8]= 25;        
	    daa[19*20+9]= 4797;      daa[19*20+10]= 865;      daa[19*20+11]= 249;      daa[19*20+12]= 475;      
	    daa[19*20+13]= 317;      daa[19*20+14]= 122;      daa[19*20+15]= 167;      daa[19*20+16]= 760;      
	    daa[19*20+17]= 10;       daa[19*20+18]= 119;      
	    
	    f[0]= 0.076;            f[1]= 0.062;            f[2]= 0.041;            f[3]= 0.037;            
	    f[4]= 0.009;            f[5]= 0.038;            f[6]= 0.049;            f[7]= 0.084;            
	    f[8]= 0.025;            f[9]= 0.081;            f[10]= 0.101;           f[11]= 0.05;            
	    f[12]= 0.022;           f[13]= 0.051;           f[14]= 0.043;           f[15]= 0.062;           
	    f[16]= 0.054;           f[17]= 0.018;           f[18]= 0.031;           f[19]= 0.066; 
	  }
	  break;
	case VT:
	  {
	    /*
	      daa[1*20+0]= 0.233108;   daa[2*20+0]= 0.199097;   daa[2*20+1]= 0.210797;   daa[3*20+0]= 0.265145;   
	      daa[3*20+1]= 0.105191;   daa[3*20+2]= 0.883422;   daa[4*20+0]= 0.227333;   daa[4*20+1]= 0.031726;   
	      daa[4*20+2]= 0.027495;   daa[4*20+3]= 0.010313;   daa[5*20+0]= 0.310084;   daa[5*20+1]= 0.493763;   
	      daa[5*20+2]= 0.2757;     daa[5*20+3]= 0.205842;   daa[5*20+4]= 0.004315;   daa[6*20+0]= 0.567957;   
	      daa[6*20+1]= 0.25524;    daa[6*20+2]= 0.270417;   daa[6*20+3]= 1.599461;   daa[6*20+4]= 0.005321;   
	      daa[6*20+5]= 0.960976;   daa[7*20+0]= 0.876213;   daa[7*20+1]= 0.156945;   daa[7*20+2]= 0.362028;   
	      daa[7*20+3]= 0.311718;   daa[7*20+4]= 0.050876;   daa[7*20+5]= 0.12866;    daa[7*20+6]= 0.250447;   
	      daa[8*20+0]= 0.078692;   daa[8*20+1]= 0.213164;   daa[8*20+2]= 0.290006;   daa[8*20+3]= 0.134252;   
	      daa[8*20+4]= 0.016695;   daa[8*20+5]= 0.315521;   daa[8*20+6]= 0.104458;   daa[8*20+7]= 0.058131;   
	      daa[9*20+0]= 0.222972;   daa[9*20+1]= 0.08151;    daa[9*20+2]= 0.087225;   daa[9*20+3]= 0.01172;    
	      daa[9*20+4]= 0.046398;   daa[9*20+5]= 0.054602;   daa[9*20+6]= 0.046589;   daa[9*20+7]= 0.051089;   
	      daa[9*20+8]= 0.020039;   daa[10*20+0]= 0.42463;   daa[10*20+1]= 0.192364;  daa[10*20+2]= 0.069245;  
	      daa[10*20+3]= 0.060863;  daa[10*20+4]= 0.091709;  daa[10*20+5]= 0.24353;   daa[10*20+6]= 0.151924;  
	      daa[10*20+7]= 0.087056;  daa[10*20+8]= 0.103552;  daa[10*20+9]= 2.08989;   daa[11*20+0]= 0.393245;  
	      daa[11*20+1]= 1.755838;  daa[11*20+2]= 0.50306;   daa[11*20+3]= 0.261101;  daa[11*20+4]= 0.004067;  
	      daa[11*20+5]= 0.738208;  daa[11*20+6]= 0.88863;   daa[11*20+7]= 0.193243;  daa[11*20+8]= 0.153323;  
	      daa[11*20+9]= 0.093181;  daa[11*20+10]= 0.201204; daa[12*20+0]= 0.21155;   daa[12*20+1]= 0.08793;   
	      daa[12*20+2]= 0.05742;   daa[12*20+3]= 0.012182;  daa[12*20+4]= 0.02369;   daa[12*20+5]= 0.120801;  
	      daa[12*20+6]= 0.058643;  daa[12*20+7]= 0.04656;   daa[12*20+8]= 0.021157;  daa[12*20+9]= 0.493845;  
	      daa[12*20+10]= 1.105667; daa[12*20+11]= 0.096474; daa[13*20+0]= 0.116646;  daa[13*20+1]= 0.042569;  
	      daa[13*20+2]= 0.039769;  daa[13*20+3]= 0.016577;  daa[13*20+4]= 0.051127;  daa[13*20+5]= 0.026235;  
	      daa[13*20+6]= 0.028168;  daa[13*20+7]= 0.050143;  daa[13*20+8]= 0.079807;  daa[13*20+9]= 0.32102;   
	      daa[13*20+10]= 0.946499; daa[13*20+11]= 0.038261; daa[13*20+12]= 0.173052; daa[14*20+0]= 0.399143;  
	      daa[14*20+1]= 0.12848;   daa[14*20+2]= 0.083956;  daa[14*20+3]= 0.160063;  daa[14*20+4]= 0.011137;  
	      daa[14*20+5]= 0.15657;   daa[14*20+6]= 0.205134;  daa[14*20+7]= 0.124492;  daa[14*20+8]= 0.078892;  
	      daa[14*20+9]= 0.054797;  daa[14*20+10]= 0.169784; daa[14*20+11]= 0.212302; daa[14*20+12]= 0.010363; 
	      daa[14*20+13]= 0.042564; daa[15*20+0]= 1.817198;  daa[15*20+1]= 0.292327;  daa[15*20+2]= 0.847049;  
	      daa[15*20+3]= 0.461519;  daa[15*20+4]= 0.17527;   daa[15*20+5]= 0.358017;  daa[15*20+6]= 0.406035;  
	      daa[15*20+7]= 0.612843;  daa[15*20+8]= 0.167406;  daa[15*20+9]= 0.081567;  daa[15*20+10]= 0.214977; 
	      daa[15*20+11]= 0.400072; daa[15*20+12]= 0.090515; daa[15*20+13]= 0.138119; daa[15*20+14]= 0.430431; 
	      daa[16*20+0]= 0.877877;  daa[16*20+1]= 0.204109;  daa[16*20+2]= 0.471268;  daa[16*20+3]= 0.178197;  
	      daa[16*20+4]= 0.079511;  daa[16*20+5]= 0.248992;  daa[16*20+6]= 0.321028;  daa[16*20+7]= 0.136266;  
	      daa[16*20+8]= 0.101117;  daa[16*20+9]= 0.376588;  daa[16*20+10]= 0.243227; daa[16*20+11]= 0.446646; 
	      daa[16*20+12]= 0.184609; daa[16*20+13]= 0.08587;  daa[16*20+14]= 0.207143; daa[16*20+15]= 1.767766; 
	      daa[17*20+0]= 0.030309;  daa[17*20+1]= 0.046417;  daa[17*20+2]= 0.010459;  daa[17*20+3]= 0.011393;  
	      daa[17*20+4]= 0.007732;  daa[17*20+5]= 0.021248;  daa[17*20+6]= 0.018844;  daa[17*20+7]= 0.02399;   
	      daa[17*20+8]= 0.020009;  daa[17*20+9]= 0.034954;  daa[17*20+10]= 0.083439; daa[17*20+11]= 0.023321; 
	      daa[17*20+12]= 0.022019; daa[17*20+13]= 0.12805;  daa[17*20+14]= 0.014584; daa[17*20+15]= 0.035933; 
	      daa[17*20+16]= 0.020437; daa[18*20+0]= 0.087061;  daa[18*20+1]= 0.09701;   daa[18*20+2]= 0.093268;  
	      daa[18*20+3]= 0.051664;  daa[18*20+4]= 0.042823;  daa[18*20+5]= 0.062544;  daa[18*20+6]= 0.0552;    
	      daa[18*20+7]= 0.037568;  daa[18*20+8]= 0.286027;  daa[18*20+9]= 0.086237;  daa[18*20+10]= 0.189842; 
	      daa[18*20+11]= 0.068689; daa[18*20+12]= 0.073223; daa[18*20+13]= 0.898663; daa[18*20+14]= 0.032043; 
	      daa[18*20+15]= 0.121979; daa[18*20+16]= 0.094617; daa[18*20+17]= 0.124746; daa[19*20+0]= 1.230985;  
	      daa[19*20+1]= 0.113146;  daa[19*20+2]= 0.049824;  daa[19*20+3]= 0.048769;  daa[19*20+4]= 0.163831;  
	      daa[19*20+5]= 0.112027;  daa[19*20+6]= 0.205868;  daa[19*20+7]= 0.082579;  daa[19*20+8]= 0.068575;  
	      daa[19*20+9]= 3.65443;   daa[19*20+10]= 1.337571; daa[19*20+11]= 0.144587; daa[19*20+12]= 0.307309; 
	      daa[19*20+13]= 0.247329; daa[19*20+14]= 0.129315; daa[19*20+15]= 0.1277;   daa[19*20+16]= 0.740372; 
	      daa[19*20+17]= 0.022134; daa[19*20+18]= 0.125733; 	    	    
	      
	      f[0]  = 0.07900;         f[1]= 0.05100;        f[2]  = 0.04200;         f[3]= 0.05300;         
	      f[4]  = 0.01500;         f[5]= 0.03700;        f[6]  = 0.06200;         f[7]= 0.07100;         
	      f[8]  = 0.02300;         f[9]= 0.06200;        f[10] = 0.09600;        f[11]= 0.05700;        
	      f[12] = 0.02400;        f[13]= 0.04300;        f[14] = 0.04400;        f[15]= 0.06400;        
	      f[16] = 0.05600;        f[17]= 0.01300;        f[18] = 0.03500;        f[19]= 0.07300; 
	    */

	    daa[1*20+0]=   1.2412691067876198;
	    daa[2*20+0]=   1.2184237953498958;
	    daa[2*20+1]=   1.5720770753326880;
	    daa[3*20+0]=   1.3759368509441177;
	    daa[3*20+1]=   0.7550654439001206;
	    daa[3*20+2]=   7.8584219153689405;
	    daa[4*20+0]=   2.4731223087544874;
	    daa[4*20+1]=   1.4414262567428417;
	    daa[4*20+2]=   0.9784679122774127;
	    daa[4*20+3]=   0.2272488448121475;
	    daa[5*20+0]=   2.2155167805137470;
	    daa[5*20+1]=   5.5120819705248678;
	    daa[5*20+2]=   3.0143201670924822;
	    daa[5*20+3]=   1.6562495638176040;
	    daa[5*20+4]=   0.4587469126746136;
	    daa[6*20+0]=   2.3379911207495061;
	    daa[6*20+1]=   1.3542404860613146;
	    daa[6*20+2]=   2.0093434778398112;
	    daa[6*20+3]=   9.6883451875685065;
	    daa[6*20+4]=   0.4519167943192672;
	    daa[6*20+5]=   6.8124601839937675;
	    daa[7*20+0]=   3.3386555146457697;
	    daa[7*20+1]=   1.3121700301622004;
	    daa[7*20+2]=   2.4117632898861809;
	    daa[7*20+3]=   1.9142079025990228;
	    daa[7*20+4]=   1.1034605684472507;
	    daa[7*20+5]=   0.8776110594765502;
	    daa[7*20+6]=   1.3860121390169038;
	    daa[8*20+0]=   0.9615841926910841;
	    daa[8*20+1]=   4.9238668283945266;
	    daa[8*20+2]=   6.1974384977884114;
	    daa[8*20+3]=   2.1459640610133781;
	    daa[8*20+4]=   1.5196756759380692;
	    daa[8*20+5]=   7.9943228564946525;
	    daa[8*20+6]=   1.6360079688522375;
	    daa[8*20+7]=   0.8561248973045037;
	    daa[9*20+0]=   0.8908203061925510;
	    daa[9*20+1]=   0.4323005487925516;
	    daa[9*20+2]=   0.9179291175331520;
	    daa[9*20+3]=   0.2161660372725585;
	    daa[9*20+4]=   0.9126668032539315;
	    daa[9*20+5]=   0.4882733432879921;
	    daa[9*20+6]=   0.4035497929633328;
	    daa[9*20+7]=   0.2888075033037488;
	    daa[9*20+8]=   0.5787937115407940;
	    daa[10*20+0]=  1.0778497408764076;
	    daa[10*20+1]=  0.8386701149158265;
	    daa[10*20+2]=  0.4098311270816011;
	    daa[10*20+3]=  0.3574207468998517;
	    daa[10*20+4]=  1.4081315998413697;
	    daa[10*20+5]=  1.3318097154194044;
	    daa[10*20+6]=  0.5610717242294755;
	    daa[10*20+7]=  0.3578662395745526;
	    daa[10*20+8]=  1.0765007949562073;
	    daa[10*20+9]=  6.0019110258426362;
	    daa[11*20+0]=  1.4932055816372476;
	    daa[11*20+1]=  10.017330817366002;
	    daa[11*20+2]=  4.4034547578962568;
	    daa[11*20+3]=  1.4521790561663968;
	    daa[11*20+4]=  0.3371091785647479;
	    daa[11*20+5]=  6.0519085243118811;
	    daa[11*20+6]=  4.3290086529582830;
	    daa[11*20+7]=  0.8945563662345198;
	    daa[11*20+8]=  1.8085136096039203;
	    daa[11*20+9]=  0.6244297525127139;
	    daa[11*20+10]= 0.5642322882556321;
	    daa[12*20+0]=  1.9006455961717605;
	    daa[12*20+1]=  1.2488638689609959;
	    daa[12*20+2]=  0.9378803706165143;
	    daa[12*20+3]=  0.4075239926000898;
	    daa[12*20+4]=  1.2213054800811556;
	    daa[12*20+5]=  1.9106190827629084;
	    daa[12*20+6]=  0.7471936218068498;
	    daa[12*20+7]=  0.5954812791740037;
	    daa[12*20+8]=  1.3808291710019667;
	    daa[12*20+9]=  6.7597899772045418;
	    daa[12*20+10]= 8.0327792947421148;
	    daa[12*20+11]= 1.7129670976916258;
	    daa[13*20+0]=  0.6883439026872615;
	    daa[13*20+1]=  0.4224945197276290;
	    daa[13*20+2]=  0.5044944273324311;
	    daa[13*20+3]=  0.1675129724559251;
	    daa[13*20+4]=  1.6953951980808002;
	    daa[13*20+5]=  0.3573432522499545;
	    daa[13*20+6]=  0.2317194387691585;
	    daa[13*20+7]=  0.3693722640980460;
	    daa[13*20+8]=  1.3629765501081097;
	    daa[13*20+9]=  2.2864286949316077;
	    daa[13*20+10]= 4.3611548063555778;
	    daa[13*20+11]= 0.3910559903834828;
	    daa[13*20+12]= 2.3201373546296349;
	    daa[14*20+0]=  2.7355620089953550;
	    daa[14*20+1]=  1.3091837782420783;
	    daa[14*20+2]=  0.7103720531974738;
	    daa[14*20+3]=  1.0714605979577547;
	    daa[14*20+4]=  0.4326227078645523;
	    daa[14*20+5]=  2.3019177728300728;
	    daa[14*20+6]=  1.5132807416252063;
	    daa[14*20+7]=  0.7744933618134962;
	    daa[14*20+8]=  1.8370555852070649;
	    daa[14*20+9]=  0.4811402387911145;
	    daa[14*20+10]= 1.0084320519837335;
	    daa[14*20+11]= 1.3918935593582853;
	    daa[14*20+12]= 0.4953193808676289;
	    daa[14*20+13]= 0.3746821107962129;
	    daa[15*20+0]=  6.4208961859142883;
	    daa[15*20+1]=  1.9202994262316166;
	    daa[15*20+2]=  6.1234512396801764;
	    daa[15*20+3]=  2.2161944596741829;
	    daa[15*20+4]=  3.6366815408744255;
	    daa[15*20+5]=  2.3193703643237220;
	    daa[15*20+6]=  1.8273535587773553;
	    daa[15*20+7]=  3.0637776193717610;
	    daa[15*20+8]=  1.9699895187387506;
	    daa[15*20+9]=  0.6047491507504744;
	    daa[15*20+10]= 0.8953754669269811;
	    daa[15*20+11]= 1.9776630140912268;
	    daa[15*20+12]= 1.0657482318076852;
	    daa[15*20+13]= 1.1079144700606407;
	    daa[15*20+14]= 3.5465914843628927;
	    daa[16*20+0]=  5.2892514169776437;
	    daa[16*20+1]=  1.3363401740560601;
	    daa[16*20+2]=  3.8852506105922231;
	    daa[16*20+3]=  1.5066839872944762;
	    daa[16*20+4]=  1.7557065205837685;
	    daa[16*20+5]=  2.1576510103471440;
	    daa[16*20+6]=  1.5839981708584689;
	    daa[16*20+7]=  0.7147489676267383;
	    daa[16*20+8]=  1.6136654573285647;
	    daa[16*20+9]=  2.6344778384442731;
	    daa[16*20+10]= 1.0192004372506540;
	    daa[16*20+11]= 2.5513781312660280;
	    daa[16*20+12]= 3.3628488360462363;
	    daa[16*20+13]= 0.6882725908872254;
	    daa[16*20+14]= 1.9485376673137556;
	    daa[16*20+15]= 8.8479984061248178;
	    daa[17*20+0]=  0.5488578478106930;
	    daa[17*20+1]=  1.5170142153962840;
	    daa[17*20+2]=  0.1808525752605976;
	    daa[17*20+3]=  0.2496584188151770;
	    daa[17*20+4]=  1.6275179891253113;
	    daa[17*20+5]=  0.8959082681546182;
	    daa[17*20+6]=  0.4198391148111098;
	    daa[17*20+7]=  0.9349753595598769;
	    daa[17*20+8]=  0.6301954684360302;
	    daa[17*20+9]=  0.5604648274060783;
	    daa[17*20+10]= 1.5183114434679339;
	    daa[17*20+11]= 0.5851920879490173;
	    daa[17*20+12]= 1.4680478689711018;
	    daa[17*20+13]= 3.3448437239772266;
	    daa[17*20+14]= 0.4326058001438786;
	    daa[17*20+15]= 0.6791126595939816;
	    daa[17*20+16]= 0.4514203099376473;
	    daa[18*20+0]=  0.5411769916657778;
	    daa[18*20+1]=  0.8912614404565405;
	    daa[18*20+2]=  1.0894926581511342;
	    daa[18*20+3]=  0.7447620891784513;
	    daa[18*20+4]=  2.1579775140421025;
	    daa[18*20+5]=  0.9183596801412757;
	    daa[18*20+6]=  0.5818111331782764;
	    daa[18*20+7]=  0.3374467649724478;
	    daa[18*20+8]=  7.7587442309146040;
	    daa[18*20+9]=  0.8626796044156272;
	    daa[18*20+10]= 1.2452243224541324;
	    daa[18*20+11]= 0.7835447533710449;
	    daa[18*20+12]= 1.0899165770956820;
	    daa[18*20+13]= 10.384852333133459;
	    daa[18*20+14]= 0.4819109019647465;
	    daa[18*20+15]= 0.9547229305958682;
	    daa[18*20+16]= 0.8564314184691215;
	    daa[18*20+17]= 4.5377235790405388;
	    daa[19*20+0]=  4.6501894691803214;
	    daa[19*20+1]=  0.7807017855806767;
	    daa[19*20+2]=  0.4586061981719967;
	    daa[19*20+3]=  0.4594535241660911;
	    daa[19*20+4]=  2.2627456996290891;
	    daa[19*20+5]=  0.6366932501396869;
	    daa[19*20+6]=  0.8940572875547330;
	    daa[19*20+7]=  0.6193321034173915;
	    daa[19*20+8]=  0.5333220944030346;
	    daa[19*20+9]=  14.872933461519061;
	    daa[19*20+10]= 3.5458093276667237;
	    daa[19*20+11]= 0.7801080335991272;
	    daa[19*20+12]= 4.0584577156753401;
	    daa[19*20+13]= 1.7039730522675411;
	    daa[19*20+14]= 0.5985498912985666;
	    daa[19*20+15]= 0.9305232113028208;
	    daa[19*20+16]= 3.4242218450865543;
	    daa[19*20+17]= 0.5658969249032649;
	    daa[19*20+18]= 1.0000000000000000;
	    
	    f[0]=  0.0770764620135024;
	    f[1]=  0.0500819370772208;
	    f[2]=  0.0462377395993731;
	    f[3]=  0.0537929860758246;
	    f[4]=  0.0144533387583345;
	    f[5]=  0.0408923608974345;
	    f[6]=  0.0633579339160905;
	    f[7]=  0.0655672355884439;
	    f[8]=  0.0218802687005936;
	    f[9]=  0.0591969699027449;
	    f[10]= 0.0976461276528445;
	    f[11]= 0.0592079410822730;
	    f[12]= 0.0220695876653368;
	    f[13]= 0.0413508521834260;
	    f[14]= 0.0476871596856874;
	    f[15]= 0.0707295165111524;
	    f[16]= 0.0567759161524817;
	    f[17]= 0.0127019797647213;
	    f[18]= 0.0323746050281867;
	    f[19]= 0.0669190817443274;
	  }
	  break;
	case BLOSUM62:
	  {
	    daa[1*20+0]= 0.735790389698;  daa[2*20+0]= 0.485391055466;  daa[2*20+1]= 1.297446705134;  
	    daa[3*20+0]= 0.543161820899;  
	    daa[3*20+1]= 0.500964408555;  daa[3*20+2]= 3.180100048216;  daa[4*20+0]= 1.45999531047;   
	    daa[4*20+1]= 0.227826574209;  
	    daa[4*20+2]= 0.397358949897;  daa[4*20+3]= 0.240836614802;  daa[5*20+0]= 1.199705704602;  
	    daa[5*20+1]= 3.020833610064;  
	    daa[5*20+2]= 1.839216146992;  daa[5*20+3]= 1.190945703396;  daa[5*20+4]= 0.32980150463;   
	    daa[6*20+0]= 1.1709490428;    
	    daa[6*20+1]= 1.36057419042;   daa[6*20+2]= 1.24048850864;   daa[6*20+3]= 3.761625208368;  
	    daa[6*20+4]= 0.140748891814;  
	    daa[6*20+5]= 5.528919177928;  daa[7*20+0]= 1.95588357496;   daa[7*20+1]= 0.418763308518;  
	    daa[7*20+2]= 1.355872344485;  
	    daa[7*20+3]= 0.798473248968;  daa[7*20+4]= 0.418203192284;  daa[7*20+5]= 0.609846305383;  
	    daa[7*20+6]= 0.423579992176;  
	    daa[8*20+0]= 0.716241444998;  daa[8*20+1]= 1.456141166336;  daa[8*20+2]= 2.414501434208;  
	    daa[8*20+3]= 0.778142664022;  
	    daa[8*20+4]= 0.354058109831;  daa[8*20+5]= 2.43534113114;   daa[8*20+6]= 1.626891056982;  
	    daa[8*20+7]= 0.539859124954;  
	    daa[9*20+0]= 0.605899003687;  daa[9*20+1]= 0.232036445142;  daa[9*20+2]= 0.283017326278;  
	    daa[9*20+3]= 0.418555732462;  
	    daa[9*20+4]= 0.774894022794;  daa[9*20+5]= 0.236202451204;  daa[9*20+6]= 0.186848046932;  
	    daa[9*20+7]= 0.189296292376;  
	    daa[9*20+8]= 0.252718447885;  daa[10*20+0]= 0.800016530518; daa[10*20+1]= 0.622711669692; 
	    daa[10*20+2]= 0.211888159615; 
	    daa[10*20+3]= 0.218131577594; daa[10*20+4]= 0.831842640142; daa[10*20+5]= 0.580737093181; 
	    daa[10*20+6]= 0.372625175087; 
	    daa[10*20+7]= 0.217721159236; daa[10*20+8]= 0.348072209797; daa[10*20+9]= 3.890963773304; 
	    daa[11*20+0]= 1.295201266783; 
	    daa[11*20+1]= 5.411115141489; daa[11*20+2]= 1.593137043457; daa[11*20+3]= 1.032447924952; 
	    daa[11*20+4]= 0.285078800906; 
	    daa[11*20+5]= 3.945277674515; daa[11*20+6]= 2.802427151679; daa[11*20+7]= 0.752042440303; 
	    daa[11*20+8]= 1.022507035889; 
	    daa[11*20+9]= 0.406193586642; daa[11*20+10]= 0.445570274261;daa[12*20+0]= 1.253758266664; 
	    daa[12*20+1]= 0.983692987457; 
	    daa[12*20+2]= 0.648441278787; daa[12*20+3]= 0.222621897958; daa[12*20+4]= 0.76768882348;  
	    daa[12*20+5]= 2.494896077113; 
	    daa[12*20+6]= 0.55541539747;  daa[12*20+7]= 0.459436173579; daa[12*20+8]= 0.984311525359; 
	    daa[12*20+9]= 3.364797763104; 
	    daa[12*20+10]= 6.030559379572;daa[12*20+11]= 1.073061184332;daa[13*20+0]= 0.492964679748; 
	    daa[13*20+1]= 0.371644693209; 
	    daa[13*20+2]= 0.354861249223; daa[13*20+3]= 0.281730694207; daa[13*20+4]= 0.441337471187; 
	    daa[13*20+5]= 0.14435695975;  
	    daa[13*20+6]= 0.291409084165; daa[13*20+7]= 0.368166464453; daa[13*20+8]= 0.714533703928; 
	    daa[13*20+9]= 1.517359325954; 
	    daa[13*20+10]= 2.064839703237;daa[13*20+11]= 0.266924750511;daa[13*20+12]= 1.77385516883; 
	    daa[14*20+0]= 1.173275900924; 
	    daa[14*20+1]= 0.448133661718; daa[14*20+2]= 0.494887043702; daa[14*20+3]= 0.730628272998; 
	    daa[14*20+4]= 0.356008498769; 
	    daa[14*20+5]= 0.858570575674; daa[14*20+6]= 0.926563934846; daa[14*20+7]= 0.504086599527; daa[14*20+8]= 0.527007339151; 
	    daa[14*20+9]= 0.388355409206; daa[14*20+10]= 0.374555687471;daa[14*20+11]= 1.047383450722;daa[14*20+12]= 0.454123625103;
	    daa[14*20+13]= 0.233597909629;daa[15*20+0]= 4.325092687057; daa[15*20+1]= 1.12278310421;  daa[15*20+2]= 2.904101656456; 
	    daa[15*20+3]= 1.582754142065; daa[15*20+4]= 1.197188415094; daa[15*20+5]= 1.934870924596; daa[15*20+6]= 1.769893238937; 
	    daa[15*20+7]= 1.509326253224; daa[15*20+8]= 1.11702976291;  daa[15*20+9]= 0.35754441246;  daa[15*20+10]= 0.352969184527;
	    daa[15*20+11]= 1.752165917819;daa[15*20+12]= 0.918723415746;daa[15*20+13]= 0.540027644824;daa[15*20+14]= 1.169129577716;
	    daa[16*20+0]= 1.729178019485; daa[16*20+1]= 0.914665954563; daa[16*20+2]= 1.898173634533; daa[16*20+3]= 0.934187509431; 
	    daa[16*20+4]= 1.119831358516; daa[16*20+5]= 1.277480294596; daa[16*20+6]= 1.071097236007; daa[16*20+7]= 0.641436011405; 
	    daa[16*20+8]= 0.585407090225; daa[16*20+9]= 1.17909119726;  daa[16*20+10]= 0.915259857694;daa[16*20+11]= 1.303875200799;
	    daa[16*20+12]= 1.488548053722;daa[16*20+13]= 0.488206118793;daa[16*20+14]= 1.005451683149;daa[16*20+15]= 5.15155629227; 
	    daa[17*20+0]= 0.465839367725; daa[17*20+1]= 0.426382310122; daa[17*20+2]= 0.191482046247; daa[17*20+3]= 0.145345046279; 
	    daa[17*20+4]= 0.527664418872; daa[17*20+5]= 0.758653808642; daa[17*20+6]= 0.407635648938; daa[17*20+7]= 0.508358924638; 
	    daa[17*20+8]= 0.30124860078;  daa[17*20+9]= 0.34198578754;  daa[17*20+10]= 0.6914746346;  daa[17*20+11]= 0.332243040634;
	    daa[17*20+12]= 0.888101098152;daa[17*20+13]= 2.074324893497;daa[17*20+14]= 0.252214830027;daa[17*20+15]= 0.387925622098;
	    daa[17*20+16]= 0.513128126891;daa[18*20+0]= 0.718206697586; daa[18*20+1]= 0.720517441216; daa[18*20+2]= 0.538222519037; 
	    daa[18*20+3]= 0.261422208965; daa[18*20+4]= 0.470237733696; daa[18*20+5]= 0.95898974285;  daa[18*20+6]= 0.596719300346; 
	    daa[18*20+7]= 0.308055737035; daa[18*20+8]= 4.218953969389; daa[18*20+9]= 0.674617093228; daa[18*20+10]= 0.811245856323;
	    daa[18*20+11]= 0.7179934869;  daa[18*20+12]= 0.951682162246;daa[18*20+13]= 6.747260430801;daa[18*20+14]= 0.369405319355;
	    daa[18*20+15]= 0.796751520761;daa[18*20+16]= 0.801010243199;daa[18*20+17]= 4.054419006558;daa[19*20+0]= 2.187774522005; 
	    daa[19*20+1]= 0.438388343772; daa[19*20+2]= 0.312858797993; daa[19*20+3]= 0.258129289418; daa[19*20+4]= 1.116352478606; 
	    daa[19*20+5]= 0.530785790125; daa[19*20+6]= 0.524253846338; daa[19*20+7]= 0.25334079019;  daa[19*20+8]= 0.20155597175;  
	    daa[19*20+9]= 8.311839405458; daa[19*20+10]= 2.231405688913;daa[19*20+11]= 0.498138475304;daa[19*20+12]= 2.575850755315;
	    daa[19*20+13]= 0.838119610178;daa[19*20+14]= 0.496908410676;daa[19*20+15]= 0.561925457442;daa[19*20+16]= 2.253074051176;
	    daa[19*20+17]= 0.266508731426;daa[19*20+18]= 1;             
	    
	    f[0]= 0.074;                 f[1]= 0.052;                 f[2]= 0.045;                 f[3]= 0.054;                 
	    f[4]= 0.025;                 f[5]= 0.034;                 f[6]= 0.054;                 f[7]= 0.074;                 
	    f[8]= 0.026;                 f[9]= 0.068;                 f[10]= 0.099;                f[11]= 0.058;                
	    f[12]= 0.025;                f[13]= 0.047;                f[14]= 0.039;                f[15]= 0.057;                
	    f[16]= 0.051;                f[17]= 0.013;                f[18]= 0.032;                f[19]= 0.073;
	  }
	  break;
	case MTMAM:
	  {
	    daa[1*20+0]= 32;              daa[2*20+0]= 2;    daa[2*20+1]= 4;               daa[3*20+0]= 11;
	    daa[3*20+1]= 0;               daa[3*20+2]= 864;  daa[4*20+0]= 0;               daa[4*20+1]= 186;
	    daa[4*20+2]= 0;               daa[4*20+3]= 0;    daa[5*20+0]= 0;               daa[5*20+1]= 246;
	    daa[5*20+2]= 8;               daa[5*20+3]= 49;   daa[5*20+4]= 0;               daa[6*20+0]= 0;
	    daa[6*20+1]= 0;               daa[6*20+2]= 0;    daa[6*20+3]= 569;             daa[6*20+4]= 0;
	    daa[6*20+5]= 274;             daa[7*20+0]= 78;   daa[7*20+1]= 18;              daa[7*20+2]= 47;
	    daa[7*20+3]= 79;              daa[7*20+4]= 0;    daa[7*20+5]= 0;               daa[7*20+6]= 22;
	    daa[8*20+0]= 8;               daa[8*20+1]= 232;  daa[8*20+2]= 458;             daa[8*20+3]= 11;
	    daa[8*20+4]= 305;             daa[8*20+5]= 550;  daa[8*20+6]= 22;              daa[8*20+7]= 0;
	    daa[9*20+0]= 75;              daa[9*20+1]= 0;    daa[9*20+2]= 19;              daa[9*20+3]= 0;
	    daa[9*20+4]= 41;              daa[9*20+5]= 0;    daa[9*20+6]= 0;               daa[9*20+7]= 0;
	    daa[9*20+8]= 0;               daa[10*20+0]= 21;  daa[10*20+1]= 6;              daa[10*20+2]= 0;
	    daa[10*20+3]= 0;              daa[10*20+4]= 27;  daa[10*20+5]= 20;             daa[10*20+6]= 0;
	    daa[10*20+7]= 0;              daa[10*20+8]= 26;  daa[10*20+9]= 232;            daa[11*20+0]= 0;
	    daa[11*20+1]= 50;             daa[11*20+2]= 408; daa[11*20+3]= 0;              daa[11*20+4]= 0;
	    daa[11*20+5]= 242;            daa[11*20+6]= 215; daa[11*20+7]= 0;              daa[11*20+8]= 0;
	    daa[11*20+9]= 6;              daa[11*20+10]= 4;  daa[12*20+0]= 76;             daa[12*20+1]= 0;
	    daa[12*20+2]= 21;             daa[12*20+3]= 0;   daa[12*20+4]= 0;              daa[12*20+5]= 22;
	    daa[12*20+6]= 0;              daa[12*20+7]= 0;   daa[12*20+8]= 0;              daa[12*20+9]= 378;
	    daa[12*20+10]= 609;           daa[12*20+11]= 59; daa[13*20+0]= 0;              daa[13*20+1]= 0;
	    daa[13*20+2]= 6;              daa[13*20+3]= 5;   daa[13*20+4]= 7;              daa[13*20+5]= 0;
	    daa[13*20+6]= 0;              daa[13*20+7]= 0;   daa[13*20+8]= 0;              daa[13*20+9]= 57;
	    daa[13*20+10]= 246;           daa[13*20+11]= 0;  daa[13*20+12]= 11;            daa[14*20+0]= 53;
	    daa[14*20+1]= 9;              daa[14*20+2]= 33;  daa[14*20+3]= 2;              daa[14*20+4]= 0;
	    daa[14*20+5]= 51;             daa[14*20+6]= 0;   daa[14*20+7]= 0;              daa[14*20+8]= 53;
	    daa[14*20+9]= 5;              daa[14*20+10]= 43; daa[14*20+11]= 18;            daa[14*20+12]= 0;
	    daa[14*20+13]= 17;            daa[15*20+0]= 342; daa[15*20+1]= 3;              daa[15*20+2]= 446;
	    daa[15*20+3]= 16;             daa[15*20+4]= 347; daa[15*20+5]= 30;             daa[15*20+6]= 21;
	    daa[15*20+7]= 112;            daa[15*20+8]= 20;  daa[15*20+9]= 0;              daa[15*20+10]= 74;
	    daa[15*20+11]= 65;            daa[15*20+12]= 47; daa[15*20+13]= 90;            daa[15*20+14]= 202;
	    daa[16*20+0]= 681;            daa[16*20+1]= 0;   daa[16*20+2]= 110;            daa[16*20+3]= 0;
	    daa[16*20+4]= 114;            daa[16*20+5]= 0;   daa[16*20+6]= 4;              daa[16*20+7]= 0;
	    daa[16*20+8]= 1;              daa[16*20+9]= 360; daa[16*20+10]= 34;            daa[16*20+11]= 50;
	    daa[16*20+12]= 691;           daa[16*20+13]= 8;  daa[16*20+14]= 78;            daa[16*20+15]= 614;
	    daa[17*20+0]= 5;              daa[17*20+1]= 16;  daa[17*20+2]= 6;              daa[17*20+3]= 0;
	    daa[17*20+4]= 65;             daa[17*20+5]= 0;   daa[17*20+6]= 0;              daa[17*20+7]= 0;
	    daa[17*20+8]= 0;              daa[17*20+9]= 0;   daa[17*20+10]= 12;            daa[17*20+11]= 0;
	    daa[17*20+12]= 13;            daa[17*20+13]= 0;  daa[17*20+14]= 7;             daa[17*20+15]= 17;
	    daa[17*20+16]= 0;             daa[18*20+0]= 0;   daa[18*20+1]= 0;              daa[18*20+2]= 156;
	    daa[18*20+3]= 0;              daa[18*20+4]= 530; daa[18*20+5]= 54;             daa[18*20+6]= 0;
	    daa[18*20+7]= 1;              daa[18*20+8]= 1525;daa[18*20+9]= 16;             daa[18*20+10]= 25;
	    daa[18*20+11]= 67;            daa[18*20+12]= 0;  daa[18*20+13]= 682;           daa[18*20+14]= 8;
	    daa[18*20+15]= 107;           daa[18*20+16]= 0;  daa[18*20+17]= 14;            daa[19*20+0]= 398;
	    daa[19*20+1]= 0;              daa[19*20+2]= 0;   daa[19*20+3]= 10;             daa[19*20+4]= 0;
	    daa[19*20+5]= 33;             daa[19*20+6]= 20;  daa[19*20+7]= 5;              daa[19*20+8]= 0;
	    daa[19*20+9]= 2220;           daa[19*20+10]= 100;daa[19*20+11]= 0;             daa[19*20+12]= 832;
	    daa[19*20+13]= 6;             daa[19*20+14]= 0;  daa[19*20+15]= 0;             daa[19*20+16]= 237;
	    daa[19*20+17]= 0;             daa[19*20+18]= 0;       
	    
	    f[0]= 0.06920;  f[1]=  0.01840;  f[2]= 0.04000;  f[3]= 0.018600;
	    f[4]= 0.00650;  f[5]=  0.02380;  f[6]= 0.02360;  f[7]= 0.055700;
	    f[8]= 0.02770;  f[9]=  0.09050;  f[10]=0.16750;  f[11]= 0.02210;
	    f[12]=0.05610;  f[13]= 0.06110;  f[14]=0.05360;  f[15]= 0.07250;
	    f[16]=0.08700;  f[17]= 0.02930;  f[18]=0.03400;  f[19]= 0.04280;
	  }
	  break;
	case LG:
	  {
	    daa[1*20+0] = 0.425093;

	    daa[2*20+0] = 0.276818; daa[2*20+1] = 0.751878;

	    daa[3*20+0] = 0.395144; daa[3*20+1] = 0.123954; daa[3*20+2] = 5.076149;
	    
	    daa[4*20+0] = 2.489084; daa[4*20+1] = 0.534551; daa[4*20+2] = 0.528768; daa[4*20+3] = 0.062556;
								 
	    daa[5*20+0] = 0.969894; daa[5*20+1] = 2.807908; daa[5*20+2] = 1.695752; daa[5*20+3] = 0.523386; daa[5*20+4] = 0.084808;

	    daa[6*20+0] = 1.038545; daa[6*20+1] = 0.363970; daa[6*20+2] = 0.541712; daa[6*20+3] = 5.243870; daa[6*20+4] = 0.003499; daa[6*20+5] = 4.128591;

	    daa[7*20+0] = 2.066040; daa[7*20+1] = 0.390192; daa[7*20+2] = 1.437645; daa[7*20+3] = 0.844926; daa[7*20+4] = 0.569265; daa[7*20+5] = 0.267959; daa[7*20+6] = 0.348847;
 
	    daa[8*20+0] = 0.358858; daa[8*20+1] = 2.426601; daa[8*20+2] = 4.509238; daa[8*20+3] = 0.927114; daa[8*20+4] = 0.640543; daa[8*20+5] = 4.813505; daa[8*20+6] = 0.423881; 
	    daa[8*20+7] = 0.311484;

	    daa[9*20+0] = 0.149830; daa[9*20+1] = 0.126991; daa[9*20+2] = 0.191503; daa[9*20+3] = 0.010690; daa[9*20+4] = 0.320627; daa[9*20+5] = 0.072854; daa[9*20+6] = 0.044265; 
	    daa[9*20+7] = 0.008705; daa[9*20+8] = 0.108882; 

	    daa[10*20+0] = 0.395337; daa[10*20+1] = 0.301848; daa[10*20+2] = 0.068427; daa[10*20+3] = 0.015076; daa[10*20+4] = 0.594007; daa[10*20+5] = 0.582457; daa[10*20+6] = 0.069673; 
	    daa[10*20+7] = 0.044261; daa[10*20+8] = 0.366317; daa[10*20+9] = 4.145067 ;

	    daa[11*20+0] = 0.536518; daa[11*20+1] = 6.326067; daa[11*20+2] = 2.145078; daa[11*20+3] = 0.282959; daa[11*20+4] = 0.013266; daa[11*20+5] = 3.234294; daa[11*20+6] = 1.807177; 
	    daa[11*20+7] = 0.296636; daa[11*20+8] = 0.697264; daa[11*20+9] = 0.159069; daa[11*20+10] = 0.137500;


	    daa[12*20+0] = 1.124035; daa[12*20+1] = 0.484133; daa[12*20+2] = 0.371004; daa[12*20+3] = 0.025548; daa[12*20+4] = 0.893680; daa[12*20+5] = 1.672569; daa[12*20+6] = 0.173735; 
	    daa[12*20+7] = 0.139538; daa[12*20+8] = 0.442472; daa[12*20+9] = 4.273607; daa[12*20+10] = 6.312358; daa[12*20+11] = 0.656604;

	    daa[13*20+0] = 0.253701; daa[13*20+1] = 0.052722;daa[13*20+2] = 0.089525; daa[13*20+3] = 0.017416; daa[13*20+4] = 1.105251; daa[13*20+5] = 0.035855; daa[13*20+6] = 0.018811; 
	    daa[13*20+7] = 0.089586; daa[13*20+8] = 0.682139; daa[13*20+9] = 1.112727; daa[13*20+10] = 2.592692; daa[13*20+11] = 0.023918; daa[13*20+12] = 1.798853;

	    daa[14*20+0] = 1.177651; daa[14*20+1] = 0.332533;daa[14*20+2] = 0.161787; daa[14*20+3] = 0.394456; daa[14*20+4] = 0.075382; daa[14*20+5] = 0.624294; daa[14*20+6] = 0.419409; 
	    daa[14*20+7] = 0.196961; daa[14*20+8] = 0.508851; daa[14*20+9] = 0.078281; daa[14*20+10] = 0.249060; daa[14*20+11] = 0.390322; daa[14*20+12] = 0.099849; 
	    daa[14*20+13] = 0.094464;
 
	    daa[15*20+0] = 4.727182; daa[15*20+1] = 0.858151;daa[15*20+2] = 4.008358; daa[15*20+3] = 1.240275; daa[15*20+4] = 2.784478; daa[15*20+5] = 1.223828; daa[15*20+6] = 0.611973; 
	    daa[15*20+7] = 1.739990; daa[15*20+8] = 0.990012; daa[15*20+9] = 0.064105; daa[15*20+10] = 0.182287; daa[15*20+11] = 0.748683; daa[15*20+12] = 0.346960; 
	    daa[15*20+13] = 0.361819; daa[15*20+14] = 1.338132;
 
	    daa[16*20+0] = 2.139501; daa[16*20+1] = 0.578987;daa[16*20+2] = 2.000679; daa[16*20+3] = 0.425860; daa[16*20+4] = 1.143480; daa[16*20+5] = 1.080136; daa[16*20+6] = 0.604545; 
	    daa[16*20+7] = 0.129836; daa[16*20+8] = 0.584262; daa[16*20+9] = 1.033739; daa[16*20+10] = 0.302936; daa[16*20+11] = 1.136863; daa[16*20+12] = 2.020366; 
	    daa[16*20+13] = 0.165001; daa[16*20+14] = 0.571468; daa[16*20+15] = 6.472279;

	    daa[17*20+0] = 0.180717; daa[17*20+1] = 0.593607;daa[17*20+2] = 0.045376; daa[17*20+3] = 0.029890; daa[17*20+4] = 0.670128; daa[17*20+5] = 0.236199; daa[17*20+6] = 0.077852; 
	    daa[17*20+7] = 0.268491; daa[17*20+8] = 0.597054; daa[17*20+9] = 0.111660; daa[17*20+10] = 0.619632; daa[17*20+11] = 0.049906; daa[17*20+12] = 0.696175; 
	    daa[17*20+13] = 2.457121; daa[17*20+14] = 0.095131; daa[17*20+15] = 0.248862; daa[17*20+16] = 0.140825;

	    daa[18*20+0] = 0.218959; daa[18*20+1] = 0.314440;daa[18*20+2] = 0.612025; daa[18*20+3] = 0.135107; daa[18*20+4] = 1.165532; daa[18*20+5] = 0.257336; daa[18*20+6] = 0.120037; 
	    daa[18*20+7] = 0.054679; daa[18*20+8] = 5.306834; daa[18*20+9] = 0.232523; daa[18*20+10] = 0.299648; daa[18*20+11] = 0.131932; daa[18*20+12] = 0.481306; 
	    daa[18*20+13] = 7.803902; daa[18*20+14] = 0.089613; daa[18*20+15] = 0.400547; daa[18*20+16] = 0.245841; daa[18*20+17] = 3.151815;

	    daa[19*20+0] = 2.547870; daa[19*20+1] = 0.170887;daa[19*20+2] = 0.083688; daa[19*20+3] = 0.037967; daa[19*20+4] = 1.959291; daa[19*20+5] = 0.210332; daa[19*20+6] = 0.245034; 
	    daa[19*20+7] = 0.076701; daa[19*20+8] = 0.119013; daa[19*20+9] = 10.649107; daa[19*20+10] = 1.702745; daa[19*20+11] = 0.185202; daa[19*20+12] = 1.898718; 
	    daa[19*20+13] = 0.654683; daa[19*20+14] = 0.296501; daa[19*20+15] = 0.098369; daa[19*20+16] = 2.188158; daa[19*20+17] = 0.189510; daa[19*20+18] = 0.249313;
	    
	    /*f[0] = 0.07906;
	    f[1] = 0.05594; 
	    f[2] = 0.04198; 
	    f[3] = 0.05305; 
	    f[4] = 0.01294; 
	    f[5] = 0.04077; 
	    f[6] = 0.07158; 
	    f[7] = 0.05734; 
	    f[8] = 0.02235; 
	    f[9] = 0.06216; 
	    f[10] = 0.09908; 
	    f[11] = 0.06460; 
	    f[12] = 0.02295; 
	    f[13] = 0.04230; 
	    f[14] = 0.04404; 
	    f[15] = 0.06120; 
	    f[16] = 0.05329; 
	    f[17] = 0.01207; 
	    f[18] = 0.03415; 
	    f[19] = 0.06915; 	   */

	    f[0] = 0.079066; f[1] = 0.055941; f[2] = 0.041977; f[3] = 0.053052;
	    f[4] = 0.012937; f[5] = 0.040767; f[6] = 0.071586; f[7] = 0.057337;
	    f[8] = 0.022355; f[9] = 0.062157; f[10] = 0.099081; f[11] = 0.064600;
	    f[12] = 0.022951; f[13] = 0.042302; f[14] = 0.044040; f[15] = 0.061197;
	    f[16] = 0.053287; f[17] = 0.012066; f[18] = 0.034155; f[19] = 0.069146;
	  }	  
	  break;
	case MTART:
	  {
	   
	    daa[1*20+0]=   0.2;
	    daa[2*20+0]=   0.2;
	    daa[2*20+1]=   0.2;
	    daa[3*20+0]=   1;
           daa[3*20+1]=   4;
           daa[3*20+2]=   500;
           daa[4*20+0]=   254;
           daa[4*20+1]=   36;
           daa[4*20+2]=   98;
           daa[4*20+3]=   11;
           daa[5*20+0]=   0.2;
           daa[5*20+1]=   154;
           daa[5*20+2]=   262;
           daa[5*20+3]=   0.2;
           daa[5*20+4]=   0.2;
           daa[6*20+0]=   0.2;
           daa[6*20+1]=   0.2;
           daa[6*20+2]=   183;
           daa[6*20+3]=   862;
           daa[6*20+4]=   0.2;
           daa[6*20+5]=   262;
           daa[7*20+0]=   200;
           daa[7*20+1]=   0.2;
           daa[7*20+2]=   121;
           daa[7*20+3]=   12;
           daa[7*20+4]=   81;
           daa[7*20+5]=   3;
           daa[7*20+6]=   44;
           daa[8*20+0]=   0.2;
           daa[8*20+1]=   41;
           daa[8*20+2]=   180;
           daa[8*20+3]=   0.2;
           daa[8*20+4]=   12;
           daa[8*20+5]=   314;
           daa[8*20+6]=   15;
           daa[8*20+7]=   0.2;
           daa[9*20+0]=   26;
           daa[9*20+1]=   2;
           daa[9*20+2]=   21;
           daa[9*20+3]=   7;
           daa[9*20+4]=   63;
           daa[9*20+5]=   11;
           daa[9*20+6]=   7;
           daa[9*20+7]=   3;
           daa[9*20+8]=   0.2;
           daa[10*20+0]=  4;
           daa[10*20+1]=  2;
           daa[10*20+2]=  13;
           daa[10*20+3]=  1;
           daa[10*20+4]=  79;
           daa[10*20+5]=  16;
           daa[10*20+6]=  2;
           daa[10*20+7]=  1;
           daa[10*20+8]=  6;
           daa[10*20+9]=  515;
           daa[11*20+0]=  0.2;
           daa[11*20+1]=  209;
           daa[11*20+2]=  467;
           daa[11*20+3]=  2;
           daa[11*20+4]=  0.2;
           daa[11*20+5]=  349;
           daa[11*20+6]=  106;
           daa[11*20+7]=  0.2;
           daa[11*20+8]=  0.2;
           daa[11*20+9]=  3;
           daa[11*20+10]= 4;
           daa[12*20+0]=  121;
           daa[12*20+1]=  5;
           daa[12*20+2]=  79;
           daa[12*20+3]=  0.2;
           daa[12*20+4]=  312;
           daa[12*20+5]=  67;
           daa[12*20+6]=  0.2;
           daa[12*20+7]=  56;
           daa[12*20+8]=  0.2;
           daa[12*20+9]=  515;
           daa[12*20+10]= 885;
           daa[12*20+11]= 106;
           daa[13*20+0]=  13;
           daa[13*20+1]=  5;
           daa[13*20+2]=  20;
           daa[13*20+3]=  0.2;
           daa[13*20+4]=  184;
           daa[13*20+5]=  0.2;
           daa[13*20+6]=  0.2;
           daa[13*20+7]=  1;
           daa[13*20+8]=  14;
           daa[13*20+9]=  118;
           daa[13*20+10]= 263;
           daa[13*20+11]= 11;
           daa[13*20+12]= 322;
           daa[14*20+0]=  49;
           daa[14*20+1]=  0.2;
           daa[14*20+2]=  17;
           daa[14*20+3]=  0.2;
           daa[14*20+4]=  0.2;
           daa[14*20+5]=  39;
           daa[14*20+6]=  8;
           daa[14*20+7]=  0.2;
           daa[14*20+8]=  1;
           daa[14*20+9]=  0.2;
           daa[14*20+10]= 12;
           daa[14*20+11]= 17;
           daa[14*20+12]= 5;
           daa[14*20+13]= 15;
           daa[15*20+0]=  673;
           daa[15*20+1]=  3;
           daa[15*20+2]=  398;
           daa[15*20+3]=  44;
           daa[15*20+4]=  664;
           daa[15*20+5]=  52;
           daa[15*20+6]=  31;
           daa[15*20+7]=  226;
           daa[15*20+8]=  11;
           daa[15*20+9]=  7;
           daa[15*20+10]= 8;
           daa[15*20+11]= 144;
           daa[15*20+12]= 112;
           daa[15*20+13]= 36;
           daa[15*20+14]= 87;
           daa[16*20+0]=  244;
           daa[16*20+1]=  0.2;
           daa[16*20+2]=  166;
           daa[16*20+3]=  0.2;
           daa[16*20+4]=  183;
           daa[16*20+5]=  44;
           daa[16*20+6]=  43;
           daa[16*20+7]=  0.2;
           daa[16*20+8]=  19;
           daa[16*20+9]=  204;
           daa[16*20+10]= 48;
           daa[16*20+11]= 70;
           daa[16*20+12]= 289;
           daa[16*20+13]= 14;
           daa[16*20+14]= 47;
           daa[16*20+15]= 660;
           daa[17*20+0]=  0.2;
           daa[17*20+1]=  0.2;
           daa[17*20+2]=  8;
           daa[17*20+3]=  0.2;
           daa[17*20+4]=  22;
           daa[17*20+5]=  7;
           daa[17*20+6]=  11;
           daa[17*20+7]=  2;
           daa[17*20+8]=  0.2;
           daa[17*20+9]=  0.2;
           daa[17*20+10]= 21;
           daa[17*20+11]= 16;
           daa[17*20+12]= 71;
           daa[17*20+13]= 54;
           daa[17*20+14]= 0.2;
           daa[17*20+15]= 2;
           daa[17*20+16]= 0.2;
           daa[18*20+0]=  1;
           daa[18*20+1]=  4;
           daa[18*20+2]=  251;
           daa[18*20+3]=  0.2;
           daa[18*20+4]=  72;
           daa[18*20+5]=  87;
           daa[18*20+6]=  8;
           daa[18*20+7]=  9;
           daa[18*20+8]=  191;
           daa[18*20+9]=  12;
           daa[18*20+10]= 20;
           daa[18*20+11]= 117;
           daa[18*20+12]= 71;
           daa[18*20+13]= 792;
           daa[18*20+14]= 18;
           daa[18*20+15]= 30;
           daa[18*20+16]= 46;
           daa[18*20+17]= 38;
           daa[19*20+0]=  340;
           daa[19*20+1]=  0.2;
           daa[19*20+2]=  23;
           daa[19*20+3]=  0.2;
           daa[19*20+4]=  350;
           daa[19*20+5]=  0.2;
           daa[19*20+6]=  14;
           daa[19*20+7]=  3;
           daa[19*20+8]=  0.2;
           daa[19*20+9]=  1855;
           daa[19*20+10]= 85;
           daa[19*20+11]= 26;
           daa[19*20+12]= 281;
           daa[19*20+13]= 52;
           daa[19*20+14]= 32;
           daa[19*20+15]= 61;
           daa[19*20+16]= 544;
           daa[19*20+17]= 0.2;
           daa[19*20+18]= 2;
           
           f[0]=  0.054116;
           f[1]=  0.018227;
           f[2]=  0.039903;
           f[3]=  0.020160;
           f[4]=  0.009709;
           f[5]=  0.018781;
           f[6]=  0.024289;
           f[7]=  0.068183;
           f[8]=  0.024518;
           f[9]=  0.092638;
           f[10]= 0.148658;
           f[11]= 0.021718;
           f[12]= 0.061453;
           f[13]= 0.088668;
           f[14]= 0.041826;
           f[15]= 0.091030;
           f[16]= 0.049194;
           f[17]= 0.029786;
           f[18]= 0.039443;
           f[19]= 0.057700;
	  }
	  break;
	case MTZOA:
	  {
           daa[1*20+0]=   3.3;
           daa[2*20+0]=   1.7;
           daa[2*20+1]=   33.6;
           daa[3*20+0]=   16.1;
           daa[3*20+1]=   3.2;
           daa[3*20+2]=   617.0;
           daa[4*20+0]=   272.5;
           daa[4*20+1]=   61.1;
           daa[4*20+2]=   94.6;
           daa[4*20+3]=   9.5;
           daa[5*20+0]=   7.3;
           daa[5*20+1]=   231.0;
           daa[5*20+2]=   190.3;
           daa[5*20+3]=   19.3;
           daa[5*20+4]=   49.1;
           daa[6*20+0]=   17.1;
           daa[6*20+1]=   6.4;
           daa[6*20+2]=   174.0;
           daa[6*20+3]=   883.6;
           daa[6*20+4]=   3.4;
           daa[6*20+5]=   349.4;
           daa[7*20+0]=   289.3;
           daa[7*20+1]=   7.2;
           daa[7*20+2]=   99.3;
           daa[7*20+3]=   26.0;
           daa[7*20+4]=   82.4;
           daa[7*20+5]=   8.9;
           daa[7*20+6]=   43.1;
           daa[8*20+0]=   2.3;
           daa[8*20+1]=   61.7;
           daa[8*20+2]=   228.9;
           daa[8*20+3]=   55.6;
           daa[8*20+4]=   37.5;
           daa[8*20+5]=   421.8;
           daa[8*20+6]=   14.9;
           daa[8*20+7]=   7.4;
           daa[9*20+0]=   33.2;
           daa[9*20+1]=   0.2;
           daa[9*20+2]=   24.3;
           daa[9*20+3]=   1.5;
           daa[9*20+4]=   48.8;
           daa[9*20+5]=   0.2;
           daa[9*20+6]=   7.3;
           daa[9*20+7]=   3.4;
           daa[9*20+8]=   1.6;
           daa[10*20+0]=  15.6;
           daa[10*20+1]=  4.1;
           daa[10*20+2]=  7.9;
           daa[10*20+3]=  0.5;
           daa[10*20+4]=  59.7;
           daa[10*20+5]=  23.0;
           daa[10*20+6]=  1.0;
           daa[10*20+7]=  3.5;
           daa[10*20+8]=  6.6;
           daa[10*20+9]=  425.2;
           daa[11*20+0]=  0.2;
           daa[11*20+1]=  292.3;
           daa[11*20+2]=  413.4;
           daa[11*20+3]=  0.2;
           daa[11*20+4]=  0.2;
           daa[11*20+5]=  334.0;
           daa[11*20+6]=  163.2;
           daa[11*20+7]=  10.1;
           daa[11*20+8]=  23.9;
           daa[11*20+9]=  8.4;
           daa[11*20+10]= 6.7;
           daa[12*20+0]=  136.5;
           daa[12*20+1]=  3.8;
           daa[12*20+2]=  73.7;
           daa[12*20+3]=  0.2;
           daa[12*20+4]=  264.8;
           daa[12*20+5]=  83.9;
           daa[12*20+6]=  0.2;
           daa[12*20+7]=  52.2;
           daa[12*20+8]=  7.1;
           daa[12*20+9]=  449.7;
           daa[12*20+10]= 636.3;
           daa[12*20+11]= 83.0;
           daa[13*20+0]=  26.5;
           daa[13*20+1]=  0.2;
           daa[13*20+2]=  12.9;
           daa[13*20+3]=  2.0;
           daa[13*20+4]=  167.8;
           daa[13*20+5]=  9.5;
           daa[13*20+6]=  0.2;
           daa[13*20+7]=  5.8;
           daa[13*20+8]=  13.1;
           daa[13*20+9]=  90.3;
           daa[13*20+10]= 234.2;
           daa[13*20+11]= 16.3;
           daa[13*20+12]= 215.6;
           daa[14*20+0]=  61.8;
           daa[14*20+1]=  7.5;
           daa[14*20+2]=  22.6;
           daa[14*20+3]=  0.2;
           daa[14*20+4]=  8.1;
           daa[14*20+5]=  52.2;
           daa[14*20+6]=  20.6;
           daa[14*20+7]=  1.3;
           daa[14*20+8]=  15.6;
           daa[14*20+9]=  2.6;
           daa[14*20+10]= 11.4;
           daa[14*20+11]= 24.3;
           daa[14*20+12]= 5.4;
           daa[14*20+13]= 10.5;
           daa[15*20+0]=  644.9;
           daa[15*20+1]=  11.8;
           daa[15*20+2]=  420.2;
           daa[15*20+3]=  51.4;
           daa[15*20+4]=  656.3;
           daa[15*20+5]=  96.4;
           daa[15*20+6]=  38.4;
           daa[15*20+7]=  257.1;
           daa[15*20+8]=  23.1;
           daa[15*20+9]=  7.2;
           daa[15*20+10]= 15.2;
           daa[15*20+11]= 144.9;
           daa[15*20+12]= 95.3;
           daa[15*20+13]= 32.2;
           daa[15*20+14]= 79.7;
           daa[16*20+0]=  378.1;
           daa[16*20+1]=  3.2;
           daa[16*20+2]=  184.6;
           daa[16*20+3]=  2.3;
           daa[16*20+4]=  199.0;
           daa[16*20+5]=  39.4;
           daa[16*20+6]=  34.5;
           daa[16*20+7]=  5.2;
           daa[16*20+8]=  19.4;
           daa[16*20+9]=  222.3;
           daa[16*20+10]= 50.0;
           daa[16*20+11]= 75.5;
           daa[16*20+12]= 305.1;
           daa[16*20+13]= 19.3;
           daa[16*20+14]= 56.9;
           daa[16*20+15]= 666.3;
           daa[17*20+0]=  3.1;
           daa[17*20+1]=  16.9;
           daa[17*20+2]=  6.4;
           daa[17*20+3]=  0.2;
           daa[17*20+4]=  36.1;
           daa[17*20+5]=  6.1;
           daa[17*20+6]=  3.5;
           daa[17*20+7]=  12.3;
           daa[17*20+8]=  4.5;
           daa[17*20+9]=  9.7;
           daa[17*20+10]= 27.2;
           daa[17*20+11]= 6.6;
           daa[17*20+12]= 48.7;
           daa[17*20+13]= 58.2;
           daa[17*20+14]= 1.3;
           daa[17*20+15]= 10.3;
           daa[17*20+16]= 3.6;
           daa[18*20+0]=  2.1;
           daa[18*20+1]=  13.8;
           daa[18*20+2]=  141.6;
           daa[18*20+3]=  13.9;
           daa[18*20+4]=  76.7;
           daa[18*20+5]=  52.3;
           daa[18*20+6]=  10.0;
           daa[18*20+7]=  4.3;
           daa[18*20+8]=  266.5;
           daa[18*20+9]=  13.1;
           daa[18*20+10]= 5.7;
           daa[18*20+11]= 45.0;
           daa[18*20+12]= 41.4;
           daa[18*20+13]= 590.5;
           daa[18*20+14]= 4.2;
           daa[18*20+15]= 29.7;
           daa[18*20+16]= 29.0;
           daa[18*20+17]= 79.8;
           daa[19*20+0]=  321.9;
           daa[19*20+1]=  5.1;
           daa[19*20+2]=  7.1;
           daa[19*20+3]=  3.7;
           daa[19*20+4]=  243.8;
           daa[19*20+5]=  9.0;
           daa[19*20+6]=  16.3;
           daa[19*20+7]=  23.7;
           daa[19*20+8]=  0.3;
           daa[19*20+9]=  1710.6;
           daa[19*20+10]= 126.1;
           daa[19*20+11]= 11.1;
           daa[19*20+12]= 279.6;
           daa[19*20+13]= 59.6;
           daa[19*20+14]= 17.9;
           daa[19*20+15]= 49.5;
           daa[19*20+16]= 396.4;
           daa[19*20+17]= 13.7;
           daa[19*20+18]= 15.6;
           
           f[0]=  0.069;
           f[1]=  0.021;
           f[2]=  0.030;
           f[3]=  0.020;
           f[4]=  0.010;
           f[5]=  0.019;
           f[6]=  0.025;
           f[7]=  0.072;
           f[8]=  0.027;
           f[9]=  0.085;
           f[10]= 0.157;
           f[11]= 0.019;
           f[12]= 0.051;
           f[13]= 0.082;
           f[14]= 0.045;
           f[15]= 0.081;
           f[16]= 0.056;
           f[17]= 0.028;
           f[18]= 0.037;
           f[19]= 0.066;
	  }
	  break;
	case PMB:
	  {
           daa[1*20+0]=   0.674995699;
           daa[2*20+0]=   0.589645178;
           daa[2*20+1]=   1.189067034;
           daa[3*20+0]=   0.462499504;
           daa[3*20+1]=   0.605460903;
           daa[3*20+2]=   3.573373315;
           daa[4*20+0]=   1.065445546;
           daa[4*20+1]=   0.31444833;
           daa[4*20+2]=   0.589852457;
           daa[4*20+3]=   0.246951424;
           daa[5*20+0]=   1.111766964;
           daa[5*20+1]=   2.967840934;
           daa[5*20+2]=   2.299755865;
           daa[5*20+3]=   1.686058219;
           daa[5*20+4]=   0.245163782;
           daa[6*20+0]=   1.046334652;
           daa[6*20+1]=   1.201770702;
           daa[6*20+2]=   1.277836748;
           daa[6*20+3]=   4.399995525;
           daa[6*20+4]=   0.091071867;
           daa[6*20+5]=   4.15967899;
           daa[7*20+0]=   1.587964372;
           daa[7*20+1]=   0.523770553;
           daa[7*20+2]=   1.374854049;
           daa[7*20+3]=   0.734992057;
           daa[7*20+4]=   0.31706632;
           daa[7*20+5]=   0.596789898;
           daa[7*20+6]=   0.463812837;
           daa[8*20+0]=   0.580830874;
           daa[8*20+1]=   1.457127446;
           daa[8*20+2]=   2.283037894;
           daa[8*20+3]=   0.839348444;
           daa[8*20+4]=   0.411543728;
           daa[8*20+5]=   1.812173605;
           daa[8*20+6]=   0.877842609;
           daa[8*20+7]=   0.476331437;
           daa[9*20+0]=   0.464590585;
           daa[9*20+1]=   0.35964586;
           daa[9*20+2]=   0.426069419;
           daa[9*20+3]=   0.266775558;
           daa[9*20+4]=   0.417547309;
           daa[9*20+5]=   0.315256838;
           daa[9*20+6]=   0.30421529;
           daa[9*20+7]=   0.180198883;
           daa[9*20+8]=   0.285186418;
           daa[10*20+0]=  0.804404505;
           daa[10*20+1]=  0.520701585;
           daa[10*20+2]=  0.41009447;
           daa[10*20+3]=  0.269124919;
           daa[10*20+4]=  0.450795211;
           daa[10*20+5]=  0.625792937;
           daa[10*20+6]=  0.32078471;
           daa[10*20+7]=  0.259854426;
           daa[10*20+8]=  0.363981358;
           daa[10*20+9]=  4.162454693;
           daa[11*20+0]=  0.831998835;
           daa[11*20+1]=  4.956476453;
           daa[11*20+2]=  2.037575629;
           daa[11*20+3]=  1.114178954;
           daa[11*20+4]=  0.274163536;
           daa[11*20+5]=  3.521346591;
           daa[11*20+6]=  2.415974716;
           daa[11*20+7]=  0.581001076;
           daa[11*20+8]=  0.985885486;
           daa[11*20+9]=  0.374784947;
           daa[11*20+10]= 0.498011337;
           daa[12*20+0]=  1.546725076;
           daa[12*20+1]=  0.81346254;
           daa[12*20+2]=  0.737846301;
           daa[12*20+3]=  0.341932741;
           daa[12*20+4]=  0.618614612;
           daa[12*20+5]=  2.067388546;
           daa[12*20+6]=  0.531773639;
           daa[12*20+7]=  0.465349326;
           daa[12*20+8]=  0.380925433;
           daa[12*20+9]=  3.65807012;
           daa[12*20+10]= 5.002338375;
           daa[12*20+11]= 0.661095832;
           daa[13*20+0]=  0.546169219;
           daa[13*20+1]=  0.303437244;
           daa[13*20+2]=  0.425193716;
           daa[13*20+3]=  0.219005213;
           daa[13*20+4]=  0.669206193;
           daa[13*20+5]=  0.406042546;
           daa[13*20+6]=  0.224154698;
           daa[13*20+7]=  0.35402891;
           daa[13*20+8]=  0.576231691;
           daa[13*20+9]=  1.495264661;
           daa[13*20+10]= 2.392638293;
           daa[13*20+11]= 0.269496317;
           daa[13*20+12]= 2.306919847;
           daa[14*20+0]=  1.241586045;
           daa[14*20+1]=  0.65577338;
           daa[14*20+2]=  0.711495595;
           daa[14*20+3]=  0.775624818;
           daa[14*20+4]=  0.198679914;
           daa[14*20+5]=  0.850116543;
           daa[14*20+6]=  0.794584081;
           daa[14*20+7]=  0.588254139;
           daa[14*20+8]=  0.456058589;
           daa[14*20+9]=  0.366232942;
           daa[14*20+10]= 0.430073179;
           daa[14*20+11]= 1.036079005;
           daa[14*20+12]= 0.337502282;
           daa[14*20+13]= 0.481144863;
           daa[15*20+0]=  3.452308792;
           daa[15*20+1]=  0.910144334;
           daa[15*20+2]=  2.572577221;
           daa[15*20+3]=  1.440896785;
           daa[15*20+4]=  0.99870098;
           daa[15*20+5]=  1.348272505;
           daa[15*20+6]=  1.205509425;
           daa[15*20+7]=  1.402122097;
           daa[15*20+8]=  0.799966711;
           daa[15*20+9]=  0.530641901;
           daa[15*20+10]= 0.402471997;
           daa[15*20+11]= 1.234648153;
           daa[15*20+12]= 0.945453716;
           daa[15*20+13]= 0.613230817;
           daa[15*20+14]= 1.217683028;
           daa[16*20+0]=  1.751412803;
           daa[16*20+1]=  0.89517149;
           daa[16*20+2]=  1.823161023;
           daa[16*20+3]=  0.994227284;
           daa[16*20+4]=  0.847312432;
           daa[16*20+5]=  1.320626678;
           daa[16*20+6]=  0.949599791;
           daa[16*20+7]=  0.542185658;
           daa[16*20+8]=  0.83039281;
           daa[16*20+9]=  1.114132523;
           daa[16*20+10]= 0.779827336;
           daa[16*20+11]= 1.290709079;
           daa[16*20+12]= 1.551488041;
           daa[16*20+13]= 0.718895136;
           daa[16*20+14]= 0.780913179;
           daa[16*20+15]= 4.448982584;
           daa[17*20+0]=  0.35011051;
           daa[17*20+1]=  0.618778365;
           daa[17*20+2]=  0.422407388;
           daa[17*20+3]=  0.362495245;
           daa[17*20+4]=  0.445669347;
           daa[17*20+5]=  0.72038474;
           daa[17*20+6]=  0.261258229;
           daa[17*20+7]=  0.37874827;
           daa[17*20+8]=  0.72436751;
           daa[17*20+9]=  0.516260502;
           daa[17*20+10]= 0.794797115;
           daa[17*20+11]= 0.43340962;
           daa[17*20+12]= 0.768395107;
           daa[17*20+13]= 3.29519344;
           daa[17*20+14]= 0.499869138;
           daa[17*20+15]= 0.496334956;
           daa[17*20+16]= 0.38372361;
           daa[18*20+0]=  0.573154753;
           daa[18*20+1]=  0.628599063;
           daa[18*20+2]=  0.720013799;
           daa[18*20+3]=  0.436220437;
           daa[18*20+4]=  0.55626163;
           daa[18*20+5]=  0.728970584;
           daa[18*20+6]=  0.50720003;
           daa[18*20+7]=  0.284727562;
           daa[18*20+8]=  2.210952064;
           daa[18*20+9]=  0.570562395;
           daa[18*20+10]= 0.811019594;
           daa[18*20+11]= 0.664884513;
           daa[18*20+12]= 0.93253606;
           daa[18*20+13]= 5.894735673;
           daa[18*20+14]= 0.433748126;
           daa[18*20+15]= 0.593795813;
           daa[18*20+16]= 0.523549536;
           daa[18*20+17]= 2.996248013;
           daa[19*20+0]=  2.063050067;
           daa[19*20+1]=  0.388680158;
           daa[19*20+2]=  0.474418852;
           daa[19*20+3]=  0.275658381;
           daa[19*20+4]=  0.998911631;
           daa[19*20+5]=  0.634408285;
           daa[19*20+6]=  0.527640634;
           daa[19*20+7]=  0.314700907;
           daa[19*20+8]=  0.305792277;
           daa[19*20+9]=  8.002789424;
           daa[19*20+10]= 2.113077156;
           daa[19*20+11]= 0.526184203;
           daa[19*20+12]= 1.737356217;
           daa[19*20+13]= 0.983844803;
           daa[19*20+14]= 0.551333603;
           daa[19*20+15]= 0.507506011;
           daa[19*20+16]= 1.89965079;
           daa[19*20+17]= 0.429570747;
           daa[19*20+18]= 0.716795463;
           
           f[0]=  0.076;
           f[1]=  0.054;
           f[2]=  0.038;
           f[3]=  0.045;
           f[4]=  0.028;
           f[5]=  0.034;
           f[6]=  0.053;
           f[7]=  0.078;
           f[8]=  0.030;
           f[9]=  0.060;
           f[10]= 0.096;
           f[11]= 0.052;
           f[12]= 0.022;
           f[13]= 0.045;
           f[14]= 0.042;
           f[15]= 0.068;
           f[16]= 0.056;
           f[17]= 0.016;
           f[18]= 0.036;
           f[19]= 0.071;
	  }
	  break;
	case HIVB:
	  {
           daa[1*20+0]=   0.30750700;
           daa[2*20+0]=   0.00500000;
           daa[2*20+1]=   0.29554300;
           daa[3*20+0]=   1.45504000;
           daa[3*20+1]=   0.00500000;
           daa[3*20+2]=   17.66120000;
           daa[4*20+0]=   0.12375800;
           daa[4*20+1]=   0.35172100;
           daa[4*20+2]=   0.08606420;
           daa[4*20+3]=   0.00500000;
           daa[5*20+0]=   0.05511280;
           daa[5*20+1]=   3.42150000;
           daa[5*20+2]=   0.67205200;
           daa[5*20+3]=   0.00500000;
           daa[5*20+4]=   0.00500000;
           daa[6*20+0]=   1.48135000;
           daa[6*20+1]=   0.07492180;
           daa[6*20+2]=   0.07926330;
           daa[6*20+3]=   10.58720000;
           daa[6*20+4]=   0.00500000;
           daa[6*20+5]=   2.56020000;
           daa[7*20+0]=   2.13536000;
           daa[7*20+1]=   3.65345000;
           daa[7*20+2]=   0.32340100;
           daa[7*20+3]=   2.83806000;
           daa[7*20+4]=   0.89787100;
           daa[7*20+5]=   0.06191370;
           daa[7*20+6]=   3.92775000;
           daa[8*20+0]=   0.08476130;
           daa[8*20+1]=   9.04044000;
           daa[8*20+2]=   7.64585000;
           daa[8*20+3]=   1.91690000;
           daa[8*20+4]=   0.24007300;
           daa[8*20+5]=   7.05545000;
           daa[8*20+6]=   0.11974000;
           daa[8*20+7]=   0.00500000;
           daa[9*20+0]=   0.00500000;
           daa[9*20+1]=   0.67728900;
           daa[9*20+2]=   0.68056500;
           daa[9*20+3]=   0.01767920;
           daa[9*20+4]=   0.00500000;
           daa[9*20+5]=   0.00500000;
           daa[9*20+6]=   0.00609079;
           daa[9*20+7]=   0.00500000;
           daa[9*20+8]=   0.10311100;
           daa[10*20+0]=  0.21525600;
           daa[10*20+1]=  0.70142700;
           daa[10*20+2]=  0.00500000;
           daa[10*20+3]=  0.00876048;
           daa[10*20+4]=  0.12977700;
           daa[10*20+5]=  1.49456000;
           daa[10*20+6]=  0.00500000;
           daa[10*20+7]=  0.00500000;
           daa[10*20+8]=  1.74171000;
           daa[10*20+9]=  5.95879000;
           daa[11*20+0]=  0.00500000;
           daa[11*20+1]=  20.45000000;
           daa[11*20+2]=  7.90443000;
           daa[11*20+3]=  0.00500000;
           daa[11*20+4]=  0.00500000;
           daa[11*20+5]=  6.54737000;
           daa[11*20+6]=  4.61482000;
           daa[11*20+7]=  0.52170500;
           daa[11*20+8]=  0.00500000;
           daa[11*20+9]=  0.32231900;
           daa[11*20+10]= 0.08149950;
           daa[12*20+0]=  0.01866430;
           daa[12*20+1]=  2.51394000;
           daa[12*20+2]=  0.00500000;
           daa[12*20+3]=  0.00500000;
           daa[12*20+4]=  0.00500000;
           daa[12*20+5]=  0.30367600;
           daa[12*20+6]=  0.17578900;
           daa[12*20+7]=  0.00500000;
           daa[12*20+8]=  0.00500000;
           daa[12*20+9]=  11.20650000;
           daa[12*20+10]= 5.31961000;
           daa[12*20+11]= 1.28246000;
           daa[13*20+0]=  0.01412690;
           daa[13*20+1]=  0.00500000;
           daa[13*20+2]=  0.00500000;
           daa[13*20+3]=  0.00500000;
           daa[13*20+4]=  9.29815000;
           daa[13*20+5]=  0.00500000;
           daa[13*20+6]=  0.00500000;
           daa[13*20+7]=  0.29156100;
           daa[13*20+8]=  0.14555800;
           daa[13*20+9]=  3.39836000;
           daa[13*20+10]= 8.52484000;
           daa[13*20+11]= 0.03426580;
           daa[13*20+12]= 0.18802500;
           daa[14*20+0]=  2.12217000;
           daa[14*20+1]=  1.28355000;
           daa[14*20+2]=  0.00739578;
           daa[14*20+3]=  0.03426580;
           daa[14*20+4]=  0.00500000;
           daa[14*20+5]=  4.47211000;
           daa[14*20+6]=  0.01202260;
           daa[14*20+7]=  0.00500000;
           daa[14*20+8]=  2.45318000;
           daa[14*20+9]=  0.04105930;
           daa[14*20+10]= 2.07757000;
           daa[14*20+11]= 0.03138620;
           daa[14*20+12]= 0.00500000;
           daa[14*20+13]= 0.00500000;
           daa[15*20+0]=  2.46633000;
           daa[15*20+1]=  3.47910000;
           daa[15*20+2]=  13.14470000;
           daa[15*20+3]=  0.52823000;
           daa[15*20+4]=  4.69314000;
           daa[15*20+5]=  0.11631100;
           daa[15*20+6]=  0.00500000;
           daa[15*20+7]=  4.38041000;
           daa[15*20+8]=  0.38274700;
           daa[15*20+9]=  1.21803000;
           daa[15*20+10]= 0.92765600;
           daa[15*20+11]= 0.50411100;
           daa[15*20+12]= 0.00500000;
           daa[15*20+13]= 0.95647200;
           daa[15*20+14]= 5.37762000;
           daa[16*20+0]=  15.91830000;
           daa[16*20+1]=  2.86868000;
           daa[16*20+2]=  6.88667000;
           daa[16*20+3]=  0.27472400;
           daa[16*20+4]=  0.73996900;
           daa[16*20+5]=  0.24358900;
           daa[16*20+6]=  0.28977400;
           daa[16*20+7]=  0.36961500;
           daa[16*20+8]=  0.71159400;
           daa[16*20+9]=  8.61217000;
           daa[16*20+10]= 0.04376730;
           daa[16*20+11]= 4.67142000;
           daa[16*20+12]= 4.94026000;
           daa[16*20+13]= 0.01412690;
           daa[16*20+14]= 2.01417000;
           daa[16*20+15]= 8.93107000;
           daa[17*20+0]=  0.00500000;
           daa[17*20+1]=  0.99133800;
           daa[17*20+2]=  0.00500000;
           daa[17*20+3]=  0.00500000;
           daa[17*20+4]=  2.63277000;
           daa[17*20+5]=  0.02665600;
           daa[17*20+6]=  0.00500000;
           daa[17*20+7]=  1.21674000;
           daa[17*20+8]=  0.06951790;
           daa[17*20+9]=  0.00500000;
           daa[17*20+10]= 0.74884300;
           daa[17*20+11]= 0.00500000;
           daa[17*20+12]= 0.08907800;
           daa[17*20+13]= 0.82934300;
           daa[17*20+14]= 0.04445060;
           daa[17*20+15]= 0.02487280;
           daa[17*20+16]= 0.00500000;
           daa[18*20+0]=  0.00500000;
           daa[18*20+1]=  0.00991826;
           daa[18*20+2]=  1.76417000;
           daa[18*20+3]=  0.67465300;
           daa[18*20+4]=  7.57932000;
           daa[18*20+5]=  0.11303300;
           daa[18*20+6]=  0.07926330;
           daa[18*20+7]=  0.00500000;
           daa[18*20+8]=  18.69430000;
           daa[18*20+9]=  0.14816800;
           daa[18*20+10]= 0.11198600;
           daa[18*20+11]= 0.00500000;
           daa[18*20+12]= 0.00500000;
           daa[18*20+13]= 15.34000000;
           daa[18*20+14]= 0.03043810;
           daa[18*20+15]= 0.64802400;
           daa[18*20+16]= 0.10565200;
           daa[18*20+17]= 1.28022000;
           daa[19*20+0]=  7.61428000;
           daa[19*20+1]=  0.08124540;
           daa[19*20+2]=  0.02665600;
           daa[19*20+3]=  1.04793000;
           daa[19*20+4]=  0.42002700;
           daa[19*20+5]=  0.02091530;
           daa[19*20+6]=  1.02847000;
           daa[19*20+7]=  0.95315500;
           daa[19*20+8]=  0.00500000;
           daa[19*20+9]=  17.73890000;
           daa[19*20+10]= 1.41036000;
           daa[19*20+11]= 0.26582900;
           daa[19*20+12]= 6.85320000;
           daa[19*20+13]= 0.72327400;
           daa[19*20+14]= 0.00500000;
           daa[19*20+15]= 0.07492180;
           daa[19*20+16]= 0.70922600;
           daa[19*20+17]= 0.00500000;
           daa[19*20+18]= 0.04105930;
           
           /*f[0]=  0.060;
           f[1]=  0.066;
           f[2]=  0.044;
           f[3]=  0.042;
           f[4]=  0.020;
           f[5]=  0.054;
           f[6]=  0.071;
           f[7]=  0.072;
           f[8]=  0.022;
           f[9]=  0.070;
           f[10]= 0.099;
           f[11]= 0.057;
           f[12]= 0.020;
           f[13]= 0.029;
           f[14]= 0.046;
           f[15]= 0.051;
           f[16]= 0.054;
           f[17]= 0.033;
           f[18]= 0.028;
           f[19]= 0.062;*/

	   f[0]= 0.060490222;           f[1]= 0.066039665;           f[2]= 0.044127815;           f[3]= 0.042109048;
           f[4]= 0.020075899;           f[5]= 0.053606488;           f[6]= 0.071567447;           f[7]= 0.072308239;
           f[8]= 0.022293943;           f[9]= 0.069730629;           f[10]= 0.098851122;          f[11]= 0.056968211;
           f[12]= 0.019768318;          f[13]= 0.028809447;          f[14]= 0.046025282;          f[15]= 0.05060433;
           f[16]= 0.053636813;          f[17]= 0.033011601;          f[18]= 0.028350243;          f[19]= 0.061625237;
	  }
	  break;
	case HIVW:
	  {
           daa[1*20+0]=   0.0744808;
           daa[2*20+0]=   0.6175090;
           daa[2*20+1]=   0.1602400;
           daa[3*20+0]=   4.4352100;
           daa[3*20+1]=   0.0674539;
           daa[3*20+2]=   29.4087000;
           daa[4*20+0]=   0.1676530;
           daa[4*20+1]=   2.8636400;
           daa[4*20+2]=   0.0604932;
           daa[4*20+3]=   0.0050000;
           daa[5*20+0]=   0.0050000;
           daa[5*20+1]=   10.6746000;
           daa[5*20+2]=   0.3420680;
           daa[5*20+3]=   0.0050000;
           daa[5*20+4]=   0.0050000;
           daa[6*20+0]=   5.5632500;
           daa[6*20+1]=   0.0251632;
           daa[6*20+2]=   0.2015260;
           daa[6*20+3]=   12.1233000;
           daa[6*20+4]=   0.0050000;
           daa[6*20+5]=   3.2065600;
           daa[7*20+0]=   1.8685000;
           daa[7*20+1]=   13.4379000;
           daa[7*20+2]=   0.0604932;
           daa[7*20+3]=   10.3969000;
           daa[7*20+4]=   0.0489798;
           daa[7*20+5]=   0.0604932;
           daa[7*20+6]=   14.7801000;
           daa[8*20+0]=   0.0050000;
           daa[8*20+1]=   6.8440500;
           daa[8*20+2]=   8.5987600;
           daa[8*20+3]=   2.3177900;
           daa[8*20+4]=   0.0050000;
           daa[8*20+5]=   18.5465000;
           daa[8*20+6]=   0.0050000;
           daa[8*20+7]=   0.0050000;
           daa[9*20+0]=   0.0050000;
           daa[9*20+1]=   1.3406900;
           daa[9*20+2]=   0.9870280;
           daa[9*20+3]=   0.1451240;
           daa[9*20+4]=   0.0050000;
           daa[9*20+5]=   0.0342252;
           daa[9*20+6]=   0.0390512;
           daa[9*20+7]=   0.0050000;
           daa[9*20+8]=   0.0050000;
           daa[10*20+0]=  0.1602400;
           daa[10*20+1]=  0.5867570;
           daa[10*20+2]=  0.0050000;
           daa[10*20+3]=  0.0050000;
           daa[10*20+4]=  0.0050000;
           daa[10*20+5]=  2.8904800;
           daa[10*20+6]=  0.1298390;
           daa[10*20+7]=  0.0489798;
           daa[10*20+8]=  1.7638200;
           daa[10*20+9]=  9.1024600;
           daa[11*20+0]=  0.5927840;
           daa[11*20+1]=  39.8897000;
           daa[11*20+2]=  10.6655000;
           daa[11*20+3]=  0.8943130;
           daa[11*20+4]=  0.0050000;
           daa[11*20+5]=  13.0705000;
           daa[11*20+6]=  23.9626000;
           daa[11*20+7]=  0.2794250;
           daa[11*20+8]=  0.2240600;
           daa[11*20+9]=  0.8174810;
           daa[11*20+10]= 0.0050000;
           daa[12*20+0]=  0.0050000;
           daa[12*20+1]=  3.2865200;
           daa[12*20+2]=  0.2015260;
           daa[12*20+3]=  0.0050000;
           daa[12*20+4]=  0.0050000;
           daa[12*20+5]=  0.0050000;
           daa[12*20+6]=  0.0050000;
           daa[12*20+7]=  0.0489798;
           daa[12*20+8]=  0.0050000;
           daa[12*20+9]=  17.3064000;
           daa[12*20+10]= 11.3839000;
           daa[12*20+11]= 4.0956400;
           daa[13*20+0]=  0.5979230;
           daa[13*20+1]=  0.0050000;
           daa[13*20+2]=  0.0050000;
           daa[13*20+3]=  0.0050000;
           daa[13*20+4]=  0.3629590;
           daa[13*20+5]=  0.0050000;
           daa[13*20+6]=  0.0050000;
           daa[13*20+7]=  0.0050000;
           daa[13*20+8]=  0.0050000;
           daa[13*20+9]=  1.4828800;
           daa[13*20+10]= 7.4878100;
           daa[13*20+11]= 0.0050000;
           daa[13*20+12]= 0.0050000;
           daa[14*20+0]=  1.0098100;
           daa[14*20+1]=  0.4047230;
           daa[14*20+2]=  0.3448480;
           daa[14*20+3]=  0.0050000;
           daa[14*20+4]=  0.0050000;
           daa[14*20+5]=  3.0450200;
           daa[14*20+6]=  0.0050000;
           daa[14*20+7]=  0.0050000;
           daa[14*20+8]=  13.9444000;
           daa[14*20+9]=  0.0050000;
           daa[14*20+10]= 9.8309500;
           daa[14*20+11]= 0.1119280;
           daa[14*20+12]= 0.0050000;
           daa[14*20+13]= 0.0342252;
           daa[15*20+0]=  8.5942000;
           daa[15*20+1]=  8.3502400;
           daa[15*20+2]=  14.5699000;
           daa[15*20+3]=  0.4278810;
           daa[15*20+4]=  1.1219500;
           daa[15*20+5]=  0.1602400;
           daa[15*20+6]=  0.0050000;
           daa[15*20+7]=  6.2796600;
           daa[15*20+8]=  0.7251570;
           daa[15*20+9]=  0.7400910;
           daa[15*20+10]= 6.1439600;
           daa[15*20+11]= 0.0050000;
           daa[15*20+12]= 0.3925750;
           daa[15*20+13]= 4.2793900;
           daa[15*20+14]= 14.2490000;
           daa[16*20+0]=  24.1422000;
           daa[16*20+1]=  0.9282030;
           daa[16*20+2]=  4.5420600;
           daa[16*20+3]=  0.6303950;
           daa[16*20+4]=  0.0050000;
           daa[16*20+5]=  0.2030910;
           daa[16*20+6]=  0.4587430;
           daa[16*20+7]=  0.0489798;
           daa[16*20+8]=  0.9595600;
           daa[16*20+9]=  9.3634500;
           daa[16*20+10]= 0.0050000;
           daa[16*20+11]= 4.0480200;
           daa[16*20+12]= 7.4131300;
           daa[16*20+13]= 0.1145120;
           daa[16*20+14]= 4.3370100;
           daa[16*20+15]= 6.3407900;
           daa[17*20+0]=  0.0050000;
           daa[17*20+1]=  5.9656400;
           daa[17*20+2]=  0.0050000;
           daa[17*20+3]=  0.0050000;
           daa[17*20+4]=  5.4989400;
           daa[17*20+5]=  0.0443298;
           daa[17*20+6]=  0.0050000;
           daa[17*20+7]=  2.8258000;
           daa[17*20+8]=  0.0050000;
           daa[17*20+9]=  0.0050000;
           daa[17*20+10]= 1.3703100;
           daa[17*20+11]= 0.0050000;
           daa[17*20+12]= 0.0050000;
           daa[17*20+13]= 0.0050000;
           daa[17*20+14]= 0.0050000;
           daa[17*20+15]= 1.1015600;
           daa[17*20+16]= 0.0050000;
           daa[18*20+0]=  0.0050000;
           daa[18*20+1]=  0.0050000;
           daa[18*20+2]=  5.0647500;
           daa[18*20+3]=  2.2815400;
           daa[18*20+4]=  8.3483500;
           daa[18*20+5]=  0.0050000;
           daa[18*20+6]=  0.0050000;
           daa[18*20+7]=  0.0050000;
           daa[18*20+8]=  47.4889000;
           daa[18*20+9]=  0.1145120;
           daa[18*20+10]= 0.0050000;
           daa[18*20+11]= 0.0050000;
           daa[18*20+12]= 0.5791980;
           daa[18*20+13]= 4.1272800;
           daa[18*20+14]= 0.0050000;
           daa[18*20+15]= 0.9331420;
           daa[18*20+16]= 0.4906080;
           daa[18*20+17]= 0.0050000;
           daa[19*20+0]=  24.8094000;
           daa[19*20+1]=  0.2794250;
           daa[19*20+2]=  0.0744808;
           daa[19*20+3]=  2.9178600;
           daa[19*20+4]=  0.0050000;
           daa[19*20+5]=  0.0050000;
           daa[19*20+6]=  2.1995200;
           daa[19*20+7]=  2.7962200;
           daa[19*20+8]=  0.8274790;
           daa[19*20+9]=  24.8231000;
           daa[19*20+10]= 2.9534400;
           daa[19*20+11]= 0.1280650;
           daa[19*20+12]= 14.7683000;
           daa[19*20+13]= 2.2800000;
           daa[19*20+14]= 0.0050000;
           daa[19*20+15]= 0.8626370;
           daa[19*20+16]= 0.0050000;
           daa[19*20+17]= 0.0050000;
           daa[19*20+18]= 1.3548200;
           
           /*f[0]=  0.038;
           f[1]=  0.057;
           f[2]=  0.089;
           f[3]=  0.034;
           f[4]=  0.024;
           f[5]=  0.044;
           f[6]=  0.062;
           f[7]=  0.084;
           f[8]=  0.016;
           f[9]=  0.098;
           f[10]= 0.058;
           f[11]= 0.064;
           f[12]= 0.016;
           f[13]= 0.042;
           f[14]= 0.046;
           f[15]= 0.055;
           f[16]= 0.081;
           f[17]= 0.020;
           f[18]= 0.021;
           f[19]= 0.051;*/

	   f[0]= 0.0377494;             f[1]= 0.057321;              f[2]= 0.0891129;             f[3]= 0.0342034;
           f[4]= 0.0240105;             f[5]= 0.0437824;             f[6]= 0.0618606;             f[7]= 0.0838496;
           f[8]= 0.0156076;             f[9]= 0.0983641;             f[10]= 0.0577867;            f[11]= 0.0641682;
           f[12]= 0.0158419;            f[13]= 0.0422741;            f[14]= 0.0458601;            f[15]= 0.0550846;
           f[16]= 0.0813774;            f[17]= 0.019597;             f[18]= 0.0205847;            f[19]= 0.0515638;

	  }
	  break;
	case JTTDCMUT:
	  {
           daa[1*20+0]=   0.531678;
           daa[2*20+0]=   0.557967;
           daa[2*20+1]=   0.451095;
           daa[3*20+0]=   0.827445;
           daa[3*20+1]=   0.154899;
           daa[3*20+2]=   5.549530;
           daa[4*20+0]=   0.574478;
           daa[4*20+1]=   1.019843;
           daa[4*20+2]=   0.313311;
           daa[4*20+3]=   0.105625;
           daa[5*20+0]=   0.556725;
           daa[5*20+1]=   3.021995;
           daa[5*20+2]=   0.768834;
           daa[5*20+3]=   0.521646;
           daa[5*20+4]=   0.091304;
           daa[6*20+0]=   1.066681;
           daa[6*20+1]=   0.318483;
           daa[6*20+2]=   0.578115;
           daa[6*20+3]=   7.766557;
           daa[6*20+4]=   0.053907;
           daa[6*20+5]=   3.417706;
           daa[7*20+0]=   1.740159;
           daa[7*20+1]=   1.359652;
           daa[7*20+2]=   0.773313;
           daa[7*20+3]=   1.272434;
           daa[7*20+4]=   0.546389;
           daa[7*20+5]=   0.231294;
           daa[7*20+6]=   1.115632;
           daa[8*20+0]=   0.219970;
           daa[8*20+1]=   3.210671;
           daa[8*20+2]=   4.025778;
           daa[8*20+3]=   1.032342;
           daa[8*20+4]=   0.724998;
           daa[8*20+5]=   5.684080;
           daa[8*20+6]=   0.243768;
           daa[8*20+7]=   0.201696;
           daa[9*20+0]=   0.361684;
           daa[9*20+1]=   0.239195;
           daa[9*20+2]=   0.491003;
           daa[9*20+3]=   0.115968;
           daa[9*20+4]=   0.150559;
           daa[9*20+5]=   0.078270;
           daa[9*20+6]=   0.111773;
           daa[9*20+7]=   0.053769;
           daa[9*20+8]=   0.181788;
           daa[10*20+0]=  0.310007;
           daa[10*20+1]=  0.372261;
           daa[10*20+2]=  0.137289;
           daa[10*20+3]=  0.061486;
           daa[10*20+4]=  0.164593;
           daa[10*20+5]=  0.709004;
           daa[10*20+6]=  0.097485;
           daa[10*20+7]=  0.069492;
           daa[10*20+8]=  0.540571;
           daa[10*20+9]=  2.335139;
           daa[11*20+0]=  0.369437;
           daa[11*20+1]=  6.529255;
           daa[11*20+2]=  2.529517;
           daa[11*20+3]=  0.282466;
           daa[11*20+4]=  0.049009;
           daa[11*20+5]=  2.966732;
           daa[11*20+6]=  1.731684;
           daa[11*20+7]=  0.269840;
           daa[11*20+8]=  0.525096;
           daa[11*20+9]=  0.202562;
           daa[11*20+10]= 0.146481;
           daa[12*20+0]=  0.469395;
           daa[12*20+1]=  0.431045;
           daa[12*20+2]=  0.330720;
           daa[12*20+3]=  0.190001;
           daa[12*20+4]=  0.409202;
           daa[12*20+5]=  0.456901;
           daa[12*20+6]=  0.175084;
           daa[12*20+7]=  0.130379;
           daa[12*20+8]=  0.329660;
           daa[12*20+9]=  4.831666;
           daa[12*20+10]= 3.856906;
           daa[12*20+11]= 0.624581;
           daa[13*20+0]=  0.138293;
           daa[13*20+1]=  0.065314;
           daa[13*20+2]=  0.073481;
           daa[13*20+3]=  0.032522;
           daa[13*20+4]=  0.678335;
           daa[13*20+5]=  0.045683;
           daa[13*20+6]=  0.043829;
           daa[13*20+7]=  0.050212;
           daa[13*20+8]=  0.453428;
           daa[13*20+9]=  0.777090;
           daa[13*20+10]= 2.500294;
           daa[13*20+11]= 0.024521;
           daa[13*20+12]= 0.436181;
           daa[14*20+0]=  1.959599;
           daa[14*20+1]=  0.710489;
           daa[14*20+2]=  0.121804;
           daa[14*20+3]=  0.127164;
           daa[14*20+4]=  0.123653;
           daa[14*20+5]=  1.608126;
           daa[14*20+6]=  0.191994;
           daa[14*20+7]=  0.208081;
           daa[14*20+8]=  1.141961;
           daa[14*20+9]=  0.098580;
           daa[14*20+10]= 1.060504;
           daa[14*20+11]= 0.216345;
           daa[14*20+12]= 0.164215;
           daa[14*20+13]= 0.148483;
           daa[15*20+0]=  3.887095;
           daa[15*20+1]=  1.001551;
           daa[15*20+2]=  5.057964;
           daa[15*20+3]=  0.589268;
           daa[15*20+4]=  2.155331;
           daa[15*20+5]=  0.548807;
           daa[15*20+6]=  0.312449;
           daa[15*20+7]=  1.874296;
           daa[15*20+8]=  0.743458;
           daa[15*20+9]=  0.405119;
           daa[15*20+10]= 0.592511;
           daa[15*20+11]= 0.474478;
           daa[15*20+12]= 0.285564;
           daa[15*20+13]= 0.943971;
           daa[15*20+14]= 2.788406;
           daa[16*20+0]=  4.582565;
           daa[16*20+1]=  0.650282;
           daa[16*20+2]=  2.351311;
           daa[16*20+3]=  0.425159;
           daa[16*20+4]=  0.469823;
           daa[16*20+5]=  0.523825;
           daa[16*20+6]=  0.331584;
           daa[16*20+7]=  0.316862;
           daa[16*20+8]=  0.477355;
           daa[16*20+9]=  2.553806;
           daa[16*20+10]= 0.272514;
           daa[16*20+11]= 0.965641;
           daa[16*20+12]= 2.114728;
           daa[16*20+13]= 0.138904;
           daa[16*20+14]= 1.176961;
           daa[16*20+15]= 4.777647;
           daa[17*20+0]=  0.084329;
           daa[17*20+1]=  1.257961;
           daa[17*20+2]=  0.027700;
           daa[17*20+3]=  0.057466;
           daa[17*20+4]=  1.104181;
           daa[17*20+5]=  0.172206;
           daa[17*20+6]=  0.114381;
           daa[17*20+7]=  0.544180;
           daa[17*20+8]=  0.128193;
           daa[17*20+9]=  0.134510;
           daa[17*20+10]= 0.530324;
           daa[17*20+11]= 0.089134;
           daa[17*20+12]= 0.201334;
           daa[17*20+13]= 0.537922;
           daa[17*20+14]= 0.069965;
           daa[17*20+15]= 0.310927;
           daa[17*20+16]= 0.080556;
           daa[18*20+0]=  0.139492;
           daa[18*20+1]=  0.235601;
           daa[18*20+2]=  0.700693;
           daa[18*20+3]=  0.453952;
           daa[18*20+4]=  2.114852;
           daa[18*20+5]=  0.254745;
           daa[18*20+6]=  0.063452;
           daa[18*20+7]=  0.052500;
           daa[18*20+8]=  5.848400;
           daa[18*20+9]=  0.303445;
           daa[18*20+10]= 0.241094;
           daa[18*20+11]= 0.087904;
           daa[18*20+12]= 0.189870;
           daa[18*20+13]= 5.484236;
           daa[18*20+14]= 0.113850;
           daa[18*20+15]= 0.628608;
           daa[18*20+16]= 0.201094;
           daa[18*20+17]= 0.747889;
           daa[19*20+0]=  2.924161;
           daa[19*20+1]=  0.171995;
           daa[19*20+2]=  0.164525;
           daa[19*20+3]=  0.315261;
           daa[19*20+4]=  0.621323;
           daa[19*20+5]=  0.179771;
           daa[19*20+6]=  0.465271;
           daa[19*20+7]=  0.470140;
           daa[19*20+8]=  0.121827;
           daa[19*20+9]=  9.533943;
           daa[19*20+10]= 1.761439;
           daa[19*20+11]= 0.124066;
           daa[19*20+12]= 3.038533;
           daa[19*20+13]= 0.593478;
           daa[19*20+14]= 0.211561;
           daa[19*20+15]= 0.408532;
           daa[19*20+16]= 1.143980;
           daa[19*20+17]= 0.239697;
           daa[19*20+18]= 0.165473;
           
           f[0]=  0.077;
           f[1]=  0.051;
           f[2]=  0.043;
           f[3]=  0.051;
           f[4]=  0.020;
           f[5]=  0.041;
           f[6]=  0.062;
           f[7]=  0.075;
           f[8]=  0.023;
           f[9]=  0.053;
           f[10]= 0.091;
           f[11]= 0.059;
           f[12]= 0.024;
           f[13]= 0.040;
           f[14]= 0.051;
           f[15]= 0.068;
           f[16]= 0.059;
           f[17]= 0.014;
           f[18]= 0.032;
           f[19]= 0.066;
	  }
	  break;
	case FLU:
	  {
	    daa[ 1*20+ 0] 	=	0.138658765	;
	    daa[ 2*20+ 0] 	=	0.053366579	;
	    daa[ 2*20+ 1] 	=	0.161000889	;
	    daa[ 3*20+ 0] 	=	0.584852306	;
	    daa[ 3*20+ 1] 	=	0.006771843	;
	    daa[ 3*20+ 2] 	=	7.737392871	;
	    daa[ 4*20+ 0] 	=	0.026447095	;
	    daa[ 4*20+ 1] 	=	0.167207008	;
	    daa[ 4*20+ 2] 	=	1.30E-05	;
	    daa[ 4*20+ 3] 	=	1.41E-02	;
	    daa[ 5*20+ 0] 	=	0.353753982	;
	    daa[ 5*20+ 1] 	=	3.292716942	;
	    daa[ 5*20+ 2] 	=	0.530642655	;
	    daa[ 5*20+ 3] 	=	0.145469388	;
	    daa[ 5*20+ 4] 	=	0.002547334	;
	    daa[ 6*20+ 0] 	=	1.484234503	;
	    daa[ 6*20+ 1] 	=	0.124897617	;
	    daa[ 6*20+ 2] 	=	0.061652192	;
	    daa[ 6*20+ 3] 	=	5.370511279	;
	    daa[ 6*20+ 4] 	=	3.91E-11	;
	    daa[ 6*20+ 5] 	=	1.195629122	;
	    daa[ 7*20+ 0] 	=	1.132313122	;
	    daa[ 7*20+ 1] 	=	1.190624465	;
	    daa[ 7*20+ 2] 	=	0.322524648	;
	    daa[ 7*20+ 3] 	=	1.934832784	;
	    daa[ 7*20+ 4] 	=	0.116941459	;
	    daa[ 7*20+ 5] 	=	0.108051341	;
	    daa[ 7*20+ 6] 	=	1.593098825	;
	    daa[ 8*20+ 0] 	=	0.214757862	;
	    daa[ 8*20+ 1] 	=	1.879569938	;
	    daa[ 8*20+ 2] 	=	1.387096032	;
	    daa[ 8*20+ 3] 	=	0.887570549	;
	    daa[ 8*20+ 4] 	=	2.18E-02	;
	    daa[ 8*20+ 5] 	=	5.330313412	;
	    daa[ 8*20+ 6] 	=	0.256491863	;
	    daa[ 8*20+ 7] 	=	0.058774527	;
	    daa[ 9*20+ 0] 	=	0.149926734	;
	    daa[ 9*20+ 1] 	=	0.246117172	;
	    daa[ 9*20+ 2] 	=	0.218571975	;
	    daa[ 9*20+ 3] 	=	0.014085917	;
	    daa[ 9*20+ 4] 	=	0.001112158	;
	    daa[ 9*20+ 5] 	=	0.02883995	;
	    daa[ 9*20+ 6] 	=	1.42E-02	;
	    daa[ 9*20+ 7] 	=	1.63E-05	;
	    daa[ 9*20+ 8] 	=	0.243190142	;
	    daa[10*20+ 0] 	=	0.023116952	;
	    daa[10*20+ 1] 	=	0.296045557	;
	    daa[10*20+ 2] 	=	8.36E-04	;
	    daa[10*20+ 3] 	=	0.005730682	;
	    daa[10*20+ 4] 	=	0.005613627	;
	    daa[10*20+ 5] 	=	1.020366955	;
	    daa[10*20+ 6] 	=	0.016499536	;
	    daa[10*20+ 7] 	=	0.006516229	;
	    daa[10*20+ 8] 	=	0.321611694	;
	    daa[10*20+ 9] 	=	3.512072282	;
	    daa[11*20+ 0] 	=	0.47433361	;
	    daa[11*20+ 1] 	=	15.30009662	;
	    daa[11*20+ 2] 	=	2.646847965	;
	    daa[11*20+ 3] 	=	0.29004298	;
	    daa[11*20+ 4] 	=	3.83E-06	;
	    daa[11*20+ 5] 	=	2.559587177	;
	    daa[11*20+ 6] 	=	3.881488809	;
	    daa[11*20+ 7] 	=	0.264148929	;
	    daa[11*20+ 8] 	=	0.347302791	;
	    daa[11*20+ 9] 	=	0.227707997	;
	    daa[11*20+10] 	=	0.129223639	;
	    daa[12*20+ 0] 	=	0.058745423	;
	    daa[12*20+ 1] 	=	0.890162346	;
	    daa[12*20+ 2] 	=	0.005251688	;
	    daa[12*20+ 3] 	=	0.041762964	;
	    daa[12*20+ 4] 	=	0.11145731	;
	    daa[12*20+ 5] 	=	0.190259181	;
	    daa[12*20+ 6] 	=	0.313974351	;
	    daa[12*20+ 7] 	=	0.001500467	;
	    daa[12*20+ 8] 	=	0.001273509	;
	    daa[12*20+ 9] 	=	9.017954203	;
	    daa[12*20+10] 	=	6.746936485	;
	    daa[12*20+11] 	=	1.331291619	;
	    daa[13*20+ 0] 	=	0.080490909	;
	    daa[13*20+ 1] 	=	1.61E-02	;
	    daa[13*20+ 2] 	=	8.36E-04	;
	    daa[13*20+ 3] 	=	1.06E-06	;
	    daa[13*20+ 4] 	=	0.104053666	;
	    daa[13*20+ 5] 	=	0.032680657	;
	    daa[13*20+ 6] 	=	0.001003501	;
	    daa[13*20+ 7] 	=	0.001236645	;
	    daa[13*20+ 8] 	=	0.119028506	;
	    daa[13*20+ 9] 	=	1.463357278	;
	    daa[13*20+10] 	=	2.986800036	;
	    daa[13*20+11] 	=	3.20E-01	;
	    daa[13*20+12] 	=	0.279910509	;
	    daa[14*20+ 0] 	=	0.659311478	;
	    daa[14*20+ 1] 	=	0.15402718	;
	    daa[14*20+ 2] 	=	3.64E-02	;
	    daa[14*20+ 3] 	=	0.188539456	;
	    daa[14*20+ 4] 	=	1.59E-13	;
	    daa[14*20+ 5] 	=	0.712769599	;
	    daa[14*20+ 6] 	=	0.319558828	;
	    daa[14*20+ 7] 	=	0.038631761	;
	    daa[14*20+ 8] 	=	0.924466914	;
	    daa[14*20+ 9] 	=	0.080543327	;
	    daa[14*20+10] 	=	0.634308521	;
	    daa[14*20+11] 	=	0.195750632	;
	    daa[14*20+12] 	=	5.69E-02	;
	    daa[14*20+13] 	=	0.00713243	;
	    daa[15*20+ 0] 	=	3.011344519	;
	    daa[15*20+ 1] 	=	0.95013841	;
	    daa[15*20+ 2] 	=	3.881310531	;
	    daa[15*20+ 3] 	=	0.338372183	;
	    daa[15*20+ 4] 	=	0.336263345	;
	    daa[15*20+ 5] 	=	0.487822499	;
	    daa[15*20+ 6] 	=	0.307140298	;
	    daa[15*20+ 7] 	=	1.585646577	;
	    daa[15*20+ 8] 	=	0.58070425	;
	    daa[15*20+ 9] 	=	0.290381075	;
	    daa[15*20+10] 	=	0.570766693	;
	    daa[15*20+11] 	=	0.283807672	;
	    daa[15*20+12] 	=	0.007026588	;
	    daa[15*20+13] 	=	0.99668567	;
	    daa[15*20+14] 	=	2.087385344	;
	    daa[16*20+ 0] 	=	5.418298175	;
	    daa[16*20+ 1] 	=	0.183076905	;
	    daa[16*20+ 2] 	=	2.140332316	;
	    daa[16*20+ 3] 	=	0.135481233	;
	    daa[16*20+ 4] 	=	0.011975266	;
	    daa[16*20+ 5] 	=	0.602340963	;
	    daa[16*20+ 6] 	=	0.280124895	;
	    daa[16*20+ 7] 	=	0.01880803	;
	    daa[16*20+ 8] 	=	0.368713573	;
	    daa[16*20+ 9] 	=	2.904052286	;
	    daa[16*20+10] 	=	0.044926357	;
	    daa[16*20+11] 	=	1.5269642	;
	    daa[16*20+12] 	=	2.031511321	;
	    daa[16*20+13] 	=	0.000134906	;
	    daa[16*20+14] 	=	0.542251094	;
	    daa[16*20+15] 	=	2.206859934	;
	    daa[17*20+ 0] 	=	1.96E-01	;
	    daa[17*20+ 1] 	=	1.369429408	;
	    daa[17*20+ 2] 	=	5.36E-04	;
	    daa[17*20+ 3] 	=	1.49E-05	;
	    daa[17*20+ 4] 	=	0.09410668	;
	    daa[17*20+ 5] 	=	4.40E-02	;
	    daa[17*20+ 6] 	=	0.155245492	;
	    daa[17*20+ 7] 	=	0.196486447	;
	    daa[17*20+ 8] 	=	2.24E-02	;
	    daa[17*20+ 9] 	=	0.03213215	;
	    daa[17*20+10] 	=	0.431277663	;
	    daa[17*20+11] 	=	4.98E-05	;
	    daa[17*20+12] 	=	0.070460039	;
	    daa[17*20+13] 	=	0.814753094	;
	    daa[17*20+14] 	=	0.000431021	;
	    daa[17*20+15] 	=	0.099835753	;
	    daa[17*20+16] 	=	0.207066206	;
	    daa[18*20+ 0] 	=	0.018289288	;
	    daa[18*20+ 1] 	=	0.099855497	;
	    daa[18*20+ 2] 	=	0.373101927	;
	    daa[18*20+ 3] 	=	0.525398543	;
	    daa[18*20+ 4] 	=	0.601692431	;
	    daa[18*20+ 5] 	=	0.072205935	;
	    daa[18*20+ 6] 	=	0.10409287	;
	    daa[18*20+ 7] 	=	0.074814997	;
	    daa[18*20+ 8] 	=	6.448954446	;
	    daa[18*20+ 9] 	=	0.273934263	;
	    daa[18*20+10] 	=	0.340058468	;
	    daa[18*20+11] 	=	0.012416222	;
	    daa[18*20+12] 	=	0.874272175	;
	    daa[18*20+13] 	=	5.393924245	;
	    daa[18*20+14] 	=	1.82E-04	;
	    daa[18*20+15] 	=	0.39255224	;
	    daa[18*20+16] 	=	0.12489802	;
	    daa[18*20+17] 	=	0.42775543	;
	    daa[19*20+ 0] 	=	3.53200527	;
	    daa[19*20+ 1] 	=	0.103964386	;
	    daa[19*20+ 2] 	=	0.010257517	;
	    daa[19*20+ 3] 	=	0.297123975	;
	    daa[19*20+ 4] 	=	0.054904564	;
	    daa[19*20+ 5] 	=	0.406697814	;
	    daa[19*20+ 6] 	=	0.285047948	;
	    daa[19*20+ 7] 	=	0.337229619	;
	    daa[19*20+ 8] 	=	0.098631355	;
	    daa[19*20+ 9] 	=	14.39405219	;
	    daa[19*20+10] 	=	0.890598579	;
	    daa[19*20+11] 	=	0.07312793	;
	    daa[19*20+12] 	=	4.904842235	;
	    daa[19*20+13] 	=	0.592587985	;
	    daa[19*20+14] 	=	0.058971975	;
	    daa[19*20+15] 	=	0.088256423	;
	    daa[19*20+16] 	=	0.654109108	;
	    daa[19*20+17] 	=	0.256900461	;
	    daa[19*20+18] 	=	0.167581647	;
	    
 
  
	    f[0]	=	0.0471	;
	    f[1]	=	0.0509	;
	    f[2]	=	0.0742	;
	    f[3]	=	0.0479	;
	    f[4]	=	0.0250	;
	    f[5]	=	0.0333	;
	    f[6]	=	0.0546	;
	    f[7]	=	0.0764	;
	    f[8]	=	0.0200	;
	    f[9]	=	0.0671	;
	    f[10]	=	0.0715	;
	    f[11]	=	0.0568	;
	    f[12]	=	0.0181	;
	    f[13]	=	0.0305	;
	    f[14]	=	0.0507	;
	    f[15]	=	0.0884	;
	    f[16]	=	0.0743	;
	    f[17]	=	0.0185	;
	    f[18]	=	0.0315	;
	    f[19]	=	0.0632	;
	  }
	  break;  
      case LG4:     
	{
	  double 
	    rates[4][190] = 
	    {
	      {
		0.269343
		, 0.254612, 0.150988
		, 0.236821, 0.031863, 0.659648
		, 2.506547, 0.938594, 0.975736, 0.175533
		, 0.359080, 0.348288, 0.697708, 0.086573, 0.095967
		, 0.304674, 0.156000, 0.377704, 0.449140, 0.064706, 4.342595
		, 1.692015, 0.286638, 0.565095, 0.380358, 0.617945, 0.202058, 0.264342
		, 0.251974, 0.921633, 1.267609, 0.309692, 0.390429, 2.344059, 0.217750, 0.104842
		, 1.085220, 0.325624, 0.818658, 0.037814, 1.144150, 0.534567, 0.222793, 0.062682, 0.567431
		, 0.676353, 0.602366, 0.217027, 0.007533, 1.595775, 0.671143, 0.158424, 0.070463, 0.764255, 8.226528
		, 0.179155, 0.971338, 1.343718, 0.133744, 0.122468, 0.983857, 0.994128, 0.220916, 0.410581, 0.387487, 0.181110
		, 1.636817, 0.515217, 0.670461, 0.071252, 1.534848, 5.288642, 0.255628, 0.094198, 0.257229, 25.667158, 6.819689, 1.591212
		, 0.235498, 0.123932, 0.099793, 0.030425, 0.897279, 0.112229, 0.022529, 0.047488, 0.762914, 1.344259, 0.865691, 0.038921, 2.030833
		, 1.265605, 0.040163, 0.173354, 0.027579, 0.259961, 0.580374, 0.088041, 0.145595, 0.143676, 0.298859, 1.020117, 0.000714, 0.190019, 0.093964
		, 5.368405, 0.470952, 5.267140, 0.780505, 4.986071, 0.890554, 0.377949, 1.755515, 0.786352, 0.527246, 0.667783, 0.659948, 0.731921, 0.837669, 1.355630
		, 1.539394, 0.326789, 1.688169, 0.283738, 1.389282, 0.329821, 0.231770, 0.117017, 0.449977, 3.531600, 0.721586, 0.497588, 2.691697, 0.152088, 0.698040, 16.321298
		, 0.140944, 0.375611, 0.025163, 0.002757, 0.801456, 0.257253, 0.103678, 0.132995, 0.345834, 0.377156, 0.839647, 0.176970, 0.505682, 1.670170, 0.091298, 0.210096, 0.013165
		, 0.199836, 0.146857, 0.806275, 0.234246, 1.436970, 0.319669, 0.010076, 0.036859, 3.503317, 0.598632, 0.738969, 0.154436, 0.579000, 4.245524, 0.074524, 0.454195, 0.232913, 1.178490
		, 9.435529, 0.285934, 0.395670, 0.130890, 6.097263, 0.516259, 0.503665, 0.222960, 0.149143, 13.666175, 2.988174, 0.162725, 5.973826, 0.843416, 0.597394, 0.701149, 4.680002, 0.300085, 0.416262
	      },
	      {
		0.133720
		, 0.337212, 0.749052
		, 0.110918, 0.105087, 4.773487
		, 3.993460, 0.188305, 1.590332, 0.304942
		, 0.412075, 2.585774, 1.906884, 0.438367, 0.242076
		, 0.435295, 0.198278, 0.296366, 7.470333, 0.008443, 3.295515
		, 7.837540, 0.164607, 0.431724, 0.153850, 1.799716, 0.269744, 0.242866
		, 0.203872, 2.130334, 9.374479, 1.080878, 0.152458, 12.299133, 0.279589, 0.089714
		, 0.039718, 0.024553, 0.135254, 0.014979, 0.147498, 0.033964, 0.005585, 0.007248, 0.022746
		, 0.075784, 0.080091, 0.084971, 0.014128, 0.308347, 0.500836, 0.022833, 0.022999, 0.161270, 1.511682
		, 0.177662, 10.373708, 1.036721, 0.038303, 0.043030, 2.181033, 0.321165, 0.103050, 0.459502, 0.021215, 0.078395
		, 0.420784, 0.192765, 0.329545, 0.008331, 0.883142, 1.403324, 0.168673, 0.160728, 0.612573, 1.520889, 7.763266, 0.307903
		, 0.071268, 0.019652, 0.088753, 0.013547, 0.566609, 0.071878, 0.020050, 0.041022, 0.625361, 0.382806, 1.763059, 0.044644, 1.551911
		, 0.959127, 1.496585, 0.377794, 0.332010, 0.318192, 1.386970, 0.915904, 0.224255, 2.611479, 0.029351, 0.068250, 1.542356, 0.047525, 0.182715
		, 11.721512, 0.359408, 2.399158, 0.219464, 9.104192, 0.767563, 0.235229, 3.621219, 0.971955, 0.033780, 0.043035, 0.236929, 0.319964, 0.124977, 0.840651
		, 2.847068, 0.218463, 1.855386, 0.109808, 4.347048, 0.765848, 0.164569, 0.312024, 0.231569, 0.356327, 0.159597, 0.403210, 1.135162, 0.106903, 0.269190, 9.816481
		, 0.030203, 0.387292, 0.118878, 0.067287, 0.190240, 0.122113, 0.007023, 0.137411, 0.585141, 0.020634, 0.228824, 0.000122, 0.474862, 3.135128, 0.030313, 0.093830, 0.119152
		, 0.067183, 0.130101, 0.348730, 0.061798, 0.301198, 0.095382, 0.095764, 0.044628, 2.107384, 0.046105, 0.100117, 0.017073, 0.192383, 8.367641, 0.000937, 0.137416, 0.044722, 4.179782
		, 0.679398, 0.041567, 0.092408, 0.023701, 1.271187, 0.115566, 0.055277, 0.086988, 0.060779, 8.235167, 0.609420, 0.061764, 0.581962, 0.184187, 0.080246, 0.098033, 1.438350, 0.023439, 0.039124
	      },	    
	      {
		0.421017
		, 0.316236, 0.693340
		, 0.285984, 0.059926, 6.158219
		, 4.034031, 1.357707, 0.708088, 0.063669
		, 0.886972, 2.791622, 1.701830, 0.484347, 0.414286
		, 0.760525, 0.233051, 0.378723, 4.032667, 0.081977, 4.940411
		, 0.754103, 0.402894, 2.227443, 1.102689, 0.416576, 0.459376, 0.508409
		, 0.571422, 2.319453, 5.579973, 0.885376, 1.439275, 4.101979, 0.576745, 0.428799
		, 0.162152, 0.085229, 0.095692, 0.006129, 0.490937, 0.104843, 0.045514, 0.004705, 0.098934
		, 0.308006, 0.287051, 0.056994, 0.007102, 0.958988, 0.578990, 0.067119, 0.024403, 0.342983, 3.805528
		, 0.390161, 7.663209, 1.663641, 0.105129, 0.135029, 3.364474, 0.652618, 0.457702, 0.823674, 0.129858, 0.145630
		, 1.042298, 0.364551, 0.293222, 0.037983, 1.486520, 1.681752, 0.192414, 0.070498, 0.222626, 4.529623, 4.781730, 0.665308
		, 0.362476, 0.073439, 0.129245, 0.020078, 1.992483, 0.114549, 0.023272, 0.064490, 1.491794, 1.113437, 2.132006, 0.041677, 1.928654
		, 1.755491, 0.087050, 0.099325, 0.163817, 0.242851, 0.322939, 0.062943, 0.198698, 0.192904, 0.062948, 0.180283, 0.059655, 0.129323, 0.065778
		, 3.975060, 0.893398, 5.496314, 1.397313, 3.575120, 1.385297, 0.576191, 1.733288, 1.021255, 0.065131, 0.129115, 0.600308, 0.387276, 0.446001, 1.298493
		, 2.565079, 0.534056, 2.143993, 0.411388, 2.279084, 0.893006, 0.528209, 0.135731, 0.518741, 0.972662, 0.280700, 0.890086, 1.828755, 0.189028, 0.563778, 7.788147
		, 0.283631, 0.497926, 0.075454, 0.043794, 1.335322, 0.308605, 0.140137, 0.150797, 1.409726, 0.119868, 0.818331, 0.080591, 1.066017, 3.754687, 0.073415, 0.435046, 0.197272
		, 0.242513, 0.199157, 0.472207, 0.085937, 2.039787, 0.262751, 0.084578, 0.032247, 7.762326, 0.153966, 0.299828, 0.117255, 0.438215, 14.506235, 0.089180, 0.352766, 0.215417, 5.054245
		, 2.795818, 0.107130, 0.060909, 0.029724, 2.986426, 0.197267, 0.196977, 0.044327, 0.116751, 7.144311, 1.848622, 0.118020, 1.999696, 0.705747, 0.272763, 0.096935, 1.820982, 0.217007, 0.172975
	      },
	      {
		0.576160
		, 0.567606, 0.498643
		, 0.824359, 0.050698, 3.301401
		, 0.822724, 4.529235, 1.291808, 0.101930
		, 1.254238, 2.169809, 1.427980, 0.449474, 0.868679
		, 1.218615, 0.154502, 0.411471, 3.172277, 0.050239, 2.138661
		, 1.803443, 0.604673, 2.125496, 1.276384, 1.598679, 0.502653, 0.479490
		, 0.516862, 2.874265, 4.845769, 0.719673, 3.825677, 4.040275, 0.292773, 0.596643
		, 0.180898, 0.444586, 0.550969, 0.023542, 2.349573, 0.370160, 0.142187, 0.016618, 0.500788
		, 0.452099, 0.866322, 0.201033, 0.026731, 2.813990, 1.645178, 0.135556, 0.072152, 1.168817, 5.696116
		, 0.664186, 2.902886, 2.101971, 0.127988, 0.200218, 2.505933, 0.759509, 0.333569, 0.623100, 0.547454, 0.363656
		, 0.864415, 0.835049, 0.632649, 0.079201, 2.105931, 1.633544, 0.216462, 0.252419, 0.665406, 7.994105, 11.751178, 1.096842
		, 0.324478, 0.208947, 0.280339, 0.041683, 4.788477, 0.107022, 0.067711, 0.171320, 3.324779, 2.965328, 5.133843, 0.084856, 4.042591
		, 1.073043, 0.173826, 0.041985, 0.270336, 0.121299, 0.351384, 0.228565, 0.225318, 0.376089, 0.058027, 0.390354, 0.214230, 0.058954, 0.126299
		, 3.837562, 0.884342, 4.571911, 0.942751, 6.592827, 1.080063, 0.465397, 3.137614, 1.119667, 0.362516, 0.602355, 0.716940, 0.506796, 1.444484, 1.432558
		, 2.106026, 0.750016, 2.323325, 0.335915, 1.654673, 1.194017, 0.617231, 0.318671, 0.801030, 4.455842, 0.580191, 1.384210, 3.522468, 0.473128, 0.432718, 5.716300
		, 0.163720, 0.818102, 0.072322, 0.068275, 3.305436, 0.373790, 0.054323, 0.476587, 1.100360, 0.392946, 1.703323, 0.085720, 1.725516, 5.436253, 0.053108, 0.498594, 0.231832
		, 0.241167, 0.302440, 1.055095, 0.246940, 9.741942, 0.249895, 0.129973, 0.052363, 11.542498, 1.047449, 1.319667, 0.139770, 1.330225, 26.562270, 0.046986, 0.737653, 0.313460, 5.165098
		, 1.824586, 0.435795, 0.179086, 0.091739, 3.609570, 0.649507, 0.656681, 0.225234, 0.473437, 19.897252, 3.001995, 0.452926, 3.929598, 1.692159, 0.370204, 0.373501, 3.329822, 0.326593, 0.860743
	      }
	    };

	  double
	    freqs[4][20] = 
	    {{0.082276,0.055172,0.043853,0.053484,0.018957,0.028152,0.046679,0.157817,0.033297,0.028284,0.054284,0.025275,0.023665,0.041874,0.063071,0.066501,0.065424,0.023837,0.038633,0.049465},
	     {0.120900,0.036460,0.026510,0.040410,0.015980,0.021132,0.025191,0.036369,0.015884,0.111029,0.162852,0.024820,0.028023,0.074058,0.012065,0.041963,0.039072,0.012666,0.040478,0.114137},
	     {0.072639,0.051691,0.038642,0.055580,0.009829,0.031374,0.048731,0.065283,0.023791,0.086640,0.120847,0.052177,0.026728,0.032589,0.039238,0.046748,0.053361,0.008024,0.037426,0.098662},
	     {0.104843,0.078835,0.043513,0.090498,0.002924,0.066163,0.151640,0.038843,0.022556,0.018383,0.038687,0.104462,0.010166,0.009089,0.066950,0.053667,0.049486,0.004409,0.012924,0.031963}};


	  makeAASubstMat(daa, f,  rates[lg4_index],  freqs[lg4_index]);

	  /*int 
	    i, 
	    j, 
	    r = 0;

	  for(i = 1; i < 20; i++)
	    for(j = 0; j < i; j++)
	      {
		daa[i * 20 + j] = rates[lg4_index][r];
		r++;
	      }
	  
	  assert(r == 190);

	  for(i = 0; i < 20; i++)
	  f[i] = freqs[lg4_index][i];	  */

	}      
	break;
      case LG4X:
	{
	   double 
	    rates[4][190] = 
	    {
	      {
		0.295719,
		0.067388, 0.448317,
		0.253712, 0.457483, 2.358429,
		1.029289, 0.576016, 0.251987, 0.189008,
		0.107964, 1.741924, 0.216561, 0.599450, 0.029955,
		0.514644, 0.736017, 0.503084, 109.901504, 0.084794, 4.117654,
		10.868848, 0.704334, 0.435271, 1.070052, 1.862626, 0.246260, 1.202023,
		0.380498, 5.658311, 4.873453, 5.229858, 0.553477, 6.508329, 1.634845, 0.404968,
		0.084223, 0.123387, 0.090748, 0.052764, 0.151733, 0.054187, 0.060194, 0.048984, 0.204296,
		0.086976, 0.221777, 0.033310, 0.021407, 0.230320, 0.195703, 0.069359, 0.069963, 0.504221, 1.495537,
		0.188789, 93.433377, 0.746537, 0.621146, 0.096955, 1.669092, 2.448827, 0.256662, 1.991533, 0.091940, 0.122332,
		0.286389, 0.382175, 0.128905, 0.081091, 0.352526, 0.810168, 0.232297, 0.228519, 0.655465, 1.994320, 3.256485, 0.457430,
		0.155567, 0.235965, 0.127321, 0.205164, 0.590018, 0.066081, 0.064822, 0.241077, 6.799829, 0.754940, 2.261319, 0.163849, 1.559944,
		1.671061, 6.535048, 0.904011, 5.164456, 0.386853, 2.437439, 3.537387, 4.320442, 11.291065, 0.170343, 0.848067, 5.260446, 0.426508, 0.438856,
		2.132922, 0.525521, 0.939733, 0.747330, 1.559564, 0.165666, 0.435384, 3.656545, 0.961142, 0.050315, 0.064441, 0.360946, 0.132547, 0.306683, 4.586081,
		0.529591, 0.303537, 0.435450, 0.308078, 0.606648, 0.106333, 0.290413, 0.290216, 0.448965, 0.372166, 0.102493, 0.389413, 0.498634, 0.109129, 2.099355, 3.634276,
		0.115551, 0.641259, 0.046646, 0.260889, 0.587531, 0.093417, 0.280695, 0.307466, 6.227274, 0.206332, 0.459041, 0.033291, 0.559069, 18.392863, 0.411347, 0.101797, 0.034710,
		0.102453, 0.289466, 0.262076, 0.185083, 0.592318, 0.035149, 0.105999, 0.096556, 20.304886, 0.097050, 0.133091, 0.115301, 0.264728, 66.647302, 0.476350, 0.148995, 0.063603, 20.561407,
		0.916683, 0.102065, 0.043986, 0.080708, 0.885230, 0.072549, 0.206603, 0.306067, 0.205944, 5.381403, 0.561215, 0.112593, 0.693307, 0.400021, 0.584622, 0.089177, 0.755865, 0.133790, 0.154902
	      },
	      {
		0.066142,
		0.590377, 0.468325,
		0.069930, 0.013688, 2.851667,
		9.850951, 0.302287, 3.932151, 0.146882,
		1.101363, 1.353957, 8.159169, 0.249672, 0.582670,
		0.150375, 0.028386, 0.219934, 0.560142, 0.005035, 3.054085,
		0.568586, 0.037750, 0.421974, 0.046719, 0.275844, 0.129551, 0.037250,
		0.051668, 0.262130, 2.468752, 0.106259, 0.098208, 4.210126, 0.029788, 0.013513,
		0.127170, 0.016923, 0.344765, 0.003656, 0.445038, 0.165753, 0.008541, 0.002533, 0.031779,
		0.292429, 0.064289, 0.210724, 0.004200, 1.217010, 1.088704, 0.014768, 0.005848, 0.064558, 7.278994,
		0.071458, 0.855973, 1.172204, 0.014189, 0.033969, 1.889645, 0.125869, 0.031390, 0.065585, 0.029917, 0.042762,
		1.218562, 0.079621, 0.763553, 0.009876, 1.988516, 3.344809, 0.056702, 0.021612, 0.079927, 7.918203, 14.799537, 0.259400,
		0.075144, 0.011169, 0.082464, 0.002656, 0.681161, 0.111063, 0.004186, 0.004854, 0.095591, 0.450964, 1.506485, 0.009457, 1.375871,
		7.169085, 0.161937, 0.726566, 0.040244, 0.825960, 2.067758, 0.110993, 0.129497, 0.196886, 0.169797, 0.637893, 0.090576, 0.457399, 0.143327,
		30.139501, 0.276530, 11.149790, 0.267322, 18.762977, 3.547017, 0.201148, 0.976631, 0.408834, 0.104288, 0.123793, 0.292108, 0.598048, 0.328689, 3.478333,
		13.461692, 0.161053, 4.782635, 0.053740, 11.949233, 2.466507, 0.139705, 0.053397, 0.126088, 1.578530, 0.641351, 0.297913, 4.418398, 0.125011, 2.984862, 13.974326,
		0.021372, 0.081472, 0.058046, 0.006597, 0.286794, 0.188236, 0.009201, 0.019475, 0.037226, 0.015909, 0.154810, 0.017172, 0.239749, 0.562720, 0.061299, 0.154326, 0.060703,
		0.045779, 0.036742, 0.498072, 0.027639, 0.534219, 0.203493, 0.012095, 0.004964, 0.452302, 0.094365, 0.140750, 0.021976, 0.168432, 1.414883, 0.077470, 0.224675, 0.123480, 0.447011,
		4.270235, 0.030342, 0.258487, 0.012745, 4.336817, 0.281953, 0.043812, 0.015539, 0.016212, 16.179952, 3.416059, 0.032578, 2.950318, 0.227807, 1.050562, 0.112000, 5.294490, 0.033381, 0.045528
	      },	    
	      {
		0.733336,
		0.558955, 0.597671,
		0.503360, 0.058964, 5.581680,
		4.149599, 2.863355, 1.279881, 0.225860,
		1.415369, 2.872594, 1.335650, 0.434096, 1.043232,
		1.367574, 0.258365, 0.397108, 2.292917, 0.209978, 4.534772,
		1.263002, 0.366868, 1.840061, 1.024707, 0.823594, 0.377181, 0.496780,
		0.994098, 2.578946, 5.739035, 0.821921, 3.039380, 4.877840, 0.532488, 0.398817,
		0.517204, 0.358350, 0.284730, 0.027824, 1.463390, 0.370939, 0.232460, 0.008940, 0.349195,
		0.775054, 0.672023, 0.109781, 0.021443, 1.983693, 1.298542, 0.169219, 0.043707, 0.838324, 5.102837,
		0.763094, 5.349861, 1.612642, 0.088850, 0.397640, 3.509873, 0.755219, 0.436013, 0.888693, 0.561690, 0.401070,
		1.890137, 0.691594, 0.466979, 0.060820, 2.831098, 2.646440, 0.379926, 0.087640, 0.488389, 7.010411, 8.929538, 1.357738,
		0.540460, 0.063347, 0.141582, 0.018288, 4.102068, 0.087872, 0.020447, 0.064863, 1.385133, 3.054968, 5.525874, 0.043394, 3.135353,
		0.200122, 0.032875, 0.019509, 0.042687, 0.059723, 0.072299, 0.023282, 0.036426, 0.050226, 0.039318, 0.067505, 0.023126, 0.012695, 0.015631,
		4.972745, 0.821562, 4.670980, 1.199607, 5.901348, 1.139018, 0.503875, 1.673207, 0.962470, 0.204155, 0.273372, 0.567639, 0.570771, 0.458799, 0.233109,
		1.825593, 0.580847, 1.967383, 0.420710, 2.034980, 0.864479, 0.577513, 0.124068, 0.502294, 2.653232, 0.437116, 1.048288, 2.319555, 0.151684, 0.077004, 8.113282,
		0.450842, 0.661866, 0.088064, 0.037642, 2.600668, 0.390688, 0.109318, 0.218118, 1.065585, 0.564368, 1.927515, 0.120994, 1.856122, 4.154750, 0.011074, 0.377578, 0.222293,
		0.526135, 0.265730, 0.581928, 0.141233, 5.413080, 0.322761, 0.153776, 0.039217, 8.351808, 0.854294, 0.940458, 0.180650, 0.975427, 11.429924, 0.026268, 0.429221, 0.273138, 4.731579,
		3.839269, 0.395134, 0.145401, 0.090101, 4.193725, 0.625409, 0.696533, 0.104335, 0.377304, 15.559906, 2.508169, 0.449074, 3.404087, 1.457957, 0.052132, 0.260296, 2.903836, 0.564762, 0.681215
	      },
	      {
		0.658412,
		0.566269, 0.540749,
		0.854111, 0.058015, 3.060574,
		0.884454, 5.851132, 1.279257, 0.160296,
		1.309554, 2.294145, 1.438430, 0.482619, 0.992259,
		1.272639, 0.182966, 0.431464, 2.992763, 0.086318, 2.130054,
		1.874713, 0.684164, 2.075952, 1.296206, 2.149634, 0.571406, 0.507160,
		0.552007, 3.192521, 4.840271, 0.841829, 5.103188, 4.137385, 0.351381, 0.679853,
		0.227683, 0.528161, 0.644656, 0.031467, 3.775817, 0.437589, 0.189152, 0.025780, 0.665865,
		0.581512, 1.128882, 0.266076, 0.048542, 3.954021, 2.071689, 0.217780, 0.082005, 1.266791, 8.904999,
		0.695190, 3.010922, 2.084975, 0.132774, 0.190734, 2.498630, 0.767361, 0.326441, 0.680174, 0.652629, 0.440178,
		0.967985, 1.012866, 0.720060, 0.133055, 1.776095, 1.763546, 0.278392, 0.343977, 0.717301, 10.091413, 14.013035, 1.082703,
		0.344015, 0.227296, 0.291854, 0.056045, 4.495841, 0.116381, 0.092075, 0.195877, 4.001286, 2.671718, 5.069337, 0.091278, 4.643214,
		0.978992, 0.156635, 0.028961, 0.209188, 0.264277, 0.296578, 0.177263, 0.217424, 0.362942, 0.086367, 0.539010, 0.172734, 0.121821, 0.161015,
		3.427163, 0.878405, 4.071574, 0.925172, 7.063879, 1.033710, 0.451893, 3.057583, 1.189259, 0.359932, 0.742569, 0.693405, 0.584083, 1.531223, 1.287474,
		2.333253, 0.802754, 2.258357, 0.360522, 2.221150, 1.283423, 0.653836, 0.377558, 0.964545, 4.797423, 0.780580, 1.422571, 4.216178, 0.599244, 0.444362, 5.231362,
		0.154701, 0.830884, 0.073037, 0.094591, 3.017954, 0.312579, 0.074620, 0.401252, 1.350568, 0.336801, 1.331875, 0.068958, 1.677263, 5.832025, 0.076328, 0.548763, 0.208791,
		0.221089, 0.431617, 1.238426, 0.313945, 8.558815, 0.305772, 0.181992, 0.072258, 12.869737, 1.021885, 1.531589, 0.163829, 1.575754, 33.873091, 0.079916, 0.831890, 0.307846, 5.910440,
		2.088785, 0.456530, 0.199728, 0.118104, 4.310199, 0.681277, 0.752277, 0.241015, 0.531100, 23.029406, 4.414850, 0.481711, 5.046403, 1.914768, 0.466823, 0.382271, 3.717971, 0.282540, 0.964421
	      }
	    };

	  double
	    freqs[4][20] = 
	    {{0.147383 , 0.017579 , 0.058208 , 0.017707 , 0.026331 , 0.041582 , 0.017494 , 0.027859 , 0.011849 , 0.076971 , 
	      0.147823 , 0.019535 , 0.037132 , 0.029940 , 0.008059 , 0.088179 , 0.089653 , 0.006477 , 0.032308 , 0.097931},
	     {0.063139 , 0.066357 , 0.011586 , 0.066571 , 0.010800 , 0.009276 , 0.053984 , 0.146986 , 0.034214 , 0.088822 , 
	      0.098196 , 0.032390 , 0.021263 , 0.072697 , 0.016761 , 0.020711 , 0.020797 , 0.025463 , 0.045615 , 0.094372},
	     {0.062457 , 0.066826 , 0.049332 , 0.065270 , 0.006513 , 0.041231 , 0.058965 , 0.080852 , 0.028024 , 0.037024 , 
	      0.075925 , 0.064131 , 0.019620 , 0.028710 , 0.104579 , 0.056388 , 0.062027 , 0.008241 , 0.033124 , 0.050760},
	     {0.106471 , 0.074171 , 0.044513 , 0.096390 , 0.002148 , 0.066733 , 0.158908 , 0.037625 , 0.020691 , 0.014608 , 
	      0.028797 , 0.105352 , 0.007864 , 0.007477 , 0.083595 , 0.055726 , 0.047711 , 0.003975 , 0.010088 , 0.027159}};

	  makeAASubstMat(daa, f,  rates[lg4_index],  freqs[lg4_index]);

	  /*int 
	    i, 
	    j, 
	    r = 0;

	  for(i = 1; i < 20; i++)
	    for(j = 0; j < i; j++)
	      {
		daa[i * 20 + j] = rates[lg4_index][r];
		r++;
	      }
	  
	  assert(r == 190);

	  for(i = 0; i < 20; i++)
	  f[i] = freqs[lg4_index][i];	  */
	  
	}
	break;
      case STMTREV:
	{
	  double rates[190] =
	    {
	      0.1159435373,  
	      0.2458816714, 0.1355713516, 
	      0.9578712472, 0.0775041665, 8.4408676914, 
	      0.2327281954, 9.1379470330, 0.1137687264, 0.0582110367, 
	      0.3309250853, 5.2854173238, 0.1727184754, 0.8191776581, 0.0009722083, 
	      0.6946680829, 0.0966719296, 0.2990806606, 7.3729791633, 0.0005604799, 3.5773486727, 
	      2.8076062202, 3.0815651393, 0.5575702616, 2.2627839242, 1.1721237455, 0.0482085663, 3.3184632572,  
	      0.2275494971, 2.8251848421, 9.5228608030, 2.3191131858, 0.0483235836, 4.4138715270, 0.0343694246, 0.0948383460, 
	      0.0627691644, 0.5712158076, 0.2238609194, 0.0205779319, 0.1527276944, 0.0206129952, 0.0328079744, 0.1239000315, 0.0802374651, 
	      0.0305818840, 0.1930408758, 0.0540967250, 0.0018843293, 0.2406073246, 0.3299454620, 0.0373753435, 0.0005918940, 0.1192904610, 1.3184058362, 
	      0.2231434272, 6.0541970908, 4.3977466558, 0.1347413792, 0.0001480536, 5.2864094506, 6.8883522181, 0.5345755286, 0.3991624551, 0.2107928508, 0.1055933141, 
	      0.1874527991, 0.2427875732, 0.0433577842, 0.0000022173, 0.0927357503, 0.0109238300, 0.0663619185, 0.0128777966, 0.0722334577, 4.3016010974, 1.1493262595, 0.4773694701, 
	      0.0458112245, 0.0310030750, 0.0233493970, 0.0000080023, 0.8419347601, 0.0027817812, 0.0361207581, 0.0490593583, 0.0197089530, 0.3634155844, 2.1032860162, 0.0861057517, 0.1735660361, 
	      1.5133910481, 0.7858555362, 0.3000131148, 0.3337627573, 0.0036260499, 1.5386413234, 0.5196922389, 0.0221252552, 1.0171151697, 0.0534088166, 6.0377879080, 0.4350064365, 0.1634497017, 
	      0.3545179411, 
	      2.3008246523, 0.7625702322, 1.9431704326, 0.6961369276, 2.3726544756, 0.1837198343, 0.9087013201, 2.5477016916, 0.3081949928, 0.1713464632, 2.7297706102, 0.3416923226, 0.0730798705, 
	      4.0107845583, 8.4630191575, 
	      4.3546170435, 1.0655012755, 1.6534489471, 0.0985354973, 0.1940108923, 0.3415280861, 0.2794040892, 0.1657005971, 0.2704552047, 2.3418182855, 0.0426297282, 1.2152488582, 4.6553742047, 
	      0.0068797851, 1.1613183519, 2.2213527952, 
	      0.0565037747, 6.7852754661, 0.0000010442, 0.0000002842, 0.9529353202, 0.0009844045, 0.0002705734, 0.5068170211, 0.0000932799, 0.0050518699, 0.3163744815, 0.0000023280, 0.1010587493, 
	      0.2890102379, 0.0041564377, 0.0495269526, 0.0002026765, 
	      0.0358664532, 0.0714121777, 0.3036789915, 1.3220740967, 1.7972997876, 0.0066458178, 0.3052655031, 0.0174305437, 21.9842817264, 0.1070890246, 0.0770894218, 0.1929529483, 0.0561599188, 
	      1.6748429971, 0.0021338646, 1.8890678523, 0.2834320440, 0.3134203648, 
	      3.2116908598, 0.0108028571, 0.0860833645, 0.0426724431, 0.3652373073, 0.0287789552, 0.1484349765, 0.5158740953, 0.0059791370, 3.3648305163, 0.8763855707, 0.0776875418, 0.9145670668, 
	      0.3963331926, 0.1080226203, 0.0640951379, 0.2278998021, 0.0388755869, 0.1836950254}; 
	    
	  double
	    freqs[20] = {0.0461811000, 0.0534080000, 0.0361971000, 0.0233326000, 0.0234170000, 0.0390397000, 0.0341284001, 0.0389164000, 0.0164640000, 0.0891534000, 
			 0.1617310001, 0.0551341000, 0.0233262000, 0.0911252000, 0.0344713001, 0.0771077000, 0.0418603001, 0.0200784000, 0.0305429000, 0.0643851996};  

	  makeAASubstMat(daa, f, rates, freqs);
	  
	}
	break;
      case DUMMY:
	  {
	    double 
	      rates[190] = {10.7,
			    0.2,   3.7,
			    4.3,   0.2, 1024.8,
			    26.1,  28.2,  46.5,   0.2,
			    0.2, 316.9,  32.4,  12.6,  13.7,
			    23.0,   0.2,  26.6, 814.8, 0.2, 159.5,
			    76.3,   9.6,  39.3,  33.1,  23.5,   2.6,  56.8,
			    0.2, 212.7, 294.3,  59.4, 131.9, 604.6,   5.1,   1.5,
			    41.9,   0.7,  15.2,   0.2,   0.2,   0.2,   0.2,   1.4,   0.2,
			    4.8,   4.9,   0.2,   0.2,  19.9,  25.3,   0.2,   0.2,  11.0, 133.7,
			    1.7,  17.2, 331.6,   5.5,   0.2, 164.7, 237.8,   5.7,  12.7,   0.2,   0.2,
			    124.0,   0.2,   0.2,   0.2,   0.2,   6.9,   1.5,   0.5,   1.2, 322.3, 416.2,  36.8,
			    2.5,   0.2,   0.4,   0.2,  77.4,   0.2,   0.2,   0.7,  15.4,  44.1, 277.9,   0.2,   8.3,
			    41.3,   8.3,   0.2,   2.5,   0.2,  75.8,   0.2,   0.4,  75.7,   2.6,  52.1,   3.9,   1.3,  11.1,
			    307.4,   7.1, 491.3,  15.0, 437.3,  25.3,   7.4, 105.9,  56.4,   2.4,  81.0,  24.2,   6.6, 122.8, 327.1,
			    782.7,   0.4,  94.0,  19.0,  24.9,  12.5,   1.6,   0.2,   9.6, 360.6,  16.6,  48.8, 651.4,   3.9,  64.7, 455.7,
			    0.2,  18.3,   0.2,   0.2, 133.5,   4.3,   0.2,   6.5,   1.2,   0.2,  10.7,   0.2,   0.6,   2.0,   0.5,   4.2,   0.2,
			    0.2,   1.0,  85.7,  29.4, 994.7,   6.6,   0.2,   0.2, 2096.0,   6.2,  15.5,   1.6,   6.5, 502.9,  11.5,  63.6,   4.3,  10.1,
			    565.4,   0.2,   1.7,  14.9,  20.9,   0.2,  13.5,  13.5,   0.2, 1987.0,  74.0,   1.7, 716.8,   5.6,   2.7,   0.9, 246.3,   3.0,   0.2};

	    double
	      freqs[20] = { 0.066446,    0.017604,    0.043105,    0.017760,    0.005969,    0.024329,    0.023622,    0.052890,    0.026973,    0.088543,    
			    0.162813,   0.025336,    0.062589,    0.061567,    0.053608,    0.074271,    0.087828,    0.027617,    0.034022,    0.043108};
	    
	    makeAASubstMat(daa, f,  rates,  freqs);

	    /* int 
	      i, j, r = 0;

	    for(i = 1; i < 20; i++)
	      for(j = 0; j < i; j++)
		{
		  daa[i * 20 + j] = rates[r];
		  r++;
		}
	       
	    assert(r == 190);

	    for(i = 0; i < 20; i++)
	    f[i] = freqs[i];	  */
	  }
	  break;
	case DUMMY2:
	  {
	    double rates[190] = { 6.5,
				  4.5,  10.6,
				  84.3,   9.5, 643.2,
				  19.5, 353.7,  10.9,  10.7,
				  6.1, 486.3,  18.0,  11.6,   0.1,
				  74.5,  21.5,  13.0, 437.4,   0.1, 342.6,
				  118.1, 183.9,  17.4, 150.3,  86.8,   7.1, 161.9,
				  2.8, 346.6, 345.3, 202.4, 111.8, 450.1,   6.2,   2.2,
				  1.5,  50.6,  25.6,   5.6,   3.4,   3.6,   4.3,   2.5,   8.4,
				  3.9,  36.9,   2.4,   5.9,  20.3,  26.1,   5.1,   3.4,  17.3, 205.0,
				  4.2, 712.1, 639.2,  10.1,   0.1, 500.5, 426.6,  29.3,   9.2,  37.9,  10.8,
				  13.4,  53.5,   9.9,   3.8,  10.5,   9.5,   9.6,   3.8,   3.6, 534.9, 142.8,  83.6,
				  4.3,   5.0,   8.7,   7.5, 238.0,   2.4,   7.7,   3.1,  11.0,  61.0, 542.3,   9.4,   3.8,
				  91.2,  69.0,   3.5,  13.4,   6.5, 145.6,   8.1,   2.6, 133.9,   2.1, 155.8,  21.2,  10.5,  12.6,
				  251.1,  82.9, 271.4,  34.8, 471.9,  10.7,  16.4, 136.7,  19.2,  36.2, 160.3,  23.9,   6.2, 249.4, 348.6,
				  467.5,  82.5, 215.5,   8.0,   7.4,   5.4,  11.6,   6.3,   3.8, 266.2,  10.7, 140.2, 295.2,   3.6, 181.2, 144.8,
				  3.4, 171.8,   6.1,   3.5, 518.6,  17.0,   9.1,  49.0,   5.7,   3.3,  98.8,   2.3,  11.1,  34.1,   1.1,  56.3,   1.5,
				  2.2,   4.3,  69.9, 202.9, 579.1,   9.4,   9.1,   2.1, 889.2,  10.8,   9.6,  20.1,   3.4, 255.9,   5.6, 264.3,   3.3,  21.7,
				  363.2,   8.4,   1.6,  10.3,  37.8,   5.1,  21.6,  76.0,   1.1, 595.0, 155.8,   9.2, 191.9, 102.2,   7.7,  10.1,  36.8,   5.0,   7.2};

	    double freqs[20] = {0.061007,    0.060799,    0.043028,    0.038515,    0.011297,    0.035406,    0.050764,    0.073749,    0.024609,    0.085629,
				0.106930,    0.046704,    0.023382,    0.056136,    0.043289,    0.073994,    0.052078,    0.018023,    0.036043,    0.058620};

	    makeAASubstMat(daa, f,  rates,  freqs);

	    /*
	    int 
	      i, j, r = 0;

	    for(i = 1; i < 20; i++)
	      for(j = 0; j < i; j++)
		{
		  daa[i * 20 + j] = rates[r];
		  r++;
		}
	       
	    assert(r == 190);

	    for(i = 0; i < 20; i++)
	    f[i] = freqs[i];	*/

	  }
	 
	  break;
	default: 
	  assert(0);
	}
    }


  /*
    
  TODO review frequency sums for fixed as well as empirical base frequencies !
  
  NUMERICAL BUG fix, rounded AA freqs in some models, such that 
  they actually really sum to 1.0 +/- epsilon 
  
  {
    double acc = 0.0;
  
    for(i = 0; i < 20; i++)
      acc += f[i];
    
    printf("%1.80f\n", acc);
    assert(acc == 1.0);  
  }
  */
 


  for (i=0; i<20; i++)  
    for (j=0; j<i; j++)               
      daa[j*20+i] = daa[i*20+j];

  
  /*
    for (i=0; i<20; i++)  
    {
    for (j=0; j<20; j++)
    {
    if(i == j)
    printf("0.0 ");
    else
    printf("%f ", daa[i * 20 + j]);
    }
    printf("\n");
    }
    
    for (i=0; i<20; i++) 
    printf("%f ", f[i]);
    printf("\n");
  */
  

  max = 0;
  
  for(i = 0; i < 19; i++)
    for(j = i + 1; j < 20; j++)
      {
	q[i][j] = temp = daa[i * 20 + j];
	if(temp > max) 
	  max = temp;
      }
 
  scaler = AA_SCALE / max;
   
  /* SCALING HAS BEEN RE-INTRODUCED TO RESOLVE NUMERICAL  PROBLEMS */   

  r = 0;
  for(i = 0; i < 19; i++)
    {      
      for(j = i + 1; j < 20; j++)
	{  
	
	  q[i][j] *= scaler;
	  
	  
	  assert(q[i][j] <= AA_SCALE_PLUS_EPSILON);
	  
	  initialRates[r++] = q[i][j];
	}
    }             
}



static void initGeneric(const int n, const unsigned int *valueVector, int valueVectorLength,		       
			double *ext_EIGN,
			double *EV,
			double *EI,
			double *frequencies,
			double *ext_initialRates,
			double *tipVector)
{
  double
    fracchange = 0.0,
    **r, 
    **a, 
    **EIGV,
    *initialRates = ext_initialRates, 
    *f, 
    *e, 
    *d, 
    *invfreq, 
    *EIGN,
    *eptr; 
  
  int 
    i, 
    j, 
    k, 
    m, 
    l;  

  r    = (double **)rax_malloc(n * sizeof(double *));
  EIGV = (double **)rax_malloc(n * sizeof(double *));  
  a    = (double **)rax_malloc(n * sizeof(double *));	  
  
  for(i = 0; i < n; i++)
    {
      a[i]    = (double*)rax_malloc(n * sizeof(double));
      EIGV[i] = (double*)rax_malloc(n * sizeof(double));
      r[i]    = (double*)rax_malloc(n * sizeof(double));
    }

  f       = (double*)rax_malloc(n * sizeof(double));
  e       = (double*)rax_malloc(n * sizeof(double));
  d       = (double*)rax_malloc(n * sizeof(double));
  invfreq = (double*)rax_malloc(n * sizeof(double));
  EIGN    = (double*)rax_malloc(n * sizeof(double));

  
  for(l = 0; l < n; l++)		 
    f[l] = frequencies[l];	
  /*assert(initialRates[numRates - 1] == 1.0);	*/
  
  i = 0;
  
  for(j = 0; j < n; j++)	 
    for(k = 0; k < n; k++)
      r[j][k] = 0.0;
  
  for(j = 0; j < n - 1; j++)
    for (k = j+1; k < n; k++)      	  
      r[j][k] = initialRates[i++];         
  
  for (j = 0; j < n; j++) 
    {
      r[j][j] = 0.0;
      for (k = 0; k < j; k++)
	r[j][k] = r[k][j];
    }                         
  
  for (j = 0; j< n; j++)
    for (k = 0; k< n; k++)
      fracchange += f[j] * r[j][k] * f[k];             
  
  m = 0;
  
  for(i=0; i< n; i++) 
    a[i][i] = 0;
  
  /*assert(r[n - 2][n - 1] == 1.0);*/
  
  for(i=0; i < n; i++) 
    {
      for(j=i+1;  j < n; j++) 
	{
	  double factor =  initialRates[m++];
	  a[i][j] = a[j][i] = factor * sqrt( f[i] * f[j]);
	  a[i][i] -= factor * f[j];
	  a[j][j] -= factor * f[i];
	}
    }             	        

  makeEigen(a, n, d, e);	   	  	    

  for(i=0; i<n; i++)     
    for(j=0; j<n; j++)       
      a[i][j] *= sqrt(f[j]);
  
  for (i=0; i<n; i++)
    {	  
      if (d[i] > -1e-8) 
	{	      
	  if (i != 0) 
	    {		    
	      double tmp = d[i], sum=0;
	      d[i] = d[0];
	      d[0] = tmp;
	      for (j=0; j < n; j++) 
		{
		  tmp = a[i][j];
		  a[i][j] = a[0][j];
		  sum += (a[0][j] = tmp);
		}
	      for (j=0; j < n; j++) 
		a[0][j] /= sum;
	    }
	  break;
	}
    }
  
  for (i=0; i< n; i++) 
    {
      EIGN[i] = -d[i];
      
      for (j=0; j<n; j++)
	EIGV[i][j] = a[j][i];
      invfreq[i] = 1 / EIGV[i][0]; 
    }                                    
  
  for(l = 1; l < n; l++)
    {
      ext_EIGN[(l - 1)] = EIGN[l] * (1.0 / fracchange); 
      assert( ext_EIGN[(l - 1)] > 0.0);
    }
  
  eptr = EV;
  
  for(i = 0; i < n; i++)		  
    for(j = 0; j < n; j++)
      {
	*eptr++ = EIGV[i][j];	             	    	     

      }
  
  for(i = 0; i < n; i++)
    for(j = 1; j < n; j++)
      EI[i * (n - 1) + (j - 1)] = EV[i * n + j] * invfreq[i];  
  
  for(i=0; i < valueVectorLength; i++)
    {
      unsigned int value = valueVector[i];
      
      for(j = 0; j < n; j++)
	tipVector[i * n + j]     = 0;	            

      if(value > 0)
	{		      
	  for (j = 0; j < n; j++) 
	    {	    
	      if ((value >> j) & 1) 
		{
		  int l;
		  for(l = 0; l < n; l++)
		    tipVector[i * n + l] += EIGV[j][l];		      		      		     		      
		}	     		  
	    }	    
	}     
    }

  for(i = 0; i < valueVectorLength; i++)
    {
       for(j = 0; j < n; j++)
	 if(tipVector[i * n + j] > MAX_TIP_EV)
	   tipVector[i * n + j] = MAX_TIP_EV;
    }


  

  for(i = 0; i < n; i++)
    {
      rax_free(EIGV[i]);
      rax_free(a[i]);
      rax_free(r[i]);
    }

  rax_free(r);
  rax_free(a);
  rax_free(EIGV);

  rax_free(f);
  rax_free(e);
  rax_free(d);
  rax_free(invfreq);
  rax_free(EIGN);
}




void initReversibleGTR(tree *tr, int model)
{ 
  double    
    *ext_EIGN         = tr->partitionData[model].EIGN,
    *EV               = tr->partitionData[model].EV,
    *EI               = tr->partitionData[model].EI,
    *frequencies      = tr->partitionData[model].frequencies,
    *ext_initialRates = tr->partitionData[model].substRates,
    *tipVector        = tr->partitionData[model].tipVector;
  
 int 
   states = tr->partitionData[model].states;

 switch(tr->partitionData[model].dataType)
   { 
   case GENERIC_32:
   case GENERIC_64:
   case SECONDARY_DATA_6:
   case SECONDARY_DATA_7: 
   case SECONDARY_DATA:
   case DNA_DATA:
   case BINARY_DATA:
     initGeneric(states, 
		 getBitVector(tr->partitionData[model].dataType), 
		 getUndetermined(tr->partitionData[model].dataType) + 1, 	       
		 ext_EIGN, 
		 EV, 
		 EI, 
		 frequencies, 
		 ext_initialRates,
		 tipVector);
#ifdef _HET
     assert(tr->partitionData[model].dataType == DNA_DATA);

      initGeneric(states, 
		  getBitVector(tr->partitionData[model].dataType), 
		  getUndetermined(tr->partitionData[model].dataType) + 1, 	       
		  tr->partitionData[model].EIGN_TIP, 
		  tr->partitionData[model].EV_TIP, 
		  tr->partitionData[model].EI_TIP, 
		  frequencies, 
		  tr->partitionData[model].substRates_TIP,
		  tr->partitionData[model].tipVector_TIP);
#endif
     break;   
   case AA_DATA: 

     assert(!(tr->partitionData[model].usePredefinedProtFreqs && tr->partitionData[model].optimizeBaseFrequencies));
     
     if(!((tr->partitionData[model].protModels == GTR) || (tr->partitionData[model].protModels == GTR_UNLINKED)))
       {
	 double 
	   f[20];
	 
	 int 
	   l;
       
	 if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
	   {
	     int 
	       i;
	    
	     //don't allow base freq opt for LG4 and LG4X models!

	     //assert(!tr->partitionData[model].optimizeBaseFrequencies);

	     for(i = 0; i < 4; i++)
	       {	
		 initProtMat(f, tr->partitionData[model].protModels, &(tr->partitionData[model].substRates_LG4[i][0]), model, tr, i);
		 
		 if(tr->partitionData[model].usePredefinedProtFreqs == TRUE)
		   {
		     for(l = 0; l < 20; l++)		
		       tr->partitionData[model].frequencies_LG4[i][l] = f[l];
		   }
		 else
		   {		    
		     memcpy(tr->partitionData[model].frequencies_LG4[i], frequencies, 20 * sizeof(double));
		   }
	       }
	   }
	 else
	   {
	     if(tr->partitionData[model].protModels == AUTO)
	       {
		 //printf("init prot mat %s partition %d\n", protModels[tr->partitionData[model].autoProtModels], model);
		 initProtMat(f, tr->partitionData[model].autoProtModels, ext_initialRates, model, tr, 0);

		
		     
		 //printf("t1\n");
		 if(tr->partitionData[model].usePredefinedProtFreqs == FALSE && !tr->partitionData[model].optimizeBaseFrequencies)		   
		   genericBaseFrequencies(tr, tr->partitionData[model].states, tr->rdta, tr->cdta, tr->partitionData[model].lower, tr->partitionData[model].upper, model, 
					  getSmoothFreqs(tr->partitionData[model].dataType),
					  getBitVector(tr->partitionData[model].dataType));
		 //printf("t2\n");
		   
	       }
	     else
	       initProtMat(f, tr->partitionData[model].protModels, ext_initialRates, model, tr, 0);
	     
	     if(tr->partitionData[model].protModels == PROT_FILE)
	       assert(tr->partitionData[model].usePredefinedProtFreqs == TRUE);  
	     
	     if(tr->partitionData[model].usePredefinedProtFreqs == TRUE)
	       {
		 for(l = 0; l < 20; l++)		
		   frequencies[l] = f[l];
	       }
	     else
	       {
		 //nothing to do here, the base freqs have already been initialized or optimized
		 //don't overwrite them!
		 
		 /*
		   if(tr->partitionData[model].optimizeBaseFrequencies)
		   {
		   for(l = 0; l < 20; l++)		
		   frequencies[l] = 1.0 / 20.0;
		   }
		 */
	       }
	     
	   }	   
       }
     else 
       {
	 assert(tr->partitionData[model].usePredefinedProtFreqs == FALSE);

	 //nothing to do here, the base freqs have already been initialized or optimized
	 //don't overwrite them!

	 /*
	   if(tr->partitionData[model].optimizeBaseFrequencies)
	   {
	   int l;
	   
	   for(l = 0; l < 20; l++)		
	   frequencies[l] = 1.0 / 20.0;
	   }
	 */
       }

     if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
       {
	 int 
	   i;     
	
	 for(i = 0; i < 4; i++)	   
	   initGeneric(states, bitVectorAA, 23,
		       tr->partitionData[model].rawEIGN_LG4[i],  tr->partitionData[model].EV_LG4[i],  tr->partitionData[model].EI_LG4[i], tr->partitionData[model].frequencies_LG4[i], 
		       tr->partitionData[model].substRates_LG4[i],
		       tr->partitionData[model].tipVector_LG4[i]);   	    	    	  
	 
	 scaleLG4X_EIGN(tr, model);	 
       }
     else
       initGeneric(states, bitVectorAA, 23, 
		   ext_EIGN, EV, EI, frequencies, ext_initialRates,
		   tipVector);                   
     break;  
   default:
     assert(0);
   }  
}


double LnGamma (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
  double x, f, z, result;

  x = alpha;
  f = 0.0;
  
  if ( x < 7.0) 
     {
       f = 1.0;  
       z = alpha - 1.0;
      
       while ((z = z + 1.0) < 7.0)  
	 {	  
	   f *= z;
	 }
       x = z;   
     
       assert(f != 0.0);
	
       f=-log(f);
     }
   
   z = 1/(x*x);
   
   result = f + (x-0.5)*log(x) - x + .918938533204673 
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;  

   return result;
}



double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, pn[6];


   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

   
   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:  
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];   
 l32:  
   a++;  
   b+=2;  
   term++;   
   an=a*term;
   for (i=0; i<2; i++) 
     pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   
   dif=fabs(gin-rn);  
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:   
   gin=rn;
 l35: 
   for (i=0; i<4; i++) 
     pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow)            
     goto l32;        
   
   for (i=0; i<4; i++) 
     pn[i]/=overflow;

   
   goto l32;
 l42:  
   gin=1-factor*gin;

 l50: 
   return (gin);
}




double PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}


double PointChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;
  
   if (p<.000002 || p>.999998 || v<=0) return (-1);
  
   g = LnGamma(v/2);
   
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3:    
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;   
   if ((t=IncompleteGamma (p1, xx, g))< 0.0) 
     {
       printf ("IncompleteGamma \n");      
       return (-1);
     }
  
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}






void makeGammaCats(int rateHetModel, double alpha, double *gammaRates, int K, boolean useMedian, double propInvariant)
{
  int 
    i;

  double 
    factor = alpha / alpha * K, 
    lnga1, 
    alfa = alpha, 
    beta = alpha,
    *gammaProbs = (double *)rax_malloc(K * sizeof(double));

  /* Note that ALPHA_MIN setting is somewhat critical due to   */
  /* numerical instability caused by very small rate[0] values */
  /* induced by low alpha values around 0.01 */

  assert(alfa >= ALPHA_MIN); 

  if(useMedian)
    {
      double  
	middle = 1.0 / (2.0*K),
	t = 0.0; 
      
      for(i = 0; i < K; i++)     
	gammaRates[i] = PointGamma((double)(i * 2 + 1) * middle, alfa, beta);
      
      for (i = 0; i < K; i++) 
	t += gammaRates[i];
       for( i = 0; i < K; i++)     
	 gammaRates[i] *= factor / t;
    }
  else
    {
      lnga1 = LnGamma(alfa + 1);

      for (i = 0; i < K - 1; i++)
	gammaProbs[i] = PointGamma((i + 1.0) / K, alfa, beta);

      for (i = 0; i < K - 1; i++)
	gammaProbs[i] = IncompleteGamma(gammaProbs[i] * beta, alfa + 1, lnga1);   

      gammaRates[0] = gammaProbs[0] * factor;
      
      gammaRates[K - 1] = (1 - gammaProbs[K - 2]) * factor;

      for (i= 1; i < K - 1; i++)  
	gammaRates[i] = (gammaProbs[i] - gammaProbs[i - 1]) * factor;      
    }
  /* assert(gammaRates[0] >= 0.00000000000000000000000000000044136090435925743185910935350715027016962154188875); */

  if(rateHetModel == GAMMA_I)
    {
      double 
	scaler = 1.0 / (1.0 - propInvariant);

      for (i = 0; i < K; i++)
        gammaRates[i] *= scaler;
    }
     

  rax_free(gammaProbs);

  return;  
}

static void genericInvariant(tree *tr, int lower, int upper, const unsigned int *bitVector, 
			     unsigned int undetermined, int states, int *numberOfInvariableColumns, int *weightOfInvariableColumns)
{
  int
    count = 0,
    sum  = 0,
    i,
    j;

  for(i = lower; i < upper; i++)	      
    {		 
      unsigned int 
	encoding = 0, 
	code;
      int 
	secSum = 0,
	position = -1;	  	

      for(j = 1; j <= tr->mxtips; j++)
	{
	  code = bitVector[tr->yVector[j][i]];
	  if(code != undetermined)
	    {
	      if(!(code & encoding))		
		encoding = encoding | code;
	      else
		encoding = encoding | (encoding & code);	     
	    }
	}

      for(j = 0; j < states; j++)
	{
	  if(encoding >> j & 1)
	    {
	      secSum++;
	      position = j;		  
	    }
	}

      if(secSum == 1)
	{
	  assert(position >= 0);
	  tr->invariant[i] = position;
	  count = count + 1;
	  sum = sum + tr->cdta->aliaswgt[i];
	}
      else
	tr->invariant[i] = states;
    }

  *numberOfInvariableColumns += count;
  *weightOfInvariableColumns += sum;
}

static void setRates(double *r, int rates, boolean JC69)
{
  int 
    i;


  if(JC69)
    for(i = 0; i < rates - 1; i++)    
      r[i] = 1.0;
  else
    for(i = 0; i < rates - 1; i++)    
      r[i] = 0.5;

  r[rates - 1] = 1.0;
}

void initRateMatrix(tree *tr)
{
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {	
      boolean 
	JC69 = FALSE;

      int 	
	i,
	states = tr->partitionData[model].states,
	rates  = (states * states - states) / 2;
      
      if(tr->partitionData[model].dataType == DNA_DATA && (tr->useJC69 || tr->useK80 || tr->useHKY85))
	JC69 = TRUE;
      
      switch(tr->partitionData[model].dataType)
	{
	case BINARY_DATA:
	case DNA_DATA:
	case SECONDARY_DATA:
	case SECONDARY_DATA_6:
	case SECONDARY_DATA_7:
	  setRates(tr->partitionData[model].substRates, rates, JC69);
#ifdef _HET
	  assert(tr->partitionData[model].dataType = DNA_DATA);
	  setRates(tr->partitionData[model].substRates_TIP, rates, JC69);
#endif
	  break;	  
	case GENERIC_32:
	case GENERIC_64:	  
	  switch(tr->multiStateModel)
	    {
	    case ORDERED_MULTI_STATE:
	      {
		int 
		  j, 
		  k, 
		  i = 0;
		
		for(j = 0; j < states; j++)
		  for(k = j + 1; k < states; k++)
		    tr->partitionData[model].substRates[i++] = (double)(k - j);			
		assert(i == rates);		
	      }
	      break;
	    case MK_MULTI_STATE:
	      for(i = 0; i < rates; i++)
		tr->partitionData[model].substRates[i] = 1.0;
	      
	      break;
	    case GTR_MULTI_STATE:
	      setRates(tr->partitionData[model].substRates, rates, JC69);
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case AA_DATA:
	  if(tr->partitionData[model].protModels == GTR || tr->partitionData[model].protModels == GTR_UNLINKED)	      
	    putWAG(tr->partitionData[model].substRates);
	  break;
	default:
	  assert(0);
	}           
      
      if(tr->partitionData[model].nonGTR)
	{
	  assert(tr->partitionData[model].dataType == SECONDARY_DATA || 
		 tr->partitionData[model].dataType == SECONDARY_DATA_6 ||
		 tr->partitionData[model].dataType == SECONDARY_DATA_7);
	  	  
	  for(i = 0; i < rates; i++)
	    {
	      if(tr->partitionData[model].symmetryVector[i] == -1)
		tr->partitionData[model].substRates[i] = 0.0;
	      else
		{
		  if(tr->partitionData[model].symmetryVector[i] == tr->partitionData[model].symmetryVector[rates - 1])
		    tr->partitionData[model].substRates[i] = 1.0;
		}
	    }
	}
    }  
}

static void setSymmetry(int *s, int *sDest, const int sCount, int *f, int *fDest, const int fCount)
{
  int i;

  for(i = 0; i < sCount; i++)
    sDest[i] = s[i];

  for(i = 0; i < fCount; i++)
    fDest[i] = f[i];
}


static void setupK80Symmetries(tree *tr)
{
  int model;

  assert(tr->useK80 || tr->useHKY85);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(tr->partitionData[model].dataType == DNA_DATA)
	{
	  int 
	    s[6] = {1, 0, 1, 1, 0, 1};

	  memcpy(tr->partitionData[model].symmetryVector, s, sizeof(int) * 6);
	}
    }
}


static void setupSecondaryStructureSymmetries(tree *tr)
{
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(tr->partitionData[model].dataType == SECONDARY_DATA || 
	 tr->partitionData[model].dataType == SECONDARY_DATA_6 || 
	 tr->partitionData[model].dataType == SECONDARY_DATA_7)
	{	
	  switch(tr->secondaryStructureModel)
	    {
	    case SEC_6_A:
	      tr->partitionData[model].nonGTR = FALSE;
	      break;
	    case SEC_6_B:
	      {
		int f[6]  = {0, 1, 2, 3, 4, 5};
		int s[15] = {2, 0, 1, 2, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2, 1};

		setSymmetry(s, tr->partitionData[model].symmetryVector, 15, f, tr->partitionData[model].frequencyGrouping, 6);
		  
		tr->partitionData[model].nonGTR = TRUE;
	      }
	      break;
	    case SEC_6_C:
	      {
		int f[6]  = {0, 2, 2, 1, 0, 1};
		int s[15] = {2, 0, 1, 2, 2, 2, 2, 0, 1, 1, 2, 2, 2, 2, 1};

		setSymmetry(s, tr->partitionData[model].symmetryVector, 15, f, tr->partitionData[model].frequencyGrouping, 6);
		
		tr->partitionData[model].nonGTR = TRUE;
	      }
	      break;
	    case SEC_6_D:
	      {
		int f[6]  = {0, 2, 2, 1, 0, 1};
		int s[15] = {2, -1, 1, 2, 2, 2, 2, -1, 1, 1, 2, 2, 2, 2, 1};

		setSymmetry(s, tr->partitionData[model].symmetryVector, 15, f, tr->partitionData[model].frequencyGrouping, 6);

		tr->partitionData[model].nonGTR = TRUE;
	      }
	      break;
	    case SEC_6_E:
	      {
		int f[6]  = {0, 1, 2, 3, 4, 5};
		int s[15] = {2, -1, 1, 2, 2, 2, 2, -1, 1, 1, 2, 2, 2, 2, 1};

		setSymmetry(s, tr->partitionData[model].symmetryVector, 15, f, tr->partitionData[model].frequencyGrouping, 6);

		tr->partitionData[model].nonGTR = TRUE;
	      }
	      break;
	    case SEC_7_A:
	      tr->partitionData[model].nonGTR = FALSE;
	      break;
	    case SEC_7_B:
	      {
	      	int f[7]  = {0, 2, 2, 1, 0, 1, 3};
		int s[21] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 21, f, tr->partitionData[model].frequencyGrouping, 7);

		tr->partitionData[model].nonGTR = TRUE;

	      }
	      break;
	    case SEC_7_C:
	      {
	      	int f[7]  = {0, 1, 2, 3, 4, 5, 6};
		int s[21] = {-1, -1, 0, -1, -1, 4, -1, -1, -1, 3, 5, 1, -1, -1, 6, -1, -1, 7, 2, 8, 9};
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 21, f, tr->partitionData[model].frequencyGrouping, 7);

		tr->partitionData[model].nonGTR = TRUE;

	      }
	      break;
	    case SEC_7_D:
	      {
	      	int f[7]  = {0, 1, 2, 3, 4, 5, 6};
		int s[21] = {2, 0, 1, 2, 2, 3, 2, 2, 0, 1, 3, 1, 2, 2, 3, 2, 2, 3, 1, 3, 3};
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 21, f, tr->partitionData[model].frequencyGrouping, 7);

		tr->partitionData[model].nonGTR = TRUE;

	      }
	      break;
	    case SEC_7_E:
	      {
	      	int f[7]  = {0, 1, 2, 3, 4, 5, 6};
		int s[21] = {-1, -1, 0, -1, -1, 1, -1, -1, -1, 0, 1, 0, -1, -1, 1, -1, -1, 1, 0, 1, 1};
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 21, f, tr->partitionData[model].frequencyGrouping, 7);

		tr->partitionData[model].nonGTR = TRUE;

	      }
	      break;
	    case SEC_7_F:
	      {
	      	int f[7]  = {0, 2, 2, 1, 0, 1, 3};
		int s[21] = {2, 0, 1, 2, 2, 3, 2, 2, 0, 1, 3, 1, 2, 2, 3, 2, 2, 3, 1, 3, 3};		
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 21, f, tr->partitionData[model].frequencyGrouping, 7);

		tr->partitionData[model].nonGTR = TRUE;

	      }
	      break;
	      
	    case SEC_16:
	      tr->partitionData[model].nonGTR = FALSE;
	      break;
	    case SEC_16_A:
	      {
	      	int f[16]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		int s[120] = {/* AA */  4,  4,  3,  4, -1, -1, -1,  4, -1, -1, -1,  3, -1, -1, -1,
			      /* AC */  4,  3, -1,  4, -1, -1, -1,  3, -1, -1, -1,  4, -1, -1,
			      /* AG */  3, -1, -1,  3, -1, -1, -1,  4, -1, -1, -1,  3, -1,
			      /* AU */ -1, -1,  2,  3, -1,  0, -1,  1,  2, -1,  2,  3,
			      /* CA */  4,  3,  4,  4, -1, -1, -1,  3, -1, -1, -1,
			      /* CC */  3,  4, -1,  3, -1, -1, -1,  4, -1, -1,
			      /* CG */  3, -1,  2,  3,  2,  0, -1,  1, -1,
			      /* CU */ -1, -1, -1,  3, -1, -1, -1,  4,
			      /* GA */  3,  4,  3,  3, -1, -1, -1,
			      /* GC */  3,  1,  2,  3,  2, -1,
			      /* GG */  3, -1, -1,  3, -1,
			      /* GU */  2, -1,  2,  3,
			      /* UA */  3,  1,  3,
			      /* UC */  3,  4,
			      /* UG */  3};
			      
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 120, f, tr->partitionData[model].frequencyGrouping, 16);
			      
		tr->partitionData[model].nonGTR = TRUE;

		}
	      break;
	    case SEC_16_B:
	      {
		int f[16]  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
		int s[120] = {/* AA */  0,  0,  0,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,
			      /* AC */  0,  0, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1, -1,
			      /* AG */  0, -1, -1,  0, -1, -1, -1,  0, -1, -1, -1,  0, -1,
			      /* AU */ -1, -1,  0,  0, -1,  0, -1,  0,  0, -1,  0,  0,
			      /* CA */  0,  0,  0,  0, -1, -1, -1,  0, -1, -1, -1,
			      /* CC */  0,  0, -1,  0, -1, -1, -1,  0, -1, -1,
			      /* CG */  0, -1,  0,  0,  0,  0, -1,  0, -1,
			      /* CU */ -1, -1, -1,  0, -1, -1, -1,  0,
			      /* GA */  0,  0,  0,  0, -1, -1, -1,
			      /* GC */  0,  0,  0,  0,  0, -1,
			      /* GG */  0, -1, -1,  0, -1,
			      /* GU */  0, -1,  0,  0,
			      /* UA */  0,  0,  0,
			      /* UC */  0,  0,
			      /* UG */  0};
			      
		
		setSymmetry(s, tr->partitionData[model].symmetryVector, 120, f, tr->partitionData[model].frequencyGrouping, 16);
			      
		tr->partitionData[model].nonGTR = TRUE;
	      }
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

    }

}


static void calculatePartitionWeights(tree *tr)
{
  if(tr->NumberOfModels > 1)    
    {
      int 
	model, 
	i;
      
      double 
	*modelWeights = (double *)rax_calloc(tr->NumberOfModels, sizeof(double)),
	wgtsum = 0.0;  
      
      for(i = 0; i < tr->cdta->endsite; i++)
	{
	  modelWeights[tr->model[i]]  += (double)tr->cdta->aliaswgt[i];
	  wgtsum                      += (double)tr->cdta->aliaswgt[i];
	}  
 	        
      for(model = 0; model < tr->NumberOfModels; model++)      		      	  	        
	tr->partitionContributions[model] = modelWeights[model] / wgtsum;             	 	    

      rax_free(modelWeights);
    }
}

void initModel(tree *tr, rawdata *rdta, cruncheddata *cdta, analdef *adef)
{  
  int 
    model, 
    j;
  
  double  
    temp;  
     
  optimizeRateCategoryInvocations = 1;      
  tr->numberOfInvariableColumns = 0;
  tr->weightOfInvariableColumns = 0;	 

  for(model = 0; model < tr->NumberOfModels; model++)
    tr->partitionData[model].numberOfCategories = 1;
  
  for (j = 0; j < tr->cdta->endsite; j++) 
    {
      tr->cdta->patrat[j] = temp = 1.0;
      tr->cdta->patratStored[j] = 1.0;
      tr->cdta->rateCategory[j] = 0;           
    } 

  for(model = 0; model < tr->NumberOfModels; model++)
    {            
      tr->partitionData[model].numberOfCategories = 1;           
      tr->partitionData[model].perSiteRates[0] = 1.0; 
      tr->partitionData[model].unscaled_perSiteRates[0] = 1.0;
    }

  updatePerSiteRates(tr, FALSE);
 

  if(tr->useK80 || tr->useHKY85)
    setupK80Symmetries(tr);

  setupSecondaryStructureSymmetries(tr);
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {               
      if(adef->useInvariant)
	{      
	  size_t
	    i;
	  
	  int 
	    count = 0,
	    total = 0,
	    states = tr->partitionData[model].states;	  


	   genericInvariant(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, 
			    getBitVector(tr->partitionData[model].dataType), 
			    getUndetermined(tr->partitionData[model].dataType), 
			    states, 
			    &(tr->numberOfInvariableColumns), 
			    &(tr->weightOfInvariableColumns)); 	  
	  
	   for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
	     {
	       if(tr->invariant[i] < states) 		
		 count += tr->cdta->aliaswgt[i];	      
	       
	       total += tr->cdta->aliaswgt[i];
	     }
	   
	  tr->partitionData[model].propInvariant = ((double)count)/((double) total);	 
	} 
    }

  initRateMatrix(tr); 

  baseFrequenciesGTR(rdta, cdta, tr);  

  calculatePartitionWeights(tr);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	k;

      tr->partitionData[model].alpha = 1.0;                
      tr->partitionData[model].brLenScaler = 1.0;

      if(tr->partitionData[model].protModels == AUTO)
	tr->partitionData[model].autoProtModels = WAG; /* initialize by WAG per default */

                     
      makeGammaCats(tr->rateHetModel, tr->partitionData[model].alpha, tr->partitionData[model].gammaRates, 4, tr->useGammaMedian, tr->partitionData[model].propInvariant);       
      
      for(k = 0; k < tr->partitionData[model].states; k++)
	tr->partitionData[model].freqExponents[k] = 0.0;

      for(j = 0; j < 4; j++)
	{
	  tr->partitionData[model].weights[j] = 0.25;
	  tr->partitionData[model].weightExponents[j] = 0.0;
	}

      initReversibleGTR(tr, model);
    }                                         

#ifdef _USE_PTHREADS
  masterBarrier(THREAD_COPY_INIT_MODEL, tr);   
#endif
}




