/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein. 
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


#ifdef __SIM_SSE3
#include <xmmintrin.h>
#include <pmmintrin.h>
/*#include <tmmintrin.h>*/
#endif

#ifdef _USE_PTHREADS
extern volatile double *reductionBuffer;
extern volatile int NumberOfThreads;
#endif

extern const unsigned int mask32[32];


static void calcDiagptableFlex(double z, int numberOfCategories, double *rptr, double *EIGN, double *diagptable, const int numStates)
{
  int 
    i, 
    l;
  
  double 
    lz, 
    lza[64];
  
  const int 
    rates = numStates - 1;
  
  assert(numStates <= 64);
  
  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);

  for(l = 0; l < rates; l++)      
    lza[l] = EIGN[l] * lz; 

  for(i = 0; i <  numberOfCategories; i++)
    {	      	       
      diagptable[i * numStates] = 1.0;

      for(l = 1; l < numStates; l++)
	diagptable[i * numStates + l] = EXP(rptr[i] * lza[l - 1]);     	          
    }        
}

static double evaluateCatFlex(int *ex1, int *ex2, int *cptr, int *wptr,
			      double *x1, double *x2, double *tipVector,
			      unsigned char *tipX1, int n, double *diagptable_start, double *vector, boolean writeVector, const boolean fastScaling, const int numStates)
{
  double   
    sum = 0.0, 
    term,
    *diagptable,  
    *left, 
    *right;
  
  int     
    i, 
    l;                           
  
  if(tipX1)
    {            
      if(writeVector)	  
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[numStates * tipX1[i]]);
	    right = &(x2[numStates * i]);
	    
	    diagptable = &diagptable_start[numStates * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < numStates; l++)
	      term += left[l] * right[l] * diagptable[l];	 	  	  
	     
	    if(fastScaling)
	      term = LOG(FABS(term));
	    else
	      term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {	       	
	    left = &(tipVector[numStates * tipX1[i]]);
	    right = &(x2[numStates * i]);
	    
	    diagptable = &diagptable_start[numStates * cptr[i]];	           	 
	    
	    for(l = 0, term = 0.0; l < numStates; l++)
	      term += left[l] * right[l] * diagptable[l];
	 	  	  
	    if(fastScaling)
	      term = LOG(FABS(term));
	    else
	      term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
	   
	    sum += wptr[i] * term;
	  }      
    }    
  else
    {
      if(writeVector)	
	for (i = 0; i < n; i++) 
	  {		       	      	      
	    left  = &x1[numStates * i];
	    right = &x2[numStates * i];
	    
	  diagptable = &diagptable_start[numStates * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < numStates; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	  
	  vector[i] = term;
	  
	  sum += wptr[i] * term;      
	}
      else
	for (i = 0; i < n; i++) 
	  {		       	      	      
	    left  = &x1[numStates * i];
	    right = &x2[numStates * i];
	    
	    diagptable = &diagptable_start[numStates * cptr[i]];	  	
	    
	    for(l = 0, term = 0.0; l < numStates; l++)
	      term += left[l] * right[l] * diagptable[l];	
	    
	    if(fastScaling)
	      term = LOG(FABS(term));
	    else
	      term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	    
	    sum += wptr[i] * term;      
	  }
    }
            
 
  
  return  sum;         
} 

static double evaluateGammaFlex_GAPPED(int *ex1, int *ex2, int *wptr,
				       double *x1_start, double *x2_start,  
				       double *tipVector, 
				       unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector, const boolean fastScaling, const int numStates,
				       double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   
    sum = 0.0, 
    term,
    *left, 
    *right,
    *x1,
    *x2;
  
  int     
    i, 
    j, 
    l; 

  const int 
    gammaStates = numStates * 4;
            
  if(tipX1)
    {          
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[numStates * tipX1[i]]);	  	  
	    
	    if(x2_gap[i / 32] & mask32[i % 32])
	      x2 = x2_gapColumn;
	    else
	      x2 = &x2_start[gammaStates * i];

	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[numStates * j]);
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	      
	      }	 
	    
	    if(fastScaling)
	      term = LOG(0.25 * FABS(term));
	    else
	      term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     
	      left = &(tipVector[numStates * tipX1[i]]);	 

	      if(x2_gap[i / 32] & mask32[i % 32])
		x2 = x2_gapColumn;
	      else
		x2 = &x2_start[gammaStates * i];

	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  right = &(x2[numStates * j]);
		  
		  for(l = 0; l < numStates; l++)
		    term += left[l] * right[l] * diagptable[j * numStates + l];	      
		}
	      
	      if(fastScaling)
		term = LOG(0.25 * FABS(term));
	      else
		term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	      
	      sum += wptr[i] * term;
	    }     	 
	}
    }              
  else
    {
      if(writeVector)
	for (i = 0; i < n; i++) 
	{	
	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1 = x1_gapColumn;
	  else
	    x1 = &x1_start[gammaStates * i];
  	 	             
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    x2 = &x2_start[gammaStates * i];


	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[numStates * j]);
	      right = &(x2[numStates * j]);	    
	      
	      for(l = 0; l < numStates; l++)
		term += left[l] * right[l] * diagptable[j * numStates + l];	
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	
	  vector[i] = term;
  
	  sum += wptr[i] * term;
	}         
      else
	for (i = 0; i < n; i++) 
	  {	  	 	             
	    if(x1_gap[i / 32] & mask32[i % 32])
	      x1 = x1_gapColumn;
	    else
	      x1 = &x1_start[gammaStates * i];
	    
	    if(x2_gap[i / 32] & mask32[i % 32])
	      x2 = x2_gapColumn;
	    else
	      x2 = &x2_start[gammaStates * i];
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		left  = &(x1[numStates * j]);
		right = &(x2[numStates * j]);	    
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	
	      }
	    
	    if(fastScaling)
	      term = LOG(0.25 * FABS(term));
	    else
	      term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	    
	    sum += wptr[i] * term;
	  }         
    }
         
  return  sum;
}



static double evaluateGammaFlex(int *ex1, int *ex2, int *wptr,
				double *x1, double *x2,  
				double *tipVector, 
				unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector, const boolean fastScaling, const int numStates)
{
  double   
    sum = 0.0, 
    term,
    *left, 
    *right;
  
  int     
    i, 
    j, 
    l; 

  const int 
    gammaStates = numStates * 4;
            
  if(tipX1)
    {          
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[numStates * tipX1[i]]);	  	  
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[gammaStates * i + numStates * j]);
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	      
	      }	 
	    
	    if(fastScaling)
	      term = LOG(0.25 * FABS(term));
	    else
	      term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     
	      left = &(tipVector[numStates * tipX1[i]]);	  	  
	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  right = &(x2[gammaStates * i + numStates * j]);
		  
		  for(l = 0; l < numStates; l++)
		    term += left[l] * right[l] * diagptable[j * numStates + l];	      
		}
	      
	      if(fastScaling)
		term = LOG(0.25 * FABS(term));
	      else
		term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	      
	      sum += wptr[i] * term;
	    }     	 
	}
    }              
  else
    {
      if(writeVector)
	for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[gammaStates * i + numStates * j]);
	      right = &(x2[gammaStates * i + numStates * j]);	    
	      
	      for(l = 0; l < numStates; l++)
		term += left[l] * right[l] * diagptable[j * numStates + l];	
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	
	  vector[i] = term;
  
	  sum += wptr[i] * term;
	}         
      else
	for (i = 0; i < n; i++) 
	  {	  	 	             
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		left  = &(x1[gammaStates * i + numStates * j]);
		right = &(x2[gammaStates * i + numStates * j]);	    
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	
	      }
	    
	    if(fastScaling)
	      term = LOG(0.25 * FABS(term));
	    else
	      term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	    
	    sum += wptr[i] * term;
	  }         
    }
         
  return  sum;
}

static double evaluateGammaFlex_perSite(int *ex1, int *ex2, int *wptr,
					double *x1, double *x2,  				       
					unsigned char *tipX1, int n, 
					int *perSiteAA,
					siteAAModels *siteProtModel,
					const boolean fastScaling, const int numStates)
{
  double   
    sum = 0.0, 
    term,
    *left, 
    *right;
  
  int     
    i, 
    j, 
    l; 

  const int 
    gammaStates = numStates * 4;
            
  if(tipX1)
    {                    
      for (i = 0; i < n; i++) 
	{
	  double 
	    *tipVector  = siteProtModel[perSiteAA[i]].tipVector,
	    *diagptable = siteProtModel[perSiteAA[i]].left;
	     
	  left = &(tipVector[numStates * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[gammaStates * i + numStates * j]);
		  
	      for(l = 0; l < numStates; l++)
		term += left[l] * right[l] * diagptable[j * numStates + l];	      
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}     	       
    }              
  else
    {
     
      for (i = 0; i < n; i++) 
	{	  
	  double 	   
	    *diagptable = siteProtModel[perSiteAA[i]].left;	 	             	  

	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[gammaStates * i + numStates * j]);
	      right = &(x2[gammaStates * i + numStates * j]);	    
	      
	      for(l = 0; l < numStates; l++)
		term += left[l] * right[l] * diagptable[j * numStates + l];	
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
  
  return  sum;
}


static double evaluateGammaInvarFlex (int *ex1, int *ex2, int *wptr, int *iptr,
				      double *x1, double *x2, 
				      double *tipVector,double *tFreqs, double invariants,
				      unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector, const boolean fastScaling,
				      const int numStates)
{
  double   
    sum = 0.0, 
    term, 
    freqs[64],
    scaler = 0.25 * (1.0 - invariants),
    *left,
    *right;
  
  int     
    i, 
    j, 
    l;     
  
  const int gammaStates = numStates * 4;
    
  for(i = 0; i < numStates; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {    
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[numStates * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[gammaStates * i + numStates * j]);
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	      
	      }
	    
	    if(iptr[i] < numStates)
	      if(fastScaling) 
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	      else
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + ex2[i] * LOG(minlikelihood);
	    else
	      if(fastScaling)
		term = LOG(scaler * FABS(term));
	      else
		term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));
	   	    	    
	    vector[i] = term;
	   
	    sum += wptr[i] * term;
	  }         
      else
	for (i = 0; i < n; i++) 
	  {
	    left = &(tipVector[numStates * tipX1[i]]);
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		right = &(x2[gammaStates * i + numStates * j]);
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	      
	      }	  
	    
	    if(iptr[i] < numStates)
	      if(fastScaling) 
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	      else
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + ex2[i] * LOG(minlikelihood);
	    else
	      if(fastScaling)
		term = LOG(scaler * FABS(term));
	      else
		term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));
	    	    
	    sum += wptr[i] * term;
	  }    	
    }                
  else
    {    
      if(writeVector)
	for (i = 0; i < n; i++) 
	  {	  	 	       	  
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		left  = &(x1[gammaStates * i + numStates * j]);
		right = &(x2[gammaStates * i + numStates * j]);	    
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	
	      }
	    
	    if(iptr[i] < numStates)
	      if(fastScaling) 
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	      else
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + (ex1[i] + ex2[i]) * LOG(minlikelihood);
	    else
	      if(fastScaling)
		term = LOG(scaler * FABS(term));
	      else
		term = LOG(scaler * FABS(term)) + (ex1[i] + ex2[i]) * LOG(minlikelihood);	  	 	
	    
	    vector[i] = term;

	    sum += wptr[i] * term;
	  }   
      else
	for (i = 0; i < n; i++) 
	  {	  	 	       	  
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		left  = &(x1[gammaStates * i + numStates * j]);
		right = &(x2[gammaStates * i + numStates * j]);	    
		
		for(l = 0; l < numStates; l++)
		  term += left[l] * right[l] * diagptable[j * numStates + l];	
	      }
	    
	    if(iptr[i] < numStates)
	      if(fastScaling) 
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	      else
		term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + (ex1[i] + ex2[i]) * LOG(minlikelihood);
	    else
	      if(fastScaling)
		term = LOG(scaler * FABS(term));
	      else
		term = LOG(scaler * FABS(term)) + (ex1[i] + ex2[i]) * LOG(minlikelihood);	  	 	
	    
	    sum += wptr[i] * term;
	  }              
    }
       
  return  sum;
}



void calcDiagptable(double z, int data, int numberOfCategories, double *rptr, double *EIGN, double *diagptable)
{
  int i, l;
  double lz;

  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);

  switch(data)
    {
    case BINARY_DATA:
       {
	double lz1;
	lz1 = EIGN[0] * lz;
	for(i = 0; i <  numberOfCategories; i++)
	  {		 
	    diagptable[2 * i] = 1.0;
	    diagptable[2 * i + 1] = EXP(rptr[i] * lz1);	   	    
	  }
      }
      break;
    case DNA_DATA:
      {
	double lz1, lz2, lz3;
	lz1 = EIGN[0] * lz;
	lz2 = EIGN[1] * lz;
	lz3 = EIGN[2] * lz;
       
	for(i = 0; i <  numberOfCategories; i++)
	  {		 	    
	    diagptable[4 * i] = 1.0;
	    diagptable[4 * i + 1] = EXP(rptr[i] * lz1);
	    diagptable[4 * i + 2] = EXP(rptr[i] * lz2);
	    diagptable[4 * i + 3] = EXP(rptr[i] * lz3);	    
	  }
      }
      break;
    case AA_DATA:
      {
	double lza[19];

	for(l = 0; l < 19; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 20] = 1.0;

	    for(l = 1; l < 20; l++)
	      diagptable[i * 20 + l] = EXP(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    case SECONDARY_DATA:
      {
	double lza[15];

	for(l = 0; l < 15; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 16] = 1.0;

	    for(l = 1; l < 16; l++)
	      diagptable[i * 16 + l] = EXP(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    case SECONDARY_DATA_6:
      {
	double lza[5];

	for(l = 0; l < 5; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 6] = 1.0;

	    for(l = 1; l < 6; l++)
	      diagptable[i * 6 + l] = EXP(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    case SECONDARY_DATA_7:
      {
	double lza[6];

	for(l = 0; l < 6; l++)      
	  lza[l] = EIGN[l] * lz; 

	for(i = 0; i <  numberOfCategories; i++)
	  {	      	       
	    diagptable[i * 7] = 1.0;

	    for(l = 1; l < 7; l++)
	      diagptable[i * 7 + l] = EXP(rptr[i] * lza[l - 1]);     	          
	  }
      }
      break;
    default:
      assert(0);
    }
}









static double evaluateGTRCATPROT (int *ex1, int *ex2, int *cptr, int *wptr,
				  double *x1, double *x2, double *tipVector,
				  unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {                 
      for (i = 0; i < n; i++) 
	{	       	
	  left = &(tipVector[20 * tipX1[i]]);
	  right = &(x2[20 * i]);
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	           	 
#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();	    
	  
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d lv = _mm_load_pd(&left[l]);
	      __m128d rv = _mm_load_pd(&right[l]);
	      __m128d mul = _mm_mul_pd(lv, rv);
	      __m128d dv = _mm_load_pd(&diagptable[l]);
	      
	      tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));		   
	    }		 		
	  
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
#else  
	  for(l = 0, term = 0.0; l < 20; l++)
	    term += left[l] * right[l] * diagptable[l];	 	  	  
#endif	    
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  
	  diagptable = &diagptable_start[20 * cptr[i]];	  	
#ifdef __SIM_SSE3
	    __m128d tv = _mm_setzero_pd();	    
	      	    
	    for(l = 0; l < 20; l+=2)
	      {
		__m128d lv = _mm_load_pd(&left[l]);
		__m128d rv = _mm_load_pd(&right[l]);
		__m128d mul = _mm_mul_pd(lv, rv);
		__m128d dv = _mm_load_pd(&diagptable[l]);
		
		tv = _mm_add_pd(tv, _mm_mul_pd(mul, dv));		   
	      }		 		
	      
	      tv = _mm_hadd_pd(tv, tv);
	      _mm_storel_pd(&term, tv);
#else  
	  for(l = 0, term = 0.0; l < 20; l++)
	    term += left[l] * right[l] * diagptable[l];	
#endif

	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 


static double evaluateGTRCATSECONDARY (int *ex1, int *ex2, int *cptr, int *wptr,
				       double *x1, double *x2, double *tipVector,
				       unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {                  
      for (i = 0; i < n; i++) 
	{	       	
	  left = &(tipVector[16 * tipX1[i]]);
	  right = &(x2[16 * i]);
	  
	  diagptable = &diagptable_start[16 * cptr[i]];	           	 
	  
	  for(l = 0, term = 0.0; l < 16; l++)
	    term += left[l] * right[l] * diagptable[l];
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[16 * i];
	  right = &x2[16 * i];
	  
	  diagptable = &diagptable_start[16 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 16; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCATSECONDARY_6 (int *ex1, int *ex2, int *cptr, int *wptr,
				       double *x1, double *x2, double *tipVector,
				       unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {                       
      for (i = 0; i < n; i++) 
	{	       	
	  left = &(tipVector[6 * tipX1[i]]);
	  right = &(x2[6 * i]);
	  
	  diagptable = &diagptable_start[6 * cptr[i]];	           	 
	  
	  for(l = 0, term = 0.0; l < 6; l++)
	    term += left[l] * right[l] * diagptable[l];
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[6 * i];
	  right = &x2[6 * i];
	  
	  diagptable = &diagptable_start[6 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 6; l++)
	    term += left[l] * right[l] * diagptable[l];	
	  
	   if(fastScaling)
	     term = LOG(FABS(term));
	   else
	     term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCATSECONDARY_7(int *ex1, int *ex2, int *cptr, int *wptr,
					double *x1, double *x2, double *tipVector,
					unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double   sum = 0.0, term;
  double  *diagptable,  *left, *right;
  int     i, l;                           
  
  if(tipX1)
    {                 
      for (i = 0; i < n; i++) 
	{	       	
	  left = &(tipVector[7 * tipX1[i]]);
	  right = &(x2[7 * i]);
	  
	  diagptable = &diagptable_start[7 * cptr[i]];	           	 
	  
	  for(l = 0, term = 0.0; l < 7; l++)
	    term += left[l] * right[l] * diagptable[l];	 	  	  
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}      
    }    
  else
    {
    
      for (i = 0; i < n; i++) 
	{		       	      	      
	  left  = &x1[7 * i];
	  right = &x2[7 * i];
	  
	  diagptable = &diagptable_start[7 * cptr[i]];	  	
	  
	  for(l = 0, term = 0.0; l < 7; l++)
	    term += left[l] * right[l] * diagptable[l];	

	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;      
	}
    }
             
  return  sum;         
} 

static double evaluateGTRCAT_BINARY (int *ex1, int *ex2, int *cptr, int *wptr,
				     double *x1_start, double *x2_start, double *tipVector, 		      
				     unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double  sum = 0.0, term;       
  int     i;
#ifndef  __SIM_SSE3
  int j;  
#endif
  double  *diagptable, *x1, *x2;                      	    
 
  if(tipX1)
    {          
      for (i = 0; i < n; i++) 
	{
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));	    		   	  
#endif
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &(x2_start[2 * i]);
	  
	  diagptable = &(diagptable_start[2 * cptr[i]]);	    	    	  
	
#ifdef __SIM_SSE3	  
	  _mm_store_pd(t, _mm_mul_pd(_mm_load_pd(x1), _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(diagptable))));
	  
	  if(fastScaling)
	    term = LOG(FABS(t[0] + t[1]));
	  else
	    term = LOG(FABS(t[0] + t[1])) + (ex2[i] * LOG(minlikelihood));			     
#else	 	    
	  for(j = 0, term = 0.0; j < 2; j++)       	   	     
	    term += x1[j] * x2[j] * diagptable[j];	      
	  	 
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));	   	    	   	 	  	  	 
#endif	  

	  sum += wptr[i] * term;
	}	
    }               
  else
    {
      for (i = 0; i < n; i++) 
	{	
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));	    		   	  
#endif	          	
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  
	  diagptable = &diagptable_start[2 * cptr[i]];		  
#ifdef __SIM_SSE3	  
	  _mm_store_pd(t, _mm_mul_pd(_mm_load_pd(x1), _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(diagptable))));
	  
	  if(fastScaling)
	    term = LOG(FABS(t[0] + t[1]));
	  else
	    term = LOG(FABS(t[0] + t[1])) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));			     
#else	  
	  for(j = 0, term = 0.0; j < 2; j++)
	    term += x1[j] * x2[j] * diagptable[j];   
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
#endif
	  
	  sum += wptr[i] * term;
	}	   
    }
       
  return  sum;         
} 


static double evaluateGTRGAMMA_BINARY(int *ex1, int *ex2, int *wptr,
				      double *x1_start, double *x2_start, 
				      double *tipVector, 
				      unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;    
  int     i, j;
#ifndef __SIM_SSE3
  int k;
#endif 
  double  *x1, *x2;             

  if(tipX1)
    {          
      for (i = 0; i < n; i++)
	{
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;
#endif
	  x1 = &(tipVector[2 * tipX1[i]]);	 
	  x2 = &x2_start[8 * i];	          	  	
#ifdef __SIM_SSE3	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[0]);
	      x2v = _mm_load_pd(&x2[j * 2]);
	      dv   = _mm_load_pd(&diagptable[j * 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);	      	      
	    }
	  
	  _mm_store_pd(t, termv);	        
	  
	  if(fastScaling)
	    term = LOG(0.25 * (FABS(t[0] + t[1])));
	  else
	    term = LOG(0.25 * (FABS(t[0] + t[1]))) + (ex2[i] * LOG(minlikelihood));	  
#else
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];	          	  	  	    	    
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ex2[i] * LOG(minlikelihood);
#endif	 
	  
	  sum += wptr[i] * term;
	}	  
    }
  else
    {         
      for (i = 0; i < n; i++) 
	{
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;
#endif	  	 	  	  
	  x1 = &x1_start[8 * i];
	  x2 = &x2_start[8 * i];
	  	  
#ifdef __SIM_SSE3	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[j * 2]);
	      x2v = _mm_load_pd(&x2[j * 2]);
	      dv   = _mm_load_pd(&diagptable[j * 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);	      	      
	    }
	  
	  _mm_store_pd(t, termv);
	  
	  
	  if(fastScaling)
	    term = LOG(0.25 * (FABS(t[0] + t[1])));
	  else
	    term = LOG(0.25 * (FABS(t[0] + t[1]))) + ((ex1[i] +ex2[i]) * LOG(minlikelihood));	  
#else	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      term += x1[j * 2 + k] * x2[j * 2 + k] * diagptable[j * 2 + k];	          	  	  	      

	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex1[i] + ex2[i]) * LOG(minlikelihood);
#endif

	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 

static double evaluateGTRGAMMAINVAR_BINARY (int *ex1, int *ex2, int *wptr, int *iptr,
					    double *x1_start, double *x2_start,
					    double *tipVector, double *tFreqs, double invariants,
					    unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{ 
  int     i, j, k;
  double  *x1, *x2; 
  double 
    freqs[2], 
    scaler = 0.25 * (1.0 - invariants),
    sum = 0.0, 
    term; 

  freqs[0] = tFreqs[0] * invariants; 
  freqs[1] = tFreqs[1] * invariants; 

  if(tipX1)
    {         
      for (i = 0; i < n; i++) 
	{
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &x2_start[8 * i];	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      term += x1[k] * x2[j * 2 + k] * diagptable[j * 2 + k];
	  
	  if(iptr[i] < 2)
	    if(fastScaling)	   
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + ex2[i] * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));	 
	  
	  sum += wptr[i] * term;
	}	  
    }
  else
    {           		

      for (i = 0; i < n; i++) 
	{	  	 	  	
	  x1 = &x1_start[8 * i];
	  x2 = &x2_start[8 * i];	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      term += x1[j * 2 + k] * x2[j * 2 + k] * diagptable[j * 2 + k];	  	 	      
	  
	  if(iptr[i] < 2)
	    if(fastScaling)	   	       
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + (ex2[i] + ex1[i]) * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));

	  sum += wptr[i] * term;
	}	  	                        
    }

  return  sum;
} 


static double evaluateGTRCAT (int *ex1, int *ex2, int *cptr, int *wptr,
			      double *x1_start, double *x2_start, double *tipVector, 		      
			      unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling)
{
  double  sum = 0.0, term;       
  int     i;
#ifndef __SIM_SSE3
  int j;  
#endif
  double  *diagptable, *x1, *x2;                      	    
 
  if(tipX1)
    {           
      for (i = 0; i < n; i++) 
	{	
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;
#endif
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[4 * cptr[i]];
	  
#ifdef __SIM_SSE3	    	  
	  x1v1 =  _mm_load_pd(&x1[0]);
	  x1v2 =  _mm_load_pd(&x1[2]);
	  x2v1 =  _mm_load_pd(&x2[0]);
	  x2v2 =  _mm_load_pd(&x2[2]);
	  dv1  =  _mm_load_pd(&diagptable[0]);
	  dv2  =  _mm_load_pd(&diagptable[2]);
	  
	  x1v1 = _mm_mul_pd(x1v1, x2v1);
	  x1v1 = _mm_mul_pd(x1v1, dv1);
	  
	  x1v2 = _mm_mul_pd(x1v2, x2v2);
	  x1v2 = _mm_mul_pd(x1v2, dv2);
	  
	  x1v1 = _mm_add_pd(x1v1, x1v2);
	  
	  _mm_store_pd(t, x1v1);
	  
	  if(fastScaling)
	    term = LOG(FABS(t[0] + t[1]));
	  else
	    term = LOG(FABS(t[0] + t[1])) + (ex2[i] * LOG(minlikelihood));
#else
	  for(j = 0, term = 0.0; j < 4; j++)
	    term += x1[j] * x2[j] * diagptable[j];
	  
	  /*{
	    double 
	      term[4],
	      sum = 0.0;
	    
	    for(j = 0; j < 4; j++)
	      {
		term[j] = ABS(x1[j] * x2[j] * diagptable[j]);
		sum += term[j];
	      }

	    printf("RRRRRRR %1.80f %1.80f %1.80f %1.80f\n", term[0]/sum, term[1]/sum, term[2]/sum, term[3]/sum);
	    }*/

	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));
#endif	    
	  sum += wptr[i] * term;
	}	
    }               
  else
    {
      for (i = 0; i < n; i++) 
	{ 
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));
	   __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;
#endif
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[4 * cptr[i]];	
	  
#ifdef __SIM_SSE3	  
	  x1v1 =  _mm_load_pd(&x1[0]);
	  x1v2 =  _mm_load_pd(&x1[2]);
	  x2v1 =  _mm_load_pd(&x2[0]);
	  x2v2 =  _mm_load_pd(&x2[2]);
	  dv1  =  _mm_load_pd(&diagptable[0]);
	  dv2  =  _mm_load_pd(&diagptable[2]);
	  
	  x1v1 = _mm_mul_pd(x1v1, x2v1);
	  x1v1 = _mm_mul_pd(x1v1, dv1);
	  
	  x1v2 = _mm_mul_pd(x1v2, x2v2);
	  x1v2 = _mm_mul_pd(x1v2, dv2);
	  
	  x1v1 = _mm_add_pd(x1v1, x1v2);
	  
	  _mm_store_pd(t, x1v1);
	  
	  if(fastScaling)
	    term = LOG(FABS(t[0] + t[1]));
	  else
	    term = LOG(FABS(t[0] + t[1])) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
#else
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    term += x1[j] * x2[j] * diagptable[j];     
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));	  
#endif
	  sum += wptr[i] * term;
	}    
    }
       
  return  sum;         
} 


#ifdef __SIM_SSE3

static double evaluateGTRGAMMA_GAPPED(int *ex1, int *ex2, int *wptr,
				      double *x1_start, double *x2_start, 
				      double *tipVector, 
				      unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling,
				      double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   sum = 0.0, term;    
  int     i, j;
  double  *x1, *x2;             

 

  if(tipX1)
    {          	
      for (i = 0; i < n; i++)
	{
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;

	  x1 = &(tipVector[4 * tipX1[i]]);	 
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    x2 = &x2_start[16 * i];	 
	  
	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[0]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);	  	 

	  if(fastScaling)
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	  else
	    term = LOG(0.25 * FABS(t[0] + t[1])) + (ex2[i] * LOG(minlikelihood));	  
	  
	  sum += wptr[i] * term;
	}     
    }
  else
    {        
      for (i = 0; i < n; i++) 
	{

	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;

	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1 = x1_gapColumn;
	  else
	    x1 = &x1_start[16 * i]; 	  	  
	 	      
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    x2 = &x2_start[16 * i];
	
	  termv = _mm_set1_pd(0.0);	  	 
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[j * 4]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[j * 4 + 2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);

	  if(fastScaling)
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	  else
	    term = LOG(0.25 * FABS(t[0] + t[1])) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));	  
	  
	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 


static double evaluateGTRGAMMA_GAPPED_SAVE(int *ex1, int *ex2, int *wptr,
					   double *x1_start, double *x2_start, 
					   double *tipVector, 
					   unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling,
					   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   sum = 0.0, term;    
  int     i, j;
  double  
    *x1, 
    *x2,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;

 

  if(tipX1)
    {        
     
      
      for (i = 0; i < n; i++)
	{
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;

	  x1 = &(tipVector[4 * tipX1[i]]);	 
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    {
	      x2 = x2_ptr;	 
	      x2_ptr += 16;
	    }
	  
	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[0]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);	  	 

	  if(fastScaling)
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	  else
	    term = LOG(0.25 * FABS(t[0] + t[1])) + (ex2[i] * LOG(minlikelihood));	  
	  
	  sum += wptr[i] * term;
	}     
    }
  else
    {        
      
      for (i = 0; i < n; i++) 
	{

	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;

	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1 = x1_gapColumn;
	  else
	    {
	      x1 = x1_ptr; 	  	  
	      x1_ptr += 16;
	    }
	 	      
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    {
	      x2 = x2_ptr;
	      x2_ptr += 16;
	    }
	
	  termv = _mm_set1_pd(0.0);	  	 
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[j * 4]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[j * 4 + 2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);

	  if(fastScaling)
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	  else
	    term = LOG(0.25 * FABS(t[0] + t[1])) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));	  
	  
	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 

#else

static double evaluateGTRGAMMA_GAPPED(int *ex1, int *ex2, int *wptr,
				      double *x1_start, double *x2_start, 
				      double *tipVector, 
				      unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling,
				      double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   sum = 0.0, term;    
  int     i, j, k;
  double  *x1, *x2;             

 

  if(tipX1)
    {          	
      for (i = 0; i < n; i++)
	{

	  x1 = &(tipVector[4 * tipX1[i]]);
	  
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    x2 = &x2_start[16 * i];	 	  

	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];	          	  	  	    	    	  
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ex2[i] * LOG(minlikelihood);	 
	  
	  sum += wptr[i] * term;
	}     
    }
  else
    {        
      for (i = 0; i < n; i++) 
	{	  	 	  	  
	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1 = x1_gapColumn;
	  else
	    x1 = &x1_start[16 * i]; 	  	  
	 	      
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    x2 = &x2_start[16 * i];
	

	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[j * 4 + k] * x2[j * 4 + k] * diagptable[j * 4 + k];
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex1[i] + ex2[i]) * LOG(minlikelihood);
	  
	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 


#endif

static double evaluateGTRGAMMA(int *ex1, int *ex2, int *wptr,
			       double *x1_start, double *x2_start, 
			       double *tipVector, 
			       unsigned char *tipX1, const int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;    
  int     i, j;
#ifndef __SIM_SSE3  
  int k;
#endif
  double  *x1, *x2;             

 

  if(tipX1)
    {          	
      for (i = 0; i < n; i++)
	{
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;
#endif
	  x1 = &(tipVector[4 * tipX1[i]]);	 
	  x2 = &x2_start[16 * i];	 
	  
#ifdef __SIM_SSE3	
	  termv = _mm_set1_pd(0.0);	    	   
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[0]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);
	  
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	  else
	    term = LOG(0.25 * FABS(t[0] + t[1])) + (ex2[i] * LOG(minlikelihood));	  
#else
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];	          	  	  	    	    	  
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ex2[i] * LOG(minlikelihood);	 
#endif
	  
	  sum += wptr[i] * term;
	}     
    }
  else
    {        
      for (i = 0; i < n; i++) 
	{
#ifdef __SIM_SSE3
	  double t[2] __attribute__ ((aligned (16)));
	  __m128d termv, x1v, x2v, dv;
#endif
	  	 	  	  
	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];	  	  
	
#ifdef __SIM_SSE3	
	  termv = _mm_set1_pd(0.0);	  	 
	  
	  for(j = 0; j < 4; j++)
	    {
	      x1v = _mm_load_pd(&x1[j * 4]);
	      x2v = _mm_load_pd(&x2[j * 4]);
	      dv   = _mm_load_pd(&diagptable[j * 4]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	      
	      x1v = _mm_load_pd(&x1[j * 4 + 2]);
	      x2v = _mm_load_pd(&x2[j * 4 + 2]);
	      dv   = _mm_load_pd(&diagptable[j * 4 + 2]);
	      
	      x1v = _mm_mul_pd(x1v, x2v);
	      x1v = _mm_mul_pd(x1v, dv);
	      
	      termv = _mm_add_pd(termv, x1v);
	    }
	  
	  _mm_store_pd(t, termv);

	  if(fastScaling)
	    term = LOG(0.25 * FABS(t[0] + t[1]));
	  else
	    term = LOG(0.25 * FABS(t[0] + t[1])) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));	  
#else 
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[j * 4 + k] * x2[j * 4 + k] * diagptable[j * 4 + k];
	          	  	  	      
	   if(fastScaling)
	     term = LOG(0.25 * FABS(term));
	    else
	      term = LOG(0.25 * FABS(term)) + (ex1[i] + ex2[i]) * LOG(minlikelihood);
#endif
	  
	  sum += wptr[i] * term;
	}                      	
    }

  return sum;
} 









static double evaluateGTRGAMMAINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
				     double *x1_start, double *x2_start,
				     double *tipVector, double *tFreqs, double invariants,
				     unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{ 
  int     i, j, k;
  double  *x1, *x2; 
  double 
    freqs[4], 
    scaler = 0.25 * (1.0 - invariants),
    sum = 0.0, 
    term; 

  freqs[0] = tFreqs[0] * invariants; 
  freqs[1] = tFreqs[1] * invariants;
  freqs[2] = tFreqs[2] * invariants;
  freqs[3] = tFreqs[3] * invariants;   

  if(tipX1)
    {         
      for (i = 0; i < n; i++) 
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[16 * i];	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[k] * x2[j * 4 + k] * diagptable[j * 4 + k];
	  
	  if(iptr[i] < 4)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + ex2[i] * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));	 
	  
	  sum += wptr[i] * term;
	}	  
    }
  else
    {           		

      for (i = 0; i < n; i++) 
	{	  	 	  	
	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      term += x1[j * 4 + k] * x2[j * 4 + k] * diagptable[j * 4 + k];	  	 	      
	  
	  if(iptr[i] < 4)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + (ex2[i] + ex1[i]) * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));

	  sum += wptr[i] * term;
	}	  	                        
    }

  return  sum;
} 




static double evaluateGTRGAMMAPROT (int *ex1, int *ex2, int *wptr,
				    double *x1, double *x2,  
				    double *tipVector, 
				    unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{
#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      right = &(x2[80 * i + 20 * j]);
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  
#else
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[80 * i + 20 * j]);
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	      
	    }	  
#endif
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}    	        
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();	 	  	  
	      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);	  
#else
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	
	    }
#endif
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}


static double evaluateGTRGAMMAPROT_GAPPED (int *ex1, int *ex2, int *wptr,
					   double *x1, double *x2,  
					   double *tipVector, 
					   unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling,
					   double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)					   
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  
    *left, 
    *right,
    *x1v,
    *x2v;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2v = x2_gapColumn;
	  else
	    x2v = &x2[80 * i];

#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      right = &(x2v[20 * j]);
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  
#else
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2v[20 * j]);
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	      
	    }	  
#endif	 
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}    	        
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{
	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1v = x1_gapColumn;
	  else
	    x1v = &x1[80 * i];

	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2v = x2_gapColumn;
	  else
	    x2v = &x2[80 * i];
	  	 	             
#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();	 	  	  
	      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      left  = &(x1v[20 * j]);
	      right = &(x2v[20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);	  
#else
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1v[20 * j]);
	      right = &(x2v[20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	
	    }
#endif
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMAPROT_GAPPED_SAVE (int *ex1, int *ex2, int *wptr,
						double *x1, double *x2,  
						double *tipVector, 
						unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling,
						double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)					   
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  
    *left, 
    *right,
    *x1v,
    *x2v;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2v = x2_gapColumn;
	  else
	    x2v = &x2[80 * i];

#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      right = &(x2v[20 * j]);
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);
	  
#else
	  left = &(tipVector[20 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2v[20 * j]);
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	      
	    }	  
#endif	 
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}    	        
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{
	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1v = x1_gapColumn;
	  else
	    x1v = &x1[80 * i];

	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2v = x2_gapColumn;
	  else
	    x2v = &x2[80 * i];
	  	 	             
#ifdef __SIM_SSE3
	  __m128d tv = _mm_setzero_pd();	 	  	  
	      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double *d = &diagptable[j * 20];
	      left  = &(x1v[20 * j]);
	      right = &(x2v[20 * j]);
	      
	      for(l = 0; l < 20; l+=2)
		{
		  __m128d mul = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
		  tv = _mm_add_pd(tv, _mm_mul_pd(mul, _mm_load_pd(&d[l])));		   
		}		 		
	    }
	  tv = _mm_hadd_pd(tv, tv);
	  _mm_storel_pd(&term, tv);	  
#else
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1v[20 * j]);
	      right = &(x2v[20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	
	    }
#endif
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}





static double evaluateGTRGAMMASECONDARY (int *ex1, int *ex2, int *wptr,
					 double *x1, double *x2,  
					 double *tipVector, 
					 unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{	     
	  left = &(tipVector[16 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[64 * i + 16 * j]);
	      
	      for(l = 0; l < 16; l++)
		term += left[l] * right[l] * diagptable[j * 16 + l];	      
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}     	     
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[64 * i + 16 * j]);
	      right = &(x2[64 * i + 16 * j]);	    
	      
	      for(l = 0; l < 16; l++)
		term += left[l] * right[l] * diagptable[j * 16 + l];	
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARY_6 (int *ex1, int *ex2, int *wptr,
					   double *x1, double *x2,  
					   double *tipVector, 
					   unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {                
      for (i = 0; i < n; i++) 
	{	     
	  left = &(tipVector[6 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[24 * i + 6 * j]);
	      
	      for(l = 0; l < 6; l++)
		term += left[l] * right[l] * diagptable[j * 6 + l];	      
	    }	
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}     	      
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[24 * i + 6 * j]);
	      right = &(x2[24 * i + 6 * j]);	    
	      
	      for(l = 0; l < 6; l++)
		term += left[l] * right[l] * diagptable[j * 6 + l];	
	    }

	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARY_7 (int *ex1, int *ex2, int *wptr,
					   double *x1, double *x2,  
					   double *tipVector, 
					   unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{	     
	  left = &(tipVector[7 * tipX1[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[28 * i + 7 * j]);
	      
	      for(l = 0; l < 7; l++)
		term += left[l] * right[l] * diagptable[j * 7 + l];	      
	    }	
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}     	    
    }              
  else
    {
      for (i = 0; i < n; i++) 
	{	  	 	             
      
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[28 * i + 7 * j]);
	      right = &(x2[28 * i + 7 * j]);	    
	      
	      for(l = 0; l < 7; l++)
		term += left[l] * right[l] * diagptable[j * 7 + l];	
	    }
	  
	  if(fastScaling)
	    term = LOG(0.25 * FABS(term));
	  else
	    term = LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}

static double evaluateGTRGAMMAPROTINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
					 double *x1, double *x2, 
					 double *tipVector,double *tFreqs, double invariants,
					 unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   
    sum = 0.0, term, freqs[20],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 20; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {         
      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[20 * tipX1[i]]);
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[80 * i + 20 * j]);
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	      
	    }	  
	  
	  if(iptr[i] < 20)	
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + ex2[i] * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		term += left[l] * right[l] * diagptable[j * 20 + l];	
	    }
	  
	  if(iptr[i] < 20)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + (ex1[i] + ex2[i]) * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARYINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
					      double *x1, double *x2, 
					      double *tipVector,double *tFreqs, double invariants,
					      unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   
    sum = 0.0, term, freqs[16],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 16; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {         
      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[16 * tipX1[i]]);
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[64 * i + 16 * j]);
	      
	      for(l = 0; l < 16; l++)
		term += left[l] * right[l] * diagptable[j * 16 + l];	      
	    }	  
	  
	  if(iptr[i] < 16)
	    if(fastScaling) 
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + ex2[i] * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[64 * i + 16 * j]);
	      right = &(x2[64 * i + 16 * j]);	    
	      
	      for(l = 0; l < 16; l++)
		term += left[l] * right[l] * diagptable[j * 16 + l];	
	    }

	  if(iptr[i] < 16)
	    if(fastScaling) 
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + (ex1[i] + ex2[i]) * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex1[i] + ex2[i]) * LOG(minlikelihood);	  	 	
	  
	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARYINVAR_6 (int *ex1, int *ex2, int *wptr, int *iptr,
						double *x1, double *x2, 
						double *tipVector,double *tFreqs, double invariants,
						unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   
    sum = 0.0, term, freqs[6],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 6; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {         
      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[6 * tipX1[i]]);
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[24 * i + 6 * j]);
	      
	      for(l = 0; l < 6; l++)
		term += left[l] * right[l] * diagptable[j * 6 + l];	      
	    }	  
	  
	  if(iptr[i] < 6)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + ex2[i] * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[24 * i + 6 * j]);
	      right = &(x2[24 * i + 6 * j]);	    
	      
	      for(l = 0; l < 6; l++)
		term += left[l] * right[l] * diagptable[j * 6 + l];	
	    }
	  
	  if(iptr[i] < 6)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + (ex2[i] + ex1[i]) * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));

	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}

static double evaluateGTRGAMMASECONDARYINVAR_7 (int *ex1, int *ex2, int *wptr, int *iptr,
						double *x1, double *x2, 
						double *tipVector,double *tFreqs, double invariants,
						unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling)
{
  double   
    sum = 0.0, term, freqs[7],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j, l;     
  double *left, *right;   
    
  for(i = 0; i < 7; i++)
    freqs[i] = tFreqs[i] * invariants;            	  
  
  if(tipX1)
    {          
      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[7 * tipX1[i]]);
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[28 * i + 7 * j]);
	      
	      for(l = 0; l < 7; l++)
		term += left[l] * right[l] * diagptable[j * 7 + l];	      
	    }	  
	  
	  if(iptr[i] < 7)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]))  + ex2[i] * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + (ex2[i] * LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}    	
    }                
  else
    {    
      for (i = 0; i < n; i++) 
	{	  	 	       	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[28 * i + 7 * j]);
	      right = &(x2[28 * i + 7 * j]);	    
	      
	      for(l = 0; l < 7; l++)
		term += left[l] * right[l] * diagptable[j * 7 + l];	
	    }
	  
	  if(iptr[i] < 7)
	    if(fastScaling)
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]]));
	    else
	      term = LOG(((scaler * FABS(term)) + freqs[iptr[i]])) + (ex2[i] + ex1[i]) * LOG(minlikelihood);
	  else
	    if(fastScaling)
	      term = LOG(scaler * FABS(term));
	    else
	      term = LOG(scaler * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));

	  sum += wptr[i] * term;
	}              
    }
       
  return  sum;
}


double evaluateIterative(tree *tr,  boolean writeVector)
{
  double 
    *pz = tr->td[0].ti[0].qz,
    result = 0.0;   

  int 
    rateHet = tr->discreteRateCategories,
    pNumber = tr->td[0].ti[0].pNumber, 
    qNumber = tr->td[0].ti[0].qNumber, 
    model;
 
  newviewIterative(tr); 

  if(writeVector)
    assert(!tr->useFastScaling);

  for(model = 0; model < tr->NumberOfModels; model++)
    {            
      if(tr->executeModel[model])
	{	
	  int 
	    width = tr->partitionData[model].width,
	    states = tr->partitionData[model].states;
	  
	  double 
	    z, 
	    partitionLikelihood = 0.0, 
	    *_vector;
	  
	  int    
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL;
	  unsigned int
	    *x1_gap = (unsigned int*)NULL,
	    *x2_gap = (unsigned int*)NULL;

	  double 
	    *x1_start   = (double*)NULL, 
	    *x2_start   = (double*)NULL,
	    *diagptable = (double*)NULL,
	    *x1_gapColumn = (double*)NULL,
	    *x2_gapColumn = (double*)NULL;;

	 

	  unsigned char 
	    *tip = (unsigned char*)NULL;

	  if(writeVector)
	    _vector = tr->partitionData[model].perSiteLL;
	  else
	    _vector = (double*)NULL;

	 
	  diagptable = tr->partitionData[model].left;


	  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
	    {	        	    
	      if(isTip(qNumber, tr->mxtips))
		{			  		  
		  x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips -1];
		  
		  if(!tr->useFastScaling)
		    ex2      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];

		  if(tr->useGappedImplementation || tr->saveMemory)
		    {
		      x2_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
		      x2_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
		    }

		  tip = tr->partitionData[model].yVector[qNumber];	 	      
		}           
	      else
		{
		  
		  x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
		  

		  if(!tr->useFastScaling)
		    ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];

		  if(tr->useGappedImplementation || tr->saveMemory)
		    {
		      x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
		      x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
		    }
		  
		  tip = tr->partitionData[model].yVector[pNumber];
		}
	    }
	  else
	    {  
	      if(tr->useGappedImplementation || tr->saveMemory)
	      {
		x1_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
		x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
		x1_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
		x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
	      }
	      
	      x1_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
	      x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
	
	      if(!tr->useFastScaling)
		{
		  ex1      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];
		  ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];     
		}
	    }


	  if(tr->multiBranch)
	    z = pz[model];
	  else
	    z = pz[0];

	  if(writeVector)
	    {
	     
	      
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    calcDiagptableFlex(z, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable, states);
		    
		    partitionLikelihood = evaluateCatFlex(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							  x1_start, x2_start, tr->partitionData[model].tipVector,
							  tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);
		  }	     	      
		  break;	      
		case GAMMA:
		  {		    
		    assert(!tr->estimatePerSiteAA);

		    calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, states);
		     
		    if(tr->useGappedImplementation)
		      partitionLikelihood = evaluateGammaFlex_GAPPED(ex1, ex2, tr->partitionData[model].wgt,
								     x1_start, x2_start, tr->partitionData[model].tipVector,
								     tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states,
								     x1_gapColumn, x2_gapColumn, x1_gap, x2_gap); 
		    else
		      partitionLikelihood = evaluateGammaFlex(ex1, ex2, tr->partitionData[model].wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);		   
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, states);
		    
		    partitionLikelihood = evaluateGammaInvarFlex(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								 x1_start, x2_start, 
								 tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								 tr->partitionData[model].propInvariant, 
								 tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);
		  }	  
		  break;
		default:
		  assert(0);
		}	   	      
	    }
	  else
	    {
	      switch(tr->partitionData[model].dataType)
		{ 
		case BINARY_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		      {		   		    
			calcDiagptable(z, BINARY_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood =  evaluateGTRCAT_BINARY(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								     x1_start, x2_start, tr->partitionData[model].tipVector, 
								     tip, width, diagptable, tr->useFastScaling);
		      }
		      break;	  	   
		    case GAMMA:	   
		      {		    		    
			calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
			
			partitionLikelihood = evaluateGTRGAMMA_BINARY(ex1, ex2, tr->partitionData[model].wgt,
								      x1_start, x2_start, tr->partitionData[model].tipVector,
								      tip, width, diagptable, tr->useFastScaling); 		   
		      }
		      break; 
		    case GAMMA_I:
		      {		    		    
			calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMAINVAR_BINARY(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									   x1_start, x2_start,
									   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									   tr->partitionData[model].propInvariant,
									   tip, width, diagptable, tr->useFastScaling);
		      }
		      break;
		    default:
		      assert(0);
		    }
		  break;	   
		case DNA_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:
		     
			  calcDiagptable(z, DNA_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
			  
			  partitionLikelihood =  evaluateGTRCAT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								x1_start, x2_start, tr->partitionData[model].tipVector, 
								tip, width, diagptable, tr->useFastScaling);			  		       
		      break;	  	   
		    case GAMMA:
		     			     	     		      
		      calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
#ifdef __SIM_SSE3
		      if(tr->saveMemory)
			partitionLikelihood = evaluateGTRGAMMA_GAPPED_SAVE(ex1, ex2, tr->partitionData[model].wgt,
									   x1_start, x2_start, tr->partitionData[model].tipVector,
									   tip, width, diagptable, tr->useFastScaling,
									   x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
		      else
#endif
			{
			  if(tr->useGappedImplementation)
			    partitionLikelihood = evaluateGTRGAMMA_GAPPED(ex1, ex2, tr->partitionData[model].wgt,
									  x1_start, x2_start, tr->partitionData[model].tipVector,
									  tip, width, diagptable, tr->useFastScaling,
									  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap); 
			  else
			    partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, tr->partitionData[model].wgt,
								   x1_start, x2_start, tr->partitionData[model].tipVector,
								   tip, width, diagptable, tr->useFastScaling); 		
			}			
		      break; 
		    case GAMMA_I:
		      {
			calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMAINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								    x1_start, x2_start,
								    tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								    tr->partitionData[model].propInvariant,
								    tip, width, diagptable, tr->useFastScaling);
		      }
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case AA_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		     	   
		      calcDiagptable(z, AA_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
			  
		      partitionLikelihood = evaluateGTRCATPROT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							       x1_start, x2_start, tr->partitionData[model].tipVector,
							       tip, width, diagptable, tr->useFastScaling);		  
				     	      
		      break;	      
		    case GAMMA:
		      if(tr->estimatePerSiteAA)
			{ 
			  int 
			    p;
			    
			  for(p = 0; p < (NUM_PROT_MODELS - 3); p++)			    
			    calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->siteProtModel[p].EIGN, tr->siteProtModel[p].left);

			  partitionLikelihood = evaluateGammaFlex_perSite(ex1, 
									  ex2, 
									  tr->partitionData[model].wgt,
									  x1_start, 
									  x2_start,
									  tip, 
									  width,
									  tr->partitionData[model].perSiteAAModel,
									  tr->siteProtModel,
									  tr->useFastScaling, 
									  20);

			}
		      else
			{
			  calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
#ifdef __SIM_SSE3
			  if(tr->saveMemory)
			    partitionLikelihood = evaluateGTRGAMMAPROT_GAPPED_SAVE(ex1, ex2, tr->partitionData[model].wgt,
										   x1_start, x2_start, tr->partitionData[model].tipVector,
										   tip, width, diagptable, tr->useFastScaling,
										   x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
			  else
#endif
			    {
			      if(tr->useGappedImplementation)
				partitionLikelihood = evaluateGTRGAMMAPROT_GAPPED(ex1, ex2, tr->partitionData[model].wgt,
										  x1_start, x2_start, tr->partitionData[model].tipVector,
										  tip, width, diagptable, tr->useFastScaling,
										  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
			      else
				partitionLikelihood = evaluateGTRGAMMAPROT(ex1, ex2, tr->partitionData[model].wgt,
									   x1_start, x2_start, tr->partitionData[model].tipVector,
									   tip, width, diagptable, tr->useFastScaling);
			    }
			}
		      break;
		    case GAMMA_I:		  	    
		      {
			calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMAPROTINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									x1_start, x2_start, 
									tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									tr->partitionData[model].propInvariant, 
									tip, width, diagptable, tr->useFastScaling);
		      }	  
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case GENERIC_32:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		      {
			calcDiagptableFlex(z, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable, states);
			
			partitionLikelihood = evaluateCatFlex(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);
		      }	     	      
		      break;	      
		    case GAMMA:
		      {
			calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, states);
			
			partitionLikelihood = evaluateGammaFlex(ex1, ex2, tr->partitionData[model].wgt,
								x1_start, x2_start, tr->partitionData[model].tipVector,
								tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);
		      }
		      break;
		    case GAMMA_I:		  	    
		      {
			calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, states);
			
			partitionLikelihood = evaluateGammaInvarFlex(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								     x1_start, x2_start, 
								     tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								     tr->partitionData[model].propInvariant, 
								     tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);
		      }	  
		      break;
		    default:
		      assert(0);
		    }
		  break;		  		 
		case SECONDARY_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		      {
			calcDiagptable(z, SECONDARY_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
			
		    partitionLikelihood = evaluateGTRCATSECONDARY(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								  x1_start, x2_start, tr->partitionData[model].tipVector,
								  tip, width, diagptable, tr->useFastScaling);
		      }	     	      
		      break;	      
		    case GAMMA:
		      {
			calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMASECONDARY(ex1, ex2, tr->partitionData[model].wgt,
									x1_start, x2_start, tr->partitionData[model].tipVector,
									tip, width, diagptable, tr->useFastScaling);
		      }
		      break;
		    case GAMMA_I:		  	    
		      {
			calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									     x1_start, x2_start, 
									     tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									     tr->partitionData[model].propInvariant, 
									     tip, width, diagptable, tr->useFastScaling);
		      }	  
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case SECONDARY_DATA_6:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		      {
			calcDiagptable(z, SECONDARY_DATA_6, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRCATSECONDARY_6(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
									x1_start, x2_start, tr->partitionData[model].tipVector,
									tip, width, diagptable, tr->useFastScaling);
		      }	     	      
		      break;	      
		    case GAMMA:
		      {
			calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMASECONDARY_6(ex1, ex2, tr->partitionData[model].wgt,
									  x1_start, x2_start, tr->partitionData[model].tipVector,
									  tip, width, diagptable, tr->useFastScaling);
		      }
		      break;
		    case GAMMA_I:		  	    
		      {
			calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_6(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									       x1_start, x2_start, 
									       tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									       tr->partitionData[model].propInvariant, 
									       tip, width, diagptable, tr->useFastScaling);
		      }	  
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case SECONDARY_DATA_7:
		  switch(tr->rateHetModel)
		    {
		    case CAT:	    
		      {
			calcDiagptable(z, SECONDARY_DATA_7, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRCATSECONDARY_7(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
									x1_start, x2_start, tr->partitionData[model].tipVector,
									tip, width, diagptable, tr->useFastScaling);		   
		      }	     	      
		      break;	      
		    case GAMMA:
		      {
			calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMASECONDARY_7(ex1, ex2, tr->partitionData[model].wgt,
									  x1_start, x2_start, tr->partitionData[model].tipVector,
									  tip, width, diagptable, tr->useFastScaling);
		      }
		      break;
		    case GAMMA_I:		  	    
		      {
			calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			
			partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_7(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									       x1_start, x2_start, 
									       tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									       tr->partitionData[model].propInvariant, 
									       tip, width, diagptable, tr->useFastScaling);
		      }	  
		      break;
		    default:
		      assert(0);
		    }
		  break;
		default:
		  assert(0);
		}
	    }	

	  if(width > 0)
	    {
	      assert(partitionLikelihood < 0.0);
	  
	      if(tr->useFastScaling)
		{	     	      		 
		    {
	      		      
		      partitionLikelihood += (tr->partitionData[model].globalScaler[pNumber] + tr->partitionData[model].globalScaler[qNumber]) * LOG(minlikelihood);
		    }
		}
	    }
	  
	  result += partitionLikelihood;	  
	  tr->perPartitionLH[model] = partitionLikelihood; 
	       
	}
    }
      
  return result;
}


double evaluateIterativeMulti(tree *tr,  boolean writeVector)
{
  double 
    result = 0.0;  
  int pNumber, qNumber, model;
  double *pz; 

  newviewIterativeMulti(tr); 

  if(writeVector)
    assert(!tr->useFastScaling);

  for(model = 0; model < tr->NumberOfModels; model++)
    {            
      if(tr->executeModel[model])
	{		  
	  int 
	    width = tr->partitionData[model].width;
	  
	  double 
	    z, 
	    partitionLikelihood, 
	    *_vector;
	  
	  int    
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL;

	  double 
	    *x1_start   = (double*)NULL, 
	    *x2_start   = (double*)NULL,
	    *diagptable = (double*)NULL;

	  


	  unsigned char 
	    *tip = (unsigned char*)NULL;

	  pNumber = tr->td[model].ti[0].pNumber;
	  qNumber = tr->td[model].ti[0].qNumber;
	  pz      = tr->td[model].ti[0].qz;

	  if(writeVector)
	    _vector = tr->partitionData[model].perSiteLL;
	  else
	    _vector = (double*)NULL;

	 
	  diagptable = tr->partitionData[model].left;


	  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
	    {	        	    
	      if(isTip(qNumber, tr->mxtips))
		{	
		  
		  
		  x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips -1];
		 

		  if(!tr->useFastScaling)
		    ex2      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];
		  
		  tip = tr->partitionData[model].yVector[qNumber];	 	      
		}           
	      else
		{
		  
		  x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
		  

		  if(!tr->useFastScaling)
		    ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];	 
		  
		  tip = tr->partitionData[model].yVector[pNumber];
		}
	    }
	  else
	    {  
	      
	      x1_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
	      x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
	
	      if(!tr->useFastScaling)
		{
		  ex1      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];
		  ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];     
		}
	    }


	  if(tr->multiBranch)
	    z = pz[model];
	  else
	    z = pz[0];

	  switch(tr->partitionData[model].dataType)
	    { 
	    case BINARY_DATA:
	       switch(tr->rateHetModel)
		{
		case CAT:	    
		  {		   		    
		    calcDiagptable(z, BINARY_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood =  evaluateGTRCAT_BINARY(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								 x1_start, x2_start, tr->partitionData[model].tipVector, 
								 tip, width, diagptable, tr->useFastScaling);
		  }
		  break;	  	   
		case GAMMA:	   
		  {		    		    
		    calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    

		    partitionLikelihood = evaluateGTRGAMMA_BINARY(ex1, ex2, tr->partitionData[model].wgt,
								  x1_start, x2_start, tr->partitionData[model].tipVector,
								  tip, width, diagptable, tr->useFastScaling); 		   
		  }
		  break; 
		case GAMMA_I:
		  {		    		    
		    calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMAINVAR_BINARY(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								       x1_start, x2_start,
								       tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								       tr->partitionData[model].propInvariant,
								       tip, width, diagptable, tr->useFastScaling);
		  }
		  break;
		default:
		  assert(0);
		}
	      break;	   
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:		  
		    {
		      calcDiagptable(z, DNA_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
		      
		      partitionLikelihood =  evaluateGTRCAT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							    x1_start, x2_start, tr->partitionData[model].tipVector, 
							    tip, width, diagptable, tr->useFastScaling);
		      		      
		    }
		  break;	  	   
		case GAMMA:		 
		    {		     		      
		      calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
		      if(tr->useGappedImplementation || tr->saveMemory)
			assert(0);
		      else		      
			partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, tr->partitionData[model].wgt,
							       x1_start, x2_start, tr->partitionData[model].tipVector,
							       tip, width, diagptable, tr->useFastScaling); 		      
		    }
		  break; 
		case GAMMA_I:
		  {
		    calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMAINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								x1_start, x2_start,
								tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								tr->partitionData[model].propInvariant,
								tip, width, diagptable, tr->useFastScaling);
		  }
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {		   
		    calcDiagptable(z, AA_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATPROT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
							     x1_start, x2_start, tr->partitionData[model].tipVector,
							     tip, width, diagptable, tr->useFastScaling);		  
		  }	     	      
		  break;	      
		case GAMMA:		  
		    {
		      calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		      
		      partitionLikelihood = evaluateGTRGAMMAPROT(ex1, ex2, tr->partitionData[model].wgt,
								 x1_start, x2_start, tr->partitionData[model].tipVector,
								 tip, width, diagptable, tr->useFastScaling);
		    }
		  break;
		case GAMMA_I:		  	    
		  {
		    calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMAPROTINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
								    x1_start, x2_start, 
								    tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								    tr->partitionData[model].propInvariant, 
								    tip, width, diagptable, tr->useFastScaling);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    calcDiagptable(z, SECONDARY_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATSECONDARY(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								  x1_start, x2_start, tr->partitionData[model].tipVector,
								  tip, width, diagptable, tr->useFastScaling);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMASECONDARY(ex1, ex2, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, tr->useFastScaling);
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									 x1_start, x2_start, 
									 tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									 tr->partitionData[model].propInvariant, 
									 tip, width, diagptable, tr->useFastScaling);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    calcDiagptable(z, SECONDARY_DATA_6, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATSECONDARY_6(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, tr->useFastScaling);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMASECONDARY_6(ex1, ex2, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, tr->useFastScaling);
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_6(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									   x1_start, x2_start, 
									   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									   tr->partitionData[model].propInvariant, 
									   tip, width, diagptable, tr->useFastScaling);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    calcDiagptable(z, SECONDARY_DATA_7, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRCATSECONDARY_7(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								    x1_start, x2_start, tr->partitionData[model].tipVector,
								    tip, width, diagptable, tr->useFastScaling);		   
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

		    partitionLikelihood = evaluateGTRGAMMASECONDARY_7(ex1, ex2, tr->partitionData[model].wgt,
								      x1_start, x2_start, tr->partitionData[model].tipVector,
								      tip, width, diagptable, tr->useFastScaling);
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
		    partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_7(ex1, ex2, tr->partitionData[model].wgt, tr->partitionData[model].invariant,
									   x1_start, x2_start, 
									   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
									   tr->partitionData[model].propInvariant, 
									   tip, width, diagptable, tr->useFastScaling);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
	  
	  if(tr->useFastScaling)
	    {	      
		{

		  partitionLikelihood += (tr->partitionData[model].globalScaler[pNumber] + tr->partitionData[model].globalScaler[qNumber]) * LOG(minlikelihood);
		}
	    }
	  
	  result += partitionLikelihood;	  
	  tr->perPartitionLH[model] = partitionLikelihood;
	}            
    }
      
  return result;
}

double evaluateGeneric (tree *tr, nodeptr p)
{
  volatile double result;
  nodeptr q = p->back; 
  int i;
  
  if(tr->multiGene)
    {     
      nodeptr startNodes[NUM_BRANCHES];  
      nodeptr q;

      findNext(p, tr, startNodes);
      
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  p = startNodes[i];
	  q = p->backs[i];

	  tr->td[i].ti[0].pNumber = p->number;
	  tr->td[i].ti[0].qNumber = q->number;          	  	 
	  tr->td[i].ti[0].qz[i] =  q->z[i];	  
	  tr->td[i].count = 1;

	  if(!p->xs[i])
	    computeTraversalInfoMulti(p, &(tr->td[i].ti[0]), &(tr->td[i].count), tr->mxtips, i);
	  if(!q->xs[i])
	    computeTraversalInfoMulti(q, &(tr->td[i].ti[0]), &(tr->td[i].count), tr->mxtips, i);
	}
      
      result = evaluateIterativeMulti(tr, FALSE);
    }
  else
    {
      tr->td[0].ti[0].pNumber = p->number;
      tr->td[0].ti[0].qNumber = q->number;          
  
      for(i = 0; i < tr->numBranches; i++)    
	tr->td[0].ti[0].qz[i] =  q->z[i];
  
      tr->td[0].count = 1;
      if(!p->x)
	computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
      if(!q->x)
	computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);  
      
#ifdef _USE_PTHREADS 
      {
	int j;
	
	masterBarrier(THREAD_EVALUATE, tr); 
	if(tr->NumberOfModels == 1)
	  {
	    for(i = 0, result = 0.0; i < NumberOfThreads; i++)          
	      result += reductionBuffer[i];  	  	     
	    
	    tr->perPartitionLH[0] = result;
	  }
	else
	  {
	    volatile double partitionResult;
	    
	    result = 0.0;
	    
	    for(j = 0; j < tr->NumberOfModels; j++)
	      {
		for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
		  partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
		result += partitionResult;
		tr->perPartitionLH[j] = partitionResult;
	      }
	  }
      }  
#else
      result = evaluateIterative(tr, FALSE);
#endif   
    }

  tr->likelihood = result;    

  

  return result;
}

double evaluateGenericMulti (tree *tr, nodeptr p, int model)
{
  volatile double result;
  nodeptr q = p->back; 
  
  if(tr->multiGene)
    {               
      int i;
      
      for(i = 0; i < tr->NumberOfModels; i++)
	tr->executeModel[i] = FALSE;
      tr->executeModel[model] = TRUE;
        
      q = p->backs[model];

      assert(q->backs[model] && p->backs[model]);
      assert(q->backs[model] == p);
      assert(p->backs[model] == q);

      tr->td[model].ti[0].pNumber = p->number;
      tr->td[model].ti[0].qNumber = q->number;          	  	 
      tr->td[model].ti[0].qz[model] =  q->z[model];	  
      tr->td[model].count = 1;

      if(!p->xs[model])
	computeTraversalInfoMulti(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->mxtips, model);
      if(!q->xs[model])
	computeTraversalInfoMulti(q, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->mxtips, model);	
      
      result = evaluateIterativeMulti(tr, FALSE);

      for(i = 0; i < tr->NumberOfModels; i++)
	tr->executeModel[i] = TRUE;      
    }
  else
    assert(0);

  return result;
}


double evaluateGenericInitrav (tree *tr, nodeptr p)
{
  volatile double result;   
  
  if(tr->multiGene)
    {
      determineFullTraversalMulti(p, tr);
      result = evaluateIterativeMulti(tr, FALSE);
    }
  else
    {
      determineFullTraversal(p, tr);
      
#ifdef _USE_PTHREADS 
      {
	int i, j;
    
	masterBarrier(THREAD_EVALUATE, tr);    

	if(tr->NumberOfModels == 1)
	  {
	    for(i = 0, result = 0.0; i < NumberOfThreads; i++)          
	      result += reductionBuffer[i];  	  	     
      
	    tr->perPartitionLH[0] = result;
	  }
	else
	  {
	    volatile double partitionResult;
	    
	    result = 0.0;
	    
	    for(j = 0; j < tr->NumberOfModels; j++)
	      {
		for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
		  partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
		result +=  partitionResult;
		tr->perPartitionLH[j] = partitionResult;
	      }
	  }
    
      }
#else
      result = evaluateIterative(tr, FALSE);
#endif

    }
 

  tr->likelihood = result;         

  return result;
}




void onlyInitrav(tree *tr, nodeptr p)
{   
  if(tr->multiGene)
    {
      determineFullTraversalMulti(p, tr);
      newviewIterativeMulti(tr); 
    }
  else
    {
      determineFullTraversal(p, tr);  

#ifdef _USE_PTHREADS  
      masterBarrier(THREAD_NEWVIEW, tr);  	 
#else
      newviewIterative(tr);   
#endif   
    }
}




static void computeFullTraversalInfoMulti(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model)
{
  if(isTip(p->number, maxTips))
    {
      assert(p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);
      return; 
    }

  {           
    if(p->backs[model])
      {
	nodeptr q = p->next->backs[model];
	nodeptr r = p->next->next->backs[model];
	assert(p == p->next->next->next);
	p->xs[model] = 1;
	p->next->xs[model] = 0;
	p->next->next->xs[model] = 0;
	
	if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
	  {
	    assert((r->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]) && (q->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]));
	  
	    ti[*counter].tipCase = TIP_TIP; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;
	    	    
	    {
	      double z;
	      z = q->z[model];
	      z = (z > zmin) ? log(z) : log(zmin);
	      ti[*counter].qz[model] = z;
	      
	      z = r->z[model];
	      z = (z > zmin) ? log(z) : log(zmin);
	      ti[*counter].rz[model] = z;	    
	    }     

	    *counter = *counter + 1;
	  }  
	else
	  {
	    if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
	      {		
		nodeptr tmp;
		 
		if(isTip(r->number, maxTips))
		  {
		    assert(r->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);
		    tmp = r;
		    r = q;
		    q = tmp;
		  }
		else
		  assert(q->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);
		
		computeFullTraversalInfoMulti(r, ti, counter, maxTips, model);	
		
		ti[*counter].tipCase = TIP_INNER; 
		ti[*counter].pNumber = p->number;
		ti[*counter].qNumber = q->number;
		ti[*counter].rNumber = r->number;
		
	      
		{
		  double z;
		  z = q->z[model];
		  z = (z > zmin) ? log(z) : log(zmin);
		  ti[*counter].qz[model] = z;
		  
		  z = r->z[model];
		  z = (z > zmin) ? log(z) : log(zmin);
		  ti[*counter].rz[model] = z;		
		}   
		
		*counter = *counter + 1;
	      }
	    else
	      {	 	  
		computeFullTraversalInfoMulti(q, ti, counter, maxTips, model);	       
		computeFullTraversalInfoMulti(r, ti, counter, maxTips, model);
		
		ti[*counter].tipCase = INNER_INNER; 
		ti[*counter].pNumber = p->number;
		ti[*counter].qNumber = q->number;
		ti[*counter].rNumber = r->number;
	
		{
		  double z;
		  z = q->z[model];
		  z = (z > zmin) ? log(z) : log(zmin);
		  ti[*counter].qz[model] = z;
		  
		  z = r->z[model];
		  z = (z > zmin) ? log(z) : log(zmin);
		  ti[*counter].rz[model] = z;		
		}   
		
		*counter = *counter + 1;
	      }
	  }          
      }
    else
      {	
	p->xs[model] = 0;
	p->next->xs[model] = 0;
	p->next->next->xs[model] = 0;
	assert(p == p->next->next->next);

	computeFullTraversalInfoMulti(p->next->back, ti, counter, maxTips, model);
	computeFullTraversalInfoMulti(p->next->next->back, ti, counter, maxTips, model);
      }
  }
}



void determineFullTraversalMulti(nodeptr p, tree *tr)
{
  int model;

  assert(p == tr->start);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      nodeptr start = tr->startVector[model];
      nodeptr q = start->backs[model];
      
      assert(start->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);
      
      tr->td[model].ti[0].pNumber = start->number;
      tr->td[model].ti[0].qNumber = q->number;
      
      tr->td[model].ti[0].qz[model] = q->z[model];    

      assert(isTip(start->number, tr->mxtips));

      /* entry number zero stores the virtual root */

      tr->td[model].count = 1; 
      computeFullTraversalInfoMulti(q, &(tr->td[model].ti[0]),  &(tr->td[model].count), tr->mxtips, model); 
      computeFullTraversalInfoMulti(start, &(tr->td[model].ti[0]),  &(tr->td[model].count), tr->mxtips, model);

      /*printf("%d %d\n", tr->td[model].count - 1, tr->mxtipsVector[model] - 2);*/
      assert(tr->td[model].count -  1 == tr->mxtipsVector[model] - 2);
    }
}



#ifdef _USE_PTHREADS

double evalCL(tree *tr, double *x2, int *_ex2, unsigned char *_tip, double *pz)
{
  double 
    *x1_start = (double*)NULL,   
    result = 0.0;

  int 
    *ex1 = (int*)NULL,
     model, 
    columnCounter, 
    offsetCounter;

    
    
  unsigned char 
    *tip = (unsigned char*)NULL;

 


  for(model = 0, columnCounter = 0, offsetCounter = 0; model < tr->NumberOfModels; model++)
    {                 	
      int 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower,	
	*ex2,
	*rateCategory, 
	*wgt,         
	*invariant;
		

      double 
	*x2_start,
	z, 
	partitionLikelihood, 	
	*diagptable = tr->partitionData[model].left;	 

      
      rateCategory = &tr->contiguousRateCategory[columnCounter];
      wgt          = &tr->contiguousWgt[columnCounter];
      invariant    = &tr->contiguousInvariant[columnCounter]; 
      tip          = &_tip[columnCounter];
      x2_start     = &x2[offsetCounter];
      ex2          = &_ex2[columnCounter];



      if(tr->multiBranch)
	z = pz[model];
      else
	z = pz[0];

      switch(tr->partitionData[model].dataType)
	{ 
	case BINARY_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    	      	   		    
	      calcDiagptable(z, BINARY_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
		
	      partitionLikelihood =  evaluateGTRCAT_BINARY(ex1, ex2, rateCategory, wgt,
							   x1_start, x2_start, tr->partitionData[model].tipVector, 
							   tip, width, diagptable, tr->useFastScaling);	      	      
	      break;	  	   
	    case GAMMA:	   	  
	      calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
		
	      partitionLikelihood = evaluateGTRGAMMA_BINARY(ex1, ex2,wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, tr->useFastScaling); 
	     
	      break; 
	    case GAMMA_I:	        		    
	      calcDiagptable(z, BINARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		
	      partitionLikelihood = evaluateGTRGAMMAINVAR_BINARY(ex1, ex2,wgt, invariant,
								 x1_start, x2_start,
								 tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								 tr->partitionData[model].propInvariant,
								 tip, width, diagptable, tr->useFastScaling);	      
	      break;
	    default:
	      assert(0);
	    }
	  break;	   
	case DNA_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	     
	      calcDiagptable(z, DNA_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
	      
	      partitionLikelihood =  evaluateGTRCAT(ex1, ex2, rateCategory,wgt,
						    x1_start, x2_start, tr->partitionData[model].tipVector, 
						    tip, width, diagptable, tr->useFastScaling);	      	   	      
	      break;	  	   
	    case GAMMA:		 	     		      
	      calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
		      
	      partitionLikelihood = evaluateGTRGAMMA(ex1, ex2,wgt,
						     x1_start, x2_start, tr->partitionData[model].tipVector,
						     tip, width, diagptable, tr->useFastScaling); 		      	      
	      break; 
	    case GAMMA_I:		  
	      calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

	      partitionLikelihood = evaluateGTRGAMMAINVAR(ex1, ex2,wgt,invariant,
							  x1_start, x2_start,
							  tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
							  tr->partitionData[model].propInvariant,
							  tip, width, diagptable, tr->useFastScaling);	    	     
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case AA_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    		 		   
	      calcDiagptable(z, AA_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
	      
	      partitionLikelihood = evaluateGTRCATPROT(ex1, ex2, rateCategory,wgt,
						       x1_start, x2_start, tr->partitionData[model].tipVector,
						       tip, width, diagptable, tr->useFastScaling);	      
	      break;	      
	    case GAMMA:		 
	      calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		      
	      partitionLikelihood = evaluateGTRGAMMAPROT(ex1, ex2,wgt,
							 x1_start, x2_start, tr->partitionData[model].tipVector,
							 tip, width, diagptable, tr->useFastScaling);	    	      
	      break;
	    case GAMMA_I:		  	    		  
	      calcDiagptable(z, AA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
	      partitionLikelihood = evaluateGTRGAMMAPROTINVAR(ex1, ex2,wgt,invariant,
							      x1_start, x2_start, 
							      tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
							      tr->partitionData[model].propInvariant, 
							      tip, width, diagptable, tr->useFastScaling);	  	      
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    	      
	      calcDiagptable(z, SECONDARY_DATA, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);

	      partitionLikelihood = evaluateGTRCATSECONDARY(ex1, ex2, rateCategory,wgt,
							    x1_start, x2_start, tr->partitionData[model].tipVector,
							    tip, width, diagptable, tr->useFastScaling);		 	      
	      break;	      
	    case GAMMA:		  
	      calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

	      partitionLikelihood = evaluateGTRGAMMASECONDARY(ex1, ex2,wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, tr->useFastScaling);		  	     
	      break;
	    case GAMMA_I:		  	    		  
	      calcDiagptable(z, SECONDARY_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
		    
	      partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR(ex1, ex2,wgt,invariant,
								   x1_start, x2_start, 
								   tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								   tr->partitionData[model].propInvariant, 
								   tip, width, diagptable, tr->useFastScaling);			     
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA_6:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    		 
	      calcDiagptable(z, SECONDARY_DATA_6, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);

	      partitionLikelihood = evaluateGTRCATSECONDARY_6(ex1, ex2, rateCategory,wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, tr->useFastScaling);		  	     	      	      
	      break;	      
	    case GAMMA:		  
	      calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

	      partitionLikelihood = evaluateGTRGAMMASECONDARY_6(ex1, ex2,wgt,
								x1_start, x2_start, tr->partitionData[model].tipVector,
								tip, width, diagptable, tr->useFastScaling);		  	      
	      break;
	    case GAMMA_I:		  	    		  
	      calcDiagptable(z, SECONDARY_DATA_6, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
	      
	      partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_6(ex1, ex2,wgt,invariant,
								     x1_start, x2_start, 
								     tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								     tr->partitionData[model].propInvariant, 
								     tip, width, diagptable, tr->useFastScaling);			      
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA_7:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	    		  
	      calcDiagptable(z, SECONDARY_DATA_7, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable);
	      
	      partitionLikelihood = evaluateGTRCATSECONDARY_7(ex1, ex2, rateCategory,wgt,
							      x1_start, x2_start, tr->partitionData[model].tipVector,
							      tip, width, diagptable, tr->useFastScaling);	      
	      break;	      
	    case GAMMA:	      
	      calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);

	      partitionLikelihood = evaluateGTRGAMMASECONDARY_7(ex1, ex2,wgt,
								x1_start, x2_start, tr->partitionData[model].tipVector,
								tip, width, diagptable, tr->useFastScaling);	      
	      break;
	    case GAMMA_I:		  	    	      
	      calcDiagptable(z, SECONDARY_DATA_7, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
	      
	      partitionLikelihood = evaluateGTRGAMMASECONDARYINVAR_7(ex1, ex2,wgt,invariant,
								     x1_start, x2_start, 
								     tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								     tr->partitionData[model].propInvariant, 
								     tip, width, diagptable, tr->useFastScaling);			     
	      break;
	    default:
	      assert(0);
	    }
	  break;
	default:
	  assert(0);
	}

      assert(!tr->useFastScaling);
      	  
      result += partitionLikelihood;	

      columnCounter += width;
      offsetCounter += width * tr->partitionData[model].states * tr->discreteRateCategories;	
    }           
      

  return result;
}







#endif




/*****************************************************************************************************/



double evaluateGenericVector (tree *tr, nodeptr p)
{
  volatile double result;
  nodeptr q = p->back; 
  int i;
  
  if(tr->multiGene)
    {   
      assert(0);
      
      /*
	nodeptr startNodes[NUM_BRANCHES];  
	nodeptr q;
	
	findNext(p, tr, startNodes);
	
	for(i = 0; i < tr->NumberOfModels; i++)
	{
	p = startNodes[i];
	q = p->backs[i];

	tr->td[i].ti[0].pNumber = p->number;
	tr->td[i].ti[0].qNumber = q->number;          	  	 
	tr->td[i].ti[0].qz[i] =  q->z[i];	  
	tr->td[i].count = 1;
	
	if(!p->xs[i])
	computeTraversalInfoMulti(p, &(tr->td[i].ti[0]), &(tr->td[i].count), tr->mxtips, i);
	if(!q->xs[i])
	    computeTraversalInfoMulti(q, &(tr->td[i].ti[0]), &(tr->td[i].count), tr->mxtips, i);
	    }
      
	    result = evaluateIterativeMulti(tr, FALSE);
      */
    }
  else
    {
      tr->td[0].ti[0].pNumber = p->number;
      tr->td[0].ti[0].qNumber = q->number;          
  
      for(i = 0; i < tr->numBranches; i++)    
	tr->td[0].ti[0].qz[i] =  q->z[i];
  
      tr->td[0].count = 1;
      if(!p->x)
	computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
      if(!q->x)
	computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);  
      
#ifdef _USE_PTHREADS 
      {
	int j;

	masterBarrier(THREAD_EVALUATE_VECTOR, tr);
	if(tr->NumberOfModels == 1)
	  {
	    for(i = 0, result = 0.0; i < NumberOfThreads; i++)	      	       
	      result += reductionBuffer[i];  	  	     	     
	    
	    tr->perPartitionLH[0] = result;
	  }
	else
	  {
	    volatile double partitionResult;
	    
	    result = 0.0;
	    
	    for(j = 0; j < tr->NumberOfModels; j++)
	      {
		for(i = 0, partitionResult = 0.0; i < NumberOfThreads; i++)          	      
		  partitionResult += reductionBuffer[i * tr->NumberOfModels + j];
		result += partitionResult;
		tr->perPartitionLH[j] = partitionResult;
	      }
	  }
      }  
#else
      result = evaluateIterative(tr, TRUE);
#endif   
    }

  tr->likelihood = result;    
  
  return result;
}
