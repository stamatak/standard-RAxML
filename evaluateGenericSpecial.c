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




void ascertainmentBiasSequence(unsigned char tip[32], int numStates, int dataType, int nodeNumber)
{ 
  assert(numStates <= 32 && numStates > 1);

  assert(nodeNumber > 0);
  
  switch(dataType)
    {
    case BINARY_DATA:          
      tip[0] = 1;
      tip[1] = 2;	
      break;
    case DNA_DATA:      
      tip[0] = 1;
      tip[1] = 2;
      tip[2] = 4;
      tip[3] = 8;	
      break;
    case AA_DATA:     
      {
	int 
	  i;
	
	for(i = 0; i < numStates; i++)	  
	  tip[i] = i;	  
      }    
      break;
    case GENERIC_32:      	
      {
	int 
	  i;
	
	for(i = 0; i < numStates; i++)	  
	  tip[i] = i;
      }
      break;
    default:
      assert(0);
    }
}

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


static void calcDiagptableFlex_LG4(double z, int numberOfCategories, double *rptr, double *EIGN[4], double *diagptable, const int numStates)
{
  int 
    i, 
    l;
  
  double 
    lz;
  
  assert(numStates <= 64);
  
  if (z < zmin) 
    lz = log(zmin);
  else
    lz = log(z);

  for(i = 0; i <  numberOfCategories; i++)
    {	      	       
      diagptable[i * numStates] = 1.0;

      for(l = 1; l < numStates; l++)
	diagptable[i * numStates + l] = EXP(rptr[i] * EIGN[i][l - 1] * lz);     	          
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

//asc

static double evaluateCatAsc(int *ex1, int *ex2,
			     double *x1, double *x2,  
			     double *tipVector, 
			     unsigned char *tipX1, int n, double *diagptable, const int numStates, 
			     double *accumulator, double *weightVector, int dataType, int nodeNumber, 
			     double *goldmanAccumulator)
{
  double
    exponent,
    sum = 0.0, 
    unobserved,
    term,
    *left, 
    *right;
  
  int     
    i,    
    l;   
         
 
   
  if(tipX1)
    {            
      unsigned char 
	tip[32];

      ascertainmentBiasSequence(tip, numStates, dataType, nodeNumber);  
      
      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[numStates * tip[i]]);	  	  
	  right = &(x2[i * numStates]);

	  term = 0.0;
	         	      
	  for(l = 0; l < numStates; l++)
	    term += left[l] * right[l] * diagptable[l];	      	 	 	  	 

	  /* assumes that pow behaves as expected/specified for underflows
	     from the man page: 
	     If  result  underflows, and is not representable, a range error occurs,
	     and 0.0 is returned.
	  */


	  exponent = pow(minlikelihood, (double)ex2[i]);

	  unobserved = FABS(term) * exponent;	    
	  
#ifdef _DEBUG_ASC
	  if(ex2[i] > 0)
	    {
	      printf("s %d\n", ex2[i]);
	      assert(0);
	    }
#endif	  
	    
	  if(weightVector)
	    *accumulator += weightVector[i] * (LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood)));

	  *goldmanAccumulator += ((LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood))) * FABS(term) * exponent);
	  
	  sum += unobserved;
	}              
    }      
  else
    {           
      for (i = 0; i < n; i++) 
	{	  	 
	  term = 0.0;
	  
	  left  = &(x1[i * numStates]);
	  right = &(x2[i * numStates]);	    
	  
	  for(l = 0; l < numStates; l++)
	    term += left[l] * right[l] * diagptable[l];		  
	  
	  /* assumes that pow behaves as expected/specified for underflows
	     from the man page: 
	     If  result  underflows, and is not representable, a range error occurs,
	     and 0.0 is returned.
	  */

	  exponent = pow(minlikelihood, (double)(ex1[i] + ex2[i]));
	  
	  unobserved = FABS(term) * exponent;
	  
#ifdef _DEBUG_ASC
	  if(ex2[i] > 0 || ex1[i] > 0)
	    {
	      printf("s %d %d\n", ex1[i], ex2[i]);
	      assert(0);
	    }
#endif
	  
	  if(weightVector)
	    *accumulator += weightVector[i] * (LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood)));

	  *goldmanAccumulator += ((LOG(FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood))) * FABS(term) * exponent);
	  
	  sum += unobserved;
	}          
    }
  
  return  sum;
}


static double evaluateGammaAsc(int *ex1, int *ex2,
			       double *x1, double *x2,  
			       double *tipVector, 
			       unsigned char *tipX1, int n, double *diagptable, const int numStates, 
			       double *accumulator, double *weightVector, int dataType, int nodeNumber,
			       double *goldmanAccumulator)
{
  double
    exponent,
    sum = 0.0, 
    unobserved,
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
      unsigned char 
	tip[32];

      ascertainmentBiasSequence(tip, numStates, dataType, nodeNumber);            

      for (i = 0; i < n; i++) 
	{
	  left = &(tipVector[numStates * tip[i]]);	  	  
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      right = &(x2[gammaStates * i + numStates * j]);
	      
	      for(l = 0; l < numStates; l++)
		term += left[l] * right[l] * diagptable[j * numStates + l];	      
	    }	 	  	 


	  /* assumes that pow behaves as expected/specified for underflows
	     from the man page: 
	     If  result  underflows, and is not representable, a range error occurs,
	     and 0.0 is returned.
	  */

	  exponent = pow(minlikelihood, (double)ex2[i]);	  

	  unobserved = 0.25 * FABS(term) * exponent;	    
	  
#ifdef _DEBUG_ASC	 	  
	  if(ex2[i] > 0)
	    {
	      printf("s %d\n", ex2[i]);
	      assert(0);
	    }
#endif	  
	  if(weightVector)
	    *accumulator += weightVector[i] * (LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood)));
	  

	  *goldmanAccumulator += ((LOG(0.25 * FABS(term)) + (ex2[i] * LOG(minlikelihood))) * 0.25 * FABS(term) * exponent);

	  sum += unobserved;
	}              
    }              
  else
    {           
      for (i = 0; i < n; i++) 
	{	  	 	             
	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      left  = &(x1[gammaStates * i + numStates * j]);
	      right = &(x2[gammaStates * i + numStates * j]);	    
	      
	      for(l = 0; l < numStates; l++)
		term += left[l] * right[l] * diagptable[j * numStates + l];	
	    }
	  
	  /* assumes that pow behaves as expected/specified for underflows
	     from the man page: 
	     If  result  underflows, and is not representable, a range error occurs,
	     and 0.0 is returned.
	  */

	  exponent = pow(minlikelihood, (double)(ex1[i] + ex2[i]));
	  
	  unobserved = 0.25 * FABS(term) * exponent;
	  
#ifdef _DEBUG_ASC
	  if(ex2[i] > 0 || ex1[i] > 0)
	    {
	      printf("s %d %d\n", ex1[i], ex2[i]);
	      assert(0);
	    }
#endif
	 	  	  	  
	   if(weightVector)
	     *accumulator += weightVector[i] * (LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood)));

	   *goldmanAccumulator += ((LOG(0.25 * FABS(term)) + ((ex1[i] + ex2[i]) * LOG(minlikelihood))) * 0.25 * FABS(term) * exponent);
	   
	  sum += unobserved;
	}             
    }        

  return  sum;
}

//asc

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

static double evaluateGammaFlex_LG4(int *ex1, int *ex2, int *wptr,
				    double *x1, double *x2,  
				    double *tipVector[4], 
				    unsigned char *tipX1, int n, double *diagptable, double *vector, boolean writeVector, const boolean fastScaling, const int numStates, double *weights)
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
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		double 
		  t = 0.0;

		left = &(tipVector[j][numStates * tipX1[i]]);
		right = &(x2[gammaStates * i + numStates * j]);
		
		for(l = 0; l < numStates; l++)
		  t += left[l] * right[l] * diagptable[j * numStates + l];	      

		term += weights[j] * t;
	      }	 
	    
	    if(fastScaling)
	      term = LOG(FABS(term));
	    else
	      term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	    	    
	    vector[i] = term;
	    
	    sum += wptr[i] * term;
	  }         
      else
	{       
	  for (i = 0; i < n; i++) 
	    {	     	      	  	  	      
	      for(j = 0, term = 0.0; j < 4; j++)
		{
		  double
		    t = 0.0;
		  
		  left = &(tipVector[j][numStates * tipX1[i]]);
		  right = &(x2[gammaStates * i + numStates * j]);
		  
		  for(l = 0; l < numStates; l++)
		    t += left[l] * right[l] * diagptable[j * numStates + l];	      

		  term += weights[j] * t;
		}
	      
	      if(fastScaling)
		term = LOG(FABS(term));
	      else
		term = LOG(FABS(term)) + (ex2[i] * LOG(minlikelihood));	   
	      
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
	      double 
		t = 0.0;
	      
	      left  = &(x1[gammaStates * i + numStates * j]);
	      right = &(x2[gammaStates * i + numStates * j]);	    
	      
	      for(l = 0; l < numStates; l++)
		t += left[l] * right[l] * diagptable[j * numStates + l];	

	      term +=  weights[j] * t;
	    }
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	
	  vector[i] = term;
  
	  sum += wptr[i] * term;
	}         
      else
	for (i = 0; i < n; i++) 
	  {	  	 	             
	    
	    for(j = 0, term = 0.0; j < 4; j++)
	      {
		double 
		  t = 0.0;
		
		left  = &(x1[gammaStates * i + numStates * j]);
		right = &(x2[gammaStates * i + numStates * j]);	    
		
		for(l = 0; l < numStates; l++)
		  t += left[l] * right[l] * diagptable[j * numStates + l];	

		term +=  weights[j] * t;
	      }
	    
	    if(fastScaling)
	      term = LOG(FABS(term));
	    else
	      term = LOG(FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	    
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


#ifdef __SIM_SSE3




static double evaluateGTRCATPROT_SAVE (int *ex1, int *ex2, int *cptr, int *wptr,
    double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling,
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double   
    sum = 0.0, 
  term,
  *diagptable,  
  *left, 
  *right,
  *left_ptr = x1,
  *right_ptr = x2;

  int     
    i, 
    l;                           

  if(tipX1)
  {                 
    for (i = 0; i < n; i++) 
    {	       	
      left = &(tipVector[20 * tipX1[i]]);

      if(isGap(x2_gap, i))
        right = x2_gapColumn;
      else
      {
        right = right_ptr;
        right_ptr += 20;
      }	  	 

      diagptable = &diagptable_start[20 * cptr[i]];	           	 

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

      if(fastScaling)
	term = LOG(term);
      else
	term = LOG(term) + (ex2[i] * LOG(minlikelihood));

      sum += wptr[i] * term;
    }      
  }    
  else
  {

    for (i = 0; i < n; i++) 
    {		       	      	      	  
      if(isGap(x1_gap, i))
        left = x1_gapColumn;
      else
      {
        left = left_ptr;
        left_ptr += 20;
      }

      if(isGap(x2_gap, i))
        right = x2_gapColumn;
      else
      {
        right = right_ptr;
        right_ptr += 20;
      }

      diagptable = &diagptable_start[20 * cptr[i]];	  	

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
      
      if(fastScaling)
	term = LOG(term);	 
      else 
	term = LOG(term) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));

      sum += wptr[i] * term;      
    }
  }

  return  sum;         
} 


static double evaluateGTRCAT_SAVE (int *ex1, int *ex2, int *cptr, int *wptr,
    double *x1_start, double *x2_start, double *tipVector, 		      
    unsigned char *tipX1, int n, double *diagptable_start, const boolean fastScaling,
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double  sum = 0.0, term;       
  int     i;

  double  *diagptable, 
          *x1, 
          *x2,
          *x1_ptr = x1_start,
          *x2_ptr = x2_start;

  if(tipX1)
  {           
    for (i = 0; i < n; i++) 
    {	
      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
      __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

      x1 = &(tipVector[4 * tipX1[i]]);

      if(isGap(x2_gap, i))
        x2 = x2_gapColumn;
      else
      {
        x2 = x2_ptr;
        x2_ptr += 4;
      }

      diagptable = &diagptable_start[4 * cptr[i]];

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
	term = LOG(t[0] + t[1]);      
      else
	term =  LOG(t[0] + t[1]) + (ex2[i] * LOG(minlikelihood));

      sum += wptr[i] * term;
    }	
  }               
  else
  {
    for (i = 0; i < n; i++) 
    { 
      double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
      __m128d x1v1, x1v2, x2v1, x2v2, dv1, dv2;

      if(isGap(x1_gap, i))
        x1 = x1_gapColumn;
      else
      {
        x1 = x1_ptr;
        x1_ptr += 4;
      }

      if(isGap(x2_gap, i))
        x2 = x2_gapColumn;
      else
      {
        x2 = x2_ptr;
        x2_ptr += 4;
      }

      diagptable = &diagptable_start[4 * cptr[i]];	

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
	term = LOG(t[0] + t[1]);
      else
	term = LOG(t[0] + t[1]) + ((ex1[i] + ex2[i]) * LOG(minlikelihood));
      
      sum += wptr[i] * term;
    }    
  }

  return  sum;         
} 

#endif




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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));	    		   	  
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));	    		   	  
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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

	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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
	  double t[2] __attribute__ ((aligned (BYTE_ALIGNMENT)));
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


static double evaluateGTRGAMMAPROT_LG4(int *ex1, int *ex2, int *wptr,
				       double *x1, double *x2,  
				       double *tipVector[4], 
				       unsigned char *tipX1, int n, double *diagptable, const boolean fastScaling, double *weights)
{
  double   sum = 0.0, term;        
  int     i, j, l;   
  double  *left, *right;              
  
  if(tipX1)
    {               
      for (i = 0; i < n; i++) 
	{	  	  	  
	  for(j = 0, term = 0.0; j < 4; j++)
	    {	   
	      double 
		t = 0.0;

	      left = &(tipVector[j][20 * tipX1[i]]);
	      right = &(x2[80 * i + 20 * j]);
	      
	      for(l = 0; l < 20; l++)	 
		t += left[l] * right[l] * diagptable[j * 20 + l];

	      term += weights[j] * t;	      
	    }	  
	  
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
	  for(j = 0, term = 0.0; j < 4; j++)
	    {
	      double 
		t = 0.0;

	      left  = &(x1[80 * i + 20 * j]);
	      right = &(x2[80 * i + 20 * j]);	    
	      
	      for(l = 0; l < 20; l++)
		t += left[l] * right[l] * diagptable[j * 20 + l];	

	      term += weights[j] * t;
	    }
	  
	  if(fastScaling)
	    term = LOG(FABS(term));
	  else
	    term = LOG(FABS(term)) + ((ex1[i] + ex2[i])*LOG(minlikelihood));
	  
	  sum += wptr[i] * term;
	}         
    }
       
  return  sum;
}




#ifdef __SIM_SSE3

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
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x1v,
    *x2v;              

  if(tipX1)
  {               
    for (i = 0; i < n; i++) 
    {
      if(x2_gap[i / 32] & mask32[i % 32])
        x2v = x2_gapColumn;
      else
      {
        x2v = x2_ptr;
        x2_ptr += 80;
      }

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


      if(fastScaling)
        term = LOG(0.25 * term);
      else
        term = LOG(0.25 * term) + (ex2[i] * LOG(minlikelihood));	   

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
      {
        x1v = x1_ptr;
        x1_ptr += 80;
      }

      if(x2_gap[i / 32] & mask32[i % 32])
        x2v = x2_gapColumn;
      else
      {
        x2v = x2_ptr;
        x2_ptr += 80;
      }

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

      if(fastScaling)
        term = LOG(0.25 * term);
      else
        term = LOG(0.25 * term) + ((ex1[i] + ex2[i])*LOG(minlikelihood));

      sum += wptr[i] * term;
    }         
  }

  return  sum;
}


#endif






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
    rateHet,
    pNumber = tr->td[0].ti[0].pNumber, 
    qNumber = tr->td[0].ti[0].qNumber, 
    model;

  if(tr->rateHetModel == CAT)
    rateHet = 1;
  else
    rateHet = 4;
  
  newviewIterative(tr); 

  if(writeVector)
    assert(!tr->useFastScaling);

#ifdef _DEBUG_MULTI_EPA  
  printf("EV: ");
#endif

  for(model = 0; model < tr->NumberOfModels; model++)
    {         
#ifdef _DEBUG_MULTI_EPA  
      printf("%d ", tr->executeModel[model]);
#endif
        
      if(tr->executeModel[model])
	{	
	  int 
	    width = tr->partitionData[model].width,
	    states = tr->partitionData[model].states,
	    ascWidth = tr->partitionData[model].states;
	  
	  double 
	    z, 
	    partitionLikelihood = 0.0, 
	    *_vector;
	  
	  int
	    tipNodeNumber = -1,
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL,
	    *ex1_asc = (int*)NULL, 
	    *ex2_asc = (int*)NULL;
	  
	  unsigned int
	    *x1_gap = (unsigned int*)NULL,
	    *x2_gap = (unsigned int*)NULL;

	  double 
	    *weights    = tr->partitionData[model].weights,
	    *x1_start   = (double*)NULL, 
	    *x2_start   = (double*)NULL,
	    *x1_start_asc   = (double*)NULL, 
	    *x2_start_asc   = (double*)NULL,
	    *diagptable = (double*)NULL,
	    *x1_gapColumn = (double*)NULL,
	    *x2_gapColumn = (double*)NULL;
	 
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
		  tipNodeNumber = qNumber;
			  		  
		  x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips -1];
		  
		  if(!tr->useFastScaling)
		    ex2      = tr->partitionData[model].expVector[pNumber - tr->mxtips - 1];

		  if(tr->saveMemory)
		    {
		      x2_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
		      x2_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];
		    }

		  tip = tr->partitionData[model].yVector[qNumber];
		  
#ifdef _USE_PTHREADS
		  if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
		    if(tr->partitionData[model].ascBias)
#endif		  
		    { 
		      x2_start_asc = &tr->partitionData[model].ascVector[(pNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];		     
		      ex2_asc = &tr->partitionData[model].ascExpVector[(pNumber - tr->mxtips - 1) * ascWidth];
		    }    
		}           
	      else
		{
		  tipNodeNumber = pNumber;
		  
		  x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
		  

		  if(!tr->useFastScaling)
		    ex2      = tr->partitionData[model].expVector[qNumber - tr->mxtips - 1];

		  if(tr->saveMemory)
		    {
		      x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
		      x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
		    }
		  
		  tip = tr->partitionData[model].yVector[pNumber];

#ifdef _USE_PTHREADS
		  if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
		    if(tr->partitionData[model].ascBias)
#endif		 
		    {		      
		      x2_start_asc = &tr->partitionData[model].ascVector[(qNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];		     
		      ex2_asc = &tr->partitionData[model].ascExpVector[(qNumber - tr->mxtips - 1) * ascWidth];
		    }		  
		}
	    }
	  else
	    {  
	      if(tr->saveMemory)
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

#ifdef _USE_PTHREADS
	      if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
		if(tr->partitionData[model].ascBias)
#endif	      
		 {
		   x1_start_asc = &tr->partitionData[model].ascVector[(pNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];
		   x2_start_asc = &tr->partitionData[model].ascVector[(qNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];

		   ex1_asc = &tr->partitionData[model].ascExpVector[(pNumber - tr->mxtips - 1) * ascWidth];
		   ex2_asc = &tr->partitionData[model].ascExpVector[(qNumber - tr->mxtips - 1) * ascWidth];
		 }
	    }


	  if(tr->multiBranch)
	    z = pz[model];
	  else
	    z = pz[0];


#ifdef _BASTIEN
	  if(tr->doBastienStuff)
	    {
	      if(tr->multiBranch)
		assert(0);
	      else
		{
		  assert(tr->td[0].ti[0].secondDerivativeP[0] == tr->td[0].ti[0].secondDerivativeQ[0]);
		  printf("\nHello, I am the second derivative at the root branch %1.40f\n", tr->td[0].ti[0].secondDerivativeQ[0]);
		}
		  
	    }
#endif

	  //printf("branch %f\n", z);

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
		    if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
			{		   			 			 
			  calcDiagptableFlex_LG4(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN_LG4, diagptable, 20);

			 
			  partitionLikelihood = evaluateGammaFlex_LG4(ex1, ex2, tr->partitionData[model].wgt,
								      x1_start, x2_start, tr->partitionData[model].tipVector_LG4,
								      tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states, weights);			  
			}
		    else
		      {
			calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, states);
		     	
			partitionLikelihood = evaluateGammaFlex(ex1, ex2, tr->partitionData[model].wgt,
								x1_start, x2_start, tr->partitionData[model].tipVector,
								tip, width, diagptable, _vector, writeVector, tr->useFastScaling, states);		   
		      }
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
#ifdef __SIM_SSE3
			  if(tr->saveMemory)
			    {			     			      
			      partitionLikelihood = evaluateGTRCAT_SAVE(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
									    x1_start, x2_start, tr->partitionData[model].tipVector,
									    tip, width, diagptable, tr->useFastScaling,  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
			    }
			  else
#endif
			    partitionLikelihood =  evaluateGTRCAT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								  x1_start, x2_start, tr->partitionData[model].tipVector, 
								  tip, width, diagptable, tr->useFastScaling);			  		       
		      break;	  	   
		    case GAMMA:
#ifdef _HET		      
		      assert(!tr->saveMemory);

		      if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
			{			  
			  calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN_TIP, diagptable);
			  
			  partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, tr->partitionData[model].wgt,
								 x1_start, x2_start, tr->partitionData[model].tipVector_TIP,
								 tip, width, diagptable, tr->useFastScaling);			
			}
		      else
			{			  
			  calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);
			  
			  partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, tr->partitionData[model].wgt,
								 x1_start, x2_start, (double*)NULL,
								 tip, width, diagptable, tr->useFastScaling);			
			}
#else
		     			     	     		      
		      calcDiagptable(z, DNA_DATA, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable);		    		    
#ifdef __SIM_SSE3
		      if(tr->saveMemory)						  			  
			partitionLikelihood = evaluateGTRGAMMA_GAPPED_SAVE(ex1, ex2, tr->partitionData[model].wgt,
									   x1_start, x2_start, tr->partitionData[model].tipVector,
									   tip, width, diagptable, tr->useFastScaling,
									   x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);			
		      else
#endif
		
			partitionLikelihood = evaluateGTRGAMMA(ex1, ex2, tr->partitionData[model].wgt,
							       x1_start, x2_start, tr->partitionData[model].tipVector,
							       tip, width, diagptable, tr->useFastScaling); 					
#endif
		      break; 
		    case GAMMA_I:
		      { 
			assert(!tr->saveMemory);
			
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
#ifdef __SIM_SSE3
		      if(tr->saveMemory)
			{			 
			  partitionLikelihood = evaluateGTRCATPROT_SAVE(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
									x1_start, x2_start, tr->partitionData[model].tipVector,
									tip, width, diagptable, tr->useFastScaling,  x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
			}
		      else
#endif			  
			partitionLikelihood = evaluateGTRCATPROT(ex1, ex2, tr->partitionData[model].rateCategory, tr->partitionData[model].wgt,
								 x1_start, x2_start, tr->partitionData[model].tipVector,
								 tip, width, diagptable, tr->useFastScaling);		  
				     	      
		      break;	      
		    case GAMMA:	
		      if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
			{						  
			  calcDiagptableFlex_LG4(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN_LG4, diagptable, 20);
			  			  
			  partitionLikelihood = evaluateGTRGAMMAPROT_LG4(ex1, ex2, tr->partitionData[model].wgt,
									 x1_start, x2_start, tr->partitionData[model].tipVector_LG4,
									 tip, width, diagptable, tr->useFastScaling, weights);			 			    
			  
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
			    partitionLikelihood = evaluateGTRGAMMAPROT(ex1, ex2, tr->partitionData[model].wgt,
								       x1_start, x2_start, tr->partitionData[model].tipVector,
								       tip, width, diagptable, tr->useFastScaling);			
			}
		      break;
		    case GAMMA_I:		  	    
		      {
			assert(!tr->saveMemory);

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
	      if(tr->useFastScaling)		    	      		      
		partitionLikelihood += (tr->partitionData[model].globalScaler[pNumber] + tr->partitionData[model].globalScaler[qNumber]) * LOG(minlikelihood);		    
	      
	      assert(partitionLikelihood < 0.0);
#ifdef _USE_PTHREADS
	      if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
		if(tr->partitionData[model].ascBias)
#endif	  	      
		  {
		    size_t
		      i;
		    
		    int	       
		      w = 0;
		    
		    double 	
		      *weightVector = (double*)NULL,
		      accumulator = 0.0,
		      goldmanAccumulator = 0.0,
		      correction;

		    if(tr->ascertainmentCorrectionType == STAMATAKIS_CORRECTION)		     
		      weightVector = tr->partitionData[model].invariableFrequencies;		      
		    //		      invariableWeight;

		    for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
		      w += tr->cdta->aliaswgt[i];	

		    switch(tr->rateHetModel)
		      {
		      case CAT:
			{
			  double 
			    rates = 1.0;
			  
			  //need to re-calculate P-matrix for the correction here assuming a rate of 1.0 
			  calcDiagptableFlex(z, 1, &rates, tr->partitionData[model].EIGN, diagptable, states);
			  
			  
			  correction = evaluateCatAsc(ex1_asc, ex2_asc, x1_start_asc, x2_start_asc, tr->partitionData[model].tipVector,
						      tip, ascWidth, diagptable, ascWidth, &accumulator, weightVector, tr->partitionData[model].dataType, 
						      tipNodeNumber, &goldmanAccumulator);     		  	 	       
			}
			break;
		      case GAMMA:			
			correction = evaluateGammaAsc(ex1_asc, ex2_asc, x1_start_asc, x2_start_asc, tr->partitionData[model].tipVector,
						      tip, ascWidth, diagptable, ascWidth, &accumulator, weightVector, tr->partitionData[model].dataType,
						      tipNodeNumber, &goldmanAccumulator);			
			break;
		      default:
			assert(0);
		      }
		    
		   
		    switch(tr->ascertainmentCorrectionType)
		      {
		      case LEWIS_CORRECTION:		    			  		  	      	     	     		   
			partitionLikelihood = partitionLikelihood - (double)w * LOG(1.0 - correction);			
			break;
		      case STAMATAKIS_CORRECTION:		       
			partitionLikelihood += accumulator;
			break;
		      case FELSENSTEIN_CORRECTION:
			partitionLikelihood += tr->partitionData[model].invariableWeight * LOG(correction);		       
			break;
		      case GOLDMAN_CORRECTION_1:			
			partitionLikelihood += ((correction * (double)w * LOG(correction)) / (1.0 - correction));					       	  
			break;
		      case GOLDMAN_CORRECTION_2:
			//printf("Goldman acc: %f\n", goldmanAccumulator);
			partitionLikelihood += (((double)w / (1.0 - correction)) * goldmanAccumulator);
			break;
		      case GOLDMAN_CORRECTION_3:
			partitionLikelihood += ((tr->partitionData[model].invariableWeight / correction) * goldmanAccumulator);
			break;
		      default:
			assert(0);
		      }
	      
#ifdef _DEBUG_ASC 	      
		    printf("E w: %f %f ARG %f ragu %f\n", partitionLikelihood, (double)w, 1.0 - correction, (((double)w / (1.0 - correction)) * goldmanAccumulator));	
		    
#endif	      		    
		  }
	      
#ifdef _USE_PTHREADS
	      if(!(tr->partitionData[model].ascBias && tr->threadID == 0))
		{
#endif
		  if(partitionLikelihood >= 0.0)
		    {
		      printf("positive log like: %f for partition %d\n", partitionLikelihood, model);
		      assert(0);
		    }
#ifdef _USE_PTHREADS	      
		}
#endif
	    }
	  result += partitionLikelihood;	  
	 
	  tr->perPartitionLH[model] = partitionLikelihood; 	       
	}
    }
#ifdef _DEBUG_MULTI_EPA
  printf("\n");
#endif 
  return result;
}




double evaluateGeneric (tree *tr, nodeptr p)
{
  volatile 
    double result;
  
  nodeptr 
    q = p->back; 
  
  int 
    i;
  
  
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;          
  
  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];
  
  tr->td[0].count = 1;
  
  if(!p->x)
    computeTraversalInfo(tr, p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
  
  if(!q->x)
    computeTraversalInfo(tr, q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);  
  
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
	volatile 
	  double partitionResult;
	
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
 
  assert(result <= 0.0);

  tr->likelihood = result;      

  return result;
}




double evaluateGenericInitrav (tree *tr, nodeptr p)
{
  volatile double 
    result;   
  
  determineFullTraversal(p, tr);
    
#ifdef _USE_PTHREADS 
  {
    int 
      i, 
      j;
      
    masterBarrier(THREAD_EVALUATE, tr);    
      
    if(tr->NumberOfModels == 1)
      {
	for(i = 0, result = 0.0; i < NumberOfThreads; i++)          
	  result += reductionBuffer[i];  	  	     
	  
	tr->perPartitionLH[0] = result;
      }
    else
      {
	volatile double 
	  partitionResult;
	  
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

  
  assert(result <= 0.0);

  tr->likelihood = result;         

  return result;
}




void onlyInitrav(tree *tr, nodeptr p)
{     
  determineFullTraversal(p, tr);  

#ifdef _USE_PTHREADS  
  masterBarrier(THREAD_NEWVIEW, tr);  	 
#else
  newviewIterative(tr);   
#endif      
}






#ifdef _USE_PTHREADS

double evalCL(tree *tr, double *x2, int *_ex2, unsigned char *_tip, double *pz, int insertion)
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

  setPartitionMask(tr, insertion, tr->executeModel);

#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("EV %s: ", tr->nameList[tr->inserts[insertion]]);
#endif

  for(model = 0, columnCounter = 0, offsetCounter = 0; model < tr->NumberOfModels; model++)
    {   
      int 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;

#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("%d", tr->executeModel[model]);
#endif    

      if(tr->executeModel[model])
	{
	  int 	 
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
	    case GENERIC_32:
	      switch(tr->rateHetModel)
		{
		case CAT:	    
		  {
		    calcDiagptableFlex(z, tr->partitionData[model].numberOfCategories, tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN, diagptable, tr->partitionData[model].states);
		    
		    partitionLikelihood = evaluateCatFlex(ex1, ex2, rateCategory, wgt,
							  x1_start, x2_start, tr->partitionData[model].tipVector,
							  tip, width, diagptable, (double*)NULL, FALSE, tr->useFastScaling, tr->partitionData[model].states);
		  }	     	      
		  break;	      
		case GAMMA:
		  {
		    calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, tr->partitionData[model].states);
		    
		    partitionLikelihood = evaluateGammaFlex(ex1, ex2, wgt,
							    x1_start, x2_start, tr->partitionData[model].tipVector,
							    tip, width, diagptable, (double*)NULL, FALSE, tr->useFastScaling, tr->partitionData[model].states);
		  }
		  break;
		case GAMMA_I:		  	    
		  {
		    calcDiagptableFlex(z, 4, tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, diagptable, tr->partitionData[model].states);
		    
		    partitionLikelihood = evaluateGammaInvarFlex(ex1, ex2, wgt, invariant,
								 x1_start, x2_start, 
								 tr->partitionData[model].tipVector, tr->partitionData[model].frequencies, 
								 tr->partitionData[model].propInvariant, 
								 tip, width, diagptable, (double*)NULL, FALSE, tr->useFastScaling, tr->partitionData[model].states);
		  }	  
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
	  
	  assert(!tr->useFastScaling);

	  /* error ? */
	  
      	  tr->perPartitionLH[model] = partitionLikelihood; 
	  
	  result += partitionLikelihood;	
	}
      
      columnCounter += width;
      offsetCounter += width * tr->partitionData[model].states * tr->discreteRateCategories;       
    }         
  
  resetPartitionMask(tr, tr->executeModel);
#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("\n");
#endif
  if(tr->perPartitionEPA)
    {
      assert(tr->perPartitionLH[tr->readPartition[insertion]] <= 0.0);
      return (tr->perPartitionLH[tr->readPartition[insertion]]);
    }
  else
    { 
      assert(result <= 0.0);
      return result;      
    }
}







#endif




/*****************************************************************************************************/



double evaluateGenericVector (tree *tr, nodeptr p)
{
  volatile double result;
  nodeptr q = p->back; 
  int i;
  
  
  {
    tr->td[0].ti[0].pNumber = p->number;
    tr->td[0].ti[0].qNumber = q->number;          
    
    for(i = 0; i < tr->numBranches; i++)    
      tr->td[0].ti[0].qz[i] =  q->z[i];
    
    tr->td[0].count = 1;
    if(!p->x)
      computeTraversalInfo(tr, p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
    if(!q->x)
      computeTraversalInfo(tr, q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);  
    
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

  assert(result <= 0.0);

  tr->likelihood = result;    
  
  return result;
}
