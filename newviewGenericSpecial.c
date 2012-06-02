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
#include <stdint.h>
#include <limits.h>
#include "axml.h"

#ifdef __SIM_SSE3

#include <stdint.h>
#include <xmmintrin.h>
#include <pmmintrin.h>

const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
       uint64_t i[2];
       __m128d m;
} absMask = {{0x7fffffffffffffffULL , 0x7fffffffffffffffULL }};




#endif


#ifdef _USE_PTHREADS
#include <pthread.h>
extern volatile int NumberOfThreads;
extern pthread_mutex_t          mutex;
#endif


extern const unsigned int mask32[32];

static inline boolean isGap(unsigned int *x, int pos)
{
  return (x[pos / 32] & mask32[pos % 32]);
}

static inline boolean noGap(unsigned int *x, int pos)
{
  return (!(x[pos / 32] & mask32[pos % 32]));
}

void makeP_Flex(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, const int numStates)
{
  int 
    i,
    j,
    k;
  
  const int
    rates = numStates - 1,
    statesSquare = numStates * numStates;

  double 
    lz1[64],
    lz2[64], 
    d1[64],  
    d2[64];

  assert(numStates <= 64);
     
  for(i = 0; i < rates; i++)
    {
      lz1[i] = EIGN[i] * z1;
      lz2[i] = EIGN[i] * z2;
    }

  for(i = 0; i < numberOfCategories; i++)
    {
      for(j = 0; j < rates; j++)
	{
	  d1[j] = EXP (rptr[i] * lz1[j]);
	  d2[j] = EXP (rptr[i] * lz2[j]);
	}

      for(j = 0; j < numStates; j++)
	{
	  left[statesSquare * i  + numStates * j] = 1.0;
	  right[statesSquare * i + numStates * j] = 1.0;

	  for(k = 1; k < numStates; k++)
	    {
	      left[statesSquare * i + numStates * j + k]  = d1[k-1] * EI[rates * j + (k-1)];
	      right[statesSquare * i + numStates * j + k] = d2[k-1] * EI[rates * j + (k-1)];
	    }
	}
    }  
}






static void newviewFlexCat(int tipCase, double *extEV,
			   int *cptr,
			   double *x1, double *x2, double *x3, double *tipVector,
			   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling, const int numStates)
{
  double
    *le, *ri, *v, *vl, *vr,
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale, addScale = 0;

  const int 
    statesSquare = numStates * numStates;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * statesSquare];
	    ri = &right[cptr[i] * statesSquare];

	    vl = &(tipVector[numStates * tipX1[i]]);
	    vr = &(tipVector[numStates * tipX2[i]]);
	    v  = &x3[numStates * i];

	    for(l = 0; l < numStates; l++)
	      v[l] = 0.0;

	    for(l = 0; l < numStates; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < numStates; j++)
		  {
		    ump_x1 += vl[j] * le[l * numStates + j];
		    ump_x2 += vr[j] * ri[l * numStates + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < numStates; j++)
		  v[j] += x1px2 * extEV[l * numStates + j];
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * statesSquare];
	    ri = &right[cptr[i] * statesSquare];

	    vl = &(tipVector[numStates * tipX1[i]]);
	    vr = &x2[numStates * i];
	    v  = &x3[numStates * i];

	    for(l = 0; l < numStates; l++)
	      v[l] = 0.0;

	    for(l = 0; l < numStates; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < numStates; j++)
		  {
		    ump_x1 += vl[j] * le[l * numStates + j];
		    ump_x2 += vr[j] * ri[l * numStates + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < numStates; j++)
		  v[j] += x1px2 * extEV[l * numStates + j];
	      }

	    scale = 1;
	    for(l = 0; scale && (l < numStates); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	    

	    if(scale)
	      {
		for(l = 0; l < numStates; l++)
		  v[l] *= twotothe256;
		
		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * statesSquare];
	  ri = &right[cptr[i] * statesSquare];

	  vl = &x1[numStates * i];
	  vr = &x2[numStates * i];
	  v = &x3[numStates * i];

	  for(l = 0; l < numStates; l++)
	    v[l] = 0.0;

	  for(l = 0; l < numStates; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < numStates; j++)
		{
		  ump_x1 += vl[j] * le[l * numStates + j];
		  ump_x2 += vr[j] * ri[l * numStates + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < numStates; j++)
		v[j] += x1px2 * extEV[l * numStates + j];
	    }

	   scale = 1;
	   for(l = 0; scale && (l < numStates); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	  
	   if(scale)
	     {
	       for(l = 0; l < numStates; l++)
		 v[l] *= twotothe256;

	       if(useFastScaling)
		 addScale += wgt[i];
	       else
		 ex3[i]  += 1;	     
	     }
	}
      break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewFlexGamma(int tipCase,
			     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			     int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			     int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling, 
			     const int numStates)
{
  double  *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr, al, ar;

  const int 
    statesSquare = numStates * numStates,
    gammaStates = 4 * numStates;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[numStates * tipX1[i]]);
		vr = &(tipVector[numStates * tipX2[i]]);
		v =  &(x3[gammaStates * i + numStates * k]);

		for(l = 0; l < numStates; l++)
		  v[l] = 0;

		for(l = 0; l < numStates; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < numStates; j++)
		      {
			al += vl[j] * left[k * statesSquare + l * numStates + j];
			ar += vr[j] * right[k * statesSquare + l * numStates + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < numStates; j++)
		      v[j] += x1px2 * extEV[numStates * l + j];
		  }
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[numStates * tipX1[i]]);
		vr = &(x2[gammaStates * i + numStates * k]);
		v =  &(x3[gammaStates * i + numStates * k]);

		for(l = 0; l < numStates; l++)
		  v[l] = 0;

		for(l = 0; l < numStates; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < numStates; j++)
		      {
			al += vl[j] * left[k * statesSquare + l * numStates + j];
			ar += vr[j] * right[k * statesSquare + l * numStates + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < numStates; j++)
		      v[j] += x1px2 * extEV[numStates * l + j];
		  }
	      }
	   
	    v = &x3[gammaStates * i];
	    scale = 1;
	    for(l = 0; scale && (l < gammaStates); l++)
	      scale = (ABS(v[l]) <  minlikelihood);

	    if(scale)
	      {
		for(l = 0; l < gammaStates; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[gammaStates * i + numStates * k]);
	     vr = &(x2[gammaStates * i + numStates * k]);
	     v =  &(x3[gammaStates * i + numStates * k]);

	     for(l = 0; l < numStates; l++)
	       v[l] = 0;

	     for(l = 0; l < numStates; l++)
	       {
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < numStates; j++)
		   {
		     al += vl[j] * left[k * statesSquare + l * numStates + j];
		     ar += vr[j] * right[k * statesSquare + l * numStates + j];
		   }

		 x1px2 = al * ar;
		 for(j = 0; j < numStates; j++)
		   v[j] += x1px2 * extEV[numStates * l + j];
	       }
	   }
	 
	 v = &(x3[gammaStates * i]);
	 scale = 1;
	 for(l = 0; scale && (l < gammaStates); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if (scale)
	   {
	     for(l = 0; l < gammaStates; l++)
	       v[l] *= twotothe256;

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	    
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}








static void newviewFlexGamma_perSite(int tipCase,
				     double *x1, double *x2, double *x3, 
				     int *perSiteAA, 
				     siteAAModels *siteProtModel,
				     int *ex3, 
				     unsigned char *tipX1, 
				     unsigned char *tipX2,
				     int n, 				     
				     int *wgt, 
				     int *scalerIncrement, 
				     const boolean useFastScaling, 
				     const int numStates)
{
  double  *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr, al, ar;

  const int     
    statesSquare = numStates * numStates,  
    gammaStates  = 4 * numStates;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    double 
	      *tipVector = siteProtModel[perSiteAA[i]].tipVector,
	      *extEV     = siteProtModel[perSiteAA[i]].EV,
	      *left      = siteProtModel[perSiteAA[i]].left,
	      *right     = siteProtModel[perSiteAA[i]].right;

	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[numStates * tipX1[i]]);
		vr = &(tipVector[numStates * tipX2[i]]);
		v =  &(x3[gammaStates * i + numStates * k]);

		for(l = 0; l < numStates; l++)
		  v[l] = 0;

		for(l = 0; l < numStates; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < numStates; j++)
		      {
			al += vl[j] * left[k * statesSquare + l * numStates + j];
			ar += vr[j] * right[k * statesSquare + l * numStates + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < numStates; j++)
		      v[j] += x1px2 * extEV[numStates * l + j];
		  }
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    double 
	      *tipVector = siteProtModel[perSiteAA[i]].tipVector,
	      *extEV     = siteProtModel[perSiteAA[i]].EV,
	      *left      = siteProtModel[perSiteAA[i]].left,
	      *right     = siteProtModel[perSiteAA[i]].right;

	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[numStates * tipX1[i]]);
		vr = &(x2[gammaStates * i + numStates * k]);
		v =  &(x3[gammaStates * i + numStates * k]);

		for(l = 0; l < numStates; l++)
		  v[l] = 0;

		for(l = 0; l < numStates; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < numStates; j++)
		      {
			al += vl[j] * left[k * statesSquare + l * numStates + j];
			ar += vr[j] * right[k * statesSquare + l * numStates + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < numStates; j++)
		      v[j] += x1px2 * extEV[numStates * l + j];
		  }
	      }
	   
	    v = &x3[gammaStates * i];
	    scale = 1;
	    for(l = 0; scale && (l < gammaStates); l++)
	      scale = (ABS(v[l]) <  minlikelihood);

	    if(scale)
	      {
		for(l = 0; l < gammaStates; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 double 	   
	   *extEV     = siteProtModel[perSiteAA[i]].EV,
	   *left      = siteProtModel[perSiteAA[i]].left,
	   *right     = siteProtModel[perSiteAA[i]].right;

	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[gammaStates * i + numStates * k]);
	     vr = &(x2[gammaStates * i + numStates * k]);
	     v =  &(x3[gammaStates * i + numStates * k]);

	     for(l = 0; l < numStates; l++)
	       v[l] = 0;

	     for(l = 0; l < numStates; l++)
	       {
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < numStates; j++)
		   {
		     al += vl[j] * left[k * statesSquare + l * numStates + j];
		     ar += vr[j] * right[k * statesSquare + l * numStates + j];
		   }

		 x1px2 = al * ar;
		 for(j = 0; j < numStates; j++)
		   v[j] += x1px2 * extEV[numStates * l + j];
	       }
	   }
	 
	 v = &(x3[gammaStates * i]);
	 scale = 1;
	 for(l = 0; scale && (l < gammaStates); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if (scale)
	   {
	     for(l = 0; l < gammaStates; l++)
	       v[l] *= twotothe256;

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	    
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


static void makeP(double z1, double z2, double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, double *right, int data, boolean saveMem, int maxCat)
{
  int i, j, k;

  switch(data)
    {
    case BINARY_DATA:
      {
	double d1, d2;

	for(i = 0; i < numberOfCategories; i++)
	  {
	    d1 = EXP(rptr[i] * EIGN[0] * z1);
	    d2 = EXP(rptr[i] * EIGN[0] * z2);
	    
	    for(j = 0; j < 2; j++)
	      {
		left[i * 4 + j * 2] = 1.0;
		right[i * 4 + j * 2] = 1.0;

		left[i * 4 + j * 2 + 1]  = d1 * EI[j];
		right[i * 4 + j * 2 + 1] = d2 * EI[j];	
	      }
	  }
      }
      break;
    case DNA_DATA:
      {
#ifdef __SIM_SSE3
	double 
	  d1[4] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
	  d2[4] __attribute__ ((aligned (BYTE_ALIGNMENT))),
	  ez1[3], 
	  ez2[3],
	  EI_16[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	
	  	  
	for(j = 0; j < 4; j++)
	  {
	    EI_16[j * 4] = 1.0;
	    for(k = 0; k < 3; k++)
	      EI_16[j * 4 + k + 1] = EI[3 * j + k];
	  }	  

	for(j = 0; j < 3; j++)
	  {
	    ez1[j] = EIGN[j] * z1;
	    ez2[j] = EIGN[j] * z2;
	  }


	for(i = 0; i < numberOfCategories; i++)
	  {	   
	    __m128d 
	      d1_0, d1_1,
	      d2_0, d2_1;
 
	    d1[0] = 1.0;
	    d2[0] = 1.0;

	    for(j = 0; j < 3; j++)
	      {
		d1[j+1] = EXP(rptr[i] * ez1[j]);
		d2[j+1] = EXP(rptr[i] * ez2[j]);
	      }

	    d1_0 = _mm_load_pd(&d1[0]);
	    d1_1 = _mm_load_pd(&d1[2]);

	    d2_0 = _mm_load_pd(&d2[0]);
	    d2_1 = _mm_load_pd(&d2[2]);
	    

	    for(j = 0; j < 4; j++)
	      {	       
		double *ll = &left[i * 16 + j * 4];
		double *rr = &right[i * 16 + j * 4];	       

		__m128d eev = _mm_load_pd(&EI_16[4 * j]);
		
		_mm_store_pd(&ll[0], _mm_mul_pd(d1_0, eev));
		_mm_store_pd(&rr[0], _mm_mul_pd(d2_0, eev));
		
		eev = _mm_load_pd(&EI_16[4 * j + 2]);
		
		_mm_store_pd(&ll[2], _mm_mul_pd(d1_1, eev));
		_mm_store_pd(&rr[2], _mm_mul_pd(d2_1, eev));

		
	      }
	  }

	if(saveMem)
	  {
	    i = maxCat;
	    	    
	    {	   
	      __m128d 
		d1_0, d1_1,
		d2_0, d2_1;
	      
	      d1[0] = 1.0;
	      d2[0] = 1.0;
	      
	      for(j = 0; j < 3; j++)
		{
		  d1[j+1] = EXP(ez1[j]);
		  d2[j+1] = EXP(ez2[j]);
		}	     
	      
	      d1_0 = _mm_load_pd(&d1[0]);
	      d1_1 = _mm_load_pd(&d1[2]);
	      
	      d2_0 = _mm_load_pd(&d2[0]);
	      d2_1 = _mm_load_pd(&d2[2]);
	      	      
	      for(j = 0; j < 4; j++)
		{	       
		  double *ll = &left[i * 16 + j * 4];
		  double *rr = &right[i * 16 + j * 4];	       
		  
		  __m128d eev = _mm_load_pd(&EI_16[4 * j]);
		  
		  _mm_store_pd(&ll[0], _mm_mul_pd(d1_0, eev));
		  _mm_store_pd(&rr[0], _mm_mul_pd(d2_0, eev));
		  
		  eev = _mm_load_pd(&EI_16[4 * j + 2]);
		  
		  _mm_store_pd(&ll[2], _mm_mul_pd(d1_1, eev));
		  _mm_store_pd(&rr[2], _mm_mul_pd(d2_1, eev));
		  
		  
		}
	    }
	  }
      
#else
	double d1[3], d2[3];

	for(i = 0; i < numberOfCategories; i++)
	  {
	    for(j = 0; j < 3; j++)
	      {
		d1[j] = EXP(rptr[i] * EIGN[j] * z1);
		d2[j] = EXP(rptr[i] * EIGN[j] * z2);
	      }

	    for(j = 0; j < 4; j++)
	      {
		left[i * 16 + j * 4] = 1.0;
		right[i * 16 + j * 4] = 1.0;

		for(k = 0; k < 3; k++)
		  {
		    left[i * 16 + j * 4 + k + 1]  = d1[k] * EI[3 * j + k];
		    right[i * 16 + j * 4 + k + 1] = d2[k] * EI[3 * j + k];
		  }
	      }
	  }
#endif
      }
      break;
    case SECONDARY_DATA:
      {
	double lz1[15], lz2[15], d1[15], d2[15];

	for(i = 0; i < 15; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;
	  }

	for(i = 0; i < numberOfCategories; i++)
	  {
	    for(j = 0; j < 15; j++)
	      {
		d1[j] = EXP (rptr[i] * lz1[j]);
		d2[j] = EXP (rptr[i] * lz2[j]);
	      }

	    for(j = 0; j < 16; j++)
	      {
		left[256 * i  + 16 * j] = 1.0;
		right[256 * i + 16 * j] = 1.0;

		for(k = 1; k < 16; k++)
		  {
		    left[256 * i + 16 * j + k]  = d1[k-1] * EI[15 * j + (k-1)];
		    right[256 * i + 16 * j + k] = d2[k-1] * EI[15 * j + (k-1)];
		  }
	      }
	  }
      }
      break;
    case SECONDARY_DATA_6:
      {
	double lz1[5], lz2[5], d1[5], d2[5];

	for(i = 0; i < 5; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;
	  }

	for(i = 0; i < numberOfCategories; i++)
	  {
	    for(j = 0; j < 5; j++)
	      {
		d1[j] = EXP (rptr[i] * lz1[j]);
		d2[j] = EXP (rptr[i] * lz2[j]);
	      }

	    for(j = 0; j < 6; j++)
	      {
		left[36 * i  + 6 * j] = 1.0;
		right[36 * i + 6 * j] = 1.0;

		for(k = 1; k < 6; k++)
		  {
		    left[36 * i + 6 * j + k]  = d1[k-1] * EI[5 * j + (k-1)];
		    right[36 * i + 6 * j + k] = d2[k-1] * EI[5 * j + (k-1)];
		  }
	      }
	  }
      }
      break;
    case SECONDARY_DATA_7:
      {
	double lz1[6], lz2[6], d1[6], d2[6];

	for(i = 0; i < 6; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;
	  }

	for(i = 0; i < numberOfCategories; i++)
	  {
	    for(j = 0; j < 6; j++)
	      {
		d1[j] = EXP (rptr[i] * lz1[j]);
		d2[j] = EXP (rptr[i] * lz2[j]);
	      }

	    for(j = 0; j < 7; j++)
	      {
		left[49 * i  + 7 * j] = 1.0;
		right[49 * i + 7 * j] = 1.0;

		for(k = 1; k < 7; k++)
		  {
		    left[49 * i + 7 * j + k]  = d1[k-1] * EI[6 * j + (k-1)];
		    right[49 * i + 7 * j + k] = d2[k-1] * EI[6 * j + (k-1)];
		  }
	      }
	  }
      }
      break;
    case AA_DATA:
      {
	double lz1[19], lz2[19], d1[19], d2[19];

	for(i = 0; i < 19; i++)
	  {
	    lz1[i] = EIGN[i] * z1;
	    lz2[i] = EIGN[i] * z2;
	  }

	for(i = 0; i < numberOfCategories; i++)
	  {
	    for(j = 0; j < 19; j++)
	      {
		d1[j] = EXP (rptr[i] * lz1[j]);
		d2[j] = EXP (rptr[i] * lz2[j]);
	      }

	    for(j = 0; j < 20; j++)
	      {
		left[400 * i  + 20 * j] = 1.0;
		right[400 * i + 20 * j] = 1.0;

		for(k = 1; k < 20; k++)
		  {
		    left[400 * i + 20 * j + k]  = d1[k-1] * EI[19 * j + (k-1)];
		    right[400 * i + 20 * j + k] = d2[k-1] * EI[19 * j + (k-1)];
		  }
	      }
	  }

	if(saveMem)
        {
          i = maxCat;

          for(j = 0; j < 19; j++)
          {
            d1[j] = EXP (lz1[j]);
            d2[j] = EXP (lz2[j]);
          }

          for(j = 0; j < 20; j++)
          {
            left[400 * i  + 20 * j] = 1.0;
            right[400 * i + 20 * j] = 1.0;

            for(k = 1; k < 20; k++)
            {
              left[400 * i + 20 * j + k]  = d1[k-1] * EI[19 * j + (k-1)];
              right[400 * i + 20 * j + k] = d2[k-1] * EI[19 * j + (k-1)];
            }
          }
        }      
      }
      break;
    default:
      assert(0);
    }

}




#ifndef __SIM_SSE3

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
				  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
				  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  double
    ump_x1, ump_x2, x1px2[2];
  int i, j, k, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &(tipVector[2 * tipX2[i]]);
	    x3 = &x3_start[2 * i];	    

	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];

	    for(j = 0; j < 2; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 2; k++)
		  {
		    ump_x1 += x1[k] * le[j * 2 + k];
		    ump_x2 += x2[k] * ri[j * 2 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }

	    for(j = 0; j < 2; j++)
	      x3[j] = 0.0;

	    for(j = 0; j < 2; j++)
	      for(k = 0; k < 2; k++)
		x3[k] += x1px2[j] * EV[j * 2 + k];	   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[2 * i];
	    x3 = &x3_start[2 * i];
	    
	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];

	    for(j = 0; j < 2; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 2; k++)
		  {
		    ump_x1 += x1[k] * le[j * 2 + k];
		    ump_x2 += x2[k] * ri[j * 2 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }

	    for(j = 0; j < 2; j++)
	      x3[j] = 0.0;

	    for(j = 0; j < 2; j++)
	      for(k = 0; k < 2; k++)
		x3[k] +=  x1px2[j] *  EV[2 * j + k];	   

	    scale = 1;
	    for(j = 0; j < 2 && scale; j++)
	      scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);

	    if(scale)
	      {
		for(j = 0; j < 2; j++)
		  x3[j] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	       
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  x3 = &x3_start[2 * i];

	  le = &left[cptr[i] * 4];
	  ri = &right[cptr[i] * 4];

	  for(j = 0; j < 2; j++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      for(k = 0; k < 2; k++)
		{
		  ump_x1 += x1[k] * le[j * 2 + k];
		  ump_x2 += x2[k] * ri[j * 2 + k];
		}
	      x1px2[j] = ump_x1 * ump_x2;
	    }

	  for(j = 0; j < 2; j++)
	    x3[j] = 0.0;

	  for(j = 0; j < 2; j++)
	    for(k = 0; k < 2; k++)
	      x3[k] +=  x1px2[j] *  EV[2 * j + k];	  

	  scale = 1;
	  for(j = 0; j < 2 && scale; j++)
	    scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);

	  if(scale)
	    {
	      for(j = 0; j < 2; j++)
		x3[j] *= twotothe256;

	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	   
	    }
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#else

static void newviewGTRCAT_BINARY( int tipCase,  double *EV,  int *cptr,
				  double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
				  int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				  int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  int i, l, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &(tipVector[2 * tipX2[i]]);
	    x3 = &x3_start[2 * i];	   

	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];

	    _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	    for(l = 0; l < 2; l++)
	      {		 		 						   		  		 		 
		__m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
		__m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
		
		al = _mm_hadd_pd(al, al);
		ar = _mm_hadd_pd(ar, ar);
		
		al = _mm_mul_pd(al, ar);
		
		__m128d vv  = _mm_load_pd(x3);
		__m128d EVV = _mm_load_pd(&EV[2 * l]);
		
		vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		
		_mm_store_pd(x3, vv);		     	  		   		  
	      }	    	   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &x2_start[2 * i];
	    x3 = &x3_start[2 * i];
	    
	    le =  &left[cptr[i] * 4];
	    ri =  &right[cptr[i] * 4];

	    _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	    for(l = 0; l < 2; l++)
	      {		 		 						   		  		 		 
		__m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
		__m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
		
		al = _mm_hadd_pd(al, al);
		ar = _mm_hadd_pd(ar, ar);
		
		al = _mm_mul_pd(al, ar);
		
		__m128d vv  = _mm_load_pd(x3);
		__m128d EVV = _mm_load_pd(&EV[2 * l]);
		
		vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		
		_mm_store_pd(x3, vv);		     	  		   		  
	      }	 
	    
	    __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	    scale = 1;
	    
	    __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
	    v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	    if(_mm_movemask_pd( v1 ) != 3)
	      scale = 0;	  	         
	    
	    if(scale)
	      {
		__m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		
		__m128d ex3v = _mm_load_pd(x3);		  
		_mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));		    		   		  
		
		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	  
	      }	 	 	   
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
	  x3 = &x3_start[2 * i];

	  le = &left[cptr[i] * 4];
	  ri = &right[cptr[i] * 4];

	  _mm_store_pd(x3, _mm_setzero_pd());	    
	  
	  for(l = 0; l < 2; l++)
	    {		 		 						   		  		 		 
	      __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&le[l * 2]));
	      __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&ri[l * 2]));
	      
	      al = _mm_hadd_pd(al, al);
	      ar = _mm_hadd_pd(ar, ar);
	      
	      al = _mm_mul_pd(al, ar);
	      
	      __m128d vv  = _mm_load_pd(x3);
	      __m128d EVV = _mm_load_pd(&EV[2 * l]);
	      
	      vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
	      
	      _mm_store_pd(x3, vv);		     	  		   		  
	    }	 	 	 	  

	  __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	  scale = 1;
	  	  
	  __m128d v1 = _mm_and_pd(_mm_load_pd(x3), absMask.m);
	  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	  if(_mm_movemask_pd( v1 ) != 3)
	    scale = 0;	  	         
	 
	  if(scale)
	    {
	      __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	      	    
	      __m128d ex3v = _mm_load_pd(x3);		  
	      _mm_store_pd(x3, _mm_mul_pd(ex3v,twoto));		    		   		  
	     
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	  
	   }	 	 
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


#endif

#ifdef __SIM_SSE3
static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
				   )
{
  double
    *x1, *x2, *x3;
 
  int i, k, l, scale, addScale = 0; 

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 x2  = &(tipVector[2 * tipX2[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     	    
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
       }
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
       {
	 x1  = &(tipVector[2 * tipX1[i]]);
	 
	 for(k = 0; k < 4; k++)
	   {	     	     
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }      
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {	 
	 for(k = 0; k < 4; k++)
	   {	     
	     x1 = &(x1_start[8 * i + 2 * k]);
	     x2 = &(x2_start[8 * i + 2 * k]);
	     x3 = &(x3_start[8 * i + 2 * k]);	     
	    	         
	     _mm_store_pd(x3, _mm_setzero_pd());	    
	    	     
	     for(l = 0; l < 2; l++)
	       {		 		 						   		  		 		 
		 __m128d al = _mm_mul_pd(_mm_load_pd(x1), _mm_load_pd(&left[k * 4 + l * 2]));
		 __m128d ar = _mm_mul_pd(_mm_load_pd(x2), _mm_load_pd(&right[k * 4 + l * 2]));
		 		       
		 al = _mm_hadd_pd(al, al);
		 ar = _mm_hadd_pd(ar, ar);
		   
		 al = _mm_mul_pd(al, ar);
		   
		 __m128d vv  = _mm_load_pd(x3);
		 __m128d EVV = _mm_load_pd(&EV[2 * l]);
		 
		 vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
		 
		 _mm_store_pd(x3, vv);		     	  		   		  
	       }	     	    
	   }
	
	 x3 = &(x3_start[8 * i]);
	 __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l += 2)
	   {
	     __m128d vv = _mm_load_pd(&x3[l]);
	     __m128d v1 = _mm_and_pd(vv, absMask.m);
	     v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	     if(_mm_movemask_pd( v1 ) != 3)
	       scale = 0;
	   }	    	         
	 
	 if(scale)
	   {
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 8; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&x3[l]);		  
		 _mm_store_pd(&x3[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }	 
       }
      break;

    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#else

static void newviewGTRGAMMA_BINARY(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
				   )
{
  double
    *x1, *x2, *x3;
  double
    ump_x1,
    ump_x2,
    x1px2[4];
  int i, j, k, l, scale, addScale = 0;


  /* C-OPT figure out if we are at an inner node who has two tips/leaves
     as descendants TIP_TIP, a tip and another inner node as descendant
     TIP_INNER, or two inner nodes as descendants INNER_INNER */

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[2 * tipX1[i]]);
	    x2 = &(tipVector[2 * tipX2[i]]);
	    x3 = &x3_start[i * 8];

	    for(j = 0; j < 8; j++)
	      x3[j] = 0.0;

	    for (j = 0; j < 4; j++)
	      {
		for (k = 0; k < 2; k++)
		  {
		    ump_x1 = 0.0;
		    ump_x2 = 0.0;

		    for (l=0; l < 2; l++)
		      {
			ump_x1 += x1[l] * left[ j*4 + k*2 + l];
			ump_x2 += x2[l] * right[j*4 + k*2 + l];
		      }

		    x1px2[k] = ump_x1 * ump_x2;
		  }

		for(k = 0; k < 2; k++)
		  for (l = 0; l < 2; l++)
		    x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	 for (i = 0; i < n; i++)
	   {
	     x1 = &(tipVector[2 * tipX1[i]]);
	     x2 = &x2_start[i * 8];
	     x3 = &x3_start[i * 8];

	     for(j = 0; j < 8; j++)
	       x3[j] = 0.0;

	     for (j = 0; j < 4; j++)
	       {
		 for (k = 0; k < 2; k++)
		   {
		     ump_x1 = 0.0;
		     ump_x2 = 0.0;

		     for (l=0; l < 2; l++)
		       {
			 ump_x1 += x1[l] * left[ j*4 + k*2 + l];
			 ump_x2 += x2[j*2 + l] * right[j*4 + k*2 + l];
		       }

		     x1px2[k] = ump_x1 * ump_x2;
		   }

		 for(k = 0; k < 2; k++)
		   for (l = 0; l < 2; l++)
		     x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

	       }	    

	     scale = 1;
	     for(l = 0; scale && (l < 8); l++)
	       scale = (ABS(x3[l]) <  minlikelihood);

	     if(scale)
	       {
		 for (l=0; l < 8; l++)
		   x3[l] *= twotothe256;
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;	       
	       }

	   }
      }
      break;
    case INNER_INNER:

      /* C-OPT here we don't do any pre-computations
	 This should be the most compute intensive loop of the three
	 cases here. If we have one or two tips as descendants
	 we can take a couple of shortcuts */


     for (i = 0; i < n; i++)
       {
	 x1 = &x1_start[i * 8];
	 x2 = &x2_start[i * 8];
	 x3 = &x3_start[i * 8];

	 for(j = 0; j < 8; j++)
	   x3[j] = 0.0;

	 for (j = 0; j < 4; j++)
	   {
	     for (k = 0; k < 2; k++)
	       {
		 ump_x1 = 0.0;
		 ump_x2 = 0.0;

		 for (l=0; l < 2; l++)
		   {
		     ump_x1 += x1[j*2 + l] * left[ j*4 + k*2 + l];
		     ump_x2 += x2[j*2 + l] * right[j*4 + k*2 + l];
		   }

		 x1px2[k] = ump_x1 * ump_x2;
	       }

	     for(k = 0; k < 2; k++)
	       for (l = 0; l < 2; l++)
		 x3[j * 2 + l] +=  x1px2[k] * EV[2 * k + l];

	   }
	 
	 scale = 1;
	 for(l = 0; scale && (l < 8); l++)
	   scale = (ABS(x3[l]) <  minlikelihood);


	 if(scale)
	   {
	     for (l=0; l<8; l++)
	       x3[l] *= twotothe256;

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }
       }
     break;

    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#endif



#ifndef __SIM_SSE3

static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
			   double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
			   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1, *x2, *x3;
  double
    ump_x1, ump_x2, x1px2[4];
  int i, j, k, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &(tipVector[4 * tipX2[i]]);
	    x3 = &x3_start[4 * i];

	    le =  &left[cptr[i] * 16];
	    ri =  &right[cptr[i] * 16];

	    for(j = 0; j < 4; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 4; k++)
		  {
		    ump_x1 += x1[k] * le[j * 4 + k];
		    ump_x2 += x2[k] * ri[j * 4 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }

	    for(j = 0; j < 4; j++)
	      x3[j] = 0.0;

	    for(j = 0; j < 4; j++)
	      for(k = 0; k < 4; k++)
		x3[k] += x1px2[j] * EV[j * 4 + k];	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    x1 = &(tipVector[4 * tipX1[i]]);
	    x2 = &x2_start[4 * i];
	    x3 = &x3_start[4 * i];

	    le =  &left[cptr[i] * 16];
	    ri =  &right[cptr[i] * 16];

	    for(j = 0; j < 4; j++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;
		for(k = 0; k < 4; k++)
		  {
		    ump_x1 += x1[k] * le[j * 4 + k];
		    ump_x2 += x2[k] * ri[j * 4 + k];
		  }
		x1px2[j] = ump_x1 * ump_x2;
	      }

	    for(j = 0; j < 4; j++)
	      x3[j] = 0.0;

	    for(j = 0; j < 4; j++)
	      for(k = 0; k < 4; k++)
		x3[k] +=  x1px2[j] *  EV[4 * j + k];	   

	    scale = 1;
	    for(j = 0; j < 4 && scale; j++)
	      scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);	    	   
	    	    
	    if(scale)
	      {		    
		for(j = 0; j < 4; j++)
		  x3[j] *= twotothe256;
		
		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;		
	      }	     
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];

	  le = &left[cptr[i] * 16];
	  ri = &right[cptr[i] * 16];

	  for(j = 0; j < 4; j++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      for(k = 0; k < 4; k++)
		{
		  ump_x1 += x1[k] * le[j * 4 + k];
		  ump_x2 += x2[k] * ri[j * 4 + k];
		}
	      x1px2[j] = ump_x1 * ump_x2;
	    }

	  for(j = 0; j < 4; j++)
	    x3[j] = 0.0;

	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      x3[k] +=  x1px2[j] *  EV[4 * j + k];
	
	  scale = 1;
	  for(j = 0; j < 4 && scale; j++)
	    scale = (x3[j] < minlikelihood && x3[j] > minusminlikelihood);

	  if(scale)
	    {		    
	      for(j = 0; j < 4; j++)
		x3[j] *= twotothe256;
	      
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;		
	    }	  
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#else


static void newviewGTRCAT_SAVE( int tipCase,  double *EV,  int *cptr,
				double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
				int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
				double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le,
  *ri,
  *x1,
  *x2,
  *x3,
  *x1_ptr = x1_start,
  *x2_ptr = x2_start, 
  *x3_ptr = x3_start, 
  EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));

  int 
    i, 
    j, 
    scale, 
    scaleGap = 0,
    addScale = 0;

  __m128d
    minlikelihood_sse = _mm_set1_pd( minlikelihood ),
                      sc = _mm_set1_pd(twotothe256),
                      EVV[8];  

  for(i = 0; i < 4; i++)
    for (j=0; j < 4; j++)
      EV_t[4 * j + i] = EV[4 * i + j];

  for(i = 0; i < 8; i++)
    EVV[i] = _mm_load_pd(&EV_t[i * 2]);

  

  {
    x1 = x1_gapColumn;	      
    x2 = x2_gapColumn;
    x3 = x3_gapColumn;

   

    le =  &left[maxCats * 16];	     	 
    ri =  &right[maxCats * 16];		   	  	  	  	         

   

    __m128d x1_0 = _mm_load_pd( &x1[0] );
    __m128d x1_2 = _mm_load_pd( &x1[2] );

    __m128d left_k0_0 = _mm_load_pd( &le[0] );
    __m128d left_k0_2 = _mm_load_pd( &le[2] );
    __m128d left_k1_0 = _mm_load_pd( &le[4] );
    __m128d left_k1_2 = _mm_load_pd( &le[6] );
    __m128d left_k2_0 = _mm_load_pd( &le[8] );
    __m128d left_k2_2 = _mm_load_pd( &le[10] );
    __m128d left_k3_0 = _mm_load_pd( &le[12] );
    __m128d left_k3_2 = _mm_load_pd( &le[14] );

    left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
    left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

    left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
    left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
    left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

    left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
    left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

    left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
    left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
    left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

    __m128d x2_0 = _mm_load_pd( &x2[0] );
    __m128d x2_2 = _mm_load_pd( &x2[2] );

    __m128d right_k0_0 = _mm_load_pd( &ri[0] );
    __m128d right_k0_2 = _mm_load_pd( &ri[2] );
    __m128d right_k1_0 = _mm_load_pd( &ri[4] );
    __m128d right_k1_2 = _mm_load_pd( &ri[6] );
    __m128d right_k2_0 = _mm_load_pd( &ri[8] );
    __m128d right_k2_2 = _mm_load_pd( &ri[10] );
    __m128d right_k3_0 = _mm_load_pd( &ri[12] );
    __m128d right_k3_2 = _mm_load_pd( &ri[14] );

    right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
    right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

    right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
    right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
    right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

    right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
    right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

    right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
    right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
    right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   

    __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
    __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

    __m128d EV_t_l0_k0 = EVV[0];
    __m128d EV_t_l0_k2 = EVV[1];
    __m128d EV_t_l1_k0 = EVV[2];
    __m128d EV_t_l1_k2 = EVV[3];
    __m128d EV_t_l2_k0 = EVV[4];
    __m128d EV_t_l2_k2 = EVV[5];
    __m128d EV_t_l3_k0 = EVV[6];
    __m128d EV_t_l3_k2 = EVV[7];

    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  

   

    if(tipCase != TIP_TIP)
    {    
      scale = 1;

      __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
      if(_mm_movemask_pd( v1 ) != 3)
        scale = 0;
      else
      {
        v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
        v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
        if(_mm_movemask_pd( v1 ) != 3)
          scale = 0;
      }

      if(scale)
      {		      
        _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
        _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      

        scaleGap = TRUE;	   
      }	
      else
      {
        _mm_store_pd(x3, EV_t_l0_k0);
        _mm_store_pd(&x3[2], EV_t_l2_k0);
      }
    }
    else
    {
      _mm_store_pd(x3, EV_t_l0_k0);
      _mm_store_pd(&x3[2], EV_t_l2_k0);
    }
  }
  
  switch(tipCase)
  {
    case TIP_TIP:      
      for (i = 0; i < n; i++)
      {
        if(noGap(x3_gap, i))
        {
          x1 = &(tipVector[4 * tipX1[i]]);
          x2 = &(tipVector[4 * tipX2[i]]);

          x3 = x3_ptr;

          if(isGap(x1_gap, i))
            le =  &left[maxCats * 16];
          else	  	  
            le =  &left[cptr[i] * 16];	  

          if(isGap(x2_gap, i))
            ri =  &right[maxCats * 16];
          else	 	  
            ri =  &right[cptr[i] * 16];

          __m128d x1_0 = _mm_load_pd( &x1[0] );
          __m128d x1_2 = _mm_load_pd( &x1[2] );

          __m128d left_k0_0 = _mm_load_pd( &le[0] );
          __m128d left_k0_2 = _mm_load_pd( &le[2] );
          __m128d left_k1_0 = _mm_load_pd( &le[4] );
          __m128d left_k1_2 = _mm_load_pd( &le[6] );
          __m128d left_k2_0 = _mm_load_pd( &le[8] );
          __m128d left_k2_2 = _mm_load_pd( &le[10] );
          __m128d left_k3_0 = _mm_load_pd( &le[12] );
          __m128d left_k3_2 = _mm_load_pd( &le[14] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

          __m128d x2_0 = _mm_load_pd( &x2[0] );
          __m128d x2_2 = _mm_load_pd( &x2[2] );

          __m128d right_k0_0 = _mm_load_pd( &ri[0] );
          __m128d right_k0_2 = _mm_load_pd( &ri[2] );
          __m128d right_k1_0 = _mm_load_pd( &ri[4] );
          __m128d right_k1_2 = _mm_load_pd( &ri[6] );
          __m128d right_k2_0 = _mm_load_pd( &ri[8] );
          __m128d right_k2_2 = _mm_load_pd( &ri[10] );
          __m128d right_k3_0 = _mm_load_pd( &ri[12] );
          __m128d right_k3_2 = _mm_load_pd( &ri[14] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );	  	  

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6];
          __m128d EV_t_l3_k2 = EVV[7];

          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	 

          _mm_store_pd(x3, EV_t_l0_k0);
          _mm_store_pd(&x3[2], EV_t_l2_k0);	  	 	   	    

          x3_ptr += 4;
        }
      }
      break;
    case TIP_INNER:      
      for (i = 0; i < n; i++)
      { 
        if(isGap(x3_gap, i))
        {
          if(scaleGap)		   		    
            addScale += wgt[i];
        }
        else
        {	      
          x1 = &(tipVector[4 * tipX1[i]]);

          x2 = x2_ptr;
          x3 = x3_ptr;

          if(isGap(x1_gap, i))
            le =  &left[maxCats * 16];
          else
            le =  &left[cptr[i] * 16];

          if(isGap(x2_gap, i))
          {		 
            ri =  &right[maxCats * 16];
            x2 = x2_gapColumn;
          }
          else
          {
            ri =  &right[cptr[i] * 16];
            x2 = x2_ptr;
            x2_ptr += 4;
          }	  	  	  	  

          __m128d x1_0 = _mm_load_pd( &x1[0] );
          __m128d x1_2 = _mm_load_pd( &x1[2] );

          __m128d left_k0_0 = _mm_load_pd( &le[0] );
          __m128d left_k0_2 = _mm_load_pd( &le[2] );
          __m128d left_k1_0 = _mm_load_pd( &le[4] );
          __m128d left_k1_2 = _mm_load_pd( &le[6] );
          __m128d left_k2_0 = _mm_load_pd( &le[8] );
          __m128d left_k2_2 = _mm_load_pd( &le[10] );
          __m128d left_k3_0 = _mm_load_pd( &le[12] );
          __m128d left_k3_2 = _mm_load_pd( &le[14] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

          __m128d x2_0 = _mm_load_pd( &x2[0] );
          __m128d x2_2 = _mm_load_pd( &x2[2] );

          __m128d right_k0_0 = _mm_load_pd( &ri[0] );
          __m128d right_k0_2 = _mm_load_pd( &ri[2] );
          __m128d right_k1_0 = _mm_load_pd( &ri[4] );
          __m128d right_k1_2 = _mm_load_pd( &ri[6] );
          __m128d right_k2_0 = _mm_load_pd( &ri[8] );
          __m128d right_k2_2 = _mm_load_pd( &ri[10] );
          __m128d right_k3_0 = _mm_load_pd( &ri[12] );
          __m128d right_k3_2 = _mm_load_pd( &ri[14] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6];
          __m128d EV_t_l3_k2 = EVV[7];


          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  

          scale = 1;

          __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;
          else
          {
            v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }

          if(scale)
          {		      
            _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
            _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      

            addScale += wgt[i];	  
          }	
          else
          {
            _mm_store_pd(x3, EV_t_l0_k0);
            _mm_store_pd(&x3[2], EV_t_l2_k0);
          }

          x3_ptr += 4;
        }

      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
      { 
        if(isGap(x3_gap, i))
        {
          if(scaleGap)		   		    
            addScale += wgt[i];
        }
        else
        {	     
          x3 = x3_ptr;

          if(isGap(x1_gap, i))
          {
            x1 = x1_gapColumn;
            le =  &left[maxCats * 16];
          }
          else
          {
            le =  &left[cptr[i] * 16];
            x1 = x1_ptr;
            x1_ptr += 4;
          }

          if(isGap(x2_gap, i))	
          {
            x2 = x2_gapColumn;
            ri =  &right[maxCats * 16];	    
          }
          else
          {
            ri =  &right[cptr[i] * 16];
            x2 = x2_ptr;
            x2_ptr += 4;
          }	 	  	  	  

          __m128d x1_0 = _mm_load_pd( &x1[0] );
          __m128d x1_2 = _mm_load_pd( &x1[2] );

          __m128d left_k0_0 = _mm_load_pd( &le[0] );
          __m128d left_k0_2 = _mm_load_pd( &le[2] );
          __m128d left_k1_0 = _mm_load_pd( &le[4] );
          __m128d left_k1_2 = _mm_load_pd( &le[6] );
          __m128d left_k2_0 = _mm_load_pd( &le[8] );
          __m128d left_k2_2 = _mm_load_pd( &le[10] );
          __m128d left_k3_0 = _mm_load_pd( &le[12] );
          __m128d left_k3_2 = _mm_load_pd( &le[14] );

          left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
          left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);

          left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
          left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);

          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
          left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
          left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);

          left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
          left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);

          left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
          left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);

          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
          left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
          left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);

          __m128d x2_0 = _mm_load_pd( &x2[0] );
          __m128d x2_2 = _mm_load_pd( &x2[2] );

          __m128d right_k0_0 = _mm_load_pd( &ri[0] );
          __m128d right_k0_2 = _mm_load_pd( &ri[2] );
          __m128d right_k1_0 = _mm_load_pd( &ri[4] );
          __m128d right_k1_2 = _mm_load_pd( &ri[6] );
          __m128d right_k2_0 = _mm_load_pd( &ri[8] );
          __m128d right_k2_2 = _mm_load_pd( &ri[10] );
          __m128d right_k3_0 = _mm_load_pd( &ri[12] );
          __m128d right_k3_2 = _mm_load_pd( &ri[14] );

          right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
          right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

          right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
          right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
          right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
          right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);

          right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
          right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);

          right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
          right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
          right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
          right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   

          __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
          __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );

          __m128d EV_t_l0_k0 = EVV[0];
          __m128d EV_t_l0_k2 = EVV[1];
          __m128d EV_t_l1_k0 = EVV[2];
          __m128d EV_t_l1_k2 = EVV[3];
          __m128d EV_t_l2_k0 = EVV[4];
          __m128d EV_t_l2_k2 = EVV[5];
          __m128d EV_t_l3_k0 = EVV[6];
          __m128d EV_t_l3_k2 = EVV[7];


          EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
          EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

          EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
          EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

          EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
          EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

          EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
          EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );

          EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
          EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
          EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

          EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  	 

          scale = 1;

          __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
          v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
          if(_mm_movemask_pd( v1 ) != 3)
            scale = 0;
          else
          {
            v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }

          if(scale)
          {		      
            _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
            _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      

            addScale += wgt[i];	  
          }	
          else
          {
            _mm_store_pd(x3, EV_t_l0_k0);
            _mm_store_pd(&x3[2], EV_t_l2_k0);
          }

          x3_ptr += 4;
        }
      }
      break;
    default:
      assert(0);
  }


  *scalerIncrement = addScale;
}


static void newviewGTRCAT( int tipCase,  double *EV,  int *cptr,
			   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
			   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1,
    *x2, 
    *x3, 
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));
    
  int 
    i, 
    j, 
    scale, 
    addScale = 0;
   
  __m128d
    minlikelihood_sse = _mm_set1_pd( minlikelihood ),
    sc = _mm_set1_pd(twotothe256),
    EVV[8];  
  
  for(i = 0; i < 4; i++)
    for (j=0; j < 4; j++)
      EV_t[4 * j + i] = EV[4 * i + j];
  
  for(i = 0; i < 8; i++)
    EVV[i] = _mm_load_pd(&EV_t[i * 2]);
  
  switch(tipCase)
    {
    case TIP_TIP:      
      for (i = 0; i < n; i++)
	{	 
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  
	  x3 = &x3_start[i * 4];
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];
	  
	  __m128d x1_0 = _mm_load_pd( &x1[0] );
	  __m128d x1_2 = _mm_load_pd( &x1[2] );
	  
	  __m128d left_k0_0 = _mm_load_pd( &le[0] );
	  __m128d left_k0_2 = _mm_load_pd( &le[2] );
	  __m128d left_k1_0 = _mm_load_pd( &le[4] );
	  __m128d left_k1_2 = _mm_load_pd( &le[6] );
	  __m128d left_k2_0 = _mm_load_pd( &le[8] );
	  __m128d left_k2_2 = _mm_load_pd( &le[10] );
	  __m128d left_k3_0 = _mm_load_pd( &le[12] );
	  __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	  left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	  left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	  
	  left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	  left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	  
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	  left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	  
	  left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	  left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	  
	  left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	  left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	  
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	  left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	  
	  __m128d x2_0 = _mm_load_pd( &x2[0] );
	  __m128d x2_2 = _mm_load_pd( &x2[2] );
	  
	  __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	  __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	  __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	  __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	  __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	  __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	  __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	  __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	  
	  right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	  right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	  right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	  right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	  
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	  right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	  
	  right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	  right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	  
	  right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	  right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	  
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	  right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	  
	  __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	  __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );	  	  

	  __m128d EV_t_l0_k0 = EVV[0];
	  __m128d EV_t_l0_k2 = EVV[1];
	  __m128d EV_t_l1_k0 = EVV[2];
	  __m128d EV_t_l1_k2 = EVV[3];
	  __m128d EV_t_l2_k0 = EVV[4];
	  __m128d EV_t_l2_k2 = EVV[5];
	  __m128d EV_t_l3_k0 = EVV[6];
	  __m128d EV_t_l3_k2 = EVV[7];
	  
	  EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	  EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	  
	  EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	  EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	  
	  EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	  
	  EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	  EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	  	  
	  EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	  EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	  EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	  
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	 
	  	  
	  _mm_store_pd(x3, EV_t_l0_k0);
	  _mm_store_pd(&x3[2], EV_t_l2_k0);	  	 	   	    
	}
      break;
    case TIP_INNER:      
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m128d x1_0 = _mm_load_pd( &x1[0] );
	  __m128d x1_2 = _mm_load_pd( &x1[2] );
	  
	  __m128d left_k0_0 = _mm_load_pd( &le[0] );
	  __m128d left_k0_2 = _mm_load_pd( &le[2] );
	  __m128d left_k1_0 = _mm_load_pd( &le[4] );
	  __m128d left_k1_2 = _mm_load_pd( &le[6] );
	  __m128d left_k2_0 = _mm_load_pd( &le[8] );
	  __m128d left_k2_2 = _mm_load_pd( &le[10] );
	  __m128d left_k3_0 = _mm_load_pd( &le[12] );
	  __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	  left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	  left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	  
	  left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	  left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	  
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	  left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	  
	  left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	  left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	  
	  left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	  left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	  
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	  left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	  
	  __m128d x2_0 = _mm_load_pd( &x2[0] );
	  __m128d x2_2 = _mm_load_pd( &x2[2] );
	  
	  __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	  __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	  __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	  __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	  __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	  __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	  __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	  __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	  
	  right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	  right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	  right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	  right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	  
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	  right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	  
	  right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	  right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	  
	  right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	  right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	  
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	  right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	  
	  __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	  __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
	  
	  __m128d EV_t_l0_k0 = EVV[0];
	  __m128d EV_t_l0_k2 = EVV[1];
	  __m128d EV_t_l1_k0 = EVV[2];
	  __m128d EV_t_l1_k2 = EVV[3];
	  __m128d EV_t_l2_k0 = EVV[4];
	  __m128d EV_t_l2_k2 = EVV[5];
	  __m128d EV_t_l3_k0 = EVV[6];
	  __m128d EV_t_l3_k2 = EVV[7];
	 
	  
	  EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	  EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	  
	  EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	  EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	  
	  EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	  
	  EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	  EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	  	  
	  EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	  EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	  EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	  
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  
	 
	  scale = 1;
	  	  	  	    
	  __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	  if(_mm_movemask_pd( v1 ) != 3)
	    scale = 0;
	  else
	    {
	      v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
	      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	      if(_mm_movemask_pd( v1 ) != 3)
		scale = 0;
	    }
	  	  
	  if(scale)
	    {		      
	      _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
	      _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
	      
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;
	    }	
	  else
	    {
	      _mm_store_pd(x3, EV_t_l0_k0);
	      _mm_store_pd(&x3[2], EV_t_l2_k0);
	    }
	 
	  	  
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m128d x1_0 = _mm_load_pd( &x1[0] );
	  __m128d x1_2 = _mm_load_pd( &x1[2] );
	  
	  __m128d left_k0_0 = _mm_load_pd( &le[0] );
	  __m128d left_k0_2 = _mm_load_pd( &le[2] );
	  __m128d left_k1_0 = _mm_load_pd( &le[4] );
	  __m128d left_k1_2 = _mm_load_pd( &le[6] );
	  __m128d left_k2_0 = _mm_load_pd( &le[8] );
	  __m128d left_k2_2 = _mm_load_pd( &le[10] );
	  __m128d left_k3_0 = _mm_load_pd( &le[12] );
	  __m128d left_k3_2 = _mm_load_pd( &le[14] );
	  
	  left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	  left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	  
	  left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	  left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	  
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	  left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	  left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	  
	  left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	  left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	  
	  left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	  left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	  
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	  left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	  left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	  
	  __m128d x2_0 = _mm_load_pd( &x2[0] );
	  __m128d x2_2 = _mm_load_pd( &x2[2] );
	  
	  __m128d right_k0_0 = _mm_load_pd( &ri[0] );
	  __m128d right_k0_2 = _mm_load_pd( &ri[2] );
	  __m128d right_k1_0 = _mm_load_pd( &ri[4] );
	  __m128d right_k1_2 = _mm_load_pd( &ri[6] );
	  __m128d right_k2_0 = _mm_load_pd( &ri[8] );
	  __m128d right_k2_2 = _mm_load_pd( &ri[10] );
	  __m128d right_k3_0 = _mm_load_pd( &ri[12] );
	  __m128d right_k3_2 = _mm_load_pd( &ri[14] );
	  
	  right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	  right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	  
	  right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	  right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	  
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	  right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	  right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	  
	  right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	  right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	  
	  right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	  right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	  
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	  right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	  right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
	  
	  __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	  __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
	  
	  __m128d EV_t_l0_k0 = EVV[0];
	  __m128d EV_t_l0_k2 = EVV[1];
	  __m128d EV_t_l1_k0 = EVV[2];
	  __m128d EV_t_l1_k2 = EVV[3];
	  __m128d EV_t_l2_k0 = EVV[4];
	  __m128d EV_t_l2_k2 = EVV[5];
	  __m128d EV_t_l3_k0 = EVV[6];
	  __m128d EV_t_l3_k2 = EVV[7];
	 
	  
	  EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	  EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	  
	  EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	  EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	  
	  EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	  EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	  
	  EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	  EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	  	  
	  EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	  EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	  EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	  
	  EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );	  	 	    		  	 

	  scale = 1;
	  	  
	  __m128d v1 = _mm_and_pd(EV_t_l0_k0, absMask.m);
	  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	  if(_mm_movemask_pd( v1 ) != 3)
	    scale = 0;
	  else
	    {
	      v1 = _mm_and_pd(EV_t_l2_k0, absMask.m);
	      v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	      if(_mm_movemask_pd( v1 ) != 3)
		scale = 0;
	    }
	  	  
	  if(scale)
	    {		      
	      _mm_store_pd(&x3[0], _mm_mul_pd(EV_t_l0_k0, sc));
	      _mm_store_pd(&x3[2], _mm_mul_pd(EV_t_l2_k0, sc));	      	      
	      
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;
	    }	
	  else
	    {
	      _mm_store_pd(x3, EV_t_l0_k0);
	      _mm_store_pd(&x3[2], EV_t_l2_k0);
	    }
	  	  
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
}



#endif





#ifdef __SIM_SSE3

static void newviewGTRGAMMA(int tipCase,
			    double *x1_start, double *x2_start, double *x3_start,
			    double *EV, double *tipVector,
			    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			    const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
			    )
{
  int 
    i, 
    j, 
    k, 
    l,
    addScale = 0;
  
  double
    *x1,
    *x2,
    *x3,
    max,
    maxima[2] __attribute__ ((aligned (BYTE_ALIGNMENT))),       
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));      
    
  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);
   
  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT))), *uX2, umpX2[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);
		 
		}
	  }   	
	  
	for (i = 0; i < n; i++)
	  {
	    x3 = &x3_start[i * 16];

	    
	    uX1 = &umpX1[16 * tipX1[i]];
	    uX2 = &umpX2[16 * tipX2[i]];	   	    	    
	    
	    for (j = 0; j < 4; j++)
	       {				 		  		  		   
		 __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		 __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		 				  
		   
		 __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
		 __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
 		 

		 //
		 // multiply left * right
		 //
		 
		 __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
		 __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );
		 
		 
		 //
		 // multiply with EV matrix (!?)
		 //
		 
		 __m128d EV_t_l0_k0 = EVV[0];
		 __m128d EV_t_l0_k2 = EVV[1];
		 __m128d EV_t_l1_k0 = EVV[2];
		 __m128d EV_t_l1_k2 = EVV[3];
		 __m128d EV_t_l2_k0 = EVV[4];
		 __m128d EV_t_l2_k2 = EVV[5];
		 __m128d EV_t_l3_k0 = EVV[6]; 
		 __m128d EV_t_l3_k2 = EVV[7];
		 
		 EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		 EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		 
		 EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		 EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		 
		 EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		 
		 EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		 EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		 
		 EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		 EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		 EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		 
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		 
		 _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
		 _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
	       }
	  }
      }
      break;
    case TIP_INNER:
      {	
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);		 
		}
	  }

	 for (i = 0; i < n; i++)
	   {
	     __m128d maxv =_mm_setzero_pd();
	     
	     x2 = &x2_start[i * 16];
	     x3 = &x3_start[i * 16];

	     uX1 = &umpX1[16 * tipX1[i]];	     

	     for (j = 0; j < 4; j++)
	       {

		 //
		 // multiply/add right side
		 //
		 double *x2_p = &x2[j*4];
		 double *right_k0_p = &right[j*16];
		 double *right_k1_p = &right[j*16 + 1*4];
		 double *right_k2_p = &right[j*16 + 2*4];
		 double *right_k3_p = &right[j*16 + 3*4];
		 __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		 __m128d x2_2 = _mm_load_pd( &x2_p[2] );

		 __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		 __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		 __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		 __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		 __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		 __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		 __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		 __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );



		 right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		 right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);

		 right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		 right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);

		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		 right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);


		 right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		 right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);


		 right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		 right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);

		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		 right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);

		 {
		   //
		   // load left side from tip vector
		   //
		   
		   __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		   __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		 
		 
		   //
		   // multiply left * right
		   //
		   
		   __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
		   __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
		   
		   
		   //
		   // multiply with EV matrix (!?)
		   //		   		  

		   __m128d EV_t_l0_k0 = EVV[0];
		   __m128d EV_t_l0_k2 = EVV[1];
		   __m128d EV_t_l1_k0 = EVV[2];
		   __m128d EV_t_l1_k2 = EVV[3];
		   __m128d EV_t_l2_k0 = EVV[4];
		   __m128d EV_t_l2_k2 = EVV[5];
		   __m128d EV_t_l3_k0 = EVV[6]; 
		   __m128d EV_t_l3_k2 = EVV[7];

		   
		   EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		   EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		   EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		   
		   EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		   EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		   
		   EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		   EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		   
		   EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		   EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		   EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		   		   
		   EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		   EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		   EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		   
		   EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		   
		   values[j * 2]     = EV_t_l0_k0;
		   values[j * 2 + 1] = EV_t_l2_k0;		   		   
		   
		   maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		   maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   
		 }
	       }

	     
	     _mm_store_pd(maxima, maxv);

	     max = MAX(maxima[0], maxima[1]);

	     if(max < minlikelihood)
	       {
		 __m128d sv = _mm_set1_pd(twotothe256);
	       		       	   	 	     
		 _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		 _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		 _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		 _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		 _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		 _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		 _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		 _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;
	       }
	     else
	       {
		 _mm_store_pd(&x3[0], values[0]);	   
		 _mm_store_pd(&x3[2], values[1]);
		 _mm_store_pd(&x3[4], values[2]);
		 _mm_store_pd(&x3[6], values[3]);
		 _mm_store_pd(&x3[8], values[4]);	   
		 _mm_store_pd(&x3[10], values[5]);
		 _mm_store_pd(&x3[12], values[6]);
		 _mm_store_pd(&x3[14], values[7]);
	       }
	   }
      }
      break;
    case INNER_INNER:     
     for (i = 0; i < n; i++)
       {
	 __m128d maxv =_mm_setzero_pd();
	 

	 x1 = &x1_start[i * 16];
	 x2 = &x2_start[i * 16];
	 x3 = &x3_start[i * 16];
	 
	 for (j = 0; j < 4; j++)
	   {
	     
	     double *x1_p = &x1[j*4];
	     double *left_k0_p = &left[j*16];
	     double *left_k1_p = &left[j*16 + 1*4];
	     double *left_k2_p = &left[j*16 + 2*4];
	     double *left_k3_p = &left[j*16 + 3*4];
	     
	     __m128d x1_0 = _mm_load_pd( &x1_p[0] );
	     __m128d x1_2 = _mm_load_pd( &x1_p[2] );
	     
	     __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
	     __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
	     __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
	     __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
	     __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
	     __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
	     __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
	     __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
	     
	     left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	     left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	     
	     left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	     left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	     
	     left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	     left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	     left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	     
	     left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	     left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	     
	     left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	     left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	     
	     left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	     left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	     left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	     
	     
	     //
	     // multiply/add right side
	     //
	     double *x2_p = &x2[j*4];
	     double *right_k0_p = &right[j*16];
	     double *right_k1_p = &right[j*16 + 1*4];
	     double *right_k2_p = &right[j*16 + 2*4];
	     double *right_k3_p = &right[j*16 + 3*4];
	     __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	     __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	     
	     __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	     __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	     __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	     __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	     __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	     __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	     __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	     __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	     
	     right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	     right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	     
	     right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	     right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	     
	     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	     right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	     
	     right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	     right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	     
	     
	     right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	     right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	     
	     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	     right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   

             //
             // multiply left * right
             //

	     __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	     __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );


             //
             // multiply with EV matrix (!?)
             //	     

	     __m128d EV_t_l0_k0 = EVV[0];
	     __m128d EV_t_l0_k2 = EVV[1];
	     __m128d EV_t_l1_k0 = EVV[2];
	     __m128d EV_t_l1_k2 = EVV[3];
	     __m128d EV_t_l2_k0 = EVV[4];
	     __m128d EV_t_l2_k2 = EVV[5];
	     __m128d EV_t_l3_k0 = EVV[6]; 
	     __m128d EV_t_l3_k2 = EVV[7];


	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );

	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );

	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );

	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );


	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
            EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
            EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );

            EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );

	    
	    values[j * 2] = EV_t_l0_k0;
	    values[j * 2 + 1] = EV_t_l2_k0;            	   	    

	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
           }
	 	 
	 
	 _mm_store_pd(maxima, maxv);
	 
	 max = MAX(maxima[0], maxima[1]);
	 
	 if(max < minlikelihood)
	   {
	     __m128d sv = _mm_set1_pd(twotothe256);
	       		       	   	 	     
	     _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	     _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	     _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	     _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	     _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	     _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	     _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	     _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
	     
	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;
	   }
	 else
	   {
	     _mm_store_pd(&x3[0], values[0]);	   
	     _mm_store_pd(&x3[2], values[1]);
	     _mm_store_pd(&x3[4], values[2]);
	     _mm_store_pd(&x3[6], values[3]);
	     _mm_store_pd(&x3[8], values[4]);	   
	     _mm_store_pd(&x3[10], values[5]);
	     _mm_store_pd(&x3[12], values[6]);
	     _mm_store_pd(&x3[14], values[7]);
	   }	 
       }
   
     break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;

}


static void newviewGTRGAMMA_GAPPED(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, int gapLength, 
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn)
{
  int 
    i, 
    j, 
    k, 
    l,
    addScale = 0, 
    scaleGap = 0;
  
  double
    *x1,
    *x2,
    *x3,   
    max,
    maxima[2] __attribute__ ((aligned (BYTE_ALIGNMENT))),        
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));      
    
  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);
   

  
  /*
    {           
    int v = gapLength - (gapLength % 4);
    
    for(i = 0; i < v; i+=4)	
    _mm_storeu_si128(&x3_gap[i], _mm_and_si128(_mm_loadu_si128(&x1_gap[i]), _mm_loadu_si128(&x2_gap[i])));
    
    for(; i < gapLength; i++)
    x3_gap[i] = x1_gap[i] & x2_gap[i];
    }
  */
 
  for(i = 0; i < gapLength; i++)
    x3_gap[i] = x1_gap[i] & x2_gap[i];    

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT))), *uX2, umpX2[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);
		 
		}
	  }   		  
		 	   
	uX1 = &umpX1[240];
	uX2 = &umpX2[240];	   	    	    
	
	for (j = 0; j < 4; j++)
	  {				 		  		  		   
	    __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	    __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	    	    
	    __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
	    __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
	    
	    __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
	    __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );		    		    		   
	    
	    __m128d EV_t_l0_k0 = EVV[0];
	    __m128d EV_t_l0_k2 = EVV[1];
	    __m128d EV_t_l1_k0 = EVV[2];
	    __m128d EV_t_l1_k2 = EVV[3];
	    __m128d EV_t_l2_k0 = EVV[4];
	    __m128d EV_t_l2_k2 = EVV[5];
	    __m128d EV_t_l3_k0 = EVV[6]; 
	    __m128d EV_t_l3_k2 = EVV[7];
	    
	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	    
	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	    
	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	    
	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	    
	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	    
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	    
	    _mm_store_pd( &x3_gapColumn[j * 4 + 0], EV_t_l0_k0 );
	    _mm_store_pd( &x3_gapColumn[j * 4 + 2], EV_t_l2_k0 );
	  }       

	for (i = 0; i < n; i++)
	  {
	    x3 = &x3_start[i * 16];

	    if(!(x3_gap[i / 32] & mask32[i % 32]))	     
	      {
		uX1 = &umpX1[16 * tipX1[i]];
		uX2 = &umpX2[16 * tipX2[i]];	   	    	    
		
		for (j = 0; j < 4; j++)
		  {				 		  		  		   
		    __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		    __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		    
		    
		    __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
		    __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
		    
		    
		    //
		    // multiply left * right
		    //
		    
		    __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
		    __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );
		    
		    
		    //
		    // multiply with EV matrix (!?)
		    //
		    
		    __m128d EV_t_l0_k0 = EVV[0];
		    __m128d EV_t_l0_k2 = EVV[1];
		    __m128d EV_t_l1_k0 = EVV[2];
		    __m128d EV_t_l1_k2 = EVV[3];
		    __m128d EV_t_l2_k0 = EVV[4];
		    __m128d EV_t_l2_k2 = EVV[5];
		    __m128d EV_t_l3_k0 = EVV[6]; 
		    __m128d EV_t_l3_k2 = EVV[7];
		    
		    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		    
		    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		    
		    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		    
		    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		    
		    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		    
		    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		    
		    _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
		    _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
		  }
	       }
	  }
      }
      break;
    case TIP_INNER:
      {	
	double 
	  *uX1, 
	  umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));		 

	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);		 
		}
	  }

	{
	  __m128d maxv =_mm_setzero_pd();
	  
	  scaleGap = 0;
	  
	  x2 = x2_gapColumn;			 
	  x3 = x3_gapColumn;
	  
	  uX1 = &umpX1[240];	     
	  
	  for (j = 0; j < 4; j++)
	    {		     		   
	      double *x2_p = &x2[j*4];
	      double *right_k0_p = &right[j*16];
	      double *right_k1_p = &right[j*16 + 1*4];
	      double *right_k2_p = &right[j*16 + 2*4];
	      double *right_k3_p = &right[j*16 + 3*4];
	      __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	      __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	      
	      __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	      __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	      __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	      __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	      __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	      __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	      __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	      __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	      	      
	      right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	      right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	      
	      right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	      right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	      
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	      right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	      	       
	      right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	      right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	      	       
	      right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	      right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	      
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	      right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);
	      
	      __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	      __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	      
	      __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
	      __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
	      
	      __m128d EV_t_l0_k0 = EVV[0];
	      __m128d EV_t_l0_k2 = EVV[1];
	      __m128d EV_t_l1_k0 = EVV[2];
	      __m128d EV_t_l1_k2 = EVV[3];
	      __m128d EV_t_l2_k0 = EVV[4];
	      __m128d EV_t_l2_k2 = EVV[5];
	      __m128d EV_t_l3_k0 = EVV[6]; 
	      __m128d EV_t_l3_k2 = EVV[7];
	      
	      EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	      EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	      
	      EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	      EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	      
	      EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	      
	      EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	      EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	      
	      EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	      EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	      EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	      
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	      
	      values[j * 2]     = EV_t_l0_k0;
	      values[j * 2 + 1] = EV_t_l2_k0;		   		   
	      
	      maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	      maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   	     		   
	    }

	  
	  _mm_store_pd(maxima, maxv);
		 
	  max = MAX(maxima[0], maxima[1]);
	  
	  if(max < minlikelihood)
	    {
	      scaleGap = 1;
	      
	      __m128d sv = _mm_set1_pd(twotothe256);
	      
	      _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	      _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	      _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	      _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	      _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	      _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	      _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	      _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     	      	     
	    }
	  else
	    {
	      _mm_store_pd(&x3[0], values[0]);	   
	      _mm_store_pd(&x3[2], values[1]);
	      _mm_store_pd(&x3[4], values[2]);
	      _mm_store_pd(&x3[6], values[3]);
	      _mm_store_pd(&x3[8], values[4]);	   
	      _mm_store_pd(&x3[10], values[5]);
	      _mm_store_pd(&x3[12], values[6]);
	      _mm_store_pd(&x3[14], values[7]);
	    }
	}
	
	for (i = 0; i < n; i++)
	   {
	     if((x3_gap[i / 32] & mask32[i % 32]))
	       {	       
		 if(scaleGap)
		   {
		     if(useFastScaling)
		       addScale += wgt[i];
		     else
		       ex3[i]  += 1;
		   }
	       }
	     else
	       {				 
		 __m128d maxv =_mm_setzero_pd();

		 if(x2_gap[i / 32] & mask32[i % 32])
		   x2 = x2_gapColumn;
		 else
		   x2 = &x2_start[i * 16];
		 
		 x3 = &x3_start[i * 16];
		 
		 uX1 = &umpX1[16 * tipX1[i]];	     
		 
		 
		 for (j = 0; j < 4; j++)
		   {		     		   
		     double *x2_p = &x2[j*4];
		     double *right_k0_p = &right[j*16];
		     double *right_k1_p = &right[j*16 + 1*4];
		     double *right_k2_p = &right[j*16 + 2*4];
		     double *right_k3_p = &right[j*16 + 3*4];
		     __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		     __m128d x2_2 = _mm_load_pd( &x2_p[2] );
		     
		     __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		     __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		     __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		     __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		     __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		     __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		     __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		     __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
		     
		     		     
		     right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		     right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
		     
		     right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		     right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
		     
		     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		     right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
		     
		     
		     right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		     right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
		     
		     
		     right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		     right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
		     
		     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		     right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);
		     
		     {
		       //
		       // load left side from tip vector
		       //
		       
		       __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		       __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		       
		       
		       //
		       // multiply left * right
			   //
		       
		       __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
		       __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
		       
		       
		       //
		       // multiply with EV matrix (!?)
		       //		   		  
		       
		       __m128d EV_t_l0_k0 = EVV[0];
		       __m128d EV_t_l0_k2 = EVV[1];
		       __m128d EV_t_l1_k0 = EVV[2];
		       __m128d EV_t_l1_k2 = EVV[3];
		       __m128d EV_t_l2_k0 = EVV[4];
		       __m128d EV_t_l2_k2 = EVV[5];
		       __m128d EV_t_l3_k0 = EVV[6]; 
		       __m128d EV_t_l3_k2 = EVV[7];
		       
		       
		       EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		       EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		       EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		       
		       EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		       EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		       
		       EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		       EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		       
		       EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		       EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		       EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		       
		       EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		       EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		       EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		       
		       EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		       
		       values[j * 2]     = EV_t_l0_k0;
		       values[j * 2 + 1] = EV_t_l2_k0;		   		   
			   
		       maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		       maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   
		     }		   
		   }

	     
		 _mm_store_pd(maxima, maxv);
		 
		 max = MAX(maxima[0], maxima[1]);
		 
		 if(max < minlikelihood)
		   {
		     __m128d sv = _mm_set1_pd(twotothe256);
		     
		     _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		     _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		     _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		     _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		     _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		     _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		     _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		     _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		     
		     if(useFastScaling)
		       addScale += wgt[i];
		     else
		       ex3[i]  += 1;
		   }
		 else
		   {
		     _mm_store_pd(&x3[0], values[0]);	   
		     _mm_store_pd(&x3[2], values[1]);
		     _mm_store_pd(&x3[4], values[2]);
		     _mm_store_pd(&x3[6], values[3]);
		     _mm_store_pd(&x3[8], values[4]);	   
		     _mm_store_pd(&x3[10], values[5]);
		     _mm_store_pd(&x3[12], values[6]);
		     _mm_store_pd(&x3[14], values[7]);
		   }
	       }
	   }
      }
      break;
    case INNER_INNER:         
      {
	__m128d maxv =_mm_setzero_pd();
	
	scaleGap = 0;
	
	x1 = x1_gapColumn;	     	    
	x2 = x2_gapColumn;	    
	x3 = x3_gapColumn;
	
	for (j = 0; j < 4; j++)
	  {
	    
	    double *x1_p = &x1[j*4];
	    double *left_k0_p = &left[j*16];
	    double *left_k1_p = &left[j*16 + 1*4];
	    double *left_k2_p = &left[j*16 + 2*4];
	    double *left_k3_p = &left[j*16 + 3*4];
	    
	    __m128d x1_0 = _mm_load_pd( &x1_p[0] );
	    __m128d x1_2 = _mm_load_pd( &x1_p[2] );
	    
	    __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
	    __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
	    __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
	    __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
	    __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
	    __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
	    __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
	    __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
	    
	    left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	    left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	    
	    left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	    left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	    
	    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	    left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	    
	    left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	    left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	    
	    left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	    left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	    
	    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	    left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	    
	    
	    double *x2_p = &x2[j*4];
	    double *right_k0_p = &right[j*16];
	    double *right_k1_p = &right[j*16 + 1*4];
	    double *right_k2_p = &right[j*16 + 2*4];
	    double *right_k3_p = &right[j*16 + 3*4];
	    __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	    __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	    
	    __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	    __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	    __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	    __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	    __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	    __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	    __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	    __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	    
	    right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	    right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	    
	    right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	    right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	    
	    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	    right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	    
	    right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	    right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	    	    
	    right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	    right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	    
	    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	    right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   		 		
	    
	    __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	    __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );		 		 	   
	    
	    __m128d EV_t_l0_k0 = EVV[0];
	    __m128d EV_t_l0_k2 = EVV[1];
	    __m128d EV_t_l1_k0 = EVV[2];
	    __m128d EV_t_l1_k2 = EVV[3];
	    __m128d EV_t_l2_k0 = EVV[4];
	    __m128d EV_t_l2_k2 = EVV[5];
	    __m128d EV_t_l3_k0 = EVV[6]; 
	    __m128d EV_t_l3_k2 = EVV[7];
	    
	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	    
	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	    
	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	    
	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	    
	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	    
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	    
	    
	    values[j * 2] = EV_t_l0_k0;
	    values[j * 2 + 1] = EV_t_l2_k0;            	   	    
	    
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
	  }
		     
	_mm_store_pd(maxima, maxv);
	
	max = MAX(maxima[0], maxima[1]);
	
	if(max < minlikelihood)
	  {
	    __m128d sv = _mm_set1_pd(twotothe256);
	    
	    scaleGap = 1;
	    
	    _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	    _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	    _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	    _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	    _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	    _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	    _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	    _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     	    	 
	  }
	else
	  {
	    _mm_store_pd(&x3[0], values[0]);	   
	    _mm_store_pd(&x3[2], values[1]);
	    _mm_store_pd(&x3[4], values[2]);
	    _mm_store_pd(&x3[6], values[3]);
	    _mm_store_pd(&x3[8], values[4]);	   
	    _mm_store_pd(&x3[10], values[5]);
	    _mm_store_pd(&x3[12], values[6]);
	    _mm_store_pd(&x3[14], values[7]);
	  }
      }

 
     for (i = 0; i < n; i++)
       { 
	 if(x3_gap[i / 32] & mask32[i % 32])
	   {	     
	     if(scaleGap)
	       {
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1; 	       
	       }
	   }
	 else
	   {
	     __m128d maxv =_mm_setzero_pd();
	     
	     if(x1_gap[i / 32] & mask32[i % 32])
	       x1 = x1_gapColumn;
	     else
	       x1 = &x1_start[i * 16];

	     if(x2_gap[i / 32] & mask32[i % 32])
	       x2 = x2_gapColumn;
	     else
	       x2 = &x2_start[i * 16];


	     x3 = &x3_start[i * 16];
	     
	     for (j = 0; j < 4; j++)
	       {
		 
		 double *x1_p = &x1[j*4];
		 double *left_k0_p = &left[j*16];
		 double *left_k1_p = &left[j*16 + 1*4];
		 double *left_k2_p = &left[j*16 + 2*4];
		 double *left_k3_p = &left[j*16 + 3*4];
		 
		 __m128d x1_0 = _mm_load_pd( &x1_p[0] );
		 __m128d x1_2 = _mm_load_pd( &x1_p[2] );
		 
		 __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
		 __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
		 __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
		 __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
		 __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
		 __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
		 __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
		 __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
		 
		 left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
		 left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
		 
		 left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
		 left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
		 
		 left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
		 left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
		 left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
		 
		 left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
		 left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
		 
		 left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
		 left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
		 
		 left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
		 left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
		 left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
		 
		 
		 //
		 // multiply/add right side
		 //
		 double *x2_p = &x2[j*4];
		 double *right_k0_p = &right[j*16];
		 double *right_k1_p = &right[j*16 + 1*4];
		 double *right_k2_p = &right[j*16 + 2*4];
		 double *right_k3_p = &right[j*16 + 3*4];
		 __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		 __m128d x2_2 = _mm_load_pd( &x2_p[2] );
		 
		 __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		 __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		 __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		 __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		 __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		 __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		 __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		 __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
		 
		 right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		 right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
		 
		 right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		 right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
		 
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		 right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
		 
		 right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		 right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
		 
		 
		 right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		 right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
		 
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		 right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
		 
		 //
		 // multiply left * right
		 //
		 
		 __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
		 __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
		 
		 
		 //
		 // multiply with EV matrix (!?)
		 //	     
		 
		 __m128d EV_t_l0_k0 = EVV[0];
		 __m128d EV_t_l0_k2 = EVV[1];
		 __m128d EV_t_l1_k0 = EVV[2];
		 __m128d EV_t_l1_k2 = EVV[3];
		 __m128d EV_t_l2_k0 = EVV[4];
		 __m128d EV_t_l2_k2 = EVV[5];
		 __m128d EV_t_l3_k0 = EVV[6]; 
		 __m128d EV_t_l3_k2 = EVV[7];
		 
		 
		 EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		 EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		 
		 EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		 EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		 
		 EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		 
		 EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		 EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		 
		 
		 EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		 EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		 EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		 
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		 
		 
		 values[j * 2] = EV_t_l0_k0;
		 values[j * 2 + 1] = EV_t_l2_k0;            	   	    
		 
		 maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		 maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
	       }
	     
	     
	     _mm_store_pd(maxima, maxv);
	     
	     max = MAX(maxima[0], maxima[1]);
	     
	     if(max < minlikelihood)
	       {
		 __m128d sv = _mm_set1_pd(twotothe256);
		 
		 _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		 _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		 _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		 _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		 _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		 _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		 _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		 _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;
	       }
	     else
	       {
		 _mm_store_pd(&x3[0], values[0]);	   
		 _mm_store_pd(&x3[2], values[1]);
		 _mm_store_pd(&x3[4], values[2]);
		 _mm_store_pd(&x3[6], values[3]);
		 _mm_store_pd(&x3[8], values[4]);	   
		 _mm_store_pd(&x3[10], values[5]);
		 _mm_store_pd(&x3[12], values[6]);
		 _mm_store_pd(&x3[14], values[7]);
	       }	 
	   }
       }    
     break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;
}


static void newviewGTRGAMMA_GAPPED_SAVE(int tipCase,
					double *x1_start, double *x2_start, double *x3_start,
					double *EV, double *tipVector,
					int *ex3, unsigned char *tipX1, unsigned char *tipX2,
					const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
					unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, 
					double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn)
{
  int     
    i, 
    j, 
    k, 
    l,
    addScale = 0, 
    scaleGap = 0;
  
  double
    *x1,
    *x2,
    *x3,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start,       
    max,
    maxima[2] __attribute__ ((aligned (BYTE_ALIGNMENT))),        
    EV_t[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));      
    
  __m128d 
    values[8],
    EVV[8];  

  for(k = 0; k < 4; k++)
    for (l=0; l < 4; l++)
      EV_t[4 * l + k] = EV[4 * k + l];

  for(k = 0; k < 8; k++)
    EVV[k] = _mm_load_pd(&EV_t[k * 2]);      
 
  

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT))), *uX2, umpX2[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));


	for (i = 1; i < 16; i++)
	  {	    
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{			 	 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  __m128d left1 = _mm_load_pd(&right[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&right[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX2[i*16 + j*4 + k], acc);
		 
		}
	  }   		  
	
	uX1 = &umpX1[240];
	uX2 = &umpX2[240];	   	    	    
	
	for (j = 0; j < 4; j++)
	  {				 		  		  		   
	    __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	    __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	    	    
	    __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
	    __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
	    
	    __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
	    __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );		    		    		   
	    
	    __m128d EV_t_l0_k0 = EVV[0];
	    __m128d EV_t_l0_k2 = EVV[1];
	    __m128d EV_t_l1_k0 = EVV[2];
	    __m128d EV_t_l1_k2 = EVV[3];
	    __m128d EV_t_l2_k0 = EVV[4];
	    __m128d EV_t_l2_k2 = EVV[5];
	    __m128d EV_t_l3_k0 = EVV[6]; 
	    __m128d EV_t_l3_k2 = EVV[7];
	    
	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	    
	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	    
	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	    
	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	    
	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	    
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	    	  
	    _mm_store_pd( &x3_gapColumn[j * 4 + 0], EV_t_l0_k0 );
	    _mm_store_pd( &x3_gapColumn[j * 4 + 2], EV_t_l2_k0 );	   
	  }  
	
       
	x3 = x3_start;
	
	for (i = 0; i < n; i++)
	  {	    
	    if(!(x3_gap[i / 32] & mask32[i % 32]))	     
	      {
		uX1 = &umpX1[16 * tipX1[i]];
		uX2 = &umpX2[16 * tipX2[i]];	   	    	    		
		
		for (j = 0; j < 4; j++)
		  {				 		  		  		   
		    __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		    __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		    
		    
		    __m128d uX2_k0_sse = _mm_load_pd( &uX2[j * 4] );
		    __m128d uX2_k2_sse = _mm_load_pd( &uX2[j * 4 + 2] );
		    
		    
		    //
		    // multiply left * right
		    //
		    
		    __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, uX2_k0_sse );
		    __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, uX2_k2_sse );
		    
		    
		    //
		    // multiply with EV matrix (!?)
		    //
		    
		    __m128d EV_t_l0_k0 = EVV[0];
		    __m128d EV_t_l0_k2 = EVV[1];
		    __m128d EV_t_l1_k0 = EVV[2];
		    __m128d EV_t_l1_k2 = EVV[3];
		    __m128d EV_t_l2_k0 = EVV[4];
		    __m128d EV_t_l2_k2 = EVV[5];
		    __m128d EV_t_l3_k0 = EVV[6]; 
		    __m128d EV_t_l3_k2 = EVV[7];
		    
		    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		    
		    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		    
		    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		    
		    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		    
		    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		    
		    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		    
		    _mm_store_pd( &x3[j * 4 + 0], EV_t_l0_k0 );
		    _mm_store_pd( &x3[j * 4 + 2], EV_t_l2_k0 );
		  }
		
		x3 += 16;
	      }
	  }
      }
      break;
    case TIP_INNER:
      {	
	double 
	  *uX1, 
	  umpX1[256] __attribute__ ((aligned (BYTE_ALIGNMENT)));		 

	for (i = 1; i < 16; i++)
	  {
	    __m128d x1_1 = _mm_load_pd(&(tipVector[i*4]));
	    __m128d x1_2 = _mm_load_pd(&(tipVector[i*4 + 2]));	   

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m128d left1 = _mm_load_pd(&left[j*16 + k*4]);
		  __m128d left2 = _mm_load_pd(&left[j*16 + k*4 + 2]);
		  
		  __m128d acc = _mm_setzero_pd();

		  acc = _mm_add_pd(acc, _mm_mul_pd(left1, x1_1));
		  acc = _mm_add_pd(acc, _mm_mul_pd(left2, x1_2));
		  		  
		  acc = _mm_hadd_pd(acc, acc);
		  _mm_storel_pd(&umpX1[i*16 + j*4 + k], acc);		 
		}
	  }

	{
	  __m128d maxv =_mm_setzero_pd();
	  
	  scaleGap = 0;
	  
	  x2 = x2_gapColumn;			 
	  x3 = x3_gapColumn;
	  
	  uX1 = &umpX1[240];	     
	  
	  for (j = 0; j < 4; j++)
	    {		     		   
	      double *x2_p = &x2[j*4];
	      double *right_k0_p = &right[j*16];
	      double *right_k1_p = &right[j*16 + 1*4];
	      double *right_k2_p = &right[j*16 + 2*4];
	      double *right_k3_p = &right[j*16 + 3*4];
	      __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	      __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	      
	      __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	      __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	      __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	      __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	      __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	      __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	      __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	      __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	      	      
	      right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	      right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	      
	      right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	      right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	      
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	      right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	      right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	      	       
	      right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	      right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	      	       
	      right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	      right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	      
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	      right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	      right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);
	      
	      __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
	      __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
	      
	      __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
	      __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
	      
	      __m128d EV_t_l0_k0 = EVV[0];
	      __m128d EV_t_l0_k2 = EVV[1];
	      __m128d EV_t_l1_k0 = EVV[2];
	      __m128d EV_t_l1_k2 = EVV[3];
	      __m128d EV_t_l2_k0 = EVV[4];
	      __m128d EV_t_l2_k2 = EVV[5];
	      __m128d EV_t_l3_k0 = EVV[6]; 
	      __m128d EV_t_l3_k2 = EVV[7];
	      
	      EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	      EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	      
	      EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	      EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	      
	      EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	      EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	      
	      EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	      EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	      
	      EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	      EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	      EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	      
	      EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	      
	      values[j * 2]     = EV_t_l0_k0;
	      values[j * 2 + 1] = EV_t_l2_k0;		   		   
	      
	      maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	      maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   	     		   
	    }

	  
	  _mm_store_pd(maxima, maxv);
		 
	  max = MAX(maxima[0], maxima[1]);
	  
	  if(max < minlikelihood)
	    {
	      scaleGap = 1;
	      
	      __m128d sv = _mm_set1_pd(twotothe256);
	      
	      _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	      _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	      _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	      _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	      _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	      _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	      _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	      _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     	      	     
	    }
	  else
	    {
	      _mm_store_pd(&x3[0], values[0]);	   
	      _mm_store_pd(&x3[2], values[1]);
	      _mm_store_pd(&x3[4], values[2]);
	      _mm_store_pd(&x3[6], values[3]);
	      _mm_store_pd(&x3[8], values[4]);	   
	      _mm_store_pd(&x3[10], values[5]);
	      _mm_store_pd(&x3[12], values[6]);
	      _mm_store_pd(&x3[14], values[7]);
	    }
	}		       	
      	
	x3 = x3_start;

	for (i = 0; i < n; i++)
	   {
	     if((x3_gap[i / 32] & mask32[i % 32]))
	       {	       
		 if(scaleGap)
		   {
		     if(useFastScaling)
		       addScale += wgt[i];
		     else
		       ex3[i]  += 1;
		   }
	       }
	     else
	       {				 
		 __m128d maxv =_mm_setzero_pd();		 
		 
		 if(x2_gap[i / 32] & mask32[i % 32])
		   x2 = x2_gapColumn;
		 else
		   {
		     x2 = x2_ptr;
		     x2_ptr += 16;
		   }
		 		 		 
		 uX1 = &umpX1[16 * tipX1[i]];	     
		 
		 
		 for (j = 0; j < 4; j++)
		   {		     		   
		     double *x2_p = &x2[j*4];
		     double *right_k0_p = &right[j*16];
		     double *right_k1_p = &right[j*16 + 1*4];
		     double *right_k2_p = &right[j*16 + 2*4];
		     double *right_k3_p = &right[j*16 + 3*4];
		     __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		     __m128d x2_2 = _mm_load_pd( &x2_p[2] );
		     
		     __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		     __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		     __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		     __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		     __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		     __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		     __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		     __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
		     
		     		     
		     right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		     right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
		     
		     right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		     right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
		     
		     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		     right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		     right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
		     
		     
		     right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		     right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
		     
		     
		     right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		     right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
		     
		     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		     right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		     right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);
		     
		     {
		       //
		       // load left side from tip vector
		       //
		       
		       __m128d uX1_k0_sse = _mm_load_pd( &uX1[j * 4] );
		       __m128d uX1_k2_sse = _mm_load_pd( &uX1[j * 4 + 2] );
		       
		       
		       //
		       // multiply left * right
			   //
		       
		       __m128d x1px2_k0 = _mm_mul_pd( uX1_k0_sse, right_k0_0 );
		       __m128d x1px2_k2 = _mm_mul_pd( uX1_k2_sse, right_k2_0 );
		       
		       
		       //
		       // multiply with EV matrix (!?)
		       //		   		  
		       
		       __m128d EV_t_l0_k0 = EVV[0];
		       __m128d EV_t_l0_k2 = EVV[1];
		       __m128d EV_t_l1_k0 = EVV[2];
		       __m128d EV_t_l1_k2 = EVV[3];
		       __m128d EV_t_l2_k0 = EVV[4];
		       __m128d EV_t_l2_k2 = EVV[5];
		       __m128d EV_t_l3_k0 = EVV[6]; 
		       __m128d EV_t_l3_k2 = EVV[7];
		       
		       
		       EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		       EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		       EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		       
		       EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		       EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		       
		       EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		       EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		       
		       EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		       EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		       EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		       
		       EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		       EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		       EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		       
		       EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		       
		       values[j * 2]     = EV_t_l0_k0;
		       values[j * 2 + 1] = EV_t_l2_k0;		   		   
			   
		       maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		       maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));		   
		     }		   
		   }

	     
		 _mm_store_pd(maxima, maxv);
		 
		 max = MAX(maxima[0], maxima[1]);
		 
		 if(max < minlikelihood)
		   {
		     __m128d sv = _mm_set1_pd(twotothe256);
		     
		     _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		     _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		     _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		     _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		     _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		     _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		     _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		     _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		     
		     if(useFastScaling)
		       addScale += wgt[i];
		     else
		       ex3[i]  += 1;
		   }
		 else
		   {
		     _mm_store_pd(&x3[0], values[0]);	   
		     _mm_store_pd(&x3[2], values[1]);
		     _mm_store_pd(&x3[4], values[2]);
		     _mm_store_pd(&x3[6], values[3]);
		     _mm_store_pd(&x3[8], values[4]);	   
		     _mm_store_pd(&x3[10], values[5]);
		     _mm_store_pd(&x3[12], values[6]);
		     _mm_store_pd(&x3[14], values[7]);
		   }		 
		 
		 x3 += 16;
	       }
	   }
      }
      break;
    case INNER_INNER:         
      {
	__m128d maxv =_mm_setzero_pd();
	
	scaleGap = 0;
	
	x1 = x1_gapColumn;	     	    
	x2 = x2_gapColumn;	    
	x3 = x3_gapColumn;
	
	for (j = 0; j < 4; j++)
	  {
	    
	    double *x1_p = &x1[j*4];
	    double *left_k0_p = &left[j*16];
	    double *left_k1_p = &left[j*16 + 1*4];
	    double *left_k2_p = &left[j*16 + 2*4];
	    double *left_k3_p = &left[j*16 + 3*4];
	    
	    __m128d x1_0 = _mm_load_pd( &x1_p[0] );
	    __m128d x1_2 = _mm_load_pd( &x1_p[2] );
	    
	    __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
	    __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
	    __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
	    __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
	    __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
	    __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
	    __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
	    __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
	    
	    left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
	    left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
	    
	    left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
	    left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
	    
	    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
	    left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
	    left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
	    
	    left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
	    left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
	    
	    left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
	    left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
	    
	    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
	    left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
	    left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
	    
	    
	    double *x2_p = &x2[j*4];
	    double *right_k0_p = &right[j*16];
	    double *right_k1_p = &right[j*16 + 1*4];
	    double *right_k2_p = &right[j*16 + 2*4];
	    double *right_k3_p = &right[j*16 + 3*4];
	    __m128d x2_0 = _mm_load_pd( &x2_p[0] );
	    __m128d x2_2 = _mm_load_pd( &x2_p[2] );
	    
	    __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
	    __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
	    __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
	    __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
	    __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
	    __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
	    __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
	    __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
	    
	    right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
	    right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
	    
	    right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
	    right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
	    
	    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
	    right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
	    right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
	    
	    right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
	    right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
	    	    
	    right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
	    right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
	    
	    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
	    right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
	    right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   		 		
	    
	    __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
	    __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );		 		 	   
	    
	    __m128d EV_t_l0_k0 = EVV[0];
	    __m128d EV_t_l0_k2 = EVV[1];
	    __m128d EV_t_l1_k0 = EVV[2];
	    __m128d EV_t_l1_k2 = EVV[3];
	    __m128d EV_t_l2_k0 = EVV[4];
	    __m128d EV_t_l2_k2 = EVV[5];
	    __m128d EV_t_l3_k0 = EVV[6]; 
	    __m128d EV_t_l3_k2 = EVV[7];
	    
	    EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
	    EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
	    
	    EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
	    EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
	    
	    EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
	    EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
	    
	    EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
	    EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
	    
	    EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
	    EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
	    EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
	    
	    EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
	    
	    
	    values[j * 2] = EV_t_l0_k0;
	    values[j * 2 + 1] = EV_t_l2_k0;            	   	    
	    
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
	    maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
	  }
		     
	_mm_store_pd(maxima, maxv);
	
	max = MAX(maxima[0], maxima[1]);
	
	if(max < minlikelihood)
	  {
	    __m128d sv = _mm_set1_pd(twotothe256);
	    
	    scaleGap = 1;
	    
	    _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
	    _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
	    _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
	    _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
	    _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
	    _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
	    _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
	    _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     	    	 
	  }
	else
	  {
	    _mm_store_pd(&x3[0], values[0]);	   
	    _mm_store_pd(&x3[2], values[1]);
	    _mm_store_pd(&x3[4], values[2]);
	    _mm_store_pd(&x3[6], values[3]);
	    _mm_store_pd(&x3[8], values[4]);	   
	    _mm_store_pd(&x3[10], values[5]);
	    _mm_store_pd(&x3[12], values[6]);
	    _mm_store_pd(&x3[14], values[7]);
	  }
      }

     
      x3 = x3_start;

     for (i = 0; i < n; i++)
       { 
	 if(x3_gap[i / 32] & mask32[i % 32])
	   {	     
	     if(scaleGap)
	       {
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1; 	       
	       }
	   }
	 else
	   {
	     __m128d maxv =_mm_setzero_pd();	     	    
	     
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
	     
	     
	     for (j = 0; j < 4; j++)
	       {
		 
		 double *x1_p = &x1[j*4];
		 double *left_k0_p = &left[j*16];
		 double *left_k1_p = &left[j*16 + 1*4];
		 double *left_k2_p = &left[j*16 + 2*4];
		 double *left_k3_p = &left[j*16 + 3*4];
		 
		 __m128d x1_0 = _mm_load_pd( &x1_p[0] );
		 __m128d x1_2 = _mm_load_pd( &x1_p[2] );
		 
		 __m128d left_k0_0 = _mm_load_pd( &left_k0_p[0] );
		 __m128d left_k0_2 = _mm_load_pd( &left_k0_p[2] );
		 __m128d left_k1_0 = _mm_load_pd( &left_k1_p[0] );
		 __m128d left_k1_2 = _mm_load_pd( &left_k1_p[2] );
		 __m128d left_k2_0 = _mm_load_pd( &left_k2_p[0] );
		 __m128d left_k2_2 = _mm_load_pd( &left_k2_p[2] );
		 __m128d left_k3_0 = _mm_load_pd( &left_k3_p[0] );
		 __m128d left_k3_2 = _mm_load_pd( &left_k3_p[2] );
		 
		 left_k0_0 = _mm_mul_pd(x1_0, left_k0_0);
		 left_k0_2 = _mm_mul_pd(x1_2, left_k0_2);
		 
		 left_k1_0 = _mm_mul_pd(x1_0, left_k1_0);
		 left_k1_2 = _mm_mul_pd(x1_2, left_k1_2);
		 
		 left_k0_0 = _mm_hadd_pd( left_k0_0, left_k0_2 );
		 left_k1_0 = _mm_hadd_pd( left_k1_0, left_k1_2);
		 left_k0_0 = _mm_hadd_pd( left_k0_0, left_k1_0);
		 
		 left_k2_0 = _mm_mul_pd(x1_0, left_k2_0);
		 left_k2_2 = _mm_mul_pd(x1_2, left_k2_2);
		 
		 left_k3_0 = _mm_mul_pd(x1_0, left_k3_0);
		 left_k3_2 = _mm_mul_pd(x1_2, left_k3_2);
		 
		 left_k2_0 = _mm_hadd_pd( left_k2_0, left_k2_2);
		 left_k3_0 = _mm_hadd_pd( left_k3_0, left_k3_2);
		 left_k2_0 = _mm_hadd_pd( left_k2_0, left_k3_0);
		 
		 
		 //
		 // multiply/add right side
		 //
		 double *x2_p = &x2[j*4];
		 double *right_k0_p = &right[j*16];
		 double *right_k1_p = &right[j*16 + 1*4];
		 double *right_k2_p = &right[j*16 + 2*4];
		 double *right_k3_p = &right[j*16 + 3*4];
		 __m128d x2_0 = _mm_load_pd( &x2_p[0] );
		 __m128d x2_2 = _mm_load_pd( &x2_p[2] );
		 
		 __m128d right_k0_0 = _mm_load_pd( &right_k0_p[0] );
		 __m128d right_k0_2 = _mm_load_pd( &right_k0_p[2] );
		 __m128d right_k1_0 = _mm_load_pd( &right_k1_p[0] );
		 __m128d right_k1_2 = _mm_load_pd( &right_k1_p[2] );
		 __m128d right_k2_0 = _mm_load_pd( &right_k2_p[0] );
		 __m128d right_k2_2 = _mm_load_pd( &right_k2_p[2] );
		 __m128d right_k3_0 = _mm_load_pd( &right_k3_p[0] );
		 __m128d right_k3_2 = _mm_load_pd( &right_k3_p[2] );
		 
		 right_k0_0 = _mm_mul_pd( x2_0, right_k0_0);
		 right_k0_2 = _mm_mul_pd( x2_2, right_k0_2);
		 
		 right_k1_0 = _mm_mul_pd( x2_0, right_k1_0);
		 right_k1_2 = _mm_mul_pd( x2_2, right_k1_2);
		 
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k0_2);
		 right_k1_0 = _mm_hadd_pd( right_k1_0, right_k1_2);
		 right_k0_0 = _mm_hadd_pd( right_k0_0, right_k1_0);
		 
		 right_k2_0 = _mm_mul_pd( x2_0, right_k2_0);
		 right_k2_2 = _mm_mul_pd( x2_2, right_k2_2);
		 
		 
		 right_k3_0 = _mm_mul_pd( x2_0, right_k3_0);
		 right_k3_2 = _mm_mul_pd( x2_2, right_k3_2);
		 
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k2_2);
		 right_k3_0 = _mm_hadd_pd( right_k3_0, right_k3_2);
		 right_k2_0 = _mm_hadd_pd( right_k2_0, right_k3_0);	   
		 
		 //
		 // multiply left * right
		 //
		 
		 __m128d x1px2_k0 = _mm_mul_pd( left_k0_0, right_k0_0 );
		 __m128d x1px2_k2 = _mm_mul_pd( left_k2_0, right_k2_0 );
		 
		 
		 //
		 // multiply with EV matrix (!?)
		 //	     
		 
		 __m128d EV_t_l0_k0 = EVV[0];
		 __m128d EV_t_l0_k2 = EVV[1];
		 __m128d EV_t_l1_k0 = EVV[2];
		 __m128d EV_t_l1_k2 = EVV[3];
		 __m128d EV_t_l2_k0 = EVV[4];
		 __m128d EV_t_l2_k2 = EVV[5];
		 __m128d EV_t_l3_k0 = EVV[6]; 
		 __m128d EV_t_l3_k2 = EVV[7];
		 
		 
		 EV_t_l0_k0 = _mm_mul_pd( x1px2_k0, EV_t_l0_k0 );
		 EV_t_l0_k2 = _mm_mul_pd( x1px2_k2, EV_t_l0_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l0_k2 );
		 
		 EV_t_l1_k0 = _mm_mul_pd( x1px2_k0, EV_t_l1_k0 );
		 EV_t_l1_k2 = _mm_mul_pd( x1px2_k2, EV_t_l1_k2 );
		 
		 EV_t_l1_k0 = _mm_hadd_pd( EV_t_l1_k0, EV_t_l1_k2 );
		 EV_t_l0_k0 = _mm_hadd_pd( EV_t_l0_k0, EV_t_l1_k0 );
		 
		 EV_t_l2_k0 = _mm_mul_pd( x1px2_k0, EV_t_l2_k0 );
		 EV_t_l2_k2 = _mm_mul_pd( x1px2_k2, EV_t_l2_k2 );
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l2_k2 );
		 
		 
		 EV_t_l3_k0 = _mm_mul_pd( x1px2_k0, EV_t_l3_k0 );
		 EV_t_l3_k2 = _mm_mul_pd( x1px2_k2, EV_t_l3_k2 );
		 EV_t_l3_k0 = _mm_hadd_pd( EV_t_l3_k0, EV_t_l3_k2 );
		 
		 EV_t_l2_k0 = _mm_hadd_pd( EV_t_l2_k0, EV_t_l3_k0 );
		 
		 
		 values[j * 2] = EV_t_l0_k0;
		 values[j * 2 + 1] = EV_t_l2_k0;            	   	    
		 
		 maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l0_k0, absMask.m));
		 maxv = _mm_max_pd(maxv, _mm_and_pd(EV_t_l2_k0, absMask.m));
	       }
	     
	     
	     _mm_store_pd(maxima, maxv);
	     
	     max = MAX(maxima[0], maxima[1]);
	     
	     if(max < minlikelihood)
	       {
		 __m128d sv = _mm_set1_pd(twotothe256);
		 
		 _mm_store_pd(&x3[0], _mm_mul_pd(values[0], sv));	   
		 _mm_store_pd(&x3[2], _mm_mul_pd(values[1], sv));
		 _mm_store_pd(&x3[4], _mm_mul_pd(values[2], sv));
		 _mm_store_pd(&x3[6], _mm_mul_pd(values[3], sv));
		 _mm_store_pd(&x3[8], _mm_mul_pd(values[4], sv));	   
		 _mm_store_pd(&x3[10], _mm_mul_pd(values[5], sv));
		 _mm_store_pd(&x3[12], _mm_mul_pd(values[6], sv));
		 _mm_store_pd(&x3[14], _mm_mul_pd(values[7], sv));	     
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;
	       }
	     else
	       {
		 _mm_store_pd(&x3[0], values[0]);	   
		 _mm_store_pd(&x3[2], values[1]);
		 _mm_store_pd(&x3[4], values[2]);
		 _mm_store_pd(&x3[6], values[3]);
		 _mm_store_pd(&x3[8], values[4]);	   
		 _mm_store_pd(&x3[10], values[5]);
		 _mm_store_pd(&x3[12], values[6]);
		 _mm_store_pd(&x3[14], values[7]);
	       }	 

	    
		 
	     x3 += 16;

	   }
       }
     break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;
}




#else

static void newviewGTRGAMMA_GAPPED(int tipCase,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *EV, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				   unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, int gapLength, 
				   double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn)
{
  int 
    i, 
    j, 
    k, 
    l, 
    scale, 
    addScale = 0, 
    scaleGap = 0;
  
  double
    *x1,
    *x2,
    *x3,
    buf,
    ump_x1,
    ump_x2,
    x1px2[4];

  for(i = 0; i < gapLength; i++)
    x3_gap[i] = x1_gap[i] & x2_gap[i];  

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256], *uX2, umpX2[256];


	for(i = 1; i < 16; i++)
	  {
	    x1 = &(tipVector[i * 4]);

	    for(j=0; j<4; j++)
	      for(k=0; k<4; k++)
		{
		  umpX1[i*16 + j*4 + k] = 0.0;
		  umpX2[i*16 + j*4 + k] = 0.0;

		  for (l=0; l < 4; l++)
		    {
		      umpX1[i*16 + j*4 + k] += x1[l] * left[j*16 + k*4 + l];
		      umpX2[i*16 + j*4 + k] += x1[l] * right[j*16 + k*4 + l];
		    }
		}
	  }

	
	x3 = x3_gapColumn;

	uX1 = &umpX1[240];
	uX2 = &umpX2[240];

	for(j = 0; j < 16; j++)
	  x3[j] = 0.0;

	for (j = 0; j < 4; j++)
	  for (k = 0; k < 4; k++)
	    {
	      buf = uX1[j*4 + k] * uX2[j*4 + k];
	      
	      for (l=0; l<4; l++)
		x3[j * 4 + l] +=  buf * EV[4 * k + l];
	    }	         

	for (i = 0; i < n; i++)
	  { 
	    if(!(x3_gap[i / 32] & mask32[i % 32]))
	      {
		x3 = &x3_start[i * 16];
		
		uX1 = &umpX1[16 * tipX1[i]];
		uX2 = &umpX2[16 * tipX2[i]];
		
		for(j = 0; j < 16; j++)
		  x3[j] = 0.0;
		
		for (j = 0; j < 4; j++)
		  for (k = 0; k < 4; k++)
		    {
		      buf = uX1[j*4 + k] * uX2[j*4 + k];
		      
		      for (l=0; l<4; l++)
			x3[j * 4 + l] +=  buf * EV[4 * k + l];
		    }	   
	      }
	  }
      }
      break;
    case TIP_INNER:
      {	
	double *uX1, umpX1[256];


	for (i = 1; i < 16; i++)
	  {
	    x1 = &(tipVector[i*4]);

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  umpX1[i*16 + j*4 + k] = 0.0;
		  for (l=0; l < 4; l++)
		    umpX1[i*16 + j*4 + k] += x1[l] * left[j*16 + k*4 + l];
		}
	  }

	x2 = x2_gapColumn;
	x3 = x3_gapColumn;

	uX1 = &umpX1[240];

	for(j = 0; j < 16; j++)
	  x3[j] = 0.0;

	for (j = 0; j < 4; j++)
	  {
	    for (k = 0; k < 4; k++)
	      {
		ump_x2 = 0.0;
		
		for (l=0; l<4; l++)
		  ump_x2 += x2[j*4 + l] * right[j* 16 + k*4 + l];
		x1px2[k] = uX1[j * 4 + k] * ump_x2;
	      }
	    
	    for(k = 0; k < 4; k++)
	      for (l=0; l<4; l++)
		x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];	    
	  }	     

	scale = 1;
	for(l = 0; scale && (l < 16); l++)
	  scale = (ABS(x3[l]) <  minlikelihood);

	if(scale)
	  {
	    scaleGap = 1;
	    for (l=0; l<16; l++)
	      x3[l] *= twotothe256;	    		 
	  }
	
	for (i = 0; i < n; i++)
	  {
	    if((x3_gap[i / 32] & mask32[i % 32]))
	      {	       
		if(scaleGap)
		  {
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;
		  }
	      }
	    else
	      {
		if(x2_gap[i / 32] & mask32[i % 32])
		  x2 = x2_gapColumn;
		else
		  x2 = &x2_start[i * 16];
		
		x3 = &x3_start[i * 16];

		uX1 = &umpX1[16 * tipX1[i]];
		
		for(j = 0; j < 16; j++)
		  x3[j] = 0.0;
		
		for (j = 0; j < 4; j++)
		  {
		    
		    
		    for (k = 0; k < 4; k++)
		      {
			ump_x2 = 0.0;
			
			for (l=0; l<4; l++)
			  ump_x2 += x2[j*4 + l] * right[j* 16 + k*4 + l];
			x1px2[k] = uX1[j * 4 + k] * ump_x2;
		      }
		    
		    for(k = 0; k < 4; k++)
		      for (l=0; l<4; l++)
			x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];
		    
		    
		  }	     
		
		scale = 1;
		for(l = 0; scale && (l < 16); l++)
		  scale = (ABS(x3[l]) <  minlikelihood);
		
		if(scale)
		  {
		    for (l=0; l<16; l++)
		      x3[l] *= twotothe256;
		    
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;		 
		  }
		
	      }
	  }
      }
      break;
    case INNER_INNER:

      x1 = x1_gapColumn;
      x2 = x2_gapColumn;
      x3 = x3_gapColumn;

      for(j = 0; j < 16; j++)
	x3[j] = 0.0;

      for (j = 0; j < 4; j++)
	{
	  for (k = 0; k < 4; k++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;
	      
	      for (l=0; l<4; l++)
		{
		  ump_x1 += x1[j*4 + l] * left[j*16 + k*4 +l];
		  ump_x2 += x2[j*4 + l] * right[j*16 + k*4 +l];
		}

	      x1px2[k] = ump_x1 * ump_x2;
	    }

	  for(k = 0; k < 4; k++)
	    for (l=0; l<4; l++)
	      x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];
	}	 

      scale = 1;
      for(l = 0; scale && (l < 16); l++)
	scale = (ABS(x3[l]) <  minlikelihood);

      if(scale)
	{
	  scaleGap = 1;
	  
	  for (l=0; l<16; l++)
	    x3[l] *= twotothe256;	 
	}		      	

     for (i = 0; i < n; i++)
       {
	 if((x3_gap[i / 32] & mask32[i % 32]))
	   {	       
	     if(scaleGap)
	       {
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;
	       }
	   }
	 else
	   {
	     if(x1_gap[i / 32] & mask32[i % 32])
	       x1 = x1_gapColumn;
	     else
	       x1 = &x1_start[i * 16];

	     if(x2_gap[i / 32] & mask32[i % 32])
	       x2 = x2_gapColumn;
	     else
	       x2 = &x2_start[i * 16];

	     x3 = &x3_start[i * 16];	 
	     
	     for(j = 0; j < 16; j++)
	       x3[j] = 0.0;
	     
	     for (j = 0; j < 4; j++)
	       {
		 for (k = 0; k < 4; k++)
		   {
		     ump_x1 = 0.0;
		     ump_x2 = 0.0;
		     
		     for (l=0; l<4; l++)
		       {
			 ump_x1 += x1[j*4 + l] * left[j*16 + k*4 +l];
			 ump_x2 += x2[j*4 + l] * right[j*16 + k*4 +l];
		       }		     		     		     

		     x1px2[k] = ump_x1 * ump_x2;
		   }
		 
		 for(k = 0; k < 4; k++)
		   for (l=0; l<4; l++)
		     x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];
	       }	 
	     
	     scale = 1;
	     for(l = 0; scale && (l < 16); l++)
	       scale = (ABS(x3[l]) <  minlikelihood);
	     
	     if(scale)
	       {
		 for (l=0; l<16; l++)
		   x3[l] *= twotothe256;
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;	   
	       }		      		
	   }
       }
     break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;

}
  



static void newviewGTRGAMMA(int tipCase,
			    double *x1_start, double *x2_start, double *x3_start,
			    double *EV, double *tipVector,
			    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			    const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
			    )
{
  int i, j, k, l, scale, addScale = 0;
  double
    *x1,
    *x2,
    *x3,
    buf,
    ump_x1,
    ump_x2;


  double x1px2[4];

  
 


  switch(tipCase)
    {
    case TIP_TIP:
      {
	double *uX1, umpX1[256], *uX2, umpX2[256];


	for(i = 1; i < 16; i++)
	  {
	    x1 = &(tipVector[i * 4]);

	    for(j=0; j<4; j++)
	      for(k=0; k<4; k++)
		{
		  umpX1[i*16 + j*4 + k] = 0.0;
		  umpX2[i*16 + j*4 + k] = 0.0;

		  for (l=0; l < 4; l++)
		    {
		      umpX1[i*16 + j*4 + k] += x1[l] * left[j*16 + k*4 + l];
		      umpX2[i*16 + j*4 + k] += x1[l] * right[j*16 + k*4 + l];
		    }
		}
	  }

	

	for (i = 0; i < n; i++)
	  {
	    x3 = &x3_start[i * 16];

	    uX1 = &umpX1[16 * tipX1[i]];
	    uX2 = &umpX2[16 * tipX2[i]];

	    for(j = 0; j < 16; j++)
	      x3[j] = 0.0;

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  buf = uX1[j*4 + k] * uX2[j*4 + k];

		  for (l=0; l<4; l++)
		    x3[j * 4 + l] +=  buf * EV[4 * k + l];
		}	   
	  }
      }
      break;
    case TIP_INNER:
      {	
	double *uX1, umpX1[256];


	for (i = 1; i < 16; i++)
	  {
	    x1 = &(tipVector[i*4]);

	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{
		  umpX1[i*16 + j*4 + k] = 0.0;
		  for (l=0; l < 4; l++)
		    umpX1[i*16 + j*4 + k] += x1[l] * left[j*16 + k*4 + l];
		}
	  }


	 for (i = 0; i < n; i++)
	   {
	     x2 = &x2_start[i * 16];
	     x3 = &x3_start[i * 16];

	     uX1 = &umpX1[16 * tipX1[i]];

	     for(j = 0; j < 16; j++)
	       x3[j] = 0.0;

	     for (j = 0; j < 4; j++)
	       {


		 for (k = 0; k < 4; k++)
		   {
		     ump_x2 = 0.0;

		     for (l=0; l<4; l++)
		       ump_x2 += x2[j*4 + l] * right[j* 16 + k*4 + l];
		     x1px2[k] = uX1[j * 4 + k] * ump_x2;
		   }

		 for(k = 0; k < 4; k++)
		   for (l=0; l<4; l++)
		     x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];


	       }	     

	     scale = 1;
	     for(l = 0; scale && (l < 16); l++)
	       scale = (ABS(x3[l]) <  minlikelihood);

	     if(scale)
	       {
		 for (l=0; l<16; l++)
		   x3[l] *= twotothe256;

		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;		 
	       }

	   }
      }
      break;
    case INNER_INNER:
     for (i = 0; i < n; i++)
       {
	 x1 = &x1_start[i * 16];
	 x2 = &x2_start[i * 16];
	 x3 = &x3_start[i * 16];

	 for(j = 0; j < 16; j++)
	   x3[j] = 0.0;

	 for (j = 0; j < 4; j++)
	   {
	     for (k = 0; k < 4; k++)
	       {
		 ump_x1 = 0.0;
		 ump_x2 = 0.0;

		 for (l=0; l<4; l++)
		   {
		     ump_x1 += x1[j*4 + l] * left[j*16 + k*4 +l];
		     ump_x2 += x2[j*4 + l] * right[j*16 + k*4 +l];
		   }




		 x1px2[k] = ump_x1 * ump_x2;
	       }

	     for(k = 0; k < 4; k++)
	       for (l=0; l<4; l++)
	         x3[j * 4 + l] +=  x1px2[k] * EV[4 * k + l];
           }	 

	 scale = 1;
	 for(l = 0; scale && (l < 16); l++)
	   scale = (ABS(x3[l]) <  minlikelihood);

	 if(scale)
	   {
	     for (l=0; l<16; l++)
	       x3[l] *= twotothe256;

	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	   
	   }		      		
       }
   
     break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;

}

#endif


#ifdef __SIM_SSE3
static void newviewGTRCATPROT_SAVE(int tipCase, double *extEV,
    int *cptr,
    double *x1, double *x2, double *x3, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
    unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,
    double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn, const int maxCats)
{
  double
    *le, 
  *ri, 
  *v, 
  *vl, 
  *vr,
  *x1_ptr = x1,
  *x2_ptr = x2, 
  *x3_ptr = x3;

  int 
    i, 
    l, 
    j, 
    scale, 
    scaleGap = 0,
    addScale = 0;

  {
    vl = x1_gapColumn;	      
    vr = x2_gapColumn;
    v = x3_gapColumn;

    le = &left[maxCats * 400];
    ri = &right[maxCats * 400];	  

    for(l = 0; l < 20; l+=2)
      _mm_store_pd(&v[l], _mm_setzero_pd());	      		

    for(l = 0; l < 20; l++)
    {
      __m128d x1v = _mm_setzero_pd();
      __m128d x2v = _mm_setzero_pd();
      double 
        *ev = &extEV[l * 20],
        *lv = &le[l * 20],
        *rv = &ri[l * 20];


      for(j = 0; j < 20; j+=2)
      {
        x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
        x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
      }

      x1v = _mm_hadd_pd(x1v, x1v);
      x2v = _mm_hadd_pd(x2v, x2v);

      x1v = _mm_mul_pd(x1v, x2v);

      for(j = 0; j < 20; j+=2)
      {
        __m128d vv = _mm_load_pd(&v[j]);
        vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
        _mm_store_pd(&v[j], vv);
      }		    	
    }

    if(tipCase != TIP_TIP)
    { 	    
      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

      scale = 1;
      for(l = 0; scale && (l < 20); l += 2)
      {
        __m128d vv = _mm_load_pd(&v[l]);
        __m128d v1 = _mm_and_pd(vv, absMask.m);
        v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
        if(_mm_movemask_pd( v1 ) != 3)
          scale = 0;
      }	    	        

      if(scale)
      {
        __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

        for(l = 0; l < 20; l+=2)
        {
          __m128d ex3v = _mm_load_pd(&v[l]);		  
          _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
        }		   		  

        scaleGap = TRUE;	   
      }
    }
  }

  switch(tipCase)
  {
    case TIP_TIP:
      {
        for (i = 0; i < n; i++)
        {
          if(noGap(x3_gap, i))
          {		
            vl = &(tipVector[20 * tipX1[i]]);
            vr = &(tipVector[20 * tipX2[i]]);
            v  = x3_ptr;

            if(isGap(x1_gap, i))
              le =  &left[maxCats * 400];
            else	  	  
              le =  &left[cptr[i] * 400];	  

            if(isGap(x2_gap, i))
              ri =  &right[maxCats * 400];
            else	 	  
              ri =  &right[cptr[i] * 400];

            for(l = 0; l < 20; l+=2)
              _mm_store_pd(&v[l], _mm_setzero_pd());	      		

            for(l = 0; l < 20; l++)
            {
              __m128d x1v = _mm_setzero_pd();
              __m128d x2v = _mm_setzero_pd();	 
              double 
                *ev = &extEV[l * 20],
                *lv = &le[l * 20],
                *rv = &ri[l * 20];

              for(j = 0; j < 20; j+=2)
              {
                x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
                x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
              }

              x1v = _mm_hadd_pd(x1v, x1v);
              x2v = _mm_hadd_pd(x2v, x2v);

              x1v = _mm_mul_pd(x1v, x2v);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
                _mm_store_pd(&v[j], vv);
              }		   
            }

            x3_ptr += 20;

          }   
        }
      }
      break;
    case TIP_INNER:
      {
        for (i = 0; i < n; i++)
        {
          if(isGap(x3_gap, i))
          {
            if(scaleGap)		   		    
              addScale += wgt[i];
          }
          else
          {	 
            vl = &(tipVector[20 * tipX1[i]]);

            vr = x2_ptr;
            v = x3_ptr;

            if(isGap(x1_gap, i))
              le =  &left[maxCats * 400];
            else
              le =  &left[cptr[i] * 400];

            if(isGap(x2_gap, i))
            {		 
              ri =  &right[maxCats * 400];
              vr = x2_gapColumn;
            }
            else
            {
              ri =  &right[cptr[i] * 400];
              vr = x2_ptr;
              x2_ptr += 20;
            }	  	  	  	  		  

            for(l = 0; l < 20; l+=2)
              _mm_store_pd(&v[l], _mm_setzero_pd());	      			   

            for(l = 0; l < 20; l++)
            {
              __m128d x1v = _mm_setzero_pd();
              __m128d x2v = _mm_setzero_pd();	
              double 
                *ev = &extEV[l * 20],
                *lv = &le[l * 20],
                *rv = &ri[l * 20];

              for(j = 0; j < 20; j+=2)
              {
                x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
                x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
              }

              x1v = _mm_hadd_pd(x1v, x1v);
              x2v = _mm_hadd_pd(x2v, x2v);

              x1v = _mm_mul_pd(x1v, x2v);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
                _mm_store_pd(&v[j], vv);
              }		    
            }

            { 	    
              __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

              scale = 1;
              for(l = 0; scale && (l < 20); l += 2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d v1 = _mm_and_pd(vv, absMask.m);
                v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
                if(_mm_movemask_pd( v1 ) != 3)
                  scale = 0;
              }	    	  
            }


            if(scale)
            {
              __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

              for(l = 0; l < 20; l+=2)
              {
                __m128d ex3v = _mm_load_pd(&v[l]);
                _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));		    
              }

              addScale += wgt[i];	  
            }
            x3_ptr += 20;
          }
        }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
      { 
        if(isGap(x3_gap, i))
        {
          if(scaleGap)		   		    
            addScale += wgt[i];
        }
        else
        {	  	     
          v = x3_ptr;

          if(isGap(x1_gap, i))
          {
            vl = x1_gapColumn;
            le =  &left[maxCats * 400];
          }
          else
          {
            le =  &left[cptr[i] * 400];
            vl = x1_ptr;
            x1_ptr += 20;
          }

          if(isGap(x2_gap, i))	
          {
            vr = x2_gapColumn;
            ri =  &right[maxCats * 400];	    
          }
          else
          {
            ri =  &right[cptr[i] * 400];
            vr = x2_ptr;
            x2_ptr += 20;
          }	 	  	  	  

          for(l = 0; l < 20; l+=2)
            _mm_store_pd(&v[l], _mm_setzero_pd());	      		

          for(l = 0; l < 20; l++)
          {
            __m128d x1v = _mm_setzero_pd();
            __m128d x2v = _mm_setzero_pd();
            double 
              *ev = &extEV[l * 20],
              *lv = &le[l * 20],
              *rv = &ri[l * 20];

            for(j = 0; j < 20; j+=2)
            {
              x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
              x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
            }

            x1v = _mm_hadd_pd(x1v, x1v);
            x2v = _mm_hadd_pd(x2v, x2v);

            x1v = _mm_mul_pd(x1v, x2v);

            for(j = 0; j < 20; j+=2)
            {
              __m128d vv = _mm_load_pd(&v[j]);
              vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
              _mm_store_pd(&v[j], vv);
            }		    

          }

          { 	    
            __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

            scale = 1;
            for(l = 0; scale && (l < 20); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }	    	  
          }

          if(scale)
          {
            __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

            for(l = 0; l < 20; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);		  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
            }		   		  

            addScale += wgt[i];	   
          }
          x3_ptr += 20;
        }
      }
      break;
    default:
      assert(0);
  }


  *scalerIncrement = addScale;

}


#endif

static void newviewGTRCATPROT(int tipCase, double *extEV,
			      int *cptr,
			      double *x1, double *x2, double *x3, double *tipVector,
			      int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			      int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;
#ifndef __SIM_SSE3
  double
    ump_x1, ump_x2, x1px2;
#endif
  int i, l, j, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 400];
	    ri = &right[cptr[i] * 400];

	    vl = &(tipVector[20 * tipX1[i]]);
	    vr = &(tipVector[20 * tipX2[i]]);
	    v  = &x3[20 * i];
#ifdef __SIM_SSE3
	    for(l = 0; l < 20; l+=2)
	      _mm_store_pd(&v[l], _mm_setzero_pd());	      		
#else
	    for(l = 0; l < 20; l++)
	      v[l] = 0.0;
#endif

	    for(l = 0; l < 20; l++)
	      {
#ifdef __SIM_SSE3
		__m128d x1v = _mm_setzero_pd();
		__m128d x2v = _mm_setzero_pd();	 
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];

		for(j = 0; j < 20; j+=2)
		  {
		    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		  }

		x1v = _mm_hadd_pd(x1v, x1v);
		x2v = _mm_hadd_pd(x2v, x2v);

		x1v = _mm_mul_pd(x1v, x2v);
		
		for(j = 0; j < 20; j+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[j]);
		    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		    _mm_store_pd(&v[j], vv);
		  }		    
#else
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 20; j++)
		  {
		    ump_x1 += vl[j] * le[l * 20 + j];
		    ump_x2 += vr[j] * ri[l * 20 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 20; j++)
		  v[j] += x1px2 * extEV[l * 20 + j];
#endif
	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 400];
	    ri = &right[cptr[i] * 400];

	    vl = &(tipVector[20 * tipX1[i]]);
	    vr = &x2[20 * i];
	    v  = &x3[20 * i];

#ifdef __SIM_SSE3
	    for(l = 0; l < 20; l+=2)
	      _mm_store_pd(&v[l], _mm_setzero_pd());	      		
#else
	    for(l = 0; l < 20; l++)
	      v[l] = 0.0;
#endif
	   

	    for(l = 0; l < 20; l++)
	      {
#ifdef __SIM_SSE3

		__m128d x1v = _mm_setzero_pd();
		__m128d x2v = _mm_setzero_pd();	
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];

		for(j = 0; j < 20; j+=2)
		  {
		    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		  }

		x1v = _mm_hadd_pd(x1v, x1v);
		x2v = _mm_hadd_pd(x2v, x2v);

		x1v = _mm_mul_pd(x1v, x2v);
		
		for(j = 0; j < 20; j+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[j]);
		    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		    _mm_store_pd(&v[j], vv);
		  }		    
#else
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 20; j++)
		  {
		    ump_x1 += vl[j] * le[l * 20 + j];
		    ump_x2 += vr[j] * ri[l * 20 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 20; j++)
		  v[j] += x1px2 * extEV[l * 20 + j];
#endif
	      }
#ifdef __SIM_SSE3
	    { 	    
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 20); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }
#else
	    scale = 1;
	    for(l = 0; scale && (l < 20); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	   
#endif 

	    if(scale)
	      {
#ifdef __SIM_SSE3
		__m128d twoto = _mm_set_pd(twotothe256, twotothe256);

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d ex3v = _mm_load_pd(&v[l]);
		    _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));		    
		  }
#else
		for(l = 0; l < 20; l++)
		  v[l] *= twotothe256;
#endif

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 400];
	  ri = &right[cptr[i] * 400];

	  vl = &x1[20 * i];
	  vr = &x2[20 * i];
	  v = &x3[20 * i];

#ifdef __SIM_SSE3
	    for(l = 0; l < 20; l+=2)
	      _mm_store_pd(&v[l], _mm_setzero_pd());	      		
#else
	    for(l = 0; l < 20; l++)
	      v[l] = 0.0;
#endif
	 
	  for(l = 0; l < 20; l++)
	    {
#ifdef __SIM_SSE3
		__m128d x1v = _mm_setzero_pd();
		__m128d x2v = _mm_setzero_pd();
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];


		for(j = 0; j < 20; j+=2)
		  {
		    x1v = _mm_add_pd(x1v, _mm_mul_pd(_mm_load_pd(&vl[j]), _mm_load_pd(&lv[j])));		    
		    x2v = _mm_add_pd(x2v, _mm_mul_pd(_mm_load_pd(&vr[j]), _mm_load_pd(&rv[j])));
		  }

		x1v = _mm_hadd_pd(x1v, x1v);
		x2v = _mm_hadd_pd(x2v, x2v);

		x1v = _mm_mul_pd(x1v, x2v);
		
		for(j = 0; j < 20; j+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[j]);
		    vv = _mm_add_pd(vv, _mm_mul_pd(x1v, _mm_load_pd(&ev[j])));
		    _mm_store_pd(&v[j], vv);
		  }		    
#else
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 20; j++)
		{
		  ump_x1 += vl[j] * le[l * 20 + j];
		  ump_x2 += vr[j] * ri[l * 20 + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 20; j++)
		v[j] += x1px2 * extEV[l * 20 + j];
#endif
	    }
#ifdef __SIM_SSE3
	    { 	    
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 20); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }
#else
	   scale = 1;
	   for(l = 0; scale && (l < 20); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
#endif	   

	   if(scale)
	     {
#ifdef __SIM_SSE3
	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 20; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  
#else
	       for(l = 0; l < 20; l++)
		 v[l] *= twotothe256;
#endif

	       if(useFastScaling)
		 addScale += wgt[i];
	       else
		 ex3[i]  += 1;	      
	     }
	}
      break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;

}






static void newviewGTRCATSECONDARY(int tipCase, double *extEV,
				   int *cptr,
				   double *x1, double *x2, double *x3, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;
  double
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 256];
	    ri = &right[cptr[i] * 256];

	    vl = &(tipVector[16 * tipX1[i]]);
	    vr = &(tipVector[16 * tipX2[i]]);
	    v  = &x3[16 * i];

	    for(l = 0; l < 16; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 16; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 16; j++)
		  {
		    ump_x1 += vl[j] * le[l * 16 + j];
		    ump_x2 += vr[j] * ri[l * 16 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 16; j++)
		  v[j] += x1px2 * extEV[l * 16 + j];
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 256];
	    ri = &right[cptr[i] * 256];

	    vl = &(tipVector[16 * tipX1[i]]);
	    vr = &x2[16 * i];
	    v  = &x3[16 * i];

	    for(l = 0; l < 16; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 16; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 16; j++)
		  {
		    ump_x1 += vl[j] * le[l * 16 + j];
		    ump_x2 += vr[j] * ri[l * 16 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 16; j++)
		  v[j] += x1px2 * extEV[l * 16 + j];
	      }

	    scale = 1;
	    for(l = 0; scale && (l < 16); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	    

	    if(scale)
	      {
		for(l = 0; l < 16; l++)
		  v[l] *= twotothe256;
		
		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 256];
	  ri = &right[cptr[i] * 256];

	  vl = &x1[16 * i];
	  vr = &x2[16 * i];
	  v = &x3[16 * i];

	  for(l = 0; l < 16; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 16; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 16; j++)
		{
		  ump_x1 += vl[j] * le[l * 16 + j];
		  ump_x2 += vr[j] * ri[l * 16 + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 16; j++)
		v[j] += x1px2 * extEV[l * 16 + j];
	    }

	   scale = 1;
	   for(l = 0; scale && (l < 16); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	  
	   if(scale)
	     {
	       for(l = 0; l < 16; l++)
		 v[l] *= twotothe256;

	       if(useFastScaling)
		 addScale += wgt[i];
	       else
		 ex3[i]  += 1;	     
	     }
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}



static void newviewGTRCATSECONDARY_6(int tipCase, double *extEV,
				   int *cptr,
				   double *x1, double *x2, double *x3, double *tipVector,
				   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				   int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;
  double
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 36];
	    ri = &right[cptr[i] * 36];

	    vl = &(tipVector[6 * tipX1[i]]);
	    vr = &(tipVector[6 * tipX2[i]]);
	    v  = &x3[6 * i];

	    for(l = 0; l < 6; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 6; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 6; j++)
		  {
		    ump_x1 += vl[j] * le[l * 6 + j];
		    ump_x2 += vr[j] * ri[l * 6 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 6; j++)
		  v[j] += x1px2 * extEV[l * 6 + j];
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 36];
	    ri = &right[cptr[i] * 36];

	    vl = &(tipVector[6 * tipX1[i]]);
	    vr = &x2[6 * i];
	    v  = &x3[6 * i];

	    for(l = 0; l < 6; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 6; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 6; j++)
		  {
		    ump_x1 += vl[j] * le[l * 6 + j];
		    ump_x2 += vr[j] * ri[l * 6 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 6; j++)
		  v[j] += x1px2 * extEV[l * 6 + j];
	      }

	    scale = 1;
	    for(l = 0; scale && (l < 6); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	   

	    if(scale)
	      {
		for(l = 0; l < 6; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	       
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 36];
	  ri = &right[cptr[i] * 36];

	  vl = &x1[6 * i];
	  vr = &x2[6 * i];
	  v = &x3[6 * i];

	  for(l = 0; l < 6; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 6; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 6; j++)
		{
		  ump_x1 += vl[j] * le[l * 6 + j];
		  ump_x2 += vr[j] * ri[l * 6 + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 6; j++)
		v[j] += x1px2 * extEV[l * 6 + j];
	    }

	   scale = 1;
	   for(l = 0; scale && (l < 6); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	  

	   if(scale)
	     {
	       for(l = 0; l < 6; l++)
		 v[l] *= twotothe256;
	       
	       if(useFastScaling)
		 addScale += wgt[i];
	       else
		 ex3[i]  += 1;
	     }
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewGTRCATSECONDARY_7(int tipCase, double *extEV,
				     int *cptr,
				     double *x1, double *x2, double *x3, double *tipVector,
				     int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				     int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;
  double
    ump_x1, ump_x2, x1px2;
  int i, l, j, scale, addScale = 0;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 49];
	    ri = &right[cptr[i] * 49];

	    vl = &(tipVector[7 * tipX1[i]]);
	    vr = &(tipVector[7 * tipX2[i]]);
	    v  = &x3[7 * i];

	    for(l = 0; l < 7; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 7; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 7; j++)
		  {
		    ump_x1 += vl[j] * le[l * 7 + j];
		    ump_x2 += vr[j] * ri[l * 7 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 7; j++)
		  v[j] += x1px2 * extEV[l * 7 + j];
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    le = &left[cptr[i] * 49];
	    ri = &right[cptr[i] * 49];

	    vl = &(tipVector[7 * tipX1[i]]);
	    vr = &x2[7 * i];
	    v  = &x3[7 * i];

	    for(l = 0; l < 7; l++)
	      v[l] = 0.0;

	    for(l = 0; l < 7; l++)
	      {
		ump_x1 = 0.0;
		ump_x2 = 0.0;

		for(j = 0; j < 7; j++)
		  {
		    ump_x1 += vl[j] * le[l * 7 + j];
		    ump_x2 += vr[j] * ri[l * 7 + j];
		  }

		x1px2 = ump_x1 * ump_x2;

		for(j = 0; j < 7; j++)
		  v[j] += x1px2 * extEV[l * 7 + j];
	      }

	    scale = 1;
	    for(l = 0; scale && (l < 7); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	    

	    if(scale)
	      {
		for(l = 0; l < 7; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		 addScale += wgt[i];
		else
		  ex3[i]  += 1;			     
	      }
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 49];
	  ri = &right[cptr[i] * 49];

	  vl = &x1[7 * i];
	  vr = &x2[7 * i];
	  v = &x3[7 * i];

	  for(l = 0; l < 7; l++)
	    v[l] = 0.0;

	  for(l = 0; l < 7; l++)
	    {
	      ump_x1 = 0.0;
	      ump_x2 = 0.0;

	      for(j = 0; j < 7; j++)
		{
		  ump_x1 += vl[j] * le[l * 7 + j];
		  ump_x2 += vr[j] * ri[l * 7 + j];
		}

	      x1px2 =  ump_x1 * ump_x2;

	      for(j = 0; j < 7; j++)
		v[j] += x1px2 * extEV[l * 7 + j];
	    }

	   scale = 1;
	   for(l = 0; scale && (l < 7); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	  

	   if(scale)
	     {
	       for(l = 0; l < 7; l++)
		 v[l] *= twotothe256;

	       if(useFastScaling)
		 addScale += wgt[i];
	       else
		 ex3[i]  += 1;	      
	     }
	}
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}



static void newviewGTRGAMMAPROT(int tipCase,
				double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double  *uX1, *uX2, *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr;
#ifndef __SIM_SSE3
  double al, ar;
#endif



  switch(tipCase)
    {
    case TIP_TIP:
      {
	double umpX1[1840], umpX2[1840];

	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
#ifdef __SIM_SSE3
		double *ll =  &left[k * 20];
		double *rr =  &right[k * 20];
		
		__m128d umpX1v = _mm_setzero_pd();
		__m128d umpX2v = _mm_setzero_pd();

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
		    umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
		umpX2v = _mm_hadd_pd(umpX2v, umpX2v);
		
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);
		_mm_storel_pd(&umpX2[80 * i + k], umpX2v);
#else
		umpX1[80 * i + k] = 0.0;
		umpX2[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  {
		    umpX1[80 * i + k] +=  v[l] *  left[k * 20 + l];
		    umpX2[80 * i + k] +=  v[l] * right[k * 20 + l];
		  }
#endif
	      }
	  }

	for(i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];

	    for(j = 0; j < 4; j++)
	      {
		v = &x3[i * 80 + j * 20];

#ifdef __SIM_SSE3
		__m128d zero =  _mm_setzero_pd();
		for(k = 0; k < 20; k+=2)		  		    
		  _mm_store_pd(&v[k], zero);

		for(k = 0; k < 20; k++)
		  { 
		    double *eev = &extEV[k * 20];
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		    __m128d x1px2v = _mm_set1_pd(x1px2);

		    for(l = 0; l < 20; l+=2)
		      {
		      	__m128d vv = _mm_load_pd(&v[l]);
			__m128d ee = _mm_load_pd(&eev[l]);

			vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			
			_mm_store_pd(&v[l], vv);
		      }
		  }

#else

		for(k = 0; k < 20; k++)
		  v[k] = 0.0;

		for(k = 0; k < 20; k++)
		  {		   
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		   
		    for(l = 0; l < 20; l++)		      					
		      v[l] += x1px2 * extEV[20 * k + l];		     
		  }
#endif
	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	double umpX1[1840], ump_x2[20];


	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
#ifdef __SIM_SSE3
		double *ll =  &left[k * 20];
				
		__m128d umpX1v = _mm_setzero_pd();
		
		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));		    					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);				
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);		
#else	    
		umpX1[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  umpX1[80 * i + k] +=  v[l] * left[k * 20 + l];
#endif

	      }
	  }

	for (i = 0; i < n; i++)
	  {
	    uX1 = &umpX1[80 * tipX1[i]];

	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[80 * i + k * 20]);
#ifdef __SIM_SSE3	       
		for(l = 0; l < 20; l++)
		  {		   
		    double *r =  &right[k * 400 + l * 20];
		    __m128d ump_x2v = _mm_setzero_pd();	    
		    
		    for(j = 0; j < 20; j+= 2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			__m128d rr = _mm_load_pd(&r[j]);
			ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
		      }
		     
		    ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
		    
		    _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		  }

		v = &(x3[80 * i + 20 * k]);

		__m128d zero =  _mm_setzero_pd();
		for(l = 0; l < 20; l+=2)		  		    
		  _mm_store_pd(&v[l], zero);
		  
		for(l = 0; l < 20; l++)
		  {
		    double *eev = &extEV[l * 20];
		    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		    __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		    for(j = 0; j < 20; j+=2)
		      {
			__m128d vv = _mm_load_pd(&v[j]);
			__m128d ee = _mm_load_pd(&eev[j]);
			
			vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			
			_mm_store_pd(&v[j], vv);
		      }		     		    
		  }			
#else
		for(l = 0; l < 20; l++)
		  {
		    ump_x2[l] = 0.0;

		    for(j = 0; j < 20; j++)
		      ump_x2[l] += v[j] * right[k * 400 + l * 20 + j];
		  }

		v = &(x3[80 * i + 20 * k]);

		for(l = 0; l < 20; l++)
		  v[l] = 0;

		for(l = 0; l < 20; l++)
		  {
		    x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		    for(j = 0; j < 20; j++)
		      v[j] += x1px2 * extEV[l * 20  + j];
		  }
#endif
	      }
	   
#ifdef __SIM_SSE3
	    { 
	      v = &(x3[80 * i]);
	      __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	      
	      scale = 1;
	      for(l = 0; scale && (l < 80); l += 2)
		{
		  __m128d vv = _mm_load_pd(&v[l]);
		  __m128d v1 = _mm_and_pd(vv, absMask.m);
		  v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		  if(_mm_movemask_pd( v1 ) != 3)
		    scale = 0;
		}	    	  
	    }
#else
	    v = &x3[80 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 80); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
#endif

	    if (scale)
	      {
#ifdef __SIM_SSE3
	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 80; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  
#else
		for(l = 0; l < 80; l++)
		  v[l] *= twotothe256;
#endif

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	       
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[80 * i + 20 * k]);
	     vr = &(x2[80 * i + 20 * k]);
	     v =  &(x3[80 * i + 20 * k]);

#ifdef __SIM_SSE3
	     __m128d zero =  _mm_setzero_pd();
	     for(l = 0; l < 20; l+=2)		  		    
	       _mm_store_pd(&v[l], zero);
#else
	     for(l = 0; l < 20; l++)
	       v[l] = 0;
#endif

	     for(l = 0; l < 20; l++)
	       {		 
#ifdef __SIM_SSE3
		 {
		   __m128d al = _mm_setzero_pd();
		   __m128d ar = _mm_setzero_pd();

		   double *ll   = &left[k * 400 + l * 20];
		   double *rr   = &right[k * 400 + l * 20];
		   double *EVEV = &extEV[20 * l];
		   
		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d lv  = _mm_load_pd(&ll[j]);
		       __m128d rv  = _mm_load_pd(&rr[j]);
		       __m128d vll = _mm_load_pd(&vl[j]);
		       __m128d vrr = _mm_load_pd(&vr[j]);
		       
		       al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
		       ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
		     }  		 
		       
		   al = _mm_hadd_pd(al, al);
		   ar = _mm_hadd_pd(ar, ar);
		   
		   al = _mm_mul_pd(al, ar);

		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d vv  = _mm_load_pd(&v[j]);
		       __m128d EVV = _mm_load_pd(&EVEV[j]);

		       vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

		       _mm_store_pd(&v[j], vv);
		     }		  		   		  
		 }		 
#else
		 al = 0.0;
		 ar = 0.0;

		 for(j = 0; j < 20; j++)
		   {
		     al += vl[j] * left[k * 400 + l * 20 + j];
		     ar += vr[j] * right[k * 400 + l * 20 + j];
		   }

		 x1px2 = al * ar;

		 for(j = 0; j < 20; j++)
		   v[j] += x1px2 * extEV[20 * l + j];
#endif
	       }
	   }
	 

#ifdef __SIM_SSE3
	 { 
	   v = &(x3[80 * i]);
	   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	   
	   scale = 1;
	   for(l = 0; scale && (l < 80); l += 2)
	     {
	       __m128d vv = _mm_load_pd(&v[l]);
	       __m128d v1 = _mm_and_pd(vv, absMask.m);
	       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	       if(_mm_movemask_pd( v1 ) != 3)
		 scale = 0;
	     }	    	  
	 }
#else
	 v = &(x3[80 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 80); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
#endif

	 if (scale)
	   {
#ifdef __SIM_SSE3
	       __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	       
	       for(l = 0; l < 80; l+=2)
		 {
		   __m128d ex3v = _mm_load_pd(&v[l]);		  
		   _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		 }		   		  
#else	     
	     for(l = 0; l < 80; l++)
	       v[l] *= twotothe256;
#endif

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	  
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#ifdef __SIM_SSE3

static void newviewGTRGAMMAPROT_GAPPED(int tipCase,
				       double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, int gapLength, 
				       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
				       )
{
  double  *uX1, *uX2, *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0,   
    gapScaling = 0;
  double *vl, *vr, *x1v, *x2v;

  for(i = 0; i < gapLength; i++)
    x3_gap[i] = x1_gap[i] & x2_gap[i];

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double umpX1[1840], umpX2[1840];

	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
		double *ll =  &left[k * 20];
		double *rr =  &right[k * 20];
		
		__m128d umpX1v = _mm_setzero_pd();
		__m128d umpX2v = _mm_setzero_pd();

		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
		    umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
		umpX2v = _mm_hadd_pd(umpX2v, umpX2v);
		
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);
		_mm_storel_pd(&umpX2[80 * i + k], umpX2v);
	      }
	  }

	{
	  uX1 = &umpX1[1760];
	  uX2 = &umpX2[1760];

	  for(j = 0; j < 4; j++)
	    {
	      v = &x3_gapColumn[j * 20];

	      __m128d zero =  _mm_setzero_pd();
	      for(k = 0; k < 20; k+=2)		  		    
		_mm_store_pd(&v[k], zero);

	      for(k = 0; k < 20; k++)
		{ 
		  double *eev = &extEV[k * 20];
		  x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		  __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		  for(l = 0; l < 20; l+=2)
		    {
		      __m128d vv = _mm_load_pd(&v[l]);
		      __m128d ee = _mm_load_pd(&eev[l]);
		      
		      vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
		      
		      _mm_store_pd(&v[l], vv);
		    }
		}
	    }	   
	}	

	for(i = 0; i < n; i++)
	  {
	    if(!(x3_gap[i / 32] & mask32[i % 32]))
	      {

		uX1 = &umpX1[80 * tipX1[i]];
		uX2 = &umpX2[80 * tipX2[i]];
		
		for(j = 0; j < 4; j++)
		  {
		    v = &x3[i * 80 + j * 20];
		    
		    
		    __m128d zero =  _mm_setzero_pd();
		    for(k = 0; k < 20; k+=2)		  		    
		      _mm_store_pd(&v[k], zero);
		    
		    for(k = 0; k < 20; k++)
		      { 
			double *eev = &extEV[k * 20];
			x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
			__m128d x1px2v = _mm_set1_pd(x1px2);
			
			for(l = 0; l < 20; l+=2)
			  {
			    __m128d vv = _mm_load_pd(&v[l]);
			    __m128d ee = _mm_load_pd(&eev[l]);
			    
			    vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			    
			    _mm_store_pd(&v[l], vv);
			  }
		      }
		  }	   
	      }
	  }
      }
      break;
    case TIP_INNER:
      {
	double umpX1[1840], ump_x2[20];


	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {
		double *ll =  &left[k * 20];
				
		__m128d umpX1v = _mm_setzero_pd();
		
		for(l = 0; l < 20; l+=2)
		  {
		    __m128d vv = _mm_load_pd(&v[l]);
		    umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));		    					
		  }
		
		umpX1v = _mm_hadd_pd(umpX1v, umpX1v);				
		_mm_storel_pd(&umpX1[80 * i + k], umpX1v);		

	      }
	  }

	{
	  uX1 = &umpX1[1760];

	  for(k = 0; k < 4; k++)
	    {
	      v = &(x2_gapColumn[k * 20]);
	       
	      for(l = 0; l < 20; l++)
		{		   
		  double *r =  &right[k * 400 + l * 20];
		  __m128d ump_x2v = _mm_setzero_pd();	    
		  
		  for(j = 0; j < 20; j+= 2)
		    {
		      __m128d vv = _mm_load_pd(&v[j]);
		      __m128d rr = _mm_load_pd(&r[j]);
		      ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
		    }
		  
		  ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
		  
		  _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		}

	      v = &(x3_gapColumn[20 * k]);

	      __m128d zero =  _mm_setzero_pd();
	      for(l = 0; l < 20; l+=2)		  		    
		_mm_store_pd(&v[l], zero);
		  
	      for(l = 0; l < 20; l++)
		{
		  double *eev = &extEV[l * 20];
		  x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		  __m128d x1px2v = _mm_set1_pd(x1px2);
		  
		  for(j = 0; j < 20; j+=2)
		    {
		      __m128d vv = _mm_load_pd(&v[j]);
		      __m128d ee = _mm_load_pd(&eev[j]);
		      
		      vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
		      
		      _mm_store_pd(&v[j], vv);
		    }		     		    
		}			
	      
	    }
	  
	  { 
	    v = x3_gapColumn;
	    __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	    
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 2)
	      {
		__m128d vv = _mm_load_pd(&v[l]);
		__m128d v1 = _mm_and_pd(vv, absMask.m);
		v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		if(_mm_movemask_pd( v1 ) != 3)
		  scale = 0;
	      }	    	  
	  }


	  if (scale)
	    {
	      gapScaling = 1;
	      __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	      
	      for(l = 0; l < 80; l+=2)
		{
		  __m128d ex3v = _mm_load_pd(&v[l]);		  
		  _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		}		   		  	      	    	       
	    }
	}

	for (i = 0; i < n; i++)
	  {	    
	    if((x3_gap[i / 32] & mask32[i % 32]))
	       {	       
		 if(gapScaling)
		   {
		     if(useFastScaling)
		       addScale += wgt[i];
		     else
		       ex3[i]  += 1;
		   }
	       }
	     else
	       {
		 uX1 = &umpX1[80 * tipX1[i]];

		  if(x2_gap[i / 32] & mask32[i % 32])
		   x2v = x2_gapColumn;
		 else
		   x2v = &x2[i * 80];
		 
		 for(k = 0; k < 4; k++)
		   {
		     v = &(x2v[k * 20]);
		     
		     for(l = 0; l < 20; l++)
		       {		   
			 double *r =  &right[k * 400 + l * 20];
			 __m128d ump_x2v = _mm_setzero_pd();	    
			 
			 for(j = 0; j < 20; j+= 2)
			   {
			     __m128d vv = _mm_load_pd(&v[j]);
			     __m128d rr = _mm_load_pd(&r[j]);
			     ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
			   }
			 
			 ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);
			 
			 _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
		       }
		     
		     v = &(x3[80 * i + 20 * k]);
		     
		     __m128d zero =  _mm_setzero_pd();
		     for(l = 0; l < 20; l+=2)		  		    
		       _mm_store_pd(&v[l], zero);
		     
		     for(l = 0; l < 20; l++)
		       {
			 double *eev = &extEV[l * 20];
			 x1px2 = uX1[k * 20 + l]  * ump_x2[l];
			 __m128d x1px2v = _mm_set1_pd(x1px2);
			 
			 for(j = 0; j < 20; j+=2)
			   {
			     __m128d vv = _mm_load_pd(&v[j]);
			     __m128d ee = _mm_load_pd(&eev[j]);
			     
			     vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));
			     
			     _mm_store_pd(&v[j], vv);
			   }		     		    
		       }			
		     
		   }
		 
		 
		 { 
		   v = &(x3[80 * i]);
		   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
		   
		   scale = 1;
		   for(l = 0; scale && (l < 80); l += 2)
		     {
		       __m128d vv = _mm_load_pd(&v[l]);
		       __m128d v1 = _mm_and_pd(vv, absMask.m);
		       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		       if(_mm_movemask_pd( v1 ) != 3)
			 scale = 0;
		     }	    	  
		 }
		 
		 
		 if (scale)
		   {
		     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		     
		     for(l = 0; l < 80; l+=2)
		       {
			 __m128d ex3v = _mm_load_pd(&v[l]);		  
			 _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		       }		   		  
		     
		     if(useFastScaling)
		       addScale += wgt[i];
		     else
		       ex3[i]  += 1;	       
		   }
	       }
	  }
      }
      break;
    case INNER_INNER:
      {
	for(k = 0; k < 4; k++)
	   {
	     vl = &(x1_gapColumn[20 * k]);
	     vr = &(x2_gapColumn[20 * k]);
	     v =  &(x3_gapColumn[20 * k]);

	     __m128d zero =  _mm_setzero_pd();
	     for(l = 0; l < 20; l+=2)		  		    
	       _mm_store_pd(&v[l], zero);
	     
	     for(l = 0; l < 20; l++)
	       {		 
		 {
		   __m128d al = _mm_setzero_pd();
		   __m128d ar = _mm_setzero_pd();

		   double *ll   = &left[k * 400 + l * 20];
		   double *rr   = &right[k * 400 + l * 20];
		   double *EVEV = &extEV[20 * l];
		   
		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d lv  = _mm_load_pd(&ll[j]);
		       __m128d rv  = _mm_load_pd(&rr[j]);
		       __m128d vll = _mm_load_pd(&vl[j]);
		       __m128d vrr = _mm_load_pd(&vr[j]);
		       
		       al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
		       ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
		     }  		 
		       
		   al = _mm_hadd_pd(al, al);
		   ar = _mm_hadd_pd(ar, ar);
		   
		   al = _mm_mul_pd(al, ar);

		   for(j = 0; j < 20; j+=2)
		     {
		       __m128d vv  = _mm_load_pd(&v[j]);
		       __m128d EVV = _mm_load_pd(&EVEV[j]);

		       vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

		       _mm_store_pd(&v[j], vv);
		     }		  		   		  
		 }		 

	       }
	   }
	 

	{ 
	   v = x3_gapColumn;
	   __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	   
	   scale = 1;
	   for(l = 0; scale && (l < 80); l += 2)
	     {
	       __m128d vv = _mm_load_pd(&v[l]);
	       __m128d v1 = _mm_and_pd(vv, absMask.m);
	       v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
	       if(_mm_movemask_pd( v1 ) != 3)
		 scale = 0;
	     }	    	  
	 }

	 if (scale)
	   {
	     gapScaling = 1;
	     __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
	     
	     for(l = 0; l < 80; l+=2)
	       {
		 __m128d ex3v = _mm_load_pd(&v[l]);		  
		 _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
	       }		   		  
	     
	    	  
	   }
      }

      for (i = 0; i < n; i++)
       {
	  if(x3_gap[i / 32] & mask32[i % 32])
	   {	     
	     if(gapScaling)
	       {
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1; 	       
	       }
	   }
	 else
	   {
	     if(x1_gap[i / 32] & mask32[i % 32])
	       x1v = x1_gapColumn;
	     else
	       x1v = &x1[i * 80];

	     if(x2_gap[i / 32] & mask32[i % 32])
	       x2v = x2_gapColumn;
	     else
	       x2v = &x2[i * 80];

	     for(k = 0; k < 4; k++)
	       {
		 vl = &(x1v[20 * k]);
		 vr = &(x2v[20 * k]);
		 v =  &(x3[80 * i + 20 * k]);
		 		 
		 __m128d zero =  _mm_setzero_pd();
		 for(l = 0; l < 20; l+=2)		  		    
		   _mm_store_pd(&v[l], zero);
		 		 
		 for(l = 0; l < 20; l++)
		   {		 
		     {
		       __m128d al = _mm_setzero_pd();
		       __m128d ar = _mm_setzero_pd();
		       
		       double *ll   = &left[k * 400 + l * 20];
		       double *rr   = &right[k * 400 + l * 20];
		       double *EVEV = &extEV[20 * l];
		       
		       for(j = 0; j < 20; j+=2)
			 {
			   __m128d lv  = _mm_load_pd(&ll[j]);
			   __m128d rv  = _mm_load_pd(&rr[j]);
			   __m128d vll = _mm_load_pd(&vl[j]);
			   __m128d vrr = _mm_load_pd(&vr[j]);
			   
			   al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
			   ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
			 }  		 
		       
		       al = _mm_hadd_pd(al, al);
		       ar = _mm_hadd_pd(ar, ar);
		       
		       al = _mm_mul_pd(al, ar);
		       
		       for(j = 0; j < 20; j+=2)
			 {
			   __m128d vv  = _mm_load_pd(&v[j]);
			   __m128d EVV = _mm_load_pd(&EVEV[j]);
			   
			   vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));
			   
			   _mm_store_pd(&v[j], vv);
			 }		  		   		  
		     }		 
		     
		   }
	       }
	     

	     
	     { 
	       v = &(x3[80 * i]);
	       __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );
	       
	       scale = 1;
	       for(l = 0; scale && (l < 80); l += 2)
		 {
		   __m128d vv = _mm_load_pd(&v[l]);
		   __m128d v1 = _mm_and_pd(vv, absMask.m);
		   v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
		   if(_mm_movemask_pd( v1 ) != 3)
		     scale = 0;
		 }	    	  
	     }
	     
	     
	     if (scale)
	       {
		 __m128d twoto = _mm_set_pd(twotothe256, twotothe256);
		 
		 for(l = 0; l < 80; l+=2)
		   {
		     __m128d ex3v = _mm_load_pd(&v[l]);		  
		     _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
		   }		   		  
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;	  
	       }
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

static void newviewGTRGAMMAPROT_GAPPED_SAVE(int tipCase,
    double *x1, double *x2, double *x3, double *extEV, double *tipVector,
    int *ex3, unsigned char *tipX1, unsigned char *tipX2,
    int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
    unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap,  
    double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
    )
{
  double  *uX1, *uX2, *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0,   
       gapScaling = 0;
  double 
    *vl, *vr, *x1v, *x2v,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x3_ptr = x3;



  switch(tipCase)
  {
    case TIP_TIP:
      {
        double umpX1[1840], umpX2[1840];

        for(i = 0; i < 23; i++)
        {
          v = &(tipVector[20 * i]);

          for(k = 0; k < 80; k++)
          {
            double *ll =  &left[k * 20];
            double *rr =  &right[k * 20];

            __m128d umpX1v = _mm_setzero_pd();
            __m128d umpX2v = _mm_setzero_pd();

            for(l = 0; l < 20; l+=2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));
              umpX2v = _mm_add_pd(umpX2v, _mm_mul_pd(vv, _mm_load_pd(&rr[l])));					
            }

            umpX1v = _mm_hadd_pd(umpX1v, umpX1v);
            umpX2v = _mm_hadd_pd(umpX2v, umpX2v);

            _mm_storel_pd(&umpX1[80 * i + k], umpX1v);
            _mm_storel_pd(&umpX2[80 * i + k], umpX2v);
          }
        }

        {
          uX1 = &umpX1[1760];
          uX2 = &umpX2[1760];

          for(j = 0; j < 4; j++)
          {
            v = &x3_gapColumn[j * 20];

            __m128d zero =  _mm_setzero_pd();
            for(k = 0; k < 20; k+=2)		  		    
              _mm_store_pd(&v[k], zero);

            for(k = 0; k < 20; k++)
            { 
              double *eev = &extEV[k * 20];
              x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
              __m128d x1px2v = _mm_set1_pd(x1px2);

              for(l = 0; l < 20; l+=2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d ee = _mm_load_pd(&eev[l]);

                vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                _mm_store_pd(&v[l], vv);
              }
            }
          }	   
        }	

        for(i = 0; i < n; i++)
        {
          if(!(x3_gap[i / 32] & mask32[i % 32]))
          {
            uX1 = &umpX1[80 * tipX1[i]];
            uX2 = &umpX2[80 * tipX2[i]];

            for(j = 0; j < 4; j++)
            {
              v = &x3_ptr[j * 20];


              __m128d zero =  _mm_setzero_pd();
              for(k = 0; k < 20; k+=2)		  		    
                _mm_store_pd(&v[k], zero);

              for(k = 0; k < 20; k++)
              { 
                double *eev = &extEV[k * 20];
                x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
                __m128d x1px2v = _mm_set1_pd(x1px2);

                for(l = 0; l < 20; l+=2)
                {
                  __m128d vv = _mm_load_pd(&v[l]);
                  __m128d ee = _mm_load_pd(&eev[l]);

                  vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                  _mm_store_pd(&v[l], vv);
                }
              }
            }	   
            x3_ptr += 80;
          }
        }
      }
      break;
    case TIP_INNER:
      {
        double umpX1[1840], ump_x2[20];


        for(i = 0; i < 23; i++)
        {
          v = &(tipVector[20 * i]);

          for(k = 0; k < 80; k++)
          {
            double *ll =  &left[k * 20];

            __m128d umpX1v = _mm_setzero_pd();

            for(l = 0; l < 20; l+=2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              umpX1v = _mm_add_pd(umpX1v, _mm_mul_pd(vv, _mm_load_pd(&ll[l])));		    					
            }

            umpX1v = _mm_hadd_pd(umpX1v, umpX1v);				
            _mm_storel_pd(&umpX1[80 * i + k], umpX1v);		

          }
        }

        {
          uX1 = &umpX1[1760];

          for(k = 0; k < 4; k++)
          {
            v = &(x2_gapColumn[k * 20]);

            for(l = 0; l < 20; l++)
            {		   
              double *r =  &right[k * 400 + l * 20];
              __m128d ump_x2v = _mm_setzero_pd();	    

              for(j = 0; j < 20; j+= 2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                __m128d rr = _mm_load_pd(&r[j]);
                ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
              }

              ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

              _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
            }

            v = &(x3_gapColumn[20 * k]);

            __m128d zero =  _mm_setzero_pd();
            for(l = 0; l < 20; l+=2)		  		    
              _mm_store_pd(&v[l], zero);

            for(l = 0; l < 20; l++)
            {
              double *eev = &extEV[l * 20];
              x1px2 = uX1[k * 20 + l]  * ump_x2[l];
              __m128d x1px2v = _mm_set1_pd(x1px2);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv = _mm_load_pd(&v[j]);
                __m128d ee = _mm_load_pd(&eev[j]);

                vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                _mm_store_pd(&v[j], vv);
              }		     		    
            }			

          }

          { 
            v = x3_gapColumn;
            __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

            scale = 1;
            for(l = 0; scale && (l < 80); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }	    	  
          }


          if (scale)
          {
            gapScaling = 1;
            __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

            for(l = 0; l < 80; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);		  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
            }		   		  	      	    	       
          }
        }

        for (i = 0; i < n; i++)
        {	    
          if((x3_gap[i / 32] & mask32[i % 32]))
          {	       
            if(gapScaling)
            {
              if(useFastScaling)
                addScale += wgt[i];
              else
                ex3[i]  += 1;
            }
          }
          else
          {
            uX1 = &umpX1[80 * tipX1[i]];

            if(x2_gap[i / 32] & mask32[i % 32])
              x2v = x2_gapColumn;
            else
            {
              x2v = x2_ptr;
              x2_ptr += 80;
            }

            for(k = 0; k < 4; k++)
            {
              v = &(x2v[k * 20]);

              for(l = 0; l < 20; l++)
              {		   
                double *r =  &right[k * 400 + l * 20];
                __m128d ump_x2v = _mm_setzero_pd();	    

                for(j = 0; j < 20; j+= 2)
                {
                  __m128d vv = _mm_load_pd(&v[j]);
                  __m128d rr = _mm_load_pd(&r[j]);
                  ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(vv, rr));
                }

                ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

                _mm_storel_pd(&ump_x2[l], ump_x2v);		   		     
              }

              v = &x3_ptr[20 * k];

              __m128d zero =  _mm_setzero_pd();
              for(l = 0; l < 20; l+=2)		  		    
                _mm_store_pd(&v[l], zero);

              for(l = 0; l < 20; l++)
              {
                double *eev = &extEV[l * 20];
                x1px2 = uX1[k * 20 + l]  * ump_x2[l];
                __m128d x1px2v = _mm_set1_pd(x1px2);

                for(j = 0; j < 20; j+=2)
                {
                  __m128d vv = _mm_load_pd(&v[j]);
                  __m128d ee = _mm_load_pd(&eev[j]);

                  vv = _mm_add_pd(vv, _mm_mul_pd(x1px2v,ee));

                  _mm_store_pd(&v[j], vv);
                }		     		    
              }			

            }


            { 
              v = x3_ptr;
              __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

              scale = 1;
              for(l = 0; scale && (l < 80); l += 2)
              {
                __m128d vv = _mm_load_pd(&v[l]);
                __m128d v1 = _mm_and_pd(vv, absMask.m);
                v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
                if(_mm_movemask_pd( v1 ) != 3)
                  scale = 0;
              }	    	  
            }


            if (scale)
            {
              __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

              for(l = 0; l < 80; l+=2)
              {
                __m128d ex3v = _mm_load_pd(&v[l]);		  
                _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
              }		   		  

              if(useFastScaling)
                addScale += wgt[i];
              else
                ex3[i]  += 1;	       
            }

            x3_ptr += 80;
          }
        }
      }
      break;
    case INNER_INNER:
      {
        for(k = 0; k < 4; k++)
        {
          vl = &(x1_gapColumn[20 * k]);
          vr = &(x2_gapColumn[20 * k]);
          v =  &(x3_gapColumn[20 * k]);

          __m128d zero =  _mm_setzero_pd();
          for(l = 0; l < 20; l+=2)		  		    
            _mm_store_pd(&v[l], zero);

          for(l = 0; l < 20; l++)
          {		 
            {
              __m128d al = _mm_setzero_pd();
              __m128d ar = _mm_setzero_pd();

              double *ll   = &left[k * 400 + l * 20];
              double *rr   = &right[k * 400 + l * 20];
              double *EVEV = &extEV[20 * l];

              for(j = 0; j < 20; j+=2)
              {
                __m128d lv  = _mm_load_pd(&ll[j]);
                __m128d rv  = _mm_load_pd(&rr[j]);
                __m128d vll = _mm_load_pd(&vl[j]);
                __m128d vrr = _mm_load_pd(&vr[j]);

                al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
                ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
              }  		 

              al = _mm_hadd_pd(al, al);
              ar = _mm_hadd_pd(ar, ar);

              al = _mm_mul_pd(al, ar);

              for(j = 0; j < 20; j+=2)
              {
                __m128d vv  = _mm_load_pd(&v[j]);
                __m128d EVV = _mm_load_pd(&EVEV[j]);

                vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

                _mm_store_pd(&v[j], vv);
              }		  		   		  
            }		 

          }
        }


        { 
          v = x3_gapColumn;
          __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

          scale = 1;
          for(l = 0; scale && (l < 80); l += 2)
          {
            __m128d vv = _mm_load_pd(&v[l]);
            __m128d v1 = _mm_and_pd(vv, absMask.m);
            v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
            if(_mm_movemask_pd( v1 ) != 3)
              scale = 0;
          }	    	  
        }

        if (scale)
        {
          gapScaling = 1;
          __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

          for(l = 0; l < 80; l+=2)
          {
            __m128d ex3v = _mm_load_pd(&v[l]);		  
            _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
          }		   		  


        }
      }

      for (i = 0; i < n; i++)
      {
        if(x3_gap[i / 32] & mask32[i % 32])
        {	     
          if(gapScaling)
          {
            if(useFastScaling)
              addScale += wgt[i];
            else
              ex3[i]  += 1; 	       
          }
        }
        else
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

          for(k = 0; k < 4; k++)
          {
            vl = &(x1v[20 * k]);
            vr = &(x2v[20 * k]);
            v =  &x3_ptr[20 * k];

            __m128d zero =  _mm_setzero_pd();
            for(l = 0; l < 20; l+=2)		  		    
              _mm_store_pd(&v[l], zero);

            for(l = 0; l < 20; l++)
            {		 
              {
                __m128d al = _mm_setzero_pd();
                __m128d ar = _mm_setzero_pd();

                double *ll   = &left[k * 400 + l * 20];
                double *rr   = &right[k * 400 + l * 20];
                double *EVEV = &extEV[20 * l];

                for(j = 0; j < 20; j+=2)
                {
                  __m128d lv  = _mm_load_pd(&ll[j]);
                  __m128d rv  = _mm_load_pd(&rr[j]);
                  __m128d vll = _mm_load_pd(&vl[j]);
                  __m128d vrr = _mm_load_pd(&vr[j]);

                  al = _mm_add_pd(al, _mm_mul_pd(vll, lv));
                  ar = _mm_add_pd(ar, _mm_mul_pd(vrr, rv));
                }  		 

                al = _mm_hadd_pd(al, al);
                ar = _mm_hadd_pd(ar, ar);

                al = _mm_mul_pd(al, ar);

                for(j = 0; j < 20; j+=2)
                {
                  __m128d vv  = _mm_load_pd(&v[j]);
                  __m128d EVV = _mm_load_pd(&EVEV[j]);

                  vv = _mm_add_pd(vv, _mm_mul_pd(al, EVV));

                  _mm_store_pd(&v[j], vv);
                }		  		   		  
              }		 

            }
          }



          { 
            v = x3_ptr;
            __m128d minlikelihood_sse = _mm_set1_pd( minlikelihood );

            scale = 1;
            for(l = 0; scale && (l < 80); l += 2)
            {
              __m128d vv = _mm_load_pd(&v[l]);
              __m128d v1 = _mm_and_pd(vv, absMask.m);
              v1 = _mm_cmplt_pd(v1,  minlikelihood_sse);
              if(_mm_movemask_pd( v1 ) != 3)
                scale = 0;
            }	    	  
          }


          if (scale)
          {
            __m128d twoto = _mm_set_pd(twotothe256, twotothe256);

            for(l = 0; l < 80; l+=2)
            {
              __m128d ex3v = _mm_load_pd(&v[l]);		  
              _mm_store_pd(&v[l], _mm_mul_pd(ex3v,twoto));	
            }		   		  

            if(useFastScaling)
              addScale += wgt[i];
            else
              ex3[i]  += 1;	  
          }
          x3_ptr += 80;
        }
      }
      break;
    default:
      assert(0);
  }

  if(useFastScaling)
    *scalerIncrement = addScale;

}

#else

static void newviewGTRGAMMAPROT_GAPPED(int tipCase,
				       double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling,
				       unsigned int *x1_gap, unsigned int *x2_gap, unsigned int *x3_gap, int gapLength, 
				       double *x1_gapColumn, double *x2_gapColumn, double *x3_gapColumn
				       )
{
  double  
    *uX1, 
    *uX2, 
    *v,
    *x1v,
    *x2v,
    *x3v,
    x1px2,
    *vl, 
    *vr, 
    al, 
    ar;
  
  int  
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0,
    scaleGap = 0;

  for(i = 0; i < gapLength; i++)
    x3_gap[i] = x1_gap[i] & x2_gap[i];

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double umpX1[1840], umpX2[1840];

	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {

		umpX1[80 * i + k] = 0.0;
		umpX2[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  {
		    umpX1[80 * i + k] +=  v[l] *  left[k * 20 + l];
		    umpX2[80 * i + k] +=  v[l] * right[k * 20 + l];
		  }
	      }
	  }

	uX1 = &umpX1[1760];
	uX2 = &umpX2[1760];
	    
	for(j = 0; j < 4; j++)
	  {
	    v = &x3_gapColumn[j * 20];

	    for(k = 0; k < 20; k++)
	      v[k] = 0.0;

	    for(k = 0; k < 20; k++)
	      {		   
		x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		   
		for(l = 0; l < 20; l++)		      					
		  v[l] += x1px2 * extEV[20 * k + l];		     
	      }
	  }	   

	for(i = 0; i < n; i++)
	  {
	    if(!(x3_gap[i / 32] & mask32[i % 32]))
	      {
		uX1 = &umpX1[80 * tipX1[i]];
		uX2 = &umpX2[80 * tipX2[i]];
		
		for(j = 0; j < 4; j++)
		  {
		    v = &x3[i * 80 + j * 20];
		    		    
		    for(k = 0; k < 20; k++)
		      v[k] = 0.0;
		    
		    for(k = 0; k < 20; k++)
		      {		   
			x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
			
			for(l = 0; l < 20; l++)		      					
			  v[l] += x1px2 * extEV[20 * k + l];		     
		      }
		  }	   
	      }
	  }
      }
      break;
    case TIP_INNER:
      {
	double umpX1[1840], ump_x2[20];

	for(i = 0; i < 23; i++)
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++)
	      {

		umpX1[80 * i + k] = 0.0;

		for(l = 0; l < 20; l++)
		  umpX1[80 * i + k] +=  v[l] * left[k * 20 + l];


	      }
	  }

	
	uX1 = &umpX1[1760];

	for(k = 0; k < 4; k++)
	  {
	    v = &(x2_gapColumn[k * 20]);

	    for(l = 0; l < 20; l++)
	      {
		ump_x2[l] = 0.0;
		
		for(j = 0; j < 20; j++)
		  ump_x2[l] += v[j] * right[k * 400 + l * 20 + j];
	      }

	    v = &(x3_gapColumn[20 * k]);

	    for(l = 0; l < 20; l++)
	      v[l] = 0;
	    
	    for(l = 0; l < 20; l++)
	      {
		x1px2 = uX1[k * 20 + l]  * ump_x2[l];
		for(j = 0; j < 20; j++)
		  v[j] += x1px2 * extEV[l * 20  + j];
	      }
	  }
	
	v = x3_gapColumn;
	scale = 1;
	for(l = 0; scale && (l < 80); l++)
	  scale = (ABS(v[l]) <  minlikelihood);

	if(scale)
	  {
	    scaleGap = 1;

	    for(l = 0; l < 80; l++)
	      v[l] *= twotothe256;	    	    	       
	  }

	for (i = 0; i < n; i++)
	  {
	    if(x3_gap[i / 32] & mask32[i % 32])
	      {	       
		if(scaleGap)
		  {
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;
		  }
	      }
	    else
	      {
		uX1 = &umpX1[80 * tipX1[i]];

		if(x2_gap[i / 32] & mask32[i % 32])
		  x2v = x2_gapColumn;
		else
		  x2v = &x2[i * 80];

		for(k = 0; k < 4; k++)
		  {
		    v = &(x2v[k * 20]);

		    for(l = 0; l < 20; l++)
		      {
			ump_x2[l] = 0.0;
			
			for(j = 0; j < 20; j++)
			  ump_x2[l] += v[j] * right[k * 400 + l * 20 + j];
		      }
		    
		    v = &(x3[80 * i + 20 * k]);
		    
		    for(l = 0; l < 20; l++)
		      v[l] = 0;
		    
		    for(l = 0; l < 20; l++)
		      {
			x1px2 = uX1[k * 20 + l]  * ump_x2[l];
			for(j = 0; j < 20; j++)
			  v[j] += x1px2 * extEV[l * 20  + j];
		      }
		  }
		

		v = &x3[80 * i];
		scale = 1;
		for(l = 0; scale && (l < 80); l++)
		  scale = (ABS(v[l]) <  minlikelihood);
		
		
		if (scale)
		  {
		    for(l = 0; l < 80; l++)
		      v[l] *= twotothe256;
		    
		    if(useFastScaling)
		      addScale += wgt[i];
		    else
		      ex3[i]  += 1;	       
		  }
	      }
	  }
      }
      break;
    case INNER_INNER:

      for(k = 0; k < 4; k++)
	{
	  vl = &(x1_gapColumn[20 * k]);
	  vr = &(x2_gapColumn[20 * k]);
	  v =  &(x3_gapColumn[20 * k]);

	  for(l = 0; l < 20; l++)
	    v[l] = 0;

	  for(l = 0; l < 20; l++)
	    {		 
	      al = 0.0;
	      ar = 0.0;
	      
	      for(j = 0; j < 20; j++)
		{
		  al += vl[j] * left[k * 400 + l * 20 + j];
		  ar += vr[j] * right[k * 400 + l * 20 + j];
		}
	      
	      x1px2 = al * ar;
	      
	      for(j = 0; j < 20; j++)
		v[j] += x1px2 * extEV[20 * l + j];
	    }
	}
	 
      v = x3_gapColumn;
      scale = 1;
      for(l = 0; scale && (l < 80); l++)
	scale = ((ABS(v[l]) <  minlikelihood));

      if(scale)
	{
	  scaleGap = 1;
     
	  for(l = 0; l < 80; l++)
	    v[l] *= twotothe256;	  	 
	}

      for (i = 0; i < n; i++)
       {
	 if(x3_gap[i / 32] & mask32[i % 32])
	   {	     
	     if(scaleGap)
	       {
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1; 	       
	       }
	   }
	 else
	   {	    
	     if(x1_gap[i / 32] & mask32[i % 32])
	       x1v = x1_gapColumn;
	     else
	       x1v = &x1[i * 80];

	     if(x2_gap[i / 32] & mask32[i % 32])
	       x2v = x2_gapColumn;
	     else
	       x2v = &x2[i * 80];

	     x3v = &x3[i * 80];

	     for(k = 0; k < 4; k++)
	       {
		 vl = &(x1v[20 * k]);
		 vr = &(x2v[20 * k]);
		 v =  &(x3v[20 * k]);

		 for(l = 0; l < 20; l++)
		   v[l] = 0;

		 for(l = 0; l < 20; l++)
		   {		 
		     al = 0.0;
		     ar = 0.0;
		     
		     for(j = 0; j < 20; j++)
		       {
			 al += vl[j] * left[k * 400 + l * 20 + j];
			 ar += vr[j] * right[k * 400 + l * 20 + j];
		       }
		     
		     x1px2 = al * ar;
		     
		     for(j = 0; j < 20; j++)
		       v[j] += x1px2 * extEV[20 * l + j];
		   }
	       }
	 
	     v = &(x3[80 * i]);
	     scale = 1;
	     for(l = 0; scale && (l < 80); l++)
	       scale = ((ABS(v[l]) <  minlikelihood));

	     if(scale)
	       {		 
		 for(l = 0; l < 80; l++)
		   v[l] *= twotothe256;
		 
		 if(useFastScaling)
		   addScale += wgt[i];
		 else
		   ex3[i]  += 1;	  
	       }
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


#endif







static void newviewGTRGAMMASECONDARY(int tipCase,
				     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				     int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				     int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double  *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr, al, ar;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[16 * tipX1[i]]);
		vr = &(tipVector[16 * tipX2[i]]);
		v =  &(x3[64 * i + 16 * k]);

		for(l = 0; l < 16; l++)
		  v[l] = 0;

		for(l = 0; l < 16; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 16; j++)
		      {
			al += vl[j] * left[k * 256 + l * 16 + j];
			ar += vr[j] * right[k * 256 + l * 16 + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < 16; j++)
		      v[j] += x1px2 * extEV[16 * l + j];
		  }
	      }	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[16 * tipX1[i]]);
		vr = &(x2[64 * i + 16 * k]);
		v =  &(x3[64 * i + 16 * k]);

		for(l = 0; l < 16; l++)
		  v[l] = 0;

		for(l = 0; l < 16; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 16; j++)
		      {
			al += vl[j] * left[k * 256 + l * 16 + j];
			ar += vr[j] * right[k * 256 + l * 16 + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < 16; j++)
		      v[j] += x1px2 * extEV[16 * l + j];
		  }
	      }
	   
	    v = &x3[64 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 64); l++)
	      scale = (ABS(v[l]) <  minlikelihood);

	    if (scale)
	      {
		for(l = 0; l < 64; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[64 * i + 16 * k]);
	     vr = &(x2[64 * i + 16 * k]);
	     v =  &(x3[64 * i + 16 * k]);

	     for(l = 0; l < 16; l++)
	       v[l] = 0;

	     for(l = 0; l < 16; l++)
	       {
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 16; j++)
		   {
		     al += vl[j] * left[k * 256 + l * 16 + j];
		     ar += vr[j] * right[k * 256 + l * 16 + j];
		   }

		 x1px2 = al * ar;
		 for(j = 0; j < 16; j++)
		   v[j] += x1px2 * extEV[16 * l + j];
	       }
	   }
	 
	 v = &(x3[64 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 64); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if (scale)
	   {
	     for(l = 0; l < 64; l++)
	       v[l] *= twotothe256;

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	    
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}



static void newviewGTRGAMMASECONDARY_6(int tipCase,
				       double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double  *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr, al, ar;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[6 * tipX1[i]]);
		vr = &(tipVector[6 * tipX2[i]]);
		v =  &(x3[24 * i + 6 * k]);

		for(l = 0; l < 6; l++)
		  v[l] = 0;

		for(l = 0; l < 6; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 6; j++)
		      {
			al += vl[j] * left[k * 36 + l * 6 + j];
			ar += vr[j] * right[k * 36 + l * 6 + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < 6; j++)
		      v[j] += x1px2 * extEV[6 * l + j];
		  }
	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[6 * tipX1[i]]);
		vr = &(x2[24 * i + 6 * k]);
		v =  &(x3[24 * i + 6 * k]);

		for(l = 0; l < 6; l++)
		  v[l] = 0;

		for(l = 0; l < 6; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 6; j++)
		      {
			al += vl[j] * left[k * 36 + l * 6 + j];
			ar += vr[j] * right[k * 36 + l * 6 + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < 6; j++)
		      v[j] += x1px2 * extEV[6 * l + j];
		  }
	      }
	   
	    v = &x3[24 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 24); l++)
	      scale = (ABS(v[l]) <  minlikelihood);

	    if(scale)
	      {
		for(l = 0; l < 24; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;		
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[24 * i + 6 * k]);
	     vr = &(x2[24 * i + 6 * k]);
	     v =  &(x3[24 * i + 6 * k]);

	     for(l = 0; l < 6; l++)
	       v[l] = 0;

	     for(l = 0; l < 6; l++)
	       {
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 6; j++)
		   {
		     al += vl[j] * left[k * 36 + l * 6 + j];
		     ar += vr[j] * right[k * 36 + l * 6 + j];
		   }

		 x1px2 = al * ar;
		 for(j = 0; j < 6; j++)
		   v[j] += x1px2 * extEV[6 * l + j];
	       }
	   }
	 
	 v = &(x3[24 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 24); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if (scale)
	   {
	     for(l = 0; l < 24; l++)
	       v[l] *= twotothe256;

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	   
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;

}


static void newviewGTRGAMMASECONDARY_7(int tipCase,
				       double *x1, double *x2, double *x3, double *extEV, double *tipVector,
				       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double  *v;
  double x1px2;
  int  i, j, l, k, scale, addScale = 0;
  double *vl, *vr, al, ar;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[7 * tipX1[i]]);
		vr = &(tipVector[7 * tipX2[i]]);
		v =  &(x3[28 * i + 7 * k]);

		for(l = 0; l < 7; l++)
		  v[l] = 0;

		for(l = 0; l < 7; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 7; j++)
		      {
			al += vl[j] * left[k * 49 + l * 7 + j];
			ar += vr[j] * right[k * 49 + l * 7 + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < 7; j++)
		      v[j] += x1px2 * extEV[7 * l + j];
		  }
	      }	   
	  }
      }
      break;
    case TIP_INNER:
      {
	for (i = 0; i < n; i++)
	  {
	    for(k = 0; k < 4; k++)
	      {
		vl = &(tipVector[7 * tipX1[i]]);
		vr = &(x2[28 * i + 7 * k]);
		v =  &(x3[28 * i + 7 * k]);

		for(l = 0; l < 7; l++)
		  v[l] = 0;

		for(l = 0; l < 7; l++)
		  {
		    al = 0.0;
		    ar = 0.0;
		    for(j = 0; j < 7; j++)
		      {
			al += vl[j] * left[k * 49 + l * 7 + j];
			ar += vr[j] * right[k * 49 + l * 7 + j];
		      }

		    x1px2 = al * ar;
		    for(j = 0; j < 7; j++)
		      v[j] += x1px2 * extEV[7 * l + j];
		  }
	      }
	   
	    v = &x3[28 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 28); l++)
	      scale = (ABS(v[l]) <  minlikelihood);

	    if (scale)
	      {
		for(l = 0; l < 28; l++)
		  v[l] *= twotothe256;

		if(useFastScaling)
		  addScale += wgt[i];
		else
		  ex3[i]  += 1;	      
	      }
	  }
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
       {
	 for(k = 0; k < 4; k++)
	   {
	     vl = &(x1[28 * i + 7 * k]);
	     vr = &(x2[28 * i + 7 * k]);
	     v =  &(x3[28 * i + 7 * k]);

	     for(l = 0; l < 7; l++)
	       v[l] = 0;

	     for(l = 0; l < 7; l++)
	       {
		 al = 0.0;
		 ar = 0.0;
		 for(j = 0; j < 7; j++)
		   {
		     al += vl[j] * left[k * 49 + l * 7 + j];
		     ar += vr[j] * right[k * 49 + l * 7 + j];
		   }

		 x1px2 = al * ar;
		 for(j = 0; j < 7; j++)
		   v[j] += x1px2 * extEV[7 * l + j];
	       }
	   }
	 
	 v = &(x3[28 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 28); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if (scale)
	   {
	     for(l = 0; l < 28; l++)
	       v[l] *= twotothe256;

	     if(useFastScaling)
	       addScale += wgt[i];
	     else
	       ex3[i]  += 1;	   
	   }
       }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
}







void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches)
{
  if(isTip(p->number, maxTips))
    return;

  {
    int i;
    nodeptr q = p->next->back;
    nodeptr r = p->next->next->back;

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {
	while (! p->x)
	 {
	   if (! p->x)
	     getxnode(p);
	 }

	ti[*counter].tipCase = TIP_TIP;
	ti[*counter].pNumber = p->number;
	ti[*counter].qNumber = q->number;
	ti[*counter].rNumber = r->number;

	for(i = 0; i < numBranches; i++)
	  {
	    double z;
	    z = q->z[i];
	  
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].qz[i] = z;

	    z = r->z[i];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].rz[i] = z;
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
		tmp = r;
		r = q;
		q = tmp;
	      }

	    while ((! p->x) || (! r->x))
	      {
		if (! r->x)
		  computeTraversalInfo(r, ti, counter, maxTips, numBranches);
		if (! p->x)
		  getxnode(p);
	      }

	    ti[*counter].tipCase = TIP_INNER;
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;

	    for(i = 0; i < numBranches; i++)
	      {
		double z;
		z = q->z[i];
			
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[i] = z;

		z = r->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[i] = z;
	      }

	    *counter = *counter + 1;
	  }
	else
	  {

	    while ((! p->x) || (! q->x) || (! r->x))
	      {
		if (! q->x)
		  computeTraversalInfo(q, ti, counter, maxTips, numBranches);
		if (! r->x)
		  computeTraversalInfo(r, ti, counter, maxTips, numBranches);
		if (! p->x)
		  getxnode(p);
	      }

	    ti[*counter].tipCase = INNER_INNER;
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;

	    for(i = 0; i < numBranches; i++)
	      {
		double z;
		z = q->z[i];	

		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[i] = z;

		z = r->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[i] = z;
	      }

	    *counter = *counter + 1;
	  }
      }
  }

  

}



void computeTraversalInfoMulti(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model)
{
  if(isTip(p->number, maxTips))
    {
      assert(p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);
      assert(p->backs[model]);
      return;
    }

 
  assert(p->backs[model]);

  {     
    nodeptr q = p->next->backs[model];
    nodeptr r = p->next->next->backs[model];
    
    assert(p == p->next->next->next);
      
    assert(q && r);

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {	  
	while (! p->xs[model])
	 {	  
	   if (! p->xs[model])
	     getxsnode(p, model); 	   
	 }

	assert(p->xs[model]);

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
		tmp = r;
		r = q;
		q = tmp;
	      }

	    while ((! p->xs[model]) || (! r->xs[model])) 
	      {	 			
		if (! r->xs[model]) 
		  computeTraversalInfoMulti(r, ti, counter, maxTips, model);
		if (! p->xs[model]) 
		  getxsnode(p, model);	
	      }
	    	   
	    assert(p->xs[model] && r->xs[model]);

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

	    while ((! p->xs[model]) || (! q->xs[model]) || (! r->xs[model])) 
	      {		
		if (! q->xs[model]) 
		  computeTraversalInfoMulti(q, ti, counter, maxTips, model);
		if (! r->xs[model]) 
		  computeTraversalInfoMulti(r, ti, counter, maxTips, model);
		if (! p->xs[model]) 
		  getxsnode(p, model);	
	      }

	    assert(p->xs[model] && r->xs[model] && q->xs[model]);

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

}



void newviewIterative (tree *tr)
{
  traversalInfo 
    *ti   = tr->td[0].ti;
  
  int 
    i, 
    model;

  for(i = 1; i < tr->td[0].count; i++)
    {
      traversalInfo 
	*tInfo = &ti[i];

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->executeModel[model])
	    {	      
	      double
		*x1_start = (double*)NULL,
		*x2_start = (double*)NULL,
		*x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1],
		*left     = (double*)NULL,
		*right    = (double*)NULL,
		*x1_gapColumn = (double*)NULL,
		*x2_gapColumn = (double*)NULL,
		*x3_gapColumn = (double*)NULL,	       
		qz, 
		rz;
	      
	      unsigned int
		*x1_gap = (unsigned int*)NULL,
		*x2_gap = (unsigned int*)NULL,
		*x3_gap = (unsigned int*)NULL;
	      	     
	      int		
		scalerIncrement = 0,
		*wgt = (int*)NULL,	       
		*ex3 = (int*)NULL;
	      
	      unsigned char
		*tipX1 = (unsigned char *)NULL,
		*tipX2 = (unsigned char *)NULL;
	    	      
	      size_t
		gapOffset,
		rateHet,
		states = (size_t)tr->partitionData[model].states,
		width = tr->partitionData[model].width,
		availableLength =  tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)],			       
		requiredLength;
	       
	      if(tr->rateHetModel == CAT)
		rateHet = 1;
	      else
		rateHet = 4;

	      if(tr->useGappedImplementation || tr->saveMemory)
		{
		  x1_gap = &(tr->partitionData[model].gapVector[tInfo->qNumber * tr->partitionData[model].gapVectorLength]);
		  x2_gap = &(tr->partitionData[model].gapVector[tInfo->rNumber * tr->partitionData[model].gapVectorLength]);
		  x3_gap = &(tr->partitionData[model].gapVector[tInfo->pNumber * tr->partitionData[model].gapVectorLength]);
		}
	      
	      if(tr->saveMemory)
		{
		   size_t
		     j,
		     setBits = 0;
		   
		   gapOffset = gapOffset = states * (size_t)getUndetermined(tr->partitionData[model].dataType);

		   availableLength = tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)];
		   
		   for(j = 0; j < (size_t)tr->partitionData[model].gapVectorLength; j++)
		     {
		       x3_gap[j] = x1_gap[j] & x2_gap[j];
		       setBits += (size_t)(precomputed16_bitcount(x3_gap[j]));
		     }
		   
		   requiredLength = (width - setBits)  * rateHet * states * sizeof(double);	
		}
	      else
		requiredLength  =  width * rateHet * states * sizeof(double);

	      if(requiredLength != availableLength)
		{		  
		  if(x3_start)
		    free(x3_start);
		 
		  x3_start = (double*)malloc_aligned(requiredLength);		 
		  
		  tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1] = x3_start;
		  
		  tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)] = requiredLength;		 
		}
	      

	      if(tr->useFastScaling)		
		wgt   =  tr->partitionData[model].wgt;	
	      else
		{
		  size_t
		    availableExpLength = tr->partitionData[model].expSpaceVector[(tInfo->pNumber - tr->mxtips - 1)],
		    requiredExpLength  = width * sizeof(int);
		  
		  ex3 = tr->partitionData[model].expVector[tInfo->pNumber - tr->mxtips - 1];
		  
		  if(requiredExpLength != availableExpLength)
		    {
		      if(ex3)
			free(ex3);
		 
		      ex3 = (int*)malloc_aligned(requiredExpLength);		 
		  
		      tr->partitionData[model].expVector[tInfo->pNumber - tr->mxtips - 1] = ex3;
		  
		      tr->partitionData[model].expSpaceVector[(tInfo->pNumber - tr->mxtips - 1)] = requiredExpLength;
		    }
		}

	      switch(tInfo->tipCase)
		{
		case TIP_TIP:		  
		  tipX1    = tr->partitionData[model].yVector[tInfo->qNumber];
		  tipX2    = tr->partitionData[model].yVector[tInfo->rNumber];		 		 		 
		  
		  if(tr->useGappedImplementation || tr->saveMemory)
		    {
		      x1_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);
		      x2_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);
		      x3_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
		    }
		  
		  if(!tr->useFastScaling)
		    {
		      size_t
			k;		      		     

		      for(k = 0; k < width; k++)
			ex3[k] = 0;
		    }
		  break;
		case TIP_INNER:		 
		  tipX1    =  tr->partitionData[model].yVector[tInfo->qNumber];		 
		  x2_start       = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];
		 		  
		  if(tr->useGappedImplementation || tr->saveMemory)
		    {
		      x1_gapColumn   = &(tr->partitionData[model].tipVector[gapOffset]);
		      x2_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->rNumber - tr->mxtips - 1) * states * rateHet];
		      x3_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
		    }		
		  
		  if(!tr->useFastScaling)
		    {
		      size_t 
			k;
		      int
			*ex2;
		      
		      ex2      = tr->partitionData[model].expVector[tInfo->rNumber - tr->mxtips - 1];		     
		      
		      for(k = 0; k < width; k++)
			ex3[k] = ex2[k];
		    }
		  break;
		case INNER_INNER:		 		 		    
		  x1_start       = tr->partitionData[model].xVector[tInfo->qNumber - tr->mxtips - 1];
		  x2_start       = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];
		      
		  if(tr->useGappedImplementation || tr->saveMemory)
		    {		      
		      x1_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->qNumber - tr->mxtips - 1) * states * rateHet];
		      x2_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->rNumber - tr->mxtips - 1) * states * rateHet];
		      x3_gapColumn   = &tr->partitionData[model].gapColumn[(tInfo->pNumber - tr->mxtips - 1) * states * rateHet];
		    }
		  
		  if(!tr->useFastScaling)
		    {
		      size_t
			k;
		      
		      int
			*ex1,
			*ex2;

		      ex1      = tr->partitionData[model].expVector[tInfo->qNumber - tr->mxtips - 1];
		      ex2      = tr->partitionData[model].expVector[tInfo->rNumber - tr->mxtips - 1];		      
		      
		      for(k = 0; k < width; k++)
			ex3[k] = ex1[k] + ex2[k];
		    }		  
		  break;
		default:
		  assert(0);
		}
	     		
	      left  = tr->partitionData[model].left;
	      right = tr->partitionData[model].right;
	      
	      if(tr->multiBranch)
		{
		  qz = tInfo->qz[model];
		  rz = tInfo->rz[model];
		}
	      else
		{
		  qz = tInfo->qz[0];
		  rz = tInfo->rz[0];
		}
	      	     	      
	      switch(tr->partitionData[model].dataType)
		{
		case BINARY_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:
		      {			
			makeP(qz, rz, tr->partitionData[model].perSiteRates,   tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN, tr->partitionData[model].numberOfCategories,
			      left, right, BINARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRCAT_BINARY(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					     x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					     ex3, tipX1, tipX2,
					     width, left, right, wgt, &scalerIncrement, tr->useFastScaling
					     );		
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI, tr->partitionData[model].EIGN,
			      4, left, right, BINARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMA_BINARY(tInfo->tipCase,
					       x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					       ex3, tipX1, tipX2,
					       width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			  makeP(qz, rz, tr->partitionData[model].perSiteRates,   tr->partitionData[model].EI,
				tr->partitionData[model].EIGN, tr->partitionData[model].numberOfCategories,
				left, right, DNA_DATA, tr->saveMemory, tr->maxCategories);

#ifdef __SIM_SSE3
			  if(tr->saveMemory)
			    {
			      assert(tr->useFastScaling);			      
			      newviewGTRCAT_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2,
						 width, left, right, wgt, &scalerIncrement, TRUE, x1_gap, x2_gap, x3_gap,
						 x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
			    }
			  else
#endif
			    {
#ifdef __AVX
			      newviewGTRCAT_AVX(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						ex3, tipX1, tipX2,
						width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else
			      newviewGTRCAT(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					    ex3, tipX1, tipX2,
					    width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#endif
			    }
			}
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      
			{			 
			  makeP(qz, rz, tr->partitionData[model].gammaRates,
				tr->partitionData[model].EI, tr->partitionData[model].EIGN,
				4, left, right, DNA_DATA, tr->saveMemory, tr->maxCategories);
		  

#ifdef __SIM_SSE3
			  if(tr->saveMemory)
			    {
			     
			    newviewGTRGAMMA_GAPPED_SAVE(tInfo->tipCase,
							x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
							ex3, tipX1, tipX2,
							width, left, right, wgt, &scalerIncrement, tr->useFastScaling,
							x1_gap, x2_gap, x3_gap, 
							x1_gapColumn, x2_gapColumn, x3_gapColumn);
			   
			    }
			  else
#endif
			    {
			      if(tr->useGappedImplementation)
				newviewGTRGAMMA_GAPPED(tInfo->tipCase,
						       x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						       ex3, tipX1, tipX2,
						       width, left, right, wgt, &scalerIncrement, tr->useFastScaling,
						       x1_gap, x2_gap, x3_gap, tr->partitionData[model].gapVectorLength, 
						       x1_gapColumn, x2_gapColumn, x3_gapColumn);
			      else
#ifdef __AVX
				newviewGTRGAMMA_AVX(tInfo->tipCase,
						    x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						    ex3, tipX1, tipX2,
						    width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else
				newviewGTRGAMMA(tInfo->tipCase,
						x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
						ex3, tipX1, tipX2,
						width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#endif
			    }
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
			  makeP(qz, rz, tr->partitionData[model].perSiteRates,
				tr->partitionData[model].EI,
				tr->partitionData[model].EIGN,
				tr->partitionData[model].numberOfCategories, left, right, AA_DATA, tr->saveMemory, tr->maxCategories);

#ifdef __SIM_SSE3
			  if(tr->saveMemory)
			    newviewGTRCATPROT_SAVE(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						   x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						   ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, TRUE, x1_gap, x2_gap, x3_gap,
						   x1_gapColumn, x2_gapColumn, x3_gapColumn, tr->maxCategories);
			  else
#endif

#ifdef __AVX
			  newviewGTRCATPROT_AVX(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else			  
			  newviewGTRCATPROT(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					    ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);			
#endif
			}
		      break;
		    case GAMMA:
		    case GAMMA_I:		     		      
			{
			  if(tr->estimatePerSiteAA)
			    {
			      int 
				p;
			   
			      for(p = 0; p < (NUM_PROT_MODELS - 3); p++)				
				  makeP(qz, rz, tr->partitionData[model].gammaRates,
					tr->siteProtModel[p].EI,
					tr->siteProtModel[p].EIGN,
					4, 
					tr->siteProtModel[p].left, 
					tr->siteProtModel[p].right, 
					AA_DATA, tr->saveMemory, tr->maxCategories);				
			      
			      newviewFlexGamma_perSite(tInfo->tipCase,
						       x1_start, x2_start, x3_start,
						       tr->partitionData[model].perSiteAAModel,
						       tr->siteProtModel,						      
						       ex3, tipX1, tipX2,
						       width, wgt, &scalerIncrement, tr->useFastScaling, 20);			      			      
			    }
			  else
			    {
			      makeP(qz, rz, tr->partitionData[model].gammaRates,
				    tr->partitionData[model].EI,
				    tr->partitionData[model].EIGN,
				    4, left, right, AA_DATA, tr->saveMemory, tr->maxCategories);
#ifdef __SIM_SSE3
			      if(tr->saveMemory)
				newviewGTRGAMMAPROT_GAPPED_SAVE(tInfo->tipCase,
								x1_start, x2_start, x3_start,
								tr->partitionData[model].EV,
								tr->partitionData[model].tipVector,
								ex3, tipX1, tipX2,
								width, left, right, wgt, &scalerIncrement, tr->useFastScaling,
								x1_gap, x2_gap, x3_gap,
								x1_gapColumn, x2_gapColumn, x3_gapColumn);
			      else
#endif
				{
				  if(tr->useGappedImplementation)
				    newviewGTRGAMMAPROT_GAPPED(tInfo->tipCase,
							       x1_start, x2_start, x3_start,
							       tr->partitionData[model].EV,
							       tr->partitionData[model].tipVector,
							       ex3, tipX1, tipX2,
							       width, left, right, wgt, &scalerIncrement, tr->useFastScaling,
							       x1_gap, x2_gap, x3_gap, tr->partitionData[model].gapVectorLength, 
							       x1_gapColumn, x2_gapColumn, x3_gapColumn);
				  else
#ifdef __AVX
				    newviewGTRGAMMAPROT_AVX(tInfo->tipCase,
							    x1_start, x2_start, x3_start,
							    tr->partitionData[model].EV,
							    tr->partitionData[model].tipVector,
							    ex3, tipX1, tipX2,
							    width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else
				    newviewGTRGAMMAPROT(tInfo->tipCase,
							x1_start, x2_start, x3_start,
							tr->partitionData[model].EV,
							tr->partitionData[model].tipVector,
							ex3, tipX1, tipX2,
							width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#endif
				}
			    }
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
			makeP(qz, rz, tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA_6, tr->saveMemory, tr->maxCategories);

			newviewGTRCATSECONDARY_6(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      4, left, right, SECONDARY_DATA_6, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMASECONDARY_6(tInfo->tipCase,
						   x1_start, x2_start, x3_start,
						   tr->partitionData[model].EV,
						   tr->partitionData[model].tipVector,
						   ex3, tipX1, tipX2,
						   width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			makeP(qz, rz, tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA_7, tr->saveMemory, tr->maxCategories);

			newviewGTRCATSECONDARY_7(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      4, left, right, SECONDARY_DATA_7, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMASECONDARY_7(tInfo->tipCase,
						   x1_start, x2_start, x3_start,
						   tr->partitionData[model].EV,
						   tr->partitionData[model].tipVector,
						   ex3, tipX1, tipX2,
						   width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			makeP(qz, rz, tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRCATSECONDARY(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					       ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      4, left, right, SECONDARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMASECONDARY(tInfo->tipCase,
						 x1_start, x2_start, x3_start,
						 tr->partitionData[model].EV,
						 tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2,
						 width, left, right, wgt, &scalerIncrement, tr->useFastScaling);		
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
			makeP_Flex(qz, rz, tr->partitionData[model].perSiteRates,
				   tr->partitionData[model].EI,
				   tr->partitionData[model].EIGN,
				   tr->partitionData[model].numberOfCategories, left, right, states);
			
			newviewFlexCat(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
				       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				       ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling, states);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP_Flex(qz, rz, tr->partitionData[model].gammaRates,
				   tr->partitionData[model].EI,
				   tr->partitionData[model].EIGN,
				   4, left, right, states);

			newviewFlexGamma(tInfo->tipCase,
					 x1_start, x2_start, x3_start,
					 tr->partitionData[model].EV,
					 tr->partitionData[model].tipVector,
					 ex3, tipX1, tipX2,
					 width, left, right, wgt, &scalerIncrement, tr->useFastScaling, states);		
		      }
		      break;
		    default:
		      assert(0);
		    }
		  break;
		case GENERIC_64:
		  break;
		default:
		  assert(0);
		}
	      if(tr->useFastScaling)
		{

		  tr->partitionData[model].globalScaler[tInfo->pNumber] = 
		    tr->partitionData[model].globalScaler[tInfo->qNumber] + 
		    tr->partitionData[model].globalScaler[tInfo->rNumber] +
		    (unsigned int)scalerIncrement;
		  assert(tr->partitionData[model].globalScaler[tInfo->pNumber] < INT_MAX);
		}
	    }
	}



    }
}



void newviewIterativeMulti (tree *tr)
{ 
  int 
    i, 
    model;

  size_t 
    rateHet;
  
  assert(tr->multiBranch);
  assert(!(tr->useGappedImplementation || tr->saveMemory));
  assert(tr->useFastScaling);
   
  if(tr->rateHetModel == CAT)
    rateHet = 1;
  else
    rateHet = 4;  

  for(model = 0; model < tr->NumberOfModels; model++)    
    {         
      if(tr->executeModel[model])
	{	      
	  double
	    *x1_start = (double*)NULL,
	    *x2_start = (double*)NULL,
	    *x3_start = (double*)NULL,
	    *left     = (double*)NULL,
	    *right    = (double*)NULL,
	    qz, 
	    rz;
	  
	  int
	    *ex3 = (int*)NULL,
	    scalerIncrement = 0,
	    *wgt = (int*)NULL;

	  unsigned char
	    *tipX1 = (unsigned char *)NULL,
	    *tipX2 = (unsigned char *)NULL;
	 
	  size_t 	 
	    states = (size_t)tr->partitionData[model].states,
	    width =  tr->partitionData[model].width,
	    requiredLength = width * rateHet * states * sizeof(double);	   
	  
	  traversalInfo 
	    *ti   = tr->td[model].ti;

	  if(tr->useFastScaling)		
	    wgt = tr->partitionData[model].wgt;		 		  					  	 	  

	  for(i = 1; i < tr->td[model].count; i++)
	    {    	      
	      traversalInfo 
		*tInfo = &ti[i];

	      size_t
		availableLength =  tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)];
		
	      x3_start = tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1];

	      if(requiredLength != availableLength)
		{		  
		  if(x3_start)
		    free(x3_start);
		 
		  x3_start = (double*)malloc_aligned(requiredLength);		 
		  
		  tr->partitionData[model].xVector[tInfo->pNumber - tr->mxtips - 1] = x3_start;
		  
		  tr->partitionData[model].xSpaceVector[(tInfo->pNumber - tr->mxtips - 1)] = requiredLength;		 
		}

	      switch(tInfo->tipCase)
		{
		case TIP_TIP:
		  tipX1    = tr->partitionData[model].yVector[tInfo->qNumber];
		  tipX2    = tr->partitionData[model].yVector[tInfo->rNumber];				 		 

		  
		  break;
		case TIP_INNER:
		  tipX1    =  tr->partitionData[model].yVector[tInfo->qNumber];		  
		  x2_start       = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];
		 
		  
		  break;
		case INNER_INNER:
		  x1_start       = tr->partitionData[model].xVector[tInfo->qNumber - tr->mxtips - 1];
		  x2_start       = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];		  

		 
		  break;
		default:
		  assert(0);
		}

	     
	      left  = tr->partitionData[model].left;
	      right = tr->partitionData[model].right;
		

	      if(tr->multiBranch)
		{
		  qz = tInfo->qz[model];
		  rz = tInfo->rz[model];
		}
	      else
		{
		  assert(0);
		  qz = tInfo->qz[0];
		  rz = tInfo->rz[0];
		}

	      switch(tr->partitionData[model].dataType)
		{
		case BINARY_DATA:
		  switch(tr->rateHetModel)
		    {
		    case CAT:
		      {			
			makeP(qz, rz, tr->partitionData[model].perSiteRates,   tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN, tr->partitionData[model].numberOfCategories,
			      left, right, BINARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRCAT_BINARY(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					     x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					     ex3, tipX1, tipX2,
					     width, left, right, wgt, &scalerIncrement, tr->useFastScaling
					     );		
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI, tr->partitionData[model].EIGN,
			      4, left, right, BINARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMA_BINARY(tInfo->tipCase,
					       x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					       ex3, tipX1, tipX2,
					       width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			  makeP(qz, rz, tr->partitionData[model].perSiteRates,   tr->partitionData[model].EI,
				tr->partitionData[model].EIGN, tr->partitionData[model].numberOfCategories,
				left, right, DNA_DATA, tr->saveMemory, tr->maxCategories);

			  newviewGTRCAT(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					ex3, tipX1, tipX2,
					width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
			}
		      break;
		    case GAMMA:
		    case GAMMA_I:
		     
			{			 
			  makeP(qz, rz, tr->partitionData[model].gammaRates,
				tr->partitionData[model].EI, tr->partitionData[model].EIGN,
				4, left, right, DNA_DATA, tr->saveMemory, tr->maxCategories);


			  if(tr->useGappedImplementation || tr->saveMemory)
			    assert(0);
			  /*newviewGTRGAMMA(tInfo->tipCase,
			    x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
			    ex3, tipX1, tipX2,
			    width, left, right, wgt, &scalerIncrement, tr->useFastScaling,
			    x1_gap, x2_gap, x3_gap, tr->partitionData[model].gapVectorLength, 
			    x1_gapColumn, x2_gapColumn, x3_gapColumn);
			  */
			  
			  else
			    newviewGTRGAMMA(tInfo->tipCase,
					    x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
					    ex3, tipX1, tipX2,
					    width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			  makeP(qz, rz, tr->partitionData[model].perSiteRates,
				tr->partitionData[model].EI,
				tr->partitionData[model].EIGN,
				tr->partitionData[model].numberOfCategories, left, right, AA_DATA, tr->saveMemory, tr->maxCategories);
			  
			  newviewGTRCATPROT(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					    ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);			
			}
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      
		      
			{
			  makeP(qz, rz, tr->partitionData[model].gammaRates,
				tr->partitionData[model].EI,
				tr->partitionData[model].EIGN,
				4, left, right, AA_DATA, tr->saveMemory, tr->maxCategories);

			  newviewGTRGAMMAPROT(tInfo->tipCase,
					      x1_start, x2_start, x3_start,
					      tr->partitionData[model].EV,
					      tr->partitionData[model].tipVector,
					      ex3, tipX1, tipX2,
					      width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			makeP(qz, rz, tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA_6, tr->saveMemory, tr->maxCategories);

			newviewGTRCATSECONDARY_6(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      4, left, right, SECONDARY_DATA_6, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMASECONDARY_6(tInfo->tipCase,
						   x1_start, x2_start, x3_start,
						   tr->partitionData[model].EV,
						   tr->partitionData[model].tipVector,
						   ex3, tipX1, tipX2,
						   width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			makeP(qz, rz, tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA_7, tr->saveMemory, tr->maxCategories);

			newviewGTRCATSECONDARY_7(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
						 x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      4, left, right, SECONDARY_DATA_7, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMASECONDARY_7(tInfo->tipCase,
						   x1_start, x2_start, x3_start,
						   tr->partitionData[model].EV,
						   tr->partitionData[model].tipVector,
						   ex3, tipX1, tipX2,
						   width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
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
			makeP(qz, rz, tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRCATSECONDARY(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
					       ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
		      }
		      break;
		    case GAMMA:
		    case GAMMA_I:
		      {
			makeP(qz, rz, tr->partitionData[model].gammaRates,
			      tr->partitionData[model].EI,
			      tr->partitionData[model].EIGN,
			      4, left, right, SECONDARY_DATA, tr->saveMemory, tr->maxCategories);

			newviewGTRGAMMASECONDARY(tInfo->tipCase,
						 x1_start, x2_start, x3_start,
						 tr->partitionData[model].EV,
						 tr->partitionData[model].tipVector,
						 ex3, tipX1, tipX2,
						 width, left, right, wgt, &scalerIncrement, tr->useFastScaling);		
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

		  tr->partitionData[model].globalScaler[tInfo->pNumber] = 
		    tr->partitionData[model].globalScaler[tInfo->qNumber] + 
		    tr->partitionData[model].globalScaler[tInfo->rNumber] +
		    (unsigned int)scalerIncrement;
		  assert(tr->partitionData[model].globalScaler[tInfo->pNumber] < INT_MAX);
		}
	    }
	}
    }
}


void newviewGeneric (tree *tr, nodeptr p)
{  
  if(isTip(p->number, tr->mxtips))
    return;

  if(tr->multiGene)
    {	           
      int i;
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  if(tr->executeModel[i])
	    {
	      tr->td[i].count = 1; 
	      computeTraversalInfoMulti(p, &(tr->td[i].ti[0]), &(tr->td[i].count), tr->mxtips, i); 
	    }
	}
      /* if(tr->td[i].count > 1)*/
      newviewIterativeMulti(tr);
    }
  else
    {
      tr->td[0].count = 1;
      computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
      
      if(tr->td[0].count > 1)
	{
#ifdef _USE_PTHREADS
	  masterBarrier(THREAD_NEWVIEW, tr);
#else
	  newviewIterative(tr);
#endif
	}
    }
}



void newviewGenericMulti (tree *tr, nodeptr p, int model)
{  
  assert(tr->multiGene);

  assert(!isTip(p->number, tr->mxtips)); 
  
  {	           
    int i;
    
    assert(p->backs[model]);

    for(i = 0; i < tr->NumberOfModels; i++)
      tr->executeModel[i] = FALSE;
    tr->executeModel[model] = TRUE;

   	   
    tr->td[model].count = 1; 
    computeTraversalInfoMulti(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->mxtips, model);    
    
    /* if(tr->td[i].count > 1)*/
    newviewIterativeMulti(tr);

   for(i = 0; i < tr->NumberOfModels; i++)
      tr->executeModel[i] = TRUE; 
  }
}

void newviewGenericMasked(tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips))
    return;

  {
    int i;

    for(i = 0; i < tr->NumberOfModels; i++)
      {
	if(tr->partitionConverged[i])
	  tr->executeModel[i] = FALSE;
	else
	  tr->executeModel[i] = TRUE;
      }
    
    if(tr->multiGene)
      {
	for(i = 0; i < tr->NumberOfModels; i++)
	  {
	    if(tr->executeModel[i])
	      {
		tr->td[i].count = 1; 
		computeTraversalInfoMulti(p, &(tr->td[i].ti[0]), &(tr->td[i].count), tr->mxtips, i); 
	      }
	    else
	      tr->td[i].count = 0; 
	  }
	/* if(tr->td[i].count > 1) ? */
	newviewIterativeMulti(tr);
      }
    else
      {
	tr->td[0].count = 1;
	computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);

	if(tr->td[0].count > 1)
	  {
#ifdef _USE_PTHREADS
	    masterBarrier(THREAD_NEWVIEW_MASKED, tr);
#else
	    newviewIterative(tr);
#endif
	  }
      }

    for(i = 0; i < tr->NumberOfModels; i++)
      tr->executeModel[i] = TRUE;
  }
}


#ifdef  _USE_PTHREADS

static void newviewMultiGrain(tree *tr,  double *x1, double *x2, double *x3, int *_ex1, int *_ex2, int *_ex3, unsigned char *_tipX1, unsigned char *_tipX2, 
			      int tipCase, double *_pz, double *_qz)
{
  int    
    scalerIncrement = 0,   
    model,         
    columnCounter = 0,
    offsetCounter = 0;

  

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      double
	*x1_start = (double*)NULL,
	*x2_start = (double*)NULL,
	*x3_start = (double*)NULL,
	*left     = tr->partitionData[model].left,
	*right    = tr->partitionData[model].right,
	pz, qz;

      int
	*wgt          = &tr->contiguousWgt[columnCounter],
	*rateCategory = &tr->contiguousRateCategory[columnCounter],
	*ex1 = (int*)NULL,
	*ex2 = (int*)NULL,
	*ex3 = (int*)NULL,
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;

      unsigned char
	*tipX1 = (unsigned char *)NULL,
	*tipX2 = (unsigned char *)NULL;     

      switch(tipCase)
	{
	case TIP_TIP:    
	  tipX1 =    &_tipX1[columnCounter];
	  tipX2 =    &_tipX2[columnCounter];
	  ex3   =    &_ex3[columnCounter];
	  x3_start = &x3[offsetCounter];
	 
	  if(!tr->useFastScaling)
	    {
	      int k;	      
	      
	      for(k = 0; k < width; k++)
		ex3[k] = 0;
	    }
	  break;
	case TIP_INNER:
	  tipX1 =    &_tipX1[columnCounter];

	  ex2   =    &_ex2[columnCounter];
	  x2_start = &x2[offsetCounter];
	  
	  ex3   =    &_ex3[columnCounter];
	  x3_start = &x3[offsetCounter];
	 
	  if(!tr->useFastScaling)
	    {
	      int k;	      
	      
	      for(k = 0; k < width; k++)
		ex3[k] = ex2[k];
	    }
	  break;
	case INNER_INNER:
	  ex1   =    &_ex1[columnCounter];
	  x1_start = &x1[offsetCounter];

	  ex2   =    &_ex2[columnCounter];
	  x2_start = &x2[offsetCounter];
	  
	  ex3   =    &_ex3[columnCounter];
	  x3_start = &x3[offsetCounter];
	 
	  if(!tr->useFastScaling)
	    {
	      int k;	      
	      
	      for(k = 0; k < width; k++)
		ex3[k] = ex1[k] + ex2[k];
	    }
	  break;
	default:
	  assert(0);
	}    

      if(tr->multiBranch)
	{
	  pz = _pz[model];
	  pz = (pz > zmin) ? log(pz) : log(zmin);
	  
	  qz = _qz[model];
	  qz = (qz > zmin) ? log(qz) : log(zmin);
	}
      else
	{	  
	  pz = _pz[0];
	  pz = (pz > zmin) ? log(pz) : log(zmin);

	  qz = _qz[0];
	  qz = (qz > zmin) ? log(qz) : log(zmin);
	}

      

      switch(tr->partitionData[model].dataType)
	{
	case BINARY_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	      
	      makeP(pz, qz, tr->partitionData[model].perSiteRates,   tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN, tr->partitionData[model].numberOfCategories,
		    left, right, BINARY_DATA, tr->saveMemory, tr->maxCategories);
	      
	      newviewGTRCAT_BINARY(tipCase,  tr->partitionData[model].EV, rateCategory,
				   x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				   ex3, tipX1, tipX2,
				   width, left, right, wgt, &scalerIncrement, tr->useFastScaling
				   );	      	      
	      break;
	    case GAMMA:
	    case GAMMA_I:	      		
	      makeP(pz, qz, tr->partitionData[model].gammaRates,
		    tr->partitionData[model].EI, tr->partitionData[model].EIGN,
		    4, left, right, BINARY_DATA, tr->saveMemory, tr->maxCategories);
	      
	      newviewGTRGAMMA_BINARY(tipCase,
				     x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
				     ex3, tipX1, tipX2,
				     width, left, right, wgt, &scalerIncrement, tr->useFastScaling);	      	      
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case DNA_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	      			
	      makeP(pz, qz, tr->partitionData[model].perSiteRates,   tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN, tr->partitionData[model].numberOfCategories,
		    left, right, DNA_DATA, tr->saveMemory, tr->maxCategories);

#ifdef __AVX
	      newviewGTRCAT_AVX(tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
				x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				ex3, tipX1, tipX2,
				width, left, right, wgt, &scalerIncrement, tr->useFastScaling); 
#else
	      
	      newviewGTRCAT(tipCase,  tr->partitionData[model].EV, rateCategory,
			    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
			    ex3, tipX1, tipX2,
			    width, left, right, wgt, &scalerIncrement, tr->useFastScaling
			    );
#endif
	      	      
	      break;
	    case GAMMA:
	    case GAMMA_I:	     		
	      makeP(pz, qz, tr->partitionData[model].gammaRates,
		    tr->partitionData[model].EI, tr->partitionData[model].EIGN,
		    4, left, right, DNA_DATA, tr->saveMemory, tr->maxCategories);
	      
#ifdef __AVX
	      newviewGTRGAMMA_AVX(tipCase,
				  x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
				  ex3, tipX1, tipX2,
				  width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else
	      newviewGTRGAMMA(tipCase,
			      x1_start, x2_start, x3_start, tr->partitionData[model].EV, tr->partitionData[model].tipVector,
			      ex3, tipX1, tipX2,
			      width, left, right, wgt, &scalerIncrement, tr->useFastScaling);			      	     
#endif
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case AA_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	      
	      makeP(pz, qz, tr->partitionData[model].perSiteRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    tr->partitionData[model].numberOfCategories, left, right, AA_DATA, tr->saveMemory, tr->maxCategories);
#ifdef __AVX
	      newviewGTRCATPROT_AVX(tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
				    x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				    ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else		      
	      newviewGTRCATPROT(tipCase,  tr->partitionData[model].EV, rateCategory,
				x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);	      	      
#endif
	      break;
	    case GAMMA:
	    case GAMMA_I:	      
	      makeP(pz, qz, tr->partitionData[model].gammaRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    4, left, right, AA_DATA, tr->saveMemory, tr->maxCategories);
#ifdef __AVX
	      newviewGTRGAMMAPROT_AVX(tipCase,
				      x1_start, x2_start, x3_start,
				      tr->partitionData[model].EV,
				      tr->partitionData[model].tipVector,
				      ex3, tipX1, tipX2,
				      width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#else	      
	      newviewGTRGAMMAPROT(tipCase,
				  x1_start, x2_start, x3_start,
				  tr->partitionData[model].EV,
				  tr->partitionData[model].tipVector,
				  ex3, tipX1, tipX2,
				  width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
#endif
	      	      
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA:
	  switch(tr->rateHetModel)
	    {
	    case CAT:	      
	      makeP(pz, qz, tr->partitionData[model].perSiteRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA, tr->saveMemory, tr->maxCategories);

	      newviewGTRCATSECONDARY(tipCase,  tr->partitionData[model].EV, rateCategory,
				     x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				     ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
	      	      
	      break;
	    case GAMMA:
	    case GAMMA_I:	      
	      makeP(pz, qz, tr->partitionData[model].gammaRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    4, left, right, SECONDARY_DATA, tr->saveMemory, tr->maxCategories);
	      
	      newviewGTRGAMMASECONDARY(tipCase,
				       x1_start, x2_start, x3_start,
				       tr->partitionData[model].EV,
				       tr->partitionData[model].tipVector,
				       ex3, tipX1, tipX2,
				       width, left, right, wgt, &scalerIncrement, tr->useFastScaling);
	    
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA_6:
	  switch(tr->rateHetModel)
	    {
	    case CAT:		      
	      makeP(pz, qz, tr->partitionData[model].perSiteRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA_6, tr->saveMemory, tr->maxCategories);
	      
	      newviewGTRCATSECONDARY_6(tipCase,  tr->partitionData[model].EV, rateCategory,
				       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				       ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);	     
	      break;
	    case GAMMA:
	    case GAMMA_I:		      
	      makeP(pz, qz, tr->partitionData[model].gammaRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    4, left, right, SECONDARY_DATA_6, tr->saveMemory, tr->maxCategories);

	      newviewGTRGAMMASECONDARY_6(tipCase,
					 x1_start, x2_start, x3_start,
					 tr->partitionData[model].EV,
					 tr->partitionData[model].tipVector,
					 ex3, tipX1, tipX2,
					 width, left, right, wgt, &scalerIncrement, tr->useFastScaling);	    
	      break;
	    default:
	      assert(0);
	    }
	  break;
	case SECONDARY_DATA_7:
	  switch(tr->rateHetModel)
	    {
	    case CAT:		      
	      makeP(pz, qz, tr->partitionData[model].perSiteRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    tr->partitionData[model].numberOfCategories, left, right, SECONDARY_DATA_7, tr->saveMemory, tr->maxCategories);
	      
	      newviewGTRCATSECONDARY_7(tipCase,  tr->partitionData[model].EV, rateCategory,
				       x1_start, x2_start, x3_start, tr->partitionData[model].tipVector,
				       ex3, tipX1, tipX2, width, left, right, wgt, &scalerIncrement, tr->useFastScaling);	      
	      break;
	    case GAMMA:
	    case GAMMA_I:	      
	      makeP(pz, qz, tr->partitionData[model].gammaRates,
		    tr->partitionData[model].EI,
		    tr->partitionData[model].EIGN,
		    4, left, right, SECONDARY_DATA_7, tr->saveMemory, tr->maxCategories);

	      newviewGTRGAMMASECONDARY_7(tipCase,
					 x1_start, x2_start, x3_start,
					 tr->partitionData[model].EV,
					 tr->partitionData[model].tipVector,
					 ex3, tipX1, tipX2,
					 width, left, right, wgt, &scalerIncrement, tr->useFastScaling);	      
	      break;
	    default:
	      assert(0);
	    }
	  break;
	default:
	  assert(0);
	}
      columnCounter += width;
      offsetCounter += width * tr->partitionData[model].states * tr->discreteRateCategories;
    }
}



void newviewClassify(tree *tr, branchInfo *b, double *z)
{
  int 
    leftNumber = b->epa->leftNodeNumber,
    rightNumber = b->epa->rightNodeNumber,
    tipCase = -1,
    *ex1 = (int*)NULL,
    *ex2 = (int*)NULL,
    *ex3 = tr->temporaryScaling;

  double
    *x1_start = (double*)NULL,
    *x2_start = (double*)NULL,
    *x3_start = tr->temporaryVector;

  unsigned char
    *tipX1 = (unsigned char*)NULL,
    *tipX2 = (unsigned char*)NULL;
  
  if (isTip(leftNumber, tr->mxtips) && isTip(rightNumber, tr->mxtips))
    {	  
      tipCase = TIP_TIP;
      
      tipX1 = tr->contiguousTips[leftNumber];
      tipX2 = tr->contiguousTips[rightNumber];
    }
  else
    {
      if (isTip(leftNumber, tr->mxtips))
	{	      
	  tipCase = TIP_INNER;
	  
	  tipX1 = tr->contiguousTips[leftNumber];

	  x2_start = b->epa->right;
	  ex2      = b->epa->rightScaling;	     	  
	}
      else
	{
	  if(isTip(rightNumber, tr->mxtips))
	    {		  
	      tipCase = TIP_INNER;
	      
	      tipX1 = tr->contiguousTips[rightNumber];
	      
	      x2_start = b->epa->left;
	      ex2      = b->epa->leftScaling;		 
	    }
	  else
	    {
	      tipCase = INNER_INNER;
	      
	      x1_start = b->epa->left;
	      ex1      = b->epa->leftScaling;
	      
	      x2_start = b->epa->right;
	      ex2      = b->epa->rightScaling;	      	    
	    }
	}
    }
    
  newviewMultiGrain(tr,  x1_start, x2_start, x3_start, ex1, ex2, ex3, tipX1, tipX2, 
		    tipCase, z, z);

}



void newviewClassifySpecial(tree *tr, double *x1_start, double *x2_start, double *x3_start, int *ex1, int *ex2, int *ex3,
			    unsigned char *tipX1,  unsigned char *tipX2, int tipCase, double *pz, double *qz)
{
  newviewMultiGrain(tr,  x1_start, x2_start, x3_start, ex1, ex2, ex3, tipX1, tipX2, 
		    tipCase, pz, qz);
  
}




#endif


