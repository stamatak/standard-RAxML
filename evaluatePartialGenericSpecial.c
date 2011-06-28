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
#endif



/********************** GTRCAT ***************************************/


#ifdef __SIM_SSE3

static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					   traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					   unsigned  char **yVector, int mxtips)
{       
  double   *x1, *x2, *x3;  
  int
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[20 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(  lVector[20 * (rNumber - mxtips)]);                    
      break;
    case INNER_INNER:            
      x1 = &(lVector[20 * (qNumber - mxtips)]);
      x2 = &(lVector[20 * (rNumber - mxtips)]);                 
      break;    
    default:
      assert(0);
    }
     
  {
    double  
      e1[20] __attribute__ ((aligned (16))),
      e2[20] __attribute__ ((aligned (16))),
      d1[20] __attribute__ ((aligned (16))), 
      d2[20] __attribute__ ((aligned (16))), 
      lz1, lz2;  
    int l, k, scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    e1[0] = 1.0;
    e2[0] = 1.0;
    
    for(l = 1; l < 20; l++)
      {
	e1[l] = EXP(EIGN[l - 1] * lz1);
	e2[l] = EXP(EIGN[l - 1] * lz2);
      }

    for(l = 0; l < 20; l+=2)
      {
	__m128d d1v = _mm_mul_pd(_mm_load_pd(&x1[l]), _mm_load_pd(&e1[l]));
	__m128d d2v = _mm_mul_pd(_mm_load_pd(&x2[l]), _mm_load_pd(&e2[l]));
	
	_mm_store_pd(&d1[l], d1v);
	_mm_store_pd(&d2[l], d2v);	
      }

    __m128d zero = _mm_setzero_pd();

    for(l = 0; l < 20; l+=2)
      _mm_store_pd(&x3[l], zero);
                
    for(l = 0; l < 20; l++)
      { 	      
	double *ev = &EV[l * 20];
	__m128d ump_x1v = _mm_setzero_pd();
	__m128d ump_x2v = _mm_setzero_pd();
	__m128d x1px2v;

	for(k = 0; k < 20; k+=2)
	  {       
	    __m128d eiv = _mm_load_pd(&EI[20 * l + k]);
	    __m128d d1v = _mm_load_pd(&d1[k]);
	    __m128d d2v = _mm_load_pd(&d2[k]);
	    ump_x1v = _mm_add_pd(ump_x1v, _mm_mul_pd(d1v, eiv));
	    ump_x2v = _mm_add_pd(ump_x2v, _mm_mul_pd(d2v, eiv));	  
	  }

	ump_x1v = _mm_hadd_pd(ump_x1v, ump_x1v);
	ump_x2v = _mm_hadd_pd(ump_x2v, ump_x2v);

	x1px2v = _mm_mul_pd(ump_x1v, ump_x2v);

	for(k = 0; k < 20; k+=2)
	  {
	    __m128d ex3v = _mm_load_pd(&x3[k]);
	    __m128d EVV  = _mm_load_pd(&ev[k]);
	    ex3v = _mm_add_pd(ex3v, _mm_mul_pd(x1px2v, EVV));
	    
	    _mm_store_pd(&x3[k], ex3v);	   	   
	  }
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < 20); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	      
	__m128d twoto = _mm_set_pd(twotothe256, twotothe256);

	for(l = 0; l < 20; l+=2)
	  {
	    __m128d ex3v = _mm_mul_pd(_mm_load_pd(&x3[l]),twoto);
	    _mm_store_pd(&x3[l],_mm_mul_pd(ex3v, twoto));	
	  }
 	
	/*
	  for(l = 0; l < 20; l++)
	  x3[l] *= twotothe256;		   
	*/

	*eVector = *eVector + 1;
      }
    
    return;      
  }
}

static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
					int w, double *EIGN, double *EI, double *EV,
					double *tipVector, unsigned char **yVector, 
					int branchReference, int mxtips)
{
  double lz, term;       
  double  d[20];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double 
    *lVector = (double *)malloc_aligned(sizeof(double) * 20 * mxtips),
    myEI[400]  __attribute__ ((aligned (16)));

  traversalInfo *trav = &ti[0];

  

  for(k = 0; k < 20; k++)
    {
      myEI[k * 20] = 1.0;
      for(l = 1; l < 20; l++)
	myEI[k * 20 + l] = EI[k * 19 + l - 1];
    }

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCATPROT(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
			    &ti[k], EIGN, myEI, EV, 
			    tipVector, yVector, mxtips);       
   
  x2 = &lVector[20 * (trav->qNumber - mxtips)];

       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < 20; l++)
    d[l] = EXP (EIGN[l-1] * lz);

  term = 0.0;
  
  for(l = 0; l < 20; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  
 
  return  term;
}


#else

static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					   traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					   unsigned  char **yVector, int mxtips)
{       
  double   *x1, *x2, *x3;  
  int
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[20 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(  lVector[20 * (rNumber - mxtips)]);                    
      break;
    case INNER_INNER:            
      x1 = &(lVector[20 * (qNumber - mxtips)]);
      x2 = &(lVector[20 * (rNumber - mxtips)]);                 
      break;    
    default:
      assert(0);
    }
     
  {
    double  d1[20], d2[20], ump_x1, ump_x2, x1px2, lz1, lz2;  
    int l, k, scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    for(l = 1; l < 20; l++)
      {
	d1[l] = x1[l] * EXP(EIGN[l - 1] * lz1);
	d2[l] = x2[l] * EXP(EIGN[l - 1] * lz2);
      }

    for(l = 0; l < 20; l++)
      x3[l] = 0.0;
           
    for(l = 0; l < 20; l++)
      {
	ump_x1 = x1[0];
	ump_x2 = x2[0];

	for(k = 1; k < 20; k++)
	  {       
	    ump_x1 += d1[k] * EI[19 * l + k-1];
	    ump_x2 += d2[k] * EI[19 * l + k-1];
	  }

	x1px2 = ump_x1 * ump_x2;

	for(k = 0; k < 20; k++)
	  x3[k] += x1px2 * EV[l * 20 + k];	
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < 20); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	
	for(l = 0; l < 20; l++)
	  x3[l] *= twotothe256;		   
	*eVector = *eVector + 1;
      }
    
    return;      
  }
}

static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
					int w, double *EIGN, double *EI, double *EV,
					double *tipVector, unsigned char **yVector, 
					int branchReference, int mxtips)
{
  double lz, term;       
  double  d[20];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double *lVector = (double *)malloc_aligned(sizeof(double) * 20 * mxtips);

  traversalInfo *trav = &ti[0];

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCATPROT(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
			    &ti[k], EIGN, EI, EV, 
			    tipVector, yVector, mxtips);       
   
  x2 = &lVector[20 * (trav->qNumber - mxtips)];

       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < 20; l++)
    d[l] = EXP (EIGN[l-1] * lz);

  term = 0.0;
  
  for(l = 0; l < 20; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  

  return  term;
}

#endif

static inline void computeVectorGTRCATSECONDARY(double *lVector, int *eVector, double ki, int i, double qz, double rz,
						traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
						unsigned  char **yVector, int mxtips)
{       
  double   *x1, *x2, *x3;  
  int 
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[16 * (pNumber  - mxtips)]);    

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[16 * yVector[qNumber][i]]);
      x2 = &(tipVector[16 * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[16 * yVector[qNumber][i]]);
      x2 = &(  lVector[16 * (rNumber - mxtips)]);                   
      break;
    case INNER_INNER:            
      x1 = &(lVector[16 * (qNumber - mxtips)]);
      x2 = &(lVector[16 * (rNumber - mxtips)]);                
      break;    
    default:
      assert(0);
    }
     
  {
    double  d1[16], d2[16], ump_x1, ump_x2, x1px2, lz1, lz2;  
    int l, k, scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    for(l = 1; l < 16; l++)
      {
	d1[l] = x1[l] * EXP(EIGN[l - 1] * lz1);
	d2[l] = x2[l] * EXP(EIGN[l - 1] * lz2);
      }

    for(l = 0; l < 16; l++)
      x3[l] = 0.0;
           
    for(l = 0; l < 16; l++)
      {
	ump_x1 = x1[0];
	ump_x2 = x2[0];

	for(k = 1; k < 16; k++)
	  {       
	    ump_x1 += d1[k] * EI[15 * l + k-1];
	    ump_x2 += d2[k] * EI[15 * l + k-1];
	  }

	x1px2 = ump_x1 * ump_x2;

	for(k = 0; k < 16; k++)
	  x3[k] += x1px2 * EV[l * 16 + k];	
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < 16); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	
	for(l = 0; l < 16; l++)
	  x3[l] *= twotothe256;		   
	*eVector = *eVector + 1;
      }
    
    return;      
  }
}


static inline void computeVectorFlex(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				     traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				     unsigned  char **yVector, int mxtips, const int numStates)
{       
  double   *x1, *x2, *x3;  
  int 
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[numStates * (pNumber  - mxtips)]);    

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[numStates * yVector[qNumber][i]]);
      x2 = &(tipVector[numStates * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[numStates * yVector[qNumber][i]]);
      x2 = &(  lVector[numStates * (rNumber - mxtips)]);                   
      break;
    case INNER_INNER:            
      x1 = &(lVector[numStates * (qNumber - mxtips)]);
      x2 = &(lVector[numStates * (rNumber - mxtips)]);                
      break;    
    default:
      assert(0);
    }
     
  {
    double  d1[64], d2[64], ump_x1, ump_x2, x1px2, lz1, lz2;  
    int l, k, scale;
    const int rates = numStates - 1;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    for(l = 1; l < numStates; l++)
      {
	d1[l] = x1[l] * EXP(EIGN[l - 1] * lz1);
	d2[l] = x2[l] * EXP(EIGN[l - 1] * lz2);
      }

    for(l = 0; l < numStates; l++)
      x3[l] = 0.0;
           
    for(l = 0; l < numStates; l++)
      {
	ump_x1 = x1[0];
	ump_x2 = x2[0];

	for(k = 1; k < numStates; k++)
	  {       
	    ump_x1 += d1[k] * EI[rates * l + k-1];
	    ump_x2 += d2[k] * EI[rates * l + k-1];
	  }

	x1px2 = ump_x1 * ump_x2;

	for(k = 0; k < numStates; k++)
	  x3[k] += x1px2 * EV[l * numStates + k];	
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < numStates); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	
	for(l = 0; l < numStates; l++)
	  x3[l] *= twotothe256;		   
	*eVector = *eVector + 1;
      }
    
    return;      
  }
}


static inline void computeVectorGTRCATSECONDARY_6(double *lVector, int *eVector, double ki, int i, double qz, double rz,
						  traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
						  unsigned  char **yVector, int mxtips)
{       
  double   *x1, *x2, *x3;  
  int 
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[6 * (pNumber  - mxtips)]);     

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[6 * yVector[qNumber][i]]);
      x2 = &(tipVector[6 * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[6 * yVector[qNumber][i]]);
      x2 = &(  lVector[6 * (rNumber - mxtips)]);                    
      break;
    case INNER_INNER:            
      x1 = &(lVector[6 * (qNumber - mxtips)]);
      x2 = &(lVector[6 * (rNumber - mxtips)]);                
      break;    
    default:
      assert(0);
    }
     
  {
    double  d1[6], d2[6], ump_x1, ump_x2, x1px2, lz1, lz2;  
    int l, k, scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    for(l = 1; l < 6; l++)
      {
	d1[l] = x1[l] * EXP(EIGN[l - 1] * lz1);
	d2[l] = x2[l] * EXP(EIGN[l - 1] * lz2);
      }

    for(l = 0; l < 6; l++)
      x3[l] = 0.0;
           
    for(l = 0; l < 6; l++)
      {
	ump_x1 = x1[0];
	ump_x2 = x2[0];

	for(k = 1; k < 6; k++)
	  {       
	    ump_x1 += d1[k] * EI[5 * l + k-1];
	    ump_x2 += d2[k] * EI[5 * l + k-1];
	  }

	x1px2 = ump_x1 * ump_x2;

	for(k = 0; k < 6; k++)
	  x3[k] += x1px2 * EV[l * 6 + k];	
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < 6); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	
	for(l = 0; l < 6; l++)
	  x3[l] *= twotothe256;		   
	*eVector = *eVector + 1;
      }
    
    return;      
  }
}


static inline void computeVectorGTRCATSECONDARY_7(double *lVector, int *eVector, double ki, int i, double qz, double rz,
						  traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
						  unsigned  char **yVector, int mxtips)
{       
  double   *x1, *x2, *x3;  
  int
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[7 * (pNumber  - mxtips)]);  

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[7 * yVector[qNumber][i]]);
      x2 = &(tipVector[7 * yVector[rNumber][i]]);     
      break;
    case TIP_INNER:     
      x1 = &(tipVector[7 * yVector[qNumber][i]]);
      x2 = &(  lVector[7 * (rNumber - mxtips)]);                    
      break;
    case INNER_INNER:            
      x1 = &(lVector[7 * (qNumber - mxtips)]);
      x2 = &(lVector[7 * (rNumber - mxtips)]);              
      break;    
    default:
      assert(0);
    }
     
  {
    double  d1[7], d2[7], ump_x1, ump_x2, x1px2, lz1, lz2;  
    int l, k, scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;        

    for(l = 1; l < 7; l++)
      {
	d1[l] = x1[l] * EXP(EIGN[l - 1] * lz1);
	d2[l] = x2[l] * EXP(EIGN[l - 1] * lz2);
      }

    for(l = 0; l < 7; l++)
      x3[l] = 0.0;
           
    for(l = 0; l < 7; l++)
      {
	ump_x1 = x1[0];
	ump_x2 = x2[0];

	for(k = 1; k < 7; k++)
	  {       
	    ump_x1 += d1[k] * EI[6 * l + k-1];
	    ump_x2 += d2[k] * EI[6 * l + k-1];
	  }

	x1px2 = ump_x1 * ump_x2;

	for(k = 0; k < 7; k++)
	  x3[k] += x1px2 * EV[l * 7 + k];	
      }                      
    
    scale = 1;
    for(l = 0; scale && (l < 7); l++)
      scale = ((x3[l] < minlikelihood) && (x3[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	
	for(l = 0; l < 7; l++)
	  x3[l] *= twotothe256;		   
	*eVector = *eVector + 1;
      }
    
    return;      
  }
}


static inline void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       unsigned char **yVector, int mxtips)
{       
  double  d1[3], d2[3],  ump_x1, ump_x2, x1px2[4], lz1, lz2; 
  double *x1, *x2, *x3;
  int j, k,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[4 * (pNumber  - mxtips)];  
 

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &(tipVector[4 * yVector[rNumber][i]]);    
      break;
    case TIP_INNER:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &lVector[4 * (rNumber - mxtips)];           
      break;
    case INNER_INNER:            
      x1 = &lVector[4 * (qNumber - mxtips)];
      x2 = &lVector[4 * (rNumber - mxtips)];     
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
  for(j = 0; j < 3; j++)
    {
      d1[j] = 
	x1[j + 1] * 
	EXP(EIGN[j] * lz1);
      d2[j] = x2[j + 1] * EXP(EIGN[j] * lz2);	    
    }
 
 
  for(j = 0; j < 4; j++)
    {     
      ump_x1 = x1[0];
      ump_x2 = x2[0];
      for(k = 0; k < 3; k++)
	{
	  ump_x1 += d1[k] * EI[j * 3 + k];
	  ump_x2 += d2[k] * EI[j * 3 + k];
	}
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < 4; j++)
    x3[j] = 0.0;

  for(j = 0; j < 4; j++)          
    for(k = 0; k < 4; k++)	
      x3[k] +=  x1px2[j] *  EV[4 * j + k];	   
      
  
  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
      x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
      x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
      x3[3] < minlikelihood && x3[3] > minusminlikelihood)
    {	     
      x3[0]   *= twotothe256;
      x3[1]   *= twotothe256;
      x3[2]   *= twotothe256;     
      x3[3]   *= twotothe256;     
      *eVector = *eVector + 1;
    }	              

  return;
}






static inline void computeVectorGTRCAT_BINARY(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					      traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					      unsigned char **yVector, int mxtips)
{       
  double  d1, d2,  ump_x1, ump_x2, x1px2[2], lz1, lz2; 
  double *x1, *x2, *x3;
  int 
    j, k,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[2 * (pNumber  - mxtips)];  

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[2 * yVector[qNumber][i]]);
      x2 = &(tipVector[2 * yVector[rNumber][i]]);   
      break;
    case TIP_INNER:     
      x1 = &(tipVector[2 * yVector[qNumber][i]]);
      x2 = &lVector[2 * (rNumber - mxtips)];                    
      break;
    case INNER_INNER:            
      x1 = &lVector[2 * (qNumber - mxtips)];
      x2 = &lVector[2 * (rNumber - mxtips)];               
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
 
  d1 = x1[1] * EXP(EIGN[0] * lz1);
  d2 = x2[1] * EXP(EIGN[0] * lz2);	        
 
  for(j = 0; j < 2; j++)
    {     
      ump_x1 = x1[0];
      ump_x2 = x2[0];
      
      ump_x1 += d1 * EI[j];
      ump_x2 += d2 * EI[j];
	
      x1px2[j] = ump_x1 * ump_x2;
    }
  
  for(j = 0; j < 2; j++)
    x3[j] = 0.0;

  for(j = 0; j < 2; j++)          
    for(k = 0; k < 2; k++)	
      x3[k] +=  x1px2[j] *  EV[2 * j + k];	   
      
  
  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
      x3[1] < minlikelihood && x3[1] > minusminlikelihood
      )
    {	     
      x3[0]   *= twotothe256;
      x3[1]   *= twotothe256;     
      *eVector = *eVector + 1;
    }	              

  return;
}

static double evaluatePartialGTRCAT_BINARY(int i, double ki, int counter,  traversalInfo *ti, double qz,
					   int w, double *EIGN, double *EI, double *EV,
					   double *tipVector, unsigned  char **yVector, 
					   int branchReference, int mxtips)
{
  double lz, term;       
  double  d;
  double   *x1, *x2; 
  int scale = 0, k;
  double *lVector = (double *)malloc(sizeof(double) * 2 * mxtips);  
  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[2 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCAT_BINARY(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], &ti[k], 
			       EIGN, EI, EV, 
			       tipVector, yVector, mxtips);       
   
  x2 = &lVector[2 * (trav->qNumber - mxtips)];

     

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d = EXP (EIGN[0] * lz);
  
  term =  x1[0] * x2[0];
  term += x1[1] * x2[1] * d; 

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  

  return  term;
}


static double evaluatePartialGTRCAT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, unsigned  char **yVector, 
				    int branchReference, int mxtips)
{
  double lz, term;       
  double  d[3];
  double   *x1, *x2; 
  int scale = 0, k;
  double *lVector = (double *)malloc_aligned(sizeof(double) * 4 * mxtips);    

  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[4 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCAT(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], &ti[k], 
			EIGN, EI, EV, 
			tipVector, yVector, mxtips);       
   
  x2 = &lVector[4 * (trav->qNumber - mxtips)]; 

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d[0] = EXP (EIGN[0] * lz);
  d[1] = EXP (EIGN[1] * lz);
  d[2] = EXP (EIGN[2] * lz);       	   
  
  term =  x1[0] * x2[0];
  term += x1[1] * x2[1] * d[0];
  term += x1[2] * x2[2] * d[1];
  term += x1[3] * x2[3] * d[2];     

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);  

  return  term;
}




static double evaluatePartialGTRCATSECONDARY(int i, double ki, int counter,  traversalInfo *ti, double qz,
					     int w, double *EIGN, double *EI, double *EV,
					     double *tipVector, unsigned char **yVector, 
					     int branchReference, int mxtips)
{
  double lz, term;       
  double  d[16];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double *lVector = (double *)malloc(sizeof(double) * 16 * mxtips);
 
  traversalInfo *trav = &ti[0];

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[16 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCATSECONDARY(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
				 &ti[k], EIGN, EI, EV, 
				 tipVector, yVector, mxtips);       
   
  x2 = &lVector[16 * (trav->qNumber - mxtips)];

      

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < 16; l++)
    d[l] = EXP (EIGN[l-1] * lz);

  term = 0.0;
  
  for(l = 0; l < 16; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  

  return  term;
}


static double evaluatePartialFlex(int i, double ki, int counter,  traversalInfo *ti, double qz,
				  int w, double *EIGN, double *EI, double *EV,
				  double *tipVector, unsigned char **yVector, 
				  int branchReference, int mxtips, const int numStates)
{
  double lz, term;       
  double  d[64];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double *lVector = (double *)malloc(sizeof(double) * numStates * mxtips);
 
  traversalInfo *trav = &ti[0];

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[numStates *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorFlex(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
		      &ti[k], EIGN, EI, EV, 
		      tipVector, yVector, mxtips, numStates);       
   
  x2 = &lVector[numStates * (trav->qNumber - mxtips)];
     
  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < numStates; l++)
    d[l] = EXP (EIGN[l-1] * lz);

  term = 0.0;
  
  for(l = 0; l < numStates; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
  

  return  term;
}


static double evaluatePartialGTRCATSECONDARY_6(int i, double ki, int counter,  traversalInfo *ti, double qz,
					       int w, double *EIGN, double *EI, double *EV,
					       double *tipVector, unsigned char **yVector, 
					       int branchReference, int mxtips)
{
  double lz, term;       
  double  d[6];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double *lVector = (double *)malloc(sizeof(double) * 6 * mxtips);
  
  traversalInfo *trav = &ti[0];

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[6 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCATSECONDARY_6(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
				   &ti[k], EIGN, EI, EV, 
				   tipVector, yVector, mxtips);       
   
  x2 = &lVector[6 * (trav->qNumber - mxtips)];

 

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < 6; l++)
    d[l] = EXP (EIGN[l-1] * lz);

  term = 0.0;
  
  for(l = 0; l < 6; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
 

  return  term;
}

static double evaluatePartialGTRCATSECONDARY_7(int i, double ki, int counter,  traversalInfo *ti, double qz,
					       int w, double *EIGN, double *EI, double *EV,
					       double *tipVector, unsigned char **yVector, 
					       int branchReference, int mxtips)
{
  double lz, term;       
  double  d[7];
  double   *x1, *x2; 
  int scale = 0, k, l;
  double *lVector = (double *)malloc(sizeof(double) * 7 * mxtips);
 
  traversalInfo *trav = &ti[0];

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[7 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCATSECONDARY_7(lVector, &scale, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
				   &ti[k], EIGN, EI, EV, 
				   tipVector, yVector, mxtips);       
   
  x2 = &lVector[7 * (trav->qNumber - mxtips)];

       

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = 1.0;
  for(l = 1; l < 7; l++)
    d[l] = EXP (EIGN[l-1] * lz);

  term = 0.0;
  
  for(l = 0; l < 7; l++)
    term += x1[l] * x2[l] * d[l];   

  term = LOG(term) + (scale * LOG(minlikelihood));   

  term = term * w;

  free(lVector);
 

  return  term;
}




/*********************************************************************************************/



void computeFullTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches)
{
  if(isTip(p->number, maxTips))
    return; 

  {     
    int i;
    nodeptr q = p->next->back;
    nodeptr r = p->next->next->back;

    /* set xnode info at this point */

    p->x = 1;
    p->next->x = 0;
    p->next->next->x = 0;     

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {	  
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
	    
	    computeFullTraversalInfo(r, ti, counter, maxTips, numBranches);	
	    	   
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
	    computeFullTraversalInfo(q, ti, counter, maxTips, numBranches);	       
	    computeFullTraversalInfo(r, ti, counter, maxTips, numBranches);
	   
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

void determineFullTraversal(nodeptr p, tree *tr)
{
  nodeptr q = p->back;
  int k;

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(k = 0; k < tr->numBranches; k++)        
    tr->td[0].ti[0].qz[k] = q->z[k];    

  assert(isTip(p->number, tr->mxtips));

  tr->td[0].count = 1; 
  computeFullTraversalInfo(q, &(tr->td[0].ti[0]),  &(tr->td[0].count), tr->mxtips, tr->numBranches); 
  computeFullTraversalInfo(p, &(tr->td[0].ti[0]),  &(tr->td[0].count), tr->mxtips, tr->numBranches);
}





double evaluatePartialGeneric (tree *tr, int i, double ki, int _model)
{
  double result;
  int 
    branchReference,
    states = tr->partitionData[_model].states;
    
#ifdef _USE_PTHREADS
  int index = i; 
#else
  int index = i - tr->partitionData[_model].lower;
#endif
  
  if(tr->multiBranch)
    branchReference = _model;
  else
    branchReference = 0;

  assert(tr->rateHetModel == CAT);
    
 

  switch(tr->partitionData[_model].dataType)
    {
    case BINARY_DATA:
      result = evaluatePartialGTRCAT_BINARY(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					    tr->partitionData[_model].wgt[index],
					    tr->partitionData[_model].EIGN, 
					    tr->partitionData[_model].EI, 
					    tr->partitionData[_model].EV,
					    tr->partitionData[_model].tipVector,
					    tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      break;
    case DNA_DATA:      
      result = evaluatePartialGTRCAT(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
				     tr->partitionData[_model].wgt[index],
				     tr->partitionData[_model].EIGN, 
				     tr->partitionData[_model].EI, 
				     tr->partitionData[_model].EV,
				     tr->partitionData[_model].tipVector,
				     tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      break;
    case AA_DATA:
      result = evaluatePartialGTRCATPROT(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					 tr->partitionData[_model].wgt[index],
					 tr->partitionData[_model].EIGN, 
					 tr->partitionData[_model].EI, 
					 tr->partitionData[_model].EV,
					 tr->partitionData[_model].tipVector, 
					 tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      break;
    case SECONDARY_DATA:
       result = evaluatePartialGTRCATSECONDARY(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					       tr->partitionData[_model].wgt[index],
					       tr->partitionData[_model].EIGN, 
					       tr->partitionData[_model].EI, 
					       tr->partitionData[_model].EV,
					       tr->partitionData[_model].tipVector, 
					       tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      break;
    case SECONDARY_DATA_6:
      result = evaluatePartialGTRCATSECONDARY_6(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
						tr->partitionData[_model].wgt[index],
						tr->partitionData[_model].EIGN, 
						tr->partitionData[_model].EI, 
						tr->partitionData[_model].EV,
						tr->partitionData[_model].tipVector, 
						tr->partitionData[_model].yVector, branchReference, tr->mxtips);
      break; 
    case SECONDARY_DATA_7:
      result = evaluatePartialGTRCATSECONDARY_7(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
						tr->partitionData[_model].wgt[index],
						tr->partitionData[_model].EIGN, 
						tr->partitionData[_model].EI, 
						tr->partitionData[_model].EV,
						tr->partitionData[_model].tipVector, 
						tr->partitionData[_model].yVector, branchReference, tr->mxtips);
    case GENERIC_32:
      result = evaluatePartialFlex(index, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
				   tr->partitionData[_model].wgt[index],
				   tr->partitionData[_model].EIGN, 
				   tr->partitionData[_model].EI, 
				   tr->partitionData[_model].EV,
				   tr->partitionData[_model].tipVector, 
				   tr->partitionData[_model].yVector, branchReference, tr->mxtips, states);
      break;
    default:
      assert(0);
    }
 

  return result;
}

