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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with
 *  thousands of taxa and mixed models".
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
extern volatile double *reductionBufferTwo;
extern volatile int NumberOfThreads;
#endif

extern const unsigned int mask32[32];


/*******************/

static void sumCAT_BINARY(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
			  unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i;
  
#ifndef __SIM_SSE3
  int j;
#endif
  double *x1, *x2;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &(tipVector[2 * tipX2[i]]);

#ifndef __SIM_SSE3
	  for(j = 0; j < 2; j++)
	    sum[i * 2 + j]     = x1[j] * x2[j];
#else
	  _mm_store_pd(&sum[i * 2], _mm_mul_pd( _mm_load_pd(x1), _mm_load_pd(x2)));
#endif
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &x2_start[2 * i];

#ifndef __SIM_SSE3
	  for(j = 0; j < 2; j++)
	    sum[i * 2 + j]     = x1[j] * x2[j];
#else
	  _mm_store_pd(&sum[i * 2], _mm_mul_pd( _mm_load_pd(x1), _mm_load_pd(x2)));  
#endif
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[2 * i];
	  x2 = &x2_start[2 * i];
#ifndef __SIM_SSE3
	  for(j = 0; j < 2; j++)
	    sum[i * 2 + j]     = x1[j] * x2[j];
#else
	  _mm_store_pd(&sum[i * 2], _mm_mul_pd( _mm_load_pd(x1), _mm_load_pd(x2)));   
#endif
	}
      break;
    default:
      assert(0);
    }
}
#ifndef __SIM_SSE3
static void sumCAT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
		   unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, j;
  double *x1, *x2;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);

	  for(j = 0; j < 4; j++)
	    sum[i * 4 + j]     = x1[j] * x2[j];

	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];

	  for(j = 0; j < 4; j++)
	    sum[i * 4 + j]     = x1[j] * x2[j];

	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];


	  for(j = 0; j < 4; j++)
	    sum[i * 4 + j]     = x1[j] * x2[j];

	}
      break;
    default:
      assert(0);
    }
}

#else

static void sumCAT_SAVE(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  int i;
  double 
    *x1, 
    *x2,    
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;

  switch(tipCase)
  {
    case TIP_TIP:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        x2 = &(tipVector[4 * tipX2[i]]);

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
      {
        x1 = &(tipVector[4 * tipX1[i]]);
        if(isGap(x2_gap, i))
          x2 = x2_gapColumn;
        else
        {
          x2 = x2_ptr;
          x2_ptr += 4;
        }

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
      }
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
      {
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

        _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
        _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));

      }    
      break;
    default:
      assert(0);
  }
}


static void sumCAT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
		   unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i;
  double 
    *x1, 
    *x2;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);

	  _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
	  _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];

	  _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
	  _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];

	  _mm_store_pd( &sum[i*4 + 0], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));
	  _mm_store_pd( &sum[i*4 + 2], _mm_mul_pd( _mm_load_pd( &x1[2] ), _mm_load_pd( &x2[2] )));

	}    
      break;
    default:
      assert(0);
    }
}

#endif




static void coreGTRCAT_BINARY(int upper, int numberOfCategories, double *sum,
			      volatile double *d1, volatile double *d2, 
			      double *rptr, double *EIGN, int *cptr, double lz, int *wgt)
{
  int i;
  double
    *d, *d_start,
    tmp_0, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e[2];
  double dd1;

  e[0] = EIGN[0];
  e[1] = EIGN[0] * EIGN[0];


  d = d_start = (double *)rax_malloc(numberOfCategories * sizeof(double));

  dd1 = e[0] * lz;

  for(i = 0; i < numberOfCategories; i++)
    d[i] = EXP(dd1 * rptr[i]);

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];
      
      d = &d_start[cptr[i]];

      inv_Li = sum[2 * i];
      inv_Li += (tmp_0 = d[0] * sum[2 * i + 1]);

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   = tmp_0 * e[0];
      d2lnLidlz2 = tmp_0 * e[1];

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(d_start);
}

#ifdef __SIM_SSE3

static void coreGTRCAT(int upper, int numberOfCategories, double *sum,
			   volatile double *d1, volatile double *d2,
		       double *rptr, double *EIGN, int *cptr, double lz, int *wgt)
{
  int i;
  double
    *d, *d_start,
    inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e1[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double e2[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double dd1, dd2, dd3;

  __m128d
    e1v[2],
    e2v[2];

  e1[0] = 0.0;
  e2[0] = 0.0;
  e1[1] = EIGN[0];
  e2[1] = EIGN[0] * EIGN[0];
  e1[2] = EIGN[1];
  e2[2] = EIGN[1] * EIGN[1];
  e1[3] = EIGN[2];
  e2[3] = EIGN[2] * EIGN[2];

  e1v[0]= _mm_load_pd(&e1[0]);
  e1v[1]= _mm_load_pd(&e1[2]);

  e2v[0]= _mm_load_pd(&e2[0]);
  e2v[1]= _mm_load_pd(&e2[2]);

  d = d_start = (double *)rax_malloc(numberOfCategories * 4 * sizeof(double));

  dd1 = EIGN[0] * lz;
  dd2 = EIGN[1] * lz;
  dd3 = EIGN[2] * lz;

  for(i = 0; i < numberOfCategories; i++)
    {
      d[i * 4 + 0] = 1.0;
      d[i * 4 + 1] = EXP(dd1 * rptr[i]);
      d[i * 4 + 2] = EXP(dd2 * rptr[i]);
      d[i * 4 + 3] = EXP(dd3 * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double 
	*s = &sum[4 * i],
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];
      
      d = &d_start[4 * cptr[i]];  
      
      __m128d tmp_0v =_mm_mul_pd(_mm_load_pd(&d[0]),_mm_load_pd(&s[0]));
      __m128d tmp_1v =_mm_mul_pd(_mm_load_pd(&d[2]),_mm_load_pd(&s[2]));

      __m128d inv_Liv    = _mm_add_pd(tmp_0v, tmp_1v);      
            	  
      __m128d dlnLidlzv   = _mm_add_pd(_mm_mul_pd(tmp_0v, e1v[0]), _mm_mul_pd(tmp_1v, e1v[1]));	  
      __m128d d2lnLidlz2v = _mm_add_pd(_mm_mul_pd(tmp_0v, e2v[0]), _mm_mul_pd(tmp_1v, e2v[1]));


      inv_Liv   = _mm_hadd_pd(inv_Liv, inv_Liv);
      dlnLidlzv = _mm_hadd_pd(dlnLidlzv, dlnLidlzv);
      d2lnLidlz2v = _mm_hadd_pd(d2lnLidlz2v, d2lnLidlz2v);                 
 
      _mm_storel_pd(&inv_Li, inv_Liv);     
      _mm_storel_pd(&dlnLidlz, dlnLidlzv);                 
      _mm_storel_pd(&d2lnLidlz2, d2lnLidlz2v);      

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(d_start);
}


#else

static void coreGTRCAT(int upper, int numberOfCategories, double *sum,
			   volatile double *d1, volatile double *d2, 
		       double *rptr, double *EIGN, int *cptr, double lz, int *wgt)
{
  int i;
  double
    *d, *d_start,
    tmp_0, tmp_1, tmp_2, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e[6];
  double dd1, dd2, dd3;

  e[0] = EIGN[0];
  e[1] = EIGN[0] * EIGN[0];
  e[2] = EIGN[1];
  e[3] = EIGN[1] * EIGN[1];
  e[4] = EIGN[2];
  e[5] = EIGN[2] * EIGN[2];

  d = d_start = (double *)rax_malloc(numberOfCategories * 4 * sizeof(double));

  dd1 = e[0] * lz;
  dd2 = e[2] * lz;
  dd3 = e[4] * lz;

  for(i = 0; i < numberOfCategories; i++)
    {
      d[i * 4] = EXP(dd1 * rptr[i]);
      d[i * 4 + 1] = EXP(dd2 * rptr[i]);
      d[i * 4 + 2] = EXP(dd3 * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];

      d = &d_start[4 * cptr[i]];

      inv_Li = sum[4 * i];
      inv_Li += (tmp_0 = d[0] * sum[4 * i + 1]);
      inv_Li += (tmp_1 = d[1] * sum[4 * i + 2]);
      inv_Li += (tmp_2 = d[2] * sum[4 * i + 3]);

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   = tmp_0 * e[0];
      d2lnLidlz2 = tmp_0 * e[1];

      dlnLidlz   += tmp_1 * e[2];
      d2lnLidlz2 += tmp_1 * e[3];

      dlnLidlz   += tmp_2 * e[4];
      d2lnLidlz2 += tmp_2 * e[5];

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;


      dlnLdlz   += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(d_start);
}

#endif





#ifdef __SIM_SSE3
static void sumGTRCATPROT_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  int 
    i, 
  l;

  double 
    *sum, 
    *left, 
    *right,
    *left_ptr = x1,
    *right_ptr = x2;

  switch(tipCase)
  {
    case TIP_TIP:
      for (i = 0; i < n; i++)
      {
        left  = &(tipVector[20 * tipX1[i]]);
        right = &(tipVector[20 * tipX2[i]]);
        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);		 
        }

      }
      break;
    case TIP_INNER:
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

        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);		 
        }

      }
      break;
    case INNER_INNER:
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

        sum = &sumtable[20 * i];

        for(l = 0; l < 20; l+=2)
        {
          __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));

          _mm_store_pd(&sum[l], sumv);		 
        }
      }
      break;
    default:
      assert(0);
  }
}


#endif


static void sumGTRCATPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			  unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  left  = &(tipVector[20 * tipX1[i]]);
	  right = &(tipVector[20 * tipX2[i]]);
	  sum = &sumtable[20 * i];
#ifdef __SIM_SSE3
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
	      
	      _mm_store_pd(&sum[l], sumv);		 
	    }
#else
	  for(l = 0; l < 20; l++)
	    sum[l] = left[l] * right[l];
#endif
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  left = &(tipVector[20 * tipX1[i]]);
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
#ifdef __SIM_SSE3
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
	      
	      _mm_store_pd(&sum[l], sumv);		 
	    }
#else
	  for(l = 0; l < 20; l++)
	    sum[l] = left[l] * right[l];
#endif
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
#ifdef __SIM_SSE3
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[l]), _mm_load_pd(&right[l]));
	      
	      _mm_store_pd(&sum[l], sumv);		 
	    }
#else
	  for(l = 0; l < 20; l++)
	    sum[l] = left[l] * right[l];
#endif
	}
      break;
    default:
      assert(0);
    }
}



static void sumGTRCATSECONDARY(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			       unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  left  = &(tipVector[16 * tipX1[i]]);
	  right = &(tipVector[16 * tipX2[i]]);
	  sum = &sumtable[16 * i];

	  for(l = 0; l < 16; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  left = &(tipVector[16 * tipX1[i]]);
	  right = &x2[16 * i];
	  sum = &sumtable[16 * i];

	  for(l = 0; l < 16; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  left  = &x1[16 * i];
	  right = &x2[16 * i];
	  sum = &sumtable[16 * i];

	  for(l = 0; l < 16; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    default:
      assert(0);
    }
}


static void sumCatFlex(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
		       unsigned char *tipX1, unsigned char *tipX2, int n, const int numStates)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  left  = &(tipVector[numStates * tipX1[i]]);
	  right = &(tipVector[numStates * tipX2[i]]);
	  sum = &sumtable[numStates * i];

	  for(l = 0; l < numStates; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  left = &(tipVector[numStates * tipX1[i]]);
	  right = &x2[numStates * i];
	  sum = &sumtable[numStates * i];

	  for(l = 0; l < numStates; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  left  = &x1[numStates * i];
	  right = &x2[numStates * i];
	  sum = &sumtable[numStates * i];

	  for(l = 0; l < numStates; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    default:
      assert(0);
    }
}


static void sumGTRCATSECONDARY_6(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				 unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  left  = &(tipVector[6 * tipX1[i]]);
	  right = &(tipVector[6 * tipX2[i]]);
	  sum = &sumtable[6 * i];

	  for(l = 0; l < 6; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  left = &(tipVector[6 * tipX1[i]]);
	  right = &x2[6 * i];
	  sum = &sumtable[6 * i];

	  for(l = 0; l < 6; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  left  = &x1[6 * i];
	  right = &x2[6 * i];
	  sum = &sumtable[6 * i];

	  for(l = 0; l < 6; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    default:
      assert(0);
    }
}

static void sumGTRCATSECONDARY_7(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				 unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = 0; i < n; i++)
	{
	  left  = &(tipVector[7 * tipX1[i]]);
	  right = &(tipVector[7 * tipX2[i]]);
	  sum = &sumtable[7 * i];

	  for(l = 0; l < 7; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  left = &(tipVector[7 * tipX1[i]]);
	  right = &x2[7 * i];
	  sum = &sumtable[7 * i];

	  for(l = 0; l < 7; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  left  = &x1[7 * i];
	  right = &x2[7 * i];
	  sum = &sumtable[7 * i];

	  for(l = 0; l < 7; l++)
	    sum[l] = left[l] * right[l];
	}
      break;
    default:
      assert(0);
    }
}


#ifdef __SIM_SSE3

static void coreGTRCATPROT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
			   volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable, int *wgt)
{
  int i, l;
  double *d1, *d_start, *sum;
  double 
    e[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
    s[20] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
    dd[20] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)rax_malloc(numberOfCategories * 20 * sizeof(double));

  e[0] = 0.0;
  s[0] = 0.0; 

  for(l = 1; l < 20; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz;
    }

  for(i = 0; i < numberOfCategories; i++)
    {      
      d1[20 * i] = 1.0;
      for(l = 1; l < 20; l++)
	d1[20 * i + l] = EXP(dd[l] * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];

      __m128d a0 = _mm_setzero_pd();
      __m128d a1 = _mm_setzero_pd();
      __m128d a2 = _mm_setzero_pd();

      d1 = &d_start[20 * cptr[i]];
      sum = &sumtable[20 * i];
          
      for(l = 0; l < 20; l+=2)
	{	  
	  __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d1[l]), _mm_load_pd(&sum[l]));
	  
	  a0 = _mm_add_pd(a0, tmpv);
	  __m128d sv = _mm_load_pd(&s[l]);	  
	  
	  a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, sv));
	  __m128d ev = _mm_load_pd(&e[l]);	  

	  a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, ev));
	}

      a0 = _mm_hadd_pd(a0, a0);
      a1 = _mm_hadd_pd(a1, a1);
      a2 = _mm_hadd_pd(a2, a2);

      _mm_storel_pd(&inv_Li, a0);     
      _mm_storel_pd(&dlnLidlz, a1);                 
      _mm_storel_pd(&d2lnLidlz2, a2);
      
      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}

#else

static void coreGTRCATPROT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
			   volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable, int *wgt)
{
  int i, l;
  double *d1, *d_start, *sum;
  double e[20], s[20], dd[20];
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)rax_malloc(numberOfCategories * 20 * sizeof(double));

  for(l = 1; l < 20; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz;
    }

  for(i = 0; i < numberOfCategories; i++)
    {
      for(l = 1; l < 20; l++)
	d1[20 * i + l] = EXP(dd[l] * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];
      
      d1 = &d_start[20 * cptr[i]];
      sum = &sumtable[20 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;

      for(l = 1; l < 20; l++)
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}

#endif



static void coreCatFlex(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
			volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable,
			const int numStates, int *wgt)
{
  int i, l;
  double 
    *d1, 
    *d_start, 
    *sum,
    e[64], 
    s[64], 
    dd[64],
    tmp,
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)rax_malloc(numberOfCategories * numStates * sizeof(double));

  for(l = 1; l < numStates; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz;
    }

  for(i = 0; i < numberOfCategories; i++)
    {
      for(l = 1; l < numStates; l++)
	d1[numStates * i + l] = EXP(dd[l] * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];

      d1 = &d_start[numStates * cptr[i]];
      sum = &sumtable[numStates * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;

      for(l = 1; l < numStates; l++)
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}

static void coreGammaFlex(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
			  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, const int numStates)
{
  double  
    *sum, 
    diagptable[1024],
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr,
    tmp,
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2;

  int     
    i, 
    j, 
    l;  

  const int 
    gammaStates = 4 * numStates;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < numStates; l++)
	{
	  diagptable[i * gammaStates + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * gammaStates + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * gammaStates + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      sum = &sumtable[i * gammaStates];
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * numStates];

	  for(l = 1; l < numStates; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * gammaStates + l * 4] * sum[j * numStates + l]);
	      dlnLidlz   +=  tmp * diagptable[j * gammaStates + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * gammaStates + l * 4 + 2];
	    }
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}

static double pow2(double x)
{
  return (x * x);
}

static double pow3(double x)
{
  return (x * x *x);
}

static double coreCatAsc(double *EIGN, double *sumtable, int upper,
			 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, const int numStates, 
			 volatile double *standard_ext_dlnLdlz,  volatile double *standard_ext_d2lnLdlz2, 
			 volatile double *felsenstein_ext_dlnLdlz,  volatile double *felsenstein_ext_d2lnLdlz2,
			 volatile double *goldman1_d1, volatile double *goldman1_d2,
			 volatile double *goldman2_d1, volatile double *goldman2_d2,
			 volatile double *goldman3_d1, volatile double *goldman3_d2,
			 double *weightVector, double *ascScaler)
{
  double  
    diagptable[1024], 
    lh = 0.0,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    standard_dlnLdlz = 0.0,
    standard_d2lnLdlz2 = 0.0,
    extended_lh = 0.0,
    extended_dlnLdlz = 0.0,
    extended_d2lnLdlz2 = 0.0,
    loglh = 0.0,
    ki, 
    kisqr;

  int     
    i,     
    l;  

 
  ki = 1.0;
  kisqr = 1.0;

  for(l = 1; l < numStates; l++)
    {
      diagptable[l * 4]     = EXP(EIGN[l-1] * ki * lz);
      diagptable[l * 4 + 1] = EIGN[l-1] * ki;
      diagptable[l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      double
	*sum = &sumtable[i * numStates],
	tmp,
	inv_Li   = 0.0,
	dlnLidlz = 0.0,
	d2lnLidlz2 = 0.0;

    
      inv_Li += sum[0];

      for(l = 1; l < numStates; l++)
	{
	  inv_Li     += (tmp = diagptable[l * 4] * sum[l]);
	  dlnLidlz   += tmp * diagptable[l * 4 + 1];
	  d2lnLidlz2 += tmp * diagptable[l * 4 + 2];
	}	            
            
      if(weightVector)
	{
	  double
	    standard_inv_Li = inv_Li,
	    standard_dlnLidlz = dlnLidlz,
	    standard_d2lnLidlz2 = d2lnLidlz2;
	    
	  standard_inv_Li = 1.0 / FABS(standard_inv_Li);

	  standard_dlnLidlz   *= standard_inv_Li;
	  standard_d2lnLidlz2 *= standard_inv_Li;

	  standard_dlnLdlz   += weightVector[i] * standard_dlnLidlz;
	  standard_d2lnLdlz2 += weightVector[i] * (standard_d2lnLidlz2 - standard_dlnLidlz * standard_dlnLidlz);
	}
      
      inv_Li = FABS(inv_Li);             
       
      //TODO scaler ???
      //g(x) := derivatives of f(x) * log(f(x)) goldman corrections 2 and 3 where f(x) is the likelihood function 
      extended_lh        += inv_Li * log(inv_Li); // Li * log(Li)
      extended_dlnLdlz   += (dlnLidlz * (log(inv_Li) + 1.0)); //Li' * (log(Li) + 1)			   
      extended_d2lnLdlz2 += ((d2lnLidlz2 * (log(inv_Li) + 1.0)) + (pow2(dlnLidlz) / inv_Li)); // (Li'' * (log(Li) + 1)) + (Li')^2 / Li

      //sum over likelihood and 1st and 2nd derivative
      //ascScaler is the numerical scaling for preventing underflow 
      lh        += inv_Li * ascScaler[i];	  	 
      dlnLdlz   += dlnLidlz * ascScaler[i];
      d2lnLdlz2 += d2lnLidlz2 * ascScaler[i];            
    } 

  //calculate the log likelihood 
  loglh = log(lh);

  //printf("loglh %f\n", loglh);

  *ext_dlnLdlz   = (dlnLdlz / (lh - 1.0));
  *ext_d2lnLdlz2 = (((lh - 1.0) * (d2lnLdlz2) - (dlnLdlz * dlnLdlz)) / ((lh - 1.0) * (lh - 1.0)));  

  *felsenstein_ext_dlnLdlz    = dlnLdlz / lh;
  *felsenstein_ext_d2lnLdlz2 = ((d2lnLdlz2 * lh) - (dlnLdlz * dlnLdlz)) / (lh * lh);

  *standard_ext_dlnLdlz   = standard_dlnLdlz;
  *standard_ext_d2lnLdlz2 = standard_d2lnLdlz2;

   //derivatives of (f(x) * log(f(x))) / (1 - f(x)) 
  //goldman 1 correction
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *goldman1_d1 = 
    (dlnLdlz * (-lh + loglh + 1.0)) / 
    pow2(1.0 - lh);
  
  *goldman1_d2 = 
    ((pow2(dlnLdlz) * (pow2(lh) - (2.0 * lh * loglh) - 1.0)) - ((lh - 1.0) * lh * d2lnLdlz2 * (lh - loglh - 1.0))) / 
    (lh * pow3(lh - 1.0));

  //derivatives of g(x) / (1 - f(x)), where the derivatives of g(x) were computed above in extended_*, i.e.,
  // g() denotes the derivatives of the sum of per-site L * log(L) over all character states
  //goldman 2 correction
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *goldman2_d1 = 
    (extended_lh * dlnLdlz - ((lh - 1.0) * extended_dlnLdlz)) 
    / pow2(1.0 - lh);
  
  *goldman2_d2 = 
    ((2.0 * dlnLdlz * extended_dlnLdlz) / pow2(1.0 - lh)) 
    + (extended_lh * ((d2lnLdlz2 / pow2(1.0 - lh))) + (2.0 * pow2(dlnLdlz) / (pow3(1.0 - lh))))
    + (extended_d2lnLdlz2 / (1.0 - lh));
  
  //derivatives of g(x) / f(x), where the derivatives of g(x) were computed above in extended_*
  //goldman 3 correction 
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *goldman3_d1 = 
    ((lh * extended_dlnLdlz) - (extended_lh * dlnLdlz)) / 
    pow2(lh);

  *goldman3_d2 = 
    ((-lh * extended_lh * d2lnLdlz2) - (2.0 * lh * dlnLdlz * extended_dlnLdlz) + (2.0 * extended_lh * pow2(dlnLdlz)) + (pow2(lh) * extended_d2lnLdlz2)) /
    pow3(lh);



  return lh;
}





static double coreGammaAsc(double *gammaRates, double *EIGN, double *sumtable, int upper,
			   volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, const int numStates, 
			   volatile double *standard_ext_dlnLdlz,  volatile double *standard_ext_d2lnLdlz2, 
			   volatile double *felsenstein_ext_dlnLdlz,  volatile double *felsenstein_ext_d2lnLdlz2,
			   volatile double *goldman1_d1, volatile double *goldman1_d2,
			   volatile double *goldman2_d1, volatile double *goldman2_d2,
			   volatile double *goldman3_d1, volatile double *goldman3_d2,
			   double *weightVector, double *ascScaler)
{
  double  
    diagptable[1024], 
    lh = 0.0,
    loglh = 0.0,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    standard_dlnLdlz = 0.0,
    standard_d2lnLdlz2 = 0.0,  
    extended_lh = 0.0,
    extended_dlnLdlz = 0.0,
    extended_d2lnLdlz2 = 0.0,
    ki, 
    kisqr;

  int     
    i, 
    j, 
    l;  

  const int 
    gammaStates = 4 * numStates;

  //claculate derivatives of P(lz) = exp(EIGN * ki * lz)

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < numStates; l++)
	{
	  diagptable[i * gammaStates + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * gammaStates + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * gammaStates + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  //loop over sites 

  for (i = 0; i < upper; i++)
    {
      double
	*sum = &sumtable[i * gammaStates],
	tmp,
	inv_Li   = 0.0, //likelihood at site i
	dlnLidlz = 0.0, //1st derivative of likelihood at site i
	d2lnLidlz2 = 0.0; //2nd derivative of likelihood at site i 


      //calc derviavtives 
      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * numStates];

	  for(l = 1; l < numStates; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * gammaStates + l * 4] * sum[j * numStates + l]);
	      dlnLidlz   += tmp * diagptable[j * gammaStates + l * 4 + 1];
	      d2lnLidlz2 += tmp * diagptable[j * gammaStates + l * 4 + 2];
	    }	  
	}    
            
      //this is for the stamatakis correction, i.e. just the derivative of log(f(x)), where f(x) is the likelihood function
      if(weightVector)
	{
	  double
	    standard_inv_Li = inv_Li,
	    standard_dlnLidlz = dlnLidlz,
	    standard_d2lnLidlz2 = d2lnLidlz2;
	  
	  standard_inv_Li = 1.0 / FABS(standard_inv_Li);

	  standard_dlnLidlz   *= standard_inv_Li;
	  standard_d2lnLidlz2 *= standard_inv_Li;

	  standard_dlnLdlz   += weightVector[i] * standard_dlnLidlz;
	  standard_d2lnLdlz2 += weightVector[i] * (standard_d2lnLidlz2 - standard_dlnLidlz * standard_dlnLidlz);	  	  
	}

     
      //multiply with 1/4 because of GAMMA 
      inv_Li = 0.25 * FABS(inv_Li);  //site likelihood       
      dlnLidlz *= 0.25;             //1st derivative 
      d2lnLidlz2 *= 0.25;           //2nd derivative

      //g(x) := derivatives of f(x) * log(f(x)) goldman corrections 2 and 3 where f(x) is the likelihood function 
      extended_lh        += inv_Li * log(inv_Li); // Li * log(Li)
      extended_dlnLdlz   += (dlnLidlz * (log(inv_Li) + 1.0)); //Li' * (log(Li) + 1)			   
      extended_d2lnLdlz2 += ((d2lnLidlz2 * (log(inv_Li) + 1.0)) + (pow2(dlnLidlz) / inv_Li)); // (Li'' * (log(Li) + 1)) + (Li')^2 / Li

      
      //sum over likelihood and 1st and 2nd derivative
      //ascScaler is the numerical scaling for preventing underflow 
      lh        += inv_Li * ascScaler[i];	  	 
      dlnLdlz   += dlnLidlz * ascScaler[i];
      d2lnLdlz2 += d2lnLidlz2 * ascScaler[i];       

      //printf("SCALE: %f\n", ascScaler[i]);
    } 

  //calculate the log likelihood 
  loglh = log(lh);

  //printf("loglh %f\n", loglh);

  //1st and 2nd derivative for the Lewis correction i.e., derivatives of log(1.0 - f(x))
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *ext_dlnLdlz   = (dlnLdlz / (lh - 1.0));
  *ext_d2lnLdlz2 = ((lh - 1.0) * (d2lnLdlz2) - pow2(dlnLdlz)) / pow2(lh - 1.0);  
  
  //1st and 2nd derivative for the Felsenstein correction, i.e., derivative of log(f(x))
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *felsenstein_ext_dlnLdlz    = dlnLdlz / lh;
  *felsenstein_ext_d2lnLdlz2 = ((d2lnLdlz2 * lh) - pow2(dlnLdlz)) / pow2(lh);
  

  //derivatives of (f(x) * log(f(x))) / (1 - f(x)) 
  //goldman 1 correction
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *goldman1_d1 = 
    (dlnLdlz * (-lh + loglh + 1.0)) / 
    pow2(1.0 - lh);
  
  *goldman1_d2 = 
    ((pow2(dlnLdlz) * (pow2(lh) - (2.0 * lh * loglh) - 1.0)) - ((lh - 1.0) * lh * d2lnLdlz2 * (lh - loglh - 1.0))) / 
    (lh * pow3(lh - 1.0));

  //derivatives of g(x) / (1 - f(x)), where the derivatives of g(x) were computed above in extended_*, i.e.,
  // g() denotes the derivatives of the sum of per-site L * log(L) over all four character states
  //goldman 2 correction
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *goldman2_d1 = 
    (extended_lh * dlnLdlz - ((lh - 1.0) * extended_dlnLdlz)) 
    / pow2(1.0 - lh);
  
  *goldman2_d2 = 
    ((2.0 * dlnLdlz * extended_dlnLdlz) / pow2(1.0 - lh)) 
    + (extended_lh * ((d2lnLdlz2 / pow2(1.0 - lh))) + (2.0 * pow2(dlnLdlz) / (pow3(1.0 - lh))))
    + (extended_d2lnLdlz2 / (1.0 - lh));
  
  //derivatives of g(x) / f(x), where the derivatives of g(x) were computed above in extended_*
  //goldman 3 correction 
  //here f(x) := L(A) + L(C) + L(G) + L(T) 
  *goldman3_d1 = 
    ((lh * extended_dlnLdlz) - (extended_lh * dlnLdlz)) / 
    pow2(lh);

  *goldman3_d2 = 
    ((-lh * extended_lh * d2lnLdlz2) - (2.0 * lh * dlnLdlz * extended_dlnLdlz) + (2.0 * extended_lh * pow2(dlnLdlz)) + (pow2(lh) * extended_d2lnLdlz2)) /
    pow3(lh);


  //normal derivatives of the log likelihood, i.e., log(f(x))
  //here f(x) := L(site_i) 
  *standard_ext_dlnLdlz   = standard_dlnLdlz;
  *standard_ext_d2lnLdlz2 = standard_d2lnLdlz2;

  return lh;
}




static void coreGammaInvarFlex(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
			       volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
			       double propInvar, int *iptr, const int numStates)
{
  double  
    *sum, 
    diagptable[1024],
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr,
    freqs[64],
    scaler =  0.25 * (1.0 - propInvar),
    tmp,
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2;
  
  int     
    i, 
    l, 
    j;
 
  const int 
    gammaStates = 4 * numStates;
  
  for(i = 0; i < numStates; i++)
    freqs[i] = frequencies[i] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < numStates; l++)
	{
	  diagptable[i * gammaStates + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * gammaStates + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * gammaStates + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }


   for(i = 0; i < upper; i++)
     {
       sum = &sumtable[i * gammaStates];
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;

       for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * numStates];

	  for(l = 1; l < numStates; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * gammaStates + l * 4] * sum[j * numStates + l]);
	      dlnLidlz   +=  tmp * diagptable[j * gammaStates + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * gammaStates + l * 4 + 2];
	    }
	}

       inv_Li = FABS(inv_Li) * scaler;

       if(iptr[i] < numStates)
	 inv_Li += freqs[iptr[i]];

       inv_Li = 1.0 / inv_Li;

       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;

       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;

       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;
}


static void coreGTRCATSECONDARY(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
				volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable, int *wgt)
{
  int i, l;
  double *d1, *d_start, *sum;
  double e[16], s[16], dd[16];
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)rax_malloc(numberOfCategories * 16 * sizeof(double));

  for(l = 1; l < 16; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz;
    }

  for(i = 0; i < numberOfCategories; i++)
    {
      for(l = 1; l < 16; l++)
	d1[16 * i + l] = EXP(dd[l] * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];
      
      d1 = &d_start[16 * cptr[i]];
      sum = &sumtable[16 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;

      for(l = 1; l < 16; l++)
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}

static void coreGTRCATSECONDARY_6(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
				  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable, int *wgt)
{
  int i, l;
  double *d1, *d_start, *sum;
  double e[6], s[6], dd[6];
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)rax_malloc(numberOfCategories * 6 * sizeof(double));

  for(l = 1; l < 6; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz;
    }

  for(i = 0; i < numberOfCategories; i++)
    {
      for(l = 1; l < 6; l++)
	d1[6 * i + l] = EXP(dd[l] * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];
      
      d1 = &d_start[6 * cptr[i]];
      sum = &sumtable[6 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;

      for(l = 1; l < 6; l++)
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}

static void coreGTRCATSECONDARY_7(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int upper,
				  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *sumtable, int *wgt)
{
  int i, l;
  double *d1, *d_start, *sum;
  double e[7], s[7], dd[7];
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double  dlnLdlz = 0.0;
  double  d2lnLdlz2 = 0.0;

  d1 = d_start = (double *)rax_malloc(numberOfCategories * 7 * sizeof(double));

  for(l = 1; l < 7; l++)
    {
      e[l]  = EIGN[l-1] * EIGN[l-1];
      s[l]  = EIGN[l-1];
      dd[l] = s[l] * lz;
    }

  for(i = 0; i < numberOfCategories; i++)
    {
      for(l = 1; l < 7; l++)
	d1[7 * i + l] = EXP(dd[l] * rptr[i]);
    }

  for (i = 0; i < upper; i++)
    {
      double
	r = rptr[cptr[i]],
	wr1 = r * wgt[i],
	wr2 = r * r * wgt[i];

      d1 = &d_start[7 * cptr[i]];
      sum = &sumtable[7 * i];

      inv_Li     = sum[0];
      dlnLidlz   = 0.0;
      d2lnLidlz2 = 0.0;

      for(l = 1; l < 7; l++)
	{
	  inv_Li     += (tmp = d1[l] * sum[l]);
	  dlnLidlz   += tmp *  s[l];
	  d2lnLidlz2 += tmp *  e[l];
	}

      inv_Li = 1.0/FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wr1 * dlnLidlz;
      d2lnLdlz2 += wr2 * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  rax_free(d_start);
}



static void coreGTRGAMMAINVAR_BINARY(double propInvar, double *frequencies, double *gammaRates, double *EIGN,
				     double *sumtable,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2,
				     int *iptr, int *wrptr, int upper, double lz)
{
  double  *sum, diagptable[12];
  int     i, j;
    double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double freqs[4];
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp_1;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  freqs[0] = frequencies[0] * propInvar;
  freqs[1] = frequencies[1] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      diagptable[i * 3]     = EXP (EIGN[0] * ki * lz);
      diagptable[i * 3 + 1] = EIGN[0] * ki;
      diagptable[i * 3 + 2] = EIGN[0] * EIGN[0] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      sum = &(sumtable[i * 8]);

      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[2 * j];


	  tmp_1      =  diagptable[3 * j] * sum[2 * j + 1];
	  inv_Li     += tmp_1;
	  dlnLidlz   += tmp_1 * diagptable[3 * j + 1];
	  d2lnLidlz2 += tmp_1 * diagptable[3 * j + 2];
	 }

      inv_Li = FABS(inv_Li) * scaler;

      if(iptr[i] < 2)
	inv_Li += freqs[iptr[i]];

      inv_Li = 1.0 / inv_Li;

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLidlz *= scaler;
      d2lnLidlz2 *= scaler;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}



static void coreGTRGAMMAINVAR(double propInvar, double *frequencies, double *gammaRates, double *EIGN,
			      double *sumtable,  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2,
			      int *iptr, int *wrptr, int upper, double lz)
{
  double  *sum, diagptable[64];
  int     i, j, k;
    double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double freqs[4];
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp_1;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  freqs[0] = frequencies[0] * propInvar;
  freqs[1] = frequencies[1] * propInvar;
  freqs[2] = frequencies[2] * propInvar;
  freqs[3] = frequencies[3] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      diagptable[i * 16]     = EXP (EIGN[0] * ki * lz);
      diagptable[i * 16 + 1] = EXP (EIGN[1] * ki * lz);
      diagptable[i * 16 + 2] = EXP (EIGN[2] * ki * lz);

      diagptable[i * 16 + 3] = EIGN[0] * ki;
      diagptable[i * 16 + 4] = EIGN[0] * EIGN[0] * kisqr;

      diagptable[i * 16 + 5] = EIGN[1] * ki;
      diagptable[i * 16 + 6] = EIGN[1] * EIGN[1] * kisqr;

      diagptable[i * 16 + 7] = EIGN[2] * ki;
      diagptable[i * 16 + 8] = EIGN[2] * EIGN[2] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      sum = &(sumtable[i * 16]);

       inv_Li      = 0.0;
       dlnLidlz    = 0.0;
       d2lnLidlz2  = 0.0;

       for(j = 0; j < 4; j++)
	 {
	   inv_Li += sum[4 * j];

	   for(k = 0; k < 3; k++)
	     {
	       tmp_1      =  diagptable[16 * j + k] * sum[4 * j + k + 1];
	       inv_Li     += tmp_1;
	       dlnLidlz   += tmp_1 * diagptable[16 * j + k * 2 + 3];
	       d2lnLidlz2 += tmp_1 * diagptable[16 * j + k * 2 + 4];
	     }
	 }

       inv_Li = FABS(inv_Li) * scaler;

      if(iptr[i] < 4)
	inv_Li += freqs[iptr[i]];

      inv_Li = 1.0 / inv_Li;

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLidlz *= scaler;
      d2lnLidlz2 *= scaler;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}


/* C-OPT the two functions below are used to optimize the branch lengths of the tree.
   together, they account for approx 25-30% of overall run time sumGAMMA requires less
   overall time than coreGAMMA though. */

static void sumGAMMA_BINARY(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
			    unsigned char *tipX1, unsigned char *tipX2, int n)
{
  double *x1, *x2, *sum;
  int i, j;
#ifndef __SIM_SSE3
  int k;
#endif

  /* C-OPT once again switch over possible configurations at inner node */

  switch(tipCase)
    {
    case TIP_TIP:
      /* C-OPT main for loop overt alignment length */
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[2 * tipX1[i]]);
	  x2 = &(tipVector[2 * tipX2[i]]);
	  sum = &sumtable[i * 8];
#ifndef __SIM_SSE3	  
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      sum[j * 2 + k] = x1[k] * x2[k];
#else
	  for(j = 0; j < 4; j++)
	    _mm_store_pd( &sum[j*2], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[0] )));	 
#endif
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  x1  = &(tipVector[2 * tipX1[i]]);
	  x2  = &x2_start[8 * i];
	  sum = &sumtable[8 * i];

#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      sum[j * 2 + k] = x1[k] * x2[j * 2 + k];
#else
	  for(j = 0; j < 4; j++)
	    _mm_store_pd( &sum[j*2], _mm_mul_pd( _mm_load_pd( &x1[0] ), _mm_load_pd( &x2[j * 2] )));
#endif
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1  = &x1_start[8 * i];
	  x2  = &x2_start[8 * i];
	  sum = &sumtable[8 * i];
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 2; k++)
	      sum[j * 2 + k] = x1[j * 2 + k] * x2[j * 2 + k];
#else
	  for(j = 0; j < 4; j++)
	    _mm_store_pd( &sum[j*2], _mm_mul_pd( _mm_load_pd( &x1[j * 2] ), _mm_load_pd( &x2[j * 2] )));
#endif
	}
      break;
    default:
      assert(0);
    }
}



static void sumGAMMA_GAPPED_SAVE(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
				 unsigned char *tipX1, unsigned char *tipX2, int n, 
				 double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  double 
    *x1, 
    *x2, 
    *sum,
    *x1_ptr = x1_start,
    *x2_ptr = x2_start;
  
  int i, j, k; 

  switch(tipCase)
    {
    case TIP_TIP:     
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  sum = &sumtable[i * 16];
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[k] * x2[k];
#else
	  for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k+=2)
	      _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[k] )));
#endif
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  x1  = &(tipVector[4 * tipX1[i]]);
	  
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    {
	      x2  = x2_ptr;
	      x2_ptr += 16;
	    }
	  
	  sum = &sumtable[16 * i];
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[k] * x2[j * 4 + k];
#else
	  for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k+=2)
	      _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[j * 4 + k] )));
#endif
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  if(x1_gap[i / 32] & mask32[i % 32])
	    x1 = x1_gapColumn;
	  else
	    {
	      x1  = x1_ptr;
	      x1_ptr += 16;
	    }
	  
	  if(x2_gap[i / 32] & mask32[i % 32])
	    x2 = x2_gapColumn;
	  else
	    {
	      x2  = x2_ptr;
	      x2_ptr += 16;
	    }

	  sum = &sumtable[16 * i];
	  
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[j * 4 + k] * x2[j * 4 + k];
#else
	   for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k+=2)
	      _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[j * 4 + k] ), _mm_load_pd( &x2[j * 4 + k] )));
#endif
	}
      break;
    default:
      assert(0);
    }
}





static void sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
		     unsigned char *tipX1, unsigned char *tipX2, int n)
{
  double *x1, *x2, *sum;
  int i, j, k;

  /* C-OPT once again switch over possible configurations at inner node */

  switch(tipCase)
    {
    case TIP_TIP:
      /* C-OPT main for loop overt alignment length */
      for (i = 0; i < n; i++)
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  sum = &sumtable[i * 16];
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[k] * x2[k];
#else
	  for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k+=2)
	      _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[k] )));
#endif
	}
      break;
    case TIP_INNER:
      for (i = 0; i < n; i++)
	{
	  x1  = &(tipVector[4 * tipX1[i]]);
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[k] * x2[j * 4 + k];
#else
	  for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k+=2)
	      _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[k] ), _mm_load_pd( &x2[j * 4 + k] )));
#endif
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  x1  = &x1_start[16 * i];
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
#ifndef __SIM_SSE3
	  for(j = 0; j < 4; j++)
	    for(k = 0; k < 4; k++)
	      sum[j * 4 + k] = x1[j * 4 + k] * x2[j * 4 + k];
#else
	   for(j = 0; j < 4; j++)	    
	    for(k = 0; k < 4; k+=2)
	      _mm_store_pd( &sum[j*4 + k], _mm_mul_pd( _mm_load_pd( &x1[j * 4 + k] ), _mm_load_pd( &x2[j * 4 + k] )));
#endif
	}
      break;
    default:
      assert(0);
    }
}




#ifndef __SIM_SSE3
static void coreGTRGAMMA_BINARY(const int upper, double *sumtable,
				volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  int i, j;
  double
    *diagptable, *diagptable_start, *sum,
    tmp_1, inv_Li, dlnLidlz, d2lnLidlz2, ki, kisqr,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  diagptable = diagptable_start = (double *)rax_malloc(sizeof(double) * 12);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      diagptable[i * 3]     = EXP (EIGN[0] * ki * lz);
      diagptable[i * 3 + 1] = EIGN[0] * ki;
      diagptable[i * 3 + 2] = EIGN[0] * EIGN[0] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      diagptable = diagptable_start;
      sum = &(sumtable[i * 8]);

      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[2 * j];

	  tmp_1      =  diagptable[3 * j] * sum[2 * j + 1];
	  inv_Li     += tmp_1;
	  dlnLidlz   += tmp_1 * diagptable[3 * j + 1];
	  d2lnLidlz2 += tmp_1 * diagptable[3 * j + 2];
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;


      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(diagptable_start);
}
#else
static void coreGTRGAMMA_BINARY(const int upper, double *sumtable,
				volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  double 
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr,  
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2,  
    *sum, 
    diagptable0[8] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable1[8] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable2[8] __attribute__ ((aligned (BYTE_ALIGNMENT)));    
    
  int     
    i, 
    j;
  
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
      
      diagptable0[i * 2] = 1.0;
      diagptable1[i * 2] = 0.0;
      diagptable2[i * 2] = 0.0;
     
      diagptable0[i * 2 + 1] = EXP(EIGN[0] * ki * lz);
      diagptable1[i * 2 + 1] = EIGN[0] * ki;
      diagptable2[i * 2 + 1] = EIGN[0] * EIGN[0] * kisqr;    
    }

  for (i = 0; i < upper; i++)
    { 
      __m128d a0 = _mm_setzero_pd();
      __m128d a1 = _mm_setzero_pd();
      __m128d a2 = _mm_setzero_pd();

      sum = &sumtable[i * 8];         

      for(j = 0; j < 4; j++)
	{	 	  	
	  double 	   
	    *d0 = &diagptable0[j * 2],
	    *d1 = &diagptable1[j * 2],
	    *d2 = &diagptable2[j * 2];
  	 	 	 
	  __m128d tmpv = _mm_mul_pd(_mm_load_pd(d0), _mm_load_pd(&sum[j * 2]));
	  a0 = _mm_add_pd(a0, tmpv);
	  a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(d1)));
	  a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(d2)));
	    	 	  
	}

      a0 = _mm_hadd_pd(a0, a0);
      a1 = _mm_hadd_pd(a1, a1);
      a2 = _mm_hadd_pd(a2, a2);

      _mm_storel_pd(&inv_Li, a0);     
      _mm_storel_pd(&dlnLidlz, a1);
      _mm_storel_pd(&d2lnLidlz2, a2); 

      inv_Li = 1.0 / FABS(inv_Li);
     
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;     

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

 
  *d1   = dlnLdlz;
  *d2 = d2lnLdlz2; 
}


#endif

#ifndef __SIM_SSE3
static void coreGTRGAMMA(const int upper, double *sumtable,
			 volatile double *d1,   volatile double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  int i, j, k;
  double
    *diagptable, *diagptable_start, *sum,
    tmp_1, inv_Li, dlnLidlz, d2lnLidlz2, ki, kisqr,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;

  diagptable = diagptable_start = (double *)rax_malloc(sizeof(double) * 64);



  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      diagptable[i * 16]     = EXP (EIGN[0] * ki * lz);
      diagptable[i * 16 + 1] = EXP (EIGN[1] * ki * lz);
      diagptable[i * 16 + 2] = EXP (EIGN[2] * ki * lz);

      diagptable[i * 16 + 3] = EIGN[0] * ki;
      diagptable[i * 16 + 4] = EIGN[0] * EIGN[0] * kisqr;

      diagptable[i * 16 + 5] = EIGN[1] * ki;
      diagptable[i * 16 + 6] = EIGN[1] * EIGN[1] * kisqr;

      diagptable[i * 16 + 7] = EIGN[2] * ki;
      diagptable[i * 16 + 8] = EIGN[2] * EIGN[2] * kisqr;
    }

  for (i = 0; i < upper; i++)
    {
      sum = &(sumtable[i * 16]);

      inv_Li      = 0.0;
      dlnLidlz    = 0.0;
      d2lnLidlz2  = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[4 * j];

	  for(k = 0; k < 3; k++)
	    {
	      tmp_1      =  diagptable[16 * j + k] * sum[4 * j + k + 1];
	      inv_Li     += tmp_1;
	      dlnLidlz   += tmp_1 * diagptable[16 * j + k * 2 + 3];
	      d2lnLidlz2 += tmp_1 * diagptable[16 * j + k * 2 + 4];
	    }
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;



      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  rax_free(diagptable_start);
}
#else

static void coreGTRGAMMA(const int upper, double *sumtable,
			 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  double 
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0,
    ki, 
    kisqr,  
    inv_Li, 
    dlnLidlz, 
    d2lnLidlz2,  
    *sum, 
    diagptable0[16] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable1[16] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable2[16] __attribute__ ((aligned (BYTE_ALIGNMENT)));    
    
  int     
    i, 
    j, 
    l;
  
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
      
      diagptable0[i * 4] = 1.0;
      diagptable1[i * 4] = 0.0;
      diagptable2[i * 4] = 0.0;

      for(l = 1; l < 4; l++)
	{
	  diagptable0[i * 4 + l] = EXP(EIGN[l-1] * ki * lz);
	  diagptable1[i * 4 + l] = EIGN[l-1] * ki;
	  diagptable2[i * 4 + l] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    { 
      __m128d a0 = _mm_setzero_pd();
      __m128d a1 = _mm_setzero_pd();
      __m128d a2 = _mm_setzero_pd();

      sum = &sumtable[i * 16];         

      for(j = 0; j < 4; j++)
	{	 	  	
	  double 	   
	    *d0 = &diagptable0[j * 4],
	    *d1 = &diagptable1[j * 4],
	    *d2 = &diagptable2[j * 4];
  	 	 
	  for(l = 0; l < 4; l+=2)
	    {
	      __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d0[l]), _mm_load_pd(&sum[j * 4 + l]));
	      a0 = _mm_add_pd(a0, tmpv);
	      a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(&d1[l])));
	      a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(&d2[l])));
	    }	 	  
	}

      a0 = _mm_hadd_pd(a0, a0);
      a1 = _mm_hadd_pd(a1, a1);
      a2 = _mm_hadd_pd(a2, a2);

      _mm_storel_pd(&inv_Li, a0);     
      _mm_storel_pd(&dlnLidlz, a1);
      _mm_storel_pd(&d2lnLidlz2, a2);       

      inv_Li = 1.0 / FABS(inv_Li);
     
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;     

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }
 
  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2; 
}


#endif









static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			 unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[20 * tipX1[i]]);
	  right = &(tipVector[20 * tipX2[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 80 + l * 20];
#ifdef __SIM_SSE3
	      for(k = 0; k < 20; k+=2)
		{
		  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
		  
		  _mm_store_pd(&sum[k], sumv);		 
		}
#else
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];
#endif
	    }
	}
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[20 * tipX1[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[80 * i + l * 20]);
	      sum = &sumtable[i * 80 + l * 20];
#ifdef __SIM_SSE3
	      for(k = 0; k < 20; k+=2)
		{
		  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
		  
		  _mm_store_pd(&sum[k], sumv);		 
		}
#else
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];
#endif
	    }
	}
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[80 * i + l * 20]);
	      right = &(x2[80 * i + l * 20]);
	      sum   = &(sumtable[i * 80 + l * 20]);

#ifdef __SIM_SSE3
	      for(k = 0; k < 20; k+=2)
		{
		  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
		  
		  _mm_store_pd(&sum[k], sumv);		 
		}
#else
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];
#endif
	    }
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMAPROT_LG4(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector[4],
			     unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:
      for(i = 0; i < n; i++)
	{	  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(tipVector[l][20 * tipX1[i]]);
	      right = &(tipVector[l][20 * tipX2[i]]);

	      sum = &sumtable[i * 80 + l * 20];
#ifdef __SIM_SSE3
	      for(k = 0; k < 20; k+=2)
		{
		  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
		  
		  _mm_store_pd(&sum[k], sumv);		 
		}
#else
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];
#endif
	    }
	}
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
	{
	 

	  for(l = 0; l < 4; l++)
	    { 
	      left = &(tipVector[l][20 * tipX1[i]]);
	      right = &(x2[80 * i + l * 20]);
	      sum = &sumtable[i * 80 + l * 20];
#ifdef __SIM_SSE3
	      for(k = 0; k < 20; k+=2)
		{
		  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
		  
		  _mm_store_pd(&sum[k], sumv);		 
		}
#else
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];
#endif
	    }
	}
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[80 * i + l * 20]);
	      right = &(x2[80 * i + l * 20]);
	      sum   = &(sumtable[i * 80 + l * 20]);

#ifdef __SIM_SSE3
	      for(k = 0; k < 20; k+=2)
		{
		  __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));
		  
		  _mm_store_pd(&sum[k], sumv);		 
		}
#else
	      for(k = 0; k < 20; k++)
		sum[k] = left[k] * right[k];
#endif
	    }
	}
      break;
    default:
      assert(0);
    }
}


#ifdef __SIM_SSE3
static void sumGAMMAPROT_GAPPED_SAVE(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
    unsigned char *tipX1, unsigned char *tipX2, int n, 
    double *x1_gapColumn, double *x2_gapColumn, unsigned int *x1_gap, unsigned int *x2_gap)
{
  int i, l, k;
  double 
    *left, 
    *right, 
    *sum,
    *x1_ptr = x1,
    *x2_ptr = x2,
    *x1v,
    *x2v;

  switch(tipCase)
  {
    case TIP_TIP:
      for(i = 0; i < n; i++)
      {
        left  = &(tipVector[20 * tipX1[i]]);
        right = &(tipVector[20 * tipX2[i]]);

        for(l = 0; l < 4; l++)
        {
          sum = &sumtable[i * 80 + l * 20];

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);		 
          }

        }
      }
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
      {
        left = &(tipVector[20 * tipX1[i]]);

        if(x2_gap[i / 32] & mask32[i % 32])
          x2v = x2_gapColumn;
        else
        {
          x2v = x2_ptr;
          x2_ptr += 80;
        }

        for(l = 0; l < 4; l++)
        {
          right = &(x2v[l * 20]);
          sum = &sumtable[i * 80 + l * 20];

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);		 
          }
        }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
      {
        if(x1_gap[i / 32] & mask32[i % 32])
          x1v = x1_gapColumn;
        else
        {
          x1v  = x1_ptr;
          x1_ptr += 80;
        }

        if(x2_gap[i / 32] & mask32[i % 32])
          x2v = x2_gapColumn;
        else
        {
          x2v  = x2_ptr;
          x2_ptr += 80;
        }

        for(l = 0; l < 4; l++)
        {
          left  = &(x1v[l * 20]);
          right = &(x2v[l * 20]);
          sum   = &(sumtable[i * 80 + l * 20]);

          for(k = 0; k < 20; k+=2)
          {
            __m128d sumv = _mm_mul_pd(_mm_load_pd(&left[k]), _mm_load_pd(&right[k]));

            _mm_store_pd(&sum[k], sumv);		 
          }
        }
      }
      break;
    default:
      assert(0);
  }
}

#endif

static void sumCatAsc(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
		      int n, const int numStates, int dataType, int qNumber, int rNumber, int maxtips)
{
  int i, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	unsigned char 
	  tip1[32],
	  tip2[32];

	ascertainmentBiasSequence(tip1, numStates, dataType, qNumber);
	ascertainmentBiasSequence(tip2, numStates, dataType, rNumber);
	
	for(i = 0; i < n; i++)
	  {
	    left  = &(tipVector[numStates * tip1[i]]);
	    right = &(tipVector[numStates * tip2[i]]);
	    	    
	    sum = &sumtable[i * numStates];
	    
	    for(k = 0; k < numStates; k++)
	      sum[k] = left[k] * right[k];	  
	  }
      }
      break;
    case TIP_INNER:
      {
	unsigned char 
	  tip[32];

	if(rNumber <= maxtips)
	  ascertainmentBiasSequence(tip, numStates, dataType, rNumber);
	else
	  ascertainmentBiasSequence(tip, numStates, dataType, qNumber);

	for(i = 0; i < n; i++)
	  {
	    left = &(tipVector[numStates * tip[i]]);
	    
	    
	    right = &(x2[i * numStates]);
	    sum = &sumtable[i * numStates];
	    
	    for(k = 0; k < numStates; k++)
	      sum[k] = left[k] * right[k];	 
	  }
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  left  = &(x1[i * numStates]);
	  right = &(x2[i * numStates]);
	  sum   = &(sumtable[i * numStates]);

	  for(k = 0; k < numStates; k++)
	    sum[k] = left[k] * right[k];	 
	}
      break;
    default:
      assert(0);
    }
}

static void sumGammaAsc(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			int n, const int numStates, int dataType, int qNumber, int rNumber, int maxtips)
{
  int i, l, k;
  double *left, *right, *sum;

  const int 
    gammaStates = numStates * 4;

  switch(tipCase)
    {
    case TIP_TIP:
      {
	unsigned char 
	  tip1[32],
	  tip2[32];

	ascertainmentBiasSequence(tip1, numStates, dataType, qNumber);
	ascertainmentBiasSequence(tip2, numStates, dataType, rNumber);

	for(i = 0; i < n; i++)
	  {
	    left  = &(tipVector[numStates * tip1[i]]);
	    right = &(tipVector[numStates * tip2[i]]);
	    
	    for(l = 0; l < 4; l++)
	      {
		sum = &sumtable[i * gammaStates + l * numStates];
		for(k = 0; k < numStates; k++)
		  sum[k] = left[k] * right[k];
	      }
	  }
      }
      break;
    case TIP_INNER:
      {
	unsigned char 
	  tip[32];

	if(rNumber < maxtips)
	  ascertainmentBiasSequence(tip, numStates, dataType, rNumber);
	else
	  ascertainmentBiasSequence(tip, numStates, dataType, qNumber);

      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[numStates * tip[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[gammaStates * i + l * numStates]);
	      sum = &sumtable[i * gammaStates + l * numStates];

	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      }
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[gammaStates * i + l * numStates]);
	      right = &(x2[gammaStates * i + l * numStates]);
	      sum   = &(sumtable[i * gammaStates + l * numStates]);

	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    default:
      assert(0);
    }
}


static void sumGammaFlex(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			 unsigned char *tipX1, unsigned char *tipX2, int n, const int numStates)
{
  int i, l, k;
  double *left, *right, *sum;

  const int gammaStates = numStates * 4;

  switch(tipCase)
    {
    case TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[numStates * tipX1[i]]);
	  right = &(tipVector[numStates * tipX2[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * gammaStates + l * numStates];
	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[numStates * tipX1[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[gammaStates * i + l * numStates]);
	      sum = &sumtable[i * gammaStates + l * numStates];

	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[gammaStates * i + l * numStates]);
	      right = &(x2[gammaStates * i + l * numStates]);
	      sum   = &(sumtable[i * gammaStates + l * numStates]);

	      for(k = 0; k < numStates; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMASECONDARY(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			      unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[16 * tipX1[i]]);
	  right = &(tipVector[16 * tipX2[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 64 + l * 16];
	      for(k = 0; k < 16; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[16 * tipX1[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[64 * i + l * 16]);
	      sum = &sumtable[i * 64 + l * 16];

	      for(k = 0; k < 16; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[64 * i + l * 16]);
	      right = &(x2[64 * i + l * 16]);
	      sum   = &(sumtable[i * 64 + l * 16]);

	      for(k = 0; k < 16; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMASECONDARY_6(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[6 * tipX1[i]]);
	  right = &(tipVector[6 * tipX2[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 24 + l * 6];
	      for(k = 0; k < 6; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[6 * tipX1[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[24 * i + l * 6]);
	      sum = &sumtable[i * 24 + l * 6];

	      for(k = 0; k < 6; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[24 * i + l * 6]);
	      right = &(x2[24 * i + l * 6]);
	      sum   = &(sumtable[i * 24 + l * 6]);

	      for(k = 0; k < 6; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    default:
      assert(0);
    }
}

static void sumGAMMASECONDARY_7(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
				unsigned char *tipX1, unsigned char *tipX2, int n)
{
  int i, l, k;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:
      for(i = 0; i < n; i++)
	{
	  left  = &(tipVector[7 * tipX1[i]]);
	  right = &(tipVector[7 * tipX2[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      sum = &sumtable[i * 28 + l * 7];
	      for(k = 0; k < 7; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case TIP_INNER:
      for(i = 0; i < n; i++)
	{
	  left = &(tipVector[7 * tipX1[i]]);

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[28 * i + l * 7]);
	      sum = &sumtable[i * 28 + l * 7];

	      for(k = 0; k < 7; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    case INNER_INNER:
      for(i = 0; i < n; i++)
	{
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[28 * i + l * 7]);
	      right = &(x2[28 * i + l * 7]);
	      sum   = &(sumtable[i * 28 + l * 7]);

	      for(k = 0; k < 7; k++)
		sum[k] = left[k] * right[k];
	    }
	}
      break;
    default:
      assert(0);
    }
}

#ifdef __SIM_SSE3

static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
			      volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, 
    diagptable0[80] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable1[80] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable2[80] __attribute__ ((aligned (BYTE_ALIGNMENT)));    
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
      
      diagptable0[i * 20] = 1.0;
      diagptable1[i * 20] = 0.0;
      diagptable2[i * 20] = 0.0;

      for(l = 1; l < 20; l++)
	{
	  diagptable0[i * 20 + l] = EXP(EIGN[l-1] * ki * lz);
	  diagptable1[i * 20 + l] = EIGN[l-1] * ki;
	  diagptable2[i * 20 + l] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    { 
      __m128d a0 = _mm_setzero_pd();
      __m128d a1 = _mm_setzero_pd();
      __m128d a2 = _mm_setzero_pd();

      sum = &sumtable[i * 80];         

      for(j = 0; j < 4; j++)
	{	 	  	
	  double 	   
	    *d0 = &diagptable0[j * 20],
	    *d1 = &diagptable1[j * 20],
	    *d2 = &diagptable2[j * 20];
  	 	 
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d tmpv = _mm_mul_pd(_mm_load_pd(&d0[l]), _mm_load_pd(&sum[j * 20 +l]));
	      a0 = _mm_add_pd(a0, tmpv);
	      a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(&d1[l])));
	      a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(&d2[l])));
	    }	 	  
	}

      a0 = _mm_hadd_pd(a0, a0);
      a1 = _mm_hadd_pd(a1, a1);
      a2 = _mm_hadd_pd(a2, a2);

      _mm_storel_pd(&inv_Li, a0);
      _mm_storel_pd(&dlnLidlz, a1);
      _mm_storel_pd(&d2lnLidlz2, a2);

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMAPROT_LG4(double *gammaRates, double *EIGN[4], double *sumtable, int upper, int *wrptr,
				 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *weights)
{
  double  *sum, 
    diagptable0[80] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable1[80] __attribute__ ((aligned (BYTE_ALIGNMENT))),
    diagptable2[80] __attribute__ ((aligned (BYTE_ALIGNMENT)));    
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;   

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
      
      diagptable0[i * 20] = 1.0;
      diagptable1[i * 20] = 0.0;
      diagptable2[i * 20] = 0.0;

      for(l = 1; l < 20; l++)
	{
	  diagptable0[i * 20 + l] = EXP(EIGN[i][l-1] * ki * lz);
	  diagptable1[i * 20 + l] = EIGN[i][l-1] * ki;
	  diagptable2[i * 20 + l] = EIGN[i][l-1] * EIGN[i][l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {       
      double 	
	inv_Li = 0.0, 
	dlnLidlz = 0.0, 
	d2lnLidlz2 = 0.0;

      sum = &sumtable[i * 80];         

      for(j = 0; j < 4; j++)
	{	 	  	
	  double
	    l0,
	    l1,
	    l2, 	   
	    *d0 = &diagptable0[j * 20],
	    *d1 = &diagptable1[j * 20],
	    *d2 = &diagptable2[j * 20];
	  
	  __m128d 
	    a0 = _mm_setzero_pd(),
	    a1 = _mm_setzero_pd(),
	    a2 = _mm_setzero_pd();
  	 	 
	  for(l = 0; l < 20; l+=2)
	    {
	      __m128d 
		tmpv = _mm_mul_pd(_mm_load_pd(&d0[l]), _mm_load_pd(&sum[j * 20 +l]));
	      
	      a0 = _mm_add_pd(a0, tmpv);
	      a1 = _mm_add_pd(a1, _mm_mul_pd(tmpv, _mm_load_pd(&d1[l])));
	      a2 = _mm_add_pd(a2, _mm_mul_pd(tmpv, _mm_load_pd(&d2[l])));
	    }
	  
	  a0 = _mm_hadd_pd(a0, a0);
	  a1 = _mm_hadd_pd(a1, a1);
	  a2 = _mm_hadd_pd(a2, a2);

	  _mm_storel_pd(&l0, a0);
	  _mm_storel_pd(&l1, a1);
	  _mm_storel_pd(&l2, a2);
	  
	  inv_Li     += weights[j] * l0;
	  dlnLidlz   += weights[j] * l1;
	  d2lnLidlz2 += weights[j] * l2; 
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}


#else

static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
			      volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[512];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 20; l++)
	{
	  diagptable[i * 128 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 128 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 128 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      sum = &sumtable[i * 80];
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 20];

	  for(l = 1; l < 20; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 128 + l * 4] * sum[j * 20 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 128 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 128 + l * 4 + 2];
	    }
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMAPROT_LG4(double *gammaRates, double *EIGN[4], double *sumtable, int upper, int *wrptr,
				 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *weights)
{
  double  *sum, diagptable[512];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 20; l++)
	{
	  diagptable[i * 128 + l * 4]     = EXP(EIGN[i][l-1] * ki * lz);
	  diagptable[i * 128 + l * 4 + 1] = EIGN[i][l-1] * ki;
	  diagptable[i * 128 + l * 4 + 2] = EIGN[i][l-1] * EIGN[i][l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      sum = &sumtable[i * 80];
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  double 
	    a1 = 0.0,
	    a2 = 0.0,
	    a3 = 0.0;

	  a1 += sum[j * 20];

	  for(l = 1; l < 20; l++)
	    {
	      a1 += (tmp = diagptable[j * 128 + l * 4] * sum[j * 20 + l]);
	      a2 +=  tmp * diagptable[j * 128 + l * 4 + 1];
	      a3 +=  tmp * diagptable[j * 128 + l * 4 + 2];
	    }

	  inv_Li     += weights[j] * a1;
	  dlnLidlz   += weights[j] * a2;
	  d2lnLidlz2 += weights[j] * a3;
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}


#endif




static void coreGTRGAMMASECONDARY(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				  volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[256];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 16; l++)
	{
	  diagptable[i * 64 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 64 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 64 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      sum = &sumtable[i * 64];
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 16];

	  for(l = 1; l < 16; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 64 + l * 4] * sum[j * 16 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 64 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 64 + l * 4 + 2];
	    }
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMASECONDARY_6(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[96];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 6; l++)
	{
	  diagptable[i * 24 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 24 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 24 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      sum = &sumtable[i * 24];
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 6];

	  for(l = 1; l < 6; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 24 + l * 4] * sum[j * 6 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 24 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 24 + l * 4 + 2];
	    }
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMASECONDARY_7(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				    volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz)
{
  double  *sum, diagptable[112];
  int     i, j, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 7; l++)
	{
	  diagptable[i * 28 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 28 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 28 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }

  for (i = 0; i < upper; i++)
    {
      sum = &sumtable[i * 28];
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;

      for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 7];

	  for(l = 1; l < 7; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 28 + l * 4] * sum[j * 7 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 28 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 28 + l * 4 + 2];
	    }
	}

      inv_Li = 1.0 / FABS(inv_Li);

      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMAPROTINVAR(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				   volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
				  double propInvar, int *iptr)
{
  double  *sum, diagptable[512];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double freqs[20];
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 20; i++)
    freqs[i] = frequencies[i] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 20; l++)
	{
	  diagptable[i * 128 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 128 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 128 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }


   for(i = 0; i < upper; i++)
     {
       sum = &sumtable[i * 80];
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;

       for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 20];

	  for(l = 1; l < 20; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 128 + l * 4] * sum[j * 20 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 128 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 128 + l * 4 + 2];
	    }
	}

       inv_Li = FABS(inv_Li) * scaler;

       if(iptr[i] < 20)
	 inv_Li += freqs[iptr[i]];

       inv_Li = 1.0 / inv_Li;

       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;

       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;

       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;
}


static void coreGTRGAMMASECONDARYINVAR(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
				       volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
				       double propInvar, int *iptr)
{
  double  *sum, diagptable[256];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double freqs[16];
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 16; i++)
    freqs[i] = frequencies[i] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 16; l++)
	{
	  diagptable[i * 64 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 64 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 64 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }


   for(i = 0; i < upper; i++)
     {
       sum = &sumtable[i * 64];
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;

       for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 16];

	  for(l = 1; l < 16; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 64 + l * 4] * sum[j * 16 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 64 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 64 + l * 4 + 2];
	    }
	}

       inv_Li = FABS(inv_Li) * scaler;

       if(iptr[i] < 16)
	 inv_Li += freqs[iptr[i]];

       inv_Li = 1.0 / inv_Li;

       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;

       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;

       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMASECONDARYINVAR_6(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
					 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
					 double propInvar, int *iptr)
{
  double  *sum, diagptable[96];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double freqs[6];
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 6; i++)
    freqs[i] = frequencies[i] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 6; l++)
	{
	  diagptable[i * 24 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 24 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 24 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }


   for(i = 0; i < upper; i++)
     {
       sum = &sumtable[i * 24];
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;

       for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 6];

	  for(l = 1; l < 6; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 24 + l * 4] * sum[j * 6 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 24 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 24 + l * 4 + 2];
	    }
	}

       inv_Li = FABS(inv_Li) * scaler;

       if(iptr[i] < 6)
	 inv_Li += freqs[iptr[i]];

       inv_Li = 1.0 / inv_Li;

       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;

       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;

       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;
}

static void coreGTRGAMMASECONDARYINVAR_7(double *gammaRates, double *EIGN, double *sumtable, int upper, int *wrptr,
					 volatile double *ext_dlnLdlz,  volatile double *ext_d2lnLdlz2, double lz, double *frequencies,
					 double propInvar, int *iptr)
{
  double  *sum, diagptable[112];
  int     i, l, j;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;
  double freqs[7];
  double scaler =  0.25 * (1.0 - propInvar);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 7; i++)
    freqs[i] = frequencies[i] * propInvar;

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;

      for(l = 1; l < 7; l++)
	{
	  diagptable[i * 28 + l * 4]     = EXP(EIGN[l-1] * ki * lz);
	  diagptable[i * 28 + l * 4 + 1] = EIGN[l-1] * ki;
	  diagptable[i * 28 + l * 4 + 2] = EIGN[l-1] * EIGN[l-1] * kisqr;
	}
    }


   for(i = 0; i < upper; i++)
     {
       sum = &sumtable[i * 28];
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;

       for(j = 0; j < 4; j++)
	{
	  inv_Li += sum[j * 7];

	  for(l = 1; l < 7; l++)
	    {
	      inv_Li     += (tmp = diagptable[j * 28 + l * 4] * sum[j * 7 + l]);
	      dlnLidlz   +=  tmp * diagptable[j * 28 + l * 4 + 1];
	      d2lnLidlz2 +=  tmp * diagptable[j * 28 + l * 4 + 2];
	    }
	}

       inv_Li = FABS(inv_Li) * scaler;

       if(iptr[i] < 7)
	 inv_Li += freqs[iptr[i]];

       inv_Li = 1.0 / inv_Li;

       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;

       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;

       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2;
}


static void getVects(tree *tr, unsigned char **tipX1, unsigned char **tipX2, double **x1_start, double **x2_start, int *tipCase, int model,
		     double **x1_gapColumn, double **x2_gapColumn, unsigned int **x1_gap, unsigned int **x2_gap,  double **x1_start_asc, double **x2_start_asc)
{
  int 
    rateHet,
    states = tr->partitionData[model].states,    
    pNumber = tr->td[0].ti[0].pNumber, 
    qNumber = tr->td[0].ti[0].qNumber;

  if(tr->rateHetModel == CAT)
    rateHet = 1;
  else
    rateHet = 4;
  
  *x1_start_asc = (double*)NULL;
  *x2_start_asc = (double*)NULL;
  
  *x1_start = (double*)NULL,
  *x2_start = (double*)NULL;
  
  *tipX1 = (unsigned char*)NULL,
  *tipX2 = (unsigned char*)NULL;

  if(isTip(pNumber, tr->mxtips) || isTip(qNumber, tr->mxtips))
    {
      if(!( isTip(pNumber, tr->mxtips) && isTip(qNumber, tr->mxtips)) )
	{
	  *tipCase = TIP_INNER;
	  if(isTip(qNumber, tr->mxtips))
	    {
	      *tipX1 = tr->partitionData[model].yVector[qNumber];
	      *x2_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
	      
#ifdef _USE_PTHREADS
	      if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
		if(tr->partitionData[model].ascBias)
#endif
		{
		  *x2_start_asc = &tr->partitionData[model].ascVector[(pNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];
		}

	      if(tr->saveMemory)
		{
		  *x2_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
		  *x2_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet];  
		}
	    }
	  else
	    {
	      *tipX1 = tr->partitionData[model].yVector[pNumber];
	      *x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];
	      
#ifdef _USE_PTHREADS
	      if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
		if(tr->partitionData[model].ascBias)
#endif	      
	     
		{
		  *x2_start_asc = &tr->partitionData[model].ascVector[(qNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];
		}

	      if(tr->saveMemory)
		{
		  *x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
		  *x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet];
		}
	    }
	}
      else
	{
	  *tipCase = TIP_TIP;	  
	  *tipX1 = tr->partitionData[model].yVector[pNumber];
	  *tipX2 = tr->partitionData[model].yVector[qNumber];
	}
    }
  else
    {
      *tipCase = INNER_INNER;

      *x1_start = tr->partitionData[model].xVector[pNumber - tr->mxtips - 1];
      *x2_start = tr->partitionData[model].xVector[qNumber - tr->mxtips - 1];


#ifdef _USE_PTHREADS
      if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
	if(tr->partitionData[model].ascBias)
#endif
	{
	  *x1_start_asc = &tr->partitionData[model].ascVector[(pNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];
	  *x2_start_asc = &tr->partitionData[model].ascVector[(qNumber - tr->mxtips - 1) * tr->partitionData[model].ascOffset];
	}           
      
      if(tr->saveMemory)
	{
	  *x1_gap = &(tr->partitionData[model].gapVector[pNumber * tr->partitionData[model].gapVectorLength]);
	  *x1_gapColumn   = &tr->partitionData[model].gapColumn[(pNumber - tr->mxtips - 1) * states * rateHet]; 
      
	  *x2_gap = &(tr->partitionData[model].gapVector[qNumber * tr->partitionData[model].gapVectorLength]);
	  *x2_gapColumn   = &tr->partitionData[model].gapColumn[(qNumber - tr->mxtips - 1) * states * rateHet]; 
	}
    }

}




void makenewzIterative(tree *tr)
{
  int 
    model, 
    tipCase;

  double
    *x1_start = (double*)NULL,
    *x2_start = (double*)NULL,
    *x1_start_asc = (double*)NULL,
    *x2_start_asc = (double*)NULL,
    *x1_gapColumn = (double*)NULL,
    *x2_gapColumn = (double*)NULL;

  unsigned char
    *tipX1,
    *tipX2; 
  
  unsigned int
    *x1_gap = (unsigned int*)NULL,
    *x2_gap = (unsigned int*)NULL;			      
 
  newviewIterative(tr);

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      /*printf("MNZIterative %d %d\n", model, tr->executeModel[model]);*/
      
      if(tr->executeModel[model])
	{
	  int 
	    
	    width = tr->partitionData[model].width,
	    states = tr->partitionData[model].states;

	 
	  getVects(tr, &tipX1, &tipX2, &x1_start, &x2_start, &tipCase, model, &x1_gapColumn, &x2_gapColumn, &x1_gap, &x2_gap, &x1_start_asc, &x2_start_asc);
	 

	  switch(tr->partitionData[model].dataType)
	    {
	    case BINARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  sumCAT_BINARY(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMA_BINARY(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				  width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:		 
#ifdef __SIM_SSE3
		  if(tr->saveMemory)
		    {
		      
		      sumCAT_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				  width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
		    }
		  else
#endif
		    sumCAT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
			   width);
		  break;
		case GAMMA:
		case GAMMA_I:		  
		      {
			if(tr->saveMemory)
			  sumGAMMA_GAPPED_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
					       width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
			else			  
#ifdef _HET
			  sumGAMMA(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector_TIP, tipX1, tipX2, width);
#else
			  sumGAMMA(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				   width);			  
#endif			  
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
#ifdef __SIM_SSE3
		  if(tr->saveMemory)	  
		    sumGTRCATPROT_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				       tipX1, tipX2, width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
		  else
#endif
		    sumGTRCATPROT(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				  tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:		  
#ifdef __SIM_SSE3
		  if(tr->saveMemory)
		    sumGAMMAPROT_GAPPED_SAVE(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
					     width, x1_gapColumn, x2_gapColumn, x1_gap, x2_gap);
		  else
#endif		  
		    if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)		      
		      {			
			sumGAMMAPROT_LG4(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector_LG4,
					 tipX1, tipX2, width);
		      }
		    else
		      sumGAMMAPROT(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				   tipX1, tipX2, width);		       
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				     tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				    tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY_6(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				       tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY_6(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				      tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY_7(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				       tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY_7(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				      tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case GENERIC_32:
	       switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumCatFlex(tipCase, tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
			     tipX1, tipX2, width, states);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGammaFlex(tipCase,  tr->partitionData[model].sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
			       tipX1, tipX2, width, states);
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
#ifdef _USE_PTHREADS
	  if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
	    if(tr->partitionData[model].ascBias)
#endif	 
	    {
	      int 
		pNumber = tr->td[0].ti[0].pNumber, 
		qNumber = tr->td[0].ti[0].qNumber,
		i,
		*ex1_asc = &tr->partitionData[model].ascExpVector[(pNumber - tr->mxtips - 1) * states],
		*ex2_asc = &tr->partitionData[model].ascExpVector[(qNumber - tr->mxtips - 1) * states];

	      switch(tipCase)
		{
		case TIP_TIP:
		   for(i = 0; i < states; i++)
		     tr->partitionData[model].ascScaler[i] = 1.0;		  
		  break;
		case TIP_INNER:
		  if(isTip(pNumber, tr->mxtips))
		    {
		      for(i = 0; i < states; i++)
			tr->partitionData[model].ascScaler[i] = pow(minlikelihood, (double)ex2_asc[i]);
		    }
		  else
		    {
		      for(i = 0; i < states; i++)
			tr->partitionData[model].ascScaler[i] = pow(minlikelihood, (double)ex1_asc[i]);
		    }
		  break;
		case INNER_INNER:
		  for(i = 0; i < states; i++)
		    tr->partitionData[model].ascScaler[i] = pow(minlikelihood, (double)(ex1_asc[i] + ex2_asc[i]));
		  break;
		default:
		  assert(0);
		}

	      /*for(i = 0; i < states; i++)
		printf("%d ", ex2_asc[i]);
	      printf("\n");
	      for(i = 0; i < states; i++)
		printf("%d ", ex1_asc[i]);
		printf("\n");*/
	      /*for(i = 0; i < states; i++)
		if(tr->partitionData[model].ascScaler[i] != 1.0)
		  printf("%f \n", tr->partitionData[model].ascScaler[i]);
	      //printf("\n");*/
	      

	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumCatAsc(tipCase, tr->partitionData[model].ascSumBuffer, x1_start_asc, x2_start_asc, tr->partitionData[model].tipVector,
			    states, states, tr->partitionData[model].dataType, pNumber, qNumber, tr->mxtips);
		  break;
		case GAMMA:
		  sumGammaAsc(tipCase, tr->partitionData[model].ascSumBuffer, x1_start_asc, x2_start_asc, tr->partitionData[model].tipVector,
			      states, states, tr->partitionData[model].dataType, pNumber, qNumber, tr->mxtips);
		  break;
		default:
		  assert(0);
		}
	    }
	}
    }
}




#ifdef _HET
void execCore(tree *tr, volatile double *_dlnLdlz, volatile double *_d2lnLdlz2, boolean isTipBranch)
#else
void execCore(tree *tr, volatile double *_dlnLdlz, volatile double *_d2lnLdlz2)
#endif
{
  int 
    model, 
    branchIndex;

  double 
    lz;  

  _dlnLdlz[0]   = 0.0;
  _d2lnLdlz2[0] = 0.0;

#ifdef _DEBUG_MULTI_EPA  
  printf("MNZC: ");
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
	    states = tr->partitionData[model].states;
	  
	  double 
	    *weights         = tr->partitionData[model].weights,
	    *sumBuffer       = (double*)NULL;
	 

	  volatile double
	    dlnLdlz   = 0.0,
	    d2lnLdlz2 = 0.0;
	  

	  sumBuffer = tr->partitionData[model].sumBuffer;


	  if(tr->multiBranch)
	    {
	      branchIndex = model;
	      lz = tr->coreLZ[model];
	      _dlnLdlz[model]   = 0.0;
	      _d2lnLdlz2[model] = 0.0;
	    }
	  else
	    {
	      branchIndex = 0;
	      lz = tr->coreLZ[0];
	    }

	  

	  switch(tr->partitionData[model].dataType)
	    {
	    case BINARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCAT_BINARY(width, tr->partitionData[model].numberOfCategories, sumBuffer,
				    &dlnLdlz, &d2lnLdlz2, 
				    tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN,  tr->partitionData[model].rateCategory, lz, tr->partitionData[model].wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMA_BINARY(width, sumBuffer,
				      &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
				      tr->partitionData[model].wgt);
		  break;
		case GAMMA_I:
		  coreGTRGAMMAINVAR_BINARY(tr->partitionData[model].propInvariant, tr->partitionData[model].frequencies,
					   tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					   sumBuffer, &dlnLdlz, &d2lnLdlz2,
					   tr->partitionData[model].invariant, tr->partitionData[model].wgt,
					   width, lz);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:	      
		{
		  switch(tr->rateHetModel)
		    {
		    case CAT:
		      coreGTRCAT(width, tr->partitionData[model].numberOfCategories, sumBuffer,
				 &dlnLdlz, &d2lnLdlz2,
				 tr->partitionData[model].perSiteRates, tr->partitionData[model].EIGN,  tr->partitionData[model].rateCategory, lz, tr->partitionData[model].wgt);
		      break;
		    case GAMMA:	
#ifdef _HET
		      if(isTipBranch)
			coreGTRGAMMA(width, sumBuffer,
				     &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN_TIP, tr->partitionData[model].gammaRates, lz,
				     tr->partitionData[model].wgt);
		      else
			coreGTRGAMMA(width, sumBuffer,
				     &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
				     tr->partitionData[model].wgt);

#else	  	     
		      coreGTRGAMMA(width, sumBuffer,
				   &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
				   tr->partitionData[model].wgt);
#endif
		      break;
		    case GAMMA_I:
		      coreGTRGAMMAINVAR(tr->partitionData[model].propInvariant, tr->partitionData[model].frequencies,
					tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					sumBuffer, &dlnLdlz, &d2lnLdlz2,
					tr->partitionData[model].invariant, tr->partitionData[model].wgt,
					width, lz);
		      break;
		    default:
		      assert(0);
		    }
		}
	      break;
	    case AA_DATA:
	     
		{
		  switch(tr->rateHetModel)
		    {
		    case CAT:
		      coreGTRCATPROT(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  tr->partitionData[model].perSiteRates,
				     tr->partitionData[model].rateCategory, width,				    
				     &dlnLdlz, &d2lnLdlz2,
				     sumBuffer, tr->partitionData[model].wgt);
		      break;
		    case GAMMA:	
		      if(tr->partitionData[model].protModels == LG4 || tr->partitionData[model].protModels == LG4X)
			{						  
			  coreGTRGAMMAPROT_LG4(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN_LG4,
					      sumBuffer, width, tr->partitionData[model].wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, weights);

			}
		      else
			coreGTRGAMMAPROT(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					 sumBuffer, width, tr->partitionData[model].wgt,
					 &dlnLdlz, &d2lnLdlz2, lz);
		      break;
		    case GAMMA_I:
		      coreGTRGAMMAPROTINVAR(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					    sumBuffer, width, tr->partitionData[model].wgt,
					    &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					    tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		      break;
		    default:
		      assert(0);
		    }
		}
	      break;
	    case SECONDARY_DATA:	      
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  tr->partitionData[model].perSiteRates,
				      tr->partitionData[model].rateCategory, width,				      
				      &dlnLdlz, &d2lnLdlz2,
				      sumBuffer, tr->partitionData[model].wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					sumBuffer, width, tr->partitionData[model].wgt,
					&dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:		 
		  coreGTRGAMMASECONDARYINVAR(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					     sumBuffer, width, tr->partitionData[model].wgt,
					     &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					     tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:	      
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY_6(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  tr->partitionData[model].perSiteRates,
					tr->partitionData[model].rateCategory, width,					
					&dlnLdlz, &d2lnLdlz2,
					sumBuffer, tr->partitionData[model].wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY_6(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					  sumBuffer, width, tr->partitionData[model].wgt,
					  &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR_6(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					       sumBuffer, width, tr->partitionData[model].wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					       tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:	    
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY_7(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  tr->partitionData[model].perSiteRates,
					tr->partitionData[model].rateCategory, width,				      
					&dlnLdlz, &d2lnLdlz2,
					sumBuffer, tr->partitionData[model].wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY_7(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					  sumBuffer, width, tr->partitionData[model].wgt,
					  &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR_7(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					       sumBuffer, width, tr->partitionData[model].wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					       tr->partitionData[model].propInvariant, tr->partitionData[model].invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case GENERIC_32:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreCatFlex(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  tr->partitionData[model].perSiteRates,
			      tr->partitionData[model].rateCategory, width,			      
			      &dlnLdlz, &d2lnLdlz2,
			      sumBuffer, states, tr->partitionData[model].wgt);
		  break;
		case GAMMA:
		  coreGammaFlex(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				sumBuffer, width, tr->partitionData[model].wgt,
				&dlnLdlz, &d2lnLdlz2, lz, states);
		  break;
		case GAMMA_I:
		  coreGammaInvarFlex(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				     sumBuffer, width, tr->partitionData[model].wgt,
				     &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
				     tr->partitionData[model].propInvariant, tr->partitionData[model].invariant, states);
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }

	 
	  
#ifdef _USE_PTHREADS
	  if(tr->partitionData[model].ascBias && tr->threadID == 0)
#else
	    if(tr->partitionData[model].ascBias)
#endif      
	    {
	      size_t
		i;

	      double 		
		*weightVector = (double*)NULL;

	      int	     
		w = 0;
	      
	      volatile double 
		d1 = 0.0,
		d2 = 0.0,
		standard_d1 = 0.0,
		standard_d2 = 0.0,
		felsenstein_d1 = 0.0,
		felsenstein_d2 = 0.0,
		goldman1_d1 = 0.0,
		goldman1_d2 = 0.0,
		goldman2_d1 = 0.0,
		goldman2_d2 = 0.0,
		goldman3_d1 = 0.0,
		goldman3_d2 = 0.0;
	      
	      if(tr->ascertainmentCorrectionType == STAMATAKIS_CORRECTION)		
		weightVector = tr->partitionData[model].invariableFrequencies;		    
	      
	       switch(tr->rateHetModel)
		 {
		 case CAT:
		   coreCatAsc(tr->partitionData[model].EIGN, tr->partitionData[model].ascSumBuffer, states,
			      &d1,  &d2, lz, states, 
			      &standard_d1, &standard_d2, 
			      &felsenstein_d1, &felsenstein_d2,
			      &goldman1_d1, &goldman1_d2,
			      &goldman2_d1, &goldman2_d2,
			      &goldman3_d1, &goldman3_d2,
			      weightVector, tr->partitionData[model].ascScaler);	  
		   break;
		 case GAMMA:
		   coreGammaAsc(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN, tr->partitionData[model].ascSumBuffer, states,
				&d1,  &d2, lz, states, 
				&standard_d1, &standard_d2, 
				&felsenstein_d1, &felsenstein_d2,
				&goldman1_d1, &goldman1_d2,
				&goldman2_d1, &goldman2_d2,
				&goldman3_d1, &goldman3_d2,
				weightVector, tr->partitionData[model].ascScaler);	  
		   break;
		 default:
		   assert(0);
		 }
	  
	       for(i = tr->partitionData[model].lower; i < tr->partitionData[model].upper; i++)
		 w += tr->cdta->aliaswgt[i];
	       
	       switch(tr->ascertainmentCorrectionType)
		 {
		 case LEWIS_CORRECTION:		   		   
		   _dlnLdlz[branchIndex]   +=  dlnLdlz - (double)w * d1;
		   _d2lnLdlz2[branchIndex] +=  d2lnLdlz2-  (double)w * d2;
		   break;
		 case STAMATAKIS_CORRECTION:
		   _dlnLdlz[branchIndex]   += dlnLdlz + standard_d1;
		   _d2lnLdlz2[branchIndex] += d2lnLdlz2 + standard_d2;
		   break;
		 case FELSENSTEIN_CORRECTION:
		   _dlnLdlz[branchIndex]   += dlnLdlz   + tr->partitionData[model].invariableWeight * felsenstein_d1;
		   _d2lnLdlz2[branchIndex] += d2lnLdlz2 + tr->partitionData[model].invariableWeight * felsenstein_d2;
		   break;
		 case GOLDMAN_CORRECTION_1:		  		 
		   _dlnLdlz[branchIndex]   += dlnLdlz   + (double)w * goldman1_d1;
		   _d2lnLdlz2[branchIndex] += d2lnLdlz2 + (double)w * goldman1_d2;
		   break;
		 case GOLDMAN_CORRECTION_2:
		   //printf("%f %f\n", goldman2_d1, goldman2_d2);
		   //exit(0);
		   _dlnLdlz[branchIndex]   += dlnLdlz   + (double)w * goldman2_d1;
		   _d2lnLdlz2[branchIndex] += d2lnLdlz2 + (double)w * goldman2_d2;
		   break;
		 case GOLDMAN_CORRECTION_3:
		   _dlnLdlz[branchIndex]   += dlnLdlz   + tr->partitionData[model].invariableWeight * goldman3_d1;
		   _d2lnLdlz2[branchIndex] += d2lnLdlz2 + tr->partitionData[model].invariableWeight * goldman3_d2;
		   break;		  
		 default:
		   assert(0);
		 }
	       
	    }  
	    else
	      {
		_dlnLdlz[branchIndex]   = _dlnLdlz[branchIndex]   + dlnLdlz;
		_d2lnLdlz2[branchIndex] = _d2lnLdlz2[branchIndex] + d2lnLdlz2;	      
	      }
	}
    }
  
#ifdef _DEBUG_MULTI_EPA
  printf("\n");
#endif 	  

#ifdef _DEBUG_ASC
  printf("\n\n");
#endif
}



#ifdef _HET
static void topLevelMakenewz(tree *tr, double *z0, int _maxiter, double *result, boolean isTipBranch)
#else
static void topLevelMakenewz(tree *tr, double *z0, int _maxiter, double *result)
#endif
{
  double   z[NUM_BRANCHES], zprev[NUM_BRANCHES], zstep[NUM_BRANCHES];
  volatile double  dlnLdlz[NUM_BRANCHES], d2lnLdlz2[NUM_BRANCHES];
  int i, maxiter[NUM_BRANCHES], numBranches, model;
  boolean 
    firstIteration = TRUE,
    outerConverged[NUM_BRANCHES],
    loopConverged;

  if(tr->multiBranch)
    numBranches = tr->NumberOfModels;
  else
    numBranches = 1;

  for(i = 0; i < numBranches; i++)
    {
      z[i] = z0[i];       

      maxiter[i] = _maxiter;
      
      outerConverged[i] = FALSE;
      tr->curvatOK[i]       = TRUE;      
    }

  if(tr->perPartitionEPA)   
    {
      if(tr->multiBranch)
	for(i = 0; i < numBranches; i++)
	  if(!(tr->executeModel[i]))
	    outerConverged[i] = TRUE;
    }

  do
    {
      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && tr->curvatOK[i] == TRUE)
	    {
	      tr->curvatOK[i] = FALSE;

	      zprev[i] = z[i];

	      zstep[i] = (1.0 - zmax) * z[i] + zmin;
	    }
	}

      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && tr->curvatOK[i] == FALSE)
	    {
	      double lz;

	      if (z[i] < zmin) 			  
		z[i] = zmin;		
	      else
		{
		  if(z[i] > zmax) 		    
		    z[i] = zmax;		      		    
		}

	      lz    = log(z[i]);

	      tr->coreLZ[i] = lz;
	    }
	}

      if(tr->multiBranch)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->executeModel[model])
		tr->executeModel[model] = !tr->curvatOK[model];

	    }
	}
      else
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->perPartitionEPA)
		{
		  if(tr->executeModel[model])
		    tr->executeModel[model] = !tr->curvatOK[0];
		}
	      else
		tr->executeModel[model] = !tr->curvatOK[0];
	    }
	}


#ifdef _USE_PTHREADS
      if(firstIteration)
	{
	  masterBarrier(THREAD_MAKENEWZ_FIRST, tr);
	  firstIteration = FALSE;
	}
      else
	masterBarrier(THREAD_MAKENEWZ, tr);

      if(!tr->multiBranch)
	{
	  dlnLdlz[0] = 0.0;
	  d2lnLdlz2[0] = 0.0;
	  for(i = 0; i < NumberOfThreads; i++)
	    {
	      dlnLdlz[0]   += reductionBuffer[i];
	      d2lnLdlz2[0] += reductionBufferTwo[i];
	    }
	}
      else
	{
	  int j;
	  for(j = 0; j < tr->NumberOfModels; j++)
	    {
	      dlnLdlz[j] = 0.0;
	      d2lnLdlz2[j] = 0.0;
	      for(i = 0; i < NumberOfThreads; i++)
		{
		  dlnLdlz[j]   += reductionBuffer[i * tr->NumberOfModels + j];
		  d2lnLdlz2[j] += reductionBufferTwo[i * tr->NumberOfModels + j];
		}
	    }
	}
#else
      if(firstIteration)
	{
	  makenewzIterative(tr);
	  firstIteration = FALSE;
	}
#ifdef _HET
      execCore(tr, dlnLdlz, d2lnLdlz2, isTipBranch);
#else
      execCore(tr, dlnLdlz, d2lnLdlz2);
#endif
#endif

      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && tr->curvatOK[i] == FALSE)
	    {
	      if ((d2lnLdlz2[i] >= 0.0) && (z[i] < zmax))
		zprev[i] = z[i] = 0.37 * z[i] + 0.63;  /*  Bad curvature, shorten branch */
	      else
		tr->curvatOK[i] = TRUE;
	    }
	}

      for(i = 0; i < numBranches; i++)
	{
	  if(tr->curvatOK[i] == TRUE && outerConverged[i] == FALSE)
	    {
	      if(d2lnLdlz2[i] < 0.0)
		{
		  double 
		    tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
		  
		  if(tantmp < 100)
		    {
		      z[i] *= EXP(tantmp);
		      if (z[i] < zmin)
			z[i] = zmin;

		      if (z[i] > 0.25 * zprev[i] + 0.75)
			z[i] = 0.25 * zprev[i] + 0.75;
		    }
		  else
		    z[i] = 0.25 * zprev[i] + 0.75;
		}
	      
	      if(z[i] > zmax) 
		z[i] = zmax;

	      maxiter[i] = maxiter[i] - 1;

	      /* the previous termination condition did not guarantee an improvement in LnL in 
		 cases where the initial branch was very long ! */
	      //if(maxiter[i] > 0 && (ABS(z[i] - zprev[i]) > zstep[i]))
	      if((ABS(z[i] - zprev[i]) > zstep[i]))
		{
		  /* We should make a more informed decision here,
		     based on the log like improvement */

		  if(maxiter[i] < -20)
		    {
		      z[i] = z0[i];
		      outerConverged[i] = TRUE;
		    }
		  else
		    outerConverged[i] = FALSE;
		}
	      else				
		outerConverged[i] = TRUE;						
	    }
	}

      loopConverged = TRUE;
      for(i = 0; i < numBranches; i++)
	loopConverged = loopConverged && outerConverged[i];
    }
  while (!loopConverged);

  for(model = 0; model < tr->NumberOfModels; model++)
    tr->executeModel[model] = TRUE;

  for(i = 0; i < numBranches; i++)
    result[i] = z[i];

#ifdef _BASTIEN
  for(i = 0; i < numBranches; i++)
    tr->secondDerivative[i] = d2lnLdlz2[i];
#endif 

}

#ifdef _USE_PTHREADS

static void sumClassify(tree *tr, int tipCase, double *_x1, double *_x2, unsigned char *_tipX1, unsigned char *_tipX2, boolean *executeModel)
{
  int 
    model = 0,
    columnCounter = 0, 
    offsetCounter = 0;
  
  double
    *x1_start = (double *)NULL,
    *x2_start = (double *)NULL;

  unsigned char 
    *tipX1 = (unsigned char*)NULL,
    *tipX2 = (unsigned char*)NULL;

#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("MNZS: ");
#endif

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;       
      
#ifdef _DEBUG_MULTI_EPA
      if(tr->threadID == THREAD_TO_DEBUG)
	printf("%d", executeModel[model]);
#endif       

      if(executeModel[model])
	{
	  double *sumBuffer = &tr->temporarySumBuffer[offsetCounter];

	  switch(tipCase)
	    {
	    case TIP_TIP:
	      tipX1 = &_tipX1[columnCounter];
	      tipX2 = &_tipX2[columnCounter];
	      break;
	    case TIP_INNER:
	      tipX1 = &_tipX1[columnCounter];
	      x2_start    = &_x2[offsetCounter];
	      break;
	    case INNER_INNER:
	      x1_start = &_x1[offsetCounter];
	      x2_start = &_x2[offsetCounter];
	      break;
	    default:
	      assert(0);	      
	    }

	  switch(tr->partitionData[model].dataType)
	    {
	    case BINARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:		
		  sumCAT_BINARY(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMA_BINARY(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
				  width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:		  
		  sumCAT(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
			 width);
		  break;
		case GAMMA:
		case GAMMA_I:		  
		  sumGAMMA(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector, tipX1, tipX2,
			   width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:		 
		  sumGTRCATPROT(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:		 
		  sumGAMMAPROT(tipCase,  sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
			       tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				     tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY(tipCase,  sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				    tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY_6(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				       tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY_6(tipCase,  sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				      tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
	      	{
		case CAT:
		  sumGTRCATSECONDARY_7(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				       tipX1, tipX2, width);
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGAMMASECONDARY_7(tipCase,  sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
				      tipX1, tipX2, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case GENERIC_32:
	      switch(tr->rateHetModel)
	      	{
		case CAT: 
		  sumCatFlex(tipCase, sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
			     tipX1, tipX2, width, tr->partitionData[model].states);		  
		  break;
		case GAMMA:
		case GAMMA_I:
		  sumGammaFlex(tipCase,  sumBuffer, x1_start, x2_start, tr->partitionData[model].tipVector,
			       tipX1, tipX2, width, tr->partitionData[model].states);		 
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
	}
      
      columnCounter += width;
      offsetCounter += width * tr->partitionData[model].states * tr->discreteRateCategories;
    } 
#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("\n");
#endif
}

static void coreClassify(tree *tr, volatile double *_dlnLdlz, volatile double *_d2lnLdlz2, double *coreLZ, boolean *executeModel)
{
  int 
    model, 
    branchIndex,
    offsetCounter = 0,
    columnCounter = 0;
  double lz;

  _dlnLdlz[0]   = 0.0;
  _d2lnLdlz2[0] = 0.0;
  
#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("MNZC: ");
#endif 
  
  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int 
	width = tr->partitionData[model].upper - tr->partitionData[model].lower;

      if(tr->multiBranch)
	branchIndex = model;
      else
	branchIndex = 0;	

     
#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("%d", executeModel[model]);
#endif    
      if(executeModel[model])
	{	  
	  double 
	    *sumBuffer       = &tr->temporarySumBuffer[offsetCounter],
	    *patrat          = tr->partitionData[model].perSiteRates;
	 
	  int 
	    *rateCategory  = &tr->contiguousRateCategory[columnCounter],
	    *wgt           = &tr->contiguousWgt[columnCounter],
	    *invariant     = &tr->contiguousInvariant[columnCounter];
	 

	  volatile double
	    dlnLdlz   = 0.0,
	    d2lnLdlz2 = 0.0;	  	  

	  if(tr->multiBranch)
	    {
	      branchIndex = model;
	      lz = coreLZ[model];
	      _dlnLdlz[model]   = 0.0;
	      _d2lnLdlz2[model] = 0.0;
	    }
	  else
	    {
	      branchIndex = 0;
	      lz = coreLZ[0];
	    }

	  switch(tr->partitionData[model].dataType)
	    {
	    case BINARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCAT_BINARY(width, tr->partitionData[model].numberOfCategories, sumBuffer,
				    &dlnLdlz, &d2lnLdlz2,
				    patrat, tr->partitionData[model].EIGN,  rateCategory, lz, wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMA_BINARY(width, sumBuffer,
				      &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
				      wgt);
		  break;
		case GAMMA_I:
		  coreGTRGAMMAINVAR_BINARY(tr->partitionData[model].propInvariant, tr->partitionData[model].frequencies,
					   tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					   sumBuffer, &dlnLdlz, &d2lnLdlz2,
					   invariant, wgt,
					   width, lz);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case DNA_DATA:	     
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCAT(width, tr->partitionData[model].numberOfCategories, sumBuffer,
			     &dlnLdlz, &d2lnLdlz2,
			     patrat, tr->partitionData[model].EIGN,  rateCategory, lz, wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMA(width, sumBuffer,
			       &dlnLdlz, &d2lnLdlz2, tr->partitionData[model].EIGN, tr->partitionData[model].gammaRates, lz,
			       wgt);
		  break;
		case GAMMA_I:
		  coreGTRGAMMAINVAR(tr->partitionData[model].propInvariant, tr->partitionData[model].frequencies,
				    tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				    sumBuffer, &dlnLdlz, &d2lnLdlz2,
				    invariant, wgt,
					width, lz);
		  break;
		default:
		  assert(0);
		}		
	      break;
	    case AA_DATA:	     
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATPROT(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  patrat,
				 rateCategory, width,				 
				 &dlnLdlz, &d2lnLdlz2,
				 sumBuffer, wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMAPROT(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				   sumBuffer, width, wgt,
				   &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMAPROTINVAR(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					sumBuffer, width, wgt,
					&dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					tr->partitionData[model].propInvariant, invariant);
		  break;
		default:
		  assert(0);
		}	    
	      break;
	    case SECONDARY_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  patrat,
				      rateCategory, width,				      
				      &dlnLdlz, &d2lnLdlz2,
				      sumBuffer, wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					sumBuffer, width, wgt,
					&dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					     sumBuffer, width, wgt,
					     &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					     tr->partitionData[model].propInvariant, invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_6:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY_6(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories, patrat,
					rateCategory, width,					
					&dlnLdlz, &d2lnLdlz2,
					sumBuffer, wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY_6(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					  sumBuffer, width, wgt,
					  &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR_6(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					       sumBuffer, width, wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					       tr->partitionData[model].propInvariant, invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case SECONDARY_DATA_7:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreGTRCATSECONDARY_7(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  patrat,
					rateCategory, width,					
					&dlnLdlz, &d2lnLdlz2,
					sumBuffer, wgt);
		  break;
		case GAMMA:
		  coreGTRGAMMASECONDARY_7(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					  sumBuffer, width, wgt,
					  &dlnLdlz, &d2lnLdlz2, lz);
		  break;
		case GAMMA_I:
		  coreGTRGAMMASECONDARYINVAR_7(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
					       sumBuffer, width, wgt,
					       &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
					       tr->partitionData[model].propInvariant, invariant);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case GENERIC_32:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  coreCatFlex(tr->partitionData[model].EIGN, lz, tr->partitionData[model].numberOfCategories,  patrat,
			      rateCategory, width,			      
			      &dlnLdlz, &d2lnLdlz2,
			      sumBuffer, tr->partitionData[model].states, wgt);
		  break;
		case GAMMA:
		  coreGammaFlex(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				sumBuffer, width, wgt,
				&dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].states);
		  break;
		case GAMMA_I:
		  coreGammaInvarFlex(tr->partitionData[model].gammaRates, tr->partitionData[model].EIGN,
				     sumBuffer, width, wgt,
				     &dlnLdlz, &d2lnLdlz2, lz, tr->partitionData[model].frequencies,
				     tr->partitionData[model].propInvariant, invariant, tr->partitionData[model].states);
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }

	  _dlnLdlz[branchIndex]   = _dlnLdlz[branchIndex]   + dlnLdlz;
	  _d2lnLdlz2[branchIndex] = _d2lnLdlz2[branchIndex] + d2lnLdlz2;
	}

      columnCounter += width;      
      offsetCounter += width * tr->partitionData[model].states * tr->discreteRateCategories;
    }

#ifdef _DEBUG_MULTI_EPA
  if(tr->threadID == THREAD_TO_DEBUG)
    printf("\n");
#endif

}

void makenewzClassify(tree *tr, int _maxiter, double *result, double *z0, double *x1_start, double *x2_start, unsigned char *tipX1,  
		      unsigned char *tipX2, int tipCase, boolean *partitionConverged, int insertion)
{
  double   
    z[NUM_BRANCHES], 
    zprev[NUM_BRANCHES], 
    zstep[NUM_BRANCHES],
    coreLZ[NUM_BRANCHES];
  volatile double  
    dlnLdlz[NUM_BRANCHES], 
    d2lnLdlz2[NUM_BRANCHES];
  
  int i, maxiter[NUM_BRANCHES], numBranches, model;
  
  boolean firstIteration = TRUE;
  boolean outerConverged[NUM_BRANCHES];
  boolean loopConverged,
    curvatOK[NUM_BRANCHES],
    executeModel[NUM_BRANCHES];

  

  if(tr->multiBranch)
    numBranches = tr->NumberOfModels;
  else
    numBranches = 1;

  

  for(i = 0; i < numBranches; i++)
    {
      z[i] = z0[i];
      maxiter[i] = _maxiter;
      outerConverged[i] = FALSE;
      curvatOK[i]       = TRUE;
      if(partitionConverged[i])
	executeModel[i] = FALSE;
      else
	executeModel[i] = TRUE;
    }
  
  if(tr->perPartitionEPA)    
    {
      setPartitionMask(tr, insertion, executeModel);

      if(tr->multiBranch)
	for(i = 0; i < numBranches; i++)
	  if(!executeModel[i])
	    outerConverged[i] = TRUE;
    }

      
  
  do
    {
      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && curvatOK[i] == TRUE)
	    {
	      curvatOK[i] = FALSE;

	      zprev[i] = z[i];

	      zstep[i] = (1.0 - zmax) * z[i] + zmin;
	    }
	}

      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && curvatOK[i] == FALSE)
	    {
	      double lz;

	      if (z[i] < zmin) z[i] = zmin;
	      else if (z[i] > zmax) z[i] = zmax;
	      lz    = log(z[i]);

	      coreLZ[i] = lz;
	    }
	}

      
      if(tr->multiBranch)
	{
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(executeModel[model])
		executeModel[model] = !curvatOK[model];	      
	    }
	}
      else
	{ 
	 
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->perPartitionEPA)
		{		 
		  if(executeModel[model])
		    {
		      executeModel[model] = !curvatOK[0];		  		    		     
		    }
		  /*else
		    executeModel[model] = FALSE;*/
		}
	      else
		executeModel[model] = !curvatOK[0];
	    }
	  
	}      

      
      if(firstIteration)
	{	  
	  sumClassify(tr, tipCase, x1_start, x2_start, tipX1, tipX2, executeModel);	  
	  firstIteration = FALSE;
	}	 
	
      coreClassify(tr, dlnLdlz, d2lnLdlz2, coreLZ, executeModel);	

      for(i = 0; i < numBranches; i++)
	{
	  if(outerConverged[i] == FALSE && curvatOK[i] == FALSE)
	    {	     
	      if ((d2lnLdlz2[i] >= 0.0) && (z[i] < zmax))
		zprev[i] = z[i] = 0.37 * z[i] + 0.63;  /*  Bad curvature, shorten branch */
	      else
		curvatOK[i] = TRUE;
	    }
	}

      for(i = 0; i < numBranches; i++)
	{
	  if(curvatOK[i] == TRUE && outerConverged[i] == FALSE)
	    {
	      if (d2lnLdlz2[i] < 0.0)
		{
		  double tantmp = -dlnLdlz[i] / d2lnLdlz2[i];
		  if (tantmp < 100)
		    {
		      z[i] *= EXP(tantmp);
		      if (z[i] < zmin)
			z[i] = zmin;

		      if (z[i] > 0.25 * zprev[i] + 0.75)
			z[i] = 0.25 * zprev[i] + 0.75;
		    }
		  else
		    z[i] = 0.25 * zprev[i] + 0.75;
		}
	      if (z[i] > zmax) z[i] = zmax;

	      maxiter[i] = maxiter[i] - 1;
	      if(maxiter[i] > 0 && (ABS(z[i] - zprev[i]) > zstep[i]))
		outerConverged[i] = FALSE;
	      else
		outerConverged[i] = TRUE;
	    }
	}

      loopConverged = TRUE;
      for(i = 0; i < numBranches; i++)
	loopConverged = loopConverged && outerConverged[i];
    }
  while (!loopConverged);   
  
  for(i = 0; i < numBranches; i++)   
    result[i] = z[i];  
}

#endif

void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result, boolean mask)
{
  int 
    i; 

#ifdef _HET
  boolean 
    isTipBranch;

  if(isTip(p->number, tr->mxtips) || isTip(q->number, tr->mxtips))
    isTipBranch = TRUE;
  else
    isTipBranch = FALSE;
#endif

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < tr->numBranches; i++)
    {     
      tr->td[0].ti[0].qz[i] =  z0[i];

      if(mask)
	{	     
	  if(tr->partitionConverged[i])
	    tr->executeModel[i] = FALSE;
	  else
	    tr->executeModel[i] = TRUE;
	}
    }
    
  tr->td[0].count = 1;
    
  if(!p->x)
    computeTraversalInfo(tr, p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);

  if(!q->x)
    computeTraversalInfo(tr, q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
   
#ifdef _HET
  topLevelMakenewz(tr, z0, maxiter, result, isTipBranch);
#else
  topLevelMakenewz(tr, z0, maxiter, result);
#endif
  
  for(i = 0; i < tr->numBranches; i++)
      tr->executeModel[i] = TRUE;
}



void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2)
{
#ifdef _HET
  assert(0);//not implemented
#else 
  int 
    i;

  assert(taxon1 != taxon2);
  assert(0 < taxon1 && taxon1 <= tr->mxtips);
  assert(0 < taxon2 && taxon2 <= tr->mxtips);

  tr->td[0].ti[0].pNumber = taxon1;
  tr->td[0].ti[0].qNumber = taxon2;
  tr->td[0].ti[0].tipCase = TIP_TIP;

  for(i = 0; i < tr->numBranches; i++)
    tr->td[0].ti[0].qz[i] =  defaultz;      /* TODO why defaultz ? */

  tr->td[0].count = 1;
 
  topLevelMakenewz(tr, z0, maxiter, result);
#endif
}



