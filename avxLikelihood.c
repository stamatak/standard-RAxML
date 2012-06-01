#include <unistd.h>

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include "axml.h"
#include <stdint.h>
#include <xmmintrin.h>
#include <pmmintrin.h>
#include <immintrin.h>

#ifdef _FMA
#include <x86intrin.h>
#define FMAMACC(a,b,c) _mm256_macc_pd(b,c,a)
#endif

const union __attribute__ ((aligned (BYTE_ALIGNMENT)))
{
  uint64_t i[4];
  __m256d m;
  
} absMask_AVX = {{0x7fffffffffffffffULL, 0x7fffffffffffffffULL, 0x7fffffffffffffffULL, 0x7fffffffffffffffULL}};



static inline __m256d hadd4(__m256d v, __m256d u)
{ 
  __m256d
    a, b;
  
  v = _mm256_hadd_pd(v, v);
  a = _mm256_permute2f128_pd(v, v, 1);
  v = _mm256_add_pd(a, v);

  u = _mm256_hadd_pd(u, u);
  b = _mm256_permute2f128_pd(u, u, 1);
  u = _mm256_add_pd(b, u);

  v = _mm256_mul_pd(v, u);	
  
  return v;
}

static inline __m256d hadd3(__m256d v)
{ 
  __m256d
    a;
  
  v = _mm256_hadd_pd(v, v);
  a = _mm256_permute2f128_pd(v, v, 1);
  v = _mm256_add_pd(a, v);
  
  return v;
}


void  newviewGTRGAMMA_AVX(int tipCase,
			 double *x1, double *x2, double *x3,
			 double *extEV, double *tipVector,
			 int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			 const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
			 )
{
 
  int  
    i, 
    k, 
    scale, 
    addScale = 0;
 
  __m256d 
    minlikelihood_avx = _mm256_set1_pd( minlikelihood ),
    twoto = _mm256_set1_pd(twotothe256);
 

  switch(tipCase)
    {
    case TIP_TIP:
      {
	double 
	  *uX1, 
	  umpX1[1024] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
	  *uX2, 
	  umpX2[1024] __attribute__ ((aligned (BYTE_ALIGNMENT)));

	for (i = 1; i < 16; i++)
	  {
	    __m256d 
	      tv = _mm256_load_pd(&(tipVector[i * 4]));

	    int 
	      j;
	    
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&left[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX1[i * 64 + j * 16 + k * 4], left1);
		}
	  
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&right[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX2[i * 64 + j * 16 + k * 4], left1);
		}	    
	  }   	
	  

	for(i = 0; i < n; i++)
	  {	    		 	    
	    uX1 = &umpX1[64 * tipX1[i]];
	    uX2 = &umpX2[64 * tipX2[i]];		  
	    
	    for(k = 0; k < 4; k++)
	      {
		__m256d	   
		  xv = _mm256_setzero_pd();
	       
		int 
		  l;
		
		for(l = 0; l < 4; l++)
		  {	       	     				      	      																	   
		    __m256d
		      x1v =  _mm256_mul_pd(_mm256_load_pd(&uX1[k * 16 + l * 4]), _mm256_load_pd(&uX2[k * 16 + l * 4]));
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
#ifdef _FMA
		    xv = FMAMACC(xv,x1v,evv);
#else						  
		    xv = _mm256_add_pd(xv, _mm256_mul_pd(x1v, evv));
#endif
		  }
		
		_mm256_store_pd(&x3[16 * i + 4 * k], xv);
	      }	         	   	    
	  }
      }
      break;
    case TIP_INNER:
      {
	double 
	  *uX1, 
	  umpX1[1024] __attribute__ ((aligned (BYTE_ALIGNMENT)));

	for (i = 1; i < 16; i++)
	  {
	    __m256d 
	      tv = _mm256_load_pd(&(tipVector[i*4]));

	    int 
	      j;
	    
	    for (j = 0; j < 4; j++)
	      for (k = 0; k < 4; k++)
		{		 
		  __m256d 
		    left1 = _mm256_load_pd(&left[j * 16 + k * 4]);		  		  		  

		  left1 = _mm256_mul_pd(left1, tv);		  
		  left1 = hadd3(left1);
		  		  		  
		  _mm256_store_pd(&umpX1[i * 64 + j * 16 + k * 4], left1);
		}	 	   
	  }   	
	
	for(i = 0; i < n; i++)
	  { 
	    __m256d
	      xv[4];	    	   
	    
	    scale = 1;
	    uX1 = &umpX1[64 * tipX1[i]];

	    for(k = 0; k < 4; k++)
	      {
		__m256d	   		 
		  xvr = _mm256_load_pd(&(x2[i * 16 + k * 4]));

		int 
		  l;

		xv[k]  = _mm256_setzero_pd();
		  
		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d  
		      x1v = _mm256_load_pd(&uX1[k * 16 + l * 4]),		     
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x2v = hadd3(x2v);
		    x1v = _mm256_mul_pd(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
			
#ifdef _FMA
		    xv[k] = FMAMACC(xv[k],x1v,evv);
#else			  
		    xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
#endif
		  }
		    
		if(scale)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }
	      }	    

	    if(scale)
	      {
		xv[0] = _mm256_mul_pd(xv[0], twoto);
		xv[1] = _mm256_mul_pd(xv[1], twoto);
		xv[2] = _mm256_mul_pd(xv[2], twoto);
		xv[3] = _mm256_mul_pd(xv[3], twoto);
		addScale += wgt[i];
	      }

	    _mm256_store_pd(&x3[16 * i],      xv[0]);
	    _mm256_store_pd(&x3[16 * i + 4],  xv[1]);
	    _mm256_store_pd(&x3[16 * i + 8],  xv[2]);
	    _mm256_store_pd(&x3[16 * i + 12], xv[3]);
	  }
      }
      break;
    case INNER_INNER:
      {
	for(i = 0; i < n; i++)
	  {	
	    __m256d
	      xv[4];
	    
	    scale = 1;

	    for(k = 0; k < 4; k++)
	      {
		__m256d	   
		 
		  xvl = _mm256_load_pd(&(x1[i * 16 + k * 4])),
		  xvr = _mm256_load_pd(&(x2[i * 16 + k * 4]));

		int 
		  l;

		xv[k] = _mm256_setzero_pd();

		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d 
		      x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x1v = hadd4(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
						  
		    xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
		  }
		
		if(scale)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }
	      }

	     if(scale)
	      {
		xv[0] = _mm256_mul_pd(xv[0], twoto);
		xv[1] = _mm256_mul_pd(xv[1], twoto);
		xv[2] = _mm256_mul_pd(xv[2], twoto);
		xv[3] = _mm256_mul_pd(xv[3], twoto);
		addScale += wgt[i];
	      }
		
	    _mm256_store_pd(&x3[16 * i],      xv[0]);
	    _mm256_store_pd(&x3[16 * i + 4],  xv[1]);
	    _mm256_store_pd(&x3[16 * i + 8],  xv[2]);
	    _mm256_store_pd(&x3[16 * i + 12], xv[3]);
	  }
      }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
  
}


/*
static void  newviewGTRGAMMA_AVX2(int tipCase,
			 double *x1, double *x2, double *x3,
			 double *extEV, double *tipVector,
			 int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			 const int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling
			 )
{
 
  int  
    i, 
    k, 
    scale, 
    addScale = 0;
 
  __m256d 
    minlikelihood_avx = _mm256_set1_pd( minlikelihood ),
    twoto = _mm256_set1_pd(twotothe256);
 

  switch(tipCase)
    {
    case TIP_TIP:
      {
	for(i = 0; i < n; i++)
	  {
	    __m256d
	      xvl = _mm256_load_pd(&(tipVector[4 * tipX1[i]])),
	      xvr = _mm256_load_pd(&(tipVector[4 * tipX2[i]]));
				  
	    for(k = 0; k < 4; k++)
	      {
		__m256d	   
		  xv = _mm256_setzero_pd();

		int 
		  l;

		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d 
		      x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x1v = hadd4(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
						  
		    xv = _mm256_add_pd(xv, _mm256_mul_pd(x1v, evv));
		  }
		
		_mm256_store_pd(&x3[16 * i + 4 * k], xv);
	      }	         	   	    
	  }
      }
      break;
    case TIP_INNER:
      {
	for(i = 0; i < n; i++)
	  { 
	    __m256d
	      xv[4];
	    
	    __m256d
	      xvl = _mm256_load_pd(&(tipVector[4 * tipX1[i]]));
	    
	    scale = 1;

	    for(k = 0; k < 4; k++)
	      {
		__m256d	   		 
		  xvr = _mm256_load_pd(&(x2[i * 16 + k * 4]));

		int 
		  l;

		xv[k]  = _mm256_setzero_pd();
		  
		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d 
		      x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x1v = hadd4(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
						  
		    xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
		  }
		    
		if(scale)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }
	      }	    

	    if(scale)
	      {
		xv[0] = _mm256_mul_pd(xv[0], twoto);
		xv[1] = _mm256_mul_pd(xv[1], twoto);
		xv[2] = _mm256_mul_pd(xv[2], twoto);
		xv[3] = _mm256_mul_pd(xv[3], twoto);
		addScale += wgt[i];
	      }

	    _mm256_store_pd(&x3[16 * i],      xv[0]);
	    _mm256_store_pd(&x3[16 * i + 4],  xv[1]);
	    _mm256_store_pd(&x3[16 * i + 8],  xv[2]);
	    _mm256_store_pd(&x3[16 * i + 12], xv[3]);
	  }
      }
      break;
    case INNER_INNER:
      {
	for(i = 0; i < n; i++)
	  {	
	    __m256d
	      xv[4];
	    
	    scale = 1;

	    for(k = 0; k < 4; k++)
	      {
		__m256d	   
		 
		  xvl = _mm256_load_pd(&(x1[i * 16 + k * 4])),
		  xvr = _mm256_load_pd(&(x2[i * 16 + k * 4]));

		int 
		  l;

		xv[k] = _mm256_setzero_pd();

		for(l = 0; l < 4; l++)
		  {	       	     				      	      															
		    __m256d 
		      x1v = _mm256_mul_pd(xvl, _mm256_load_pd(&left[k * 16 + l * 4])),
		      x2v = _mm256_mul_pd(xvr, _mm256_load_pd(&right[k * 16 + l * 4]));			    
			
		    x1v = hadd4(x1v, x2v);			
		
		    __m256d 
		      evv = _mm256_load_pd(&extEV[l * 4]);
						  
		    xv[k] = _mm256_add_pd(xv[k], _mm256_mul_pd(x1v, evv));
		  }
		
		if(scale)
		  {
		    __m256d 	     
		      v1 = _mm256_and_pd(xv[k], absMask_AVX.m);

		    v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
		    
		    if(_mm256_movemask_pd( v1 ) != 15)
		      scale = 0;
		  }
	      }

	     if(scale)
	      {
		xv[0] = _mm256_mul_pd(xv[0], twoto);
		xv[1] = _mm256_mul_pd(xv[1], twoto);
		xv[2] = _mm256_mul_pd(xv[2], twoto);
		xv[3] = _mm256_mul_pd(xv[3], twoto);
		addScale += wgt[i];
	      }
		
	    _mm256_store_pd(&x3[16 * i],      xv[0]);
	    _mm256_store_pd(&x3[16 * i + 4],  xv[1]);
	    _mm256_store_pd(&x3[16 * i + 8],  xv[2]);
	    _mm256_store_pd(&x3[16 * i + 12], xv[3]);
	  }
      }
      break;
    default:
      assert(0);
    }

  if(useFastScaling)
    *scalerIncrement = addScale;
  
}
*/

void newviewGTRCAT_AVX(int tipCase,  double *EV,  int *cptr,
			   double *x1_start, double *x2_start,  double *x3_start, double *tipVector,
			   int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			   int n,  double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le,
    *ri,
    *x1,
    *x2, 
    *x3;
    
  int 
    i, 
    j, 
    scale, 
    addScale = 0;
   
  __m256d 
    minlikelihood_avx = _mm256_set1_pd( minlikelihood ),
    twoto = _mm256_set1_pd(twotothe256);
  
  switch(tipCase)
    {
    case TIP_TIP:      
      for (i = 0; i < n; i++)
	{	 
	  int 
	    l;
	  
	  le = &left[cptr[i] * 16];
	  ri = &right[cptr[i] * 16];

	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  
	  __m256d	   
	    vv = _mm256_setzero_pd();
	   	   	    
	  for(l = 0; l < 4; l++)
	    {	       	     				      	      															
	      __m256d 
		x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
			
	      x1v = hadd4(x1v, x2v);			
		
	      __m256d 
		evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
	      vv = FMAMACC(vv,x1v,evv);
#else				
	      vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
	    }	  		  

	  _mm256_store_pd(&x3_start[4 * i], vv);	    	   	    
	}
      break;
    case TIP_INNER:      
      for (i = 0; i < n; i++)
	{
	  int 
	    l;

	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];	 
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m256d	   
	    vv = _mm256_setzero_pd();
	  
	  for(l = 0; l < 4; l++)
	    {	       	     				      	      															
	      __m256d 
		x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
			
	      x1v = hadd4(x1v, x2v);			
		
	      __m256d 
		evv = _mm256_load_pd(&EV[l * 4]);
				
#ifdef _FMA
	      vv = FMAMACC(vv,x1v,evv);
#else	      
	      vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));
#endif
	    }	  		  
	  
	  
	  __m256d 	     
	    v1 = _mm256_and_pd(vv, absMask_AVX.m);

	  v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	    
	  if(_mm256_movemask_pd( v1 ) == 15)
	    {	     	      
	      vv = _mm256_mul_pd(vv, twoto);	      
	      addScale += wgt[i];
	    }       
	  
	  _mm256_store_pd(&x3_start[4 * i], vv);	 	  	  
	}
      break;
    case INNER_INNER:
      for (i = 0; i < n; i++)
	{
	  int 
	    l;

	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  
	  le =  &left[cptr[i] * 16];
	  ri =  &right[cptr[i] * 16];

	  __m256d	   
	    vv = _mm256_setzero_pd();
	  
	  for(l = 0; l < 4; l++)
	    {	       	     				      	      															
	      __m256d 
		x1v = _mm256_mul_pd(_mm256_load_pd(x1), _mm256_load_pd(&le[l * 4])),
		x2v = _mm256_mul_pd(_mm256_load_pd(x2), _mm256_load_pd(&ri[l * 4]));			    
			
	      x1v = hadd4(x1v, x2v);			
		
	      __m256d 
		evv = _mm256_load_pd(&EV[l * 4]);
#ifdef _FMA
	      vv = FMAMACC(vv,x1v,evv);
#else						
	      vv = _mm256_add_pd(vv, _mm256_mul_pd(x1v, evv));						      	
#endif
	    }	  		  

	 
	  __m256d 	     
	    v1 = _mm256_and_pd(vv, absMask_AVX.m);

	  v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	    
	  if(_mm256_movemask_pd( v1 ) == 15)
	    {	
	      vv = _mm256_mul_pd(vv, twoto);	      
	      addScale += wgt[i];
	    }	

	  _mm256_store_pd(&x3_start[4 * i], vv);
	  	  
	}
      break;
    default:
      assert(0);
    }

  
  *scalerIncrement = addScale;
}

void newviewGTRCATPROT_AVX(int tipCase, double *extEV,
			       int *cptr,
			       double *x1, double *x2, double *x3, double *tipVector,
			       int *ex3, unsigned char *tipX1, unsigned char *tipX2,
			       int n, double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling)
{
  double
    *le, *ri, *v, *vl, *vr;

  int i, l, j, scale, addScale = 0;

#ifdef _FMA
  int k;
#endif

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

	    __m256d vv[5];
	    
	    vv[0] = _mm256_setzero_pd();
	    vv[1] = _mm256_setzero_pd();
	    vv[2] = _mm256_setzero_pd();
	    vv[3] = _mm256_setzero_pd();
	    vv[4] = _mm256_setzero_pd();	   	    

	    for(l = 0; l < 20; l++)
	      {	       
		__m256d 
		  x1v = _mm256_setzero_pd(),
		  x2v = _mm256_setzero_pd();	
				
		double 
		  *ev = &extEV[l * 20],
		  *lv = &le[l * 20],
		  *rv = &ri[l * 20];														

#ifdef _FMA		
		for(k = 0; k < 20; k += 4) 
		  {
		    __m256d vlv = _mm256_load_pd(&vl[k]);
		    __m256d lvv = _mm256_load_pd(&lv[k]);
		    x1v = FMAMACC(x1v,vlv,lvv);
		    __m256d vrv = _mm256_load_pd(&vr[k]);
		    __m256d rvv = _mm256_load_pd(&rv[k]);
		    x2v = FMAMACC(x2v,vrv,rvv);
		  }
#else		
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
		x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));

		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
		x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));	
#endif

		x1v = hadd4(x1v, x2v);			
#ifdef _FMA
		for(k = 0; k < 5; k++) 
		  {
		    __m256d evv = _mm256_load_pd(&ev[k*4]);
		    vv[k] = FMAMACC(vv[k],x1v,evv);
		  }	  
#else		
		__m256d 
		  evv[5];
	    	
		evv[0] = _mm256_load_pd(&ev[0]);
		evv[1] = _mm256_load_pd(&ev[4]);
		evv[2] = _mm256_load_pd(&ev[8]);
		evv[3] = _mm256_load_pd(&ev[12]);
		evv[4] = _mm256_load_pd(&ev[16]);		
		
		vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
		vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
		vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
		vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
		vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      		      	  
#endif
	      }
	    _mm256_store_pd(&v[0], vv[0]);
	    _mm256_store_pd(&v[4], vv[1]);
	    _mm256_store_pd(&v[8], vv[2]);
	    _mm256_store_pd(&v[12], vv[3]);
	    _mm256_store_pd(&v[16], vv[4]);
	  }
      }
      break;
    case TIP_INNER:      	
      for (i = 0; i < n; i++)
	{
	  le = &left[cptr[i] * 400];
	  ri = &right[cptr[i] * 400];
	  
	  vl = &(tipVector[20 * tipX1[i]]);
	  vr = &x2[20 * i];
	  v  = &x3[20 * i];	   
	  
	  __m256d vv[5];
	  
	  vv[0] = _mm256_setzero_pd();
	  vv[1] = _mm256_setzero_pd();
	  vv[2] = _mm256_setzero_pd();
	  vv[3] = _mm256_setzero_pd();
	  vv[4] = _mm256_setzero_pd();
	  
	 

	  for(l = 0; l < 20; l++)
	    {	       
	      __m256d 
		x1v = _mm256_setzero_pd(),
		x2v = _mm256_setzero_pd();	
	      
	      double 
		*ev = &extEV[l * 20],
		*lv = &le[l * 20],
		*rv = &ri[l * 20];														
#ifdef _FMA
	      for(k = 0; k < 20; k += 4) 
		{
		  __m256d vlv = _mm256_load_pd(&vl[k]);
		  __m256d lvv = _mm256_load_pd(&lv[k]);
		  x1v = FMAMACC(x1v,vlv,lvv);
		  __m256d vrv = _mm256_load_pd(&vr[k]);
		  __m256d rvv = _mm256_load_pd(&rv[k]);
		  x2v = FMAMACC(x2v,vrv,rvv);
		}
#else	      
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
	      
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));
#endif

	      x1v = hadd4(x1v, x2v);			
	      
	      __m256d 
		evv[5];
	      
	      evv[0] = _mm256_load_pd(&ev[0]);
	      evv[1] = _mm256_load_pd(&ev[4]);
	      evv[2] = _mm256_load_pd(&ev[8]);
	      evv[3] = _mm256_load_pd(&ev[12]);
	      evv[4] = _mm256_load_pd(&ev[16]);		

#ifdef _FMA
	      for(k = 0; k < 5; k++)
		vv[k] = FMAMACC(vv[k],x1v,evv[k]);		 
#else	      
	      vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
	      vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
	      vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
	      vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
	      vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
	    }	  

	   	     
	  __m256d minlikelihood_avx = _mm256_set1_pd( minlikelihood );
	  
	  scale = 1;
	  
	  for(l = 0; scale && (l < 20); l += 4)
	    {	       
	      __m256d 
		v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
	      v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	      
	      if(_mm256_movemask_pd( v1 ) != 15)
		scale = 0;
	    }	    	  	  
	 

	  if(scale)
	    {
	      __m256d 
		twoto = _mm256_set1_pd(twotothe256);
	      
	      for(l = 0; l < 20; l += 4)
		vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 
	  
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	      
	    }

	  _mm256_store_pd(&v[0], vv[0]);
	  _mm256_store_pd(&v[4], vv[1]);
	  _mm256_store_pd(&v[8], vv[2]);
	  _mm256_store_pd(&v[12], vv[3]);
	  _mm256_store_pd(&v[16], vv[4]);	       
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

	  __m256d vv[5];
	  
	  vv[0] = _mm256_setzero_pd();
	  vv[1] = _mm256_setzero_pd();
	  vv[2] = _mm256_setzero_pd();
	  vv[3] = _mm256_setzero_pd();
	  vv[4] = _mm256_setzero_pd();
	  
	  for(l = 0; l < 20; l++)
	    {	       
	      __m256d 
		x1v = _mm256_setzero_pd(),
		x2v = _mm256_setzero_pd();	
	      
	      double 
		*ev = &extEV[l * 20],
		*lv = &le[l * 20],
		*rv = &ri[l * 20];														
	      
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[0]), _mm256_load_pd(&lv[0])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[4]), _mm256_load_pd(&lv[4])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[8]), _mm256_load_pd(&lv[8])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[12]), _mm256_load_pd(&lv[12])));
	      x1v = _mm256_add_pd(x1v, _mm256_mul_pd(_mm256_load_pd(&vl[16]), _mm256_load_pd(&lv[16])));
	      
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[0]), _mm256_load_pd(&rv[0])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[4]), _mm256_load_pd(&rv[4])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[8]), _mm256_load_pd(&rv[8])));			    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[12]), _mm256_load_pd(&rv[12])));				    
	      x2v = _mm256_add_pd(x2v,  _mm256_mul_pd(_mm256_load_pd(&vr[16]), _mm256_load_pd(&rv[16])));

	      x1v = hadd4(x1v, x2v);			
#ifdef _FMA
	       for(k = 0; k < 5; k++) 
		 {
		   __m256d evv = _mm256_load_pd(&ev[k*4]);
		   vv[k] = FMAMACC(vv[k],x1v,evv);
		 }
#else	      
	      __m256d 
		evv[5];
	      
	      evv[0] = _mm256_load_pd(&ev[0]);
	      evv[1] = _mm256_load_pd(&ev[4]);
	      evv[2] = _mm256_load_pd(&ev[8]);
	      evv[3] = _mm256_load_pd(&ev[12]);
	      evv[4] = _mm256_load_pd(&ev[16]);		
	      
	      vv[0] = _mm256_add_pd(vv[0], _mm256_mul_pd(x1v, evv[0]));
	      vv[1] = _mm256_add_pd(vv[1], _mm256_mul_pd(x1v, evv[1]));
	      vv[2] = _mm256_add_pd(vv[2], _mm256_mul_pd(x1v, evv[2]));
	      vv[3] = _mm256_add_pd(vv[3], _mm256_mul_pd(x1v, evv[3]));
	      vv[4] = _mm256_add_pd(vv[4], _mm256_mul_pd(x1v, evv[4]));				      	
#endif
	    }	  

	   	     
	  __m256d minlikelihood_avx = _mm256_set1_pd( minlikelihood );
	  
	  scale = 1;
	  
	  for(l = 0; scale && (l < 20); l += 4)
	    {	       
	      __m256d 
		v1 = _mm256_and_pd(vv[l / 4], absMask_AVX.m);
	      v1 = _mm256_cmp_pd(v1,  minlikelihood_avx, _CMP_LT_OS);
	      
	      if(_mm256_movemask_pd( v1 ) != 15)
		scale = 0;
	    }	    	  	  

	  if(scale)
	    {
	      __m256d 
		twoto = _mm256_set1_pd(twotothe256);
	      
	      for(l = 0; l < 20; l += 4)
		vv[l / 4] = _mm256_mul_pd(vv[l / 4] , twoto);		    		 
	  
	      if(useFastScaling)
		addScale += wgt[i];
	      else
		ex3[i]  += 1;	      
	    }

	  _mm256_store_pd(&v[0], vv[0]);
	  _mm256_store_pd(&v[4], vv[1]);
	  _mm256_store_pd(&v[8], vv[2]);
	  _mm256_store_pd(&v[12], vv[3]);
	  _mm256_store_pd(&v[16], vv[4]);
	 
	}
      break;
    default:
      assert(0);
    }
  
  if(useFastScaling)
    *scalerIncrement = addScale;
}

static void newviewGTRGAMMAPROT_AVX2(int tipCase,
			     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			     int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
			     double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling) 
{
  double	
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr;
  
  int	
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0;

 
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif


#if GCC_VERSION < 40500
   __m256d
    bitmask = _mm256_set_pd(0,0,0,-1);
#else
  __m256i
    bitmask = _mm256_set_epi32(0, 0, 0, 0, 0, 0, -1, -1);
#endif



  /* this is required for doing some pre-computations that help to save 
     numerical operations. What we are actually computing here are additional lookup tables 
     for each possible state a certain data-type can assume.
     for DNA with ambuguity coding this is 15, for proteins this is 22 or 23, since there 
     also exist one or two amibguity codes for protein data.
     Essentially this is very similar to the tip vectors which we also use as lookup tables */
  
  switch(tipCase) 
    {
    case TIP_TIP: 
      {
	/* allocate pre-compute memory space */
	double umpX1[1840], umpX2[1840];
	/* multiply all possible tip state vectors with the respective P-matrices 
	 * maxStateValue=23,span=80,states=20  for the protein model (see above) */
	for(i = 0; i < 23 ; i++) 
	  {
	    v = &(tipVector[20 * i]);
	    for(k = 0; k < 80; k++) 
	      {
		double *ll =  &left[k * 20];
		double *rr =  &right[k * 20];
		__m256d umpX1v = _mm256_setzero_pd();
		__m256d umpX2v = _mm256_setzero_pd();
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
#ifdef _FMA
		    __m256d llv = _mm256_load_pd(&ll[l]);
		    umpX1v = FMAMACC(umpX1v,vv,llv);
		    __m256d rrv = _mm256_load_pd(&rr[l]);
		    umpX2v = FMAMACC(umpX2v,vv,rrv);
#else		    
		    umpX1v = _mm256_add_pd(umpX1v,_mm256_mul_pd(vv,_mm256_load_pd(&ll[l])));
		    umpX2v = _mm256_add_pd(umpX2v,_mm256_mul_pd(vv,_mm256_load_pd(&rr[l])));
#endif
		  }
		
		umpX1v = hadd3(umpX1v);
		umpX2v = hadd3(umpX2v);
		_mm256_maskstore_pd(&umpX1[80 * i + k],bitmask, umpX1v);
		_mm256_maskstore_pd(&umpX2[80 * i + k],bitmask, umpX2v);
	      } 
	  }
	for(i = 0; i < n; i++) 
	  {
	    /* access the precomputed arrays (pre-computed multiplication of conditional with the tip state) */
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];
	    /* loop over discrete GAMMA rates */
	    for(j = 0; j < 4; j++) 
	      {
		/* the rest is the same as for CAT */
		v = &x3[i * 80 + j * 20];
		__m256d zero = _mm256_setzero_pd();
		for(k = 0; k < 20; k+=4)
		  _mm256_store_pd(&v[k],zero);
		
		for(k = 0; k < 20; k++) 
		  {			 
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];
		    __m256d x1px2v = _mm256_set1_pd(x1px2);
		    
		    for(l = 0; l < 20; l+=4) 
		      {
			__m256d vv = _mm256_load_pd(&v[l]);
			__m256d extEvv = _mm256_load_pd(&extEV[20 * k + l]);
#ifdef _FMA
			vv = FMAMACC(vv,x1px2v,extEvv);
#else
			vv = _mm256_add_pd(vv,_mm256_mul_pd(x1px2v,extEvv));
#endif
			_mm256_store_pd(&v[l],vv);
		      } 
		  } 
	      } 
	  } 
      } 
      break;
    case TIP_INNER: 
      {
	/* we do analogous pre-computations as above, with the only difference that we now do them only for one tip vector 
	   double *umpX1 = (double*)malloc(sizeof(double) * precomputeLength),
	   *ump_x2 = (double*)malloc(sizeof(double) * states); */
	double 
	  umpX1[1840], 
	  ump_x2[20];
	/* precompute P and left tip vector product */
	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);
	    for(k = 0; k < 80; k++) 
	      {
		__m256d umpX1v = _mm256_setzero_pd();
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);

#ifdef _FMA
		    __m256d leftv = _mm256_load_pd(&left[k * 20 + l]);
		    umpX1v = FMAMACC(umpX1v,vv,leftv);
#else
		    umpX1v = _mm256_add_pd(umpX1v,_mm256_mul_pd(vv,_mm256_load_pd(&left[k * 20 + l])));
#endif
		  }
		umpX1v = hadd3(umpX1v);
		_mm256_maskstore_pd(&umpX1[80 * i + k],bitmask,umpX1v);
	      } 
	  }
	
	for (i = 0; i < n; i++) 
	  {
	    /* access pre-computed value based on the raw sequence data tipX1 that is used as an index */
	    uX1 = &umpX1[80 * tipX1[i]];
	    /* loop over discrete GAMMA rates */
	    for(k = 0; k < 4; k++) 
	      {
		v = &(x2[80 * i + k * 20]);
		for(l = 0; l < 20; l++) 
		  {
		    __m256d ump_x2v = _mm256_setzero_pd();
		    for(j = 0; j < 20; j += 4) 
		      {
#ifdef _FMA
			__m256d vv = _mm256_load_pd(&v[j]);
			__m256d rightv = _mm256_load_pd(&right[k*400+l*20+j]);
			ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
			ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(_mm256_load_pd(&v[j]),_mm256_load_pd(&right[k*400+l*20+j])));
#endif
		      } 
		    ump_x2v = hadd3(ump_x2v);
		    _mm256_maskstore_pd(&ump_x2[l],bitmask,ump_x2v);
		  }
		v = &(x3[80 * i + 20 * k]);
		__m256d zero = _mm256_setzero_pd();
		for(l = 0; l < 20; l+=4)
		  _mm256_store_pd(&v[l],zero);
		
		for(l = 0; l < 20; l++) 
		  {
		    x1px2 = uX1[k * 20 + l]	* ump_x2[l];
		    __m256d x1px2v = _mm256_set1_pd(x1px2);
		    
		    for(j = 0; j < 20; j+=4) 
		      {
			/*v[j] += x1px2 * extEV[l * states	+ j];*/
			__m256d vv = _mm256_load_pd(&v[j]);
#ifdef _FMA
			__m256d extEVv =_mm256_load_pd(&extEV[l * 20 + j]);
			vv = FMAMACC(vv,x1px2v,extEVv);
#else
			vv = _mm256_add_pd(vv,_mm256_mul_pd(x1px2v,_mm256_load_pd(&extEV[l * 20 + j])));
#endif

			_mm256_store_pd(&v[j],vv);
		      }
		  } 
	      }
	    /* also do numerical scaling as above. Note that here we need to scale 
	       4 * 4 values for DNA or 4 * 20 values for protein data.
	       If they are ALL smaller than our threshold, we scale. Note that,
	       this can cause numerical problems with GAMMA, if the values generated 
	       by the four discrete GAMMA rates are too different.
	       For details, see: 
	       F. Izquierdo-Carrasco, S.A. Smith, A. Stamatakis: "Algorithms, Data Structures, and Numerics for Likelihood-based Phylogenetic Inference of Huge Trees"
	    */
	    v = &x3[80 * i];
	    __m256d minlikelihood_avx = _mm256_set1_pd(minlikelihood);
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 4) 
	      {
		__m256d vv = _mm256_load_pd(&v[l]);
		__m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		if(_mm256_movemask_pd(vv_abs) != 15)
		  scale = 0;
	      }
	    
	    if (scale) 
	      {
		__m256d twotothe256v = _mm256_set_pd(twotothe256,twotothe256,twotothe256,twotothe256);
		for(l = 0; l < 80; l += 4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    _mm256_store_pd(&v[l],_mm256_mul_pd(vv,twotothe256v));
		  }
		if(useFastScaling)
		  addScale += wgt[i];				
		else
		  ex3[i] += 1;
	      } 
	  } 
      } 
      break;
    case INNER_INNER:
      /* same as above, without pre-computations */
      for(i = 0; i < n; i++) 
	{
	  for(k = 0; k < 4; k++) 
	    {
	      vl = &(x1[80 * i + 20 * k]);
	      vr = &(x2[80 * i + 20 * k]);
	      v  = &(x3[80 * i + 20 * k]);
	      __m256d zero = _mm256_setzero_pd();
	      
	      for(l = 0; l < 20; l+=4)
		_mm256_store_pd(&v[l],zero);
	      
	      for(l = 0; l < 20; l++) 
		{		 
		  __m256d al = _mm256_setzero_pd();
		  __m256d ar = _mm256_setzero_pd();
		  for(j = 0; j < 20; j += 4) 
		    {
		      __m256d leftv = _mm256_load_pd(&left[k * 400 + l * 20 + j]);
		      __m256d rightv = _mm256_load_pd(&right[k * 400 + l * 20 + j]);
		      __m256d vlv = _mm256_load_pd(&vl[j]);
		      __m256d vrv = _mm256_load_pd(&vr[j]);
#ifdef _FMA
		      al = FMAMACC(al,vlv,leftv);
		      ar = FMAMACC(ar,vrv,rightv);
#else
		      al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		      ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif
		    }
		  al = hadd3(al);
		  ar = hadd3(ar);
		  al = _mm256_mul_pd(ar,al);
		  for(j = 0; j < 20; j += 4) 
		    {
		      __m256d vv = _mm256_load_pd(&v[j]);
		      __m256d extEvv = _mm256_load_pd(&extEV[20 * l + j]);
#ifdef _FMA
		      vv = FMAMACC(vv,al,extEvv);
#else
		      vv = _mm256_add_pd(vv,_mm256_mul_pd(al,extEvv));
#endif
		      _mm256_store_pd(&v[j],vv);
		    } 
		} 
	    }
	  v = &(x3[80 * i]);
	  scale = 1;
	  __m256d minlikelihood_avx = _mm256_set1_pd(minlikelihood);
	  for(l = 0; scale && (l < 80); l += 4) 
	    {
	      __m256d vv = _mm256_load_pd(&v[l]);
	      __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
	      vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
	      if(_mm256_movemask_pd(vv_abs) != 15)
		scale = 0;
	    }
	  if(scale) 
	    {	
	      __m256d twotothe256v = _mm256_set_pd(twotothe256,twotothe256,twotothe256,twotothe256);
	      for(l = 0; l < 80; l += 4) 
		{
		  __m256d vv = _mm256_load_pd(&v[l]);
		  _mm256_store_pd(&v[l],_mm256_mul_pd(vv,twotothe256v));
		}
	      if(useFastScaling)
		addScale += wgt[i];					
	      else
		ex3[i] += 1;
	    } 
	}
      break;
    default:
      assert(0);
    }
  /* as above, increment the global counter that counts scaling multiplications by the scaling multiplications carried out for computing the likelihood array at node p */
  if(useFastScaling)
    *scalerIncrement = addScale;
}



 

void newviewGTRGAMMAPROT_AVX(int tipCase,
			     double *x1, double *x2, double *x3, double *extEV, double *tipVector,
			     int *ex3, unsigned char *tipX1, unsigned char *tipX2, int n, 
			     double *left, double *right, int *wgt, int *scalerIncrement, const boolean useFastScaling) 
{
  double	
    *uX1, 
    *uX2, 
    *v, 
    x1px2, 
    *vl, 
    *vr;
  
  int	
    i, 
    j, 
    l, 
    k, 
    scale, 
    addScale = 0;

 
#ifndef GCC_VERSION
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#endif


#if GCC_VERSION < 40500
   __m256d
    bitmask = _mm256_set_pd(0,0,0,-1);
#else
  __m256i
    bitmask = _mm256_set_epi32(0, 0, 0, 0, 0, 0, -1, -1);
#endif 
  
  switch(tipCase) 
    {
    case TIP_TIP: 
      {
       
	double 
	  umpX1[1840] __attribute__ ((aligned (BYTE_ALIGNMENT))), 
	  umpX2[1840] __attribute__ ((aligned (BYTE_ALIGNMENT)));

	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);
	    
	    for(k = 0; k < 80; k++) 
	      {
		double 
		  *ll =  &left[k * 20],
		  *rr =  &right[k * 20];
		
		__m256d 
		  umpX1v = _mm256_setzero_pd(),
		  umpX2v = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
#ifdef _FMA
		    __m256d llv = _mm256_load_pd(&ll[l]);
		    umpX1v = FMAMACC(umpX1v,vv,llv);
		    __m256d rrv = _mm256_load_pd(&rr[l]);
		    umpX2v = FMAMACC(umpX2v,vv,rrv);
#else		    
		    umpX1v = _mm256_add_pd(umpX1v,_mm256_mul_pd(vv,_mm256_load_pd(&ll[l])));
		    umpX2v = _mm256_add_pd(umpX2v,_mm256_mul_pd(vv,_mm256_load_pd(&rr[l])));
#endif
		  }
		
		umpX1v = hadd3(umpX1v);
		umpX2v = hadd3(umpX2v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
		_mm256_maskstore_pd(&umpX2[80 * i + k], bitmask, umpX2v);
	      } 
	  }

	for(i = 0; i < n; i++) 
	  {	    
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];
	   
	    for(j = 0; j < 4; j++) 
	      {     	
		__m256d vv[5];  

		v = &x3[i * 80 + j * 20];
			
		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();

		for(k = 0; k < 20; k++) 
		  {			 
		    x1px2 = uX1[j * 20 + k] * uX2[j * 20 + k];

		    __m256d x1px2v = _mm256_set1_pd(x1px2);		    
		    
		    __m256d extEvv = _mm256_load_pd(&extEV[20 * k]);
#ifdef _FMA
		    vv[0] = FMAMACC(vv[0],x1px2v,extEvv);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[0],vv[0]);
		    
		    extEvv = _mm256_load_pd(&extEV[20 * k + 4]);
#ifdef _FMA
		    vv[1] = FMAMACC(vv[1],x1px2v,extEvv);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

		    extEvv = _mm256_load_pd(&extEV[20 * k + 8]);
#ifdef _FMA
		    vv[2] = FMAMACC(vv[2],x1px2v,extEvv);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[8],vv[2]);

		    extEvv = _mm256_load_pd(&extEV[20 * k + 12]);
#ifdef _FMA
		    vv[3] = FMAMACC(vv[3],x1px2v,extEvv);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[12],vv[3]);

		    extEvv = _mm256_load_pd(&extEV[20 * k + 16]);
#ifdef _FMA
		    vv[4] = FMAMACC(vv[4],x1px2v,extEvv);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v,extEvv));
#endif
		    _mm256_store_pd(&v[16],vv[4]);
		  } 
	      } 
	  } 
      } 
      break;
    case TIP_INNER: 
      {

	double 
	  umpX1[1840] __attribute__ ((aligned (BYTE_ALIGNMENT))),
	  ump_x2[20] __attribute__ ((aligned (BYTE_ALIGNMENT)));

	for(i = 0; i < 23; i++) 
	  {
	    v = &(tipVector[20 * i]);

	    for(k = 0; k < 80; k++) 
	      {
		__m256d umpX1v = _mm256_setzero_pd();
		for(l = 0; l < 20; l+=4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    __m256d leftv = _mm256_load_pd(&left[k * 20 + l]);
#ifdef _FMA
		   
		    umpX1v = FMAMACC(umpX1v, vv, leftv);
#else
		    umpX1v = _mm256_add_pd(umpX1v, _mm256_mul_pd(vv, leftv));
#endif
		  }
		umpX1v = hadd3(umpX1v);
		_mm256_maskstore_pd(&umpX1[80 * i + k], bitmask, umpX1v);
	      } 
	  }
	
	for (i = 0; i < n; i++) 
	  {	   
	    uX1 = &umpX1[80 * tipX1[i]];
	   	    
	    for(k = 0; k < 4; k++) 
	      {
		v = &(x2[80 * i + k * 20]);
		
		for(l = 0; l < 20; l++) 
		  {
		    __m256d ump_x2v = _mm256_setzero_pd();
		    		  
		    __m256d vv = _mm256_load_pd(&v[0]);
		    __m256d rightv = _mm256_load_pd(&right[k*400+l*20+0]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    vv = _mm256_load_pd(&v[4]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+4]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[8]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+8]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[12]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+12]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif

		    vv = _mm256_load_pd(&v[16]);
		    rightv = _mm256_load_pd(&right[k*400+l*20+16]);
#ifdef _FMA
		    ump_x2v = FMAMACC(ump_x2v,vv,rightv);
#else
		    ump_x2v = _mm256_add_pd(ump_x2v, _mm256_mul_pd(vv, rightv));
#endif
		    
		    ump_x2v = hadd3(ump_x2v);
		    _mm256_maskstore_pd(&ump_x2[l], bitmask, ump_x2v);
		  }
		
		v = &(x3[80 * i + 20 * k]);
	

		__m256d vv[5]; 

		vv[0] = _mm256_setzero_pd();
		vv[1] = _mm256_setzero_pd();
		vv[2] = _mm256_setzero_pd();
		vv[3] = _mm256_setzero_pd();
		vv[4] = _mm256_setzero_pd();
		
		for(l = 0; l < 20; l++) 
		  {
		    x1px2 = uX1[k * 20 + l]	* ump_x2[l];
		    __m256d x1px2v = _mm256_set1_pd(x1px2);	
	    		 
#ifdef _FMA
		    __m256d ev = _mm256_load_pd(&extEV[l * 20 + 0]);
		    vv[0] = FMAMACC(vv[0],x1px2v, ev);
#else
		    vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 0])));
#endif
		    _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 4]);
		    vv[1] = FMAMACC(vv[1],x1px2v, ev);
#else
		    vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 4])));
#endif
		    _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 8]);
		    vv[2] = FMAMACC(vv[2],x1px2v, ev);
#else
		    vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 8])));
#endif
		    _mm256_store_pd(&v[8],vv[2]);
		    
#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 12]);
		    vv[3] = FMAMACC(vv[3],x1px2v, ev);
#else
		    vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 12])));
#endif
		    _mm256_store_pd(&v[12],vv[3]);


#ifdef _FMA
		    ev = _mm256_load_pd(&extEV[l * 20 + 16]);
		    vv[4] = FMAMACC(vv[4],x1px2v, ev);
#else
		    vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(x1px2v, _mm256_load_pd(&extEV[l * 20 + 16])));
#endif
		    _mm256_store_pd(&v[16],vv[4]);

		  } 
	      }
	   
	    v = &x3[80 * i];
	    __m256d minlikelihood_avx = _mm256_set1_pd(minlikelihood);
	    scale = 1;
	    for(l = 0; scale && (l < 80); l += 4) 
	      {
		__m256d vv = _mm256_load_pd(&v[l]);
		__m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
		vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
		if(_mm256_movemask_pd(vv_abs) != 15)
		  scale = 0;
	      }
	    
	    if(scale) 
	      {		
		__m256d twotothe256v = _mm256_set_pd(twotothe256,twotothe256,twotothe256,twotothe256);
		for(l = 0; l < 80; l += 4) 
		  {
		    __m256d vv = _mm256_load_pd(&v[l]);
		    _mm256_store_pd(&v[l],_mm256_mul_pd(vv,twotothe256v));
		  }
		if(useFastScaling)
		  addScale += wgt[i];				
		else
		  ex3[i] += 1;
	      } 
	  } 
      } 
      break;
    case INNER_INNER:      
      for(i = 0; i < n; i++) 
	{ 
	  scale = 1;
	  
	  for(k = 0; k < 4; k++) 
	    {
	      vl = &(x1[80 * i + 20 * k]);
	      vr = &(x2[80 * i + 20 * k]);
	      v  = &(x3[80 * i + 20 * k]);	      	   

	      __m256d vv[5]; 
	      
	      vv[0] = _mm256_setzero_pd();
	      vv[1] = _mm256_setzero_pd();
	      vv[2] = _mm256_setzero_pd();
	      vv[3] = _mm256_setzero_pd();
	      vv[4] = _mm256_setzero_pd();
	      
	      for(l = 0; l < 20; l++) 
		{		  
		  __m256d al = _mm256_setzero_pd();
		  __m256d ar = _mm256_setzero_pd();
       		  
		  __m256d leftv  = _mm256_load_pd(&left[k * 400 + l * 20 + 0]);
		  __m256d rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 0]);
		  __m256d vlv = _mm256_load_pd(&vl[0]);
		  __m256d vrv = _mm256_load_pd(&vr[0]);
		  
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));		  
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 4]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 4]);
		  vlv = _mm256_load_pd(&vl[4]);
		  vrv = _mm256_load_pd(&vr[4]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 8]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 8]);
		  vlv = _mm256_load_pd(&vl[8]);
		  vrv = _mm256_load_pd(&vr[8]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 12]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 12]);
		  vlv = _mm256_load_pd(&vl[12]);
		  vrv = _mm256_load_pd(&vr[12]);
#ifdef _FMA
		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  leftv = _mm256_load_pd(&left[k * 400 + l * 20 + 16]);
		  rightv = _mm256_load_pd(&right[k * 400 + l * 20 + 16]);
		  vlv = _mm256_load_pd(&vl[16]);
		  vrv = _mm256_load_pd(&vr[16]);

#ifdef _FMA		    
		  al = FMAMACC(al, vlv, leftv);
		  ar = FMAMACC(ar, vrv, rightv);
#else
		  al = _mm256_add_pd(al,_mm256_mul_pd(vlv,leftv));
		  ar = _mm256_add_pd(ar,_mm256_mul_pd(vrv,rightv));
#endif

		  /**************************************************************************************************************/

		  al = hadd3(al);
		  ar = hadd3(ar);
		  al = _mm256_mul_pd(ar,al);
		  
		  /************************************************************************************************************/
#ifdef _FMA		    
		  __m256d ev =  _mm256_load_pd(&extEV[20 * l + 0]);
		  vv[0] = FMAMACC(vv[0], al, ev);		 
#else
		  vv[0] = _mm256_add_pd(vv[0],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 0])));			  		 		  
#endif
		  _mm256_store_pd(&v[0],vv[0]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 4]);
		  vv[1] = FMAMACC(vv[1], al, ev);		 
#else
		  vv[1] = _mm256_add_pd(vv[1],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 4])));		  		 
#endif
		  _mm256_store_pd(&v[4],vv[1]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 8]);
		  vv[2] = FMAMACC(vv[2], al, ev);		 
#else
		  vv[2] = _mm256_add_pd(vv[2],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 8])));		  		 
#endif
		  _mm256_store_pd(&v[8],vv[2]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 12]);
		  vv[3] = FMAMACC(vv[3], al, ev);		 
#else
		  vv[3] = _mm256_add_pd(vv[3],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 12])));		  		 
#endif
		  _mm256_store_pd(&v[12],vv[3]);

#ifdef _FMA		    
		  ev =  _mm256_load_pd(&extEV[20 * l + 16]);
		  vv[4] = FMAMACC(vv[4], al, ev);		 
#else
		  vv[4] = _mm256_add_pd(vv[4],_mm256_mul_pd(al, _mm256_load_pd(&extEV[20 * l + 16])));			 	  
#endif
		  _mm256_store_pd(&v[16],vv[4]);		 
		} 
	    }
	  v = &(x3[80 * i]);
	  scale = 1;
	  __m256d minlikelihood_avx = _mm256_set1_pd(minlikelihood);	 

	  for(l = 0; scale && (l < 80); l += 4) 
	    {
	      __m256d vv = _mm256_load_pd(&v[l]);
	      __m256d vv_abs = _mm256_and_pd(vv,absMask_AVX.m);
	      vv_abs = _mm256_cmp_pd(vv_abs,minlikelihood_avx,_CMP_LT_OS);
	      if(_mm256_movemask_pd(vv_abs) != 15)
		scale = 0;	     
	    }

	  if(scale) 
	    {		     	      
	      __m256d twotothe256v = _mm256_set_pd(twotothe256,twotothe256,twotothe256,twotothe256);
	      for(l = 0; l < 80; l += 4) 
		{
		  __m256d vv = _mm256_load_pd(&v[l]);
		  _mm256_store_pd(&v[l],_mm256_mul_pd(vv,twotothe256v));
		}
	      if(useFastScaling)
		addScale += wgt[i];					
	      else
		ex3[i] += 1;
	    } 
	}
      break;
    default:
      assert(0);
    }
 
  if(useFastScaling)
    *scalerIncrement = addScale;
}
