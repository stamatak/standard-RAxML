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
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdint.h>

#ifdef __AVX

#ifdef __SIM_SSE3

#define _SSE3_WAS_DEFINED

#undef __SIM_SSE3

#endif

#endif


#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#endif

#ifdef __AVX

#include <xmmintrin.h>
#include <immintrin.h>

#endif


#include "axml.h"



extern const unsigned int mask32[32]; 
/* vector-specific stuff */

extern char **globalArgv;
extern int globalArgc;

#ifdef __SIM_SSE3

#define INTS_PER_VECTOR 4
#define INT_TYPE __m128i
#define CAST __m128i*
#define SET_ALL_BITS_ONE _mm_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO _mm_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm_load_si128
#define VECTOR_BIT_AND _mm_and_si128
#define VECTOR_BIT_OR  _mm_or_si128
#define VECTOR_STORE  _mm_store_si128
#define VECTOR_AND_NOT _mm_andnot_si128

#endif

#ifdef __AVX

#define INTS_PER_VECTOR 8
#define INT_TYPE __m256d
#define CAST double*
#define SET_ALL_BITS_ONE (__m256d)_mm256_set_epi32(0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF)
#define SET_ALL_BITS_ZERO (__m256d)_mm256_set_epi32(0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000, 0x00000000)
#define VECTOR_LOAD _mm256_load_pd
#define VECTOR_BIT_AND _mm256_and_pd
#define VECTOR_BIT_OR  _mm256_or_pd
#define VECTOR_STORE  _mm256_store_pd
#define VECTOR_AND_NOT _mm256_andnot_pd

#endif

extern double masterTime;
extern char  workdir[1024];
extern char run_id[128];

/********************************DNA FUNCTIONS *****************************************************************/



static void checkSeed(analdef *adef)
{ 
  static boolean seedChecked = FALSE;

  if(!seedChecked) 
    {
      /*printf("Checking seed\n");*/

      if(adef->parsimonySeed <= 0)
	{
	  printf("Error: you need to specify a random number seed with \"-p\" for the randomized stepwise addition\n");
	  printf("parsimony algorithm or random tree building algorithm such that runs can be reproduced and debugged ... exiting\n");      
	}
  
      assert(adef->parsimonySeed > 0);
      seedChecked = TRUE;
    }
}

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }
}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips, boolean full)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->x)
    getxnodeLocal(p);  
  
  if(full)
    {
       if(q->number > maxTips) 
	 computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  else
    {
      if(q->number > maxTips && !q->x) 
	computeTraversalInfoParsimony(q, ti, counter, maxTips, full);
      
      if(r->number > maxTips && !r->x) 
	computeTraversalInfoParsimony(r, ti, counter, maxTips, full);
    }
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}



#if (defined(__SIM_SSE3) || defined(__AVX))
#define BIT_COUNT(x)  precomputed16_bitcount(x)

/* 
   The critical speed of this function
   is mostly a machine/architecure issue
   Nontheless, I decided not to use 
   __builtin_popcount(x)
   for better general portability, albeit it's worth testing
   on x86 64 bit architectures
*/
#else
#define BIT_COUNT(x)  precomputed16_bitcount(x)
#endif



#if (defined(__SIM_SSE3) || defined(__AVX))

#ifdef __SIM_SSE3

static unsigned int vectorCount(__m128i b)
{
  const unsigned int 
    mu1 = 0x55555555,
    mu2 = 0x33333333,
    mu3 = 0x0F0F0F0F,
    mu4 = 0x0000003F;
  
  unsigned int 
    tcnt[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
  
  __m128i 
    m1 = _mm_set_epi32 (mu1, mu1, mu1, mu1),
    m2 = _mm_set_epi32 (mu2, mu2, mu2, mu2),
    m3 = _mm_set_epi32 (mu3, mu3, mu3, mu3),
    m4 = _mm_set_epi32 (mu4, mu4, mu4, mu4),
    tmp1, 
    tmp2;
 

  /* b = (b & 0x55555555) + (b >> 1 & 0x55555555); */
  tmp1 = _mm_srli_epi32(b, 1);                    /* tmp1 = (b >> 1 & 0x55555555)*/
  tmp1 = _mm_and_si128(tmp1, m1); 
  tmp2 = _mm_and_si128(b, m1);                    /* tmp2 = (b & 0x55555555) */
  b    = _mm_add_epi32(tmp1, tmp2);               /*  b = tmp1 + tmp2 */

  /* b = (b & 0x33333333) + (b >> 2 & 0x33333333); */
  tmp1 = _mm_srli_epi32(b, 2);                    /* (b >> 2 & 0x33333333) */
  tmp1 = _mm_and_si128(tmp1, m2); 
  tmp2 = _mm_and_si128(b, m2);                    /* (b & 0x33333333) */
  b    = _mm_add_epi32(tmp1, tmp2);               /* b = tmp1 + tmp2 */

  /* b = (b + (b >> 4)) & 0x0F0F0F0F; */
  tmp1 = _mm_srli_epi32(b, 4);                    /* tmp1 = b >> 4 */
  b = _mm_add_epi32(b, tmp1);                     /* b = b + (b >> 4) */
  b = _mm_and_si128(b, m3);                       /*           & 0x0F0F0F0F */

  /* b = b + (b >> 8); */
  tmp1 = _mm_srli_epi32 (b, 8);                   /* tmp1 = b >> 8 */
  b = _mm_add_epi32(b, tmp1);                     /* b = b + (b >> 8) */
  
  /* b = (b + (b >> 16)) & 0x0000003F; */
  tmp1 = _mm_srli_epi32 (b, 16);                  /* b >> 16 */
  b = _mm_add_epi32(b, tmp1);                     /* b + (b >> 16) */
  b = _mm_and_si128(b, m4);                       /* (b >> 16) & 0x0000003F; */
   
  _mm_store_si128((__m128i *)tcnt, b);

  return tcnt[0] + tcnt[1] + tcnt[2] + tcnt[3];
}

#endif

static inline unsigned int populationCount(INT_TYPE v_N)
{
#ifdef __AVX
  {
    unsigned long int
      res[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
    unsigned int a, b;
    
    _mm256_store_pd((double*)res, v_N);
    
    a = __builtin_popcountl(res[0]) + __builtin_popcountl(res[1]);
    b = __builtin_popcountl(res[2]) + __builtin_popcountl(res[3]);
    
    return (a + b);	   
  }
#else	  
  return (vectorCount(v_N)); 
#endif
}

void newviewParsimonyIterativeFast(tree *tr)
{    
  INT_TYPE
    allOne = SET_ALL_BITS_ONE;

  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *this[2];

		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C,
		      v_A, v_C;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);		  		  		  		  
		    
		    v_N = VECTOR_BIT_OR(l_A, l_C);
		    
		    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);		  
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *this[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    INT_TYPE
		      s_r, s_l, v_N,
		      l_A, l_C, l_G, l_T,
		      v_A, v_C, v_G, v_T;	    	 
		    
		    s_l = VECTOR_LOAD((CAST)(&left[0][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[0][i]));
		    l_A = VECTOR_BIT_AND(s_l, s_r);
		    v_A = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[1][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[1][i]));
		    l_C = VECTOR_BIT_AND(s_l, s_r);
		    v_C = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[2][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[2][i]));
		    l_G = VECTOR_BIT_AND(s_l, s_r);
		    v_G = VECTOR_BIT_OR(s_l, s_r);
		    
		    s_l = VECTOR_LOAD((CAST)(&left[3][i]));
		    s_r = VECTOR_LOAD((CAST)(&right[3][i]));
		    l_T = VECTOR_BIT_AND(s_l, s_r);
		    v_T = VECTOR_BIT_OR(s_l, s_r);
		    
		    v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));	  	 	    	  
		    
		    VECTOR_STORE((CAST)(&this[0][i]), VECTOR_BIT_OR(l_A, VECTOR_AND_NOT(v_N, v_A)));
		    VECTOR_STORE((CAST)(&this[1][i]), VECTOR_BIT_OR(l_C, VECTOR_AND_NOT(v_N, v_C)));
		    VECTOR_STORE((CAST)(&this[2][i]), VECTOR_BIT_OR(l_G, VECTOR_AND_NOT(v_N, v_G)));
		    VECTOR_STORE((CAST)(&this[3][i]), VECTOR_BIT_OR(l_T, VECTOR_AND_NOT(v_N, v_T)));	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);	
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *this[20];

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[20], 
		      v_A[20];	    	 
		    
		    for(j = 0; j < 20; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < 20; j++)		    
		      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);
		  }
	      }
	      break;
	    default:
	      {
		parsimonyNumber
		  *left[32], 
		  *right[32],
		  *this[32];

		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i += INTS_PER_VECTOR)
		  {	 	  
		    size_t j;
		    
		    INT_TYPE
		      s_r, s_l, 
		      v_N = SET_ALL_BITS_ZERO,
		      l_A[32], 
		      v_A[32];	    	 
		    
		    for(j = 0; j < states; j++)
		      {
			s_l = VECTOR_LOAD((CAST)(&left[j][i]));
			s_r = VECTOR_LOAD((CAST)(&right[j][i]));
			l_A[j] = VECTOR_BIT_AND(s_l, s_r);
			v_A[j] = VECTOR_BIT_OR(s_l, s_r);
			
			v_N = VECTOR_BIT_OR(v_N, l_A[j]);
		      }
		    
		    for(j = 0; j < states; j++)		    
		      VECTOR_STORE((CAST)(&this[j][i]), VECTOR_BIT_OR(l_A[j], VECTOR_AND_NOT(v_N, v_A[j])));		 	  	 	 	  	  	  	
		    
		    v_N = VECTOR_AND_NOT(v_N, allOne);
		    
		    totalScore += populationCount(v_N);
		  }	  			
	      }
	    }	  	 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}

static inline unsigned int evaluatePopcount(INT_TYPE v_N)
{
#ifdef __AVX            	       	   	      
  unsigned long int
    res[4] __attribute__ ((aligned (BYTE_ALIGNMENT)));
	     
  unsigned int a, b;
	     
  _mm256_store_pd((double*)res, v_N);
  
  a = __builtin_popcountl(res[0]) + __builtin_popcountl(res[1]);
  b = __builtin_popcountl(res[2]) + __builtin_popcountl(res[3]);
	     
  return (a + b);	            
#else      	       
  unsigned int
    sum = 0,
    counts[INTS_PER_VECTOR] __attribute__ ((aligned (BYTE_ALIGNMENT)));

  VECTOR_STORE((CAST)counts, v_N);

  sum += BIT_COUNT(counts[0]) + BIT_COUNT(counts[1]);
  sum += BIT_COUNT(counts[2]) + BIT_COUNT(counts[3]);          

  return sum;
#endif
}

unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  INT_TYPE 
    allOne = SET_ALL_BITS_ONE;

  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                       
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),		 
		   v_N = VECTOR_BIT_OR(l_A, l_C);
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += evaluatePopcount(v_N);
		 
		 if(sum >= bestScore)
		   return sum;		   	       
	       }
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	                        
		 INT_TYPE      
		   l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[0][i])), VECTOR_LOAD((CAST)(&right[0][i]))),
		   l_C = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[1][i])), VECTOR_LOAD((CAST)(&right[1][i]))),
		   l_G = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[2][i])), VECTOR_LOAD((CAST)(&right[2][i]))),
		   l_T = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[3][i])), VECTOR_LOAD((CAST)(&right[3][i]))),
		   v_N = VECTOR_BIT_OR(VECTOR_BIT_OR(l_A, l_C), VECTOR_BIT_OR(l_G, l_T));     
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += evaluatePopcount(v_N);
		 
		 if(sum >= bestScore)		 
		   return sum;	        
	       }	   	 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i += INTS_PER_VECTOR)
		{                	       
		  int 
		    j;
		  
		  INT_TYPE      
		    l_A,
		    v_N = SET_ALL_BITS_ZERO;     
		  
		  for(j = 0; j < 20; j++)
		    {
		      l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		      v_N = VECTOR_BIT_OR(l_A, v_N);
		    }
		  
		  v_N = VECTOR_AND_NOT(v_N, allOne);
		  
		  sum += evaluatePopcount(v_N);	       
		  
		  if(sum >= bestScore)	    
		    return sum;		    	       
		}
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       *left[32],  
	       *right[32]; 

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i += INTS_PER_VECTOR)
	       {                	       
		 size_t
		   j;
		 
		 INT_TYPE      
		   l_A,
		   v_N = SET_ALL_BITS_ZERO;     
		 
		 for(j = 0; j < states; j++)
		   {
		     l_A = VECTOR_BIT_AND(VECTOR_LOAD((CAST)(&left[j][i])), VECTOR_LOAD((CAST)(&right[j][i])));
		     v_N = VECTOR_BIT_OR(l_A, v_N);
		   }
		 
		 v_N = VECTOR_AND_NOT(v_N, allOne);
		 
		 sum += evaluatePopcount(v_N);	       
		 
		 if(sum >= bestScore)	      
		   return sum;		       
	       }
	   }
	 }
    }
  
  return sum;
}


#else

void newviewParsimonyIterativeFast(tree *tr)
{    
  int 
    model,
    *ti = tr->ti,
    count = ti[0],
    index; 

  for(index = 4; index < count; index += 4)
    {      
      unsigned int
	totalScore = 0;

      size_t
	pNumber = (size_t)ti[index],
	qNumber = (size_t)ti[index + 1],
	rNumber = (size_t)ti[index + 2];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  size_t
	    k,
	    states = tr->partitionData[model].states,
	    width = tr->partitionData[model].parsimonyLength;	 
            
	  unsigned int	
	    i;      
                 
	  switch(states)
	    {
	    case 2:       
	      {
		parsimonyNumber
		  *left[2],
		  *right[2],
		  *this[2];
		
		parsimonyNumber
		   o_A,
		   o_C,
		   t_A,
		   t_C,	
		   t_N;
		
		for(k = 0; k < 2; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 2 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];		   

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		  
		    t_N = ~(t_A | t_C);	  

		    this[0][i] = t_A | (t_N & o_A);
		    this[1][i] = t_C | (t_N & o_C);		   
		    
		    totalScore += BIT_COUNT(t_N);   
		  }
	      }
	      break;
	    case 4:
	      {
		parsimonyNumber
		  *left[4],
		  *right[4],
		  *this[4];

		for(k = 0; k < 4; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 4 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
		  }

		parsimonyNumber
		   o_A,
		   o_C,
		   o_G,
		   o_T,
		   t_A,
		   t_C,
		   t_G,
		   t_T,	
		   t_N;

		for(i = 0; i < width; i++)
		  {	 	  
		    t_A = left[0][i] & right[0][i];
		    t_C = left[1][i] & right[1][i];
		    t_G = left[2][i] & right[2][i];	  
		    t_T = left[3][i] & right[3][i];

		    o_A = left[0][i] | right[0][i];
		    o_C = left[1][i] | right[1][i];
		    o_G = left[2][i] | right[2][i];	  
		    o_T = left[3][i] | right[3][i];

		    t_N = ~(t_A | t_C | t_G | t_T);	  

		    this[0][i] = t_A | (t_N & o_A);
		    this[1][i] = t_C | (t_N & o_C);
		    this[2][i] = t_G | (t_N & o_G);
		    this[3][i] = t_T | (t_N & o_T); 
		    
		    totalScore += BIT_COUNT(t_N);   
		  }
	      }
	      break;
	    case 20:
	      {
		parsimonyNumber
		  *left[20],
		  *right[20],
		  *this[20];

		parsimonyNumber
		  o_A[20],
		  t_A[20],	  
		  t_N;

		for(k = 0; k < 20; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * 20 * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		  }

		for(i = 0; i < width; i++)
		  {	 	  
		    size_t k;
		    
		    t_N = 0;

		    for(k = 0; k < 20; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;

		    for(k = 0; k < 20; k++)		      
		      this[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += BIT_COUNT(t_N); 
		  }
	      }
	      break;
	    default:
	      {		
		parsimonyNumber
		  *left[32],
		  *right[32],
		  *this[32];
		
		parsimonyNumber
		  o_A[32],
		  t_A[32],	  
		  t_N;
		
		assert(states <= 32);
		
		for(k = 0; k < states; k++)
		  {
		    left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		    right[k] = &(tr->partitionData[model].parsVect[(width * states * rNumber) + width * k]);
		    this[k]  = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
		  }
		
		for(i = 0; i < width; i++)
		  {	 	  
		    t_N = 0;
		    
		    for(k = 0; k < states; k++)
		      {
			t_A[k] = left[k][i] & right[k][i];
			o_A[k] = left[k][i] | right[k][i];
			t_N = t_N | t_A[k];
		      }
		    
		    t_N = ~t_N;
		    
		    for(k = 0; k < states; k++)		      
		      this[k][i] = t_A[k] | (t_N & o_A[k]);		   
		    
		    totalScore += BIT_COUNT(t_N); 
		  }
	      }			      
	    } 
	}

      tr->parsimonyScore[pNumber] = totalScore + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];      
    }
}



unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  size_t 
    pNumber = (size_t)tr->ti[1],
    qNumber = (size_t)tr->ti[2];

  int
    model;

  unsigned int 
    bestScore = tr->bestParsimony,    
    sum;

  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr); 

  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = tr->partitionData[model].states,
	width = tr->partitionData[model].parsimonyLength, 
	i;

       switch(states)
	 {
	 case 2:
	   {
	     parsimonyNumber 
	       t_A,
	       t_C,	      
	       t_N,
	       *left[2],
	       *right[2];
	     
	     for(k = 0; k < 2; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 2 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 2 * pNumber) + width * k]);
	       }     
	     
	     for(i = 0; i < width; i++)
	       {                	                       
		 t_A = left[0][i] & right[0][i];
		 t_C = left[1][i] & right[1][i];
		 
		  t_N = ~(t_A | t_C);

		  sum += BIT_COUNT(t_N);    
		 
		 if(sum >= bestScore)
		   return sum;		   	       
	       }
	   }
	   break;
	 case 4:
	   {
	     parsimonyNumber
	       t_A,
	       t_C,
	       t_G,
	       t_T,
	       t_N,
	       *left[4],
	       *right[4];
      
	     for(k = 0; k < 4; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * 4 * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * 4 * pNumber) + width * k]);
	       }        

	     for(i = 0; i < width; i++)
	       {                	                        
		  t_A = left[0][i] & right[0][i];
		  t_C = left[1][i] & right[1][i];
		  t_G = left[2][i] & right[2][i];	  
		  t_T = left[3][i] & right[3][i];

		  t_N = ~(t_A | t_C | t_G | t_T);

		  sum += BIT_COUNT(t_N);     
		 
		 if(sum >= bestScore)		 
		   return sum;	        
	       }	   	 
	   }
	   break;
	 case 20:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[20],
	       *right[20];
	     
	      for(k = 0; k < 20; k++)
		{
		  left[k]  = &(tr->partitionData[model].parsVect[(width * 20 * qNumber) + width * k]);
		  right[k] = &(tr->partitionData[model].parsVect[(width * 20 * pNumber) + width * k]);
		}  
	   
	      for(i = 0; i < width; i++)
		{ 
		  t_N = 0;
		  
		  for(k = 0; k < 20; k++)
		    {
		      t_A = left[k][i] & right[k][i];
		      t_N = t_N | t_A;
		    }
  	       
		  t_N = ~t_N;

		  sum += BIT_COUNT(t_N);      
		  
		  if(sum >= bestScore)	    
		    return sum;		    	       
		}
	   }
	   break;
	 default:
	   {
	     parsimonyNumber
	       t_A,
	       t_N,
	       *left[32], 
	       *right[32];  

	     assert(states <= 32);

	     for(k = 0; k < states; k++)
	       {
		 left[k]  = &(tr->partitionData[model].parsVect[(width * states * qNumber) + width * k]);
		 right[k] = &(tr->partitionData[model].parsVect[(width * states * pNumber) + width * k]);
	       }  
	   
	     for(i = 0; i < width; i++)
	       {                	       
		 t_N = 0;
		  
		 for(k = 0; k < states; k++)
		   {
		     t_A = left[k][i] & right[k][i];
		     t_N = t_N | t_A;
		   }
  	       
		  t_N = ~t_N;

		  sum += BIT_COUNT(t_N);      
		  		  		 
		 if(sum >= bestScore)			  
		   return sum;			   
	       }	     	     
	   }
	 }
    }
  
  return sum;
}

#endif






static unsigned int evaluateParsimony(tree *tr, nodeptr p, boolean full)
{
  volatile unsigned int result;
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(full)
    {
      if(p->number > tr->mxtips)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }
  else
    {
      if(p->number > tr->mxtips && !p->x)
	computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips, full);
      if(q->number > tr->mxtips && !q->x)
	computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips, full); 
    }

  ti[0] = counter;

  result = evaluateParsimonyIterativeFast(tr);

  return result;
}


static void newviewParsimony(tree *tr, nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
    tr->ti[0] = counter;            
    
    newviewParsimonyIterativeFast(tr);      
  }
}





/****************************************************************************************************************************************/

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches); 
   
  newviewParsimony(tr, p);     
} 

/*
  static nodeptr buildNewTip (tree *tr, nodeptr p)
  { 
  nodeptr  q;
  
  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
  assert(q == q->next->next->next);
  assert(q->x || q->next->x || q->next->next->x);
  return  q;
  } 
*/

static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  q->next->back = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;
 
  return  q;
} 

static void buildSimpleTree (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, s, p);
}


static void testInsertParsimony (tree *tr, nodeptr p, nodeptr q)
{ 
  unsigned int 
    mp;
 
  nodeptr  
    r = q->back;   

  boolean 
    doIt = TRUE;
    
  if(tr->grouped)
    {
      int 
	rNumber = tr->constraintVector[r->number],
	qNumber = tr->constraintVector[q->number],
	pNumber = tr->constraintVector[p->number];

      doIt = FALSE;
     
      if(pNumber == -9)
	pNumber = checker(tr, p->back);
      if(pNumber == -9)
	doIt = TRUE;
      else
	{
	  if(qNumber == -9)
	    qNumber = checker(tr, q);

	  if(rNumber == -9)
	    rNumber = checker(tr, r);

	  if(pNumber == rNumber || pNumber == qNumber)
	    doIt = TRUE;       
	}
    }

  if(doIt)
    {
      insertParsimony(tr, p, q);   
  
      mp = evaluateParsimony(tr, p->next->next, FALSE);          
      
      if(mp < tr->bestParsimony)
	{
	  tr->bestParsimony = mp;
	  tr->insertNode = q;
	  tr->removeNode = p;
	}
  
      hookupDefault(q, r, tr->numBranches);
      p->next->next->back = p->next->back = (nodeptr) NULL;
    }
       
  return;
} 


static void restoreTreeParsimony(tree *tr, nodeptr p, nodeptr q)
{ 
  nodeptr
    r = q->back;
  
  int counter = 4;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches);
  
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr); 
}


static void addTraverseParsimony (tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav, boolean doAll)
{        
  if (doAll || (--mintrav <= 0))               
    testInsertParsimony(tr, p, q);	                 

  if (((q->number > tr->mxtips)) && ((--maxtrav > 0) || doAll))
    {	      
      addTraverseParsimony(tr, p, q->next->back, mintrav, maxtrav, doAll);	      
      addTraverseParsimony(tr, p, q->next->next->back, mintrav, maxtrav, doAll);              	     
    }
}


static nodeptr findAnyTipFast(nodeptr p, int numsp)
{ 
  return  (p->number <= numsp)? p : findAnyTipFast(p->next->back, numsp);
} 


static void makePermutationFast(int *perm, int n, analdef *adef)
{    
  int  
    i, 
    j, 
    k;

  checkSeed(adef);      

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {      
      double d =  randum(&adef->parsimonySeed);

      k =  (int)((double)(n + 1 - i) * d);
      
      j        = perm[i];

      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}

static nodeptr  removeNodeParsimony (nodeptr p, tree *tr)
{ 
  nodeptr  q, r;         

  q = p->next->back;
  r = p->next->next->back;   
    
  hookupDefault(q, r, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
  
  return  q;
}

static int rearrangeParsimony(tree *tr, nodeptr p, int mintrav, int maxtrav, boolean doAll)  
{   
  nodeptr  
    p1, 
    p2, 
    q, 
    q1, 
    q2;
  
  int      
    mintrav2; 

  boolean 
    doP = TRUE,
    doQ = TRUE;
           
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3; 

  assert(mintrav == 1);

  if(maxtrav < mintrav)
    return 0;

  q = p->back;

  if(tr->constrained)
    {    
      if(! tipHomogeneityChecker(tr, p->back, 0))
	doP = FALSE;
	
      if(! tipHomogeneityChecker(tr, q->back, 0))
	doQ = FALSE;
		        
      if(doQ == FALSE && doP == FALSE)
	return 0;
    }  

  if((p->number > tr->mxtips) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      if ((p1->number > tr->mxtips) || (p2->number > tr->mxtips)) 
	{	  	  
	  removeNodeParsimony(p, tr);	  	 

	  if ((p1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p1->next->back, mintrav, maxtrav, doAll);         
	      addTraverseParsimony(tr, p, p1->next->next->back, mintrav, maxtrav, doAll);          
	    }
	 
	  if ((p2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, p, p2->next->back, mintrav, maxtrav, doAll);
	      addTraverseParsimony(tr, p, p2->next->next->back, mintrav, maxtrav, doAll);          
	    }
	    
	   
	  hookupDefault(p->next,       p1, tr->numBranches); 
	  hookupDefault(p->next->next, p2, tr->numBranches);	   	    	    

	  newviewParsimony(tr, p);
	}
    }  
       
  if ((q->number > tr->mxtips) && (maxtrav > 0) && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;

      if (
	  (
	   (q1->number > tr->mxtips) && 
	   ((q1->next->back->number > tr->mxtips) || (q1->next->next->back->number > tr->mxtips))
	   )
	  ||
	  (
	   (q2->number > tr->mxtips) && 
	   ((q2->next->back->number > tr->mxtips) || (q2->next->next->back->number > tr->mxtips))
	   )
	  )
	{	   

	  removeNodeParsimony(q, tr);
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if ((q1->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q1->next->back, mintrav2 , maxtrav, doAll);
	      addTraverseParsimony(tr, q, q1->next->next->back, mintrav2 , maxtrav, doAll);         
	    }
	 
	  if ((q2->number > tr->mxtips)) 
	    {
	      addTraverseParsimony(tr, q, q2->next->back, mintrav2 , maxtrav, doAll);
	      addTraverseParsimony(tr, q, q2->next->next->back, mintrav2 , maxtrav, doAll);          
	    }	   
	   
	  hookupDefault(q->next,       q1, tr->numBranches); 
	  hookupDefault(q->next->next, q2, tr->numBranches);
	   
	  newviewParsimony(tr, q);
	}
    }

  return 1;
} 


static void restoreTreeRearrangeParsimony(tree *tr)
{    
  removeNodeParsimony(tr->removeNode, tr);  
  restoreTreeParsimony(tr, tr->removeNode, tr->insertNode);  
}

/*
static boolean isInformative2(tree *tr, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = 15;

  unsigned char
    nucleotide,
    target = 0;
  	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;      	           
    }
  
  
  if(check[1] > 1)
    {
      informativeCounter++;    
      target = target | 1;
    }
  if(check[2] > 1)
    {
      informativeCounter++; 
      target = target | 2;
    }
  if(check[4] > 1)
    {
      informativeCounter++; 
      target = target | 4;
    }
  if(check[8] > 1)
    {
      informativeCounter++; 
      target = target | 8;
    }
	  
  if(informativeCounter >= 2)
    return TRUE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(j == 3 || j == 5 || j == 6 || j == 7 || j == 9 || j == 10 || j == 11 || 
	     j == 12 || j == 13 || j == 14)
	    {
	      if(check[j] > 1)
		{
		  if(!(target & j))
		    return TRUE;
		}
	    }
	} 
    }
     
  return FALSE;	     
}
*/

static boolean isInformative(tree *tr, int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

  const unsigned int
    *bitVector = getBitVector(dataType);

  unsigned char
    nucleotide;
  
	
  for(j = 0; j < 256; j++)
    check[j] = 0;
  
  for(j = 1; j <= tr->mxtips; j++)
    {	   
      nucleotide = tr->yVector[j][site];	    
      check[nucleotide] =  check[nucleotide] + 1;
      assert(bitVector[nucleotide] > 0);	           
    }
  
  for(j = 0; j < undetermined; j++)
    {
      if(check[j] > 0)
	informativeCounter++;    
    } 
	  
  if(informativeCounter <= 1)
    return FALSE;    
  else
    {        
      for(j = 0; j < undetermined; j++)
	{
	  if(check[j] > 1)
	    return TRUE;
	} 
    }
     
  return FALSE;	     
}


static void determineUninformativeSites(tree *tr, int *informative)
{
  int 
    i,
    number = 0;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.

     Alternative definition: informative columns contain at least two types
     of nucleotides, and each nucleotide must appear at least twice in each 
     column. Kind of a pain if we intend to check for this when using, e.g.,
     amibiguous DNA encoding.
  */

  for(i = 0; i < tr->cdta->endsite; i++)
    {
      if(isInformative(tr, tr->dataVector[i], i))
	informative[i] = 1;
      else
	{
	  informative[i] = 0;
	  number++;
	}            
    }
  
 
  /* printf("Uninformative Patterns: %d\n", number); */
}


static void reorderNodes(tree *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if((p->number <= tr->mxtips))    
    return;
  else
    {              
      for(i = tr->mxtips + 1; (i <= (tr->mxtips + tr->mxtips - 1)) && (found == 0); i++)
	{
	  if (p == np[i] || p == np[i]->next || p == np[i]->next->next)
	    {
	      if(p == np[i])			       
		tr->nodep[*count + tr->mxtips + 1] = np[i];		 		
	      else
		{
		  if(p == np[i]->next)		  
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next;		     	   
		  else		   
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next->next;		    		    
		}

	      found = 1;	      	     
	      *count = *count + 1;
	    }
	}            
     
      assert(found != 0);

      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}





  
static void compressDNA(tree *tr, int *informative, boolean saveMemory)
{
  size_t
    totalNodes,
    i,
    model;
  
  if(saveMemory)
    totalNodes = (size_t)tr->innerNodes + 1 + (size_t)tr->mxtips;
  else
    totalNodes = 2 * (size_t)tr->mxtips;

 

  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    {
      size_t
	k,
	states = (size_t)tr->partitionData[model].states,       
	compressedEntries,
	compressedEntriesPadded,
	entries = 0, 
	lower = tr->partitionData[model].lower,
	upper = tr->partitionData[model].upper;

      parsimonyNumber 
	**compressedTips = (parsimonyNumber **)malloc(states * sizeof(parsimonyNumber*)),
	*compressedValues = (parsimonyNumber *)malloc(states * sizeof(parsimonyNumber));
      
      for(i = lower; i < upper; i++)    
	if(informative[i])
	  entries += (size_t)tr->cdta->aliaswgt[i];     
  
      compressedEntries = entries / PCF;

      if(entries % PCF != 0)
	compressedEntries++;

#if (defined(__SIM_SSE3) || defined(__AVX))
      if(compressedEntries % INTS_PER_VECTOR != 0)
	compressedEntriesPadded = compressedEntries + (INTS_PER_VECTOR - (compressedEntries % INTS_PER_VECTOR));
      else
	compressedEntriesPadded = compressedEntries;
#else
      compressedEntriesPadded = compressedEntries;
#endif     

      
      tr->partitionData[model].parsVect = (parsimonyNumber *)malloc_aligned((size_t)compressedEntriesPadded * states * totalNodes * sizeof(parsimonyNumber));
     
      for(i = 0; i < compressedEntriesPadded * states * totalNodes; i++)      
	tr->partitionData[model].parsVect[i] = 0;          

      for(i = 0; i < (size_t)tr->mxtips; i++)
	{
	  size_t
	    w = 0,
	    compressedIndex = 0,
	    compressedCounter = 0,
	    index = 0;

	  for(k = 0; k < states; k++)
	    {
	      compressedTips[k] = &(tr->partitionData[model].parsVect[(compressedEntriesPadded * states * (i + 1)) + (compressedEntriesPadded * k)]);
	      compressedValues[k] = 0;
	    }                
	      
	  for(index = lower; index < (size_t)upper; index++)
	    {
	      if(informative[index])
		{
		  const unsigned int 
		    *bitValue = getBitVector(tr->partitionData[model].dataType);

		  parsimonyNumber 
		    value = bitValue[tr->yVector[i + 1][index]];	  
	      
		  for(w = 0; w < (size_t)tr->cdta->aliaswgt[index]; w++)
		    {	   
		      for(k = 0; k < states; k++)
			{
			  if(value & mask32[k])
			    compressedValues[k] |= mask32[compressedCounter];
			}
		     
		      compressedCounter++;
		  
		      if(compressedCounter == PCF)
			{
			  for(k = 0; k < states; k++)
			    {
			      compressedTips[k][compressedIndex] = compressedValues[k];
			      compressedValues[k] = 0;
			    }			 
			  
			  compressedCounter = 0;
			  compressedIndex++;
			}
		    }
		}
	    }
                           
	  for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	    {	
	      for(;compressedCounter < PCF; compressedCounter++)	      
		for(k = 0; k < states; k++)
		  compressedValues[k] |= mask32[compressedCounter];		  
	  
	      for(k = 0; k < states; k++)
		{
		  compressedTips[k][compressedIndex] = compressedValues[k];
		  compressedValues[k] = 0;
		}	      	      
	      
	      compressedCounter = 0;
	    }	 	
	}               
  
      tr->partitionData[model].parsimonyLength = compressedEntriesPadded;   

      free(compressedTips);
      free(compressedValues);
    }
  
  tr->parsimonyScore = (unsigned int*)malloc_aligned(sizeof(unsigned int) * totalNodes);  
          
  for(i = 0; i < totalNodes; i++) 
    tr->parsimonyScore[i] = 0;
}



static void stepwiseAddition(tree *tr, nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;

  unsigned int 
    mp;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips, FALSE);              
  tr->ti[0] = counter;
  tr->ti[1] = p->number;
  tr->ti[2] = p->back->number;
    
  mp = evaluateParsimonyIterativeFast(tr);
  
  if(mp < tr->bestParsimony)
    {    
      tr->bestParsimony = mp;
      tr->insertNode = q;     
    }
 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {	      
      stepwiseAddition(tr, p, q->next->back);	      
      stepwiseAddition(tr, p, q->next->next->back);              	     
    }
}

static void markNodesInTree(nodeptr p, tree *tr, unsigned char *nodesInTree)
{
  if(isTip(p->number, tr->mxtips))
    nodesInTree[p->number] = 1;
  else
    {
      markNodesInTree(p->next->back, tr, nodesInTree);
      markNodesInTree(p->next->next->back, tr, nodesInTree);
    }

}

void makeParsimonyTreeFast(tree *tr, analdef *adef, boolean full)
{   
  nodeptr  
    p, 
    f;    

  size_t
    model;

  int 
    i, 
    nextsp,
    *perm        = (int *)malloc((size_t)(tr->mxtips + 1) * sizeof(int)),
    *informative = (int *)malloc(sizeof(int) * (size_t)tr->cdta->endsite);  

  unsigned int 
    randomMP, 
    startMP;        

  /* double t; */

  determineUninformativeSites(tr, informative);     

  compressDNA(tr, informative, FALSE);

  free(informative); 

  tr->ti = (int*)malloc(sizeof(int) * 4 * (size_t)tr->mxtips);  
 
  /*t = gettime();*/

  if(!full)
    {           
      unsigned int 
	score;
	
      int 
	j = 0;

      unsigned char
	*nodesInTree = (unsigned char*)calloc((size_t)(tr->mxtips + 1), sizeof(unsigned char));	      
	
      tr->start = findAnyTipFast(tr->start, tr->rdta->numsp);
	
      tr->bestParsimony = INT_MAX;

      score = evaluateParsimony(tr, tr->start->back, TRUE);
		
      assert(tr->start);
      
      checkSeed(adef);

      markNodesInTree(tr->start, tr, nodesInTree);
      markNodesInTree(tr->start->back, tr, nodesInTree);
	
      j = tr->ntips + 1;
	
      if(tr->grouped)
	{
	  for(i = 1; i <= tr->mxtips; i++)      
	    {
	      if(tr->constraintVector[i] == -1) 
		{
		  perm[j++] = i;		
		  tr->constraintVector[i] = -9;
		}
	    }
	}
      else
	{
	  if(tr->constrained)
	    { 
	      for(i = 1; i <= tr->mxtips; i++)
		tr->constraintVector[i] = 0;
	      
	      for(i = 1; i <= tr->mxtips; i++)
		{		  
		  if(nodesInTree[i] == 0) 	      
		    perm[j++] = i;
		  else
		    tr->constraintVector[i] = 1;		    
		}
	    }
	  else
	    {
	      for(i = 1; i <= tr->mxtips; i++)      
		if(nodesInTree[i] == 0) 
		  perm[j++] = i;	  
	    }
	}
	
      for(i = tr->ntips + 1; i <= tr->mxtips; i++) 
	{	     
	  int k, j;
	  
	  k =  (int)((double)(tr->mxtips + 1 - i) * randum(&adef->parsimonySeed));
	  
	  assert(i + k <= tr->mxtips);
	  j        = perm[i];
	  perm[i]     = perm[i + k];
	  perm[i + k] = j;
	}    
	
      f = tr->start;     

      free(nodesInTree);
    }
  else
    {
      assert(!tr->constrained);

      makePermutationFast(perm, tr->mxtips, adef);
      
      tr->ntips = 0;    
      
      tr->nextnode = tr->mxtips + 1;       
      
      buildSimpleTree(tr, perm[1], perm[2], perm[3]);      
      
      f = tr->start;
    }     
  
  while(tr->ntips < tr->mxtips) 
    {	
      nodeptr q;
      
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
        
      if(tr->grouped && !full)
	{
	  int 
	    number = p->back->number;	  	 

	  tr->constraintVector[number] = -9;
	}
          
      stepwiseAddition(tr, q, f->back);      	  	 
      
      {
	nodeptr	  
	  r = tr->insertNode->back;
	
	int counter = 4;
	
	hookupDefault(q->next,       tr->insertNode, tr->numBranches);
	hookupDefault(q->next->next, r, tr->numBranches);
	
	computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips, FALSE);              
	tr->ti[0] = counter;
	
	newviewParsimonyIterativeFast(tr);	
      }
    }    
  
  /*printf("ADD: %d\n", tr->bestParsimony); */
  
  nodeRectifier(tr);
  
  randomMP = tr->bestParsimony;        
  
  do
    {
      startMP = randomMP;
      nodeRectifier(tr);
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{
	  rearrangeParsimony(tr, tr->nodep[i], 1, 20, FALSE);
	  if(tr->bestParsimony < randomMP)
	    {		
	      restoreTreeRearrangeParsimony(tr);
	      randomMP = tr->bestParsimony;
	    }
	}      		  	   
    }
  while(randomMP < startMP);
  
  /*printf("OPT: %d %f\n", tr->bestParsimony, gettime() - t);*/

  
     
  free(perm);  
  free(tr->parsimonyScore);
  
  for(model = 0; model < (size_t) tr->NumberOfModels; model++)
    free(tr->partitionData[model].parsVect);
  
  free(tr->ti);
} 



 
static void insertRandom (nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, numBranches);
  hookupDefault(p->next->next, r, numBranches); 
} 







static void buildSimpleTreeRandom (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertRandom(s, p, tr->numBranches);
}

int checker(tree *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->mxtips))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
	return group;

      group = checker(tr, p->next->back);
      if(group != -9) 
	return group;

      group = checker(tr, p->next->next->back);
      if(group != -9) 
	return group;

      return -9;
    }
}







static int markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp)
{
  if(isTip(p->number, numsp))
    return 0;
  else
    {
      branches[*counter] = p->next;
      branches[*counter + 1] = p->next->next;
      
      *counter = *counter + 2;
      
      return ((2 + markBranches(branches, p->next->back, counter, numsp) + 
	       markBranches(branches, p->next->next->back, counter, numsp)));
    }
}



nodeptr findAnyTip(nodeptr p, int numsp)
{ 
  return  isTip(p->number, numsp) ? p : findAnyTip(p->next->back, numsp);
} 


int randomInt(int n)
{
  return rand() %n;
}

void makePermutation(int *perm, int n, analdef *adef)
{    
  int  i, j, k;

  checkSeed(adef);          

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {    
      k =  (int)((double)(n + 1 - i) * randum(&adef->parsimonySeed));

      assert(i + k <= n);
      
      j        = perm[i];
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }
}








boolean tipHomogeneityChecker(tree *tr, nodeptr p, int grouping)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(tr->constraintVector[p->number] != grouping) 
	return FALSE;
      else 
	return TRUE;
    }
  else
    {   
      return  (tipHomogeneityChecker(tr, p->next->back, grouping) && tipHomogeneityChecker(tr, p->next->next->back,grouping));      
    }
}







void makeRandomTree(tree *tr, analdef *adef)
{  
  nodeptr p, f, randomBranch;    
  int nextsp;
  int *perm, branchCounter;
  nodeptr *branches;
  
  branches = (nodeptr *)malloc(sizeof(nodeptr) * (2 * tr->mxtips));
  perm = (int *)malloc((tr->mxtips + 1) * sizeof(int));                         
  
  makePermutation(perm, tr->mxtips, adef);              
  
  tr->ntips = 0;       	       
  tr->nextnode = tr->mxtips + 1;    
  
  buildSimpleTreeRandom(tr, perm[1], perm[2], perm[3]);
  
  while (tr->ntips < tr->mxtips) 
    {	       
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];
      
      /*printf("ADDING SPECIES %d\n", nextsp);*/
      
      buildNewTip(tr, p);  	
      
      f = findAnyTip(tr->start, tr->mxtips);
      f = f->back;
      
      branchCounter = 1;
      branches[0] = f;
      markBranches(branches, f, &branchCounter, tr->mxtips);

      assert(branchCounter == ((2 * (tr->ntips - 1)) - 3));
      
      randomBranch = branches[randomInt(branchCounter)];
      
      insertRandom(p->back, randomBranch, tr->numBranches);
      
    }
  free(perm);            
  free(branches);
}



void nodeRectifier(tree *tr)
{
  nodeptr *np = (nodeptr *)malloc(2 * tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1];
  tr->rooted      = FALSE;

  /* TODO why is tr->rooted set to FALSE here ?*/
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count); 

 
  free(np);
}


void makeParsimonyTree(tree *tr, analdef *adef)
{  
  makeParsimonyTreeFast(tr, adef, TRUE);        
}

void makeParsimonyTreeIncomplete(tree *tr, analdef *adef)
{    
  makeParsimonyTreeFast(tr, adef, FALSE);        
}
 

static void setupBranchMetaInfo(tree *tr, nodeptr p, int nTips, branchInfo *bInf)
{
  int 
    countBranches = tr->branchCounter;

  if(isTip(p->number, tr->mxtips))    
    {      
      p->bInf       = &bInf[countBranches];
      p->back->bInf = &bInf[countBranches];               	      

      bInf[countBranches].oP = p;
      bInf[countBranches].oQ = p->back;
      
      bInf[countBranches].epa->leftNodeNumber = p->number;
      bInf[countBranches].epa->rightNodeNumber = p->back->number;
         
      bInf[countBranches].epa->branchNumber = countBranches;	                 
      bInf[countBranches].epa->originalBranchLength = p->z[0];

      tr->branchCounter =  tr->branchCounter + 1;
      return;
    }
  else
    {
      nodeptr q;
      assert(p == p->next->next->next);

      p->bInf       = &bInf[countBranches];
      p->back->bInf = &bInf[countBranches];

      bInf[countBranches].oP = p;
      bInf[countBranches].oQ = p->back;

      bInf[countBranches].epa->leftNodeNumber = p->number;
      bInf[countBranches].epa->rightNodeNumber = p->back->number;

               
      bInf[countBranches].epa->branchNumber = countBranches;
      bInf[countBranches].epa->originalBranchLength = p->z[0];

      tr->branchCounter =  tr->branchCounter + 1;      

      q = p->next;

      while(q != p)
	{
	  setupBranchMetaInfo(tr, q->back, nTips, bInf);	
	  q = q->next;
	}
     
      return;
    }
}
 


static void setupJointFormat(tree *tr, nodeptr p, int ntips, branchInfo *bInf, int *count)
{
  if(isTip(p->number, tr->mxtips))    
    {      
      p->bInf->epa->jointLabel = *count;
      *count = *count + 1;
           
      return;
    }
  else
    {                           
      setupJointFormat(tr, p->next->back, ntips, bInf, count);            
      setupJointFormat(tr, p->next->next->back, ntips, bInf, count);     
      
      p->bInf->epa->jointLabel = *count;
      *count = *count + 1; 
      
      return;
    }
}
 





static void setupBranchInfo(tree *tr, nodeptr q)
{
  nodeptr 
    originalNode = tr->nodep[tr->mxtips + 1];

  int 
    count = 0;

  tr->branchCounter = 0;

  setupBranchMetaInfo(tr, q, tr->ntips, tr->bInf);
    
  assert(tr->branchCounter == tr->numberOfBranches);

  if(tr->wasRooted)
    {
      assert(tr->leftRootNode->back == tr->rightRootNode);
      assert(tr->leftRootNode       == tr->rightRootNode->back);      

      if(!isTip(tr->leftRootNode->number, tr->mxtips))
	{
	  setupJointFormat(tr,  tr->leftRootNode->next->back, tr->ntips, tr->bInf, &count);
	  setupJointFormat(tr,  tr->leftRootNode->next->next->back, tr->ntips, tr->bInf, &count);
	}
      
       tr->leftRootNode->bInf->epa->jointLabel = count;
       tr->rootLabel = count;
       count = count + 1;

       if(!isTip(tr->rightRootNode->number, tr->mxtips))
	 {
	  setupJointFormat(tr,  tr->rightRootNode->next->back, tr->ntips, tr->bInf, &count);
	  setupJointFormat(tr,  tr->rightRootNode->next->next->back, tr->ntips, tr->bInf, &count);
	}	       
    }
  else
    {
      setupJointFormat(tr, originalNode->back, tr->ntips, tr->bInf, &count);
      setupJointFormat(tr, originalNode->next->back, tr->ntips, tr->bInf, &count);
      setupJointFormat(tr, originalNode->next->next->back, tr->ntips, tr->bInf, &count);      
    }  

  assert(count == tr->numberOfBranches);
}

static void testInsertFast(tree *tr, nodeptr r, nodeptr q)
{
  unsigned int
    result;
  
  nodeptr  
    x = q->back;      
  
  int 
    i,
    *inserts = tr->inserts;
    	           
  assert(!tr->grouped);                             
 
  hookupDefault(r->next,       q, tr->numBranches);
  hookupDefault(r->next->next, x, tr->numBranches);	                         
   
  newviewParsimony(tr, r);   
    
  for(i = 0; i < tr->numberOfTipsForInsertion; i++)
    {           	  	    
      hookupDefault(r, tr->nodep[inserts[i]], tr->numBranches);
      
      tr->bestParsimony = INT_MAX;

      result = evaluateParsimony(tr, r, FALSE);            
      
      r->back = (nodeptr) NULL;
      tr->nodep[inserts[i]]->back = (nodeptr) NULL;
      
      tr->bInf[q->bInf->epa->branchNumber].epa->parsimonyScore[i] = result;	  	         
    }
 
  hookupDefault(q, x, tr->numBranches);
  
  r->next->next->back = r->next->back = (nodeptr) NULL;
 
}


static void traverseTree(tree *tr, nodeptr r, nodeptr q)
{       
  testInsertFast(tr, r, q);

  if(!isTip(q->number, tr->rdta->numsp))
    {   
      nodeptr 
	a = q->next;

      while(a != q)
	{
	  traverseTree(tr, r, a->back);
	  a = a->next;
	}      
    }
} 



typedef struct
  {
    unsigned int parsimonyScore;  
    int number;
  }
  infoMP;


static int infoCompare(const void *p1, const void *p2)
{
  infoMP *rc1 = (infoMP *)p1;
  infoMP *rc2 = (infoMP *)p2;

  unsigned int i = rc1->parsimonyScore;
  unsigned int j = rc2->parsimonyScore;

  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}

void classifyMP(tree *tr, analdef *adef)
{    
  int  
    *informative = (int *)malloc(sizeof(int) * (size_t)tr->cdta->endsite),
    i, 
    j,  
    *perm;    
  
  infoMP
    *inf;
    
  nodeptr     
    r, 
    q;    

  char   
    jointFormatTreeFileName[1024],  
    originalLabelledTreeFileName[1024],
    labelledTreeFileName[1024],
    likelihoodWeightsFileName[1024];
              
  FILE    
    *likelihoodWeightsFile,
    *treeFile;
  
  unsigned int 
    score;        

  assert(adef->restart);

  determineUninformativeSites(tr, informative);     

  compressDNA(tr, informative, TRUE);

  free(informative); 

  tr->ti = (int*)malloc(sizeof(int) * 4 * (size_t)tr->mxtips);    
  
  tr->numberOfBranches = 2 * tr->ntips - 3;

  printBothOpen("\nRAxML Evolutionary Placement Algorithm using parsimony\n"); 

  tr->bestParsimony = INT_MAX;

  score = evaluateParsimony(tr, tr->start->back, TRUE);

  printBothOpen("\nparsimony score of reference tree: %u\n\n", score);

  perm        = (int *)calloc(((size_t)tr->mxtips) + 1, sizeof(int));
  tr->inserts = (int *)calloc((size_t)tr->mxtips,   sizeof(int));

  markTips(tr->start,       perm, tr->mxtips);
  markTips(tr->start->back, perm ,tr->mxtips);

  tr->numberOfTipsForInsertion = 0;

  for(i = 1; i <= tr->mxtips; i++)
    {
      if(perm[i] == 0)
	{
	  tr->inserts[tr->numberOfTipsForInsertion] = i;
	  tr->numberOfTipsForInsertion = tr->numberOfTipsForInsertion + 1;
	}
    }    

  free(perm);
  
  printBothOpen("RAxML will place %d Query Sequences into the %d branches of the reference tree with %d taxa\n\n",  tr->numberOfTipsForInsertion, (2 * tr->ntips - 3), tr->ntips);  

  assert(tr->numberOfTipsForInsertion == (tr->mxtips - tr->ntips));      

  tr->bInf              = (branchInfo*)malloc(tr->numberOfBranches * sizeof(branchInfo)); 

  for(i = 0; i < tr->numberOfBranches; i++)
    {      
      tr->bInf[i].epa = (epaBranchData*)malloc(sizeof(epaBranchData));                    

      tr->bInf[i].epa->parsimonyScore = (unsigned int*)malloc(tr->numberOfTipsForInsertion * sizeof(unsigned int));

      for(j = 0; j < tr->numberOfTipsForInsertion; j++)
	tr->bInf[i].epa->parsimonyScore[j] = INT_MAX;
                
      tr->bInf[i].epa->branchNumber = i;
      tr->bInf[i].epa->countThem = (int*)calloc(tr->numberOfTipsForInsertion, sizeof(int));
      
      sprintf(tr->bInf[i].epa->branchLabel, "I%d", i);     
    } 

  r = tr->nodep[(tr->nextnode)++]; 
    

  q = findAnyTip(tr->start, tr->rdta->numsp);

  assert(isTip(q->number, tr->rdta->numsp));
  assert(!isTip(q->back->number, tr->rdta->numsp));
	 
  q = q->back; 
  
  setupBranchInfo(tr, q);   
 
  traverseTree(tr, r, q);
                   
  printBothOpen("Overall Classification time: %f\n\n", gettime() - masterTime);	
  
 
  strcpy(jointFormatTreeFileName,      workdir);  
  strcpy(originalLabelledTreeFileName, workdir);  
  strcpy(labelledTreeFileName,         workdir);  
  strcpy(likelihoodWeightsFileName,    workdir);

  strcat(jointFormatTreeFileName,      "RAxML_portableTree."); 
  strcat(originalLabelledTreeFileName, "RAxML_originalLabelledTree.");
  strcat(labelledTreeFileName,         "RAxML_labelledTree.");
  strcat(likelihoodWeightsFileName,    "RAxML_equallyParsimoniousPlacements.");
  
  strcat(jointFormatTreeFileName,      run_id); 
  strcat(originalLabelledTreeFileName, run_id);
  strcat(labelledTreeFileName,         run_id);
  strcat(likelihoodWeightsFileName,    run_id);

  strcat(jointFormatTreeFileName,      ".jplace");
  
  free(tr->tree_string);

  tr->treeStringLength *= 16;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char));  
  
 

  treeFile = myfopen(originalLabelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, FALSE, FALSE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile); 

  treeFile = myfopen(jointFormatTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, TRUE, FALSE);
  
  fprintf(treeFile, "{\n");
  fprintf(treeFile, "\t\"tree\": \"%s\", \n", tr->tree_string);
  fprintf(treeFile, "\t\"placements\": [\n");
      

  inf = (infoMP*)malloc(sizeof(infoMP) * tr->numberOfBranches); 
  
  /* joint format */
        
  

  for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
    {
      unsigned int
	lmax;
      
      int 	   
	validEntries = tr->numberOfBranches;
      
      for(j =  0; j < tr->numberOfBranches; j++) 
	{
	  inf[j].parsimonyScore = tr->bInf[j].epa->parsimonyScore[i];	  
	  inf[j].number         = tr->bInf[j].epa->jointLabel;
	}
      
      qsort(inf, tr->numberOfBranches, sizeof(infoMP), infoCompare);	            
      
      j = 0;
      
      lmax = inf[0].parsimonyScore;
      
      fprintf(treeFile, "\t{\"p\":[");
      
      while(j < validEntries && inf[j].parsimonyScore == lmax)	  
	{ 	    
	  if(j > 0)
	    {
	      if(tr->wasRooted && inf[j].number == tr->rootLabel)		  
		assert(0);		 
	      else
		fprintf(treeFile, ",[%d, %u]", inf[j].number, inf[j].parsimonyScore);
	    }
	  else
	    {
	      if(tr->wasRooted && inf[j].number == tr->rootLabel)		  
		assert(0);		  
	      else
		fprintf(treeFile, "[%d, %u]", inf[j].number, inf[j].parsimonyScore);
	    }	  
	  
	  j++;
	}
      
      if(i == tr->numberOfTipsForInsertion - 1)
	fprintf(treeFile, "], \"n\":[\"%s\"]}\n", tr->nameList[tr->inserts[i]]);
      else
	fprintf(treeFile, "], \"n\":[\"%s\"]},\n", tr->nameList[tr->inserts[i]]);
    }  
    
  fprintf(treeFile, "\t ],\n");
  /*  fprintf(treeFile, "\t\"metadata\": {\"invocation\": \"RAxML EPA parsimony\"},\n");*/
  fprintf(treeFile, "\t\"metadata\": {\"invocation\": ");

  fprintf(treeFile, "\"");
  
  {
    int i;
    
    for(i = 0; i < globalArgc; i++)
      fprintf(treeFile,"%s ", globalArgv[i]);
  }
  fprintf(treeFile, "\", \"raxml_version\": \"%s\"", programVersion);
  fprintf(treeFile,"},\n");

  fprintf(treeFile, "\t\"version\": 2,\n");
  fprintf(treeFile, "\t\"fields\": [\n");
  fprintf(treeFile, "\t\"edge_num\", \"parsimony\"\n");
  fprintf(treeFile, "\t]\n");
  fprintf(treeFile, "}\n");
  
  fclose(treeFile);

  /* JSON format end */ 
           
  likelihoodWeightsFile = myfopen(likelihoodWeightsFileName, "wb");
 
        
  for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
    {
      unsigned int	 
	lmax;
	
      int 
	validEntries = tr->numberOfBranches;
	
      for(j =  0; j < tr->numberOfBranches; j++) 
	{
	  inf[j].parsimonyScore = tr->bInf[j].epa->parsimonyScore[i];
	  inf[j].number = j;
	}
	
      qsort(inf, tr->numberOfBranches, sizeof(infoMP), infoCompare);	 	     
		
      j = 0;
      
      lmax = inf[0].parsimonyScore;
      
      while(j < validEntries && inf[j].parsimonyScore == lmax)	  
	{ 	   	    
	  fprintf(likelihoodWeightsFile, "%s I%d %u\n", tr->nameList[tr->inserts[i]], inf[j].number, inf[j].parsimonyScore);
	  tr->bInf[inf[j].number].epa->countThem[i] = 1;
	  j++;
	}			      	   
    }
      
  free(inf);
   
  fclose(likelihoodWeightsFile); 
    
  
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, FALSE, FALSE, FALSE);
  treeFile = fopen(labelledTreeFileName, "wb");
  fprintf(treeFile, "%s\n", tr->tree_string);
  fclose(treeFile);
  

  printBothOpen("Equally parsimonious read placements written to file: %s\n\n", likelihoodWeightsFileName);
  printBothOpen("Labelled reference tree with branch labels (without query sequences) written to file: %s\n\n", originalLabelledTreeFileName); 
  printBothOpen("Labelled reference tree with branch labels in portable pplacer/EPA format (without query sequences) written to file: %s\n\n", jointFormatTreeFileName);
  printBothOpen("Labelled reference tree including branch labels and query sequences written to file: %s\n\n", labelledTreeFileName); 
  
  exit(0);
}

#ifdef __AVX

#ifdef _SSE3_WAS_DEFINED

#define __SIM_SSE3

#undef _SSE3_WAS_DEFINED

#endif

#endif
