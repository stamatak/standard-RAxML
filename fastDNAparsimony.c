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

#ifdef __SIM_SSE3

#include <xmmintrin.h>
#include <pmmintrin.h>
  
#endif

#include "axml.h"

extern double masterTime;
extern char  permFileName[1024];
extern char  infoFileName[1024];
extern const unsigned int mask32[32];

/*
  #ifdef _USE_PTHREADS
  extern volatile int NumberOfThreads;
  extern volatile int *reductionBufferParsimony;
  #endif
*/


/********************************DNA FUNCTIONS *****************************************************************/


/*static void hookupParsimony(nodeptr p, nodeptr q)
{
  p->back = q;
  q->back = p;
  }*/

static void getxnodeLocal (nodeptr p)
{
  nodeptr  s;

  if((s = p->next)->x || (s = s->next)->x)
    {
      p->x = s->x;
      s->x = 0;
    }
}

static void computeTraversalInfoParsimony(nodeptr p, int *ti, int *counter, int maxTips)
{        
  nodeptr 
    q = p->next->back,
    r = p->next->next->back;
  
  if(! p->x)
    getxnodeLocal(p);  
  
  if(q->number > maxTips && !q->x) 
    computeTraversalInfoParsimony(q, ti, counter, maxTips);
  
  if(r->number > maxTips && !r->x) 
    computeTraversalInfoParsimony(r, ti, counter, maxTips);
  
  
  ti[*counter]     = p->number;
  ti[*counter + 1] = q->number;
  ti[*counter + 2] = r->number;
  *counter = *counter + 4;
}

#define BIT_COUNT(x)  precomputed16_bitcount(x)

/*
  __builtin_popcount(x)
*/



#ifdef __SIM_SSE3

static unsigned int vectorCount(__m128i b)
{
  const unsigned int 
    mu1 = 0x55555555,
    mu2 = 0x33333333,
    mu3 = 0x0F0F0F0F,
    mu4 = 0x0000003F;
  
  unsigned int 
    tcnt[4] __attribute__ ((aligned (16)));
  
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



#ifdef __SIM_SSE3

void newviewParsimonyIterativeFast(tree *tr)
{    
  __m128i 
    allOne = _mm_set_epi32 (0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
  
  int     
    *ti = tr->ti,
    count = ti[0],
    index;

  unsigned int  
    width = tr->compressedWidth;

  /*
    #ifdef _USE_PTHREADS
    const int 
    tid = tr->threadID,
    n = NumberOfThreads;    
    #endif
  */

  for(index = 4; index < count; index += 4)
    {        
      int 
	pNumber = ti[index],
	qNumber = ti[index + 1],
	rNumber = ti[index + 2];
      
      parsimonyNumber      
	*leftState_A  = tr->parsimonyState_A[qNumber],
	*rightState_A = tr->parsimonyState_A[rNumber],
	*thisState_A  = tr->parsimonyState_A[pNumber],
	
	*leftState_C  = tr->parsimonyState_C[qNumber],
	*rightState_C = tr->parsimonyState_C[rNumber],
	*thisState_C  = tr->parsimonyState_C[pNumber],

	*leftState_G  = tr->parsimonyState_G[qNumber],
	*rightState_G = tr->parsimonyState_G[rNumber],
	*thisState_G  = tr->parsimonyState_G[pNumber],

	*leftState_T  = tr->parsimonyState_T[qNumber],
	*rightState_T = tr->parsimonyState_T[rNumber],
	*thisState_T  = tr->parsimonyState_T[pNumber]; 	           
      
      
      
      unsigned int	
	i,            
	ts = 0;      
                 
       for(i = 0; i < width; i+=4)
	{
	  /*
	    #ifdef _USE_PTHREADS
	    if((i / 4) % n == tid)
	    {
	    #endif	 
	  */
	  
	  __m128i
	    s_r, s_l, v_N,
	    l_A, l_C, l_G, l_T,
	    v_A, v_C, v_G, v_T;
	    
	  s_l = _mm_load_si128((__m128i *)(&leftState_A[i]));
	  s_r = _mm_load_si128((__m128i *)(&rightState_A[i]));
	  l_A = _mm_and_si128(s_l, s_r);
	  v_A = _mm_or_si128(s_l, s_r);

	  s_l = _mm_load_si128((__m128i *)(&leftState_C[i]));
	  s_r = _mm_load_si128((__m128i *)(&rightState_C[i]));
	  l_C = _mm_and_si128(s_l, s_r);
	  v_C = _mm_or_si128(s_l, s_r);
	  
	  s_l = _mm_load_si128((__m128i *)(&leftState_G[i]));
	  s_r = _mm_load_si128((__m128i *)(&rightState_G[i]));
	  l_G = _mm_and_si128(s_l, s_r);
	  v_G = _mm_or_si128(s_l, s_r);

	  s_l = _mm_load_si128((__m128i *)(&leftState_T[i]));
	  s_r = _mm_load_si128((__m128i *)(&rightState_T[i]));
	  l_T = _mm_and_si128(s_l, s_r);
	  v_T = _mm_or_si128(s_l, s_r);
	  	   	    
	  v_N = _mm_or_si128(_mm_or_si128(l_A, l_C), _mm_or_si128(l_G, l_T));	  	 	    	  
	 
	  _mm_store_si128((__m128i *)(&thisState_A[i]), _mm_or_si128(l_A, _mm_andnot_si128(v_N, v_A)));
	  _mm_store_si128((__m128i *)(&thisState_C[i]), _mm_or_si128(l_C, _mm_andnot_si128(v_N, v_C)));
	  _mm_store_si128((__m128i *)(&thisState_G[i]), _mm_or_si128(l_G, _mm_andnot_si128(v_N, v_G)));
	  _mm_store_si128((__m128i *)(&thisState_T[i]), _mm_or_si128(l_T, _mm_andnot_si128(v_N, v_T)));	  	 	 	  
	  	  
	  v_N = _mm_andnot_si128(v_N, allOne);
	  
	  ts += vectorCount(v_N);
	  
	  /*
	    #ifdef _USE_PTHREADS
	    }
	    #endif
	  */
	}	     
	  
      tr->parsimonyScore[pNumber] = ts + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];          	 	    	      	  	
    }
}


unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  __m128i 
    allOne = _mm_set_epi32 (0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);

  int   
    pNumber = tr->ti[1],
    qNumber = tr->ti[2];
  
  unsigned int 
    bestScore = tr->bestParsimony,
    width = tr->compressedWidth,
    i,   
    sum;

  parsimonyNumber
    *rightState_A = tr->parsimonyState_A[pNumber], 
    *leftState_A  = tr->parsimonyState_A[qNumber],
    
    *rightState_C = tr->parsimonyState_C[pNumber], 
    *leftState_C  = tr->parsimonyState_C[qNumber],

    *rightState_G = tr->parsimonyState_G[pNumber], 
    *leftState_G  = tr->parsimonyState_G[qNumber],

    *rightState_T = tr->parsimonyState_T[pNumber], 
    *leftState_T  = tr->parsimonyState_T[qNumber];
  
  /*
    #ifdef _USE_PTHREADS
    const int 
    tid = tr->threadID,
    n = NumberOfThreads;  
    #endif
  */
            
  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr);       
  
  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];
  
  for(i = 0; i < width; i+=4)
    {     
      /*
	#ifdef _USE_PTHREADS
	if((i / 4) % n == tid)
	{
	#endif	
      */
      
      unsigned int 
	counts[4] __attribute__ ((aligned (16)));
                 
      __m128i      
	l_A = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_A[i])), _mm_load_si128((__m128i *)(&rightState_A[i]))),
	l_C = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_C[i])), _mm_load_si128((__m128i *)(&rightState_C[i]))),
	l_G = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_G[i])), _mm_load_si128((__m128i *)(&rightState_G[i]))),
	l_T = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_T[i])), _mm_load_si128((__m128i *)(&rightState_T[i]))),
	v_N = _mm_or_si128(_mm_or_si128(l_A, l_C), _mm_or_si128(l_G, l_T));     

      _mm_store_si128((__m128i *)counts, _mm_andnot_si128(v_N, allOne));
                      
      sum += BIT_COUNT(counts[0]) + BIT_COUNT(counts[1]);
      sum += BIT_COUNT(counts[2]) + BIT_COUNT(counts[3]);

      if(sum >= bestScore)
	break;
      
      /*
	#ifdef _USE_PTHREADS
	}
	#endif  
      */
    }   
  
  return sum;      
}


#else



void newviewParsimonyIterativeFast(tree *tr)
{  
  int
    count = tr->ti[0],
    *ti   = tr->ti,
    index;

  unsigned int
    width = tr->compressedWidth;

  for(index = 4; index < count; index += 4)
    {        
      int 
	pNumber = ti[index],
	qNumber = ti[index + 1],
	rNumber = ti[index + 2];
      
      parsimonyNumber      
	*leftState_A  = tr->parsimonyState_A[qNumber],
	*rightState_A = tr->parsimonyState_A[rNumber],
	*thisState_A  = tr->parsimonyState_A[pNumber],
	
	*leftState_C  = tr->parsimonyState_C[qNumber],
	*rightState_C = tr->parsimonyState_C[rNumber],
	*thisState_C  = tr->parsimonyState_C[pNumber],

	*leftState_G  = tr->parsimonyState_G[qNumber],
	*rightState_G = tr->parsimonyState_G[rNumber],
	*thisState_G  = tr->parsimonyState_G[pNumber],

	*leftState_T  = tr->parsimonyState_T[qNumber],
	*rightState_T = tr->parsimonyState_T[rNumber],
	*thisState_T  = tr->parsimonyState_T[pNumber];
      
      register parsimonyNumber
	o_A,
	o_C,
	o_G,
	o_T,
	t_A,
	t_C,
	t_G,
	t_T,	
	t_N;
      
      unsigned int 
	i,	   
	ts = 0;             	         	  

      for(i = 0; i < width; i++)
	{		 	    
	  t_A = leftState_A[i] & rightState_A[i];
	  t_C = leftState_C[i] & rightState_C[i];
	  t_G = leftState_G[i] & rightState_G[i];	  
	  t_T = leftState_T[i] & rightState_T[i];

	  o_A = leftState_A[i] | rightState_A[i];
	  o_C = leftState_C[i] | rightState_C[i];
	  o_G = leftState_G[i] | rightState_G[i];	  
	  o_T = leftState_T[i] | rightState_T[i];

	  t_N = ~(t_A | t_C | t_G | t_T);	  

	  thisState_A[i] = t_A | (t_N & o_A);
	  thisState_C[i] = t_C | (t_N & o_C);
	  thisState_G[i] = t_G | (t_N & o_G);
	  thisState_T[i] = t_T | (t_N & o_T);	 	 	  
	  
	  ts += BIT_COUNT(t_N);   	   		      
	}	
	  
      tr->parsimonyScore[pNumber] = ts + tr->parsimonyScore[rNumber] + tr->parsimonyScore[qNumber];          	 	    	      	  	
    }
}


unsigned int evaluateParsimonyIterativeFast(tree *tr)
{
  int   
    pNumber = tr->ti[1],
    qNumber = tr->ti[2];  
  
  unsigned int 
    bestScore = tr->bestParsimony,
    width = tr->compressedWidth,
    i,   
    sum, 
    t_A,
    t_C,
    t_G,
    t_T,
    t_N;

  parsimonyNumber
    *rightState_A = tr->parsimonyState_A[pNumber], 
    *leftState_A  = tr->parsimonyState_A[qNumber],
    
    *rightState_C = tr->parsimonyState_C[pNumber], 
    *leftState_C  = tr->parsimonyState_C[qNumber],

    *rightState_G = tr->parsimonyState_G[pNumber], 
    *leftState_G  = tr->parsimonyState_G[qNumber],

    *rightState_T = tr->parsimonyState_T[pNumber], 
    *leftState_T  = tr->parsimonyState_T[qNumber];
            
  if(tr->ti[0] > 4)
    newviewParsimonyIterativeFast(tr);       
  
  sum = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];
  
  for(i = 0; i < width; i++)
    {                
       t_A = leftState_A[i] & rightState_A[i];
       t_C = leftState_C[i] & rightState_C[i];
       t_G = leftState_G[i] & rightState_G[i];	  
       t_T = leftState_T[i] & rightState_T[i];

       t_N = ~(t_A | t_C | t_G | t_T);

       sum += BIT_COUNT(t_N);     
       if(sum >= bestScore)
	 break;
    }         

  return sum;      
}

#endif






static unsigned int evaluateParsimony(tree *tr, nodeptr p)
{
  volatile unsigned int result;
  nodeptr q = p->back;
  int
    *ti = tr->ti,
    counter = 4;
  
  ti[1] = p->number;
  ti[2] = q->number;

  if(p->number > tr->mxtips && !p->x)
    computeTraversalInfoParsimony(p, ti, &counter, tr->mxtips);
  if(q->number > tr->mxtips && !q->x)
    computeTraversalInfoParsimony(q, ti, &counter, tr->mxtips); 
  
  ti[0] = counter;

  /*
    #ifdef _USE_PTHREADS
    {
    int i;    
    
    masterBarrier(THREAD_FAST_EVALUATE_PARSIMONY, tr);    
    
    for(i = 0, result = 0; i < NumberOfThreads; i++)
    result += reductionBufferParsimony[i]; 
    }
    #else
  */
  
  result = evaluateParsimonyIterativeFast(tr);

  /* #endif */

  return result;
}


static void newviewParsimony(tree *tr, nodeptr  p)
{     
  if(p->number <= tr->mxtips)
    return;

  {
    int 
      counter = 4;     
           
    computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips);              
    tr->ti[0] = counter;
         
    /*
      #ifdef _USE_PTHREADS
      masterBarrier(THREAD_FAST_NEWVIEW_PARSIMONY, tr); 
      #else
    */
    
    newviewParsimonyIterativeFast(tr);    

    /* #endif */
  }
}





/****************************************************************************************************************************************/

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q,  tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches); 
   
  newviewParsimony(tr, p);     
} 

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
  unsigned int mp;
 
  nodeptr  r = q->back;   
    
  insertParsimony(tr, p, q);   
  
  mp = evaluateParsimony(tr, p->next->next);          
      
  if(mp < tr->bestParsimony)
    {
      tr->bestParsimony = mp;
      tr->insertNode = q;
      tr->removeNode = p;
    }
  
  hookupDefault(q, r, tr->numBranches);
  p->next->next->back = p->next->back = (nodeptr) NULL;
       
  return;
} 


static void restoreTreeParsimony(tree *tr, nodeptr p, nodeptr q)
{ 
  nodeptr
    r = q->back;
  
  int counter = 4;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches);
  
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr); 
  
  /*insertParsimony(tr, p, q);  

  if((p->number > tr->mxtips) && (q->number > tr->mxtips))
    {
      while ((! p->x)) 
	{
	  if (! (p->x))
	    newviewParsimony(tr, p);		     
	}
    }
  if((p->number <= tr->mxtips) && (q->number > tr->mxtips))
    {
      while ((! q->x)) 
	{		  
	  if (! (q->x)) 
	    newviewParsimony(tr, q);
	}
    }
  if((p->number > tr->mxtips) && (q->number > tr->mxtips))
    {
      while ((! p->x) || (! q->x)) 
	{
	  if (! (p->x))
	    newviewParsimony(tr, p);
	  if (! (q->x))
	    newviewParsimony(tr, q);
	}
	}*/	
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

static int randomIntFast(int n)
{
  return rand() % n;
}

static void makePermutationFast(int *perm, int n, analdef *adef)
{    
  int  i, j, k;


  if(adef->parsimonySeed == 0)   
    srand((unsigned int) gettimeSrand());          


  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {
      if(adef->parsimonySeed == 0) 
	k        = randomIntFast(n + 1 - i);
      else
	k =  (int)((double)(n + 1 - i) * randum(&adef->parsimonySeed));
      
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
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2; 
           
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3; 

  
  if(maxtrav < mintrav)
    return 0;

  q = p->back;

  if(p->number > tr->mxtips) 
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
       
  if ((q->number > tr->mxtips) && maxtrav > 0) 
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

static boolean isInformative(tree *tr, int dataType, int site)
{
  int
    informativeCounter = 0,
    check[256],   
    j,   
    undetermined = getUndetermined(dataType);

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
  
  /*
    printf("Uninformative Patterns: %d\n", number);
  */
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
     
      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}

static void nodeRectifierFast(tree *tr)
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
  
static parsimonyNumber *compressDNA(tree *tr, int *informative)
{
  int
    i,
    entries = 0,   
    compressedEntries,
    compressedEntriesPadded;

  parsimonyNumber
    *compressedScratch;
  
  for(i = 0; i < tr->cdta->endsite; i++)    
    if(informative[i])
      entries += tr->cdta->aliaswgt[i];     
  
  compressedEntries = entries / PCF;

  if(entries % PCF != 0)
    compressedEntries++;
  
  
  /*printf("compression %d -> %d\n", entries, compressedEntries);  */

#ifdef __SIM_SSE3
  if(compressedEntries % 4 != 0)
    compressedEntriesPadded = compressedEntries + (4 - (compressedEntries % 4));
  else
    compressedEntriesPadded = compressedEntries;
#else
  compressedEntriesPadded = compressedEntries;
#endif

  /* printf("padded %d\n", compressedEntriesPadded); */

  compressedScratch = (parsimonyNumber *)malloc_aligned(compressedEntriesPadded * 8 * tr->mxtips * sizeof(parsimonyNumber));
     
  for(i = 0; i < compressedEntriesPadded * 8 * tr->mxtips; i++)      
    compressedScratch[i] = 0;    
      

  for(i = 0; i < tr->mxtips; i++)
    {
      parsimonyNumber
	*compressedTip_A = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1))],
	*compressedTip_C = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1)) + compressedEntriesPadded * 1],
	*compressedTip_G = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1)) + compressedEntriesPadded * 2],
	*compressedTip_T = &compressedScratch[(compressedEntriesPadded * 4 * (i + 1)) + compressedEntriesPadded * 3],
	compressedValue_A = 0,
	compressedValue_C = 0,
	compressedValue_G = 0,
	compressedValue_T = 0;
      
      int       
	w = 0,
	compressedIndex = 0,
	compressedCounter = 0,
	index = 0;
	      
      for(index = 0; index < tr->cdta->endsite; index++)
	{
	  if(informative[index])
	    {
	      parsimonyNumber value = (parsimonyNumber)(tr->yVector[i + 1][index]);	  
	      
	      for(w = 0; w < tr->cdta->aliaswgt[index]; w++)
		{	    	     	    
		  if(value & 1)
		    compressedValue_A |= mask32[compressedCounter];
		  if(value & 2)
		    compressedValue_C |= mask32[compressedCounter];
		  if(value & 4)
		    compressedValue_G |= mask32[compressedCounter];
		  if(value & 8)
		    compressedValue_T |= mask32[compressedCounter];
	          
		  compressedCounter++;
		  
		  if(compressedCounter == PCF)
		    {
		      compressedTip_A[compressedIndex] = compressedValue_A;
		      compressedTip_C[compressedIndex] = compressedValue_C;
		      compressedTip_G[compressedIndex] = compressedValue_G;
		      compressedTip_T[compressedIndex] = compressedValue_T;
		      
		      compressedValue_A = 0;
		      compressedValue_C = 0;
		      compressedValue_G = 0;
		      compressedValue_T = 0;
		      
		      compressedCounter = 0;
		      compressedIndex++;
		    }
		}
	    }
	}
                           
      for(;compressedIndex < compressedEntriesPadded; compressedIndex++)
	{	
	  for(;compressedCounter < PCF; compressedCounter++)
	    {	        	      	      
	      compressedValue_A |= mask32[compressedCounter];		 
	      compressedValue_C |= mask32[compressedCounter];		 
	      compressedValue_G |= mask32[compressedCounter];		  
	      compressedValue_T |= mask32[compressedCounter];	      	     
	    }
	  
	   compressedTip_A[compressedIndex] = compressedValue_A;
	   compressedTip_C[compressedIndex] = compressedValue_C;
	   compressedTip_G[compressedIndex] = compressedValue_G;
	   compressedTip_T[compressedIndex] = compressedValue_T;
	   
	   compressedValue_A = 0;
	   compressedValue_C = 0;
	   compressedValue_G = 0;
	   compressedValue_T = 0;
	   
	   compressedCounter = 0;
	}	 	
    }
    
  tr->parsimonyState_A = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * tr->mxtips);
  tr->parsimonyState_C = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * tr->mxtips);
  tr->parsimonyState_G = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * tr->mxtips);
  tr->parsimonyState_T = (parsimonyNumber**)malloc_aligned(sizeof(parsimonyNumber*) * 2 * tr->mxtips);
  
  tr->parsimonyScore = (unsigned int*)malloc_aligned(sizeof(unsigned int) * 2 * tr->mxtips);  
    
  for(i = 0; i < 2 * tr->mxtips; i++) 
    {
      tr->parsimonyState_A[i] = &compressedScratch[i * 4 * compressedEntriesPadded];       
      tr->parsimonyState_C[i] = &compressedScratch[i * 4 * compressedEntriesPadded + compressedEntriesPadded * 1]; 
      tr->parsimonyState_G[i] = &compressedScratch[i * 4 * compressedEntriesPadded + compressedEntriesPadded * 2]; 
      tr->parsimonyState_T[i] = &compressedScratch[i * 4 * compressedEntriesPadded + compressedEntriesPadded * 3];          
 
      tr->parsimonyScore[i] = 0;
    }
  
  tr->compressedWidth = (unsigned int)compressedEntriesPadded;
  
  return compressedScratch; 
}



static void stepwiseAddition(tree *tr, nodeptr p, nodeptr q)
{            
  nodeptr 
    r = q->back;

  unsigned int 
    mp,
    bestParsimony = tr->bestParsimony;
  
  int 
    counter = 4;
  
  p->next->back = q;
  q->back = p->next;

  p->next->next->back = r;
  r->back = p->next->next;
   
  computeTraversalInfoParsimony(p, tr->ti, &counter, tr->mxtips);              
  tr->ti[0] = counter;
    
  newviewParsimonyIterativeFast(tr);      
  
#ifdef __SIM_SSE3
  {
    __m128i 
      allOne = _mm_set_epi32 (0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF, 0xFFFFFFFF);
    
    int   
      pNumber = p->number,
      qNumber = p->back->number;
  
    unsigned int       
      width = tr->compressedWidth,
      i;
    
    parsimonyNumber
      *rightState_A = tr->parsimonyState_A[pNumber], 
      *leftState_A  = tr->parsimonyState_A[qNumber],
      
      *rightState_C = tr->parsimonyState_C[pNumber], 
      *leftState_C  = tr->parsimonyState_C[qNumber],
      
      *rightState_G = tr->parsimonyState_G[pNumber], 
      *leftState_G  = tr->parsimonyState_G[qNumber],
      
      *rightState_T = tr->parsimonyState_T[pNumber], 
      *leftState_T  = tr->parsimonyState_T[qNumber];
        
    mp = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];

    for(i = 0; i < width; i+=4)
      {     
	unsigned int 
	  counts[4] __attribute__ ((aligned (16)));
	
	__m128i      
	  l_A = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_A[i])), _mm_load_si128((__m128i *)(&rightState_A[i]))),
	  l_C = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_C[i])), _mm_load_si128((__m128i *)(&rightState_C[i]))),
	  l_G = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_G[i])), _mm_load_si128((__m128i *)(&rightState_G[i]))),
	  l_T = _mm_and_si128(_mm_load_si128((__m128i *)(&leftState_T[i])), _mm_load_si128((__m128i *)(&rightState_T[i]))),
	  v_N = _mm_or_si128(_mm_or_si128(l_A, l_C), _mm_or_si128(l_G, l_T));     
	
	_mm_store_si128((__m128i *)counts, _mm_andnot_si128(v_N, allOne));
	
	mp += BIT_COUNT(counts[0]) + BIT_COUNT(counts[1]);
	mp += BIT_COUNT(counts[2]) + BIT_COUNT(counts[3]);  
	if(mp >= bestParsimony)
	  goto SKIP;
      }           
  }
#else
  {
    int   
      pNumber = p->number,
      qNumber = p->back->number;  
    
    unsigned int 
      width = tr->compressedWidth,
      i,         
      t_A,
      t_C,
      t_G,
      t_T,
      t_N;
    
    parsimonyNumber
      *rightState_A = tr->parsimonyState_A[pNumber], 
      *leftState_A  = tr->parsimonyState_A[qNumber],
      
      *rightState_C = tr->parsimonyState_C[pNumber], 
      *leftState_C  = tr->parsimonyState_C[qNumber],
      
      *rightState_G = tr->parsimonyState_G[pNumber], 
      *leftState_G  = tr->parsimonyState_G[qNumber],
      
      *rightState_T = tr->parsimonyState_T[pNumber], 
      *leftState_T  = tr->parsimonyState_T[qNumber];
    
    mp = tr->parsimonyScore[pNumber] + tr->parsimonyScore[qNumber];
    
    for(i = 0; i < width; i++)
      {                
	t_A = leftState_A[i] & rightState_A[i];
	t_C = leftState_C[i] & rightState_C[i];
	t_G = leftState_G[i] & rightState_G[i];	  
	t_T = leftState_T[i] & rightState_T[i];
	
	t_N = ~(t_A | t_C | t_G | t_T);
	
	mp += BIT_COUNT(t_N);
	
	if(mp >= bestParsimony)
	  goto SKIP;
      }            
  }
#endif

  tr->bestParsimony = mp;
  tr->insertNode = q;     
  
 SKIP: 
  q->back = r;
  r->back = q;
   
  if(q->number > tr->mxtips && tr->parsimonyScore[q->number] > 0)
    {	      
      stepwiseAddition(tr, p, q->next->back);	      
      stepwiseAddition(tr, p, q->next->next->back);              	     
    }
}


void makeParsimonyTreeFastDNA(tree *tr, analdef *adef)
{   
  nodeptr  
    p, 
    f;    
  
  int  
    i, 
    nextsp,
    *perm, 
    *informative;
  
  unsigned int 
    randomMP, 
    startMP;        

  parsimonyNumber 
    *parsimonySpace;  
  
  /* stuff for informative sites */
 
  informative = (int *)malloc(sizeof(int) * tr->cdta->endsite);   

  perm = (int *)malloc((tr->mxtips + 1) * sizeof(int));   
  
  /* end */
       
  determineUninformativeSites(tr, informative);     

  makePermutationFast(perm, tr->mxtips, adef);

  parsimonySpace = compressDNA(tr, informative);

  tr->ti = (int*)malloc(sizeof(int) * 4 * tr->mxtips);

  /*
    #ifdef _USE_PTHREADS
    masterBarrier(THREAD_INIT_FAST_PARSIMONY, tr); 
    #endif
  */
  
  tr->ntips = 0;    
  
  tr->nextnode = tr->mxtips + 1;       
 
  buildSimpleTree(tr, perm[1], perm[2], perm[3]);      
 

  /* t = gettime(); */

  f = tr->start;

  while(tr->ntips < tr->mxtips) 
    {	
      nodeptr q;

      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];                 
      q = tr->nodep[(tr->nextnode)++];
      p->back = q;
      q->back = p;
                   
     
      stepwiseAddition(tr, q, f->back);      
     
      /*restoreTreeParsimony(tr, q, tr->insertNode);*/

      {
	nodeptr	  
	  r = tr->insertNode->back;
  
	int counter = 4;
	
	hookupDefault(q->next,       tr->insertNode, tr->numBranches);
	hookupDefault(q->next->next, r, tr->numBranches);
	
	computeTraversalInfoParsimony(q, tr->ti, &counter, tr->mxtips);              
	tr->ti[0] = counter;
	
	newviewParsimonyIterativeFast(tr);	
      }
    }    

  /* printf("ADD: %d\n", tr->bestParsimony); */

  
  /*
    printf("ADD: %f %d\n", gettime() - t, tr->bestParsimony); 

    t = gettime();
  */
  
  nodeRectifierFast(tr);
 
  randomMP = tr->bestParsimony;        
  
  do
    {
      startMP = randomMP;
      nodeRectifierFast(tr);
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
 
  
  /* printf("OPT: %d\n", tr->bestParsimony); */
    /*          
    printf("REARRANGEMENT MP Score %d Time %f\n", tr->bestParsimony, gettime() - masterTime); 
  */
 
  nodeRectifierFast(tr);             

  free(informative);     
  free(perm);    
  free(tr->parsimonyState_A);
  free(tr->parsimonyState_C);
  free(tr->parsimonyState_G);
  free(tr->parsimonyState_T);
  free(parsimonySpace);
 

  free(tr->parsimonyScore);
  free(tr->ti);
} 



 
