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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
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

#include "axml.h"


extern char  workdir[1024];
extern char run_id[128];

extern const char binaryStateNames[2];   
extern const char dnaStateNames[4];
extern const char protStateNames[20];
extern const char genericStateNames[32];

static char getStateCharacter(int dataType, int state)
{
  char 
    result;

  assert(MIN_MODEL < dataType && dataType < MAX_MODEL);

  switch(dataType)
    {
    case BINARY_DATA:
      result = binaryStateNames[state];
      break;
    case DNA_DATA:
       result = dnaStateNames[state];
      break;
    case AA_DATA:
      result =  protStateNames[state];
      break; 
    case GENERIC_32:
      result = genericStateNames[state];
      break;
    default:
      assert(0);
    }

  return  result;
}

static void makeP_Flex_Ancestral(double *rptr, double *EI,  double *EIGN, int numberOfCategories, double *left, const int numStates)
{
  int 
    i,
    j,
    k;
  
  const int
    rates = numStates - 1,
    statesSquare = numStates * numStates;

  double 
    z1 = 0.0,
    lz1[64],
    d1[64];

  assert(numStates <= 64);
     
  for(i = 0; i < rates; i++)    
    lz1[i] = EIGN[i] * z1;
     

  for(i = 0; i < numberOfCategories; i++)
    {
      for(j = 0; j < rates; j++)	
	d1[j] = EXP (rptr[i] * lz1[j]);
	 
      for(j = 0; j < numStates; j++)
	{
	  left[statesSquare * i  + numStates * j] = 1.0;	 

	  for(k = 1; k < numStates; k++)	    
	    left[statesSquare * i + numStates * j + k]  = d1[k-1] * EI[rates * j + (k-1)];	     
	}
    }  
}

static void ancestralCat(double *v, double *sumBuffer, double *diagptable, int i, int numStates)
{
  int 
    l,
    j;
  
  double 
    *ancestral = &sumBuffer[numStates * i],
    sum = 0.0,
    *term = (double*)malloc(sizeof(double) * numStates);

  for(l = 0; l < numStates; l++)
    {
      double 
	ump_x1 = 0.0;
      
      for(j = 0; j < numStates; j++)	
	ump_x1 += v[j] * diagptable[l * numStates + j];

      sum += ump_x1;
      term[l] = ump_x1;
      
    }
		
  for(l = 0; l < numStates; l++)          
    ancestral[l] = term[l] / sum;	
   
  free(term);
}

static void newviewFlexCat_Ancestral(int tipCase, double *extEV,
				     int *cptr,
				     double *x1, double *x2, double *tipVector,
				     unsigned char *tipX1, unsigned char *tipX2,
				     int n, double *left, double *right, const int numStates, double *diagptable, double *sumBuffer)
{
  double   
    *x3 = (double *)malloc(sizeof(double) * numStates),
    *le, *ri, *v, *vl, *vr,
    ump_x1, ump_x2, x1px2;
  int 
    i, l, j, scale;

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

	    v  = x3;

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
	   	    
	    ancestralCat(v, sumBuffer, &diagptable[cptr[i] * statesSquare], i, numStates);
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
	    v  = x3;

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
	      }
	   
	    ancestralCat(v, sumBuffer, &diagptable[cptr[i] * statesSquare], i, numStates);	    
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
	  v = x3;

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
	     }
	  
	   ancestralCat(v, sumBuffer, &diagptable[cptr[i] * statesSquare], i, numStates);	
	}
      break;
    default:
      assert(0);
    }

  free(x3);

 

}



static void ancestralGamma(double *_v, double *sumBuffer, double *diagptable, int i, int numStates, int gammaStates)
{
  int 
    statesSquare = numStates * numStates,
    k,
    j,
    l;

  double 
    *v,
    *ancestral = &sumBuffer[gammaStates * i],
    sum = 0.0,
    *term = (double*)malloc(sizeof(double) * numStates);	      	      

  for(l = 0; l < numStates; l++)
    term[l] = 0.0;

  for(k = 0; k < 4; k++)
    {	  
      v =  &(_v[numStates * k]);

      for(l = 0; l < numStates; l++)
	{
	  double
	    al = 0.0;
	 
	  for(j = 0; j < numStates; j++)	    
	    al += v[j] * diagptable[k * statesSquare + l * numStates + j];
	  
	  term[l] += al;
	  sum += al;
	}
    }
  
  for(l = 0; l < numStates; l++)        
    ancestral[l] = term[l] / sum;       
   
  free(term);
}


static void newviewFlexGamma_Ancestral(int tipCase,
				       double *x1, double *x2, double *extEV, double *tipVector,
				       unsigned char *tipX1, unsigned char *tipX2,
				       int n, double *left, double *right,
				       const int numStates, double *diagptable, double *sumBuffer)
{
  double  
    *v,
    *x3 = (double*)malloc(sizeof(double) * 4 * numStates);
  double x1px2;
  int  i, j, l, k, scale;
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
		v =  &(x3[numStates * k]);

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
	    
	    ancestralGamma(x3, sumBuffer, diagptable, i, numStates, gammaStates);	   
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
		v =  &(x3[numStates * k]);

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
	   
	    v = x3;
	    scale = 1;
	    for(l = 0; scale && (l < gammaStates); l++)
	      scale = (ABS(v[l]) <  minlikelihood);

	    if(scale)
	      {
		for(l = 0; l < gammaStates; l++)
		  v[l] *= twotothe256;	    
	      }

	    ancestralGamma(x3, sumBuffer, diagptable, i, numStates, gammaStates);	    	    
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
	     v =  &(x3[numStates * k]);

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
	 
	 v = x3;
	 scale = 1;
	 for(l = 0; scale && (l < gammaStates); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));

	 if (scale)
	   {
	     for(l = 0; l < gammaStates; l++)
	       v[l] *= twotothe256;	     	    
	   }
	 
	 ancestralGamma(x3, sumBuffer, diagptable, i, numStates, gammaStates);	 	 
       }
      break;
    default:
      assert(0);
    }

 

  free(x3);

}
void newviewIterativeAncestral(tree *tr)
{
  traversalInfo 
    *ti   = tr->td[0].ti;
  
  int 
    i, 
    model;
  
  assert(!tr->useGappedImplementation);
  assert(!tr->saveMemory);
  assert(!tr->estimatePerSiteAA);
  
  for(i = 1; i < tr->td[0].count; i++)
    {
      traversalInfo 
	*tInfo = &ti[i];

      for(model = 0; model < tr->NumberOfModels; model++)
	{	  	    
	  double
	    *x1_start = (double*)NULL,
	    *x2_start = (double*)NULL,	   
	    *left     = (double*)NULL,
	    *right    = (double*)NULL,	   	   	       
	    qz, 
	    rz;
	      	  	  	  	  
	  unsigned char
	    *tipX1 = (unsigned char *)NULL,
	    *tipX2 = (unsigned char *)NULL;
	  
	  size_t
	    rateHet,
	    states = (size_t)tr->partitionData[model].states,
	    width = tr->partitionData[model].width;
	    	       
	  if(tr->rateHetModel == CAT)
	    rateHet = 1;
	  else
	    rateHet = 4;	  	  
	 	  
	  switch(tInfo->tipCase)
	    {
	    case TIP_TIP:		  
	      tipX1    = tr->partitionData[model].yVector[tInfo->qNumber];
	      tipX2    = tr->partitionData[model].yVector[tInfo->rNumber];		 		 		 	      	     	      	      
	      break;
	    case TIP_INNER:		 
	      tipX1    = tr->partitionData[model].yVector[tInfo->qNumber];		 
	      x2_start = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];	            
	      break;
	    case INNER_INNER:		 		 		    
	      x1_start = tr->partitionData[model].xVector[tInfo->qNumber - tr->mxtips - 1];
	      x2_start = tr->partitionData[model].xVector[tInfo->rNumber - tr->mxtips - 1];	      	      	  
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
	  
	  switch(tr->rateHetModel)
	    {
	    case CAT:
	      {	
		double
		  *diagptable = (double*)malloc_aligned(tr->partitionData[model].numberOfCategories * states * states * sizeof(double));
		
		makeP_Flex(qz, rz, tr->partitionData[model].perSiteRates,
			   tr->partitionData[model].EI,
			   tr->partitionData[model].EIGN,
			   tr->partitionData[model].numberOfCategories, left, right, states);
		

		makeP_Flex_Ancestral(tr->partitionData[model].perSiteRates,
				     tr->partitionData[model].EI,
				     tr->partitionData[model].EIGN,
				     tr->partitionData[model].numberOfCategories, diagptable, states);
				     

		newviewFlexCat_Ancestral(tInfo->tipCase,  tr->partitionData[model].EV, tr->partitionData[model].rateCategory,
					 x1_start, x2_start, tr->partitionData[model].tipVector,
					 tipX1, tipX2, width, left, right, states, diagptable, 
					 tr->partitionData[model].sumBuffer);

		free(diagptable);
	      }
	      break;
	    case GAMMA:
	    case GAMMA_I:
	      {	
		double
		  *diagptable = (double*)malloc_aligned(4 * states * states * sizeof(double));
		
		makeP_Flex(qz, rz, tr->partitionData[model].gammaRates,
			   tr->partitionData[model].EI,
			   tr->partitionData[model].EIGN,
			   4, left, right, states);
		
		makeP_Flex_Ancestral(tr->partitionData[model].gammaRates,
				     tr->partitionData[model].EI,
				     tr->partitionData[model].EIGN,
				     4, diagptable, states);
		

		newviewFlexGamma_Ancestral(tInfo->tipCase,
					   x1_start, x2_start,
					   tr->partitionData[model].EV,
					   tr->partitionData[model].tipVector,
					   tipX1, tipX2,
					   width, left, right, states, diagptable, tr->partitionData[model].sumBuffer);
		
		free(diagptable);
	      }
	      break;
	    default:
	      assert(0);
	    }	  		 			  			 
	}    
    }

}
static void traversalInfoAncestralRoot(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches)
{
  int 
    i;
  
  nodeptr 
    q = p,
    r = p->back;

  for(i = 0; i < numBranches; i++)
    {
      double 
	z = sqrt(q->z[i]);
       if(z < zmin) 
	 z = zmin;
       if(z > zmax)
	 z = zmax;
       
       z = log(z);
       
       ti[*counter].qz[i] = z;
       ti[*counter].rz[i] = z;     
    }

  if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
    {
      ti[*counter].tipCase = TIP_TIP;      
      ti[*counter].qNumber = q->number;
      ti[*counter].rNumber = r->number;
     
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
	  
	  ti[*counter].tipCase = TIP_INNER;	   
	  ti[*counter].qNumber = q->number;
	  ti[*counter].rNumber = r->number;	  	  
	  
	  *counter = *counter + 1;
	}
      else
	{	   
	  ti[*counter].tipCase = INNER_INNER;	    
	  ti[*counter].qNumber = q->number;
	  ti[*counter].rNumber = r->number;	  	 
	  
	  *counter = *counter + 1;
	}
    }
}


void newviewGenericAncestral (tree *tr, nodeptr p, boolean atRoot)
{  
  if(atRoot)
    {
      assert(!tr->multiGene);
      tr->td[0].count = 1;
      traversalInfoAncestralRoot(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
      
      if(tr->td[0].count > 1)
	{
#ifdef _USE_PTHREADS
	  masterBarrier(THREAD_NEWVIEW_ANCESTRAL, tr);
#else
	  newviewIterativeAncestral(tr);
#endif
	}
    }   
  else
    {
      if(isTip(p->number, tr->mxtips))
	return;
      
      if(tr->multiGene)       
	assert(0);     
      else
	{
	  tr->td[0].count = 1;
	  computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->mxtips, tr->numBranches);
	  
	  if(tr->td[0].count > 1)
	    {
#ifdef _USE_PTHREADS
	      masterBarrier(THREAD_NEWVIEW_ANCESTRAL, tr);
#else
	      newviewIterativeAncestral(tr);
#endif
	    }
	}
    }
}

typedef struct {
  double *probs;
  char c;
  int states;
} ancestralState;


static void computeAncestralRec(tree *tr, nodeptr p, int *counter, FILE *probsFile, FILE *statesFile, boolean atRoot)
{
#ifdef _USE_PTHREADS
  size_t 
    accumulatedOffset = 0;
#endif

  int 
    model,
    globalIndex = 0;
  
  ancestralState 
    *a = (ancestralState *)malloc(sizeof(ancestralState) * tr->cdta->endsite),
    *unsortedA = (ancestralState *)malloc(sizeof(ancestralState) * tr->rdta->sites);
  
  if(!atRoot)
    {
      if(isTip(p->number, tr->mxtips))
	return;
  
      computeAncestralRec(tr, p->next->back,       counter, probsFile, statesFile, atRoot);
      computeAncestralRec(tr, p->next->next->back, counter, probsFile, statesFile, atRoot);

      newviewGeneric(tr, p);
    }

  newviewGenericAncestral(tr, p, atRoot);

#ifdef _USE_PTHREADS
  masterBarrier(THREAD_GATHER_ANCESTRAL, tr);
#endif

  if(atRoot)
    {
      fprintf(probsFile, "ROOT\n");
      fprintf(statesFile, "ROOT ");
    }
  else
    {
      fprintf(probsFile, "%d\n", p->number);
      fprintf(statesFile, "%d ", p->number);
    }

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      int	
	offset,
	i,
	width = tr->partitionData[model].upper - tr->partitionData[model].lower,	
	states = tr->partitionData[model].states;
#ifdef _USE_PTHREADS
      double
	*ancestral = tr->ancestralStates;
#else
      double 
	*ancestral = tr->partitionData[model].sumBuffer;
#endif

      if(tr->rateHetModel == CAT)
	offset = 1;
      else
	offset = 4;            

      for(i = 0; i < width; i++, globalIndex++)
	{
	  double
	    equal = 1.0 / (double)states,
	    max = -1.0;
	    
	  boolean
	    approximatelyEqual = TRUE;

	  int
	    max_l = -1,
	    l;
	  
	  char 
	    c;

	  a[globalIndex].states = states;
	  a[globalIndex].probs = (double *)malloc(sizeof(double) * states);
	  
	  for(l = 0; l < states; l++)
	    {
	      double 
		value = ancestral[offset * states * i + l];

	      if(value > max)
		{
		  max = value;
		  max_l = l;
		}
	      
	      approximatelyEqual = approximatelyEqual && (ABS(equal - value) < 0.000001);
	      
	      a[globalIndex].probs[l] = value;	      	      
	    }

	  
	  if(approximatelyEqual)
	    c = '?';	  
	  else
	    c = getStateCharacter(tr->partitionData[model].dataType, max_l);
	  
	  a[globalIndex].c = c;	  
	}

#ifdef _USE_PTHREADS
      accumulatedOffset += width * offset * states;
#endif            
    }

  {
    int 
      j, 
      k;
    
    for(j = 0; j < tr->cdta->endsite; j++)
      {
	for(k = 0; k < tr->rdta->sites; k++)	    
	  if(j == tr->patternPosition[k])		
	    {
	      int 
		sorted = j,
		unsorted = tr->columnPosition[k] - 1;
	      
	      unsortedA[unsorted].states = a[sorted].states;
	      unsortedA[unsorted].c = a[sorted].c;
	      unsortedA[unsorted].probs = (double*)malloc(sizeof(double) * unsortedA[unsorted].states);
	      memcpy(unsortedA[unsorted].probs,  a[sorted].probs, sizeof(double) * a[sorted].states);	      
	    }	   
	}  

    for(k = 0; k < tr->rdta->sites; k++)
      {
	for(j = 0; j < unsortedA[k].states; j++)
	  fprintf(probsFile, "%f ", unsortedA[k].probs[j]);
	fprintf(probsFile, "\n");
	fprintf(statesFile, "%c", unsortedA[k].c);
      }
    fprintf(probsFile, "\n");
    fprintf(statesFile, "\n");
  }


  *counter = *counter + 1;

  {
    int j;

    for(j = 0; j < tr->rdta->sites; j++)
      free(unsortedA[j].probs);
    for(j = 0; j < tr->cdta->endsite; j++)
      free(a[j].probs);
  }

  free(a);
  free(unsortedA);
}

static char *ancestralTreeRec(char *treestr, tree *tr, nodeptr p)
{        
  if(isTip(p->number, tr->rdta->numsp)) 
    {         
      sprintf(treestr, "%s", tr->nameList[p->number]);    
      
      while (*treestr) 
	treestr++;
    }
  else 
    {                    
      *treestr++ = '(';     
      treestr = ancestralTreeRec(treestr, tr, p->next->back);     
      *treestr++ = ',';
      treestr = ancestralTreeRec(treestr, tr, p->next->next->back);          
      *treestr++ = ')'; 
      sprintf(treestr, "%d", p->number);      
    }
        
  while (*treestr) 
    treestr++;
  
  return  treestr;
}


static char *ancestralTree(char *treestr, tree *tr)
{
  *treestr++ = '(';
  treestr = ancestralTreeRec(treestr, tr, tr->leftRootNode);
  *treestr++ = ',';
  treestr = ancestralTreeRec(treestr, tr, tr->rightRootNode);	 
  *treestr++ = ')';
  sprintf(treestr, "ROOT");
  while (*treestr) 
    treestr++;  
  
  *treestr++ = ';';
  
  *treestr++ = '\0';
  while (*treestr) treestr++;     
  return  treestr;
}

void computeAncestralStates(tree *tr, double referenceLikelihood, analdef *adef)
{
  int 
    counter = 0;
  
  char 
    treeFileName[2048],
    ancestralProbsFileName[2048],
    ancestralStatesFileName[2048];

  FILE
    *treeFile,
    *probsFile,
    *statesFile;

#ifdef _USE_PTHREADS
  tr->ancestralStates = (double*)malloc(getContiguousVectorLength(tr) * sizeof(double));
#endif

  /*  assert(!adef->compressPatterns);*/

  strcpy(ancestralProbsFileName,         workdir);
  strcpy(ancestralStatesFileName,         workdir);
  strcpy(treeFileName,         workdir);

  strcat(ancestralProbsFileName,         "RAxML_marginalAncestralProbabilities.");
  strcat(ancestralStatesFileName,        "RAxML_marginalAncestralStates.");
  strcat(treeFileName,                   "RAxML_nodeLabelledRootedTree.");

  strcat(ancestralProbsFileName,         run_id);
  strcat(ancestralStatesFileName,        run_id);
  strcat(treeFileName,                   run_id);
  
  probsFile = myfopen(ancestralProbsFileName,   "w");
  statesFile = myfopen(ancestralStatesFileName, "w");
  treeFile   = myfopen(treeFileName,            "w");

  assert(tr->leftRootNode == tr->rightRootNode->back);
  
  computeAncestralRec(tr, tr->leftRootNode, &counter, probsFile, statesFile, FALSE);
  computeAncestralRec(tr, tr->rightRootNode, &counter, probsFile, statesFile, FALSE);
  computeAncestralRec(tr, tr->rightRootNode, &counter, probsFile, statesFile, TRUE);
  
  evaluateGeneric(tr, tr->rightRootNode);

  if(fabs(tr->likelihood - referenceLikelihood) > 0.5)
    {
      printf("Something suspiciuous is going on with the marginal ancestral probability computations\n");
      assert(0);
    } 
  
  assert(counter == tr->mxtips - 1);
    
  ancestralTree(tr->tree_string, tr);

  fprintf(treeFile, "%s\n", tr->tree_string);

  fclose(probsFile);
  fclose(statesFile);
  fclose(treeFile);

  printBothOpen("Marginal Ancestral Probabilities written to file:\n%s\n\n", ancestralProbsFileName);
  printBothOpen("Ancestral Sequences based on Marginal Ancestral Probabilities written to file:\n%s\n\n", ancestralStatesFileName); 
  printBothOpen("Node-laballed ROOTED tree written to file:\n%s\n", treeFileName);
}


