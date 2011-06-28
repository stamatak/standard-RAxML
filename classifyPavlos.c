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
extern double masterTime;


#ifdef _USE_PTHREADS
extern volatile int NumberOfThreads;
extern volatile int NumberOfJobs;
#endif



static double getBranch(tree *tr, double *b, double *bb)
{
  double z = 0.0;

  if(!tr->multiBranch)
    {
      assert(tr->fracchange != -1.0);
      assert(b[0] == bb[0]);
      z = b[0];
      if (z < zmin) 
	z = zmin;      	 
      if(z > zmax)
	z = zmax;
      z = -log(z) * tr->fracchange;
      return z;	
    }
  else
    {       
      int i;
      double x = 0;
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  assert(b[i] == bb[i]);
	  assert(tr->partitionContributions[i] != -1.0);
	  assert(tr->fracchanges[i] != -1.0);
	  x = b[i];
	  if (x < zmin) 
	    x = zmin;      	 
	  if(x > zmax)
	    x = zmax;
	  x = -log(x) * tr->fracchanges[i];
	  z += x * tr->partitionContributions[i];
	}	
      
      return z;
    } 

}


static char *Tree2StringClassifyRec(char *treestr, tree *tr, nodeptr p, int *countBranches, 
				    int *inserts, boolean subtreeSummary, boolean *foundTheSubtree, int CUTOFF, boolean originalTree)
{        
  branchInfo *bInf = p->bInf;
  int        i, countQuery = 0;   

  *countBranches = *countBranches + 1;

  if(!subtreeSummary && !originalTree)
    {
      for(i = 0; i < tr->numberOfTipsForInsertion; i++)
	if(bInf->epa->countThem[i] > 0)
	  countQuery++;  
      
      if(countQuery > 0)
	{
	  *treestr++ = '(';   
	  for(i = 0; i <  tr->numberOfTipsForInsertion; i++)
	    {
	      if(bInf->epa->countThem[i] > 0)
		{	      
		  sprintf(treestr,"QUERY___%s___%d", tr->nameList[inserts[i]], bInf->epa->countThem[i]);
		  while (*treestr) treestr++;	
		  *treestr++ = ',';
		}
	    }
	}
    }

  if(isTip(p->number, tr->rdta->numsp)) 
    {
      char *nameptr = tr->nameList[p->number];  

      if(subtreeSummary)
	{
	  for(i = 0; i <  tr->numberOfTipsForInsertion; i++)
	    {
	      if(p->bInf->epa->executeThem[i] >= CUTOFF)
		foundTheSubtree[i] = TRUE;
	      else
		foundTheSubtree[i] = FALSE;
	    }	  
	}

        
      sprintf(treestr, "%s", nameptr);    
      while (*treestr) treestr++;
    }
  else 
    {                    
      *treestr++ = '(';     
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, 
				       countBranches, inserts,subtreeSummary,foundTheSubtree, CUTOFF, originalTree);     
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, 
				       countBranches, inserts, subtreeSummary,foundTheSubtree, CUTOFF, originalTree);          
      *treestr++ = ')';                    

      if(subtreeSummary)
	{  
	  for(i = 0; i < tr->numberOfTipsForInsertion; i++)
	    {
	      if(p->bInf->epa->executeThem[i] >= CUTOFF && 
		 p->next->back->bInf->epa->executeThem[i] < CUTOFF && 
		 p->next->next->back->bInf->epa->executeThem[i] < CUTOFF)
		foundTheSubtree[i] = TRUE;
	      else
		foundTheSubtree[i] = FALSE;
	    }
	}
    }
   
  if(!subtreeSummary && countQuery > 0)
    {
      sprintf(treestr, ":1.0[%s]", p->bInf->epa->branchLabel);
      while (*treestr) treestr++;
      *treestr++ = ')'; 
    }
    
  if(originalTree)
    sprintf(treestr, ":%8.20f[%s", p->bInf->epa->originalBranchLength, p->bInf->epa->branchLabel);
  else
    sprintf(treestr, ":1.0[%s", p->bInf->epa->branchLabel);
  while (*treestr) treestr++;

  assert(!(subtreeSummary == TRUE && originalTree == TRUE));
  assert(!(countQuery > 0 &&  originalTree == TRUE));

  if(subtreeSummary)
    {
      for(i = 0; i < tr->numberOfTipsForInsertion; i++)
	{
	  if(foundTheSubtree[i])
	    {
	      sprintf(treestr," %s", tr->nameList[inserts[i]]);
	      while (*treestr) treestr++;
	    }
	}
    }

  sprintf(treestr, "]");            	 	        
  while (*treestr) treestr++;

  return  treestr;
}


static int findRoot(nodeptr p,  int numberOfTipsForInsertion, int ntips)
{
  if(isTip(p->number, ntips))
    {
      int i;
      for(i = 0; i < numberOfTipsForInsertion; i++)
	if(p->bInf->epa->countThem[i] > 0)
	  return 0;

      return p->number;
    }
  else
    {
      int left;
      assert(p == p->next->next->next);

      left = findRoot(p->next->back, numberOfTipsForInsertion, ntips);

      if(left > 0)
	return left;
      else
	return findRoot(p->next->next->back, numberOfTipsForInsertion, ntips);
    }
}

static void calcSubtree(nodeptr p, int nTips, int numberOfTipsForInsertion)
{
  int i;

  if(isTip(p->number, nTips))    
    {
      for(i = 0; i < numberOfTipsForInsertion; i++)
	p->bInf->epa->executeThem[i] = p->bInf->epa->countThem[i];
    
      return;
    }
  else
    {
      nodeptr q;
      assert(p == p->next->next->next);
      
      q = p->next;

      while(q != p)
	{
	  calcSubtree(q->back, nTips, numberOfTipsForInsertion);	
	  q = q->next;
	}

     
      for(i = 0; i < numberOfTipsForInsertion; i++)
	p->bInf->epa->executeThem[i] = p->next->back->bInf->epa->executeThem[i] + 
	  p->next->next->back->bInf->epa->executeThem[i] + p->bInf->epa->countThem[i];
          
      return;
    }
}


static char *Tree2StringClassify(char *treestr, tree *tr, int *inserts, 
				 boolean subtreeSummary, int CUTOFF, branchInfo *bInf, int root, 
				 boolean  originalTree)
{
  nodeptr p;
  int countBranches = 0; 

  if(!subtreeSummary)
    {           
      p = tr->start->back;
      
      assert(!isTip(p->number, tr->mxtips));
      
      *treestr++ = '(';
      treestr = Tree2StringClassifyRec(treestr, tr, p->back, &countBranches, 
				       inserts, FALSE, (boolean *)NULL, -1, originalTree);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, &countBranches, 
				       inserts, FALSE, (boolean *)NULL, -1, originalTree);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, &countBranches, 
				       inserts, FALSE, (boolean *)NULL, -1, originalTree);
      *treestr++ = ')';
      *treestr++ = ';';
      
      assert(countBranches == 2 * tr->ntips - 3);
      
      *treestr++ = '\0';
      while (*treestr) treestr++;     
      return  treestr;
    }
  else
    {
      int i, j;
      boolean *foundTheSubtree = (boolean*)malloc(sizeof(boolean) * tr->numberOfTipsForInsertion);

      assert(root > 0);
      assert(originalTree == FALSE);
      
     
      for(i = 0; i < tr->numberOfBranches; i++)
	for(j = 0; j < tr->numberOfTipsForInsertion; j++)
	  bInf[i].epa->executeThem[j] = 0;
           
      p = tr->nodep[root];
	  
      calcSubtree(p, tr->mxtips, tr->numberOfTipsForInsertion);
      calcSubtree(p->back, tr->mxtips, tr->numberOfTipsForInsertion);     	      

      assert(isTip(p->number, tr->mxtips));           
      assert(!isTip(p->back->number, tr->mxtips));
	  
      *treestr++ = '(';
      *treestr++ = '(';	          
      sprintf(treestr, "%s", tr->nameList[p->number]);
      countBranches++;
      while (*treestr) treestr++;
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->back, &countBranches, 
				       inserts, TRUE, foundTheSubtree, CUTOFF, originalTree);      
      *treestr++ = ')';
      sprintf(treestr,"[ROOT]");
      while (*treestr) treestr++;
      *treestr++ = ')';
      *treestr++ = ';';
	  
      assert(countBranches == 2 * tr->ntips - 2);
      *treestr++ = '\0';
      while (*treestr) treestr++;
      free(foundTheSubtree);
      return  treestr;   
    }
}




static void markTips(nodeptr p, int *perm, int maxTips)
{
  if(isTip(p->number, maxTips))
    {
      perm[p->number] = 1;
      return;
    }
  else
    {
      nodeptr q = p->next;
      while(q != p)
	{
	  markTips(q->back, perm, maxTips);
	  q = q->next;
	}      
    }
}

static void testInsertThorough(tree *tr, nodeptr r, nodeptr q, boolean bootstrap)
{
  double 
    result,           
    qz[NUM_BRANCHES],
    z[NUM_BRANCHES];
  
  nodeptr  
    x = q->back;      
  
  int 
    *inserts = tr->inserts,    
    j;   

  for(j = 0; j < tr->numBranches; j++)    
    {
      qz[j] = q->z[j];
      z[j] = sqrt(qz[j]); 

      if(z[j] < zmin) 
	z[j] = zmin;
      
      if(z[j] > zmax)
	z[j] = zmax;
    }
  
  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {                
      if(q->bInf->epa->executeThem[j])
	{	 
	  nodeptr s = tr->nodep[inserts[j]];	  	 	    	  
	  
	  hookup(r->next,       q, z, tr->numBranches);
	  hookup(r->next->next, x, z, tr->numBranches);
	  hookupDefault(r, s, tr->numBranches);      		     
	  
	  newviewGeneric(tr, r);	     
	  
	  localSmooth(tr, r, smoothings);
	  
	  result = evaluateGeneric(tr, r);	 	       
	  	  
	  if(bootstrap)
	    tr->bInf[q->bInf->epa->branchNumber].epa->bootstrapBranches[j] = getBranch(tr, r->z, r->back->z);
	  else
	    tr->bInf[q->bInf->epa->branchNumber].epa->branches[j] = getBranch(tr, r->z, r->back->z);
	 
	  tr->bInf[q->bInf->epa->branchNumber].epa->likelihoods[j] = result;	      	 				              
	}
    }

  hookup(q, x, qz, tr->numBranches);
   
  r->next->next->back = r->next->back = (nodeptr) NULL; 
}



static void testInsertFast(tree *tr, nodeptr r, nodeptr q)
{
  double
    result,
    qz[NUM_BRANCHES], 
    z[NUM_BRANCHES];
  
  nodeptr  
    x = q->back;      
  
  int 
    i,
    *inserts = tr->inserts;
    	           

  assert(!tr->grouped);                    
  
  for(i = 0; i < tr->numBranches; i++)
    {	
      qz[i] = q->z[i];
      z[i] = sqrt(q->z[i]);      
      
      if(z[i] < zmin) 
	z[i] = zmin;
      if(z[i] > zmax)
	z[i] = zmax;
    }        
  
  hookup(r->next,       q, z, tr->numBranches);
  hookup(r->next->next, x, z, tr->numBranches);	                         
  
  newviewGeneric(tr, r);   
  
  for(i = 0; i < tr->numberOfTipsForInsertion; i++)
    {
      if(q->bInf->epa->executeThem[i])
	{	  	    
	  hookupDefault(r, tr->nodep[inserts[i]], tr->numBranches);
	  result = evaluateGeneric (tr, r);	     
	      
	  r->back = (nodeptr) NULL;
	  tr->nodep[inserts[i]]->back = (nodeptr) NULL;
	  	  
	  tr->bInf[q->bInf->epa->branchNumber].epa->likelihoods[i] = result;	  	     
	}    
    }
 
  hookup(q, x, qz, tr->numBranches);
  
  r->next->next->back = r->next->back = (nodeptr) NULL;
}

static void addTraverseRob(tree *tr, nodeptr r, nodeptr q,
			   boolean thorough, boolean bootstrap)
{       
  if(thorough)
    testInsertThorough(tr, r, q, bootstrap);
  else
    testInsertFast(tr, r, q);

  if(!isTip(q->number, tr->rdta->numsp))
    {   
      nodeptr a = q->next;

      while(a != q)
	{
	  addTraverseRob(tr, r, a->back, thorough, bootstrap);
	  a = a->next;
	}      
    }
} 

#ifdef _USE_PTHREADS

static size_t getContiguousVectorLength(tree *tr)
{
  size_t length = 0;
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      size_t realWidth = tr->partitionData[model].upper - tr->partitionData[model].lower;
      int states = tr->partitionData[model].states;

      length += (realWidth * (size_t)states * (size_t)(tr->discreteRateCategories));      	
    }

  return length;
}

static size_t getContiguousScalingLength(tree *tr)
{
  size_t length = 0;
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)    
    length += tr->partitionData[model].upper - tr->partitionData[model].lower;

  return length;
}

static void allocBranchX(tree *tr)
{
  int i;

  for(i = 0; i < tr->numberOfBranches; i++)
    {
      branchInfo *b = &(tr->bInf[i]);

      b->epa->left  = (double*)malloc_aligned(sizeof(double) * tr->contiguousVectorLength);
      b->epa->leftScaling = (int*)malloc(sizeof(int) * tr->contiguousScalingLength);

      b->epa->right = (double*)malloc_aligned(sizeof(double)  * tr->contiguousVectorLength);
      b->epa->rightScaling = (int*)malloc(sizeof(int) * tr->contiguousScalingLength);

      b->epa->leftParsimony = (parsimonyVector *)b->epa->left;
      b->epa->rightParsimony = (parsimonyVector *)b->epa->right;
    }

}

static void updateClassify(tree *tr, double *z, boolean *partitionSmoothed, boolean *partitionConverged, double *x1, double *x2, unsigned char *tipX1, unsigned char *tipX2, int tipCase)
{   
  int i;

  boolean smoothedPartitions[NUM_BRANCHES];

  double 
    newz[NUM_BRANCHES], 
    z0[NUM_BRANCHES];
  
  for(i = 0; i < tr->numBranches; i++)   
    z0[i] = z[i];          

  makenewzClassify(tr, newzpercycle, newz, z0, x1, x2, tipX1, tipX2, tipCase, partitionConverged); 

  for(i = 0; i < tr->numBranches; i++)    
    smoothedPartitions[i]  = partitionSmoothed[i];
 
  for(i = 0; i < tr->numBranches; i++)
    {         
      if(!partitionConverged[i])
	{		    
	  if(ABS(newz[i] - z0[i]) > deltaz)  	     
	    smoothedPartitions[i] = FALSE;       
	    	             
	  z[i] = newz[i];	 
	}
    }


  for(i = 0; i < tr->numBranches; i++)    
    partitionSmoothed[i]  = smoothedPartitions[i];  
}


static double localSmoothClassify (tree *tr, int maxtimes, int leftNodeNumber, int rightNodeNumber, int insertionNodeNumber, double *e1, double *e2, double *e3, 
				   branchInfo *b)
{ 
  int tipCase;
  
  boolean 
    allSmoothed,
    partitionSmoothed[NUM_BRANCHES],
    partitionConverged[NUM_BRANCHES];

  double 
    result,
    *x1 = (double*)NULL,
    *x2 = (double*)NULL,
    *x3 = (double*)NULL;
	  
  int
    i,
    *ex1 = (int*)NULL,
    *ex2 = (int*)NULL,
    *ex3 = (int*)NULL; 

  unsigned char 
    *tipX1 = (unsigned char *)NULL,
    *tipX2 = (unsigned char *)NULL;
  
  x3  = tr->temporaryVector;
  ex3 = tr->temporaryScaling;
    
  
  for(i = 0; i < tr->numBranches; i++)	
    partitionConverged[i] = FALSE;	

  while (--maxtimes >= 0) 
    {     
      for(i = 0; i < tr->numBranches; i++)	
	partitionSmoothed[i] = TRUE;     
	 	
      /* e3 */
      if(isTip(leftNodeNumber, tr->mxtips) && isTip(rightNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_TIP;
	  
	  tipX1 = tr->contiguousTips[leftNodeNumber];
	  tipX2 = tr->contiguousTips[rightNodeNumber];

	  newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);
	}
      else
	{
	  if (isTip(leftNodeNumber, tr->mxtips))
	    {
	      tipCase = TIP_INNER;

	      tipX1 = tr->contiguousTips[leftNodeNumber];
	      
	      x2  = b->epa->right;	     
	      ex2 = b->epa->rightScaling; 	  
	      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);
	    }
	  else 
	    {
	      if(isTip(rightNodeNumber, tr->mxtips))
		{
		  tipCase = TIP_INNER;

		  tipX1 = tr->contiguousTips[rightNodeNumber];
	  
		  x2  = b->epa->left;	 
		  ex2 = b->epa->leftScaling; 
		  newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1);
		}       
	      else
		{
		  tipCase = INNER_INNER;
		  
		  x1  = b->epa->left;
		  ex1 = b->epa->leftScaling;
		  
		  x2  = b->epa->right;
		  ex2 = b->epa->rightScaling;
		  newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);
		}
	    }
	}
	
     

      tipCase = TIP_INNER;
      
      x2    = tr->temporaryVector;
      tipX1 = tr->contiguousTips[insertionNodeNumber];

      updateClassify(tr, e3, partitionSmoothed, partitionConverged, x1, x2, tipX1, tipX2, tipCase);

      /* e1 **********************************************************/

      tipX1 = tr->contiguousTips[insertionNodeNumber];

      if(isTip(rightNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_TIP;

	  tipX2 = tr->contiguousTips[rightNodeNumber];
	}
      else
	{	 
	  tipCase = TIP_INNER;
		  
	  x2  = b->epa->right;
	  ex2 = b->epa->rightScaling;		  	
	}
	
      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e3, e2);

      if(isTip(leftNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_INNER;

	  tipX1 = tr->contiguousTips[leftNodeNumber];
	}
      else
	{
	  tipCase = INNER_INNER;

	  x1      =  b->epa->left;
	}

      updateClassify(tr, e1, partitionSmoothed, partitionConverged, x1, x3, tipX1, (unsigned char*)NULL, tipCase);

      /* e2 *********************************************************/

      tipX1 = tr->contiguousTips[insertionNodeNumber];

      if(isTip(leftNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_TIP;
	  
	  tipX2 = tr->contiguousTips[leftNodeNumber];
	}
      else
	{	 
	  tipCase = TIP_INNER;
		  
	  x2  = b->epa->left;
	  ex2 = b->epa->leftScaling;		  	
	}
	
      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e3, e1);

      if(isTip(rightNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_INNER;

	  tipX1 = tr->contiguousTips[rightNodeNumber];
	}
      else
	{
	  tipCase = INNER_INNER;

	  x1      =  b->epa->right;
	}

      updateClassify(tr, e2, partitionSmoothed, partitionConverged, x1, x3, tipX1, (unsigned char*)NULL, tipCase);


      allSmoothed = TRUE;
      for(i = 0; i < tr->numBranches; i++)
	{
	  if(partitionSmoothed[i] == FALSE)
	    allSmoothed = FALSE;
	  else
	    partitionConverged[i] = TRUE;
	}
      
      if(allSmoothed)
	break;
    }

  
  if(isTip(leftNodeNumber, tr->mxtips) && isTip(rightNodeNumber, tr->mxtips))
    {
      tipCase = TIP_TIP;

      tipX1 = tr->contiguousTips[leftNodeNumber];
      tipX2 = tr->contiguousTips[rightNodeNumber];

      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);
    }
  else
    {
      if (isTip(leftNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_INNER;

	  tipX1 = tr->contiguousTips[leftNodeNumber];	  

	  x2  = b->epa->right;	     
	  ex2 = b->epa->rightScaling; 	  
	  newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);
	}
      else 
	{
	  if(isTip(rightNodeNumber, tr->mxtips))
	    {
	      tipCase = TIP_INNER;

	      tipX1 = tr->contiguousTips[rightNodeNumber];
	      
	      x2  = b->epa->left;	 
	      ex2 = b->epa->leftScaling; 
	      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1);
	    }       
	  else
	    {
	      tipCase = INNER_INNER;
	      
	      x1  = b->epa->left;
	      ex1 = b->epa->leftScaling;
	      
	      x2  = b->epa->right;
	      ex2 = b->epa->rightScaling;
	      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);
	    }
	}
    }
  
     
  tipCase = TIP_INNER;
  
  x2    = tr->temporaryVector;
  tipX1 = tr->contiguousTips[insertionNodeNumber];    
 
  result = evalCL(tr, x2, ex3, tipX1, e3);

  return result;
}

/**
   void testInsertThoroughIterativeOld(tree *tr, int branchNumber, boolean bootstrap)
   {
   int 
   i;
   
   double 
   defaultArray[NUM_BRANCHES];
   
   boolean 
   executeAll[NUM_BRANCHES];
   
   for(i = 0; i < tr->numBranches; i++)
   {
   defaultArray[i] = defaultz;
   executeAll[i] = TRUE;
   }
   
   {
   double	    
   result, 
   lzqr, 
   lzqs, 
   lzrs, 
   lzsum, 
   lzq, 
   lzr, 
   lzs, 
   lzmax, 
   zqs[NUM_BRANCHES],
   zrs[NUM_BRANCHES],
   *qz,
   e1[NUM_BRANCHES], 
   e2[NUM_BRANCHES], 
   e3[NUM_BRANCHES];  
   
   int
   tipCase,	   
   insertions,
   leftNodeNumber,
   rightNodeNumber;
   
   
   branchInfo 
   *b = &(tr->bInf[branchNumber]);
   
   leftNodeNumber  = b->epa->leftNodeNumber;
   rightNodeNumber = b->epa->rightNodeNumber;
   
   for(insertions = 0; insertions < tr->numberOfTipsForInsertion; insertions++)
   { 
   if(b->epa->executeThem[insertions])
   {	     
   double 
   *x1 = (double*)NULL,
   *x2 = (double*)NULL,
   *x3 = (double*)NULL;
   
   int
   model = 0,
   *ex1 = (int*)NULL,
   *ex2 = (int*)NULL,
   *ex3 = (int*)NULL; 
   
   unsigned char 
   *tipX1 = (unsigned char *)NULL,
   *tipX2 = (unsigned char *)NULL;	  	  
   
   qz = b->epa->branchLengths;      
   
   tipX1 = tr->contiguousTips[tr->inserts[insertions]];
   
   if(isTip(leftNodeNumber, tr->mxtips))
   {
   tipX2 = tr->contiguousTips[leftNodeNumber];
   x2 = (double*)NULL;
   tipCase = TIP_TIP;
   }
   else
   {
   tipX2 = (unsigned char*)NULL;
   x2 = b->epa->left;
   tipCase = TIP_INNER;
   }
   
   makenewzClassify(tr, iterations, zqs, defaultArray, (double*)NULL, x2, tipX1, tipX2, tipCase, executeAll);  
   
   
   if(isTip(rightNodeNumber, tr->mxtips))
   { 
   x2 = (double*)NULL;
   
   tipX2 = tr->contiguousTips[rightNodeNumber];
   
   tipCase = TIP_TIP;
   }
   else
   {
   tipX2 = (unsigned char*)NULL;
   x2 = b->epa->right;
   tipCase = TIP_INNER;
   }
   
   makenewzClassify(tr, iterations, zrs, defaultArray, (double*)NULL, x2, tipX1, tipX2, tipCase, executeAll);      
   
   
   for(model = 0; model < tr->numBranches; model++)
   {
   lzqr = (qz[model] > zmin) ? log(qz[model]) : log(zmin);		  
   lzqs = (zqs[model] > zmin) ? log(zqs[model]) : log(zmin);
   lzrs = (zrs[model] > zmin) ? log(zrs[model]) : log(zmin);
   lzsum = 0.5 * (lzqr + lzqs + lzrs);
   
   lzq = lzsum - lzrs;
   lzr = lzsum - lzqs;
   lzs = lzsum - lzqr;
   lzmax = log(zmax);
   
   if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
   else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
   else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
   
   e1[model] = exp(lzq);
   e2[model] = exp(lzr);
   e3[model] = exp(lzs);
   }
   
   x3  = tr->temporaryVector;
   ex3 = tr->temporaryScaling;
   
   if(isTip(leftNodeNumber, tr->mxtips) && isTip(rightNodeNumber, tr->mxtips))
   {
   tipCase = TIP_TIP;
   
   tipX1 = tr->contiguousTips[leftNodeNumber];
   tipX2 = tr->contiguousTips[rightNodeNumber];
   
   newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);	
   }
   else
   {
   if (isTip(leftNodeNumber, tr->mxtips))
   {
   tipCase = TIP_INNER;
   
   tipX1 = tr->contiguousTips[leftNodeNumber];
   
   x2  = b->epa->right;	     
   ex2 = b->epa->rightScaling; 
   newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);	
   }
   else 
   {
   if(isTip(rightNodeNumber, tr->mxtips))
   {
   tipCase = TIP_INNER;
   
   tipX1 = tr->contiguousTips[rightNodeNumber];
   
   x2  = b->epa->left;	 
   ex2 = b->epa->leftScaling;
   newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1);	
   }       
   else
   {
   tipCase = INNER_INNER;
   
   x1  = b->epa->left;
   ex1 = b->epa->leftScaling;
   
   x2  = b->epa->right;
   ex2 = b->epa->rightScaling;
   newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);	
   }
   }	      
   }
   
   result = localSmoothClassify(tr, smoothings, leftNodeNumber, rightNodeNumber, tr->inserts[insertions], e1, e2, e3, b);
   
   b->epa->likelihoods[insertions] = result;	      			      
   
   if(bootstrap)
   b->epa->bootstrapBranches[insertions] = getBranch(tr, e3, e3);
   else
   b->epa->branches[insertions] = getBranch(tr, e3, e3);
   }	  
   }
   }
   }
*/




void testInsertThoroughIterative(tree *tr, int branchNumber, boolean bootstrap)
{    
  double	    
    result, 
    z,
    e1[NUM_BRANCHES],
    e2[NUM_BRANCHES],
    e3[NUM_BRANCHES],      
    *x3  = tr->temporaryVector;
     
  branchInfo 
    *b = &(tr->bInf[branchNumber]); 
  
  int   
    tipCase,
    model,
    insertions,
    leftNodeNumber = b->epa->leftNodeNumber,
    rightNodeNumber = b->epa->rightNodeNumber,
    *ex3 = tr->temporaryScaling;	         
   	     	  
  for(insertions = 0; insertions < tr->numberOfTipsForInsertion; insertions++)
    { 
      if(b->epa->executeThem[insertions])
	{
	  double 	    
	    *x1 = (double*)NULL,
	    *x2 = (double*)NULL;
	    
	  int	   
	    *ex1 = (int*)NULL,
	    *ex2 = (int*)NULL; 
	  
	  unsigned char 
	    *tipX1 = (unsigned char *)NULL,
	    *tipX2 = (unsigned char *)NULL;	     	    	  	  	    	  	  
	  
	  for(model = 0; model < tr->numBranches; model++)
	    {
	      z = sqrt(b->epa->branchLengths[model]);
	      if(z < zmin) 
		z = zmin;
	      if(z > zmax)
		z = zmax;

	      e1[model] = z;
	      e2[model] = z;
	      e3[model] = defaultz;
	    }
	      	  	  	    
	  if(isTip(leftNodeNumber, tr->mxtips) && isTip(rightNodeNumber, tr->mxtips))
	    {
	      tipCase = TIP_TIP;
	      
	      tipX1 = tr->contiguousTips[leftNodeNumber];
	      tipX2 = tr->contiguousTips[rightNodeNumber];
	      
	      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);	
	    }
	  else
	    {
	      if (isTip(leftNodeNumber, tr->mxtips))
		{
		  tipCase = TIP_INNER;
		  
		  tipX1 = tr->contiguousTips[leftNodeNumber];
		  
		  x2  = b->epa->right;	     
		  ex2 = b->epa->rightScaling; 
		  newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);	
		}
	      else 
		{
		  if(isTip(rightNodeNumber, tr->mxtips))
		    {
		      tipCase = TIP_INNER;
		      
		      tipX1 = tr->contiguousTips[rightNodeNumber];
		      
		      x2  = b->epa->left;	 
		      ex2 = b->epa->leftScaling;
		      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1);	
		    }       
		  else
		    {
		      tipCase = INNER_INNER;
		      
		      x1  = b->epa->left;
		      ex1 = b->epa->leftScaling;
		      
		      x2  = b->epa->right;
		      ex2 = b->epa->rightScaling;
		      newviewClassifySpecial(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2);	
		    }
		}	      
	    }
	  
	  result = localSmoothClassify(tr, smoothings, leftNodeNumber, rightNodeNumber, tr->inserts[insertions], e1, e2, e3, b);
	      	    
	  b->epa->likelihoods[insertions] = result;	      			      
	  
	  if(bootstrap)
	    b->epa->bootstrapBranches[insertions] = getBranch(tr, e3, e3);
	  else
	    b->epa->branches[insertions] = getBranch(tr, e3, e3);
	}	  
    }
}







void addTraverseRobIterative(tree *tr, int branchNumber)
{      
  int 
    inserts;

  branchInfo 
    *b = &(tr->bInf[branchNumber]);

  double 
    result,
    z[NUM_BRANCHES],
    defaultArray[NUM_BRANCHES];       

      		 	    	        
  assert(!tr->useFastScaling);
     
  for(inserts = 0; inserts < tr->numBranches; inserts++) 
    {	      
      z[inserts] = sqrt(b->epa->branchLengths[inserts]);      
      defaultArray[inserts] = defaultz;
      
      if(z[inserts] < zmin) 
	z[inserts] = zmin;
      if(z[inserts] > zmax)
	z[inserts] = zmax;
    }        
                     
  newviewClassify(tr, b, z);   
        
  for(inserts = 0; inserts < tr->numberOfTipsForInsertion; inserts++) 
    { 	       		     
      result = evalCL(tr, tr->temporaryVector, tr->temporaryScaling, tr->contiguousTips[tr->inserts[inserts]], defaultArray);
	  	  
      b->epa->likelihoods[inserts] = result;	  	 			        
    } 
} 





#endif






static void setupBranchMetaInfo(tree *tr, nodeptr p, int nTips, branchInfo *bInf)
{
  int 
    i,
    countBranches = tr->branchCounter;

  if(isTip(p->number, tr->mxtips))    
    {      
      p->bInf       = &bInf[countBranches];
      p->back->bInf = &bInf[countBranches];               	      

      bInf[countBranches].oP = p;
      bInf[countBranches].oQ = p->back;
      
      bInf[countBranches].epa->leftNodeNumber = p->number;
      bInf[countBranches].epa->rightNodeNumber = p->back->number;

      bInf[countBranches].epa->originalBranchLength = getBranch(tr, p->z, p->back->z);      
      bInf[countBranches].epa->branchNumber = countBranches;	
      
      for(i = 0; i < tr->numBranches; i++)
	bInf[countBranches].epa->branchLengths[i] = p->z[i];
      
#ifdef _USE_PTHREADS
      if(!p->back->x)
	newviewGeneric(tr, p->back);
      masterBarrier(THREAD_GATHER_LIKELIHOOD, tr);
#endif

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

      bInf[countBranches].epa->originalBranchLength = getBranch(tr, p->z, p->back->z);           
      bInf[countBranches].epa->branchNumber = countBranches;
      
      for(i = 0; i < tr->numBranches; i++)
	bInf[countBranches].epa->branchLengths[i] = p->z[i];


#ifdef _USE_PTHREADS
      if(!p->x)
	newviewGeneric(tr, p);            
	 
      if(!isTip(p->back->number, tr->mxtips))
	{
	  if(!p->back->x)
	    newviewGeneric(tr, p->back);	 
	}	      

      masterBarrier(THREAD_GATHER_LIKELIHOOD, tr);
#endif      

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
 





static void setupBranchInfo(tree *tr, nodeptr q)
{
  tr->branchCounter = 0;

  setupBranchMetaInfo(tr, q, tr->ntips, tr->bInf);
    
  assert(tr->branchCounter == tr->numberOfBranches);
}

static int infoCompare(const void *p1, const void *p2)
{
  info *rc1 = (info *)p1;
  info *rc2 = (info *)p2;

  double i = rc1->lh;
  double j = rc2->lh;

  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}

static void consolidateInfoMLHeuristics(tree *tr, boolean EPA)
{
  int 
    throwAwayStart,
    i, 
    j;

  info 
    *inf = (info*)malloc(sizeof(info) * tr->numberOfBranches);

  if(EPA)
    throwAwayStart = MAX(5, (int)(0.5 + (double)(tr->numberOfBranches) * tr->fastEPAthreshold));
  else
    throwAwayStart = MAX(5, (int)(tr->numberOfBranches/10));  

  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {     
      for(i = 0; i < tr->numberOfBranches; i++)
	{      
	  inf[i].lh = tr->bInf[i].epa->likelihoods[j];
	  inf[i].number = i;
	}
      
      qsort(inf, tr->numberOfBranches, sizeof(info), infoCompare);

      for(i = throwAwayStart; i < tr->numberOfBranches; i++)       
	tr->bInf[inf[i].number].epa->executeThem[j] = 0;	   	
    }
  
   for(i = 0; i < tr->numberOfBranches; i++)
     for(j = 0; j < tr->numberOfTipsForInsertion; j++)    
       tr->bInf[i].epa->likelihoods[j] = unlikely;     
  

  free(inf);
}


typedef struct
{
  int number;
  unsigned int score;
} infoMP;

static int infoCompareMP(const void *p1, const void *p2)
{
  infoMP *rc1 = (infoMP *)p1;
  infoMP *rc2 = (infoMP *)p2;

  double i = rc1->score;
  double j = rc2->score;

  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}

static void consolidateInfoMPHeuristics(tree *tr)
{
  int 
    throwAwayStart,
    i, 
    j;

  infoMP
    *inf = (infoMP*)malloc(sizeof(infoMP) * tr->numberOfBranches);
 
  throwAwayStart = MAX(5, (int)(0.5 + (double)(tr->numberOfBranches) * tr->fastEPAthreshold));
   

  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {     
      for(i = 0; i < tr->numberOfBranches; i++)
	{      
	  inf[i].score = tr->bInf[i].epa->parsimonyScores[j];
	  inf[i].number = i;
	}
      
      qsort(inf, tr->numberOfBranches, sizeof(infoMP), infoCompareMP);

      for(i = throwAwayStart; i < tr->numberOfBranches; i++)       
	tr->bInf[inf[i].number].epa->executeThem[j] = 0;	   	
    }
  
   for(i = 0; i < tr->numberOfBranches; i++)
     for(j = 0; j < tr->numberOfTipsForInsertion; j++)    
       tr->bInf[i].epa->likelihoods[j] = unlikely;     
   
  free(inf);
}

static void consolidateInfo(tree *tr)
{
  int i, j;

  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {
      double
	max = unlikely;

      int 
	max_index = -1;

      for(i = 0; i < tr->numberOfBranches; i++)
	{      
	  if(tr->bInf[i].epa->likelihoods[j] > max)
	    {
	      max = tr->bInf[i].epa->likelihoods[j];
	      max_index = i;
	    }
	}

      assert(max_index >= 0);

      tr->bInf[max_index].epa->countThem[j] = tr->bInf[max_index].epa->countThem[j] + 1;
    }
}

static void consolidateInfoBootstrap(tree *tr)
{
  int i, j;

  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {
      double
	max = unlikely;

      int 
	max_index = -1;

      for(i = 0; i < tr->numberOfBranches; i++)
	{      
	  if(tr->bInf[i].epa->likelihoods[j] > max)
	    {
	      max = tr->bInf[i].epa->likelihoods[j];
	      max_index = i;
	    }
	}

      assert(max_index >= 0);

      tr->bInf[max_index].epa->countThem[j] = tr->bInf[max_index].epa->countThem[j] + 1;
      tr->bInf[max_index].epa->branches[j] += tr->bInf[max_index].epa->bootstrapBranches[j];
    }
}


/*#define _ERICK */

#ifdef _ERICK

/* function to lexicographically compare strings in C using qsort */

static int compareStrings(const void *p1, const void *p2)
{
  char **rc1 = (char **)p1;
  char **rc2 = (char **)p2;  
  
  return strcmp(*rc1, *rc2);
}


/* 
   fix order of subtree by looking for the taxon with the minimum rank according to
   the lexicographic order stored in order.
*/
   
   


static int fixOrder(tree *tr, nodeptr p, int *order)
{


  if(isTip(p->number, tr->mxtips))
    {
      /* 
	 if node p is a tip just return the rank according to the lexicographic sorting 
      */

      assert(order[p->number] > 0 && order[p->number] <= tr->mxtips);

      return (order[p->number]);
    }
  else
    {
      /* 
	 p is an inner node, so first get the minimum rank 
	 of the taxa contained in the left and right subtree,
	 given by p->next->back and 
	 p->next->next->back
      */

      int
	leftRank = fixOrder(tr, p->next->back, order),
	rightRank = fixOrder(tr, p->next->next->back, order);

      /* 
	 now if the minimum lexikographic rank of a taxon on the 
	 left subtree is larger than the minimum lexicographic rank 
	 of a taxon in the right subtree we need to exchange them, i.e.,
	 the right subtree becomes the left subtree and the 
	 left subtree becomes the right one 
      */

      if(leftRank >= rightRank)
	{
	  nodeptr 	    
	    left = p->next->back,
	    right = p->next->next->back;	 	     	  
	  
	  hookup(p->next, right, right->z, tr->numBranches); 
	  hookup(p->next->next, left, left->z, tr->numBranches); 	 
	}

      /* 
	 now just return the minimum of the left and right subtree ranks 
	 for the recursion
      */

      return (MIN(leftRank, rightRank));
    }

}


static nodeptr findCanonicTip(tree *tr, int *perm)
{
  char 
    **names = (char **)malloc(sizeof(char *) * ((size_t)tr->ntips));

  int
    *order = (int *)calloc(((size_t)(tr->mxtips + 1)), sizeof(int)),
    lookup,
    i,
    j;
  

  /*
    add those taxa in the alignment that form part of the reference 
    tree to the list that shall be sorted 
  */

  for(i = 1, j = 0; i <= tr->mxtips; i++)
    {
      if(perm[i] == 1)
	{
	  int 
	    len = strlen(tr->nameList[i]);
	  
	  names[j] = (char*)malloc(len * sizeof(char));

	  strcpy(names[j], tr->nameList[i]);

	  j++;
	}
    }
  
  /*
    sort the list lexicographically
  */

  qsort(names, tr->ntips, sizeof(char *), compareStrings);


  /* 
     now store the lexicographic order/rank for every taxon 
     that is in the reference tree in an array.
     if a taxon is not in the ref tree the respective entry 
     in array order will be zero 
  */

  for(i = 0; i < tr->ntips; i++)
    {
      lookup =  lookupWord(names[i], tr->nameHash);

      assert(lookup > 0 && lookup <= tr->mxtips);
      
      order[lookup] = i;
    }


  /*
    now set the canonic starting taxon from where we
    start traversing the tree to the lexicographically 
    smallest tip
  */

  lookup = lookupWord(names[0], tr->nameHash);

  /*
    make sure that the starting taxon is part of the reference 
    tree to feel better 
  */

  assert(perm[lookup] == 1);
  
  printBothOpen("Re-rooting to %s %s\n", names[0], tr->nameList[lookup]);
  
  /* 
     free string arrays
  */

  for(i = 0; i < tr->ntips; i++)
    free(names[i]);

  free(names);

  /*
    call recursive function to fix ordering of subtrees such 
    that the left subtree always contains the lexicographically 
    smaller taxon 
  */

  fixOrder(tr, tr->nodep[lookup]->back, order);
  
  /* 
     free rank/order array
  */

  free(order);
  
  /* return canonic starting node */

  return (tr->nodep[lookup]);
}

#endif


#ifdef _PAVLOS

static void printPerBranchReadAlignments(tree *tr)
{
  int 
    i, 
    j;
  
  for(j = 0; j < tr->numberOfBranches; j++) 
    {
      int 
	readCounter = 0;
      
      for(i = 0; i < tr->numberOfTipsForInsertion; i++)        	
	if(tr->bInf[j].epa->countThem[i] > 0)	    
	  readCounter++;

       if(readCounter > 0)
	{
	  int l, w;

	  char 
	    alignmentFileName[2048] = "",
	    buf[64] = "";

	  FILE 
	    *af;

	  strcat(alignmentFileName, workdir);
	  strcat(alignmentFileName, "RAxML_BranchAlignment.I");

	  sprintf(buf, "%d", j);

	  strcat(alignmentFileName, buf);
	  
	  af = myfopen(alignmentFileName, "wb");
	  

	  fprintf(af, "%d %d\n", readCounter, tr->rdta->sites);
	  
	  for(i = 0; i < tr->numberOfTipsForInsertion; i++)        	
	    {
	      if(tr->bInf[j].epa->countThem[i] > 0)	    	      	
		{		 		 
		  unsigned char *tip   =  tr->yVector[tr->inserts[i]];
		   
		  fprintf(af, "%s ", tr->nameList[tr->inserts[i]]);

		  for(l = 0; l < tr->cdta->endsite; l++)
		    {
		      for(w = 0; w < tr->cdta->aliaswgt[l]; w++)
			fprintf(af, "%c", getInverseMeaning(tr->dataVector[l], tip[l]));	      
		    }

		  fprintf(af, "\n");
		}	    		  
	    }
	  fclose(af);

	  printBothOpen("Branch read alignment for branch %d written to file %s\n", j, alignmentFileName);
	}
    }	  
}


#endif

void classifyML(tree *tr, analdef *adef)
{
  int
    root = -1,
    i, 
    j,  
    *perm;    
  
#ifdef _ERICK
  nodeptr
    canonicRoot;
#endif

  nodeptr     
    r, 
    q;    
     
  boolean     
    thorough = adef->thoroughInsertion;

  char 
    labelledTreeFileName[1024],
    originalLabelledTreeFileName[1024],
    classificationFileName[1024];

  FILE 
    *treeFile, 
    *classificationFile;

  adef->outgroup = FALSE;
  tr->doCutoff   = FALSE;

  assert(adef->restart);
  
  tr->numberOfBranches = 2 * tr->ntips - 3;

  printBothOpen("\nRAxML Evolutionary Placement Algorithm\n");

  printBothOpen("Thorough  Insertion Method with branch length optimization: %s\n", thorough?"YES":"NO");  

  evaluateGenericInitrav(tr, tr->start); 
  
  modOpt(tr, adef, TRUE, 1.0, FALSE);

  perm    = (int *)calloc(tr->mxtips + 1, sizeof(int));
  tr->inserts = (int *)calloc(tr->mxtips, sizeof(int));

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

#ifdef _ERICK  
  canonicRoot = findCanonicTip(tr, perm);
  tr->start = canonicRoot;
#endif

  free(perm);
  
  printBothOpen("RAxML will place %d Query Sequences into the %d branches of the reference tree with %d taxa\n\n",  tr->numberOfTipsForInsertion, (2 * tr->ntips - 3), tr->ntips);

  if(adef->rapidBoot)
    printBothOpen("RAxML will execute %d Bootstrap replicates to infer classification support\n\n",  adef->multipleRuns);

  assert(tr->numberOfTipsForInsertion == (tr->mxtips - tr->ntips));      

  tr->bInf              = (branchInfo*)malloc(tr->numberOfBranches * sizeof(branchInfo)); 

  for(i = 0; i < tr->numberOfBranches; i++)
    {      
      tr->bInf[i].epa = (epaBranchData*)malloc(sizeof(epaBranchData));

      tr->bInf[i].epa->countThem   = (int*)calloc(tr->numberOfTipsForInsertion, sizeof(int));      
      
      tr->bInf[i].epa->executeThem = (int*)calloc(tr->numberOfTipsForInsertion, sizeof(int));
      for(j = 0; j < tr->numberOfTipsForInsertion; j++)
	tr->bInf[i].epa->executeThem[j] = 1;

      tr->bInf[i].epa->parsimonyScores = (unsigned int*)calloc(tr->numberOfTipsForInsertion, sizeof(unsigned int));

      tr->bInf[i].epa->branches    = (double*)calloc(tr->numberOfTipsForInsertion, sizeof(double));     
      tr->bInf[i].epa->bootstrapBranches    = (double*)calloc(tr->numberOfTipsForInsertion, sizeof(double));     
      
      tr->bInf[i].epa->likelihoods = (double*)calloc(tr->numberOfTipsForInsertion, sizeof(double));      
      tr->bInf[i].epa->branchNumber = i;     
      sprintf(tr->bInf[i].epa->branchLabel, "I%d", i);     
    } 

  r = tr->nodep[(tr->nextnode)++]; 
    
#ifdef _ERICK
  q = canonicRoot;
#else
  q = findAnyTip(tr->start, tr->rdta->numsp);
#endif

  assert(isTip(q->number, tr->rdta->numsp));
  assert(!isTip(q->back->number, tr->rdta->numsp));
	 
  q = q->back;   
     
         
  
#ifdef _USE_PTHREADS 
  tr->contiguousVectorLength = getContiguousVectorLength(tr);
  tr->contiguousScalingLength = getContiguousScalingLength(tr);
  allocBranchX(tr);
  masterBarrier(THREAD_INIT_EPA, tr); 
#endif 
  
  setupBranchInfo(tr, q);   

 
  
  if(tr->fastEPA_MP || tr->fastEPA_ML)
    {	 
      int heuristicInsertions =  MAX(5, (int)(0.5 + (double)(tr->numberOfBranches) * tr->fastEPAthreshold));	  	 
      
      thorough = FALSE;   	         	 
	        
      if(tr->fastEPA_MP)
	{ 
	  printBothOpen("Searching for %d out of %d most promising branches with MP heuristics\n", heuristicInsertions, tr->numberOfBranches);	 
   
	  makeParsimonyInsertions(tr, q, r);
	  
	  evaluateGenericInitrav(tr, q->back);

#ifdef _USE_PTHREADS
	  setupBranchInfo(tr, q);    	
#endif	  

	  consolidateInfoMPHeuristics(tr);
	}     
	  
      if(tr->fastEPA_ML)
	{
	  printBothOpen("Searching for %d out of %d most promising branches with ML heuristics\n", heuristicInsertions, tr->numberOfBranches);	      

#ifdef _USE_PTHREADS	 
	  NumberOfJobs = tr->numberOfBranches;
	  masterBarrier(THREAD_INSERT_CLASSIFY, tr);               			 
#else  		
	  addTraverseRob(tr, r, q, thorough, FALSE);
#endif

	  consolidateInfoMLHeuristics(tr, TRUE);
	}           

     
      
      thorough = TRUE;
            
      

#ifdef _USE_PTHREADS
      NumberOfJobs = tr->numberOfBranches;
      masterBarrier(THREAD_INSERT_CLASSIFY_THOROUGH, tr);	      
#else     
      addTraverseRob(tr, r, q, thorough, FALSE);
#endif
      consolidateInfo(tr);              
    }
  else
    {            
	               	                 
#ifdef _USE_PTHREADS
      NumberOfJobs = tr->numberOfBranches;
      if(thorough)   
	masterBarrier(THREAD_INSERT_CLASSIFY_THOROUGH, tr);	           
      else
	masterBarrier(THREAD_INSERT_CLASSIFY, tr);
#else                
      addTraverseRob(tr, r, q, thorough, FALSE);	         	       
#endif
      
      if(adef->rapidBoot)
	{
	  int 
	    replicates,
	    *originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int)),
	    *originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));
       	 
	  memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
	  memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite); 	 	  	  	 		  	 

	  consolidateInfoMLHeuristics(tr, FALSE);
	  
	  for(i = 0; i < tr->numberOfBranches; i++)
	    for(j = 0; j < tr->numberOfTipsForInsertion; j++)         
	      tr->bInf[i].epa->branches[j] = 0.0;   
	  
	  for(replicates = 0; replicates < adef->multipleRuns; replicates++)
	    {	     	      	      
	      computeNextReplicate(tr, &adef->rapidBoot, originalRateCategories, originalInvariant, TRUE, TRUE);
	      	 
	      resetBranches(tr);
	      evaluateGenericInitrav(tr, q->back);	  
	      treeEvaluate(tr, 1);
	      
#ifdef _USE_PTHREADS
	      NumberOfJobs = tr->numberOfBranches;
	      masterBarrier(THREAD_CONTIGUOUS_REPLICATE, tr); 
	      setupBranchInfo(tr, q); 
	      if(thorough)	   	    
		masterBarrier(THREAD_INSERT_CLASSIFY_THOROUGH_BS, tr);	    
	      else
		masterBarrier(THREAD_INSERT_CLASSIFY, tr);	     	     	  
#else
	      addTraverseRob(tr, r, q, thorough, TRUE); 
#endif  
	     	      
	      consolidateInfoBootstrap(tr);	      
	    }	  
	}
      else
	consolidateInfo(tr);	
    }
      
  printBothOpen("Overall Classification time: %f\n\n", gettime() - masterTime);			               	


#ifdef _PAVLOS
  assert(adef->compressPatterns  == FALSE);
  printPerBranchReadAlignments(tr);
#endif

  strcpy(labelledTreeFileName,         workdir);
  strcpy(originalLabelledTreeFileName, workdir);
  strcpy(classificationFileName,       workdir);
   
  strcat(labelledTreeFileName,         "RAxML_labelledTree.");
  strcat(originalLabelledTreeFileName, "RAxML_originalLabelledTree.");
  strcat(classificationFileName,       "RAxML_classification.");
  
  
  strcat(labelledTreeFileName,         run_id);
  strcat(originalLabelledTreeFileName, run_id);
  strcat(classificationFileName,       run_id);
 
  free(tr->tree_string);
  tr->treeStringLength *= 16;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char));
 
 
  treeFile = myfopen(labelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, FALSE, -1, (branchInfo*)NULL, -1, FALSE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);
  
  treeFile = myfopen(originalLabelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, FALSE, -1, (branchInfo*)NULL, -1, TRUE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);
    
  if(adef->rapidBoot)
    {
      root = findRoot(tr->start,  tr->numberOfTipsForInsertion, tr->mxtips);           
      
      if(root == 0)
	root = findRoot(tr->start->back,  tr->numberOfTipsForInsertion, tr->mxtips);
      
      if(root == 0)    
	printf("We are in deep shit, the tree can't be rooted\n");	      
      else
	{
	  char extendedTreeFileName[1024];
	  int cutoff, k;	     
	  
	  for(k = 20; k < 100; k += 5)
	    {
	      cutoff = (int)(((double)adef->multipleRuns * (double)k / 100.0) + 0.5);
	      sprintf(extendedTreeFileName, "%s.%d", labelledTreeFileName, k);	     
	      
	      treeFile = myfopen(extendedTreeFileName, "wb");	  
	      Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, cutoff, tr->bInf, root, FALSE);	     
	      fprintf(treeFile, "%s\n", tr->tree_string);
	      fclose(treeFile);
	      
	      printBothOpen("Least common ancestor file for cutoff %d (%d) written to file %s\n", cutoff, k, extendedTreeFileName);	      
	    }
	}
    }

  classificationFile = myfopen(classificationFileName, "wb");
  
  for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
    for(j = 0; j < tr->numberOfBranches; j++) 
      {       
	if(tr->bInf[j].epa->countThem[i] > 0)	    
	  {
	    if(thorough)
	      fprintf(classificationFile, "%s I%d %d %8.20f\n", tr->nameList[tr->inserts[i]], j, tr->bInf[j].epa->countThem[i], 
		      tr->bInf[j].epa->branches[i] / (double)(tr->bInf[j].epa->countThem[i]));
	    else
	      fprintf(classificationFile, "%s I%d %d\n", tr->nameList[tr->inserts[i]], j, tr->bInf[j].epa->countThem[i]);
	  }
      }
  
  fclose(classificationFile);  

  printBothOpen("\n\nLabelled reference tree including branch labels and query sequences written to file: %s\n\n", labelledTreeFileName); 
  printBothOpen("Labelled reference tree with branch labels (without query sequences) written to file: %s\n\n", originalLabelledTreeFileName); 
  printBothOpen("Classification result file written to file: %s\n\n", classificationFileName);
   

  if(!adef->rapidBoot)
    {
      info 
	*inf = (info*)malloc(sizeof(info) * tr->numberOfBranches);

      FILE 
	*likelihoodWeightsFile;
      FILE
	*entropyFile; /*pp 20110531 */

      char 
	likelihoodWeightsFileName[1024];
      char
	entropyFileName[1024]; /*pp 20110531 */

      strcpy(likelihoodWeightsFileName,       workdir);
      strcat(likelihoodWeightsFileName,       "RAxML_classificationLikelihoodWeights.");
      strcat(likelihoodWeightsFileName,       run_id);

      strcpy(entropyFileName, workdir); /*pp 20110531 */
      strcat(entropyFileName, "RAxML_entropy."); /*pp 20110531 */
      strcat(entropyFileName, run_id); /*pp 20110531 */

      likelihoodWeightsFile = myfopen(likelihoodWeightsFileName, "wb");
      entropyFile = myfopen(entropyFileName, "wb"); /*pp 20110531 */

      for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
	{
	  double
	    lmax = 0.0,
	    acc = 0.0;

	  int 
	    validEntries = 0;

	  for(j =  0; j < tr->numberOfBranches; j++) 
	    {
	      inf[j].lh = tr->bInf[j].epa->likelihoods[i];
	      inf[j].number = j;
	    }

	  qsort(inf, tr->numberOfBranches, sizeof(info), infoCompare);	 

	  for(j =  0; j < tr->numberOfBranches; j++) 
	    if(inf[j].lh == unlikely)	     
	      break;
	    else	     
	      validEntries++;	     	      

	  assert(validEntries > 0);

	  j = 0;
	  
	  lmax = inf[0].lh;


	  double entropy = 0.; /*pp 20110531 */

	  while(/*acc <= 0.95 &&*/ j < validEntries)	  /*pp 20110531 */
	    { 
	      int 
		k;
	      
	      double 
		all = 0.0,
		prob = 0.0;

	      for(k =  0; k < validEntries; k++) 	   
		all += exp(inf[k].lh - lmax);	     
	      
	      acc += (prob = (exp(inf[j].lh - lmax) / all));
	      if(prob > 0)
		entropy -= ( prob * log(prob) ); /*pp 20110531 */
	      /* else if prob == 0 then lim[p->0] p*logp = 0, i.e. does not contribute anyway to the entropy */
		
	      fprintf(likelihoodWeightsFile, "%s I%d %f %f\n", tr->nameList[tr->inserts[i]], inf[j].number, prob, acc);
	      
	      j++;
	    }
	  /* normalize entropy by dividing with the log(validEntries) which is the maximum Entropy possible */
	  fprintf(entropyFile, "%s\t%f\n", tr->nameList[tr->inserts[i]], entropy/log(validEntries) ); /*pp 20110531 */
	}      
      
      free(inf);

      fclose(likelihoodWeightsFile); 
      fclose(entropyFile); /*pp 20110531 */
      printBothOpen("Classification result file using likelihood wieghts written to file: %s\n\n", likelihoodWeightsFileName);
    }          

 
    
  exit(0);
}


