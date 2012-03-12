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

extern char **globalArgv;
extern int globalArgc;
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
				    int *inserts, boolean originalTree, boolean jointLabels, boolean likelihood)
{        
  branchInfo *bInf = p->bInf;
  int        i, countQuery = 0;   

  *countBranches = *countBranches + 1;

  

  if(!originalTree)
    {
      for(i = 0; i < tr->numberOfTipsForInsertion; i++)
	if(bInf->epa->countThem[i] > 0)
	  countQuery++;  
      
      if(countQuery > 0)
	{
	  int 
	    localCounter = 0;
	  
	  *treestr++ = '(';

	  if(countQuery > 1)
	    *treestr++ = '(';

	  for(i = 0; i <  tr->numberOfTipsForInsertion; i++)
	    {
	      if(bInf->epa->countThem[i] > 0)
		{	      		  
		  if(likelihood)
		    {
		      char 
			branchLength[128];
		      
		      sprintf(branchLength, "%f", bInf->epa->branches[i]);		  
		      sprintf(treestr,"QUERY___%s:%s", tr->nameList[inserts[i]], branchLength);
		    }
		  else
		    sprintf(treestr,"QUERY___%s", tr->nameList[inserts[i]]);
		  
		  while (*treestr) treestr++;
		  if(localCounter < countQuery - 1)
		    *treestr++ = ',';
		  localCounter++;
		}	      
	    }

	  if(countQuery > 1)
	    {
	      sprintf(treestr,"):0.0,");
	      while (*treestr) treestr++;
	    }
	  else
	    *treestr++ = ',';
	  
	}
    }

  if(isTip(p->number, tr->rdta->numsp)) 
    {
      char *nameptr = tr->nameList[p->number];  
        
      sprintf(treestr, "%s", nameptr);    
      while (*treestr) treestr++;
    }
  else 
    {                    
      *treestr++ = '(';     
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, 
				       countBranches, inserts, originalTree, jointLabels, likelihood);     
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, 
				       countBranches, inserts, originalTree, jointLabels, likelihood);          
      *treestr++ = ')';                         
    }
   
  if(countQuery > 0)
    {
      sprintf(treestr, ":%8.20f[%s]", 0.5 * p->bInf->epa->originalBranchLength, p->bInf->epa->branchLabel);
      while (*treestr) treestr++;
      *treestr++ = ')'; 
    }
    
  if(originalTree)
    {
      if(jointLabels)
	{
	  if(tr->wasRooted)
	    {	      
	      if(p == tr->leftRootNode)
		{
		  sprintf(treestr, ":%8.20f{%d", p->bInf->epa->originalBranchLength * 0.5, p->bInf->epa->jointLabel);  
		  assert(tr->rootLabel == p->bInf->epa->jointLabel);
		}
	      else
		{
		  if(p == tr->rightRootNode)
		    {
		      sprintf(treestr, ":%8.20f{%d", p->bInf->epa->originalBranchLength * 0.5, tr->numberOfBranches);  
		      assert(tr->rootLabel == p->bInf->epa->jointLabel);
		    }
		  else
		    sprintf(treestr, ":%8.20f{%d", p->bInf->epa->originalBranchLength, p->bInf->epa->jointLabel);       
		}
	    }
	  else
	    sprintf(treestr, ":%8.20f{%d", p->bInf->epa->originalBranchLength, p->bInf->epa->jointLabel);  
	}
      else
	sprintf(treestr, ":%8.20f[%s", p->bInf->epa->originalBranchLength, p->bInf->epa->branchLabel);	
    }
  else    
    {
      if(countQuery > 0)
	sprintf(treestr, ":%8.20f[%s", 0.5 * p->bInf->epa->originalBranchLength, p->bInf->epa->branchLabel);
      else
	sprintf(treestr, ":%8.20f[%s", p->bInf->epa->originalBranchLength, p->bInf->epa->branchLabel);
    }
     
  while (*treestr) treestr++;
  
  assert(!(countQuery > 0 &&  originalTree == TRUE));

  if(jointLabels) 
    sprintf(treestr, "}");   
  else
    sprintf(treestr, "]");            	 	        
  
  while (*treestr) treestr++;

  return  treestr;
}




char *Tree2StringClassify(char *treestr, tree *tr, int *inserts, 
			  boolean  originalTree, boolean jointLabels, boolean likelihood)
{
  nodeptr 
    p;
  
  int 
    countBranches = 0; 

      
  if(jointLabels && tr->wasRooted)
    { 
      assert(originalTree);
      
      *treestr++ = '(';
      treestr = Tree2StringClassifyRec(treestr, tr, tr->leftRootNode, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, tr->rightRootNode, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood);	 
      *treestr++ = ')';
      *treestr++ = ';';
      
      assert(countBranches == 2 * tr->ntips - 2);
      
      *treestr++ = '\0';
      while (*treestr) treestr++;     
      return  treestr;
    }
  else
    {
      if(jointLabels)
	p = tr->nodep[tr->mxtips + 1];
      else
	p = tr->start->back;
      
      assert(!isTip(p->number, tr->mxtips));
      
      *treestr++ = '(';
      treestr = Tree2StringClassifyRec(treestr, tr, p->back, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood);
      *treestr++ = ')';
      *treestr++ = ';';
      
      assert(countBranches == 2 * tr->ntips - 3);
      
      *treestr++ = '\0';
      while (*treestr) treestr++;     
      return  treestr;
    }
}




void markTips(nodeptr p, int *perm, int maxTips)
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


static boolean containsRoot(nodeptr p, tree *tr, int rootNumber)
{

  if(isTip(p->number, tr->mxtips))
    {
      if(p->number == rootNumber)
	return TRUE;
      else
	return FALSE;
    }
  else
    {
      if(p->number == rootNumber)
	return TRUE;
      else
	{
	  if(containsRoot(p->next->back, tr, rootNumber))	    
	    return TRUE;	    
	  else
	    {
	      if(containsRoot(p->next->next->back, tr, rootNumber))
		return TRUE;
	      else
		return FALSE;
	    }
	}

    }
}

static nodeptr findRootDirection(nodeptr p, tree *tr, int rootNumber)
{  
  if(containsRoot(p, tr, rootNumber))
    return p;
  
  if(containsRoot(p->back, tr, rootNumber))
    return (p->back);
    
  /* one of the two subtrees must contain the root */

  assert(0);
}


static void testInsertThorough(tree *tr, nodeptr r, nodeptr q)
{
  double 
    originalBranchLength = getBranch(tr, q->z, q->back->z),
    result,           
    qz[NUM_BRANCHES],
    z[NUM_BRANCHES];
  
  nodeptr      
    root = (nodeptr)NULL,
    x = q->back;      
  
  int 
    *inserts = tr->inserts,    
    j;     

  boolean 
    atRoot = FALSE;

  if(!tr->wasRooted)
    root = findRootDirection(q, tr, tr->mxtips + 1);
  else
    {
      if((q == tr->leftRootNode && x == tr->rightRootNode) ||
	 (x == tr->leftRootNode && q == tr->rightRootNode))
	atRoot = TRUE;
      else
	{
	  nodeptr 
	    r1 = findRootDirection(q, tr, tr->leftRootNode->number),
	    r2 = findRootDirection(q, tr, tr->rightRootNode->number);

	  assert(r1 == r2);

	  root = r1;
	}
    }
  
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
	  nodeptr 
	    s = tr->nodep[inserts[j]];	  	 	    	  

	  double 
	    ratio,
	    modifiedBranchLength,
	    distalLength;
	  
	  hookup(r->next,       q, z, tr->numBranches);
	  hookup(r->next->next, x, z, tr->numBranches);
	  hookupDefault(r, s, tr->numBranches);      		     
	  
	  newviewGeneric(tr, r);	     
	  
	  localSmooth(tr, r, smoothings);
	  
	  result = evaluateGeneric(tr, r);	 	       
	  	  
	 
	  tr->bInf[q->bInf->epa->branchNumber].epa->branches[j] = getBranch(tr, r->z, r->back->z);
	 
	  tr->bInf[q->bInf->epa->branchNumber].epa->likelihoods[j] = result;	 

	  modifiedBranchLength = getBranch(tr, q->z, q->back->z) + getBranch(tr, x->z, x->back->z);

	  ratio = originalBranchLength / modifiedBranchLength;

	  if(tr->wasRooted && atRoot)
	    {	     
	      /* always take distal length from left root node and then fix this later */

	      if(x == tr->leftRootNode)
		distalLength = getBranch(tr, x->z, x->back->z);
	      else
		{
		  assert(x == tr->rightRootNode);
		  distalLength = getBranch(tr, q->z, q->back->z);
		}
	    }
	  else
	    {
	      if(root == x)
		distalLength = getBranch(tr, x->z, x->back->z);
	      else
		{
		  assert(root == q);
		  distalLength = getBranch(tr, q->z, q->back->z);
		}	      	      
	    }

	  distalLength *= ratio;
          
	  assert(distalLength <= originalBranchLength);
	     
	  tr->bInf[q->bInf->epa->branchNumber].epa->distalBranches[j] = distalLength;	  
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


#ifndef _USE_PTHREADS

#ifdef _SPECIES_STUFF

static double testInsertSpecies(tree *tr, nodeptr r, nodeptr attachmentBranch, int insertNumber, boolean optimizeOtherBranch)
{
  double
    z[NUM_BRANCHES],
    result,
    qz[NUM_BRANCHES],
    zmins[NUM_BRANCHES];
  
  nodeptr  
    x, q;      
  
  int 
    i;
    	
  if(isTip(attachmentBranch->number, tr->mxtips))
    x = attachmentBranch;
  else
    x = attachmentBranch->back;

  q = x->back;

  assert(isTip(x->number, tr->mxtips));

  assert(!tr->grouped);                    
  
  for(i = 0; i < tr->numBranches; i++)    	
    {
      qz[i] = q->z[i];
      zmins[i] = zmax;
    }
     
  initrav(tr, q);

  hookup(r->next,       q, qz, tr->numBranches);
  hookup(r->next->next, x, zmins, tr->numBranches);	                         
  hookup(r, tr->nodep[insertNumber], zmins, tr->numBranches);

  printf("read name: %s\n", tr->nameList[insertNumber]);    

  newviewGeneric(tr, r);   
  
  if(optimizeOtherBranch)
    {     
      makenewzGeneric(tr, r->next, q, qz, 10, z, FALSE);
      hookup(q, r->next, z, tr->numBranches); 
      newviewGeneric(tr, r);
    }
 
  result = evaluateGeneric (tr, r);	     
	      
  r->back = (nodeptr) NULL;
  tr->nodep[insertNumber]->back = (nodeptr) NULL;
	  	   
  hookup(q, x, qz, tr->numBranches);
  
  r->next->next->back = r->next->back = (nodeptr) NULL;

  return result;
}

#endif

#endif

static void addTraverseRob(tree *tr, nodeptr r, nodeptr q,
			   boolean thorough)
{       
  if(thorough)
    testInsertThorough(tr, r, q);
  else
    testInsertFast(tr, r, q);

  if(!isTip(q->number, tr->rdta->numsp))
    {   
      nodeptr a = q->next;

      while(a != q)
	{
	  addTraverseRob(tr, r, a->back, thorough);
	  a = a->next;
	}      
    }
} 

#ifdef _USE_PTHREADS

size_t getContiguousVectorLength(tree *tr)
{
  size_t length = 0;
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      size_t 
	realWidth = tr->partitionData[model].upper - tr->partitionData[model].lower;
      
      int 
	states = tr->partitionData[model].states;

      length += (realWidth * (size_t)states * (size_t)(tr->discreteRateCategories));      	
    }

  return length;
}

static size_t getContiguousScalingLength(tree *tr)
{
  size_t 
    length = 0;
  
  int 
    model;

  for(model = 0; model < tr->NumberOfModels; model++)    
    length += tr->partitionData[model].upper - tr->partitionData[model].lower;

  return length;
}

static void allocBranchX(tree *tr)
{
  int 
    i = 0;

  for(i = 0; i < tr->numberOfBranches; i++)
    {
      branchInfo 
	*b = &(tr->bInf[i]);

      b->epa->left  = (double*)malloc_aligned(sizeof(double) * tr->contiguousVectorLength, 16);
      b->epa->leftScaling = (int*)malloc(sizeof(int) * tr->contiguousScalingLength);

      b->epa->right = (double*)malloc_aligned(sizeof(double)  * tr->contiguousVectorLength, 16);
      b->epa->rightScaling = (int*)malloc(sizeof(int) * tr->contiguousScalingLength);     
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





void testInsertThoroughIterative(tree *tr, int branchNumber)
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
  
  nodeptr      
    root = (nodeptr)NULL,
    x = b->oP,
    q = b->oQ;

  boolean 
    atRoot = FALSE;

  int   
    tipCase,
    model,
    insertions,
    leftNodeNumber = b->epa->leftNodeNumber,
    rightNodeNumber = b->epa->rightNodeNumber,
    *ex3 = tr->temporaryScaling;	         

  assert(x->number == leftNodeNumber);
  assert(q->number == rightNodeNumber);

  if(!tr->wasRooted)
    root = findRootDirection(q, tr, tr->mxtips + 1);
  else
    {
      if((q == tr->leftRootNode && x == tr->rightRootNode) ||
	 (x == tr->leftRootNode && q == tr->rightRootNode))
	atRoot = TRUE;
      else
	{
	  nodeptr 
	    r1 = findRootDirection(q, tr, tr->leftRootNode->number),
	    r2 = findRootDirection(q, tr, tr->rightRootNode->number);

	  assert(r1 == r2);

	  root = r1;
	}
    }
	     	  
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
	  
	   double 
	    ratio,
	    modifiedBranchLength,
	    distalLength;

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
	  	
	  b->epa->branches[insertions] = getBranch(tr, e3, e3);	  

	  modifiedBranchLength = getBranch(tr, e1, e1) + getBranch(tr, e2, e2);

	  ratio = b->epa->originalBranchLength / modifiedBranchLength;

	  if(tr->wasRooted && atRoot)
	    {	     
	      /* always take distal length from left root node and then fix this later */

	      if(x == tr->leftRootNode)
		distalLength = getBranch(tr, e1, e1);
	      else
		{
		  assert(x == tr->rightRootNode);
		  distalLength = getBranch(tr, e2, e2);
		}
	    }
	  else
	    {
	      if(root == x)
		distalLength = getBranch(tr, e1, e1);
	      else
		{
		  assert(root == q);
		  distalLength = getBranch(tr, e2, e2);
		}	      	      
	    }

	  distalLength *= ratio;
          
	  assert(distalLength <= b->epa->originalBranchLength);
	     
	  b->epa->distalBranches[insertions] = distalLength;	  	
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

static void consolidateInfoMLHeuristics(tree *tr, int throwAwayStart)
{
  int 
    i, 
    j;

  info 
    *inf = (info*)malloc(sizeof(info) * tr->numberOfBranches);

  assert(tr->useEpaHeuristics);

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

	  strcpy(alignmentFileName, workdir);
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
    i, 
    j,  
    *perm;    
  

  nodeptr     
    r, 
    q;    

  char
    entropyFileName[1024],
    jointFormatTreeFileName[1024],
    labelledTreeFileName[1024],
    originalLabelledTreeFileName[1024],
    classificationFileName[1024];

  FILE 
    *entropyFile,
    *treeFile, 
    *classificationFile;

  adef->outgroup = FALSE;
  tr->doCutoff   = FALSE;

  assert(adef->restart);
  
  tr->numberOfBranches = 2 * tr->ntips - 3;

  printBothOpen("\nRAxML Evolutionary Placement Algorithm\n"); 

  evaluateGenericInitrav(tr, tr->start); 
  
  modOpt(tr, adef, TRUE, 1.0, FALSE);

  printBothOpen("\nLikelihood of reference tree: %f\n\n", tr->likelihood);

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

  free(perm);
  
  printBothOpen("RAxML will place %d Query Sequences into the %d branches of the reference tree with %d taxa\n\n",  tr->numberOfTipsForInsertion, (2 * tr->ntips - 3), tr->ntips);  

  assert(tr->numberOfTipsForInsertion == (tr->mxtips - tr->ntips));      

  tr->bInf              = (branchInfo*)malloc(tr->numberOfBranches * sizeof(branchInfo)); 

  for(i = 0; i < tr->numberOfBranches; i++)
    {      
      tr->bInf[i].epa = (epaBranchData*)malloc(sizeof(epaBranchData));

      tr->bInf[i].epa->countThem   = (int*)calloc(tr->numberOfTipsForInsertion, sizeof(int));      
      
      tr->bInf[i].epa->executeThem = (int*)calloc(tr->numberOfTipsForInsertion, sizeof(int));
      
      for(j = 0; j < tr->numberOfTipsForInsertion; j++)
	tr->bInf[i].epa->executeThem[j] = 1;

      tr->bInf[i].epa->branches          = (double*)calloc(tr->numberOfTipsForInsertion, sizeof(double));   
      tr->bInf[i].epa->distalBranches    = (double*)calloc(tr->numberOfTipsForInsertion, sizeof(double)); 
         
      tr->bInf[i].epa->likelihoods = (double*)calloc(tr->numberOfTipsForInsertion, sizeof(double));      
      tr->bInf[i].epa->branchNumber = i;
      
      sprintf(tr->bInf[i].epa->branchLabel, "I%d", i);     
    } 

  r = tr->nodep[(tr->nextnode)++]; 
    

  q = findAnyTip(tr->start, tr->rdta->numsp);

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
  
  if(tr->useEpaHeuristics)
    {	 
      int 
	heuristicInsertions =  MAX(5, (int)(0.5 + (double)(tr->numberOfBranches) * tr->fastEPAthreshold));	  	    	         	 	             
	        
      printBothOpen("EPA heuristics: determining %d out of %d most promising insertion branches\n", heuristicInsertions, tr->numberOfBranches);	      

#ifdef _USE_PTHREADS	 
      NumberOfJobs = tr->numberOfBranches;
      masterBarrier(THREAD_INSERT_CLASSIFY, tr);               			 
#else  		
      addTraverseRob(tr, r, q, FALSE);
#endif
      
      consolidateInfoMLHeuristics(tr, heuristicInsertions);
    }           
           
#ifdef _USE_PTHREADS
  NumberOfJobs = tr->numberOfBranches;
  masterBarrier(THREAD_INSERT_CLASSIFY_THOROUGH, tr);	      
#else     
  addTraverseRob(tr, r, q, TRUE);
#endif
  consolidateInfo(tr);                
      
  printBothOpen("Overall Classification time: %f\n\n", gettime() - masterTime);			               	


#ifdef _PAVLOS
  assert(adef->compressPatterns  == FALSE);
  printPerBranchReadAlignments(tr);
#endif
  
  strcpy(entropyFileName,              workdir);
  strcpy(jointFormatTreeFileName,      workdir);
  strcpy(labelledTreeFileName,         workdir);
  strcpy(originalLabelledTreeFileName, workdir);
  strcpy(classificationFileName,       workdir);
  
  strcat(entropyFileName,              "RAxML_entropy.");
  strcat(jointFormatTreeFileName,      "RAxML_portableTree.");
  strcat(labelledTreeFileName,         "RAxML_labelledTree.");
  strcat(originalLabelledTreeFileName, "RAxML_originalLabelledTree.");
  strcat(classificationFileName,       "RAxML_classification.");
  
  strcat(entropyFileName,              run_id);
  strcat(jointFormatTreeFileName,      run_id);
  strcat(labelledTreeFileName,         run_id);
  strcat(originalLabelledTreeFileName, run_id);
  strcat(classificationFileName,       run_id);
 
  strcat(jointFormatTreeFileName,      ".jplace");
  
  free(tr->tree_string);
  tr->treeStringLength *= 16;

  tr->tree_string  = (char*)calloc(tr->treeStringLength, sizeof(char));
 
 
  treeFile = myfopen(labelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, FALSE, FALSE, TRUE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);
  
 
  treeFile = myfopen(originalLabelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, FALSE, TRUE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);

  /* 
     JSON format only works for sequential version 
     porting this to Pthreads will be a pain in the ass
     
  */

  treeFile = myfopen(jointFormatTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, TRUE, TRUE);
  
  fprintf(treeFile, "{\n");
  fprintf(treeFile, "\t\"tree\": \"%s\", \n", tr->tree_string);
  fprintf(treeFile, "\t\"placements\": [\n");
      
  {
    info 
      *inf = (info*)malloc(sizeof(info) * tr->numberOfBranches);
        
    for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
      {
	double
	  maxprob = 0.0,
	  lmax = 0.0,
	  acc = 0.0;
	
	int 	   
	  validEntries = 0;

	for(j =  0; j < tr->numberOfBranches; j++) 
	  {
	    inf[j].lh            = tr->bInf[j].epa->likelihoods[i];
	    inf[j].pendantBranch = tr->bInf[j].epa->branches[i];
	    inf[j].distalBranch  = tr->bInf[j].epa->distalBranches[i];
	    inf[j].number        = tr->bInf[j].epa->jointLabel;
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

	fprintf(treeFile, "\t{\"p\":[");

	/* 
	   Erick's cutoff:
	   
	   I keep at most 7 placements and throw away anything that has less than
	   0.01*best_ml_ratio.
	   
	   my old cutoff was at 0.95 accumulated likelihood weight:
	   	     
	   while(acc <= 0.95)
	*/

	/*#define _ALL_ENTRIES*/
	  
#ifdef _ALL_ENTRIES
	assert(validEntries == tr->numberOfBranches);
	while(j < validEntries)	  
#else
	while(j < validEntries && j < 7)	  
#endif
	  { 
	    int 
	      k;
	    
	    double 
	      all = 0.0,
	      prob = 0.0;

	    for(k =  0; k < validEntries; k++) 	   
	      all += exp(inf[k].lh - lmax);	     
	      
	    acc += (prob = (exp(inf[j].lh - lmax) / all));
	      
	    if(j == 0)
	      maxprob = prob;
#ifndef _ALL_ENTRIES
	    if(prob >= maxprob * 0.01)
#endif
	      {
		if(j > 0)
		  {
		    if(tr->wasRooted && inf[j].number == tr->rootLabel)
		      {
			double 
			  b = getBranch(tr, tr->leftRootNode->z, tr->rightRootNode->z);
			
			if(inf[j].distalBranch > 0.5 * b)
			  fprintf(treeFile, ",[%d, %f, %f, %f, %f]", tr->numberOfBranches, inf[j].lh, prob, inf[j].distalBranch - 0.5 * b, inf[j].pendantBranch);
			else
			  fprintf(treeFile, ",[%d, %f, %f, %f, %f]", inf[j].number, inf[j].lh, prob, 0.5 * b - inf[j].distalBranch, inf[j].pendantBranch); 
		      }
		    else
		      fprintf(treeFile, ",[%d, %f, %f, %f, %f]", inf[j].number, inf[j].lh, prob, inf[j].distalBranch, inf[j].pendantBranch);
		  }
		else
		  {
		    if(tr->wasRooted && inf[j].number == tr->rootLabel)
		      {
			double 
			  b = getBranch(tr, tr->leftRootNode->z, tr->rightRootNode->z);
			
			if(inf[j].distalBranch > 0.5 * b)
			  fprintf(treeFile, "[%d, %f, %f, %f, %f]", tr->numberOfBranches, inf[j].lh, prob, inf[j].distalBranch - 0.5 * b, inf[j].pendantBranch);
			else
			  fprintf(treeFile, "[%d, %f, %f, %f, %f]", inf[j].number, inf[j].lh, prob, 0.5 * b - inf[j].distalBranch, inf[j].pendantBranch); 
		      }
		    else
		      fprintf(treeFile, "[%d, %f, %f, %f, %f]", inf[j].number, inf[j].lh, prob,  inf[j].distalBranch, inf[j].pendantBranch);
		  }
	      }
	    	      
	    j++;
	  }
	  
	if(i == tr->numberOfTipsForInsertion - 1)
	  fprintf(treeFile, "], \"n\":[\"%s\"]}\n", tr->nameList[tr->inserts[i]]);
	else
	  fprintf(treeFile, "], \"n\":[\"%s\"]},\n", tr->nameList[tr->inserts[i]]);
      }      
    
    free(inf);      
  }

#ifdef _ALL_ENTRIES
  assert(j == tr->numberOfBranches);
#endif

  fprintf(treeFile, "\t ],\n");
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
  fprintf(treeFile, "\t\"edge_num\", \"likelihood\", \"like_weight_ratio\", \"distal_length\", \n");
  fprintf(treeFile, "\t\"pendant_length\"\n");
  fprintf(treeFile, "\t]\n");
  fprintf(treeFile, "}\n");
  
  fclose(treeFile);

  /* JSON format end */

  classificationFile = myfopen(classificationFileName, "wb");
  
  for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
    for(j = 0; j < tr->numberOfBranches; j++) 
      {       
	if(tr->bInf[j].epa->countThem[i] > 0)    	  
	  fprintf(classificationFile, "%s I%d %d %8.20f\n", tr->nameList[tr->inserts[i]], j, tr->bInf[j].epa->countThem[i], 
		  tr->bInf[j].epa->branches[i] / (double)(tr->bInf[j].epa->countThem[i]));	    
      }
  
  fclose(classificationFile);  

  
  printBothOpen("\n\nLabelled reference tree including branch labels and query sequences written to file: %s\n\n", labelledTreeFileName); 
  printBothOpen("Labelled reference tree with branch labels (without query sequences) written to file: %s\n\n", originalLabelledTreeFileName); 
  printBothOpen("Labelled reference tree with branch labels in portable pplacer/EPA format (without query sequences) written to file: %s\n\n", jointFormatTreeFileName);
  printBothOpen("Classification result file written to file: %s\n\n", classificationFileName);
   

  
  {
    info 
      *inf = (info*)malloc(sizeof(info) * tr->numberOfBranches);
    
    FILE 
      *likelihoodWeightsFile;
    
    char 
      likelihoodWeightsFileName[1024];
    
    strcpy(likelihoodWeightsFileName,       workdir);
    strcat(likelihoodWeightsFileName,       "RAxML_classificationLikelihoodWeights.");
    strcat(likelihoodWeightsFileName,       run_id);
    
    likelihoodWeightsFile = myfopen(likelihoodWeightsFileName, "wb");
    entropyFile           = myfopen(entropyFileName, "wb");
        
    for(i = 0; i < tr->numberOfTipsForInsertion; i++)    
      {
	double
	  entropy = 0.0,
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
	
	while(acc <= 0.95 && j < validEntries)	  
	  { 
	    int 
	      k;
	    
	    double 
	      all = 0.0,
	      prob = 0.0;
	    
	    for(k =  0; k < validEntries; k++) 	   
	      all += exp(inf[k].lh - lmax);	     
	    
	    acc += (prob = (exp(inf[j].lh - lmax) / all));
	    
	    fprintf(likelihoodWeightsFile, "%s I%d %f %f\n", tr->nameList[tr->inserts[i]], inf[j].number, prob, acc);
	    
#ifndef _USE_PTHREADS
	    

#ifdef _SPECIES_STUFF	    
	    /* Species test */
	    
	    if(j == 0)
	      {
		boolean 
		  tipBranch;
		
		nodeptr 
		  attachmentBranch = tr->bInf[inf[j].number].oP,
		  check =  tr->bInf[inf[j].number].oQ;
		
		double 
		  speciesLikelihood;
		
		assert(attachmentBranch == check->back);
		
		tipBranch = (isTip(attachmentBranch->number, tr->mxtips) || isTip(check->number, tr->mxtips));
		
		if(tipBranch)
		  {
		    evaluateGeneric(tr, q);
		    printf("Tree likelihood: %f\n", tr->likelihood);
		    
		    printf("Best placement for read: %s lh: %f \n", tr->nameList[tr->inserts[i]], inf[j].lh);		     
		    
		    speciesLikelihood = testInsertSpecies(tr, r, attachmentBranch, tr->inserts[i], FALSE);
		    
		    printf("species likelihood: %f\n", speciesLikelihood);
		    
		    speciesLikelihood = testInsertSpecies(tr, r, attachmentBranch, tr->inserts[i], TRUE);
		    
		    printf("species likelihood with re-opt: %f\n", speciesLikelihood);
		    
		    printf("Like ratio: %f\n\n", -2.0 * inf[j].lh + 2.0 * speciesLikelihood);		    
		  }
	      }

#endif
#endif 

	    j++;
	  }
	
	j = 0;
	
	while(j < validEntries)
	  { 
	    int 
	      k;
	    
	    double 
	      all = 0.0,
	      prob = 0.0;
	    
	    for(k =  0; k < validEntries; k++) 	   
	      all += exp(inf[k].lh - lmax);	     
	    
	    prob = exp(inf[j].lh - lmax) / all;	      	    
	    
	    if(prob > 0)
	      entropy -= ( prob * log(prob) ); /*pp 20110531 */	      			     
	    
	    j++;
	  }
	
	/* normalize entropy by dividing with the log(validEntries) which is the maximum Entropy possible */
	
	fprintf(entropyFile, "%s\t%f\n", tr->nameList[tr->inserts[i]], entropy / log((double)validEntries));	      	   
      }     
      
    free(inf);

    fclose(entropyFile);
    fclose(likelihoodWeightsFile); 
    
    printBothOpen("Classification result file using likelihood wieghts written to file: %s\n\n", likelihoodWeightsFileName);
    printBothOpen("Classification entropy result file using likelihood wieghts written to file: %s\n\n", entropyFileName);
  }          
     
  exit(0);
}


