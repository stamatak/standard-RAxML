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
extern char bootStrapFile[1024];

#ifdef _USE_PTHREADS
extern volatile int NumberOfThreads;
extern volatile int NumberOfJobs;
#endif

static double getBranch(tree *tr, double *b, double *bb)
{
  double z = 0.0;

  if(!tr->multiBranch)
    {     
      assert(b[0] == bb[0]);
      z = b[0];
      if (z < zmin) 
	z = zmin;      	 
      if(z > zmax)
	z = zmax;
      z = -log(z);
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
	  x = b[i];
	  if (x < zmin) 
	    x = zmin;      	 
	  if(x > zmax)
	    x = zmax;
	  x = -log(x);
	  
	  z += x * tr->partitionContributions[i];
	}	
      
      return z;
    } 

}

static double getBranchPerPartition(tree *tr, double *b, double *bb, int j)
{
  double z = 0.0;

  if(!tr->multiBranch)
    {    
      assert(b[0] == bb[0]);
      z = b[0];
      if (z < zmin) 
	z = zmin;      	 
      if(z > zmax)
	z = zmax;
      z = -log(z);
      return z;	
    }
  else
    {                 
      int 
	i = tr->readPartition[j];
    
      assert(b[i] == bb[i]);    
      z = b[i];
      if (z < zmin) 
	z = zmin;      	 
      if(z > zmax)
	z = zmax;
      z = -log(z);
      
      return z;
    } 

}


static char *Tree2StringClassifyRec(char *treestr, tree *tr, nodeptr p, int *countBranches, 
				    int *inserts, boolean originalTree, boolean jointLabels, boolean likelihood, boolean subtreePlacement)
{        
  branchInfo 
    *bInf = p->bInf;

  //for subtree placements only, this variable tells us if we are 
  //at the branch from which we pruned the subtree
  boolean
    atPruningBranch = FALSE;
  
  int        
    i, 
    countQuery = 0;   

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

  //normally p->bInf and p->back->bInf that define a branch 
  //must point to the same branch meta-data entry, except for one branch 
  //in the subtree placement procedure, this is the branch from which we pruned the 
  //subtree

  if(p->bInf != p->back->bInf)
    {
      //if this was a subtree placement set pruning branch to true
      if(subtreePlacement)
	{
	  assert(originalTree);
	  atPruningBranch = TRUE;	  
	}
      else
	//otherwise exit with an assertion 
	assert(p->bInf == p->back->bInf);
    }
     

  if(isTip(p->number, tr->rdta->numsp)) 
    {
      char 
	*nameptr = tr->nameList[p->number];  
        
      sprintf(treestr, "%s", nameptr);    
      
      while(*treestr) 
	treestr++;
    }
  else 
    {          
      *treestr++ = '(';     
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, 
				       countBranches, inserts, originalTree, jointLabels, likelihood, subtreePlacement);     
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, 
				       countBranches, inserts, originalTree, jointLabels, likelihood, subtreePlacement);          
      *treestr++ = ')';                         
    }
   
  if(countQuery > 0)
    {
      sprintf(treestr, ":%8.20f[%s]", 0.5 * p->bInf->epa->originalBranchLength, p->bInf->epa->branchLabel);
      while (*treestr) 
	treestr++;
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
	    {	     
	      //if we are doing subtree placements and this is the branch from which we pruned the tree
	      //we actually output both branch labels that were connected to the subtree root node that 
	      //we pruned
	      if(subtreePlacement && atPruningBranch)
		sprintf(treestr, ":%8.20f{%d,%d", p->bInf->epa->originalBranchLength, 
			p->bInf->epa->jointLabel, p->back->bInf->epa->jointLabel);
	      else
		sprintf(treestr, ":%8.20f{%d", p->bInf->epa->originalBranchLength, 
			p->bInf->epa->jointLabel);  
	    }
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
			  boolean  originalTree, boolean jointLabels, boolean likelihood, int rootNumber, boolean subtreePlacement)
{
  nodeptr 
    p;
  
  int 
    countBranches = 0; 

      
  if(jointLabels && tr->wasRooted)
    { 
      assert(originalTree);
      assert(!subtreePlacement);
      
      *treestr++ = '(';
      treestr = Tree2StringClassifyRec(treestr, tr, tr->leftRootNode, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood, subtreePlacement);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, tr->rightRootNode, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood, subtreePlacement);	 
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
	p = tr->nodep[rootNumber];
      else
	p = tr->start->back;
      
      assert(!isTip(p->number, tr->mxtips));
      
      *treestr++ = '(';
      treestr = Tree2StringClassifyRec(treestr, tr, p->back, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood, subtreePlacement);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->back, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood, subtreePlacement);
      *treestr++ = ',';
      treestr = Tree2StringClassifyRec(treestr, tr, p->next->next->back, &countBranches, 
				       inserts, originalTree, jointLabels, likelihood, subtreePlacement);
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


void setPartitionMask(tree *tr, int i, boolean *executeModel)
{
  int 
    model = 0;

  if(tr->perPartitionEPA)
    {
      for(model = 0; model < tr->NumberOfModels; model++)   
	executeModel[model] = FALSE;

      executeModel[tr->readPartition[i]] = TRUE;  
    }
  else
    {
      for(model = 0; model < tr->NumberOfModels; model++)   
	executeModel[model] = TRUE;
    }
}

void resetPartitionMask(tree *tr, boolean *executeModel)
{
  int 
    model = 0;
  
  for(model = 0; model < tr->NumberOfModels; model++)
    executeModel[model] = TRUE;
}



static boolean allSmoothedEPA(tree *tr)
{
  int i;
  boolean result = TRUE;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      if(tr->partitionSmoothed[i] == FALSE)
	result = FALSE;
      else
	tr->partitionConverged[i] = TRUE;
    }

  return result;
}

static boolean updateEPA(tree *tr, nodeptr p, int j)
{       
  nodeptr  q; 
  boolean smoothedPartitions[NUM_BRANCHES];
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  double _deltaz;

  q = p->back;   

  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];    
  
  setPartitionMask(tr, j, tr->executeModel);
  makenewzGeneric(tr, p, q, z0, newzpercycle, z, FALSE);
  
#ifdef _BASTIEN
  assert(0);
#endif 

  for(i = 0; i < tr->numBranches; i++)    
    smoothedPartitions[i]  = tr->partitionSmoothed[i];
      
  for(i = 0; i < tr->numBranches; i++)
    {         
      if(!tr->partitionConverged[i])
	{	  
	    _deltaz = deltaz;
	    
	  if(ABS(z[i] - z0[i]) > _deltaz)  
	    {	      
	      smoothedPartitions[i] = FALSE;       
	    }	 	  
	  
	  p->z[i] = q->z[i] = z[i];	 
	}
    }
  
  for(i = 0; i < tr->numBranches; i++)    
    tr->partitionSmoothed[i]  = smoothedPartitions[i];
  
  return TRUE;
}

static boolean localSmoothEPA(tree *tr, nodeptr p, int maxtimes, int j)
{ 
  nodeptr  q;
  int i;
  
  if (isTip(p->number, tr->rdta->numsp)) return FALSE;
  
   for(i = 0; i < tr->numBranches; i++)	
     tr->partitionConverged[i] = FALSE;	

  while (--maxtimes >= 0) 
    {     
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;
	 	
      q = p;
      do 
	{
	  if (! updateEPA(tr, q, j)) return FALSE;
	  q = q->next;
        } 
      while (q != p);
      
      if (allSmoothedEPA(tr)) 
	break;
    }

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }

  return TRUE;
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
	   

	  if(tr->perPartitionEPA)
	    {
	      setPartitionMask(tr, j, tr->executeModel);	     
	      newviewGenericMasked(tr, r);	     

	      setPartitionMask(tr, j, tr->executeModel);
	      localSmoothEPA(tr, r, smoothings, j);

	      setPartitionMask(tr, j, tr->executeModel);
	      evaluateGeneric(tr, r);

	      result = tr->perPartitionLH[tr->readPartition[j]];
	      	      
	      resetPartitionMask(tr, tr->executeModel);
	    }
	  else
	    {
	      newviewGeneric(tr, r); 
	      localSmooth(tr, r, smoothings);
	      result = evaluateGeneric(tr, r);	     
	    }
	  	  

	  if(tr->perPartitionEPA)
	    tr->bInf[q->bInf->epa->branchNumber].epa->branches[j] = getBranchPerPartition(tr, r->z, r->back->z, j);
	  else
	    tr->bInf[q->bInf->epa->branchNumber].epa->branches[j] = getBranch(tr, r->z, r->back->z);	  
	 
	  tr->bInf[q->bInf->epa->branchNumber].epa->likelihoods[j] = result;	 

	  if(tr->perPartitionEPA)
	    modifiedBranchLength = getBranchPerPartition(tr, q->z, q->back->z, j) + getBranchPerPartition(tr, x->z, x->back->z, j);
	  else	      
	    modifiedBranchLength = getBranch(tr, q->z, q->back->z) + getBranch(tr, x->z, x->back->z);

	  ratio = originalBranchLength / modifiedBranchLength;

	  if(tr->wasRooted && atRoot)
	    {	     
	      /* always take distal length from left root node and then fix this later */

	      if(x == tr->leftRootNode)
		{
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, x->z, x->back->z, j);
		  else
		    distalLength = getBranch(tr, x->z, x->back->z);
		}
	      else
		{
		  assert(x == tr->rightRootNode);
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, q->z, q->back->z, j);
		  else
		    distalLength = getBranch(tr, q->z, q->back->z);
		}
	    }
	  else
	    {
	      if(root == x)
		{
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, x->z, x->back->z, j);
		  else
		    distalLength = getBranch(tr, x->z, x->back->z);
		}
	      else
		{
		  assert(root == q); 
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, q->z, q->back->z, j);
		  else
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

	  if(!tr->perPartitionEPA)
	    {
	      result = evaluateGeneric (tr, r);	     	      
	    }
	  else
	    {
	      setPartitionMask(tr, i, tr->executeModel);
	      evaluateGeneric(tr, r);
	      
	      result = tr->perPartitionLH[tr->readPartition[i]];

	      resetPartitionMask(tr, tr->executeModel);	     
	    }

	  
	  r->back = (nodeptr) NULL;
	  tr->nodep[inserts[i]]->back = (nodeptr) NULL;
	  	  
	  tr->bInf[q->bInf->epa->branchNumber].epa->likelihoods[i] = result;	  	     
	}    
    }
 
  hookup(q, x, qz, tr->numBranches);
  
  r->next->next->back = r->next->back = (nodeptr) NULL;
}




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

      b->epa->left  = (double*)rax_malloc(sizeof(double) * tr->contiguousVectorLength);
      b->epa->leftScaling = (int*)rax_malloc(sizeof(int) * tr->contiguousScalingLength);

      b->epa->right = (double*)rax_malloc(sizeof(double)  * tr->contiguousVectorLength);
      b->epa->rightScaling = (int*)rax_malloc(sizeof(int) * tr->contiguousScalingLength);     
    }
}

static void updateClassify(tree *tr, double *z, boolean *partitionSmoothed, boolean *partitionConverged, double *x1, double *x2, unsigned char *tipX1, unsigned char *tipX2, int tipCase, int insertions)
{   
  int i;

  boolean smoothedPartitions[NUM_BRANCHES];

  double 
    newz[NUM_BRANCHES], 
    z0[NUM_BRANCHES];
  
  for(i = 0; i < tr->numBranches; i++)   
    z0[i] = z[i];          

  makenewzClassify(tr, newzpercycle, newz, z0, x1, x2, tipX1, tipX2, tipCase, partitionConverged, insertions); 

#ifdef _BASTIEN
  assert(0);
#endif 

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
				   branchInfo *b, int insertions)
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

	  newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);
	}
      else
	{
	  if (isTip(leftNodeNumber, tr->mxtips))
	    {
	      tipCase = TIP_INNER;

	      tipX1 = tr->contiguousTips[leftNodeNumber];
	      
	      x2  = b->epa->right;	     
	      ex2 = b->epa->rightScaling; 	  
	      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);
	    }
	  else 
	    {
	      if(isTip(rightNodeNumber, tr->mxtips))
		{
		  tipCase = TIP_INNER;

		  tipX1 = tr->contiguousTips[rightNodeNumber];
	  
		  x2  = b->epa->left;	 
		  ex2 = b->epa->leftScaling; 
		  newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1, insertions);
		}       
	      else
		{
		  tipCase = INNER_INNER;
		  
		  x1  = b->epa->left;
		  ex1 = b->epa->leftScaling;
		  
		  x2  = b->epa->right;
		  ex2 = b->epa->rightScaling;
		  newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);
		}
	    }
	}
	
     

      tipCase = TIP_INNER;
      
      x2    = tr->temporaryVector;
      tipX1 = tr->contiguousTips[insertionNodeNumber];

      updateClassify(tr, e3, partitionSmoothed, partitionConverged, x1, x2, tipX1, tipX2, tipCase, insertions);

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
	
      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e3, e2, insertions);

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

      updateClassify(tr, e1, partitionSmoothed, partitionConverged, x1, x3, tipX1, (unsigned char*)NULL, tipCase, insertions);

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
	
      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e3, e1, insertions);

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

      updateClassify(tr, e2, partitionSmoothed, partitionConverged, x1, x3, tipX1, (unsigned char*)NULL, tipCase, insertions);


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

      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);
    }
  else
    {
      if (isTip(leftNodeNumber, tr->mxtips))
	{
	  tipCase = TIP_INNER;

	  tipX1 = tr->contiguousTips[leftNodeNumber];	  

	  x2  = b->epa->right;	     
	  ex2 = b->epa->rightScaling; 	  
	  newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);
	}
      else 
	{
	  if(isTip(rightNodeNumber, tr->mxtips))
	    {
	      tipCase = TIP_INNER;

	      tipX1 = tr->contiguousTips[rightNodeNumber];
	      
	      x2  = b->epa->left;	 
	      ex2 = b->epa->leftScaling; 
	      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1, insertions);
	    }       
	  else
	    {
	      tipCase = INNER_INNER;
	      
	      x1  = b->epa->left;
	      ex1 = b->epa->leftScaling;
	      
	      x2  = b->epa->right;
	      ex2 = b->epa->rightScaling;
	      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);
	    }
	}
    }
  
     
  tipCase = TIP_INNER;
  
  x2    = tr->temporaryVector;
  tipX1 = tr->contiguousTips[insertionNodeNumber];    
 
  result = evalCL(tr, x2, ex3, tipX1, e3, insertions);

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
	      	     
	      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);	
	    }
	  else
	    {
	      if (isTip(leftNodeNumber, tr->mxtips))
		{
		  tipCase = TIP_INNER;
		  
		  tipX1 = tr->contiguousTips[leftNodeNumber];
		  
		  x2  = b->epa->right;	     
		  ex2 = b->epa->rightScaling; 
		  newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);	
		}
	      else 
		{
		  if(isTip(rightNodeNumber, tr->mxtips))
		    {
		      tipCase = TIP_INNER;
		      
		      tipX1 = tr->contiguousTips[rightNodeNumber];
		      
		      x2  = b->epa->left;	 
		      ex2 = b->epa->leftScaling;
		      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e2, e1, insertions);	
		    }       
		  else
		    {
		      tipCase = INNER_INNER;
		      
		      x1  = b->epa->left;
		      ex1 = b->epa->leftScaling;
		      
		      x2  = b->epa->right;
		      ex2 = b->epa->rightScaling;
		      newviewMultiGrain(tr, x1, x2, x3, ex1, ex2, ex3, tipX1, tipX2, tipCase, e1, e2, insertions);	
		    }
		}	      
	    }
	  
	  result = localSmoothClassify(tr, smoothings, leftNodeNumber, rightNodeNumber, tr->inserts[insertions], e1, e2, e3, b, insertions);
	      	    
	  b->epa->likelihoods[insertions] = result;	      			      
	  	
	  if(tr->perPartitionEPA)
	    b->epa->branches[insertions] = getBranchPerPartition(tr, e3, e3, insertions);	
	  else
	    b->epa->branches[insertions] = getBranch(tr, e3, e3);	  

	  if(tr->perPartitionEPA)
	    modifiedBranchLength = getBranchPerPartition(tr, e1, e1, insertions) + getBranchPerPartition(tr, e2, e2, insertions);
	  else
	    modifiedBranchLength = getBranch(tr, e1, e1) + getBranch(tr, e2, e2);

	  ratio = b->epa->originalBranchLength / modifiedBranchLength;

	  if(tr->wasRooted && atRoot)
	    {	     
	      /* always take distal length from left root node and then fix this later */

	      if(x == tr->leftRootNode)
		{
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, e1, e1, insertions);
		  else
		    distalLength = getBranch(tr, e1, e1);
		}
	      else
		{
		  assert(x == tr->rightRootNode);
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, e2, e2, insertions);
		  else
		    distalLength = getBranch(tr, e2, e2);
		}
	    }
	  else
	    {
	      if(root == x)
		{
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, e1, e1, insertions);
		  else
		    distalLength = getBranch(tr, e1, e1);
		}
	      else
		{
		  assert(root == q);
		  if(tr->perPartitionEPA)
		    distalLength = getBranchPerPartition(tr, e2, e2, insertions);
		  else
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
                     
  newviewClassify(tr, b, z, inserts);   
        
  for(inserts = 0; inserts < tr->numberOfTipsForInsertion; inserts++) 
    { 	       		     
      result = evalCL(tr, tr->temporaryVector, tr->temporaryScaling, tr->contiguousTips[tr->inserts[inserts]], defaultArray, inserts);
	  	  
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
      if(!tr->doSubtreeEPA)
	{
	  if(!p->back->x)
	    newviewGeneric(tr, p->back);
	  masterBarrier(THREAD_GATHER_LIKELIHOOD, tr);
	}
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
      if(!tr->doSubtreeEPA)
	{
	  if(!p->x)
	    newviewGeneric(tr, p);            
	  
	  if(!isTip(p->back->number, tr->mxtips))
	    {
	      if(!p->back->x)
		newviewGeneric(tr, p->back);	 
	    }	      
	  
	  masterBarrier(THREAD_GATHER_LIKELIHOOD, tr);
	}
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
    *inf = (info*)rax_malloc(sizeof(info) * (size_t)tr->numberOfBranches);

  assert(tr->useEpaHeuristics);

  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {     
      for(i = 0; i < tr->numberOfBranches; i++)
	{      
	  inf[i].lh = tr->bInf[i].epa->likelihoods[j];
	  inf[i].number = i;
	}
      
      qsort(inf, (size_t)tr->numberOfBranches, sizeof(info), infoCompare);

      for(i = throwAwayStart; i < tr->numberOfBranches; i++)       
	tr->bInf[inf[i].number].epa->executeThem[j] = 0;	   	
    }
  
   for(i = 0; i < tr->numberOfBranches; i++)
     for(j = 0; j < tr->numberOfTipsForInsertion; j++)    
       tr->bInf[i].epa->likelihoods[j] = unlikely;     
  

  rax_free(inf);
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







static void analyzeReads(tree *tr)
{ 
  int 
    i,
    *inserts = tr->inserts;

  tr->readPartition = (int *)rax_malloc(sizeof(int) * (size_t)tr->numberOfTipsForInsertion);
  
  for(i = 0; i < tr->numberOfTipsForInsertion; i++)
    {
      size_t
	j;
      
      int
	whichPartition = -1,
	partitionCount = 0,
	model,
	nodeNumber = tr->nodep[inserts[i]]->number;

       unsigned char 	
	 *tipX1 = tr->yVector[nodeNumber];
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int 	   
	    nonGap = 0;
	  
	  unsigned char
	    undetermined = getUndetermined(tr->partitionData[model].dataType);

	  for(j = tr->partitionData[model].lower; j < tr->partitionData[model].upper; j++)
	    {	      
	      if(tipX1[j] != undetermined)
		nonGap++;
	    }
       
	  if(nonGap > 0)
	    {
	      partitionCount++;
	      whichPartition = model;
	    }	  	  
	}
      
      assert(partitionCount == 1);
      assert(whichPartition >= 0 && whichPartition < tr->NumberOfModels);    

      tr->readPartition[i] = whichPartition; 
    }
}

static void printStandardFormat(tree *tr, char *jointFormatTreeFileName, int rootNumber, int numberOfTipsForInsertion, boolean subtreePlacement)
{
  FILE
    *treeFile = myfopen(jointFormatTreeFileName, "wb"); 

  info 
    *inf = (info*)rax_malloc(sizeof(info) * (size_t)tr->numberOfBranches);

  int
    i,
    j;
  
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, TRUE, TRUE, rootNumber, subtreePlacement);
  
  fprintf(treeFile, "{\n");
  fprintf(treeFile, "\t\"tree\": \"%s\", \n", tr->tree_string);
  fprintf(treeFile, "\t\"placements\": [\n");
                 
  for(i = 0; i < numberOfTipsForInsertion; i++)    
    {
      double
	all,
	maxprob = 0.0,
	lmax = 0.0,
	acc = 0.0;
      
      int 
	k,
	validEntries = 0;
      
      for(j =  0; j < tr->numberOfBranches; j++) 
	{
	  inf[j].lh            = tr->bInf[j].epa->likelihoods[i];
	  inf[j].pendantBranch = tr->bInf[j].epa->branches[i];
	  inf[j].distalBranch  = tr->bInf[j].epa->distalBranches[i];
	  inf[j].number        = tr->bInf[j].epa->jointLabel;
	}
      
      qsort(inf, (size_t)tr->numberOfBranches, sizeof(info), infoCompare);	 
      
      for(j =  0; j < tr->numberOfBranches; j++) 
	if(inf[j].lh == unlikely)
	  break;
	else	     
	  validEntries++;	  

      //for subtree placements make sure we have as many valid entries in the list 
      //as we have branches in the reference tree (after sutree pruning!)
      if(subtreePlacement)
	assert(validEntries == (2 * tr->ntips) - 3);
           
      //otherwise also make sure that there is at least one valid entry!
      assert(validEntries > 0);
      
      j = 0;
      
      lmax = inf[0].lh;
      
      for(k = 0, all = 0.0; k < validEntries; k++)
	all += exp(inf[k].lh - lmax);
      
      fprintf(treeFile, "\t{\"p\":[");
      
      /* 
	 Erick's cutoff:
	 
	 I keep at most 7 placements and throw away anything that has less than
	 0.01*best_ml_ratio.
	 
	 my old cutoff was at 0.95 accumulated likelihood weight:
	 
	 while(acc <= 0.95)
      */            
      
      assert(tr->numberOfEPAEntries > 0);
      tr->numberOfEPAEntries = MIN(tr->numberOfEPAEntries, (unsigned int)validEntries);

      while(j < validEntries)	  
	  { 
	    if((!tr->useAccumulatedEPACutoff && (unsigned int)j < tr->numberOfEPAEntries) || (tr->useAccumulatedEPACutoff && acc <= tr->accumulatedEPACutoff))
	      {
		double 
		  prob = 0.0;
	    
		acc += (prob = (exp(inf[j].lh - lmax) / all));
	    
		if(j == 0)
		  maxprob = prob;

		if((!tr->useAccumulatedEPACutoff && prob >= maxprob * tr->probThresholdEPA) || (tr->useAccumulatedEPACutoff))
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
	      }
	    	      
	    j++;
	  }
      
      if(subtreePlacement)
	{
	  if(i == numberOfTipsForInsertion - 1)
	    fprintf(treeFile, "], \"n\":[\"%s\"]}\n", "subtree");
	  else
	    fprintf(treeFile, "], \"n\":[\"%s\"]},\n", "subtree");
	}
      else
	{
	  if(i == numberOfTipsForInsertion - 1)
	    fprintf(treeFile, "], \"n\":[\"%s\"]}\n", tr->nameList[tr->inserts[i]]);
	  else
	    fprintf(treeFile, "], \"n\":[\"%s\"]},\n", tr->nameList[tr->inserts[i]]);
	}
    }      
    
  rax_free(inf);      

#ifdef _ALL_ENTRIES
  assert(j == tr->numberOfBranches);
#endif

  fprintf(treeFile, "\t ],\n");
  fprintf(treeFile, "\t\"metadata\": {\"invocation\": ");

  fprintf(treeFile, "\"");
  
    
  for(i = 0; i < globalArgc; i++)
    fprintf(treeFile,"%s ", globalArgv[i]);
  
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
}

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

  if(tr->perPartitionEPA)
    printBothOpen("\nRAxML Evolutionary Placement Algorithm for partitioned multi-gene datasets (experimental version)\n"); 
  else
    printBothOpen("\nRAxML Evolutionary Placement Algorithm\n"); 

  if(adef->useBinaryModelFile)
    {      
      if(tr->ntips != tr->binaryFile_ntips)
	{
	  printf("\nError: number of %d tips in reference tree and %d tips in the tree used to estimate\n", tr->ntips, tr->binaryFile_ntips);
	  printf("the binary model parameter file do not match, RAxML will exit now!\n\n");
	}
      assert(tr->ntips == tr->binaryFile_ntips);
      evaluateGenericInitrav(tr, tr->start);
      //treeEvaluate(tr, 2);
    }
  else
    {
      evaluateGenericInitrav(tr, tr->start); 
  
      modOpt(tr, adef, TRUE, 1.0);
    }

  printBothOpen("\nLikelihood of reference tree: %f\n\n", tr->likelihood);

  perm    = (int *)rax_calloc((size_t)tr->mxtips + 1, sizeof(int));
  tr->inserts = (int *)rax_calloc((size_t)tr->mxtips, sizeof(int));

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

  rax_free(perm);
  
  printBothOpen("RAxML will place %d Query Sequences into the %d branches of the reference tree with %d taxa\n\n",  tr->numberOfTipsForInsertion, (2 * tr->ntips - 3), tr->ntips);  

  assert(tr->numberOfTipsForInsertion == (tr->mxtips - tr->ntips));      

  tr->bInf              = (branchInfo*)rax_malloc((size_t)tr->numberOfBranches * sizeof(branchInfo)); 

  for(i = 0; i < tr->numberOfBranches; i++)
    {      
      tr->bInf[i].epa = (epaBranchData*)rax_malloc(sizeof(epaBranchData));

      tr->bInf[i].epa->countThem   = (int*)rax_calloc((size_t)tr->numberOfTipsForInsertion, sizeof(int));      
      
      tr->bInf[i].epa->executeThem = (int*)rax_calloc((size_t)tr->numberOfTipsForInsertion, sizeof(int));
      
      for(j = 0; j < tr->numberOfTipsForInsertion; j++)
	tr->bInf[i].epa->executeThem[j] = 1;

      tr->bInf[i].epa->branches          = (double*)rax_calloc((size_t)tr->numberOfTipsForInsertion, sizeof(double));   
      tr->bInf[i].epa->distalBranches    = (double*)rax_calloc((size_t)tr->numberOfTipsForInsertion, sizeof(double)); 
         
      tr->bInf[i].epa->likelihoods = (double*)rax_calloc((size_t)tr->numberOfTipsForInsertion, sizeof(double));      
      tr->bInf[i].epa->branchNumber = i;
      
      sprintf(tr->bInf[i].epa->branchLabel, "I%d", i);     
    } 

  r = tr->nodep[(tr->nextnode)++]; 
    

  q = findAnyTip(tr->start, tr->rdta->numsp);

  assert(isTip(q->number, tr->rdta->numsp));
  assert(!isTip(q->back->number, tr->rdta->numsp));
	 
  q = q->back;       
  
  if(tr->perPartitionEPA)
    analyzeReads(tr);

 

#ifdef _USE_PTHREADS
  tr->contiguousVectorLength = getContiguousVectorLength(tr);
  tr->contiguousScalingLength = getContiguousScalingLength(tr);
  allocBranchX(tr);
  masterBarrier(THREAD_INIT_EPA, tr);
#endif 
 
  setupBranchInfo(tr, q);   
  
#ifdef _USE_PTHREADS  	
  masterBarrier(THREAD_FREE_VECTORS, tr); 
#endif

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
  
  rax_free(tr->tree_string);
  tr->treeStringLength *= 16;

  tr->tree_string  = (char*)rax_calloc((size_t)tr->treeStringLength, sizeof(char));
 
 
  treeFile = myfopen(labelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, FALSE, FALSE, TRUE, tr->mxtips + 1, FALSE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);
  
 
  treeFile = myfopen(originalLabelledTreeFileName, "wb");
  Tree2StringClassify(tr->tree_string, tr, tr->inserts, TRUE, FALSE, TRUE, tr->mxtips + 1, FALSE);
  fprintf(treeFile, "%s\n", tr->tree_string);    
  fclose(treeFile);

  /* 
     JSON format only works for sequential version 
     porting this to Pthreads will be a pain in the ass     
  */
 
  printStandardFormat(tr, jointFormatTreeFileName, tr->mxtips + 1, tr->numberOfTipsForInsertion, FALSE);
  
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
      *inf = (info*)rax_malloc(sizeof(info) * (size_t)tr->numberOfBranches);
    
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
	  all,
	  entropy = 0.0,
	  lmax = 0.0,
	  acc = 0.0;
	
	int
	  k,
	  validEntries = 0;
	
	for(j =  0; j < tr->numberOfBranches; j++) 
	  {
	    inf[j].lh = tr->bInf[j].epa->likelihoods[i];
	    inf[j].number = j;
	  }
	
	qsort(inf, (size_t)tr->numberOfBranches, sizeof(info), infoCompare);	 
	
	for(j =  0; j < tr->numberOfBranches; j++) 
	  if(inf[j].lh == unlikely)
	    break;
	  else	     
	    validEntries++;	     	      
	
	assert(validEntries > 0);
	
	j = 0;
	
	lmax = inf[0].lh;
	   
	for(k = 0, all = 0.0; k < validEntries; k++)
	  all += exp(inf[k].lh - lmax);

	while(acc <= 0.95 && j < validEntries)	  
	  { 
	    double 
	      prob = 0.0;
	    
	    acc += (prob = (exp(inf[j].lh - lmax) / all));
	    
	    fprintf(likelihoodWeightsFile, "%s I%d %f %f\n", tr->nameList[tr->inserts[i]], inf[j].number, prob, acc);
	   	    
	    j++;
	  }
	
	j = 0;
	
	while(j < validEntries)
	  { 
	    double 
	      prob = 0.0;

	    prob = exp(inf[j].lh - lmax) / all;	      	    
	    
	    if(prob > 0)
	      entropy -= ( prob * log(prob) ); /*pp 20110531 */	      			     
	    
	    j++;
	  }
	
	/* normalize entropy by dividing with the log(validEntries) which is the maximum Entropy possible */
	
	fprintf(entropyFile, "%s\t%f\n", tr->nameList[tr->inserts[i]], entropy / log((double)validEntries));	      	   
      }     
      
    rax_free(inf);

    fclose(entropyFile);
    fclose(likelihoodWeightsFile); 
    
    printBothOpen("Classification result file using likelihood weights written to file: %s\n\n", likelihoodWeightsFileName);
    printBothOpen("Classification entropy result file using likelihood weights written to file: %s\n\n", entropyFileName);
  }          
     
  exit(0);
}

/***** subtree EPA ****/



//#define _DEBUG_SUBTREE_EPA

//function to determine if nodeptr p is the root of a subtree 
//that contains all taxa of the subtree specified in a line of the 
//subtree specification file passed via -z

static boolean isRoot(nodeptr p, tree *tr, int tipsFound, int *subtreeTips, int *counter)
{
  if(isTip(p->number, tr->rdta->numsp))
    {
      int 
	i;

      //if it is a tip, check if the tip number is contained in the list of tip numbers
      //stored in subtreeTips that span the subtree we want to place

      for(i = 0; i < tipsFound; i++)
	{
	  if(subtreeTips[i] == p->number)
	    {
	      //if this tip forms part of the subtree we are looking for 
	      //increment the counter of subtree tips we have found and return TRUE
	      *counter = *counter + 1;
	      return TRUE;
	    }
	}
      return FALSE;
      
    }
  else
    {     
      //if our node p is an annier node 
      //check if the left and right subtrees only contain tips that belong to the subtree we want to place
      //if this is the case return TRUE

      nodeptr 
	q = p->next;

      while(q != p)
	{
	  if(isRoot(q->back, tr, tipsFound, subtreeTips, counter) == FALSE)
	    return FALSE;
	  q = q->next;
	}

      return TRUE;
    }
}


//function that checks if the subtree that shall be placed as specified in the 
//file passed via -z is monophyletic, if this is the case it returns the root of that subtree
//otherwise it will return NULL
static nodeptr findRoot(tree *tr, int tipsFound, int *subtreeTips)
{
  nodeptr 
    p = (nodeptr)NULL,
    *subtrees = (nodeptr *)rax_malloc(sizeof(nodeptr) * 3 * (size_t)tr->mxtips);

  boolean 
    monophyletic = FALSE;

  int 
    i,
    count = 0,
    tipCounter = 0;
    
  //function that finds the roots of all subtrees containing tipsFound taxa 
  //in the comprehensive reference tree 
  //it stores the subtree roots in an array called subtrees and 
  //returns the number of subtrees it found in count
  collectSubtrees(tr, subtrees, &count, tipsFound);
  
  //now loop over all subtrees that have the size of the subtree we are looking for 
  //and determine if one of them contains all taxa of the subtree we actually want to place
  
  for(i = 0; (i < count) && (!monophyletic); i++)
    if(isRoot(subtrees[i], tr, tipsFound, subtreeTips, &tipCounter))
      {
	//set the return value if this is the subtree we have been looking for
	p = subtrees[i];
	
	//break the loop
	monophyletic = TRUE;
	
	//make sure that the subtree really contains all taxa of the subtree we want to place
	assert(tipCounter == tipsFound);
      }	    

  rax_free(subtrees);

  //return the subtree root 

  return p;
}

#ifdef _DEBUG_SUBTREE_EPA
//debugging function for printing the subtree we want to place, once we have found his root

static void printRec(nodeptr p, tree *tr)
{
  if(isTip(p->number, tr->mxtips))
    printf("%s", tr->nameList[p->number]);
  else
    {
      nodeptr
	q = p->next;

      printf("(");
      while(q != p)
	{
	  printRec(q->back, tr);	  

	  if(q->next != p)
	    printf(",");
	  q = q->next;
	}
      printf(")");
    }
}
#endif
  
//function that computes the insertion likelihood of our subtree rooted at subTreeRoot 
//into a branch q <-> q->back of the remaining reference tree
//the rootNumber parameter is required for the jplace standard output format to determine 
//the direction of the root for defining pendant and distal branch lengths

static void testInsertSubtree(tree *tr, nodeptr subTreeRoot, nodeptr q, int rootNumber)
{
  double 
    //store the original branch length of the reference tree
    originalBranchLength = getBranch(tr, q->z, q->back->z),
    result,           
    qz[NUM_BRANCHES],
    z[NUM_BRANCHES];
  
  nodeptr      
    root = (nodeptr)NULL,  
    x = q->back; 

  int 
    j;
  
  //rooted trees and the EPA extension for mult-gene datasets are currently not allowed!
  assert(!tr->wasRooted);
  assert(!tr->perPartitionEPA);
  
  //determine if the virtual root for the jplace output format lies in 
  //the direction of or q->back, that is, on which end of our insertion branch it lies
  root = findRootDirection(q, tr, rootNumber);

  //make sure that we found it!
  assert(root);

  //store the original branch lengths between 
  //q and q->back (x)
  for(j = 0; j < tr->numBranches; j++)    
    {
      qz[j] = q->z[j];
      z[j] = sqrt(qz[j]); 

      if(z[j] < zmin) 
	z[j] = zmin;
      
      if(z[j] > zmax)
	z[j] = zmax;
    }  
    
  //insert the subtree into the branch connecting q and x 
  hookup(subTreeRoot->next,       q, z, tr->numBranches);
  hookup(subTreeRoot->next->next, x, z, tr->numBranches);
  
  //set a default value for the branch to which the subtree is attached
  hookupDefault(subTreeRoot, subTreeRoot->back, tr->numBranches);      		     
  
  //re-calculate the cond likelihood vector at the subtree root
  newviewGeneric(tr, subTreeRoot);

  //optimize the three branch lengths around the insertion position
  localSmooth(tr, subTreeRoot, smoothings);

  //calculate the insertion likelihood
  result = evaluateGeneric(tr, subTreeRoot);
  

  //now store data for generating jplace output later-on
  {
    double 
      modifiedBranchLength = 0.0,
      distalLength = 0.0,
      ratio = 0.0;
    

    //get the pendant branch for this placement
    tr->bInf[q->bInf->epa->branchNumber].epa->branches[0] = getBranch(tr, subTreeRoot->z, subTreeRoot->back->z);

    //store the likelihood for placing the subtree in here
    tr->bInf[q->bInf->epa->branchNumber].epa->likelihoods[0] = result;
    
    //get the sum of branch lengths on the path from q to x on which we inserted our subtree
    modifiedBranchLength = getBranch(tr, q->z, q->back->z) + getBranch(tr, x->z, x->back->z);
    
    //calculate the ration between the original branch connecting q and x 
    //and the current pathe length between q and x 
    ratio = originalBranchLength / modifiedBranchLength;
    
    //if the root lies in the direction of x 
    //the branch leading from the subtree insertion point to x is the distal branch
    if(root == x)		
      distalLength = getBranch(tr, x->z, x->back->z);	       
    //otherwise it's the branch leading from the subtree insertion point to q
    else
      {
	assert(root == q); 	
	distalLength = getBranch(tr, q->z, q->back->z);
      }
    
    //now scale the distal length by the ratio 
    distalLength *= ratio;
    
    //must be shorter than the originalBranch 
    assert(distalLength <= originalBranchLength);
    
    //store the distal branch
    tr->bInf[q->bInf->epa->branchNumber].epa->distalBranches[0] = distalLength;
  }
  
#ifdef _DEBUG_SUBTREE_EPA
  printf("insertion likelihood %f\n", result);
#endif
  
  //repair the original branch between q and x such that we can move on 
  hookup(q, x, qz, tr->numBranches);
  
  //set the attachment points of the subtree root to NULL for a more clean implementation
  //the actual subtree we are placing in the tree is attached to subTreeRoot->back!
  subTreeRoot->next->next->back = subTreeRoot->next->back = (nodeptr) NULL; 
}

//function to recursively place the subtree rooted at subtreeRoot into all branches of the reference 
//tree
static void placeSubtreeRec(tree *tr, nodeptr subtreeRoot, nodeptr q, int *insertionCount, int rootNumber)
{     
  //insert subtree and calculate likelihood into the branch q and q->back
  testInsertSubtree(tr, subtreeRoot, q, rootNumber);
  
  //increment the number of subtree insertions we have conducted
  *insertionCount = *insertionCount + 1;

  //if q is not a tip descend into the left and right subtree of q
  if(!isTip(q->number, tr->rdta->numsp))
    {   
      nodeptr 
	a = q->next;

      while(a != q)
	{
	  placeSubtreeRec(tr, subtreeRoot, a->back, insertionCount, rootNumber);
	  a = a->next;
	}      
    }
} 

//function to calculate all placements for a specific subtree into the remaining tree
//subTreeRoot is the root of the subtree to be placed
//tipsFound the number of tips this subtree contains
//subTreeIndex the index of the subtree in the input file

static void placeSubtree(nodeptr subTreeRoot, tree *tr, int tipsFound, int subTreeIndex)
{
  int 
    rootNumber = 0,
    //number of branches in the reference tree
    branchesInReferenceTree = 2 * (tr->mxtips - tipsFound) - 3,
    //number of tips in the reference tree
    taxaInReferenceTree = tr->mxtips - tipsFound,
    insertionCount = 0, 
    i;

  //branch length buffers to store the branches such as to be able to repair the 
  //original tree
  double
    z1[NUM_BRANCHES],
    z2[NUM_BRANCHES],
    z3[NUM_BRANCHES];
  
  nodeptr 
    //nodes to repair the pruning of the subtree to be placed 
    //and restore the original comprehensive tree once we are done with placing
    //thus subtree
    p1 = subTreeRoot->back,
    p2 = subTreeRoot->next->back,
    p3 = subTreeRoot->next->next->back,
    q,
    start;

  char 
    buf[64] = "",
    subTreeFileName[1024] = "";
  //first generate a file name for the jplace output file for the placements of the current subtree

  sprintf(buf, "%d", subTreeIndex);

  strcpy(subTreeFileName, workdir);
  strcat(subTreeFileName, "RAxML_subtreePlacement.");
  strcat(subTreeFileName, run_id);
  strcat(subTreeFileName, ".");
  strcat(subTreeFileName, buf);
  strcat(subTreeFileName, ".jplace");  

  //store the branch lengths of the initial comprehensive tree
  for(i = 0; i < tr->numBranches; i++)
    {
      z1[i] = p1->z[i];
      z2[i] = p2->z[i];
      z3[i] = p3->z[i];
    }
	
  //prune the subtree and return the node q that is located on 
  //one end of the branch from which we pruned the subtree
  q = removeNodeBIG(tr, subTreeRoot, tr->numBranches);

  //find an arbitrary tip in the remaining tree into which we are going to place the subtree
  start = findAnyTip(q, tr->mxtips);
  
  //set an ARBITRARY inner node of the remaining subtree as root, this si required 
  //for the jplace format!
  rootNumber = start->back->number;          

  //make sure that the rootNumber node for the jplace format is an inner node and not a tip
  assert(rootNumber > tr->mxtips);

  //do a full traversal of the reference tree, might not be required, but feels safer
  evaluateGenericInitrav(tr, start);

  //TODO maybe re-optimize branch lengths as well?
  //Alexey could you please test ?
  //Note that, for this we need to save all branch lengths prior to pruning the subtree and 
  //restore them again before we exit the function again!

  printBothOpen("\nLikelihood of reference tree for subtree %d: %f\n\n", subTreeIndex, tr->likelihood);

  //now call the function that calculates the likelihoods for 
  //placing the subtree into all branches of the reference tree
  placeSubtreeRec(tr, subTreeRoot, start->back, &insertionCount, rootNumber);

  //make sure that we have placed the subtree into all branches
  assert(insertionCount == branchesInReferenceTree);

 
  //set the actual number of tips in the reference tree, required such that 
  //an assertion in printStandardFormat does not fail
  tr->ntips =  taxaInReferenceTree;

  //print jplace-formatted placement result to file
  printStandardFormat(tr, subTreeFileName, rootNumber, 1, TRUE);


  //repair the original comprehensive tree by re-inserting the subtree 
  //into its original position with original branch length values
  hookup(subTreeRoot,             p1, z1, tr->numBranches);
  hookup(subTreeRoot->next,       p2, z2, tr->numBranches);
  hookup(subTreeRoot->next->next, p3, z3, tr->numBranches);
  
  //do a full traversal of the now, again comprehensive tree, maybe not required but feels 
  //safer
  evaluateGenericInitrav(tr, tr->start);

#ifdef _DEBUG_SUBTREE_EPA
  printf("repaired like %f\n\n", tr->likelihood);
#endif

  printBothOpen("Subtree placements for subtree %d written to file %s\n\n", subTreeIndex, subTreeFileName);

  printBothOpen("****************************************************************\n");
}

//top level function that computes subtree placements for a list of subtrees
void subtreeEPA(tree *tr, analdef *adef) 
{
  int  
    i = 0,
    //input file line count
    lineCount = 0,
    //count of subtrees in that file
    subTreeCount = 0,
    //data structure for storing tip node numbers of the subtree under consideration
    *subtreeTips = (int *)rax_malloc(sizeof(int) * (size_t)(tr->mxtips)); 
  
  char 
    *line = (char*)NULL;

  size_t 
    len = 0;

  ssize_t 
    read;     
  
  FILE 
    *f = myfopen(bootStrapFile, "rb");

  nodeptr 
    q;
  
  //no outgroups and hauristic likelihood evaluations allowd
  adef->outgroup = FALSE;
  tr->doCutoff   = FALSE;
  
  //make sure that an input tree was provided via -t 
  //and the the input tree contains all taxa!
  assert(adef->restart);
  assert(tr->ntips == tr->mxtips);
  
  //set tree was rooted to FALSE,
  //since the reference trees will always be subtrees of the input tree
  //the information on the root is potentially useless
  //note that, the order of branch labels is maintained, however the way the output reference trees are actually rooted may change!
  if(tr->wasRooted)
    {      
      printBothOpen("WARNING: setting rooted to FALSE!\n");
      tr->wasRooted = FALSE;
    }

  //set up variables required for storing branch meta-data

  tr->branchCounter = 0;
  tr->numberOfBranches = (2 * tr->mxtips) - 3;

  //initialize branch meta-data structure

  tr->bInf              = (branchInfo*)rax_malloc((size_t)tr->numberOfBranches * sizeof(branchInfo));

  for(i = 0; i < tr->numberOfBranches; i++)
    {      
      tr->bInf[i].epa = (epaBranchData*)rax_malloc(sizeof(epaBranchData));      

      tr->bInf[i].epa->branches          = (double*)rax_calloc(1, sizeof(double));   
      tr->bInf[i].epa->distalBranches    = (double*)rax_calloc(1, sizeof(double)); 
         
      tr->bInf[i].epa->likelihoods = (double*)rax_calloc(1, sizeof(double));      
      tr->bInf[i].epa->branchNumber = i;
      
      sprintf(tr->bInf[i].epa->branchLabel, "I%d", i);     
    } 

  //link branch meta-data information to the full, comprehensive reference tree
  //just like in standard EPA

  q = findAnyTip(tr->start, tr->rdta->numsp);

  setupBranchInfo(tr, q->back);
    
  //make sure we linked branch meta-data to as many branches as the tree actually has
  assert(tr->branchCounter == tr->numberOfBranches);

  //use binary model file for avoiding re-optimization of model parameters 
  //for faster production deployement 
  //Attention: not tested yet!
  if(adef->useBinaryModelFile)
    {      
      if(tr->ntips != tr->binaryFile_ntips)
	{
	  printf("\nError: number of %d tips in reference tree and %d tips in the tree used to estimate\n", tr->ntips, tr->binaryFile_ntips);
	  printf("the binary model parameter file do not match, RAxML will exit now!\n\n");
	}
      assert(tr->ntips == tr->binaryFile_ntips);
      evaluateGenericInitrav(tr, tr->start);
      //treeEvaluate(tr, 2);
    }
  //just re-optimize all parameters 
  else
    {
      evaluateGenericInitrav(tr, tr->start); 
  
      modOpt(tr, adef, TRUE, 0.1);
    }

  printBothOpen("\nLikelihood of comprehensive tree: %f\n\n", tr->likelihood);
  
  printBothOpen("****************************************************************\n");

  //loop over the lines of the input file 

  while((read = rax_getline(&line, &len, f)) != -1)
    {
      nodeptr
	p,
	subTreeRoot;

      int 	
	tipsFound = 0;

      i = 0;         

      while(i < read - 1)
	{	  	  	  
	  //if there is not a whitechar in the current line start putting together the taxon name 
	  //it contains until we meet the next whitespace
	  if(!whitechar(line[i]))
	    {
	      int	       
		tipNumber = -1,
		nameLength = 0;
	      
	      char  
		str[nmlngth + 2]="";      	      

	      //copy taxon name from the line of the file into a string buffer
	      do
		{
		  str[nameLength] = line[i];
		  i++;
		  nameLength++;
		}
	      while(!whitechar(line[i]) && i < read - 1);
	 
	      //make sure the taxon name is not too long
	      assert(nameLength < nmlngth + 2);

	      //look up if the taxon name specified in the file actually exists 
	      tipNumber = treeFindTipByLabelString(str, tr, FALSE);
	      
	      if(tipNumber == 0)
		{
		  printBothOpen("\nCouldn't find tip with name %s, exiting\n", str);
		  exit(-1);
		}

#ifdef _DEBUG_SUBTREE_EPA
	      printf("tipnumber %d\n", tipNumber);
#endif
	      //if the taxon name exists add it to our list of taxa 
	      //that span the current subtree

	      subtreeTips[tipsFound] = tipNumber;
	      
#ifdef _DEBUG_SUBTREE_EPA	      
	      printf("found tip %s\n", tr->nameList[tipNumber]);
#endif
	      
	      tipsFound++;
	    }
	  else
	    i++;
	}       
      
      lineCount++;

      //now that we have processed the entire line and determined 
      //that all taxon names in that line exist determine if they actually span
      //a subtree

      if(tipsFound > 0)	
	{
	  printBothOpen("\nFound a subtree with %d tips for placement\n", tipsFound);
          	
	  //determine if the taxa in the list span a subtree of the current tree
	  p = findRoot(tr, tipsFound, subtreeTips); 
	  
	  //if they don't exit 
	  if(p == (nodeptr)NULL)
	    {
	      printBothOpen("\nSubtree in line %d of input file is not monophyletic as it ought to be, exiting ...\n", lineCount);
	      exit(-1);
	    }

#ifdef _DEBUG_SUBTREE_EPA
	  printRec(p, tr);
	  printf("\n\n");
#endif

	  //p is the root node of the subtree but we actually need to go back 
	  //to p->back which is the node that defines the branch on which the 
	  //subtree root is hanging, it's this node we actually need to place 
	  //into the remaining reference tree
	  subTreeRoot = p->back;

	  //make sure the node from which the subtree is hanging is not a tip
	  assert(!isTip(subTreeRoot->number, tr->mxtips));

	  //make sure the remaining tree has at least 3 taxa, otherwise this doesn't make sense
	  if(tr->mxtips - tipsFound > 2)
	    {
	      int 
		j;

	      subTreeCount++;

	      //make sure to reset all entries in the branch meta-data structure to default values 
	      //because we will be re-using them for every subtree specified in the file
	      for(j = 0; j < tr->numberOfBranches; j++) 
		{
		  tr->bInf[j].epa->likelihoods[0]    = unlikely;  
		  tr->bInf[j].epa->branches[0]       = -1.0;
		  tr->bInf[j].epa->distalBranches[0] = -1.0;		 
		}

	      placeSubtree(subTreeRoot, tr, tipsFound, subTreeCount);	   
	    }
	  else
	    {	 
	      //just skip the subtree
	      printBothOpen("The subtree specified in line %d of the input file is so large that\n", lineCount);
	      printBothOpen("the remaining tree into which it shall be placed only has %d taxa\n", tr->mxtips - tipsFound);
	      printBothOpen("... skipping this subtree ... \n");
	    }
	}         

    }
   
  //free some data structures and we are done

  if(line)
    rax_free(line);   
  
  fclose(f);

  rax_free(subtreeTips);

  //free the branch length meta-data 
  
  for(i = 0; i < tr->numberOfBranches; i++)
    {            
      rax_free(tr->bInf[i].epa->branches);
      rax_free(tr->bInf[i].epa->distalBranches);         
      rax_free(tr->bInf[i].epa->likelihoods);      
      rax_free(tr->bInf[i].epa);                
    } 

  rax_free(tr->bInf);
  
  printBothOpen("\n\nRAxML computed placements for %d subtrees in the given tree and will exit now, na zdorovja!\n\n", subTreeCount);

  exit(0);
}
