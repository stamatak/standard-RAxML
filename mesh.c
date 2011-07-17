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

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"

extern const unsigned int mask32[32];
extern int Thorough;
extern double masterTime;
extern char run_id[128];
extern char  workdir[1024];

int terraceCounter = 0;


static void findNextRec(nodeptr p, tree *tr, int *count, nodeptr *nodes, int depth)
{
  int 
    model,
    notFound;

  if(isTip(p->number, 
	   tr->mxtips))
    {
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(!nodes[model] && (p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]))
	    {
	      assert(p->backs[model]);
	      nodes[model] = p;
	      count[model] = depth;
	    }
	}
    }
  else
    {
      assert(p == p->next->next->next);

      for(model = 0, notFound = 0; model < tr->NumberOfModels; model++)
	{
	  if(!nodes[model])
	    {
	      if(p->backs[model])
		{
		  nodes[model] = p;
		  count[model] = depth;	       
		}
	      else
		notFound++;
	    }
	  
	}

      if(notFound == 0)
	return;

      findNextRec(p->next->back, tr, count, nodes, depth + 1);
      findNextRec(p->next->next->back, tr, count, nodes, depth + 1);
    }

}

void findNext(nodeptr p, tree *tr, nodeptr *result)
{
  int 
    depth,
    model,
    leftDepth[NUM_BRANCHES],
    rightDepth[NUM_BRANCHES];
  
  nodeptr
    leftNodes[NUM_BRANCHES],
    rightNodes[NUM_BRANCHES];

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      leftDepth[model]  = -1;
      rightDepth[model] = -1;

      result[model]     = (nodeptr)NULL;      
      leftNodes[model]  = (nodeptr)NULL;
      rightNodes[model] = (nodeptr)NULL;
    }

  depth = 0;
  findNextRec(p, tr, leftDepth, leftNodes, 0);

  depth = 0; 
  findNextRec(p->back, tr, rightDepth, rightNodes, 0);


  for(model = 0; model < tr->NumberOfModels; model++)
    {
      assert(rightDepth[model] >= 0 || leftDepth[model] >= 0);
      assert(rightNodes[model] || leftNodes[model]);

      if(rightDepth[model] >= 0 && leftDepth[model] >= 0)
	{
	  if(rightDepth[model] < leftDepth[model])
	    result[model] = rightNodes[model];
	  else 
	    result[model] = leftNodes[model];
	}
      else
	{
	  if(rightDepth[model] >= 0)
	    result[model] = rightNodes[model];
	  else
	    result[model] = leftNodes[model];
	}
	 
      assert(result[model]->backs[model] && result[model]->next->backs[model] && result[model]->next->next->backs[model]);
      /* printf("Partition %d node %d\n", model, result[model]->number); */
    
    }
}


static int containsModel(nodeptr p, int model, tree *tr)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH])
	return 1;
      else
	return 0;
    }
  else
    {
      assert(p = p->next->next->next);
      return (containsModel(p->next->back, model, tr) || containsModel(p->next->next->back, model, tr));
    }
}

static void reduceTreeModelREC(nodeptr p, nodeptr reference, int model, tree *tr, int *branchCounter)
{

  if(isTip(p->number, tr->mxtips))
    {
      assert(p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);

      p->backs[model]         = reference;
      reference->backs[model] = p;
      p->z[model] = reference->z[model] = defaultz;

      *branchCounter = *branchCounter + 1;

    }
  else
    {
      nodeptr q = p->next;
      nodeptr r = p->next->next;
      
      int left  = containsModel(q->back, model, tr);      
      int right = containsModel(r->back, model, tr);

      assert(p = p->next->next->next);

      if(left && right)
	{
	  p->backs[model]         = reference;
	  reference->backs[model] = p;
	  p->z[model] = reference->z[model] = defaultz;

	  *branchCounter = *branchCounter + 1;

	  reduceTreeModelREC(q->back, q,  model, tr, branchCounter);   	 	
	  reduceTreeModelREC(r->back, r,  model, tr, branchCounter);
	} 
      else
	{
	  if(left || right)
	    {
	      if(left)		
		reduceTreeModelREC(q->back, reference,  model, tr, branchCounter);
	      else	
		reduceTreeModelREC(r->back, reference,  model, tr, branchCounter);		 
	    }
	  else	    
	    assert(0);	      	    
	}     
    }
}





static void initravPresence(nodeptr p, int ntips)
{
  if(!(isTip(p->number, ntips)))
    {
      int i;

      nodeptr 
	l = p->next->back,
	r = p->next->next->back;

      assert(p->next->next->next == p);

      initravPresence(l, ntips);
      initravPresence(r, ntips);
	
      for(i = 0; i < VECTOR_LENGTH; i++)
	p->isPresent[i] = l->isPresent[i] | r->isPresent[i];
    }

}


static int sanityCheckRec(tree *tr, nodeptr p, int model)
{
  
  if(isTip(p->number, tr->ntips))
    {
      assert(p->backs[model]);
      assert(p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);     
      return 1;
    }
  else
    {          
      assert(p->backs[model] && p->next->backs[model] && p->next->next->backs[model]);

      return 
	(sanityCheckRec(tr, p->next->backs[model], model) + sanityCheckRec(tr, p->next->next->backs[model], model));
      
    }


}








void setupPointerMesh(tree *tr)
{
  if(tr->multiGene)
    {
      int 
	model;
      
      assert(isTip(tr->start->number, tr->mxtips));
      
      /*printf("Global start at tip %d\n", tr->start->number);*/      
      /*initravPresence(tr->start->back, tr->mxtips);*/
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{ 
	  int 
	    branchCounter = 0;
	  
	  assert(isTip(tr->startVector[model]->number, tr->mxtips));
	  
	  reduceTreeModelREC(tr->startVector[model]->back, tr->startVector[model], model, tr, &branchCounter);
	  
	  /*printf("Partition %d has %d branches\n", model, 2 * tr->mxtipsVector[model] - 3);*/
	  
	  assert(branchCounter == 2 * tr->mxtipsVector[model] - 3);	  
	}            
    }
}




static nodeptr getPartition(int model, tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH])
	{
	  assert(p->backs[model]);	  	  	  
	  return p;
	} 
      else
	return ((nodeptr)NULL);
    }
  else
    {
      if(p->backs[model])
	{
	  assert(p->next->backs[model] && p->next->next->backs[model]);	
	  return p;
	}
      else
	{	 
	  nodeptr 
	    r1 = (nodeptr)NULL, 
	    r2 = (nodeptr)NULL;

	  r1 = getPartition(model, tr, p->next->back); 	  	 
	  r2 = getPartition(model, tr, p->next->next->back);
	  
	  assert(!(r1 && r2));

	  if(r1)
	    {
	      assert(r1->backs[model] != p->next);
	      return r1;
	    }
	  else
	    {
	      if(r2)
		{
		  assert(r2->backs[model] != p->next->next);
		  return r2;
		}
	      else
		return ((nodeptr)NULL);
	    }
	}

    }

}


static nodeptr getPartitionOneSide(int model, tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->mxtips))
    {      
      if((p->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]))
	{
	  assert(p->backs[model]);	  	  	  
	  return p;
	} 
      else
	return ((nodeptr)NULL);      
    }
  else
    {
      if(p->backs[model])
	{	 
	  assert(p->next->backs[model] && p->next->next->backs[model]);	
	  return p;
	}
      else
	{	 
	  nodeptr 
	    r1 = (nodeptr)NULL, 
	    r2 = (nodeptr)NULL;

	  r1 = getPartitionOneSide(model, tr, p->next->back); 	  	 
	  r2 = getPartitionOneSide(model, tr, p->next->next->back);
	  
	  if(r1 && r2)
	    {
	      assert(r1->backs[model] = r2);
	      assert(r2->backs[model] = r1);
	    }

	  if(r1)
	    return r1;
	  if(r2)
	    return r2;

	  return ((nodeptr)NULL);
	}

    }

}


static void hookup_MULTI(nodeptr p, nodeptr q)
{
  p->back = q;
  q->back = p;
}

static void hookupMESH(nodeptr p,  nodeptr q, double z, int model)
{
  p->backs[model] = q;
  q->backs[model] = p;
  q->z[model] = p->z[model] = z;
}








static void removeNodeBIG_MULTI(tree *tr, nodeptr p)
{  
  double   zqr[NUM_BRANCHES], result[NUM_BRANCHES];
  nodeptr qGlobal, rGlobal;
  int model;
        
  /* q---p---r  */
  /*     |      */
  /*   p->back  */

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      int 
	down = 0,
	left = 0, 
	right = 0;

      nodeptr 
	downPtr = (nodeptr)NULL;

      tr->removeNodes[model] = (nodeptr)NULL;
      tr->leftNodes[model]   = (nodeptr)NULL;
      tr->rightNodes[model]  = (nodeptr)NULL;  
      tr->storedBacks[model]  = (nodeptr)NULL;

      if(p->backs[model])
	{	 
	  /* 
	     this is the easy case all departing branches of node 
	     p contain taxa for partition model 
	  */

	  assert(p->next->backs[model] && p->next->next->backs[model]);

	  /* root of the subtree to be removed */
	  tr->removeNodes[model] = p;	 
	  tr->leftNodes[model]   = p->next->backs[model];
	  tr->rightNodes[model]  = p->next->next->backs[model];
	  
	  if(isTip(tr->leftNodes[model]->number, tr->mxtips) && isTip(tr->rightNodes[model]->number, tr->mxtips))
	    tr->removeNodes[model] = (nodeptr)NULL;

	  assert(p->backs[model]->backs[model] == p);
	  
	  /*
	    {	   
	    
	    int 
	    down = containsModel(p->back, model, tr),
	    left = containsModel(p->next->back, model, tr),
	    right = containsModel(p->next->next->back, model, tr);
	    
	    assert(down && left && right);
	    }
	  */
	}
      else
	{ 	 
	  assert(!(p->next->backs[model]) && !(p->next->next->backs[model]));
	  down  = containsModel(p->back, model, tr);
	  left  = containsModel(p->next->back, model, tr);
	  right = containsModel(p->next->next->back, model, tr);
	 	  	  
	  if((left || right) && down)
	    {	      
	      downPtr = getPartition(model, tr, p->back);	 
	     	  
	      assert(downPtr);    	     
	      assert(downPtr->next->backs[model] && downPtr->next->next->backs[model]);
	      assert(downPtr->backs[model]->backs[model] == downPtr);

	      if(((!left && right) || (!right && left)))	      	      
		{
		  if(!isTip(downPtr->backs[model]->number, tr->mxtips))		
		    {
		      /* here we need to go back one step */
		      tr->removeNodes[model] = downPtr->backs[model];		
		      tr->leftNodes[model]   = downPtr->backs[model]->next->backs[model];
		      tr->rightNodes[model]  = downPtr->backs[model]->next->next->backs[model];
		      
		      if(isTip(tr->leftNodes[model]->number, tr->mxtips) && isTip(tr->rightNodes[model]->number, tr->mxtips))
			tr->removeNodes[model] = (nodeptr)NULL;

		      assert(downPtr->backs[model]->backs[model] = downPtr);
		      assert(downPtr->backs[model]->z[model] == downPtr->z[model]);
		    }
		}
	      else
		{		 
		  assert(0);		 
		}
	    }	  	  		 
	}	  
      

      if(tr->removeNodes[model])
	{
	  int i;
	  
	  tr->zDown[model]  = tr->removeNodes[model]->z[model];
	  tr->zLeft[model]  = tr->leftNodes[model]->z[model];
	  tr->zRight[model] = tr->rightNodes[model]->z[model];
	  
	  assert(!isTip(tr->removeNodes[model]->number, tr->mxtips));
	  assert(!(isTip(tr->leftNodes[model]->number, tr->mxtips) && isTip(tr->rightNodes[model]->number, tr->mxtips)));

	  for(i = 0; i < tr->NumberOfModels; i++)
	    tr->partitionConverged[i] = TRUE;

	  tr->partitionConverged[model] = FALSE;

	  zqr[model] = tr->leftNodes[model]->z[model] * tr->rightNodes[model]->z[model];

	  hookupMESH(tr->leftNodes[model], tr->rightNodes[model], zqr[model], model);	 
	  
	  makenewzGeneric(tr, tr->leftNodes[model], tr->rightNodes[model], zqr, iterations, result, TRUE);      
	  
	  hookupMESH(tr->leftNodes[model], tr->rightNodes[model], result[model], model);
	  
	  tr->removeNodes[model]->next->backs[model] = (nodeptr)NULL;
	  tr->removeNodes[model]->next->next->backs[model] = (nodeptr)NULL;
	  	 	  
	  tr->storedBacks[model] = tr->removeNodes[model]->backs[model];

	  tr->removeNodes[model]->backs[model] = (nodeptr)NULL;
	}
    }

 

  qGlobal = p->next->back;
  rGlobal = p->next->next->back;

  hookup_MULTI(qGlobal, rGlobal);
    
  p->next->next->back = p->next->back = (nodeptr)NULL;  
}


/************************************************************************/


static void resetLikelihoodList(tree *tr)
{
  int model;
  
  for(model = 0; model < tr->NumberOfModels; model++)
    tr->likelihoodList[model]->count = 0; 
}

static void insertLikelihoodList(tree *tr, int model, int left, int right, double likelihood)
{
  lhList *list = tr->likelihoodList[model];

  if(list->count == list->size)
    {
      list->size *= 2;
      list->entries = realloc(list->entries, list->size * sizeof(lhEntry));
      assert(list->entries);
    }

  list->entries[list->count].left = left;
  list->entries[list->count].right = right;
  list->entries[list->count].likelihood = likelihood;
  
  list->count = list->count + 1;
}

static double lookupLikelihoodList(tree *tr, int model, int left, int right)
{
  lhList *list = tr->likelihoodList[model];

  int i;

  for(i = 0; i < list->count; i++)
    {
      lhEntry *entry = &(list->entries[i]);
      if((entry->left == left && entry->right == right) || (entry->left == right && entry->right == left))
	return entry->likelihood;
    }

  return 0.0;    
}

static void testInsertBIG_MULTI (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  
    r;
  
  double 
    likelihood = 0.0;
  
  int
    evalCount = 0,
    model;
 
  r = q->back;                   

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      if(tr->removeNodes[model])
	{
	  int left = containsModel(q, model, tr);
	  int right = containsModel(r, model, tr);
	      
	  if(left || right)
	    {	    
	      nodeptr 
		qModel = (nodeptr)NULL,
		rModel = (nodeptr)NULL;

	      if(left && right)
		{	      
		  qModel = getPartition(model, tr, q);
		  rModel = getPartition(model, tr, r);	 	 	
		}
	      else
		{
		  if(left)
		    {
		      qModel = getPartitionOneSide(model, tr, q);
		      rModel = qModel->backs[model];
		    }
		  else
		    {
		      assert(right);
		      rModel = getPartitionOneSide(model, tr, r);
		      qModel = rModel->backs[model];
		    }
		  assert(!(isTip(rModel->number, tr->mxtips) && isTip(qModel->number, tr->mxtips)));

		}
	
	      assert(!isTip(tr->removeNodes[model]->number, tr->mxtips));	  	  
	  
	      assert(qModel && rModel);
	  	 	 	  	            
	      if(!isTip(rModel->number, tr->mxtips))
		assert(rModel->next->backs[model] && rModel->next->next->backs[model]);
	      else	      	   
		assert(rModel->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);
		  
	      if(!isTip(qModel->number, tr->mxtips))
		assert(qModel->next->backs[model] && qModel->next->next->backs[model]);
	      else	     
		assert(qModel->isPresent[model / MASK_LENGTH] & mask32[model % MASK_LENGTH]);	   
	      	       
	      assert(qModel->backs[model] == rModel);
	      assert(rModel->backs[model] == qModel);
	      
	      if(!((rModel == tr->leftNodes[model] && qModel == tr->rightNodes[model]) ||
		   (rModel == tr->rightNodes[model] && qModel == tr->leftNodes[model])))
		{
		  double storedLH = lookupLikelihoodList(tr, model, rModel->number, qModel->number);

		  if(storedLH == 0.0)
		    {		    
		      double z, znew, partLH, qz[NUM_BRANCHES];
		      
		      assert(rModel->backs[model] == qModel);
		      assert(qModel->backs[model] == rModel);
		      assert(qModel->z[model] == rModel->z[model]);
		      z = qModel->z[model];
		      qz[model] = z;

		      znew = sqrt(z);
		      if(znew < zmin)
			znew = zmin;
		      if(znew > zmax)
			znew = zmax;
		      
		      hookupMESH(qModel, tr->removeNodes[model]->next, znew, model);	 
		      hookupMESH(rModel, tr->removeNodes[model]->next->next, znew, model);
		      hookupMESH(tr->removeNodes[model], tr->storedBacks[model], tr->zDown[model], model);
		      
		      if(Thorough)
			{
			  int i;
			  
			  double  
			    zqr[NUM_BRANCHES], 
			    zqs[NUM_BRANCHES], 
			    zrs[NUM_BRANCHES],
			    e1[NUM_BRANCHES], 
			    e2[NUM_BRANCHES], 
			    e3[NUM_BRANCHES],
			    defaultArray[NUM_BRANCHES];

			  double 
			    lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;

			  for(i = 0; i < tr->NumberOfModels; i++)
			    {
			      defaultArray[i] = defaultz;
			      tr->partitionConverged[i] = TRUE;
			    }

			  tr->partitionConverged[model] = FALSE;

			  makenewzGeneric(tr, qModel, rModel, qz, iterations, zqr, TRUE);           
			  makenewzGeneric(tr, qModel, tr->storedBacks[model], defaultArray, iterations, zqs, TRUE);                  
			  makenewzGeneric(tr, rModel, tr->storedBacks[model], defaultArray, iterations, zrs, TRUE);
			  			  
			  lzqr = (zqr[model] > zmin) ? log(zqr[model]) : log(zmin); 
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

			  hookupMESH(qModel, tr->removeNodes[model]->next, e1[model], model);	 
			  hookupMESH(rModel, tr->removeNodes[model]->next->next, e2[model], model);
			  hookupMESH(tr->removeNodes[model], tr->storedBacks[model], e3[model], model);

			  assert(!isTip(qModel->backs[model]->number, tr->mxtips));
			  newviewGenericMulti(tr, qModel->backs[model], model);

			  localSmoothMulti(tr, tr->removeNodes[model], smoothings, model);
			  partLH = evaluateGenericMulti(tr, qModel, model);
			  evalCount++;
			}
		      else
			{
			  assert(!isTip(qModel->backs[model]->number, tr->mxtips));
			  newviewGenericMulti(tr, qModel->backs[model], model);			  
		      
			  partLH = evaluateGenericMulti(tr, qModel, model);
			  evalCount++;
			}
		      
		      insertLikelihoodList(tr, model, rModel->number, qModel->number, partLH);
		      likelihood += partLH; 
		      
		      
		      hookupMESH(qModel, rModel, z, model);
		      
		      tr->removeNodes[model]->backs[model]             = (nodeptr)NULL;	
		      tr->removeNodes[model]->next->backs[model]       =  (nodeptr)NULL;
		      tr->removeNodes[model]->next->next->backs[model] = (nodeptr)NULL;		  		  		 
		    }
		  else
		    likelihood += storedLH;
		}
	      else
		likelihood += tr->storedPerPartitionLH[model];
	    }
	  else
	    likelihood += tr->storedPerPartitionLH[model];
	}
      else	   
	likelihood += tr->storedPerPartitionLH[model];	      	    
    }
  
  tr->likelihood = likelihood;
     
  if(tr->likelihood > tr->endLH)
    {			  
      tr->insertNode = q;
      tr->removeNode = p;               
      tr->endLH = tr->likelihood;
      printBothOpen("%f %d %d\n", likelihood, tr->insertNode->number,  tr->removeNode->number);
    }        

  if(evalCount == 0)
    terraceCounter++;

  return;   
}


static boolean testInsertNoMesh (tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;  
  int i;
  
  r = q->back; 
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
       
  if (! insertBIG(tr, p, q, tr->numBranches))       
    return FALSE;         
      
  evaluateGeneric(tr, p->next->next);           
      
  if(tr->likelihood > tr->endLH)
    {			  
      tr->insertNode = q;
      tr->removeNode = p;               
      tr->endLH = tr->likelihood;   
      printBothOpen("%f %d %d\n", tr->likelihood, tr->insertNode->number,  tr->removeNode->number);
    }         

  hookup(q, r, qz, tr->numBranches); 
      
  p->next->next->back = p->next->back = (nodeptr) NULL;

  if(Thorough)
    {
       nodeptr s = p->back;
       hookup(p, s, pz, tr->numBranches);
    }
  
  return TRUE;
}


/*********************************************************************************************************************************/

static void addTraverseBIG_MULTI(tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  
  if (--mintrav <= 0)     
    testInsertBIG_MULTI(tr, p, q);
  
  if ((!isTip(q->number, tr->rdta->numsp)) && (--maxtrav > 0)) 
    {    
      addTraverseBIG_MULTI(tr, p, q->next->back, mintrav, maxtrav);
      addTraverseBIG_MULTI(tr, p, q->next->next->back, mintrav, maxtrav);    
    }
} 


static void addTraverseNoMesh(tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
    {              
      if (! testInsertNoMesh(tr, p, q))  return;

    }
  
  if ((!isTip(q->number, tr->rdta->numsp)) && (--maxtrav > 0)) 
    {    
      addTraverseNoMesh(tr, p, q->next->back, mintrav, maxtrav);
      addTraverseNoMesh(tr, p, q->next->next->back, mintrav, maxtrav);    
    }
} 

/*********************************************************************************************************************************/



static int rearrangeBIG_MULTI(tree *tr, nodeptr pGlobal, int mintrav, int maxtrav)   
{   
  nodeptr  p1, p2, qGlobal, q1, q2;
  int      mintrav2, model;  
  
  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  qGlobal = pGlobal->back;  

  if (!isTip(pGlobal->number, tr->rdta->numsp)) 
    {     
      p1 = pGlobal->next->back;
      p2 = pGlobal->next->next->back;
           
      if(!isTip(p1->number, tr->rdta->numsp) || !isTip(p2->number, tr->rdta->numsp))
	{
	  resetLikelihoodList(tr);
	  removeNodeBIG_MULTI(tr, pGlobal);
	    
	  if(!isTip(p1->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG_MULTI(tr, pGlobal, p1->next->back,
				   mintrav, maxtrav);         
	      addTraverseBIG_MULTI(tr, pGlobal, p1->next->next->back,
				   mintrav, maxtrav);          
	    }
	  
	  if(!isTip(p2->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG_MULTI(tr, pGlobal, p2->next->back,
				   mintrav, maxtrav);
	      addTraverseBIG_MULTI(tr, pGlobal, p2->next->next->back,
				   mintrav, maxtrav);          
	    }	  	  	  	 	   	    	    	  
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->removeNodes[model])
		{	
		  
		  hookupMESH(tr->removeNodes[model], tr->storedBacks[model], tr->zDown[model], model);		
		  hookupMESH(tr->removeNodes[model]->next,       tr->leftNodes[model], tr->zLeft[model], model);
		  hookupMESH(tr->removeNodes[model]->next->next, tr->rightNodes[model], tr->zRight[model], model);		  		 
		  
		  newviewGenericMulti(tr, tr->removeNodes[model], model);	       
		}		
	    }	    	  	  
	  
	  hookup_MULTI(pGlobal->next,       p1); 
	  hookup_MULTI(pGlobal->next->next, p2);	  	  
	}
    }
  
  if (!isTip(qGlobal->number, tr->rdta->numsp) && maxtrav > 0) 
    {
      q1 = qGlobal->next->back;
      q2 = qGlobal->next->next->back;
      
          
      if (
	  (
	   ! isTip(q1->number, tr->rdta->numsp) && 
	   (! isTip(q1->next->back->number, tr->rdta->numsp) || ! isTip(q1->next->next->back->number, tr->rdta->numsp))
	   )
	  ||
	  (
	   ! isTip(q2->number, tr->rdta->numsp) && 
	   (! isTip(q2->next->back->number, tr->rdta->numsp) || ! isTip(q2->next->next->back->number, tr->rdta->numsp))
	   )
	  )
	{
	  resetLikelihoodList(tr);
	  removeNodeBIG_MULTI(tr, qGlobal);
	   
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if ( !isTip(q1->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG_MULTI(tr, qGlobal, q1->next->back,
				   mintrav2 , maxtrav);
	      addTraverseBIG_MULTI(tr, qGlobal, q1->next->next->back,
				   mintrav2 , maxtrav);         
	    }
	  
	  if ( ! isTip(q2->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG_MULTI(tr, qGlobal, q2->next->back,
				   mintrav2 , maxtrav);
	      addTraverseBIG_MULTI(tr, qGlobal, q2->next->next->back,
				   mintrav2 , maxtrav);          
	    }
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {
	      if(tr->removeNodes[model])
		{				  		  
		  hookupMESH(tr->removeNodes[model], tr->storedBacks[model], tr->zDown[model], model); 
		  hookupMESH(tr->removeNodes[model]->next,       tr->leftNodes[model], tr->zLeft[model], model);
		  hookupMESH(tr->removeNodes[model]->next->next, tr->rightNodes[model], tr->zRight[model], model);
		  		 
		  newviewGenericMulti(tr, tr->removeNodes[model], model);
		}	     	
	    }	    	
	  
	  hookup_MULTI(qGlobal->next,       q1); 
	  hookup_MULTI(qGlobal->next->next, q2);	    	 	 	 
	}
    }
  
  return  1;
} 


static int rearrangeNoMesh(tree *tr, nodeptr p, int mintrav, int maxtrav)   
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;   
  
  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  q = p->back;
  
   
  if (!isTip(p->number, tr->rdta->numsp)) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
     
      if(!isTip(p1->number, tr->rdta->numsp) || !isTip(p2->number, tr->rdta->numsp))
	{
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      p1z[i] = p1->z[i];
	      p2z[i] = p2->z[i];	   	   
	    }
	  
	  if (! removeNodeBIG(tr, p,  tr->numBranches)) return badRear;
	  
	  if (!isTip(p1->number, tr->rdta->numsp)) 
	    {
	      addTraverseNoMesh(tr, p, p1->next->back,
				mintrav, maxtrav);         
	      addTraverseNoMesh(tr, p, p1->next->next->back,
				mintrav, maxtrav);          
	    }
	  
	  if (!isTip(p2->number, tr->rdta->numsp)) 
	    {
	      addTraverseNoMesh(tr, p, p2->next->back,
				mintrav, maxtrav);
	      addTraverseNoMesh(tr, p, p2->next->next->back,
				mintrav, maxtrav);          
	    }
	  	  
	  hookup(p->next,       p1, p1z, tr->numBranches); 
	  hookup(p->next->next, p2, p2z, tr->numBranches);

	  newviewGeneric(tr, p);	   	    	    	
	}
    }  
  
  if (!isTip(q->number, tr->rdta->numsp) && maxtrav > 0) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;
           
      if (
	  (
	   ! isTip(q1->number, tr->rdta->numsp) && 
	   (! isTip(q1->next->back->number, tr->rdta->numsp) || ! isTip(q1->next->next->back->number, tr->rdta->numsp))
	   )
	  ||
	  (
	   ! isTip(q2->number, tr->rdta->numsp) && 
	   (! isTip(q2->next->back->number, tr->rdta->numsp) || ! isTip(q2->next->next->back->number, tr->rdta->numsp))
	   )
	  )
	{
	  
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      q1z[i] = q1->z[i];
	      q2z[i] = q2->z[i];
	    }
	  
	  if (! removeNodeBIG(tr, q, tr->numBranches)) return badRear;
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if (!isTip(q1->number, tr->rdta->numsp)) 
	    {
	      addTraverseNoMesh(tr, q, q1->next->back,
				mintrav2 , maxtrav);
	      addTraverseNoMesh(tr, q, q1->next->next->back,
				mintrav2 , maxtrav);         
	    }
	  
	  if (!isTip(q2->number, tr->rdta->numsp)) 
	    {
	      addTraverseNoMesh(tr, q, q2->next->back,
				mintrav2 , maxtrav);
	      addTraverseNoMesh(tr, q, q2->next->next->back,
				mintrav2 , maxtrav);          
	    }	   
	  
	  hookup(q->next,       q1, q1z, tr->numBranches); 
	  hookup(q->next->next, q2, q2z, tr->numBranches);
	  

	  newviewGeneric(tr, q);	 
	}
    } 
  
  return  1;
} 

/*********************************************************************************************************************************/

static void executeInsert(tree *tr)
{
  nodeptr 
    p = tr->removeNode,
    r = p->next->back,
    q = p->next->next->back;  
 
  hookup_MULTI(r, q);
   
  r = tr->insertNode;
  q = r->back;
 
  hookup_MULTI(q, p->next);
  hookup_MULTI(r, p->next->next); 
}



/*********************************************************************************************************************************/


static void treeOptimizeRapidMesh(tree *tr, int mintrav, int maxtrav, analdef *adef, char fileName[1024])
{
  int i; 
  FILE *f;

  nodeRectifier(tr);

  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3;     
 
  tr->startLH = tr->endLH = tr->likelihood;
   
  tr->endLH = tr->likelihood;
  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {   

      tr->bestOfNode = unlikely;              
          
      rearrangeBIG_MULTI(tr, tr->nodep[i], mintrav, maxtrav);     
    }     

  evaluateGenericInitrav(tr, tr->start);

  if(tr->insertNode && tr->removeNode)
    {
      printBothOpen("Done %f -> %f best insert at %d remove at %d\n", tr->likelihood, tr->endLH, tr->insertNode->number, tr->removeNode->number);
      executeInsert(tr);
    }
  else
    printBothOpen("SPR moves did not yield an imporved tree\n");
  
  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, 
	      FALSE, FALSE, adef, NO_BRANCHES, FALSE, FALSE);

  f = myfopen(fileName, "wb");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);
}

static void treeOptimizeRapidNoMesh(tree *tr, int mintrav, int maxtrav, analdef *adef, char fileName[1024])
{
  int i;   

  FILE *f;

  nodeRectifier(tr);

  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3;     
 
  tr->startLH = tr->endLH = tr->likelihood;
 
  printf("Enter %f\n", tr->likelihood);

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {           

      tr->bestOfNode = unlikely;                          
      rearrangeNoMesh(tr, tr->nodep[i], mintrav, maxtrav);
    }     

  evaluateGenericInitrav(tr, tr->start);
  
  if(tr->insertNode && tr->removeNode)
    {
      printBothOpen("Done %f -> %f best insert at %d remove at %d\n", tr->likelihood, tr->endLH, tr->insertNode->number, tr->removeNode->number);
      executeInsert(tr);
    }
  else
    printBothOpen("SPR moves did not yield an imporved tree\n");
 
  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, 
	      FALSE, FALSE, adef, NO_BRANCHES, FALSE, FALSE); 

  f = myfopen(fileName, "wb");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);  
}

/*********************************************************************************************************************************/






void meshTreeSearch(tree *tr, analdef *adef, int thorough)
{
  double 
    t, 
    lh = 0.0;

  int 
    model;

  char fileName[1024];

  strcpy(fileName,         workdir);
  strcat(fileName,         "RAxML_sprTree.");
  strcat(fileName,         run_id);

  Thorough = thorough;  

  t = gettime();
  modOpt(tr, adef, FALSE, 5.0, FALSE);
  printBothOpen("Model optimization time: %f %f\n", gettime() - t, tr->likelihood);

  t = gettime();
  
  evaluateGenericInitrav(tr, tr->start);

  for(model = 0; model < tr->NumberOfModels; model++)
    {            
      tr->likelihoodList[model] = (lhList *)malloc(sizeof(lhList));
      tr->likelihoodList[model]->count = 0;
      tr->likelihoodList[model]->size = 2048;
      tr->likelihoodList[model]->entries = (lhEntry*)malloc(2048 * sizeof(lhEntry));

      tr->storedPerPartitionLH[model] = tr->perPartitionLH[model];
      lh += tr->storedPerPartitionLH[model];
    } 
      
 

  if(tr->multiGene)
    {
      treeOptimizeRapidMesh(tr, 1, 10, adef, fileName);
      printBothOpen("Number of terrace moves encountered during SPR: %d\n", terraceCounter);
    }
  else
    treeOptimizeRapidNoMesh(tr, 1, 10, adef, fileName);

  printBothOpen("SPR time: %f\n", gettime() - t);

  printBothOpen("Total time: %f\n", gettime() - masterTime);

  exit(0);
}
