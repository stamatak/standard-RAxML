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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 *  of taxa and mixed models". 
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


extern char run_id[128];
extern char workdir[1024];
extern double masterTime;

/*
  function below not very interesting, standard RAxML stuff 
*/

static double testInsertThorough(tree *tr, nodeptr r, nodeptr q, boolean useVector)
{
  double 
    result,           
    qz[NUM_BRANCHES],
    z[NUM_BRANCHES];
  
  nodeptr  
    x = q->back,
    s = r->back;
  
  int     
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
    
  hookup(r->next,       q, z, tr->numBranches);
  hookup(r->next->next, x, z, tr->numBranches);
  hookupDefault(r, s, tr->numBranches);      		     
    
  newviewGeneric(tr, r);	     
    
  localSmooth(tr, r, smoothings);
	  
  if(useVector)
    result = evaluateGenericVector(tr, r);
  else
    result = evaluateGeneric(tr, r);	 	       	  	   

  hookup(q, x, qz, tr->numBranches);
      
  r->next->next->back = r->next->back = (nodeptr) NULL; 

  return result;
}


/*
  structure to store likelihood and insertion node 
  for each window position 
*/

typedef struct 
{
  double lh;
  nodeptr p;
} positionData;


/*
  traverse the tree recursively and insert taxon into each position
*/

static void traverseBias(nodeptr p, nodeptr q, tree *tr, int *branchCounter, positionData *pd, int windowSize)
{
  double 
    sum = 0.0,
    result = 0.0;

  int 
    i;


  /* 
     compute the likelihood of inserting the tip attached to 
     p between q and q->back 

     Actually in testInsertThorough() we compute the per site 
     log likes which are then stored in an array tr->perSiteLL[i]
  */
  
  result = testInsertThorough(tr, p, q, TRUE); 


  /* 
     stuff below can be removed at some point 
     it just makes me feel better
  */

  for(i = 0; i < tr->cdta->endsite; i++)
    sum += tr->perSiteLL[i];

  assert(ABS(sum - result) < 0.001);

  /*************************************/

  /* 
     for each window position just compute the likelihood over 
     the window.

     If its is better than the current best one, store the likelihood
     and the insertion node in the data structure 
  */

  for(i = 0; i < tr->cdta->endsite - windowSize; i++)
    {     
      int 
	j;

      for(j = i, sum = 0.0; j < i + windowSize; j++)
	sum += tr->perSiteLL[j];

      if(sum > pd[i].lh)
	{      
	  pd[i].lh = sum;
	  pd[i].p = q;
	}	    	
    }

  *branchCounter = *branchCounter + 1;
  

  /* and here comes the recursion */
 
  if(!isTip(q->number, tr->rdta->numsp))
    {    
      traverseBias(p, q->next->back,       tr, branchCounter, pd, windowSize);
      traverseBias(p, q->next->next->back, tr, branchCounter, pd, windowSize);    
    }    
}


/*
  functions to compute the node distances between inferred and true 
  placement positions. I think that they are correct.
*/

static int findRec(nodeptr ref, nodeptr best, int mxtips, int value)
{
  if(isTip(ref->number, mxtips))
    {
      if(ref == best || ref->back == best)
	return value;
      else
	return 0;
    }
  else
    {
      int 
	d1,
	d2;

      if(ref == best || ref->back == best)
	return value;

      d1 = findRec(ref->next->back, best, mxtips, value + 1);
      d2 = findRec(ref->next->next->back, best, mxtips, value + 1);
      
      assert((d1 > 0 && d2 == 0) || (d2 > 0 && d1 == 0) || (d1 == 0 && d2 == 0)); 
      
      return (d1 + d2);
    }
}

static double getNodeDistance(nodeptr ref, nodeptr best, int mxtips)
{  
  int 
    d1 = findRec(ref, best, mxtips, 0),
    d2 = findRec(ref->back, best, mxtips, 0);

  assert((d1 > 0 && d2 == 0) || (d2 > 0 && d1 == 0) || (d1 == 0 && d2 == 0));    

  return ((double)(d1 + d2));
}

void computePlacementBias(tree *tr, analdef *adef)
{
  int 
    windowSize = adef->slidingWindowSize,
    k,   
    i, 
    tips,
    numTraversalBranches = (2 * (tr->mxtips - 1)) - 3; /* compute number of branches into which we need to insert once we have removed a taxon */    
    
  char 
    fileName[1024];

  FILE 
    *outFile;

  /* data for each sliding window starting position */

  positionData
    *pd = (positionData *)malloc(sizeof(positionData) * (tr->cdta->endsite - windowSize));

  double 
    *nodeDistances = (double*)calloc(tr->cdta->endsite, sizeof(double)), /* array to store node distnces ND for every sliding window position */
    *distances = (double*)calloc(tr->cdta->endsite, sizeof(double)); /* array to store avg distances for every site */

  strcpy(fileName,         workdir);
  strcat(fileName,         "RAxML_SiteSpecificPlacementBias.");
  strcat(fileName,         run_id);

  outFile = myfopen(fileName, "w");

  printBothOpen("Likelihood of comprehensive tree %f\n\n", tr->likelihood);

  if(windowSize > tr->cdta->endsite)
    {
      printBothOpen("The size of your sliding window is %d while the number of sites in the alignment is %d\n\n", windowSize, tr->cdta->endsite);
      exit(-1);
    }

  if(windowSize >=  (int)(0.9 * tr->cdta->endsite))    
    printBothOpen("WARNING: your sliding window of size %d is only slightly smaller than you alignment that has %d sites\n\n",  windowSize, tr->cdta->endsite);

  printBothOpen("Sliding window size: %d\n\n", windowSize);

  /* prune and re-insert on tip at a time into all branches of the remaining tree */

  for(tips = 1; tips <= tr->mxtips; tips++)    
    {
      nodeptr 
	myStart,
	p = tr->nodep[tips]->back, /* this is the node at which we are prunung */
	p1 =  p->next->back,
	p2 =  p->next->next->back;
      
      double
	pz[NUM_BRANCHES],
	p1z[NUM_BRANCHES], 
	p2z[NUM_BRANCHES];

      int 
	branchCounter = 0;
      
      /* reset array values for this tip */

      for(i = 0; i < tr->cdta->endsite; i++)
	{
	  pd[i].lh = unlikely;
	  pd[i].p = (nodeptr)NULL;
	}

      /* store the three branch lengths adjacent to the position at which we prune */

      for(i = 0; i < tr->numBranches; i++)
	{
	  p1z[i] = p1->z[i];
	  p2z[i] = p2->z[i];	   	   
	  pz[i] = p->z[i];
	}           

      /* prune the taxon, optimizing the branch between p1 and p2 */

      removeNodeBIG(tr, p,  tr->numBranches);  

      printBothOpen("Pruning taxon Number %d [%s]\n", tips, tr->nameList[tips]);
      
      /* find any tip to start traversing the tree */

      myStart = findAnyTip(p1, tr->mxtips);      

      /* insert taxon, compute likelihood and remove taxon again from all branches */

      traverseBias(p, myStart->back, tr, &branchCounter, pd, windowSize);     

      assert(branchCounter == numTraversalBranches);                    

      /* for every sliding window position calc ND to the true/correct position at p */

      for(i = 0; i < tr->cdta->endsite - windowSize; i++)		
	nodeDistances[i] = getNodeDistance(p1, pd[i].p, tr->mxtips); 

      /* now analyze */

      for(i = 0; i < tr->cdta->endsite; i++)
	{
	  double 
	    d = 0.0;

	  int 
	    s = 0;	  	  
	  
	  /* 
	     check site position, i.e., doe we have windowSize data points available 
	     or fewer because we are at the start or the end of the alignment 
	  */

	  /*
	    for each site just accumulate the node distances we have for all 
	    sliding windows that passed over this site 
	  */

	  if(i < windowSize)
	    {
	      for(k = 0; k < i + 1; k++, s++)	    		
		d += nodeDistances[k];		 		
	    }
	  else
	    {
	      if(i < tr->cdta->endsite - windowSize)
		{
		  for(k = i - windowSize + 1; k <= i; k++, s++)	    		    
		    d += nodeDistances[k];		      		   
		}	    
	      else
		{				  
		  for(k = i - windowSize; k < (tr->cdta->endsite - windowSize); k++, s++)	    		    
		    d += nodeDistances[k + 1];		     		 
		}
	    }
	   	  
	  /* 
	     now just divide the accumultaed ND distance by 
	     the number of distances we have for this position and then add it to the acc 
	     distances over all taxa.
	     I just realized that the version on which I did the tests 
	     I sent to Simon I used 
	     
	     distances[i] = d / ((double)s);

	     instead of 
	     
	     distances[i] += d / ((double)s);

	     gamo tin poutana mou
	  */
	     
	  distances[i] += (d / ((double)s));	  
	}
      

      /* 
	 re-connect taxon to its original position 
      */

      hookup(p->next,       p1,      p1z, tr->numBranches); 
      hookup(p->next->next, p2,      p2z, tr->numBranches);
      hookup(p,             p->back, pz, tr->numBranches);      

      /*
	fix likelihood vectors 
      */

      newviewGeneric(tr, p);          
    }

  /* 
     now just compute the average ND over all taxa 
  */
    

  for(i = 0; i < tr->cdta->endsite; i++)
    {
      double 
	avg = distances[i] / ((double)tr->mxtips);

      fprintf(outFile, "%d %f\n", i, avg);
    }

  printBothOpen("\nTime for EPA-based site-specific placement bias calculation: %f\n", gettime() - masterTime); 
  printBothOpen("Site-specific placement bias statistics written to file %s\n", fileName);

  fclose(outFile);

  exit(0);
}
