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

extern int Thorough;
extern infoList iList;
extern char seq_file[1024];











void computeBOOTRAPID (tree *tr, analdef *adef, long *radiusSeed) 
{ 
  int i, bestTrav, impr;                
  double lh, previousLh, difference, epsilon;              
  bestlist *bestT, *bt;  
  int countIT;
  
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
  saveBestTree(bestT, tr);

  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;  
  initBestTree(bt, 5, tr->mxtips);

  initInfoList(10);
 
  difference = 10.0;
  epsilon = 0.01;         

  bestTrav = adef->bestTrav = 5 + 11 * randum(radiusSeed);
   
  Thorough = 1;  
 
  impr = 1;

  if(tr->doCutoff)
    tr->itCount = 0;

  tr->bigCutoff = TRUE;
  
  for(countIT = 0; countIT < 2 && impr; countIT++) 
    {              
      recallBestTree(bestT, 1, tr);       
      treeEvaluate(tr, 1);	 	                    
      saveBestTree(bestT, tr); 
           
      lh = previousLh = tr->likelihood;
         
      treeOptimizeRapid(tr, 1, bestTrav, adef, bt);       

      impr = 0;
	  
      for(i = 1; i <= bt->nvalid; i++)
	{	    		  	   
	  recallBestTree(bt, i, tr);	    	  	  
	 

	  treeEvaluate(tr, 0.25);	  	  
	      
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	       	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}

    } 
  
  tr->bigCutoff = FALSE;

  recallBestTree(bestT, 1, tr);   
  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList(); 
}







void optimizeRAPID(tree *tr, analdef *adef) 
{ 
  int i,  impr, bestTrav;
   
  double lh, previousLh, difference, epsilon;              
  bestlist *bestT, *bt;   
  
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  initInfoList(50);
 
  difference = 10.0;
  epsilon = 0.01;    
    
  Thorough = 0; 

 
         
  saveBestTree(bestT, tr);  
  bestTrav = adef->bestTrav = determineRearrangementSetting(tr, adef, bestT, bt);                   
  saveBestTree(bestT, tr); 

  impr = 1;
  if(tr->doCutoff)
    tr->itCount = 0;

  while(impr)
    {              
      recallBestTree(bestT, 1, tr);     
      treeEvaluate(tr, 1);


      saveBestTree(bestT, tr);     
         
      lh = previousLh = tr->likelihood;
         
      treeOptimizeRapid(tr, 1, bestTrav, adef, bt);   
      
      impr = 0;
	  
      for(i = 1; i <= bt->nvalid; i++)
	{	    		  	   
	  recallBestTree(bt, i, tr);	    
	  treeEvaluate(tr, 0.25);	    	 	

	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	       	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}          
    }

  recallBestTree(bestT, 1, tr);
  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList(); 
}



void thoroughOptimization(tree *tr, analdef *adef, topolRELL_LIST *rl, int index) 
{ 
  int i, impr;                
  int rearrangementsMin = 1, rearrangementsMax = adef->stepwidth;
  double lh, previousLh, difference, epsilon;              
  bestlist *bestT, *bt;  
    
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;   
  initBestTree(bt, 20, tr->mxtips);

  initInfoList(50);
 
  difference = 10.0;
  epsilon = 0.01;       
 
  
  saveBestTree(bestT, tr);  
 
  impr = 1;
  if(tr->doCutoff)
    tr->itCount = 0;
  
  Thorough = 1;
  impr = 1;

  while(1)
    {	     	
      recallBestTree(bestT, 1, tr);    
      if(impr)
	{	    	
	  rearrangementsMin = 1;
	  rearrangementsMax = adef->stepwidth;	    	  
	}			  			
      else
	{		       	   
	  rearrangementsMax += adef->stepwidth;
	  rearrangementsMin += adef->stepwidth; 	        	      
	  if(rearrangementsMax > adef->max_rearrange)	     	     	 
	    goto cleanup; 	   
	}
          
      treeEvaluate(tr, 1.0);	      
      previousLh = lh = tr->likelihood;	      
      saveBestTree(bestT, tr);                           
      
      treeOptimizeRapid(tr, rearrangementsMin, rearrangementsMax, adef, bt);
      
      impr = 0;			      	
		
      for(i = 1; i <= bt->nvalid; i++)
	{	
	  recallBestTree(bt, i, tr);	    
	 	  
	  treeEvaluate(tr, 0.25);	    	 
	  
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	  	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}

	
    }

 cleanup:  
  saveTL(rl, tr, index);  
  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList();
}




/*********************************************************************************************************************/










static boolean qupdate (tree *tr, nodeptr p)
{    
  nodeptr  q;
  double   z0[NUM_BRANCHES], z[NUM_BRANCHES];	
  int i;
      
  q = p->back;
  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];
      
  makenewzGeneric(tr, p, q, z0, 1, z, FALSE); 	 
      
  for(i = 0; i < tr->numBranches; i++)
    p->z[i] = q->z[i] = z[i]; 
        
  return TRUE;    
} 





static boolean qsmoothLocal(tree *tr, nodeptr p, int n)
{
  nodeptr  q;
  
  if(n == 0)
    return TRUE;
  else
    {
      if (! qupdate(tr, p))               return FALSE; /*  Adjust branch */
      if (!isTip(p->number, tr->rdta->numsp)) 
	{                                  /*  Adjust descendants */
	  q = p->next;
	  while (q != p) 
	    {
	      if (! qsmoothLocal(tr, q->back, n - 1))   return FALSE;
	      q = q->next;
	    }	
	  
	  newviewGeneric(tr, p);
	}
      
      return TRUE;
    }
} 

static void quickSmoothLocal(tree *tr, int n)
{
  nodeptr p = tr->insertNode;
  nodeptr q;

  if(n == 0)
    {
      evaluateGeneric(tr, p);
    }
  else
    {
      qsmoothLocal(tr, p->back, n - 1);
      if(!isTip(p->number, tr->rdta->numsp))
	{
	  q = p->next;
	  while(q != p)
	    {
	      qsmoothLocal(tr, q->back, n - 1);
	      q = q->next;
	    }
	}
      evaluateGeneric(tr, p);
    }
}




int treeOptimizeThorough(tree *tr, int mintrav, int maxtrav)
{
  int i;   
  bestlist *bestT;  

  nodeRectifier(tr);

  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
  
  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3;      
 
  tr->startLH = tr->endLH = tr->likelihood;  

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {     

      
      tr->bestOfNode = unlikely;     
      if(rearrangeBIG(tr, tr->nodep[i], mintrav, maxtrav))
	{          
	 
	  if((tr->endLH > tr->startLH) && (tr->bestOfNode != unlikely))
	    {			    
	      restoreTreeFast(tr);	     
	      quickSmoothLocal(tr, 3);
	      tr->startLH = tr->endLH = tr->likelihood;	 		     
	    }	 
	  else
	    {		 
	      if(tr->bestOfNode != unlikely)
		{		     
		  resetBestTree(bestT);		  		  		  
		  saveBestTree(bestT, tr);
		  restoreTreeFast(tr);		  
		  quickSmoothLocal(tr, 3);		  		    
		  if(tr->likelihood < tr->startLH)		    		    
		    {
		      int res;
		      res = recallBestTree(bestT, 1, tr);		      		    
		      assert(res > 0);
		    }
		  else				    
		    tr->startLH = tr->endLH = tr->likelihood;		  
		}
	    }
	}

    	
    }    

  freeBestTree(bestT);
  free(bestT);

  return 1;     
}
