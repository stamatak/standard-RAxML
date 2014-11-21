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
extern char resultFileName[1024];
extern char tree_file[1024];
extern char workdir[1024];
extern char run_id[128];
extern FILE *INFILE;
extern double masterTime;



boolean initrav (tree *tr, nodeptr p)
{ 
  nodeptr  q;
  
  if (!isTip(p->number, tr->rdta->numsp)) 
    {      
      q = p->next;
      
      do 
	{	   
	  if (! initrav(tr, q->back))  return FALSE;		   
	  q = q->next;	
	} 
      while (q != p);  
      
      newviewGeneric(tr, p);
    }
  
  return TRUE;
} 









//#define _DEBUG_UPDATE

boolean update(tree *tr, nodeptr p)
{       
  nodeptr  q; 
  boolean smoothedPartitions[NUM_BRANCHES];
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  double _deltaz;

#ifdef _DEBUG_UPDATE
  double 
    startLH;

  evaluateGeneric(tr, p);

  startLH = tr->likelihood;
#endif

  q = p->back;   

  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];    

  if(tr->numBranches > 1)
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, TRUE);  
  else
    makenewzGeneric(tr, p, q, z0, newzpercycle, z, FALSE);
  
  for(i = 0; i < tr->numBranches; i++)    
    smoothedPartitions[i]  = tr->partitionSmoothed[i];
      
  for(i = 0; i < tr->numBranches; i++)
    {         
      if(!tr->partitionConverged[i])
	{	  
	    _deltaz = deltaz;
	    
	  if(ABS(z[i] - z0[i]) > _deltaz)  
	    {	    
	      /*
		#ifdef _DEBUG_UPDATE  
		printf("accept\n");
		#endif
	      */
	      smoothedPartitions[i] = FALSE;       
	    }	 
	  /*
	    #ifdef _DEBUG_UPDATE
	    else
	    printf("reject\n");
	    #endif
	  */
	  	  
	  p->z[i] = q->z[i] = z[i];	 
#ifdef _BASTIEN
	  //printf("update %f\n", tr->secondDerivative[i]);
	  p->secondDerivative[i] = q->secondDerivative[i] = tr->secondDerivative[i];
	  p->secondDerivativeValid[i] = q->secondDerivativeValid[i] = TRUE;
#endif
	}
    }
  
#ifdef _DEBUG_UPDATE
  evaluateGeneric(tr, p);

  if(tr->likelihood <= startLH)
    {
      if(fabs(tr->likelihood - startLH) > 0.01)
	{
	  printf("%f %f\n", startLH, tr->likelihood);	  
	  assert(0);             
	}
    }
#endif

  for(i = 0; i < tr->numBranches; i++)    
    tr->partitionSmoothed[i]  = smoothedPartitions[i];
  
  return TRUE;
}




boolean smooth (tree *tr, nodeptr p)
{
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */
  if (! isTip(p->number, tr->rdta->numsp)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  if (! smooth(tr, q->back))   return FALSE;
	  q = q->next;
	}	
      
      if(tr->multiBranch)		  
	newviewGenericMasked(tr, p);	
      else
	newviewGeneric(tr, p);     
    }
  
  return TRUE;
} 

static boolean allSmoothed(tree *tr)
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



boolean smoothTree (tree *tr, int maxtimes)
{
  nodeptr  p, q;   
  int i, count = 0;
   
  p = tr->start;
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  while (--maxtimes >= 0) 
    {    
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;		

      if (! smooth(tr, p->back))       return FALSE;
      if (!isTip(p->number, tr->rdta->numsp)) 
	{
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smooth(tr, q->back))   return FALSE;
	      q = q->next;
	    }
	}
         
      count++;

      if (allSmoothed(tr)) 
	break;      
    }

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;



  return TRUE;
} 



boolean localSmooth (tree *tr, nodeptr p, int maxtimes)
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
	  if (! update(tr, q)) return FALSE;
	  q = q->next;
        } 
      while (q != p);
      
      if (allSmoothed(tr)) 
	break;
    }

  for(i = 0; i < tr->numBranches; i++)
    {
      tr->partitionSmoothed[i] = FALSE; 
      tr->partitionConverged[i] = FALSE;
    }

  return TRUE;
}





static void resetInfoList(void)
{
  int i;

  iList.valid = 0;

  for(i = 0; i < iList.n; i++)    
    {
      iList.list[i].node = (nodeptr)NULL;
      iList.list[i].likelihood = unlikely;
    }    
}

void initInfoList(int n)
{
  int i;

  iList.n = n;
  iList.valid = 0;
  iList.list = (bestInfo *)rax_malloc(sizeof(bestInfo) * n);

  for(i = 0; i < n; i++)
    {
      iList.list[i].node = (nodeptr)NULL;
      iList.list[i].likelihood = unlikely;
    }
}

void freeInfoList(void)
{ 
  rax_free(iList.list);   
}


void insertInfoList(nodeptr node, double likelihood)
{
  int i;
  int min = 0;
  double min_l =  iList.list[0].likelihood;

  for(i = 1; i < iList.n; i++)
    {
      if(iList.list[i].likelihood < min_l)
	{
	  min = i;
	  min_l = iList.list[i].likelihood;
	}
    }

  if(likelihood > min_l)
    {
      iList.list[min].likelihood = likelihood;
      iList.list[min].node = node;
      iList.valid += 1;
    }

  if(iList.valid > iList.n)
    iList.valid = iList.n;
}


boolean smoothRegion (tree *tr, nodeptr p, int region)
{ 
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */

  if(region > 0)
    {
      if (!isTip(p->number, tr->rdta->numsp)) 
	{                                 
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smoothRegion(tr, q->back, --region))   return FALSE;
	      q = q->next;
	    }	
	  
	    newviewGeneric(tr, p);
	}
    }
  
  return TRUE;
}

boolean regionalSmooth (tree *tr, nodeptr p, int maxtimes, int region)
  {
    nodeptr  q;
    int i;

    if (isTip(p->number, tr->rdta->numsp)) return FALSE;            /* Should be an error */

    for(i = 0; i < tr->numBranches; i++)
      tr->partitionConverged[i] = FALSE;

    while (--maxtimes >= 0) 
      {	
	for(i = 0; i < tr->numBranches; i++)	  
	  tr->partitionSmoothed[i] = TRUE;
	  
	q = p;
	do 
	  {
	    if (! smoothRegion(tr, q, region)) return FALSE;
	    q = q->next;
	  } 
	while (q != p);
	
	if (allSmoothed(tr)) 
	  break;
      }

    for(i = 0; i < tr->numBranches; i++)
      tr->partitionSmoothed[i] = FALSE;
    for(i = 0; i < tr->numBranches; i++)
      tr->partitionConverged[i] = FALSE;
   
    return TRUE;
  } /* localSmooth */





nodeptr  removeNodeBIG (tree *tr, nodeptr p, int numBranches)
{  
  double   zqr[NUM_BRANCHES], result[NUM_BRANCHES];
  nodeptr  q, r;
  int i;
        
  q = p->next->back;
  r = p->next->next->back;
  
  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        
   
  makenewzGeneric(tr, q, r, zqr, iterations, result, FALSE);   

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 
      
  p->next->next->back = p->next->back = (node *) NULL;



  return  q; 
}

nodeptr  removeNodeRestoreBIG (tree *tr, nodeptr p)
{
  nodeptr  q, r;
        
  q = p->next->back;
  r = p->next->next->back;  

  newviewGeneric(tr, q);
  newviewGeneric(tr, r);
  
  hookup(q, r, tr->currentZQR, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
     
  return  q;
}


boolean insertBIG (tree *tr, nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r, s;
  int i;
  
  r = q->back;
  s = p->back;
      
  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];
  
#ifndef _HET
  if(Thorough)
    { 
      double  zqr[NUM_BRANCHES], zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;      
      double defaultArray[NUM_BRANCHES];	
      double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
      double *qz;
      
      qz = q->z;
      
      for(i = 0; i < numBranches; i++)
	defaultArray[i] = defaultz;
      
      makenewzGeneric(tr, q, r, qz, iterations, zqr, FALSE);           
      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs, FALSE);                  
      makenewzGeneric(tr, r, s, defaultArray, iterations, zrs, FALSE);
      
      
      for(i = 0; i < numBranches; i++)
	{
	  lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
	  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
	  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
	  lzsum = 0.5 * (lzqr + lzqs + lzrs);
	  
	  lzq = lzsum - lzrs;
	  lzr = lzsum - lzqs;
	  lzs = lzsum - lzqr;
	  lzmax = log(zmax);
	  
	  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
	  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
	  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
	  
	  e1[i] = exp(lzq);
	  e2[i] = exp(lzr);
	  e3[i] = exp(lzs);
	}
      hookup(p->next,       q, e1, numBranches);
      hookup(p->next->next, r, e2, numBranches);
      hookup(p,             s, e3, numBranches);      		  
    }
  else
#endif
    {       
      double  z[NUM_BRANCHES]; 
      
      for(i = 0; i < numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  
	  if(z[i] < zmin) 
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);	                         
    }
  
  newviewGeneric(tr, p);
  
  if(Thorough)
    {     
      localSmooth(tr, p, smoothings);   
      for(i = 0; i < numBranches; i++)
	{
	  tr->lzq[i] = p->next->z[i];
	  tr->lzr[i] = p->next->next->z[i];
	  tr->lzs[i] = p->z[i];            
	}
    }           
  
  return  TRUE;
}

boolean insertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  
  r = q->back;
  s = p->back;

  if(Thorough)
    {                        
      hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
      hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
      hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES];
      int i;
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  double zz;
	  zz = sqrt(q->z[i]);     
	  if(zz < zmin) 
	    zz = zmin;
	  if(zz > zmax)
	    zz = zmax;
  	  z[i] = zz;
	}

      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);
    }   
    
  newviewGeneric(tr, p);
       
  return  TRUE;
}


static void restoreTopologyOnly(tree *tr, bestlist *bt)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[NUM_BRANCHES], pz[NUM_BRANCHES], p1z[NUM_BRANCHES], p2z[NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;
      
  p1 = p->next->back;
  p2 = p->next->next->back;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      p1z[i] = p1->z[i];
      p2z[i] = p2->z[i];
    }
  
  hookup(p1, p2, tr->currentZQR, tr->numBranches);
  
  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];           
    }
  
  r = q->back;
  s = p->back;
  
  if(Thorough)
    {                        
      hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
      hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
      hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
    }
  else
    { 	
      double  z[NUM_BRANCHES];	
      for(i = 0; i < tr->numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  if(z[i] < zmin)
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);
    }     
  
  tr->likelihood = tr->bestOfNode;
    
  saveBestTree(bt, tr);
  
  tr->likelihood = currentLH;
  
  hookup(q, r, qz, tr->numBranches);
  
  p->next->next->back = p->next->back = (nodeptr) NULL;
  
  if(Thorough)    
    hookup(p, s, pz, tr->numBranches);          
      
  hookup(p->next,       p1, p1z, tr->numBranches); 
  hookup(p->next->next, p2, p2z, tr->numBranches);      
}


boolean testInsertBIG (tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = TRUE;
  double startLH = tr->endLH;
  int i;
  
  r = q->back; 
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
  
  if(tr->grouped)
    {
      int rNumber, qNumber, pNumber;
      
      doIt = FALSE;
      
      rNumber = tr->constraintVector[r->number];
      qNumber = tr->constraintVector[q->number];
      pNumber = tr->constraintVector[p->number];
      
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
      if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;         
      
      evaluateGeneric(tr, p->next->next);   
       
      if(tr->likelihood > tr->bestOfNode)
	{
	  tr->bestOfNode = tr->likelihood;
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      tr->currentZQR[i] = tr->zqr[i];           
	      tr->currentLZR[i] = tr->lzr[i];
	      tr->currentLZQ[i] = tr->lzq[i];
	      tr->currentLZS[i] = tr->lzs[i];      
	    }
	}
      
      if(tr->likelihood > tr->endLH)
	{			  
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    tr->currentZQR[i] = tr->zqr[i];      
	  tr->endLH = tr->likelihood;                      
	}        
      
      hookup(q, r, qz, tr->numBranches);
      
      p->next->next->back = p->next->back = (nodeptr) NULL;
      
      if(Thorough)
	{
	  nodeptr s = p->back;
	  hookup(p, s, pz, tr->numBranches);      
	} 
      
      if((tr->doCutoff) && (tr->likelihood < startLH))
	{
	  tr->lhAVG += (startLH - tr->likelihood);
	  tr->lhDEC++;
	  if((startLH - tr->likelihood) >= tr->lhCutoff)
	    return FALSE;	    
	  else
	    return TRUE;
	}
      else
	return TRUE;
    }
  else
    return TRUE;  
}



 
void addTraverseBIG(tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
    {              
      if (! testInsertBIG(tr, p, q))  return;

    }
  
  if ((!isTip(q->number, tr->rdta->numsp)) && (--maxtrav > 0)) 
    {    
      addTraverseBIG(tr, p, q->next->back, mintrav, maxtrav);
      addTraverseBIG(tr, p, q->next->next->back, mintrav, maxtrav);    
    }
} 





int rearrangeBIG(tree *tr, nodeptr p, int mintrav, int maxtrav)   
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  boolean doP = TRUE, doQ = TRUE;
  
  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
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
  
  if (!isTip(p->number, tr->rdta->numsp) && doP) 
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
	      addTraverseBIG(tr, p, p1->next->back,
			     mintrav, maxtrav);         
	      addTraverseBIG(tr, p, p1->next->next->back,
			     mintrav, maxtrav);          
	    }
	  
	  if (!isTip(p2->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, p, p2->next->back,
			     mintrav, maxtrav);
	      addTraverseBIG(tr, p, p2->next->next->back,
			     mintrav, maxtrav);          
	    }
	  	  
	  hookup(p->next,       p1, p1z, tr->numBranches); 
	  hookup(p->next->next, p2, p2z, tr->numBranches);	   	    	    
	  newviewGeneric(tr, p);	   	    
	}
    }  
  
  if (!isTip(q->number, tr->rdta->numsp) && maxtrav > 0 && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;
      
      /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
	((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
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
	  
	  if (/*! q1->tip*/ !isTip(q1->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, q, q1->next->back,
			     mintrav2 , maxtrav);
	      addTraverseBIG(tr, q, q1->next->next->back,
			     mintrav2 , maxtrav);         
	    }
	  
	  if (/*! q2->tip*/ ! isTip(q2->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, q, q2->next->back,
			     mintrav2 , maxtrav);
	      addTraverseBIG(tr, q, q2->next->next->back,
			     mintrav2 , maxtrav);          
	    }	   
	  
	  hookup(q->next,       q1, q1z, tr->numBranches); 
	  hookup(q->next->next, q2, q2z, tr->numBranches);
	  
	  newviewGeneric(tr, q); 	   
	}
    } 
  
  return  1;
} 





double treeOptimizeRapid(tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt)
{
  int i, index,
    *perm = (int*)NULL;   

  nodeRectifier(tr);

  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3;  
    
  resetInfoList();
  
  resetBestTree(bt);
 
  tr->startLH = tr->endLH = tr->likelihood;
 
  if(tr->doCutoff)
    {
      if(tr->bigCutoff)
	{	  
	  if(tr->itCount == 0)    
	    tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
	  else    		 
	    tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
	}
      else
	{
	  if(tr->itCount == 0)    
	    tr->lhCutoff = tr->likelihood / -1000.0;    
	  else    		 
	    tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
	}    

      tr->itCount = tr->itCount + 1;
      tr->lhAVG = 0;
      tr->lhDEC = 0;
    }
  
  if(adef->permuteTreeoptimize)
    {
      int n = tr->mxtips + tr->mxtips - 2;   
      perm = (int *)rax_malloc(sizeof(int) * (n + 1));
      makePermutation(perm, 1, n, adef);
    }

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {           
      tr->bestOfNode = unlikely;          

      if(adef->permuteTreeoptimize)
	index = perm[i];
      else
	index = i;
      


      if(rearrangeBIG(tr, tr->nodep[index], mintrav, maxtrav))
	{    
	  if(Thorough)
	    {
	      if(tr->endLH > tr->startLH)                 	
		{			   	     
		  restoreTreeFast(tr);	 	 
		  tr->startLH = tr->endLH = tr->likelihood;	 
		  saveBestTree(bt, tr);
		}
	      else
		{ 		  
		  if(tr->bestOfNode != unlikely)		    	     
		    restoreTopologyOnly(tr, bt);		    
		}	   
	    }
	  else
	    {
	      insertInfoList(tr->nodep[index], tr->bestOfNode);	    
	      if(tr->endLH > tr->startLH)                 	
		{		      
		  restoreTreeFast(tr);	  	      
		  tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
		}	    	  
	    }
	}     
    }     

  if(!Thorough)
    {           
      Thorough = 1;  
      
      for(i = 0; i < iList.valid; i++)
	{ 
	  

	  tr->bestOfNode = unlikely;
	  
	  if(rearrangeBIG(tr, iList.list[i].node, mintrav, maxtrav))
	    {	  
	      if(tr->endLH > tr->startLH)                 	
		{	 	     
		  restoreTreeFast(tr);	 	 
		  tr->startLH = tr->endLH = tr->likelihood;	 
		  saveBestTree(bt, tr);
		}
	      else
		{ 
	      
		  if(tr->bestOfNode != unlikely)
		    {	     
		      restoreTopologyOnly(tr, bt);
		    }	
		}      
	    }
	}       
          
      Thorough = 0;
    }

  if(adef->permuteTreeoptimize)
    rax_free(perm);

  return tr->startLH;     
}




boolean testInsertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{    
  if(Thorough)
    {
      if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;    
      
      evaluateGeneric(tr, p->next->next);               
    }
  else
    {
      if (! insertRestoreBIG(tr, p, q))       return FALSE;
      
      {
	nodeptr x, y;
	x = p->next->next;
	y = p->back;
			
	if(! isTip(x->number, tr->rdta->numsp) && isTip(y->number, tr->rdta->numsp))
	  {
	    while ((! x->x)) 
	      {
		if (! (x->x))
		  newviewGeneric(tr, x);		     
	      }
	  }
	
	if(isTip(x->number, tr->rdta->numsp) && !isTip(y->number, tr->rdta->numsp))
	  {
	    while ((! y->x)) 
	      {		  
		if (! (y->x))
		  newviewGeneric(tr, y);
	      }
	  }
	
	if(!isTip(x->number, tr->rdta->numsp) && !isTip(y->number, tr->rdta->numsp))
	  {
	    while ((! x->x) || (! y->x)) 
	      {
		if (! (x->x))
		  newviewGeneric(tr, x);
		if (! (y->x))
		  newviewGeneric(tr, y);
	      }
	  }				      	
	
      }
	
      tr->likelihood = tr->endLH;
    }
     
  return TRUE;
} 

void restoreTreeFast(tree *tr)
{
  removeNodeRestoreBIG(tr, tr->removeNode);    
  testInsertRestoreBIG(tr, tr->removeNode, tr->insertNode);
}

//#define _DEBUG_RELL

static void updateRellTrees(tree *tr, int numberOfSamples)
{
  int i, j;

#ifdef _DEBUG_RELL
  int 
    impr = 0;
#endif

  evaluateGenericVector(tr, tr->start); 

  for(i = 0; i < numberOfSamples; i++)
    {
      double 
	like = 0.0;

      int 
	w = 0,
	*weights = &(tr->resample[i * tr->originalCrunchedLength]);

      for(j = 0; j < tr->originalCrunchedLength; j++)
	{
	  w += weights[j];
	  like += (double)weights[j] * tr->perSiteLL[j];
	}

      if(like > tr->rellTrees->t[i]->likelihood)
	{
#ifdef _DEBUG_RELL
	  impr = 1;
	  if(tr->rellTrees->t[i]->likelihood != unlikely)
	    printf("[%d] better tree found %f -> %f\n", i, tr->rellTrees->t[i]->likelihood, like);
	  else
	    printf("[%d] better tree found unlikely -> %f\n", i, like);
#endif	  
	  saveTreeList(tr->rellTrees, tr, i, like);	
	}      

    }

#ifdef _DEBUG_RELL
  if(impr)
    printf("\n\n");
#endif
  
}

int determineRearrangementSetting(tree *tr,  analdef *adef, bestlist *bestT, bestlist *bt)
{
  int i, mintrav, maxtrav, bestTrav, impr, index, MaxFast,
    *perm = (int*)NULL;
  double startLH; 
  boolean cutoff;  

  MaxFast = 26;

  startLH = tr->likelihood;

  cutoff = tr->doCutoff;
  tr->doCutoff = FALSE;
 
    
  mintrav = 1;
  maxtrav = 5;

  bestTrav = maxtrav = 5;

  impr = 1;

  resetBestTree(bt);

  if(adef->permuteTreeoptimize)
    {
      int n = tr->mxtips + tr->mxtips - 2;   
      perm = (int *)rax_malloc(sizeof(int) * (n + 1));
      makePermutation(perm, 1, n, adef);
    }
  

  while(impr && maxtrav < MaxFast)
    {	
      recallBestTree(bestT, 1, tr);     
      nodeRectifier(tr);
     
      
      if (maxtrav > tr->ntips - 3)  
	maxtrav = tr->ntips - 3;    
 
      tr->startLH = tr->endLH = tr->likelihood;
          
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{                



	  if(adef->permuteTreeoptimize)
	    index = perm[i];
	  else
	    index = i;	 	 

	  tr->bestOfNode = unlikely;
	  if(rearrangeBIG(tr, tr->nodep[index], mintrav, maxtrav))
	    {	     
	      if(tr->endLH > tr->startLH)                 	
		{		 	 	      
		  restoreTreeFast(tr);	        	  	 	  	      
		  tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  	      
		}	         	       	
	    }
	}
      
      treeEvaluate(tr, 0.25);
      saveBestTree(bt, tr);                                    

      /*printf("DETERMINE_BEST: %d %f\n", maxtrav, tr->likelihood);*/

      if(tr->likelihood > startLH)
	{	 
	  startLH = tr->likelihood; 	  	  	  
	  printLog(tr, adef, FALSE);	  
	  bestTrav = maxtrav;	 
	  impr = 1;
	}
      else
	{
	  impr = 0;
	}
      maxtrav += 5;
      
      if(tr->doCutoff)
	{
	  tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));       
  
	  tr->itCount =  tr->itCount + 1;
	  tr->lhAVG = 0;
	  tr->lhDEC = 0;
	}
    }

  recallBestTree(bt, 1, tr);   
  tr->doCutoff = cutoff;

  if(adef->permuteTreeoptimize)
    rax_free(perm);

  
  return bestTrav;     
}



/*
  #define _TERRACES
  
  define _TERRACES to better explore trees on terraces
*/

#ifdef _TERRACES
/* count the number of BS or ML tree search replicates */

int bCount = 0;
#endif

void computeBIGRAPID (tree *tr, analdef *adef, boolean estimateModel) 
{ 
  unsigned int
    vLength = 0;
  int
    i,
    impr, 
    bestTrav,
    rearrangementsMax = 0, 
    rearrangementsMin = 0,    
    thoroughIterations = 0,
    fastIterations = 0;
   
  double lh, previousLh, difference, epsilon;              
  bestlist *bestT, *bt;  
    
#ifdef _TERRACES
  /* store the 20 best trees found in a dedicated list */

  bestlist
    *terrace;
  
  /* output file names */

  char 
    terraceFileName[1024],
    buf[64];
#endif

  hashtable *h = (hashtable*)NULL;
  unsigned int **bitVectors = (unsigned int**)NULL;
  
 
  if(tr->searchConvergenceCriterion)
    {          
      bitVectors = initBitVector(tr, &vLength);
      h = initHashTable(tr->mxtips * 4);   
    }

  bestT = (bestlist *) rax_malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) rax_malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

#ifdef _TERRACES 
  /* initialize the tree list and the output file name for the current tree search/replicate */


  terrace = (bestlist *) rax_malloc(sizeof(bestlist));      
  terrace->ninit = 0;
  initBestTree(terrace, 20, tr->mxtips); 
  
  sprintf(buf, "%d", bCount);
  
  strcpy(terraceFileName,         workdir);
  strcat(terraceFileName,         "RAxML_terrace.");
  strcat(terraceFileName,         run_id);
  strcat(terraceFileName,         ".BS.");
  strcat(terraceFileName,         buf);
  
  printf("%s\n", terraceFileName);
#endif

  initInfoList(50);
 
  difference = 10.0;
  epsilon = 0.01;    
    
  Thorough = 0;     
  
  if(estimateModel)
    {
      if(adef->useBinaryModelFile)
	{
	  readBinaryModel(tr, adef);
	  evaluateGenericInitrav(tr, tr->start);
	  treeEvaluate(tr, 2);
	}
      else
	{
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, 10.0);
	}
    }
  else
    treeEvaluate(tr, 2);  

 



  printLog(tr, adef, FALSE); 

  saveBestTree(bestT, tr);
  
  if(!adef->initialSet)   
    bestTrav = adef->bestTrav = determineRearrangementSetting(tr, adef, bestT, bt);                   
  else
    bestTrav = adef->bestTrav = adef->initial;

  if(estimateModel)
    {
      if(adef->useBinaryModelFile)	
	treeEvaluate(tr, 2);
      else
	{
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, 5.0);
	}
    }
  else
    treeEvaluate(tr, 1);
  
  saveBestTree(bestT, tr); 
  impr = 1;
  if(tr->doCutoff)
    tr->itCount = 0;

 

  while(impr)
    {              
      recallBestTree(bestT, 1, tr); 

      if(tr->searchConvergenceCriterion)
	{
	  int bCounter = 0;	      
	  
	  if(fastIterations > 1)
	    cleanupHashTable(h, (fastIterations % 2));		
	  
	  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, fastIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
				  &bCounter, 1, FALSE, FALSE);	    
	  
	  assert(bCounter == tr->mxtips - 3);	    	   
	  
	  if(fastIterations > 0)
	    {
	      double rrf = convergenceCriterion(h, tr->mxtips);
	      
	      if(rrf <= 0.01) /* 1% cutoff */
		{
		  printBothOpen("ML fast search converged at fast SPR cycle %d with stopping criterion\n", fastIterations);
		  printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
		  cleanupHashTable(h, 0);
		  cleanupHashTable(h, 1);
		  goto cleanup_fast;
		}
	      else		    
		printBothOpen("ML search convergence criterion fast cycle %d->%d Relative Robinson-Foulds %f\n", fastIterations - 1, fastIterations, rrf);
	    }
	}

	 
      fastIterations++;	


      treeEvaluate(tr, 1.0);
      
      
      saveBestTree(bestT, tr);           
      printLog(tr, adef, FALSE);         
      printResult(tr, adef, FALSE);    
      lh = previousLh = tr->likelihood;
   
     
      treeOptimizeRapid(tr, 1, bestTrav, adef, bt);   
      
      impr = 0;
	  
      for(i = 1; i <= bt->nvalid; i++)
	{	    		  	   
	  recallBestTree(bt, i, tr);
	  
	  treeEvaluate(tr, 0.25);	    	 		      	 

	  if(adef->rellBootstrap)
	    updateRellTrees(tr, NUM_RELL_BOOTSTRAPS);

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

 

  if(tr->searchConvergenceCriterion)
    {
      cleanupHashTable(h, 0);
      cleanupHashTable(h, 1);
    }

 cleanup_fast:

  Thorough = 1;
  impr = 1;
  
  recallBestTree(bestT, 1, tr); 
  if(estimateModel)
    {
      if(adef->useBinaryModelFile)	
	treeEvaluate(tr, 2);
      else
	{
	  evaluateGenericInitrav(tr, tr->start);
	  modOpt(tr, adef, FALSE, 1.0);
	}
    }
  else
    treeEvaluate(tr, 1.0);

  while(1)
    {	
      recallBestTree(bestT, 1, tr);    
      if(impr)
	{	    
	  printResult(tr, adef, FALSE);
	  rearrangementsMin = 1;
	  rearrangementsMax = adef->stepwidth;	

	 

	  if(tr->searchConvergenceCriterion)
	    {
	      int bCounter = 0;	      

	      if(thoroughIterations > 1)
		cleanupHashTable(h, (thoroughIterations % 2));		
		
	      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, thoroughIterations % 2, BIPARTITIONS_RF, (branchInfo *)NULL,
				      &bCounter, 1, FALSE, FALSE);	    
	      
	      assert(bCounter == tr->mxtips - 3);	    	   
	      
	      if(thoroughIterations > 0)
		{
		  double rrf = convergenceCriterion(h, tr->mxtips);
		  
		  if(rrf <= 0.01) /* 1% cutoff */
		    {
		      printBothOpen("ML search converged at thorough SPR cycle %d with stopping criterion\n", thoroughIterations);
		      printBothOpen("Relative Robinson-Foulds (RF) distance between respective best trees after one succseful SPR cycle: %f%s\n", rrf, "%");
		      goto cleanup;
		    }
		  else		    
		    printBothOpen("ML search convergence criterion thorough cycle %d->%d Relative Robinson-Foulds %f\n", thoroughIterations - 1, thoroughIterations, rrf);
		}
	    }

	 
	   	  
	  thoroughIterations++;	  
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
      
      printLog(tr, adef, FALSE);
      treeOptimizeRapid(tr, rearrangementsMin, rearrangementsMax, adef, bt);
	
      impr = 0;			      		            

      for(i = 1; i <= bt->nvalid; i++)
	{		 
	  recallBestTree(bt, i, tr);	 	    	    	
	  
	  treeEvaluate(tr, 0.25);
	    	 
	  if(adef->rellBootstrap)
	    updateRellTrees(tr, NUM_RELL_BOOTSTRAPS);
	  
#ifdef _TERRACES
	  /* save all 20 best trees in the terrace tree list */
	  saveBestTree(terrace, tr);
#endif

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

#ifdef _TERRACES
  {
    double
      bestLH = tr->likelihood;
    FILE 
      *f = myfopen(terraceFileName, "w");
    
    /* print out likelihood of best tree found */

    printf("best tree: %f\n", tr->likelihood);

    /* print out likelihoods of 20 best trees found during the tree search */

    for(i = 1; i <= terrace->nvalid; i++)
      {
	recallBestTree(terrace, i, tr);
	
	/* if the likelihood scores are smaller than some epsilon 0.000001
	   print the tree to file */
	   
	if(ABS(bestLH - tr->likelihood) < 0.000001)
	  {
	    printf("%d %f\n", i, tr->likelihood);
	    Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, adef, NO_BRANCHES, FALSE, FALSE, FALSE, FALSE);
	    
	    fprintf(f, "%s\n", tr->tree_string); 
	  }
      }

    fclose(f);
    /* increment tree search counter */
    bCount++;
  }
#endif
 
  
  if(tr->searchConvergenceCriterion)
    {
      freeBitVectors(bitVectors, 2 * tr->mxtips);
      rax_free(bitVectors);
      freeHashTable(h);
      rax_free(h);
    }
  
  freeBestTree(bestT);
  rax_free(bestT);
  freeBestTree(bt);
  rax_free(bt);

#ifdef _TERRACES
  /* free terrace tree list */

  freeBestTree(terrace);
  rax_free(terrace);
#endif

  freeInfoList();  
  printLog(tr, adef, FALSE);
  printResult(tr, adef, FALSE);
}




boolean treeEvaluate (tree *tr, double smoothFactor)       /* Evaluate a user tree */
{
  boolean result;
 
  if(tr->useBrLenScaler)
    assert(0);
 
  result = smoothTree(tr, (int)((double)smoothings * smoothFactor));
  
  assert(result); 

  evaluateGeneric(tr, tr->start);         

  return TRUE;
}

static void setupBranches(tree *tr, nodeptr p,  branchInfo *bInf)
{
  int    
    countBranches = tr->branchCounter;

  if(isTip(p->number, tr->mxtips))    
    {      
      p->bInf       = &bInf[countBranches];
      p->back->bInf = &bInf[countBranches];               	      

      bInf[countBranches].oP = p;
      bInf[countBranches].oQ = p->back;
           
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

      tr->branchCounter =  tr->branchCounter + 1;      

      q = p->next;

      while(q != p)
	{
	  setupBranches(tr, q->back, bInf);	
	  q = q->next;
	}
     
      return;
    }
}
 


static void smoothTreeRandom (tree *tr, int maxtimes, analdef *adef)
{
  int i, k, *perm;

  tr->branchCounter = 0;
  tr->numberOfBranches = 2 * tr->mxtips - 3;
  
  tr->bInf = (branchInfo*)rax_malloc(tr->numberOfBranches * sizeof(branchInfo)); 

  perm = (int*)rax_malloc(sizeof(int) * tr->numberOfBranches);

  setupBranches(tr, tr->start->back, tr->bInf);



  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  /*printf("%d \n", maxtimes);*/

  while (--maxtimes >= 0) 
    {    
      for(i = 0; i < tr->numBranches; i++)	
	tr->partitionSmoothed[i] = TRUE;		

      /*printf("%d %d\n", maxtimes, tr->numberOfBranches);*/

      makePermutation(perm, 0, tr->numberOfBranches - 1, adef);

      for(k = 0; k < tr->numberOfBranches; k++)
	{
	  /*printf("%d Node %d\n", k, tr->bInf[k].oP->number);*/
	  update(tr, tr->bInf[perm[k]].oP);
	  newviewGeneric(tr, tr->bInf[perm[k]].oP);     
	}

      /*if (! smooth(tr, p->back))       return FALSE;
      if (!isTip(p->number, tr->rdta->numsp)) 
	{
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smooth(tr, q->back))   return FALSE;
	      q = q->next;
	    }
	}             
      */

      if (allSmoothed(tr)) 
	break;      
    }

  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;

  rax_free(tr->bInf);
  rax_free(perm);
} 


void treeEvaluateRandom (tree *tr, double smoothFactor, analdef *adef)       
{
 
  smoothTreeRandom(tr, (int)((double)smoothings * smoothFactor), adef);
  

  evaluateGeneric(tr, tr->start);       
}

void treeEvaluateProgressive(tree *tr)
{  
  int 
    i, k;

  tr->branchCounter = 0;
  tr->numberOfBranches = 2 * tr->mxtips - 3;
  
  tr->bInf = (branchInfo*)rax_malloc(tr->numberOfBranches * sizeof(branchInfo));  

  setupBranches(tr, tr->start->back, tr->bInf);

  assert(tr->branchCounter == tr->numberOfBranches);
  
  for(i = 0; i < tr->numBranches; i++)
    tr->partitionConverged[i] = FALSE;
  
  for(i = 0; i < 10; i++)
    {
      for(k = 0; k < tr->numberOfBranches; k++)
	{      
	  update(tr, tr->bInf[k].oP);
	  newviewGeneric(tr, tr->bInf[k].oP);     
	}
      evaluateGenericInitrav(tr, tr->start);
      printf("It %d %f \n", i, tr->likelihood);
    }
}
   

  



