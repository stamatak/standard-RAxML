// this file contains legacy code of functions that are not being used any more but that could
// become useful in the future again!


//old mygetopt, basically an adaption of the GNU implementation 

static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static 
    int sp = 1;

  register 
    int c;
  
  register 
    char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else
    {
      if(strcmp(argv[*optind], "--") == 0)
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0)
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0')
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':')
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc)
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    }
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    }
  else
    {
      if(argv[*optind][++sp] == '\0')
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
 
  return(c);
}


//old MRE implementation 


#else




static void mre(hashtable *h, boolean icp, entry*** sbi, int* len, int which, int n, unsigned int vectorLength, boolean sortp, tree *tr, boolean bootStopping)
{
  entry **sbw;
  unsigned int 
    i = 0,
    j = 0,
    k = 0;  
 
  sbw = (entry **) rax_calloc(h->entryCount, sizeof(entry *));  

  for(i = 0; i < h->tableSize; i++)
    {		
      if(h->table[i] != NULL)
	{
	  entry *e = h->table[i];
	  do
	    {
	      sbw[j] = e;
	      j++;
	      e = e->next;
	    }
	  while(e != NULL);
	}
    }

  assert(j == h->entryCount);

  if(which == 0)    
    qsort(sbw, h->entryCount, sizeof(entry *), _sortByWeight0);      
  else    
    qsort(sbw, h->entryCount, sizeof(entry *), _sortByWeight1);      

  /* ***********************************          */
  /* SOS SBI is never rax_freed ********************* */
  /* ******************************************** */
  /**** this will cause problems for repeated invocations */
  /**** with the bootstopping MRE VERSION !!!!!!        ***/
      


  *sbi = (entry **)rax_calloc(n - 3, sizeof(entry *));

  *len = 0;

  if(icp == FALSE)
    {     
      for(i = 0; i < h->entryCount && (*len) < n-3; i++)
	{	
	  boolean compatflag = TRUE;

	  assert(*len < n-3);		  
	
	  /*  for(k = 0; k < (unsigned int)(*len); k++)	  */
	  /*if(sbw[i]->supportFromTreeset[which] <= mr_thresh) */
	    for(k = ((unsigned int)(*len)); k > 0; k--)
	      {
		/*
		  k indexes sbi
		  j indexes sbw
		  need to compare the two
		*/
		
		if(!compatible((*sbi)[k-1], sbw[i], vectorLength))		
		  {
		    compatflag = FALSE;	    
		    break;
		  }
	      }
	  
	  if(compatflag)
	    {	      
	      (*sbi)[*len] = sbw[i];
	      (*len)++;
	    }         
	}
    }
  else 
    {      
      for(i = 0; i < (unsigned int)(n-3); i++)
	{
	  (*sbi)[i] = sbw[i];
	  (*len)++;
	}
    }

  rax_free(sbw);

  if (sortp == TRUE)
    qsort(*sbi, (*len), sizeof(entry *), sortByIndex);    

  return;
}







static void printBip(entry *curBip, entry **consensusBips, const unsigned int consensusBipLen, const int numtips, const unsigned int vectorLen, 
		     boolean *processed, tree *tr, FILE *outf, const int numberOfTrees, boolean topLevel, unsigned int *printCounter)
{
  int
    branchLabel,     
    printed = 0;

  unsigned int 
    i,
    j;

  unsigned int *subBip = (unsigned int *)rax_calloc(vectorLen, sizeof(unsigned int));

  double 
    support = 0.0;

  for(i = 0; i < consensusBipLen; i++)
    {
      if((!processed[i]) && issubset(consensusBips[i]->bitVector, curBip->bitVector, vectorLen))
	{
	  boolean processThisRound = TRUE;
	  
	  for (j = 0; j < consensusBipLen; j++)	    
	    if(j != i && !processed[j] && issubset(consensusBips[i]->bitVector, consensusBips[j]->bitVector, vectorLen))		
	      processThisRound = FALSE;			   
	  
	  if(processThisRound == TRUE)
	    {
	      processed[i] = TRUE;

	      for(j = 0; j < vectorLen; j++)
		subBip[j] |= consensusBips[i]->bitVector[j];
	      
	      if(printed == 0 && !topLevel)		
		fprintf(outf, "(");		
	      else		
		fprintf(outf, ",");
		
	      printBip(consensusBips[i], consensusBips, consensusBipLen, numtips, vectorLen, processed, tr, outf, numberOfTrees, FALSE, printCounter);

	      printed += 1;
	    }
	}
    }	
  
  for(i = 0; i < ((unsigned int)numtips); i++)
    {
      if((((curBip->bitVector[i/MASK_LENGTH] & mask32[i%MASK_LENGTH]) > 0) && ((subBip[i/MASK_LENGTH] & mask32[i%MASK_LENGTH]) == 0) ) == TRUE)
	{
	  if(printed == 0 && !topLevel)	    
	    fprintf(outf,"(");	   
	  else	    
	    fprintf(outf,",");
	   
	  fprintf(outf,"%s", tr->nameList[i+1]);
	  printed += 1;
	}
    }

  rax_free(subBip);

  support = ((double)(curBip->supportFromTreeset[0])) / ((double) (numberOfTrees));
  branchLabel = (int)(0.5 + support * 100.0);
  
  if(!topLevel)
    {
      *printCounter = *printCounter + 1;
      fprintf(outf,"):1.0[%d]", branchLabel);
    }
}

void computeConsensusOnly(tree *tr, char *treeSetFileName, analdef *adef)
{        
  hashtable 
    *h = initHashTable(tr->mxtips * FC_INIT * 10); 

  hashNumberType
    entries = 0;
  
  int  
    numberOfTrees = 0, 
    i, 
    j, 
    l,
    treeVectorLength,     
    consensusBipsLen,
    mr_thresh;

  unsigned int
    printCounter = 0,
    vectorLength,
    **bitVectors = initBitVector(tr, &vectorLength),
    *topBip;

  entry 
    topBipE,
    **consensusBips;
 
  boolean  
    *processed;

  char 
    consensusFileName[1024];
  
  FILE 
    *outf,
    *treeFile = getNumberOfTrees(tr, treeSetFileName, adef);
     
  numberOfTrees = tr->numberOfTrees; 

  checkTreeNumber(numberOfTrees, treeSetFileName);

  mr_thresh = ((double)numberOfTrees / 2.0);  
  
  assert(sizeof(unsigned char) == 1);
  
  if(numberOfTrees % MASK_LENGTH == 0)
    treeVectorLength = numberOfTrees / MASK_LENGTH;
  else
    treeVectorLength = 1 + (numberOfTrees / MASK_LENGTH);  

  for(i = 1; i <= numberOfTrees; i++)
    {                  
      int 
	bCount = 0;
      
      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE, FALSE);               
      
      assert(tr->mxtips == tr->ntips);
      
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			      &bCount, treeVectorLength, FALSE, FALSE);
      
      assert(bCount == tr->mxtips - 3);                     
    } 

  if(tr->consensusType == MR_CONSENSUS || tr->consensusType == STRICT_CONSENSUS)
    {
      consensusBips = (entry **)rax_calloc(tr->mxtips - 3, sizeof(entry *));
      consensusBipsLen = 0;
    }

  for(j = 0; j < (int)h->tableSize; j++)
    {		
      if(h->table[j] != NULL)
	{
	  entry *e = h->table[j];
	  
	  do
	    {	
	      int cnt = 0;
	      
	      unsigned int 
		*set = e->treeVector;
	      
	      for(l = 0; l < numberOfTrees; l++)					     
		if((set[l / MASK_LENGTH] != 0) && (set[l / MASK_LENGTH] & mask32[l % MASK_LENGTH]))		    			    
		  cnt++;		    		     		
	      
	      if(tr->consensusType == MR_CONSENSUS)
		{
		  if(cnt > mr_thresh)			      
		    {
		      consensusBips[consensusBipsLen] = e;
		      consensusBipsLen++;
		    }
		}

	      if(tr->consensusType == STRICT_CONSENSUS)
		{
		  if(cnt == numberOfTrees)
		    {
		      consensusBips[consensusBipsLen] = e;
		      consensusBipsLen++;
		    }
		}

	      e->supportFromTreeset[0] = cnt;
	      e = e->next;
	      entries++;
	    }
	  while(e != NULL);		
	}	  	        
    }	

  fclose(treeFile); 
  assert(entries == h->entryCount);
  
  if(tr->consensusType == MR_CONSENSUS || tr->consensusType == STRICT_CONSENSUS)
    assert(consensusBipsLen <= (tr->mxtips - 3));

  if(tr->consensusType == MRE_CONSENSUS)    
    mre(h, FALSE, &consensusBips, &consensusBipsLen, 0, tr->mxtips, vectorLength, FALSE, tr);    

  
  /* printf("Bips OLD %d\n", consensusBipsLen); */

  processed = (boolean *) rax_calloc(consensusBipsLen, sizeof(boolean));

  topBip = (unsigned int *) rax_calloc(vectorLength, sizeof(unsigned int));
  
  for(i = 1; i < tr->mxtips; i++)
    topBip[i / MASK_LENGTH] |= mask32[i % MASK_LENGTH];  

  topBipE.bitVector = topBip;
  topBipE.supportFromTreeset[0] = numberOfTrees;

  strcpy(consensusFileName,         workdir);  
  
  switch(tr->consensusType)
    {
    case MR_CONSENSUS:
      strcat(consensusFileName,         "RAxML_MajorityRuleConsensusTree.");
      break;
    case MRE_CONSENSUS:
      strcat(consensusFileName,         "RAxML_MajorityRuleExtendedConsensusTree.");
      break;
    case STRICT_CONSENSUS:
      strcat(consensusFileName,         "RAxML_StrictConsensusTree.");
      break;
    default:
      assert(0);
    }
  
  strcat(consensusFileName,         run_id);

  outf = myfopen(consensusFileName, "wb");

  fprintf(outf, "(%s", tr->nameList[1]);
  printBip(&topBipE, consensusBips, consensusBipsLen, tr->mxtips, vectorLength, processed, tr, outf, numberOfTrees, TRUE, &printCounter);  
  fprintf(outf, ");\n");

  assert(consensusBipsLen == (int)printCounter);

  fclose(outf);
  
  switch(tr->consensusType)
    {
    case MR_CONSENSUS:
      printBothOpen("RAxML Majority Rule consensus tree written to file: %s\n", consensusFileName);
      break;
    case MRE_CONSENSUS:
      printBothOpen("RAxML extended Majority Rule consensus tree written to file: %s\n", consensusFileName);
      break;
    case STRICT_CONSENSUS:
      printBothOpen("RAxML strict consensus tree written to file: %s\n", consensusFileName);
      break;
    default:
      assert(0);
    }
  
  
  rax_free(topBip);
  rax_free(processed);

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  freeHashTable(h);
  rax_free(h);
  rax_free(consensusBips);  
  
  exit(0);   
}


//old slow initial code for tree plausibility checking 

/* function to extract the bit mask for the taxa that are present in the small tree */

static void setupMask(unsigned int *smallTreeMask, nodeptr p, int numsp)
{
  if(isTip(p->number, numsp))
    smallTreeMask[(p->number - 1) / MASK_LENGTH] |= mask32[(p->number - 1) % MASK_LENGTH];
  else
    {    
      nodeptr 
	q = p->next;

      /* I had to change this function to account for mult-furcating trees.
	 In this case an inner node can have more than 3 cyclically linked 
	 elements, because there might be more than 3 outgoing branches 
	 from an inner node */

      while(q != p)
	{
	  setupMask(smallTreeMask, q->back, numsp);
	  q = q->next;
	}
      
      //old code below 
      //setupMask(smallTreeMask, p->next->back, numsp);	  
      //setupMask(smallTreeMask, p->next->next->back, numsp);      
    }
}


static void newviewBipartitions(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength)
{
  if(isTip(p->number, numsp))
    return;
  {
    nodeptr 
      q = p->next->back, 
      r = p->next->next->back;
    unsigned int       
      *vector = bitVectors[p->number],
      *left  = bitVectors[q->number],
      *right = bitVectors[r->number];
    unsigned 
      int i;           

    while(!p->x)
      {	
	if(!p->x)
	  getxnode(p);
      }

    p->hash = q->hash ^ r->hash;

    if(isTip(q->number, numsp) && isTip(r->number, numsp))
      {		
	for(i = 0; i < vectorLength; i++)
	  vector[i] = left[i] | right[i];	  	
      }
    else
      {	
	if(isTip(q->number, numsp) || isTip(r->number, numsp))
	  {
	    if(isTip(r->number, numsp))
	      {	
		nodeptr tmp = r;
		r = q;
		q = tmp;
	      }	   
	    	    
	    while(!r->x)
	      {
		if(!r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	    	 
	  }
	else
	  {	    
	    while((!r->x) || (!q->x))
	      {
		if(!q->x)
		  newviewBipartitions(bitVectors, q, numsp, vectorLength);
		if(!r->x)
		  newviewBipartitions(bitVectors, r, numsp, vectorLength);
	      }	   	    	    	    	   

	    for(i = 0; i < vectorLength; i++)
	      vector[i] = left[i] | right[i];	 
	  }

      }     
  }     
}


/* this function actually traverses the small tree, generates the bit vectors for all 
   non-trivial bipartitions and simultaneously counts how many bipartitions (already stored in the has table) are shared with the big tree
*/

static int bitVectorTraversePlausibility(unsigned int **bitVectors, nodeptr p, int numsp, unsigned int vectorLength, hashtable *h,
					 int *countBranches, int firstTaxon, tree *tr, boolean multifurcating)
{

  /* trivial bipartition */

  if(isTip(p->number, numsp))
    return 0;
  else
    {
      int 
	found = 0;

      nodeptr 
	q = p->next;          

      /* recursively descend into the tree and get the bips of all subtrees first */

      do 
	{
	  found = found + bitVectorTraversePlausibility(bitVectors, q->back, numsp, vectorLength, h, countBranches, firstTaxon, tr, multifurcating);
	  q = q->next;
	}
      while(q != p);
           
      /* compute the bipartition induced by the current branch p, p->back,
	 here we invoke two different functions, depending on whether we are dealing with 
	 a multi-furcating or bifurcating tree.
      */

      if(multifurcating)
	newviewBipartitionsMultifurcating(bitVectors, p, numsp, vectorLength);
      else
	newviewBipartitions(bitVectors, p, numsp, vectorLength);
      
      assert(p->x);      

      /* if p->back does not lead to a tip this is an inner branch that induces a non-trivial bipartition.
	 in this case we need to lookup if the induced bipartition is already contained in the hash table 
      */

      if(!(isTip(p->back->number, numsp)))
	{	
	  /* this is the bit vector to insert into the hash table */
	  unsigned int 
	    *toInsert = bitVectors[p->number];
	  
	  /* compute the hash number on that bit vector */
	  hashNumberType 
	    position = oat_hash((unsigned char *)toInsert, sizeof(unsigned int) * vectorLength) % h->tableSize;	 	 

	  /* each bipartition can be stored in two forms (the two bit-wise complements
	     we always canonically store that version of the bit-vector that does not contain the 
	     first taxon of the small tree, we use an assertion to make sure that all is correct */

	  assert(!(toInsert[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH]));	 
	  	      
	  /* increment the branch counter to assure that all inner branches are traversed */
	  
	  *countBranches =  *countBranches + 1;	
	  	 	
	  /* now look up this bipartition in the hash table, If it is present the number of 
	     shared bipartitions between the small and the big tree is incremented by 1 */
	   
	  found = found + findHash(toInsert, h, vectorLength, position);	
	}
      return found;
    }
}



////multifurcating trees 

void plausibilityChecker(tree *tr, analdef *adef)
{
  FILE 
    *treeFile,
    *rfFile;
  
  tree 
    *smallTree = (tree *)rax_malloc(sizeof(tree));

  char 
    rfFileName[1024];
 
  /* init hash table for big reference tree */
  
  hashtable
    *h      = initHashTable(tr->mxtips * 2 * 2);
  
  /* init the bit vectors we need for computing and storing bipartitions during 
     the tree traversal */
  unsigned int 
    vLength, 
    **bitVectors = initBitVector(tr, &vLength);
   
  int
    numberOfTreesAnalyzed = 0,
    branchCounter = 0,
    i;

  double 
    avgRF = 0.0;

  /* set up an output file name */

  strcpy(rfFileName,         workdir);  
  strcat(rfFileName,         "RAxML_RF-Distances.");
  strcat(rfFileName,         run_id);

  rfFile = myfopen(rfFileName, "wb");  

  assert(adef->mode ==  PLAUSIBILITY_CHECKER);

  /* open the big reference tree file and parse it */

  treeFile = myfopen(tree_file, "r");

  printBothOpen("Parsing reference tree %s\n", tree_file);

  treeReadLen(treeFile, tr, FALSE, TRUE, TRUE, adef, TRUE, FALSE);

  assert(tr->mxtips == tr->ntips);

  printBothOpen("The reference tree has %d tips\n", tr->ntips);

  fclose(treeFile);
  
  /* extract all induced bipartitions from the big tree and store them in the hastable */
  
  bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vLength, h, 0, BIPARTITIONS_RF, (branchInfo *)NULL,
			  &branchCounter, 1, FALSE, FALSE);
     
  assert(branchCounter == tr->mxtips - 3);   
  
  /* now see how many small trees we have */

  treeFile = getNumberOfTrees(tr, bootStrapFile, adef);

  checkTreeNumber(tr->numberOfTrees, bootStrapFile);

  /* allocate a data structure for parsing the potentially mult-furcating tree */

  allocateMultifurcations(tr, smallTree);

  /* loop over all small trees */

  for(i = 0; i < tr->numberOfTrees;  i++)
    {          
      int           
	numberOfSplits = readMultifurcatingTree(treeFile, smallTree, adef, TRUE);

      if(numberOfSplits > 0)
	{
	  unsigned int
	    entryCount = 0,
	    k,
	    j,	
	    *masked    = (unsigned int *)rax_calloc(vLength, sizeof(unsigned int)),
	    *smallTreeMask = (unsigned int *)rax_calloc(vLength, sizeof(unsigned int));

	  hashtable
	    *rehash = initHashTable(tr->mxtips * 2 * 2);

	  double
	    rf,
	    maxRF;

	  int 
	    bCounter = 0,  
	    bips,
	    firstTaxon,
	    taxa = 0;

	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Small tree %d has %d tips and %d bipartitions\n", i, smallTree->ntips, numberOfSplits);    

	  /* compute the maximum RF distance for computing the relative RF distance later-on */
	  
	  /* note that here we need to pay attention, since the RF distance is not normalized 
	     by 2 * (n-3) but we need to account for the fact that the multifurcating small tree 
	     will potentially contain less bipartitions. 
	     Hence the normalization factor is obtained as 2 * numberOfSplits, where numberOfSplits is the number of bipartitions
	     in the small tree.
	  */
	  
	  maxRF = (double)(2 * numberOfSplits);
	  
	  /* now set up a bit mask where only the bits are set to one for those 
	     taxa that are actually present in the small tree we just read */
	  
	  /* note that I had to apply some small changes to this function to make it work for 
	     multi-furcating trees ! */

	  setupMask(smallTreeMask, smallTree->start,       smallTree->mxtips);
	  setupMask(smallTreeMask, smallTree->start->back, smallTree->mxtips);

	  /* now get the index of the first taxon of the small tree.
	     we will use this to unambiguously store the bipartitions 
	  */
	  
	  firstTaxon = smallTree->start->number;
	  
	  /* make sure that this bit vector is set up correctly, i.e., that 
	     it contains as many non-zero bits as there are taxa in this small tree 
	  */
	  
	  for(j = 0; j < vLength; j++)
	    taxa += BIT_COUNT(smallTreeMask[j]);
	  assert(taxa == smallTree->ntips);
	  
	  /* now re-hash the big tree by applying the above bit mask */
	  
	  
	  /* loop over hash table */
	  
	  for(k = 0, entryCount = 0; k < h->tableSize; k++)	     
	    {    
	      if(h->table[k] != NULL)
		{
		  entry *e = h->table[k];
		  
		  /* we resolve collisions by chaining, hence the loop here */
		  
		  do
		    {
		      unsigned int 
			*bitVector = e->bitVector; 
		      
		      hashNumberType 
			position;
		      
		      int 
			count = 0;
		      
		      /* double check that our tree mask contains the first taxon of the small tree */
		      
		      assert(smallTreeMask[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH]);
		      
		      /* if the first taxon is set then we will re-hash the bit-wise complement of the 
			 bit vector.
			 The count variable is used for a small optimization */
		      
		      if(bitVector[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH])		    
			{
			  //hash complement
			  
			  for(j = 0; j < vLength; j++)
			    {
			      masked[j] = (~bitVector[j]) & smallTreeMask[j];			     
			      count += BIT_COUNT(masked[j]);
			    }
			}
		      else
			{
			  //hash this vector 
			  
			  for(j = 0; j < vLength; j++)
			    {
			      masked[j] = bitVector[j] & smallTreeMask[j];  
			      count += BIT_COUNT(masked[j]);      
			    }
			}
		      
		      /* note that padding the last bits is not required because they are set to 0 automatically by smallTreeMask */	
		      
		      /* make sure that we will re-hash  the canonic representation of the bipartition 
			 where the bit for firstTaxon is set to 0!
		      */
		      
		      assert(!(masked[(firstTaxon - 1) / MASK_LENGTH] & mask32[(firstTaxon - 1) % MASK_LENGTH]));
		      
		      /* only if the masked bipartition of the large tree is a non-trivial bipartition (two or more bits set to 1 
			 will we re-hash it */
		      
		      if(count > 1)
			{
			  /* compute hash */
			  position = oat_hash((unsigned char *)masked, sizeof(unsigned int) * vLength);
			  position = position % rehash->tableSize;
			  
			  /* re-hash to the new hash table that contains the bips of the large tree, pruned down 
			     to the taxa contained in the small tree
			  */
			  insertHashPlausibility(masked, rehash, vLength, position);
			}		
		      
		      entryCount++;
		      
		      e = e->next;
		    }
		  while(e != NULL);
		}
	    }
	  
	  /* make sure that we tried to re-hash all bipartitions of the original tree */
	  
	  assert(entryCount == (unsigned int)(tr->mxtips - 3));
	  
	  /* now traverse the small tree and count how many bipartitions it shares 
	     with the corresponding induced tree from the large tree */
	  
	  /* the following function also had to be modified to account for multi-furcating trees ! */
	  
	  bips = bitVectorTraversePlausibility(bitVectors, smallTree->start->back, smallTree->mxtips, vLength, rehash, &bCounter, firstTaxon, smallTree, TRUE);
	  
	  /* compute the relative RF */
	  
	  rf = (double)(2 * (numberOfSplits - bips)) / maxRF;           
	  
	  assert(numberOfSplits >= bips);

	  assert(rf <= 1.0);
	  
	  avgRF += rf;
	  
	  if(numberOfTreesAnalyzed % 100 == 0)
	    printBothOpen("Relative RF tree %d: %f\n\n", i, rf);

	  fprintf(rfFile, "%d %f\n", i, rf);
	  
	  /* I also modified this assertion, we nee to make sure here that we checked all non-trivial splits/bipartitions 
	     in the multi-furcating tree whech can be less than n - 3 ! */
	  
	  assert(bCounter == numberOfSplits);         
	  
	  /* free masks and hast table for this iteration */
	  
	  rax_free(smallTreeMask);
	  rax_free(masked);
	  freeHashTable(rehash);
	  rax_free(rehash);
	  numberOfTreesAnalyzed++;
	}
    }

  printBothOpen("Number of small trees skipped: %d\n\n", tr->numberOfTrees - numberOfTreesAnalyzed);
  
  printBothOpen("Average RF distance %f\n\n", avgRF / (double)numberOfTreesAnalyzed);

  printBothOpen("Total execution time: %f secs\n\n", gettime() - masterTime);

  printBothOpen("\nFile containing all %d pair-wise RF distances written to file %s\n\n", numberOfTreesAnalyzed, rfFileName);

  

  fclose(treeFile);
  fclose(rfFile);    
  
  /* free the data structure used for parsing the potentially multi-furcating tree */

  freeMultifurcations(smallTree);
  rax_free(smallTree);

  freeBitVectors(bitVectors, 2 * tr->mxtips);
  rax_free(bitVectors);
  
  freeHashTable(h);
  rax_free(h);
}

