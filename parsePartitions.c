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
#include <strings.h>




#include "axml.h"

/*****************************FUNCTIONS FOR READING MULTIPLE MODEL SPECIFICATIONS************************************************/


extern char modelFileName[1024];
extern char excludeFileName[1024];
extern char secondaryStructureFileName[1024];


extern char seq_file[1024];

extern char *protModels[NUM_PROT_MODELS];

static boolean lineContainsOnlyWhiteChars(char *line)
{
  int i, n = strlen(line);

  if(n == 0)
    return TRUE;

  for(i = 0; i < n; i++)
    {
      if(!whitechar(line[i]))
	return FALSE;
    }
  return TRUE;
}


static int isNum(char c)
{
  
  return (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
	  c == '5' || c == '6' || c == '7' || c == '8' || c == '9');
}


static void skipWhites(char **ch)
{
  while(**ch == ' ' || **ch == '\t')
    *ch = *ch + 1;
}

static void analyzeIdentifier(char **ch, int modelNumber, tree *tr)
{
  char 
    ident[2048] = "",
    model[2048] = "",
    thisModel[2048] = "";
  
  int 
    i = 0, 
    n, 
    j,
    containsComma = 0;

  while(**ch != '=')
    {
      if(**ch != ' ' && **ch != '\t')
	{
	  ident[i] = **ch;      
	  i++;
	}
      *ch = *ch + 1;
    }
  
  n = i;
  i = 0;
  
  for(i = 0; i < n; i++)
    if(ident[i] == ',') 
      containsComma = 1;

  if(!containsComma)
    {
      printf("Error, model file must have format: DNA or AA model, then a comma, and then the partition name\n");
      exit(-1);
    }
  else
    {
      boolean 
	useProteinSubstitutionFile = FALSE,
	found = FALSE;
      
      int 
	openBracket = 0,
	closeBracket = 0,
	openPos = 0,
	closePos = 0;
            
      i = 0;
      
      while(ident[i] != ',')
	{	 
	  if(ident[i] == '[')
	    {
	      openPos = i;
	      openBracket++;
	    }
	  if(ident[i] == ']')
	    {
	      closePos = i;
	      closeBracket++;
	    }
	  model[i] = ident[i];
	  i++;
	}     

      if(closeBracket > 0 || openBracket > 0)
	{
	  if((closeBracket == 1) && (openBracket == 1) && (openPos < closePos))
	    useProteinSubstitutionFile = TRUE;
	  else
	    {
	      printf("\nError: Apparently you want to specify a user-defined protein substitution model that shall be read from file\n");
	      printf("It must be enclosed in opening and closing bracktes like this: [fileName]\n\n");
	      printf("you specified: %s\n\n", model);
	      exit(-1);
	    }
	}
      
      if(useProteinSubstitutionFile)
	{
	  char 
	    protFileName[2048] = "";

	  int 
	    pos,
	    k,
	    lower = 0,
	    upper = i - 1;

	  while(model[lower] == '[' || model[lower] == ' ')
	    lower++;

	  while(model[upper] == ']' || model[upper] == ' ')
	    upper--;
	  
	  assert(lower < upper);

	  for(k = lower, pos = 0; k <= upper; k++, pos++)
	    protFileName[pos] = model[k];
	  
	  protFileName[pos] = '\0';

	  

	  if(!filexists(protFileName))
	    {
	      printf("\n\ncustom protein substitution file [%s] you want to use does not exist!\n", protFileName);
	      printf("you need to specify the full path\n");
	      printf("the file name shall not contain blanks!\n\n");
	      exit(0);
	    }
	  

	  strcpy(tr->initialPartitionData[modelNumber].proteinSubstitutionFileName, protFileName);
	  /*printf("%s \n", tr->initialPartitionData[modelNumber].proteinSubstitutionFileName);*/
	  
	  tr->initialPartitionData[modelNumber].protModels = PROT_FILE;		  
	  tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = TRUE;
	  tr->initialPartitionData[modelNumber].dataType   = AA_DATA;
	}
      else
	{
	  /* AA */
	  
	  for(i = 0; i < NUM_PROT_MODELS && !found; i++)
	    {	
	      strcpy(thisModel, protModels[i]);
	      
	      if(strcasecmp(model, thisModel) == 0)
		{	      	      
		  tr->initialPartitionData[modelNumber].protModels = i;		  
		  tr->initialPartitionData[modelNumber].usePredefinedProtFreqs = TRUE;
		  tr->initialPartitionData[modelNumber].dataType   = AA_DATA;		  
		  found = TRUE;
		}
	      

	      if(i != GTR && i != GTR_UNLINKED)
		{
		  strcpy(thisModel, protModels[i]);
		  strcat(thisModel, "F");
		  
		  if(strcasecmp(model, thisModel) == 0)
		    {	      
		      tr->initialPartitionData[modelNumber].protModels = i;		  
		      tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = FALSE;
		      tr->initialPartitionData[modelNumber].dataType   = AA_DATA;		     
		      found = TRUE;
		    }
		}
	      
	      if(found && (tr->initialPartitionData[modelNumber].protModels == GTR || tr->initialPartitionData[modelNumber].protModels == GTR_UNLINKED))
		tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = FALSE;		    		
	    }
	  
	  if(!found)
	    {		  	  
	      if(strcasecmp(model, "DNA") == 0)
		{	     	      
		  tr->initialPartitionData[modelNumber].protModels = -1;		  
		  tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = FALSE;
		  tr->initialPartitionData[modelNumber].dataType   = DNA_DATA;
		  
		  found = TRUE;
		}	  
	      else
		{	    	  
		  if(strcasecmp(model, "BIN") == 0)
		    {	     	      
		      tr->initialPartitionData[modelNumber].protModels = -1;		  
		      tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = FALSE;
		      tr->initialPartitionData[modelNumber].dataType   = BINARY_DATA;
		      
		      found = TRUE;
		    }
		  else
		    {
		      if(strcasecmp(model, "MULTI") == 0)
			{	     	      
			  tr->initialPartitionData[modelNumber].protModels = -1;		  
			  tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = FALSE;
			  tr->initialPartitionData[modelNumber].dataType   = GENERIC_32;
			  
			  found = TRUE;
			}
		      else
			{
			  if(strcasecmp(model, "CODON") == 0)
			    {	     	      
			      tr->initialPartitionData[modelNumber].protModels = -1;		  
			      tr->initialPartitionData[modelNumber].usePredefinedProtFreqs  = FALSE;
			      tr->initialPartitionData[modelNumber].dataType   = GENERIC_64;
			      
			      found = TRUE;
			    }
			}
		    }
		  
		}
	    }

	  if(!found)
	    {
	      printf("ERROR: you specified the unknown model %s for partition %d\n", model, modelNumber);
	      exit(-1);
	    }
	}
	  
      i = 0;
      while(ident[i++] != ',');      
	  
      tr->initialPartitionData[modelNumber].partitionName = (char*)malloc((n - i + 1) * sizeof(char));          
	  
      j = 0;
      while(i < n)	
	tr->initialPartitionData[modelNumber].partitionName[j++] =  ident[i++];
      
      tr->initialPartitionData[modelNumber].partitionName[j] = '\0';                      
    }
}



static void setModel(int model, int position, int *a)
{
  if(a[position] == -1)
    a[position] = model;
  else
    {
      printf("ERROR trying to assign model %d to position %d \n", model, position);
      printf("while already model %d has been assigned to this position\n", a[position]);
      exit(-1);
    }      
}


static int myGetline(char **lineptr, int *n, FILE *stream)
{
  char *line, *p;
  int size, copy, len;
  int chunkSize = 256 * sizeof(char);

   if (*lineptr == NULL || *n < 2) 
    {
      line = (char *)realloc(*lineptr, chunkSize);
      if (line == NULL)
	return -1;
      *lineptr = line;
      *n = chunkSize;
    }

   line = *lineptr;
   size = *n;
  
   copy = size;
   p = line;
   
   while(1)
     {
       while (--copy > 0)
	 {
	   register int c = getc(stream);
	   if (c == EOF)
	     goto lose;
	   else
	     {
	       *p++ = c;
	       if(c == '\n' || c == '\r')	
		 goto win;
	     }
	 }

       /* Need to enlarge the line buffer.  */
       len = p - line;
       size *= 2;
       line = realloc (line, size);
       if (line == NULL)
	 goto lose;
       *lineptr = line;
       *n = size;
       p = line + len;
       copy = size - len;
     }
   
 lose:
  if (p == *lineptr)
    return -1;
  /* Return a partial line since we got an error in the middle.  */
 win:
  *p = '\0';
  return p - *lineptr;
}



void parsePartitions(analdef *adef, rawdata *rdta, tree *tr)
{
  FILE *f; 
  int numberOfModels = 0; 
  int nbytes = 0;
  char *ch;
  char *cc = (char *)NULL;
  char **p_names;
  int n, i, l;
  int lower, upper, modulo;
  char buf[256];
  int **partitions;
  int pairsCount;
  int as, j;
  int k; 

  f = myfopen(modelFileName, "rb");   

 
  while(myGetline(&cc, &nbytes, f) > -1)
    {     
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  numberOfModels++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }     
  
  rewind(f);
      
  p_names = (char **)malloc(sizeof(char *) * numberOfModels);
  partitions = (int **)malloc(sizeof(int *) * numberOfModels);
      
 
  
  tr->initialPartitionData = (pInfo*)malloc(sizeof(pInfo) * numberOfModels);

      
  for(i = 0; i < numberOfModels; i++) 
    {     
      tr->initialPartitionData[i].protModels = adef->proteinMatrix;
      tr->initialPartitionData[i].usePredefinedProtFreqs  = adef->protEmpiricalFreqs;
      tr->initialPartitionData[i].dataType   = -1;
    }

  for(i = 0; i < numberOfModels; i++)    
    partitions[i] = (int *)NULL;
    
  i = 0;
  while(myGetline(&cc, &nbytes, f) > -1)
    {          
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  n = strlen(cc);	 
	  p_names[i] = (char *)malloc(sizeof(char) * (n + 1));
	  strcpy(&(p_names[i][0]), cc);
	  i++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }         

  for(i = 0; i < numberOfModels; i++)
    {           
      ch = p_names[i];     
      pairsCount = 0;
      skipWhites(&ch);
      
      if(*ch == '=')
	{
	  printf("Identifier missing prior to '=' in %s\n", p_names[i]);
	  exit(-1);
	}
      
      analyzeIdentifier(&ch, i, tr);
      ch++;
            
    numberPairs:
      pairsCount++;
      partitions[i] = (int *)realloc((void *)partitions[i], (1 + 3 * pairsCount) * sizeof(int));
      partitions[i][0] = pairsCount;
      partitions[i][3 + 3 * (pairsCount - 1)] = -1; 	
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  exit(-1);
	}   
      
      l = 0;
      while(isNum(*ch))		 
	{
	  /*printf("%c", *ch);*/
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      lower = atoi(buf);
      partitions[i][1 + 3 * (pairsCount - 1)] = lower;   
      
      skipWhites(&ch);
      
      /* NEW */
      
      if((*ch != '-') && (*ch != ','))
	{
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {
	      upper = lower;
	      goto SINGLE_NUMBER;
	    }
	  else
	    {
	      printf("'-' or ',' expected in %s\n", p_names[i]);
	      exit(-1);
	    }
	}	 
      
      if(*ch == ',')
	{	     
	  upper = lower;
	  goto SINGLE_NUMBER;
	}
      
      /* END NEW */
      
      ch++;   
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  exit(-1);
	}    
      
      l = 0;
      while(isNum(*ch))
	{    
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      upper = atoi(buf);     
    SINGLE_NUMBER:
      partitions[i][2 + 3 * (pairsCount - 1)] = upper;        	  
      
      if(upper < lower)
	{
	  printf("Upper bound %d smaller than lower bound %d for this partition: %s\n", upper, lower,  p_names[i]);
	  exit(-1);
	}
      
      skipWhites(&ch);
      
      if(*ch == '\0' || *ch == '\n' || *ch == '\r') /* PC-LINEBREAK*/
	{    
	  goto parsed;
	}
      
      if(*ch == ',')
	{	 
	  ch++;
	  goto numberPairs;
	}
      
      if(*ch == '\\')
	{
	  ch++;
	  skipWhites(&ch);
	  
	  if(!isNum(*ch))
	    {
	      printf("%c Number expected in %s\n", *ch, p_names[i]);
	      exit(-1);
	    }     
	  
	  l = 0;
	  while(isNum(*ch))
	    {
	      buf[l] = *ch;
	      ch++;	
	      l++;
	    }
	  buf[l] = '\0';
	  modulo = atoi(buf);      
	  partitions[i][3 + 3 * (pairsCount - 1)] = modulo; 	
	  
	  skipWhites(&ch);
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {	     
	      goto parsed;
	    }
	  if(*ch == ',')
	    {	       
	      ch++;
	      goto numberPairs;
	    }
	}  
      
      assert(0);
       
    parsed:
      i = i;
    }
  
  fclose(f);
 
  /*********************************************************************************************************************/ 

  for(i = 0; i <= rdta->sites; i++)
    tr->model[i] = -1;
  
  for(i = 0; i < numberOfModels; i++)
    {   
      as = partitions[i][0];     
      
      for(j = 0; j < as; j++)
	{
	  lower = partitions[i][1 + j * 3];
	  upper = partitions[i][2 + j * 3]; 
	  modulo = partitions[i][3 + j * 3];	
	 
	  if(modulo == -1)
	    {
	      for(k = lower; k <= upper; k++)
		setModel(i, k, tr->model);
	    }
	  else
	    {
	      for(k = lower; k <= upper; k += modulo)
		{
		  if(k <= rdta->sites)
		    setModel(i, k, tr->model);	      
		}
	    }
	}        
    }


  for(i = 1; i < rdta->sites + 1; i++)
    {
      
      if(tr->model[i] == -1)
	{
	  printf("ERROR: Alignment Position %d has not been assigned any model\n", i);
	  exit(-1);
	}      
    }  

  for(i = 0; i < numberOfModels; i++)
    {
      free(partitions[i]);
      free(p_names[i]);
    }
  
  free(partitions);
  free(p_names);    
    
  tr->NumberOfModels = numberOfModels;     
  
  if(adef->perGeneBranchLengths)
    {
      if(tr->NumberOfModels > NUM_BRANCHES)
	{
	  printf("You are trying to use %d partitioned models for an individual per-gene branch length estimate.\n", tr->NumberOfModels);
	  printf("Currently only %d are allowed to improve efficiency.\n", NUM_BRANCHES);
	  printf("\n");
	  printf("In order to change this please replace the line \"#define NUM_BRANCHES   %d\" in file \"axml.h\" \n", NUM_BRANCHES);
	  printf("by \"#define NUM_BRANCHES   %d\" and then re-compile RAxML.\n", tr->NumberOfModels);
	  exit(-1);
	}
      else
	{
	  tr->multiBranch = 1;
	  tr->numBranches = tr->NumberOfModels;
	}
    }
}

/*******************************************************************************************************************************/

void handleExcludeFile(tree *tr, analdef *adef, rawdata *rdta)
{
  FILE *f;  
  char buf[256];
  int
    ch,
    j, value, i,
    state = 0,
    numberOfModels = 0,
    l = -1,
    excludeRegion   = 0,
    excludedColumns = 0,
    modelCounter    = 1;
  int
    *excludeArray, *countArray, *modelList;
  int
    **partitions;

  printf("\n\n");

  f = myfopen(excludeFileName, "rb");    

  while((ch = getc(f)) != EOF)
    {
      if(ch == '-')
	numberOfModels++;
    } 

  excludeArray = (int*)malloc(sizeof(int) * (rdta->sites + 1));
  countArray   = (int*)malloc(sizeof(int) * (rdta->sites + 1));
  modelList    = (int *)malloc((rdta->sites + 1)* sizeof(int));

  partitions = (int **)malloc(sizeof(int *) * numberOfModels);  
  for(i = 0; i < numberOfModels; i++)
    partitions[i] = (int *)malloc(sizeof(int) * 2);

  rewind(f);
  
  while((ch = getc(f)) != EOF)
    {     
      switch(state)
	{
	case 0: /* get first number */
	  if(!whitechar(ch))
	    {
	      if(!isNum(ch))
		{
		  printf("exclude file must have format: number-number [number-number]*\n");
		  exit(-1);
		}
	      l = 0;
	      buf[l++] = ch;
	      state = 1;
	    }
	  break;
	case 1: /*get the number or detect - */
	  if(!isNum(ch) && ch != '-')
	    {
	      printf("exclude file must have format: number-number [number-number]*\n");
	      exit(-1);
	    }
	  if(isNum(ch))
	    {
	      buf[l++] = ch;
	    }
	  else
	    {
	      buf[l++] = '\0';	     
	      value = atoi(buf);
	      partitions[excludeRegion][0] = value;
	      state = 2;
	    }
	  break;
	case 2: /*get second number */
	  if(!isNum(ch))
	    {
	      printf("exclude file must have format: number-number [number-number]*\n");
	      exit(-1);
	    }
	  l = 0;
	  buf[l++] = ch;
	  state = 3;
	  break;
	case 3: /* continue second number or find end */	 
	  if(!isNum(ch) && !whitechar(ch))
	    {
	      printf("exclude file must have format: number-number [number-number]*\n");
	      exit(-1);
	    }
	  if(isNum(ch))
	    {
	      buf[l++] = ch;
	    }
	  else
	    {	      
	      buf[l++] = '\0';	     
	      value = atoi(buf);
	      partitions[excludeRegion][1] = value;
	      excludeRegion++;
	      state = 0;
	    }
	  break;
	default:
	  assert(0);
	}
    }
     
  if(state == 3)
    {
      buf[l++] = '\0';     
      value = atoi(buf);
      partitions[excludeRegion][1] = value;
      excludeRegion++;
    }
  
  assert(excludeRegion == numberOfModels);

  for(i = 0; i <= rdta->sites; i++)
    {
      excludeArray[i] = -1;
      countArray[i] = 0;      
      modelList[i] = -1;
    }  

  for(i = 0; i < numberOfModels; i++)
    {
      int lower = partitions[i][0];
      int upper = partitions[i][1];

      if(lower > upper)
	{
	  printf("Misspecified exclude region %d\n", i);
	  printf("lower bound %d is greater than upper bound %d\n", lower, upper);
	  exit(-1);
	}

      if(lower == 0)
	{
	  printf("Misspecified exclude region %d\n", i);
	  printf("lower bound must be greater than 0\n");
	  exit(-1);
	}

      if(upper > rdta->sites)
	{
	  printf("Misspecified exclude region %d\n", i);
	  printf("upper bound %d must be smaller than %d\n", upper, (rdta->sites + 1));
	  exit(-1);
	}	
      for(j = lower; j <= upper; j++)
	{
	  if(excludeArray[j] != -1)
	    {
	      printf("WARNING: Exclude regions %d and %d overlap at position %d (already excluded %d times)\n", 
		     excludeArray[j], i, j, countArray[j]);
	    }
	  excludeArray[j] = i;
	  countArray[j]   =  countArray[j] + 1;	 
	}
    }

  for(i = 1; i <= rdta->sites; i++)
    {
      if(excludeArray[i] != -1)
	excludedColumns++;
      else
	{
	  modelList[modelCounter] = tr->model[i];
	  modelCounter++;
	}
    }

  printf("You have excluded %d out of %d columns\n", excludedColumns, rdta->sites);

  if(excludedColumns == rdta->sites)
    {
      printf("Error: You have excluded all sites\n");
      exit(-1);
    }

  if(adef->useSecondaryStructure && (excludedColumns > 0))
    {
      char mfn[2048];
      int countColumns;
      FILE *newFile;

      assert(adef->useMultipleModel);

      strcpy(mfn, secondaryStructureFileName);
      strcat(mfn, ".");
      strcat(mfn, excludeFileName);

      newFile = myfopen(mfn, "wb");

      printBothOpen("\nA secondary structure file with analogous structure assignments for non-excluded columns is printed to file %s\n", mfn);        	     	    
		  
      for(i = 1, countColumns = 0; i <= rdta->sites; i++)
	{		  
	  if(excludeArray[i] == -1)
	    fprintf(newFile, "%c", tr->secondaryStructureInput[i - 1]);
	  else
	    countColumns++;
	}
		  
      assert(countColumns == excludedColumns);
		  
      fprintf(newFile,"\n");
		  
      fclose(newFile);
    }


  if(adef->useMultipleModel && (excludedColumns > 0))
    {      
      char mfn[2048];
      FILE *newFile;

      strcpy(mfn, modelFileName);
      strcat(mfn, ".");
      strcat(mfn, excludeFileName);

      newFile = myfopen(mfn, "wb");

      printf("\nA partition file with analogous model assignments for non-excluded columns is printed to file %s\n", mfn);     
	      
      for(i = 0; i < tr->NumberOfModels; i++)
	{
	  boolean modelStillExists = FALSE;
		  
	  for(j = 1; (j <= rdta->sites) && (!modelStillExists); j++)
	    {
	      if(modelList[j] == i)
		modelStillExists = TRUE;
	    }

	  if(modelStillExists)
	    {	  	      
	      int k = 1;
	      int lower, upper;
	      int parts = 0;

	      switch(tr->partitionData[i].dataType)
		{
		case AA_DATA:		      		     
		  {
		    char AAmodel[1024];
		    
		    strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
		    if(tr->partitionData[i].usePredefinedProtFreqs == FALSE)
		      strcat(AAmodel, "F");		  
		    
		    fprintf(newFile, "%s, ", AAmodel);
		  }
		  break;
		case DNA_DATA:
		  fprintf(newFile, "DNA, ");
		  break;
		case BINARY_DATA:
		  fprintf(newFile, "BIN, ");
		  break;
		case GENERIC_32:
		  fprintf(newFile, "MULTI, ");
		  break;
		case GENERIC_64:
		  fprintf(newFile, "CODON, ");
		  break;
		default:
		  assert(0);
		}

	      fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);
	      
	      while(k <= rdta->sites)
		{
		  if(modelList[k] == i)
		    {
		      lower = k;
		      while((modelList[k + 1] == i) && (k <= rdta->sites))		      			
			k++;
		      upper = k;
		      
		      if(lower == upper)		  
			{
			  if(parts == 0)
			    fprintf(newFile, "%d", lower);
			  else
			    fprintf(newFile, ",%d", lower);
			}
		      else
			{
			  if(parts == 0)
			    fprintf(newFile, "%d-%d", lower, upper);
			  else
			    fprintf(newFile, ",%d-%d", lower, upper);
			}		  
		      parts++;
		    }
		  k++;
		}
	      fprintf(newFile, "\n");
	    }		  
	}	
      fclose(newFile);
    }

  
  {
    FILE *newFile;
    char mfn[2048];
   

    strcpy(mfn, seq_file);
    strcat(mfn, ".");
    strcat(mfn, excludeFileName);
    
    newFile = myfopen(mfn, "wb");
    
    printf("\nAn alignment file with excluded columns is printed to file %s\n\n\n", mfn);
    
    fprintf(newFile, "%d %d\n", tr->mxtips, rdta->sites - excludedColumns);
    
    for(i = 1; i <= tr->mxtips; i++)
      {   
	unsigned char *tipI =  &(rdta->y[i][1]);
	fprintf(newFile, "%s ", tr->nameList[i]);
	
	for(j = 0; j < rdta->sites; j++)
	  {
	    if(excludeArray[j + 1] == -1)	      
	      fprintf(newFile, "%c", getInverseMeaning(tr->dataVector[j + 1], tipI[j]));	       	  
	  }
	
	fprintf(newFile, "\n");
      }
    
    fclose(newFile);
  }

  
  fclose(f);
  for(i = 0; i < numberOfModels; i++)
    free(partitions[i]);
  free(partitions);  
  free(excludeArray);
  free(countArray);
  free(modelList);
}


void parseProteinModel(double *externalAAMatrix, char *fileName)
{
  FILE 
    *f = myfopen(fileName, "rb"); 
  
  int 
    doublesRead = 0,
    result = 0,
    i, 
    j;
  
  double 
    acc = 0.0;

  printf("\nParsing user-defined protein model from file %s\n", fileName);  

  while(doublesRead < 420)
    {     
      result = fscanf(f, "%lf", &(externalAAMatrix[doublesRead++]));           

      if(result == EOF)
	{
	  printf("Error protein model file must consist of exactly 420 entries \n");
	  printf("The first 400 entries are for the rates of the AA matrix, while the\n");
	  printf("last 20 should contain the empirical base frequencies\n");
	  printf("Reached End of File after %d entries\n", (doublesRead - 1));
	  exit(-1);
	}    
    }
       
  fclose(f);

  /* CHECKS */
  for(i = 0; i < 20; i++)
    for(j = 0; j < 20; j++)
      {
	if(i != j)
	  {
	    if(externalAAMatrix[i * 20 + j] != externalAAMatrix[j * 20 + i])
	      {
		printf("Error user-defined Protein model matrix must be symmetric\n");
		printf("Entry P[%d][%d]=%f at position %d is not equal to P[%d][%d]=%f at position %d\n", 
		       i, j,  externalAAMatrix[i * 20 + j], (i * 20 + j),
		       j, i,  externalAAMatrix[j * 20 + i], (j * 20 + i));
		exit(-1);
	      }
	  }
      }

  acc = 0.0;

  for(i = 400; i < 420; i++)    
    acc += externalAAMatrix[i];         

  if((acc > 1.0 + 1.0E-6) || (acc <  1.0 - 1.0E-6))
    {
      printf("Base frequencies in user-defined AA substitution matrix do not sum to 1.0\n");
      printf("the sum is %1.80f\n", acc);
      exit(-1);
    }
}




void parseSecondaryStructure(tree *tr, analdef *adef, int sites)
{
  if(adef->useSecondaryStructure)
    {
      FILE *f = myfopen(secondaryStructureFileName, "rb");

      int
	i,
	k,
	countCharacters = 0,
	ch,
	*characters,
	**brackets,
	opening,
	closing,
	depth,
	numberOfSymbols,
	numSecondaryColumns;      

      unsigned char bracketTypes[4][2] = {{'(', ')'}, {'<', '>'},  {'[', ']'},  {'{', '}'}};          

      numberOfSymbols = 4;     

      tr->secondaryStructureInput = (char*)malloc(sizeof(char) * sites);

      while((ch = fgetc(f)) != EOF)
	{
	  if(ch == '(' || ch == ')' || ch == '<' || ch == '>' || ch == '[' || ch == ']' || ch == '{' || ch == '}' || ch == '.')
	    countCharacters++;
	  else
	    {
	      if(!whitechar(ch))
		{
		  printf("Secondary Structure file %s contains character %c at position %d\n", secondaryStructureFileName, ch, countCharacters + 1);
		  printf("Allowed Characters are \"( ) < > [ ] { } \" and \".\" \n");
		  errorExit(-1);
		}
	    }
	}
      
      if(countCharacters != sites)
	{
	  printf("Error: Alignment length is: %d, secondary structure file has length %d\n", sites, countCharacters);
	  errorExit(-1);
	}
    
      characters = (int*)malloc(sizeof(int) * countCharacters); 

      brackets = (int **)malloc(sizeof(int*) * numberOfSymbols);
      
      for(k = 0; k < numberOfSymbols; k++)	  
	brackets[k]   = (int*)calloc(countCharacters, sizeof(int));

      rewind(f);

      countCharacters = 0;
      while((ch = fgetc(f)) != EOF)
	{
	  if(!whitechar(ch)) 
	    {
	      tr->secondaryStructureInput[countCharacters] = ch;
	      characters[countCharacters++] = ch;
	    }
	}
      
      assert(countCharacters == sites);

      for(k = 0; k < numberOfSymbols; k++)
	{
	  for(i = 0, opening = 0, closing = 0, depth = 0; i < countCharacters; i++)
	    {
	      if((characters[i] == bracketTypes[k][0] || characters[i] == bracketTypes[k][1]) && 
		 (tr->extendedDataVector[i+1] == AA_DATA || tr->extendedDataVector[i+1] == BINARY_DATA ||
		  tr->extendedDataVector[i+1] == GENERIC_32 || tr->extendedDataVector[i+1] == GENERIC_64))
		{
		  printf("Secondary Structure only for DNA character positions \n");
		  printf("I am at position %d of the secondary structure file and this is not part of a DNA partition\n", i+1);
		  errorExit(-1);
		}
	      
	      if(characters[i] == bracketTypes[k][0])
		{	      
		  depth++;
		  /*printf("%d %d\n", depth, i);*/
		  brackets[k][i] = depth;
		  opening++;
		}
	      if(characters[i] == bracketTypes[k][1])
		{	  
		  brackets[k][i] = depth; 
		  /*printf("%d %d\n", depth, i);  */
		  depth--;	
		  
		  closing++;
		}	  	  	  
	      
	      if(closing > opening)
		{
		  printf("at position %d there is a closing bracket too much\n", i+1);
		  errorExit(-1);
		}
	    }	

	  if(depth != 0)
	    {
	      printf("Problem: Depth: %d\n", depth);
	      printf("Your secondary structure file may be missing a closing or opening paraenthesis!\n");
	    }
	  assert(depth == 0);
	  
	  if(countCharacters != sites)
	    {
	      printf("Problem: sec chars: %d sites: %d\n",countCharacters, sites);
	      printf("The number of sites in the alignment does not match the length of the secondary structure file\n");
	    }
	  assert(countCharacters == sites);
	
      
	  if(closing != opening)
	    {
	      printf("Number of opening brackets %d should be equal to number of closing brackets %d\n", opening, closing);
	      errorExit(-1);
	    }
	}
      
      for(i = 0, numSecondaryColumns = 0; i < countCharacters; i++)
	{
	  int checkSum = 0;

	  for(k = 0; k < numberOfSymbols; k++)
	    {
	      if(brackets[k][i] > 0)
		{
		  checkSum++;
		  
		  switch(tr->secondaryStructureModel)
		    {
		    case SEC_16:
		    case SEC_16_A:
		    case SEC_16_B:
		    case SEC_16_C:
		    case SEC_16_D:
		    case SEC_16_E:
		    case SEC_16_F:
		    case SEC_16_I:
		    case SEC_16_J:
		    case SEC_16_K:
		      tr->extendedDataVector[i+1] = SECONDARY_DATA;
		      break;
		    case SEC_6_A:
		    case SEC_6_B:
		    case SEC_6_C:
		    case SEC_6_D:
		    case SEC_6_E:
		      tr->extendedDataVector[i+1] = SECONDARY_DATA_6;
		      break;
		    case SEC_7_A:
		    case SEC_7_B:
		    case SEC_7_C:
		    case SEC_7_D:
		    case SEC_7_E:
		    case SEC_7_F:
		      tr->extendedDataVector[i+1] = SECONDARY_DATA_7;
		      break;
		    default:
		      assert(0);
		    }
		 
		  numSecondaryColumns++;
		}
	    }
	  assert(checkSum <= 1);
	}
      
      assert(numSecondaryColumns % 2 == 0);
      
      /*printf("Number of secondary columns: %d merged columns: %d\n", numSecondaryColumns, numSecondaryColumns / 2);*/

      tr->numberOfSecondaryColumns = numSecondaryColumns;
      if(numSecondaryColumns > 0)
	{
	  int model = tr->NumberOfModels;
	  int countPairs;
	  pInfo *partBuffer = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);

	  for(i = 1; i <= sites; i++)
	    {
	      for(k = 0; k < numberOfSymbols; k++)
		{
		  if(brackets[k][i-1] > 0)
		    tr->model[i] = model;
		}

	    }

	  /* now make a copy of partition data */

	 
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      partBuffer[i].partitionName = (char*)malloc((strlen(tr->extendedPartitionData[i].partitionName) + 1) * sizeof(char));
	      strcpy(partBuffer[i].partitionName, tr->extendedPartitionData[i].partitionName);
	      strcpy(partBuffer[i].proteinSubstitutionFileName, tr->extendedPartitionData[i].proteinSubstitutionFileName);
	      partBuffer[i].dataType =  tr->extendedPartitionData[i].dataType;
	      partBuffer[i].protModels=  tr->extendedPartitionData[i].protModels;
	      partBuffer[i].usePredefinedProtFreqs=  tr->extendedPartitionData[i].usePredefinedProtFreqs;	      
	    }

	  for(i = 0; i < tr->NumberOfModels; i++)
	    free(tr->extendedPartitionData[i].partitionName);
	  free(tr->extendedPartitionData);
	 
	  tr->extendedPartitionData = (pInfo*)malloc(sizeof(pInfo) * (tr->NumberOfModels + 1));
	  
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      tr->extendedPartitionData[i].partitionName = (char*)malloc((strlen(partBuffer[i].partitionName) + 1) * sizeof(char));
	      strcpy(tr->extendedPartitionData[i].partitionName, partBuffer[i].partitionName);
	      strcpy(tr->extendedPartitionData[i].proteinSubstitutionFileName, partBuffer[i].proteinSubstitutionFileName);
	      tr->extendedPartitionData[i].dataType =  partBuffer[i].dataType;
	      tr->extendedPartitionData[i].protModels= partBuffer[i].protModels;
	      tr->extendedPartitionData[i].usePredefinedProtFreqs=  partBuffer[i].usePredefinedProtFreqs;	      
	      free(partBuffer[i].partitionName);
	    }
	  free(partBuffer);

	  tr->extendedPartitionData[i].partitionName = (char*)malloc(64 * sizeof(char));

	  switch(tr->secondaryStructureModel)
	    {
	    case SEC_16:
	    case SEC_16_A:
	    case SEC_16_B:
	    case SEC_16_C:
	    case SEC_16_D:
	    case SEC_16_E:
	    case SEC_16_F:
	    case SEC_16_I:
	    case SEC_16_J:
	    case SEC_16_K:
	      strcpy(tr->extendedPartitionData[i].partitionName, "SECONDARY STRUCTURE 16 STATE MODEL");
	      tr->extendedPartitionData[i].dataType = SECONDARY_DATA;
	      break;
	    case SEC_6_A:
	    case SEC_6_B:
	    case SEC_6_C:
	    case SEC_6_D:
	    case SEC_6_E:
	      strcpy(tr->extendedPartitionData[i].partitionName, "SECONDARY STRUCTURE 6 STATE MODEL");
	      tr->extendedPartitionData[i].dataType = SECONDARY_DATA_6;
	      break;
	    case SEC_7_A:
	    case SEC_7_B:
	    case SEC_7_C:
	    case SEC_7_D:
	    case SEC_7_E:
	    case SEC_7_F:
	      strcpy(tr->extendedPartitionData[i].partitionName, "SECONDARY STRUCTURE 7 STATE MODEL");
	      tr->extendedPartitionData[i].dataType = SECONDARY_DATA_7;
	      break;
	    default:
	      assert(0);
	    }

	  tr->extendedPartitionData[i].protModels= -1;
	  tr->extendedPartitionData[i].usePredefinedProtFreqs=  FALSE;	 

	  tr->NumberOfModels++;	 
	  
	  if(adef->perGeneBranchLengths)
	    {
	      if(tr->NumberOfModels > NUM_BRANCHES)
		{
		  printf("You are trying to use %d partitioned models for an individual per-gene branch length estimate.\n", tr->NumberOfModels);
		  printf("Currently only %d are allowed to improve efficiency.\n", NUM_BRANCHES);
		  printf("Note that the number of partitions has automatically been incremented by one to accomodate secondary structure models\n");
		  printf("\n");
		  printf("In order to change this please replace the line \"#define NUM_BRANCHES   %d\" in file \"axml.h\" \n", NUM_BRANCHES);
		  printf("by \"#define NUM_BRANCHES   %d\" and then re-compile RAxML.\n", tr->NumberOfModels);
		  exit(-1);
		}
	      else
		{
		  tr->multiBranch = 1;
		  tr->numBranches = tr->NumberOfModels;
		}
	    }
	  
	  assert(countCharacters == sites);

	  tr->secondaryStructurePairs = (int*)malloc(sizeof(int) * countCharacters);
	  for(i = 0; i < countCharacters; i++)
	    tr->secondaryStructurePairs[i] = -1;
	  /*
	    for(i = 0; i < countCharacters; i++)
	    printf("%d", brackets[i]);
	    printf("\n");
	  */
	  countPairs = 0;

	  for(k = 0; k < numberOfSymbols; k++)
	    {
	      i = 0;
	      
	  
	      while(i < countCharacters)
		{
		  int 
		    j = i,
		    bracket = -1,
		    openBracket,
		    closeBracket;
		  
		  while(j < countCharacters && ((bracket = brackets[k][j]) == 0))
		    {
		      i++;
		      j++;
		    }

		  assert(bracket >= 0);

		  if(j == countCharacters)
		    {
		      assert(bracket == 0);
		      break;
		    }
	    
		  openBracket = j;
		  j++;
		  
		  while(bracket != brackets[k][j] && j < countCharacters)
		    j++;
		  assert(j < countCharacters);
		  closeBracket = j;
		  
		  assert(closeBracket < countCharacters && openBracket < countCharacters);

		  assert(brackets[k][closeBracket] > 0 && brackets[k][openBracket] > 0);
		  
		  /*printf("%d %d %d\n", openBracket, closeBracket, bracket);*/
		  brackets[k][closeBracket] = 0;
		  brackets[k][openBracket]  = 0;	      	      
		  countPairs++;
		  
		  tr->secondaryStructurePairs[closeBracket] = openBracket;
		  tr->secondaryStructurePairs[openBracket] = closeBracket;
		}

	      assert(i == countCharacters);
	    }
	     
	  assert(countPairs == numSecondaryColumns / 2);

	  
	  /*for(i = 0; i < countCharacters; i++)
	    printf("%d ", tr->secondaryStructurePairs[i]);
	    printf("\n");*/
	  

	  adef->useMultipleModel = TRUE;

	}

      
      for(k = 0; k < numberOfSymbols; k++)	  
	free(brackets[k]);
      free(brackets);
      free(characters);

      fclose(f);
    }
}
