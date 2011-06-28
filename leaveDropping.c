#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "axml.h"


/* the define below makes the algorithm run much faster and should therefore always be used :-) */

#define _SORTED

/*
  The define below is usually not required, Andre just 
  used it, to have a one to one comparison/verification with Nick's original 
  python script. 
  The only thing this changes is that taxa are sorted in exactly the same way as 
  in Nick's script. In some cases, since we are dealing with a greedy algorithm 
  the taxon order has an impact on the dropsets size and composition.


  #define _COMPARE_TO_NICK

*/

extern char run_id[128];
extern char workdir[1024];
extern double masterTime;
extern const unsigned int mask32[];



typedef struct hash_el
{
  unsigned int fullKey;
  void *value;
  struct hash_el *next;
} HashElem;

typedef struct hash_table
{
  unsigned int tableSize;
  unsigned int entryCount;
  void *commonAttributes;
  unsigned int (*hashFunction)(struct hash_table *h, void *value);
  boolean (*equalFunction)(struct hash_table *hashtable, void *entrA, void *entryB);
  HashElem **table;
} HashTable;

typedef struct 
{
  HashTable *hashTable;
  HashElem *hashElem;
  unsigned int index;
} HashTableIterator;



typedef struct profile_elem
{
  /* TODO should better be named taxaSet and treeSet */
  unsigned int *bitVector;
  unsigned int *treeVector;
  unsigned int treeVectorSupport;
  unsigned int id;
} ProfileElem;

typedef struct
{
  unsigned int *dropset;
  unsigned int sizeOfDropset;
  unsigned int numberOfMergingPairs; 
} Dropset;

typedef struct 
{
  unsigned int dropsetVectorLength;
  unsigned int bipartitionVetorLength;
  unsigned int *randForTaxa;
} DropsetAttr;

typedef struct 
{
  unsigned int bitVectorLength; 
  unsigned int treeVectorLength;  
  unsigned int *randForTaxa;	/* random numbers to hash the vectors */
  unsigned int lastByte;		/* the padding bits */
} ProfileElemAttr;


typedef struct 
{
  void *arrayTable;
  void *commonAttributes; 
  unsigned int length;
} Array;




#define BIT_COUNT(x)  precomputed16_bitcount(x)

#define FLIP_NTH_BIT(bitVector,n) (bitVector[n / MASK_LENGTH] |= mask32[ n % MASK_LENGTH ])
#define UNFLIP_NTH_BIT(bitVeoctor,n) (bitVector[n / MASK_LENGTH] &= ~mask32[ n % MASK_LENGTH ])
#define NTH_BIT_IS_SET(bitVector,n) (bitVector[n / MASK_LENGTH] & mask32[n % MASK_LENGTH])


typedef struct 
{
  List *elems;
  unsigned int length;
} SparseBitVector;

/* SWITCHES */

#define PROFILE_TABLE(x,y) createHashTable(x,y,randomProfileElemHashFunction, bitVectorEqual)
#define DROPSET_TABLE(x,y) createHashTable(x,y,randomDropsetHashFunctionSparse, dropsetEqualFunctionSparse)
#define SORT_ARRAY(array, elemType, sortFunction) qsort(array->arrayTable, array->length, sizeof(elemType*), sortFunction)
#define HASH_TABLE_SIZE(x) (x * LOG(x))

/* END */

unsigned int genericBitCount(unsigned int* bitVector, unsigned int bitVectorLength)
{
  unsigned int 
    i, 
    result = 0;

  for(i = 0; i < bitVectorLength; i++)
    result += BIT_COUNT(bitVector[i]);
  
  return result; 
}


#ifdef _ANDRE_UNUSED_FUNCTIONS 

/* stuff added by Andre that is not used however */

static unsigned int naiveHashFunctionAdapterVersion(unsigned int *bv, unsigned int length)
{
  unsigned int 
    i, 
    result = 0; 
  
  for( i = 0; i < length; ++i)
    result ^= bv[i];
  
  /* return enhanceHash(result); */
  
  return result; 
}



static void *searchHashTableSeq(HashTable *hashtable, void *value, unsigned int hashValue)
{
  unsigned int 
    position = hashValue % hashtable->tableSize;
  
  HashElem 
    *elem, 
    *result = (HashElem*)NULL; 

  for(elem = hashtable->table[position]; 
      elem; 
      elem = elem->next)
    {
      if(elem->fullKey  == hashValue && 
	 hashtable->equalFunction(hashtable, elem->value, value))
	{
	  result = elem->value; 
	  break;
	}
    }

  return result;
}

static Array* dropsetHashToArray(HashTable *dropsets)
{
  unsigned int 
    count = 0; 
  
  HashTableIterator 
    *hashTableIterator = createHashTableIterator(dropsets);
  
  Array 
    *array = calloc(1, sizeof(Array));
  
  array->arrayTable = calloc(dropsets->entryCount, sizeof(Dropset*));
  array->length = dropsets->entryCount;
  
  do
    {      
      ((Dropset**)array->arrayTable)[count++] = getCurrentValueFromHashTableIterator(hashTableIterator);
    } 
  while(hashTableIteratorNext(hashTableIterator));
  
  assert(count == dropsets->entryCount);
  
  free(hashTableIterator);
  
  return array;
}

static void printDropset(unsigned int *bitVector, unsigned int numberOfTaxa, char **nameList)
{
  unsigned int 
    i;
  
  printBothOpen("DROPSET: ");
  
  for(i = 0; i < numberOfTaxa; ++i)
    if(NTH_BIT_IS_SET(bitVector, i))
      printBothOpen("%s,", nameList[i+1]);
  
  printBothOpen("\t");
}

static unsigned int naiveDropsetHashFunction(HashTable *hashTable, void *value)
{
  unsigned int 
    i, 
    result = 0,
    length = ((DropsetAttr*)hashTable->commonAttributes)->dropsetVectorLength;
  
  Dropset 
    *dropset = (Dropset*)value;
  
  for(i = 0; i < length; ++i )
    result ^= dropset->dropset[i];
  
  return result;
}

static List* appendToList(void *value, List *list)
{  
  List 
    *listElem = calloc(1, sizeof(List));
  
  listElem->value = value;
  listElem->next = list;
  
  return listElem;
}

static void insertIntoHashTableSeq(HashTable *hashTable, void *value, unsigned int index)
{  
  HashElem
    *hashElem = calloc(1, sizeof(HashElem));
  
  hashElem->fullKey = index;
  
  index =  hashElem->fullKey % hashTable->tableSize;
  
  hashElem->value = value;
  hashElem->next = hashTable->table[index];
  hashTable->table[index] = hashElem;
}
#endif

static HashTable *createHashTable(unsigned int size, 
				  void *commonAttr, 
				  unsigned int (*hashFunction)(HashTable *hash_table, void *value),
				  boolean (*equalFunction)(HashTable *hash_table, void *entryA, void *entryB))
{  
  static const unsigned int 
    initTable[] = {64, 128, 256, 512, 1024, 2048, 4096, 
		   8192, 16384, 32768, 65536, 131072, 
		   262144, 524288, 1048576, 2097152,
		   4194304, 8388608, 16777216, 33554432, 
		   67108864, 134217728, 268435456, 
		   536870912, 1073741824, 2147483648U};
  
  HashTable 
    *hashTable = calloc(1, sizeof(HashTable));
  
  unsigned int
    tableSize,
    i,
    primeTableLength = sizeof(initTable)/sizeof(initTable[0]),
    maxSize = (unsigned int)-1;    
  
  hashTable->hashFunction = hashFunction;
  hashTable->equalFunction = equalFunction;
  hashTable->commonAttributes = commonAttr;
 
  assert(size <= maxSize);  
  for(i = 0; initTable[i] < size && i < primeTableLength; ++i);
  assert(i < primeTableLength);

  tableSize = initTable[i];

  hashTable->table = calloc(tableSize, sizeof(HashElem*));
  hashTable->tableSize = tableSize;  
  hashTable->entryCount = 0;  

  return hashTable;
}


/* NOTE: computing hashvalue outside...cannot afford to compute it
   twice for bit vectors */

static void *searchHashTable(HashTable *hashtable, void *value, unsigned int hashValue)
{
  unsigned int 
    position = hashValue % hashtable->tableSize;
  
  HashElem 
    *elem, 
    *result = (HashElem*)NULL; 


  
  for(elem = hashtable->table[position]; 
      elem; 
      elem = elem->next)
    {
      if(elem->fullKey  == hashValue && 
	 hashtable->equalFunction(hashtable, elem->value, value))
	{
	  result = elem->value; 
	  break;
	}
    }
  


  return result;
}







static void insertIntoHashTable(HashTable *hashTable, void *value, unsigned int index)
{  
  HashElem 
    *hashElem = calloc(1, sizeof(HashElem));
  
  hashElem->fullKey = index;
  
  index =  hashElem->fullKey % hashTable->tableSize;
    
  hashElem->value = value;
  hashElem->next = hashTable->table[index];
  hashTable->table[index] = hashElem;
}


static void destroyHashTable(HashTable *hashTable, void (*freeValue)(void *value))
{
  unsigned 
    int i; 
  
  HashElem 
    *elemA, 
    *elemB,
    **table = hashTable->table;
  
  for(i = 0; i < hashTable->tableSize; ++i)
    {
      elemA = table[i];
      while(elemA != NULL)
	{
	  elemB = elemA; 
	  elemA = elemA->next; 
	  if(freeValue)
	    freeValue(elemB->value);
	  free(elemB);
	}


      
    }

  free(hashTable->commonAttributes);  
  free(hashTable->table);
  free(hashTable);
}


static void updateEntryCount(HashTable *hashTable)
{
  unsigned int 
    i, 
    result = 0;

  for(i = 0; i < hashTable->tableSize; ++i)
    {
      HashElem 
	*elem = ((HashElem**)hashTable->table)[i];
      
      while(elem)
	{
	  result++;
	  elem = elem->next;
	}
    }

  hashTable->entryCount = result;
}


static HashTableIterator *createHashTableIterator(HashTable *hashTable)
{
  unsigned 
    int i; 
  
  HashTableIterator 
    *hashTableIterator = calloc(1, sizeof(HashTableIterator));
  
  hashTableIterator->hashTable = hashTable;
  hashTableIterator->hashElem = NULL;
  hashTableIterator->index = hashTable->tableSize;
  
  if(!hashTable->entryCount)
    return hashTableIterator;

  for(i = 0; i < hashTable->tableSize; ++i)
    {
      if(hashTable->table[i])
	{
	  hashTableIterator->hashElem = hashTable->table[i];
	  hashTableIterator->index = i;
	  break;
	}
    }
  
  return hashTableIterator;
}

static boolean hashTableIteratorNext(HashTableIterator *hashTableIterator)
{
  unsigned int 
    i, 
    tableSize = hashTableIterator->hashTable->tableSize;
  
  HashElem 
    *next = hashTableIterator->hashElem->next;
  
  if(next)
    {
      hashTableIterator->hashElem = next;
      return TRUE;
    }
  
  /* TODO test case for this! */
  /* TODO looks like it should be optimised ... =/ */
  
  i = hashTableIterator->index + 1;
  
  if(i > tableSize)
    {
      hashTableIterator->index = i;
      return FALSE;
    }
  
  while((next = hashTableIterator->hashTable->table[i]) == ((HashElem*)NULL))
    {
      if( ++i >= tableSize) 
	{
	  hashTableIterator->index = i; 
	  return FALSE;
	}
    }
  
  hashTableIterator->index = i;
  hashTableIterator->hashElem = next;

  return (next != ((HashElem*)NULL));
}


/* TODO what about performance of this? */

static void *getCurrentValueFromHashTableIterator(HashTableIterator *hashTableIterator)
{  
  return ((hashTableIterator->hashElem) 
	  ?  hashTableIterator->hashElem->value
	  : NULL);
}

/* TODO implement remove */
/* TODO hash table change */


static unsigned int randomProfileElemHashFunction(HashTable *hashTable, void *value)
{
  unsigned int 
    i, 
    result = 0,
    length = ((ProfileElemAttr*)hashTable->commonAttributes)->bitVectorLength;
     
  ProfileElem 
    *profileElem = (ProfileElem*)value;

  for(i = 0; i < length * MASK_LENGTH; ++i) 
    if(NTH_BIT_IS_SET(profileElem->bitVector, i))
      result ^= ((ProfileElemAttr*)hashTable->commonAttributes)->randForTaxa[i];
  
  return result; 
}


static boolean bitVectorEqual(HashTable *hashtable, void *entryA, void *entryB)
{
  unsigned int 
    i, 
    bitVectorLength = ((ProfileElemAttr*)hashtable->commonAttributes)->bitVectorLength; 
  
  unsigned int
    *a = ((ProfileElem*)entryA)->bitVector,
    *b = ((ProfileElem*)entryB)->bitVector;
  
  for(i = 0; i < bitVectorLength; ++i)
    {
      if(a[i] != b[i])
	return FALSE;
    }
  
  return TRUE;
}




static Array* profileToArray(HashTable *profile, boolean updateFrequencyCount, boolean assignIds)
{
  HashTableIterator* 
    hashTableIterator = createHashTableIterator(profile);
  
  Array 
    *result = calloc(1, sizeof(Array));
  
  unsigned int 
    count = 0;

  /* remember to always copy s.t. free() runs w/o problems */
  
  ProfileElemAttr 
    *profileElemAttr;
  
  result->commonAttributes = calloc(1, sizeof(ProfileElemAttr));
  result->commonAttributes = memcpy(result->commonAttributes, profile->commonAttributes, sizeof(ProfileElemAttr));
  profileElemAttr = result->commonAttributes;

  result->length = profile->entryCount;
  result->arrayTable = calloc(profile->entryCount, sizeof(ProfileElem*));

  if(!hashTableIterator)
    return result;
  
  do
    {
      ProfileElem *profileElem = getCurrentValueFromHashTableIterator(hashTableIterator);

      if(updateFrequencyCount)
	profileElem->treeVectorSupport = genericBitCount(profileElem->treeVector, profileElemAttr->treeVectorLength);

      if(assignIds)
	profileElem->id = count;
      
      ((ProfileElem**)result->arrayTable)[count] = profileElem;
      assert(profileElem->bitVector && profileElem->treeVector);
      count++;
    }
  while(hashTableIteratorNext(hashTableIterator));
  
  assert(count == profile->entryCount);
  
  free(hashTableIterator);
  
  return result;
}


/* THIS IS AN ADAPTER */

static Array *convertHashtableToArray(hashtable *oldHashtable, unsigned int bitVectorLength, unsigned int treeVectorLength, 
				      unsigned int *randForTaxa, unsigned int lastByte)
{ 
  unsigned int
    count = 0,
    i;
  
  entry 
    *ent, 
    *entTmp; 
  
  Array 
    *result; 
  
  ProfileElemAttr 
    *attr = calloc(1, sizeof(ProfileElemAttr));
  
  attr->bitVectorLength = bitVectorLength; 
  attr->treeVectorLength = treeVectorLength;  
  attr->randForTaxa = randForTaxa;
  attr->lastByte = lastByte;
  
  HashTable 
    *hashTable = PROFILE_TABLE(HASH_TABLE_SIZE(oldHashtable->entryCount), attr);
  
  for(i = 0; i < oldHashtable->tableSize; ++i)
    {
      ent = oldHashtable->table[i];
      while(ent)
	{
	  ProfileElem 
	    *profileElem = calloc(1, sizeof(ProfileElem));
	  
	  profileElem->bitVector = ent->bitVector;
	  profileElem->treeVector = ent->treeVector;

	  insertIntoHashTable(hashTable, profileElem, 
			      hashTable->hashFunction(hashTable, profileElem));
	  
	  entTmp = ent->next;
	  free(ent);
	  ent = entTmp;
	}
    }
  
  free(oldHashtable);

  updateEntryCount(hashTable);
   
  
  for(i = 0; i < hashTable->tableSize; ++i)
    {	
      HashElem 
	*elem = ((HashElem**)hashTable->table)[i];
      
      while(elem)
	{
	  ((ProfileElem*)elem->value)->id = count++;
	  elem = elem->next;
	}
    }
  
  assert(count == hashTable->entryCount);

  result = profileToArray(hashTable, TRUE, TRUE);
  destroyHashTable(hashTable, NULL);
  
  return result;
}


#ifdef _SORTED

static int sortProfileElems(const void *a, const void *b)
{
  unsigned int aFreq = (*(ProfileElem**)a)->treeVectorSupport,
    bFreq = (*(ProfileElem**)b)->treeVectorSupport;
  
  if(aFreq < bFreq)
    return 1;
  else if(bFreq < aFreq)
    return -1;
  else
    return 0;
}

#endif


static Array *getInfrequentBipartitions(Array *oldArray, unsigned int threshold) 
{
  Array 
    *array = calloc(1, sizeof(Array)); 
  
  ProfileElem** 
    oldArrayTable = oldArray->arrayTable;

  unsigned int 
    i, 
    numberInfrequent = 0, 
    count = 0; 

  for(i = 0; i < oldArray->length; ++i)
    numberInfrequent += (oldArrayTable[i]->treeVectorSupport <= threshold);

  array->commonAttributes = calloc(1, sizeof(ProfileElemAttr));
  memcpy(array->commonAttributes, oldArray->commonAttributes, sizeof(ProfileElemAttr));
  array->length = numberInfrequent;
  array->arrayTable = calloc(numberInfrequent, sizeof(ProfileElem*));
  
  for(i = 0; i < oldArray->length; ++i)
    {
      ProfileElem *elem = ((ProfileElem**)oldArray->arrayTable)[i];
      assert(elem->bitVector && elem->treeVector); 
      if( elem->treeVectorSupport <= threshold)
	((ProfileElem**)array->arrayTable)[count++] = elem; 
    }

  assert(count == numberInfrequent);  
  
#ifdef _SORTED
  SORT_ARRAY(array, ProfileElem, sortProfileElems);
#endif

  return array;
}


static boolean isUnionOfTreesAboveThreshold( const ProfileElemAttr *profileElemAttr, 
					     const ProfileElem *elemA, const ProfileElem *elemB, 
					     unsigned int frequencyThreshold)
{  
  unsigned int 
    i, 
    count = 0,
    length = profileElemAttr->treeVectorLength;  

  for(i = 0; i < length; ++i )
    count += BIT_COUNT(elemA->treeVector[i] | elemB->treeVector[i]);
  
  return (count > frequencyThreshold);
}





static Dropset *getBestDropset(HashTable *dropsets

#ifdef _COMPARE_TO_NICK
			       ,char **nameList, unsigned int numberOfTaxa
#endif
			       
			       )
{
  
#ifdef _COMPARE_TO_NICK
  unsigned int j;   
#endif
  
  int 
    bestImpact = 0;
  
  unsigned int 
    i;  
  
  Dropset 
    *bestDropset = (Dropset *)NULL;
  
  HashTableIterator 
    *hashTableIterator = createHashTableIterator(dropsets);
  
  if(!getCurrentValueFromHashTableIterator(hashTableIterator))
    return NULL;

  do
    {
      Dropset 
	*dropset = getCurrentValueFromHashTableIterator(hashTableIterator);
      
      unsigned int 
	droppedTaxa = dropset->sizeOfDropset,
	gainedBips = dropset->numberOfMergingPairs;
      
      int 
	impact = gainedBips - droppedTaxa;
      
      if(impact > bestImpact)
	{
	  printBothOpen("found better dropset [gained-dropped=impact] %d\t-\t%d\t=%d\n", gainedBips, droppedTaxa, impact);
	  bestImpact = impact;
	  bestDropset = dropset;
	}
      
#ifdef _COMPARE_TO_NICK
      else if(impact > 0 && impact == bestImpact)
	{	  
	  /* resolve by number of taxa */
	  unsigned int 
	    fromBest = bestDropset->sizeOfDropset,
	    fromCurrent = dropset->sizeOfDropset;
	  
	  if(fromCurrent < fromBest)
	    {
	      printBothOpen("[SIZE] ");
	      printDropset(dropset->dropset, numberOfTaxa, nameList);
	      printBothOpen(" is shorter than ");
	      printDropset(bestDropset->dropset, numberOfTaxa, nameList);
	      bestDropset = dropset;
	      printBothOpen("\n");
	    }
	  else if(fromBest == fromCurrent)
	    {
	      for(j = 0; j < fromCurrent; ++j)
		{		  
		  if( bestDropset->dropset[j] < dropset->dropset[j] )
		    {
		      printBothOpen("[LEX] ");
		      printDropset(bestDropset->dropset, numberOfTaxa, nameList);
		      printBothOpen("  less than ");
		      printDropset(dropset->dropset, numberOfTaxa, nameList);
		      printBothOpen("\n");
		      break;
		    }
		  else if( dropset->dropset[j] < bestDropset->dropset[j] )
		    {
		      printBothOpen("[LEX] ");
		      printDropset(dropset->dropset, numberOfTaxa, nameList);
		      printBothOpen("  less than ");
		      printDropset(bestDropset->dropset, numberOfTaxa, nameList);
		      bestDropset = dropset;
		      printBothOpen("\n");
		      break;
		    }
		}
	    }
	}
#endif
    } 
  while(hashTableIteratorNext(hashTableIterator));

  free(hashTableIterator);  
  
  /* we did not find anything */
  if(!bestDropset)
    return bestDropset;

  /* else convert to non-sparse representation */
  unsigned int 
    *bvPtr = bestDropset->dropset; 
  
  bestDropset->dropset = calloc(((DropsetAttr*)dropsets->commonAttributes)->dropsetVectorLength, sizeof(unsigned int));

  for(i = 0; i < bestDropset->sizeOfDropset; ++i)
    FLIP_NTH_BIT(bestDropset->dropset, bvPtr[i]);
  
  free(bvPtr);
  
  return bestDropset;
}
 
 
static void insertBipartitionPairDropset_helper(HashTable *dropsets,
						unsigned int *diff, unsigned int diffCount)
{
  Dropset 
    *dropset = calloc(1, sizeof(Dropset)),
    *foundInHashTable;
  unsigned int 
    i, 
    ctr = 0;
  
  dropset->dropset = calloc(diffCount, sizeof(unsigned int));
 
  for(i = 0; i < ((DropsetAttr*)dropsets->commonAttributes)->dropsetVectorLength * MASK_LENGTH; ++i)
    {
      if(NTH_BIT_IS_SET(diff, i))
	{
	  dropset->dropset[ctr] = i;
	  ctr++;
	  if(ctr == diffCount)
	    break;
	}
    }
  
  free(diff);

  assert(ctr == diffCount);

  dropset->numberOfMergingPairs = 1;
  dropset->sizeOfDropset = diffCount;
  unsigned int hashValue = dropsets->hashFunction(dropsets, dropset);  
  
  foundInHashTable = (Dropset*)searchHashTable(dropsets, dropset, hashValue);
  
  if(!foundInHashTable)
    insertIntoHashTable(dropsets, dropset, hashValue);
  else
    {
      foundInHashTable->numberOfMergingPairs++;
            
      free(dropset->dropset);
      free(dropset);
    }
}


static void insertBipartitionPairDropset(HashTable *dropsets,
					 const ProfileElemAttr  *profileElemAttr, 
					 const ProfileElem *elemA, const ProfileElem *elemB,
					 const unsigned int *droppedTaxa )
{
  unsigned int 
    i, 
    diffCount1 = 0, 
    diffCount2 = 0,
    length = profileElemAttr->bitVectorLength,
    *diff1 = calloc(length, sizeof(unsigned int)), 
    *diff2 = calloc(length, sizeof(unsigned int)),
    lastByte = profileElemAttr->lastByte;

  
  for(i = 0; i < length; ++i )
    {
      diff1[i] = ( elemA->bitVector[i] ^ elemB->bitVector[i] ) ;

      if(i == length - 1)
	diff2[i] = (elemA->bitVector[i] & ~ droppedTaxa[i])  ^ ~ (elemB->bitVector[i] | droppedTaxa[i] | lastByte); 
      else
	diff2[i] = (elemA->bitVector[i] & ~ droppedTaxa[i]) ^ ~ (elemB->bitVector[i] | droppedTaxa[i]); 
    }

  
  diffCount1 = genericBitCount(diff1,length);
  diffCount2 = genericBitCount(diff2,length);
  
  if (!(diffCount1 && diffCount2)) 
    printf("problem with bip %d and %d\n", elemA->id, elemB->id);
  
  assert( diffCount1 && diffCount2 ) ;  

  if(diffCount1 < diffCount2)
    {
      insertBipartitionPairDropset_helper(dropsets, diff1, diffCount1);
      free(diff2);
    }
  else if(diffCount1 > diffCount2)
    {
      insertBipartitionPairDropset_helper(dropsets, diff2, diffCount2);
      free(diff1);
    }
  else
    {
      insertBipartitionPairDropset_helper(dropsets, diff1, diffCount1);
      insertBipartitionPairDropset_helper(dropsets, diff2, diffCount2);
    }
}




static unsigned int randomDropsetHashFunctionSparse(HashTable *hashTable, void *value)
{
  Dropset 
    *dropset = (Dropset*)value;
  
  unsigned int 
    i, 
    result = 0,
    length = dropset->sizeOfDropset,
    *randForTaxa = ((DropsetAttr*)hashTable->commonAttributes)->randForTaxa;

  for(i = 0; i < length; ++i)
    result ^= randForTaxa[dropset->dropset[i]];

  return result; 
}


static boolean dropsetEqualFunctionSparse(HashTable *hashTable, void *entryA, void *entryB)
{
  unsigned int i,
    aLength = ((Dropset*)entryA)->sizeOfDropset,
    bLength = ((Dropset*)entryB)->sizeOfDropset;
  
  unsigned int
    *a = ((Dropset*)entryA)->dropset,
    *b = ((Dropset*)entryB)->dropset;
  
  if(aLength != bLength)
    return FALSE;
  
  for(i = 0; i < aLength; i++)
    if(a[i] != b[i])
      return FALSE;

  return TRUE;
}


static HashTable* potentialProfileDropsets(Array *infrequentBipartitions, 		
					   unsigned int frequencyThreshold,
					   unsigned int *droppedTaxa )
{
  DropsetAttr 
    *dropsetAttr = calloc(1, sizeof(DropsetAttr));
  
  dropsetAttr->dropsetVectorLength = ((ProfileElemAttr*)infrequentBipartitions->commonAttributes)->bitVectorLength;
  dropsetAttr->bipartitionVetorLength = infrequentBipartitions->length;
  dropsetAttr->randForTaxa = ((ProfileElemAttr*)infrequentBipartitions->commonAttributes)->randForTaxa;

  HashTable 
    *dropsets = DROPSET_TABLE(HASH_TABLE_SIZE(infrequentBipartitions->length), dropsetAttr);
  
  ProfileElem 
    *elemA, 
    *elemB;   
  
  unsigned int 
    firstCount, 
    secondCount; 
  
  for(firstCount = 0; firstCount < infrequentBipartitions->length; ++firstCount)
    {
      elemA = ((ProfileElem**)(infrequentBipartitions->arrayTable))[firstCount];
      for(secondCount = firstCount + 1; secondCount < infrequentBipartitions->length; ++secondCount)
	{
	  elemB = ((ProfileElem**)(infrequentBipartitions->arrayTable))[secondCount];
	  assert(elemB->treeVector && elemB->bitVector);

	  if( elemA->treeVectorSupport + elemB->treeVectorSupport > frequencyThreshold )
	    {
	      if( isUnionOfTreesAboveThreshold(infrequentBipartitions->commonAttributes, elemA, elemB, frequencyThreshold))
		insertBipartitionPairDropset(dropsets, infrequentBipartitions->commonAttributes, elemA, elemB, droppedTaxa);
	    }
#ifdef _SORTED
	  else
	    break;
#endif
	}
    }
  
  /* due to parallelisation, we have to compute the entry count a posteriori */
  /* TODO rewrite this step for the sequential code */

  updateEntryCount(dropsets);

  return dropsets; 
}


static void destroyProfileElem(ProfileElem *profileElem)
{
  free(profileElem->bitVector);
  free(profileElem->treeVector);
  free(profileElem);
}


static void destroyDropset(void *dropset_) 
{
  Dropset 
    *dropset = (Dropset*)dropset_;
  
  free(dropset->dropset);
  free(dropset);
}


static List* getListOfConsensusBips(Array *allBipartitions, unsigned int frequencyThreshold)
{
  List* 
    result = (List*)NULL; 	
  
  unsigned int 
    i;
  
  ProfileElem 
    *elem; 

  for(i = 0; i < allBipartitions->length; ++i)
    {
      elem = ((ProfileElem**)allBipartitions->arrayTable)[i];
      if(elem->treeVectorSupport > frequencyThreshold)
	{
	  List 
	    *tmp = calloc(1, sizeof(List));
	  
	  tmp->value = elem;
	  tmp->next = result;	
	  result = tmp;
	}
    }
  
  return result;
}
  

static unsigned int getLengthOfList(List* list)
{ 
  unsigned int 
    result = 0; 

  for( ; list ; list = list->next)
    result++;
    
  return result;
}





static HashTable *updateAndInsertElem(ProfileElem *elem, HashTable *hashTable, unsigned int *droppedTaxa, ProfileElemAttr* profileElemAttr)
{
  unsigned int 
    j; 
  
  ProfileElem 
    *foundInHashTable, *foundComplementInHashTable,
    *complement = calloc(1, sizeof(ProfileElem));
  
  unsigned int 
    lastByte = ((ProfileElemAttr*)hashTable->commonAttributes)->lastByte;
  
  complement->bitVector = calloc(profileElemAttr->bitVectorLength, sizeof(unsigned int));

  for(j = 0; j < profileElemAttr->bitVectorLength; ++j )
    {
      if(j == profileElemAttr->bitVectorLength - 1)
	{
	  elem->bitVector[j] &= ~ ( droppedTaxa[j] | lastByte );
	  complement->bitVector[j] = ~ (elem->bitVector[j] | droppedTaxa[j] | lastByte);
	}
      else
	{
	  elem->bitVector[j] &= ~ droppedTaxa[j];
	  complement->bitVector[j] = ~ (elem->bitVector[j] | droppedTaxa[j] );
	}
    }
  
  /* bipartition vanishes */
  
  unsigned int 
    numberOfElements = genericBitCount(elem->bitVector, profileElemAttr->bitVectorLength),
    numberElementsInComplement = genericBitCount(complement->bitVector, profileElemAttr->bitVectorLength); 
  
  if(numberOfElements < 2 || numberElementsInComplement < 2)
    {
      destroyProfileElem(elem);      
      free(complement->bitVector);
      free(complement);
      return hashTable;
    }
  
  unsigned int 
    hashValue = hashTable->hashFunction(hashTable, elem);

  foundInHashTable = searchHashTable(hashTable, elem, hashValue);
  
  if(!foundInHashTable)
    {
      unsigned int 
	complementHashValue = hashTable->hashFunction(hashTable, complement);
      
      foundComplementInHashTable = searchHashTable(hashTable, complement, complementHashValue);
      
      if(foundComplementInHashTable)
	foundInHashTable = foundComplementInHashTable;
    }

  if(foundInHashTable)
    {      
      for (j = 0; j < profileElemAttr->treeVectorLength; ++j )
	foundInHashTable->treeVector[j] |= elem->treeVector[j];      
      
      foundInHashTable->treeVectorSupport = genericBitCount(foundInHashTable->treeVector, profileElemAttr->treeVectorLength);
      
      destroyProfileElem(elem);
      free(complement->bitVector);
      free(complement);
    }
  else
    {
      insertIntoHashTable(hashTable, elem, hashValue);
      free(complement->bitVector);
      free(complement);
    }  

  return hashTable;
}


static Array* restrictProfile( Array *infrequentBipartitions, List *consensusBipartitions, unsigned int *droppedTaxa )
{
  unsigned int 
    i; 
  
  ProfileElemAttr 
    *profileElemAttr = calloc(1, sizeof(ProfileElemAttr));
  
  profileElemAttr = memcpy(profileElemAttr, infrequentBipartitions->commonAttributes, sizeof(ProfileElemAttr));
  
  HashTable 
    *hashTable = PROFILE_TABLE(HASH_TABLE_SIZE(infrequentBipartitions->length + getLengthOfList(consensusBipartitions)),
			       profileElemAttr);
  
  unsigned int 
    lengthOfConsensus = getLengthOfList(consensusBipartitions);
  
  Array 
    *result; 
  
  List 
    *listElem, 
    *listTmp;

  Array 
    *tmpArray = calloc(1, sizeof(Array));

  
  tmpArray->arrayTable = calloc(infrequentBipartitions->length + lengthOfConsensus, sizeof(ProfileElem*));
  tmpArray->commonAttributes = calloc(1, sizeof(ProfileElemAttr));
  tmpArray->commonAttributes = memcpy(tmpArray->commonAttributes, infrequentBipartitions->commonAttributes, sizeof(ProfileElemAttr));

  for(i = 0; i < infrequentBipartitions->length; ++i)
    ((ProfileElem**)tmpArray->arrayTable)[i] = ((ProfileElem**)infrequentBipartitions->arrayTable)[i];

  tmpArray->length = infrequentBipartitions->length;

  free(infrequentBipartitions->commonAttributes);
  free((ProfileElem**)infrequentBipartitions->arrayTable);
  free(infrequentBipartitions);

  listElem = consensusBipartitions;
  
  while(listElem)
    {
      ((ProfileElem**)tmpArray->arrayTable)[tmpArray->length] = (ProfileElem*)listElem->value;
      tmpArray->length++;
      listTmp = listElem;
      listElem = listElem->next;
      free(listTmp);
    }  
  
  for(i = 0; i < tmpArray->length; ++i)
    {
      hashTable = updateAndInsertElem(((ProfileElem**)tmpArray->arrayTable)[i],
  				      hashTable,
  				      droppedTaxa,
  				      profileElemAttr);
    }

  free(tmpArray->arrayTable);
  free(tmpArray->commonAttributes);
  free(tmpArray);

  updateEntryCount(hashTable);

  result = profileToArray(hashTable, TRUE, FALSE);

  destroyHashTable(hashTable, NULL);
  return result;
}


/**
   Computes in a greedy manner a set of rogue taxa. 

  @return the dropset
  @param  profile            -- the hashtable used in bipartitionList.c
  @param  tree               -- the standard tree
  @param resultBipartitions  -- bipartitions after dropping the rogue
                                taxa. necessary part of the output!
 */

static unsigned int* determineGreedyDropset(hashtable *profile, tree *tree, Array **resultBipartitions)
{    

  unsigned int 
    numberOfTaxa = tree->mxtips, 
    numberOfTrees = tree->numberOfTrees, 
    frequencyThreshold,    
    i,
    numberBipartitions,
    treeVectorLength = GET_BITVECTOR_LENGTH(numberOfTrees),
    bitVectorLength = GET_BITVECTOR_LENGTH(numberOfTaxa),
    *droppedTaxa = calloc(bitVectorLength, sizeof(unsigned int)); 
  
  boolean 
    droppedTaxaThisRound = FALSE;
  
  HashTable 
    *dropsets;  
  
  Dropset 
    *bestDropset; 
  
  Array 
    *allBipartitions;

 

  /* initialise random numbers: one for each taxon  */
  unsigned int 
    *randForTaxa = calloc(numberOfTaxa, sizeof(unsigned int));
  
  for(i = 0; i < numberOfTaxa; ++i)  
    randForTaxa[i] = rand();
  
  switch(tree->consensusType)
    {           
    case MR_CONSENSUS: 
      frequencyThreshold =  numberOfTrees / 2;
      break;
    case STRICT_CONSENSUS:
      frequencyThreshold =  numberOfTrees - 1;
      break;
    case MRE_CONSENSUS:
    default:    
      assert(0);
    }

  /* 
     prepare a compensator for the padding bits (bits that are
     necessarily in the bitvector but do not represent taxa)  
  */
  
  unsigned int 
    lastByte = 0;
  
  for(i = numberOfTaxa; i < MASK_LENGTH * bitVectorLength; ++i)
    lastByte |= mask32[i % MASK_LENGTH];

  assert(numberOfTrees > 0 && numberOfTaxa > 0 &&  bitVectorLength > 0 && treeVectorLength > 0 && frequencyThreshold > 0);
  
  /* 
     caling adapter: adapts the old hashtable to the new implementation 
  */  
  
  allBipartitions = convertHashtableToArray(profile, bitVectorLength, treeVectorLength, randForTaxa, lastByte);

  printBothOpen("tree->consensusType=%d, numberOfTaxa=%d, bitVectorLength=%d, numberOfTrees=%d, treeVectorLength=%d, frequencyThreshold=%d, profile->entryCount=%d\n\n",
		tree->consensusType, numberOfTaxa, bitVectorLength,
		numberOfTrees, treeVectorLength, frequencyThreshold,allBipartitions->length); 

 

  do 
    {
      printBothOpen("================================================================\n");

      
      Array 
	*infrequentBipartitions = getInfrequentBipartitions(allBipartitions, frequencyThreshold);
	  
      List 
	*consensusBipartitions = getListOfConsensusBips(allBipartitions, frequencyThreshold);
      
      assert(allBipartitions->length == infrequentBipartitions->length + getLengthOfList(consensusBipartitions)); 
    
      numberBipartitions = allBipartitions->length;
    
      free((ProfileElem**)allBipartitions->arrayTable);
      free(allBipartitions->commonAttributes);
      free(allBipartitions);
      
      /* 
	 NOTE: this below would be possible and might make sense in some
	 situations. But as there is no equivalent in Nick's script and
	 this somehow contradicts the philosophy stated in the paper (to
	 remove all the rogues), I comment it out.
	 It would of advantage, if the algorithm is used to get a better
	 resolved tree. However, when we stop via this criteria, rogues
	 might remain in the tree. If it is the target to yield a high
	 number of branches with the highest possible support, then it
	 would be better, if those rogues were dropped as well.
      */
      
      /* 
	 if(getLengthOfList(consensusBipartitions) >= numberOfTaxa - 3 - genericBitCount(droppedTaxa))
         {
	 printBothOpen("we have enough bips: %d / %d\n", getLengthOfList(consensusBipartitions), numberOfTaxa-3); 	 
	 break; 
         } 
	 else 
	 printBothOpen("we have %d bips, we need %d.\n", getLengthOfList(consensusBipartitions), numberOfTaxa-3); 
      */

       
      printBothOpen("divided bips: %d = %d infreq + %d consensus\n", allBipartitions->length, infrequentBipartitions->length, getLengthOfList(consensusBipartitions));

     
      dropsets =  potentialProfileDropsets(infrequentBipartitions, frequencyThreshold, droppedTaxa);
    
      

#ifdef _COMPARE_TO_NICK
      bestDropset = getBestDropset(dropsets , tree->nameList, numberOfTaxa );
#else
      bestDropset = getBestDropset(dropsets);
#endif
    
      droppedTaxaThisRound = (bestDropset != NULL);      
      
    
      if(droppedTaxaThisRound)
	{		
	  /* 
	     update dropped taxa 
	  */
	  
	  for(i = 0 ; i < ((DropsetAttr*)dropsets->commonAttributes)->dropsetVectorLength; ++i)
	    droppedTaxa[i] |= bestDropset->dropset[i];
	  
	  /* 
	     update profile array 
	  */
	  
	  allBipartitions = restrictProfile(infrequentBipartitions, consensusBipartitions, droppedTaxa );
	}
      else
	{
	  allBipartitions = restrictProfile( infrequentBipartitions, consensusBipartitions, droppedTaxa );
	  *resultBipartitions = allBipartitions;
	}
      
      /* just debug */
      
      if(droppedTaxaThisRound)
	{
	  printBothOpen("dropping taxa: " );
	  for(i = 0; i < numberOfTaxa; ++i)
	    if(NTH_BIT_IS_SET(bestDropset->dropset, i))
	      printBothOpen("%s,", tree->nameList[i+1]);
	  printBothOpen("\n");
	}	  
      
      destroyHashTable(dropsets, destroyDropset);            
    }      
  while(droppedTaxaThisRound);
  
  free(randForTaxa);

  return droppedTaxa;
}

static unsigned int renormalizeBipartitions(Array *bipartitions, unsigned int *droppedTaxa, unsigned int numberOfTaxa)
{
  unsigned 
    i, 
    resultIndex = 0,
    bitVectorLength = ((ProfileElemAttr*)bipartitions->commonAttributes)->bitVectorLength,
    lastByte = 0;
  

  for(i = numberOfTaxa; i < MASK_LENGTH * bitVectorLength; ++i)
    lastByte |= mask32[i % MASK_LENGTH];
  
  /* if the reference leave did not change, nothing needs to be done */
  
  if(!NTH_BIT_IS_SET(droppedTaxa, 0))
    return 0; 

  /* find adequate new leave */
  
  for(i = 0; i < numberOfTaxa; ++i)
    {
      if(!NTH_BIT_IS_SET(droppedTaxa,i))
	{
	  resultIndex = i;
	  break;
	}
    }
  
  /* is only the case if everything was dropped  */
  
  assert(resultIndex);
 
  /* invert bit-vectors according to the reference leave */
  
  for(i = 0; i < bipartitions->length; ++i)
    {
      ProfileElem 
	*elem =  ((ProfileElem**)bipartitions->arrayTable)[i];
      
      if(NTH_BIT_IS_SET(elem->bitVector,resultIndex))
	{
	  unsigned int 
	    j;
	  
	  for(j = 0; j < bitVectorLength; ++j)
	    {	      
	      elem->bitVector[j] = 
		(j == bitVectorLength - 1)
		? ~ (elem->bitVector[j] | droppedTaxa[j] | lastByte)
		: ~ (elem->bitVector[j] | droppedTaxa[j]);
	    }
	}
    }

  return resultIndex;
}

static hashtable *reconvertHashtable(Array *bipartitionArray)
{
  unsigned int
    i, 
    count = 0; 
  
  hashtable 
    *htable = initHashTable(bipartitionArray->length); 
  
  assert(htable->tableSize >= bipartitionArray->length );
  
  for(i = 0; i < bipartitionArray->length; ++i)
    {
      ProfileElem 
	*elem = ((ProfileElem**)bipartitionArray->arrayTable)[i];
      
      entry 
	*ent = initEntry();
      
      ent->bitVector = elem->bitVector;
      ent->treeVector = elem->treeVector;
      ent->bipNumber = count++;
      htable->table[i] = ent; 
      free(elem);
    }
  
  htable->entryCount = count;
  
  return htable;
}

static void pruneTaxon(tree *tr, unsigned int k)
{

  assert(k > 0 && k <= ((unsigned int)(tr->mxtips)));

  nodeptr 
    p = tr->nodep[k],    
    q = p->back,
    q1 = q->next->back,
    q2 = q->next->next->back;


  hookupDefault(q1, q2, tr->numBranches);

  tr->start = findAnyTip(q1, tr->mxtips);

  assert(p != tr->start && q != tr->start);
}
  



void computeRogueTaxa(tree *tr,  char* treeSetFileName, analdef *adef)
{
  char 
    dropFileName[1024];

  hashtable 
    *h = initHashTable(tr->mxtips * FC_INIT * 10);

  Array 
    **resultBipartitions = calloc(1, sizeof(Array*));

  unsigned int 
    droppedTaxaNum = 0,
    numberOfTrees = 0, 
    i,    
    treeVectorLength, 
    vectorLength,
    **bitVectors = initBitVector(tr, &vectorLength),
    *droppedTaxa,
    referenceLeaveIndex = 0;  

  FILE
    *treeFile = getNumberOfTrees(tr, treeSetFileName, adef),
    *outf;    

  numberOfTrees = tr->numberOfTrees;      

  assert(adef->leaveDropMode);
  assert(sizeof(unsigned char) == 1); 
  
  treeVectorLength = GET_BITVECTOR_LENGTH(numberOfTrees);

  /* read the trees and process the bipartitions */

  for(i = 1; i <= numberOfTrees; i++)
    {                  
      int 
	bCount = 0;
      
      treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE);
      
      assert(tr->mxtips == tr->ntips);
     
      bitVectorInitravSpecial(bitVectors, tr->nodep[1]->back, tr->mxtips, vectorLength, h, (i - 1), BIPARTITIONS_BOOTSTOP, (branchInfo *)NULL,
			      &bCount, treeVectorLength, FALSE, FALSE);
      
      assert(bCount == tr->mxtips - 3); 
    }  


  rewind(treeFile); 
               
  droppedTaxa = determineGreedyDropset(h, tr, resultBipartitions);
       

  referenceLeaveIndex = renormalizeBipartitions(*resultBipartitions, droppedTaxa, tr->mxtips);

  h = reconvertHashtable(*resultBipartitions);

  /* TODO output for egrep */
    
  printBothOpen("> ALL dropped taxa: ");
 
  for(i = 0; i < (unsigned int)tr->mxtips; ++i)
    {
      if(NTH_BIT_IS_SET(droppedTaxa, i))
	{
	  printBothOpen("%d %s,", (i+1), tr->nameList[i + 1]);
	  droppedTaxaNum++;
	}
    }
  printBothOpen("\n"); 
    

  printBothOpen("\nDropping %u taxa\n", droppedTaxaNum);

  printBothOpen("\nTime for dropset calculation: %f seconds\n", gettime() - masterTime );
    
  if(droppedTaxaNum > 0)
    {
      strcpy(dropFileName, workdir);
      strcat(dropFileName, "RAxML_prunedTrees.");
      strcat(dropFileName, run_id);
      
      outf = myfopen(dropFileName, "w");
      
      for(i = 1; i <= numberOfTrees; i++)
	{                  
	  unsigned int
	    k;
	  
	  int
	    tips = 0;
	  
	  /*printf("Tree %d\n", i);*/
	  
	  treeReadLen(treeFile, tr, FALSE, FALSE, TRUE, adef, TRUE);
	  
	  assert(tr->mxtips == tr->ntips);
	  
	  for(k = 0; k < (unsigned int)tr->mxtips; k++)
	    {
	      if(NTH_BIT_IS_SET(droppedTaxa, k))
		pruneTaxon(tr, (k+1));
	    }
	  
	  tips = countTips(tr->start, tr->mxtips) + countTips(tr->start->back, tr->mxtips);
	  
	  assert((unsigned)tips == ((unsigned)tr->mxtips - droppedTaxaNum));
	  
	  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, FALSE, adef, NO_BRANCHES, FALSE, FALSE);
	  fprintf(outf, "%s", tr->tree_string);
	  
	  /*printf("%u %u\n", (unsigned)tips, ((unsigned)tr->mxtips - droppedTaxaNum));*/
	}
      
      
      printBothOpen("\nA tree collection, where the taxa from the dropset have been pruned in each tree\n");
      printBothOpen("has been written to file %s\n\n", dropFileName);
      fclose(outf);
    }
      
  printBothOpen("Total execution time: %f\n", gettime() - masterTime);
  
  fclose(treeFile); 

  freeBitVectors(bitVectors, 2 * tr->mxtips);   
  free(droppedTaxa);
  free(bitVectors); 
  freeHashTable(h);
  free(h);

  exit(0);
}
