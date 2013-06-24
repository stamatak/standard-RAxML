
#define MEM_ALLOC_NO_GUARDS 1

#include "axml.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if !defined(__APPLE__)
#include <malloc.h>
#endif


#ifdef _RAXML_USE_LLALLOC

// the llalloc library implementation in lockless_alloc/ll_alloc.c exports the alloction functions prefixed
// with 'llalloc'. The following are the forward declarations of the llalloc* functions 

#define PREFIX(X)   llalloc##X

void *PREFIX(memalign)(size_t align, size_t size);
void *PREFIX(malloc)(size_t size);
void *PREFIX(realloc)(void *p, size_t size);
int PREFIX(posix_memalign)(void **p, size_t align, size_t size);
void *PREFIX(calloc)(size_t n, size_t size);
void PREFIX(free)(void *p);


// wrappers that forward the rax_* functions to the corresponding llalloc* functions


void *rax_memalign(size_t align, size_t size) 
{
  return PREFIX(memalign)(align, size);
}

void *rax_malloc( size_t size ) 
{
  return PREFIX(malloc)(size);
}
void *rax_realloc( void *p, size_t size ) {
  return PREFIX(realloc)(p, size);
}


void rax_free(void *p) {
  PREFIX(free)(p);
}

int rax_posix_memalign(void **p, size_t align, size_t size) {
  return PREFIX(posix_memalign)(p, align, size);
}
void *rax_calloc(size_t n, size_t size) {
  return PREFIX(calloc)(n,size);
}

void *rax_malloc_aligned(size_t size) 
{
  const size_t BYTE_ALIGNMENT = 32;
  return rax_memalign(BYTE_ALIGNMENT, size);
  
}

#else 
// if llalloc should not be used, forward the rax_* functions to the corresponding standard function



void *rax_malloc( size_t size ) 
{
  void 
    *ptr = (void *)NULL;
  
  int 
    res = posix_memalign(&ptr, BYTE_ALIGNMENT, size);

  if(res != 0) 
    assert(0);
  
  return ptr;
}

void *rax_realloc(void *p, size_t size, boolean needsMemoryAlignment) 
{  
  //it's actually not that easy to implement an aligned realloc
  //because we need to know the size of the array pointed to by 
  //the pointer passed as argument
  //hence I added this boolean flag that should increase programmer 
  //awareness about this issue
  
  if(needsMemoryAlignment)
    {
      assert(0);
      return (void*)NULL;
    }
  else
    return realloc(p, size);
}

void rax_free(void *p) 
{
  free(p);
}

void *rax_calloc(size_t n, size_t size) 
{
   void 
     *ptr = rax_malloc(size * n);

   if (ptr == (void*)NULL) 
     return (void*)NULL;

   memset(ptr, 0, size * n);

   return ptr;
}



#endif



