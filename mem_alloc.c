
#define MEM_ALLOC_NO_GUARDS 1

#include "mem_alloc.h"

#include <stdio.h>
#include <stdlib.h>

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

void *rax_memalign(size_t align, size_t size) 
{
#if defined (__APPLE__)
  return malloc(size); // apple has no memalign, but seem to return 16byte (32byte?) aligned blocks by default
#else
  return memalign(align, size);
#endif
    
}

void *rax_malloc( size_t size ) 
{
  return malloc(size);
}
void *rax_realloc( void *p, size_t size ) 
{
  return realloc(p, size);
}


void rax_free(void *p) 
{
  free(p);
}

int rax_posix_memalign(void **p, size_t align, size_t size) 
{
  return posix_memalign(p, align, size);
}
void *rax_calloc(size_t n, size_t size) 
{
  return calloc(n,size);
}

void *rax_malloc_aligned(size_t size) 
{
  const size_t BYTE_ALIGNMENT = 32;
  return rax_memalign(BYTE_ALIGNMENT, size);  
}

#endif



