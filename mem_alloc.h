#ifndef __mem_alloc_h
#define __mem_alloc_h
#include <stddef.h>


void *rax_memalign(size_t align, size_t size);
void *rax_malloc(size_t size);
void *rax_realloc(void *p, size_t size);
void  rax_free(void *p);
int   rax_posix_memalign(void **p, size_t align, size_t size);
void *rax_calloc(size_t n, size_t size);
void *rax_malloc_aligned(size_t size);




#endif
