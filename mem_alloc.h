#ifndef __mem_alloc_h
#define __mem_alloc_h
#include <stddef.h>

typedef  int boolean;

void *rax_malloc(size_t size);
void *rax_realloc(void *p, size_t size, boolean needsMemoryAlignment);
void  rax_free(void *p);
void *rax_calloc(size_t n, size_t size);




#endif
