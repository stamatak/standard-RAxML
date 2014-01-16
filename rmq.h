#ifndef _rmq_h_
#define _rmq_h_

#include <math.h>

#define false 0
#define true 1

typedef int DT;                 // use long for 64bit-version (but take care of fast log!)
typedef unsigned int DTidx;     // for indexing in arrays
typedef unsigned char DTsucc;
typedef unsigned short DTsucc2;
typedef int bool;
DTidx query(DTidx, DTidx);
void RMQ_succinct(DT* a, DTidx n);
void RMQ_succinct_destroy(void);

#endif
