#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>

#include <omp.h>

#ifndef MAT_PTR_TYPE
#define MAT_PTR_TYPE int
#endif
#ifndef MAT_VAL_TYPE
#define MAT_VAL_TYPE double
#endif

#ifndef PRE_BENCH_REPEAT
#define PRE_BENCH_REPEAT 5
#endif

#ifndef BENCH_REPEAT
#define BENCH_REPEAT 20
#endif

#define SUBSTITUTION_FORWARD 0
#define SUBSTITUTION_BACKWARD 1

#ifndef BENCH_REPEAT
#define BENCH_REPEAT 1
#endif

#ifndef BENCH_SM_REPEAT
#define BENCH_SM_REPEAT 1
#endif

#ifndef MAT_PTR_TYPE
#define MAT_PTR_TYPE int
#endif

#ifndef BENCH_REPEAT_TRSV
#define BENCH_REPEAT_TRSV 1
#endif

#ifndef BENCH_REPEAT_TRSM
#define BENCH_REPEAT_TRSM 1
#endif

#ifndef BENCH_REPEAT_TRSV_COO
#define BENCH_REPEAT_TRSV_COO 1
#endif
#ifndef BENCH_REPEAT_MV
#define BENCH_REPEAT_MV 1
#endif

#ifndef BENCH_REPEAT_MM
#define BENCH_REPEAT_MM 1
#endif

#ifndef THREAD_MAX
#define THREAD_MAX 64
#endif

#ifndef THREAD_MAX_TEST
#define THREAD_MAX_TEST 24
#endif

#ifndef REPEAT_GEMM
#define REPEAT_GEMM 1
#endif

#ifndef REPEAT_GEMM_INSP
#define REPEAT_GEMM_INSP 1
#endif

#ifndef repeat
#define repeat 1
#endif

#ifndef repeat_inp
#define repeat_inp 1
#endif

// struct timeval t1, t2;
