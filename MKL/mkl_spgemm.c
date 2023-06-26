#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "mkl.h"
#include "../include/biio.h"
#include "../include/common.h"

#ifdef MKL_ILP64
#define INT_PRINT_FORMAT "%lld"
#else
#define INT_PRINT_FORMAT "%d"
#endif

#define REPEAT 5
#define ALIGN 128

/* To avoid constantly repeating the part of code that checks inbound SparseBLAS functions' status,
   use macro CALL_AND_CHECK_STATUS */
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
          if(function != SPARSE_STATUS_SUCCESS)             \
          {                                                 \
          printf(error_message); fflush(0);                 \
          status = 1;                                       \
          goto memory_free;                                 \
          }                                                 \
} while(0)

int main(int argc,char **argv) {

    if (argc < 2)
    {
        printf("Run the code by './mkl_spmv matrix.cbd'.\n");
        return 0;
    }
    struct timeval t1, t2, t3, t4;
    char  *filename;
    filename = argv[1];
    char *name = (char*)malloc(30*sizeof(char));
    for(int i = strlen(filename);i > 0;i--){
        if(filename[i] == '/'){
            // printf("%s\n",filename + i);
            strncpy(name,filename + i + 1,strlen(filename) - i - 5);
            name[strlen(filename) - i - 5] = '\0';
            // name = filename + i + 1;
            break;
        }
    }
    int m;
    int n;
    int nnz;
    int isSymmetric;
    MAT_PTR_TYPE *csrRowPtr;
    int *csrColIdx;
    MAT_VAL_TYPE *csrVal;
    read_Dmatrix_32(&m, &n, &nnz, &csrRowPtr, &csrColIdx, &csrVal,&isSymmetric, filename);
//    printf("input matrix A: ( %i, %i ) nnz = %i\n", m, n, nnz);

/* Declaration of values */
    double  *values_C = NULL;
    MKL_INT *columns_C = NULL;
    MKL_INT *pointerB_C = NULL;
    MKL_INT *pointerE_C = NULL;
    double  *rslt_mv = NULL, *rslt_mv_trans = NULL, *x = NULL, *y = NULL;

    double   left, right, residual;
    MKL_INT  rows, cols, i, j, ii, status;

    sparse_index_base_t    indexing;
    struct matrix_descr    descr_type_gen;
    sparse_matrix_t        csrA = NULL, csrB = NULL, csrC = NULL;

/* Allocation of memory */
    x = (double *)mkl_malloc(sizeof(double) * n, ALIGN);
    y = (double *)mkl_malloc(sizeof(double) * m, ALIGN);
    rslt_mv = (double *)mkl_malloc(sizeof(double) * n, ALIGN);
    rslt_mv_trans = (double *)mkl_malloc(sizeof(double) * m, ALIGN);

/* Set values of the variables*/
    descr_type_gen.type = SPARSE_MATRIX_TYPE_GENERAL;
    /* The following decriptor fields are not applicable for matrix type
     * SPARSE_MATRIX_TYPE_GENERAL but must be specified for other matrix types:
     */
    // descr_type_gen.mode = SPARSE_FILL_MODE_FULL;
    // descr_type_gen.diag = SPARSE_DIAG_NON_UNIT;
    status = 0, ii = 0;

 //Vectors x and y
    // for( i = 0; i < m; i++ )
    // {
    //     y[i] = 1.0;
    // }
    // for( i = 0; i < m; i++ )
    // {
    //     x[i] = 1.0;
    // }

/* Prepare arrays, which are related to matrices.
   Create handles for matrices A and B stored in CSR format */
    CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &csrA, SPARSE_INDEX_BASE_ZERO, m, n, csrRowPtr, csrRowPtr+1, csrColIdx, csrVal),
                          "Error after MKL_SPARSE_D_CREATE_CSR, csrA \n");
    CALL_AND_CHECK_STATUS(mkl_sparse_convert_csr(csrA,SPARSE_OPERATION_TRANSPOSE, &csrB),

                           "Error after MKL_SPARSE_CONVERT_CSR, csrB \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &csrB, SPARSE_INDEX_BASE_ZERO, m, n, csrRowPtr, csrRowPtr+1, csrColIdx, csrVal),
    //                       "Error after MKL_SPARSE_D_CREATE_CSR, csrB \n");

    CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrA ),
                          "Error after MKL_SPARSE_OPTIMIZE, csrA \n");
    CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrB ),
                          "Error after MKL_SPARSE_OPTIMIZE, csrB \n");

/* Compute C = A * B  */
    double time_sum = 0;
    double time_min = 0.0;
    double time;

for(int i = 0;i < REPEAT;i++){
    gettimeofday(&t1, NULL);
    CALL_AND_CHECK_STATUS(mkl_sparse_spmm( SPARSE_OPERATION_NON_TRANSPOSE, csrA, csrB, &csrC ),
                          "Error after MKL_SPARSE_SPMM \n");
    gettimeofday(&t2, NULL);
    time = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
    time_sum += time;
    if(i == 0) time_min = time;
    else time_min = time_min > time ? time : time_min;
}

    CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( csrC, &indexing, &rows, &cols, &pointerB_C, &pointerE_C, &columns_C, &values_C ),
                          "Error after MKL_SPARSE_D_EXPORT_CSR  \n");

    int nnzC = 0;
    for( i = 0; i < m; i++ )
    {
        nnzC += pointerE_C[i] - pointerB_C[i];
    }
    time_sum /= REPEAT;
    //printf("%.8lf,%.8lf\n",time_sum,time_min);
    double GFlops_avg = (nnzC * 2)/(time_sum * 1000000);
    double GFlops_max = (nnzC * 2)/(time_min * 1000000);
    printf("mkl,%s,%d,0,0,%.3lf,%.3lf\n",name,nnz,GFlops_avg,GFlops_max);

/* Analytic Routines for MKL_SPARSE_D_MV.
   HINTS: provides estimate of number and type of upcoming matrix-vector operations
   OPTIMIZE: analyze sparse matrix; choose proper kernels and workload balancing strategy */
    // CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrA, SPARSE_OPERATION_TRANSPOSE,     descr_type_gen, 1 ),
    //                       "Error after MKL_SPARSE_SET_MV_HINT, csrA \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrB, SPARSE_OPERATION_NON_TRANSPOSE, descr_type_gen, 1 ),
    //                       "Error after MKL_SPARSE_SET_MV_HINT, csrB \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_set_mv_hint( csrC, SPARSE_OPERATION_NON_TRANSPOSE, descr_type_gen, 1 ),
    //                       "Error after MKL_SPARSE_SET_MV_HINT, csrC \n");

    // CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrA ),
    //                       "Error after MKL_SPARSE_OPTIMIZE, csrA \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrB ),
    //                       "Error after MKL_SPARSE_OPTIMIZE, csrB \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_optimize( csrC ),
    //                       "Error after MKL_SPARSE_OPTIMIZE, csrC \n");

/* Execution Routines */
/* Step 1:
          Need to compute the following variables:
                 rslt_mv = C * x
                    left = <rslt_mv, y>              */
    // CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrC, descr_type_gen, x, 0.0, rslt_mv ),
    //                       "Error after MKL_SPARSE_D_MV, csrC*x  \n");
    // left = cblas_ddot( m, rslt_mv, 1, y, 1 );

/* Step 2:
          Need to compute the following variables:
           rslt_mv       =     B * x
           rslt_mv_trans = (A)^t * y
                   right = <rslt_mv, rslt_mv_trans>  */

    // CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_NON_TRANSPOSE, 1.0, csrB, descr_type_gen, x, 0.0, rslt_mv ),
    //                       "Error after MKL_SPARSE_D_MV, csrB*x  \n");
    // CALL_AND_CHECK_STATUS(mkl_sparse_d_mv( SPARSE_OPERATION_TRANSPOSE,     1.0, csrA, descr_type_gen, y, 0.0, rslt_mv_trans),
    //                       "Error after MKL_SPARSE_D_MV, csrA*y  \n");
    // right = cblas_ddot( m, rslt_mv, 1, rslt_mv_trans, 1);

/* Step 3:
    //       Compare values obtained for left and right  */
    // residual = fabs(left - right)/(fabs(left)+1);

    // printf( "\n The difference between < C*x , y > and < B*x , (A^t)*y > = %g,\n", residual );
    // printf( " which means that MKL_SPARSE_SPMM arrived correct at a solution.\n" );


/* Deallocate memory */
memory_free:
 //Release matrix handle. Not necessary to deallocate arrays for which we don't allocate memory: values_C, columns_C, pointerB_C, and pointerE_C.
 //These arrays will be deallocated together with csrC structure.
    // if( mkl_sparse_destroy( csrC ) != SPARSE_STATUS_SUCCESS)
    // { printf(" Error after MKL_SPARSE_DESTROY, csrC \n");fflush(0); status = 1; }

 //Deallocate arrays for which we allocate memory ourselves.
    mkl_free(rslt_mv_trans); mkl_free(rslt_mv); mkl_free(x); mkl_free(y);

 //Release matrix handle and deallocate arrays for which we allocate memory ourselves.
    if( mkl_sparse_destroy( csrA ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, csrA \n");fflush(0); status = 1; }
    // mkl_free(values_A); mkl_free(columns_A); mkl_free(rowIndex_A);

    if( mkl_sparse_destroy( csrB ) != SPARSE_STATUS_SUCCESS)
    { printf(" Error after MKL_SPARSE_DESTROY, csrB \n");fflush(0); status = 1; }
    // mkl_free(values_B); mkl_free(columns_B); mkl_free(rowIndex_B);


    return status;
}
