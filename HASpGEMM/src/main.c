#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "common.h"
#include "biio.h"
#include "Test_HASpGEMM.h"


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        printf("Run the code by './thmkl_spgemm matrixA.mtx matrixB.mtx'.\n");
        return 0;
    }

    char *filenameA, *filenameB;
    filenameA = argv[1]; // matrixA
    filenameB = argv[1]; // matrixB

    int mA;
    int nA;
    int nnzA;
    int isSymmetricA;
    MAT_PTR_TYPE *csrRowPtrA_0based;
    int *csrColIdxA_0based;
    MAT_VAL_TYPE *csrValA;
    // get matrix A (support cbd file or mtx file)
    read_Dmatrix_32(&mA, &nA, &nnzA, &csrRowPtrA_0based, &csrColIdxA_0based, &csrValA, &isSymmetricA, filenameA);
    srand(1);


    int mB;
    int nB;
    int nnzB;
    int isSymmetricB;
    MAT_PTR_TYPE *csrRowPtrB_0based;
    int *csrColIdxB_0based;
    MAT_VAL_TYPE *csrValB;

    mB = nA;
    nB = mA;
    nnzB = nnzA;
    csrRowPtrB_0based = (MAT_PTR_TYPE *)malloc(sizeof(MAT_PTR_TYPE) * (mB + 1));
    csrColIdxB_0based = (int *)malloc(sizeof(int) * nnzB);
    csrValB = (MAT_VAL_TYPE *)malloc(sizeof(MAT_VAL_TYPE) * nnzB);
    // Get matrix B 
    matrix_transposition(mA, nA, nnzA,
                         csrRowPtrA_0based, csrColIdxA_0based, csrValA,
                         csrColIdxB_0based, csrRowPtrB_0based, csrValB);

    float *csrValBs = (float *)malloc(sizeof(float) * nnzB);
    for (int i = 0; i < nnzB; i++)
        csrValBs[i] = csrValB[i];

    // Gets the name of the matrix for recording matrix information
    char *tempA = strtok(filenameA, "/");   
    char *mat_nameA;
    char *mat_group_nameA;
    int t = 0;
    int t_for = 4;
    while (tempA)
    {
        tempA = strtok(NULL, "/");
        if (tempA[strlen(tempA) - 4] == '.')
        {
            mat_nameA = strtok(tempA, ".");
            break;
        }
        if(t == (t_for-2))
        {
            mat_group_nameA = tempA;
        }
        t++;
    }
    char mat_name[100] = "";
    char mat_group_name[100] = "";
    strcat(mat_name, mat_nameA);
    strcat(mat_group_name, mat_group_nameA);

    // create matrixA-1based csr
    MAT_PTR_TYPE *csrRowPtrA_1based = (MAT_PTR_TYPE *)malloc((mA + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxA_1based = (int *)malloc(nnzA * sizeof(int));
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < nnzA; i++)
        csrColIdxA_1based[i] = csrColIdxA_0based[i] + 1;
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < mA + 1; i++)
        csrRowPtrA_1based[i] = csrRowPtrA_0based[i] + 1;

    // create matrixB-1based csr
    MAT_PTR_TYPE *csrRowPtrB_1based = (MAT_PTR_TYPE *)malloc((mB + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxB_1based = (int *)malloc(nnzB * sizeof(int));
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < nnzA; i++)
        csrColIdxB_1based[i] = csrColIdxB_0based[i] + 1;
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < mB + 1; i++)
        csrRowPtrB_1based[i] = csrRowPtrB_0based[i] + 1;


    int nnzCub;
    double avg_nnzCub;
    int max_nnzCub, min_nnzCub;
    double avg_time, min_time, avg_gflops, max_gflops, analyze_time;
    int repeat_value;
    double C_row_variation;

    // test 
    Test_HA_d_csrmultcsr(mat_name, nA, mA, nB, mB, csrRowPtrA_1based, csrColIdxA_1based, csrValA, csrRowPtrB_1based, csrColIdxB_1based, csrValB, &nnzCub, &avg_nnzCub, &max_nnzCub, &min_nnzCub, &avg_time, &min_time, &avg_gflops, &max_gflops, &analyze_time, &C_row_variation, &repeat_value);

    //如果需要输出其他参数修改printf函数即可
    printf("HASpGEMM,%s,%d,%d,%.3lf,%.3lf,%.3lf,%.3lf,%.3lf\n",mat_name,nnzA,mA,avg_nnzCub,avg_gflops,max_gflops,avg_time,analyze_time);

    free(csrRowPtrA_0based);
    free(csrColIdxA_0based);
    free(csrValA);

    free(csrRowPtrB_0based);
    free(csrColIdxB_0based);
    free(csrValB);
    free(csrValBs);

    free(csrRowPtrA_1based);
    free(csrColIdxA_1based);
    free(csrRowPtrB_1based);
    free(csrColIdxB_1based);

    return 0;
}
