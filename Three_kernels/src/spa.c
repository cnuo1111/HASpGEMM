#include <time.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>

#include "common.h"
#include "biio.h"
#include "Test_d_csrmultcsr.h"

int main(int argc, char **argv)
{

    if (argc < 2)
    {
        printf("Run the code by './thmkl_spgemm matrixA.mtx matrixB.mtx option'.\n");
        return 0;
    }

    char *filenameA, *filenameB;
    filenameA = argv[1]; 
    filenameB = argv[1]; 
    int option;
    option = 1; 

    int mA;
    int nA;
    int nnzA;
    int isSymmetricA;
    MAT_PTR_TYPE *csrRowPtrA_0based;
    int *csrColIdxA_0based;
    MAT_VAL_TYPE *csrValA;
    read_Dmatrix_32(&mA, &nA, &nnzA, &csrRowPtrA_0based, &csrColIdxA_0based, &csrValA, &isSymmetricA, filenameA);
    srand(1);

    for (int i = 0; i < nnzA; i++)
        csrValA[i] = rand() % 10 + 1;
    float *csrValAs = (float *)malloc(sizeof(float) * nnzA);
    for (int i = 0; i < nnzA; i++)
        csrValAs[i] = csrValA[i];
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
    matrix_transposition(mA, nA, nnzA,
                         csrRowPtrA_0based, csrColIdxA_0based, csrValA,
                         csrColIdxB_0based, csrRowPtrB_0based, csrValB);


    float *csrValBs = (float *)malloc(sizeof(float) * nnzB);
    for (int i = 0; i < nnzB; i++)
        csrValBs[i] = csrValB[i];

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

    double avg_row_A = 0;
    int min_row_A = csrRowPtrA_0based[1] - csrRowPtrA_0based[0], max_row_A = 0;
    double avg_row_B = 0;

    int min_row_B, max_row_B;
    double avg_col_B = 0;
    int min_col_B = csrRowPtrB_0based[1] - csrRowPtrB_0based[0], max_col_B = 0;
    double row, col;
    double *rowA = (double*)malloc(mA * sizeof(double));
    double *colB = (double*)malloc(mB * sizeof(double));
    for (MAT_PTR_TYPE i = 0; i < mA; i++)
    {
        row = csrRowPtrA_0based[i + 1] - csrRowPtrA_0based[i];
        rowA[i] = row;
        if (row < min_row_A)
            min_row_B = min_row_A = row;
        if (row > max_row_A)
            max_row_B = max_row_A = row;
        avg_row_A += row;
    }
    avg_row_A = avg_row_A / (double)mA;
    avg_row_B = avg_row_A;
    for (MAT_PTR_TYPE i = 0; i < mB; i++)
    {
        col = csrRowPtrB_0based[i + 1] - csrRowPtrB_0based[i];
        colB[i] = col;
        if (col < min_col_B)
            min_col_B = col;
        if (col > max_col_B)
            max_col_B = col;
        avg_col_B += col;
    }
    avg_col_B = avg_col_B / (double)mB;

    double A_row_variation = matrix_variance(rowA, mA);
    double B1_row_variation = A_row_variation;
    double B_row_variation = matrix_variance(colB, mB);
    double C_row_variation;
    free(rowA);
    free(colB);

    // create matrixA-1based csr
    MAT_PTR_TYPE *csrRowPtrA_1based = (MAT_PTR_TYPE *)malloc((mA + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxA_1based = (int *)malloc(nnzA * sizeof(int));
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < nnzA; i++)
        csrColIdxA_1based[i] = csrColIdxA_0based[i] + 1;

    for (MAT_PTR_TYPE i = 0; i < mA + 1; i++)
        csrRowPtrA_1based[i] = csrRowPtrA_0based[i] + 1;

    // create matrixB-1based csr
    MAT_PTR_TYPE *csrRowPtrB_1based = (MAT_PTR_TYPE *)malloc((mB + 1) * sizeof(MAT_PTR_TYPE));
    int *csrColIdxB_1based = (int *)malloc(nnzB * sizeof(int));
#pragma omp parallel for
    for (MAT_PTR_TYPE i = 0; i < nnzA; i++)
        csrColIdxB_1based[i] = csrColIdxB_0based[i] + 1;

    for (MAT_PTR_TYPE i = 0; i < mB + 1; i++)
        csrRowPtrB_1based[i] = csrRowPtrB_0based[i] + 1;

    int nnzCub;
    double avg_nnzCub;
    int max_nnzCub, min_nnzCub;
    double avg_time_esc, min_time_esc, avg_gflops_esc, max_gflops_esc;
    double avg_time_hash, min_time_hash, avg_gflops_hash, max_gflops_hash;
    double avg_time_spa, min_time_spa, avg_gflops_spa, max_gflops_spa;
    double sum_time;
    int repeat_value;

    Test_d_csrmultcsr(option, mat_name, nA, mA, nB, mB, csrRowPtrA_1based, csrColIdxA_1based, csrValA, csrRowPtrB_1based, csrColIdxB_1based, csrValB, 3, &nnzCub, &avg_nnzCub, &max_nnzCub, &min_nnzCub, &avg_time_spa, &min_time_spa, &avg_gflops_spa, &max_gflops_spa, &C_row_variation, &repeat_value);

    printf("SPA,%s,%d,%d,%.3lf,%.3lf,%.3lf\n",mat_name,nnzA,mA,avg_nnzCub,avg_gflops_spa,max_gflops_spa);

    free(csrRowPtrA_0based);
    free(csrColIdxA_0based);
    free(csrValA);
    free(csrValAs);

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
