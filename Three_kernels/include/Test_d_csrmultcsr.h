#ifndef _CSRMULTCSR_
#define _CSRMULTCSR_
#include "thmkl_dcsrmultcsr_esc.h"
#include "thmkl_dcsrmultcsr_hash.h"
#include "thmkl_dcsrmultcsr_spa.h"
void Test_d_csrmultcsr(int option,
                       char *mat_name,
                       int nA,
                       int mA,
                       int nB,
                       int mB,
                       int *csrRowPtrA_1based,
                       int *csrColIdxA_1based,
                       MAT_VAL_TYPE *csrValA,
                       int *csrRowPtrB_1based,
                       int *csrColIdxB_1based,
                       MAT_VAL_TYPE *csrValB,
                       int method,
                       int *nnzCub, double *avg_nnzCub, int *max_nnzCub, int *min_nnzCub,
                       MAT_VAL_TYPE *avg_time, MAT_VAL_TYPE *min_time,
                       MAT_VAL_TYPE *avg_gflops, MAT_VAL_TYPE *max_gflops,
                       MAT_VAL_TYPE *C_row_variation, int *repeat_value)
{
    if (option == 1 || option == 0)
    {
        struct timeval t1, t2;
        MAT_PTR_TYPE nnzc = 0;
        MAT_PTR_TYPE *csrRowPtrC_1based = (MAT_PTR_TYPE *)malloc((mA + 1) * sizeof(MAT_PTR_TYPE));
        int *csrColIdxC_1based;
        MAT_VAL_TYPE *csrValC;

        MAT_PTR_TYPE nnzA = csrRowPtrA_1based[mA] - 1;
        MAT_PTR_TYPE nnzB = csrRowPtrB_1based[mB] - 1;

        if(nnzA < 100000)
            *repeat_value = 10;
        else if(nnzA < 1000000)
            *repeat_value = 3;
        else 
            *repeat_value = 1;

        // printf("%d %d\n",*repeat_value,method);

        char transa = 'N';
        int request;
        int sort = 3;
        int info;
        double sum_time = 0;
        *avg_time = 0, *min_time;
        double max_time = 0;

        // omp_set_num_threads(THREAD_MAX_TEST);
        for (int i = 0; i < *repeat_value; i++)
        {
            memset(csrRowPtrC_1based, 0, (mA + 1) * sizeof(MAT_PTR_TYPE));

            request = 1;
            gettimeofday(&t1, NULL);
            // 第一步获取结构
            if (method == 1)//esc
            {
                thmkl_dcsrmultcsr_esc(&transa, &request, &sort, &mA, &nA, &nB,
                                      (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                                      (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                                      (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                                      &nnzc, &info);
            }
            else if (method == 2)//hash
            {
                thmkl_dcsrmultcsr_hash(&transa, &request, &sort, &mA, &nA, &nB,
                                       (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                                       (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                                       (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                                       &nnzc, &info);
            }
            else if (method == 3)//spa
            {
                thmkl_dcsrmultcsr_spa(&transa, &request, &sort, &mA, &nA, &nB,
                                      (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                                      (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                                      (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                                      &nnzc, &info);
            }
            gettimeofday(&t2, NULL);
            MAT_VAL_TYPE time1_dcsrmultcsr = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;

            // 给数组开空间
            nnzc = csrRowPtrC_1based[mA] - 1;
            csrColIdxC_1based = (int *)malloc((nnzc) * sizeof(int));
            csrValC = (MAT_VAL_TYPE *)malloc((nnzc) * sizeof(MAT_VAL_TYPE));
            memset(csrColIdxC_1based, 0, (nnzc) * sizeof(int));
            memset(csrValC, 0, (nnzc) * sizeof(MAT_VAL_TYPE));
            request = 2;

            gettimeofday(&t1, NULL);
            // 第二步进行求解
            if (method == 1)//esc
            {
                thmkl_dcsrmultcsr_esc(&transa, &request, &sort, &mA, &nA, &nB,
                                      (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                                      (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                                      (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                                      &nnzc, &info);
            }
            else if (method == 2)//hash
            {
                thmkl_dcsrmultcsr_hash(&transa, &request, &sort, &mA, &nA, &nB,
                                       (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                                       (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                                       (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                                       &nnzc, &info);
            }
            else if (method == 3)//spa
            {
                thmkl_dcsrmultcsr_spa(&transa, &request, &sort, &mA, &nA, &nB,
                                      (MAT_VAL_TYPE *)csrValA, csrColIdxA_1based, csrRowPtrA_1based,
                                      (MAT_VAL_TYPE *)csrValB, csrColIdxB_1based, csrRowPtrB_1based,
                                      (MAT_VAL_TYPE *)csrValC, csrColIdxC_1based, csrRowPtrC_1based,
                                      &nnzc, &info);
            }

            gettimeofday(&t2, NULL);
            MAT_VAL_TYPE time2_dcsrmultcsr = (t2.tv_sec - t1.tv_sec) * 1000.0 + (t2.tv_usec - t1.tv_usec) / 1000.0;
            MAT_VAL_TYPE time = time1_dcsrmultcsr + time2_dcsrmultcsr;
            if (i == 0) // 获取最短时间
                *min_time = time;
            sum_time += time;
            *min_time = time < *min_time ? time : *min_time;
        }
        // printf("time:%lf\n",sum_time);

        // 计算C总计算量（nnzCub）、C平均行计算量、C最大行计算量、C最小行计算量、C行计算量variation
        double *rowC = (double*)malloc(mA * sizeof(double));
        *nnzCub = 0;
        int cub = 0, cubnum = 0;
        *min_nnzCub = (csrRowPtrB_1based[csrColIdxA_1based[0]] - csrRowPtrB_1based[csrColIdxA_1based[0] - 1]) * nnzA;
        *max_nnzCub = 0;
        int m = 0;
        int i = 0;
        int j = 0;
        while (i < nnzA)
        {
            cubnum = 0;
            while (i < csrRowPtrA_1based[m + 1] - 1)
            {
                int rowB = csrColIdxA_1based[i] - 1;
                cub = csrRowPtrB_1based[rowB + 1] - csrRowPtrB_1based[rowB];
                cubnum += cub;
                i++;
            }
            *nnzCub += cubnum;
            rowC[j] = (double)cubnum;
            if (cubnum > *max_nnzCub)
                *max_nnzCub = cubnum;
            if (cubnum < *min_nnzCub)
                *min_nnzCub = cubnum;
            j++;
            m++;
        }
        *avg_nnzCub = (double)*nnzCub / (double)mA;
        *C_row_variation = matrix_variance(rowC, mA);

        // 计算平均时间
        *avg_time = sum_time / (double)*repeat_value;
        // 计算平均性能
        *avg_gflops = (2 * (*nnzCub)) / ((*avg_time) * 1000000);
        // 计算最高性能
        *max_gflops = ((*nnzCub) * 2) / ((*min_time) * 1000000);
    }
}
#endif
