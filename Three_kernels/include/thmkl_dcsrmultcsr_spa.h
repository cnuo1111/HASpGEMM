#include <stdlib.h>
#include <string.h>
#include "until.h"

#define setbit(x, y) x |= (1 << y)  // set the yth bit of x is 1
#define clrbit(x, y) x &= ~(1 << y) // set the yth bit of x is 0
#define getbit(x, y) ((x) >> (y)&1) // get the yth bit of x

#define Bound 5000
#define binbd1 64
#define binbd2 4096

void thmkl_dcsrmultcsr_spa(const char *trans, const int *request, const int *sort,
                           const int *m, const int *n, const int *k,
                           double *a, int *ja, int *ia,
                           double *b, int *jb, int *ib,
                           double *c, int *jc, int *ic,
                           const int *nzmax, int *info)
{
    int *iCub = (int *)malloc(*m * sizeof(int));
    memset(iCub, 0, sizeof(int) * *m);

    for (int i = 0; i < *m; i++)
    {
        for (int j = ia[i] - ia[0]; j < ia[i + 1] - ia[0]; j++)
        {
            int rowAtop = ja[j] - ia[0];
            iCub[i] += ib[rowAtop + 1] - ib[rowAtop];
        }
        // if(i==0) printf("nnztop_row0 = %d\n", iCub[i]);
    }

    // printf("ASD\n");
    // step 1: compute and get the number of nnzC
    if (*request == 1)
    {
#pragma omp parallel for
        for (int i = 0; i < *m; i++)
        {
            int pid = 0;
            if (*m > 100000)
                pid = i + 1;
            else
                pid = i;
            int rowid = i;
            int rowsize = iCub[rowid];
            if (rowsize == 0)
                continue;
            int len = (*k + 31) / 32;
            unsigned int *mask = (unsigned int *)malloc(sizeof(unsigned int) * len);
            memset(mask, 0, sizeof(unsigned int) * len);
            for (int offsetA = ia[rowid] - ia[0]; offsetA < ia[rowid + 1] - ia[0]; offsetA++)
            {
                int col = ja[offsetA] - ia[0];
                for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                {
                    const int key = jb[l] - ib[0];
                    setbit(mask[key / 32], key % 32);
                }
            }
            int nnzr = 0;
            for (int cid = 0; cid < *k; cid++)
            {
                if (getbit(mask[cid / 32], cid % 32) == 1)
                {
                    nnzr++;
                }
            }
            ic[pid] = nnzr;
            free(mask);
        }
        if (*m > 100000)
            scan_par(ic, *m + 1);
        else
            exclusive_scan(ic, *m + 1);
        for (int i = 0; i < *m + 1; i++)
        {
            ic[i] += ia[0];
        }
    }

    // step2: compute the value of C
    else
    {
        // printf("2\n");
#pragma omp parallel for
        for (int i = 0; i < *m; i++)
        {
            int rowid = i;
            int rowsize = iCub[rowid];
            if (rowsize == 0)
                continue;
            if (rowsize == 0)
                continue;
            int len = (*k + 31) / 32;
            unsigned int *mask = (unsigned int *)malloc(sizeof(unsigned int) * len);
            memset(mask, 0, sizeof(unsigned int) * len);
            double *d_dense_row_value = (double *)malloc((*k) * sizeof(double));
            memset(d_dense_row_value, 0, (*k) * sizeof(double));

            for (int offsetA = ia[rowid] - ia[0]; offsetA < ia[rowid + 1] - ia[0]; offsetA++)
            {
                int col = ja[offsetA] - ia[0];
                for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                {
                    const int key = jb[l] - ib[0];
                    setbit(mask[key / 32], key % 32);
                    d_dense_row_value[key] += b[l] * a[offsetA];
                }
            }

            int nnzr = ic[rowid] - ic[0];
            for (int cid = 0; cid < *k; cid++)
            {
                if (getbit(mask[cid / 32], cid % 32) == 1)
                {
                    c[nnzr] = d_dense_row_value[cid];
                    jc[nnzr] = cid + ic[0];
                    nnzr++;
                }
            }

            free(mask);
            free(d_dense_row_value);
        }
    }
    // free(binPtr);
    // free(binIdx);
    free(iCub);
}