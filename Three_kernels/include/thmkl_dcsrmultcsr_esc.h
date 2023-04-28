#include <stdlib.h>
#include <string.h>
#include "until.h"

#define setbit(x, y) x |= (1 << y)  // set the yth bit of x is 1
#define clrbit(x, y) x &= ~(1 << y) // set the yth bit of x is 0
#define getbit(x, y) ((x) >> (y)&1) // get the yth bit of x

#define Bound 5000
#define binbd1 64
#define binbd2 4096

void thmkl_dcsrmultcsr_esc(const char *trans, const int *request, const int *sort,
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
            int *jCub = (int *)malloc(rowsize * sizeof(int));
            memset(jCub, 0, rowsize * sizeof(int));

            int incr = 0;
            for (int l = ia[rowid] - ia[0]; l < ia[rowid + 1] - ia[0]; l++)
            {
                int rowB = ja[l] - ia[0];
                for (int k = ib[rowB] - ib[0]; k < ib[rowB + 1] - ib[0]; k++)
                {
                    jCub[incr] = jb[k] - ib[0];
                    incr++;
                }
            }
            quick_sort_key1(&jCub[0], rowsize);
            // merge_sort_key1(&jCub[0], rowsize);

            int nnzr = rowsize > 0 ? 1 : 0;
            for (int idx = 1; idx < rowsize; idx++)
            {
                nnzr = jCub[idx] == jCub[idx - 1] ? nnzr : nnzr + 1;
            }
            ic[pid] = nnzr;
            free(jCub);
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
#pragma omp parallel for
        for (int i = 0; i < *m; i++)
        {
            int rowid = i;
            int rowsize = iCub[rowid];
            if (rowsize == 0)
                continue;
            int *jCub = (int *)malloc(rowsize * sizeof(int));
            double *valCub = (double *)malloc(rowsize * sizeof(double));
            int *d_flagCub = (int *)malloc(rowsize * sizeof(int));
            memset(jCub, 0, rowsize * sizeof(int));
            memset(valCub, 0, rowsize * sizeof(double));
            memset(d_flagCub, 0, rowsize * sizeof(int));

            int incr = 0;
            for (int l = ia[rowid] - ia[0]; l < ia[rowid + 1] - ia[0]; l++)
            {
                int rowB = ja[l] - ia[0];
                double val = a[l];
                for (int k = ib[rowB] - ib[0]; k < ib[rowB + 1] - ib[0]; k++)
                {
                    jCub[incr] = jb[k];
                    valCub[incr] = val * b[k];
                    incr++;
                }
            }

            // sort
            quick_sort_key_val_pair1(jCub, valCub, rowsize);
            // merge_sort_key_val_pair1(jCub, valCub, rowsize);

            // compress
            d_flagCub[0] = 1;
            for (int idx = 0; idx < rowsize - 1; idx++)
            {
                d_flagCub[idx + 1] = jCub[idx + 1] == jCub[idx] ? 0 : 1;
            }
            segmented_sum(valCub, d_flagCub, rowsize);

            int incrn = 0;
            for (int idx = 0; idx < rowsize; idx++)
            {
                if (d_flagCub[idx] == 1)
                {
                    jc[ic[rowid] - ic[0] + incrn] = jCub[idx];
                    c[ic[rowid] - ic[0] + incrn] = valCub[idx];
                    incrn++;
                }
            }

            free(jCub);
            free(valCub);
            free(d_flagCub);
        }
    }
    // free(binPtr);
    // free(binIdx);
    free(iCub);
}