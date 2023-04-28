#include <stdlib.h>
#include <string.h>
#include "until.h"

#define setbit(x, y) x |= (1 << y)  // set the yth bit of x is 1
#define clrbit(x, y) x &= ~(1 << y) // set the yth bit of x is 0
#define getbit(x, y) ((x) >> (y)&1) // get the yth bit of x

#define Bound 5000
#define binbd1 64
#define binbd2 4096

void thmkl_dcsrmultcsr_hash(const char *trans, const int *request, const int *sort,
                            const int *m, const int *n, const int *k,
                            double *a, int *ja, int *ia,
                            double *b, int *jb, int *ib,
                            double *c, int *jc, int *ic,
                            const int *nzmax, int *info)
{
    // if (*m != *n)
    // {
    //     printf("Cannot multiply matrix A of size %i x %i and matrix B of size %i x %i, return.\n",
    //            *m, *k, mB, *n);
    //     return;
    // }
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
            int hashsize_full_reg = rowsize / 0.75;
            int *tmpHashtable = (int *)malloc(hashsize_full_reg * sizeof(int));
            memset(tmpHashtable, -1, sizeof(int) * hashsize_full_reg);
            for (int blkj = ia[rowid] - ia[0]; blkj < ia[rowid + 1] - ia[0]; blkj++)
            {
                int col = ja[blkj] - ia[0];
                for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                {
                    const int key = jb[l] - ib[0];
                    int hashadr = (key) % hashsize_full_reg;
                    while (1)
                    {
                        const int keyexist = tmpHashtable[hashadr];
                        if (keyexist == key)
                        {
                            break;
                        }
                        else if (keyexist == -1)
                        {
                            tmpHashtable[hashadr] = key;
                            ic[pid]++;
                            break;
                        }
                        else
                        {
                            hashadr = (hashadr + 1) % hashsize_full_reg;
                        }
                    }
                }
            }
            free(tmpHashtable);
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
            int hashsize_full_reg = (ic[rowid + 1] - ic[rowid]) / 0.75;
            int *tmpHashtable = (int *)malloc(hashsize_full_reg * sizeof(int));
            memset(tmpHashtable, -1, sizeof(int) * hashsize_full_reg);
            double *tmpValue = (double *)malloc(hashsize_full_reg * sizeof(double));
            memset(tmpValue, 0, sizeof(double) * hashsize_full_reg);
            for (int blkj = ia[rowid] - ia[0]; blkj < ia[rowid + 1] - ia[0]; blkj++)
            {
                int col = ja[blkj] - ia[0];
                for (int l = ib[col] - ib[0]; l < ib[col + 1] - ib[0]; l++)
                {
                    const int key = jb[l] - ib[0];
                    int hashadr = (key * 107) % hashsize_full_reg;
                    // int hashadr = key % hashsize_full_reg;
                    while (1)
                    {
                        const int keyexist = tmpHashtable[hashadr];
                        if (keyexist == key)
                        {
                            tmpValue[hashadr] += b[l] * a[blkj];
                            break;
                        }
                        else if (keyexist == -1)
                        {
                            tmpHashtable[hashadr] = key;
                            tmpValue[hashadr] = b[l] * a[blkj];
                            break;
                        }
                        else
                        {
                            hashadr = (hashadr + 1) % hashsize_full_reg;
                        }
                    }
                }
            }
            int cptr = ic[rowid] - ic[0];
            for (int k = 0; k < hashsize_full_reg; k++)
            {
                if (tmpHashtable[k] != -1)
                {
                    jc[cptr] = tmpHashtable[k] + ib[0];
                    c[cptr] = tmpValue[k];
                    cptr++;
                }
            }
            free(tmpHashtable);
            free(tmpValue);
            int nnzcnt = ic[rowid + 1] - ic[rowid];
            quick_sort_key_val_pair1(jc + ic[rowid] - ic[0], c + ic[rowid] - ic[0], nnzcnt);
            // merge_sort_key_val_pair1(jc + ic[rowid] - ic[0], c + ic[rowid] - ic[0], nnzcnt);
        }
    }
    // free(binPtr);
    // free(binIdx);
    free(iCub);
}