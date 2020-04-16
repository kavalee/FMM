/* Written by Aaron Kavaler aakavalee@gmail.com,
 * algorithm from https://math.berkeley.edu/~strain/128b.S20/fmm1.pdf */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fmm.h"
#include "util.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))




int inline getSCell(CauchyMultiplier* cm, double s, int Q) {
    cm->muls += 1;
    return (int) (s * Q);
}
int inline getTCell(CauchyMultiplier* cm, double t, int Q) {
    cm->muls += 1;
    return (int) (t * Q);

}

double* slowMultiply(CauchyMultiplier* cm) {
    double* f = calloc(cm->n, sizeof(double));
    for (int i = 0; i < cm->n; i++) {
        for (int j = 0; j < cm->n; j++) {
            f[i] += cm->g[j] / ( cm->t[i] - cm->s[j] );
        }
    }
    return f;
}

void computeTMoment(CauchyMultiplier* cm, double** tMoments, int ti, int si, double* tCen, double* sCen, double** sMoments) {
    double tcsI = 1.0 / (tCen[ti] - sCen[si]);
    double powM = 1;
    for (int m = 0; m < cm->p; m++) {
        double powK = tcsI;
        for (int k = 0; k < cm->p; k++) {
            tMoments[ti][m] += (cm->BINOMIAL_CACHE[m + k][k] * powK ) * ( powM * sMoments[si][k] ) ;
            cm->adds += 2;
            cm->muls += 2 + m + k - 1;
            powK *= tcsI;
        }
        powM *= tcsI;
    }
}
void computeDirect(CauchyMultiplier* cm, double* f, int sSize, int* sIndex, int tSize, int* tIndex) {
    cm->adds += tSize * sSize * 2;
    cm->muls += tSize * sSize * 2;
    for (int i = 0; i < tSize; i++) {
        for (int j = 0; j < sSize; j++) {
            f[tIndex[i]] += cm->g[sIndex[j]] / (cm->t[tIndex[i]] - cm->s[sIndex[j]]);
        }
    }
}
int** createReverseIndex(CauchyMultiplier* cm, double* s, int Q,  int* sCellSizes, int sFlag) {
    for (int i = 0; i < cm->n; i++) {
        if (sFlag == 0) {
            sCellSizes[getSCell(cm, s[i], Q)]++;
        } else {
            sCellSizes[getTCell(cm, s[i], Q)]++;
        }
    }
    int** reverseIndex = calloc(Q, sizeof(int*));
    for (int si = 0; si < Q; si++) {
        reverseIndex[si] = calloc(sCellSizes[si], sizeof(int));
    }
    int* tempCellSizes = calloc(Q, sizeof(int));
    for (int i = 0; i < cm->n; i++) {
        int cell = (sFlag == 0) ? getSCell(cm, s[i], Q) : getTCell(cm , s[i], Q);
        reverseIndex[cell][tempCellSizes[cell]] = i;
        tempCellSizes[cell]++;
    }
    free(tempCellSizes);
    return reverseIndex;
}


double* fastMultiply(CauchyMultiplier* cm) {
    if (cm->n < 4) {
        return slowMultiply(cm);
    }

    double *f = calloc(cm->n, sizeof(double));
    cm->sMin = cm->s[0];
    cm->sMax = cm->s[0];
    cm->tMin = cm->t[0];
    cm->tMax = cm->t[0];
    for (int i = 1; i < cm->n; i++) {
        cm->sMin = MIN(cm->sMin, cm->s[i]);
        cm->tMin = MIN(cm->tMin, cm->t[i]);
        cm->sMax = MAX(cm->sMax, cm->s[i]);
        cm->tMax = MAX(cm->tMax, cm->t[i]);
    }
    cm->sH = cm->sMax - cm->sMin;
    cm->tH = cm->tMax - cm->tMin;

    double **sMoments;
    double* prevSCen;
    int L = (int) MAX(log2(cm->n) - 6, 2);
    int QF = 2 << (L - 1);
    for (int l = L; l >= 2; l--) {
        int Q = 2 << (l - 1);

        int *sCellSizes = calloc(Q, sizeof(int));
        int **sReverseIndex = createReverseIndex(cm, cm->s, Q, sCellSizes, 0);
        int *tCellSizes = calloc(Q, sizeof(int));
        int **tReverseIndex = createReverseIndex(cm, cm->t, Q, tCellSizes, 1);

        double *sCen = calloc(Q, sizeof(double));
        cm->adds += Q;
        cm->muls += 3 * Q;
        for (int si = 0; si < Q; si++) {
            sCen[si] = ((double) si / ((double) Q) + 1.0 / (double) (2.0 * Q));
        }
        double *tCen = calloc(Q, sizeof(double));
        cm->adds += Q;
        cm->muls += 3 * Q;
        for (int ti = 0; ti < Q; ti++) {
            tCen[ti] = ((double) ti / ((double) Q) + 1.0 / (double) (2.0 * Q));
        }

        if(L == l) {
            sMoments = calloc(Q, sizeof(double));
            for (int si = 0; si < Q; si++) {
                sMoments[si] = calloc(cm->p, sizeof(double));
            }
            cm->adds += (cm->n + 1) * cm->p;
            cm->muls += cm->n * cm->p;
            for (int i = 0; i < cm->n; i++) {
                int si = getSCell(cm, cm->s[i], Q);
                double smc = cm->s[i] - sCen[si];
                double powK = cm->g[i];
                for (int k = 0; k < cm->p; k++) {
                    sMoments[si][k] += powK;
                    powK *= smc;
                }
            }
        } else {
            double** tempSMoments = calloc(Q, sizeof(double*));
            for (int si = 0; si < Q; si++) {
                double dcA = prevSCen[2 * si] - sCen[si];
                double dcB = prevSCen[(2 * si) + 1] - sCen[si];
                double* powACache = calloc(cm->p, sizeof(double));
                double* powBCache = calloc(cm->p, sizeof(double));
                powACache[0] = 1;
                powBCache[0] = 1;
                cm->muls += 2 * cm->p;
                for (int i = 1; i < cm->p; i++) {
                    powACache[i] = powACache[i - 1] * dcA;
                    powBCache[i] = powBCache[i - 1] * dcB;
                }
                tempSMoments[si] = calloc(cm->p, sizeof(double));
                for (int k = 0; k < cm->p; k++) {
                    double sum = 0;
                    for (int m = 0; m <= k; m++) {
                        cm->adds += 2;
                        cm->muls += 3;
                        sum += cm->BINOMIAL_CACHE[k][m] * (  powACache[k - m] * sMoments[2 * si][m]
                                                           + powBCache[k - m] * sMoments[(2 * si) + 1][m]);
                    }
                    tempSMoments[si][k] = sum;
                }
                free(powACache);
                free(powBCache);
            }

            for (int i = 0; i < (2 * Q); i++) {
                free(sMoments[i]);
            }
            free(sMoments);
            sMoments = calloc(Q, sizeof(double*));
            for (int i = 0; i < Q; i++) {
                sMoments[i] = calloc(cm->p, sizeof(double));
                memcpy(sMoments[i], tempSMoments[i], cm->p * sizeof(double));
                free(tempSMoments[i]);
            }
            free(tempSMoments);
        }

        double **tMoments = calloc(Q, sizeof(double *));
        for (int ti = 0; ti < Q; ti++) {
            tMoments[ti] = calloc(cm->p, sizeof(double));
        }
        for (int si = 0; si < Q - 2; si++) {
            computeTMoment(cm, tMoments, si + 2, si, tCen, sCen, sMoments);
            computeTMoment(cm, tMoments, si, si + 2, tCen, sCen, sMoments);
        }
        for (int si = 0; si < Q - 2; si += 2) {
            computeTMoment(cm, tMoments, si + 3, si, tCen, sCen, sMoments);
            computeTMoment(cm, tMoments, si, si + 3, tCen, sCen, sMoments);
        }

        if (L == l) {
            for (int si = 0; si < Q; si++) {
                computeDirect(cm, f, sCellSizes[si], sReverseIndex[si], tCellSizes[si], tReverseIndex[si]);
            }
            for (int si = 0; si < Q - 1; si++) {
                computeDirect(cm, f, sCellSizes[si + 1], sReverseIndex[si + 1], tCellSizes[si], tReverseIndex[si]);
                computeDirect(cm, f, sCellSizes[si], sReverseIndex[si], tCellSizes[si + 1], tReverseIndex[si + 1]);
            }
        }

        cm->adds += cm->n * cm->p;
        cm->muls += 2 * cm->n * cm->p;
        for (int i = 0; i < cm->n; i++) {
            double tPow = 1;
            int ti = getTCell(cm, cm->t[i], Q);
            double cmt = tCen[ti] - cm->t[i];
            for (int m = 0; m < cm->p; m++) {
                f[i] += tMoments[ti][m] * tPow;
                tPow *= cmt;
            }
        }

        if (L != l){
            free(prevSCen);
        }
        prevSCen = calloc(Q, sizeof(double));
        memcpy(prevSCen, sCen, sizeof(double) * Q);

        for (int i = 0; i < Q; i++) {
            free(tMoments[i]);
            free(sReverseIndex[i]);
            free(tReverseIndex[i]);
        }

        free(tMoments);
        free(sCellSizes);
        free(sReverseIndex);
        free(tCellSizes);
        free(tReverseIndex);
        free(sCen);
        free(tCen);
    }
    for (int i = 0; i < 4; i++) {
        free(sMoments[i]);
    }
    free(sMoments);
    free(prevSCen);
    return f;
}
void initializeBinomialCache(CauchyMultiplier* cm, int N) {
    cm->BINOMIAL_CACHE = calloc(N, sizeof(double*));
    cm->BINOMIAL_CACHE[0] = calloc(N, sizeof(double*));
    cm->BINOMIAL_CACHE[0][0] = 1;
    for (int r = 1; r < N; r++) {
        cm->BINOMIAL_CACHE[r] = calloc(N, sizeof(double));
        cm->BINOMIAL_CACHE[r][0] = 1;
        for (int c = 1; c < N; c++) {
            cm->BINOMIAL_CACHE[r][c] = cm->BINOMIAL_CACHE[r - 1][c - 1] + cm->BINOMIAL_CACHE[r - 1][c];
        }
    }
}

CauchyMultiplier* newCauchyMultiplier(double* s, double* t, double* g, int n, int p) {
    CauchyMultiplier* cm = calloc(1, sizeof(CauchyMultiplier));
    cm->s = s;
    cm->t = t;
    cm->g = g;
    cm->n = n;
    cm->p = p;
    cm->adds = 0;
    cm->muls = 0;
    initializeBinomialCache(cm, 2 * p);
    return cm;
}
void freeCauchyMultiplier(CauchyMultiplier* cm) {
    free(cm->s);
    free(cm->t);
    free(cm->g);
    for (int i = 0; i < 2 * cm->p; i++) {
        free(cm->BINOMIAL_CACHE[i]);
    }
    free(cm->BINOMIAL_CACHE);
    free(cm);
}


