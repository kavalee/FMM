#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))


typedef struct CauchyMultiplier{
    double* s;
    double* t;
    double* g;
    int n;
    int p;
    double sMin;
    double sMax;
    double tMin;
    double tMax;
    double sH;
    double tH;
    double** BINOMIAL_CACHE;
    int adds;
    int muls;


} CauchyMultiplier;

int inline getSCell(CauchyMultiplier* cm, double s, int Q) {
    cm->adds++;
    cm->muls += 2;
    double x = (s - cm->sMin)/cm->sH;
    return MIN( (int) (x * Q), (Q - 1) );
}
int inline getTCell(CauchyMultiplier* cm, double t, int Q) {
    cm->adds++;
    cm->muls += 2;
    double x = (t - cm->tMin)/cm->tH;
    return MIN( (int) (x * Q), (Q - 1) );

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


void printMatrix(double** m, int rows, int cols) {
    printf("\n");
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < cols; c++) {
            printf("%.4e ", m[r][c]);
        }
        printf("\n");
    }
}
void computeTMoment(CauchyMultiplier* cm, double** tMoments, int ti, int si, double* cen, double** sMoments) {
    for (int m = 0; m < cm->p; m++) {
        for (int k = 0; k < cm->p; k++) {
            tMoments[ti][m] += cm->BINOMIAL_CACHE[m + k][k] * pow(cen[ti] - cen[si], -m - k - 1) * sMoments[si][k];
            cm->adds += 2;
            cm->muls += 2 + m + k - 1;
        }
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
int** createReverseIndex(CauchyMultiplier* cm, double* s, int Q,  int* sCellSizes) {
    for (int i = 0; i < cm->n; i++) {
        sCellSizes[getSCell(cm, s[i], Q)]++;
    }
    int** reverseIndex = calloc(Q, sizeof(int*));
    for (int si = 0; si < Q; si++) {
        reverseIndex[si] = calloc(sCellSizes[si], sizeof(int));
    }
    int* tempCellSizes = calloc(Q, sizeof(int));
    for (int i = 0; i < cm->n; i++) {
        int cell = getSCell(cm, s[i], Q);
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
    printf("%lf\n", cm->sH);


    int L = (int) log2(cm->n);
    for (int l = 2; l <= L; l++) {
        int Q = 2 << (l - 1);

        int *sCellSizes = calloc(Q, sizeof(int));
        int **sReverseIndex = createReverseIndex(cm, cm->s, Q, sCellSizes);
        int *tCellSizes = calloc(Q, sizeof(int));
        int **tReverseIndex = createReverseIndex(cm, cm->t, Q, tCellSizes);

        double *cen = calloc(Q, sizeof(double));
        cm->adds += Q;
        cm->muls += 3 * Q;
        for (int si = 0; si < Q; si++) {
            cen[si] = (double) si / ((double) Q) + 1.0 / (double) (2.0 * Q);
        }

        double **sPowers = calloc(cm->n, sizeof(double *));
        for (int i = 0; i < cm->n; i++) {
            sPowers[i] = calloc(cm->p, sizeof(double));
            sPowers[i][0] = cm->g[i];
            double smc = cm->s[i] - cen[getSCell(cm, cm->s[i], Q)];
            cm->adds++;
            cm->muls += cm->p;
            for (int k = 1; k < cm->p; k++) {
                sPowers[i][k] = sPowers[i][k - 1] * smc;
            }
        }
        double **sMoments = calloc(Q, sizeof(double *));
        cm->adds += cm->n * cm->p;
        for (int si = 0; si < Q; si++) {
            sMoments[si] = calloc(cm->p, sizeof(double));
            for (int i = 0; i < sCellSizes[si]; i++) {
                for (int k = 0; k < cm->p; k++) {
                    sMoments[si][k] += sPowers[sReverseIndex[si][i]][k];
                }
            }
        }
        for (int i = 0; i < cm->n; i++) {
            free(sPowers[i]);
        }
        free(sPowers);

        double **tMoments = calloc(Q, sizeof(double *));
        for (int ti = 0; ti < Q; ti++) {
            tMoments[ti] = calloc(cm->p, sizeof(double));
        }


        for (int si = 0; si < Q - 2; si++) {
            computeTMoment(cm, tMoments, si + 2, si, cen, sMoments);
            computeTMoment(cm, tMoments, si, si + 2, cen, sMoments);
        }
        for (int si = 0; si < Q - 2; si += 2) {
            computeTMoment(cm, tMoments, si + 3, si, cen, sMoments);
            computeTMoment(cm, tMoments, si, si + 3, cen, sMoments);
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
            double cmt = cen[ti] - cm->t[i];
            for (int m = 0; m < cm->p; m++) {
                f[i] += tMoments[ti][m] * tPow;
                tPow *= cmt;
            }
        }
        for (int i = 0; i < Q; i++) {
            free(sMoments[i]);
            free(tMoments[i]);
            free(sReverseIndex[i]);
            free(tReverseIndex[i]);
        }
        free(sMoments);
        free(tMoments);
        free(sCellSizes);
        free(sReverseIndex);
        free(tCellSizes);
        free(tReverseIndex);
        free(cen);
    }
    return f;
}
void updatePowerCache() {}
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
}

void printVector(double* vec, int len) {
    printf("[");
    for (int i = 0; i < len; i++) {
        printf("%.4e ," , vec[i]);
    }
    printf("]\n");
}

double* generateRandomNVector(int N) {
    double* f = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        f[i] = (double) rand() / RAND_MAX;
    }
    return f;
}
void testMultiplyAccuracy(int N) {
    double* s = generateRandomNVector(N);
    double* t = generateRandomNVector(N);
    double* g = generateRandomNVector(N);
    CauchyMultiplier* cm = newCauchyMultiplier(s, t, g, N, 10);
    double* fast = fastMultiply(cm);
    double* slow = slowMultiply(cm);
    printVector(fast, N);
    printVector(slow, N);
    freeCauchyMultiplier(cm);
}
void testMultiplySpeed(int N) {
    double* s = generateRandomNVector(N);
    double* t = generateRandomNVector(N);
    double* g = generateRandomNVector(N);


    CauchyMultiplier* cm = newCauchyMultiplier(s, t, g, N, 10);
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double* fast = fastMultiply(cm);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    long delta = (end.tv_sec - start.tv_sec) * pow(10,9)  +  (end.tv_nsec - start.tv_nsec);
    double sec = ((double) delta ) / pow(10,9);
    printf("Took %lf,  \n", sec);
    printf("Adds:%d \nMuls:%d \n", cm->adds, cm->muls);
    //double* slow = slowMultiply(cm);
    //printVector(fast, N);
    //printVector(slow, N);
    freeCauchyMultiplier(cm);
}

int main() {
    srand(time(NULL));
    double s[] = {0.0450129590219045, 0.14058320819439485, 0.1523630580433034,
                  0.3334969196167895, 0.7051533278481288, 0.7657209605120057, 0.9009994392272345, 0.9335046628841073};
    double t[] = {0.06572442852925431, 0.3853804136623742, 0.4612045175524898,
                  0.6612689052573647, 0.6922738788833437, 0.7465256805234415, 0.8149811433845008, 0.8530371620489543};
    double g[] = {0.05297700176067044, 0.18096566993453178, 0.23747267040559616,
                  0.568132463298329, 0.6383567455540913, 0.7756066277976567, 0.8870273842314716, 0.9643387582918744};
    int n = 8;
    int p = 10;
    CauchyMultiplier* cm = newCauchyMultiplier(s, t, g, n, p);
    //testMultiply(8);
    testMultiplySpeed(100000);
    printVector(fastMultiply(cm), cm->n);
    return 0;
}
