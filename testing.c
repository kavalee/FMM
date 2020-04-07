//
// Created by kavalee on 4/6/20.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "util.h"
#include "fmm.h"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

int PRECISIONS[] = {10,20,30};
#define NUM_PRECISIONS 3

double* generateRandomNVector(int N, double size) {
    double* f = calloc(N, sizeof(double));
    for (int i = 0; i < N; i++) {
        f[i] = size * (double) rand() / RAND_MAX;
    }
    return f;
}
void testMultiplyAccuracy(int N) {
    double* s = generateRandomNVector(N, 1);
    double* t = generateRandomNVector(N, 1);
    double* g = generateRandomNVector(N, 1);
    CauchyMultiplier* cm = newCauchyMultiplier(s, t, g, N, 10);
    double* fast = fastMultiply(cm);
    double* slow = slowMultiply(cm);
    printVector(fast, N);
    printVector(slow, N);
    freeCauchyMultiplier(cm);
}
void testMultiplySpeed(int N) {
    double* s = generateRandomNVector(N, 1);
    double* t = generateRandomNVector(N, 1);
    double* g = generateRandomNVector(N, 1);


    FILE* file = fopen("output/test.txt", "w");
    fprintf(file, "10");

    CauchyMultiplier* cm = newCauchyMultiplier(s, t, g, N, 10);
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    double* fast = fastMultiply(cm);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    long delta = (end.tv_sec - start.tv_sec) * pow(10,9)  +  (end.tv_nsec - start.tv_nsec);
    double sec = ((double) delta ) / pow(10,9);
    printf("Took %lf,  \n", sec);
    printf("Adds:%u \nMuls:%u \n", cm->adds, cm->muls);
    //double* slow = slowMultiply(cm);
    //printVector(fast, N);
    //printVector(slow, N);
    freeCauchyMultiplier(cm);
    fclose(file);
}
double maxErrorAbs(double *fast, double *slow, int n) {
    double max = 0;
    //printf("%d\n",n );
    for (int i = 0; i < n; i++) {
        if (slow[i] != 0) {
            double errAbs = fabs(fast[i] - slow[i]);
            //printf("%lf, ", errAbs);
            max = MAX(errAbs, max);
        }
    }
    //printf("\n");
    return max;
}
void generateErrorData(int N, int num_trials) {
    char buffer[200];
    sprintf(buffer, "output/speed-N%d-trials%d.csv", N, num_trials);
    FILE *speed = fopen(buffer, "w");
    sprintf(buffer, "output/error-N%d-trials%d.csv", N, num_trials);
    FILE *error = fopen(buffer, "w");
    for (int n = 0; n < N; n++) {
        printf("%d\n", n);
        double *max = calloc(NUM_PRECISIONS, sizeof(double));
        for (int t = 0; t < num_trials; t++) {
            double *s = generateRandomNVector(n, 1);
            double *t = generateRandomNVector(n, 1);
            double *g = generateRandomNVector(n, 1);

            for (int p = 0; p < NUM_PRECISIONS; p++) {
                CauchyMultiplier *cm = newCauchyMultiplier(s, t, g, n, PRECISIONS[p]);
                double *slow = slowMultiply(cm);
                double *fast = fastMultiply(cm);
                max[p] += maxErrorAbs(slow, fast, n);
                free(slow);
                free(fast);
                free(cm);
            }
            free(s);
            free(t);
            free(g);
        }
        fprintf(error, "%d, ", n);
        for (int i = 0; i < NUM_PRECISIONS; i++) {
            fprintf(error, "%.4e", max[i] / num_trials);
            if (i != NUM_PRECISIONS - 1) {
                fprintf(error, ", ");
            }
        }
        fprintf(error, "\n");
        free(max);

    }
    fclose(speed);
    fclose(error);
}
void generateFlopsData(int N, int num_trials) {
    char buffer[200];
    sprintf(buffer, "output/flops-N%d-trials%d.csv", N, num_trials);
    FILE *speed = fopen(buffer, "w");
    for (int n = 0; n < N; n++) {
        printf("%d\n", n);
        double *flops = calloc(NUM_PRECISIONS, sizeof(double));
        for (int t = 0; t < num_trials; t++) {
            double *s = generateRandomNVector(n, 1);
            double *t = generateRandomNVector(n, 1);
            double *g = generateRandomNVector(n, 1);
            for (int p = 0; p < NUM_PRECISIONS; p++) {
                CauchyMultiplier *cm = newCauchyMultiplier(s, t, g, n, PRECISIONS[p]);
                double *fast = fastMultiply(cm);
                flops[p] += cm->adds + cm->muls;
                free(fast);
                free(cm);
            }
            free(s);
            free(t);
            free(g);
        }
        fprintf(speed, "%d, ", n);
        for (int i = 0; i < NUM_PRECISIONS; i++) {
            fprintf(speed, "%0.4e", flops[i] / num_trials);
            if (i != NUM_PRECISIONS - 1) {
                fprintf(speed, ", ");
            }
        }
        fprintf(speed, "\n");
        free(flops);

    }
    fclose(speed);
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
    maxErrorAbs(fastMultiply(cm), slowMultiply(cm), cm->n);
    generateFlopsData(10000,10);
    printVector(fastMultiply(cm), cm->n);
    return 0;
}