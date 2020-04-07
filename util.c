//
// Created by kavalee on 4/6/20.
//

#include "util.h"
#include <stdio.h>



void printVector(double* vec, int len) {
    printf("[");
    for (int i = 0; i < len; i++) {
        printf("%.4e ," , vec[i]);
    }
    printf("]\n");
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