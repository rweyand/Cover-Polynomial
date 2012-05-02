#include "matrix.h"
#include "typedefs.h"
#include <stdlib.h>

extern int n;

void matrixMul(matrix a, matrix b, matrix result, int size, scalar module) { //basic matrix multiplication - FFT would indeed be better
    matrix buf = (matrix) malloc(sizeof (scalar) * size * size);
    int i, j, k;
       

    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            matrix(buf, i, j, size) = 0;
            for (k = 0; k < size; k++) {
                matrix(buf, i, j, size) = MYMOD(matrix(buf, i, j, size) + MYMOD(matrix(a, i, k, size) * matrix(b, k, j, size), module), module);
            }
        }
    }

    for (i = 0; i < size * size; i++) {
        result[i] = buf[i];
    }


    free(buf);
}

void matrixPow(matrix a, int k, matrix result, int size, scalar module) { // determination of the number of directed walks from one node to the other by adjacency matrix potentiation
    int i, exp;
    for (i = 0; i < size * size; i++) {
        result[i] = a[i];
    }
    exp = k - 1;
    while (exp != 0) {
        if ((exp % 2) == 1) {
            matrixMul(a, result, result, size, module);
            exp--;
        }
        matrixMul(a, a, a, size, module);
        exp = exp / 2;
    }
    //printf("done matrix pow \n");
}
