#ifndef MATRIX_H
#define	MATRIX_H
#include "typedefs.h"

#define matrix(m,i,j,n) m[(n)*(i)+(j)]


void matrixMul(matrix, matrix, matrix, int, scalar);

void matrixPow(matrix, int, matrix, int, scalar);


#endif
