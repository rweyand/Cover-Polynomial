/* 
 * File:   main.c
 * Author: Reinhard Weyand
 *
 * Implementation of the fast evaluation algorithm for the coefficients of the cover polynomial
 * as presented in "Computing The Tutte Polynomial In Vertex-Exponential Time" by Björklund et. al. (Appendix p.3-4).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <unistd.h>
#include <limits.h>

#include "remaindering.h"
#include "tools.h"
#include "matrix.h"
#include "poly.h"
#include "typedefs.h"
scalar* P_U_I;
scalar* C_U_J; //C-polynomial to the power of J
scalar* PC_IJ; //(P^i)*(C*j)
polynom P_U_small;
polynom C_U_small;

void freeAll() { //cascade for cleaning up Memory   
    int i, j;
    for (i = 0; i < n + 1; i++)
        for (j = 0; j < n + 1; j++)
            free(VALS(i, j));
    free(vals);
    free(PC_IJ);
    free(P_U_I);
    free(C_U_J);
    free(modules);
    free(coefficients);

}

void init() {
    P_U_I = (scalar*) malloc(((n + 1)*(n + 1) + 1) * sizeof (scalar));
    P_U_I[0] = 0;
    C_U_J = (scalar*) malloc(((n + 1)*(n + 1) + 1) * sizeof (scalar));
    C_U_J[0] = 0;
    PC_IJ = (scalar*) malloc(((n + 1)*(n + 1)*(n + 1) + 1) * sizeof (scalar));
    PC_IJ[0] = 0;
}

/*
 * Compute the number of directed walks of length l in G[S] by polynomial exponentiation
 * S is the subset of vertices under consideration
 * s and t are the start- and end-vertex for the walk
 * l is the length of the walks to count, i.e. the exponent for matrix exponentiation
 */
scalar w(set S, int s, int t, int l, int m) {
    scalar module = modules[m];
    int i, j; //positions in submatrix
    int size = getNumOfNodes(S);
    i = getPos(s, S);
    j = getPos(t, S);
    if (!(ELEM(s, S) && ELEM(t, S)))
        return 0;
    if ((s == t) && (l == 0))
        return 1;
    scalar* subMatrix = (scalar*) malloc(size * size * sizeof (scalar));
    scalar* result = (scalar*) malloc(size * size * sizeof (scalar));
    getSubMatrix(S, subMatrix);
    matrixPow(subMatrix, l, result, size, module);
    scalar res = matrix(result, i, j, size);
    //  printf("freeing matrices of size %d x %d \n computed result: %d \n",size,size,res);	
    free(result);
    //	printf("freeing submatrix \n");
    free(subMatrix);
    return res;
}

void computeP_U_S(set U, set S, polynom* subsetPolynom, int m) {
    scalar module = modules[m];
    int s, t, k, i, sign;
    int size = popCnt(U^S);
    int size_S = popCnt(S);
    polynom poly = *subsetPolynom;
    for (i = 0; i < n + 3; i++)
        poly[i] = 0;
    if (U == 0 || S == 0)return;
    for (s = 0; s < n; s++) {
        for (t = 0; t < n; t++) {
            for (k = 0; k <= size; k++) {
                sign = (k % 2 == 0) ? 1 : -1;
                poly[size_S + k + 1] = MYMOD(poly[size_S + k + 1] + MYMOD(sign * (w(S, s, t, size_S + k - 1, m) * binCoeff(size, k, module)), module), module);
                if (poly[0] < (size_S + k + 1))
                    poly[0] = size_S + k + 1;
            }
        }
    }
}

void computeP_U(set U, int m) {
    scalar module = modules[m];
    int i, S, size;
    size = popCnt(U);
    polynom subsetPolynom;
    P_U_small = (polynom) malloc((n + 3) * sizeof (scalar));
    for (i = 0; i < n + 3; i++)
        P_U_small[i] = 0;
    for (S = 0; S <= U; S++) {
        if ((U | S) == U) {
            subsetPolynom = (scalar*) malloc((n + 3) * sizeof (scalar));
            computeP_U_S(U, S, &subsetPolynom, m);
            polyAdd(P_U_small, subsetPolynom, P_U_small, module);
            free(subsetPolynom);
        }
    }
}

void computeC_U_S(set U, set S, polynom* subsetPolynom, int m) {
    scalar module = modules[m];
    int s, k, i, sign;
    int size = popCnt(U^S);
    int size_S = popCnt(S);
    polynom poly = *subsetPolynom;
    for (i = 0; i < n + 2; i++)
        poly[i] = 0;
    if (S == 0)
        return;
    for (k = 0; k <= size; k++) {
        for (s = 0; s < n; s++) {
            sign = (k % 2 == 0) ? 1 : -1;
            poly[size_S + k + 1] = MYMOD(poly[size_S + k + 1] + MYMOD(sign * (w(S, s, s, size_S + k, m) * binCoeff(size, k, module)), module), module);
            if (poly[0] < (size_S + k + 1))
                poly[0] = size_S + k + 1;
        }
    }
}

void computeC_U(set U, int m) {

    scalar module = modules[m];
    int i, S, size;
    polynom subsetPolynom;
    C_U_small = (polynom) malloc((n + 2) * sizeof (polynom));
    for (i = 0; i < n + 2; i++)
        C_U_small[i] = 0;
    size = popCnt(U);
    for (S = 0; S <= U; S++) {
        if ((S | U) == U) {
            subsetPolynom = (polynom) malloc((n + 2) * sizeof (scalar));
            computeC_U_S(U, S, &subsetPolynom, m);
            polyAdd(C_U_small, subsetPolynom, C_U_small, module);
            free(subsetPolynom);
        }
    }
    for (i = 2; i < n + 2; i++) {
        C_U_small[i] /= (i - 1);
    }

}

void computeC_D(int i, int j, int m) {
    scalar module = modules[m];
    VAL(i, j, m) = 0;
    int sign, size;
    set S;
    for (S = 1; S < pow2; S++) {
        computeP_U(S, m);
        computeC_U(S, m);
        size = popCnt(S^(pow2 - 1));
        sign = size & 1 ? -1 : 1;
        polyPow(P_U_small, i, P_U_I, module);
        polyPow(C_U_small, j, C_U_J, module);
        polyMul(P_U_I, C_U_J, PC_IJ, module);
        VAL(i, j, m) = MYMOD(VAL(i, j, m) + ((PC_IJ[0] < n + 1) ? 0 : MYMOD(sign * PC_IJ[n + 1], module)), module);
        free(P_U_small);
        free(C_U_small);
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        fprintf(stderr, "no file specified\n");
        exit(EXIT_FAILURE);
    }
    int opt = -1;
    char* filename = (char*) malloc(50 * sizeof (char));
    while ((opt = getopt(argc, argv, "f:")) != -1) {
        switch (opt) {
            case 'f':
                strcpy(filename, optarg);
                parseGraph(filename);
                break;
            default:
                fprintf(stderr, "no file specified \n");
        }
    }
    init();
    pow2 = mypow(n);
    modules = getModules(n, m, &numModules);
    int i, j, m;
    vals = (scalar**) malloc((n + 1)*(n + 1) * sizeof (scalar*));
    for (i = 0; i < n + 1; i++)
        for (j = 0; j < n + 1; j++)
            VALS(i, j) = (scalar*) malloc(numModules * sizeof (scalar));
    for (m = 0; m < numModules; m++) {
        for (i = 0; i <= n; i++) {
            for (j = 0; j <= n; j++) {
                computeC_D(i, j, m);
            }
        }
    }
    mpz_t i_j_fac, i_fac, j_fac;
    mpz_init(i_fac);
    mpz_init(j_fac);
    mpz_init(i_j_fac);
    coefficients = (mpz_t*) malloc((n + 1)*(n + 1) * sizeof (mpz_t));
    for (i = 0; i < n + 1; i++)
        for (j = 0; j < n + 1; j++) {
            mpz_fac_ui(i_fac, (unsigned long) i);
            mpz_fac_ui(j_fac, (unsigned long) j);
            mpz_mul(i_j_fac, i_fac, j_fac);
            chineseRemainder(modules, VALS(i, j), numModules, COEFFS(i, j));
            mpz_div(COEFFS(i, j), COEFFS(i, j), i_j_fac);
        }
    printCoefficients();
    free(filename);
    mpz_clears(i_fac, j_fac, i_j_fac, NULL);
    freeAll();
    return (EXIT_SUCCESS);
}
