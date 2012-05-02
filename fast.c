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
#include <gmp.h>

#include "remaindering.h"
#include "tools.h"
#include "matrix.h"
#include "poly.h"
#include "typedefs.h"
scalar** P_U_I;
scalar** C_U_J; //C-polynomial to the power of J
scalar** PC_IJ; //(P^i)*(C*j)

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

/*
 * compute the polynomial
 
 
 */
void computeC_X(int m) { //compute number of spanning directed cycles for X according to björklund et.al. by fast moebius inversi
    printf("computing C_X for module %ld \n", modules[m]);
    scalar module = modules[m];
    int S, j, k, s;
    scalar* b_flat = (scalar*) malloc(pow2 * (n + 1) * sizeof (scalar));
    int size;
    for (s = 0; s < n; s++)
        for (k = 0; k <= n; k++) {
            for (S = 0; S < pow2; S++) {
                //  smallC(S, m) = 0;
                size = popCnt(S);
                B(S, 0) = w(S, s, s, k, m);
                for (j = 1; j <= n; j++) {
                    if ((1 << (j - 1)) & S) {
                        B(S, j) = MYMOD(B(S, j - 1) - B(S^(1 << (j - 1)), j - 1), module);
                    } else {
                        B(S, j) = B(S, j - 1);
                    }
                }
                if (k == size) {
                    smallC(S, m) = MYMOD(smallC(S, m) + B(S, n), module);
                }
            }
        }
    for (S = 1; S < pow2; S++) {
        if (mpz_cmp_si(maxCoeff, LONG_MAX / n) > 0) {
            scalar intInverse = getIntInverse(popCnt(S), module);
            smallC(S, m) = MYMOD(smallC(S, m) * intInverse, module);
        } else
            smallC(S, m) /= popCnt(S);
    }
    free(b_flat);
    printf("done for %ld \n", modules[m]);
}

void computeP_X(int m) { //compute number of spanning directed paths for X according to björklund et.al. by fast moebius inversion
    printf("computing P_X for module %ld \n", modules[m]);
    scalar module = modules[m];
    int k, j, s, t;
    set S;
    int size;
    scalar* a_flat = (scalar*) malloc(pow2 * (n + 1) * sizeof (scalar));
    for (s = 0; s < n; s++) {
        for (t = 0; t < n; t++) {
            for (k = 0; k <= n; k++) {
                for (S = 0; S < pow2; S++) {
                    size = popCnt(S);
                    A(S, 0) = w(S, s, t, k, m);
                    for (j = 1; j <= n; j++) {
                        if ((1 << (j - 1)) & S)
                            A(S, j) = MYMOD(A(S, j - 1) - A(S^(1 << (j - 1)), j - 1), module); //f[j](X) = - f[j-1](X\{j}) + f[j-	 
                        else
                            A(S, j) = A(S, j - 1); //f[j](X) = f[j-1](X)
                    }
                    if (k == (size - 1)) {
                        smallP(S, m) = MYMOD(smallP(S, m) + A(S, n), module);
                    }
                }

            }
        }
    }
    free(a_flat);
    printf("done for %ld \n", modules[m]);
}

void computeC_U_moebius(int m) {
    computeC_X(m);
    printf("computing generating function for spanning directed cycles with module %ld \n", modules[m]);
    scalar module = modules[m];
    int S, i, size, j;
    scalar** c_flat = (scalar**) malloc(pow2 * (n + 1) * sizeof (scalar*));
    for (S = 0; S < pow2; S++) {
        for (j = 0; j < (n + 1); j++) {
            C_POLY(S, j) = (scalar*) malloc((n + 2) * sizeof (scalar));
            C_POLY(S, j)[0] = 0;
        }
        C__U(S, m) = (scalar*) malloc((n + 2) * sizeof (scalar));
        size = popCnt(S);
        for (i = 0; i < n + 2; i++)
            C_POLY(S, 0)[i] = 0;
        C_POLY(S, 0)[size + 1] = smallC(S, m);
        if (smallC(S, m) != 0)
            C_POLY(S, 0)[0] = size + 1;
        for (j = 1; j <= n; j++) {
            if (ELEM(j - 1, S)) {
                polyAdd(C_POLY(S, j - 1), C_POLY(S^(1 << (j - 1)), j - 1), C_POLY(S, j), module);
            } else {
                polyCopy(C_POLY(S, j - 1), C_POLY(S, j));
            }
        }
        for (i = 0; i <= C_POLY(S, n)[0]; i++) {
            C__U(S, m)[i] = C_POLY(S, n)[i];
        }
    }
    for (i = 0; i < (n + 1); i++) {
        for (S = 0; S < pow2; S++) {
            free(C_POLY(S, i));
        }
    }
    free(c_flat);
    printf("done with cycles for %ld \n", modules[m]);
}

void computeP_U_moebius(int m) {
    computeP_X(m);
    printf("computing generating function for spanning directed paths with module %ld \n", modules[m]);
    scalar module = modules[m];
    int S, i, size, j;

    scalar** c_flat = (scalar**) malloc(pow2 * (n + 1) * sizeof (scalar*));
    for (S = 0; S < pow2; S++) {
        for (j = 0; j < (n + 1); j++) {
            C_POLY(S, j) = (scalar*) malloc((n + 3) * sizeof (scalar));
            C_POLY(S, j)[0] = 0;
        }
        P__U(S, m) = (scalar*) malloc((n + 3) * sizeof (scalar));
        size = popCnt(S);
        for (i = 0; i < n + 1; i++)
            C_POLY(S, 0)[i] = 0;
        C_POLY(S, 0)[size + 1] = smallP(S, m);
        if (smallP(S, m) != 0)
            C_POLY(S, 0)[0] = size + 1;
        for (j = 1; j <= n; j++) {
            //		printf("S = %d, j = %d!\n",S,j);
            if (ELEM(j - 1, S)) {
                polyAdd(C_POLY(S, j - 1), C_POLY(S^(1 << (j - 1)), j - 1), C_POLY(S, j), module);
            } else {
                polyCopy(C_POLY(S, j - 1), C_POLY(S, j));

            }
        }
        for (i = 0; i <= C_POLY(S, n)[0]; i++) {
            P__U(S, m)[i] = C_POLY(S, n)[i];
        }
    }
    for (i = 0; i < (n + 1); i++) {
        for (S = 0; S < pow2; S++) {
            free(C_POLY(S, i));
        }
    }

    free(c_flat);
    printf("done with paths for %ld \n", modules[m]);
}

void computeC_D(int m) {
    int i, j, S, sign, size;
    scalar module = modules[m];
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
            VAL(i, j, m) = 0;
            for (S = 0; S < pow2; S++) {
                size = popCnt(S^(pow2 - 1));
                sign = size & 1 ? -1 : 1;
                polyPow(P__U(S, m), i, P_U_I[m], module);
                polyPow(C__U(S, m), j, C_U_J[m], module);
                polyMul(P_U_I[m], C_U_J[m], PC_IJ[m], module);
                VAL(i, j, m) = MYMOD(VAL(i, j, m) + ((PC_IJ[m][0] < n + 1) ? 0 : MYMOD(sign * PC_IJ[m][n + 1], module)), module);
            }
        }
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

    pow2 = mypow(n);
    modules = getModules(n, m, &numModules);

    int i, j, m;
    P_U_I = (scalar**) malloc(numModules * sizeof (polynom));
    C_U_J = (scalar**) malloc(numModules * sizeof (polynom));
    PC_IJ = (scalar**) malloc(numModules * sizeof (polynom));
    for (m = 0; m < numModules; m++) {
        P_U_I[m] = (polynom) malloc(((n + 1)*(n + 1) + 1) * sizeof (scalar));
        C_U_J[m] = (polynom) malloc(((n + 1)*(n + 1) + 1) * sizeof (scalar));
        PC_IJ[m] = (polynom) malloc(((n + 1)*(n + 1)*(n + 1) + 1) * sizeof (scalar));
        P_U_I[m][0] = 0;
        C_U_J[m][0] = 0;
        PC_IJ[m][0] = 0;
    }
    vals = (scalar**) malloc((n + 1)*(n + 1) * sizeof (scalar*));
    for (i = 0; i < n + 1; i++)
        for (j = 0; j < n + 1; j++)
            VALS(i, j) = (scalar*) malloc(numModules * sizeof (scalar));
    c = (scalar*) malloc(pow2 * numModules * sizeof (scalar));
    p = (scalar*) malloc(pow2 * numModules * sizeof (scalar));
    P_U = (scalar**) malloc(pow2 * numModules * sizeof (scalar*));
    C_U = (scalar**) malloc(pow2 * numModules * sizeof (scalar*));

#pragma omp parallel
    {
#pragma omp for private(i,j)
        for (m = 0; m < numModules; m++) {
            computeP_U_moebius(m);
            computeC_U_moebius(m);
            computeC_D(m);
        }
    }

    printf("done with main loop \n");
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
    set S;
    for (m = 0; m < numModules; m++)
        for (S = 0; S < pow2; S++) {
            free(C__U(S, m));
            free(P__U(S, m));
        }
    free(p);
    free(c);
    free(C_U);
    free(P_U);
    return (EXIT_SUCCESS);
}
