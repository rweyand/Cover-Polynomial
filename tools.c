#include "tools.h"
#include "matrix.h"
#include "poly.h"
#include "remaindering.h"
#include <gmp.h>
#include <stdlib.h>
#include <stdio.h>
#include "typedefs.h"

void printCoefficients() { //print coefficients in rectengular format
    int i = 0, j = 0;
    printf("\n ------------------------------------------------------------- \n coefficients in rectengular format:\n(row -> exponent of y, column -> exponent of falling powers of x) \n\n\n");
    printf("   ");
    for (i = 0; i <= n; i++)
        printf("%d ", i);

    printf("\n");
    printf("  "
            " ");
    for (i = 0; i <= 2 * n; i++)
        printf("-");
    printf("\n");

    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {

            if (j != 0) printf(" ");
            else printf("%d |", i);
            if(i+j<=n)gmp_printf("%Zd", COEFFS(i, j));
            else printf(" ");
            mpz_clear(COEFFS(i,j));
            if (j == n)printf("\n");
        }
    }
    printf("\n");
}

void parseGraph(char* filename) { //routine to parse graph from file (specified format in wiki)
    FILE* f = fopen(filename, "r");
    if (f == NULL) {
        fprintf(stderr, "could not open file %s\n", filename);
        exit(EXIT_FAILURE);
    }
    int e;
    int i = 0;
    int j = 0;
    m = 0;
    fscanf(f, "%d", &n);
    if (n > MAXN)fprintf(stderr, "exceeded MAXN, please adjust MAXN in tools.h\n");
    while (!feof(f)) {
        fscanf(f, "%d", &e);
        if (!isspace(e)) {
            if (e < 0 || e > 1) {
                fprintf(stderr, "value different from 1 or 0 in adjacency matrix\n");
                exit(EXIT_FAILURE);
            } else {
                if (j < n && i < n) {
                    adj[j][i] = e;
                    m = m + adj[j][i];
                    j = j + (i == (n - 1) ? 1 : 0);
                    i = (i + 1) % n;
                }
            }
        }
    }
    fclose(f);
}





int getPos(int node, set S) { //get Position of node in subset S if reduced to induced graph
    if (S == 0)return 0;
    int count = 0;
    int mask = 1;
    while (mask != (1 << node)) {
        if (mask & S)
            count++;
        mask = mask << 1;
    }
    return count;
}

void polyCopy(polynom a, polynom r) { //copy polynomial a to memory location r
    int i;
    for (i = 0; i <= a[0]; i++)
        r[i] = a[i];
}

void getSubsets(set *result, set superSet, int subSetSize) { // avoid iterating over non-subsets of the superset by shrinking, iterating over all subsets and memorizing positions of elements
    int i = 0, j = 0;
    int shrunkSet = myPow(subSetSize);
    set* elements = (set*) malloc(n * sizeof (set));
    while (superSet > 0) {
        if (1 & superSet) {
            elements[j++] = i;
        }
        superSet = superSet >> 1;
        i++;
    }
    set S;
    for (S = 0; S < shrunkSet; S++) {
        result[S] = expandSubset(S, n, elements);
    }
    free(elements);
}

int expandSubset(set subSet, int n, int* positions) { // from shrunk representation and positions, reconstruct corresponding subsets produced by getSubsets()   
    int set = 0;
    int i;
    for (i = 0; i < n; i++) {
        set = set | (((1 << i) & subSet ? 1 : 0) << positions[i]);
    }
    return set;
}

void getSubMatrix(int nodes, scalar* result) { //extract submatrix from given matrix to result
    int i, j, k, l;
    k = 0;
    l = 0;
    int size = getNumOfNodes(nodes);
    for (i = 0; i < n; i++) {
        if (ELEM(i, nodes)) {
            for (j = 0; j < n; j++) {
                if (ELEM(j, nodes)) {
                    matrix(result, k, l, size) = adj[i][j];
                    l++;
                }
            }
            k++;
        }
        l = 0;
    }
}

int getNumOfNodes(set U) { //get size of vertex subset by adding number of 1s in binary representation
    int n = 0;
    set S = U;
    while (S > 0) {
        n += (S & 1);
        S = S >> 1;
    }
    return n;
}

int popCnt(set U) {
    switch (sizeof (set)) {
        case 4:
            U = U - ((U >> 1) & 0x55555555);
            U = (U & 0x33333333) + ((U >> 2) & 0x33333333);
            U = (U & 0x0F0F0F0F) + ((U >> 4) & 0x0F0F0F0F);
            U = (U & 0x00FF00FF) + ((U >> 8) & 0x00FF00FF);
            return (U & 0x0000FFFF) + ((U >> 16) & 0x0000FFFF);
        default:
            return getNumOfNodes(U);
    }
}

scalar mypow(int x) { //return 2 to the power of x
    return (scalar) (1 << x);
}

scalar fac(scalar x, scalar module) { //straight forward computation of faculty
    int i;
    scalar fac = 1;
    for (i = 1; i <= x; i++)
        fac = MYMOD((fac * i), module);
    return fac;
}

int binCoeff(scalar n, scalar k, scalar module) {
    return fac(n, module) / (fac(k, module) * fac(n - k, module));
}
