/* 
 * File:   tools.h
 * Author: reinhard
 *
 * Created on 30. März 2011, 13:30
 */

#ifndef TOOLS_H
#define	TOOLS_H


#define MAXN 30
#include "typedefs.h"
#include <gmp.h>
int n; //number of nodes
int m; // number of edges
int numModules;
int adj[MAXN][MAXN];
int pow2; //2 to the power of n
scalar *p; //p value from bhkk-paper
scalar *c; //c value from bhkk-paper
scalar** P_U; //polynomial according to paper of björklund et.al. : size at first position, exponents in increasing order
scalar** C_U;
scalar* modules;
scalar** vals;
mpz_t* coefficients;
mpz_t maxCoeff;


// definitions for more dimensional arrays
#define A(S,j) a_flat[(S)*(n+1)  + (j)]
#define B(S,j) b_flat[(S)*(n+1) + (j)]
#define ELEM(i,S) ((1 << (i)) & S)
#define C_FLAT(S,j,i) c_flat[(S)*(n+1)  + (j)][(i)]
#define C_POLY(S,j) c_flat[(S)*(n+1) + (j)]
#define C__U(S,m) C_U[(S)*numModules+(m)]
#define P__U(S,m) P_U[(S)*numModules+(m)]
#define myPow(x) (1 << (x))
#define VALS(i,j) vals[(i)*(n+1) + (j)]
#define VAL(i,j,m) vals[(i)*(n+1) + (j)][(m)]
#define COEFFS(i,j) coefficients[(i)*(n+1)+(j)]
#define smallC(S,m) c[(S)*numModules+(m)]
#define smallP(S,m) p[(S)*numModules+(m)]


void printCoefficients();

void parseGraph(char*);

void freeAll();

void init();

int getPos(int, set);

void polyCopy(polynom, polynom);

void getSubsets(set*, set, int);

int expandSubset(set, int, int*);

void getSubMatrix(int, scalar*);

int getNumOfNodes(set);

int popCnt(set);

scalar mypow(int);

scalar fac(scalar, scalar);

int binCoeff(scalar, scalar, scalar);


#endif	/* TOOLS_H */

