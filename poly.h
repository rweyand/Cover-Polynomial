#ifndef POLY_H
#define	POLY_H
#include "typedefs.h"
extern int n;

void polyAdd(polynom, polynom, polynom, int);

void polyMul(polynom, polynom, polynom, int);

void polyPow(polynom, int, polynom, scalar);

void printPoly(polynom);

#endif	/* POLY_H */
