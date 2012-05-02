#ifndef REM_H
#define	REM_H
#include "typedefs.h"
#include <gmp.h>

#define maxModules 10		//define max number of modules in advance - now workaround known for me
extern int n;
scalar nextPrime(scalar);

int isPrime(scalar);

scalar* getModules(int, int, int*);

void getInverse(mpz_t, mpz_t, mpz_t*);

scalar getIntInverse(scalar, scalar);

void chineseRemainder(scalar*, scalar*, int, mpz_t);

#endif	/* REM_H */
