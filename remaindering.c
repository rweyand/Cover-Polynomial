#include "remaindering.h"
#include "typedefs.h"
#include <limits.h>
#include <stdio.h>
#include "stdlib.h"
#include <gmp.h>
#include "tools.h"


scalar nextPrime(scalar number) {
    if (number < 0) printf("ERROR!\n");
    scalar i;
    for (i = number - 1; i > n; i--)
        if (isPrime(i)) return i;
    exit(1);
}

int isPrime(scalar number) {
    scalar i;
    for (i = 2; i * i <= number; i++) {
        if (number % i == 0)
            return 0;
    }
    return 1;
}

scalar* getModules(int n, int m, int* numModules) {//initializes module array and returns number of used modules;
    if (n < 10) {
        scalar* mods = (scalar*) malloc(sizeof (scalar));
        mods[0] = (scalar) INT_MAX / 2;
        *numModules = 1;
        return mods;
    }
    mpz_t m_choose_n, maxModule;
    scalar modules[maxModules];
    mpz_init(m_choose_n);
    mpz_init(maxModule);
    mpz_init(maxCoeff);
    int n_ = (unsigned long int) n;
    int m_ = (unsigned long int) m;
    mpz_bin_uiui(m_choose_n, m_, n_);
    mpz_set(maxCoeff, m_choose_n);
    mpz_sqrt(maxModule, maxCoeff);
    if (mpz_cmp_si(maxCoeff,LONG_MAX/n) <= 0) {
        scalar* mods = (scalar*) malloc(sizeof (scalar));
        mods[0] = (scalar) INT_MAX;
        *numModules = 1;
        return mods;
    }
    gmp_printf("Value of maxModule: %Zd \n", maxModule);
    if (mpz_cmp_si(maxModule, LONG_MAX) > 0) {
        fprintf(stderr, "not enough memory to compute coefficients!\n");
        exit(1);
    }
    scalar run = nextPrime((scalar) mpz_get_si(maxModule));
    printf("SCALAR run is %ld \n", run);
    mpz_set_si(maxCoeff, run);
    modules[0] = run;
    int modCount = 1;
    while (mpz_cmp(maxCoeff, m_choose_n) < 0) {
        run = nextPrime(run);
        modules[modCount] = run;
        mpz_mul_si(maxCoeff, maxCoeff, run);
        modCount++;
    }
    scalar* mods = (scalar*) malloc(modCount * sizeof (scalar));
    int i;
    printf("**using the following modules for remaindering: ");
    for (i = 0; i < modCount; i++) {
        mods[i] = modules[i];
        printf("%ld ", mods[i]);
    }
    printf("\n");
    *numModules = modCount;
    mpz_clears(m_choose_n, maxModule, maxCoeff, NULL);
    return mods;
}

void getInverse(mpz_t elem, mpz_t module, mpz_t* inverse) {
    mpz_t y, gcd;
    mpz_inits(y, gcd, *inverse, NULL);
    mpz_gcdext(gcd, y, *inverse, elem, module);
    mpz_clears(y, gcd, NULL);
}

scalar getIntInverse(scalar a_, scalar m) {
    mpz_t elem, module, inverse;
    mpz_inits(elem, module, NULL);
    mpz_set_si(module, m);
    mpz_set_si(elem, a_);
    getInverse(elem, module, &inverse);
    scalar ret = (scalar) mpz_get_si(inverse);
    mpz_clears(elem, module, inverse, NULL);
    return ret;
}

void chineseRemainder(scalar* modules, scalar* values, int numModules, mpz_t number) {//given the used modules and the computed values - determine exact value
    int i;
    mpz_t N;
    mpz_inits(N, number, NULL);
    mpz_set_si(N, 1L);
    for (i = 0; i < numModules; i++) {
        mpz_mul_si(N, N, (long) modules[i]);
    }
    for (i = 0; i < numModules; i++) {
        mpz_t N_i, n_i, inverse;
        mpz_inits(N_i, n_i, NULL);
        mpz_set_si(n_i, (long) modules[i]);
        mpz_cdiv_q(N_i, N, n_i);
        getInverse(n_i, N_i, &inverse);
        mpz_mul(N_i, N_i, inverse);
        mpz_mul_si(N_i, N_i, (long) values[i]);
        mpz_add(number, number, N_i);
        mpz_clears(N_i, n_i, inverse, NULL);
    }
    mpz_fdiv_r(number, number, N);
    mpz_clears(N,NULL);
}

