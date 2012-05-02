#include "poly.h"
#include <stdio.h>

void polyAdd(polynom a, polynom b, polynom r, int module) { //straight forward addition of polynoms
        int i;
        int sizea = a[0];
        int sizeb = b[0];
        polynom x, y;
        if (sizea != sizeb) {
            x = sizea < sizeb ? b : a;
            y = sizea <= sizeb ? a : b;
        } else {
            x = a;
            y = b;
        }
        sizea = x[0]; //bigger size
        sizeb = y[0]; //smaller size
        r[0] = sizea > 0 ? sizea : 0;
        for (i = 1; i <= sizea; i++) {
            r[i] = (i <= sizeb) ? (MYMOD(x[i] + y[i],module)) : x[i];
        }
    }

    void polyMul(polynom a, polynom b, polynom r, int module) { //multiplication of given polynomials - FFT would do better
        int i, j, sizer;
        int sizea = a[0];
        int sizeb = b[0];
        if (sizea == 0 || sizeb == 0) {
            r[0] = 0;
            return; //if one is empty set to zero
        }

        for (i = 0; i <= sizea + sizeb - 1; i++) {
            r[i] = 0;
        }
        for (i = 1; i <= sizea; i++) {
            for (j = 1; j <= sizeb; j++) {
                r[i + j - 1] = MYMOD(r[i + j - 1] + MYMOD(a[i] * b[j],module),module);
            }
        }
        sizer = (sizea + sizeb - 1);
        r[0] = sizer > 0 ? sizer : 0;

    }

    void polyPow(polynom a, int exp, polynom r, scalar module) {//repeated multiplication to compute polynomial power
        int i, j, size;
        if (exp == 0) { //set polynomial to 1
            r[0] = 1;
            r[1] = 1;
            for (i = 2; i <= n; i++)
                r[i] = 0;
            return;
        } else if (a[0] == 0) { //if empty, set to zero
            r[0] = 0;
            for (i = 1; i <= n; i++)
                r[i] = 0;
        }

        r[0] = 0;
        scalar buf[a[0] * exp + 1]; //standard repeated polynomial multiplication
        for (i = 0; i <= a[0]; i++)
            buf[i] = a[i];
        for (i = 1; i < exp; i++) {
            polyMul(buf, a, r, module);
            size = r[0];
            for (j = 0; j <= size; j++) {
                buf[j] = r[j];
                r[j] = 0;
            }

        }
        for (i = 0; i <= buf[0]; i++) {
            r[i] = buf[i];
        }
    }

    void printPoly(polynom poly) {
        int i;
        for (i = 0; i <= poly[0] + 1; i++)
            printf("%ld, ", poly[i]);
        printf("\n");
    }
