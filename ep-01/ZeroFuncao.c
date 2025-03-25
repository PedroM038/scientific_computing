#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utils.h"
#include "ZeroFuncao.h"
#include "DoubleType.h"

#include <stdint.h>
#include <math.h>

int AlmostEqualUlps(Double_t a, Double_t b, int maxULPs) {
   if (a.parts.sign != b.parts.sign) {
        return (a.f == 0.0 && b.f == 0.0); 
    }
   int64_t ulpsDiff = llabs(a.i - b.i);
   return ulpsDiff <= maxULPs;
}

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    real_t b = 0;
    real_t c = 0;
    for (int i = p.grau; i > 0; i--) {
        b = p.p[i] + x*b;
        c = b + x*c;
    }
    b = b*x + p.p[0];
    *px = b;
    if (dpx) *dpx = c;
}

void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    *px = 0;
    if (dpx) *dpx = 0;

    for (int i = 0; i <= p.grau; i++) {
        *px += p.p[i] * pow(x, i);
    }

    if (dpx) {
        for (int i = 1; i <= p.grau; i++) {
            *dpx += i * p.p[i] * pow(x, i-1);
        }
    }
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao(Polinomio p, real_t xl, real_t xu, int criterioParada, int *it, real_t *raiz) {
    real_t xm_old, xm_new;
    real_t fa, fb, dx;
    real_t erro;

    fa = 0.0;
    fb = 0.0;
    dx = 0.0;
    xm_old = 0.0;
    erro = 0.0;
    xm_new = (xl + xu) / 2;
    *it = 0;

    calcPolinomio_rapido(p, xl, &fa, &dx);
    calcPolinomio_rapido(p, xm_new, &fb, &dx);
    if (fa * fb < 0) xu = xm_new;
    else if (fa * fb > 0) xl = xm_new;
    else return xm_new;

    do {
        (*it)++;
        xm_old = xm_new;
        xm_new = (xl + xu) / 2;
        calcPolinomio_rapido(p, xl, &fa, &dx);
        calcPolinomio_rapido(p, xm_new, &fb, &dx);
        erro = fabs(xm_new - xm_old)/fabs(xm_new);

        if (fa * fb < 0) xu = xm_new;
        else if (fa * fb > 0) xl = xm_new;

        Double_t a, b;
        a.f = xm_new;
        b.f = xm_old;

        if ((criterioParada == 1 && fabs(erro) < EPS) ||
            (criterioParada == 2 && fabs(fb) <= ZERO) ||
            (criterioParada == 3 && AlmostEqualUlps(a, b, ULPS))) {
            break;
        }

    } while (*it < MAXIT);

    *raiz = xm_new;
    return erro;
}

real_t bisseccao_lento(Polinomio p, real_t xl, real_t xu, int criterioParada, int *it, real_t *raiz) {
    real_t xm_old, xm_new;
    real_t fa, fb, dx;
    real_t erro;

    fa = 0.0;
    fb = 0.0;
    dx = 0.0;
    xm_old = 0.0;
    erro = 0.0;
    xm_new = (xl + xu) / 2;
    *it = 0;

    calcPolinomio_lento(p, xl, &fa, &dx);
    calcPolinomio_lento(p, xm_new, &fb, &dx);
    if (fa * fb < 0) xu = xm_new;
    else if (fa * fb > 0) xl = xm_new;
    else return xm_new;

    do {
        (*it)++;
        xm_old = xm_new;
        xm_new = (xl + xu) / 2;
        calcPolinomio_lento(p, xl, &fa, &dx);
        calcPolinomio_lento(p, xm_new, &fb, &dx);
        erro = fabs(xm_new - xm_old)/fabs(xm_new);

        if (fa * fb < 0) xu = xm_new;
        else if (fa * fb > 0) xl = xm_new;

        Double_t a, b;
        a.f = xm_new;
        b.f = xm_old;

        if ((criterioParada == 1 && fabs(erro) < EPS) ||
            (criterioParada == 2 && fabs(fb) <= ZERO) ||
            (criterioParada == 3 && AlmostEqualUlps(a, b, ULPS))) {
            break;
        }

    } while (*it < MAXIT);

    *raiz = xm_new;
    return erro;
}

real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz){
    real_t x_old, x_new;
    real_t f, df;
    real_t erro;

    x_old = 0.0;
    x_new = x0;
    erro = 0.0;
    *it = 0;

    do {
        (*it)++;
        x_old = x_new;
        calcPolinomio_rapido(p, x_old, &f, &df);
        if (df == 0) {
            x_new = x_old;
            break;
        };
        x_new = x_old - f/df;
        erro = fabs(x_new - x_old);

        Double_t a, b;
        a.f = x_new;
        b.f = x_old;

        if (criterioParada == 1 && fabs(erro) < EPS) break;
        if (criterioParada == 2 && fabs(f) <= ZERO) break;
        if (criterioParada == 3 && AlmostEqualUlps(a, b, ULPS)) break;

    } while (*it < MAXIT);

    *raiz = x_new;
    return erro;
}

real_t newtonRaphson_lento (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz){
    real_t x_old, x_new;
    real_t f, df;
    real_t erro;

    x_old = 0.0;
    x_new = x0;
    erro = 0.0;
    *it = 0;

    do {
        (*it)++;
        x_old = x_new;
        calcPolinomio_lento(p, x_old, &f, &df);
        if (df == 0) {
            x_new = x_old;
            break;
        };
        x_new = x_old - f/df;
        erro = fabs(x_new - x_old)/fabs(x_new);

        Double_t a, b;
        a.f = x_new;
        b.f = x_old;

        if (criterioParada == 1 && fabs(erro) < EPS) break;
        if (criterioParada == 2 && fabs(f) <= ZERO) break;
        if (criterioParada == 3 && AlmostEqualUlps(a, b, ULPS)) break;

    } while (*it < MAXIT);

    *raiz = x_new;
    return erro;
}
