#include <stdio.h>
#include <math.h>
#include <string.h>

#include "utils.h"
#include "ZeroFuncao.h"
#include "DoubleType.h"

#include <stdint.h>
#include <math.h>

int AlmostEqualUlps(Double_t a, Double_t b, int maxULPs) {
    // Se um é negativo e outro é positivo, só podem ser iguais se forem +0.0 e -0.0
    if (a.parts.sign != b.parts.sign) {
        return (a.f == 0.0 && b.f == 0.0); 
    }

    // Calcula a diferença em ULPs corretamente
    int64_t ulpsDiff = llabs(a.i - b.i);

    // Verifica se a diferença está dentro do limite aceitável
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
        erro = fabs(xm_new - xm_old);

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
        erro = fabs(xm_new - xm_old);

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


