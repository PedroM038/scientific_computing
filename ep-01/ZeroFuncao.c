#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz) {
    real_t x_new = x0, x_old, erro, fx, dfx;
    *it = 0;

    do {
        x_old = x_new;

        // Calcula f(x) e f'(x)
        calcPolinomio_rapido(p, x_old, &fx, &dfx);

        // Verifica se a derivada é muito pequena (evita divisão por zero)
        if (fabs(dfx) < ZERO) {
            fprintf(stderr, "Erro: Derivada próxima de zero.\n");
            return (real_t)-1;
        }

        // Atualiza x usando a fórmula de Newton-Raphson
        x_new = x_old - fx / dfx;

        // Calcula erro relativo
        erro = fabs(x_new - x_old) / fabs(x_new);
        (*it)++;
    } while (erro > criterioParada && (*it) < MAXIT);

    *raiz = x_new;
    return erro;
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz)
{
    real_t xm_new, xm_old, erro, fa, fb, fxm;
    real_t verificador;

    calcPolinomio_rapido(p, a, &fa, NULL);
    calcPolinomio_rapido(p, b, &fb, NULL);
    *it = 0;

    if (fa*fb > 0) {
        fprintf(stderr, "Erro: f(a) e f(b) têm mesmo sinal\n");
        return (real_t) -1;
    } else {
        do {
            xm_old = xm_new;
            xm_new = (a + b)/2;
            calcPolinomio_rapido(p, xm_new, &fxm, NULL);
            verificador = fa*fxm;
            if (verificador < 0) {
                b = xm_new;
            } else if (verificador > 0) {
                a = xm_new;
                fa = fxm;
            } else {
                erro = 0.0;
                break;
            }
            erro = fabs(xm_new - xm_old)/fabs(xm_new);
            (*it)++;
        } while (erro > criterioParada && (*it) < MAXIT);
    }
    *raiz = xm_new;
    return erro;
}

real_t newtonRaphson_lento(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz) {
    real_t x_new = x0, x_old, erro, fx, dfx;
    *it = 0;

    do {
        x_old = x_new;

        // Calcula f(x) e f'(x)
        calcPolinomio_lento(p, x_old, &fx, &dfx);

        // Verifica se a derivada é muito pequena (evita divisão por zero)
        if (fabs(dfx) < ZERO) {
            fprintf(stderr, "Erro: Derivada próxima de zero.\n");
            return (real_t)-1;
        }

        // Atualiza x usando a fórmula de Newton-Raphson
        x_new = x_old - fx / dfx;

        // Calcula erro relativo
        erro = fabs(x_new - x_old) / fabs(x_new);
        (*it)++;
    } while (erro > criterioParada && (*it) < MAXIT);

    *raiz = x_new;
    return erro;
}

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao_lento (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz)
{
    real_t xm_new, xm_old, erro, fa, fb, fxm;
    real_t verificador;

    calcPolinomio_lento(p, a, &fa, NULL);
    calcPolinomio_lento(p, b, &fb, NULL);
    *it = 0;

    if (fa*fb > 0) {
        fprintf(stderr, "Erro: f(a) e f(b) têm mesmo sinal\n");
        return (real_t) -1;
    } else {
        do {
            xm_old = xm_new;
            xm_new = (a + b)/2;
            calcPolinomio_lento(p, xm_new, &fxm, NULL);
            verificador = fa*fxm;
            if (verificador < 0) {
                b = xm_new;
            } else if (verificador > 0) {
                a = xm_new;
                fa = fxm;
            } else {
                erro = 0.0;
                break;
            }
            erro = fabs(xm_new - xm_old)/fabs(xm_new);
            (*it)++;
        } while (erro > criterioParada && (*it) < MAXIT);
    }
    *raiz = xm_new;
    return erro;
}

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    if (p.grau < 0) {
        fprintf(stderr, "Erro: Polinômio de grau negativo\n");
        return;
    }

    // Caso especial para polinômios de grau 0
    *px = p.p[0];
    if (p.grau == 0) {
        if (dpx) *dpx = 0.0;
        return;
    }

    // Aplicação do método de Horner para cálculo eficiente
    real_t soma = p.p[p.grau];
    real_t soma_derivada = (dpx) ? p.p[p.grau] * p.grau : 0.0;

    for (int i = p.grau - 1; i >= 0; --i) {
        soma = soma * x + p.p[i];
        if (dpx && i > 0) {
            soma_derivada = soma_derivada * x + p.p[i] * i;
        }
    }

    *px = soma;
    if (dpx) *dpx = soma_derivada;
}

void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx)
{
    *px = 0.0;
    if (dpx) *dpx = 0.0;

    real_t pot_x = 1.0;  // Inicializa x^i com x^0 = 1

    for (int i = 0; i <= p.grau; ++i) {
        *px += p.p[i] * pot_x;
        if (dpx && i > 0) {
            *dpx += p.p[i] * i * (pot_x / x);  // Usa divisão para evitar recalcular potências
        }
        pot_x *= x;  // Atualiza x^i iterativamente
    }
}
