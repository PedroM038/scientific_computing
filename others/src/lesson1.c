#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "../includes/lesson1.h"

void bhaskara(float a, float b, float c, float *x1, float *x2)
{
    float delta = b * b - 4 * a * c;
    if (delta < 0)
    {
        printf("Nao existem raizes reais\n");
        return;
    }
    *x1 = (-b + sqrt(delta)) / (2 * a);
    *x2 = (-b - sqrt(delta)) / (2 * a);
}

double bissecao(double xl, double xu, double f(double x), double epsilon)
{
    double xm_old, xm_new;

    xm_new = (xl + xu) / 2;
    if (f(xl)*f(xm_new) < 0)
        xu = xm_new;
    else if (f(xl)*f(xm_new) > 0)
        xl = xm_new;
    else
        return xm_new;

    do{
        xm_old = xm_new;
        xm_new = (xl + xu) / 2;
        if (f(xl)*f(xm_new) < 0)
            xu = xm_new;
        else if (f(xl)*f(xm_new) > 0)
            xl = xm_new;
        else
            return xm_new;
    } while (fabs((xm_new - xm_old) / xm_new)*100 > epsilon);

    return xm_new;
}

double iteracao(double x0, double f(double x), double phi(double x),
                double epsilon_x, double epsilon_f, int max_iter)
{
    double x_new = x0, x_old;
    int criterio1, criterio2, criterio3;
    int iter = 1;

    do{
        x_old = x_new;
        x_new = phi(x_old);
        criterio1 = (fabs(x_old - x_new) < epsilon_x);
        criterio2 = (fabs(f(x_new)) < epsilon_f);
        criterio3 = (iter == max_iter);
        ++iter;
    } while (!criterio1 && !criterio2 && !criterio3);

    return x_new;
}

void polinomio_e_derivada(double *a, int n, double x, double *Px, double *Dx) {
    double b = 0;
    double c = 0;
    for (int i = 0; i < n; --i){
        b = b*x + a[i];
        c = c*x + b;
    }
    b = b*x + a[0];
    *Px = b;
    *Dx = c;
}