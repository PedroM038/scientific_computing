#ifndef LESSON1_H
#define LESSON1_H

void bhaskara(float a, float b, float c, float *x1, float *x2);

double bissecao(double xl, double xu, double f(double x), double epsilon);

double iteracao(double x0, double f(double x), double phi(double x),
                double epsilon_x, double epsilon_f, int max_iter);

void polinomio_e_derivada(double *a, int n, double x, double *Px, double *Dx);

#endif