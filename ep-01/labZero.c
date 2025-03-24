#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "ZeroFuncao.h"

int main ()
{
  real_t a, b;
  Polinomio pol;

  scanf("%d", &pol.grau);

  pol.p = (real_t *) malloc((pol.grau+1)*sizeof(real_t));
  if (pol.p == NULL) {
    fprintf(stderr, "Erro de alocacao de memoria\n");
    return 1;
  }

  for (int i = pol.grau; i >= 0; --i)
    scanf("%lf", &pol.p[i]);

  scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.

  //============================ Calculo de Polinomio Rápido ============================

  printf("RAPIDO\n");

  // Teste bisseccao
  int it;
  real_t raiz;
  
  //criterio de parada 1
  rtime_t t0 = timestamp();
  real_t erro = bisseccao(pol, a, b, EPS, &it, &raiz);
  rtime_t t1 = timestamp() - t0;
  printf("bissec: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = bisseccao(pol, a, b, ULPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 3
  t0 = timestamp();
  erro = bisseccao(pol, a, b, ZERO, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  // Teste newtonRaphson, chute inicial = a
  //criterio de parada 1
  t0 = timestamp();
  erro = newtonRaphson(pol, a, EPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = newtonRaphson(pol, a, ULPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 3
  t0 = timestamp();
  erro = newtonRaphson(pol, a, ZERO, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);
  
  //============================ Calculo de Polinomio Lento ============================


  printf("LENTO\n");
  // Teste bisseccao
  
  //criterio de parada 1
  t0 = timestamp();
  erro = bisseccao_lento(pol, a, b, EPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = bisseccao_lento(pol, a, b, ULPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 3
  t0 = timestamp();
  erro = bisseccao_lento(pol, a, b, ZERO, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  // Teste newtonRaphson, chute inicial = a
  //criterio de parada 1
  t0 = timestamp();
  erro = newtonRaphson_lento(pol, a, EPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = newtonRaphson_lento(pol, a, ULPS, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 3
  t0 = timestamp();
  erro = newtonRaphson_lento(pol, a, ZERO, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  free(pol.p);
  return 0;
}
