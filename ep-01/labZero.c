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

  //============================ Calculo de Polinomio Lento ============================

  printf("\n\nLENTO\n\n");

  // Teste bisseccao
  int it;
  real_t raiz;
  
  //criterio de parada 1
  rtime_t t0 = timestamp();
  real_t erro = bisseccao_lento(pol, a, b, 1, &it, &raiz);
  rtime_t t1 = timestamp() - t0;
  printf("bissec: %.15e %.15e %d %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = bisseccao_lento(pol, a, b, 2, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: %.15e %.15e %d %.8e\n", raiz, erro, it, t1);

  //criterio de parada 3
  t0 = timestamp();
  erro = bisseccao_lento(pol, a, b, 3, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: %.15e %.15e %d %.8e\n", raiz, erro, it, t1);

  /* Teste newtonRaphson

  //criterio de parada 1
  t0 = timestamp();
  erro = newtonRaphson_lento(pol, a, 1, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = newtonRaphson_lento(pol, a, 2, &it, &raiz);
  t1 = timestamp() - t0;
  printf("newton: raiz = %.15e, erro = %.15e, it = %d, tempo = %.8e\n", raiz, erro, it, t1);

  */
  //============================ Calculo de Polinomio Rápido ============================


  printf("\n\nRAPIDO\n\n");
  // Teste bisseccao
  
  //criterio de parada 1
  t0 = timestamp();
  erro = bisseccao(pol, a, b, 1, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: %.15e %.15e %d %.8e\n", raiz, erro, it, t1);

  //criterio de parada 2
  t0 = timestamp();
  erro = bisseccao(pol, a, b, 2, &it, &raiz);
  t1 = timestamp() - t0;
  printf("bissec: %.15e %.15e %d %.8e\n", raiz, erro, it, t1);

   //criterio de parada 3
   t0 = timestamp();
   erro = bisseccao(pol, a, b, 3, &it, &raiz);
   t1 = timestamp() - t0;
   printf("bissec: %.15e %.15e %d %.8e\n", raiz, erro, it, t1);

  /* Teste newtonRaphson
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
  */
  free(pol.p);
  return 0;
}
