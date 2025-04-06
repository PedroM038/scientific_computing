#include <stdio.h>
#include <float.h>

#include "includes/utils.h"
#include "includes/ZeroFuncao.h"

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

  
  int it;
  real_t raiz;
  real_t erro;
  
  //============================ Calculo de Polinomio Rápido ============================

  printf("\nRAPIDO\n");
    
  // Teste bisseccao
  //criterio de parada 1
  rtime_t t0 = timestamp();
  erro = bisseccao(pol, a, b, 1, &it, &raiz, RAPIDO);
  rtime_t t0f = timestamp() - t0;
  printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, t0f);
  
  //criterio de parada 2
  rtime_t t1 = timestamp();
  erro = bisseccao(pol, a, b, 2, &it, &raiz, RAPIDO);
  rtime_t t1f = timestamp() - t1;
  printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, t1f);
  
  //criterio de parada 3
  rtime_t t2 = timestamp();
  erro = bisseccao(pol, a, b, 3, &it, &raiz, RAPIDO);
  rtime_t t2f = timestamp() - t2;
  printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, t2f);
  
  //Teste newtonRaphson
  //criterio de parada 1
  rtime_t t3 = timestamp();
  erro = newtonRaphson(pol, a, 1, &it, &raiz, RAPIDO);
  rtime_t t3f = timestamp() - t3;
  printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, t3f);
  
  //criterio de parada 2
  rtime_t t4 = timestamp();
  erro = newtonRaphson(pol, a, 2, &it, &raiz, RAPIDO);
  rtime_t t4f = timestamp() - t4;
  printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, t4f);
  
  //criterio de parada 3
  rtime_t t5 = timestamp();
  erro = newtonRaphson(pol, a, 3, &it, &raiz, RAPIDO);
  rtime_t t5f = timestamp() - t5;
  printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, t5f);
    
  //============================ Calculo de Polinomio Lento ============================

  printf("\nLENTO\n");
  
  //criterio de parada 1
  t0 = timestamp();
  erro = bisseccao(pol, a, b, 1, &it, &raiz, LENTO);
  t0f = timestamp() - t0;
  printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, t0f);

  //criterio de parada 2
  t1 = timestamp();
  erro = bisseccao(pol, a, b, 2, &it, &raiz, LENTO);
  t1f = timestamp() - t1;
  printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, t1f);

  //criterio de parada 3
  t2 = timestamp();
  erro = bisseccao(pol, a, b, 3, &it, &raiz, LENTO);
  t2f = timestamp() - t2;
  printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, t2f);

  // Teste newtonRaphson

  //Teste newtonRaphson
  //criterio de parada 1
  t3 = timestamp();
  erro = newtonRaphson(pol, a, 1, &it, &raiz, LENTO);
  t3f = timestamp() - t3;
  printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, t3f);

  //criterio de parada 2
  t4 = timestamp();
  erro = newtonRaphson(pol, a, 2, &it, &raiz, LENTO);
  t4f = timestamp() - t4;
  printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, t4f);

  //criterio de parada 3
  t5 = timestamp();
  erro = newtonRaphson(pol, a, 3, &it, &raiz, LENTO);
  t5f = timestamp() - t5;
  printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, t5f);
  
  free(pol.p);
  return 0;
}
