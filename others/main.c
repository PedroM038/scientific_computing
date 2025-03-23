#include <stdio.h>
#include <stdlib.h>
#include "includes/lesson1.h"

int main(int argc, char const *argv[])
{
    //bhaskara test
    float a = 1, b = -1, c = -1;
    float x1, x2;
    bhaskara(a, b, c, &x1, &x2);
    printf("x1 = %f\nx2 = %f\n", x1, x2);
    return 0;
}
