#include "simpi.h"
#define MATRIX_DIMENSION_XY 10

void set_matrix_elem(double* M, int x, int y, double f)
{
  M[x + y * MATRIX_DIMENSION_XY] = f;
}
void quadratic_matrix_print(double* C)
{
  printf("\n");
  for (int a = 0; a < MATRIX_DIMENSION_XY; a++) {
    printf("\n");
    for (int b = 0; b < MATRIX_DIMENSION_XY; b++)
      printf("%.2f,", C[a + b * MATRIX_DIMENSION_XY]);
  }
  printf("\n");
}

int main(int argc, char* argv[])
{
  printf("im alive\n");
  int par_id = atoi(argv[1]);
  printf("par_id %d\n", par_id);
  int synch_size = atoi(argv[2]);
  printf("synch_size %d\n", synch_size);
  simpi my_simpi(par_id, synch_size);
  printf("simpi worked\n");
  matrix A(my_simpi, 10, 10);
  printf("matrix worked\n");
  for (int x = 0; x < MATRIX_DIMENSION_XY; x++) {
    for (int y = 0; y < MATRIX_DIMENSION_XY; y++) {
      A.get(x, y) = (double)x * y;
    }
  }
  printf("matrix init worked\n");

  quadratic_matrix_print(A.arr);
  printf("matrix init worked\n");
}