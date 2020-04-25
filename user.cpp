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
  int par_id = atoi(argv[1]);
  printf("par_id %d\n", par_id);
  int synch_size = atoi(argv[2]);
  simpi my_simpi(par_id, synch_size);
  printf("%d: simpi initialized\n", par_id);
  matrix A(my_simpi, 10, 10);
  printf("%d: matrix initialized\n", par_id);
  for (int x = 0; x < MATRIX_DIMENSION_XY; x++) {
    for (int y = 0; y < MATRIX_DIMENSION_XY; y++) {
      A.get(x, y) = (double)x * y;
    }
  }
  printf("%d: matrix values initialized\n", par_id);

  quadratic_matrix_print(A.arr);
}