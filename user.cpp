#include <signal.h>
#include <string.h>

#include "simpi.h"
#define MATRIX_DIMENSION_XY 10

int par_id;

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

void segfault_printer(int dummy)
{
  char buf[20];
  sprintf(buf, "%d: segfaulted\n", par_id);
  write(STDOUT_FILENO, buf, strlen(buf));
  exit(1);
}

int main(int argc, char* argv[])
{
  par_id = atoi(argv[1]);
  int synch_size = atoi(argv[2]);
  simpi my_simpi(par_id, synch_size);
  matrix A(my_simpi, 10, 10);

  for (int x = 0; x < MATRIX_DIMENSION_XY; x++) {
    for (int y = 0; y < MATRIX_DIMENSION_XY; y++) {
      A.get(x, y) = (double)x * y;
    }
  }
  quadratic_matrix_print(A.arr);
}