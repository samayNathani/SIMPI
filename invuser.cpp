#include <signal.h>
#include <string.h>

#include "simpi.h"
#define MATRIX_DIMENSION_X 4
#define MATRIX_DIMENSION_Y 4

int par_id;

void quadratic_matrix_print(double* C)
{
  printf("\n");
  for (int a = 0; a < MATRIX_DIMENSION_Y; a++) {
    printf("\n");
    for (int b = 0; b < MATRIX_DIMENSION_X; b++)
      printf("%.2f,", C[b + a * MATRIX_DIMENSION_X]);
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
  matrix A(my_simpi, MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  srand(time(NULL));
  for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
    for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
      A.get(x, y) = rand()%10 + 1;
    }
  }
  matrix B = A.inverse();
  if (par_id == 0) {
    printf("\nInput matrix: \n");
    quadratic_matrix_print(A.arr);
    printf("\nInverse matrix: \n");
    quadratic_matrix_print(B.arr);

  }
}