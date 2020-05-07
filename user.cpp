#include <signal.h>
#include <string.h>

#include "alg.cpp"
#define MATRIX_DIMENSION_X 2
#define MATRIX_DIMENSION_Y 3

int par_id;

// void quadratic_matrix_print(double* C)
// {
//   printf("\n");
//   for (int a = 0; a < MATRIX_DIMENSION_Y; a++) {
//     printf("\n");
//     for (int b = 0; b < MATRIX_DIMENSION_X; b++)
//       printf("%.2f,", C[b + a * MATRIX_DIMENSION_X]);
//   }
//   printf("\n");
// }

void segfault_printer(int dummy)
{
  char buf[20];
  sprintf(buf, "%d: segfaulted\n", par_id);
  write(STDOUT_FILENO, buf, strlen(buf));
  exit(1);
}

int main(int argc, char* argv[])
{
  signal(SIGSEGV, segfault_printer);
  par_id = atoi(argv[1]);
  int synch_size = atoi(argv[2]);
  simpi my_simpi(par_id, synch_size);
  matrix A(my_simpi, 4, 4);
  matrix B(my_simpi, 4, 4);
  matrix C(my_simpi, 4, 4);

  my_simpi.synch();

  for (int y = 0; y < B.get_y(); y++) {
    for (int x = 0; x < B.get_x(); x++) {
      B.set(x + y * B.get_x(), 3);
    }
  }

  for (int y = 0; y < A.get_y(); y++) {
    for (int x = 0; x < A.get_x(); x++) {
      A.set(x + y * A.get_x(), 3);
    }
  }

  // double num = A.get(0, 0);
  // printf("%f\n", num);

  // double num2 = A.get(0, 0);
  // printf("%f\n", num2);

  get_print_matrixval(A, 1, 2);
  my_simpi.synch();

  get_print_matrixval(A, 0, 2);

  // quadratic_matrix_print(B);
  // quadratic_matrix_print(B);
  // quadratic_matrix_print(A);
  // quadratic_matrix_print(C);
  // matrixmul(A,B,C,0,1);
  // quadratic_matrix_print(C);
}