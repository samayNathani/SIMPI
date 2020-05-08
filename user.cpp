#include <signal.h>
#include <string.h>

// #include "alg.cpp"
#include "simpi.h"

// #define MATRIX_DIMENSION_X 2
// #define MATRIX_DIMENSION_Y 3

int par_id;

void quadratic_matrix_print(matrix C)
{
  printf("\n");
  for (int a = 0; a < C.get_x(); a++) {
    printf("\n");
    for (int b = 0; b < C.get_y(); b++)
      printf("%.2f,", C.get_algbera(a + b * C.get_x()));
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
  signal(SIGSEGV, segfault_printer);
  par_id = atoi(argv[1]);
  int synch_size = atoi(argv[2]);
  simpi my_simpi(par_id, synch_size);
  matrix A(my_simpi, 4, 4);
  matrix B(my_simpi, 2, 4);
  matrix C(my_simpi, 4, 2);

  my_simpi.synch();
  int c = 100;
  for (int y = 0; y < A.get_y(); y++) {
    for (int x = 0; x < A.get_x(); x++) {
      A.set(x + y * A.get_x(), c++);
    }
  }


  c = 0;
  for (int y = 0; y < B.get_y(); y++) {
    for (int x = 0; x < B.get_x(); x++) {
      B.set(x + y * B.get_x(), c++);
    }
  }
 
  // for(int i = 0; i < 4; i++){
  C = B.transpose();

  // }
  quadratic_matrix_print(C);
  my_simpi.synch();
  // quadratic_matrix_print(A);
  // matrix C = matrixmul(A, B);
  // quadratic_matrix_print(C);
}