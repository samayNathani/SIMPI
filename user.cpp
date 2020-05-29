#include <signal.h>
#include <string.h>

#include "simpi.h"
#define MATRIX_DIMENSION_X 2
#define MATRIX_DIMENSION_Y 3

int par_id;

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
  int num_workers = atoi(argv[2]);
  SIMPI_INIT(par_id, num_workers);
  matrix A(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix C(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  vector D(10);
  SIMPI_SYNCH();

  for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
    for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
      A.get(x, y) = (double)x + y;
    }
  }
  SIMPI_SYNCH();

  //C = A.inverse();
  std::cout << A;
  
  std::cout << A;

  //C= A.inverse();

  //std::cout << C; 
  //SIMPI_FINALIZE();
}