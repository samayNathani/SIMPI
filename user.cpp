#include <signal.h>
#include <string.h>

#include "simpi.h"
#define MATRIX_DIMENSION_X 4
#define MATRIX_DIMENSION_Y 4

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
  matrix upper(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  matrix lower(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);



  //matrix C(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  // vector D(10);
  SIMPI_SYNCH();

  for (int y = 0; y < MATRIX_DIMENSION_Y; y++) {
    for (int x = 0; x < MATRIX_DIMENSION_X; x++) {
      A.get(y,x) = rand()%10 + 1;
    }
  }
  // A.get(0,0) = 1;
  // A.get(0,1) = 0;
  // A.get(0,2) = 0;
  // A.get(0,3) = 0;
  // A.get(1,0) = 0.125;
  // A.get(1,1) = 1;
  // A.get(1,2) = 0;
  // A.get(1,3) = 0;
  // A.get(3,0) = 0.375;
  // A.get(3,1) = -0.428571;
  // A.get(3,2) = -0.608247;
  // A.get(3,3) = 1;
  // A.get(2,0) = 0.5;
  // A.get(2,1) = 2.85714;
  // A.get(2,2) = 1;
  // A.get(2,3) = 0;

  SIMPI_SYNCH();
  //A.luDecomposition(&upper,&lower);
  std::cout<<A;
  matrix inv(MATRIX_DIMENSION_X, MATRIX_DIMENSION_Y);
  A.inverseLU(&inv);
  SIMPI_SYNCH();
  std::cout<<inv;
  SIMPI_SYNCH();



  //matrix C = A.inverse();

  // std::cout << "Input is: " << '\n' << A;

  // std::cout << "Lower is: " << '\n' << lower;

  // std::cout << "Upper is: " << '\n' << upper;

  // std::cout << A;
  // std::cout << lower;
  // std::cout << upper;



  //C= A.inverse();

  //std::cout << C; 
  SIMPI_FINALIZE();
}