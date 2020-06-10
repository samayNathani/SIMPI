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
  //  printf("YES I AM STUPID\n");
  //  printf("ProcessId = %d\n", par_id);
  //  printf("YES I AM STUPID\n");
  matrix B(3, 8);
  matrix C(8, 8);
  matrix E(8, 8);
  vector D(10);
  vector F(10);
  vector G(10);
  SIMPI_SYNCH();
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
   c = 0;
  for (int y = 0; y < D.get_size(); y++) {
    D.set(y, c++);
  }
   c = 10;
  for (int y = 0; y < F.get_size(); y++) {
    F.set(y, c++);
  }
  
  std::cout << D;
   std::cout << F;
 
  SIMPI_SYNCH();

  //C = A.inverse();
  // printf("YES I AM STUPID\n");
  SIMPI_SYNCH();
  int i = D * F;
  if (par_id == 0){
    printf("GOT HERE\n");
    printf("%d\n", i);
  }
  SIMPI_SYNCH();
  // std::cout << G;
  SIMPI_SYNCH();
  // std::cout << A;
  // std::cout << B;
  // C= A.inverse();

  // std::cout << A; 
  //SIMPI_FINALIZE();
}