#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>
#include "simpi.h"
#define MATRIX_DIMENSION_XY 10
// double get_algbera(int pos){return arr[pos];}
// void set(int pos, int val){arr[pos] = val;}
// int get_size(){return dim;}
// void set(int pos,  double val){arr[pos] = val;}
// using namespace std;
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

matrix matrixtranspose(matrix A, matrix C, int par_id, int par_count)
{
  int Arow = A.get_x();
  int Acol = A.get_y();

  int rows = Arow;
  int rpp = rows / par_count;
  int start = rpp * par_id;
  int end = start + rpp;
  for (int i = start; i < end; i++) {
    for (int j = 0; j < Acol; j++) {
      // if(A[j*Arow + i] == 4){
      //     printf("LESSS GOOO = %d\n", j+i*Arow);
      // }
      C.set(j + i * Acol, A.get_algbera(j * Arow + i));
    }
  }
  return C;
}

matrix matrixadd(matrix A, matrix B, matrix C, int par_id, int par_count)
{
  // printf("YESS\n");
  printf("YEs\n");
  int Arow = A.get_x();
  int Acol = A.get_y();
  int Brow = B.get_x();
  int Bcol = B.get_y();
  // printf("YEs\n");

  if (Arow != Brow || Acol != Bcol) {
    // raise error
    return C;
  }
  int rows = Arow;
  int rpp = rows / par_count;
  int start = rpp * par_id;
  int end = start + rpp;
  for (int i = start; i < end; i++) {
    for (int j = 0; j < Acol; j++) {
      // if(A[j*Arow + i] == 4){
      //     printf("LESSS GOOO = %d\n", j+i*Arow);
      // }
      // printf("YESS\n");
      C.set(j + i * Acol,
            (B.get_algbera(j + i * Acol) + A.get_algbera(j + i * Acol)));
      // printf("YE\n");
    }
  }
  // printf("End of alg\n");
  return C;
}

matrix matrixsubtract(matrix A, matrix B, matrix C, int par_id, int par_count)
{
  int Arow = A.get_x();
  int Acol = A.get_y();
  int Brow = B.get_x();
  int Bcol = B.get_y();
  if (Arow != Brow || Acol != Bcol) {
    return C;
  }
  int rows = Arow;
  int rpp = rows / par_count;
  int start = rpp * par_id;
  int end = start + rpp;
  for (int i = start; i < end; i++) {
    for (int j = 0; j < Acol; j++) {
      // if(A[j*Arow + i] == 4){
      //     printf("LESSS GOOO = %d\n", j+i*Arow);
      // }
      C.set(j + i * Acol,
            (A.get_algbera(j + i * Acol) - B.get_algbera(j + i * Acol)));
    }
  }
  return C;
}

matrix matrixmul(matrix A, matrix B, matrix C, int par_id, int par_count)
{
  int last_ndx = A.get_x() * A.get_y() - 1;

  int Arow = A.get_x();
  int Acol = A.get_y();
  int Brow = B.get_x();
  int Bcol = B.get_y();
  int rpp = Bcol / par_count;
  int start = rpp * par_id;
  int end = start + rpp;

  if (Acol != Brow) {
    // Send error
    printf("error");
  }
  std::cout << "B:" << &B << std::endl;
  std::cout << "B.arr:" << B.arr << std::endl;
  std::cout << "A:" << &A << std::endl;
  std::cout << "A.arr:" << A.arr << std::endl;
  std::cout << "C:" << &C << std::endl;
  std::cout << "C.arr:" << C.arr << std::endl;
  for (int a = start; a < end; a++) {
    for (int b = 0; b < Arow; b++) {
      int sum = 0;
      for (int c = 0; c < Brow; c++) {
        // new_matrix.arr[a+b*Bcol] +=
        std::cout << "test\n" << std::endl;
        std::cout << sum << std::endl;
        std::cout << c + a * Brow << std::endl;
        std::cout << B.get_algbera(c + a * Brow) << std::endl;
        if (c + a * Brow > last_ndx || c * Arow + b > last_ndx ||
            a * Arow + b) {
          printf("out of bounds!");
          exit(1);
        }
        sum = sum + B.get_algbera(c + a * Brow) * A.get_algbera(c * Arow + b);
        std::cout << sum << std::endl;
      }
      C.set(a * Arow + b, sum);
      // new_matrix.arr[a+b*Bcol] +=
    }
  }
  return C;
}

matrix scaler_matrixmult(matrix A, int scaler, int par_id, int par_count)
{
  int size = A.get_x();
  int rpp = size / par_count;
  int start = rpp * par_id;
  int end = start + rpp;
  for (int i = start; i < end; i++) {
    for (int j = 0; j < A.get_y(); j++) {
      int pos = (A.get_y() * i + j);
      A.set(pos, A.get_algbera(pos) * scaler);
    }
  }
  return A;
}

vector scaler_vectormult(vector A, int scaler, int par_id, int par_count)
{
  int size = A.get_size();
  int rpp = size / par_count;
  int start = rpp * par_id;
  int end = start + rpp;
  for (int i = start; i < end; i++) {
    A.set(i, A.get(i) * scaler);
  }
  return A;
}
void get_print_matrixval(matrix A)
{
  printf("%f\n", A.get(0, 0));
}
int dot_product(vector A, vector B, int par_id, int par_count)
{
  int size = A.get_size();
  if (par_count > 3) {
    par_count = 3;
  }
  int rpp = size / par_count;
  int start = rpp * par_id;
  int end = start + rpp;
  int sum = 0;
  for (int i = start; i < end; i++) {
    sum += A.get(i) * B.get(0);
  }
  return sum;
}

// int cross_product(vector A, vector B, int par_id, int par_count)
// {
//     if ((A.get_size() == 2) && (B.get_size() == 2))
//         {
//             //return determinant
//             return 0;
//         }
//     if ((A.get_size() == 3) && (B.get_size() == 3))
//     {
//         if (par_count > 3){
//             par_count = 3;
//         }
//         int rpp = 3/par_count;
//         int start = rpp*par_id;
//         int end = start + rpp;
//         int sum = 0;
//         for (int i = start; i< end; i++)
//         {
//             int j = i;
//             if (i == 2)cout << "B:" << &B << endl;
// std::cout << "B.arr:" << B.arr<< std::endl;
// A.get(j + 1), B.get(i), B.get(j + 1)};
//             // sum += determinant(tempmatrix);
//         }

//     }
// }

// int main(int argc, char *argv[])
// {// void quadratic_matrix_print(double* C)
// // {
// //   printf("\n");
// //   for (int a = 0; a < MATRIX_DIMENSION_Y; a++) {
// //     printf("\n");
// //     for (int b = 0; b < MATRIX_DIMENSION_X; b++)
// //       printf("%.2f,", C[b + a * MATRIX_DIMENSION_X]);
// //   }
// //   printf("\n");
// // }
// // // int a = 0;int main(int argc, char *argv[])

// // // simpi* poin = &simp;
// // // int a = 0;
// // // long b = 100;
// // // ``
// // // printf("HERE\n");
// // // simpi* addy = &simp;
// // // matrix A(addy, 4, 4);
// // // printf("NOT HERE\n");
// // // matrix B(simp, 4, 4);
// // // matrix C(simp, 4, 4);
// //  //matrices A,B and C
// // // int *ready; //needed for synch

// // // A = (matrix)mmap(0, sizeof(matrix) * 16, PROT_READ | PROT_WRITE,
// MAP_SHARED | MAP_ANONYMOUS, -1, 0);
// // // B = (matrix)mmap(0, sizeof(matrix) * 16, PROT_READ | PROT_WRITE,
// MAP_SHARED | MAP_ANONYMOUS, -1, 0);
// // // C = (matrix)mmap(0, sizeof(matrix) * 16, PROT_READ | PROT_WRITE,
// MAP_SHARED | MAP_ANONYMOUS, -1, 0);

// // 	//TODO: initialize the matrices A and B
// // // int count = 0;
// // // for (int i = 0; i < 4; i++)
// // //     for (int j = 0; j< 4; j++)
// // //     {
// // //         A.set(i + j*2, count);
// // //         count++;
// // //     }
// // // count = 0;
// // // for (int i = 0; i < 4; i++)
// // // {
// // //     for (int j = 0; j< 4; j++)
// // //     {
// // //         B.set(i + j * 4, count);
// // //         count++;
// // //     }
// // // }
// // // printf("\n");
// // // quadratic_matrix_print(A);
// // // quadratic_matrix_print(B);

// // // matrixtranspose(A, C);
// // // matrixadd(A,B,C);

// // // for (int k = 0; k< 8; k++){
// // //     printf(" %f", A[k]);
// // // }

// // quadratic_matrix_print(C);
// // //TODO: quadratic_matrix_multiplication_parallel(par_id, par_count,A,B,C,
// ...);
// // //rows = ID;
// // //rpp = rows/count
// // //start = id * rpp
// //end = start * rpp

// // quadratic_matrix_multiplication(A,B,C);

// // if(par_id==0)
// //     quadratic_matrix_print(C);

// //lets test the result:

// return 0;
// }

// matrix B(simp, 4, 4);
// matrix C(simp, 4, 4);
// // matrices A,B and C
// // int *ready; //needed for synch

// // A = (matrix)mmap(0, sizeof(matrix) * 16, PROT_READ | PROT_WRITE,
// MAP_SHARED | MAP_ANONYMOUS, -1, 0);
// // B = (matrix)mmap(0, sizeof(matrix) * 16, PROT_READ | PROT_WRITE,
// MAP_SHARED | MAP_ANONYMOUS, -1, 0);
// // C = (matrix)mmap(0, sizeof(matrix) * 16, PROT_READ | PROT_WRITE,
// MAP_SHARED | MAP_ANONYMOUS, -1, 0);

// 	// TODO: initialize the matrices A and B
// int count = 0;
// for (int i = 0; i < 4; i++)
//     for (int j = 0; j< 4; j++)
//     {
//         A.set(i + j*2, count);
//         count++;
//     }
// count = 0;
// for (int i = 0; i < 4; i++)
// {
//     for (int j = 0; j< 4; j++)
//     {
//         B.set(i + j * 4, count);
//         count++;
//     }
// }
// printf("\n");
// quadratic_matrix_print(A);
// quadratic_matrix_print(B);

// matrixtranspose(A, C);
// matrixadd(A,B,C);

// for (int k = 0; k< 8; k++){
//     printf(" %f", A[k]);
// }

// quadratic_matrix_print(C);
// //TODO: quadratic_matrix_multiplication_parallel(par_id, par_count,A,B,C,
// ...); rows = ID; rpp = rows/count start = id * rpp end = start * rpp

// quadratic_matrix_multiplication(A,B,C);
// par_id
// lets test the result:

// return 0;
// }
