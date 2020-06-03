#include <stdio.h>
#include "simpi.h"

matrix &multiply(matrix A, matrix B)
{
    matrix *result = new matrix(A.get_x(), B.get_y());
    int number_of_processes = A.get_simpi()->get_synch_info()->par_count;
    int parId = A.get_simpi()->get_id();

    if (number_of_processes > A.get_x())
    {
        number_of_processes = A.get_x();
    }
    int Arow = A.get_x();
    int Acol = A.get_y();
    int Brow = B.get_x();
    int rpp = Arow / number_of_processes;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;
    if (A.get_x() % number_of_processes != 0)
    {
        int leftover = Arow % number_of_processes;
        if (parId <= leftover)
        {
            parId += (Arow - leftover);
            int start = parId;
            int end = start + 1;
            for (int a = start; a < end; a++)
            {
                for (int b = 0; b < Arow; b++)
                {
                    int sum = 0;
                    for (int c = 0; c < Brow; c++)
                    {

                        sum = sum + B.get_algbera(c + a * Brow) * A.get_algbera(c * Arow + b);
                    }
                    result->set(a * Arow + b, sum);
                }
            }
        }
    }
    if (Acol != Brow)
    {
        // Send error
        printf("error");
    }

    for (int a = start; a < end; a++)
    {
        for (int b = 0; b < Arow; b++)
        {
            int sum = 0;
            for (int c = 0; c < Brow; c++)
            {
                sum = sum + B.get_algbera(c + a * Brow) * A.get_algbera(c * Arow + b);
            }
            result->set(a * Arow + b, sum);
        }
    }
    A.get_simpi()->synch();
    return *result;
}

matrix &transpose(matrix A)
{
    matrix *result = new matrix(A.get_y(), A.get_x());
    int number_of_processes = A.get_simpi()->get_synch_info()->par_count;

    if (number_of_processes > A.get_x())
    {
        number_of_processes = A.get_x();
    }

    int Arow = A.get_x();
    int Acol = A.get_y();

    int rpp = Arow / number_of_processes;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Acol; j++)
        {
            // if(A[j*Arow + i] == 4){
            //     printf("LESSS GOOO = %d\n", j+i*Arow);
            // }
            result->set(j + i * Acol, A.get_algbera(j * Arow + i));
        }
    }
    return *result;
}

matrix &add(matrix A, matrix B)
{
    matrix *result = new matrix(A.get_x(), A.get_y());

    int number_of_processes = A.get_simpi()->get_synch_info()->par_count;
    int parId = A.get_simpi()->get_id();

    if (number_of_processes > A.get_x())
    {
        number_of_processes = A.get_x();
    }

    int Arow = A.get_x();
    int Acol = A.get_y();
    int Brow = B.get_x();
    int Bcol = B.get_y();
    int rpp = Arow / number_of_processes;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;

    if (Arow != Brow || Acol != Bcol)
    {
        // raise error
    }
    if (A.get_x() % number_of_processes != 0)
    {
        int leftover = Arow % number_of_processes;
        if (parId <= leftover)
        {
            parId += (Arow - leftover);
            int start = parId;
            int end = start + 1;
            for (int i = start; i < end; i++)
            {
                for (int j = 0; j < Arow; j++)
                {
                    result->set(j + i * Arow,
                                (B.get_algbera(j + i * Arow) + A.get_algbera(j + i * Arow)));
                }
            }
        }
    }
    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Arow; j++)
        {
            result->set(j + i * Arow,
                        (B.get_algbera(j + i * Arow) + A.get_algbera(j + i * Arow)));
        }
    }
    // printf("End of alg\n");
    return *result;
}

matrix& subtract(matrix A, matrix B)
{
  matrix* result = new matrix(A.get_x(), A.get_y());

  int number_of_processes = A.get_simpi()->get_synch_info()->par_count;
  int parId = A.get_simpi()->get_id();
  
  if (number_of_processes > A.get_x()){
    number_of_processes = A.get_x();
  }

  int Arow = A.get_x();
//   int Acol = A.get_y();
//   int Brow = B.get_x();
//   int Bcol = B.get_y();
  int rpp = Arow / number_of_processes;
  int start = rpp * A.get_simpi()->get_id();
  int end = start + rpp;
  

	if (A.get_x() % number_of_processes != 0){
      int leftover = Arow % number_of_processes;
      if (parId <= leftover){
		            parId += (Arow - leftover); 
                int start = parId;
                int end = start + 1;
                for (int i = start; i < end; i++) {
                  for (int j = 0; j < Arow; j++) {
					result->set(j + i * Arow,
           				 (A.get_algbera(j + i * Arow) - B.get_algbera(j + i *Arow)));
                  }
                }              
                    
      }
  }
  for (int i = start; i < end; i++) {
    for (int j = 0; j < Arow; j++) {
      result->set(j + i * Arow,
            (A.get_algbera(j + i * Arow) - B.get_algbera(j + i * Arow)));
      
    }
  }
  return *result;
}

matrix &scalar_matrix_mult(matrix A, int scaler)
{
    matrix *result = new matrix(A.get_x(), A.get_y());

    int size = A.get_x();
    int number_of_processes = A.get_simpi()->get_synch_info()->par_count;
    // int parId = A.get_simpi()->get_id();

    if (number_of_processes > A.get_x())
    {
        number_of_processes = A.get_x();
    }

    int rpp = size / number_of_processes;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < A.get_y(); j++)
        {
            int pos = (A.get_y() * i + j);
            result->set(pos, A.get_algbera(pos) * scaler);
        }
    }
    return *result;
}

bool matrix_is_equal(matrix A, matrix B)
{
    int number_of_processes = A.get_simpi()->get_synch_info()->par_count;
    // int parId = A.get_simpi()->get_id();

    if (number_of_processes > A.get_x())
    {
        number_of_processes = A.get_x();
    }

    int Arow = A.get_x();
    int Acol = A.get_y();
    int Brow = B.get_x();
    int Bcol = B.get_y();
    int rpp = Arow / number_of_processes;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;

    if (Arow != Brow || Acol != Bcol)
    {
        // error
    }

    for (int i = start; i < end; i++)
    {
        for (int j = 0; j < Acol; j++)
        {

            if (A.get_algbera(j + i * Acol) != B.get_algbera(j + i * Acol))
            {
                return false;
            }
        }
    }
    return true;
}

vector &scalar_vector_mult(vector A, int scaler)
{
    vector *result = new vector(A.get_size());
    int size = A.get_size();
    int rpp = size / A.get_simpi()->get_synch_info()->par_count;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        result->set(i, A.get(i) * scaler);
    }
    return *result;
}

vector &add(vector A, vector B)
{
    vector *result = new vector(A.get_size());
    int size = A.get_size();
    int rpp = size / A.get_simpi()->get_synch_info()->par_count;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        result->set(i, A.get(i) + B.get(i));
    }
    return *result;
}

vector &subtract(vector A, vector B)
{
    vector *result = new vector( A.get_size());
    int size = A.get_size();
    int rpp = size / A.get_simpi()->get_synch_info()->par_count;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        result->set(i, A.get(i) - B.get(i));
    }
    return *result;
}

bool vector_is_equal(vector A, vector B)
{
    int size = A.get_size();
    int rpp = size / A.get_simpi()->get_synch_info()->par_count;
    int start = rpp * A.get_simpi()->get_id();
    int end = start + rpp;
    for (int i = start; i < end; i++)
    {
        if (A.get(i) != B.get(i))
        {
            return false;
        }
    }
    return true;
}

//Multiplication

matrix &operator*(matrix &lhs, matrix &rhs)
{
    return multiply(lhs, rhs);
}
matrix &operator*(matrix &lhs, int rhs)
{
    return scalar_matrix_mult(lhs, rhs);
}
matrix &operator*(int lhs, matrix &rhs)
{
    return scalar_matrix_mult(rhs, lhs);
}
vector &operator*(vector &lhs, int rhs)
{
    return scalar_vector_mult(lhs, rhs);
}
vector &operator*(int lhs, vector &rhs)
{
    return scalar_vector_mult(rhs, lhs);
}

//*=
void operator*=(matrix &lhs, matrix &rhs)
{
    lhs = multiply(lhs, rhs);
}
void operator*=(matrix &lhs, int rhs)
{
    lhs = scalar_matrix_mult(lhs, rhs);
}
void operator*=(int lhs, matrix &rhs)
{
    rhs = scalar_matrix_mult(rhs, lhs);
}
void operator*=(int lhs, vector &rhs)
{
    rhs = scalar_vector_mult(rhs, lhs);
}

//Adition
matrix &operator+(matrix &lhs, matrix &rhs)
{
    return add(lhs, rhs);
}
vector &operator+(vector &lhs, vector &rhs)
{
    return add(lhs, rhs);
}
void operator+=(matrix &lhs, matrix &rhs)
{
    lhs = add(lhs, rhs);
}
void operator+=(vector &lhs, vector &rhs)
{
    lhs = add(lhs, rhs);
}
void operator-=(matrix &lhs, matrix &rhs)
{
    lhs = subtract(lhs, rhs);
}
void operator-=(vector &lhs, vector &rhs)
{
    lhs = subtract(lhs, rhs);
}
matrix &operator-(matrix &lhs, matrix &rhs)
{
    return subtract(lhs, rhs);
}
vector &operator-(vector &lhs, vector &rhs)
{
    return subtract(lhs, rhs);
}
bool operator==(matrix &lhs, matrix &rhs)
{
    return matrix_is_equal(lhs, rhs);
}
bool operator==(vector &lhs, vector &rhs)
{
    return vector_is_equal(lhs, rhs);
}
