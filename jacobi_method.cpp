#include <iostream>
#include "matrix.h"
#include "vector.h"

vector* jacobi_method(matrix* equations, vector* constants, int processCount, int id);

int main()
{
    std::cout << "Sample Test" << std::endl;
    matrix* eq = new matrix(2,2);
    vector* constants = new vector(2);
    eq->get(0,0) = 2.0;
    eq->get(0,1) = 1.0;
    eq->get(1,0) = 2.0;
    eq->get(1,1) = 3.0;
    constants->get(0) = 5.0;
    constants->get(1) = 7.0;
    vector* solution = jacobi_method(eq, constants, 1, 0);

    std::cout << "x: " << solution->get(0) << std::endl;
    std::cout << "y: " << solution->get(1) << std::endl;

    return 0;

}

vector* jacobi_method(matrix* equations, vector* constants, int processCount, int id) {
    vector *prev = new vector(constants->getDim()); // shared mem containing a copy of values
    vector *solution = new vector(constants->getDim()); // shared mem containing actual calculated values

    int n = constants->getDim();

    int work = n / processCount;
    int i, j, k;
    int start = id * n;
    int end = start + work;

    //setup
    for (i = start; i < end; i++) {
        int temp = equations->get(i, i);
        equations->get(i, i) = constants->get(i);
        constants->get(i) = temp;
        for (j = 0; j < equations->gety(); j++) {
            if (j != i) {
                equations->get(i, j) *= -1;
            }
            equations->get(i, j) /= constants->get(i);
        }
        prev->get(i) = 0;
        solution->get(i) = constants->get(i);
    }

    // first iteration with all 0s
    //lock mutex
    for (i = start; i < end; i++) {
        float rowSum = 0;
        for (j = 0; j < equations->gety(); j++) {
            if (j == i) {
                rowSum += equations->get(i, j);
            } else {
                rowSum += (equations->get(i, j) * prev->get(j));
            }
        }
        solution->get(i) = rowSum;
    }
    //unlock mutex

    //wait for all processes here

    //lock mutex
    //loop until calculated difference is extremely small or upper bound is reached ** consult prof eckhardt
    for (k = 0; k < 100; k++)
    {
        //lock mutex
        for (i = start; i < end; i++)
        {
            //save prev value for comparision
            prev->get(i) = solution->get(i);
            float rowSum = 0;
            for (j = 0; j < equations->gety(); j++)
            {
                if (j == i) {
                    rowSum += equations->get(i, j);
                } else {
                    rowSum += (equations->get(i, j) * prev->get(j));
                }
            }
            solution->get(i) = rowSum;
        }
        //unlock mutex
        //wait for all processes
    }
    //unlock mutex

    return solution;

}