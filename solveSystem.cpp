#include <iostream>
#include "matrix.h"
#include "vector.h"

vector* solveSystem(matrix* equations, vector* constants, int processCount, int id);
int test();

int main()
{
    return test();
}
/*
 * Solves a system of linear equations
 * matrix*<float>[n][n], vector*<float>[n][n] -> vector*<float>[n][n]
 * No error checking yet! consult if needed
 */
vector* solveSystem(matrix* equations, vector* constants, int processCount, int id) {

    //TODO: initialize shared memory if id is 0
    vector *prev = new vector(constants->getDim()); // shared mem containing a copy of values
    vector *solution = new vector(constants->getDim()); // shared mem containing actual calculated values
    matrix* saveEq = new matrix(equations->getx(),equations->gety());
    vector* saveConst = new vector(constants->getDim());

    int n = constants->getDim();

    int work = n / processCount;
    /*
     * implement remainder - multi-core processing
     */
    int i, j, k;
    int start = id * n;
    int end = start + work;

    //Save Matrix and Vector
    for (i = start; i < end; i++) {
        for (j = 0; j < equations->gety(); j++) {
            saveEq->get(i,j) = equations->get(i, j);
        }
        saveConst->get(i) = constants->get(i);
    }

    //wait for all processes

    //setup
    for (i = start; i < end; i++) {
        float temp = equations->get(i, i);
        equations->get(i, i) = constants->get(i);
        constants->get(i) = temp;
        for (j = 0; j < equations->gety(); j++) {
            if (j != i) {
                equations->get(i, j) *= -1;
            }
            equations->get(i, j) /= constants->get(i);
        }
        prev->get(i) = 1;
        solution->get(i) = constants->get(i);
    }

    // first iteration with all 0s
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

    //TODO: synchObject()
    //wait for all processes here

    for (k = 0; k < 100; k++)
    {
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
        //TODO: synchObject()
        //wait for all processes
    }

    //restore original matrix and vector
    for (i = start; i < end; i++) {
        for (j = 0; j < equations->gety(); j++) {
            equations->get(i, j) = saveEq->get(i,j);
        }
        constants->get(i) = saveConst->get(i);
    }
    //wait for all processes

    if (id == 0)
    {
        return solution;
    }
    return nullptr;
}

/*
 * 3 X 3 System of Linear Equations example, x = 1, y = 2, z = 3
 */
int test()
{
    std::cout << "Sample Test" << std::endl;
    matrix* eq = new matrix(2,2);
    vector* constants = new vector(2);
    eq->get(0,0) = 9.0;
    eq->get(0,1) = 3.0;
    eq->get(1,0) = 3.0;
    eq->get(1,1) = 1.0;
    constants->get(0) = 21.0;
    constants->get(1) = 7.0;
    vector* solution = solveSystem(eq, constants, 1, 0);

    std::cout << "x: " << solution->get(0) << std::endl;
    std::cout << "y: " << solution->get(1) << std::endl;

    return 0;
}