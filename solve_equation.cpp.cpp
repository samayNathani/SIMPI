#include <iostream>
#include "matrix.h"
#include "vector.h"

int main()
{
    std::cout << "Test work" << std::endl;
    return 0;

}

int solve_system(matrix* equations, vector* constants, vector* solution)
{
    //**********************************************************************************************
    //Simple method - Single process, depends on Inversion
    // equations * solution = constants
    // matrix inverted = inverse(equations)
    // solution = inverted * constants
    // DOT PRODUCT
    // store values in solution vector
    //**********************************************************************************************
    equations->invert(); // gap filler - chilam

    for (int i = 0; i < equations->getx(); i++)
    {
        float row_sum = 0;
        for (int j = 0; j < equations->gety(); j++)
        {
            row_sum += (equations->get(i,j) * constants->arr[j]);
        }
        solution->arr[i] = row_sum;
    }
    //**********************************************************************************************
    //Gaussian Elimination method - Single process, depends on the existence of row echelon form
    // equations | constants = solution (augmented matrix)
    // convert to reduced row echelon form, lower-left triangular form
    // solve with back-substitution
    // check for already parallel option,  look on youtube for an explanation,
    //**********************************************************************************************
    matrix* augmented = new matrix(equations->getx(), equations->gety() + 1);

    for (int i = 0; i < equations->getx(); i++)
    {
        for (int j = 0; j < equations->gety(); j++)
        {
            augmented->get(i,j) = equations->get(i, j);
        }
    }
    for (int i = 0; i < equations->getx(); i++)
    {
        augmented->get(i, equations->gety() + 1) = solution->arr[i];
    }
    int n = constants->getDim();
    int rows = equations->getx();
    for(int i = 0; i < n; i++)
    {
        for(int j = i + 1; j < n; j++)
        {
            if( abs(equations->arr[i * i + i]) < abs(equations->arr[j * i + i]))
            {
                for(int k = 0; k < n + 1; k++)
                {
                    /* swapping mat[i][k] and mat[j][k] */
                    equations->arr[i * rows + k] = equations->arr[i * rows + k] + equations->arr[j * rows + k];
                    equations->arr[j * rows + k] = equations->arr[i * rows + k] - equations->arr[j * rows + k];
                    equations->arr[i * rows + k] = equations->arr[i * rows + k] - equations->arr[j * rows + k];
                }
            }
        }
    }
    /* performing Gaussian elimination */
    for(int i = 0; i < n - 1; i++)
    {
        for(int j = i + 1; j < n; j++)
        {
            float f = equations->arr[j * rows + i] / equations->arr[i * rows + i];
            for(int k = 0; k < n + 1; k++)
            {
                equations->arr[j * rows + k] = equations->arr[j * rows + k] - f * equations->arr[i * rows + k];
            }
        }
    }
    /* Backward substitution for discovering values of unknowns */
    for(int i = n-1; i>=0; i--)
    {
        solution->arr[i] = equations->arr[i * rows + n];
        for(int j = i+1; j < n; j++)
        {
            if(i != j)
            {
                solution->arr[i] = solution->arr[i] - equations->arr[i * rows + j] * solution->arr[j];
            }
        }
        solution->arr[i] = solution->arr[i] / equations->arr[i * rows + i];
    }
    //**********************************************************************************************
    return 0;

}

