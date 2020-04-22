#include <iostream>


// GLM library for math (3x3 and 4x4)
// lapack for linear algebra package
// return an int instead of the matrix -- put result matrix as an argument

// C++ program to find adjoint and inverse of a matrix 
//#include<bits/stdc++.h> 
using namespace std; 
#define N 4

void getCofactor(float *A, float* temp, int p, int q, int n, int order);
int determinant(float* A, int n, int order);
void adjoint(float* A,int* adj, int order);
bool inverse(float* A, float* inverse, int order); 

  
// Function to get cofactor of A[p][q] in temp[][]. n is current 
// dimension of A[][]
void getCofactor(float* A, float *temp, int p, int q, int n, int order) 
{ 
    int i = 0, j = 0; 
  
    // Looping for each element of the matrix 
    for (int row = 0; row < n; row++) 
    { 
        for (int col = 0; col < n; col++) 
        { 
            //  Copying into temporary matrix only those element 
            //  which are not in given row and column 
            if (row != p && col != q) 
            { 
                temp[(i) + (j++)*order] = A[row + col*order]; 
  
                // Row is filled, so increase row index and 
                // reset col index 
                if (j == n - 1) 
                { 
                    j = 0; 
                    i++; 
                } 
            } 
        } 
    } 
} 
  
/* Recursive function for finding determinant of matrix. 
   n is current dimension of A[][]. */
int determinant(float* A, int n, int order) 
{ 
    int D = 0; // Initialize result 
  
    //  Base case : if matrix contains single element 
    if (n == 1) 
        return A[0]; 
  
    float temp[order*order]; // To store cofactors 
  
    int sign = 1;  // To store sign multiplier 
  
     // Iterate for each element of first row 
    for (int f = 0; f < n; f++) 
    { 
        // Getting Cofactor of A[0][f] 
        getCofactor(A, temp, 0, f, n, order); 
        D += sign * A[0+f*order] * determinant(temp, n - 1, order); 
  
        // terms are to be added with alternate sign 
        sign = -sign; 
    } 
  
    return D; 
} 
  
// Function to get adjoint of A[N][N] in adj[N][N]. 
void adjoint(float* A, int* adj, int order) 
{ 
    if (order == 1) 
    { 
        adj[0] = 1; 
        return; 
    } 
  
    // temp is used to store cofactors of A[][] 
    int sign = 1;
    float temp[order * order]; 
  
    for (int i=0; i<order; i++) 
    { 
        for (int j=0; j<order; j++) 
        { 
            // Get cofactor of A[i][j] 
            getCofactor(A, temp, i, j, order, order); 
  
            // sign of adj[j][i] positive if sum of row 
            // and column indexes is even. 
            sign = ((i+j)%2==0)? 1: -1; 
  
            // Interchanging rows and columns to get the 
            // transpose of the cofactor matrix 
            adj[j+i*order] = (sign)*(determinant(temp, order-1, order)); //THIS CHANGE MIGHT BE WRONG
        } 
    } 
} 
  
// Function to calculate and store inverse, returns false if 
// matrix is singular 
bool inverse(float* A, float* inverse, int order) 
{ 
    // Find determinant of A[][] 
    int det = determinant(A, order, order); 
    if (det == 0) 
    { 
        cout << "Singular matrix, can't find its inverse"; 
        return false; 
    } 
  
    // Find adjoint 
    int adj[order*order]; 
    adjoint(A, adj, order); 
  
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)" 
    for (int i=0; i<order; i++) 
        for (int j=0; j<order; j++) 
            inverse[i + j*order] = adj[i + j*order]/float(det); 
  
    return true; 
} 
  
// Generic function to display the matrix.  We use it to display 
// both adjoin and inverse. adjoin is integer matrix and inverse 
// is a float. 
template<class T> 
void display(T A, int order) 
{ 
    for (int i=0; i<order; i++) 
    { 
        for (int j=0; j<order; j++) 
            cout << A[i + j*order] << " "; 
        cout << endl; 
    } 
} 
  
// Driver program 
int main() 
{ 
    // int A[N][N] = { {5, -2, 2, 7}, 
    //                 {1, 0, 0, 3}, 
    //                 {-3, 1, 5, 0}, 
    //                 {3, -1, -9, 4}};
    //col + row*order
    float A[16] = {5, 1, -3, 3, -2, 0, 1, -1, 2, 0, 5, -9, 7, 3, 0 , 4};
    int order = 4;
  
    //int adj[N][N];  // To store adjoint of A[][] 
    // float inv[N][N]; // To store inverse of A[][] 
    
    int adj[order*order];
    float inv[order*order]; 
    cout << "Input matrix is :\n"; 
    display(A, 4); 
  
    cout << "\nThe Adjoint is :\n"; 
    adjoint(A, adj, 4); 
    display(adj, 4); 
  
    cout << "\nThe Inverse is :\n"; 
    if (inverse(A, inv, 4)) 
        display(inv, 4); 
  
    return 0; 
} 