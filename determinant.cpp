#include <stio.h>
#include <stdlib.h>
#inlucde <math.h>

float determinant(float *A,int dim){
    if (dim == 1){
        return A[0];
    }
    // base case
    if (dim==2){
        float ret{A[0]*A[3]-A[1]*A[2]};
        return ret;
    }
    int sign = 0; //0 to add 1 to subtract
    int det = 0;
    // for each value in the top row
    for (int x=0;x<dim;x++){
        // make new matrix
        float subA[(dim-1)*(dim-1)];
        int idx = 0;
        for (int i=1;i<dim;i++){
            for (int j=0;j<dim;j++){
                if (j!=x){
                    subA[idx] = A[i+dim*j];
                    idx++;
                }
            }
        }
        if (sign==0){
            det += A[x]*determinant(subA,dim-1);
            sign = 1;
        }
        else{
            det -= A[x]*determinant(subA,dim-1);
            sign = 0;
        }
    }
    return det;
}


int main(){
    float A[] = {1,0,0,0,0,4,0,0,0,0,9,0,0,0,0,1};
    float det = determinant(A,4);
    // 36
    printf("deta: %.0f\n",det);
    float B[] = {1,0,0,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,1,0,0,0,0,0,1};
    float detb = determinant(B,5);
    // 2
    printf("detb: %.0f\n",detb);
   return 0;
}
