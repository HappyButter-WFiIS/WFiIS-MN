#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/tred2.c"
#include "/opt/NR/numerical_recipes.c/tqli.c"
#include "/opt/NR/numerical_recipes.c/pythag.c"

#define n 7

/* My functions */

// fill functions
void fillAMatrix(float**);
void fillIdentityMatrix(float**);
void fillVector(float*);
// operation functions
void multiplyMatrixByVector(float**, float*, float*);
float dotProduct(float*, float*);
void HotellingReduction(float**, float, float*);
void normalize(float*);
void copyVector(float*, float*);
// printing functions
void printMatrix(float**);
void printVector(float *);

int main()
{
	// set variables
	float **A, **W, **Z;
	float *d, *e;
	
	float *x_previous;
	float *x_current;

	float lambda = 0;

	A = matrix(1, n, 1, n);
	W = matrix(1, n, 1, n);
	Z = matrix(1, n, 1, n);
	
	d = vector(1,n);
	e = vector(1,n);
	x_previous = vector(1,n);
	x_current = vector(1,n);

	
	// fill A, W & Z matrixes
	fillAMatrix(A);
	fillAMatrix(W);
	fillIdentityMatrix(Z);

	// printMatrix(A);			// uncomment to print A matrix

	// A -> T   	| Numerical Recipes libary usage
	tred2(A, n, d, e);
	tqli(d, e, n, Z);

  	printf("==========================================================");
	printf("\n\t* Results after numerical_recipes libary useage:\n\n");
	printVector(d);

  // iterative method
	printf("\n==========================================================");
	printf("\n\t* Results after iterative method useage:\n\n");
	
  
  for(int k = 1; k<=n; ++k)
	{
		fillVector(x_previous);
		for(int i = 1; i<=8; ++i)
		{
			multiplyMatrixByVector(W, x_previous, x_current);
			lambda = (dotProduct(x_current, x_previous) / dotProduct(x_previous, x_previous));
			normalize(x_current);
			copyVector(x_previous, x_current);
		}
		HotellingReduction(W, lambda, x_previous);
    	printf("%d iteration last lambda value = %g\n", k, lambda);
	}


	// free the memory
	free_matrix(A, 1, n, 1, n);
	free_matrix(W, 1, n, 1, n);
	free_matrix(Z, 1, n, 1, n);

	free_vector(d, 1, n);    
	free_vector(e, 1, n);
	free_vector(x_previous, 1, n);
	free_vector(x_current, 1, n);

	return 0;
}


///////////////////////////////////////////////////////
void fillAMatrix(float** matrix)
{
 	for(int i = 1; i <=n; ++i)
 		for(int j = 1; j<=n; ++j)
    		matrix[i][j] = sqrt(i+j);
}
///////////////////////////////////////////////////////
void fillIdentityMatrix(float ** matrix)
{
 	for(int i = 1; i <=n; ++i)
 		for(int j = 1; j<=n; ++j){
 			if(i == j)
    			matrix[i][j] = 1;
    		else
    			matrix[i][j] = 0;
 		}
}
///////////////////////////////////////////////////////
void fillVector(float* vec)
{
	for(int a = 1; a<=n; ++a)
		vec[a] = 1;
}
///////////////////////////////////////////////////////
void multiplyMatrixByVector(float** matrix, float* vec, float* result)
{
  float sum = 0.0;
	for(int a = 1; a<=n; ++a)
	{
     sum = 0.0;
		for(int b = 1; b<=n; ++b)
		{	sum += matrix[a][b] * vec[b];}
    result[a] = sum;
  }
}
///////////////////////////////////////////////////////
float dotProduct(float* vec_1, float* vec_2)
{
	float sum = 0.0;
	for(int a = 1; a<=n; ++a)
		sum += (vec_1[a] * vec_2[a]);
	return sum; 
}
///////////////////////////////////////////////////////
void HotellingReduction(float** WResult, float lam, float* x)
{
	for(int a = 1; a <=n; ++a)
 		for(int b = 1; b<=n; ++b)
    		WResult[a][b] -= (lam * x[a] * x[b]); 
}

///////////////////////////////////////////////////////
void normalize(float* x)
{
	float normValue = sqrt(dotProduct(x,x));
	for(int a = 1; a<=n; ++a)
		x[a] = (x[a] / normValue);
}
///////////////////////////////////////////////////////
void copyVector(float* dest, float* source)
{
	for(int a = 1; a<=n; ++a)
		dest[a] = source[a];
}
///////////////////////////////////////////////////////
void printMatrix(float** toPrintMatrix)
{
 	for(int i = 1; i <=n; ++i){
 		for(int j = 1; j<=n; ++j)
    		printf("%g\t", toPrintMatrix[i][j]);
    	printf("\n");
    }
}	
///////////////////////////////////////////////////////
void printVector(float * vec)
{
	for(int i = 1; i<=n; ++i)
		printf("%g\t", vec[i]);
}