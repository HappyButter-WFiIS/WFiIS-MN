#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>

#define N 4

void print(const gsl_matrix * mat )
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%g\t", gsl_matrix_get(mat, i, j));
		}
		printf("\n");
	}
}


int main()
{
// dane wejsciowe
	gsl_matrix *A;
	A = gsl_matrix_calloc(N,N);

	for(int i = 0; i<N; ++i)
	{
		for (int j = 0; j<N; ++j)
		{
			double value = 1.0/(i+j+2.0);
			gsl_matrix_set(A, i, j, value);
		}
	}

	gsl_matrix *AAA;
	AAA = gsl_matrix_calloc(N,N);
	gsl_matrix_memcpy(AAA,A);
	
	int signum = 1;
	gsl_permutation * p;
	p = gsl_permutation_alloc(N);

// Ad 1
	gsl_linalg_LU_decomp(A, p, &signum);
// Ad 2 
	double detA = signum;
	for (int i = 0; i<N; i++)
	{
		for(int j = 0; j<N; j++)
		{		
			if (i==j)
			{
				detA = detA * gsl_matrix_get(A, i, j);	
				printf("%g\t", gsl_matrix_get(A, i, j));			
			}
		}
	}
	printf("\n\n\ndet(A')= %g\n", detA );

// Ad 3
	gsl_matrix *A_prim;
	A_prim = gsl_matrix_calloc(N,N);


	for(int i = 0; i < N; i++)
	{
		gsl_vector *b;
		b = gsl_vector_calloc(N);
		gsl_vector_set(b,i,1.0);
		gsl_vector *x;
		x = gsl_vector_calloc(N);

		gsl_linalg_LU_solve(A, p, b, x); 

		for (int j = 0; j<N; j++)
		{
			double value = gsl_vector_get(x,j);
			gsl_matrix_set(A_prim, j, i, value);
		}
		// free x, b 
		gsl_vector_free(b);
		gsl_vector_free(x);
	}
	printf("==================A==================\n");
	print(A);
	printf("==============A_odwrotne=============\n");
	print(A_prim);
//	Ad 4
	gsl_matrix *C;
	C = gsl_matrix_calloc(N,N);

	for (int i = 0; i<N; i++)
	{
		for (int j = 0; j<N; j++)
		{
			double result = 0.0;
			for (int k = 0; k<N; k++)
			{
				double a_1 = gsl_matrix_get(AAA, i, k);
				double a_2 = gsl_matrix_get(A_prim, k, j);
				result += a_1 * a_2;
			}
			gsl_matrix_set(C,i,j,result);
		}
	}
	printf("=============A * A_odwrotne==============\n");
	print(C);

// Ad 5
	double max_AAA = gsl_matrix_max(AAA);
	double max_A_prim = gsl_matrix_max(A_prim);
	printf("=================max====================\n");
	printf("max = %g\n", max_AAA * max_A_prim);

// zwalnianie pamieci;
	gsl_matrix_free(A);
	gsl_matrix_free(AAA);
	gsl_matrix_free(A_prim);
	gsl_matrix_free(C);
	gsl_permutation_free(p);


return 0;
}
