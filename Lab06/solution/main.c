#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define N 5
#define I_MAX 30


double licz_r(double*, double* , int , double);


int main()
{
	// set variables
	double* a = malloc((N+1) * sizeof(double));
	a[5] = 1.0;
	a[4] = 14.0;
	a[3] = 33.0;
	a[2] = -92.0;
	a[1] = -196.0;
	a[0] = 240.0;
	double* b = malloc((N+1) * sizeof(double));
	double* c = malloc(N * sizeof(double));
	double Rj, RjPrim;
	double x0, x1;
	int n;

	// calculate
	for(int L = 1; L<=N; ++L) // L - kolejne zero wielominu
	{
		n=N-L+1;
		x0 = 0.0;

		for(int i = 1; i<I_MAX; ++i) // przybliżanie zera
		{
			Rj = licz_r(a, b, n, x0);
			RjPrim = licz_r(b, c, n-1, x0);
			
			x1 = x0 - Rj/RjPrim;
			printf("zero[%d] w iteracji[%d] = %g\n", L, i, x1);

			if(fabs(x1 - x0) < 1.0E-7)
				break;

			x0 = x1; 
		}
		
		for(int i = 0; i<=n-1; ++i)
			a[i] = b[i];
	}


	// free memory
	free(a);
	free(b);
  	free(c);
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
double licz_r(double* a, double* b, int n, double x0)
{
	b[n] = 0.0;
	for(int k = n - 1; k>=0; --k)
		b[k] = a[k+1] + x0*b[k+1];

	return (double)(a[0] + x0*b[0]);  
}
///////////////////////////////////////////////////////////////////////////////