#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double* fillYArray(double*, int);
double* fillXArray(int size, double min, double max);
double* fillXArrayWithCzybyszew(int size, double min, double max);
void fillFArray(double**, double*, int);
double getNewtonInterpolation(double**F, double* xm, double x, int size);


int main()
{
	// set variables for Newton interpolation 
	int n[] = {5,10,15,20};
	double xmin = -5.0;
	double xmax = 5.0;
	double W = 0;

	FILE *file_1;
	file_1 = fopen("zad_1.dat", "w+");

	double** x = malloc( 4 * sizeof(double*) );
	double** y = malloc( 4 * sizeof(double*) );
	double***f = malloc( 4 * sizeof(double**) );

	for(int i = 0; i < 4; ++i)
	{
		x[i] = fillXArray(n[i]+1, xmin, xmax);
		y[i] = fillYArray(x[i], n[i]+1);

		f[i] = malloc((n[i]+1) * sizeof(double*));
		for(int j = 0; j < n[i]+1; ++j)
		{
			f[i][j] = malloc((n[i]+1) * sizeof(double));
			f[i][j][0] = y[i][j];
		}
	}

	for(int i = 0; i < 4; ++i)
		fillFArray(f[i], x[i], n[i] + 1);


	// main part
	for(int i = 0; i<4; ++i)
	{
		for(float j = xmin; j<=xmax + 0.001; j+=0.01)
		{
			W = getNewtonInterpolation(f[i], x[i], j, n[i]+1);
			// printf("%f  %g\n", j, W); 						// checking if filled right 
			fprintf(file_1, "%f  %g\n", j, W);
		}
		fprintf(file_1,"\n\n");
	}


	// free
	for(int i = 0 ; i<4; ++i)
	{
		free(x[i]);
		free(y[i]);
	
		for(int j = 0; j < n[i]+1; ++j)
			free(f[i][j]);
		free(f[i]);
	}

	free(x);
	free(y);
	free(f);
	fclose(file_1);

//==========================================================================================================
	// set variables for Newton interpolation with Czybyszew zeros
	FILE *file_2;
	file_2 = fopen("zad_2.dat", "w+");

	x = malloc( 4 * sizeof(double*) );
	y = malloc( 4 * sizeof(double*) );
	f = malloc( 4 * sizeof(double**) );

	for(int i = 0; i < 4; ++i)
	{
		x[i] = fillXArrayWithCzybyszew(n[i]+1, xmin, xmax);
		y[i] = fillYArray(x[i], n[i]+1);

		f[i] = malloc((n[i]+1) * sizeof(double*));
		for(int j = 0; j < n[i]+1; ++j)
		{
			f[i][j] = malloc((n[i]+1) * sizeof(double));
			f[i][j][0] = y[i][j];
		}
	}

	for(int i = 0; i < 4; ++i)
		fillFArray(f[i], x[i], n[i] + 1);

	// main part
	for(int i = 0; i<4; ++i)
	{
		for(float j = xmin; j<=xmax + 0.001; j+=0.01)
		{
			W = getNewtonInterpolation(f[i], x[i], j, n[i]+1);
			// printf("%f  %g\n", j, W); 						// checking if filled right 
			fprintf(file_2, "%f  %g\n", j, W);
		}
		fprintf(file_2,"\n\n");
	}


	// free
	for(int i = 0 ; i<4; ++i)
	{
		free(x[i]);
		free(y[i]);
	
		for(int j = 0; j < n[i]+1; ++j)
			free(f[i][j]);
		free(f[i]);
	}

	free(x);
	free(y);
	free(f);
	fclose(file_2);

	//-----------------------------------------------------------
	// nie do konca przemyslalem zarzadzanie pamiecia 
	// na poczatku zadania. Mimo to program dziala zgodnie
	// z oczekiwaniami :) 
	//-----------------------------------------------------------

	return 0;	
}

////////////////////////////////////////////////////////////////////////////////////////////
double* fillYArray(double* points, int size)
{
	double* results = malloc(size*sizeof(double));

	for(int i = 0; i < size; ++i)
	{
		results[i] = 1/(points[i]*points[i] + 1);
		// printf("Y[%d] = %g\n", i, results[i]);  // checking if filled right  
	}
	return results;
}

////////////////////////////////////////////////////////////////////////////////////////////
double* fillXArray(int size, double min, double max)
{
	double* result = malloc(size*sizeof(double));
	double 	step = (max - min) / (size - 1);
	
	for(int i = 0; i < size; ++i)
	{
		result[i] = min + (step * i);
	  // printf("x[%d] = %g\n", i, result[i]);    // checking if filled right
	}
	return result;
}
////////////////////////////////////////////////////////////////////////////////////////////
double* fillXArrayWithCzybyszew(int size, double min, double max)
{
	double* result = malloc(size*sizeof(double));

	for(int i = 0; i < size; ++i)
	{
		result[i] = 0.5*((min-max)*cos(M_PI*(2*i + 1)/(2*(size-1)+2)));
	  	// printf("x[%d] = %g\n", i, result[i]);    // checking if filled right
	}
	return result;
}
////////////////////////////////////////////////////////////////////////////////////////////
void fillFArray(double**F, double*X, int size)
{
	for(int j = 1; j < size; ++j)
	{
		for(int i = j; i < size; ++i)
		{
			F[i][j] = (F[i][j-1] - F[i-1][j-1]) / (X[i] - X[i-j]);
			// if(i==j)								
			// 	printf("%g\t", F[i][j]);			// checking if filled right
		}
	}
}
////////////////////////////////////////////////////////////////////////////////////////////
double getNewtonInterpolation(double**F, double* xm, double x, int size)
{
	double multi;
	double result = 0.0;

	for(int j = 0; j<size; ++j)
	{
		multi = 1.0;
		for(int i = 0; i<=j-1; ++i)
			multi *= (x-xm[i]);
		result += F[j][j] * multi; 
	}
	return result;
}