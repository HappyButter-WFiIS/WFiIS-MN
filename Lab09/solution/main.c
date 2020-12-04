#include <stdio.h>
#include <stdlib.h> 
#include <math.h>
#include <time.h>

#define xmax 4.0
#define xmin -4.0
#define x0 2.0
#define N 201
#define frand() ((double)rand())/(RAND_MAX+1.0)

// functions to set points to approximation
double fun(double x);
double Crand();
double fNoise(double x);
void fillX(double* x);
void fillY(double* y, double *x);

// set Gram matrix
double getAlpha(double *phi, double *x);
double getBeta(double *phi_prev, double *phi, double *x);
void fillPhi(double** phi, double* x, int m);

// calculate aproximation
double getC(double* y, double* phi);
double getS(double* phi);
void getApproximation(double* y, double** phi, int m, double** F);

// save to file
void saveToFileGram(double** phi, double* x);
void saveToFilePoints(double* x, double* y);
void saveToFileApprox(double* x, double** F);


int main()
{
	// Use current time as  
    // seed for random generator 
    srand(time(0)); 

	// set variables
	int m = 51;

	double* x = malloc(sizeof(double) * N);
	double* y = malloc(sizeof(double) * N);
	double** F = malloc(sizeof(double*) * 3);
	double** phi = malloc(sizeof(double*) * m); 
	
	for (int i = 0; i<m; ++i)
		phi[i] = malloc(sizeof(double) * N); 

	for (int i = 0; i<3; ++i)
		F[i] = malloc(sizeof(double) * N); 


	clock_t t; 
    t = clock(); 

	// main part 
	fillX(x);
	fillY(y, x);
/*
	// uncomment to find out what  
	// happens if there's no noise
	
	for(int i = 0; i<N; ++i)		
		y[i] = fun(x[i]);
*/	
	fillPhi(phi, x, m);
	getApproximation(y, phi, m, F);

    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
  
    printf("fun() took %f seconds to execute \n", time_taken); 

	saveToFileGram(phi, x);
	saveToFilePoints(x, y);
	saveToFileApprox(x, F);


	// free memory
	for (int i = 0; i<m; ++i)
		free(phi[i]);
	free(phi);

	for (int i = 0; i<3; ++i)
		free(F[i]);
	free(F);

	free(x);
	free(y);
	return 0;
}


/////////////////////////////////////////////////////////////////
double fun(double x)
{
	double sigma = (xmax-xmin) / 16;	
	double result;
	result = sin( 14*M_PI*x / (xmax - xmin) );
	result *= (exp(-(x-x0)*(x-x0) / (2*sigma*sigma)) + exp(-(x+x0)*(x+x0) / (2*sigma*sigma)));
	return result;
}
/////////////////////////////////////////////////////////////////
double Crand()
{
	return (frand() - 0.5) / 5;
}
/////////////////////////////////////////////////////////////////
double fNoise(double x)
{
	return fun(x) + Crand();
}
/////////////////////////////////////////////////////////////////
void fillX(double* x)
{
	double step = (xmax-xmin) / (N-1);
	for(int i = 0; i<N; ++i)
		x[i] = xmin + step * i;
}
/////////////////////////////////////////////////////////////////
void fillY(double* y, double *x)
{
	for(int i = 0; i<N; ++i)
		y[i] = fNoise(x[i]);	
}
/////////////////////////////////////////////////////////////////
double getAlpha(double *phi, double *x)
{
	double counter = 0.0, denominator = 0.0;

	for(int i = 0; i < N; ++i)
	{
		counter += x[i] * phi[i] * phi[i];
		denominator += phi[i] * phi[i];
	}
	
	return (double) counter/denominator;
}
/////////////////////////////////////////////////////////////////
double getBeta(double *phi_prev, double *phi, double *x)
{
	double counter = 0.0, denominator = 0.0;
	for (int i = 0; i < N; ++i)
	{
		counter += x[i] * phi_prev[i] * phi[i];
		denominator += phi_prev[i] * phi_prev[i];
	}

	return (double) counter/denominator;
}
/////////////////////////////////////////////////////////////////
void fillPhi(double** phi, double* x, int m)
{
	double alpha, beta;

	for (int i = 0; i < N; ++i)
		phi[0][i] = 1; 
	
	alpha = getAlpha(phi[0], x);
	for (int i = 0; i < N; ++i)
		phi[1][i] = (x[i] - alpha) * phi[0][i]; 

	for (int i = 2; i<m; ++i)
	{
		alpha = getAlpha(phi[i-1], x);
		beta = getBeta(phi[i-2], phi[i-1], x);

		for(int j = 0; j<N; ++j)
		{
			phi[i][j] = (x[j] - alpha) * phi[i-1][j];
			phi[i][j] -= beta * phi[i-2][j];
		}
	}
}
/////////////////////////////////////////////////////////////////
double getC(double* y, double* phi)
{
	double sum = 0.0;
	for(int i = 0; i < N; ++i)
		sum += y[i] * phi[i];
	return sum;
}
/////////////////////////////////////////////////////////////////
double getS(double* phi)
{
	double sum = 0.0;
	for(int i = 0; i < N; ++i)
		sum += phi[i] * phi[i];	
	return sum;
}
/////////////////////////////////////////////////////////////////
void getApproximation(double* y, double** phi, int m, double** F)
{
	double sum;
	for(int i = 0; i<N; ++i)
	{
		sum = 0.0;
		for(int j = 0; j < m; ++j)
		{
			sum += (getC(y, phi[j]) / getS(phi[j])) * phi[j][i];
			if(j == 10)
				F[0][i] = sum;
			if(j == 30)
				F[1][i] = sum;
			if(j == 50)
				F[2][i] = sum;
		}
	}	 
}
/////////////////////////////////////////////////////////////////
void saveToFileGram(double** phi, double* x)
{
	FILE * file_Gram;
	file_Gram = fopen("Gram.dat", "w+");

	for(int i = 0; i<N; ++i)
	{
		fprintf(file_Gram, "%g  ", x[i]);
		for(int j = 0; j<7; ++j)
			fprintf(file_Gram, "%g  ", phi[j][i] / phi[j][0]);
		fprintf(file_Gram, "\n");
	}

	fclose(file_Gram);	
}
/////////////////////////////////////////////////////////////////
void saveToFilePoints(double* x, double* y)
{
	FILE * file_points;
	file_points = fopen("pkt.dat", "w+");

	for(int i = 0; i < N; ++i)
		fprintf(file_points, "%g  %g\n", x[i], y[i]);

	fclose(file_points);
}
/////////////////////////////////////////////////////////////////
void saveToFileApprox(double* x, double** F)
{
	FILE * file_approx;
	file_approx = fopen("approx.dat", "w+");

	for (int i = 0; i<3; ++i)
	{
		for(int j = 0; j < N; ++j)
			fprintf(file_approx, "%g  %g\n", x[j], F[i][j]);
		fprintf(file_approx, "\n\n");
	}
	
	fclose(file_approx);
}