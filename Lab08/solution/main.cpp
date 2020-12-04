#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/gaussj.c"


// helper functions 
float f1(float x);
float f2(float x);
float getStep(float n, float min, float max);

void fillXVector(float *x, float n, float min, float max);
void fillYVector(float *result, float *x, int n, float (*f)(float));
void fillVectors(float *x, float *y, float n, float min, float max, float (*f)(float));

void fillAMartrix(float **A, int n, float min, float max);
void fillDMatrix(float **d, int n, float *y, float min, float max);

// main functions
void wyzM(float **A,float **d, float *y, int n, float min, float max);
float wyzSx(float *xw, float *yw, float **m, int n, float x, float min, float max);

int main()
{
	// set variables
	float xmin = -5.0;
	float xmax = 5.0;
	float delta = 0.01;
	float realValue, Sx;
	int N = 10;
	int n[] = {5, 8, 21};

//===========================================
// by alokowac i zwalniac pamiec jeden raz 
// tworze na poczatku tablice wymiaru dla
// najwiekszego n (= 21)
//===========================================
	float *x = vector(1, n[2]);
	float *y = vector(1, n[2]);

	float **A = matrix(1, n[2], 1, n[2]);
	float **d = matrix(1, n[2], 1, 1);

	FILE *file_pochodne, *file_f1, *file_f2;
	file_pochodne  = fopen("pochodne.dat", "w+");
	file_f1  = fopen("f1.dat", "w+");
	file_f2  = fopen("f2.dat", "w+");


	// main part 1 (pochodne)
	fillVectors(x, y, N, xmin, xmax, f1);
	wyzM(A, d, y, N, xmin, xmax);

	for (int i = 1; i <= N ; ++i)
	{
		realValue = ( f1(x[i]-delta)-2*f1(x[i])+f1(x[i]+delta) )/ (delta*delta);
		fprintf(file_pochodne , "%g  %g  %g\n", x[i], d[i][1], realValue);
    }


	// main part 2 (f1)
    for(int i = 0; i < 3; ++i)
    {
		fillVectors(x, y, n[i], xmin, xmax, f1);
		wyzM(A, d, y, n[i], xmin, xmax);
			
		for(float j = xmin; j<=xmax + 0.01; j+=0.1)
		{
			Sx = wyzSx(x, y, d, n[i], j, xmin, xmax);
			fprintf(file_f1 , "%g  %g\n", j, Sx);
		}
		fprintf(file_f1, "\n\n");
    }

	
	// main part 3 (f2)
    for(int i = 0; i < 3; ++i)
    {
		fillVectors(x, y, n[i], xmin, xmax, f2);
		wyzM(A, d, y, n[i], xmin, xmax);
			
		for(float j = xmin; j<=xmax + 0.01; j+=0.1)
		{
			Sx = wyzSx(x, y, d, n[i], j, xmin, xmax);
			fprintf(file_f2 , "%g  %g\n", j, Sx);
		}
		fprintf(file_f2, "\n\n");
    }


	// free memory
	free_vector(x, 1, n[2]);
	free_vector(y, 1, n[2]);

	free_matrix(A, 1, n[2], 1, n[2]);
	free_matrix(d, 1, n[2], 1, 1);
	
	fclose(file_pochodne );
	fclose(file_f1);
	fclose(file_f2);
	return 0;
}


//////////////////////////////////////////////////////////////////////////////////
float f1(float x)
{
	return 1/(1+(x*x));
}
//////////////////////////////////////////////////////////////////////////////////
float f2(float x)
{
	return cos(2*x);
}
//////////////////////////////////////////////////////////////////////////////////
float getStep(float n, float min, float max)
{
	return (max - min) / (n-1);
} 
//////////////////////////////////////////////////////////////////////////////////
void fillXVector(float *x, float n, float min, float max)
{
	float step = getStep(n, min, max);
	for (int i=1; i<=n; ++i)
		x[i] = min + ((i-1) * step);
}
//////////////////////////////////////////////////////////////////////////////////
void fillYVector(float *result, float *x, int n, float (*f)(float))
{
	for(int i = 1; i<=n; ++i)
		result[i] = f(x[i]);
}
//////////////////////////////////////////////////////////////////////////////////
void fillVectors(float *x, float *y, float n, float min, float max, float (*f)(float))
{
	fillXVector(x, n, min, max);
	fillYVector(y, x, n, f);
}
//////////////////////////////////////////////////////////////////////////////////
void fillAMartrix(float **A, int n, float min, float max)
{
	float milambda = 0.5;

	for(int i = 1; i<=n; ++i)
	{
		for(int j = 1; j<=n; ++j)
		{
			if(i==j)
				A[i][j] = 2;
			else if(i == j + 1 || i+1 == j)
				A[i][j] = milambda;
			else
				A[i][j] = 0;
		}
	}
	A[1][2] = A[n][n-1] = 0;
	A[1][1] = A[n][n] = 1;
}
//////////////////////////////////////////////////////////////////////////////////
void fillDMatrix(float **d, int n, float *y, float min, float max)
{
	float step = getStep(n, min, max);
	d[1][1] = d[n][1] = 0;

	for(int i = 2; i<n; ++i)
		d[i][1] = (6/(step*2))*((y[i+1] + y[i-1] - 2*y[i])/step);
}
//////////////////////////////////////////////////////////////////////////////////
void wyzM(float **A,float **d, float *y, int n, float min, float max)
{
	fillAMartrix(A, n, min, max);
	fillDMatrix(d, n, y, min, max);
	gaussj(A, n, d, 1);
}
//////////////////////////////////////////////////////////////////////////////////
float wyzSx(float *xw, float *yw, float **m, int n, float x, float min, float max)
{
	int i = 2;
	while(i <= n)
	{
		if(x <= xw[i] && x >= xw[i-1])
			break;
		++i;
	}

	float step = getStep(n, min, max);
	float Ai = (yw[i] - yw[i-1])/step - (step/6)*(m[i][1] - m[i-1][1]);
	float Bi = yw[i-1] - m[i-1][1] * ( (step*step) / 6);

	float Sx = ( m[i-1][1]*(xw[i]-x)*(xw[i]-x)*(xw[i]-x) )/ (6*step);
	Sx += ( m[i][1]*(x - xw[i-1])*(x - xw[i-1])*(x - xw[i-1]) )/ (6*step);
	Sx += Ai * (x - xw[i-1]);
	Sx += Bi;

	return Sx;
}