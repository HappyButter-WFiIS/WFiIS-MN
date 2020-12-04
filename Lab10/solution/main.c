#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

#define N 200
#define xmin -10.0
#define ymin -10.0
#define xmax 10.0
#define ymax 10.0


double d_rand ( const double min, const double max);
double function(double x, double y);


int main()
{
	// Use current time as  
    // seed for random generator 
    srand(time(0)); 


	// set variables
	double* x = malloc(N * sizeof(double));
	double* y = malloc(N * sizeof(double));

	for(int i = 0; i<N; ++i)
		x[i] = y[i] = 5.0;

	double T;
	double deltaX, deltaY;

	FILE * file_w0;
	file_w0 = fopen("w0.dat", "w+");

	FILE * file_T;
	file_T = fopen("T.dat", "w+");
	

	// main part
	for(int it = 0; it <= 20; ++it)
	{
		T = 10.0/pow(2, it);

		for(int k = 0; k<100; ++k)
		{
			for(int i = 0; i<N; ++i)
			{
				deltaX = d_rand(-1, 1);
				deltaY = d_rand(-1, 1);

				if(x[i] + deltaX > xmax)
					deltaX = 10.0 - x[i];
				if(x[i] + deltaX < xmin)
					deltaX = -10.0 - x[i];
				if(y[i] + deltaY > ymax)
					deltaY = 10.0 - y[i];
				if(y[i] + deltaY < ymin)
					deltaY = -10.0 - y[i];

				if(function(x[i] + deltaX, y[i] + deltaY) < function(x[i], y[i]))
				{
					x[i] = x[i] + deltaX;
					y[i] = y[i] + deltaY;	 
				}
				else if(d_rand(0, 1) < exp(-(function(x[i] + deltaX, y[i] + deltaY) - function(x[i], y[i])) / T ))
				{
					x[i] = x[i] + deltaX;
					y[i] = y[i] + deltaY;
				}
			}
			fprintf(file_w0, "%g\n", function(x[0], y[0]));
		}
		if(it == 0 || it == 7 || it == 20)
		{	
			for (int i = 0; i < N; ++i)
				fprintf(file_T, "%g  %g\n", x[i], y[i]);
			fprintf(file_T, "\n\n");
		}
	}

	int minIndex = 0;
	double minValue = function(x[0], y[0]);
	for(int i = 1; i<N; ++i)
	{
		if(function(x[i], y[i]) < minValue)
		{
			minIndex = i;
			minValue = function(x[i], y[i]);
		}
	}

	printf("\nMin value found at point: (%g, %g)\n", x[minIndex], y[minIndex]);
	printf("And the min value is... %g\n", function(x[minIndex], y[minIndex]));
	
	// free
	free(x);
	free(y);
	fclose(file_w0);
	fclose(file_T);

	return 0;
}


/////////////////////////////////////////////////////////////////////////////////
double d_rand ( const double min, const double max)
{
	double r = (double)rand() / RAND_MAX; 
	r = r * (max - min) + min;
	return r;
}
/////////////////////////////////////////////////////////////////////////////////
double function(double x, double y)
{
	return (double) sin(x)*sin(y)\
					- exp( -(x+(M_PI/2))*(x+(M_PI/2))\
			      		   -(y-(M_PI/2))*(y-(M_PI/2)));
}
/////////////////////////////////////////////////////////////////////////////////


