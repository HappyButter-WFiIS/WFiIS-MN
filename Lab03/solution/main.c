#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 2001

int main()
{
	float *b, *d_0, *d_1, *d_2, *x_s, *x_n;
	float beta = 0.4;
	float F_0 = 0.1;
	float omega = 0.8; 
	float V_0 = 0.0;
	float x_0 = 1.0;
	float w = 1.0;
	float h = 0.02;

	float a1 = 1;
	float a2 = w*w*h*h - 2.0 - beta*h;
	float a3 = 1 + beta*h;


	// uzupełniamy b
	b = malloc(N * sizeof(float));
	b[0] = 1.0;
	b[1] = 0.0;

	for (int i = 2; i < N; ++i)
	{
		b[i] = F_0 * sin(omega*i*h) * h * h;
 	}

 	// uzupełniam d_0
 	d_0 = malloc(N * sizeof(float));
 	d_0[0] = d_0[1] = 1.0;

 	for (int i = 2; i < N; ++i)
 		d_0[i] = a3;

 	// uzupełniam d_1
 	d_1 = malloc(N * sizeof(float));
 	d_1[0] = 0.0; 
 	d_1[1] = -1.0;

	for (int i = 2; i < N; ++i)
 		d_1[i] = a2;

 	// uzupełniam d_2
 	d_2 = malloc(N * sizeof(float));
 	d_2[0] = 0.0; 
 	d_2[1] = 0.0;

	for (int i = 2; i < N; ++i)
 		d_2[i] = a1;


 	// alokuje pamięć dla x_s i x_n
 	x_s = calloc(N , sizeof(float));
 	x_n = calloc(N , sizeof(float));

 	// obiczenia
 	float sum_s = 0.0;
 	float sum_n = 0.0;

 	int it = 0;
 	while (it < 100000)
 	{
 		++it;

 		x_n[0] = (1/d_0[0]) * b[0];
 		x_n[1] = (1/d_0[1]) * (b[1] - d_1[1] * x_s[0]);

 		sum_n += x_n[0] * x_n[0] + x_n[1] * x_n[1];

	 	for (int i = 2; i < N; ++i)
	 	{
	 		x_n[i] = (1/d_0[i]) * (b[i] - d_1[i] * x_s[i-1] - d_2[i] * x_s[i-2]); 
	 		sum_n += x_n[i] * x_n[i];
	 	}

	 	if (fabs(sum_n - sum_s) < 0.000001)
	 		break;

	 	for (int i = 0; i < N; ++i)
	 		x_s[i] = x_n[i];
	 	
	 	sum_s = sum_n;
	 	sum_n = 0.0;
 	}

	for(int i = 0; i < N; ++i)
		printf("%.2f %g\n", i*h, x_n[i]);

  printf("\n\n");
	


// ================================================================
// 							free
// ================================================================
 	free(b);
 	free(d_0);
 	free(d_1);
 	free(d_2);
 	free(x_s);
 	free(x_n);
	return 0;
}
