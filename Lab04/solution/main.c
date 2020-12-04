#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "/opt/NR/numerical_recipes.c/nrutil.h"
#include "/opt/NR/numerical_recipes.c/nrutil.c"
#include "/opt/NR/numerical_recipes.c/tred2.c"
#include "/opt/NR/numerical_recipes.c/tqli.c"
#include "/opt/NR/numerical_recipes.c/pythag.c"

#define nx 20
#define ny 20
#define n nx*ny
#define m 10
#define t -0.021


void print(float **);
void multiply(float **, float **, float **); // a x b = result

int main()
{
  // set variables
  float **H, **Y, **X;
  H = matrix(1, n, 1, n);
  Y = matrix(1, n, 1, n);
  X = matrix(1, n, 1, n);
  
  float *d, *e;
  d = vector(1,n);
  e = vector(1,n);
  int *indx = ivector(1, n);
  
  
  // fill H matrix
  int p;
  
  for(int i = 1; i <= nx; ++i)
  {
    for(int j = 1; j <= ny; ++j)
    {
      p = j+(i-1)*ny;
      
      for(int k=1; k <= n; ++k)
        H[p][k] = 0.0;
      if(i>1) H[p][p-ny] = t;
      if(i<nx) H[p][p+ny] = t;
      H[p][p] = -4*t;
      if(j>1) H[p][p-1] = t;
      if(j<ny) H[p][p+1] = t;
    }
  }
      
  //print(H); sprawdza³em wynik poœredni dla nx=3 i ny=3 - zgadza sie
  // fill Y matrix
  for (int i = 1; i <= n; ++i) 
  {
    for (int j = 1; j <= n; ++j) 
    {
      if (i == j) 
        Y[i][j] = 1;
      else
         Y[i][j] = 0;
    }
  }

  tred2(H, n, d, e);
  tqli(d, e, n, Y);
  
  multiply(H,Y,X);   // H*Y=X

  // sort energies
  for(int i = 1; i<=n; ++i) 
    indx[i] = i; // inicjalizacja
  
  for(int i = 1; i<=n-1; ++i)
  {
    for(int k = n; k>=i+1; --k)
    {
      float e1 = d[k-1];
      float e2 = d[k];
      float f1 = indx[k-1];
      float f2 = indx[k];
      if(e2<e1)   //wymieniamy energie i indeksy wektorów miejscami
      { 
        d[k] = e1;
        d[k-1] = e2;
        indx[k] = f1;
        indx[k-1] = f2;
      }
    }
  }
  
  // save to file 
  FILE *fp;
  fp=fopen("dane.dat","w");
  for(int i = 1; i<=nx; ++i)
  {
    for(int j = 1; j<=ny; ++j)
    {
      int z = j+(i-1)*ny;
      fprintf(fp,"%6d %6d ",i,j);
      for(int k = 1; k<=m; ++k) 
        fprintf(fp," %12.6f ",X[z][indx[k]]);
      fprintf(fp,"\n");
    }
    fprintf(fp,"\n");
  }
  fclose(fp);

  // free memory 
  free_matrix(H, 1, n, 1, n);
  free_matrix(Y, 1, n, 1, n);
  free_matrix(X, 1, n, 1, n);
  
  free_vector(d,1,n);    
  free_vector(e,1,n);  
  
  free_ivector(indx, 1, n);
  
  // end
  return 0;
}

//////////////////////////////////////////////////////////////////////////////////////////
void print(float **toPrint)
{
  for(int i = 1; i<=n; ++i)
  {
    for(int j = 1; j<=n; ++j)
      printf("%g  ", toPrint[i][j]);
    printf("\n");
  }
}
//////////////////////////////////////////////////////////////////////////////////////////
void multiply(float **a, float ** b, float **result)
{
  float rowSum;
  
  for(int i = 1; i<=n; ++i)
  {
    for(int j = 1; j<=n; ++j)
    {
      rowSum = 0.0;
      for(int k = 1; k<=n; ++k)
        rowSum += a[i][k] * b[k][j];
      result[i][j] = rowSum;
    }
  }  
}

//////////////////////////////////////////////////////////////////////////////////////////