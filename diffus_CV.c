//
// Copyright (c) 2018, Vishal_S
// All rights reserved. Please read the "license.txt" for license terms.
//
// Project Title: CFD
//
// Developer: Vishal S
//
// Contact Info: vishalsubbu97@gmail.com
//
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>


double domain[2];
int no_of_CV;
double factor;
double temp_domain0;
double temp_domain1;
double conductivity;
double area;
double source_p,source_u;
double **coeff_matrix;
double *temp;
double *X;
double *xcv;
double delx;


void write_to_file(double *func,double *X,int no_of_CV);
void  initialise(double *X,double *xcv,int no_of_location,double factor,double *domain);
void TDMA(double *temp,double **coeff_matrix,int no_of_CV);
void matrix(double *temp,double **coeff_matrix,int no_of_CV);

//////////////////////////////////////////////////////////////////////
void main()
{
  ///Input parameters ////////
  domain[0] = 0.0;
  domain[1] = 1.0;
  factor = 1.0;
  temp_domain0 = 100.0;
  temp_domain1 = 400.0;
  area = 1.0;
  conductivity = 10.0;
  source_p = 0.0;
  source_u = 0.0;
  /////////////////////////


  int option =0;
  int i=0;
  printf("\n........ Welcome........\n");
  // the domain boundaries are initialised to 0 and 1
  char choice = 'y';
  while ( (choice =='y')||(choice == 'Y'))
    {
      printf(" Enter no of CVs :");
      scanf (" %d",&no_of_CV);
      /*printf(" Enter multiplcation factor :");
	scanf (" %lf",&factor);*/

      //Allocating memroy to variables
      printf("Allocating memory to varibles.....");
      X  = (double *) malloc((no_of_CV+2) * sizeof(double));
      temp =(double *) malloc((no_of_CV+2) * sizeof(double));
      xcv = (double *) malloc((no_of_CV+1) * sizeof(double));
      coeff_matrix = (double **) malloc ((no_of_CV)*sizeof(double *));
      for (i=0;i<no_of_CV;i++)
	{
	  coeff_matrix[i] = (double *) malloc (4 * sizeof(double));
	}
      printf("done\n");

      //Initialising the grid points
      printf("Initialising the grid points....");
      initialise(X,xcv,no_of_CV,factor,domain);
      printf("done\n");
      //Computing the coefficient matrix
      printf("Computing the matrix  .....");
      matrix(temp,coeff_matrix,no_of_CV);
      printf("done\n");

      //Inverting the matrix
      printf("Inverting the matrix  .....");
      TDMA(temp,coeff_matrix,no_of_CV);
      printf("done\n");

      //Writing into file
      printf("Writing into file....");
      write_to_file(temp,X,no_of_CV);
      printf("done\n");
      
      ///Freeing Memory
      printf("Freeing Memory.... ");
      free(X);
      free(temp);
      free(xcv);
      free(coeff_matrix);
      printf("done\n");
      printf("Do you want to continue (y/n) ");
      scanf(" %c",&choice);
    }
  printf(" Program terminated.\n ");
}
/////////////////////////////////////////////////////////////////////////////////////////////


void  initialise(double *X,double *xcv,int no_of_location,double factor,double *domain)
{
  if(factor == 1.0)
    {
      int i =0;
      delx = ( ( domain[1]  - domain[0]) / no_of_location );
      X[0] = domain[0];
      X[1] = domain [0] + (delx/2);
      xcv[0] =delx/2;
      for (i=2;i<=no_of_location;i++)
	{
	  X[i] = X[i-1] + delx ;
	  xcv[i-1] = delx;
	}
      X[no_of_location+1] = X[no_of_location] + (delx/2);
      xcv[no_of_location] = delx/2;
    }
}
	  
void matrix(double *temp, double **coeff_matrix, int no_of_CV)
{
  temp[0] =temp_domain0;
  temp[no_of_CV +1 ] =temp_domain1;
  int i=0,j=0;
  double coeff;
  coeff_matrix[0][0] = 0.0;
  coeff_matrix[0][1] = 3*coeff - source_p;
  coeff_matrix[0][2] = -1*coeff;
  coeff_matrix[0][3] = 2*coeff*temp_domain0 + source_u;
  
  for (i=1;i<no_of_CV-1;i++)
    {
        coeff_matrix[i][0] = -1*coeff;
	coeff_matrix[i][1] = 2*coeff - source_p;
	coeff_matrix[i][2] = -1*coeff;
	coeff_matrix[i][3] = source_u;
    }
  coeff_matrix[no_of_CV-1][0] = -1* coeff;
  coeff_matrix[no_of_CV-1][1] = 3*coeff - source_p;
  coeff_matrix[no_of_CV-1][2] = 0.0;
  coeff_matrix[no_of_CV-1][3] = 2*coeff*temp_domain1 + source_u;
}

void TDMA(double *temp,double **coeff_matrix,int no_of_CV)
{
  int i=0;
  int j=0;
  double denom;
  for (j=3;j>=0;j--)
    {
      coeff_matrix[0][j]/=coeff_matrix[0][1];
    }
  for( i=1;i<no_of_CV-1;i++)
    {
      denom =(coeff_matrix[i][1] -(coeff_matrix[i][0]*coeff_matrix[i-1][2]));
      coeff_matrix[i][3] = (coeff_matrix[i][3] - (coeff_matrix[i-1][3]*coeff_matrix[i][0]))/ denom;
      coeff_matrix[i][2] = coeff_matrix[i][2]/denom;
      coeff_matrix[i][1] = 1.0;
      coeff_matrix[i][0] = 0;;
    }
  denom = coeff_matrix[no_of_CV-1][1] - (coeff_matrix[no_of_CV-1][0]*coeff_matrix[no_of_CV-2][2]);
  coeff_matrix[no_of_CV-1][3] = (coeff_matrix[no_of_CV-1][3] - coeff_matrix[no_of_CV-1][0]*coeff_matrix[no_of_CV-2][3])/denom;
  coeff_matrix[no_of_CV-1][2] = 0.0;
  coeff_matrix[no_of_CV-1][1] = 1.0;
  coeff_matrix[no_of_CV-1][0] = 0.0;

  // Solving for temp by back substitution

  temp[no_of_CV] = coeff_matrix[no_of_CV-1][3];
  for(j=no_of_CV-1;j>0;j--)
    {
      temp[j] = coeff_matrix[j-1][3] - coeff_matrix[j-1][2]*temp[j+1] ;
    }
}


void write_to_file(double *func,double *X,int no_of_CV)
{
	FILE  *fpw;
	char filename[20];
	char no [20];
	strcpy(no,"output_");
	sprintf(no,"%d",no_of_CV);
	strcat(no,".out");
	strcpy(filename,no);
	fpw =fopen(filename,"w");
	fprintf(fpw,"X              	func\n");
	int i=0;
	for(i=0;i<no_of_CV+2;i++)
	{
	  fprintf(fpw,"%lf	%lf\n",X[i],func[i]);	
	}
	fclose(fpw);
	
}
