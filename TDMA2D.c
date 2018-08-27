#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>




void TDMA2D (double **field,int no_of_points_x,int no_of_points_y);
void TDMA(double *field,double *aw,double *ae,double *ap,double *su,int no_of_points);
double **doublealloc2D (int nx, int ny);
void display(double **variable,int nx,int ny);
void difference (double **variable1,double **variable2,double **error,int x,int y);
double *doublealloc1D (int nx);
double max1( double variable1, double variable2);
double max ( double **variable ,int x, int y);
void write_to_file(double **func,int nx,int ny);
void display1D(double * filed, int no_of_points);



double **aw ;
double **ae ;
double **ap ;
double **as ;
double **an ;
double **Su ;
double *aw1,*ae1,*ap1;
double *su;
double *field;
double convergeTDMA = 0.0001;

int main()
{
  
  /*aw = doublealloc2D(10,10);
  ae = doublealloc2D(10,10);
  ap = doublealloc2D(10,10);
  as = doublealloc2D(10,10);
  an = doublealloc2D(10,10);
  Su = doublealloc2D(10,10);*/

 

  ae1 = doublealloc1D(10);
  ap1 = doublealloc1D(10);
  aw1 = doublealloc1D(10);
  su = doublealloc1D(10);
  field =  doublealloc1D(10);
  int i=0,j=0;
  /*for(i=0;i<10;i++)
    {
      for(j=0;j<10;j++)
	{
	  aw[i][j] = 0.0;
	  an[i][j] = 1.0;
	  ap[i][j] = 1.0;
	  as[i][j] = 2.0;
	  ae[i][j] = 0.0;
	  Su[i][j] = 10.0;
	}
    }*/
  for(i=1;i<9;i++)
    {
      aw1[i] = -0.5;
      ap1[i] = 1.0;
      ae1[i] = -0.5;
      su[i]  = 1.0;
    }
  su[1] = su[1] -(aw1[1]*field[0]);
  aw1[1] = 0.0;
  su[8] = su[8] -(ae1[8]*field[9]);
  ae1[8] = 0.0;
  TDMA(field,aw1,ae1,ap1,su,10);
  
  //display(ap,10,10);
  //double **field = doublealloc2D(10,10);
  //display(field,10,10);
  //write_to_file(field,10,10);
  //TDMA2D (field,10,10);
  //display(field,10,10);

  free(aw);
  free(ap);
  free(an);
  free(as);
  free(ae);
  free(Su);
  free(aw1);
  free(ap1);
  free(ae1);
  free(su);
  free(field);

  return 0;
}
  
  







double max ( double **variable ,int x, int y)
{
  int i=0;int j=0;
  double max1 = variable[0][0];
  for(i=0;i<x;i++)
    {
      for(j=0;j<y;j++)
	{
	  if(max1<variable[i][j])
	    max1 = variable[i][j];
	}
    }
  return max1;
}

double max1( double variable1, double variable2)
{
  if(variable1>variable2)
    return variable1;
  else
    return variable2;
}


double **doublealloc2D (int nx, int ny)
{

  double *space;
  double **arr2d;
  int i, j;

  /* first we set aside space for the array itself */

  space = (double *) malloc (nx * ny * sizeof (double));

  /* next we allocate space of an array of pointers, each
     to eventually point to the first element of a
     2 dimensional array of pointers to pointers */

  arr2d = (double **) malloc (nx * sizeof (double *));

  /* and for each of these we assign a pointer to a newly
     allocated array of pointers to a row */

  for (i = 0; i < nx; i++) {
      arr2d[i] = space + i * ny;
    }

  /*  initialising all elements to 0.0 to be on safer side */
  for (i = 0; i < nx; i++) {
      for (j = 0; j < ny; j++) {
	  arr2d[i][j] = 0.0;
	}
    }

  return (arr2d);
}


double *doublealloc1D (int nx)
{

  double *space;
  int i;

  space = (double *) malloc (nx * sizeof (double));

  /*  initialising all elements to 0.0 to be on safer side */
  for (i = 0; i < nx; i++)
    {
      space[i] = 0.0;
    }

  return (space);
}







void difference (double **variable1,double **variable2,double **error,int x,int y)
{
  int i=0,j=0;
  double diff=0.0;
  for( i=0;i<x;i++)
    {
      for(j=0;j<y;j++)
	{
	  diff = variable1[i][j] - variable2[i][j];
	  if(diff>=0.00)
	    {
	      error[i][j] = diff;
	    }
	  else
	    {
	      error[i][j] = -1*diff;
	    }
	}
    }
}

void display(double **variable,int nx,int ny)
{
  int i=0,j=0;
  for(i=0;i<nx;i++)
    {
      for(j=0;j<ny;j++)
	printf( "   %lf",variable[i][j]);
      printf("\n");
    }
  printf("\n\n\n");
}









void TDMA2D (double **field,int no_of_points_x,int no_of_points_y)
{

  double **error_field = doublealloc2D(no_of_points_x,no_of_points_y);
  double **field_old = doublealloc2D(no_of_points_x,no_of_points_y);
  double **field_new = doublealloc2D(no_of_points_x,no_of_points_y);
  double *field1 = doublealloc1D(no_of_points_y);
  double *aw1 = doublealloc1D(no_of_points_y);
  double *ae1 = doublealloc1D(no_of_points_y);
  double *ap1 = doublealloc1D(no_of_points_y);
  double *an1 = doublealloc1D(no_of_points_y);
  double *as1 = doublealloc1D(no_of_points_y);
  double *su1 = doublealloc1D(no_of_points_y);
  /*
  Equation trying to be solved here is :
  a_w T_w + a_e T_e + a_p T_p + a_n T_n + a_s T_s = Su
  */

  int i=0,j=0;
  for(i=0;i<no_of_points_x;i++)
    for(j=0;j<no_of_points_y;j++)
      {
	field_old[i][j] = field[i][j];
	field_new[i][j] = field[i][j];
      }
  int REPEAT = 1;
  while(REPEAT)
    {
      printf(" Enetering loop for %d time ......\n",REPEAT);
      REPEAT++;
      //For thw first sweep since aw will not be defined//
      for(j=0;j<no_of_points_y;j++)
	{
	  aw1[j]  = as[1][j];
	  ap1[j]  = ap[1][j];
	  ae1[j]  = an[1][j];
	  su1[j]  = Su[1][j] -ae[1][j] * field_new[2][j];
	}
      TDMA(field1,aw1,ae1,ap1,su1,no_of_points_y);
      for(j=1;j<no_of_points_y-1;j++)
	{
	  field_new[1][j] = field1[j];
	}
      if(REPEAT ==2)
	{
	  write_to_file(field_new,no_of_points_x,no_of_points_y);
	}
      
      //For the all the interior nodes
      for(i=2;i<no_of_points_x-2;i++)
	{
	  for(j=0;j<no_of_points_y;j++)
	    {
	      aw1[j]  = as[i][j];
	      ap1[j]  = ap[i][j];
	      ae1[j]  = an[i][j];
	      su1[j]  = Su[i][j] - ae[i][j] * field_new[i+1][j] - aw[i][j] *field_new[i-1][j];
	    }
	  TDMA(field1,aw1,ae1,ap1,su1,no_of_points_y);
	  
      
	  for(j=1;j<no_of_points_y-1;j++)
	    {
	      field_new[i][j] = field1[j];
	    }
	  if(REPEAT ==2)
	    {
	      write_to_file(field_new,no_of_points_x,no_of_points_y);
	    }
	}
      //For the final sweep as ae is not defined//
      for(j=0;j<no_of_points_y;j++)
	{
	  aw1[j]  = an[i][j];
	  ap1[j]  = ap[i][j];
	  ae1[j]  = as[i][j];
	  su1[j]  = Su[i][j] - aw[i][j] *field_new[i-1][j];
	}
      TDMA(field1,aw1,ae1,ap1,su1,no_of_points_y);
      for(j=1;j<no_of_points_y-1;j++)
	{
	  field_new[i][j] = field1[j];
	}
      if(REPEAT ==2)
	{
	  write_to_file(field_new,no_of_points_x,no_of_points_y);
	}
      difference(field_old,field_new,error_field,no_of_points_x,no_of_points_y);
      for(i=0;i<no_of_points_x;i++)
	{
	  for(j=0;j<no_of_points_y;j++)
	    {
	      field_old[i][j] = field_new[i][j];
	    }
	}
      if(max(error_field,no_of_points_x,no_of_points_y)<convergeTDMA)
	 REPEAT = 0;
    }
  for(i=0;i<no_of_points_x;i++)
    for(j=0;j<no_of_points_y;j++)
       field[i][j] = field_old[i][j] ;

  free(error_field);
  free(field_old);
  free(field_new);
  free(field1);
}


void TDMA(double *field,double *aw,double *ae,double *ap,double *su,int no_of_points)
{
  int i=0;
  int j=0;
  double denom;
  display1D(ap,no_of_points);
  display1D(aw,no_of_points);
  display1D(ae,no_of_points);
  display1D(su,no_of_points);
  su[1]/= ap[1];
  ae[1]/= ap[1];
  ap[1]/= ap[1];
  aw[1]/= ap[1];
  for( i=2;i<no_of_points-2;i++)
    {
      denom =(ap[i] -(aw[i]*ae[i-1]));
      su[i] = (su[i] - (su[i-1]*aw[i]))/ denom;
      ae[i] = ae[i]/denom;
      ap[i] = 1.0;
      aw[i] = 0;;
    }
  denom = ap[no_of_points-2] - (aw[no_of_points-2]*ae[no_of_points-3]);
  su[no_of_points-2] = (su[no_of_points-2] - aw[no_of_points-2]*su[no_of_points-3])/denom;
  ae[no_of_points-2] = 0.0;
  ap[no_of_points-2] = 1.0;
  aw[no_of_points-2] = 0.0;

  // Solving for field by back substitution

  field[no_of_points-2] = su[no_of_points-2];
  for(j=no_of_points-3;j>0;j--)
    {
      field[j] = su[j] - ae[j]*field[j+1] ;
    }
  display1D(field,no_of_points);
}


void write_to_file(double **func,int nx,int ny)
{
	FILE  *fpw;
	char filename[20];
	char no [20];
	strcpy(no,"iter1");
	sprintf(no,"%d",nx*ny);
	strcat(no,".out");
	strcpy(filename,no);
	fpw =fopen(filename,"a");
	int i=0;
	int j =0;
	for(i=0;i<nx;i++)
	{
	  for(j=0;j<ny;j++)
	    {
	      fprintf(fpw,"%.2lf    ",func[i][j]);
	    }
	  fprintf(fpw,"\n");
	}
	fprintf(fpw,"\n\n\n");
	fclose(fpw);
	
}

void display1D(double * field, int no_of_points)
{
  int i =0;
  for(i=0;i<no_of_points;i++)
    {
      printf("%lf   ",field[i]);
    }
  printf("\n");
}

