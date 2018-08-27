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


double domain[2][2];
int no_of_CV_x,no_of_CV_y;
double temp_domain0;
double temp_domain1;
double tau;
double area;
double **temp,**tempold;
double **Pnew,**Pold;
double **Unew,**Uold;
double **Vnew,**Vold;
double **Ustar,**Vstar;
double *X,*deltaX;
double *Y,*deltaY;
double initialtemp,initialP;
double delx;
double dely;
double delt;
double stop_time;
double rho;
double convergence;
double convergeP,convergeU,convergeV,convergeTDMA;
double Prelax,Urelax,Vrelax;

double Ttop,Tbottom,Tleft,Tright;
double Utop,Ubottom,Vleft,Vright;

double **ap,**ae,**aw,**an,**as,**Su;
int loop_no=0;
int loopno=0;


double max ( double **variable ,int x, int y);
double max1 ( double variable1, double variable2);
double **doublealloc2D(int rows,int columns);
void difference (double **variable1,double **variable2,double **error,int x,int y);
double *doublealloc1D (int nx);
int UPWIND(int field_type);
int CENTRAL(int field_type);
void  initialise(int no_of_CV_x,int no_of_CV_y,double domain[2][2]);
void boundary_velocity();
void boundary_temperature();
void TDMA2D (double **field,int no_of_points_x,int no_of_points_y);
void guass_jacobi (double **field,int no_of_points_x,int no_of_points_y);
void TDMA(double *field,double *aw,double *ae,double *ap,double *su,int no_of_CV);
void write_to_file(double **func,int nx,int ny);
void freecoeff();
void display(double **variable,int nx,int ny);
void write_to_file1(double **variable,int nx,int ny);
void display1D( double *variable, int np);
void amax(double **variable, int nx,int ny);
void TEMP();
double max1d(double *variable,  int nx);
int  max1di(double *variable,int nx);
double diverge();
void piso();
void u_momentum();
void v_momentum();
void  plapalacian();
void write_all(int t);






int main()
{
  ///Input parameters ////////
  domain[0][0] = 0.0;
  domain[0][1] = 1.0;
  domain[1][0] = 0.0;
  domain[1][1] = 1.0;
  stop_time = 0.5;
  delt = 0.005;
  convergeP = 0.1;
  convergence = 0.01;
  convergeU = 0.1;
  convergeV = 0.1;
  convergeTDMA = 0.01;
  Ttop = 0.0;
  Tbottom = 100.0;
  Tleft = 100;
  Tright = 0.0;
  Utop = 2.0;
  Ubottom = 0.0;
  Vleft = 0.0;
  Vright = 0.0;
  area = 1.0;
  tau = 0.01;
  rho  = 1.0;
  no_of_CV_x = 33;
  no_of_CV_y = 33;
  initialP = 0.0;
  initialtemp = 50.0;
  Prelax = 0.5;
  Urelax = 0.5;
  Vrelax = 0.5;
  /////////////////////////

  printf("no_of_CV_x : %d\nno_ofCV_y : %d\n",no_of_CV_x,no_of_CV_y);

  printf("Initialising variables....\n");
  initialise(no_of_CV_x,no_of_CV_y,domain);
  printf("Exiting Initialise....\n");
  printf("Calling Boundary variable....\n");
  boundary_velocity();
  printf("Boundary values initialised...\n");
  write_all(0);
  printf("Calling PISO....\n");
  piso();
  //TEMP();
  printf("PISO over .....\n");
  printf("writing into data file ....\n");
  write_all(10);
  //write_to_file();
  printf("written into file...\n");
  printf("Clearing memory.....\n");


  //freeing variables
  free(X);
  free(Y);
  free(deltaX);
  free(deltaY);
  free(Pold);
  free(Pnew);
  free(Vold);
  free(Vnew);
  free(Vstar);
  free(Ustar);
  free(Unew);
  free(Uold);
  free(tempold);
  printf("Memory cleared ..... Terminating program\n");
  return 0;
}



void u_momentum()
{
  int i=0,j=0;
  double convec_term,diffus_term;
  double Unorth,Usouth,Ueast,Uwest,Vnorth,Vsouth;
  for(i=1;i<no_of_CV_y+1;i++)
    for(j=1;j<no_of_CV_x;j++)
      {
	diffus_term = (Uold[i+1][j] + Uold[i-1][j] - 2.0*Uold[i][j] )/(dely*dely);
	diffus_term += (Uold[i][j+1] + Uold[i][j-1] - 2.0*Uold[i][j] )/(delx*delx);
	Ueast = ( Uold[i][j] + Uold[i][j+1] )/2.0;
	Uwest = ( Uold[i][j] + Uold[i][j-1] )/2.0;
	Unorth = (Uold[i][j] + Uold[i+1][j] )/2.0;
	Usouth = (Uold[i][j] + Uold[i-1][j] )/2.0;
	Vnorth = (Vold[i][j] + Vold[i][j+1] )/2.0;
	Vsouth = (Vold[i-1][j] +Vold[i-1][j+1])/2.0; 
	convec_term = ((( Ueast *Ueast) - (Uwest*Uwest))/(delx)) + ((Unorth*Vnorth - Usouth*Vsouth)/dely);
	Ustar[i][j] = Uold[i][j]  + delt*(-convec_term + (tau*diffus_term));
      }
}


void v_momentum()
{
  int i=0,j=0;
  double convec_term,diffus_term;
  double Vnorth,Vsouth,Veast,Vwest,Ueast,Uwest;
  for(i=1;i<no_of_CV_y;i++)
    for(j=1;j<no_of_CV_x+1;j++)
      {
	diffus_term = (Vold[i+1][j] + Vold[i-1][j] - 2.0*Vold[i][j] )/(dely*dely);
	diffus_term += (Vold[i][j+1] + Vold[i][j-1] - 2.0*Vold[i][j] )/(delx*delx);
	Veast = ( Vold[i][j] + Vold[i][j+1] )/2.0;
	Vwest = ( Vold[i][j] + Vold[i][j-1] )/2.0;
	Vnorth = (Vold[i][j] + Vold[i+1][j] )/2.0;
	Vsouth = (Vold[i][j] + Vold[i-1][j] )/2.0;
	Ueast = (Uold[i][j] + Uold[i+1][j] )/2.0;
	Uwest = (Uold[i][j-1] +Uold[i+1][j-1])/2.0; 
	convec_term = ((Ueast*Veast - Uwest*Vwest)/delx) + ((Vnorth*Vnorth - Vsouth*Vsouth)/dely);
	Vstar[i][j] = Vold[i][j]  + delt*(-convec_term + (tau*diffus_term));
      }
}


void  plapalacian()
{
  int i=0,j=0;
  int i1 = 0 ,i2 = 0 , j1 = 0 , j2 = 0 ;
  an = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  as = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  ae = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  aw = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  ap = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  Su = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
	
  for(i=2;i<no_of_CV_y;i++)
    for(j=2;j<no_of_CV_x;j++)
      {
	Su[i][j] = -(((Ustar[i][j-1] - Ustar[i][j])/(delx*delt)) + ((Vstar[i-1][j] - Vstar[i][j])/(dely*delt)));
	ae[i][j] = (1.0/delx)*(1.0/delx);
	aw[i][j] = (1.0/delx)*(1.0/delx);
	an[i][j] = (1.0/dely)*(1.0/dely);
	as[i][j] = (1.0/dely)*(1.0/dely);
	ap[i][j] = -(aw[i][j] + ae[i][j]+ an[i][j] + as[i][j] );
      }
  i1 = 1;
  i2 = no_of_CV_y;
  for(j=2;j<no_of_CV_x;j++)
    {
      Su[i1][j] = -(((Ustar[i1][j-1] - Ustar[i1][j])/(delx*delt)) + ((Vstar[i1-1][j] - Vstar[i1][j])/(dely*delt)));
      ae[i1][j] = (1.0/delx)*(1.0/delx);
      aw[i1][j] = (1.0/delx)*(1.0/delx);
      an[i1][j] = (1.0/dely)*(1.0/dely);
      as[i1][j] = 0.0;
      ap[i1][j] = -(aw[i1][j] + ae[i1][j]+ an[i1][j] + as[i1][j] );


      Su[i2][j] = -(((Ustar[i2][j-1] - Ustar[i2][j])/(delx*delt)) + ((Vstar[i2-1][j] - Vstar[i2][j])/(dely*delt)));
      ae[i2][j] = (1.0/delx)*(1.0/delx);
      aw[i2][j] = (1.0/delx)*(1.0/delx);
      an[i2][j] = 0.0;
      as[i2][j] =(1.0/dely)*(1.0/dely);
      ap[i2][j] = -(aw[i2][j] + ae[i2][j]+ an[i2][j] + as[i2][j] );
    }

  j1 = 1;
  j2 = no_of_CV_x;
  for(i=2;i<no_of_CV_y;i++)
    {
      Su[i][j1] = -(((Ustar[i][j1-1] - Ustar[i][j1])/(delx*delt)) + ((Vstar[i-1][j1] - Vstar[i][j1])/(dely*delt)));
      ae[i][j1] =(1.0/delx)*(1.0/delx);
      aw[i][j1] = 0.0;
      an[i][j1] = (1.0/dely)*(1.0/dely);
      as[i][j1] = (1.0/dely)*(1.0/dely);
      ap[i][j1] = -(aw[i][j1] + ae[i][j1]+ an[i][j1] + as[i][j1] );


      Su[i][j2] = -(((Ustar[i][j2-1] - Ustar[i][j2])/(delx*delt)) + ((Vstar[i-1][j2] - Vstar[i][j2])/(dely*delt)));
      ae[i][j2] = 0.0;
      aw[i][j2] = (1.0/delx)*(1.0/delx);
      an[i][j2] = (1.0/dely)*(1.0/dely);
      as[i][j2] = (1.0/dely)*(1.0/dely);
      ap[i][j2] = -(aw[i][j2]+ an[i][j2] + as[i][j2] );
    }

  //fixing i1,j1 node to be zero
  /*Su[i1][j1] = Uold[i1][j1-1] - Uold[i1][j1] + Vold[i1-1][j1] - Vold[i1][j1];
    aw[i1][j1] = 0.0;
    ae[i1][j1] = DU[i1][j1];
    an[i1][j1] = DV[i1][j1];
    as[i1][j1] = 0.0;
    ap[i1][j1] = -(aw[i1][j1] + ae[i1][j1]+ an[i1][j1] + as[i1][j1] );*/

  Su[i1][j1] = 0.0;
  aw[i1][j1] = 0.0;
  ae[i1][j1] = 0.0;
  an[i1][j1] = 0.0;
  as[i1][j1] = 0.0;
  ap[i1][j1] = 1.0; 
	

  //for i1+1,j1 node

  Su[i1+1][j1] = -(((Ustar[i1+1][j1-1] - Ustar[i1+1][j1])/(delx*delt)) + ((Vstar[i1][j1] - Vstar[i1+1][j1])/(dely*delt)));
  aw[i1+1][j1] = 0.0;
  ae[i1+1][j1] = (1.0/delx)*(1.0/delx);
  an[i1+1][j1] = (1.0/dely)*(1.0/dely);
  as[i1+1][j1] = 0.0;
  ap[i1+1][j1] = -(aw[i1+1][j1] + ae[i1+1][j1]+ an[i1+1][j1] + as[i1+1][j1] );

  //for i1,j1+1 node
  Su[i1][j1+1] = -(((Ustar[i1][j1] - Ustar[i1][j1+1])/(delx*delt)) +(( Vstar[i1-1][j1+1] - Vstar[i1][j1+1])/(dely*delt)));
  aw[i1][j1+1] = 0.0;
  ae[i1][j1+1] = (1.0/delx)*(1.0/delx);
  an[i1][j1+1] = (1.0/dely)*(1.0/dely);
  as[i1][j1+1] = 0.0;
  ap[i1][j1+1] = -(aw[i1][j1+1] + ae[i1][j1+1]+ an[i1][j1+1] + as[i1][j1+1] );

  //for i1,j2 node

  Su[i1][j2] = -((( Ustar[i1][j2-1] - Ustar[i1][j2])/(delx*delt)) + ((Vstar[i1-1][j2] - Vstar[i1][j2])/(dely*delt)));
  ae[i1][j2] = 0.0;
  aw[i1][j2] = (1.0/delx)*(1.0/delx);
  an[i1][j2] = (1.0/dely)*(1.0/dely);
  as[i1][j2] = 0.0;
  ap[i1][j2] = -(aw[i1][j2] + ae[i1][j2]+ an[i1][j2] + as[i1][j2] );

  //for i2,j1 node

  Su[i2][j1] = -(((Ustar[i2][j1-1] - Ustar[i2][j1])/(delx*delt)) + ((Vstar[i2-1][j1] - Vstar[i2][j1])/(dely*delt)));
  ae[i2][j1] = (1.0/delx)*(1.0/delx);
  aw[i2][j1] = 0.0;
  an[i2][j1] = 0.0;
  as[i2][j1] = (1.0/dely)*(1.0/dely);
  ap[i2][j1] = -(aw[i2][j1] + ae[i2][j1]+ an[i2][j1] + as[i2][j1] );

  //for i2,j2 node

  Su[i2][j2] = -(((Ustar[i2][j2-1] - Ustar[i2][j2])/(delx*delt)) + ((Vstar[i2-1][j2] - Vstar[i2][j2])/(dely*delt)));
  ae[i2][j2] = 0.0;
  aw[i2][j2] = (1.0/delx)*(1.0/delx);
  an[i2][j2] = 0.0;
  as[i2][j2] = (1.0/dely)*(1.0/dely);
  ap[i2][j2] = -(aw[i2][j2] + ae[i2][j2]+ an[i2][j2] + as[i2][j2] );

	      
	    
  /*printf("Printing for Pressur\n");
  printf("Printing for ap:\n");
  amax(ap,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for ae:\n");
  amax(ae,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for aw:\n");
  amax(aw,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for an:\n");
  amax(an,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for as:\n");
  amax(as,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for Su:\n");
  amax(Su,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(ap,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(ae,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(aw,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(an,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(as,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(Su,no_of_CV_y+2,no_of_CV_x+2);*/
	

}


double diverge()
{
  int i=0,j=0;
  double max =0 ;
  double difference;
  for(i =1;i<no_of_CV_y+1;i++)
    for(j=1;j<no_of_CV_x+1;j++)
      {
	difference = (Uold[i][j] - Uold[i][j-1])/delx;
	difference += (Vold[i][j] - Vold[i-1][j])/dely;
	if(difference>max)
	  max = difference;
	if(-difference>max)
	  max = -difference;
	
      }

  return max;
}


void piso()
{
  int i=0,j=0;
  int ui,uj;
  int I ,J ;
  int vi,vj;
  int loop=1;
  int t = 1;
  printf("Entered piso \n");
  double time = 0.0;
  while(time < stop_time)
    {
      printf("Time : %.2lf\n",time);
      
      loop = 0;
      //solve for u momentum
      printf("Solving U momentum.....");
      u_momentum();
      for (i =1;i<no_of_CV_x;i++)
	{
	  Ustar[no_of_CV_y+1][i] = ((2.0*Utop) - Ustar[no_of_CV_y][i] ) ;
	  Ustar[0][i] = ((2.0*Ubottom) - Ustar[1][i]) ;
	}
      printf("solved\n");
      //solve for v momentum
      printf("Solving V momentum.....");
      v_momentum();
      for (i =1;i<no_of_CV_y;i++)
	{
	  Vstar[i][0] = ((2.0*Vleft) -  Vstar[i][1] );
	  Vstar[i][no_of_CV_x+1] = ((2.0* Vright) - Vstar[i][no_of_CV_x] );
	}
      printf("solved\n");
      //solve for pressure
      // do
      //	{
	  printf("Entering loop for %d time\n",loop);
	  loop++;
	  printf("Solving P momentum.....");
	  plapalacian();
	  //guass_jacobi (Pold,no_of_CV_y+2,no_of_CV_x+2);
	  TDMA2D(Pold,no_of_CV_y+2,no_of_CV_x+2);
	  //boundary for pressure
	  for(I=0;I<no_of_CV_y+2;I++)
	    {
	      Pold[I][0] = Pold[I][1];
	      Pold[I][no_of_CV_x+1] = Pold[I][no_of_CV_x];
	    }
	  for(J=0;J<no_of_CV_x+2;J++)
	    {
	      Pold[0][J] = Pold[1][J];
	      Pold[no_of_CV_y+1][J] = Pold[no_of_CV_y][J];
	    }
	  printf("solved\n");

	  //updating Unew
	  printf("Updating U.....\n");
	  //printf("1");
	  for (ui=1;ui<no_of_CV_y+1;ui++)
	    {
	      for(uj=1;uj<no_of_CV_x;uj++)
		{
		  //printf(" ui : %d " , ui);
		  //printf(" uj : %d\n " , uj);
		  Unew[ui][uj] = Ustar[ui][uj] + (delt* (Pold[ui][uj] - Pold[ui][uj+1])/delx);
		  Uold[ui][uj] = Unew[ui][uj];
		  Ustar[ui][uj] = Uold[ui][uj];
		}
	    }
	  //printf("2");
	  for (i =1;i<no_of_CV_x;i++)
	    {
	      //printf(" i : %d " , i);
	      Uold[no_of_CV_y+1][i] = ((2.0*Utop) - Uold[no_of_CV_y][i] ) ;
	      Uold[0][i] = ((2.0*Ubottom) - Uold[1][i]) ;
	      Unew[no_of_CV_y+1][i] = ((2.0*Utop) - Unew[no_of_CV_y][i] ) ;
	      Unew[0][i] = ((2.0*Ubottom) - Unew[1][i]) ;
	      Ustar[no_of_CV_y+1][i] = ((2.0*Utop) - Ustar[no_of_CV_y][i] ) ;
	      Ustar[0][i] = ((2.0*Ubottom) - Ustar[1][i]) ;
	      //printf("doop");
	    }
	  //printf("klad");

	  //updating Vnew
	  /*printf("Updating V.....\n");
	  printf("sasdf");
	  printf("3");*/
	  for (vi=1;vi<no_of_CV_y;vi++)
	    {
	      for(vj=1;vj<no_of_CV_x+1;vj++)
		{
       		  //printf(" vi : %d " , vi);
		  //printf(" vj : %d\n " , vj);
		  Vnew[vi][vj] = Vstar[vi][vj] + (delt* ( Pold[vi][vj] - Pold[vi+1][vj])/dely);
		  Vold[vi][vj] = Vnew[vi][vj];
		  Vstar[vi][vj] = Vold[vi][vj];
		  //printf("Outofloop");
	      
		}
	      //printf("out of loop ");
	    }
	  //printf("Before 4 ");
	  //printf("4");
	  for (i =1;i<no_of_CV_y;i++)
	    {
	      //printf(" i : %d " , i);
	      Vold[i][0] = ((2.0*Vleft) -  Vold[i][1] );
	      Vold[i][no_of_CV_x+1] = ((2.0* Vright) - Vold[i][no_of_CV_x] );
	      Vnew[i][0] = ((2.0*Vleft) -  Vnew[i][1] );
	      Vnew[i][no_of_CV_x+1] = ((2.0* Vright) - Vnew[i][no_of_CV_x] );
	      Vstar[i][0] = ((2.0*Vleft) -  Vstar[i][1] );
	      Vstar[i][no_of_CV_x+1] = ((2.0* Vright) - Vstar[i][no_of_CV_x] );
	    }
	  
	  freecoeff();
	  //printf("Diveregnce  : %.2lf\n ", diverge());
	  //	}while(diverge() > convergence);
      //write_all(t);
      t++;
      time  = time + delt;
    }
}
  
  




	      
void  initialise(int no_of_CV_x,int no_of_CV_y,double domain[2][2])
{
  
  int i =0;int j=0;
  X  =doublealloc1D(no_of_CV_x+2);
  Y  = doublealloc1D(no_of_CV_y+2);
  deltaX =doublealloc1D(no_of_CV_x+1);
  deltaY = doublealloc1D(no_of_CV_y+1);
  delx = ( ( domain[0][1]  - domain[0][0]) / no_of_CV_x );
  dely =  ( ( domain[1][1]  - domain[1][0]) / no_of_CV_y );
  X[0] = domain[0][0];
  Y[0] = domain[1][0];
  X[1] = domain [0][0] + (delx/2);
  Y[1] = domain[1][0] + (dely/2);
  deltaX[0] = delx/2;
  deltaY[0] = dely/2;
  for (i=2;i<=no_of_CV_x;i++)
    {
      X[i] = X[i-1] + delx ;
      deltaX[i-1] = delx;
    }
  X[no_of_CV_x+1] = X[no_of_CV_x] + (delx/2);
  deltaX[no_of_CV_x] = delx/2;
  for (i=2;i<=no_of_CV_y;i++)
    {
      Y[i] = Y[i-1] + dely ;
      deltaY[i-1] = dely;
    }
  Y[no_of_CV_y+1] = Y[no_of_CV_y] + (dely/2);
  deltaY[no_of_CV_y] = dely/2;

  printf("delx : %.2lf\ndely : %.2lf\n",delx,dely);

  Uold = doublealloc2D(no_of_CV_y+2,no_of_CV_x+1);
  Unew = doublealloc2D(no_of_CV_y+2,no_of_CV_x+1);
  Ustar = doublealloc2D(no_of_CV_y+2,no_of_CV_x+1);
  Vold = doublealloc2D(no_of_CV_y+1,no_of_CV_x+2);
  Vnew = doublealloc2D(no_of_CV_y+1,no_of_CV_x+2);
  Vstar = doublealloc2D(no_of_CV_y+1,no_of_CV_x+2);
  tempold = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  Pold = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  Pnew = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  for(i=0;i<no_of_CV_y+2;i++)
    for(j=0;j<no_of_CV_x+2;j++)
      {
	Pold[i][j] = initialP;
	tempold[i][j] = initialtemp;
      }
}

void boundary_velocity()
{
  //Top surface has some velocty
  //all the other surfaces have no slip condition

  //All the sides are mailtained at fixed temperaure
  //top surface

  int i=0,j=0;
  //for u velocity

  //top walll and bottom wall
  for (i =0;i<no_of_CV_x+1;i++)
    {
      Uold[no_of_CV_y+1][i] = Utop;
      Uold[no_of_CV_y][i]   = Utop;
      Uold[0][i] = Ubottom;
      Uold[1][i] = Ubottom;
    }
  //for the side walls
  for (i =1;i<no_of_CV_y+1;i++)
    {
      Uold[i][0] =0.0;
      Uold[i][no_of_CV_x] = 0.0;
    }
  
  //for v velocity

  //top walll and bottom wall
  for (i =0;i<no_of_CV_x+2;i++)
    {
      Vold[no_of_CV_y][i] =0.0;
      Vold[0][i] = 0.0;
    }
  //for the side walls
  for (i =1;i<no_of_CV_y;i++)
    {
      Vold[i][0] = Vleft;
      Vold[i][1] = Vleft;
      Vold[i][no_of_CV_x+1] = Vright;
      Vold[i][no_of_CV_x]   = Vright;
    }
  
}

void boundary_temperature()
{
  int i=0,j=0;
  //For the top and bottom wall
  for(i=0;i<no_of_CV_x+2;i++)
    {
      tempold[no_of_CV_y+1][i] =Ttop;
      tempold[0][i] = Tbottom;
    }
  for(j=0;j<no_of_CV_y+2;j++)
    {
      tempold[j][0] = Tleft;
      tempold[j][no_of_CV_x+1]  =Tright;
    }

}

void DIFFUS()
{
  int i=0;int j=0;
  an = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  as = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  ae = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  aw = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  ap = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  Su = doublealloc2D(no_of_CV_y+2,no_of_CV_x+2);
  double De,Dn,Ds,Dw;
  double Fe,Fw,Fn,Fs;

  for(i=1;i<no_of_CV_y+1;i++)
    for(j=1;j<no_of_CV_x+1;j++)
      {
	De = tau/deltaX[j];
	Dw = tau/deltaX[j-1];
	Dn = tau/deltaY[i];
	Ds = tau/deltaY[i-1];
	Fe = rho*Uold[i][j];
	Fw = rho*Uold[i][j-1] ;
	Fn = rho*Vold[i][j] ;
	Fs = rho*Vold[i-1][j];
	ae[i][j] = -(De - (Fe/2.0));
	aw[i][j] = -(Dw + (Fw/2.0));
	as[i][j] = -(Ds - (Fs/2.0));
	an[i][j] = -(Dn + (Fn/2.0));
	ap[i][j] = -(aw[i][j] + ae[i][j]+ an[i][j] + as[i][j] ) + (Fe - Fw) + (Fn - Fs);
	if(i==1)
	  {
	    if(j==1)
	      {
		Su[i][j] = 0.0 - as[i][j]*tempold[i-1][j] - aw[i][j]*tempold[i][j-1] ;
		aw[i][j] = 0.0;
		as[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	    else if(j==no_of_CV_x)
	      {
		Su[i][j] = 0.0 - ae[i][j]*tempold[i][j+1]- as[i][j]*tempold[i-1][j] ;
		ae[i][j] = 0.0;
		as[i][j] = 0.0;

		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	    else
	      {
		Su[i][j] = 0.0 - as[i][j]*tempold[i-1][j] ;
		as[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	  }
	else if(i==no_of_CV_x)
	  {
	    if(j==1)
	      {
		Su[i][j] = 0.0 - aw[i][j]*tempold[i][j-1] - an[i][j]*tempold[i+1][j];
		aw[i][j] = 0.0;
		an[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	    else if(j==no_of_CV_y)
	      {
		Su[i][j] = 0.0 - an[i][j]*tempold[i+1][j] - ae[i][j]*tempold[i][j+1] ;
		ae[i][j] = 0.0;
		an[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	    else
	      {
		Su[i][j] = 0.0 - an[i][j]*tempold[i+1][j] ;
		an[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	  }
	else
	  {
	    if(j==1)
	      {
		Su[i][j] = 0.0 - aw[i][j]*tempold[i][j-1];
		aw[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	    else if(j==no_of_CV_y)
	      {
		Su[i][j] =0.0 - ae[i][j]*tempold[i][j+1];
		ae[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
	    else
	      {
		Su[i][j] = 0.0;
		if((ap[i][j] < 0.00001)&&(ap[i][j] > -  0.00001))
		  printf("Zero encountered here\n");
	      }
		 
	  }
	      
      }
  write_to_file1(ap,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(ae,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(aw,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(an,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(as,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file1(Su,no_of_CV_y+2,no_of_CV_x+2);
	
  
      
  
  printf("Printing for Temperature\n");
  printf("Printing for ap:\n");
  amax(ap,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for ae:\n");
  amax(ae,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for aw:\n");
  amax(aw,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for an:\n");
  amax(an,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for as:\n");
  amax(as,no_of_CV_y+2,no_of_CV_x+2);
  printf("Printing for Su:\n");
  amax(Su,no_of_CV_y+2,no_of_CV_x+2);
	

}


void TEMP()
{
  int i=0,j=0;
  for(i=0;i<no_of_CV_y+2;i++)
    for(j=0;j<no_of_CV_x+1;j++)
      Uold[i][j] = 1.0;
  for(i=0;i<no_of_CV_y+1;i++)
    for(j=0;j<no_of_CV_x+2;j++)
      Vold[i][j] = 1.0;
  boundary_temperature();
  write_to_file(tempold,no_of_CV_y+2,no_of_CV_x+2);
  DIFFUS();
  TDMA2D(tempold,no_of_CV_y+2,no_of_CV_x+2);
  write_to_file(tempold,no_of_CV_y+2,no_of_CV_x+2);
}

void guass_jacobi (double **field,int no_of_points_x,int no_of_points_y)
{

  double **error_field = doublealloc2D(no_of_points_x,no_of_points_y);
  double **field_old = doublealloc2D(no_of_points_x,no_of_points_y);
  double **field_new = doublealloc2D(no_of_points_x,no_of_points_y);
  double *field1 = doublealloc1D(no_of_points_y);
  double diff = 0.0;
  double max1 = 0.0;
  int I,J;
  int g;
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
  do
    {
      printf(" Enetering GS for %d time ......\n",REPEAT);
      REPEAT++;
      //For thw first sweep since aw will not be defined//

      for(i=1;i<no_of_points_x-1;i++)
	for(j=1;j<no_of_points_y-1;j++)
	  {
	    field_new[i][j] = Su[i][j] - aw[i][j]*field_old[i][j-1] - ae[i][j]*field_old[i][j+1] - an[i][j]*field_old[i+1][j] - as[i][j]*field_old[i-1][j];
	    field_new[i][j] = field_new[i][j]/ap[i][j];
	  }
     

       difference(field_old,field_new,error_field,no_of_points_x,no_of_points_y);

       for(i=0;i<no_of_points_x;i++)
	 for(j=0;j<no_of_points_y;j++)
	   field_old[i][j] = field_new[i][j] ;

      /*if(REPEAT >60000)
	REPEAT = 0;*/
      printf(" MAx error : %lf\n",max(error_field,no_of_points_x,no_of_points_y));
    }while(max(error_field,no_of_points_x,no_of_points_y)>convergeTDMA);
  for(i=0;i<no_of_points_x;i++)
    for(j=0;j<no_of_points_y;j++)
       field[i][j] = field_old[i][j] ;

  free(error_field);
  free(field_old);
  free(field_new);
  free(field1);
}

  

  


void TDMA2D (double **field,int no_of_points_x,int no_of_points_y)
{

  double **error_field = doublealloc2D(no_of_points_x,no_of_points_y);
  double **field_old = doublealloc2D(no_of_points_x,no_of_points_y);
  double **field_new = doublealloc2D(no_of_points_x,no_of_points_y);
  double *field1 = doublealloc1D(no_of_points_y);
  double *aw1 = doublealloc1D(no_of_points_x);
  double *ae1 = doublealloc1D(no_of_points_x);
  double *ap1 = doublealloc1D(no_of_points_x);
  double *an1 = doublealloc1D(no_of_points_x);
  double *as1 = doublealloc1D(no_of_points_x);
  double *su1 = doublealloc1D(no_of_points_x);
  double diff = 0.0;
  double max1 = 0.0;
  int I,J;
  int g;
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
      printf(" Enetering TDMA for %d time ......\n",REPEAT);
      REPEAT++;
      //For thw first sweep since aw will not be defined//
      for(j=0;j<no_of_points_x;j++)
	{
	  aw1[j]  = as[j][1];
	  ap1[j]  = ap[j][1];
	  ae1[j]  = an[j][1];
	  su1[j]  = Su[j][1] -ae[j][1] * field_old[j][2];
	}
      /* max1  = max1d(su1,no_of_points_x);
      if(max1>200)
	{
	  J = max1di(su1,no_of_points_x);
	  I = 0;
	  printf("Su1 is large here : %d  %d",I,J);
	  printf("Value of aw : %lf, field : %lf\n",aw[J][I],field_old[J][I-1]);
	  printf("Value of ae : %lf, field : %lf\n",ae[J][I],field_old[J][I+1]);
	  //scanf("%d",&g);
	  }*/
      /* printf("ap : \n");
      display1D(ap1,no_of_points_x);
      printf("ae : \n");
      display1D(ae1,no_of_points_x);
      printf("aw : \n");
      display1D(aw1,no_of_points_x);
      printf("su : \n");
      display1D(su1,no_of_points_x);*/
      
      TDMA(field1,aw1,ae1,ap1,su1,no_of_points_x);
      for(j=1;j<no_of_points_x-1;j++)
	{
	  field_new[j][1] = field1[j];
	}
      
      //For the all the interior nodes
      for(i=2;i<no_of_points_y-2;i++)
	{
	  for(j=0;j<no_of_points_x;j++)
	    {
	      aw1[j]  = as[j][i];
	      ap1[j]  = ap[j][i];
	      ae1[j]  = an[j][i];
	      su1[j]  = Su[j][i] - ae[j][i] * field_old[j][i+1] - aw[j][i] *field_new[j][i-1];
	     
	    }
	  /*printf("ap : \n");
	  display1D(ap1,no_of_points_x);
	  printf("ae : \n");
	  display1D(ae1,no_of_points_x);
	  printf("aw : \n");
	  display1D(aw1,no_of_points_x);
	  printf("su : \n");
	  display1D(su1,no_of_points_x);*/
	  /* max1  = max1d(su1,no_of_points_x);
	  if(max1>20000)
	    {
	      J = max1di(su1,no_of_points_x);
	      I = i;
	      printf("su2 is large here : %d  %d   %lf\n",I,J,su1[J]);
	      printf("Value of aw : %lf, field : %lf\n",aw[J][I],field_old[J][I-1]);
	      printf("Value of ae : %lf, field : %lf\n",ae[J][I],field_old[J][I+1]);
	      //scanf("%d",&g);
	      }*/
	  TDMA(field1,aw1,ae1,ap1,su1,no_of_points_x);
	  
      
	  for(j=1;j<no_of_points_x-1;j++)
	    {
	      field_new[j][i] = field1[j];
	    }
	}
      //For the final sweep as ae is not defined//
      for(j=0;j<no_of_points_x;j++)
	{
	  aw1[j]  = as[j][i];
	  ap1[j]  = ap[j][i];
	  ae1[j]  = an[j][i];
	  su1[j]  = Su[j][i] - aw[j][i] *field_new[j][i-1];
	}
      /*max1  = max1d(su1,no_of_points_x);
      if(max1>20000)
	{
	  J = max1di(su1,no_of_points_x);
	  I = i;
	  printf("Su3 is large here : %d  %d",I,J);
	  printf("Value of aw : %lf, field : %lf\n",aw[J][I],field_old[J][I-1]);
	  printf("Value of ae : %lf, field : %lf\n",ae[J][I],field_old[J][I+1]);
	  //scanf("%d",&g);
	  }*/
      /* printf("ap : \n");
      display1D(ap1,no_of_points_x);
      printf("ae : \n");
      display1D(ae1,no_of_points_x);
      printf("aw : \n");
      display1D(aw1,no_of_points_x);
      printf("su : \n");
      display1D(su1,no_of_points_x);*/
      TDMA(field1,aw1,ae1,ap1,su1,no_of_points_x);
      for(j=1;j<no_of_points_x-1;j++)
	{
	  field_new[j][i] = field1[j];
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
      /*if(REPEAT >60000)
	REPEAT = 0;*/
      printf(" MAx error : %lf\n",max(error_field,no_of_points_x,no_of_points_y));
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
  su[1]/= ap[1];
  ae[1]/= ap[1];
  ap[1]/= ap[1];
  aw[1]/= ap[1];
  for( i=2;i<no_of_points-2;i++)
    {
      denom =(ap[i] -(aw[i]*ae[i-1]));
      if((denom<0.0001)&&(denom>(-0.0001)))
	printf(".......OOOOPPPPPPPS  DENOM BECAME ZERO HERE VERY PATHETIC ............\n");
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
}

void freecoeff()
{
  free(an);
  free(as);
  free(ae);
  free(aw);
  free(ap);
  free(Su);
}


void write_to_file(double **func,int nx,int ny)
{
  FILE  *fpw;
  char filename[20];
  char no [20];
  strcpy(no,"out");
  strcat(no,".out");
  strcpy(filename,no);
  fpw =fopen(filename,"a");
  fprintf(fpw,"i      j    func\n");
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


double max ( double **variable ,int x, int y)
{
  int i=0;int j=0;
  double max1 = 0.0;
  for(i=0;i<x;i++)
    {
      for(j=0;j<y;j++)
	{
	  max1 += variable[i][j];
	}
    }
  max1 = max1/(x*y);
  /*printf(" MAx error : %lf     ",max1);
  for(i=0;i<x;i++)
    {
      for(j=0;j<y;j++)
	{
	  if(max1 ==variable[i][j])
	    printf( " i = %d|| j = %d\n ",i,j);
	}
	}*/
  
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


void write_to_file1(double **variable,int nx,int ny)
{
  int i=0,j=0;
  FILE  *fpw;
  char filename[20];
  char no [20];
  strcpy(no,"var");
  strcat(no,".out");
  strcpy(filename,no);
  fpw =fopen(filename,"a");
  for(i=0;i<nx;i++)
    {
      for(j=0;j<ny;j++)
	{
	  fprintf(fpw,"%.2lf    ",variable[i][j]);
	}
      fprintf(fpw,"\n");
    }
  fprintf(fpw,"\n\n\n");
  fclose(fpw);
}

void display1D( double *variable, int np)
{
  int i=0;
  for(i=0;i<np;i++)
    printf("%.2lf   ",variable[i]);
  printf("\n\n");
}


void amax(double **variable, int nx,int ny)
{
  double max;
  double min;
  max = variable[1][1];
  min = variable[1][1];
  int i=0,j=0;
  int zero =0;
  for(i=1;i<nx-1;i++)
    {
    for(j=1;j<ny-1;j++)
      {
	if(max<variable[i][j])
	  max = variable[i][j];
	if(min>variable[i][j])
	  min = variable[i][j];
	if((variable[i][j]>(-0.0001))&&(variable[i][j]<0.0001))
	  zero = 1;
	
      }
    }
    
  printf(" Max value = %.3lf       Min value = %.3lf\n",max,min);
  if(zero)
    printf("Zero encounteres\n");
  else
    printf("No zeroes found\n");
  
}

int  max1di(double *variable,int nx)
{
  double max;
  double min;
  int I = 0;
  max = variable[1];
  min = variable[1];
  int i=0;
  int zero =0;
  for(i=1;i<nx-1;i++)
    {
      if(max < variable[i])
	{
	  max = variable[i];
	  I = i;
	}
    }
  return I;
	

}

double max1d(double *variable,  int nx)
{
  double max = variable[1];
  int i=0;
  for(i=1;i<nx-1;i++)
    {
      if(max<variable[i])
	max  = variable[i];
    }
  return max;
  
}

void write_all(int t)
{
  int i=0,j=0;
  FILE  *fpw;
  char filename[20];
  char no [20];
  char time[10];
  sprintf (filename, "velocity%d.txt", t);
  fpw =fopen(filename,"w");

  //for U
  fprintf(fpw,"U momentum\n");
  for(i=0;i<no_of_CV_y+2;i++)
    {
      for(j=0;j<no_of_CV_x+1;j++)
	{
	  fprintf(fpw,"%.2lf    ",Uold[i][j]);
	}
      fprintf(fpw,"\n");
    }


  fprintf(fpw,"\n");fprintf(fpw,"\n");fprintf(fpw,"\n");
  //for V
  fprintf(fpw,"V momentum\n");
  for(i=0;i<no_of_CV_y+1;i++)
    {
      for(j=0;j<no_of_CV_x+2;j++)
	{
	  fprintf(fpw,"%.2lf    ",Vold[i][j]);
	}
      fprintf(fpw,"\n");
    }

  //For P
  fprintf(fpw,"P momentum\n");
  fprintf(fpw,"\n");fprintf(fpw,"\n");fprintf(fpw,"\n");

  for(i=0;i<no_of_CV_y+2;i++)
    {
      for(j=0;j<no_of_CV_x+2;j++)
	{
	  fprintf(fpw,"%.2lf    ",Pold[i][j]);
	}
      fprintf(fpw,"\n");
    }
  fprintf(fpw,"\n\n\n");
  fclose(fpw);
}
  

  
