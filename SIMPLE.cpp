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

#include<iostream>
#include<cmath>
#include <fstream>
using namespace std;


const int gx = 100;
const int gy = 1000;
double rho=1.0;
double delx=.01;
double mu=.001;
double max_err=0;
double Po=10;
double Uo=.1;
int I1 = 35,I2 = 65;
int J1 = 220,J2 = 280 ;


double u[gx+2][gy+2],v[gx+2][gy+2],su[gx+2][gy+2],p[gx+2][gy+2],Fe,Fw,Fs,Fn,aw[gx+2][gy+2],ap[gx+2][gy+2],ae[gx+2][gy+2],an[gx+2][gy+2],as[gx+2][gy+2],D;
double cP[gx+2][gy+2],cW[gx+2][gy+2],cE[gx+2][gy+2],cN[gx+2][gy+2],cS[gx+2][gy+2],cSU[gx+2][gy+2],Du[gx+2][gy+2],Dv[gx+2][gy+2],p_err[gx+2][gy+2];



int TDMA2D (int option);
void TDMA(double *field,double *aw,double *ae,double *ap,double *su,int start,int end);
void gauss_siedel(int option);
double max(double variable1,double variable2)
{
  double max1;
  if(variable1>variable2)
    max1  = variable1;
  else
    max1 = variable2;
  return max1;
}
void write_to_vtk();

float convection_diffusive_flux_u(int i,int j)// to assign values to corresponding fluxes accepts corresponding values of i,j
{

 
  Fe=(u[i][j]+u[i][j+1])*0.5*rho;

  Fw=((u[i][j]+u[i][j-1]))*0.5*rho;
   
  Fs=((v[i][j-1]+v[i][j]))*0.5*rho;
     
  Fn=((v[i+1][j-1]+v[i+1][j]))*0.5*rho;

  D=mu/delx;

      



  return(0);
}

float convection_diffusive_flux_v(int i,int j)// to assign values to corresponding fluxes accepts corresponding values of i,j
{

 
  Fe=(u[i][j+1]+u[i][j+1])*0.5*rho;

  Fw=((u[i][j]+u[i-1][j]))*0.5*rho;
   
  Fs=((v[i-1][j]+v[i][j]))*0.5*rho;
     
  Fn=((v[i][j]+v[i+1][j]))*0.5*rho;

  D=mu/delx;

      



  return(0);
}

int coeff_U()
{
  int i=0,j=0;
  for(i=1;i<=gx;i++)
    {
      for(j=2;j<=gy;j++)
	{
	  convection_diffusive_flux_u(i,j);
	  ae[i][j] = -max(0.0,-Fe) - D ;
	  aw[i][j] = -max(0.0,Fw) - D;
	  an[i][j] = Fn/2.0 - D;
	  as[i][j] = -Fs/2.0 - D;
	  ap[i][j] = -(ae[i][j] + aw[i][j] + an[i][j] + as[i][j]) + (Fe - Fw) + (Fn - Fs);
	  su[i][j] = (p[i][j-1]-p[i][j]);

	} 
    }
  // for right boundary
  for (i = 1;i<=gx;i++)
    {
      // right boundary
      convection_diffusive_flux_u(i,gy);
      ae[i][gy] = -max(0.0,-Fe) - D ;
      aw[i][gy] = -max(0.0,Fw) - D;
      an[i][gy] = Fn/2.0 - D;
      as[i][gy] = -Fs/2.0 - D;
      ap[i][gy] = -( aw[i][gy] + an[i][gy] + as[i][gy]) + (Fe - Fw) + (Fn - Fs);
      su[i][gy] = (p[i][gy-1]-p[i][gy]);
      ae[i][gy] = 0.0;


      //  left boundary
      convection_diffusive_flux_u(i,2);
      ae[i][2] = -max(0.0,-Fe) - D ;
      aw[i][2] = -max(0.0,Fw) - D;
      an[i][2] = Fn/2.0 - D;
      as[i][2] = -Fs/2.0 - D;
      ap[i][2] = -(aw[i][2]+ ae[i][2] + an[i][2] + as[i][2]) + (Fe - Fw) + (Fn - Fs);
      su[i][2] = (p[i][1]-p[i][2]) - (aw[i][2]*Uo);
           
      aw[i][2] = 0.0;


    }

  for (j = 2;j<=gy;j++)
    {
      // top wall
      convection_diffusive_flux_u(gx,j);
      ae[gx][j] = -max(0.0,-Fe) - D ;
      aw[gx][j] = -max(0.0,Fw) - D;
      an[gx][j] = Fn/2.0 - D;
      as[gx][j] = -Fs/2.0 - D;
      ap[gx][j] = -( 2.0*an[gx][j] + aw[gx][j] + ae[gx][j] + as[gx][j]) + (Fe - Fw) + (Fn - Fs);
      su[gx][j] = (p[gx][j-1]-p[gx][j]);
      an[gx][j] = 0.0;


      //  bottom boundary
      convection_diffusive_flux_u(1,j);
      ae[1][j] = -max(0.0,-Fe) - D ;
      aw[1][j] = -max(0.0,Fw) - D;
      an[1][j] =  Fn/2.0 - D;
      as[1][j] = -Fs/2.0 - D;
      ap[1][j] = -(aw[1][j] + ae[1][j] + an[1][j] + 2.0*as[1][j]) + (Fe - Fw) + (Fn - Fs);
      su[1][j] = (p[1][j-1]-p[1][j]) ;
      as[1][j] = 0.0;


    }

  // for 1,2  node
  convection_diffusive_flux_u(1,2);
  ae[1][2] = -max(0.0,-Fe) - D ;
  aw[1][2] = -max(0.0,Fw) - D;
  an[1][2] =  Fn/2.0 - D;
  as[1][2] = -Fs/2.0 - D;
  ap[1][2] = -(ae[1][2]+aw[1][2]+ an[1][2] + 2.0*as[1][2]) + (Fe - Fw) + (Fn - Fs);
  su[1][2] = (p[1][1]-p[1][2])-(Uo*aw[1][2]) ;
  as[1][2] = 0.0;
  aw[1][2] = 0.0;

  // for gx,2 node
  convection_diffusive_flux_u(gx,2);
  ae[gx][2] = -max(0.0,-Fe) - D ;
  aw[gx][2] = -max(0.0,Fw) - D;
  an[gx][2] = Fn/2.0 - D;
  as[gx][2] = -Fs/2.0 - D;
  ap[gx][2] = -( 2.0*an[gx][2]+aw[gx][2]+ ae[gx][2] + as[gx][2]) + (Fe - Fw) + (Fn - Fs);
  su[gx][2] = (p[gx][1]-p[gx][2])-(Uo*aw[gx][2]);
  an[gx][2] = 0.0;
  aw[gx][2] = 0.0;

  // for 1,gy node
  convection_diffusive_flux_u(1,gy);
  ae[1][gy] = -max(0.0,-Fe) - D ;
  aw[1][gy] = -max(0.0,Fw) - D;
  an[1][gy] = Fn/2.0 - D;
  as[1][gy] = -Fs/2.0 - D;
  ap[1][gy] = -( aw[1][gy] + an[1][gy] + 2.0*as[1][gy]) + (Fe - Fw) + (Fn - Fs);
  su[1][gy] = (p[1][gy-1]-p[1][gy]);
  ae[1][gy] = 0.0;
  as[1][gy] = 0.0;


  // for gx,gy node
  convection_diffusive_flux_u(gx,gy);
  ae[gx][gy] = -max(0.0,-Fe) - D ;
  aw[gx][gy] = -max(0.0,Fw) - D;
  an[gx][gy] = Fn/2.0 - D;
  as[gx][gy] = -Fs/2.0 - D;
  ap[gx][gy] = -( aw[gx][gy] + as[gx][gy] + 2.0*an[gx][gy]) + (Fe - Fw) + (Fn - Fs);
  su[gx][gy] = (p[gx][gy-1]-p[gx][gy]);
  ae[gx][gy] = 0.0;
  an[gx][gy] = 0.0;

  
  for(i=I1-1;i<=I2+1;i++)
    {
      for(j=J1;j<=J2+1;j++)
	{
	  ae[i][j] = 0.0;
	  aw[i][j] = 0.0;
	  an[i][j] = 0.0;
	  as[i][j] = 0.0;
	  ap[i][j] = 1.0;
	  su[i][j] = 0.0;
	}
    }	   

  for(i=1;i<=gx;i++)
    {
      for(j=2;j<=gy;j++)
	{
           
	  Du[i][j] = (1.0/ap[i][j]);// d for correction equation
          

	}
    }
   
  return (0);
}

int coeff_v()
{
  int i=0,j=0;
  for(i=2;i<=gx;i++)
    for(j=1;j<=gy;j++)
      {
	convection_diffusive_flux_v(i,j);
	ae[i][j] = (Fe/2) - D ;
	aw[i][j] = (-Fw/2) - D;
	an[i][j] = (-max(0,-Fn)) - D;
	as[i][j] = -max(0,Fs) - D;
	ap[i][j] = -(ae[i][j] + aw[i][j] + an[i][j] + as[i][j]) + (Fe - Fw) + (Fn - Fs);
	su[i][j] = (p[i-1][j]-p[i][j]);

      } 
  // for right boundary
  for (i = 2;i<=gx;i++)
    {
      // right boundary
      convection_diffusive_flux_v(i,gy);
      ae[i][gy] = (Fe/2) - D ;
      aw[i][gy] = (-Fw/2) - D;
      an[i][gy] = (-max(0,-Fn)) - D;
      as[i][gy] = -max(0,Fs) - D;
      ap[i][gy] = -(aw[i][gy] + an[i][gy] + as[i][gy]) + (Fe - Fw) + (Fn - Fs);
      su[i][gy] = (p[i-1][gy]-p[i][gy]);
      ae[i][gy] = 0.0;


      //  left boundary
      convection_diffusive_flux_v(i,1);
      ae[i][1] = (Fe/2) - D ;
      aw[i][1] = (-Fw/2) - D;
      an[i][1] = (-max(0,-Fn)) - D;
      as[i][1] = -max(0,Fs) - D;
      ap[i][1] = -(ae[i][1] + (2*aw[i][1]) + an[i][1] + as[i][1]) + (Fe - Fw) + (Fn - Fs);
      su[i][1] = (p[i-1][1]-p[i][1]);
      aw[i][1] = 0.0;


    }

  for (j = 2;j<=gy;j++)
    {
      // top wall
      convection_diffusive_flux_v(gx,j);
      ae[gx][j] = (Fe/2) - D ;
      aw[gx][j] = (-Fw/2) - D;
      an[gx][j] = (-max(0,-Fn)) - D;
      as[gx][j] = -max(0,Fs) - D;
      ap[gx][j] = -(ae[gx][j] + aw[gx][j] + an[gx][j] + as[gx][j]) + (Fe - Fw) + (Fn - Fs);
      su[gx][j] = (p[gx-1][j]-p[gx][j]);
      an[gx][j] = 0.0;


      //  bottom boundary
      convection_diffusive_flux_v(2,j);
      ae[2][j] = (Fe/2) - D ;
      aw[2][j] = (-Fw/2) - D;
      an[2][j] = (-max(0,-Fn)) - D;
      as[2][j] = -max(0,Fs) - D;
      ap[2][j] = -(ae[2][j] + aw[2][j] + an[2][j] + as[2][j]) + (Fe - Fw) + (Fn - Fs);
      su[2][j] = (p[2-1][j]-p[2][j]);
      as[2][j] = 0.0;


    }

  // for 2,1  node
  convection_diffusive_flux_v(2,1);
  ae[2][1] = (Fe/2) - D ;
  aw[2][1] = (-Fw/2) - D;
  an[2][1] =  (-max(0,-Fn)) - D;
  as[2][1] =  -max(0,Fs) - D;
  ap[2][1] = -(ae[2][1]+(2*aw[2][1])+ an[2][1] + as[2][1]) + (Fe - Fw) + (Fn - Fs);
  su[2][1] = (p[1][1]-p[2][1]) ;
  as[2][1] = 0.0;
  aw[2][1] = 0.0;

  // for gx,1 node
  convection_diffusive_flux_v(gx,1);
  ae[gx][1] = (Fe/2) - D ;
  aw[gx][1] = (-Fw/2) - D;
  an[gx][1] =  (-max(0,-Fn)) - D;
  as[gx][1] =  -max(0,Fs) - D;
  ap[gx][1] = -( 2.0*an[gx][1]+(2*aw[gx][1])+ ae[gx][1] + as[gx][1]) + (Fe - Fw) + (Fn - Fs);
  su[gx][1] = (p[gx-1][1]-p[gx][1]);
  an[gx][1] = 0.0;
  aw[gx][1] = 0.0;

  // for 2,gy node
  convection_diffusive_flux_v(2,gy);
  ae[2][gy] = (Fe/2) - D ;
  aw[2][gy] = (-Fw/2) - D;
  an[2][gy] = (-max(0,-Fn)) - D;
  as[2][gy] = -max(0,Fs) - D;
  ap[2][gy] = -( aw[2][gy] + an[2][gy] + as[2][gy]) + (Fe - Fw) + (Fn - Fs);
  su[2][gy] = (p[1][gy]-p[2][gy]);
  ae[2][gy] = 0.0;
  as[2][gy] = 0.0;


  // for gx,gy node
  convection_diffusive_flux_v(gx,gy);
  ae[gx][gy] =(Fe/2) - D ;
  aw[gx][gy] = (-Fw/2) - D;
  an[gx][gy] =(-max(0,-Fn)) - D;
  as[gx][gy] = -max(0,Fs) - D;
  ap[gx][gy] = -( aw[gx][gy] + as[gx][gy] + an[gx][gy]) + (Fe - Fw) + (Fn - Fs);
  su[gx][gy] = (p[gx-1][gy]-p[gx][gy]);
  ae[gx][gy] = 0.0;
  an[gx][gy] = 0.0;

  for(i=I1;i<=I2+1;i++)
    {
      for(j=J1-1;j<=J2+1;j++)
	{
	  ae[i][j] = 0.0;
	  aw[i][j] = 0.0;
	  an[i][j] = 0.0;
	  as[i][j] = 0.0;
	  ap[i][j] = 1.0;
	  su[i][j] = 0.0;
	}
    }
 
  for(i=2;i<=gx;i++)
    {
      for(j=1;j<=gy;j++)
	{
           
	  Dv[i][j] = (1.0/ap[i][j]);// d for correction equation
          

	}
    }
   
  return (0);
}





int velocity_pressure_initialisation()
{ 
  int i,j;

  for(i=0;i<=gx+1;i++)
    {   for(j=0;j<=gy+1;j++)     
	{
      
	  u[i][j]=0.1;//guess values
	  v[i][j]=0.0 ;
	  p[i][j]=0.1;// linear guess

	}
    }
 
  return(0);
}


int Boundary_condition_u()
{ int i;
  int j;
  for(i=0;i<=gx;i++)//inlet velocity
    {   
                 
      u[i][1]=Uo;//inlet velocty
      //u[i][gy+1]=u[i][gy];//outlet grad=0
    }
  for(i=1;i<=gy;i++)//wall
    {   
                  
                  
      u[0][i]=-u[1][i];//no slip lower wall
      u[gx+1][i]=-u[gx][i];//no slip upper wall
    }

  for(i=I1-1;i<=I2+1;i++)
    {
      for(j=J1;j<=J2+1;j++)
	{
	  u[i][j] = 0.0;
	}
    }
  return(0);
}


int Boundary_condition_v()
{ int i; int j;
  for(i=1;i<=gx;i++)//inlet velocity
    {    
      v[i][0]=0;
                   
      // v[i][0]=-v[i][1];//inlet velocity v=0
      v[i][gy+1]=v[i][gy];//outlet grad=0
    }
  for(i=1;i<=gy;i++)//inlet velocity
    {   
                  
      v[0][i] = 0.0;
      v[1][i]=0.0;//no slip lower wall
      v[gx+1][i]=0.0;//no slip upper wall
    }

  for(i=I1;i<=I2+1;i++)
    {
      for(j=J1-1;j<=J2+1;j++)
	{
	  v[i][j] = 0.0;
	}
    }
  return(0);
}

    
                 
int Boundary_condition_p()
{ int i;int j;
  for(i=1;i<=gx;i++)
    {   
      p[i-1][0]=p[i][0];//grad at input
      p[i][gy+1]=0;
      // p[i][0]=(2.0*Po)-p[i][1];// pressure grad.
      //p[i][gy+1]=-p[i][gy];//pressure=0  outlet
    }
  for(i=1;i<=gy;i++)
    {   
                  
                  
      p[0][i]=p[1][i];// wall grad=0
      p[gx+1][i]=p[gx][i];// wall grad=0
    }

  for(i=I1+1;i<I2;i++)
    {
      for(j=J1+1;j<J2;j++)
	{
	  p[i][j] =0.0;
	}
    }
  return(0);
}





int Boundary_condition_all()
{    
            
  Boundary_condition_u();
  Boundary_condition_v();
           
  Boundary_condition_p();
  return(0);
}

int continuity_eq()
{

  int i,j;

  for(i=1;i<=gx;i++)
    {  for(j=1;j<=gy;j++)
	{  
	  cE[i][j]=-Du[i][j+1];
	  cW[i][j]=-Du[i][j];
	  cN[i][j]=-Dv[i+1][j];
	  cS[i][j]=-Dv[i][j];
	  cP[i][j]=-(cE[i][j]+cW[i][j]+cN[i][j]+cS[i][j]);
      
	  cSU[i][j]=(u[i][j]-u[i][j+1])+(v[i][j]-v[i+1][j]);
        }
    }
         
  //left boundary (i,1) and right boundary(i,gy)
     
  for(i=1;i<=gx;i++)
    {   cW[i][1]=0.0;
            
      cE[i][gy]=0.0;

    }
           
  for(j=1;j<=gy;j++)// upper and lower wall
    {   
      cP[1][j]+=cS[1][j];
      cS[1][j] = 0.0; 

      cP[gx][j] +=cN[gx][j];        
      cN[gx][j]=0.0;

    }
   for(i=I1;i<=I2;i++)
    {
      cE[i][J1-1]= 0.0;
      cW[i][J1-1]=-Du[i][J1-1];
      cN[i][J1-1]=-Dv[i+1][J1-1];
      cS[i][J1-1]=-Dv[i][J1-1];
      cP[i][J1-1]= -(cE[i][J1-1]+cW[i][J1-1]+cN[i][J1-1]+cS[i][J1-1]);
      cSU[i][J1-1]= (u[i][J1-1]-u[i][J1])+(v[i][J1-1]-v[i+1][J1-1]);

      
      cE[i][J2+1]= -Du[i][J2+2];
      cW[i][J2+1]= 0.0;
      cN[i][J2+1]=-Dv[i+1][J2+1];
      cS[i][J2+1]=-Dv[i][J2+1];
      cP[i][J2+1]= -(cE[i][J2+1]+cW[i][J2+1]+cN[i][J2+1]+cS[i][J2+1]);
      cSU[i][J2+1]= (u[i][J2+1]-u[i][J2+2])+(v[i][J2+1]-v[i+1][J2+1]);
      
      
    }
    for(j=J1;j<=J2;j++)
    {
      cE[I1-1][j]=-Du[I1-1][j+1];
      cW[I1-1][j]=-Du[I1-1][j];
      cN[I1-1][j]= 0.0;
      cS[I1-1][j]=-Dv[I1-1][j];
      cP[I1-1][j]=-(cE[I1-1][j]+cW[I1-1][j]+cN[I1-1][j]+cS[I1-1][j]);
      cSU[I1-1][j]=(u[I1-1][j]-u[I1-1][j+1])+(v[I1-1][j]-v[I1][j]);

      cE[I2+1][j]=-Du[I2+1][j+1];
      cW[I2+1][j]=-Du[I2+1][j];
      cN[I2+1][j]=-Dv[I2+2][j];
      cS[I2+1][j]= 0.0;
      cP[I2+1][j]=-(cE[I2+1][j]+cW[I2+1][j]+cN[I2+1][j]+cS[I2+1][j]);
      cSU[I2+1][j]=(u[I2+1][j]-u[I2+1][j+1])+(v[I2+1][j]-v[I2+2][j]);
      

      
      
    }

    for(i=I1;i<=I2;i++)
      {
	for(j=J1;j<=J2;j++)
	  {
	    cE[i][j]= 0.0;
	    cW[i][j]= 0.0 ;
	    cN[i][j]= 0.0;
	    cS[i][j]= 0.0;
	    cP[i][j]= 1.0;
      
	    cSU[i][j]= 0.0;
	  }
      }
       

      

      
  return(0);
}

int simple()
{   int i,j;
  do
    { 
      cout<<"Entering SIMPLE";
      coeff_U();
      TDMA2D (1);
      // gauss_siedel(1);
      Boundary_condition_all();
      coeff_v();
      TDMA2D (2);
      //gauss_siedel(2);
      Boundary_condition_all();
      continuity_eq();
      TDMA2D (3);
      //gauss_siedel(3);

      for(i=I1;i<=I2;i++)
	{
	  p_err[i][J1] = p_err[i][J1-1];
	  p_err[i][J2] = p_err[i][J2+1];
	}
      for(j=J1;j<=J2;j++)
	{
	  p_err[I1][j] = p_err[I1-1][j];
	  p_err[I2][j] = p_err[I2+1][j];

	}

      for(i=2;i<=gx;i++)// updating v velocity
	{
	  for(j=1;j<=gy;j++)
	    {
           
	      v[i][j]=v[i][j]+(0.3*(Dv[i][j]*(p_err[i-1][j]-p_err[i][j])));
          

	    }
       
	}
      for(i=1;i<=gx;i++)//updating u velocity
	{
	  for(j=2;j<=gy;j++)
	    {
           
	      u[i][j]=u[i][j]+(0.3*(Du[i][j]*(p_err[i][j-1]-p_err[i][j])));
          

	    }
	}
     
      max_err = p_err[1][1];
      cout<<max_err;
      for(i=1;i<=gx;i++)
	{
	  for(j=1;j<=gy;j++)
	    {  
	      if(p_err[i][j]>max_err)
		{
		  cout<<"P_err[i][j]" <<p_err[i][j]<<"\n";
		  max_err=p_err[i][j];
		}
            
	      p[i][j]=p[i][j]+ (0.7*(p_err[i][j]));//updating pressure
          

	    }
	}
      cout<<" MAX ERROR : "<<max_err<<"\n";
      Boundary_condition_all();

    }while(max_err>.07);


    

  return(0);
}
/*
void write_all()
{
  int t=0;
  int i=0,j=0;
  FILE  *fpw;
  char filename[20];
  char no [20];
  char time[10];
  sprintf (filename, "velocity%d.csv", t);
  fpw =fopen(filename,"w");

  //for U
  fprintf(fpw, "nx = %d",gx+2);
    fprintf(fpw, "ny = %d",gy+2);
    fprintf(fpw, "\n");
  fprintf(fpw, " x coord , y coord , z coord , U \n");
  for(i=0;i<=gx+1;i++)
    {
      for(j=0;j<=gy+1;j++)
	{
	  fprintf(fpw,"%d,%d,%d,%lf\n",i,j,0,u[i][j]);
	}
    }


  fprintf(fpw,"\n");
  //for V
  fprintf(fpw, " x coord , y coord , z coord , V \n");
  for(i=0;i<=gx+1;i++)
  {
  for(j=0;j<=gy+1;j++)
  {
  fprintf(fpw,"%d,%d,%d,%lf, ",i,j,0,v[i][j]);
  }
  }

  //For P
  fprintf(fpw, " x coord , y coord , z coord , P \n");

  for(i=0;i<=gx+1;i++)
  {
  for(j=0;j<=gy+1;j++)
  {
  fprintf(fpw,"%d,%d,%d,%lf\n ",i,j,0,p[i][j]);
  }
  }
  fprintf(fpw,"\n");
  fclose(fpw);
}
*/


void write_all()
{
  int t=0;
  int i=0,j=0;
  FILE  *fpw;
  char filename[20];
  char no [20];
  char time[10];
  sprintf (filename, "U.txt");
  fpw =fopen(filename,"w");

  //for U
  for(i=0;i<=gx+1;i++)
    {
      for(j=0;j<=gy+1;j++)
	{
	  fprintf(fpw,"%.2lf    ",u[i][j]);
	}
      fprintf(fpw,"\n");
    }


  fprintf(fpw,"\n");fprintf(fpw,"\n");fprintf(fpw,"\n");
  fclose(fpw);
  //for V
  sprintf (filename, "V.txt");
  fpw =fopen(filename,"w");
  for(i=0;i<=gx+1;i++)
    {
      for(j=0;j<=gy+1;j++)
	{
	  fprintf(fpw,"%.2lf    ",v[i][j]);
	}
      fprintf(fpw,"\n");
    }
   fclose(fpw);

  //For P

  sprintf (filename, "P.txt");
  fpw =fopen(filename,"w");

  for(i=0;i<=gx+1;i++)
    {
      for(j=0;j<=gy+1;j++)
	{
	  fprintf(fpw,"%.2lf    ",p[i][j]);
	}
      fprintf(fpw,"\n");
    }
  fprintf(fpw,"\n\n\n");
  fclose(fpw);
  }

void write_to_vtk()
{
    int i = 0,j = 0 ;
    double X = 0.0, Y = 0.0;
    double velocity_u = 0.0, velocity_v = 0.0;
    ofstream fp;
    fp.open("U.vtk");
    if(fp.is_open())
    {
        fp<<"# vtk DataFile Version 2.0 \n";
        fp<<"Temp profile\n";
        fp<<"ASCII\n";
        fp<<"DATASET RECTILINEAR_GRID\n";
        fp<<"DIMENSIONS "<<gx+1<<" "<<gy+1<<" "<<"1"<<"\n";
        fp<<"X_COORDINATES"<<" "<<gx+1<<" "<<"float\n";
        for(i=0;i<=gx;i++)
        {
            fp<<X<<" ";
            X = X + delx;
        }
        fp<<"\n";
        fp<<"Y_COORDINATES"<<" "<<gy+1<<" "<<"float\n";
        for(j=0;j<=gy;j++)
        {
            fp<<Y<<" ";
            Y = Y + delx;
        }
        fp<<"\n";
        fp<<"Z_COORDINATES"<<" "<<"1"<<" "<<"float\n";
        fp<<"0\n";
        fp<<"\n";

        fp<<"\n";
        fp<<"\n";
        fp<<"POINT_DATA "<<(gx+1)*(gy+1)<<"\n";
        fp<<"SCALARS U double  1\n";
        fp<<"LOOKUP_TABLE default\n";
        for(j=0;j<=gy;j++)
        {
            for(i=0;i<=gx;i++)
            {
                velocity_u = (u[i][j] + u[i][j+1] )/2.0;
                fp<<velocity_u<<" ";
            }
        }
        fp<<"\n";
        fp<<"\n";
        fp<<"POINT_DATA "<<(gx+1)*(gy+1)<<"\n";
        fp<<"SCALARS V double  1\n";
        fp<<"LOOKUP_TABLE default\n";
        for(j=0;j<=gy;j++)
        {
            for(i=0;i<=gx;i++)
            {
                velocity_v = (v[i][j] + v[i+1][j])/2.0;
                fp<<velocity_v<<" ";
            }
        }
        fp<<"\n";
        fp<<"\n";
        fp<<"POINT_DATA "<<(gx+1)*(gy+1)<<"\n";
        fp<<"SCALARS Pressure double  1\n";
        fp<<"LOOKUP_TABLE default\n";
        for(j=0;j<=gy;j++)
        {
            for(i=0;i<=gx;i++)
            {
                fp<<p[i][j]<<" ";
            }
        }

        
    }
    else
    {
        printf("Cannot open file\n");
    }
    fp.close();
    
    
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////       
int main()
{   int i,j;
  velocity_pressure_initialisation();
  //cout<<"111111";
  Boundary_condition_all();

    
  cout<<"222222";
  simple();
  cout<<"33333";
  write_all();
  cout<<"444444";
    write_to_vtk();

  return(0);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int TDMA2D (int option)
{
  double error_max =0.0;
  double difference = 0.0;
  double field_old[gx+2][gy+2];
  double field_new[gx+2][gy+2];
  double field1[gx+2];
  double aw1[gx+2];
  double ae1[gx+2];
  double ap1[gx+2];
  double an1[gx+2];
  double as1[gx+2];
  double su1[gx+2];
  double max1 = 0.0;
  int I,J;
  int g;
  /*
    Equation trying to be solved here is :
    a_w T_w + a_e T_e + a_p T_p + a_n T_n + a_s T_s = Su
  */




  if(option== 1)
    {
      int i=0,j=0;
      for(i=0;i<=gx+1;i++)
	for(j=0;j<=gy+1;j++)
	  {
	    field_old[i][j] = u[i][j];
	    field_new[i][j] = u[i][j];
	  }
      int REPEAT = 1;
      while(REPEAT)
    	{
	  printf(" Enetering TDMA for %d time ......\n",REPEAT);
	  REPEAT++;
   	
   
	  for (j=2;j<=gy;j++)
	    {
      	      for(i=1;i<=gx;i++)
		{
		  aw1[i]  = as[i][j];
		  ap1[i]  = ap[i][j];
		  ae1[i]  = an[i][j];
		  su1[i]  = su[i][j] - ae[i][j] * field_old[i][j+1] - aw[i][j] *field_new[i][j-1];
	     
	        }  

	      TDMA(field1,aw1,ae1,ap1,su1,1,gx);
	      for(i=1;i<=gx;i++)
		{
	          field_new[i][j] = field1[i];
		}
	    }
	  error_max = field_new[1][2] - field_old[1][2];

	  for(i=1;i<=gx;i++)
	    {
	      for(j=2;j<=gy;j++)
		{
		  difference = field_new[i][j] - field_old[i][j];
		  if(difference<0.0)
		    {
		      difference = -1.0*difference;
		    }
		  if(error_max<difference)
		    error_max = difference;
		}
	    }
               
              
      
	  for(i=1;i<=gx;i++)
	    {
	      for(j=2;j<=gy;j++)
		{
		  field_old[i][j] = field_new[i][j];
		}
	    }

	  if(error_max<0.25)
	    REPEAT = 0;
	  printf(" MAx error : %lf\n",error_max);

	}
      for(i=1;i<=gx;i++)
	for(j=2;j<=gy;j++)
	  u[i][j] = field_old[i][j] ;
    
    }


  if(option== 2)
    {
      int i=0,j=0;
      for(i=2;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  {
	    field_old[i][j] = v[i][j];
	    field_new[i][j] = v[i][j];
	  }
      int REPEAT = 1;
      while(REPEAT)
    	{
	  printf(" Enetering TDMA for %d time ......\n",REPEAT);
	  REPEAT++;
   	
   
	  for (j=1;j<=gy;j++)
	    {
      	      for(i=2;i<=gx;i++)
		{
		  aw1[i]  = as[i][j];
		  ap1[i]  = ap[i][j];
		  ae1[i]  = an[i][j];
		  su1[i]  = su[i][j] - ae[i][j] * field_old[i][j+1] - aw[i][j] *field_new[i][j-1];
	     
	        }  

	      TDMA(field1,aw1,ae1,ap1,su1,2,gx);
	      for(i=1;i<=gx;i++)
		{
	          field_new[i][j] = field1[i];
		}
	    }
	  error_max = field_new[2][1] - field_old[2][1];

	  for(i=2;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  difference = field_new[i][j] - field_old[i][j];
		  if(difference<0.0)
		    {
		      difference = -1.0*difference;
		    }
		  if(error_max<difference)
		    error_max = difference;
		}
	    }
               
              
      
	  for(i=2;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  field_old[i][j] = field_new[i][j];
		}
	    }

	  if(error_max<0.25)
	    REPEAT = 0;
	  printf(" MAx error : %lf\n",error_max);

	}
      for(i=2;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  v[i][j] = field_old[i][j] ;
    
    }


  if(option== 3)
    {
      int i=0,j=0;
      for(i=1;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  {
	    field_old[i][j] = p_err[i][j];
	    field_new[i][j] = p_err[i][j];
	  }
      int REPEAT = 1;
      while(REPEAT)
    	{
	  printf(" Enetering TDMA for %d time ......\n",REPEAT);
	  REPEAT++;
   	
   
	  for (j=1;j<=gy;j++)
	    {
      	      for(i=1;i<=gx;i++)
		{
		  aw1[i]  = cS[i][j];
		  ap1[i]  = cP[i][j];
		  ae1[i]  = cN[i][j];
		  su1[i]  = cSU[i][j] - cE[i][j] * field_old[i][j+1] - cW[i][j] *field_new[i][j-1];
	     
	        }  

	      TDMA(field1,aw1,ae1,ap1,su1,1,gx);
	      for(i=1;i<=gx;i++)
		{
	          field_new[i][j] = field1[i];
		}
	    }
	  error_max = field_new[1][1] - field_old[1][1];

	  for(i=1;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  difference = field_new[i][j] - field_old[i][j];
		  if(difference<0.0)
		    {
		      difference = -1.0*difference;
		    }
		  if(error_max<difference)
		    error_max = difference;
		}
	    }
               
              
      
	  for(i=1;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  field_old[i][j] = field_new[i][j];
		}
	    }

	  if(error_max<0.25)
	    REPEAT = 0;
	  printf(" MAx error : %lf\n",error_max);

	}
      for(i=1;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  p_err[i][j] = field_old[i][j] ;
    
    }

  return(0);
}


void TDMA(double *field,double *aw,double *ae,double *ap,double *su,int start,int end)
{
  int i=0;
  int j=0;
  double denom;
  su[start]/= ap[start];
  ae[start]/= ap[start];
  ap[start]/= ap[start];
  aw[start]/= ap[start];
  for( i=start+1;i<=end;i++)
    {
      denom =(ap[i] -(aw[i]*ae[i-1]));
      if((denom<0.0001)&&(denom>(-0.0001)))
	printf(".......OOOOPPPPPPPS  DENOM BECAME ZERO HERE VERY PATHETIC ............\n");
      su[i] = (su[i] - (su[i-1]*aw[i]))/ denom;
      ae[i] = ae[i]/denom;
      ap[i] = 1.0;
      aw[i] = 0;;
    }

  // Solving for field by back substitution

  for(j=end;j>=start;j--)
    {
      field[j] = su[j] - ae[j]*field[j+1] ;
    }
}


void gauss_siedel(int option)
{
  double field_old[gx+2][gy+2];
  double field_new[gx+2][gy+2];
  double difference =0.0;
  double error_max = 0.0;

  if( option ==1)
    {
      int i=0,j=0;
      for(i=1;i<=gx;i++)
	for(j=2;j<=gy;j++)
	  {
	    field_old[i][j] = u[i][j];
	    field_new[i][j] = u[i][j];
	  }
      int REPEAT = 1;
      while(REPEAT)
    	{
	  printf(" Enetering GS1 for %d time ......\n",REPEAT);
	  REPEAT++;
   	
   
	  for (i=1;i<=gx;i++)
	    {
      	      for(j=2;j<=gy;j++)
		{
		  field_new[i][j] = su[i][j] - ae[i][j]*field_old[i][j+1] - aw[i][j]*field_old[i][j-1] - an[i][j]*field_old[i+1][j] - as[i][j]*field_old[i-1][j];
		  field_new[i][j] = field_new[i][j]/ap[i][j];
	        }  

	    }
	  error_max = field_new[1][2] - field_old[1][2];

	  for(i=1;i<=gx;i++)
	    {
	      for(j=2;j<=gy;j++)
		{
		  difference = field_new[i][j] - field_old[i][j];
		  if(difference<0.0)
		    {
		      difference = -1.0*difference;
		    }
		  if(error_max<difference)
		    error_max = difference;
		}
	    }
               
              
      
	  for(i=1;i<=gx;i++)
	    {
	      for(j=2;j<=gy;j++)
		{
		  field_old[i][j] = field_new[i][j];
		}
	    }

	  if(error_max<0.25)
	    REPEAT = 0;
	  printf(" MAx error : %lf\n",error_max);

	}
      for(i=1;i<=gx;i++)
	for(j=2;j<=gy;j++)
	  u[i][j] = field_old[i][j] ;
    
    }

  if(option== 2)
    {
      int i=0,j=0;
      for(i=2;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  {
	    field_old[i][j] = v[i][j];
	    field_new[i][j] = v[i][j];
	  }
      int REPEAT = 1;
      while(REPEAT)
    	{
	  printf(" Enetering GS2 for %d time ......\n",REPEAT);
	  REPEAT++;
   	
   
	  for (j=1;j<=gy;j++)
	    {
      	      for(i=2;i<=gx;i++)
		{
		  field_new[i][j] = su[i][j] - ae[i][j]*field_old[i][j+1] - aw[i][j]*field_old[i][j-1] - an[i][j]*field_old[i+1][j] - as[i][j]*field_old[i-1][j];
		  field_new[i][j] = field_new[i][j]/ap[i][j];
	     
	        }
	    }
	  error_max = field_new[2][1] - field_old[2][1];

	  for(i=2;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  difference = field_new[i][j] - field_old[i][j];
		  if(difference<0.0)
		    {
		      difference = -1.0*difference;
		    }
		  if(error_max<difference)
		    error_max = difference;
		}
	    }
               
              
      
	  for(i=2;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  field_old[i][j] = field_new[i][j];
		}
	    }

	  if(error_max<0.25)
	    REPEAT = 0;
	  printf(" MAx error : %lf\n",error_max);

	}
      for(i=2;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  v[i][j] = field_old[i][j] ;
    
    }

  if(option== 3)
    {
      int i=0,j=0;
      for(i=1;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  {
	    field_old[i][j] = p_err[i][j];
	    field_new[i][j] = p_err[i][j];
	  }
      int REPEAT = 1;
      while(REPEAT)
    	{
	  printf(" Enetering GS3 for %d time ......\n",REPEAT);
	  REPEAT++;
   	
   
	  for (j=1;j<=gy;j++)
	    {
      	      for(i=1;i<=gx;i++)
		{
		  field_new[i][j] = cSU[i][j] - cE[i][j]*field_old[i][j+1] - cW[i][j]*field_old[i][j-1] - cN[i][j]*field_old[i+1][j] - cS[i][j]*field_old[i-1][j];
		  field_new[i][j] = field_new[i][j]/cP[i][j];
                
	        }
	    }
	  error_max = field_new[1][1] - field_old[1][1];

	  for(i=1;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  difference = field_new[i][j] - field_old[i][j];
		  if(difference<0.0)
		    {
		      difference = -1.0*difference;
		    }
		  if(error_max<difference)
		    error_max = difference;
		}
	    }
               
              
      
	  for(i=1;i<=gx;i++)
	    {
	      for(j=1;j<=gy;j++)
		{
		  field_old[i][j] = field_new[i][j];
		}
	    }

	  if(error_max<0.25)
	    REPEAT = 0;
	  printf(" MAx error : %lf\n",error_max);

	}
      for(i=1;i<=gx;i++)
	for(j=1;j<=gy;j++)
	  p_err[i][j] = field_old[i][j] ;
    
    }

}








    






 
       
































