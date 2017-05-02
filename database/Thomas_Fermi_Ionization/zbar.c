#Code written by Jim Glosli, lightly edited by Jeff Haack

#include "zBar.h"
#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#define SQ(x)   ( (x) * (x))
#define CUBE(x)   ((x) * (x) * (x))
#define MAX(A,B) ((A) > (B) ? (A) : (B))

void fttfqSetV(ION_STATE *ion)
{
   double a1 = 0.003323;
   double a2 = 0.9718;
   double a3 = 9.26148e-5;
   double a4 = 3.10165;
   double b0 = -1.7630;
   double b1 = 1.43175;
   double b2 = 0.31546;
   double c1 = -0.366667;
   double c2 = 0.983333;

   double T = ion->T + 1e-100; 
   double T0 = T/pow(ion->Z,(4./3.)); 
   double Tf = T0/(1.+T0);
   ion->A = a1*pow(T0,a2)+a3*pow(T0,a4);
   ion->B = -exp(b0+b1*Tf+b2*pow(Tf,7.));
   ion->C = c1*Tf+c2;
}
double fttfqFASTV(ION_STATE *ion)
{
   double A = ion->A; 
   double B = ion->B; 
   double C = ion->C; 
   double Z = ion->Z; 
   double vol = ion->vol; 

   double alpha = 14.3139;
   double beta = 0.6624;
   double R = (Z*vol);
   double Q1 = A*pow(R,-B);
   double Y1 = pow(R,-C) ; 
   double Y2 = pow(Q1,C); 
   double Y = Y1 + Y2; 
   double Q = pow(Y,(1./C));
   //double Q = pow(pow(R,-C)+pow(Q1,C),(1./C));
   double x = alpha*pow(Q,beta);
   double zbar = Z*x/(1.+x+sqrt(1.+2.*x));

   double dR = Z; 
   double dQ1 = -B*Q1/R*dR; 
   double dY = C*(-Y1/R*dR + Y2/Q1*dQ1); 
   double dQ = Q/(C*Y)*dY; 
   double dx = x * beta/Q*dQ; 

   double t1 = sqrt(1.+2.*x);
   double dt1 = 1/t1 * dx; 
   double t2 = 1 + x + t1; 
   double dt2  = (dx + dt1); 
   double zBar =   Z*x/t2; 
   double dzBar = zBar/x*dx - zBar/t2*dt2; 
   double rho = zbar/vol; 
   double drho = (dzBar - rho)/vol; 
   ion->zBar = zBar; 
   ion->nf = rho; 
   ion->dnf = drho; 

   return rho;
}
double fttfq(double Z, double rho, double T)
{
   double alpha = 14.3139;
   double beta = 0.6624;
   double R = rho/Z;
   double a1 = 0.003323;
   double a2 = 0.9718;
   double a3 = 9.26148e-5;
   double a4 = 3.10165;
   double b0 = -1.7630;
   double b1 = 1.43175;
   double b2 = 0.31546;
   double c1 = -0.366667;
   double c2 = 0.983333;

   double T0 = T/pow(Z,(4./3.)); 
   double Tf = T0/(1.+T0);
   double A = a1*pow(T0,a2)+a3*pow(T0,a4);
   double B = -exp(b0+b1*Tf+b2*pow(Tf,7.));
   double C = c1*Tf+c2;
   double Q1 = A*pow(R,B);
   double Q = pow(pow(R,C)+pow(Q1,C),(1./C));
   double x = alpha*pow(Q,beta);
   double zbar = Z*x/(1.+x+sqrt(1.+2.*x));
   printf ("Z=%e r=%e x=%e zbar=%e\n",Z,rho,x,zbar/Z); 
   return zbar;
}


// Zero temperature version

double fttfq0(double Z,double rho)
{
   double alpha = 14.3139;
   double beta = 0.6624;

   double x = alpha*pow((rho/(Z)),beta);
   double f = x/(1. + x + sqrt(1.+2.*x));
   double zbar = f*Z;
   return zbar; 
}
void zBarFunc(int nIonType,ION_STATE *ion)
{
   double ntotal =0; 
   double n_max = 0.0; 
   int i;
   for (i=0;i<nIonType;i++)
   {
     fttfqSetV(ion+i);
     ntotal += ion[i].n; 
     n_max = MAX(n_max,ion[i].n); 
   }
   if (ntotal == 0.0)
   {
      ion->zBar = ion->Z; 
      ion->nf = 0.0; 
      ion->dnf = 0.0;  
      return; 
   }
   double vol = 1/ntotal; 
   double nf=0.0; 
   int flag = 0; 
   int cnt =0; 
   for (i=0;i<nIonType;i++)
   {
      ion[i].vol = vol; 
      fttfqFASTV(ion+i); // free electron number density ; 
      nf += ion[i].nf*ion[i].n; 
      if ( ion[i].n/n_max > 1e-8) {flag += 1; cnt++;}
      flag *= 2; 
   }
   nf *= vol; 
   int loop =0; 
   double error=0; 
   double delta, delta_nf, f, df;
   for (loop =0;loop<50;loop++) 
   {
      f  =0; 
      df =0; 
      error=0.0; 
      for (i=0;i<nIonType;i++) 
      {
         fttfqFASTV(ion+i);
         delta=-(ion[i].nf - nf)/ion[i].dnf; 
         if (ion[i].vol+delta> 1.0/ion[i].n) delta = 0.5*(1.0/ion[i].n-ion[i].vol);
         if (ion[i].vol+delta< 0.0)          delta = -0.5*ion[i].vol;
         ion[i].vol  += delta;     
         f  += ion[i].n*ion[i].vol; 
         df += ion[i].n/ion[i].dnf; 
         error = SQ(ion[i].nf/nf-1.0);
      }

      delta_nf = -(f-1.0)/df; 
      if (delta_nf> 0.2*nf) delta_nf=0.2*nf; 
      if (delta_nf< -0.2*nf) delta_nf=-0.2*nf; 
      nf   += delta_nf;

      error += SQ(f-1); 
      error = sqrt(error); 

      if (error < 1e-12) break;
   }

}


#define Na  6.0221415e23

void zBarFunc2(int nspec, double T, double *Z, double *n, double *zBar)
{
  int i;
  ION_STATE *ion = (ION_STATE*) malloc(nspec*sizeof(ION_STATE)); 
  
  for (i=0;i<nspec;i++)    {
    ion[i].T = T;   //eV
    ion[i].n = n[i]/Na; //moles/cc
    ion[i].Z = Z[i];
  }

  zBarFunc(nspec,ion);

  for(i=0;i<nspec;i++) {
    zBar[i] = ion[i].zBar;
  }

}
