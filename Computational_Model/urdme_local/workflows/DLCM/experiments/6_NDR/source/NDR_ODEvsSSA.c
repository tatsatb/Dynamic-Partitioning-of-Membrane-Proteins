/* [Remove/modify this line not to overwrite this file] */
/* Generated by RPARSE 2018-08-30 12:20 */

/* Reactions:
     @ > betaN*vol > N
     N > N*(wa*D_a+wb*D_b)/kt > @
     N+D > N*D/(kc*vol) > @
     N > N > @
     @ > betaD*vol/(1.0+pow(R/(VolR*vol),2.0)) > D
     D > D*(wa*N_a+wb*N_b)/kt > @
     D > D > @
     @ > betaR*vol/(kRS/pow((qa*D_a+qb*D_b)/VolND*(N/(vol*VolND)),2.0)+1.0) > R
     R > R > @
*/

#include <math.h>
#include "propensities.h"
#include "report.h"

enum Species {
  N,
  D,
  R,
  N_a,
  N_b,
  D_a,
  D_b
};

const int NR = 9; /* number of reactions */

/* rate constants */
const double VolND = 4.000000000000000000000e+02;
const double VolR = 4.000000000000000000000e+02;
const double betaN = 4.000000000000000000000e+04;
const double betaD = 2.000000000000000000000e+05;
const double betaR = 1.200000000000000000000e+08;
const double kt = 8.000000000000000000000e+02;
const double kc = 2.000000000000000000000e+02;
const double kRS = 1.000000000000000000000e+07;
const double wa = 1.000000000000000000000e+00;
const double wb = 1.000000000000000000000e+00;
const double qa = 1.000000000000000000000e+00;
const double qb = 1.000000000000000000000e+00;

/* forward declaration */
double rFun1(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun2(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun3(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun4(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun5(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun6(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun7(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun8(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);
double rFun9(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd);

/* static propensity vector */
static PropensityFun ptr[] = {rFun1,rFun2,rFun3,rFun4,rFun5,rFun6,rFun7,rFun8,rFun9};

/* propensity definitions */
double rFun1(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return betaN*vol;
}

double rFun2(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return xstate[N]*(wa*xstate[D_a]+wb*xstate[D_b])/kt;
}

double rFun3(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return xstate[N]*xstate[D]/(kc*vol);
}

double rFun4(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return xstate[N];
}

double rFun5(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return betaD*vol/(1.0+pow(xstate[R]/(VolR*vol),2.0));
}

double rFun6(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return xstate[D]*(wa*xstate[N_a]+wb*xstate[N_b])/kt;
}

double rFun7(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return xstate[D];
}

double rFun8(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return betaR*vol/(kRS/pow((qa*xstate[D_a]+qb*xstate[D_b])/VolND*(xstate[N]/(vol*VolND)),2.0)+1.0);
}

double rFun9(const int *xstate,double time,double vol,
             const double *ldata,const double *gdata,int sd)
{
  return xstate[R];
}

/* URDME solver interface */
PropensityFun *ALLOC_propensities(size_t Mreactions)
{
  if (Mreactions > NR) PERROR("Wrong number of reactions.");
  return ptr;
}

void FREE_propensities(PropensityFun *ptr)
{ /* do nothing since a static array was used */ }

