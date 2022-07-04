#include<stdio.h>
#include<math.h>
#include<stdlib.h>
/*--------------------------------------------------------*/
#define sqr(x) ((x)*(x))
#define fabs(x) ((x)>0.0 ? (x) : -(x))
#define max(a,b) ((a)>(b)? (a):(b))
#define sgn(x) ((x)>=0.0 ? 1.0 : -1.0)
#define pi M_PI

int n=3;		/* number of equations */
double tt,x[3];
double work[9];

double dt;	/* time-step */

/*==========================*/
double B0,Bi;
double K;
double Delta;
double P;
double D;
double Tau=7.0;
/*==========================*/
double TD[100];
double LD[100];
int N;


int F();
FILE *out;

main()
{
	int i;
double JJ;
double Jmin;
int k;
	out=fopen("output.dat","w");
	dt=1e-5;
	dt=1e-4;
/*-------------------------------------*/
TD[0]=13.0;
LD[0]=4.037824750588342;
TD[1]=16.0;
LD[1]=5.03382569395331;
TD[2]=20.0;
LD[2]=5.903089986991944;
TD[3]=23.0;
LD[3]=5.964858081947589;
TD[4]=27.0;
LD[4]=5.5319895514125506;
TD[5]=30.0;
LD[5]=4.8909795969896885;
TD[6]=36.0;
LD[6]=4.439963935920905;
TD[7]=38.0;
LD[7]=4.1673173347481764;
TD[8]=43.0;
LD[8]=4.263399331334003;
N=9;
srand48(1);
/*-------------------------------------*/
Jmin=1e9;

Jmin=60.0;

for(k=0;k<50000;++k)
{
Bi=drand48()*1e-5;
B0=Bi+drand48()*1e-5;
K=20*drand48();
Delta=2*drand48();
P=drand48()*1e5;
D=2*drand48();
/*-------------------------------------*/
	i=0;
	tt=0.0;
	x[0]=1.0e4;
	x[1]=0.0;
	x[2]=1.0e-3;
	JJ=0.0;

	while(tt<50.0 && i<N)
	{
		rk(x,tt,dt,3,F,work);
		tt+=dt;
		if(tt>TD[i])
		{
			JJ+=sqr(log(x[2])-LD[i]);
			i++;
		}
	}
	JJ/=(1.0*N);
if(JJ<Jmin)
{
	//Jmin=JJ;
/* choose to print only the lowest or all options below a J_threshold */
printf("%g \t %g %g %g %g %g %g\n",JJ,B0,Bi,K,Delta,P,D);
fprintf(out,"%g \t %g %g %g %g %g %g\n",JJ,B0,Bi,K,Delta,P,D);
fflush(out);
}
}
}
/*--------------------------------------------------------------------------*/
rk(x,tt,h,n,F,work)/* fourth order general Runge-Kutta for a system of ODEs */
double *x;	/* solves dx/dt=F(x,t) */
int n;	/* n equations/unknowns */
double *work; /* workspace array of 3*n doubles */
double tt,h;
int (*F)(); /* subroutine calculating RHS's of ODEs */
{
	int i;
	double t;
	double h2=h/2.0;
	double h6=h/6.0;

	double *r1,*r0;
	double *r2,*r3,*r4;
	double *x1;
	r0=work;
	r1=r0+n;
	x1=r0+2*n;
	r4=r3=r2=r1;

	t=tt;
	F(x,r1,t);

	for(i=0;i<n;++i)
	{
		x1[i]=x[i]+h2*r1[i];
		r0[i]=r1[i];
	}
	t=tt+h2;
	F(x1,r2,t);

	for(i=0;i<n;++i)
	{
		x1[i]=x[i]+h2*r2[i];
		r0[i]+=2.0*r2[i];
	}
	t=tt+h2;
	F(x1,r3,t);

	for(i=0;i<n;++i)
	{
		x1[i]=x[i]+h*r3[i];
		r0[i]+=2.0*r3[i];
	}
	t=tt+h;
	F(x1,r4,t);

	for(i=0;i<n;++i)
		x[i]+=h6*(r0[i]+r4[i]);
  /* dont update time */
}
/*--------------------------------------------------------------------------*/
F(X,R,t)	/* calculates ODE rhs's, returns them in R[] */
double *X,*R,t;
{
	double T,I,V;
	double Tp,Ip,Vp;

	double beta;

if(t<Tau)
beta=B0;
else
beta=Bi+(B0-Bi)*exp(-K*(t-Tau));

	T=X[0];
	I=X[1];
	V=X[2];

	Tp=D*(1.0e4-T)-beta*T*V;
	Ip=beta*T*V-Delta*I;
	Vp=P*I-23*V;

	R[0]=Tp;
	R[1]=Ip;
	R[2]=Vp;
}

/*----------------------------------------------------*/

