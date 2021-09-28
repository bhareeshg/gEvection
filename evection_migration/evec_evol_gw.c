#include<stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
#include<math.h>
#include <unistd.h>
#include<string.h>
#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286 
#define G 4*PI*PI

double a2;
double m0;
double m1;
double m2;
double J2;
double R0;
double n2;
double c=63197.639580028;
double rate_a1,rate_e1,rate_a2;
double havg(double a1,double e1,double sigma)
{
	double csigma=sqrt(G*(m0+m1)*a1)*(1-sqrt(1-(e1*e1)));
	return n2*csigma - (G*J2*m0*pow(R0,2))/(2.*pow(a1,3)*pow(pow(-1 + csigma/sqrt(a1*G*(m0+m1)),2),1.5)) - (a1*m2*(2*a1*G*m0 + 6*sqrt(a1*G*m0)*csigma-3*pow(csigma,2) + 15*(2*sqrt(a1*G*m0)-csigma)*csigma*cos(2*sigma)))/(8.*pow(a2,3)*m0);
}

int func (double t, const double org[], double f[],
      void *params)
{
  	(void)(t); /* avoid unused parameter warning */
	
	// get orbital elements from the canonical variables
	double sigma=org[0];
	double csigma=org[1];
	double a1=org[2];
	double e1=sqrt(1-pow(1-(csigma/sqrt(G*(m0+m1)*a1)),2));

	// orbital evolution due to the companion as well as J2 term. See Touma and Wisdom 1998.
	f[0]=n2+(3*G*J2*m0*pow(R0,2)*(-1+csigma/sqrt(a1*G*(m0+m1))))/(2.*pow(a1,3)*sqrt(a1*G*(m0+m1))*pow(pow(-1 + csigma/sqrt(a1*G*(m0+m1)),2),2.5))-(a1*m2*(6*sqrt(a1*G*(m0+m1))-6*csigma + 15*(2*sqrt(a1*G*(m0+m1))-csigma)*cos(2*sigma)-15*csigma*cos(2*sigma)))/(8.*pow(a2,3)*(m0+m1));
	f[1]=-(15*a1*m2*(2*sqrt(a1*G*(m0+m1))-csigma)*csigma*sin(2*sigma))/(4.*pow(a2,3)*(m0+m1));

	double fgw,dedt,dadt,fprec;
	// add post newtonian precession. See Eggleton & Kiseleva–Eggleton 2001.
	fprec=(3*G*m0*m1*sqrt(a1*G*(m0 + m1))*(-1+csigma/sqrt(a1*G*(m0 + m1))))/(pow(a1,3)*pow(c,2)*pow(pow(-1+csigma/sqrt(a1*G*(m0 + m1)),2),1.5));
	f[0]+=fprec;

	// add GR radiation. See Peters 1964.
	fgw=-(304.0/15.0)*((pow(G,3)*m0*m1*(m0+m1))/(pow(c,5)*pow(a1,4)*pow(1-e1*e1,2.5)))*(1+(121.0*e1*e1/304.0));
	dedt=fgw*e1;
	if(e1>1e-2)
	{
		dedt+=(rate_e1*e1);
	}
	dadt=-(64.0/5.0)*(pow(G,3)*m0*m1*(m0+m1)/(pow(c,5)*pow(a1,3)*pow(1-e1*e1,3.5)))*(1+(73.0*e1*e1/24.0)+(37.0*pow(e1,4)/96.0));
	dadt+=(rate_a1*a1);
	f[1]+=sqrt(G*(m0+m1))*(((1-sqrt(1-(e1*e1)))*dadt/(2*sqrt(a1)))+(a1*e1*dedt/sqrt(a1*(1-(e1*e1)))));
	f[2]=dadt;
	
	return GSL_SUCCESS;
}

int integrateGSL(double a10,double e10,double sigma10,double tmax, FILE *fp)
{
	gsl_odeiv2_system sys = {func, NULL, 3, NULL};
	gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rkf45, 1e-8, 1e-10, 1e-6);
	
	//define variables which are evolved using GSL library:
	//resonance angle sigma, it's conjugate variable csigma and semi-major axis axis of the binary
	double org[3];
	org[0]=sigma10;
	org[1]=sqrt(G*(m0+m1)*a10)*(1-sqrt(1-e10*e10));
	org[2]=a10;

	// define initial timestep 
	double dt=1000;
	//define header
	fprintf(fp,"time,a1,sigma,e,a2,J2,a1rate,e1rate,a2rate\n");

	double t=0,ti;
	int i=0;
	double et,ch;
	int status;
	//define initial and final semi-major axis a2.
	double a2init=a2;
	double a2finl=pow(3*m2/m0,0.333)*a10;

	double sbeta,sdelta,salpha,eta,sigma,ean,a1;
	double porg[3];
	double pt,pet,pa;

	ti=t+dt;
	double fgw,dedt,dadt;
	while(t<tmax)
	{
		// record previous orbital variables
		porg[0]=org[0]; porg[1]=org[1]; porg[2]=org[2];
		pt=t;
		pa=org[2];
		pet=sqrt(1-pow(1-(org[1]/sqrt(G*(m0+m1)*a1)),2));
		
		// evolve the system
		status = gsl_odeiv2_driver_apply (d, &t, ti, org);
		if (status != GSL_SUCCESS)
		  break;
		a1=org[2];
		et=sqrt(1-pow(1-(org[1]/sqrt(G*(m0+m1)*a1)),2));

		// check if the timestep is too long
		if(fabs((pet-et)/et)>1e-1 || fabs((pa-a1)/a1)>1e-1)
		{
			dt=dt/2;
			t=pt;
			ti=t+dt;
			org[0]=porg[0]; org[1]=porg[1]; org[2]=porg[2]; 
			gsl_odeiv2_driver_reset(d);
			continue;
		}

		// output every 10th timestep
		if(i%10==0)
		{
			fprintf(fp,"%e,%e,%e,%e,%e,%e,%e,%e,%e\n",t,a1,org[0],et,a2,J2,rate_a1,rate_e1,rate_a2);
		}
		ti=t+dt;

		// if semi-major axis/eccentricity is too small, end the simulation.
		if(a1<5e-3 && et<1e-3)
		{
			break;
		}
		
		// update the a2: a2dot/a2=rate_a2 -> exponential evolution
		a2=(a2init*exp(rate_a2*t));
		if(a2<a2finl)
		{
			break;
		}
		n2=sqrt(G*(m2+m1+m0)/(a2*a2*a2));
		i++;
	}
	gsl_odeiv2_driver_free (d);
}

// function returns rate of change of resonance angle: used to calculate the location of resonance
double sigma_dot(double a2l, void *params)
{
	double e1=((double *)params)[0];
	double a1=((double *)params)[1];
	n2 = sqrt(G*(m2+m0+m1)/(a2l*a2l*a2l));
	double csigma=sqrt(G*(m0+m1)*a1)*(1-sqrt(1-e1*e1));
	double sigma=M_PI/2;
	double prec1=n2+(3*G*J2*m0*pow(R0,2)*(-1+csigma/sqrt(a1*G*(m0+m1))))/(2.*pow(a1,3)*sqrt(a1*G*(m0+m1))*pow(pow(-1 + csigma/sqrt(a1*G*(m0+m1)),2),2.5))-(a1*m2*(6*sqrt(a1*G*(m0+m1))-6*csigma + 15*(2*sqrt(a1*G*(m0+m1))-csigma)*cos(2*sigma)-15*csigma*cos(2*sigma)))/(8.*pow(a2l,3)*(m0+m1));
	double prec2=(3*G*m0*m1*sqrt(a1*G*(m0 + m1))*(-1+csigma/sqrt(a1*G*(m0 + m1))))/(pow(a1,3)*pow(c,2)*pow(pow(-1+csigma/sqrt(a1*G*(m0 + m1)),2),1.5));
	return prec1+prec2;


}

// calculate the a2 at which the evection resonance occurs i.e. solve the equation:
// sigma_dot(a2) = 0
double get_start_a2(double a1,double e1)
{
	int status;
	int iter = 0, max_iter = 4000;
	double r;
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double x_lo = 1e-5; 
	double x_hi = 1e10;
	gsl_function F;
	F.function = &sigma_dot;
	double param_arr[2];
	param_arr[0]=e1;
	param_arr[1]=a1;
	F.params = param_arr;
	
	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc (T);
	gsl_root_fsolver_set (s, &F, x_lo, x_hi);
	do
	{
		iter++;
		status = gsl_root_fsolver_iterate (s);
		r = gsl_root_fsolver_root (s);
		x_lo = gsl_root_fsolver_x_lower (s);
		x_hi = gsl_root_fsolver_x_upper (s);
		status = gsl_root_test_interval (x_lo, x_hi,0, 0.0001);
	}
	while (status == GSL_CONTINUE && iter < max_iter);
	gsl_root_fsolver_free (s);
	return r;
}

int main(int argn,char *argv[])
{

	double a1,tmax,alpha1,e1;	

	// read input
	char *finp=argv[1];
	FILE *fpip = fopen(finp,"r");
	char line[1265];
	fgets(line, sizeof(line), fpip);
	fgets(line, sizeof(line), fpip);
	sscanf(line,"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",&a1,&e1,&alpha1,&a2,&m0,&m1,&m2,&tmax);
	fgets(line, sizeof(line), fpip);
	fgets(line, sizeof(line), fpip);
	sscanf(line,"%lf,%lf\n",&J2,&R0);
	fgets(line, sizeof(line), fpip);
	fgets(line, sizeof(line), fpip);
	sscanf(line,"%lf,%lf,%lf\n",&rate_e1,&rate_a1,&rate_a2);
	fgets(line, sizeof(line), fpip);
	fgets(line, sizeof(line), fpip);
	line[strlen(line)-1]='\0';
	FILE *fpout=fopen(line,"w");
	fclose(fpip);

	//a2=get_start_a2(a1,e1);
	//printf("a2: %e\n",a2);
	//exit(0);
	n2 = sqrt(G*(m2+m0+m1)/(a2*a2*a2));
	integrateGSL(a1,e1,alpha1,tmax,fpout);
	fclose(fpout);

}
