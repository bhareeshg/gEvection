/**
 * J2 precession
 * 
 * This example presents an implementation of the J2 gravitational moment. 
 * The equation of motions are integrated with the 15th order IAS15 
 * integrator. The parameters in this example have been chosen to 
 * represent those of Saturn, but one can easily change them or even 
 * include higher order terms in the multipole expansion.
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

const double J2planet			= 0.001; 	// J2 of Saturn (Murray and Dermott p 531) 
const double Mplanet			= 10.0; 	// mass of Saturn in solar masses 
const double Rplanet			= 1.0; 	// radius of Saturn in AU
const double ObliquityPlanet		= 0.;			// obliquity of the planet

const double tmax			= 10000000;			// Maximum integration time

void heartbeat(struct reb_simulation* r);
void force_J2(struct reb_simulation* r);
FILE* fp; 
int main(int argc, char* argv[]){
	struct reb_simulation* r = reb_create_simulation();
	// Setup constants
	r->integrator			= REB_INTEGRATOR_WHFAST;
	r->dt 				= 1e-6;			// initial timestep
	r->N_active			= 3; 			// only the star and the planet are massive.
	
	// Planet
	struct reb_particle planet = {0};
	planet.m  = Mplanet;
	reb_add(r, planet);

	double a1,e1,i1,omega1,Omega1,f1,m1;
	double a2,e2,i2,omega2,Omega2,f2,m2;
	
	int err;
	a1=4.6356718543987965;
	e1=0.3940047725058190429; i1=0.00*M_PI/180.0;omega1=0.00*M_PI/180.0;Omega1=0.00;f1=0.01*M_PI/180;m1=0.000;
    	struct reb_particle p = reb_tools_orbit_to_particle_err(r->G,planet,m1,a1,e1,i1,Omega1,omega1,f1,&err); 
	reb_add(r, p); 
	r->dt=2*M_PI*sqrt(a1*a1*a1/(r->G*Mplanet))/100;
	printf("DT: %e\n",r->dt);

	a2=100000; e2=0.0001; i2=0.00*M_PI/180.0;omega2=0.00*M_PI/180.0;Omega2=0.0;f2=90.0*M_PI/180;m2=1e6;
    	struct reb_particle star = reb_tools_orbit_to_particle_err(r->G,planet,m2,a2,e2,i2,Omega2,omega2,f2,&err); 
	reb_add(r, star); 
	
	reb_move_to_com(r);

	system("rm -v a.txt");					// delete previous output

	// Setup callback functions
	r->heartbeat 		= heartbeat;
	r->additional_forces 	= force_J2;

	FILE* f = fopen("a.txt","w");
	fprintf(f,"t,a,e,omega,Omega,res,M2\n");
        fclose(f);
	fp = fopen("a.txt","a");
	reb_integrate(r, 2*M_PI*tmax);
	fclose(fp);

}

void force_J2(struct reb_simulation* r){
	if (J2planet==0) return;
	// Star 
	const struct reb_particle planet = r->particles[0];		// cache
	const int N = r->N;
	for (int i=1;i<2;i++){
		const struct reb_particle p = r->particles[i]; 	// cache
		const double sprx = p.x-planet.x;
		const double spry = p.y-planet.y;
		const double sprz = p.z-planet.z;
		const double prx  = sprx*cos(-ObliquityPlanet) + sprz*sin(-ObliquityPlanet);
		const double pry  = spry;
		const double prz  =-sprx*sin(-ObliquityPlanet) + sprz*cos(-ObliquityPlanet);
		const double pr2  = prx*prx + pry*pry + prz*prz; 		// distance^2 relative to planet
		const double fac  = -3.*r->G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);

		const double pax  = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
		const double pay  = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
		const double paz  = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);
		
		r->particles[i].ax += pax*cos(ObliquityPlanet) + paz*sin(ObliquityPlanet);
		r->particles[i].ay += pay;
		r->particles[i].az +=-pax*sin(ObliquityPlanet) + paz*cos(ObliquityPlanet);
	}
}

void heartbeat(struct reb_simulation* r){
	
	
	/*struct reb_particle planet=r->particles[0];
	int err;
    	struct reb_particle pert = r->particles[1]; 
	struct reb_orbit orbitp=reb_tools_particle_to_orbit_err(r->G,pert,planet,&err);
	double ap=orbitp.a,ep=orbitp.e,ip=orbitp.inc,omegap=orbitp.omega,Omegap=orbitp.Omega,f1p=orbitp.f;
	double mp=0;
	double vmin=4.64*Rplanet;
	double vmax=13.36*Rplanet;
	ap=vmin+(vmax-vmin)*(r->t/(2*M_PI)/tmax);
	//printf("AIMMEDIATE: %e\n",ap/Rplanet);
	r->particles[1]  = reb_tools_orbit_to_particle_err(r->G,planet,mp,ap,ep,ip,Omegap,omegap,f1p,&err); 
	*/

	if(reb_output_check(r, 4000.*r->dt)){				// output something to screen	
		reb_output_timing(r, tmax*2*M_PI);
	}
	if(reb_output_check(r,r->dt*40000.0)){				// output semimajor axis to file
		const struct reb_particle planet = r->particles[0];
		const struct reb_particle star = r->particles[2];
		struct reb_orbit ostar = reb_tools_particle_to_orbit(r->G, planet,star);
		//const int N = r->N;
		int i=1;
		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],planet);
		double res=ostar.M-o.omega-o.Omega;
		fprintf(fp,"%.15e,%.15e,%.15e,%.15e,%e,%e,%e\n",r->t,o.a,o.e,o.omega,o.Omega,res,ostar.M);	
	}
}

