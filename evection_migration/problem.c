/**
 * Planetary migration in the GJ876 system
 *
 * This example applies dissipative forces to two
 * bodies orbiting a central object. The forces are specified
 * in terms of damping timescales for the semi-major axis and
 * eccentricity. This mimics planetary migration in a protostellar disc. 
 * The example reproduces the study of Lee & Peale (2002) on the 
 * formation of the planetary system GJ876. For a comparison, 
 * see figure 4 in their paper. The IAS15 or WHFAST integrators
 * can be used. Note that the forces are velocity dependent.
 * Special thanks goes to Willy Kley for helping me to implement
 * the damping terms as actual forces. 
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "rebound.h"

double* tau_a;     /**< Migration timescale in years for all particles */
double* tau_e;     /**< Eccentricity damping timescale in years for all particles */

const double J2planet			= 0.001; 	// J2 of Saturn (Murray and Dermott p 531) 
const double Mplanet			= 10.0; 	// mass of Saturn in solar masses 
const double Rplanet			= 1.0; 	// radius of Saturn in AU
const double ObliquityPlanet		= 0.;			// obliquity of the planet

const double tmax			= 100000000;			// Maximum integration time

void migration_forces(struct reb_simulation* r);
void heartbeat(struct reb_simulation* r);
void force_J2(struct reb_simulation* r);
FILE *fp;

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_create_simulation();
    // Setup constants
    r->integrator    = REB_INTEGRATOR_WHFAST;
    r->dt         = 1e-2*2.*M_PI;        // in year/(2*pi)
    r->additional_forces = force_J2;     //Set function pointer to add dissipative forces.
    r->heartbeat = heartbeat;  
    r->force_is_velocity_dependent = 1;

    // Initial conditions
    // Parameters are those of Lee & Peale 2002, Figure 4. 
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
    r->dt=2*M_PI*sqrt(a1*a1*a1/(r->G*Mplanet))/50;
    printf("DT: %e\n",r->dt);

    a2=100000; e2=0.0001; i2=0.00*M_PI/180.0;omega2=0.00*M_PI/180.0;Omega2=0.0;f2=90.0*M_PI/180;m2=1e6;
    struct reb_particle star = reb_tools_orbit_to_particle_err(r->G,planet,m2,a2,e2,i2,Omega2,omega2,f2,&err); 
    reb_add(r, star); 
	
    reb_move_to_com(r);

    system("rm -v a.txt");					// delete previous output
    fp = fopen("a.txt","w");
    fprintf(fp,"t,a,e,omega,Omega,res,M2\n");
    fclose(fp);
    
    tau_a = calloc(sizeof(double),r->N);
    tau_e = calloc(sizeof(double),r->N);

    tau_a[1] = 2.*M_PI*1e9;    // Migration timescale of planet 2 is 20000 years.
    tau_e[1] = 2.*M_PI*0.0;     // Eccentricity damping timescale is 200 years (K=100). 

    reb_move_to_com(r);          

    fp = fopen("a.txt","a");
    reb_integrate(r, 2*M_PI*tmax);
    fclose(fp);
}

void migration_forces(struct reb_simulation* r){
    const double G = r->G;
    const int N = r->N;
    struct reb_particle* const particles = r->particles;
    struct reb_particle com = particles[0]; // calculate migration forces with respect to center of mass;
    for(int i=1;i<N;i++){
        if (tau_e[i]!=0||tau_a[i]!=0){
            struct reb_particle* p = &(particles[i]);
            const double dvx = p->vx-com.vx;
            const double dvy = p->vy-com.vy;
            const double dvz = p->vz-com.vz;

            if (tau_a[i]!=0){     // Migration
                p->ax +=  dvx/(2.*tau_a[i]);
                p->ay +=  dvy/(2.*tau_a[i]);
                p->az +=  dvz/(2.*tau_a[i]);
            }
            if (tau_e[i]!=0){     // Eccentricity damping
                const double mu = G*(com.m + p->m);
                const double dx = p->x-com.x;
                const double dy = p->y-com.y;
                const double dz = p->z-com.z;

                const double hx = dy*dvz - dz*dvy; 
                const double hy = dz*dvx - dx*dvz;
                const double hz = dx*dvy - dy*dvx;
                const double h = sqrt ( hx*hx + hy*hy + hz*hz );
                const double v = sqrt ( dvx*dvx + dvy*dvy + dvz*dvz );
                const double r = sqrt ( dx*dx + dy*dy + dz*dz );
                const double vr = (dx*dvx + dy*dvy + dz*dvz)/r;
                const double ex = 1./mu*( (v*v-mu/r)*dx - r*vr*dvx );
                const double ey = 1./mu*( (v*v-mu/r)*dy - r*vr*dvy );
                const double ez = 1./mu*( (v*v-mu/r)*dz - r*vr*dvz );
                const double e = sqrt( ex*ex + ey*ey + ez*ez );        // eccentricity
                const double a = -mu/( v*v - 2.*mu/r );            // semi major axis
                const double prefac1 = 1./(1.-e*e) /tau_e[i]/1.5;
                const double prefac2 = 1./(r*h) * sqrt(mu/a/(1.-e*e))  /tau_e[i]/1.5;
                p->ax += -dvx*prefac1 + (hy*dz-hz*dy)*prefac2;
                p->ay += -dvy*prefac1 + (hz*dx-hx*dz)*prefac2;
                p->az += -dvz*prefac1 + (hx*dy-hy*dx)*prefac2;
            }
        }
        com = reb_get_com_of_pair(com,particles[i]);
    }
}

void force_J2(struct reb_simulation* r){
	if (J2planet==0) 
	{
		migration_forces(r);
		return;
	}
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
	migration_forces(r);

}

void heartbeat(struct reb_simulation* r){
	
	if(reb_output_check(r, 4000.*r->dt)){				// output something to screen	
		reb_output_timing(r, tmax*2*M_PI);
	}
	if(reb_output_check(r,r->dt*40000.0)){				// output semimajor axis to file
		const struct reb_particle planet = r->particles[0];
		const struct reb_particle star = r->particles[2];
		struct reb_orbit ostar = reb_tools_particle_to_orbit(r->G, planet,star);
		int i=1;
		struct reb_orbit o = reb_tools_particle_to_orbit(r->G, r->particles[i],planet);
		double res=ostar.M-o.omega-o.Omega;
		fprintf(fp,"%.15e,%.15e,%.15e,%.15e,%e,%e,%e\n",r->t,o.a,o.e,o.omega,o.Omega,res,ostar.M);	
	}
}

