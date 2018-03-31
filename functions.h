#include "param.h"

//
//  interact two particles
//
extern inline void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    // this is a little stupid  since it applies the force in only one direction
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoffSQ )
        return;

    double r = sqrt( r2 );

	if (r2 != 0)
    {
	   if (r2/(cutoffSQ) < *dmin * (*dmin))
       {
	      *dmin = r/cutoff;
       }
           (*davg) += r/cutoff;
           (*navg) ++;
    }

    r2 = fmax( r2, min_r_SQ );
    r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}