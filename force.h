#include "param.h"
//
//  interact two particles
//
extern inline void apply_force_SOA( particle_SOA_t *p,int I, int J, double *dmin, double *davg, int *navg)
{

    // this is a little stupid  since it applies the force in only one direction
    // double dx = neighbor.x - particle.x;
    // double dy = neighbor.y - particle.y;

    double dx = p->x[J] - p->x[I];
    double dy = p->y[J] - p->y[I];

    double r2 = dx * dx + dy * dy;

    if( r2 > cutoffSQ )
    {
        //printf(" R2 = %f, dx = %f, dy = %f, I= %d, J= %d \n", r2, dx, dy, I, J);
        return;
    }

    double r = sqrt( r2 );
    double dist = r/cutoff;

    if (dist < *dmin)
    {
      *dmin = dist;
    }

    (*davg) += dist;
    (*navg) ++;

    r2 = fmax( r2, min_r_SQ);
    r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    // but do both at the same time!!!!
    //double coef = ( 1 - cutoff / r ) / r2 / mass;
    double numer = r - cutoff;
    double denom = r * r2 * mass;
    double coef = numer / denom;

    double accelX = coef * dx;
    double accelY = coef * dy;

    p->ax[I] += accelX;
    p->ay[I] += accelY;
    p->ax[J] -= accelX;  // force applied in opposite direction
    p->ay[J] -= accelY;  // force applied in opposite direction
}
