#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define cutoffSQ (cutoff*cutoff)
#define INVCutoff (1/cutoff)
#define min_r   (cutoff/100)
#define min_r_SQ (min_r*min_r)
#define dt      0.0005

//
//  saving parameters
//
const int NSTEPS = 1000;
const int SAVEFREQ = 10;

//
// particle data structure
//
typedef struct 
{
  double x;    /// position X
  double y;    /// position y 
  double vx;   /// velocity x
  double vy;   /// velocity y
  double ax;   /// accel in x
  double ay;   // accel in y
} particle_t;

//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
//void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );

inline void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;

    // this takes 17 percent of the time. 
    double r2 = dx * dx + dy * dy;   
    double r = sqrt( r2 );
    double rInvCutoff = r*INVCutoff;

    if( r2 > cutoffSQ )
        return;
	if (r2 != 0)
    {
	   if (r2/(cutoffSQ) < *dmin * (*dmin))
       {
	      *dmin = rInvCutoff;
       }
           (*davg) += rInvCutoff;
           (*navg) ++;
    }
		
    r2 = fmax( r2, min_r_SQ);
    r = sqrt( r2 );
 
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;


    particle.ax += coef * dx;
    particle.ay += coef * dy;
}


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );


#endif
