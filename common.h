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
#include <immintrin.h>
#include <cmath> 
#include <vector>
#include <algorithm>


inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

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
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

//
// particle data structure
// This is a structure of Array. This was done since the access stride was not 1.
//By converting to Structure of Arrays we can get an access stride of 1.
typedef struct
{
	double * x;
	double * y;
	double * vx;
	double * vy;
	double * ax;
	double * ay;

} particle_SOA_t;
//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
int getNumberofBins( int size);
double getSize();
double getBinSize();
int getRowsPerProc(int BinsEachSide, int NumberofProcessors );
std::vector<int> getBoarderPeers(int rank);
int MapBinToProc(int Bin, int NumberOfProcs);
int MapParticleToBin(particle_t &particle, const int NumofBinsEachSide);
int MapParticleToProc(particle_t &particle, const int NumofBinsEachSide, const  int NumberofProcessors );
void init_particles( int n, particle_t *p );
void init_particles_SOA( int n, particle_SOA_t *p );
void apply_force_SOA( particle_SOA_t &p,int I, int J, double *dmin, double *davg, int *navg);
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );
void move_SOA( particle_SOA_t &p,int I);


//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );
void save_SOA( FILE *f, int n, particle_SOA_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );

#endif
