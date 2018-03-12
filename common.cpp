#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"
#include "param.h"
#include <immintrin.h>

double size;

//
//  timer
//
double read_timer( )
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//
void init_particles_SOA( int n, particle_SOA_t *p )
{
    srand48( time( NULL ) );

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;

    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;

    for( int i = 0; i < n; i++ )
    {
        //
        //  make sure particles are not spatially sorted
        //
        int j = lrand48()%(n-i);
        int k = shuffle[j];
        shuffle[j] = shuffle[n-i-1];

        //
        //  distribute particles evenly to ensure proper spacing
        //
        // p[i].x = size*(1.+(k%sx))/(1+sx);
        // p[i].y = size*(1.+(k/sx))/(1+sy);

        p->x[i] = size*(1.+(k%sx))/(1+sx);
        p->y[i] = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        // p[i].vx = drand48()*2-1;
        // p[i].vy = drand48()*2-1;

        p->vx[i] = drand48()*2-1;
        p->vy[i] = drand48()*2-1;
    }


    free( shuffle );
}



/*

void apply_force_vector_4( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{



    const double SubCutoff =  1 - cutoff;
    const double MassConst = mass;

    __m256d ParticleX4, ParticleY4, NeighborX4, NeighborY4, dx4, dy4, r2_4, r_4, SubCutoff4, Mass4, ParticleAx4,ParticleAy4;

    SubCutoff4 = _mm256_broadcast_sd(&SubCutoff);
    Mass4  = _mm256_broadcast_sd(&MassConst);


    //    double dx = neighbor.x - particle.x;
    //    double dy = neighbor.y - particle.y;

    ParticleX4 = _mm256_set_pd((particle).x,(particle+1).x,(particle+2).x,(particle+3).x);
    ParticleY4 = _mm256_set_pd(particle[0].y,particle[1].y,particle[2].y,particle[3].y);
    NeighborX4 = _mm256_set_pd(neighbor[0].x,neighbor[1].x,neighbor[2].x,neighbor[3].x);
    NeighborY4 = _mm256_set_pd(neighbor[0].y,neighbor[1].y,neighbor[2].y,neighbor[3].y);

    dx4 = _mm256_sub_pd(NeighborX4, ParticleX4);
    dy4 = _mm256_sub_pd(NeighborY4, ParticleY4);

    // calculate r squared
    r2_4 = _mm256_sub_pd(_mm256_mul_pd(dx4,dx4),_mm256_mul_pd(dy4,dy4));
    // get r
    r_4 = _mm256_sqrt_pd(r2_4);


    // if( r2 > cutoffSQ )
    //     return;

    // double r = sqrt( r2 );

    // if (r2 != 0)
    // {
    //    if (r2/(cutoffSQ) < *dmin * (*dmin))
    //    {
    //       *dmin = r/cutoff;
    //    }
    //        (*davg) += r/cutoff;
    //        (*navg) ++;
    // }

    // r2 = fmax( r2, min_r_SQ );


    // r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    // double coef = ( SubCutoff / r ) / r2 / mass;
    coef4 = _mm256_div_pd(_mm256_div_pd(_mm256_div_pd(SubCutoff4,r_4),r2_4),Mass4);

    // particle.ax += coef * dx;
    // particle.ay += coef * dy;
    ParticleAx4 = _mm256_add_pd(ParticleAx4,_mm256_mul_pd(coef4, dx4));
    ParticleAy4 = _mm256_add_pd(ParticleAy4,_mm256_mul_pd(coef4, dy4));
}*/


void move_SOA( particle_SOA_t *p,int I)
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p->vx[I] += p->ax[I] * dt;
    p->vy[I] += p->ay[I] * dt;
    p->x[I]  += p->vx[I] * dt;
    p->y[I]  += p->vy[I] * dt;

    // once the force is applied awesome the accel is zero.
    p->ax[I] = 0;
    p->ay[I] = 0;




    // p.vx[J] += p.ax[J] * dt;
    // p.vy[J] += p.ay[J] * dt;
    // p.x[J]  += p.vx[J] * dt;
    // p.y[J]  += p.vy[J] * dt;

    //
    //  bounce from walls
    //
    while( p->x[I] < 0 || p->x[I] > size )
    {
        p->x[I]  = p->x[I] < 0 ? -p->x[I] : 2*size-p->x[I];
        p->vx[I] = -p->vx[I];
    }
    while( p->y[I] < 0 || p->y[I] > size )
    {
        p->y[I]  = p->y[I] < 0 ? -p->y[I] : 2*size-p->y[I];
        p->vy[I] = -p->vy[I];
    }
    // while( p.x[J] < 0 || p.x[J] > size )
    // {
    //     p.x[J]  = p.x[J] < 0 ? -p.x[J] : 2*size-p.x[J];
    //     p.vx[J] = -p.vx[J];
    // }
    // while( p.y[J] < 0 || p.y[J] > size )
    // {
    //     p.y[J]  = p.y[J] < 0 ? -p.y[J] : 2*size-p.y[J];
    //     p.vy[J] = -p.vy[J];
    // }

}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  I/O routines
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void save_SOA( FILE *f, int n, particle_SOA_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
    {
        fprintf( f, "%g %g\n", p->x[i], p->y[i]);
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  command line option processing
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
