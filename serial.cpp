#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "force.h"

#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)

//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    ///////////////////////////////////////////////////////////////// End of input processing ///////////////////////////////////////////////////

    set_size( n );

    particle_SOA_t *particlesSOA =(particle_SOA_t*)malloc( sizeof(particle_SOA_t));
    particlesSOA->x = (double*)malloc( n* sizeof(double));
    particlesSOA->y = (double*)malloc( n* sizeof(double));
    particlesSOA->vx = (double*)malloc( n* sizeof(double));
    particlesSOA->vy = (double*)malloc( n* sizeof(double));
    particlesSOA->ax = (double*)malloc( n* sizeof(double));
    particlesSOA->ay = (double*)malloc( n* sizeof(double));

    init_particles_SOA( n, particlesSOA );

    //
    //  simulate a number of time steps
    //

    int loopcount = 0;
    int count = 0;
    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        davg = 0.0;
        dmin = 1.0;

        //
        //  compute forces
        //
        for( int i = 0; i < n; i++ )
        {
            for (int j = i+1; j < n; j++ )
            {
                loopcount++;
                apply_force_SOA( particlesSOA,i, j, &dmin, &davg, &navg);
    	    }
            move_SOA( particlesSOA,i);
    	}

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;

          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
            //save( fsave, n, particles );
            save_SOA( fsave, n, particlesSOA );
        }
    }

    //printf( "dmin = %f,davg = %f,navg = %f, loops = %d, applyforcescount = %d \n", dmin, davg, navg, loopcount, count );

    simulation_time = read_timer( ) - simulation_time;

    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    //
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");

    //
    // Printing summary data
    //
    if( fsum)
        fprintf(fsum,"%d %g\n",n,simulation_time);

    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );

    free(particlesSOA->x);
    free(particlesSOA->y);
    free(particlesSOA->vx);
    free(particlesSOA->vy);
    free(particlesSOA->ax);
    free(particlesSOA->ay);

    free(particlesSOA);

    if( fsave )
        fclose( fsave );

    return 0;
}
