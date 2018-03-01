#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "omp.h"

particle_t *particles;
 //setup simulation 
int navg,nabsavg=0,numthreads; 
double dmin, absmin=1.0,davg,absavg=0.0;
int MasterSave = 0;
FILE *fsave;
FILE *fsum;


void RunSimulation(particle_t * particles,int n);
void SetStats(particle_t * particles, const int n,const int step);

//
//  benchmarking program
//
int main( int argc, char **argv )
{   
    ///////////////////// PROCESS USER INPUT /////////////////////////////////////////////////////////////////////////////////////////////////
    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" ); 
        printf( "-no turns off all correctness checks and particle output\n");   
        return 0;
    }


    if( find_option( argc, argv, "-no" ) != -1 ) 
    {     
        MasterSave = 1;
    }

    int n = read_int( argc, argv, "-n", 1000 );
    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );

    fsave = savename ? fopen( savename, "w" ) : NULL;
    fsum = sumname ? fopen ( sumname, "a" ) : NULL; 

    ////////////////////////////////////////////////////////SIMULATION SETUP ////////////////////////////////////////////////////////////////////////////////////////////

    // allocate memory for our wonderful particles
    //particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    particles = (particle_t*) malloc( n * sizeof(particle_t) );

    ////////////////////////////////////////////////////// INIT SIMULATION /////////////////////////////////////////////////////////////////////////////////////////////

    set_size( n );
    init_particles( n, particles );

    ////////////////////////////////////////////////////// RUN SIMULATION ////////////////////////////////////////////////////////////////////////////////////////////////

    //  simulate a number of time steps
    //double simulation_time = read_timer( );
    double simulation_time = omp_get_wtime();

    numthreads = omp_get_num_threads();
    // Simulation code was in here 
    RunSimulation(particles,n);

    //////////////////////////////////////////////////////// END OF SIMULATION //////////////////////////////////////////////////////////////////////////////////////////////
    //simulation_time = read_timer( ) - simulation_time;
    simulation_time = omp_get_wtime() - simulation_time;

    ///////////////////////////////////////////////////////FREE SIMULATION MEMORY ////////////////////////////////////////////////////////////////////////////////////////////
    free( particles );

    //////////////////////////////////////////////////////  STATS ////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    printf( "n = %d,threads = %d, simulation time = %g seconds", n,numthreads, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg)
      { 
        absavg /= nabsavg;
      }
    // 
    //  -the minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation were particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
        printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
        if (absmin < 0.4)
        { 
            printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
        }
        if (absavg < 0.8) 
        {
                printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
        }
    }
    printf("\n");
    
    //
    // Printing summary data
    //
    if(fsum)
    {
        fprintf(fsum,"%d %d %g\n",n,numthreads,simulation_time);
    }

    //
    // Clearing space
    //
    if( fsum )
    {
        fclose( fsum );
    }

    if( fsave )
    {
        fclose( fsave );
    }
    
    return 0;
}


void RunSimulation(particle_t * particles,int n)
{
    #pragma omp parallel private(dmin) 
    {
        for( int step = 0; step < 1000; step++ )
        {
            // global vars for some reason 
            navg = 0;
            davg = 0.0;
            dmin = 1.0;
            //shared(particles) 
            //  compute all forces
            #pragma omp for reduction (+:navg) reduction(+:davg) schedule(static) 
            for( int i = 0; i < n; ++i )
            {
                particles[i].ax = particles[i].ay = 0;

                for (int j = 0; j < n; ++j )
                {
                    apply_force( particles[i], particles[j],&dmin,&davg,&navg);
                }
            }
            
            //// need a barrier here to wait for all threads to complete the force updates
            #pragma omp barrier

            // once the forces have been applied 
            #pragma omp for
            for( int i = 0; i < n; ++i ) 
            {
                move( particles[i] );
            }

            //////////////////////////////////////////////////////// set stats ////////////////////////////////////////////////////////
            #pragma omp master
            SetStats(particles,n, step);
        }
    }
}


void SetStats(particle_t * particles, const int n,const int step)
{

      //
      //  compute statistical data
      //
      //#pragma omp master
      if (navg)
      { 
        absavg += davg/navg;
        nabsavg++;
      }

      //#pragma omp critical
      if (dmin < absmin) 
      {
         absmin = dmin; 
      }
      //
      //  save if necessary
      //
      //#pragma omp master
      if( fsave && (step%SAVEFREQ) == 0 )
      {
          save( fsave, n, particles );
      } // if( find_option( argc, argv, "-no" ) == -1 ) 
}
