#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "common.h"
#include <vector>
//#include "mpi_tools.h"

MPI_Datatype PARTICLE;
void ScatterParticlesToProcs(particle_t *particles, const int NumofParticles, const int NumofBinsEachSide, const int NumberofProcessors, const int rank, particle_t * local ,int * nlocal);

//
//  benchmarking program
//
int main(int argc, char **argv)
{
    int navg, nabsavg=0;
    double dmin, absmin=1.0,davg,absavg=0.0;
    double rdavg,rdmin;
    int rnavg; 

    //
    //  process command line parameters
    //
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

    //This will probably stay the same. 
    set_size( n );
    
    //
    //  set up MPI
    //
    int n_proc, rank;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &n_proc );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    
    //
    //  allocate generic resources
    //
    FILE *fsave = savename && rank == 0 ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname && rank == 0 ? fopen ( sumname, "a" ) : NULL;

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// probably have to set stuff up in here 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // //
    // //  set up the data partitioning across processors
    // //
    // int particle_per_proc = (n + n_proc - 1) / n_proc;
    // int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    // for( int i = 0; i < n_proc+1; i++ )
    //     partition_offsets[i] = min( i * particle_per_proc, n );
    
    // int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    // for( int i = 0; i < n_proc; i++ )
    //     partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //


    if (rank == 0)
    { // if we are the master node. 
        init_particles(n, particles);

    }

    int NumofBinsEachSide = getbinNumber();
    int NumofBins = NumofBinsEachSide*NumofBinsEachSide;
    // allocate storage for local partition
    //
    int *nlocal;
    particle_t *local; //  = (particle_t*) malloc( nlocal * sizeof(particle_t) );

    /// Make the vector of vectors. The bins are vectors
    std::vector< std::vector<int> > Bins(NumofBins, std::vector<int>(0));

    // send the assign the particles to each procssor based on which bin they are located in. local particles will be populated from this array. 
    ScatterParticlesToProcs(particles, n, NumofBinsEachSide, n_proc,rank, local ,nlocal);


    //MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 // probaly going to have to change this as well 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer();
    for (int step = 0; step < NSTEPS; step++)
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;


        //  // 
        // //  collect all global data locally (not good idea to do)
        // //
        // MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        
        // //
        // //  save current step if necessary (slightly different semantics than in other codes)
        // //
        // if( find_option( argc, argv, "-no" ) == -1 )
        //   if( fsave && (step%SAVEFREQ) == 0 )
        //     save( fsave, n, particles );
        
        // //
        // //  compute all forces
        // //
        // for( int i = 0; i < nlocal; i++ )
        // {
        //     local[i].ax = local[i].ay = 0;
        //     for (int j = 0; j < n; j++ )
        //         apply_force( local[i], particles[j], &dmin, &davg, &navg );
        // }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// new apply forces function that works across processors. 

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        if( find_option( argc, argv, "-no" ) == -1 )
        { //This should stay the same 
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

          if (rank == 0)
          {
            //
            // Computing statistical data
            //
            if (rnavg) 
            {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin)
            { 
                absmin = rdmin;
            }
          }
        }

    //             //
    //     //  move particles
    //     //
    //     for( int i = 0; i < nlocal; i++ )
    //         move( local[i] );
    // }

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        ////new move function that handles cross processor particles 
}
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    simulation_time = read_timer() - simulation_time;

    if (rank == 0) 
    {
        printf("n = %d, simulation time = %g seconds", n, simulation_time);

        if (find_option(argc, argv, "-no") == -1) 
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
            printf(", absmin = %lf, absavg = %lf", absmin, absavg);
            if (absmin < 0.4) printf("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
            if (absavg < 0.8) printf("\nThe average distance is below 0.8 meaning that most particles are not interacting");
        }
        printf("\n");

        //
        // Printing summary data
        //
        if (fsum)
        {
            fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
        }
    }

    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    //free( partition_offsets );
    //free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}




// would be much cleaner to use multiple returns in C++ 14 but I am not placing faith in compiler support. It's a very new feature
void ScatterParticlesToProcs(particle_t *particles, const int NumofParticles, const int NumofBinsEachSide, const int NumberofProcessors, const int rank, particle_t * local ,int * nlocal)
{ // send particles to the correct processor based on which bin the particle is in and what processor the bin belongs to. 


    //int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    std::vector<int> partition_offsets(NumberofProcessors);
    
    //int *partition_sizes = (int*) malloc( NumberofProcessors * sizeof(int) );
    
    std::vector<int> partition_sizes(NumberofProcessors);

    // temp storage for particles in correct order
    nlocal = (int *) malloc( NumberofProcessors * sizeof(int) );
    *nlocal = partition_sizes[rank];
    local = (particle_t*) malloc( *nlocal * sizeof(particle_t) );


    std::vector< std::vector<particle_t> > ParticlesPerProcecesor(NumberofProcessors,std::vector<particle_t>() );

    if(rank == 0) // only the master node should do this 
    {   //sort the particles based on the bin location 
        //returns a proc number 
        for(int particleIndex = 0; particleIndex < NumofParticles; particleIndex++)
        {
            int ProcNumber = MapParticleToProc(particles[particleIndex],NumofBinsEachSide,NumberofProcessors);
            ParticlesPerProcecesor[ProcNumber].push_back(particles[particleIndex]); // add the the vector at each index

        }

    }

    std::vector<particle_t>  particlesassignedtoproc;

    // append all the particles to particlesassignedtoproc 0 to NumberofProcessors
    int ProcNumCount = 0;
    int runningoffset = 0;

    for(auto ParticleVector : ParticlesPerProcecesor)
    {
        particlesassignedtoproc.insert(particlesassignedtoproc.end(),ParticleVector.begin(),ParticleVector.end());
        partition_sizes[ProcNumCount] = ParticleVector.size();
        partition_offsets[ProcNumCount] = runningoffset; 
        runningoffset += ParticleVector.size();
        ProcNumCount +=1; 
    }

    *nlocal = partition_sizes[rank];

    local = (particle_t*) malloc( *nlocal * sizeof(particle_t) );

    // scatter the particles to the processors. More scattered than the programmer's brain. 
    MPI_Scatterv( particlesassignedtoproc.data(), partition_sizes.data(), partition_offsets.data(), PARTICLE, local, *nlocal, PARTICLE, 0, MPI_COMM_WORLD );

    // pharaoh, let my particles go!!!!
    // free(partition_offsets);
    // free(partition_sizes);
    //free(particlesassignedtoproc);
}


