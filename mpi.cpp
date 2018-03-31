#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <vector>
#include <set>

#include "common.h"
#include "functions.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
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
    
    MPI_Datatype PARTICLE;
    MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
    MPI_Type_commit( &PARTICLE );
    
    //
    //  set up the data partitioning across processors
    //
    int particle_per_proc = (n + n_proc - 1) / n_proc;
    int *partition_offsets = (int*) malloc( (n_proc+1) * sizeof(int) );
    for( int i = 0; i < n_proc+1; i++ )
        partition_offsets[i] = min( i * particle_per_proc, n );
    
    int *partition_sizes = (int*) malloc( n_proc * sizeof(int) );
    for( int i = 0; i < n_proc; i++ )
        partition_sizes[i] = partition_offsets[i+1] - partition_offsets[i];
    
    //
    //  allocate storage for local partition
    //
    int nlocal = partition_sizes[rank];
    particle_t *local = (particle_t*) malloc( nlocal * sizeof(particle_t) );
    
    //
    //  initialize and distribute the particles (that's fine to leave it unoptimized)
    //
    set_size( n );
    if( rank == 0 )
        init_particles( n, particles );
    //MPI_Scatterv( particles, partition_sizes, partition_offsets, PARTICLE, local, nlocal, PARTICLE, 0, MPI_COMM_WORLD );
    
    int NumofBinsEachSide = getbinNumber();
    int NumofBins = NumofBinsEachSide*NumofBinsEachSide;

    /// Make the vector of vectors. The bins are vectors
    std::vector< std::vector<int> > Bins(NumofBins, std::vector<int>(0));

    //printf("Num of bins each side: %d size is: %f BinSize: %f \n", getbinNumber(), getSize(), getBinSize());
    
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
    for( int step = 0; step < NSTEPS; step++ )
    {
        navg = 0;
        dmin = 1.0;
        davg = 0.0;
        
        // Clear bins
        for(int clear = 0; clear < NumofBins; clear++ )
        {
            Bins[clear].clear();
        }

        std::set<int> BinsWithParticles;

        for(int particleIndex = 0; particleIndex < n; ++particleIndex)
        {
            double binsize = getBinSize();

            int BinX = (int)(particles[particleIndex].x/binsize);
            int BinY = (int)(particles[particleIndex].y/binsize);

            //printf("Adding particle\n");
            int BinNum = BinX + NumofBinsEachSide*BinY;
            //printf("Particle added to Bin %d", BinNum);
            Bins[BinNum].push_back(particleIndex);

            // store the bin which contain a particle. We will ignore the empty ones
            BinsWithParticles.insert(BinNum);
        }
        
        // store the particle indices from each surrounding bin.
        std::vector<int> BinMembers;

        // for each bin apply forces from particles within the cutoff distance.
        std::set<int>::iterator it;
        for( it = BinsWithParticles.begin(); it != BinsWithParticles.end(); it++ )
        {
            // this is a vector of the bins with particles. We don't care about empty bins.
            int BinIndex = *it;

            //determine if the bin is in the four courners or edges.
            bool Left = ((BinIndex%NumofBinsEachSide) == 0) ? true : false;
            bool Right = ((BinIndex%NumofBinsEachSide) == (NumofBinsEachSide-1) ) ? true : false;
            bool Top =  ((BinIndex < NumofBinsEachSide) )? true : false;
            bool Bottom = ((BinIndex > (NumofBins - NumofBinsEachSide - 1) ) )? true : false;

            // calculate the indexes of surrounding bins.
            /// <<West    North ^  East >>
            int North = BinIndex - NumofBinsEachSide;
            int NorthEast = North + 1;
            int NorthWest = North -1;
            int East = BinIndex + 1;
            int West = BinIndex -1;
            int South = BinIndex + NumofBinsEachSide;
            int SouthEast = South +1;
            int SouthWest = South -1;

            // we are not at an edge or a corner. -- Most common case
            if( (Left | Right | Top | Bottom) == false)
            {
                //NumofPeerBins = 8; N NE NW E W S SE SW
                //printf("Test1 %d ", BinIndex);
                BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                // printf("Test2");
                BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                // printf("Test3");
                BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                // printf("Test4");
                BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                // printf("Test5");
                BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                // printf("Test6");
                BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin() ,Bins[SouthWest].end());
                //printf("Test7");
                BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                //printf("Test8");
                BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
            }

            // Top Row
            else if( Top )
            {
                // most common case for the top row -- Not in a corner.
                if( (Left | Right) == false)
                {
                    //printf("Top Row %d \n", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
                    //printVector(BinMembers);
                }
                else if( (!Left) && Right ) // Yes this would be called a corner case!!!
                {
                    // Right == East
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                }
                else if ( Left && (!Right) )// left corner
                {
                    //printf("Top Row Left %d \n", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
                }

            }
            else if( Bottom )
            {
                // most common case for the top row -- Not in a corner.
                if((Left | Right) == false)
                {
                    //printf("Bottom Row %d \n", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                    BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                    //printVector(BinMembers);
                }
                else if( (!Left) && Right ) // Yes this would be called a corner case!!!
                {
                    // Right == East
                    //printf("Bottom Row Right %d \n", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    //printVector(BinMembers);
                }
                else if( Left && (!Right) )// left corner
                {
                    //printf("Bottom Row Left %d ", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                    BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                    //printVector(BinMembers);
                }

            }
            else if(Left)  // not in a corner on the left side
            {
                //printf("Left %d \n", BinIndex );
                BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
                //printVector(BinMembers);
            }
            else if(Right) // must be the right side
            {
                //printf("Right %d \n", BinIndex );
                BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                //printVector(BinMembers);
            }
            else
            {
              printf("Getting another bin case\n");
            }

            // apply force within the bin
            for(int Inside_Bin = 0; Inside_Bin < Bins[BinIndex].size(); Inside_Bin++)
            {
                int Index = Bins[BinIndex][Inside_Bin];

                for(int Inside_BinJ =0; Inside_BinJ < Bins[BinIndex].size(); Inside_BinJ++)
                {
                    int Index2 = Bins[BinIndex][Inside_BinJ];
                    apply_force( particles[Index], particles[Index2], &dmin, &davg, &navg);
                }

            }
        
            // force to neighbor bins
            for(int calcForceindexI = 0; calcForceindexI < Bins[BinIndex].size(); calcForceindexI++ )
            {
                int ParticleThisBin = Bins[BinIndex][calcForceindexI];

                for (int calcForceindexJ = 0; calcForceindexJ < BinMembers.size(); calcForceindexJ++ )
                {
                    apply_force( particles[ParticleThisBin], particles[ BinMembers[calcForceindexJ] ], &dmin, &davg, &navg);
                }
            }
        
            BinMembers.clear();
        }   // ends BinsWithParticles iterator loop
        
        BinsWithParticles.clear();

        // 
        //  collect all global data locally (not good idea to do)
        //
        //MPI_Allgatherv( local, nlocal, PARTICLE, particles, partition_sizes, partition_offsets, PARTICLE, MPI_COMM_WORLD );
        
        //
        //  save current step if necessary (slightly different semantics than in other codes)
        //
        if( find_option( argc, argv, "-no" ) == -1 )
          if( fsave && (step%SAVEFREQ) == 0 )
            save( fsave, n, particles );
        
        //
        //  compute all forces
        //
        /*
        for( int i = 0; i < nlocal; i++ )
        {
            local[i].ax = local[i].ay = 0;
            for (int j = 0; j < n; j++ )
                apply_force( local[i], particles[j], &dmin, &davg, &navg );
        }
     */
        if( find_option( argc, argv, "-no" ) == -1 )
        {
          
          MPI_Reduce(&davg,&rdavg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&navg,&rnavg,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
          MPI_Reduce(&dmin,&rdmin,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);

 
          if (rank == 0){
            //
            // Computing statistical data
            //
            if (rnavg) {
              absavg +=  rdavg/rnavg;
              nabsavg++;
            }
            if (rdmin < absmin) absmin = rdmin;
          }
        }

        //
        //  move particles
        //
        for( int i = 0; i < nlocal; i++ )
            move( local[i] );
    }
    simulation_time = read_timer( ) - simulation_time;
  
    if (rank == 0) {  
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
        fprintf(fsum,"%d %d %g\n",n,n_proc,simulation_time);
    }
  
    //
    //  release resources
    //
    if ( fsum )
        fclose( fsum );
    free( partition_offsets );
    free( partition_sizes );
    free( local );
    free( particles );
    if( fsave )
        fclose( fsave );
    
    MPI_Finalize( );
    
    return 0;
}
