#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>
#include <set>

// this is for the viz script
//#define VIZ

#ifdef VIZ
  #include <sys/stat.h>
#endif

void printVector(std::vector<int> A)
{
    printf("Vector: ");
    for(int y=0; y< A.size(); y++)
    {
        printf("%d ",A[y]);
    }
    printf("\n");

}
//
//  benchmarking program
//
int main( int argc, char **argv )
{
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    #ifdef VIZ
        char *outname = (char*) malloc(32* sizeof(char));
        const char* dirName = "out";
        mkdir(dirName, 0777);
    #endif

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

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );

    // particle_SOA_t *particlesSOA =(particle_SOA_t*)malloc( sizeof(particle_SOA_t));
    // particlesSOA->x = (double*)malloc( n* sizeof(double));
    // particlesSOA->y = (double*)malloc( n* sizeof(double));
    // particlesSOA->vx = (double*)malloc( n* sizeof(double));
    // particlesSOA->vy = (double*)malloc( n* sizeof(double));
    // particlesSOA->ax = (double*)malloc( n* sizeof(double));
    // particlesSOA->ay = (double*)malloc( n* sizeof(double));

    set_size( n );
    init_particles( n, particles );

    //init_particles_SOA( n, particlesSOA );

    double size = getSize(); 
    int NumofBinsEachSide = getNumberofBins(size);
    int NumofBins = NumofBinsEachSide*NumofBinsEachSide;
    // //std::vector<int> * Bins =( std::vector<int> *)malloc( (NumofBins^2) * sizeof( std::vector<int> ));
    // //std::vector<std:vector<int>> Bins;int
    // printf("NumofBinsEachSide %d \n", NumofBinsEachSide);
    // printf("Num of bins %d \n",NumofBins );

    /// Make the vector of vectors. The bins are vectors
    std::vector< std::vector<int> > Bins(NumofBins, std::vector<int>(0));

    //
    //  simulate a number of time steps
    //

    printf("Num of bins each side: %d size is: %f BinSize: %f \n", getNumberofBins(size), getSize(), getBinSize());

    double simulation_time = read_timer( );

    for( int step = 0; step < NSTEPS; step++ )
    {
       navg = 0;
       davg = 0.0;
       dmin = 1.0;

        for(int clear = 0; clear < NumofBins; clear++ )
        {
            Bins[clear].clear();
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // re populate bins after the move update.

        std::set<int> BinsWithParticles;

        for(int particleIndex = 0; particleIndex < n; ++particleIndex)
        {
            // CHECKED///////////////////////
              double binsize = getBinSize();
              //printf("Test4\n");
              // get the bin index


               int BinX = (int)(particles[particleIndex].x/binsize);
               int BinY = (int)(particles[particleIndex].y/binsize);

               // int BinX = (int)(particlesSOA->x[particle]/binsize);
               // int BinY = (int)(particlesSOA->y[particle]/binsize);

               //printf("Adding particle\n");
               int BinNum = BinX + NumofBinsEachSide*BinY;
               //printf("Particle added to Bin %d", BinNum);
               Bins[BinNum].push_back(particleIndex);

               // store the bin which contain a particle. We will ignore the empty ones
               BinsWithParticles.insert(BinNum);

               // int BinX = (int)(particlesSOA->x[particle]/binsize);
               // int BinY = (int)(particlesSOA->y[particle]/binsize);

               //particles[particleIndex].ax = particles[particleIndex].ay = 0;

               //particlesSOA->ax[particle] = 0;
               //particlesSOA->ay[particle] = 0;


               //printf("P %d X: %f Y: %f BinNum: %d \n", particle,particles[particle].x, particles[particle].y,BinNum );
               //printf("There are %d particles in bin %d\n",Bins[BinNum].size(),BinNum );
        //     //addParticleToBin(n,BinX,BinY);
        }


        //
        //  compute forces
        //
    //     for( int i = 0; i < n; i++ )
    //     {
    //         particles[i].ax = particles[i].ay = 0;
    //         for (int j = 0; j < n; j++ )
                // apply_force( particles[i], particles[j],&dmin,&davg,&navg);
    //     }


        // store the particle indices from each surrounding bin.
        std::vector<int> BinMembers;

        // for each bin apply forces from particles within the cutoff distance.
        //for(int BinIndex = 0; BinIndex < NumofBins; BinIndex++ )
        std::set<int>::iterator it;
        for( it = BinsWithParticles.begin(); it != BinsWithParticles.end(); it++ )
        {

            // this is a vector of the bins with particles. We don't care about empty bins.
            int BinIndex = *it;

         // if(Bins[BinIndex].size() > 0)
         // {
            //printf("#########################################################\n");
            //printf("Bin number %d \n",BinIndex );
            //printf("#########################################################\n");
            /////////CHECKED !!!
            //determine if the bin is in the four courners or edges.
            bool Left = ((BinIndex%NumofBinsEachSide) == 0) ? true : false;
            bool Right = ((BinIndex%NumofBinsEachSide) == (NumofBinsEachSide-1) ) ? true : false;
            bool Top =  ((BinIndex < NumofBinsEachSide) )? true : false;
            bool Bottom = ((BinIndex > (NumofBins - NumofBinsEachSide - 1) ) )? true : false;

            //printf("BinIndex %d Left %d Right %d Top %d Bottom %d\n",BinIndex,Left,Right,Top,Bottom);

            ///////////CHECKED
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

            //printf(" BinIndex %d North %d NorthEast %d NorthWest %d, East %d, West %d, South %d, SouthEast %d, SouthWest %d\n",BinIndex, North,NorthEast,NorthWest,East, West,South, SouthEast,SouthWest  );

            // Note: Left = West
            // Right = East

            /// NOW! Take all the particles in each of the bins and apply forces across all 9 bins.
            // We are always going to apply forces within a bin we are in
            //BinMembers.insert(BinMembers.end(),Bins[BinIndex].begin(),Bins[BinIndex].end());


            // we are not at an edge or a corner. -- Most common case
            if( (Left | Right | Top | Bottom) == false)
            {
                ///////////////// CHECKED
                //printf("ALL %d \n", BinIndex );
                //NumofPeerBins = 8;
                //N NE NW E W S SE SW
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
                //printVector(BinMembers);
                //printf("Testlast");
            }
            //printf("After");
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
                    //printf("Top Row Right %d \n", BinIndex );
                    // Right == East
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                    //printVector(BinMembers);
                }
                else if ( Left && (!Right) )// left corner
                {
                    //printf("Top Row Left %d \n", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
                    //printVector(BinMembers);
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

            for(int Inside_Bin = 0; Inside_Bin < Bins[BinIndex].size(); Inside_Bin++)
            {
                int Index = Bins[BinIndex][Inside_Bin];
                // particles[Index].ax = particles[Index].ay = 0;

                for(int Inside_BinJ =0; Inside_BinJ < Bins[BinIndex].size(); Inside_BinJ++)
                {
                    int Index2 = Bins[BinIndex][Inside_BinJ];
                    
                    apply_force( particles[Index], particles[Index2], &dmin, &davg, &navg);
                    //apply_force_SOA( particlesSOA,Index, Inside_BinJ, &dmin, &davg, &navg);
                }

            }

            //printf(" There will be %d x %d interactions\n",Bins[BinIndex].size(), BinMembers.size());
            // force to neighbors
            for(int calcForceindexI = 0; calcForceindexI < Bins[BinIndex].size(); calcForceindexI++ )
            {
                // apply forces

                int ParticleThisBin = Bins[BinIndex][calcForceindexI];    ////// PERFECT
                    // once the force is applied awesome the accel is zero.
                // particlesSOA->ax[ParticleThisBin] = 0;
                // particlesSOA->ay[ParticleThisBin] = 0;

                for (int calcForceindexJ = 0; calcForceindexJ < BinMembers.size(); calcForceindexJ++ )
                {
                    //printf("Interaction!\n");
                    //apply_force_SOA( particlesSOA,ParticleThisBin, BinMembers[calcForceindexJ], &dmin, &davg, &navg);
                    apply_force( particles[ParticleThisBin], particles[ BinMembers[calcForceindexJ] ], &dmin, &davg, &navg);
                    //apply_force_SOA( particlesSOA,ParticleThisBin, BinMembers[calcForceindexJ], &dmin, &davg, &navg);
                 }
            }

            BinMembers.clear();

            // printf("##########################################################################################\n");
            // printf("##########################################################################################\n");
            // printf("##########################################################################################\n");
            // printf("\n");

        //}// Bins[BinNum].size() > 0

    } // end of bin calcs

    BinsWithParticles.clear();

        //
        //  move particles
        //
        for( int i = 0; i < n; i++ )
            move( particles[i] );

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
              save( fsave, n, particles );
        }


        #ifdef VIZ
            sprintf( outname, "out/fout-%05d.txt", step );

            FILE *fout = fopen( outname, "w");
            save( fout, n, particles);
            fclose( fout );
        #endif


    }
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
    free( particles );
    if( fsave )
        fclose( fsave );

    return 0;
}
