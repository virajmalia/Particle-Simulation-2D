#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include <vector>

#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)

//
//  benchmarking program
//
int main( int argc, char **argv )
{    

    // debugging output
    char *outname = (char*) malloc(32* sizeof(char));

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

    printf("Num of bins: %d size is: %f BinSize: %f \n", getbinNumber(), getSize(), getBinSize());

    // particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    // init_particles( n, particles );

    //printf("Test -1\n");
    particle_SOA_t *particlesSOA =(particle_SOA_t*)malloc( sizeof(particle_SOA_t));
    particlesSOA->x = (double*)malloc( n* sizeof(double));
    particlesSOA->y = (double*)malloc( n* sizeof(double)); 
    particlesSOA->vx = (double*)malloc( n* sizeof(double));
    particlesSOA->vy = (double*)malloc( n* sizeof(double));
    particlesSOA->ax = (double*)malloc( n* sizeof(double));
    particlesSOA->ay = (double*)malloc( n* sizeof(double));

    /// make bins 

    //printf("Test 0\n");
    init_particles_SOA( n, particlesSOA );

    //vector <int> * v 

    // //matrix of bins 
    
    int NumofBinsEachSide = getbinNumber();
    int NumofBins = NumofBinsEachSide*NumofBinsEachSide;
    // //std::vector<int> * Bins =( std::vector<int> *)malloc( (NumofBins^2) * sizeof( std::vector<int> ));
    // //std::vector<std:vector<int>> Bins; 
    // printf("NumofBinsEachSide %d \n", NumofBinsEachSide);
    // printf("Num of bins %d \n",NumofBins );

    /// Make the vector of vectors. The bins are vectors
    std::vector<std::vector<int> > Bins(NumofBins, std::vector<int>(0));

   // printf("Vector is %d ",Bins.size() );

    //printf("Test 3\n");

    //printf("BinSize is: %d \n ", Bins[1001].size());

	// Custom initializations to reduce arithmetic in apply_force()
	//double register mass_inv = 1/mass;
	//double register cutoff_sq = cutoff * cutoff;
	//double register min_r_sq = min_r * min_r;
    
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

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // clear out the bins since the particles have moved
        for(int clear = 0; clear < NumofBins; clear++ )
        {
            Bins[clear].clear();
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // re populate bins after the move update. 
        for(int particle = 0; particle < n; ++particle)
        {
              double binsize = getBinSize();
              //printf("Test4\n");
        //     // get the bin index
               int BinX = (int)(particlesSOA->x[particle]/binsize);
               int BinY = (int)(particlesSOA->y[particle]/binsize);
            
               //printf("Adding particle\n");
               int BinNum = BinX + NumofBinsEachSide*BinY;
               //printf("Particle added to Bin %d", BinNum);
               Bins[BinNum].push_back(particle);
               //printf("There are %d particles in this bin\n",Bins[BinNum].size());
        //     //addParticleToBin(n,BinX,BinY);
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        // store the particle indexs from each surrounding bin.
        std::vector<int> BinMembers; 

        // for each bin apply forces fromparticles within the cutoff distance.
        for(int BinIndex = 0; BinIndex < NumofBins; BinIndex++ )
        {
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

            //printf("North %d NorthEast %d NorthWest %d, East %d, West %d, South %d, SouthEast %d, SouthWest %d",North,NorthEast,NorthWest,East, West,South, SouthEast,SouthWest  );

            // Note: Left = West 
            // Right = East 

            /// NOW! Take all the particles in each of the bins and apply forces across all 9 bins.
            // We are always going to apply forces within a bin we are intrested in
            //BinMembers.insert(BinMembers.end(),Bins[BinIndex].begin(),Bins[BinIndex].end());

            // we are not at an edge or a corner. -- Most common case
            if( (Left | Right | Top | Bottom) == false)
            {

                //printf("ALL %d ", BinIndex );
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
                //printf("Testlast");
            }
            //printf("After");
            // Top Row 
            else if( Top )
            {
                // most common case for the top row -- Not in a corner. 
                if( (Left| Right) == false)
                {
                    //printf("Top Row %d ", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
                }
                else if( (not Left) & Right) // Yes this would be called a corner case!!!
                {
                    //printf("Top Row Right %d ", BinIndex );
                    // Right == East
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                }
                else // left corner 
                {
                    //printf("Top Row Left %d ", BinIndex );
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
                    //printf("Bottom Row %d ", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                    BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                }
                else if( (not Left) & Right) // Yes this would be called a corner case!!!
                {
                    // Right == East
                    //printf("Bottom Row Right %d ", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                    BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                    BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                }
                else // left corner 
                {
                    //printf("Bottom Row Left %d ", BinIndex );
                    BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                    BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                    BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());

                }

            }
            else if(Left)  // not in a corner on the left side
            {
                //printf("Left %d ", BinIndex );
                BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                BinMembers.insert(BinMembers.end(),Bins[NorthEast].begin(),Bins[NorthEast].end());
                BinMembers.insert(BinMembers.end(),Bins[East].begin(),Bins[East].end());
                BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
                BinMembers.insert(BinMembers.end(),Bins[SouthEast].begin(),Bins[SouthEast].end());
            }
            else // must be the right side 
            {
                //printf("Right %d ", BinIndex );
                BinMembers.insert(BinMembers.end(),Bins[NorthWest].begin(),Bins[NorthWest].end());
                BinMembers.insert(BinMembers.end(),Bins[North].begin(),Bins[North].end());
                BinMembers.insert(BinMembers.end(),Bins[West].begin(),Bins[West].end());
                BinMembers.insert(BinMembers.end(),Bins[SouthWest].begin(),Bins[SouthWest].end());
                BinMembers.insert(BinMembers.end(),Bins[South].begin(),Bins[South].end());
            }




            for(int calcForceindexI = 0; calcForceindexI < Bins[BinIndex].size(); calcForceindexI++ )
            {
                // apply forces 
                int ParticleThisBin = Bins[BinIndex][calcForceindexI];


                for (int calcForceindexJ = 0; calcForceindexJ < BinMembers.size(); calcForceindexJ++ )
                {
                      printf("Interaction");
                    apply_force_SOA( particlesSOA,ParticleThisBin, BinMembers[calcForceindexJ], &dmin, &davg, &navg);

                 }
            }

            BinMembers.clear();

        } // end of bin calcs

        //move particles
        
        for( int i = 0; i < n; i++ ) 
        {
            //move( particles[i] );  
            move_SOA( particlesSOA,i);
        }


        //
        //  compute forces
        //
     //    for( int i = 0; i < n; i++ )
     //    {
     //        //particles[i].ax = particles[i].ay = 0;
     //        //particlesSOA->ax[i] = 0;
     //        //particlesSOA->ay[i] = 0;


     //        for (int j = i+1; j < n; j++ )
     //        {

     //            loopcount++;
     //        //particlesSOA->ax[j] = 0;
     //        //particlesSOA->ay[j] = 0;

     //            apply_force_SOA( particlesSOA,i, j, &dmin, &davg, &navg);


		   //    //apply_force( particles[i], particles[j],&dmin,&davg,&navg);
    	//      }
     //         move_SOA( particlesSOA,i);
    	// }
 
 
        //

        sprintf( outname, "out/fout-%05d.txt", step );

        FILE *fout = fopen( outname, "w");
        save_SOA( fout, n, particlesSOA);
        fclose( fout );


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

    printf( "dmin = %f,davg = %f,navg = %f, loops = %d, applyforcescount = %d \n", dmin, davg, navg, loopcount, count );

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
    //free( particles );
    
    
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
