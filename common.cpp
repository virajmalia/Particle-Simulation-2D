
#include "common.h"


double size;

Local_Space_t LocalSpaceInfo; 


#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define cutoffSQ (cutoff*cutoff)
#define INVCutoff (1/cutoff)
#define min_r   (cutoff/100)
#define min_r_SQ (min_r*min_r)
#define dt      0.0005


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
    gettimeofday( &end, NULL );//
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//
void set_size( int n )
{
    size = sqrt( density * n );
}

int getNumberofBins( double size)
{
    // need to round up for partial bins
    return (int)ceil( size/cutoff );
}

double getSize()
{
    return size;
}

double getBinSize()
{
    return cutoff;
}

void init_particles( int n, particle_t *p )
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
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;
    }
    free( shuffle );
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

//
//  interact two particles
//
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{

    // this is a little stupid  since it applies the force in only one direction
    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;
    double r2 = dx * dx + dy * dy;
    if( r2 > cutoffSQ )
        return;

    double r = sqrt( r2 );

	if (r2 != 0)
    {
	   if (r2/(cutoffSQ) < *dmin * (*dmin))
       {
	      *dmin = r/cutoff;
       }
           (*davg) += r/cutoff;
           (*navg) ++;
    }

    r2 = fmax( r2, min_r_SQ );
    r = sqrt( r2 );

    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}


//
//  interact two particles
//
void apply_force_SOA( particle_SOA_t &p, int I, int J, double *dmin, double *davg, int *navg)
{

    // this is a little stupid  since it applies the force in only one direction
    // double dx = neighbor.x - particle.x;
    // double dy = neighbor.y - particle.y;

    double dx = p.x[J] - p.x[I];
    double dy = p.y[J] - p.y[I];

    double r2 = dx * dx + dy * dy;

    if( r2 > cutoffSQ )
    {
        //printf(" R2 = %f, dx = %f, dy = %f, I= %d, J= %d \n", r2, dx, dy, I, J);
        return;
    }

    double r = sqrt( r2 );
    double dist = r/cutoff;

    if (r2 != 0)
    {
       if (dist < *dmin)
       {
          *dmin = dist;
       }
        (*davg) += dist;
        // since we are doing pairs of particles at the same time.
        (*navg) ++;
    }

    r2 = fmax( r2, min_r_SQ);
    r = sqrt( r2 );
    //
    //  very simple short-range repulsive force
    // but do both at the same time!!!!
    double coef = ( 1 - cutoff / r ) / r2 / mass;

    double accelX = coef * dx;
    double accelY = coef * dy;

    p.ax[I] += accelX;
    p.ay[I] += accelY;
    // p->ax[J] -= accelX;  // force applied in opposite direction
    // p->ay[J] -= accelY;  // force applied in opposite direction
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //

    //printf("Particle Vel is %f", p.vy);
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt; // this line will crash the MPI program on odd values of n in of 2^n processors
    // double temp = p.vy * dt;
    // temp = temp + p.y;
    // double temp2 = temp; 
    // p.y = temp2;
    //printf("Particle Y: %f Vel is %f",p.y , p.vy);
    // // once the force is applied awesome the accel is zero.
    p.ax = 0;
    p.ay = 0;


    // //
    // //  bounce from walls
    // //
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

void move_SOA( particle_SOA_t &p,int I)
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx[I] += p.ax[I] * dt;
    p.vy[I] += p.ay[I] * dt;
    p.x[I]  += p.vx[I] * dt;
    p.y[I]  += p.vy[I] * dt;




    // once the force is applied awesome the accel is zero.
    p.ax[I] = 0;
    p.ay[I] = 0;

    //
    //  bounce from walls
    //
    while( p.x[I] < 0 || p.x[I] > size )
    {
        p.x[I]  = p.x[I] < 0 ? -p.x[I] : 2*size-p.x[I];
        p.vx[I] = -p.vx[I];
    }
    while( p.y[I] < 0 || p.y[I] > size )
    {
        p.y[I]  = p.y[I] < 0 ? -p.y[I] : 2*size-p.y[I];
        p.vy[I] = -p.vy[I];
    }

}

Bin_Location_t GetBinLocation(const int BinIndex, const int NumofBinsEachSide,const int NumofBins )
{ // assumes square geomerty!

    Bin_Location_t Temp; 

    Temp.Left = ((BinIndex%NumofBinsEachSide) == 0) ? true : false;
    Temp.Right = ((BinIndex%NumofBinsEachSide) == (NumofBinsEachSide-1) ) ? true : false;
    Temp.Top =  ((BinIndex < NumofBinsEachSide) )? true : false;
    Temp.Bottom = ((BinIndex > (NumofBins - NumofBinsEachSide - 1) ) )? true : false;

    return Temp; 
}

Neighbor_Indexes_t GetNeighborBinIndexes(const int BinIndex, const int NumofBinsEachSide)
{
    Neighbor_Indexes_t Temp;

        Temp.North = BinIndex - NumofBinsEachSide;
        Temp.NorthEast = Temp.North + 1;
        Temp.NorthWest = Temp.North -1;
        Temp.East = BinIndex + 1;
        Temp.West = BinIndex -1;
        Temp.South = BinIndex + NumofBinsEachSide;
        Temp.SouthEast = Temp.South +1;
        Temp.SouthWest = Temp.South -1;

    return Temp; 
};



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////// MPI /////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Neighbor_Indexes_t GetGhostBinLocations(const int BinIndex)
{// assumes only the top or bottom row!
    Neighbor_Indexes_t Temp;
        // use for top host row
        Temp.North = BinIndex;
        Temp.NorthEast = Temp.North + 1;
        Temp.NorthWest = Temp.North -1;
        Temp.East = -1;
        Temp.West = -1;
        // use for bottom  ghost row 
        Temp.South = BinIndex;
        Temp.SouthEast = Temp.South +1;
        Temp.SouthWest = Temp.South -1;

    return Temp; 
}


/// warning !! dependent on row divion!. Blocking to come later 
int getRowsPerProc(int NumberOfBinsperSide, int NumberofProcessors)
{  // we are just dividing the space up by row 
    // need to convert to a float first otherwise the value will be rounded down. 
    return floor((float)NumberOfBinsperSide/NumberofProcessors);
}


int getNumberofBinsLocal(int GlobalNumberOfBinsEachSide, int NumofBins, int rank, int NumberofProcessors)
{
    // this is for Rank 0 to NumberofProcessors -1. The last Proc gets the remainder 
    int RowEachProc = getRowsPerProc(GlobalNumberOfBinsEachSide, NumberofProcessors);

    if((rank < (NumberofProcessors -1)) || (rank == 0))
    {
        return GlobalNumberOfBinsEachSide * RowEachProc;
    }
    else
    {   // the last proc takes care of the remainder 

        return (NumofBins) - ( (NumberofProcessors -1)  * GlobalNumberOfBinsEachSide * RowEachProc );
        //return (GlobalNumberOfBinsEachSide % RowEachProc) * GlobalNumberOfBinsEachSide;
    }


}

void set_local_space(double size, int rank, int GlobalNumberOfBinsEachSide, int NumberofProcessors)
{

    int RowsPerProc = getRowsPerProc(GlobalNumberOfBinsEachSide, NumberofProcessors);
    double BinSize = getBinSize();

    if((rank < NumberofProcessors -1) || (rank == 0))
    {   
        LocalSpaceInfo.Floor = ( (rank + 1) * ( BinSize *  RowsPerProc) );
    }
    else // we are the remainder 
    {
        LocalSpaceInfo.Floor = size; 
    }

    LocalSpaceInfo.localSizeX = size; // always the same in the X
    LocalSpaceInfo.Ceiling = ( rank * (BinSize *  RowsPerProc) );
    LocalSpaceInfo.localSizeY = LocalSpaceInfo.Floor - LocalSpaceInfo.Ceiling;
}

double getLocalYSize()
{
    return LocalSpaceInfo.localSizeY; 
}

double getLocalXSize()
{
    return LocalSpaceInfo.localSizeX;
}

// binning must have already occured 
std::vector<particle_t> getGhostParticlesTop(const int rank, const int LocalNumofBinsEachSide, const int NumberofProcessors, const std::vector< std::vector<int> > & LocalBins, const std::vector <particle_t> & localParticleVec)
{ // need to finish 
    // lowest number is top left 
    // highest number in bottom right. 
    // int topleft = std::min_element(BinsByProc[rank]);
    // int bottonright = std::max_element(BinsByProc[rank]);
    std::vector<particle_t> GhostParticlesTop;

    if( (rank > 0 ) )  // if rank 0 and numprocessors is 1 we don't have any peers
    {
        for(int BinNum = 0; BinNum < LocalBins.size(); BinNum++) 
        {
            //printf("Local size is : %d, BinsEachSide is %d \n", LocalBins.size(),LocalNumofBinsEachSide );
            for(int Ghostparticle = 0; Ghostparticle < LocalBins[BinNum].size();Ghostparticle++)
            {
                 GhostParticlesTop.push_back(localParticleVec[ LocalBins[BinNum][Ghostparticle] ]);
            }  
        }

    }
        
    return GhostParticlesTop;

}


// binning must have already occured 
std::vector<particle_t> getGhostParticlesBottom(const int rank, const int LocalNumofBinsEachSide, const int NumberoflocalBins, const int NumberofProcessors, const std::vector< std::vector<int> > & LocalBins, const std::vector <particle_t> & localParticleVec)
{
    std::vector<particle_t> GhostParticlesBottom;

    if( (rank < (NumberofProcessors - 1) ))  // if rank 0 and numprocessors is 1 we don't have any peers  // if rank 0 and numprocessors is 1 we don't have any peers
    {
        for(int BinNum = (NumberoflocalBins - LocalNumofBinsEachSide); BinNum < NumberoflocalBins; BinNum++) 
        {
            for(int Ghostparticle = 0; Ghostparticle < LocalBins[BinNum].size();Ghostparticle++)
            {
                 GhostParticlesBottom.push_back(localParticleVec[ LocalBins[BinNum][Ghostparticle] ]);
            }  
        }

    }

    // PROC OF RANK = NUMBEROFPROC -1 DOES NOT HAVE A BOTTOM PEER
        
    return GhostParticlesBottom;
}


std::vector<int> getBoarderPeers(int rank, int NumberofProcessors)
{ ///FIXME!
    std::vector<int> Peers;

    if(NumberofProcessors > 1) // if it's not we are the only process 
    {
        if(rank == 0)
        { // we can only have a 
            Peers.push_back(rank+1);
        }
        else if((rank > 0) && (rank <(NumberofProcessors -1) ) )
        {
            Peers.push_back(rank-1);
            Peers.push_back(rank+1);
        }
        else // rank is equal to the numberofprc -1 
        {
            Peers.push_back(rank-1); 
        }
    }

    return Peers;
}

int MaplocalBinToGlobalBin(int rank, int localbinNumber, int NumberOfBinsperSide,int NumberofProcessors)
{
    int RowsPerProc = getRowsPerProc(NumberOfBinsperSide,NumberofProcessors); 

    return ((rank * NumberOfBinsperSide * RowsPerProc) + localbinNumber); // assuming localbins are the same lenght as the global bins. 
}

int MapGlobalBinToLocalBin(int rank, int GlobalBinNumber, int NumberOfBinsperSide,int NumberofProcessors)
{
    int RowsPerProc = getRowsPerProc(NumberOfBinsperSide,NumberofProcessors); 

    return (GlobalBinNumber - (rank * NumberOfBinsperSide * RowsPerProc)); // assuming localbins are the same lenght as the global bins. 
}




std::vector< std::vector<int> > PopulateProcBinVector(int NumberOfBinsperSide,int NumberofProcessors)
{ 

    std::vector< std::vector<int> > MapOfBinsToProcs(NumberofProcessors, std::vector<int>(0));

    // for(int ProcId = 0; ProcId < NumberofProcessors; ProcId)
    // {
    //     MapOfBinsToProcs[ProcId].push_back();
    // }

    //// do this later

    // int LogOfProcs = std::log2(NumberofProcessors);
    // int BinsPerProc = 0;

    // if(LogOfProcs % 2 == 0)
    // {
    //     //BinsPerProc = (NumberOfBins*NumberOfBinsperSide)/NumberofProcessors; 

    //     // terrible n^3 bu it's always a small number and only run once. 
    //     for(int ProcNum = 0; ProcNum < NumberofProcessors; NumberofProcessors++)
    //     {
    //         int offset = ProcNum * LogOfProcs; 

    //         //MapOfBinsToProcs.push_back(ProcNum);

    //         for (int col = 0; col < LogOfProcs; col++)
    //         {
    //             for(int Row = 0; Row < LogOfProcs; LogOfProcs++)
    //             {
    //                 int BlockNum = offset + Row + (col* size);
    //                 MapOfBinsToProcs[ProcNum].push_back(BlockNum); 
    //             } 
    //         }

    //     } //  for(int ProcNum = 0; ProcNum < NumberofProcessors; NumberofProcessors++;)

    // }
    // else // we have an odd power of two 
    // {
    //     // divide bins into 2 sets 
    //     int topsection = getRowsPerProc(NumberOfBinsperSide,NumberofProcessors); // rounds down 
    //     int bottomsection = NumberOfBinsperSide - topsection; 

    //     // do top 
    //     // for(int ProcNum = 0; ProcNum < NumberofProcessors/2; NumberofProcessors++)
    //     // {
    //     //     int offset = ProcNum * LogOfProcs; 

    //     //     //MapOfBinsToProcs.push_back(ProcNum);

    //     //     for (int col = 0; col < LogOfProcs; col++)
    //     //     {
    //     //         for(int Row = 0; Row < LogOfProcs; LogOfProcs++)
    //     //         {
    //     //             int BlockNum = offset + Row + (col* size);
    //     //             MapOfBinsToProcs[ProcNum].push_back(BlockNum); 
    //     //         } 
    //     //     }

    //     // } //  for(int ProcNum = 0; ProcNum < NumberofProcessors; NumberofProcessors++;)

    //     // do bottom 


    //     //divide each set of bins by LogOfProcs -1
    //     // same algorithm as above. 
 
    // }

    return MapOfBinsToProcs;
}

// bins will not not move. Once set, their location is static 
int MapBinToProc(const int GlobalBin, const int NumberofProcessors, const int NumberOfBinsperSide)
{ // untested this needs to be figured out

    int RowsPerProc = getRowsPerProc(NumberOfBinsperSide,NumberofProcessors); 

    int MostBins = (NumberofProcessors-1) * NumberOfBinsperSide * RowsPerProc;

    if(GlobalBin < MostBins)
    {
        return floor(GlobalBin / (NumberOfBinsperSide * RowsPerProc) ); // assuming localbins are the same size as the global bins. 
    }
    else
    { // the remainder
        return (NumberofProcessors -1);
    }
}

int MapParticleToBin(particle_t &particle, const int NumofBinsEachSide)
{
    // this will give you which bin the particle must be locatd in 
    double binsize = getBinSize();
              //printf("Test4\n");
              // get the bin index

   int BinX = (int)(particle.x/binsize);
   int BinY = (int)(particle.y/binsize);

   // int BinX = (int)(particlesSOA->x[particle]/binsize);
   // int BinY = (int)(particlesSOA->y[particle]/binsize);

   //printf("Adding particle\n");
   return  (BinX + NumofBinsEachSide*BinY);
}

int MapParticleToProc(particle_t &particle, const int NumofBinsEachSide, const  int NumberofProcessors )
{   // this will retun the processor to which a given particle belongs. 
    int BinNum = MapParticleToBin(particle,NumofBinsEachSide); // global map!!
    return MapBinToProc(BinNum,NumberofProcessors,NumofBinsEachSide);
}
// int getProcessorForBin(int Bin)
// {
//     return 
// }

//
//  Initialize the particle positions and velocities
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//  I/O routines
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void save( FILE *f, int n, particle_t *p )
{
    static bool first = true;
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

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
