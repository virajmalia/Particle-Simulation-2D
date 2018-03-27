
#include "common.h"



double size;

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

int getbinNumber()
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

int getRowsPerProc(int NumberOfBinsperSide, int NumberofProcessors, int size)
{ 
    // need to convert to a float first otherwise the value will be rounded down. 
    return ceil((float)NumberOfBinsperSide/NumberofProcessors);
}

std::vector<int> getBoarderPeers(int rank)
{ ///FIXME!

    return std::vector<int> (rank);
}

std::vector< std::vector<int> > PopulateProcBinVector(int NumberOfBins,int NumberofProcessors, int size)
{ 

    std::vector< std::vector<int> > MapOfBinsToProcs(NumberofProcessors, std::vector<int>(0));

    int LogOfProcs = std::log2(NumberofProcessors);
    int BinsPerProc = 0;

    if(LogOfProcs % 2 == 0)
    {
        BinsPerProc = NumberOfBins/NumberofProcessors; 

        // terrible n^3 bu it's always a small number and only run once. 
        for(int ProcNum = 0; ProcNum < NumberofProcessors; NumberofProcessors++)
        {
            int offset = ProcNum * LogOfProcs; 

            //MapOfBinsToProcs.push_back(ProcNum);

            for (int col = 0; col < LogOfProcs; col++)
            {
                for(int Row = 0; Row < LogOfProcs; LogOfProcs++)
                {
                    int BlockNum = offset + Row + (col* size);
                    MapOfBinsToProcs[ProcNum].push_back(BlockNum); 
                } 
            }

        } //  for(int ProcNum = 0; ProcNum < NumberofProcessors; NumberofProcessors++;)

    }
    else // we have an odd power of two 
    {
        // divide bins into 2 sets 

        //divide each set of bins by LogOfProcs -1
        // same algorithm as above. 
 
    }

    return MapOfBinsToProcs;
}


std::vector<int > getGhostbins(int size, int rank, const std::vector< std::vector<int> > & BinsByProc)
{ // need to finish 


    // lowest number is top left 
    // highest number in bottom right. 
    // int topleft = std::min_element(BinsByProc[rank]);
    // int bottonright = std::max_element(BinsByProc[rank]);

    std::vector<int> ghostbins;

    return ghostbins;

}

// bins will not not move. Once set, their location is static 
int MapBinToProc(int Bin, int NumberofProcessors)
{ // untested this needs to be figured out
    int Processor = 0; 
    return Processor;

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
    int BinNum = MapParticleToBin(particle,NumofBinsEachSide);
    return MapBinToProc(BinNum,NumberofProcessors);
}
// int getProcessorForBin(int Bin)
// {
//     return 
// }

//
//  Initialize the particle positions and velocities


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
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    // once the force is applied awesome the accel is zero.
    p.ax = 0;
    p.ay = 0;

    //
    //  bounce from walls
    //
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
