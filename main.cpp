#include "pngwriter.h"
#include "particle.h"
#include <sys/time.h>
#include <mpi.h>

using namespace std;

const int MASTER_NODE = 0;
const int TAG_PARTICLES = 45;
const int TAG_KILL = 101;
const int MAX_PARTICLE_SIZE = 8;
const int MAX_MASS = 655355555;
const int VELOCITY_SCALE = 1000;

const float G = 6.674 * pow(10.0,-11);

particle* generateParticles( int imageSize , int num );

void updatePositions( particle* objects , int numParticles );

void calculatePhysics( particle* objects , int numParticles, int index);

void sequential( particle* p , int imageSize , int numParticles, int timesteps , float maxMass, pngwriter* png );

int getParticleColor( particle* p , int numParticles );

float getMaxMass( particle* objects , int numParticles );

void drawParticles( particle* objects , int numParticles , pngwriter* png , float maxMass );

void distributeWork( particle* p , int numProcesses, int numParticles );

void gatherWork( particle* p , int numProcesses,  int numParticles);

void killSlaves( int numProcesses );

void printParticles( particle* p , int numparticles , int id );

unsigned long long int sequentialTime = 0;

unsigned long long int parallelTime = 0;

MPI_Datatype NBODY_PARTICLE_TYPE;

int main( int argc, char *argv[] ) {

	int id;

	int cores;

	MPI_Init( &argc , &argv );

	MPI_Comm_size(MPI_COMM_WORLD,&cores);

	MPI_Comm_rank(MPI_COMM_WORLD,&id);

	MPI_Type_contiguous(6, MPI_FLOAT, &NBODY_PARTICLE_TYPE);

	MPI_Type_commit(&NBODY_PARTICLE_TYPE);

	int imageSize = 500;

	int timeSteps = 300;

	int numParticles = 80;

	if( numParticles % (cores-1) != 0 )

		exit(1);

	//make sure this divides evenly
	int particlesPerNode = numParticles / ( cores - 1 );

	if( id == MASTER_NODE ) {

		particle* p = generateParticles( imageSize , 20 );

		pngwriter* png = new pngwriter( imageSize , imageSize , 0 , "" );

		float maxMass = getMaxMass( p , numParticles );

		//sequential( p , imageSize , numParticles , timeSteps , maxMass , png );

		//send it to all the nodes
		for( int i = 0 ; i < timeSteps ; i++ ) {

			distributeWork( p , cores , particlesPerNode );

			gatherWork( p , cores , particlesPerNode );

			updatePositions( p , numParticles );

			png->pngwriter_rename( i );

			drawParticles( p , numParticles , png , maxMass );

			}

		killSlaves( cores );

		exit(0);

		//system("convert -delay 2 -loop 0 *.png par_anim.gif && rm *.png");

		}

	else {

		particle* p = new particle[ numParticles ];

		MPI_Status stat;

		while( true ) {

			MPI_Recv( p , numParticles , NBODY_PARTICLE_TYPE , MASTER_NODE , MPI_ANY_TAG , MPI_COMM_WORLD, &stat );

			//if( id == 1 )

				//printParticles(p , particlesPerNode , id);

			if( stat.MPI_TAG == TAG_KILL )

				exit(0);

			for( int i = 0 ; i < particlesPerNode ; i++ )

				calculatePhysics( p , particlesPerNode , (id - 1) * particlesPerNode + i );

			MPI_Send( p , particlesPerNode , NBODY_PARTICLE_TYPE , MASTER_NODE , TAG_PARTICLES , MPI_COMM_WORLD );

			}

		}
	//cout << "proc " << id << "exiting";
	//MPI_Finalize();

	return 0;

	}

particle* generateParticles( int imageSize , int num ) {

	particle* p = new particle[num];

	srand( 548268452 );

	for( int i = 0 ; i < num ; i++ ) {

		p[i].id = i;

		p[i].x = rand() % imageSize;

		p[i].y = rand() % imageSize;

		p[i].vx = 0.0;

		p[i].vy = 0.0;

		p[i].m = rand() % MAX_MASS;


		}

	return p;

	}

void sequential( particle* p , int imageSize , int numParticles, int timesteps , float maxMass, pngwriter* png ) {

	for( int i = 0 ; i < timesteps ; i++ ) {

		png->pngwriter_rename( i );

		for( int j = 0 ; j < numParticles ; j++ ) {

			calculatePhysics( p , numParticles , j );

			}

		updatePositions( p , numParticles );

		drawParticles( p , numParticles , png , maxMass );

		}

	system("convert -delay 2 -loop 0 *.png animation.gif && rm *.png");

	}

void calculatePhysics( particle* objects , int numParticles, int index ) {

	for( int i = 0 ; i < numParticles ; i++ ) {

		if( i == index )

			continue;

		float r = sqrt( pow( objects[i].x - objects[index].x , 2 ) + pow( objects[i].y - objects[index].y , 2 ) ) + 0.1;

		float g = (G * objects[index].m * objects[i].m) / pow( r , 2 );

		float Fx = g * ( (objects[i].x - objects[index].x) / r );

		float Fy = g * ( (objects[i].y - objects[index].y) / r );

		objects[index].vx += (Fx / objects[index].m) * VELOCITY_SCALE;

		objects[index].vy += (Fy / objects[index].m) * VELOCITY_SCALE;

		}

	}

void updatePositions( particle* objects , int numParticles ) {

	for( int i = 0 ; i < numParticles ; i++ ) {

		objects[i].x += objects[i].vx;

		objects[i].y += objects[i].vy;

		}

	}

float getMaxMass( particle* objects , int numParticles ) {

	float maxMass = 0;

	for( int i = 0 ; i < numParticles ; i++ )

		maxMass = ( objects[i].m > maxMass ) ? objects[i].m : maxMass;

	return maxMass;

	}

void drawParticles( particle* objects , int numParticles , pngwriter* png , float maxMass ) {

	png->clear();

	for( int i = 0 ; i < numParticles ; i++ ) {

		int radius = round( (objects[i].m / maxMass) * MAX_PARTICLE_SIZE );

		int x = round( objects[i].x );

		int y = round( objects[i].y );

		int color = getParticleColor( &objects[i] , numParticles );

		png->filledcircle( x , y , radius , color , color , color );

		}

	png->write_png();

	}

int getParticleColor( particle* p , int numParticles ) {

	return round( 1000.0 + (p->id / numParticles * 1.0) * 64535.0 );

	}

/*long long int getTimeInMicroseconds() {

	struct timeval tv;
	struct timezone tz;
	struct tm *tm;

	gettimeofday(&tv, &tz);
	tm=localtime(&tv.tv_sec);

	return tm->tm_hour * 60 * 60 * 1000000 + tm->tm_min * 60 * 1000000 + tm->tm_sec * 1000000 + tv.tv_usec;

	}
*/

void killSlaves( int numProcesses ) {

	for( int i = 0 ; i < numProcesses ; i++ )

		MPI_Send( 0 , 0 , MPI_INT, i , TAG_KILL , MPI_COMM_WORLD );

	}

void distributeWork( particle* p , int numProcesses, int numParticles ) {

	for( int i = 1 ; i < numProcesses ; i++ )

		MPI_Send( p , numParticles , NBODY_PARTICLE_TYPE, i , TAG_PARTICLES , MPI_COMM_WORLD );

	}

void gatherWork( particle* p , int numProcesses,  int numParticles) {

	particle* buf = new particle[numParticles];

	MPI_Status stat;

	int particlesPerNode = numParticles / (numProcesses - 1);

	for( int i = 1 ; i < numProcesses ; i++ ) {

		MPI_Recv( buf , numParticles , NBODY_PARTICLE_TYPE , MPI_ANY_SOURCE , MPI_ANY_TAG , MPI_COMM_WORLD, &stat );

		for( int j = 0 ; j < particlesPerNode ; j++ )

			p[ (stat.MPI_SOURCE - 1) * particlesPerNode + j ] = buf[ (stat.MPI_SOURCE - 1) * particlesPerNode + j ];

		}

	}
void printParticles( particle* p , int numparticles , int id ) {

	cout << "Particles for process " << id << endl;

	for( int i = 0 ; i < numparticles ; i++ ) {

		cout << "id: " << p[i].id << endl << "x: " << p[i].x << endl << "y: " << p[i].y << endl << endl;

	}
}
