#include <iostream>
#include <mpi.h>
#include "pngwriter.h"
#include <vector>
#include <sys/time.h>
#include "particle.h"
#include <time.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int MASTER_NODE = 0;
int TAG_PARTICLES = 45;
int TAG_KILL = 101;

particle* generateParticles( int imageSize , int num );

float* particlesToFloats( particle* p , int num );

int main( int argc, char *argv[] ) {

	int id;

	int cores;

	MPI_Status stat;

	MPI_Init( &argc , &argv );

	MPI_Comm_size(MPI_COMM_WORLD,&cores);

	MPI_Comm_rank(MPI_COMM_WORLD,&id);

	int imageSize = 500;

	int timeSteps = 300;

	int numParticles = 80;

	//make sure this divides evenly
	int particlesPerNode = numParticles / cores * 1.0;

	if( id == MASTER_NODE ) {

		//generate some random particles
		particle* p = generateParticles( imageSize , 20 );

		//send it to all the nodes
		for( int i = 0 ; i < timeSteps ; i++ ) {

			for( int j = 0 ; j < cores ; j++ ) {

				MPI_Send( particlesToFloats(p,numParticles) , numParticles , MPI_FLOAT , j , TAG_PARTICLES , MPI_COMM_WORLD );

				}

			float* f = new float[ 6 * particlesPerNode ];

			for( int j = 0 ; j < cores ; j++ ) {

				MPI_Recv( f , 6 * particlesPerNode , MPI_FLOAT , j , TAG_PARTICLES , MPI_COMM_WORLD , &stat );

				}

			pngwriter* png = new pngwriter( imageSize , imageSize , 0 , "output.png" );

			}

		}

	else {



		}

	MPI_Finalize();

	return 0;

	}

particle* generateParticles( int imageSize , int num ) {

	particle* p = new particle[num];



	for( int i = 0 ; i < num ; i++ ) {
		srand( i );
		p[i].id = i;

		p[i].x = rand() % imageSize;

		p[i].y = rand() % imageSize;

		p[i].vx = 0.0;

		p[i].vy = 0.0;

		p[i].m = rand() % 10;


		}

	return p;

	}



float* particlesToFloats( particle* p , int num ) {

		float* f = new float[ 6 * num ];

		for( int i = 0 ; i < num ; i++ ) {

			f[i * 6] = p[i].id;
			f[i * 6 + 1] = p[i].x;
			f[i * 6 + 2] = p[i].y;
			f[i * 6 + 3] = p[i].vx;
			f[i * 6 + 4] = p[i].vy;
			f[i * 6 + 5] = p[i].m;

			}

		return f;

		}



