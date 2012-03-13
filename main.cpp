/*
 * main.cpp
 *
 *  Created on: Mar 12, 2012
 *      Author: njordan
 */



#include <iostream>
#include <mpi/mpi.h>
#include "pngwriter.h"
#include <vector>
#include <sys/time.h>

using namespace std;

int main( int argc, char *argv[] ) {

	int id;

	int cores;

	MPI_Status stat;

	MPI_Init( &argc , &argv );

	MPI_Comm_size(MPI_COMM_WORLD,&cores);

	MPI_Comm_rank(MPI_COMM_WORLD,&id);

	return 0;

	}
