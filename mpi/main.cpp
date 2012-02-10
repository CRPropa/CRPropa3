#include "Slave.h"
#include "Master.h"

#include <mpi.h>

#include <iostream>

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: mpc-mpi CONFIGFILE" << std::endl;
		return 1;
	}

	MPI_Init(&argc, &argv);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if (rank == 0) {
		Master master;
		master.run();
	} else {
		Slave slave;
		slave.load(argv[1]);
		slave.run();
	}

	MPI_Finalize();

	return 0;
}
