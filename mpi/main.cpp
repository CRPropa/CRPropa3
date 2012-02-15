#include "Slave.h"
#include "Master.h"

#include "kisslog.h"
#include "kissconvert.h"

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <sstream>

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "Usage: mpc-mpi <CONFIGFILE>" << std::endl;
		return 1;
	}

	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	std::string filename = "log-" + kiss::str(rank) + ".txt";
	std::ofstream log_file(filename.c_str());
	kiss::Log::setLogStream(log_file);

	KISS_LOG_INFO
		<< "Starting " << rank;

	if (rank == 0) {
		Master master;
		master.run();
	} else {
		Slave slave;
		slave.load(argv[1]);
		slave.run();
	}

	MPI_Finalize();

	KISS_LOG_INFO
		<< "Done.";

	return 0;
}
