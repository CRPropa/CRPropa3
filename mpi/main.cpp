#include "Slave.h"
#include "Master.h"

#include <kiss/logger.h>
#include <kiss/convert.h>

#include <mpi.h>

#include <iostream>
#include <fstream>
#include <sstream>

using namespace kiss;
using namespace std;

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "Usage: mpc-mpi <CONFIGFILE>" << endl;
		return 1;
	}

	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	string filename = "log-" + str(rank) + ".txt";
	ofstream log_file(filename.c_str());
	Logger::setLogStream(log_file);

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
