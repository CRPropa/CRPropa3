#include "Slave.h"

#include "mpc/XMLImport.h"

Slave::Slave() {
}

Slave::~Slave() {
}

void Slave::load(const std::string &filename) {
	mpc::XMLImport import(&chain);
	import.import(filename);
}

void Slave::run() {
	int work;
	double result;
	MPI_Status status;

	while (1) {

		/* Receive a message from the master */

		MPI_Recv(&work, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		/* Check the tag of the received message. */

		if (status.MPI_TAG == TAG_STOP) {
			return;
		} else if (status.MPI_TAG == TAG_WORK) {

			/* Do the work */

			result = work * work;

			/* Send the result back */

			MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		} else {
			printf("unknown tag");
		}
	}
}

