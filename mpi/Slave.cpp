#include "Slave.h"

#include "mpc/XMLImport.h"
#include "mpc/Vector3.h"

#include "kissconvert.h"

#include <fstream>
#include <sstream>

using namespace mpc;
using namespace std;
using namespace kiss;

Slave::Slave() {
}

Slave::~Slave() {
}

void Slave::load(const string &filename) {
	XMLImport import(&chain);
	import.import(filename);
}

void Slave::run() {
	job_t job;
	double result;
	MPI_Status status;

	while (1) {

		/* Receive a message from the master */

		MPI_Recv(&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		/* Check the tag of the received message. */

		if (status.MPI_TAG == TAG_STOP) {
			return;
		} else if (status.MPI_TAG == TAG_WORK) {
			processJob(job);
			/* Send the result back */

			MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		} else {
			printf("unknown tag");
		}
	}
}

struct JobInfo {
	Vector3 partition_offset;
	double partition_size, partition_margin;
};

void Slave::processJob(job_t job) {
	string filename = "job_" + str(job) + ".dat";
	ifstream in(filename.c_str(), ios::binary);

//	JobInfo info;
//	in.read((char *) &info, sizeof(info));
//	chain.setSimulationPartition(info.partition_offset, info.partition_size,
//			info.partition_margin);

	Candidate candidate;

}

