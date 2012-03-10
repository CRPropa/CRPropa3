#include "Slave.h"

#include "mpc/IO.h"
#include "mpc/XMLImport.h"
#include "mpc/ModuleChain.h"
#include "mpc/Vector3.h"
#include "mpc/Units.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/GlutDisplay.h"
#include "mpc/module/Output.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/magneticField/turbulentMagneticFieldGrid.h"

#include <kiss/convert.h>

#include <fstream>
#include <sstream>

using namespace mpc;
using namespace std;
using namespace kiss;

Slave::Slave() :
		running(true) {
}

Slave::~Slave() {
}

void Slave::load(const string &filename) {
//	XMLImport import(&chain);
//	import.import(filename);
	UniformMagneticField *field = new UniformMagneticField(
			Vector3(0., 0., 1e-20));
//		TurbulentMagneticFieldGrid *field = new TurbulentMagneticFieldGrid(Vector3(0, 0, 0), 64, 0.1 * Mpc, 1e-12, 1, 8, -11. / 3., 10);
	chain.add(25, new DeflectionCK(field));

	// interactions -------------------------------------------------------
	chain.add(30, new NuclearDecay());
	chain.add(31, new PhotoDisintegration());
	chain.add(32, new ElectronPairProduction(ElectronPairProduction::CMB));
	chain.add(33, new PhotoPionProduction(PhotoPionProduction::CMBIR));

	// break conditions ---------------------------------------------------
//		chain.add(new MinimumEnergy(5 * EeV), 50);
//		chain.add(new MaximumTrajectoryLength(100 * Mpc), 51);
	chain.add(
			52,
			new SphericalBoundary(Vector3(0, 0, 0) * Mpc, 20 * Mpc, 0.1 * Mpc,
					Candidate::Detected));
//		chain.add(new SmallObserverSphere(Vector3(0, 0, 0) * Mpc, 1 * Mpc), 53);

// output -------------------------------------------------------------
	chain.add(79, new ShellOutput());
//		chain.add(new TrajectoryOutput("trajectories.csv"), 80);
//		chain.add(new GlutDisplay(), 80);
	chain.add(100, new FlaggedOutput("final.txt", Candidate::Detected));
}

void Slave::acquireJob() {
	job_t job;
	double result;
	MPI_Status status;

	MPI_Send(0, 0, MPI_DATATYPE_NULL, 0, TAG_REQUEST_JOB, MPI_COMM_WORLD);
	MPI_Recv(&job, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

	if (status.MPI_TAG == TAG_STOP) {
		running = false;
	} else if (status.MPI_TAG == TAG_WORK) {
		processJob(job);
		MPI_Send(&result, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	} else {
		printf("unknown tag");
	}
}

void Slave::run() {
	while (running)
		acquireJob();
}

struct JobInfo {
	Vector3 partition_offset;
	double partition_size, partition_margin;
};

void Slave::processJob(job_t job) {

	// any local jobs available?
	// otherwise ask master/other ranks

	// collect particles
	// process job
	// store particles for each block in files
	// register job at master

	string in_filename = "job_" + str(job) + ".dat";
	ifstream in_stream(in_filename.c_str(), ios::binary);
	StreamInput in(in_stream);

	string out_filename = "job_" + str(job) + "_out.dat";
	ofstream out_stream(out_filename.c_str(), ios::binary);
	StreamOutput out(out_stream);

//	JobInfo info;
//	in.read((char *) &info, sizeof(info));
//	chain.setSimulationPartition(info.partition_offset, info.partition_size,
//			info.partition_margin);
	while (true) {
		Candidate candidate;
		if (!read(in, candidate))
			break;
		chain.process(&candidate);
		write(out, candidate);
//		for (size_t iSec = 0; iSec < secondaries.size(); iSec++) {
//			write(out, *secondaries[iSec]);
//			delete secondaries[iSec];
//		}
	}

}
