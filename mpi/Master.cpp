#include "Master.h"

#include "mpc/IO.h"
#include "mpc/Candidate.h"

#include "kissLog.h"
#include "kissConvert.h"

#include <fstream>

using namespace std;
using namespace kiss;
using namespace mpc;

Master::Master() :
		_NextJob(0) {
	MPI_Comm_size(MPI_COMM_WORLD, &_RankCount);
	_Slaves.reserve(_RankCount);
}

Master::~Master() {
}

void Master::run() {
	sendFirstWorkItem();

	job_t work = getNextJob();
	size_t total = 1;
	while (work != 0) {

		unit_result_t result;
		MPI_Status status;

		receiveResult(result, status);

		sendJob(work, status.MPI_SOURCE);

		work = getNextJob();
		total--;
		if (total == 0)
			break;
	}

	receiveOutstandingResults();

	stopSlaves();
}

void Master::receiveResult(unit_result_t &result, MPI_Status &status) {
	MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
			MPI_COMM_WORLD, &status);
	KISS_LOG_DEBUG
		<< "receive result '" << result << "' from '" << status.MPI_SOURCE
				<< "'";
}

void Master::sendJob(job_t &job, int rank) {
	KISS_LOG_DEBUG
		<< "sent job '" << job << "' to '" << rank << "'";
	_Slaves[rank].job = job;
	Candidate candidate;
	ParticleState initial;
	initial.setId(getNucleusId(56, 26));
	initial.setEnergy(200 * EeV);
	initial.setPosition(Vector3(0, 0, 0));
	initial.setDirection(Vector3(1, 0, 0));

	candidate.current = initial;
	candidate.initial = initial;
	candidate.setNextStep(0.01 * Mpc);

	string out_filename = "job_" + str(job) + ".dat";
	ofstream out(out_filename.c_str(), ios::binary);
	write(out, candidate);

	MPI_Send(&job, 1, MPI_INT, rank, TAG_WORK, MPI_COMM_WORLD);
}

int Master::getNextJob() {
	job_t next = 0;
	if (_FreeJobs.size()) {
		next = _FreeJobs.back();
		_FreeJobs.pop_back();
	} else {
		_NextJob++;
		return _NextJob;
	}

	return next;
}

void Master::sendFirstWorkItem() {
	for (size_t rank = 1; rank < _RankCount; ++rank) {
		job_t work = getNextJob();
		sendJob(work, rank);
	}
}

void Master::receiveOutstandingResults() {
	for (size_t rank = 1; rank < _RankCount; ++rank) {
		unit_result_t result;
		MPI_Status status;
		MPI_Recv(&result, 1, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG,
				MPI_COMM_WORLD, &status);
	}

}
void Master::stopSlaves() {
	for (size_t rank = 1; rank < _RankCount; ++rank) {
		MPI_Send(0, 0, MPI_INT, rank, TAG_STOP, MPI_COMM_WORLD);
	}
}
