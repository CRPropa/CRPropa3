#include "Master.h"

#include "kisslog.h"

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
	size_t total = 100;
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
