#ifndef MASTER_H_
#define MASTER_H_

#include "common.h"

#include <mpi.h>

#include <vector>

class Master {
	struct SlaveInfo {
		job_t job;
	};
	std::vector<job_t> _FreeJobs;
	std::vector<SlaveInfo> _Slaves;
	job_t _NextJob;
	int _RankCount;

	void sendFirstWorkItem();
	int getNextJob();
	void stopSlaves();
	void receiveOutstandingResults();
	void receiveResult(unit_result_t &result, MPI_Status &status);
	void sendJob(job_t &work, int rank);
public:
	Master();
	virtual ~Master();
	void run();
};

#endif /* MASTER_H_ */
