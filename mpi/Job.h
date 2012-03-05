#ifndef MPC_JOB_H_
#define MPC_JOB_H_

#include "common.h"

#include "kiss/io.h"

#include <fstream>

class Job {
	std::ifstream in_stream;
	kiss::StreamInput in;

	std::ofstream out_stream;
	kiss::StreamOutput out;

public:
	Job(job_t job);
	size_t processSome();
};

#endif /* MPC_JOB_H_ */
