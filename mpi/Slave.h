#ifndef SLAVE_H_
#define SLAVE_H_

#include "common.h"

#include "mpc/ModuleChain.h"

#include <mpi.h>

class Slave {
	bool running;
	mpc::ModuleChain chain;
	void processJob(job_t job);
	void acquireJob();
public:
	Slave();
	virtual ~Slave();

	void load(const std::string &filename);
	void run();
};

#endif /* SLAVE_H_ */
