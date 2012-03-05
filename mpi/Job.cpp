#include "Job.h"

#include "kiss/convert.h"

using namespace std;
using namespace kiss;

Job::Job(job_t job) :
		in(in_stream), out(out_stream) {
	string in_filename = "job_" + str(job) + ".dat";
	in_stream.open(in_filename.c_str(), ios::binary);

	string out_filename = "job_" + str(job) + "_out.dat";
	out_stream.open(out_filename.c_str(), ios::binary);
}

size_t processSome(size_t count) {
	size_t processed = 0;
	while (true) {
//		Candidate candidate;
//		if (!read(in, candidate))
//			break;
//		chain.process(&candidate);
//		write(out, candidate);
//		processed++;
//		for (size_t iSec = 0; iSec < secondaries.size(); iSec++) {
//			write(out, *secondaries[iSec]);
//			delete secondaries[iSec];
//		}
	}
	return processed;
}
