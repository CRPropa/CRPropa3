#ifndef MPC_SPATIAL_PARTITIONING_H
#define MPC_SPATIAL_PARTITIONING_H

#include "mpc/ModuleList.h"
#include "mpc/Candidate.h"

#include <vector>
#include <iostream>

namespace mpc {

struct Count {
	size_t count;
	Count();
//	operator size_t();
	void operator +=(size_t v);
};

struct Index {
	int x, y, z;
	Index();
	Index(Vector3 v);
	bool operator<(const Index &rhs) const;
};

class SpatialPartitioning: public Referenced {

public:

	typedef std::vector<ref_ptr<Candidate> > candidate_vector_t;

	SpatialPartitioning(ModuleList *moduleList);

	void run(candidate_vector_t &candidates, bool recursive);
	void run(Source *source, size_t count, bool recursive);

	void setPartitionOrigin(const Vector3 &origin);
	void setPartitionSize(double size);
	void setCurrentPartition(const Vector3 &offset);

	void setVerbose(bool verbose);
private:
	void run(Candidate *candidate, bool recursive,
			Loki::AssocVector<Index, Count> &partitions);
	Vector3 partitionOrigin, currentPartition, currentPartitionMargin;
	double partitionSize, partitionSizeMargin;
	bool verbose;
	ref_ptr<ModuleList> moduleList;
};

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::Index &idx);

#endif // MPC_SPATIAL_PARTITIONING_H
