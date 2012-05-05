#include "mpc/SpatialPartitioning.h"

namespace mpc {

Count::Count() :
		count(0) {

}

//operator Count::size_t () {
//	return count;
//}

void Count::operator +=(size_t v) {
	count += v;
}

Index::Index() :
		x(0), y(0), z(0) {

}
Index::Index(Vector3d v) :
		x(floor(v.x)), y(floor(v.y)), z(floor(v.z)) {

}

bool Index::operator<(const Index &rhs) const {
	if (x < rhs.x)
		return true;
	else if (x > rhs.x)
		return false;

	if (y < rhs.y)
		return true;
	else if (y > rhs.y)
		return false;

	if (z < rhs.z)
		return true;
	else
		return false;

}

SpatialPartitioning::SpatialPartitioning(ModuleList *moduleList) :
		partitionOrigin(0, 0, 0), partitionSize(1 * Mpc), verbose(false), moduleList(
				moduleList) {
}

bool findActivePosition(Candidate *candidate, bool recursive,
		Vector3d &position) {
	if (candidate->isActive()) {
		position = candidate->current.getPosition();
		return true;
	}

	if (recursive) {
		for (size_t i = 0; i < candidate->secondaries.size(); i++) {
			if (findActivePosition(candidate->secondaries[i], recursive,
					position))
				return true;
		}
	}

	return false;
}

void SpatialPartitioning::run(candidate_vector_t &candidates, bool recursive) {
	size_t count = candidates.size();
	Loki::AssocVector<Index, Count> partitions;

	for (size_t i = 0; i < count; i++) {
		if (candidates[i]->isActive()) {
			Vector3d pos = candidates[i]->current.getPosition()
					- partitionOrigin;
			Index idx(pos / partitionSize);
			partitions[idx] += 1;
		}
	}

	while (partitions.size() > 0) {

		Loki::AssocVector<Index, Count>::iterator iPartition =
				partitions.begin();
		Index nextPartition = iPartition->first;
		size_t nextCount = iPartition->second.count;
		for (iPartition = partitions.begin(); iPartition != partitions.end();
				iPartition++) {
//			std::cout << "testing: " << iPartition->second << std::endl;
			if (iPartition->second.count > nextCount)
				nextPartition = iPartition->first;
		}

		Vector3d o(nextPartition.x, nextPartition.y, nextPartition.z);
		setCurrentPartition(o * partitionSize);

		// do the loop
		partitions.clear();
#pragma omp parallel shared(partitions)
		{
			Loki::AssocVector<Index, Count> _partitions;
#pragma omp for schedule(dynamic, 1000)
			for (size_t i = 0; i < count; i++) {
				run(candidates[i], recursive, _partitions);
			}

			Loki::AssocVector<Index, Count>::iterator i;
#pragma omp critical
			for (i = _partitions.begin(); i != _partitions.end(); i++) {
//				std::cout << "cp ";
//				std::cout << i->first.x << ", " << i->first.y << ", "
//						<< i->first.z;
//				std::cout << " -> " << i->second << std::endl;
				partitions[i->first] += i->second.count;
			}
		}
	}
}

void SpatialPartitioning::run(Candidate *candidate, bool recursive,
		Loki::AssocVector<Index, Count> &partitions) {
	size_t active = 0;

	while (candidate->isActive()) {
		Vector3d relPos = candidate->current.getPosition() - currentPartition;
		double lo = std::min(relPos.x, std::min(relPos.y, relPos.z));
		double hi = std::max(relPos.x, std::max(relPos.y, relPos.z));
		if ((lo < 0.) || (hi > partitionSize)) {
			break;
		}
		candidate->limitNextStep(lo + 0.1 * partitionSize);
		candidate->limitNextStep(partitionSize - hi + 0.1 * partitionSize);

		moduleList->process(candidate);
	}

	if (candidate->isActive()) {
		Vector3d pos = candidate->current.getPosition() - partitionOrigin;
		Index idx(pos / partitionSize);
//		std::cout << pos / Mpc;
//		std::cout << " -> ";
//		std::cout << idx.x << ", " << idx.y << ", " << idx.z << std::endl;
		partitions[idx] += 1;
	}

// propagate secondaries
	if (recursive) {
		for (size_t i = 0; i < candidate->secondaries.size(); i++)
			run(candidate->secondaries[i], recursive, partitions);
	}
}

void SpatialPartitioning::run(Source *source, size_t count, bool recursive) {
	candidate_vector_t candidates;
	candidates.reserve(count);
	for (size_t i = 0; i < count; i++) {
		ParticleState state;
		source->prepare(state);
		ref_ptr<Candidate> candidate = new Candidate(state);
		candidates.push_back(candidate);
	}

	run(candidates, recursive);
}

void SpatialPartitioning::setPartitionOrigin(const Vector3d &origin) {
	partitionOrigin = origin;
}

void SpatialPartitioning::setPartitionSize(double size) {
	partitionSize = size;
	partitionSizeMargin = size * 1.2;
}

struct _UpdateSimulationVolume {
	Vector3d currentPartition;
	double partitionSize;
	_UpdateSimulationVolume(Vector3d currentPartition, double partitionSize) :
			currentPartition(currentPartition), partitionSize(partitionSize) {

	}

	void operator()(ref_ptr<Module> &module) {
		SimulationVolumeDependentModule *svm =
				dynamic_cast<SimulationVolumeDependentModule *>(module.get());
		if (svm)
			svm->updateSimulationVolume(currentPartition, partitionSize);
	}
};

void SpatialPartitioning::setCurrentPartition(const Vector3d &offset) {
	currentPartition = offset;
	currentPartitionMargin = offset
			- Vector3d(partitionSize * 0.1, partitionSize * 0.1,
					partitionSize * 0.1);
	if (verbose) {
		std::cout << "mpc::SpatialPartitioningExecutor::setCurrentPartition -> "
				<< offset / Mpc << std::endl;
	}
	std::for_each(moduleList->getModules().begin(),
			moduleList->getModules().end(),
			_UpdateSimulationVolume(currentPartitionMargin,
					partitionSizeMargin));
}

void SpatialPartitioning::setVerbose(bool verbose) {
	this->verbose = verbose;
}

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::Index &idx) {
	out << idx.x << ", " << idx.y << ", " << idx.z;
	return out;
}
