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
		partitionOrigin(0, 0, 0), partitionSize(10 * Mpc), partitionMarginInner(
				0.5 * Mpc), partitionMarginOuter(1 * Mpc), verbose(false), moduleList(
				moduleList) {
	updateMargins();

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

bool toBeRemoved(ref_ptr<Candidate> &candidate) {
	bool inactive = (candidate->isActive() == false);
	bool nochild = (candidate->secondaries.size() == 0);
	return (inactive && nochild);
}

void SpatialPartitioning::run(candidate_vector_t &candidates, bool recursive,
		bool deleteInactive) {
	Loki::AssocVector<Index, Count> partitions;

	for (size_t i = 0; i < candidates.size(); i++) {
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
			if (iPartition->second.count > nextCount) {
				nextPartition = iPartition->first;
				nextCount = iPartition->second.count;
			}
		}

		Vector3d o(nextPartition.x, nextPartition.y, nextPartition.z);
		setCurrentPartition(o * partitionSize);

		// do the loop
		partitions.clear();
		size_t cent = std::max(1ul, candidates.size() / 100), pc = 0;

#pragma omp parallel shared(partitions, pc)
		{
			Loki::AssocVector<Index, Count> _partitions;
#pragma omp for schedule(dynamic, 1000)
			for (size_t i = 0; i < candidates.size(); i++) {
				if (verbose && (i % cent == 0)) {
					std::cout << pc << "% - " << i << std::endl;
					pc++;
				}
				run(candidates[i], recursive, _partitions, deleteInactive);
			}

			Loki::AssocVector<Index, Count>::iterator i;
#pragma omp critical
			for (i = _partitions.begin(); i != _partitions.end(); i++) {
				partitions[i->first] += i->second.count;
			}
		}

		if (deleteInactive) {
			candidates.erase(
					remove_if(candidates.begin(), candidates.end(),
							toBeRemoved), candidates.end());
		}
	}
}

void SpatialPartitioning::run(Candidate *candidate, bool recursive,
		Loki::AssocVector<Index, Count> &partitions, bool deleteInactive) {
	size_t active = 0;

	while (candidate->isActive()) {
		Vector3d relPos = candidate->current.getPosition()
				- currentPartitionInner;
		double lo = std::min(relPos.x, std::min(relPos.y, relPos.z));
		double hi = std::max(relPos.x, std::max(relPos.y, relPos.z));
		if ((lo <= 0.) || (hi >= partitionSizeInner)) {
			break;
		}
		candidate->limitNextStep(lo);
		candidate->limitNextStep(partitionSizeInner - hi);

		moduleList->process(candidate);
	}

	if (candidate->isActive()) {
		Vector3d pos = candidate->current.getPosition() - partitionOrigin;
		Index idx(pos / partitionSize);
		partitions[idx] += 1;
	}

// propagate secondaries
	if (recursive) {
		for (size_t i = 0; i < candidate->secondaries.size(); i++) {
			run(candidate->secondaries[i], recursive, partitions,
					deleteInactive);
		}
		if (deleteInactive) {
			candidate->secondaries.erase(
					remove_if(candidate->secondaries.begin(),
							candidate->secondaries.end(), toBeRemoved),
					candidate->secondaries.end());
		}
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

	run(candidates, recursive, true);
}

void SpatialPartitioning::setPartitionOrigin(const Vector3d &origin) {
	partitionOrigin = origin;
}

void SpatialPartitioning::setPartitionSize(double size) {
	partitionSize = size;
	updateMargins();

}

void SpatialPartitioning::setPartitionMargin(double inner, double outer) {
	partitionMarginInner = inner;
	partitionMarginOuter = outer;
	updateMargins();
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
	updateMargins();
	if (verbose) {
		std::cout << "mpc::SpatialPartitioningExecutor::setCurrentPartition -> "
				<< offset / Mpc << std::endl;
	}

	std::for_each(moduleList->getModules().begin(),
			moduleList->getModules().end(),
			_UpdateSimulationVolume(currentPartitionOuter, partitionSizeOuter));
}

void SpatialPartitioning::updateMargins() {
	currentPartitionInner = currentPartition
			- Vector3d(partitionMarginInner, partitionMarginInner,
					partitionMarginInner);
	partitionSizeInner = partitionSize + 2 * partitionMarginInner;

	currentPartitionOuter = currentPartition
			- Vector3d(partitionMarginOuter, partitionMarginOuter,
					partitionMarginOuter);
	partitionSizeOuter = partitionSize + 2 * partitionMarginOuter;
}

void SpatialPartitioning::setVerbose(bool verbose) {
	this->verbose = verbose;
}

} // namespace mpc

std::ostream &operator<<(std::ostream &out, const mpc::Index &idx) {
	out << idx.x << ", " << idx.y << ", " << idx.z;
	return out;
}
