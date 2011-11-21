/*
 * Propagator.h
 *
 *  Created on: 18.11.2011
 *      Author: gmueller
 */

#ifndef PROPAGATOR_H_
#define PROPAGATOR_H_

#include <list>
#include <typeinfo>

namespace mpc {

struct Priority {
	enum Enum {
		Start = 0,
		BeforeInteractions = 20,
		Interactions = 25,
		AfterInteractions = 30,
		BeforeIntegration = 45,
		Integration = 50,
		AfterIntegration = 55,
		BeforeCommit = 70,
		Commit = 75,
		AfterCommit = 80,
		End = 100
	};
};

class Feature {
public:
	virtual ~Feature() {
	}

	virtual void apply(Candidate &c, size_t priority) = 0;
};

class Propagator {
	struct FeatureEntry {
		size_t priority;
		Feature *feature;
		bool operator<(FeatureEntry const& rhs) const {
			return (priority < rhs.priority);
		}
	};
	std::list<FeatureEntry> startFeatures;
	std::list<FeatureEntry> mainFeatures;
	std::list<FeatureEntry> endFeatures;

	void check() {
		size_t integratorCount = 0;

		std::list<FeatureEntry>::iterator iEntry = mainFeatures.begin();
		while (iEntry != mainFeatures.end()) {
			FeatureEntry &entry = *iEntry;
			iEntry++;

			if (entry.priority == Priority::Integration) {
				integratorCount++;
			}
		}

		if (integratorCount > 0) {
			std::cerr << "Warning: more than one integration feature present."
					<< std::endl;
		}
	}
public:

	void add(size_t priority, Feature *feature) {
		FeatureEntry entry;
		entry.priority = priority;
		entry.feature = feature;

		if (priority == Priority::Start) {
			startFeatures.push_back(entry);
		} else if (priority == Priority::End) {
			endFeatures.push_back(entry);
		} else {
			mainFeatures.push_back(entry);
			mainFeatures.sort();
		}

		check();
	}

	void apply(Candidate &c) {
		std::list<FeatureEntry>::iterator iStartEntry = startFeatures.begin();
		while (iStartEntry != startFeatures.end()) {
			FeatureEntry &entry = *iStartEntry;
			iStartEntry++;

			entry.feature->apply(c, entry.priority);
		}

		while (c.getStatus() == Candidate::Active) {
			std::list<FeatureEntry>::iterator iEntry = mainFeatures.begin();
			while (iEntry != mainFeatures.end()) {
				FeatureEntry &entry = *iEntry;
				iEntry++;

				entry.feature->apply(c, entry.priority);
			}
		}

		std::list<FeatureEntry>::iterator iEndEntry = endFeatures.begin();
		while (iEndEntry != startFeatures.end()) {
			FeatureEntry &entry = *iEndEntry;
			iEndEntry++;

			entry.feature->apply(c, entry.priority);
		}
	}

	void print() {
		std::list<FeatureEntry>::iterator iEntry = mainFeatures.begin();
		while (iEntry != mainFeatures.end()) {
			FeatureEntry &entry = *iEntry;
			iEntry++;

			std::cout << entry.priority << " -> "
					<< typeid(*entry.feature).name() << std::endl;
		}
	}
};

} // namespace mpc

#endif /* PROPAGATOR_H_ */
