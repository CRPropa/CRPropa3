#include "crpropa/module/Tools.h"
#include "crpropa/Clock.h"

#include <iostream>
#include <sstream>

using namespace std;

namespace crpropa {

PerformanceModule::~PerformanceModule() {
	double total = 0;
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		total += m.time;
	}
	cout << "Performance for " << calls << " calls:" << endl;
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		cout << " - " << floor((1000 * m.time / total) + 0.5) / 10 << "% -> "
				<< m.module->getDescription() << ": " << (m.time / calls)
				<< endl;
	}
}

void PerformanceModule::add(Module *module) {
	_module_info info;
	info.module = module;
	info.time = 0;
	modules.push_back(info);
}

void PerformanceModule::process(Candidate *candidate) const {
	vector<double> times(modules.size());
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		double start = Clock::getInstance().getMillisecond();
		m.module->process(candidate);
		double end = Clock::getInstance().getMillisecond();
		times[i] = end - start;
	}

#pragma omp critical
	{
		for (size_t i = 0; i < modules.size(); i++) {
			_module_info &m = modules[i];
			m.time += times[i];
		}
		calls++;
	}
}

string PerformanceModule::getDescription() const {
	stringstream sstr;
	sstr << "PerformanceModule (";
	for (size_t i = 0; i < modules.size(); i++) {
		_module_info &m = modules[i];
		if (i > 0)
			sstr << ", ";
		sstr << m.module->getDescription();
	}
	sstr << ")";
	return sstr.str();
}

PropertyStatistics::PropertyStatistics(const std::string &key) :
		key(key) {
	setDescription("PropertyStatistics: " + key);
}

PropertyStatistics::~PropertyStatistics() {
	std::cout << "Property Statistics for " << key << ":" << std::endl;
	Loki::AssocVector<std::string, size_t>::iterator i;
	for (i = properties.begin(); i != properties.end(); i++) {
		std::cout << i->first << "\t\t-> " << i->second << std::endl;
	}
}

void PropertyStatistics::process(Candidate *candidate) const {
	std::string property;
	if (candidate->getProperty(key, property)) {
#pragma omp critical
		{
			Loki::AssocVector<std::string, size_t>::iterator i =
					properties.find(property);
			if (i == properties.end()) {
				properties[property] = 1;
			} else {
				i->second++;
			}
		}
	}
}

ParticleSelector::ParticleSelector(Module *m) :
		module(m) {
}

ParticleSelector::ParticleSelector(Module *m, std::vector<int> particleIDs) :
		module(m), ids(particleIDs) {
}

void ParticleSelector::add(int particleID) {
	ids.push_back(particleID);
}

void ParticleSelector::process(Candidate* candidate) const {
	int id = candidate->current.getId();
	if (std::find(ids.begin(), ids.end(), id) != ids.end())
		module->process(candidate);
}

std::string ParticleSelector::getDescription() const {
	return "ParticleSelector";
}

} // namespace crpropa
