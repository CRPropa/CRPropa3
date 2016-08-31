#include "crpropa/module/ParticleCollector.h"
#include "crpropa/Units.h"

namespace crpropa {

ParticleCollector::ParticleCollector(std::size_t size) {
        nBuffer = size;
        container.reserve(nBuffer); // for 1e6 candidates ~ 500MB of RAM
}

void ParticleCollector::process(Candidate* c) const {
#pragma omp critical
        {
                if (container.size() < nBuffer)
	        	container.push_back(c);
        }
}

ParticleCollector::~ParticleCollector() {
        clearContainer();
}

std::size_t ParticleCollector::getCount() const {
        return container.size();
}

ref_ptr<Candidate> ParticleCollector::operator[](const std::size_t i) const {
	return container[i];
}

void ParticleCollector::clearContainer() {
        container.clear();
}

std::vector<ref_ptr<Candidate> > ParticleCollector::getAll() const {
        return container;
}
std::string ParticleCollector::getDescription() const {
        return "ParticleCollector";
}

ParticleCollector::iterator ParticleCollector::begin()
{
  return container.begin();
}

ParticleCollector::const_iterator ParticleCollector::begin() const
{
  return container.begin();
}

ParticleCollector::iterator ParticleCollector::end()
{
  return container.end();
}

ParticleCollector::const_iterator ParticleCollector::end() const
{
  return container.end();
}

} // namespace crpropa
