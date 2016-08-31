#include "crpropa/module/ParticleContainerOutput.h"
#include "crpropa/Units.h"

namespace crpropa {

ParticleContainerOutput::ParticleContainerOutput(std::size_t size) {
        nBuffer = size;
        container.reserve(nBuffer); // for 1e6 candidates ~ 500MB of RAM
}

void ParticleContainerOutput::process(Candidate* c) const {
#pragma omp critical
        {
                if (container.size() < nBuffer)
	        	container.push_back(c);
        }
}

ParticleContainerOutput::~ParticleContainerOutput() {
        clearContainer();
}

std::size_t ParticleContainerOutput::getCount() const {
        return container.size();
}

ref_ptr<Candidate> ParticleContainerOutput::operator[](const std::size_t i) const {
	return container[i];
}

void ParticleContainerOutput::clearContainer() {
        container.clear();
}

std::vector<ref_ptr<Candidate> > ParticleContainerOutput::getAll() const {
        return container;
}
std::string ParticleContainerOutput::getDescription() const {
        return "ParticleContainerOutput";
}

ParticleContainerOutput::iterator ParticleContainerOutput::begin()
{
  return container.begin();
}

ParticleContainerOutput::const_iterator ParticleContainerOutput::begin() const
{
  return container.begin();
}

ParticleContainerOutput::iterator ParticleContainerOutput::end()
{
  return container.end();
}

ParticleContainerOutput::const_iterator ParticleContainerOutput::end() const
{
  return container.end();
}

} // namespace crpropa
