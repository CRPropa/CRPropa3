#include "crpropa/module/ParticleCollector.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/Units.h"

#include <fstream>

namespace crpropa {

ParticleCollector::ParticleCollector() : nBuffer(10e6), clone(false), recursive(false)  {
        container.reserve(nBuffer); // for 1e6 candidates ~ 500MB of RAM
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer) : clone(false), recursive(false)  {
	container.reserve(nBuffer);
}

ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone) : recursive(false) {
	container.reserve(nBuffer);
}


ParticleCollector::ParticleCollector(const std::size_t nBuffer, const bool clone, const bool recursive) {
	container.reserve(nBuffer);
}


void ParticleCollector::process(Candidate* c) const {
#pragma omp critical
        {
                if (container.size() < nBuffer)
			if(clone)
		        	container.push_back(c->clone(recursive));
			else
				container.push_back(c);
        }
}

void ParticleCollector::reprocess(Module *action) const {
	for (ParticleCollector::iterator itr = container.begin(); itr != container.end(); ++itr){
		if (clone)
			action->process((*(itr->get())).clone(false));
	       	else
        	        action->process(itr->get());
	}
}

void ParticleCollector::dump(const std::string &filename) const {
	TextOutput output(filename, Output::Everything);
	reprocess(&output);
	output.close();
}

void ParticleCollector::load(const std::string &filename){

        std::ifstream infile(filename.c_str());
        if (!infile.good())
                throw std::runtime_error(
                                "crpropa::ParticleCollector: could not open file " + filename);

        std::string line;

	double lengthScale = Mpc; // default Mpc
	double energyScale = EeV; // default EeV

        while (std::getline(infile,line)) {
                std::stringstream stream(line);
                if (stream.peek() == '#')
                        continue;

		ref_ptr<Candidate> c = new Candidate(); 
		double val_d; int val_i;
		double x, y, z;
		stream >> val_d;
		c->setTrajectoryLength(val_d*lengthScale); // D
		stream >> val_d;
		c->setRedshift(val_d); // z
		stream >> val_i;
		c->setSerialNumber(val_i); // SN
		stream >> val_i;
        	c->current.setId(val_i); // ID
		stream >> val_d;
		c->current.setEnergy(val_d*energyScale); // E
		stream >> x >> y >> z;
		c->current.setPosition(Vector3d(x, y, z)*lengthScale); // X, Y, Z
		stream >> x >> y >> z;
		c->current.setDirection(Vector3d(x, y, z)*lengthScale); // Px, Py, Pz
		stream >> val_i; // SN0
		stream >> val_i;
		c->source.setId(val_i); // ID0
		stream >> val_d;
		c->source.setEnergy(val_d*energyScale);	// E0
		stream >> x >> y >> z;
		c->source.setPosition(Vector3d(x, y, z)*lengthScale); // X0, Y0, Z0
		stream >> x >> y >> z;
		c->source.setDirection(Vector3d(x, y, z)*lengthScale); // P0x, P0y, P0z
		stream >> val_i; // SN1
		stream >> val_i;
		c->created.setId(val_i); // ID1
		stream >> val_d;
		c->created.setEnergy(val_d*energyScale); // E1
	        stream >> x >> y >> z;
                c->created.setPosition(Vector3d(x, y, z)*lengthScale); // X1, Y1, Z1
                stream >> x >> y >> z;
                c->created.setDirection(Vector3d(x, y, z)*lengthScale); // P1x, P1y, P1z

		process(c);
        }
        infile.close();
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
