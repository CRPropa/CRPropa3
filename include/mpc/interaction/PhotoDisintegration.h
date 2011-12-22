#ifndef PHOTODISINTEGRATION_H_
#define PHOTODISINTEGRATION_H_

#include "mpc/Module.h"
#include "mpc/ParticleState.h"

#include <fstream>

namespace mpc {

struct IndexStruct {
	int id;
	size_t start;
	size_t end;
};

class PhotoDisintegration: public Module {
private:
	int n;
	double dx;
	double xMin;
	std::vector<IndexStruct> nucleusIndex;
	std::vector<IndexStruct> channelIndex;
	std::vector<double> tabulatedMeanFreePath;

public:
	PhotoDisintegration();

	~PhotoDisintegration();

	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	double getLambdaSum(Candidate candidate);

	double interpolate(size_t start, double x);
	void disintegrate(Candidate *candidate,
			std::vector<Candidate *> &secondaries, size_t channel);

	void createSecondary(Candidate *candidate,
			std::vector<Candidate *> &secondaries, int id, double energy);

};

} // namespace mpc

#endif /* PHOTODISINTEGRATION_H_ */
