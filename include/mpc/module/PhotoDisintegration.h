#ifndef PHOTODISINTEGRATION_H_
#define PHOTODISINTEGRATION_H_

#include "mpc/Module.h"
#include "mpc/MersenneTwister.h"

#include <vector>
#include <map>

namespace mpc {

struct DisintegrationMode {
	int channel; // #n #p #H2 #H3 #He3 #He4 emitted
	std::vector<double> y; // mean free path in [m]
};

class PhotoDisintegration: public Module {
private:
	MTRand mtrand;
	std::map<int, std::vector<DisintegrationMode> > modeMap;
	std::vector<double> x; // log10(gamma)
	int cached_id;
	int cached_channel;
	double cached_distance;

public:
	PhotoDisintegration();
	std::string getDescription() const;
	void process(Candidate *candidate, std::vector<Candidate *> &secondaries);
	double getMeanFreePath(std::vector<double> y, double lg);
	void disintegrate(Candidate *candidate,
			std::vector<Candidate *> &secondaries);
	void createSecondary(Candidate *candidate,
			std::vector<Candidate *> &secondaries, int id, double energy);
};

} // namespace mpc

#endif /* PHOTODISINTEGRATION_H_ */
