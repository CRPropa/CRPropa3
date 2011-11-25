#ifndef GLUTDISPLAY_H_
#define GLUTDISPLAY_H_

#include "mpc/Propagator.h"

namespace mpc {

class GlutDisplay: public Feature {
public:
	int counter;
	int refresh;

	GlutDisplay(double r);
	~GlutDisplay();

	void apply(Candidate &candidate);
	std::string description() const;
};

} // namspace mpc

#endif /* GLUTDISPLAY_H_ */

