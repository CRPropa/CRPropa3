#ifndef GLUTDISPLAY_H_
#define GLUTDISPLAY_H_

#include "mpc/ModuleChain.h"

namespace mpc {

/**
 @class GlutDisplay
 @brief Visualization of the particle trajectory.

 This module displays the particle trajectory in an interactive window.
 */
class GlutDisplay: public Module {
public:
	int counter;

	GlutDisplay();
	~GlutDisplay();

	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namspace mpc

#endif /* GLUTDISPLAY_H_ */

