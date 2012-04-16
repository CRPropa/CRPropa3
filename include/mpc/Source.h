#ifndef MPC_SOURCE_H
#define MPC_SOURCE_H

#include "mpc/ParticleState.h"
#include "mpc/Referenced.h"

#include <vector>

namespace mpc {

/**
 @class Source
 @brief Base cosmic ray source
 */
class Source: public Referenced {
public:
	virtual ~Source() {
	}

	virtual void prepare(ParticleState &state) const = 0;
};

/**
 @class BasicSource
 @brief Simple cosmic ray source
 */
class BasicSource: public Source {
public:
	Vector3 position;
	double index1, index2, breakpoint, Emin, Emax;
	int id;
	BasicSource();
	BasicSource(const Vector3 &sposition, int type, double Emin = 5 * EeV,
			double Emax = 1000 * EeV, double index1 = -1, double index2 = -1,
			double breakpoint = 1);

	virtual void prepare(ParticleState &state) const;
};

} // namespace mpc

#endif /* MPC_SOURCE_H */
