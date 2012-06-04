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
 @brief Cosmic ray source with (broken) power law spectrum
 */
class BasicSource: public Source {
public:
	Vector3d position;
	double index1, index2, breakpoint, Emin, Emax;
	int id;
	BasicSource();
	BasicSource(const Vector3d &sposition, int id, double Emin = 5 * EeV,
			double Emax = 1000 * EeV, double index1 = -1, double index2 = -1,
			double breakpoint = 1);

	virtual void prepare(ParticleState &state) const;
};

/**
 @class CompositeSource
 @brief Cosmic ray source with composition and power law spectrum
 */
class CompositeSource: public Source {
public:
	/**
	 @class Isotope
	 @brief Source isotope
	 */
	struct Isotope {
	public:
		int id;
		double abundance;
		double probability;
		Isotope();
		Isotope(int id, double abundance);
	};

	Vector3d position;
	double index, Emin, Emax;
	std::vector<Isotope> composition;
	CompositeSource();
	CompositeSource(const Vector3d &sposition, double Emin = 5 * EeV,
			double Emax = 1000 * EeV, double index = -1);
	void addToComposition(int id, double abundance);
	void normalize();
	virtual void prepare(ParticleState &state) const;
};

} // namespace mpc

#endif /* MPC_SOURCE_H */
