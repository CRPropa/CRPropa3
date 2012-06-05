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
	/** source position */
	Vector3d position;
	/** differential spectral index, minimum and maximum energy */
	double index, Emin, Emax;
	/** source isotopes ids */
	std::vector<int> isotope;
	/** relative abundance of source isotopes at equal energies */
	std::vector<double> abundance;
	/** cumulative probability of source isotopes */
	std::vector<double> probability;

	CompositeSource();
	CompositeSource(const Vector3d &sposition, double Emin = 5 * EeV,
			double Emax = 1000 * EeV, double index = -1);

	/** add an isotope the source */
	void addToComposition(int id, double abundance);
	/** calculate the probability for each isotope */
	void normalize();
	/** prepare a random particle */
	virtual void prepare(ParticleState &state) const;
};

} // namespace mpc

#endif /* MPC_SOURCE_H */
