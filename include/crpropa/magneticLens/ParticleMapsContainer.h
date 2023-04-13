#ifndef CRPROPA_PARTICLEMAPSCONTAINER_HH
#define CRPROPA_PARTICLEMAPSCONTAINER_HH

#include <map>
#include <vector>
#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/magneticLens/MagneticLens.h"

#include "crpropa/Vector3.h"

namespace crpropa {
/**
 * \addtogroup MagneticLenses 
 * @{
 */

/** 
 @class ParticleMapsContainer
 @brief A container for particle maps

 The maps are stored with discrete energies on a logarithmic scale. The
 default energy width is 0.02 with an energy bin from 10**17.99 - 10**18.01 eV.
 */
class ParticleMapsContainer {
private:
	std::map<int, std::map <int, double*> > _data;
	Pixelization _pixelization;
	double _deltaLogE;
	double _bin0lowerEdge;

	// get the bin number of the energy
	int energy2Idx(double energy) const;
	double idx2Energy(int idx) const;

	// weights of the particles
	double _sumOfWeights;
	std::map<int, double > _weightsPID;
	std::map<int, map<int, double> > _weights_pidEnergy;

	// lazy update of weights
	bool _weightsUpToDate;
	void _updateWeights();

public:
	/** Constructor.
	 @param deltaLogE		width of logarithmic energy bin [in eV]
	 @param bin0lowerEdge	logarithm of energy of the lower edge of first bin [in log(eV)]
	 */
	ParticleMapsContainer(double deltaLogE = 0.02, double bin0lowerEdge = 17.99) : _deltaLogE(deltaLogE), _bin0lowerEdge(bin0lowerEdge), _pixelization(6), _weightsUpToDate(false), _sumOfWeights(0) {
	}
	/** Destructor.
	 */
	~ParticleMapsContainer();

	size_t getNumberOfPixels() {
		return _pixelization.getNumberOfPixels();
	}

	/** Get the map for the particleId with the given energy.
	 @param particleId		id of the particle following the PDG numbering scheme
	 @param energy			the energy of the particle [in Joules]
	 @returns The map for a given particleId with a given energy
	 */
	double *getMap(const int particleId, double energy);

	/** Adds a particle to the map container.
	 @param particleId			id of the particle following the PDG numbering scheme
	 @param energy				the energy of the particle [in Joules]
	 @param galacticLongitude	galactic longitude [radians]
	 @param galacticLatitude	galactic latitude [radians]
	 @param weight				relative weight for the specific particle
	*/
	void addParticle(const int particleId, double energy, double galacticLongitude, double galacticLatitude, double weight = 1);
	/** Adds a particle to the map container.
	 @param particleId			id of the particle following the PDG numbering scheme
	 @param energy				the energy of the particle [in Joules]
	 @param v					vector containing the arrival directions of a particle
	 @param weight				relative weight for the specific particle
	*/
	void addParticle(const int particleId, double energy, const Vector3d &v, double weight = 1);

	/** Get all particle ids in the map.
	 @returns Vector of all ids.
	 */
	std::vector<int> getParticleIds();

	/** Get energies in map.
	 @param pid	id of the particle following the PDG numbering scheme
	 @returns Energies are returned in units of eV (unlike in other CRPropa modules) for performance reasons
	 */
	std::vector<double> getEnergies(int pid);

	void applyLens(MagneticLens &lens);;

	/** Get random particles from map.
	 The arguments are the vectors where the information will be stored.
	 @param N					number of particles to be selected
	 @param particleId			id of the particle following the PDG numbering scheme
	 @param energy				energy of interest [in eV]
	 @param galacticLongitudes	longitude in the interval [-pi, pi] [in radians]
	 @param galacticLatitudes	latitude in the interval [-pi/2, pi/2] [in radians]
	 */
	void getRandomParticles(size_t N, vector<int> &particleId,
		vector<double> &energy, vector<double> &galacticLongitudes,
		vector<double> &galacticLatitudes);

	/** Places a particle with given id and energy according to the  probability maps. 
	 @param pid					id of the particle following the PDG numbering scheme
	 @param energy				energy of interest [in eV]
	 @param galacticLongitude	longitude in the interval [-pi, pi] [in radians]
	 @param galacticLatitude	latitude in the interval [-pi/2, pi/2] [in radians]
	 @returns Returns false if operation not possible; true otherwise.
	 */
	bool placeOnMap(int pid, double energy, double &galacticLongitude, double &galacticLatitude);

	/** Force weight update prior to getting random particles. 
	 Only necessary when reusing pointer to maps after calculating weights.
	*/
	void forceWeightUpdate();

	/** Get sum of weights of the maps.
	 @returns Sum of all weights.
	 */
	double getSumOfWeights() {
		if (!_weightsUpToDate)
			_updateWeights();
		return _sumOfWeights;
	}

	/** Get weight for a given particle and energy
	 @param pid					id of the particle following the PDG numbering scheme
	 @param energy				energy of interest [in eV]
	 @returns Weight for the chosen particle and energy.
	 */
	double getWeight(int pid, double energy) {
		if (!_weightsUpToDate)
			_updateWeights();
		return _weights_pidEnergy[pid][energy2Idx(energy)];
	}
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PARTICLEMAPSCONTAINER_HH
