#ifndef PARTICLEMAPSCONTAINER_HH
#define PARTICLEMAPSCONTAINER_HH 

#include <map>
#include <vector>
#include "crpropa/magneticLens/Pixelization.h"
#include "crpropa/magneticLens/MagneticLens.h"

namespace crpropa{

/// Container for particlemaps
/// The maps are stored with discrete energies on a logarithmic scale. The
/// default energy width is 0.02 with an energy bin from 10**17.99 - 10**18.01
/// eV
class ParticleMapsContainer 
{
	private:
    std::map< int , std::map <int , double*> > _data;				
		Pixelization _pixelization;
    double _deltaLogE;
    double _bin0lowerEdge;

		// get the bin number of the energy
		int energy2Idx(double energy) const;
		double idx2Energy(int idx) const;

  public:

		ParticleMapsContainer(double deltaLogE = 0.02, double bin0lowerEdge = 17.99) : _deltaLogE(deltaLogE), _bin0lowerEdge(bin0lowerEdge), _pixelization(6)
		{
		}

    ~ParticleMapsContainer();

    size_t getNumberOfPixels()
    {
      return _pixelization.getNumberOfPixels();
    }

		/// returns the map for the particleId with the given energy
		double *getMap(const int particleId, double energy);
			
		/// adds a particle to the map container
    /// particleId is HEP particleId, energy [eV], galacticLongitude and
    /// galacticLatitude in [rad]
		void addParticle(const int particleId, double energy, double galacticLongitude, double galacticLatitude, double weight=1);
			
		/// reads the particles from the given file and weights the particles with
		/// the energy of the mother particle at the source.
		/// this assumes that the particles in the fiels are observed at earth and the galactic
		/// center is in directon (0,1,0) and the galactic north-pole is in
		/// direction (0,0,1).
		void addParticlesFromFile(const std::string inputFileName, double sourceEnergyWeightExponent=0);

		// returns a vector of all particle ids in th
		std::vector<int> getParticleIds();

		std::vector<double> getEnergies(int pid);

		void applyLens(MagneticLens &lens);;

		void getRandomParticles(size_t N, vector<int> &particleId, 
			vector<double> &energy, vector<double> &galacticLongitudes,
			vector<double> &galacticLatitudes);
			
};


/*
Python code that we have to write







*/




} // namespace parsec

#endif // PARTICLEMAPSCONTAINER_HH
