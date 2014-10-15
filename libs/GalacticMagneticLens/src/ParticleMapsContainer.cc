#include "HepPID/ParticleIDMethods.hh"
#include "parsec/Random.hh"
#include <fstream>
#include "parsec/ParticleMapsContainer.h"
#include <iostream>
namespace parsec 
{

ParticleMapsContainer::~ParticleMapsContainer()
{
  for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) {
		for(std::map<int, double*>::iterator energy_iter = pid_iter->second.begin();
			energy_iter != pid_iter->second.end(); ++energy_iter) {
			  delete[] (energy_iter->second);
		}
	}
}

int ParticleMapsContainer::energy2Idx(double energy) const
{
  double lE = log10(energy);
  return int((lE - _bin0lowerEdge) / _deltaLogE);
}

double ParticleMapsContainer::idx2Energy(int idx) const
{
  return pow(10, idx * _deltaLogE + _bin0lowerEdge + _deltaLogE / 2);
}

		
double* ParticleMapsContainer::getMap(const int particleId, double energy)
{
  if (_data.find(particleId) == _data.end())
  {
    return NULL;
  }
  int energyIdx  = energy2Idx(energy);
  if (_data[particleId].find(energyIdx) == _data[particleId].end())
  {
    return NULL;
  }
  return _data[particleId][energy2Idx(energy)];
}
			
			
void ParticleMapsContainer::addParticle(const int particleId, double energy, double galacticLongitude, double galacticLatitude, double weight)
{
  if (_data.find(particleId) == _data.end())
  {
    map<int, double*> M;
    _data[particleId] = M;
    _weights[particleId] = 0;
  }

  int energyIdx  = energy2Idx(energy);
  if (_data[particleId].find(energyIdx) == _data[particleId].end())
  {
    _data[particleId][energyIdx] = new double[_pixelization.getNumberOfPixels()];
    std::fill(_data[particleId][energyIdx], _data[particleId][energyIdx] + _pixelization.getNumberOfPixels(), 0);
  }

  uint32_t pixel = _pixelization.direction2Pix(galacticLongitude, galacticLatitude);
  _data[particleId][energyIdx][pixel] +=weight;
  _weights[particleId] +=weight;
}
 

void ParticleMapsContainer::addParticlesFromFile(const std::string inputFileName, double sourceEnergyWeightExponent)
{
  std::ifstream infile(inputFileName.c_str());

  while (infile.good()) {
    if (infile.peek() != '#') {
      int particleId, sourceId; 
      double trajectoryLength, energy, sourceEnergy, pX, pY, pZ, x, y, z, pX0, pY0, pZ0, x0, y0, z0, redShift;

      infile >> trajectoryLength 
        >> particleId 
        >> sourceId
        >> energy 
        >> sourceEnergy
        >> x >> y >> z 
        >> x0 >> y0 >> z0
        >> pX >> pY >> pZ 
        >> pX0 >> pY0 >> pZ0 
        >> redShift;
      
      double weight = pow(sourceEnergy, sourceEnergyWeightExponent);
      double galacticLongitude = atan2(-x, -y);
      double galacticLatitude =  acos(z) - M_PI /2;
      addParticle(particleId, energy, galacticLongitude, galacticLatitude, weight);
    }
  }
}


// returns a vector of all particle ids in th
std::vector<int> ParticleMapsContainer::getParticleIds()
{
	std::vector<int> ids;
	for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) 
	{
		ids.push_back(pid_iter->first);
	}
	return ids;
}


std::vector<double> ParticleMapsContainer::getEnergies(int pid)
{
	std::vector<double> energies;
	if (_data.find(pid) != _data.end())
	{
		for(std::map<int, double*>::iterator iter = _data[pid].begin(); 
			iter != _data[pid].end(); ++iter) 
		{
			energies.push_back( idx2Energy(iter->first) );
		}
	}
	return energies;
}


void ParticleMapsContainer::applyLens(MagneticLens &lens)
{
	for(std::map<int, std::map<int, double*> >::iterator pid_iter = _data.begin(); 
			pid_iter != _data.end(); ++pid_iter) {
		for(std::map<int, double*>::iterator energy_iter = pid_iter->second.begin();
			energy_iter != pid_iter->second.end(); ++energy_iter) {
		//	// transform only nuclei
			double energy = idx2Energy(energy_iter->first) / 1E18;
			double charge = HepPID::charge(pid_iter->first);
      if (fabs(charge) > DBL_EPSILON)
      {
			  lens.transformModelVector(energy_iter->second, energy / charge);
      }
			
		}
	}
}


void ParticleMapsContainer::getRandomParticles(size_t N, vector<int> &particleId, 
	vector<double> &energy, vector<double> &galacticLongitudes,
	vector<double> &galacticLatitudes)
{
	double sumOfWeights = 0;
	for(std::map<int, double>::iterator iter = _weights.begin(); iter != _weights.end(); ++iter)
  {
		sumOfWeights += iter->second;
  }

	particleId.resize(N);
	energy.resize(N);
	galacticLongitudes.resize(N);
	galacticLatitudes.resize(N);

	for(size_t i=0; i< N; i++)
	{
		double r = Random::instance().rand() * sumOfWeights;
    std::map<int, double>::iterator iter = _weights.begin();
   
    while ((r-= iter->second) > 0)
    {
      ++iter; 
    }

    particleId[i] = iter->first;
		bool foundParticle = false;

	//	// loop over maps
    for(std::map<int, double*>::iterator energy_iter =
        _data[iter->first].begin(); !foundParticle; ++energy_iter )
		{
			for(size_t j=0; j< _pixelization.getNumberOfPixels() && !foundParticle; j++)
			{
				// if this is too slow I need to store the weights for the energies
				// alongside the maps
				r+= energy_iter->second[j];
				if (r >=0)
				{
					foundParticle = true;
					energy[i] = idx2Energy(energy_iter->first);
					_pixelization.pix2Direction (j, galacticLongitudes[i], galacticLatitudes[i]);
				}
			}
		}
	}
}

} // namespace parsec
