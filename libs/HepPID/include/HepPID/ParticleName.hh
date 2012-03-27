#ifndef PARTICLENAME_HH
#define PARTICLENAME_HH
// ----------------------------------------------------------------------
//
// ParticleName.hh
// Author: Lynn Garren and Walter Brown
//
//  Create a map that gives a standard name for each pre-defined 
//  particle ID number.  This map is initialized if and only if 
//  the public functions are called. Because the map is static, 
//  the initialization happens only once.
//
//
// ----------------------------------------------------------------------

#include <string>
#include <map>
#include <iostream>

namespace HepPID {

/// get a known HepPID Particle name
std::string  particleName( const int & );
/// lookup a known ID
int          particleName( const std::string & );

/// list all known names
void  listParticleNames( std::ostream & os );

/// verify that this number has a valid name
bool validParticleName( const int & );
/// verify that this string has a valid id
bool validParticleName( const std::string & );

// forward definition of ParticleNameMap class
class ParticleNameMap;
/// access the ParticleNameMap for other purposes
ParticleNameMap const &  getParticleNameMap();

}  // namespace HepPID

#endif // PARTICLENAME_HH
