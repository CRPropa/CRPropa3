// ----------------------------------------------------------------------
//
// translateGeanttoPDT.cc
// Author: Lynn Garren
//
// translate a Geant ID number to the standard numbering scheme
//
// ----------------------------------------------------------------------

#include <iostream>

#include "HepPID/Version.hh"
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"

#define IDMAX 49

namespace HepPID {

int translateGeanttoPDT( const int id )
{
    static int IDG2H[IDMAX] = {
          22,  // gamma
         -11,  // e
          11,  // e-
          12,  // nu (nu_e)
         -13,  // mu
          13,  // mu-
         111,  // pi0
         211,  // pi
        -211,  // pi-
         130,  // K0L
         321,  // K
        -321,  // K-
        2112,  // n
        2212,  // p
       -2212,  // p~
         310,  // K0s
         221,  // eta
        3122,  // Lambda
        3222,  // Sigma
        3212,  // Sigma0
        3112,  // Sigma-
        3322,  // Xi0
        3312,  // Xi-
        3334,  // Omega-
       -2112,  // n~
       -3122,  // Lambda~
       -3112,  // Sigma~
       -3212,  // Sigma0~
       -3222,  // sigma-~
       -3322,  // Xi0~
       -3312,  // Xi-~
       -3334,  // Omega-~
         -15,  // tau
          15,  // tau-
         411,  // D
        -411,  // D-
         421,  // D0
        -421,  // D0~
         431,  // Ds
        -431,  // Ds-
        4122,  // Lambda_
          24,  // W
         -24,  // W-
          23,  // Z
      1000010020,    // deuteron
      1000010030,    // tritium
      1000020040,    // alpha
           0,   // geantino
      1000020030 };   // He3
    static int nwhine_max = 10;
    int nwhine = 0;
    int gtran = 0;

//..............geantino
    if( id == 48 ) {
      nwhine = nwhine + 1;
      if( nwhine <= nwhine_max ) {
	std::cout << "GTRAN: geantino " << id
	          << " not known to HEP (set to 0)" << std::endl;
      }
//...............normal translation
    } else if( id <= IDMAX ) {
      gtran = IDG2H[id-1];
//............anything } else {
    } else {
      nwhine = nwhine + 1;
      if( nwhine <= nwhine_max ) {
	std::cout << "GTRAN: GEANT particle " << id 
	          << " not known to HEP (set to 0)" << std::endl;
      }
    }
    return gtran;
}

}	// HepPID
