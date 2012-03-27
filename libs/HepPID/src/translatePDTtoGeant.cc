// ----------------------------------------------------------------------
//
// translatePDTtoGeant.cc
// Author: Lynn Garren
//
//     Translate from the standard numbering scheme to GEANT.
//     Particles unknown to GEANT are entered as geantino's.
//     A warning message is also printed.
//
// ----------------------------------------------------------------------

#include <iostream>

#include "HepPID/Version.hh"
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"

#define IDMAX 49

namespace HepPID {

int translatePDTtoGeant( const int id )
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
      1000010020,   	// deuteron
      1000010030,    	// tritium
      1000020040,  	// alpha
           0,   	// geantino
      1000020030 };   	// He3
    static int nwhine_max = 10;
    int nwhine = 0;
    int gtran, i;

//...........first check if its an e, mu or tau neutrino
    if( abspid(id) == 12 || abspid(id) == 14 || abspid(id) == 16 ) {
       gtran = 4;
    } else {
//..............loop over all GEANT particles, see if it matches
//..............the current HEP particle
       gtran = IDMAX + 1;
       for( i=0; i < IDMAX; ++i ) {
           if( IDG2H[i] == id ) { gtran = i; }
       }
    }
    if( gtran > IDMAX || gtran == 48 ) {
       gtran = 48;	// geantino
       nwhine = nwhine + 1;
       if( nwhine <= nwhine_max ) {
	 std::cout << "translatePDTtoGeant: HEP particle " << id 
	           <<  " not known to GEANT (converted to geantino)" << std::endl;
       }
    }
    return gtran;
}

}	// HepPID
