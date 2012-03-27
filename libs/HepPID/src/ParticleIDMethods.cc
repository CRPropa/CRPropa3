// ----------------------------------------------------------------------
//
// ParticleIDMethods.cc
//
// ----------------------------------------------------------------------

#include <cmath>	// for pow()

#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"

namespace HepPID {

namespace {

// internal function used by hasXXX methods
bool findQ( const int & pid, const int & q )
{
    if( isDyon(pid) ) { return false; }
    if( isRhadron(pid) ) {
        int iz = 7;
        for( int i=6; i > 1; --i ) {
	   if( digit(location(i),pid) == 0 ) {
	       iz = i;
	   } else if ( i == iz-1 ) {
	       // ignore squark or gluino
	   } else {
	       if( digit(location(i),pid) == q ) { return true; }
	   }
	}
        return false;
    }
    if( digit(nq3,pid) == q || digit(nq2,pid) == q || digit(nq1,pid) == q ) { return true; }
    if( isPentaquark(pid) ) { 
        if( digit(nl,pid) == q || digit(nr,pid) == q ) { return true; }
    }
    return false;
}

}

// absolute value
int abspid( const int & pid )
{
  return (pid < 0) ? -pid : pid;
}

// returns everything beyond the 7th digit (e.g. outside the numbering scheme)
int extraBits( const int & pid )
{
    return abspid(pid)/10000000;
}

//  split the PID into constituent integers
unsigned short digit( location loc, const int & pid )
{
    //  PID digits (base 10) are: n nr nl nq1 nq2 nq3 nj
    //  the location enum provides a convenient index into the PID
    int numerator = (int) std::pow(10.0,(loc-1));
    return (abspid(pid)/numerator)%10;
}

//  return the first two digits if this is a "fundamental" particle
//  ID = 100 is a special case (internal generator ID's are 81-100)
int fundamentalID( const int & pid )
{
    if( extraBits(pid) > 0 ) return 0;
    if( digit(nq2,pid) == 0 && digit(nq1,pid) == 0) {
        return abspid(pid)%10000;
    } else if( abspid(pid) <= 100 ) {
        return abspid(pid);
    } else {
        return 0;
    }
}

// Ion numbers are +/- 10LZZZAAAI. 
int Z( const int & pid )
{
    // a proton can also be a Hydrogen nucleus
    if( abspid(pid) == 2212 ) { return 1; }
    if( isNucleus(pid) ) return (abspid(pid)/10000)%1000;
    return 0;
}

// Ion numbers are +/- 10LZZZAAAI. 
int A( const int & pid )
{
    // a proton can also be a Hydrogen nucleus
    if( abspid(pid) == 2212 ) { return 1; }
    if( isNucleus(pid) ) return (abspid(pid)/10)%1000;
    return 0;
}

// if this is a nucleus (ion), get nLambda
// Ion numbers are +/- 10LZZZAAAI. 
int lambda( const int & pid )
{
    // a proton can also be a Hydrogen nucleus
    if( abspid(pid) == 2212 ) { return 0; }
    if( isNucleus(pid) ) return digit(n8,pid);
    return 0;
}


// ---  boolean methods:
//

//  check to see if this is a valid PID
bool isValid( const int & pid )
{
    if( extraBits(pid) > 0 ) {
        if( isNucleus(pid) )   { return true; }
        if( isQBall(pid) )   { return true; }
        return false; 
    }
    if( isSUSY(pid) ) { return true; }
    if( isRhadron(pid) ) { return true; }
    if( isDyon(pid) ) { return true; }
    // Meson signature
    if( isMeson(pid) )   { return true; }
    // Baryon signature
    if( isBaryon(pid) )  { return true; }
    // DiQuark signature
    if( isDiQuark(pid) ) { return true; }
    // fundamental particle
    if( fundamentalID(pid) > 0 ) { 
      if(pid > 0 ) { 
        return true; 
      } else {
        if( hasFundamentalAnti(pid) ) { return true; }
        return false;
      }
    }
    // pentaquark
    if( isPentaquark(pid) ) { return true; }
    // don't recognize this number
    return false;
}

// if this is a fundamental particle, does it have a valid antiparticle?
bool hasFundamentalAnti( const int & pid )
{
    // these are defined by the generator and therefore are always valid
    if( fundamentalID(pid) <= 100 && fundamentalID(pid) >= 80 ) { return true; }
    // check id's from 1 to 79
    if( fundamentalID(pid) > 0 && fundamentalID(pid) < 80 ) { 
       if( validParticleName(-pid) ) { return true; }
    }
    return false;
}

//  check to see if this is a valid meson
bool isMeson( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( abspid(pid) <= 100 ) { return false; }
    if( fundamentalID(pid) <= 100 && fundamentalID(pid) > 0 ) { return false; }
    if( isRhadron(pid) ) { return false; }
    int aid = abspid(pid);
    if( aid == 130 || aid == 310 || aid == 210 ) { return true; }
    // EvtGen uses some odd numbers
    if( aid == 150 || aid == 350 || aid == 510 || aid == 530 ) { return true; }
    // pomeron, etc.
    if( pid == 110 || pid == 990 || pid == 9990 ) { return true; }
    if(    digit(nj,pid) > 0 && digit(nq3,pid) > 0 
        && digit(nq2,pid) > 0 && digit(nq1,pid) == 0 ) {
        // check for illegal antiparticles
        if( digit(nq3,pid) == digit(nq2,pid) && pid < 0 ) {
            return false;
        } else {
            return true;
        }
    }
    return false;
}

//  check to see if this is a valid baryon
bool isBaryon( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( abspid(pid) <= 100 ) { return false; }
    if( fundamentalID(pid) <= 100 && fundamentalID(pid) > 0 ) { return false; }
    if( isRhadron(pid) ) { return false; }
    if( isPentaquark(pid) ) { return false; }
    if( abspid(pid) == 2110 || abspid(pid) == 2210 ) { return true; }
    if(    digit(nj,pid) > 0  && digit(nq3,pid) > 0 
        && digit(nq2,pid) > 0 && digit(nq1,pid) > 0 ) { return true; }
    return false;
}

//  check to see if this is a valid diquark
bool isDiQuark( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( abspid(pid) <= 100 ) { return false; }
    if( fundamentalID(pid) <= 100 && fundamentalID(pid) > 0 ) { return false; }
    if(    digit(nj,pid) > 0  && digit(nq3,pid) == 0 
        && digit(nq2,pid) > 0 && digit(nq1,pid) > 0 ) {  // diquark signature
       // EvtGen uses the diquarks for quark pairs, so, for instance, 
       //   5501 is a valid "diquark" for EvtGen
       //if( digit(nj) == 1 && digit(nq2) == digit(nq1) ) { 	// illegal
       //   return false; 
       //} else {
          return true;
       //}
    }
    return false;
}

// is this a valid hadron ID?
bool isHadron( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( isMeson(pid) )   { return true; }
    if( isBaryon(pid) )  { return true; }
    if( isPentaquark(pid) ) { return true; }
    if( isRhadron(pid) ) { return true; }
    return false;
}
// is this a valid lepton ID?
bool isLepton( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) >= 11 && fundamentalID(pid) <= 18 ) { return true; }
    return false;
}

//
// This implements the 2006 Monte Carlo nuclear code scheme.
// Ion numbers are +/- 10LZZZAAAI. 
// AAA is A - total baryon number
// ZZZ is Z - total charge
// L is the total number of strange quarks.
// I is the isomer number, with I=0 corresponding to the ground state.
bool isNucleus( const int & pid )
{
     // a proton can also be a Hydrogen nucleus
     if( abspid(pid) == 2212 ) { return true; }
     // new standard: +/- 10LZZZAAAI
     if( ( digit(n10,pid) == 1 ) && ( digit(n9,pid) == 0 ) ) {
        // charge should always be less than or equal to baryon number
	// the following line is A >= Z
        if( (abspid(pid)/10)%1000 >= (abspid(pid)/10000)%1000 ) { return true; }
     }
     return false;
}

//  check to see if this is a valid pentaquark
bool isPentaquark( const int & pid )
{
    // a pentaquark is of the form 9abcdej,
    // where j is the spin and a, b, c, d, and e are quarks
    if( extraBits(pid) > 0 ) { return false; }
    if( digit(n,pid) != 9 )  { return false; }
    if( digit(nr,pid) == 9 || digit(nr,pid) == 0 )  { return false; }
    if( digit(nj,pid) == 9 || digit(nl,pid) == 0 )  { return false; }
    if( digit(nq1,pid) == 0 )  { return false; }
    if( digit(nq2,pid) == 0 )  { return false; }
    if( digit(nq3,pid) == 0 )  { return false; }
    if( digit(nj,pid) == 0 )  { return false; }
    // check ordering
    if( digit(nq2,pid) > digit(nq1,pid) )  { return false; }
    if( digit(nq1,pid) > digit(nl,pid) )  { return false; }
    if( digit(nl,pid) > digit(nr,pid) )  { return false; }
    return true;
}

// is this a SUSY?
bool isSUSY( const int & pid )
{
    // fundamental SUSY particles have n = 1 or 2
    if( extraBits(pid) > 0 ) { return false; }
    if( digit(n,pid) != 1 && digit(n,pid) != 2 )  { return false; }
    if( digit(nr,pid) != 0 )  { return false; }
    // check fundamental part
    if( fundamentalID(pid) == 0 )  { return false; }
    return true;
}

// is this an R-hadron?
bool isRhadron( const int & pid )
{
    // an R-hadron is of the form 10abcdj, 100abcj, or 1000abj
    // where j is the spin, b, c, and d are quarks or gluons,
    // and a (the digit following the zero's) is a SUSY particle
    if( extraBits(pid) > 0 ) { return false; }
    if( digit(n,pid) != 1 )  { return false; }
    if( digit(nr,pid) != 0 )  { return false; }
    // make sure this isn't a SUSY particle
    if( isSUSY(pid) ) { return false; }
    // All R-hadrons have at least 3 core digits
    if( digit(nq2,pid) == 0 )  { return false; }
    if( digit(nq3,pid) == 0 )  { return false; }
    if( digit(nj,pid) == 0 )  { return false; }
    return true;
}

// is this a Dyon (magnetic monopole)?
bool isDyon( const int & pid )
{
    ///Magnetic monopoles and Dyons are assumed to have one unit of 
    ///Dirac monopole charge and a variable integer number xyz units 
    ///of electric charge.  
    ///
    ///Codes 411xyz0 are then used when the magnetic and electrical 
    ///charge sign agree and 412xyz0 when they disagree, 
    ///with the overall sign of the particle set by the magnetic charge.  
    ///For now no spin information is provided.
    ///
    if( extraBits(pid) > 0 ) { return false; }
    if( digit(n,pid) != 4 )  { return false; }
    if( digit(nr,pid) != 1 )  { return false; }
    if( (digit(nl,pid) != 1) && (digit(nl,pid) != 2) )  { return false; }
    // All Dyons have at least 1 core digit
    if( digit(nq3,pid) == 0 )  { return false; }
    // Dyons have spin zero for now
    if( digit(nj,pid) != 0 )  { return false; }
    return true;
}

// Check for QBalls
// Ad-hoc numbering for such particles is 100xxxx0, 
// where xxxx is the charge in tenths. 
bool isQBall( const int & pid )
{
    if( extraBits(pid) != 1 ) { return false; }
    if( digit(n,pid) != 0 )  { return false; }
    if( digit(nr,pid) != 0 )  { return false; }
    // check the core number
    if( (abspid(pid)/10)%10000 == 0 )  { return false; }
    // these particles have spin zero for now
    if( digit(nj,pid) != 0 )  { return false; }
    return true;
}

// does this particle contain an up quark?
bool hasUp( const int & pid)
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) > 0 ) { return false; }
    return findQ(pid,2);
}
// does this particle contain a down quark?
bool hasDown( const int & pid)
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) > 0 ) { return false; }
    return findQ(pid,1);
}
// does this particle contain a strange quark?
bool hasStrange( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) > 0 ) { return false; }
    return findQ(pid,3);
}
// does this particle contain a charm quark?
bool hasCharm( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) > 0 ) { return false; }
    return findQ(pid,4);
}
// does this particle contain a bottom quark?
bool hasBottom( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) > 0 ) { return false; }
    return findQ(pid,5);
}
// does this particle contain a top quark?
bool hasTop( const int & pid )
{
    if( extraBits(pid) > 0 ) { return false; }
    if( fundamentalID(pid) > 0 ) { return false; }
    return findQ(pid,6);
}

// ---  other information
//
// jSpin returns 2J+1, where J is the total spin
int  jSpin( const int & pid )
{
    if( fundamentalID(pid) > 0 ) { 
	// some of these are known
	int fund = fundamentalID(pid);
	if( fund > 0 && fund < 7 ) return 2;
	if( fund == 9 ) return 3; 
	if( fund > 10 && fund < 17 ) return 2;
	if( fund > 20 && fund < 25 ) return 3;
        return 0; 
    } else if( extraBits(pid) > 0 ) { 
        return 0; 
    }
    return abspid(pid)%10;
}
// sSpin returns 2S+1, where S is the spin
int  sSpin( const int & pid )
{
    if( !isMeson(pid) ) { return 0; }
    int inl = digit(nl,pid);
    //int tent = digit(n,pid);
    int js = digit(nj,pid);
    if( digit(n,pid) == 9 ) { return 0; }	// tentative ID
    //if( tent == 9 ) { return 0; }	// tentative assignment
    if( inl == 0 && js >= 3 ) { 
        return 1;
    } else if( inl == 0  && js == 1 ) {
        return 0;
    } else if( inl == 1  && js >= 3 ) {
        return 0;
    } else if( inl == 2  && js >= 3 ) {
        return 1;
    } else if( inl == 1  && js == 1 ) {
        return 1;
    } else if( inl == 3  && js >= 3 ) {
        return 1;
    }
    // default to zero
    return 0;
}
// lSpin returns 2L+1, where L is the orbital angular momentum
int  lSpin( const int & pid )
{
    if( !isMeson(pid) ) { return 0; }
    int inl = digit(nl,pid);
    //int tent = digit(n,pid);
    int js = digit(nj,pid);
    if( digit(n,pid) == 9 ) { return 0; }	// tentative ID
    if( inl == 0 && js == 3 ) { 
        return 0;
    } else if( inl == 0 && js == 5 ) {
        return 1;
    } else if( inl == 0 && js == 7 ) {
        return 2;
    } else if( inl == 0 && js == 9 ) {
        return 3;
    } else if( inl == 0  && js == 1 ) {
        return 0;
    } else if( inl == 1  && js == 3 ) {
        return 1;
    } else if( inl == 1  && js == 5 ) {
        return 2;
    } else if( inl == 1  && js == 7 ) {
        return 3;
    } else if( inl == 1  && js == 9 ) {
        return 4;
    } else if( inl == 2  && js == 3 ) {
        return 1;
    } else if( inl == 2  && js == 5 ) {
        return 2;
    } else if( inl == 2  && js == 7 ) {
        return 3;
    } else if( inl == 2  && js == 9 ) {
        return 4;
    } else if( inl == 1  && js == 1 ) {
        return 1;
    } else if( inl == 3  && js == 3 ) {
        return 2;
    } else if( inl == 3  && js == 5 ) {
        return 3;
    } else if( inl == 3  && js == 7 ) {
        return 4;
    } else if( inl == 3  && js == 9 ) {
        return 5;
    }
    // default to zero
    return 0;
}

// 3 times the charge
int threeCharge( const int & pid )
{
    int charge=0;
    int ida, sid;
    unsigned short q1, q2, q3, ql;
    static int ch100[100] = { -1, 2,-1, 2,-1, 2,-1, 2, 0, 0,
                       -3, 0,-3, 0,-3, 0,-3, 0, 0, 0,
                        0, 0, 0, 3, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 3, 0, 0, 3, 0, 0, 0,
                        0, -1, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 6, 3, 6, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    q1 = digit(nq1,pid);
    q2 = digit(nq2,pid);
    q3 = digit(nq3,pid);
    ql = digit(nl,pid);
    ida = abspid(pid);
    sid = fundamentalID(pid);
    if( ida == 0 ) {      // illegal
        return 0;
    } else if( isQBall(pid) ) {		// QBall
       charge = 3*((abspid(pid)/10)%10000);
    } else if( extraBits(pid) > 0 ) {		// ion
        return 0;
    } else if( isDyon(pid) ) {		// Dyon
        charge = 3*( (abspid(pid)/10)%1000 );
	// this is half right
	// the charge sign will be changed below if pid < 0
	if( ql == 2 ) {
            charge = -charge; 
	}
    } else if( sid > 0 && sid <= 100 ) {	// use table
        charge = ch100[sid-1];
        if(ida==1000017 || ida==1000018) { charge = 0; }
        if(ida==1000034 || ida==1000052) { charge = 0; }
        if(ida==1000053 || ida==1000054) { charge = 0; }
        if(ida==5100061 || ida==5100062) { charge = 6; }
    } else if( digit(nj,pid) == 0 ) { 		// KL, Ks, or undefined
        return 0;
    } else if( isMeson(pid) ) {			// mesons
            if( q2 == 3 || q2 == 5 ) {
                charge = ch100[q3-1] - ch100[q2-1];
            } else {
                charge = ch100[q2-1] - ch100[q3-1];
            }
    } else if( isRhadron(pid) ) {		// Rhadron
        if (( q1 == 0 ) || ( q1 == 9 )) {
            if( q2 == 3 || q2 == 5 ) {
                charge = ch100[q3-1] - ch100[q2-1];
            } else {
                charge = ch100[q2-1] - ch100[q3-1];
            }
        } else if( ql == 0 ) {
            charge = ch100[q3-1] + ch100[q2-1] + ch100[q1-1];
        } else if ( digit(nr,pid) == 0 ) {
            charge = ch100[q3-1] + ch100[q2-1] + ch100[q1-1] + ch100[ql-1];
        }
    } else if( isDiQuark(pid) ) {			// diquarks
        charge = ch100[q2-1] + ch100[q1-1];
    } else if( isBaryon(pid) ) { 			// baryons
        charge = ch100[q3-1] + ch100[q2-1] + ch100[q1-1];
    } else {		// unknown
        return 0;
    }
    if( charge == 0 ) {
        return 0;
    } else if( pid < 0 ) {
        charge = -charge; 
    }
    return charge;
}

// the actual charge
double charge( const int & pid )
{
   int tc = threeCharge(pid);
   if( isQBall(pid) ) {
       return double(tc)/30.;
   } else {
       return double(tc)/3.;
   }
}

} // HepPID
