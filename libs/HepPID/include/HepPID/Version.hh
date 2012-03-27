// ----------------------------------------------------------------------
//
// Version.hh
// Author:  Lynn Garren
//
//  for now, this is a free function
//
// ----------------------------------------------------------------------
#ifndef HepPIDVERSION_HH
#define HepPIDVERSION_HH

#include <string>
#include <iostream>

namespace HepPID {

void version( );			//!< print HepPID version
void writeVersion( std::ostream & os );	//!< write HepPID version to os
std::string versionName( );	//!< return HepPID version

}	// HepPID

#endif // HepPIDVERSION_HH
