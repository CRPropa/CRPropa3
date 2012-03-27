// ----------------------------------------------------------------------
//
// version.cc
// Author: Lynn Garren
//
//  for now, this is a free function
//
// ----------------------------------------------------------------------

#include "HepPID/Version.hh"

namespace HepPID {

std::string versionName( )
{
    return "3.04.01";
}

void version( )
{
    std::cout << " --------------- HepPID Version " << versionName()
              << " --------------- " << std::endl;
}

void writeVersion( std::ostream & os )
{
    os << "             HepPID Version: " << versionName() << std::endl;
}

}	// HepPID
