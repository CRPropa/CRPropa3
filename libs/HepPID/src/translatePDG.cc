// ------------------------------------
//
// translatePDG.cc
// Author: Lynn Garren
//
// translate an ID number to or from the standard numbering scheme and the PDG table
// use static maps
//
//  The maps are initialized if and only if the public functions are called.
//  Because the maps are static, the initialization happens only once.
//
//  The user NEVER calls PDGtoPDTMapInit()
//  We use a data table (struct SList) so that compile time is not impacted.
//
//  public functions:
//        int translatePDGtabletoPDT( const int id )
//        int translatePDTtoPDGtable( const int id )
//        PDGtoPDTMap const & getPDGtoPDTMap()
//        PDTtoPDGMap const & getPDTtoPDGMap()
//
// ------------------------------------

#include <map>
#include <utility>	// make_pair

#include "HepPID/Version.hh"
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"

namespace HepPID {

 typedef  std::map< int, int >  PDGtoPDTMap;
 typedef  std::map< int, int >  PDTtoPDGMap;

namespace {	// PDGtoPDTMapInit is private

 PDGtoPDTMap const & getPDGtoPDTMap();
 PDTtoPDGMap const & getPDTtoPDGMap();

PDGtoPDTMap const & PDGtoPDTMapInit()
{

  static PDGtoPDTMap  m;

  static const struct {
      int hid;	// PDG
      int pid;	// PDT
  } SList[] = {
    {            1,            1 },
    {           -1,           -1 },
    {            2,            2 },
    {           -2,           -2 },
    {            3,            3 },
    {           -3,           -3 },
    {            4,            4 },
    {           -4,           -4 },
    {            5,            5 },
    {           -5,           -5 },
    {            6,            6 },
    {           -6,           -6 },
    {            7,            7 },
    {           -7,           -7 },
    {            8,            8 },
    {           -8,           -8 },
    {           11,           11 },
    {          -11,          -11 },
    {           12,           12 },
    {          -12,          -12 },
    {           13,           13 },
    {          -13,          -13 },
    {           14,           14 },
    {          -14,          -14 },
    {           15,           15 },
    {          -15,          -15 },
    {           16,           16 },
    {          -16,          -16 },
    {           17,           17 },
    {          -17,          -17 },
    {           18,           18 },
    {          -18,          -18 },
    {           21,           21 },
    {           22,           22 },
    {           23,           23 },
    {           24,           24 },
    {          -24,          -24 },
    {           25,           25 },
    {           32,           32 },
    {           33,           33 },
    {           34,           34 },
    {          -34,          -34 },
    {           35,           35 },
    {           36,           36 },
    {           37,           37 },
    {          -37,          -37 },
    {           39,           39 },
    {           41,           41 },
    {          -41,          -41 },
    {           42,           42 },
    {          -42,          -42 },
    {          110,          110 },
    {          990,          990 },
    {         9990,         9990 },
    {          130,          130 },
    {          310,          310 },
    {          211,          211},
    {          111,          111},
    {          221,          221},
    {          113,          113},
    {          213,          213},
    {          223,          223},
    {          331,          331},
    {      9010221,      9010221},
    {      9000111,      9000111},
    {      9000211,      9000211},
    {          333,          333},
    {        10223,        10223},
    {        10113,        10113},
    {        10213,        10213},
    {        20113,        20113},
    {        20213,        20213},
    {          225,          225},
    {        20223,        20223},
    {       100221,       100221},
    {       100111,       100111},
    {       100211,       100211},
    {          115,          115},
    {          215,          215},
    {      9000113,      9000113},
    {      9000213,      9000213},
    {      9020221,      9020221},
    {        20333,        20333},
    {        10111,        10111},
    {        10211,        10211},
    {       100113,       100113},
    {       100213,       100213},
    {       100331,       100331},
    {      9030221,      9030221},
    {          335,          335},
    {      9010113,      9010113},
    {      9010213,      9010213},
    {        10225,        10225},
    {        30223,        30223},
    {          227,          227},
    {        10115,        10115},
    {        10215,        10215},
    {       100333,       100333},
    {          117,          117},
    {          217,          217},
    {        30113,        30113},
    {        30213,        30213},
    {        10331,        10331},
    {      9010111,      9010111},
    {      9010211,      9010211},
    {          337,          337},
    {      9050225,      9050225},
    {      9060225,      9060225},
    {          119,          119},
    {          219,          219},
    {          229,          229},
    {      9080225,      9080225},
    {      9090225,      9090225},
    {          321,          321},
    {          311,          311},
    {          323,          323},
    {          313,          313},
    {        10313,        10313},
    {        10323,        10323},
    {        20313,        20313},
    {        20323,        20323},
    {       100313,       100313},
    {       100323,       100323},
    {        10311,        10311},
    {        10321,        10321},
    {          325,          325},
    {          315,          315},
    {        30313,        30313},
    {        30323,        30323},
    {        10315,        10315},
    {        10325,        10325},
    {          317,          317},
    {          327,          327},
    {        20315,        20315},
    {        20325,        20325},
    {          319,          319},
    {          329,          329},
    {          411,          411},
    {          421,          421},
    {          423,          423},
    {          413,          413},
    {        10423,        10423},
    {          425,          425},
    {          415,          415},
    {          431,          431},
    {          433,          433},
    {        10431,        10431},
    {        20433,        20433},
    {        10433,        10433},
    {          435,          435},
    {          521,          521},
    {          511,          511},
    {          513,          513},
    {          523,          523},
    {          531,          531},
    {          541,          541},
    {          441,          441},
    {          443,          443},
    {        10441,        10441},
    {        20443,        20443},
    {          445,          445},
    {       100443,       100443},
    {        30443,        30443},
    {      9000443,      9000443},
    {      9010443,      9010443},
    {      9020443,      9020443},
    {          553,          553},
    {        10551,        10551},
    {        20553,        20553},
    {          555,          555},
    {       100553,       100553},
    {       110551,       110551},
    {       120553,       120553},
    {       100555,       100555},
    {       200553,       200553},
    {       300553,       300553},
    {      9000553,      9000553},
    {      9010553,      9010553},
    {         2212,         2212},
    {        12212,        12212},
    {        22212,        22212},
    {        32212,        32212},
    {        42212,        42212},
    {         2112,         2112},
    {        12112,        12112},
    {        22112,        22112},
    {        32112,        32112},
    {        42112,        42112},
    {         1114,         1114},
    {         2114,         2114},
    {         2214,         2214},
    {         2224,         2224},
    {        31114,         31114},
    {        32114,         32114},
    {        32214,         32214},
    {        32224,         32224},
    {         1112,         1112},
    {         1212,         1212},
    {         2122,         2122},
    {         2222,         2222},
    {        11114,         11114},
    {        12114,         12114},
    {        12214,         12214},
    {        12224,         12224},
    {         1116,         1116},
    {         1216,         1216},
    {         2126,         2126},
    {         2226,         2226},
    {        21112,        21112},
    {        21212,        21212},
    {        22122,        22122},
    {        22222,        22222},
    {        21114,        21114},
    {        22114,        22114},
    {        22214,        22214},
    {        22224,        22224},
    {        11116,        11116},
    {        11216,        11216},
    {        12126,        12126},
    {        12226,        12226},
    {         1118,         1118},
    {         2118,         2118},
    {         2218,         2218},
    {         2228,         2228},
    {         1214,         1214},
    {         2124,         2124},
    {        21214,        21214},
    {        22124,        22124},
    {        31214,        31214},
    {        32124,        32124},
    {         2116,         2116},
    {         2216,         2216},
    {        12116,        12116},
    {        12216,        12216},
    {         1218,         1218},
    {         2128,         2128},
    {         3122,         3122},
    {    13122,    13122 },
    {    23122,    23122 },
    {    33122,    33122 },
    {    43122,    43122 },
    {    53122,    53122 },
    {     3124,     3124 },
    {    13124,    13124 },
    {    23124,    23124 },
    {     3126,     3126 },
    {    13126,    13126 },
    {     3128,     3128 },
    {    23126,    23126 },
    {         3222,         3222},
    {         3212,         3212},
    {         3112,         3112},
    {    13112,    13112 },
    {    13212,    13212 },
    {    13222,    13222 },
    {    23112,    23112 },
    {    23212,    23212 },
    {    23222,    23222 },
    {         3224,         3224},
    {         3214,         3214},
    {         3114,         3114},
    {    13224,    13224 },
    {    13214,    13214 },
    {    13114,    13114 },
    {    23224,    23224 },
    {    23214,    23214 },
    {    23114,    23114 },
    {     3226,     3226 },
    {     3216,     3216 },
    {     3116,     3116 },
    {    13226,    13226 },
    {    13216,    13216 },
    {    13116,    13116 },
    {     3228,     3228 },
    {     3218,     3218 },
    {     3118,     3118 },
    {         3322,         3322},
    {         3312,         3312},
    {         3324,         3324},
    {        13314,        13314},
    {        13324,        13324},
    {         3314,         3314},
    {         3334,         3334},
    {         4122,         4122},
    {        14122,        14122},
    {         4222,         4222},
    {         4212,         4212},
    {         4112,         4112},
    {         4232,         4232},
    {         4132,         4132},
    {         4322,         4322},
    {         4312,         4312},
    {         4332,         4332},
    {         5122,         5122},
    {      9221132,      9221132},
    {      9331122,      9331122}
  };

  int listSize = sizeof(SList)/sizeof(SList[0]);
  for( int k=0; k!=listSize; ++k) {
      m.insert( std::make_pair( SList[k].hid, SList[k].pid) );
  }
  return m;
}  // PDGtoPDTMapInit()

PDTtoPDGMap const & PDTtoPDGMapInit()
{
    static PDTtoPDGMap m;
    static PDGtoPDTMap const & hmap = getPDGtoPDTMap();
    
    for(PDGtoPDTMap::const_iterator cit=hmap.begin(), mend=hmap.end(); cit!=mend; ++cit ) {
	m.insert( std::make_pair( cit->second, cit->first ));
    }
    return m;
}
  
//
// getPDGtoPDTMap is the ONLY function allowed to call PDGtoPDTMapInit
//
PDGtoPDTMap const & getPDGtoPDTMap()
{
  static PDGtoPDTMap const & hmap = PDGtoPDTMapInit();
  return hmap;
}  // getPDGtoPDTMap()

//
// getPDTtoPDGMap is the ONLY function allowed to call PDGtoPDTMapInit
//
PDTtoPDGMap const & getPDTtoPDGMap()
{
  static PDTtoPDGMap const & hmap = PDTtoPDGMapInit();
  return hmap;
}  // getPDTtoPDGMap()
 
} // unnamed namespace
  
int translatePDGtabletoPDT( const int id )
{
    static PDGtoPDTMap const & hmap = getPDGtoPDTMap();

    PDGtoPDTMap::const_iterator const cit = hmap.find( id );
    // found it in the map
    if ( cit != hmap.end() ) { return cit->second; }
    // check to see if someone has defined a valid particle type
    // that isn't in the map
    if( isValid(id) ) { return id; }
    return 0;
}

int translatePDTtoPDGtable( const int id )
{
    static PDTtoPDGMap const & pmap = getPDTtoPDGMap();

    PDTtoPDGMap::const_iterator const cit = pmap.find( id );
    // found it in the map
    if ( cit != pmap.end() ) { return cit->second; }
    // check to see if someone has defined a valid particle type
    // that isn't in the map
    if( isValid(id) ) { return id; }
    return 0;
}

void writePDGTranslationLine( int i, std::ostream & os  )
{
    // only write map entries
    static PDGtoPDTMap const & hmap = getPDGtoPDTMap();

    PDGtoPDTMap::const_iterator const cit = hmap.find( i );
    // found it in the map
    if ( cit != hmap.end() ) { 
        int id = cit->second;
	os << " PDG number: " ;
	os.width(10);
	os << i << "  HepPID number: " ;
	os.width(11);
	os << id << "  " << particleName(id) << std::endl;
	// check reverse translation
	int iback =  translatePDTtoPDGtable(id);
	if(iback != i) {
	    os << " WARNING: " << id << " translates back to " 
	       << iback << " not to " << i << std::endl;
	}
    }
    return;
}  // writePDGTranslationLine()

void  writePDGTranslation( std::ostream & os )
{
    writeVersion( os );
    os << "     HepPID Particle List" << std::endl;
    os << std::endl;

    int id, i, j, q1, q2, q3, l, m;
    // special cases
    for( id=1; id<102; ++id) {
        writePDGTranslationLine(  id, os );
        writePDGTranslationLine( -id, os );
    }
    for( i=11; i<1000; ++i) {
        id = i*10;
        writePDGTranslationLine(  id, os );
        writePDGTranslationLine( -id, os );
    }
    // diquark
    for( q2=1; q2<7; ++q2) {
	for( q1=1; q1<7; ++q1) {
            for( j=1; j<4; ++j) {
        	id = 1000*q2+100*q1+j;
        	writePDGTranslationLine(  id, os );
        	writePDGTranslationLine( -id, os );
	    }
	}
    }
    // mesons
    for( q2=1; q2<10; ++q2) {
	for( q1=1; q1<10; ++q1) {
	    for( j=1; j<10; ++j) {
		for( m=0; m<9; ++m) {
		    for( l=0; l<10; ++l) {
			id = 100000*m+10000*l+100*q2+10*q1+j;
        		writePDGTranslationLine(  id, os );
        		writePDGTranslationLine( -id, os );
		    }
		}
	    }
	}
    }
    // baryons
    for( q3=1; q3<10; ++q3) {
	for( q2=1; q2<10; ++q2) {
	    for( q1=1; q1<10; ++q1) {
		for( j=1; j<10; ++j) {
		    for( m=0; m<9; ++m) {
			id = 10000*m+1000*q3+100*q2+10*q1+j;
        		writePDGTranslationLine(  id, os );
        		writePDGTranslationLine( -id, os );
		    }
		}
	    }
	}
    }
    // ions
    for( i=1; i<3; ++i) {
	for( m=1; m<3; ++m) {
	    id = 1000+10*i+m;
            writePDGTranslationLine(  id, os );
            writePDGTranslationLine( -id, os );
	}
    }

    return;
}  // writePDGTranslation()

}	// HepPID
