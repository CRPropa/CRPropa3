// ------------------------------------
//
// translateHerwig.cc
// Author: Lynn Garren
//
// ..convert from HERWIG numbering scheme to PDT numbering scheme
//  use a static map for both translateHerwigtoPDT and translatePDTtoHerwig
//
//  The maps are initialized if and only if the public functions are called.
//  Because the maps are static, the initialization happens only once.
//
//  The user NEVER calls HerwigPDTMapInit()
//  We use a data table (struct SList) so that compile time is not impacted.
//
//  public functions:
//        int translateHerwigtoPDT( const int id )
//        int translatePDTtoHerwig( const int id )
//        HerwigPDTMap const & getHerwigPDTMap()
//        PDTHerwigMap const & getPDTHerwigMap()
//
// ------------------------------------

#include <map>
#include <utility>	// make_pair

#include "HepPID/Version.hh"
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"

namespace HepPID {

 typedef  std::map< int, int >  HerwigPDTMap;
 typedef  std::map< int, int >  PDTHerwigMap;

namespace {	// HerwigPDTMapInit is private

 HerwigPDTMap const & getHerwigPDTMap();
 PDTHerwigMap const & getPDTHerwigMap();

HerwigPDTMap const & HerwigPDTMapInit()
{

  static HerwigPDTMap  m;

  static const struct {
      int hid;	// Herwig
      int pid;	// PDT
  } SList[] = {
     {            1,            1 },
     {            2,            2 },
     {            3,            3 },
     {            4,            4 },
     {            5,            5 },
     {            6,            6 },
     {            7,            7 },
     {            8,            8 },
     {           -1,           -1 },
     {           -2,           -2 },
     {           -3,           -3 },
     {           -4,           -4 },
     {           -5,           -5 },
     {           -6,           -6 },
     {           -7,           -7 },
     {           -8,           -8 },
     {           11,           11 },
     {           12,           12 },
     {           13,           13 },
     {           14,           14 },
     {           15,           15 },
     {           16,           16 },
     {          -11,          -11 },
     {          -12,          -12 },
     {          -13,          -13 },
     {          -14,          -14 },
     {          -15,          -15 },
     {          -16,          -16 },
     {           21,           21 },
     {           22,           22 },
     {           23,           23 },
     {           24,           24 },
     {          -24,          -24 },
     {           25,           25 },
     {           26,           51 },
     {           32,           32 },
     {           35,           35 },
     {           36,           36 },
     {           37,           37 },
     {          -37,          -37 },
     {           39,           39 },
     {           91,           91 },
     {           98,      9920022 },
     {           99,      9922212 },
     {          -99,     -9922212 },
     {          111,          111 },
     {          221,          221 },
     {          113,          113 },
     {          223,          223 },
     {          331,          331 },
     {          225,          225 },
     {        20113,        20113 },
     {        20223,        20223 },
     {          115,          115 },
     {         -211,         -211 },
     {         -213,         -213 },
     {       -20213,       -20213 },
     {         -215,         -215 },
     {         -321,         -321 },
     {         -323,         -323 },
     {       -20323,       -20323 },
     {         -325,         -325 },
     {          211,          211 },
     {          213,          213 },
     {        20213,        20213 },
     {          215,          215 },
     {         -311,         -311 },
     {         -313,         -313 },
     {       -20313,       -20313 },
     {         -315,         -315 },
     {          321,          321 },
     {          323,          323 },
     {        20323,        20323 },
     {          325,          325 },
     {          311,          311 },
     {          313,          313 },
     {        20313,        20313 },
     {          315,          315 },
     {          333,          333 },
     {        20333,        20333 },
     {          335,          335 },
     {          310,          310 },
     {          130,          130 },
     {        10111,        10111 },
     {        10211,        10211 },
     {       -10211,       -10211 },
     {         2212,         2212 },
     {         2214,         2214 },
     {         2112,         2112 },
     {         2114,         2114 },
     {         1114,         1114 },
     {         3122,         3122 },
     {         3212,         3212 },
     {         3214,         3214 },
     {         3112,         3112 },
     {         3114,         3114 },
     {         3312,         3312 },
     {         3314,         3314 },
     {         2224,         2224 },
     {         3222,         3222 },
     {         3224,         3224 },
     {         3322,         3322 },
     {         3324,         3324 },
     {         3334,         3334 },
     {        -2212,        -2212 },
     {        -2214,        -2214 },
     {        -2112,        -2112 },
     {        -2114,        -2114 },
     {        -1114,        -1114 },
     {        -3122,        -3122 },
     {        -3212,        -3212 },
     {        -3214,        -3214 },
     {        -3112,        -3112 },
     {        -3114,        -3114 },
     {        -3312,        -3312 },
     {        -3314,        -3314 },
     {        -2224,        -2224 },
     {        -3222,        -3222 },
     {        -3224,        -3224 },
     {        -3322,        -3322 },
     {        -3324,        -3324 },
     {        -3334,        -3334 },
     {         2203,         2203 },
     {         2101,         2101 },
     {         1103,         1103 },
     {         3201,         3201 },
     {         3101,         3101 },
     {         3303,         3303 },
     {        -2203,        -2203 },
     {        -2101,        -2101 },
     {        -1103,        -1103 },
     {        -3201,        -3201 },
     {        -3101,        -3101 },
     {        -3303,        -3303 },
     {          411,          411 },
     {          413,          413 },
     {        20413,        20413 },
     {          415,          415 },
     {          421,          421 },
     {          423,          423 },
     {        20423,        20423 },
     {          425,          425 },
     {          431,          431 },
     {          433,          433 },
     {        20433,        20433 },
     {          435,          435 },
     {         4222,         4222 },
     {         4224,         4224 },
     {         4122,         4122 },
     {         4212,         4212 },
     {         4214,         4214 },
     {         4112,         4112 },
     {         4114,         4114 },
     {         4232,         4232 },
     {         4322,         4322 },
     {         4324,         4324 },
     {         4132,         4132 },
     {         4312,         4312 },
     {         4314,         4314 },
     {         4332,         4332 },
     {         4334,         4334 },
     {          441,          441 },
     {          443,          443 },
     {        10441,        10441 },
     {       100443,       100443 },
     {        30443,        30443 },
     {         -411,         -411 },
     {         -413,         -413 },
     {       -20413,       -20413 },
     {         -415,         -415 },
     {         -421,         -421 },
     {         -423,         -423 },
     {       -20423,       -20423 },
     {         -425,         -425 },
     {         -431,         -431 },
     {         -433,         -433 },
     {       -20433,       -20433 },
     {         -435,         -435 },
     {        -4222,        -4222 },
     {        -4224,        -4224 },
     {        -4122,        -4122 },
     {        -4212,        -4212 },
     {        -4214,        -4214 },
     {        -4112,        -4112 },
     {        -4114,        -4114 },
     {        -4232,        -4232 },
     {        -4322,        -4322 },
     {        -4324,        -4324 },
     {        -4132,        -4132 },
     {        -4312,        -4312 },
     {        -4314,        -4314 },
     {        -4332,        -4332 },
     {        -4334,        -4334 },
     {         -511,         -511 },
     {         -521,         -521 },
     {         -531,         -531 },
     {         5222,         5222 },
     {         5122,         5122 },
     {         5112,         5112 },
     {         5232,         5232 },
     {         5132,         5132 },
     {         5332,         5332 },
     {         -541,         -541 },
     {          553,          553 },
     {         -651,         -651 },
     {          611,          611 },
     {          621,          621 },
     {          631,          631 },
     {         6222,         6222 },
     {         6122,         6122 },
     {         6112,         6112 },
     {         6232,         6232 },
     {         6132,         6132 },
     {         6332,         6332 },
     {          641,          641 },
     {          651,          651 },
     {          663,          663 },
     {          511,          511 },
     {          521,          521 },
     {          531,          531 },
     {        -5222,        -5222 },
     {        -5122,        -5122 },
     {        -5112,        -5112 },
     {        -5232,        -5232 },
     {        -5132,        -5132 },
     {        -5332,        -5332 },
     {          541,          541 },
     {         -611,         -611 },
     {         -621,         -621 },
     {         -631,         -631 },
     {        -6222,        -6222 },
     {        -6122,        -6122 },
     {        -6112,        -6112 },
     {        -6232,        -6232 },
     {        -6132,        -6132 },
     {        -6332,        -6332 },
     {         -641,         -641 },
     {         -513,         -513 },
     {         -523,         -523 },
     {         -533,         -533 },
     {       -20513,       -20513 },
     {       -20523,       -20523 },
     {       -20533,       -20533 },
     {         -515,         -515 },
     {         -525,         -525 },
     {         -535,         -535 },
     {          513,          513 },
     {          523,          523 },
     {          533,          533 },
     {        20513,        20513 },
     {        20523,        20523 },
     {        20533,        20533 },
     {          515,          515 },
     {          525,          525 },
     {          535,          535 },
     {        10113,        10113 },
     {        10213,        10213 },
     {       -10213,       -10213 },
     {        10223,        10223 },
     {        10333,        10333 },
     {      9000111,      9000111 },
     {      9000211,      9000211 },
     {     -9000211,     -9000211 },
     {      9010221,      9010221 },
     {        10221,        10221 },
     {          543,          543 },
     {         -543,         -543 },
     {        20543,        20543 },
     {       -20543,       -20543 },
     {          545,          545 },
     {         -545,         -545 },
     {        10443,        10443 },
     {        20443,        20443 },
     {          445,          445 },
     {          551,          551 },
     {        10553,        10553 },
     {        10551,        10551 },
     {        20553,        20553 },
     {          555,          555 },
     {        10313,        10313 },
     {        10323,        10323 },
     {       -10313,       -10313 },
     {       -10323,       -10323 },
     {        10413,        10413 },
     {        10423,        10423 },
     {        10433,        10433 },
     {       -10413,       -10413 },
     {       -10423,       -10423 },
     {       -10433,       -10433 },
     {        10513,        10513 },
     {        10523,        10523 },
     {        10533,        10533 },
     {        10543,        10543 },
     {       -10513,       -10513 },
     {       -10523,       -10523 },
     {       -10533,       -10533 },
     {       -10543,       -10543 },
     {        10321,        10321 },
     {        10311,        10311 },
     {       -10311,       -10311 },
     {       -10321,       -10321 },
     {        10411,        10411 },
     {        10421,        10421 },
     {        10431,        10431 },
     {       -10411,       -10411 },
     {       -10421,       -10421 },
     {       -10431,       -10431 },
     {        10511,        10511 },
     {        10521,        10521 },
     {        10531,        10531 },
     {        10541,        10541 },
     {       -10511,       -10511 },
     {       -10521,       -10521 },
     {       -10531,       -10531 },
     {       -10541,       -10541 },
     {         5114,         5114 },
     {         5212,         5212 },
     {         5214,         5214 },
     {         5224,         5224 },
     {         5322,         5322 },
     {         5324,         5324 },
     {         5312,         5312 },
     {         5314,         5314 },
     {         5334,         5334 },
     {        -5114,        -5114 },
     {        -5212,        -5212 },
     {        -5214,        -5214 },
     {        -5224,        -5224 },
     {        -5322,        -5322 },
     {        -5324,        -5324 },
     {        -5312,        -5312 },
     {        -5314,        -5314 },
     {        -5334,        -5334 },
     {        10325,        10325 },
     {        10315,        10315 },
     {       -10315,       -10315 },
     {       -10325,       -10325 },
     {        30323,        30323 },
     {        30313,        30313 },
     {       -30313,       -30313 },
     {       -30323,       -30323 },
     {        20325,        20325 },
     {        20315,        20315 },
     {       -20315,       -20315 },
     {       -20325,       -20325 },
     {          327,          327 },
     {          317,          317 },
     {         -317,         -317 },
     {         -327,         -327 },
     {        10215,        10215 },
     {        10115,        10115 },
     {       -10215,       -10215 },
     {        30213,        30213 },
     {        30113,        30113 },
     {       -30213,       -30213 },
     {          217,          217 },
     {          117,          117 },
     {         -217,         -217 },
     {       100553,       100553 },
     {       110551,       110551 },
     {       120553,       120553 },
     {       100555,       100555 },
     {       200553,       200553 },
     {       300553,       300553 },
     {          227,          227 },
     {          337,          337 },
     {        10225,        10225 },
     {        10335,        10335 },
     {        30223,        30223 },
     {      1000001,      1000001 },
     {      1000002,      1000002 },
     {      1000003,      1000003 },
     {      1000004,      1000004 },
     {      1000005,      1000005 },
     {      1000006,      1000006 },
     {     -1000001,     -1000001 },
     {     -1000002,     -1000002 },
     {     -1000003,     -1000003 },
     {     -1000004,     -1000004 },
     {     -1000005,     -1000005 },
     {     -1000006,     -1000006 },
     {      2000001,      2000001 },
     {      2000002,      2000002 },
     {      2000003,      2000003 },
     {      2000004,      2000004 },
     {      2000005,      2000005 },
     {      2000006,      2000006 },
     {     -2000001,     -2000001 },
     {     -2000002,     -2000002 },
     {     -2000003,     -2000003 },
     {     -2000004,     -2000004 },
     {     -2000005,     -2000005 },
     {     -2000006,     -2000006 },
     {      1000011,      1000011 },
     {      1000012,      1000012 },
     {      1000013,      1000013 },
     {      1000014,      1000014 },
     {      1000015,      1000015 },
     {      1000016,      1000016 },
     {     -1000011,     -1000011 },
     {     -1000012,     -1000012 },
     {     -1000013,     -1000013 },
     {     -1000014,     -1000014 },
     {     -1000015,     -1000015 },
     {     -1000016,     -1000016 },
     {      2000011,      2000011 },
     {      2000012,      2000012 },
     {      2000013,      2000013 },
     {      2000014,      2000014 },
     {      2000015,      2000015 },
     {      2000016,      2000016 },
     {     -2000011,     -2000011 },
     {     -2000012,     -2000012 },
     {     -2000013,     -2000013 },
     {     -2000014,     -2000014 },
     {     -2000015,     -2000015 },
     {     -2000016,     -2000016 },
     {      1000021,      1000021 },
     {      1000022,      1000022 },
     {      1000023,      1000023 },
     {      1000025,      1000025 },
     {      1000035,      1000035 },
     {      1000024,      1000024 },
     {      1000037,      1000037 },
     {     -1000024,     -1000024 },
     {     -1000037,     -1000037 },
     {      1000039,      1000039 }
  };

  int listSize = sizeof(SList)/sizeof(SList[0]);
  for( int k=0; k!=listSize; ++k) {
      m.insert( std::make_pair( SList[k].hid, SList[k].pid) );
  }
  return m;
}  // HerwigPDTMapInit()

PDTHerwigMap const & PDTHerwigMapInit()
{
    static PDTHerwigMap m;
    static HerwigPDTMap const & hmap = getHerwigPDTMap();
    
    for(HerwigPDTMap::const_iterator cit=hmap.begin(), mend=hmap.end(); cit!=mend; ++cit ) {
	m.insert( std::make_pair( cit->second, cit->first ));
    }
    return m;
}
  
//
// getHerwigPDTMap is the ONLY function allowed to call HerwigPDTMapInit
//
HerwigPDTMap const & getHerwigPDTMap()
{
  static HerwigPDTMap const & hmap = HerwigPDTMapInit();
  return hmap;
}  // getHerwigPDTMap()

//
// getPDTHerwigMap is the ONLY function allowed to call HerwigPDTMapInit
//
PDTHerwigMap const & getPDTHerwigMap()
{
  static PDTHerwigMap const & hmap = PDTHerwigMapInit();
  return hmap;
}  // getPDTHerwigMap()
 
} // unnamed namespace
  
int translateHerwigtoPDT( const int id )
{
    static HerwigPDTMap const & hmap = getHerwigPDTMap();
    
    HerwigPDTMap::const_iterator const cit = hmap.find( id );
    // found it in the map
    if ( cit != hmap.end() ) { return cit->second; }
    // check to see if someone has defined a valid particle type
    // that isn't in the map
    if( isValid(id) ) { return id; }
    return 0;
}

int translatePDTtoHerwig( const int id )
{
    static PDTHerwigMap const & pmap = getPDTHerwigMap();

    PDTHerwigMap::const_iterator const cit = pmap.find( id );
    // found it in the map
    if ( cit != pmap.end() ) { return cit->second; }
    // check to see if someone has defined a valid particle type
    // that isn't in the map
    if( isValid(id) ) { return id; }
    return 0;
}

void writeHerwigTranslationLine( int i, std::ostream & os  )
{
    // only write map entries
    static HerwigPDTMap const & hmap = getHerwigPDTMap();

    HerwigPDTMap::const_iterator const cit = hmap.find( i );
    // found it in the map
    if ( cit != hmap.end() ) { 
        int id = cit->second;
	os << " Herwig number: " ;
	os.width(10);
	os << i << "  HepPID number: " ;
	os.width(10);
	os << id << "  " << particleName(id) << std::endl;
	// check reverse translation
	int iback =  translatePDTtoHerwig(id);
	if(iback != i) {
	    os << " WARNING: " << id << " translates back to " 
	       << iback << " not to " << i << std::endl;
	}
    }
    return;
}  // writeHerwigTranslationLine()

void  writeHerwigTranslation( std::ostream & os )
{
    writeVersion( os );
    os << "     HepPID Particle List" << std::endl;
    os << std::endl;

    int id, i, j, q1, q2, q3, l, m, n;
    // special cases
    for( id=1; id<101; ++id) {
        writeHerwigTranslationLine(  id, os );
        writeHerwigTranslationLine( -id, os );
    }
    // SUSY
    for( n=1; n<3; ++n) {
        for( i=1; i<40; ++i) {
	     id = 1000000*n+i;
             writeHerwigTranslationLine(  id, os );
             writeHerwigTranslationLine( -id, os );
	}
    }
    // diquark
    for( i=11; i<100; ++i) {
        for( j=0; j<10; ++j) {
            id = 100*i+j;
            writeHerwigTranslationLine(  id, os );
            writeHerwigTranslationLine( -id, os );
	}
    }
    // mesons
    for( q2=1; q2<10; ++q2) {
	for( q1=1; q1<10; ++q1) {
	    for( j=1; j<10; ++j) {
		for( m=0; m<9; ++m) {
		    for( l=0; l<10; ++l) {
			id = 100000*m+10000*l+100*q2+10*q1+j;
        		writeHerwigTranslationLine(  id, os );
        		writeHerwigTranslationLine( -id, os );
			id = 9000000+100000*m+10000*l+100*q2+10*q1+j;
        		writeHerwigTranslationLine(  id, os );
        		writeHerwigTranslationLine( -id, os );
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
		    id = 1000*q3+100*q2+10*q1+j;
        	    writeHerwigTranslationLine(  id, os );
        	    writeHerwigTranslationLine( -id, os );
		}
	    }
	}
    }
    // pentaquarks
    for( l=1; l<9; ++l ) {
        for ( m=1; m<9; ++m ) {
	    for( q3=1; q3<9; ++q3) {
		for( q2=1; q2<9; ++q2) {
		    for( q1=1; q1<9; ++q1) {
			id = 9*1000000+l*100000+m*10000+1000*q3+100*q2+10*q1+2;
        		writeHerwigTranslationLine(  id, os );
        		writeHerwigTranslationLine( -id, os );
		    }
		}
	    }
        }
    }
    return;
}  // writeHerwigTranslation()

}	// HepPID
