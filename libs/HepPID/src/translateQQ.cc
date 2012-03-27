// ------------------------------------
//
// translateQQ.cc
// Author: Lynn Garren
//
// translate an ID number to or from the standard numbering scheme and QQ
// use static maps
//
//  The maps are initialized if and only if the public functions are called.
//  Because the maps are static, the initialization happens only once.
//
//  The user NEVER calls QQPDTMapInit()
//  We use a data table (struct SList) so that compile time is not impacted.
//
//  public functions:
//        int translateQQtoPDT( const int id )
//        int translatePDTtoQQ( const int id )
//        QQPDTMap const & getQQPDTMap()
//        PDTQQMap const & getPDTQQMap()
//        int translateQQbar( const int id )
//        int translateInverseQQbar( const int id )
//        QQbarMap const & getQQbarMap()
//        InverseQQbarMap const & getInverseQQbarMap()
//
// ------------------------------------

#include <map>
#include <utility>	// make_pair

#include "HepPID/Version.hh"
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"

namespace HepPID {

 typedef  std::map< int, int >  QQPDTMap;
 typedef  std::map< int, int >  PDTQQMap;
 typedef  std::map< int, int >  QQbarMap;
 typedef  std::map< int, int >  InverseQQbarMap;

namespace {	// QQPDTMapInit is private

 QQPDTMap const & getQQPDTMap();
 PDTQQMap const & getPDTQQMap();
 QQbarMap const & getQQbarMap();
 InverseQQbarMap const & getInverseQQbarMap();

QQPDTMap const & QQPDTMapInit()
{

  static QQPDTMap  m;

  static const struct {
      int hid;	// Pythia
      int pid;	// PDT
  } SList[] = {
    {    -13,         21 },
    {    -12,         -6 },
    {    -11,         -5 },
    {    -10,         -4 },
    {     -9,         -3 },
    {     -8,         -1 },
    {     -7,         -2 },
    {     -6,          6 },
    {     -5,          5 },
    {     -4,          4 },
    {     -3,          3 },
    {     -2,          1 },
    {     -1,          2 },
    {      0,       10022},
    {      1,         22 },
    {      2,         23 },
    {      3,         24 },
    {      4,        -24 },
    {      5,         82 },
    {      7,         11 },
    {      8,        -11 },
    {      9,         12 },
    {     10,        -12 },
    {     11,         13 },
    {     12,        -13 },
    {     13,         14 },
    {     14,        -14 },
    {     15,         15 },
    {     16,        -15 },
    {     17,         16 },
    {     18,        -16 },
    {     19,      20313 },
    {     20,     -20313 },
    {     21,        211 },
    {     22,       -211 },
    {     23,        321 },
    {     24,       -321 },
    {     25,        311 },
    {     26,       -311 },
    {     27,        421 },
    {     28,       -421 },
    {     29,        411 },
    {     30,       -411 },
    {     31,        431 },
    {     32,       -431 },
    {     33,       -521 },
    {     34,        521 },
    {     35,       -511 },
    {     36,        511 },
    {     37,       -531 },
    {     38,        531 },
    {     39,       -541 },
    {     40,        541 },
    {     41,        621 },
    {     42,       -621 },
    {     43,        611 },
    {     44,       -611 },
    {     45,        631 },
    {     46,       -631 },
    {     47,        641 },
    {     48,       -641 },
    {     49,        651 },
    {     50,       -651 },
    {     51,        111 },
    {     52,        221 },
    {     53,        331 },
    {     54,        441 },
    {     55,        551 },
    {     56,        661 },
    {     57,        310 },
    {     58,        130 },
    {     59,      10313 },
    {     60,     -10313 },
    {     61,        213 },
    {     62,       -213 },
    {     63,        323 },
    {     64,       -323 },
    {     65,        313 },
    {     66,       -313 },
    {     67,        423 },
    {     68,       -423 },
    {     69,        413 },
    {     70,       -413 },
    {     71,        433 },
    {     72,       -433 },
    {     73,       -523 },
    {     74,        523 },
    {     75,       -513 },
    {     76,        513 },
    {     77,       -533 },
    {     78,        533 },
    {     79,       -543 },
    {     80,        543 },
    {     81,        623 },
    {     82,       -623 },
    {     83,        613 },
    {     84,       -613 },
    {     85,        633 },
    {     86,       -633 },
    {     87,        643 },
    {     88,       -643 },
    {     89,        653 },
    {     90,       -653 },
    {     91,        113 },
    {     92,        223 },
    {     93,        333 },
    {     94,        443 },
    {     95,        553 },
    {     96,        663 },
    {     97,     100553 },
    {     98,     200553 },
    {     99,     300553 },
    {    100,      10551 },
    {    101,      20553 },
    {    102,        555 },
    {    103,     110551 },
    {    104,     120553 },
    {    105,     100555 },
    {    106,      30113 },
    {    107,      20213 },
    {    108,      20113 },
    {    109,     -20213 },
    {    110,      10441 },
    {    111,      20443 },
    {    112,        445 },
    {    121,       3122 },
    {    122,      -3122 },
    {    123,       4122 },
    {    124,      -4122 },
    {    125,       4232 },
    {    126,      -4232 },
    {    127,       4132 },
    {    128,      -4132 },
    {    129,       3212 },
    {    130,      -3212 },
    {    131,       4212 },
    {    132,      -4212 },
    {    133,       4322 },
    {    134,      -4322 },
    {    135,       4312 },
    {    136,      -4312 },
    {    137,       2212 },
    {    138,      -2212 },
    {    139,       3222 },
    {    140,      -3222 },
    {    141,       4222 },
    {    142,      -4222 },
    {    143,       2112 },
    {    144,      -2112 },
    {    145,       3112 },
    {    146,      -3112 },
    {    147,       4112 },
    {    148,      -4112 },
    {    149,       3322 },
    {    150,      -3322 },
    {    151,       3312 },
    {    152,      -3312 },
    {    153,       4332 },
    {    154,      -4332 },
    {    155,       4422 },
    {    156,      -4422 },
    {    157,       4412 },
    {    158,      -4412 },
    {    159,       4432 },
    {    160,      -4432 },
    {    161,       3214 },
    {    162,      -3214 },
    {    163,       4214 },
    {    164,      -4214 },
    {    165,       4324 },
    {    166,      -4324 },
    {    167,       4314 },
    {    168,      -4314 },
    {    169,       2214 },
    {    170,      -2214 },
    {    171,       3224 },
    {    172,      -3224 },
    {    173,       4224 },
    {    174,      -4224 },
    {    175,       2114 },
    {    176,      -2114 },
    {    177,       3114 },
    {    178,      -3114 },
    {    179,       4114 },
    {    180,      -4114 },
    {    181,       3324 },
    {    182,      -3324 },
    {    183,       3314 },
    {    184,      -3314 },
    {    185,       4334 },
    {    186,      -4334 },
    {    187,       4424 },
    {    188,      -4424 },
    {    189,       4414 },
    {    190,      -4414 },
    {    191,       4434 },
    {    192,      -4434 },
    {    193,       2224 },
    {    194,      -2224 },
    {    195,       1114 },
    {    196,      -1114 },
    {    197,       3334 },
    {    198,      -3334 },
    {    199,       4444 },
    {    200,      -4444 },
    {    201,      10323 },
    {    202,     -10323 },
    {    203,      20323 },
    {    204,     -20323 },
    {    205,      30213 },
    {    206,     -30213 },
    {    207,         84 },
    {    208,        -84 },
    {    209,         85 },
    {    210,        -85 },
    {    211,      30443 },
    {    212,    9000443 },
    {    213,    9010443 },
    {    214,    9020443 },
    {    215,      10443 },
    {    216,    9000553 },
    {    217,    9010553 },
    {    218,      10553 },
    {    219,     100443 },
    {    220,    9020553 },
    {    221,      10411 },
    {    222,      20413 },
    {    223,      10413 },
    {    224,        415 },
    {    225,     -10411 },
    {    226,     -20413 },
    {    227,     -10413 },
    {    228,       -415 },
    {    229,      10421 },
    {    230,      20423 },
    {    231,      10423 },
    {    232,        425 },
    {    233,     -10421 },
    {    234,     -20423 },
    {    235,     -10423 },
    {    236,       -425 },
    {    237,      10431 },
    {    238,      20433 },
    {    239,      10433 },
    {    240,        435 },
    {    241,     -10431 },
    {    242,     -20433 },
    {    243,     -10433 },
    {    244,       -435 },
    {    251,    9000111 },
    {    252,    9000211 },
    {    253,   -9000211 },
    {    254,        115 },
    {    255,        215 },
    {    256,       -215 },
    {    257,    9010221 },
    {    258,      10221 },
    {    259,      20223 },
    {    260,      20333 },
    {    261,        225 },
    {    262,        335 },
    {    263,      10223 },
    {    264,      10333 },
    {    265,      10113 },
    {    266,      10213 },
    {    267,     -10213 },
    {    268,      10311 },
    {    269,     -10311 },
    {    270,      10321 },
    {    271,     -10321 },
    {    272,        315 },
    {    273,       -315 },
    {    274,        325 },
    {    275,       -325 },
    {    276,         86 },
    {    277,        -86 },
    {    278,        317 },
    {    279,       -317 },
    {    280,        327 },
    {    281,       -327 },
    {    291,         87 },
    {    292,        -87 },
    {    293,         88 },
    {    294,        -88 },
    {    295,         89 },
    {    296,        -89 },
    {    297,         90 },
    {    298,        -90 },
    {    401,       5122 },
    {    402,      -5122 },
    {    403,       5232 },
    {    404,      -5232 },
    {    405,       5132 },
    {    406,      -5132 },
    {    407,       5242 },
    {    408,      -5242 },
    {    409,       5142 },
    {    410,      -5142 },
    {    411,       5342 },
    {    412,      -5342 },
    {    413,       5212 },
    {    414,      -5212 },
    {    415,       5322 },
    {    416,      -5322 },
    {    417,       5312 },
    {    418,      -5312 },
    {    419,       5422 },
    {    420,      -5422 },
    {    421,       5412 },
    {    422,      -5412 },
    {    423,       5432 },
    {    424,      -5432 },
    {    425,       5222 },
    {    426,      -5222 },
    {    427,       5112 },
    {    428,      -5112 },
    {    429,       5332 },
    {    430,      -5332 },
    {    431,       5442 },
    {    432,      -5442 },
    {    433,       5522 },
    {    434,      -5522 },
    {    435,       5512 },
    {    436,      -5512 },
    {    437,       5532 },
    {    438,      -5532 },
    {    439,       5542 },
    {    440,      -5542 },
    {    441,       5214 },
    {    442,      -5214 },
    {    443,       5324 },
    {    444,      -5324 },
    {    445,       5314 },
    {    446,      -5314 },
    {    447,       5424 },
    {    448,      -5424 },
    {    449,       5414 },
    {    450,      -5414 },
    {    451,       5434 },
    {    452,      -5434 },
    {    453,       5224 },
    {    454,      -5224 },
    {    455,       5114 },
    {    456,      -5114 },
    {    457,       5334 },
    {    458,      -5334 },
    {    459,       5444 },
    {    460,      -5444 },
    {    461,       5524 },
    {    462,      -5524 },
    {    463,       5514 },
    {    464,      -5514 },
    {    465,       5534 },
    {    466,      -5534 },
    {    467,       5544 },
    {    468,      -5544 },
    {    469,       5554 },
    {    470,      -5554 },
    {    471,      10521 },
    {    472,      20523 },
    {    473,      10523 },
    {    474,        525 },
    {    475,     -10521 },
    {    476,     -20523 },
    {    477,     -10523 },
    {    478,       -525 },
    {    479,      10511 },
    {    480,      20513 },
    {    481,      10513 },
    {    482,        515 },
    {    483,     -10511 },
    {    484,     -20513 },
    {    485,     -10513 },
    {    486,       -515 },
    {    487,      10531 },
    {    488,      20533 },
    {    489,      10533 },
    {    490,        535 },
    {    491,     -10531 },
    {    492,     -20533 },
    {    493,     -10533 },
    {    494,       -535 },
    {    495,         92 },
    {    496,        -92 },
    {    497,         93 },
    {    498,        -93 }
  };

  int listSize = sizeof(SList)/sizeof(SList[0]);
  for( int k=0; k!=listSize; ++k) {
      m.insert( std::make_pair( SList[k].hid, SList[k].pid) );
  }
  return m;
}  // QQPDTMapInit()

// we need a separate map for the QQ quark pair pseudo-particles
// use diquark particle ID numbers
QQbarMap const & QQbarMapInit()
{
  static QQbarMap  m;

  static const struct {
      int hid;	// Pythia
      int pid;	// PDT
  } SList[] = {
    {      1, 2203 },
    {      2, 2101 },
    {      3, 3203 },
    {      4, 4203 },
    {      5, 5203 },
    {      6, 6203 },
    {      7, 2103 },
    {      8, 1103 },
    {      9, 3103 },
    {     10, 4103 },
    {     11, 5103 },
    {     12, 6103 },
    {     13, 3201 },
    {     14, 3101 },
    {     15, 3303 },
    {     16, 4303 },
    {     17, 5303 },
    {     18, 6303 },
    {     19, 4201 },
    {     20, 4101 },
    {     21, 4301 },
    {     22, 4403 },
    {     23, 5403 },
    {     24, 6403 },
    {     25, 5201 },
    {     26, 5101 },
    {     27, 5301 },
    {     28, 5401 },
    {     29, 5503 },
    {     30, 6503 },
    {     31, 6201 },
    {     32, 6101 },
    {     33, 6301 },
    {     34, 6401 },
    {     35, 6501 },
    {     36, 6603 },
    {     37,   81 }
  };

  int listSize = sizeof(SList)/sizeof(SList[0]);
  for( int k=0; k!=listSize; ++k) {
      m.insert( std::make_pair( SList[k].hid, SList[k].pid) );
  }
  return m;
}  // QQbarMapInit()

PDTQQMap const & PDTQQMapInit()
{
    static PDTQQMap m;
    static QQPDTMap const & hmap = getQQPDTMap();
    
    for(QQPDTMap::const_iterator cit=hmap.begin(), mend=hmap.end(); cit!=mend; ++cit ) {
	m.insert( std::make_pair( cit->second, cit->first ));
    }
    return m;
}

InverseQQbarMap const & InverseQQbarMapInit()
{
    static InverseQQbarMap m;
    static QQbarMap const & hmap = getQQbarMap();
    
    for(QQbarMap::const_iterator cit=hmap.begin(), mend=hmap.end(); cit!=mend; ++cit ) {
	m.insert( std::make_pair( cit->second, cit->first ));
    }
    return m;
}
  
//
// getQQPDTMap is the ONLY function allowed to call QQPDTMapInit
//
QQPDTMap const & getQQPDTMap()
{
  static QQPDTMap const & hmap = QQPDTMapInit();
  return hmap;
}  // getQQPDTMap()

//
// getPDTQQMap is the ONLY function allowed to call QQPDTMapInit
//
PDTQQMap const & getPDTQQMap()
{
  static PDTQQMap const & hmap = PDTQQMapInit();
  return hmap;
}  // getPDTQQMap()
//
// getQQbarMap is the ONLY function allowed to call QQbarMapInit
//
QQbarMap const & getQQbarMap()
{
  static QQbarMap const & hmap = QQbarMapInit();
  return hmap;
}  // getQQbarMap()

//
// getInverseQQbarMap is the ONLY function allowed to call QQbarMapInit
//
InverseQQbarMap const & getInverseQQbarMap()
{
  static InverseQQbarMap const & hmap = InverseQQbarMapInit();
  return hmap;
}  // getInverseQQbarMap()

} // unnamed namespace
  
int translateQQbar( const int id )
{
    static QQbarMap const & hmap = getQQbarMap();

    QQbarMap::const_iterator const cit = hmap.find( id );
    // found it in the map
    if ( cit != hmap.end() ) { return cit->second; }
    // for QQ, you can only use the map
    return 0;
}

int translateInverseQQbar( const int id )
{
    static InverseQQbarMap const & pmap = getInverseQQbarMap();

    InverseQQbarMap::const_iterator const cit = pmap.find( id );
    // found it in the map
    if ( cit != pmap.end() ) { return cit->second; }
    // for QQ, you can only use the map
    return 0;
}

int translateQQtoPDT( const int id )
{
    static QQPDTMap const & hmap = getQQPDTMap();

    QQPDTMap::const_iterator const cit = hmap.find( id );
    // found it in the map
    if ( cit != hmap.end() ) { return cit->second; }
    // for QQ, you can only use the map
    return 0;
}

int translatePDTtoQQ( const int id )
{
    static PDTQQMap const & pmap = getPDTQQMap();

    PDTQQMap::const_iterator const cit = pmap.find( id );
    // found it in the map
    if ( cit != pmap.end() ) { return cit->second; }
    // for QQ, you can only use the map
    return 0;
}

void  writeQQTranslation ( std::ostream & os )
{
    int id, iq, iback;
    writeVersion( os );
    os << "     HepPID Particle List" << std::endl;
    os << std::endl;

    // quark pairs have overlapping QQ ID numbers
    for( iq=1; iq<40; ++iq) {
	id = translateQQbar( iq );
	if ( id != 0 ) {
	    os << " QQ number: " ;
	    os.width(10);
	    os << iq << "  HepPID number: " ;
	    os.width(10);
	    os << id << "  " << particleName(id) << std::endl;
	    // check reverse translation
	    iback =  translateInverseQQbar(id);
	    if(iback != iq) {
	        os << " WARNING: " << id << " translates back to " 
		   << iback << " not to " << iq << std::endl;
	    }
	}
    }
    // regular QQ particles
    for( iq=-13; iq<501; ++iq) {
	id = translateQQtoPDT( iq );
	if ( id != 0 ) {
	    os << " QQ number: " ;
	    os.width(10);
	    os << iq << "  HepPID number: " ;
	    os.width(10);
	    os << id << "  " << particleName(id) << std::endl;
	    // check reverse translation
	    iback =  translatePDTtoQQ(id);
	    if(iback != iq) {
	        os << " WARNING: " << id << " translates back to " 
		   << iback << " not to " << iq << std::endl;
	    }
	}
    }
    return;
}  // writeQQTranslation()

}	// HepPID
