// ----------------------------------------------------------------------
//
// translateIsajet.cc
// Author: Lynn Garren
//
// translate an ID number to or from the standard numbering scheme and Isajet
// use static maps
//
// Isajet uses a different numbering scheme
// Private methods will attempt to convert mesons and baryons not in the map
//
//  The maps are initialized if and only if the public functions are called.
//  Because the maps are static, the initialization happens only once.
//
//  The user NEVER calls IsajetPDTMapInit()
//  We use a data table (struct SList) so that compile time is not impacted.
//
//  public functions:
//        int translateIsajettoPDT( const int id )
//        int translatePDTtoIsajet( const int id )
//        IsajetPDTMap const & getIsajetPDTMap()
//        PDTIsajetMap const & getPDTIsajetMap()
//
// ----------------------------------------------------------------------

#include <map>
#include <utility>	// make_pair

#include "HepPID/Version.hh"
#include "HepPID/ParticleIDTranslations.hh"
#include "HepPID/ParticleIDMethods.hh"
#include "HepPID/ParticleName.hh"

namespace HepPID {

 typedef  std::map< int, int >  IsajetPDTMap;
 typedef  std::map< int, int >  PDTIsajetMap;

namespace {	// IsajetPDTMapInit is private

 IsajetPDTMap const & getIsajetPDTMap();
 PDTIsajetMap const & getPDTIsajetMap();

IsajetPDTMap const & IsajetPDTMapInit()
{

  static IsajetPDTMap  m;

  static const struct {
      int id;	// Isajet
      int pid;	// PDT
  } SList[] = {
     {          1,            2 },
     {         -1,           -2 },
     {          2,            1 },
     {         -2,           -1 },
     {          3,            3 },
     {         -3,           -3 },
     {          4,            4 },
     {         -4,           -4 },
     {          5,            5 },
     {         -5,           -5 },
     {          6,            6 },
     {         -6,           -6 },
     {          7,            7 },
     {         -7,           -7 },
     {          8,            8 },
     {         -8,           -8 },
     {          9,           21 },
     {         10,           22 },
     {         11,           12 },
     {        -11,          -12 },
     {         12,           11 },
     {        -12,          -11 },
     {         13,           14 },
     {        -13,          -14 },
     {         14,           13 },
     {        -14,          -13 },
     {         15,           16 },
     {        -15,          -16 },
     {         16,           15 },
     {        -16,          -15 },
     {         20,          310 },
     {        -20,          130 },
     {         21,      1000002 },
     {        -21,     -1000002 },
     {         22,      1000001 },
     {        -22,     -1000001 },
     {         23,      1000003 },
     {        -23,     -1000003 },
     {         24,      1000004 },
     {        -24,     -1000004 },
     {         25,      1000005 },
     {        -25,     -1000005 },
     {         26,      1000006 },
     {        -26,     -1000006 },
     {         29,      1000021 },
     {         30,      1000022 },
     {         31,      1000012 },
     {        -31,     -1000012 },
     {         32,      1000011 },
     {        -32,     -1000011 },
     {         33,      1000014 },
     {        -33,     -1000014 },
     {         34,      1000013 },
     {        -34,     -1000013 },
     {         35,      1000016 },
     {        -35,     -1000016 },
     {         36,      1000015 },
     {        -36,     -1000015 },
     {         39,      1000024 },
     {        -39,     -1000024 },
     {         40,      1000023 },
     {         41,      2000002 },
     {        -41,     -2000002 },
     {         42,      2000001 },
     {        -42,     -2000001 },
     {         43,      2000003 },
     {        -43,     -2000003 },
     {         44,      2000004 },
     {        -44,     -2000004 },
     {         45,      2000005 },
     {        -45,     -2000005 },
     {         46,      2000006 },
     {        -46,     -2000006 },
     {         49,      1000037 },
     {        -49,     -1000037 },
     {         50,      1000025 },
     {         51,      2000012 },
     {        -51,     -2000012 },
     {         52,      2000011 },
     {        -52,     -2000011 },
     {         53,      2000014 },
     {        -53,     -2000014 },
     {         54,      2000013 },
     {        -54,     -2000013 },
     {         55,      2000016 },
     {        -55,     -2000016 },
     {         56,      2000015 },
     {        -56,     -2000015 },
     {         60,      1000035 },
     {         80,           24 },
     {        -80,          -24 },
     {         81,           25 },
     {         82,           51 },
     {         83,           35 },
     {         84,           36 },
     {         85,           55 },
     {        -85,          -55 },
     {         86,           37 },
     {        -86,          -37 },
     {         87,           53 },
     {        -87,          -53 },
     {         88,           52 },
     {        -88,          -52 },
     {         89,           54 },
     {        -89,          -54 },
     {         90,           23 },
     {         91,      1000039 },
     {         92,           39 },
     {        110,          111 },
     {        111,          113 },
     {        112,          225 },
     {        120,          211 },
     {       -120,         -211 },
     {        121,          213 },
     {       -121,         -213 },
     {        130,          321 },
     {       -130,         -321 },
     {        131,          323 },
     {       -131,         -323 },
     {        132,          325 },
     {       -132,         -325 },
     {        140,         -421 },
     {       -140,          421 },
     {        141,         -423 },
     {       -141,          423 },
     {        150,          521 },
     {       -150,         -521 },
     {        151,          523 },
     {       -151,         -523 },
     {        160,         -621 },
     {       -160,          621 },
     {        161,         -623 },
     {       -161,          623 },
     {        170,          721 },
     {       -170,         -721 },
     {        171,          723 },
     {       -171,         -723 },
     {        180,          821 },
     {       -180,         -821 },
     {        181,          823 },
     {       -181,         -823 },
     {        220,          221 },
     {        221,          223 },
     {        230,          311 },
     {       -230,         -311 },
     {        231,          313 },
     {       -231,         -313 },
     {        232,          315 },
     {       -232,         -315 },
     {        240,         -411 },
     {       -240,          411 },
     {        241,         -413 },
     {       -241,          413 },
     {        250,          511 },
     {       -250,         -511 },
     {        251,          513 },
     {       -251,         -513 },
     {        260,         -611 },
     {       -260,          611 },
     {        261,         -613 },
     {       -261,          613 },
     {        270,          711 },
     {       -270,         -711 },
     {        271,          713 },
     {       -271,         -713 },
     {        280,          811 },
     {       -280,         -811 },
     {        281,          813 },
     {       -281,         -813 },
     {        330,          331 },
     {        331,          333 },
     {        340,         -431 },
     {       -340,          431 },
     {        341,         -433 },
     {       -341,          433 },
     {        350,          531 },
     {       -350,         -531 },
     {        351,          533 },
     {       -351,         -533 },
     {        360,         -631 },
     {       -360,          631 },
     {        361,         -633 },
     {       -361,          633 },
     {        370,          731 },
     {       -370,         -731 },
     {        371,          733 },
     {       -371,         -733 },
     {        380,          831 },
     {       -380,         -831 },
     {        381,          833 },
     {       -381,         -833 },
     {        440,          441 },
     {        441,          443 },
     {        450,          541 },
     {       -450,         -541 },
     {        451,          543 },
     {       -451,         -543 },
     {        460,          641 },
     {       -460,         -641 },
     {        461,          643 },
     {       -461,         -643 },
     {        470,          741 },
     {       -470,         -741 },
     {        471,          743 },
     {       -471,         -743 },
     {        480,          841 },
     {       -480,         -841 },
     {        481,          843 },
     {       -481,         -843 },
     {        550,          551 },
     {        551,          553 },
     {        560,         -651 },
     {       -560,          651 },
     {        561,         -653 },
     {       -561,          653 },
     {        570,          751 },
     {       -570,         -751 },
     {        571,          753 },
     {       -571,         -753 },
     {        580,          851 },
     {       -580,         -851 },
     {        581,          853 },
     {       -581,         -853 },
     {        660,          661 },
     {        661,          663 },
     {        670,          761 },
     {       -670,         -761 },
     {        671,          763 },
     {       -671,         -763 },
     {        680,          861 },
     {       -680,         -861 },
     {        681,          863 },
     {       -681,         -863 },
     {        770,          771 },
     {        771,          773 },
     {        780,          871 },
     {       -780,         -871 },
     {        781,          873 },
     {       -781,         -873 },
     {        880,          881 },
     {        881,          883 },
     {       1100,         2203 },
     {      -1100,        -2203 },
     {       1111,         2224 },
     {      -1111,        -2224 },
     {       1120,         2212 },
     {      -1120,        -2212 },
     {       1121,         2214 },
     {      -1121,        -2214 },
     {       1130,         3222 },
     {      -1130,        -3222 },
     {       1131,         3224 },
     {      -1131,        -3224 },
     {       1140,         4222 },
     {      -1140,        -4222 },
     {       1141,         4224 },
     {      -1141,        -4224 },
     {       1150,         5222 },
     {      -1150,        -5222 },
     {       1151,         5224 },
     {      -1151,        -5224 },
     {       1160,         6222 },
     {      -1160,        -6222 },
     {       1161,         6224 },
     {      -1161,        -6224 },
     {       1170,         7222 },
     {      -1170,        -7222 },
     {       1171,         7224 },
     {      -1171,        -7224 },
     {       1180,         8222 },
     {      -1180,        -8222 },
     {       1181,         8224 },
     {      -1181,        -8224 },
     {       1200,         2101 },
     {      -1200,        -2101 },
     {       1220,         2112 },
     {      -1220,        -2112 },
     {       1221,         2114 },
     {      -1221,        -2114 },
     {       1230,         3212 },
     {      -1230,        -3212 },
     {       1231,         3214 },
     {      -1231,        -3214 },
     {       1240,         4212 },
     {      -1240,        -4212 },
     {       1241,         4214 },
     {      -1241,        -4214 },
     {       1250,         5212 },
     {      -1250,        -5212 },
     {       1251,         5214 },
     {      -1251,        -5214 },
     {       1260,         6212 },
     {      -1260,        -6212 },
     {       1261,         6214 },
     {      -1261,        -6214 },
     {       1270,         7212 },
     {      -1270,        -7212 },
     {       1271,         7214 },
     {      -1271,        -7214 },
     {       1280,         8212 },
     {      -1280,        -8212 },
     {       1281,         8214 },
     {      -1281,        -8214 },
     {       1300,         3201 },
     {      -1300,        -3201 },
     {       1330,         3322 },
     {      -1330,        -3322 },
     {       1331,         3324 },
     {      -1331,        -3324 },
     {       1340,         4322 },
     {      -1340,        -4322 },
     {       1341,         4324 },
     {      -1341,        -4324 },
     {       1350,         5322 },
     {      -1350,        -5322 },
     {       1351,         5324 },
     {      -1351,        -5324 },
     {       1360,         6322 },
     {      -1360,        -6322 },
     {       1361,         6324 },
     {      -1361,        -6324 },
     {       1370,         7322 },
     {      -1370,        -7322 },
     {       1371,         7324 },
     {      -1371,        -7324 },
     {       1380,         8322 },
     {      -1380,        -8322 },
     {       1381,         8324 },
     {      -1381,        -8324 },
     {       1400,        -4201 },
     {      -1400,         4201 },
     {       1440,         4422 },
     {      -1440,        -4422 },
     {       1441,         4424 },
     {      -1441,        -4424 },
     {       1450,         5422 },
     {      -1450,        -5422 },
     {       1451,         5424 },
     {      -1451,        -5424 },
     {       1460,         6422 },
     {      -1460,        -6422 },
     {       1461,         6424 },
     {      -1461,        -6424 },
     {       1500,         5201 },
     {      -1500,        -5201 },
     {       1550,         5522 },
     {      -1550,        -5522 },
     {       1551,         5524 },
     {      -1551,        -5524 },
     {       1560,         6522 },
     {      -1560,        -6522 },
     {       1561,         6524 },
     {      -1561,        -6524 },
     {       1600,        -6201 },
     {      -1600,         6201 },
     {       1660,         6622 },
     {      -1660,        -6622 },
     {       1661,         6624 },
     {      -1661,        -6624 },
     {       2130,         3122 },
     {      -2130,        -3122 },
     {       2140,         4122 },
     {      -2140,        -4122 },
     {       2150,         5122 },
     {      -2150,        -5122 },
     {       2160,         6122 },
     {      -2160,        -6122 },
     {       2170,         7122 },
     {      -2170,        -7122 },
     {       2180,         8122 },
     {      -2180,        -8122 },
     {       2200,         1103 },
     {      -2200,        -1103 },
     {       2221,         1114 },
     {      -2221,        -1114 },
     {       2230,         3112 },
     {      -2230,        -3112 },
     {       2231,         3114 },
     {      -2231,        -3114 },
     {       2240,         4112 },
     {      -2240,        -4112 },
     {       2241,         4114 },
     {      -2241,        -4114 },
     {       2250,         5112 },
     {      -2250,        -5112 },
     {       2251,         5114 },
     {      -2251,        -5114 },
     {       2260,         6112 },
     {      -2260,        -6112 },
     {       2261,         6114 },
     {      -2261,        -6114 },
     {       2270,         7112 },
     {      -2270,        -7112 },
     {       2271,         7114 },
     {      -2271,        -7114 },
     {       2280,         8112 },
     {      -2280,        -8112 },
     {       2281,         8114 },
     {      -2281,        -8114 },
     {       2300,         3101 },
     {      -2300,        -3101 },
     {       2330,         3312 },
     {      -2330,        -3312 },
     {       2331,         3314 },
     {      -2331,        -3314 },
     {       2340,         4312 },
     {      -2340,        -4312 },
     {       2341,         4314 },
     {      -2341,        -4314 },
     {       2350,         5312 },
     {      -2350,        -5312 },
     {       2351,         5314 },
     {      -2351,        -5314 },
     {       2360,         6312 },
     {      -2360,        -6312 },
     {       2361,         6314 },
     {      -2361,        -6314 },
     {       2370,         7312 },
     {      -2370,        -7312 },
     {       2371,         7314 },
     {      -2371,        -7314 },
     {       2380,         8312 },
     {      -2380,        -8312 },
     {       2381,         8314 },
     {      -2381,        -8314 },
     {       2400,        -4101 },
     {      -2400,         4101 },
     {       2440,         4412 },
     {      -2440,        -4412 },
     {       2441,         4414 },
     {      -2441,        -4414 },
     {       2450,         5412 },
     {      -2450,        -5412 },
     {       2451,         5414 },
     {      -2451,        -5414 },
     {       2460,         6412 },
     {      -2460,        -6412 },
     {       2461,         6414 },
     {      -2461,        -6414 },
     {       2500,         5101 },
     {      -2500,        -5101 },
     {       2550,         5512 },
     {      -2550,        -5512 },
     {       2551,         5514 },
     {      -2551,        -5514 },
     {       2560,         6512 },
     {      -2560,        -6512 },
     {       2561,         6514 },
     {      -2561,        -6514 },
     {       2600,        -6101 },
     {      -2600,         6101 },
     {       2660,         6612 },
     {      -2660,        -6612 },
     {       2661,         6614 },
     {      -2661,        -6614 },
     {       3140,         4232 },
     {      -3140,        -4232 },
     {       3150,         5232 },
     {      -3150,        -5232 },
     {       3160,         6232 },
     {      -3160,        -6232 },
     {       3170,         7232 },
     {      -3170,        -7232 },
     {       3180,         8232 },
     {      -3180,        -8232 },
     {       3240,         4132 },
     {      -3240,        -4132 },
     {       3250,         5132 },
     {      -3250,        -5132 },
     {       3260,         6132 },
     {      -3260,        -6132 },
     {       3270,         7132 },
     {      -3270,        -7132 },
     {       3280,         8132 },
     {      -3280,        -8132 },
     {       3300,         3303 },
     {      -3300,        -3303 },
     {       3331,         3334 },
     {      -3331,        -3334 },
     {       3340,         4332 },
     {      -3340,        -4332 },
     {       3341,         4334 },
     {      -3341,        -4334 },
     {       3350,         5332 },
     {      -3350,        -5332 },
     {       3351,         5334 },
     {      -3351,        -5334 },
     {       3360,         6332 },
     {      -3360,        -6332 },
     {       3361,         6334 },
     {      -3361,        -6334 },
     {       3370,         7332 },
     {      -3370,        -7332 },
     {       3371,         7334 },
     {      -3371,        -7334 },
     {       3380,         8332 },
     {      -3380,        -8332 },
     {       3381,         8334 },
     {      -3381,        -8334 },
     {       3400,        -4301 },
     {      -3400,         4301 },
     {       3440,         4432 },
     {      -3440,        -4432 },
     {       3441,         4434 },
     {      -3441,        -4434 },
     {       3450,         5432 },
     {      -3450,        -5432 },
     {       3451,         5434 },
     {      -3451,        -5434 },
     {       3460,         6432 },
     {      -3460,        -6432 },
     {       3461,         6434 },
     {      -3461,        -6434 },
     {       3500,         5301 },
     {      -3500,        -5301 },
     {       3550,         5532 },
     {      -3550,        -5532 },
     {       3551,         5534 },
     {      -3551,        -5534 },
     {       3560,         6532 },
     {      -3560,        -6532 },
     {       3561,         6534 },
     {      -3561,        -6534 },
     {       3600,        -6301 },
     {      -3600,         6301 },
     {       3660,         6632 },
     {      -3660,        -6632 },
     {       3661,         6634 },
     {      -3661,        -6634 },
     {       4150,         5242 },
     {      -4150,        -5242 },
     {       4160,         6242 },
     {      -4160,        -6242 },
     {       4250,         5142 },
     {      -4250,        -5142 },
     {       4260,         6142 },
     {      -4260,        -6142 },
     {       4350,         5342 },
     {      -4350,        -5342 },
     {       4360,         6342 },
     {      -4360,        -6342 },
     {       4400,         4403 },
     {      -4400,        -4403 },
     {       4441,         4444 },
     {      -4441,        -4444 },
     {       4450,         5442 },
     {      -4450,        -5442 },
     {       4451,         5444 },
     {      -4451,        -5444 },
     {       4460,         6442 },
     {      -4460,        -6442 },
     {       4461,         6444 },
     {      -4461,        -6444 },
     {       4500,         5401 },
     {      -4500,        -5401 },
     {       4550,         5542 },
     {      -4550,        -5542 },
     {       4551,         5544 },
     {      -4551,        -5544 },
     {       4560,         6542 },
     {      -4560,        -6542 },
     {       4561,         6544 },
     {      -4561,        -6544 },
     {       4600,         6401 },
     {      -4600,        -6401 },
     {       4660,         6642 },
     {      -4660,        -6642 },
     {       4661,         6644 },
     {      -4661,        -6644 },
     {       5160,         6252 },
     {      -5160,        -6252 },
     {       5260,         6152 },
     {      -5260,        -6152 },
     {       5360,         6352 },
     {      -5360,        -6352 },
     {       5460,         6452 },
     {      -5460,        -6452 },
     {       5500,         5503 },
     {      -5500,        -5503 },
     {       5551,         5554 },
     {      -5551,        -5554 },
     {       5560,         6552 },
     {      -5560,        -6552 },
     {       5561,         6554 },
     {      -5561,        -6554 },
     {       5600,        -6501 },
     {      -5600,         6501 },
     {       5660,         6652 },
     {      -5660,        -6652 },
     {       5661,         6654 },
     {      -5661,        -6654 },
     {       6600,         6603 },
     {      -6600,        -6603 },
     {       6661,         6664 },
     {      -6661,        -6664 },
     {      10016,           93 },
     {     -10016,          -93 },
     {      20016,           94 },
     {     -20016,          -94 },
     {      10110,      9010221 },
     {      10111,        20113 },
     {      10121,        20213 },
     {     -10121,       -20213 },
     {      10131,        10323 },
     {     -10131,       -10323 },
     {      10231,        10313 },
     {     -10231,       -10313 },
     {      10441,       100443 },
     {      20440,        10441 },
     {      20441,        20443 },
     {      20442,          445 },
     {      30131,       100323 },
     {     -30131,      -100323 },
     {      30231,       100313 },
     {     -30231,      -100313 }
  };

  int listSize = sizeof(SList)/sizeof(SList[0]);
  for( int k=0; k!=listSize; ++k) {
      m.insert( std::make_pair( SList[k].id, SList[k].pid) );
  }
  return m;
}  // IsajetPDTMapInit()

PDTIsajetMap const & PDTIsajetMapInit()
{
    static PDTIsajetMap m;
    static IsajetPDTMap const & hmap = getIsajetPDTMap();
    
    for(IsajetPDTMap::const_iterator cit=hmap.begin(), mend=hmap.end(); cit!=mend; ++cit ) {
	m.insert( std::make_pair( cit->second, cit->first ));
    }
    return m;
}

// if a number isn't in the map, we try to convert it
int convIsajettoPDT( const int id )
{
    // we have no idea what to do, these numbers must be in the map
    if( abspid(id) <= 100 ) { return 0; }
    if( abspid(id) > 99999 ) { return 0; }

    // find constituents
    int istran;
    unsigned short js = digit(nj,id);
    unsigned short i1 = digit(nq3,id);
    unsigned short i2 = digit(nq2,id);
    unsigned short i3 = digit(nq1,id);
    unsigned short i4 = digit(nl,id);
    
    // mesons
    if(i1 != 0 && i2 != 0 && i3 == 0) {
          //   u and d have opposite definitions - sometimes
          if(i2 <= 2 && i1 <= 2){
              //     don't change
          } else {
	      if(i1 == 2) { 
	          i1 = 1; 
	      } else if(i1 == 1) { 
	          i1 = 2; 
	      }
	      if(i2 == 2) { 
	          i2 = 1; 
	      } else if(i2 == 1) { 
	          i2 = 2; 
	      }
          }
          istran=i1*100 + i2*10 + 2*js+1 + i4*10000;
          if( id < 0 ) { istran = -istran; }
          //  charmed and top mesons have wrong sign
          if(i1 == 4 && i2 != 4) { istran = -istran; }
          if(i1 == 6 && i2 != 6 && i2 != 4) { istran = -istran; }
          // ...check for illegal antiparticles
          if(i2 == i1 && id < 0) { istran=0; }
	  return istran;
    }
    // diquarks
    if(i2 != 0 && i3 != 0 && i1 == 0) {
          // ...         u and d have opposite definitions
	  if(i3 == 2) { 
	      i3 = 1; 
	  } else if(i3 == 1) { 
	      i3 = 2; 
	  }
	  if(i2 == 2) { 
	      i2 = 1; 
	  } else if(i2 == 1) { 
	      i2 = 2; 
	  }
	  istran = 0;
          if(i2 < i3){
            istran=i3*1000 + i2*100 + 1;
          } else if(i2 == i3){
            istran=i2*1000 + i3*100 + 3;
          } else {
            istran=i2*1000 + i3*100 + 1;
          }
          if( id < 0 ) { istran = -istran; }
          // ...         charmed and top quarks have wrong sign
          if(i2 == 4 && i3 != 4) { istran=-istran; }
          if(i2 == 6 && i3 != 6 && i3 != 4) { istran=-istran; }
	  return istran;
    }
    // baryons
    if( i1 != 0 && i3 != 0 && i2 != 0 ) {
          //   u and d have opposite definitions
	  if(i3 == 2) { 
	      i3 = 1; 
	  } else if(i3 == 1) { 
	      i3 = 2; 
	  }
	  if(i2 == 2) { 
	      i2 = 1; 
	  } else if(i2 == 1) { 
	      i2 = 2; 
	  }
	  if(i1 == 2) { 
	      i1 = 1; 
	  } else if(i1 == 1) { 
	      i1 = 2; 
	  }
	  istran = 0;
          if(i1 <= 2){
            istran=i3*1000 + i2*100 + i1*10 + 2*js+2;
          } else if(i3 <= 2 && i2 <= 2){
            istran=i1*1000 + i3*100 + i2*10 + 2*js+2;
          } else {
            istran=i1*1000 + i2*100 + i3*10 + 2*js+2;
          }
          if( id < 0 ) { istran = -istran; }
    }
    // unknown
    return 0;

}

// if a number isn't in the map, we try to convert it
int convPDTtoIsajet( const int id )
{
    // we have no idea what to do, these numbers must be in the map
    if( fundamentalID(id) != 0 ) { return 0; }
    if( abspid(id) > 99999 ) { return 0; }

    // find constituents
    int istran;
    unsigned short js = digit(nj,id);
    unsigned short i1 = digit(nq3,id);
    unsigned short i2 = digit(nq2,id);
    unsigned short i3 = digit(nq1,id);
    unsigned short i4 = digit(nl,id);

    // mesons
    if(i1 != 0 && i2 != 0 && i3 == 0) {
          //   u and d have opposite definitions - sometimes
          if(i2 <= 2 && i1 <= 2){
              //     don't change
          } else {
	      if(i1 == 2) { 
	          i1 = 1; 
	      } else if(i1 == 1) { 
	          i1 = 2; 
	      }
	      if(i2 == 2) { 
	          i2 = 1; 
	      } else if(i2 == 1) { 
	          i2 = 2; 
	      }
          }
          istran=i1*100 + i2*10 + (js-1)/2 + i4*10000;
          if( id < 0 ) { istran = -istran; }
          // ...         charmed and top mesons have wrong sign
          if(i2 == 4 && i1 != 4) { istran = -istran; }
          if(i2 == 6 && i1 != 6 && i1 != 4) { istran = -istran; }
          // ...check for illegal antiparticles
          if(i2 == i1 && id < 0) { istran=0; }
	  return istran;
    }
    // diquarks
    if(i1 == 0){
          // ...         u and d have opposite definitions
	  if(i3 == 2) { 
	      i3 = 1; 
	  } else if(i3 == 1) { 
	      i3 = 2; 
	  }
	  if(i2 == 2) { 
	      i2 = 1; 
	  } else if(i2 == 1) { 
	      i2 = 2; 
	  }
	  istran = 0;
          if(i3 < i2){
            istran=i3*1000 + i2*100 + (js-1)/2;
          } else {
            istran=i2*1000 + i3*100 + (js-1)/2;
          }
          if( id < 0 ) { istran = -istran; }
          // ...         charmed and top mesons have wrong sign
          if(i2 == 4 && i3 != 4) { istran=-istran; }
          if(i2 == 6 && i3 != 6 && i3 != 4) { istran=-istran; }
	  return istran;
    }
    // ...spin 1/2 or spin 3/2 baryons
    if( i1 != 0 && i3 != 0 && i2 != 0 && ( js == 2 || js == 4) && i4 == 0 ) {
          //   u and d have opposite definitions
	  if(i3 == 2) { 
	      i3 = 1; 
	  } else if(i3 == 1) { 
	      i3 = 2; 
	  }
	  if(i2 == 2) { 
	      i2 = 1; 
	  } else if(i2 == 1) { 
	      i2 = 2; 
	  }
	  if(i1 == 2) { 
	      i1 = 1; 
	  } else if(i1 == 1) { 
	      i1 = 2; 
	  }
	  istran = 0;
          if(i3 <= 2){
            istran=i3*1000 + i2*100 + i1*10 + (js-2)/2;
          } else if(i1 <= 2 && i2 <= 2){
            istran=i2*1000 + i1*100 + i3*10 + (js-2)/2;
          } else {
            istran=i1*1000 + i2*100 + i3*10 + (js-2)/2;
          }
          if( id < 0 ) { istran = -istran; }
	  return istran;
    }
    // unknown
    return 0;
}
  
//
// getIsajetPDTMap is the ONLY function allowed to call IsajetPDTMapInit
//
IsajetPDTMap const & getIsajetPDTMap()
{
  static IsajetPDTMap const & hmap = IsajetPDTMapInit();
  return hmap;
}  // getIsajetPDTMap()

//
// getPDTIsajetMap is the ONLY function allowed to call IsajetPDTMapInit
//
PDTIsajetMap const & getPDTIsajetMap()
{
  static PDTIsajetMap const & hmap = PDTIsajetMapInit();
  return hmap;
}  // getPDTIsajetMap()
 
} // unnamed namespace

int translateIsajettoPDT( const int id )
{
    static IsajetPDTMap const & hmap = getIsajetPDTMap();

    IsajetPDTMap::const_iterator const cit = hmap.find( id );
    // found it in the map
    if ( cit != hmap.end() ) { return cit->second; }
    // try converting anyway
    return convIsajettoPDT(id);
}

int translatePDTtoIsajet( const int id )
{
    static PDTIsajetMap const & pmap = getPDTIsajetMap();

    PDTIsajetMap::const_iterator const cit = pmap.find( id );
    // found it in the map
    if ( cit != pmap.end() ) { return cit->second; }
    // try converting anyway
    return convPDTtoIsajet(id);
}

void writeIsajetTranslationLine( int i, std::ostream & os  )
{
    // only write map entries
    static IsajetPDTMap const & hmap = getIsajetPDTMap();

    IsajetPDTMap::const_iterator const cit = hmap.find( i );
    // found it in the map
    if ( cit != hmap.end() ) { 
        int id = cit->second;
	os << " Isajet number: " ;
	os.width(10);
	os << i << "  HepPID number: " ;
	os.width(10);
	os << id << "  " << particleName(id) << std::endl;
	// check reverse translation
	int iback =  translatePDTtoIsajet(id);
	if(iback != i) {
	    os << " WARNING: " << id << " translates back to " 
	       << iback << " not to " << i << std::endl;
	}
    }
    return;
}  // writeIsajetTranslationLine()

void  writeIsajetTranslation( std::ostream & os )
{
    writeVersion( os );
    os << "     HepPID Particle List" << std::endl;
    os << std::endl;

    int id, j, q1, q2, q3, m;
    // special cases
    for( id=1; id<101; ++id) {
        writeIsajetTranslationLine(  id, os );
        writeIsajetTranslationLine( -id, os );
    }
    // diquark
    for( q2=1; q2<7; ++q2) {
	for( q1=1; q1<7; ++q1) {
	    for( j=0; j<2; ++j) {
        	id = 1000*q2+100*q1+j;
        	writeIsajetTranslationLine(  id, os );
        	writeIsajetTranslationLine( -id, os );
	    }
	}
    }
    // mesons
    for( q2=1; q2<9; ++q2) {
	for( q1=1; q1<9; ++q1) {
	    for( j=0; j<3; ++j) {
		for( m=0; m<4; ++m) {
		    id = 10000*m+100*q2+10*q1+j;
        	    writeIsajetTranslationLine(  id, os );
        	    writeIsajetTranslationLine( -id, os );
		}
	    }
	}
    }
    // baryons
    for( q3=1; q3<7; ++q3) {
	for( q2=1; q2<7; ++q2) {
	    for( q1=1; q1<7; ++q1) {
		for( j=1; j<2; ++j) {
		    id = 1000*q3+100*q2+10*q1+j;
        	    writeIsajetTranslationLine(  id, os );
        	    writeIsajetTranslationLine( -id, os );
		}
	    }
	}
    }
    return;
}  // writeIsajetTranslation()

}	// HepPID
