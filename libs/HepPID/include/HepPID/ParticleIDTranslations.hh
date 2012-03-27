#ifndef PARTICLE_ID_TRANSLATIONS_HH
#define PARTICLE_ID_TRANSLATIONS_HH
// ----------------------------------------------------------------------
//
// ParticleIDTranslations.hh
// Author: Lynn Garren
//
// ..convert between various numbering implementations
//
// ----------------------------------------------------------------------

#include <iostream>

//! The HepPID namespace has independent particle ID translation methods

///
/// \namespace HepPID
/// The HepPID namespace contains a set of independent 
/// particle ID translation methods
///
namespace HepPID {

// translate between generator ID's and standard numbering scheme

// Herwig translations
/// translate Herwig to PDG standard
int   translateHerwigtoPDT( const int herwigID);
/// translate PDG standard to Herwig
int   translatePDTtoHerwig( const int pid );
/// output the translation list
void  writeHerwigTranslation( std::ostream & os );

// Isajet translations
/// translate Isajet to PDG standard
int   translateIsajettoPDT( const int isajetID );
/// translate PDG standard to Isajet
int   translatePDTtoIsajet( const int pid );
/// output the translation list
void  writeIsajetTranslation( std::ostream & os );

// Pythia translations
/// translate Pythia to PDG standard
int   translatePythiatoPDT( const int pythiaID );
/// translate PDG standard to Pythia
int   translatePDTtoPythia( const int pid );
/// output the translation list
void  writePythiaTranslation( std::ostream & os );

// EvtGen translations
/// translate EvtGen to PDG standard
int   translateEvtGentoPDT( const int evtGenID );
/// translate PDG standard to EvtGen
int   translatePDTtoEvtGen( const int pid );
/// output the translation list
void  writeEvtGenTranslation( std::ostream & os );

// PDG table translations (yes,there can be differences)
/// translate PDG table to PDG standard
int   translatePDGtabletoPDT( const int pdgID);
/// translate PDG standard to PDG table
int   translatePDTtoPDGtable( const int pid );
/// output the translation list
void  writePDGTranslation( std::ostream & os );

// QQ translations
/// translate QQ to PDG standard
int   translateQQtoPDT( const int qqID);
/// translate PDG standard to QQ
int   translatePDTtoQQ( const int pid );
/// QQ helper function
int   translateQQbar( const int id );
/// QQ helper function
int   translateInverseQQbar( const int id );
/// output the translation list
void  writeQQTranslation( std::ostream & os );

// Geant3 translations
/// translate Geant3 to PDG standard
int translateGeanttoPDT( const int geantID);
/// translate PDG standard to Geant3
int translatePDTtoGeant( const int pid );

}  // namespace HepPID

#endif // PARTICLE_ID_TRANSLATIONS_HH
