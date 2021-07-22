#include "crpropa/DiffusionTensor.h"
#include "crpropa/Units.h"
#include "kiss/logger.h"


using namespace crpropa;

/*
void DiffusionTensor::setDescription(const std::string &d) {
    description = d;
}

*/


// QLTDiffusion ------------------------------------------------------------------
QLTDiffusion::QLTDiffusion(double eps , double kap, double alp ){
    epsilon = eps;
    kappa0 = kap;
    alpha = alp;
}

void QLTDiffusion::setEpsilon(double eps){
    epsilon = eps;
}
void QLTDiffusion::setKappa0(double kap){
    kappa0 = kap;
}
void QLTDiffusion::setAlpha(double alph){
    alpha = alph;
}

double QLTDiffusion::getAlpha() const{
    return alpha;
}
double QLTDiffusion::getKappa0() const{
    return kappa0;
}
double QLTDiffusion::getEpsilon() const{
    return epsilon;
}

double QLTDiffusion::getKappaParallel(Candidate *cand){
    double rig = cand-> current.getRigidity();
    return kappa0*pow(std::abs(rig)/(4.0e9*volt),alpha);
}

double QLTDiffusion::getKappaPerpendicular(Candidate *cand){
    double rig = cand-> current.getRigidity();
    return epsilon*kappa0*pow(std::abs(rig)/(4.0e9*volt),alpha);
}

double QLTDiffusion::getKappaPerpendicular2(Candidate *cand){
    return getKappaPerpendicular(cand); 
}

std::string QLTDiffusion::getDescription() const{
    std::stringstream s;
    s << "Diffusion tensor for quasi-linear-theory (QLT)  with: \n"
    << "epsilon: \t" << epsilon <<"\n"
    << "kappa0: \t"<< kappa0 << "m2/s \n"
    << "alpha: \t" << alpha << "\n";
    return s.str();
}

// QLTTurbulent ------------------------------------------------------------------

QLTTurbulent::QLTTurbulent(ref_ptr<MagneticField> background, ref_ptr<MagneticField> turbulent, double kap, double aPara, double aPerp){
    backgroundField = background;
    turbulentField = turbulent;
    kappa0 = kap;
    alphaPara = aPara;
    alphaPerp = aPerp;
    normToEarthPosition(); 
}

double QLTTurbulent::getKappaParallel(Candidate *cand) const {
	ParticleState &current = cand->current;
    double rig = current.getRigidity();
    Vector3d pos = current.getPosition();
    double eta = turbulentField->getField(pos).getR() / backgroundField -> getField(pos).getR();
    double k = kappa0 * eta*eta/normTurbulence/normTurbulence * pow(rig/(4e9*volt), alphaPara);
    return k;
}

double QLTTurbulent::getKappaPerpendicular(Candidate *cand) const{
	ParticleState &current = cand->current;
    double rig = current.getRigidity();
    Vector3d pos = current.getPosition();
    double eta = turbulentField->getField(pos).getR() / backgroundField -> getField(pos).getR();
    double k = kappa0 * eta*eta*normTurbulence*normTurbulence* pow(rig/(4e9*volt), alphaPerp);
    return k;
}

double QLTTurbulent::getKappaPerpendicular2(Candidate *cand) const{
    return getKappaPerpendicular(cand);
}

double QLTTurbulent::getKappa0() const{
    return kappa0;
}
double QLTTurbulent::getAlphaPara() const{
    return alphaPara;
}
double QLTTurbulent::getAlphaPerp() const{
    return alphaPerp;
}
double QLTTurbulent::getNormTurbulence() const{
    return normTurbulence;
}

void QLTTurbulent::setKappa0(double kap){
    kappa0 = kap;
}
void QLTTurbulent::setAlphaPara(double a){
    alphaPara = a;
}
void QLTTurbulent::setAlphaPerp(double a){
    alphaPerp = a;
}
void QLTTurbulent::setAlpha(double a){
    alphaPara = a;
    alphaPerp = a;
}
void QLTTurbulent::setNormTurbulence(double eta){
    normTurbulence = eta;
}

void QLTTurbulent::normToEarthPosition(Vector3d posEarth){
    double b = turbulentField -> getField(posEarth).getR();
    double B = backgroundField ->getField(posEarth).getR();
    setNormTurbulence(b/B);
}

std::string QLTTurbulent::getDescription() const{
    std::stringstream ss;
    ss << "Diffusion Tensor for the quasi-linear-theory (QLT) with turbulence dependence \n \n"
    << "Energyscaling with spectral index:\n"
    << "alpha parallel \t" << alphaPara << "\n"
    << "alpha perpendicular \t" << alphaPerp << "\n \n"
    << "norming the turbulence at value: \t" << normTurbulence << "\n"
    << "norming the diffusion coefficent to: \t" << kappa0 << "m2/s \n"; 
    return ss.str();
}

// QLTRigidity ------------------------------------------------------------------

double QLTRigidity::calculateLamorRadius(ParticleState &state) const {
    double fieldStrength;
    Vector3d pos = state.getPosition();

    if(hasTurbulentField){
        fieldStrength = (backgroundField->getField(pos) + turbulentField ->getField(pos)).getR();
    }
    else{ 
        // no turbulent field: total field from background
        fieldStrength = backgroundField -> getField(pos).getR(); 
    }

    return state.getMomentum().getR()/std::abs(state.getCharge())/fieldStrength;
}

QLTRigidity::QLTRigidity(ref_ptr<MagneticField> magField, ref_ptr<TurbulentField> turbField, double kap, double aPara, double aPerp){
    hasTurbulentField = true;
    setMagneticField(magField);
    setTurbulentField(turbField);
    kappa0 = kap;
    alphaPara = aPara;
    alphaPerp = aPerp;
    normToPosition(); // default norming to earth, change later if neccesarry
}

QLTRigidity::QLTRigidity(ref_ptr<MagneticField> magField, double kap, double aPara, double aPerp){
    hasTurbulentField = false;
    setMagneticField(magField);
    kappa0 = kap;
    alphaPara = aPara;
    alphaPerp = aPerp;
    normToPosition(); // default norming to earth, change later if neccesarry
}

void QLTRigidity::setMagneticField(ref_ptr<MagneticField> field){
    backgroundField = field;
    normToPosition(normPos);
}

void QLTRigidity::setTurbulentField(ref_ptr<TurbulentField> field){
    turbulentField = field;
    hasTurbulentField = true;
    correlationLength = field -> getCorrelationLength();
    normToPosition(normPos);
}

void QLTRigidity::setKappa0(double kap){
    kappa0 = kap;
}
void QLTRigidity::setAlphaPara(double a){
    alphaPara = a;
}
void QLTRigidity::setAlphaPerp(double a){
    alphaPerp = a;
}
void QLTRigidity::setAlpha(double a){
    alphaPara = a;
    alphaPerp = a;
}

void QLTRigidity::normToPosition(const Vector3d &pos){
    normPos = pos;
    if(hasTurbulentField) {
        Vector3d b = turbulentField -> getField(pos);
        Vector3d B = backgroundField -> getField(pos);
        normEta = b.getR()/B.getR();
        normB = (b+B).getR();
    }
    else{
        try
        {
            double b = backgroundField -> getTurbulentField(pos).getR();
            double B = backgroundField -> getRegularField(pos).getR();
            normEta = b/B;
            normB = backgroundField -> getField(pos).getR();
        }
        catch(const std::exception& e)
        {
            KISS_LOG_ERROR << "Diffusiontensor QLTRigidty: Exception in calculating turbulence level \n"
            << "pleas use a magnatic field with function getTurbulenceStrength, getRegularField and getField or provide a seperate turbulenceField\n"
            << "Following Eception occures:\n" << e.what();
        }
        
    }
}

void QLTRigidity::setNormEta(double eta){
    normEta = eta;
}

void QLTRigidity::setNormB(double B){
    normB = B;
}

void QLTRigidity::setHasTurbulentField(bool use){
    hasTurbulentField = use;
}

