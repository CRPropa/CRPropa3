#include "crpropa/DiffusionTensor.h"
#include "crpropa/Units.h"

using namespace crpropa;

/*
void DiffusionTensor::setDescription(const std::string &d) {
    description = d;
}

*/
QLTDiffusion::QLTDiffusion(double eps , double kap, double alp ){
    epsilon = eps;
    kappa0 = kap;
    alpha = alp;
}

void QLTDiffusion::setEpsilon(double epsilon){
    epsilon = epsilon;
}
void QLTDiffusion::setKappa0(double kappa0){
    kappa0 = kappa0;
}
void QLTDiffusion::setAlpha(double alpha){
    alpha = alpha;
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
    double rig = cand-> current.getEnergy()/cand->current.getCharge();
    return kappa0*pow(std::abs(rig)/(4.0e9*volt),alpha);
}

double QLTDiffusion::getKappaPerpendicular(Candidate *cand){
    double rig = cand-> current.getEnergy()/cand->current.getCharge();
    return epsilon*kappa0*pow(std::abs(rig)/(4.0e9*volt),alpha);
}

double QLTDiffusion::getKappaPerpendicular2(Candidate *cand){
    double rig = cand-> current.getEnergy()/cand->current.getCharge();
    return epsilon*kappa0*pow(std::abs(rig)/(4.0e9*volt),alpha);    
}

std::string QLTDiffusion::getDescription() const{
    std::stringstream s;
    s << "Diffusion tensor for quasi-linear-theory (QLT)  with: \n"
    << "epsilon: \t" << epsilon <<"\n"
    << "kappa0: \t"<< kappa0 << "m2/s \n"
    << "alpha: \t" << alpha << "\n";
    return s.str();
}