#include "crpropa/DiffusionTensor.h"

using namespace crpropa;


void DiffusionTensor::setDescription(const std::string &d) {
    description = d;
}


/*QLTDiffusion::QLTDiffusion(){
    setEpsilon(0.1);
    setKappa0(6.1e24);
    setAlpha(1./3.);
}*/

QLTDiffusion::QLTDiffusion(double epsilon , double kappa0, double alpha ){
    setEpsilon(epsilon);
    setKappa0(kappa0);
    setAlpha(alpha);
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

void QLTDiffusion::calculateBTensor(double BTen[], Candidate *cand){
    double rig = cand-> current.getEnergy() / cand->current.getCharge();
    double DifCoeff = kappa0*pow(std::abs(rig)/ 4.0e9, alpha);
    BTen[0] = pow( 2  * DifCoeff, 0.5);
    BTen[4] = pow(2 * epsilon * DifCoeff, 0.5);
    BTen[8] = pow(2 * epsilon * DifCoeff, 0.5);
    return;
}

std::string QLTDiffusion::getDescription() const{
    std::stringstream s;
    s << "Diffusion tensor for quasi-linear-theory (QLT)  with: \n"
    << "epsilon: \t" << epsilon <<"\n"
    << "kappa0: \t"<< kappa0 << "m2/s \n"
    << "alpha: \t" << alpha << "\n";
    return s.str();
}