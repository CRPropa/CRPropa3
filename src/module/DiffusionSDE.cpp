#include "crpropa/module/DiffusionSDE.h"

#include <crpropa/Random.h>

#include <iostream>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdexcept>

using namespace crpropa;

// Defining Cash-Karp coefficients
const double a[] = { 0., 0., 0., 0., 0., 0., 1. / 5., 0., 0., 0., 0.,
		0., 3. / 40., 9. / 40., 0., 0., 0., 0., 3. / 10., -9. / 10., 6. / 5.,
		0., 0., 0., -11. / 54., 5. / 2., -70. / 27., 35. / 27., 0., 0., 1631.
				/ 55296., 175. / 512., 575. / 13824., 44275. / 110592., 253.
				/ 4096., 0. };

const double b[] = { 37. / 378., 0, 250. / 621., 125. / 594., 0., 512.
		/ 1771. };

const double bs[] = { 2825. / 27648., 0., 18575. / 48384., 13525.
		/ 55296., 277. / 14336., 1. / 4. };



DiffusionSDE::DiffusionSDE(ref_ptr<MagneticField> field, double tolerance, 
				 double minStep, double maxStep, double epsilon) :
  minStep(0)
{
  setField(field);
  setMaximumStep(maxStep);
  setMinimumStep(minStep);
  setTolerance(tolerance);
  setEpsilon(epsilon);
  setScale(1.);
  setAlpha(1./3.);
  }


void DiffusionSDE::process(Candidate *candidate) const {
	// save the new previous particle state
	ParticleState &current = candidate->current;
	candidate->previous = current;
	
	double step = clip(candidate->getNextStep(), minStep, maxStep);

	// rectilinear propagation for neutral particles
	if (current.getCharge() == 0) {
		Vector3d dir = current.getDirection();
		current.setPosition(current.getPosition() + dir * step);
		candidate->setCurrentStep(step);
		candidate->setNextStep(maxStep);
		return;
	}

	Vector3d PosIn = current.getPosition();
	Vector3d DirIn = current.getDirection();
	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();

    // Calculate the Diffusion tensor
	double BTensor[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	calculateBTensor(rig, BTensor, PosIn, DirIn, z);
    
    // Generate random numers
	double eta[] = {0., 0., 0.};
	for(size_t i=0; i < 3; i++) {
	  eta[i] =  Random::instance().randNorm();
	}

	double TStep = BTensor[0] * eta[0];
	double NStep = BTensor[4] * eta[1];
	double BStep = BTensor[8] * eta[2];

	double h = step / c_light;
	double hTry, r;

	Vector3d TVec(0.);
	Vector3d NVec(0.);
	Vector3d BVec(0.);
	Vector3d PosOut = Vector3d(0.);
	Vector3d DirOut = Vector3d(0.);
	Vector3d PosErr = Vector3d(0.);
	Vector3d PosTest = Vector3d(0.);



	do {
	  hTry = h;
	  double propStep =  TStep * pow(hTry, 0.5) / c_light;
	  tryStep(PosIn, PosOut, PosErr,PosTest, TVec, NVec, BVec, z, propStep);
	  // calculate the relative position error r and the next time step h
	  r = PosErr.getR() / tolerance;
	  h *= 0.95 * pow(r, -0.2);
	  // prevent h from too strong variations
	  h = clip(h, 0.1 * hTry, 5 * hTry);

	} while (r > 1 && hTry >= minStep / c_light && TVec.getR()==TVec.getR());
	
	// Exception: Rectilinear propagation in case of vanishing magnetic field.
	if (TVec.getR() != TVec.getR()) {
	  Vector3d dir = current.getDirection();
      current.setPosition(current.getPosition() + dir * step);
	  candidate->setCurrentStep(step);
	  candidate->setNextStep(step);
	  return;
	}
	// Integration of the SDE with a Mayorama-Euler-method
	Vector3d PO = PosOut + (NVec * NStep + BVec * BStep) * pow(hTry, 0.5) ;
    
    // Throw error message if somethin went wrong with propagation.
    // Deactivate candidate.
	bool NaN = std::isnan(PO.getR());
	if (NaN == true){
	  std::cout << "\nCandidate with 'nan'-position occured: \n";
	  std::cout << "position = " << PO << "\n";
	  std::cout << "PosIn = " << PosIn << "\n";
	  std::cout << "TVec = " << TVec << "\n";
	  std::cout << "TStep = " << std::abs(TStep) << "\n";
	  std::cout << "NVec = " << NVec << "\n";
	  std::cout << "NStep = " << NStep << "\n";
	  std::cout << "BVec = " << BVec << "\n";
	  std::cout << "BStep = " << BStep << "\n";
	  candidate->setActive(false);
	  std::cout << "Candidate is deactivated!\n";
	  std::cout << "-------------------------\n";
	  return;
	}
	
	DirOut = (PO -PosIn).getUnitVector();
	current.setPosition(PO);
	current.setDirection(DirOut);
	candidate->setCurrentStep(hTry * c_light);
	candidate->setNextStep(h * c_light);
    // Debugging and Testing
    // Delete comments if additional information should be stored in candidate
    /*
	std::stringstream s;
	const std::string AL = "arcLength";
	if (candidate->hasProperty(AL) == false){
	  s << (TStep + NStep + BStep) * pow(hTry, 0.5);
	  const std::string value = s.str();
	  candidate->setProperty(AL, value);
	  return;
	}
	else {
	  std::string arcLenString;
	  candidate->getProperty(AL, arcLenString);
	  double arcLen = ::atof(arcLenString.c_str());
	  arcLen += (TStep + NStep + BStep) * pow(hTry, 0.5);
	  s << arcLen;
	  const std::string value = s.str();
	  candidate->setProperty(AL, value);
	}
*/
}


void DiffusionSDE::tryStep(const Vector3d &PosIn, Vector3d &POut, Vector3d &PosErr,Vector3d &PosTest, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec,double z, double propStep) const {

	Vector3d k[] = {Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.)};
	POut = PosIn;
	PosTest = PosIn;
	//calculate the sume k_i * b_i
	for (size_t i = 0; i < 6; i++) {

		Vector3d y_n = PosIn;
		for (size_t j = 0; j < i; j++)
		  y_n += k[j] * a[i * 6 + j] * propStep;
		
		// update k_i = direction of the regular magnetic mean field
		Vector3d BField(0.);
		try {
		  BField = field->getField(y_n, z);
		} 
		catch (std::exception &e) {
		  std::cerr << "DiffusionSDE: Exception in getField." << std::endl;
		  std::cerr << e.what() << std::endl;
		}
		
		k[i] = BField.getUnitVector() * c_light;

		POut += k[i] * b[i] * propStep;
		PosTest += k[i] * bs[i] * propStep;

		PosErr +=  (k[i] * (b[i] - bs[i])) / c_light;
	}
	TVec = (POut-PosIn).getUnitVector();
	// Choose a random perpendicular vector as the Normal-vector.
	// Prevent 'nan's in the NVec-vector in the case of <TVec, NVec> = 0.
	while (NVec.getR()==0.){
	  Vector3d RandomVector = Random::instance().randVector();
	  NVec = TVec.cross( RandomVector );
	}
	NVec = NVec.getUnitVector();

	// Calculate the Binormal-vector
	BVec = (TVec.cross(NVec)).getUnitVector();
}

void DiffusionSDE::calculateBTensor(double r, double BTen[], Vector3d pos, Vector3d dir, double z) const {
 
    double DifCoeff = scale * 6.1e24 * pow((std::abs(r) / 4.0e9), alpha);
    BTen[0] = pow( 2  * DifCoeff, 0.5);
    BTen[4] = pow(2 * epsilon * DifCoeff, 0.5);
    BTen[8] = pow(2 * epsilon * DifCoeff, 0.5);
    return;

}


void DiffusionSDE::setMinimumStep(double min) {
	if (min < 0)
		throw std::runtime_error("DiffusionSDE: minStep < 0 ");
	if (min > maxStep)
		throw std::runtime_error("DiffusionSDE: minStep > maxStep");
	minStep = min;
}

void DiffusionSDE::setMaximumStep(double max) {
	if (max < minStep)
		throw std::runtime_error("DiffusionSDE: maxStep < minStep");
	maxStep = max;
}


void DiffusionSDE::setTolerance(double tol) {
	if ((tol > 1) or (tol < 0))
		throw std::runtime_error(
				"DiffusionSDE: tolerance error not in range 0-1");
	tolerance = tol;
}

void DiffusionSDE::setEpsilon(double e) {
	if ((e > 1) or (e < 0))
		throw std::runtime_error(
				"DiffusionSDE: epsilon not in range 0-1");
	epsilon = e;
}


void DiffusionSDE::setAlpha(double a) {
	if ((a > 2.) or (a < 0))
		throw std::runtime_error(
				"DiffusionSDE: alpha not in range 0-2");
	alpha = a;
}

void DiffusionSDE::setScale(double s) {
	if (s < 0)
		throw std::runtime_error(
				"DiffusionSDE: Scale error: Scale < 0");
	scale = s;
}

void DiffusionSDE::setField(ref_ptr<MagneticField> f) {
	field = f;
}

double DiffusionSDE::getMinimumStep() const {
	return minStep;
}

double DiffusionSDE::getMaximumStep() const {
	return maxStep;
}

double DiffusionSDE::getTolerance() const {
	return tolerance;
}

double DiffusionSDE::getEpsilon() const {
	return epsilon;
}

double DiffusionSDE::getAlpha() const {
	return alpha;
}

double DiffusionSDE::getScale() const {
	return scale;
}



std::string DiffusionSDE::getDescription() const {
	std::stringstream s;
	s << "minStep: " << minStep / kpc  << " kpc, ";
	s << "maxStep: " << maxStep / kpc  << " kpc, ";
	s << "tolerance: " << tolerance << "\n";
	
	if (epsilon != 0.1) {
	  s << "epsilon: " << epsilon << ", ";
	  }
	
	if (alpha != 1./3.) {
	  s << "alpha: " << alpha << "\n";
	  }

	if (scale != 1.) {
	  s << "scale: " << scale << "\n";
	  }

	return s.str();
}
