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
	
	double h = clip(candidate->getNextStep(), minStep, maxStep) / c_light;
	Vector3d PosIn = current.getPosition();
	Vector3d DirIn = current.getDirection();

    // rectilinear propagation for neutral particles
	if (current.getCharge() == 0) {
		Vector3d dir = current.getDirection();
		current.setPosition(current.getPosition() + dir * h * c_light);
		candidate->setCurrentStep(h * c_light);
		candidate->setNextStep(maxStep);
		return;
	}
	
	double z = candidate->getRedshift();
	double rig = current.getEnergy() / current.getCharge();

    // Calculate the Diffusion tensor
	double BTensor[] = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
	calculateBTensor(rig, BTensor, PosIn, DirIn, z);
    
    // Generate random numbers
	double eta[] = {0., 0., 0.};
	for(size_t i=0; i < 3; i++) {
	  	eta[i] =  Random::instance().randNorm();
	}
	
	double TStep = BTensor[0] * eta[0];
	double NStep = BTensor[4] * eta[1];
	double BStep = BTensor[8] * eta[2];

	Vector3d TVec(0.);
	Vector3d NVec(0.);
	Vector3d BVec(0.);
	
	Vector3d DirOut = Vector3d(0.);


	double propTime = TStep * pow(h, 0.5) / c_light;
	size_t counter = 0;
	double r=42.; //arbitrary number larger than one


	do {
		Vector3d PosOut = Vector3d(0.);
		Vector3d PosErr = Vector3d(0.);
	  	tryStep(PosIn, PosOut, PosErr, z, propTime);
	    // calculate the relative position error r and the next time step h
	  	r = PosErr.getR() / tolerance;
	  	propTime *= 0.5;
		counter += 1;

    //Check for better break condition
    // Check zero magnetic field
	} while (r > 1 && fabs(propTime) >= minStep/c_light); 
	//std::cout << "r = " << r << "\n";
	//std::cout << "propTime = " << propTime*c_light/kpc << "\n";
	//std::cout << "counter = " << counter << "\n";
	
	size_t stepNumber = pow(2, counter-1);
	double allowedTime = TStep * pow(h, 0.5) / c_light / stepNumber ;
	//std::cout << "allowed Time = " << allowedTime * c_light / kpc << "\n";
	Vector3d Start = PosIn;
	Vector3d PosOut = Vector3d(0.);
	Vector3d PosErr = Vector3d(0.);
	for (size_t j=0; j<stepNumber; j++) {
		tryStep(Start, PosOut, PosErr, z, allowedTime);
		Start = PosOut;
	}		
	
    // Normalize the tangent vector
	TVec = (PosOut-PosIn).getUnitVector();

    // Exception: Rectilinear propagation in case of vanishing magnetic field.
	if (TVec.getR() != TVec.getR()) {
	  	Vector3d dir = current.getDirection();
      		current.setPosition(current.getPosition() + dir * h *c_light);
	 	candidate->setCurrentStep(h *c_light);
	  	candidate->setNextStep(h *c_light);
	  	return;
	}
	
    // Choose a random perpendicular vector as the Normal-vector.
    // Prevent 'nan's in the NVec-vector in the case of <TVec, NVec> = 0.
	while (NVec.getR()==0.){
	  	Vector3d RandomVector = Random::instance().randVector();
	  	NVec = TVec.cross( RandomVector );
	}
	NVec = NVec.getUnitVector();

   // Calculate the Binormal-vector
	BVec = (TVec.cross(NVec)).getUnitVector();
   

    // Integration of the SDE with a Mayorama-Euler-method
	Vector3d PO = PosOut + (NVec * NStep + BVec * BStep) * pow(h, 0.5) ;
    
    // Throw error message if something went wrong with propagation.
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
	
	DirOut = (PO - PosIn).getUnitVector();
	current.setPosition(PO);
	current.setDirection(DirOut);
	candidate->setCurrentStep(h * c_light);
/*
	double nextStep = allowedTime * 0.95 * pow(r, -0.2);
	nextStep = clip(nextStep, 0.1*allowedTime, 5.*allowedTime)*c_light;
	nextStep = pow(nextStep/(BTensor[0]*0.8), 2.);
*/
	//double nextTime = pow((allowedTime*c_light)/(BTensor[0]*0.8), 2.);
	//std::cout << "nextTime = " << nextTime*c_light/kpc << "\n";
	//nextTime = nextTime * 0.95 *pow(r, -0.2);
	//std::cout << "nextTime = " << nextTime*c_light/kpc << "\n";
	//double nextStep = clip(nextTime, 0.1*h, 5.*h)*c_light;
	double nextStep;
	if (stepNumber>1){
		nextStep = h*pow(stepNumber, -2.)*c_light;
	}
	else {
		nextStep = 4 * h*c_light;
	}

	candidate->setNextStep(nextStep);
	//candidate->setNextStep(maxStep);
    // Debugging and Testing
    // Delete comments if additional information should be stored in candidate
    
	std::stringstream s;
	const std::string AL = "arcLength";
	if (candidate->hasProperty(AL) == false){
	  s << (TStep + NStep + BStep) * pow(h, 0.5);
	  const std::string value = s.str();
	  candidate->setProperty(AL, value);
	  return;
	}
	else {
	  std::string arcLenString;
	  candidate->getProperty(AL, arcLenString);
	  double arcLen = ::atof(arcLenString.c_str());
	  arcLen += (TStep + NStep + BStep) * pow(h, 0.5);
	  s << arcLen;
	  const std::string value = s.str();
	  candidate->setProperty(AL, value);
	}

}


void DiffusionSDE::tryStep(const Vector3d &PosIn, Vector3d &POut, Vector3d &PosErr,double z, double propStep) const {

	Vector3d k[] = {Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.),Vector3d(0.)};
	POut = PosIn;
	//calculate the sum k_i * b_i
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
		//PosErr +=  (k[i] * (b[i] - bs[i])) / c_light;
		PosErr +=  (k[i] * (b[i] - bs[i])) * propStep / kpc;
		
	}
	//std::cout << "PosErr = " << PosErr <<"\n";
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
