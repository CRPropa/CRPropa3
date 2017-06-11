#include "crpropa/advectionField/AdvectionField.h"


namespace crpropa {



void AdvectionFieldList::addField(ref_ptr<AdvectionField> field) {
	fields.push_back(field);
}

Vector3d AdvectionFieldList::getField(const Vector3d &position) const {
	Vector3d b(0.);
	for (int i = 0; i < fields.size(); i++)
		b += fields[i]->getField(position);
	return b;
}

double AdvectionFieldList::getDivergence(const Vector3d &position) const {
	double D=0.;
	// Work on default values for divergence or an error handling
	for (int i = 0; i < fields.size(); i++)
		D += fields[i]->getDivergence(position);
	return D;
}


//----------------------------------------------------------------
UniformAdvectionField::UniformAdvectionField(const Vector3d &value) :
			value(value) {
	}

Vector3d UniformAdvectionField::getField(const Vector3d &position) const {
	return value;
	}

double UniformAdvectionField::getDivergence(const Vector3d &position) const {
	return 0.;
	}

//----------------------------------------------------------------

ConstantSphericalAdvectionField::ConstantSphericalAdvectionField(Vector3d origin, double vWind) {
	setOrigin(origin);
	setVWind(vWind);
}

Vector3d ConstantSphericalAdvectionField::getField(const Vector3d &position) const {
	Vector3d Pos = position-origin;
	return vWind * Pos.getUnitVector();
}

double ConstantSphericalAdvectionField::getDivergence(const Vector3d &position) const {
	double R = (position-origin).getR();	
	return 2*vWind/R;
}

void ConstantSphericalAdvectionField::setOrigin(Vector3d o) {
	origin=o;
	return;
}

void ConstantSphericalAdvectionField::setVWind(double v) {
	vWind = v;
	return;
}

Vector3d ConstantSphericalAdvectionField::getOrigin() const {
	return origin;
}

double ConstantSphericalAdvectionField::getVWind() const {
	return vWind;
}



//----------------------------------------------------------------

SphericalAdvectionField::SphericalAdvectionField(Vector3d origin, double radius, double vMax, double tau, double alpha) {
	setOrigin(origin);
	setRadius(radius);
	setVMax(vMax);
	setTau(tau);
	setAlpha(alpha);
}

Vector3d SphericalAdvectionField::getField(const Vector3d &position) const {
	Vector3d Pos = position-origin;
	double R = Pos.getR();
	if (R>radius) {
		return Vector3d(0.);
	}
	double v_R = getV(R);
	return v_R * Pos.getUnitVector();
}

double SphericalAdvectionField::getDivergence(const Vector3d &position) const {
	double R = (position-origin).getR();
	if (R>radius) {
		return 0.;
	}
	double D = 2*vMax/R * ( 1-( 1-alpha*(pow(R, alpha)/(2*tau)) )*exp(-( pow(R, alpha)/tau )) );
	return D;
}

double SphericalAdvectionField::getV(const double &r) const {
	double f = vMax * (1-exp(-(pow(r, alpha)/tau)));
	return f;
}

void SphericalAdvectionField::setOrigin(Vector3d o) {
	origin = o;
	return;
}

void SphericalAdvectionField::setRadius(double r) {
	radius = r;
	return;
}

void SphericalAdvectionField::setVMax(double v) {
	vMax = v;
	return;
}

void SphericalAdvectionField::setTau(double t) {
	tau = t;
	return;
}

void SphericalAdvectionField::setAlpha(double a) {
	alpha = a;
	return;
}

Vector3d SphericalAdvectionField::getOrigin() const {
	return origin;
}

double SphericalAdvectionField::getRadius() const {
	return radius;
}

double SphericalAdvectionField::getVMax() const {
	return vMax;
}

double SphericalAdvectionField::getTau() const {
	return tau;
}

double SphericalAdvectionField::getAlpha() const {
	return alpha;
}

std::string SphericalAdvectionField::getDescription() const {
	std::stringstream s;
	s << "Origin: " << origin / kpc  << " kpc, ";
	s << "Radius: " << radius / kpc  << " kpc, ";
	s << "vMax: " << vMax / km * sec << " km/s, ";
	s << "tau: " << tau << ", ";
	s << "alpha: " << alpha << "\n";
	return s.str();
}


//-----------------------------------------------------------------

SphericalAdvectionShock::SphericalAdvectionShock(Vector3d origin, double r_0, double v_0, double l) {
	setOrigin(origin);
	setR0(r_0);
	setV0(v_0);
	setLambda(l);
}


Vector3d SphericalAdvectionShock::getField(const Vector3d &pos) const {
	Vector3d R = pos-origin;
	Vector3d e_r = R.getUnitVector();
	double r = R.getR();

	double v = v_0 * ( 1 + (pow(r_0/(2*r), 2.) -1 ) * g(r));

	return v * e_r;
}


double SphericalAdvectionShock::getDivergence(const Vector3d &pos) const {
	double r = (pos-origin).getR();

	double d1 = 2./r*(1-g(r));
	double d2 = (pow(r_0/(2*r), 2.)-1)*g_prime(r);

	return v_0 * (d1+d2);
}


double SphericalAdvectionShock::g(double r) const {
	double a = (r-r_0)/lambda;
	return 1. / (1+exp(-a));
}

double SphericalAdvectionShock::g_prime(double r) const {
	double a = (r-r_0)/lambda;
	return 1. / (2*lambda*(1+cosh(-a)));
}	


void SphericalAdvectionShock::setOrigin(Vector3d o) {
	origin = o;
}

void SphericalAdvectionShock::setR0(double r) {
	r_0 = r;
}

void SphericalAdvectionShock::setV0(double v) {
	v_0 = v;
}

void SphericalAdvectionShock::setLambda(double l) {
	lambda = l;
}

Vector3d SphericalAdvectionShock::getOrigin() const {
	return origin;
}

double SphericalAdvectionShock::getR0() const {
	return r_0;
}

double SphericalAdvectionShock::getV0() const {
	return v_0;
}

double SphericalAdvectionShock::getLambda() const {
	return lambda;
}


} // namespace crpropa
