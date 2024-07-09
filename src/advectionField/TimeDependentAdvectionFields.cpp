#include "crpropa/advectionField/TimeDependentAdvectionFields.h"


namespace crpropa {

OneDimensionalTimeDependentShock::OneDimensionalTimeDependentShock(double v_sh, double v1, double v0, double l_sh){

	setVshock(v_sh);
	setSpeeds(v1, v0);
	setShockWidth(l_sh);
	setShockPosition(0.); 
	setShockTime(0.); 
	}

Vector3d OneDimensionalTimeDependentShock::getField(const Vector3d &position, const double &time) const {

    double A = 0.5 * ( v1 + v0 ); 
    double B = 0.5 * ( v1 - v0 );

	double x = position.x;
    double xsh = v_sh * (time-t_sh0) + x_sh0; // shock position 

	Vector3d v(0.);
    v.x = A - B * tanh( (x - xsh) / l_sh );

	return v;
	}

double OneDimensionalTimeDependentShock::getDivergence(const Vector3d &position, const double &time) const {

	double dvdx = 0.;
	double x = position.x;

	if (time < t_sh0){
        // no shock
		dvdx = 0;
	}
    else{

        double B = 0.5 * ( v1 - v0 );
	    double xsh = v_sh * (time-t_sh0) + x_sh0; // shock position 
        dvdx = - B / l_sh * ( 1 - tanh( (x - xsh) / l_sh ) * tanh( (x - xsh) / l_sh ) );
    }
	
	return dvdx;
	}


void OneDimensionalTimeDependentShock::setVshock(double v) {
	v_sh = v;
	return;
}
void OneDimensionalTimeDependentShock::setSpeeds(double u1, double u0) {
	v1 = u1;
	v0 = u0;
	return;
}
void OneDimensionalTimeDependentShock::setShockWidth(double w) {
	l_sh = w;
	return;
}

void OneDimensionalTimeDependentShock::setShockPosition(double x) {
	x_sh0 = x;
	return;
}

void OneDimensionalTimeDependentShock::setShockTime(double t) {
	t_sh0 = t;
	return;
}

double OneDimensionalTimeDependentShock::getVshock() const {
    return v_sh;
}

double OneDimensionalTimeDependentShock::getV1() const {
    return v1;
}

double OneDimensionalTimeDependentShock::getV0() const {
    return v0;
}

double OneDimensionalTimeDependentShock::getShockWidth() const {
    return l_sh;
}

double OneDimensionalTimeDependentShock::getShockPosition() const {
    return x_sh0;
}

double OneDimensionalTimeDependentShock::getShockTime() const {
    return t_sh0;
}

//----------------------------------------------------------------

SedovTaylorBlastWave::SedovTaylorBlastWave(double E0, double rho0, double l_sh){

	setEnergy(E0);
	setDensity(rho0);
	setShockWidth(l_sh);
	
	}

Vector3d SedovTaylorBlastWave::getField(const Vector3d &position, const double &time) const {

	double r = position.getR();
	Vector3d e_r = position.getUnitVector();

	if (time == 0)
		return 0. * e_r ;
	double A = rho0;

    double R = pow(E0 / A, 1./5.) * pow(time, 2. / 5. ); // position of shock front 
    double vs = 2. / 5. * pow(E0 / A, 1. / 5.) * pow(time, -3. / 5.); // shock speed

    double xi = r / R; // dimensionless radius
    double V = 3 * (pow(xi, 8) + 1) / (3 * pow(xi, 8) + 5); 
    double u = 0.5 * xi * vs * V * ( 1 - tanh( (xi - 1) * R / l_sh) ); 

	return u * e_r;

	}

double SedovTaylorBlastWave::getDivergence(const Vector3d &position, const double &time) const {

	if (time == 0)
		return 0;
	double r = position.getR();
	double A = rho0; 

	double R = pow(E0 / A, 1./5.) * pow(time, 2. / 5.); 
    double vs = 2. / 5. * pow(E0 / A, 1. / 5.) * pow(time, -3. / 5.); 

	double a = (-3 * pow(r,3) * (1 + pow(r/R, 8)) * 1./pow(cosh((r - R)/l_sh),2)) / (l_sh * (5 + (3 * pow(r,8))/pow(R,8)))
     + (9 * pow(r,2) * (1 + pow(r/R, 8)) * (1 - tanh((r - R)/l_sh))) / (5 + (3 * pow(r,8))/pow(R,8))
     - (72 *pow(r,10) * (1 + pow(r/R,8)) * (1 - tanh((r - R)/l_sh)))/(pow(5 + (3 * pow(r,8))/pow(R,8),2) * pow(R,8)) 
     + (24 * pow(r,10) * (1 - tanh((r - R)/l_sh)))/((5 + (3 * pow(r,8))/pow(R,8)) * pow(R,8));
	
    double dudr = 0.5 * vs * pow(r, -2) * 1./ R * a;

	return dudr;

}


void SedovTaylorBlastWave::setShockWidth(double w) {
	l_sh = w;
	return;
}

void SedovTaylorBlastWave::setDensity(double rho) {
	rho0 = rho;
	return;
}

void SedovTaylorBlastWave::setEnergy(double e) {
	E0 = e;
	return;
}

double SedovTaylorBlastWave::getShockRadius(double time) const{
	double A = rho0;
	return  pow(E0 / A, 1./5.) * pow(time, 2./5.); 
}

double SedovTaylorBlastWave::getShockSpeed(double time) const{
    double A = rho0;
    return 2. / 5. * pow(E0 / A, 1. / 5.) * pow(time, -3. / 5.);
}

double SedovTaylorBlastWave::getShockWidth() const {
    return l_sh;
}

double SedovTaylorBlastWave::getDensity() const {
    return rho0;
}

double SedovTaylorBlastWave::getEnergy() const {
    return E0;
}

//----------------------------------------------------------------

} // namespace crpropa
