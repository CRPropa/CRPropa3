#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

PeriodicMagneticField::PeriodicMagneticField(ref_ptr<MagneticField> field,
		const Vector3d &extends) :
		field(field), extends(extends), origin(0, 0, 0), reflective(false) {

}

PeriodicMagneticField::PeriodicMagneticField(ref_ptr<MagneticField> field,
		const Vector3d &extends, const Vector3d &origin, bool reflective) :
		field(field), extends(extends), origin(origin), reflective(reflective) {

}

Vector3d &PeriodicMagneticField::getOrigin() {
	return origin;
}

void PeriodicMagneticField::setOrigin(const Vector3d &origin) {
	this->origin = origin;
}

Vector3d &PeriodicMagneticField::getExtends() {
	return extends;
}

void PeriodicMagneticField::setExtends(const Vector3d &origin) {
	this->extends = extends;
}

bool PeriodicMagneticField::isReflective() {
	return reflective;
}

void PeriodicMagneticField::setReflective(bool reflective) {
	this->reflective = reflective;
}

Vector3d PeriodicMagneticField::getField(const Vector3d &position) const {
	Vector3d n = ((position - origin) / extends).floor();
	Vector3d p = position - origin - n * extends;

	if (reflective) {
		long mx = (long) ::fabs(n.x) % 2;
		if (mx == 1)
			p.x = extends.x - p.x;
		long my = (long) ::fabs(n.y) % 2;
		if (my == 1)
			p.y = extends.y - p.y;
		long mz = (long) ::fabs(n.z) % 2;
		if (mz == 1)
			p.z = extends.z - p.z;
	}

	return field->getField(p);
}

void MagneticFieldList::addField(ref_ptr<MagneticField> field) {
	fields.push_back(field);
}

Vector3d MagneticFieldList::getField(const Vector3d &position) const {
	Vector3d b;
	for (int i = 0; i < fields.size(); i++)
		b += fields[i]->getField(position);
	return b;
}

MagneticFieldEvolution::MagneticFieldEvolution(ref_ptr<MagneticField> field,
	double m) :
	field(field), m(m) {
}

Vector3d MagneticFieldEvolution::getField(const Vector3d &position,
	double z) const {
	return field->getField(position) * pow(1+z, m);
}

Vector3d MagneticDipoleField::getField(const Vector3d &position) const {
		Vector3d r = (position - origin);
		Vector3d unit_r = r.getUnitVector();
		
		if (r.getR() == 0) { // skip singularity
			return Vector3d(0,0,0);
		}
		return unit_r * moment.dot(unit_r) / pow((r.getR()/radius), 3) * mu0 / (4*M_PI);
}

Vector3d MagneticBottle::getField(const Vector3d &position) const { // Ulrich Stroth

	double x = position.x;
	double y = position.y;
	double z = position.z;

	if(fabs(z) <= L/2){
	
	double r = sqrt( x*x + y*y );
	double phi = atan2(y,x);
	
	
	double B0 = (Bz_max + Bz_min)/2;
	double b = Bz_max/ B0 -1;
	
	double k = 2 * M_PI /L;
	
	double Bz =  B0 * ( 1 - b * ( 1 + k*k * r*r / 4 ) * cos (k*z) );
	double Br = -B0 * ( b * k*r/2 * sin(k*z) );
	
	double Bx = Br * cos(phi);
	double By = Br * sin(phi);
	
	return Vector3d(Bx, By, Bz);
		
	}else{
	
	return Vector3d(0., 0., Bz_max);
		
	}
	
}

Vector3d GyroField::getField(const Vector3d &position) const {

	double x = position.x;
	double y = position.y;
	double z = position.z;
	
	double r = sqrt( x*x + y*y );
	
	double B;
	
	if (r < radius){
		B = B_max + r/radius * (B_min - B_max);
	}	
	else{
		B = B_min;
	}
	
	return Vector3d(0., 0., B);
	
	
}

Vector3d LongConductorField::getField(const Vector3d &position) const {

	double x = position.x;
	double y = position.y;
	double z = position.z;
	
	double r = sqrt( x*x + y*y );
	double phi = atan2(y,x);
	
	double B = radius/r * (B_radius);
	
	
	return Vector3d(-B*sin(phi), B*cos(phi), 0);
	
	
}

Vector3d CircleField::getField(const Vector3d &position) const {

	double x = position.x;
	double y = position.y;
	double z = position.z;
	
	double phi = atan2(y,x);
	
	return Vector3d(-B*sin(phi), B*cos(phi), 0);
	
	
}

Vector3d HongQinField::getField(const Vector3d &position) const {

	double x = position.x;
	double y = position.y;
	double z = position.z;
	
	double Bxy = 1. + delta * (x*x/4 +  y*y)/ (x*x+y*y);
	
	
	return Vector3d(0, 0, B * Bxy);
	
	
}


#ifdef CRPROPA_HAVE_MUPARSER
RenormalizeMagneticField::RenormalizeMagneticField(ref_ptr<MagneticField> field,
		std::string expression) :
		field(field), expression(expression) {

	p =  new mu::Parser();
	p->DefineVar("B", &Bmag);
	p->DefineConst("tesla", tesla);
	p->DefineConst("gauss", gauss);
	p->DefineConst("muG", muG);
	p->DefineConst("nG", nG);
	p->SetExpr(expression);
}

Vector3d RenormalizeMagneticField::getField(const Vector3d &position) {
	Vector3d B = field->getField(position);
	Bmag = B.getR();
	return B * p->Eval();
}
#endif

} // namespace crpropa
