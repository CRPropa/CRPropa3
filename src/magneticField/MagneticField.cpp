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

JF12Field_disk::JF12Field_disk() {

	// spiral arm parameters
	pitch = 11.5 * M_PI / 180;
	sinPitch = sin(pitch);
	cosPitch = cos(pitch);
	tan90MinusPitch = tan(M_PI / 2 - pitch);

	rArms[0] = 5.1 * kpc;
	rArms[1] = 6.3 * kpc;
	rArms[2] = 7.1 * kpc;
	rArms[3] = 8.3 * kpc;
	rArms[4] = 9.8 * kpc;
	rArms[5] = 11.4 * kpc;
	rArms[6] = 12.7 * kpc;
	rArms[7] = 15.5 * kpc;

	// regular field parameters

	bRing = 0.1 * muG;
	hDisk = 0.40 * kpc;
	wDisk = 0.27 * kpc;

	bDisk[0] = 0.1 * muG;
	bDisk[1] = 3.0 * muG;
	bDisk[2] = -0.9 * muG;
	bDisk[3] = -0.8 * muG;
	bDisk[4] = -2.0 * muG;
	bDisk[5] = -4.2 * muG;
	bDisk[6] = 0.0 * muG;
	bDisk[7] = 2.7 * muG;
	
}

double logisticFunction_disk(double x, double x0, double w) {
	return 1. / (1. + exp(-2. * (fabs(x) - x0) / w));
}

Vector3d JF12Field_disk::getField(const Vector3d &pos) const {
	
	Vector3d b(0.);

	double r = sqrt(pos.x * pos.x + pos.y * pos.y); // in-plane radius
	double d = pos.getR(); // distance to galactic center
	if ((d < 1 * kpc) or (d > 20 * kpc))
		return b; // 0 field for d < 1 kpc or d > 20 kpc

	double phi = pos.getPhi(); // azimuth
	double sinPhi = sin(phi);
	double cosPhi = cos(phi);

	double lfDisk = logisticFunction_disk(pos.z, hDisk, wDisk);
	
	if (r > 3 * kpc) {
		double bMag;
		if (r < 5 * kpc) {
			// molecular ring
			bMag = bRing * (5 * kpc / r) * (1 - lfDisk);
			b.x += -bMag * sinPhi;
			b.y += bMag * cosPhi;

		} else {
			// spiral region
			double r_negx = r * exp(-(phi - M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + M_PI) / tan90MinusPitch);
			if (r_negx > rArms[7])
				r_negx = r * exp(-(phi + 3 * M_PI) / tan90MinusPitch);

			for (int i = 7; i >= 0; i--)
				if (r_negx < rArms[i])
					bMag = bDisk[i];

			bMag *= (5 * kpc / r) * (1 - lfDisk);
			b.x += bMag * (sinPitch * cosPhi - cosPitch * sinPhi);
			b.y += bMag * (sinPitch * sinPhi + cosPitch * cosPhi);
		}
	}

	return b;
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
