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


		if (r.getR() == 0) { // singularity
			return moment * 2 * mu0 / 3;
		}
		return (unit_r * (unit_r.dot(moment)) * 3 - moment) / pow(r.getR() / radius, 3) * mu0 / (4*M_PI);
}

Vector3d SingleModeHelicalMagneticField::getField(const Vector3d &position) const {
		Vector3d r = (position - origin);

		Vector3d e1 = unitVectorOrigin; // first of the unit vectors which span the polarization plane
		if ((e1.x != 0.0) && (e1.y != 0.0) && (e1.z != 0.0)) {
			e1 /= e1.getR();
		}

		Vector3d e2 = unitVector2; // second of the unit vectors which span the polarization plane
		if ((e2.x != 0.0) && (e2.y != 0.0) && (e2.z != 0.0)) {
			e2 /= e2.getR();
		}

		Vector3d wavevector = e1.cross(e2);
		if ((wavevector.x != 0.0) && (wavevector.y != 0.0) && (wavevector.z != 0.0)) {
			wavevector = wavevector/wavevector.getR() * 2 * M_PI / wavelength;
		}

		return amplitudeOrigin * (e1 * cos(wavevector.dot(r)) + e2 * handedness * sin(wavevector.dot(r)));
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
