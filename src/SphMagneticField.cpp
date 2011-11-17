#include "MagneticField.h"
#include <stdint.h>
#include <sstream>
#include <unistd.h>

SampledSphMagneticField::SampledSphMagneticField(Hep3Vector &origin, double &size,
		double &stepSize, std::string &sphFileName) {
	this->origin = origin;
	this->size = size;
	this->gridStepSize = gridStepSize;
	initialize(sphFileName);
}

SampledSphMagneticField::~SampledSphMagneticField() {
}

void SampledSphMagneticField::initialize(std::string sphFileName) {
	std::cout << "[TSphMagField] initializing field" << std::endl;
	Vector3f originKpc = Vector3f(_fOriginKpc.x(), _fOriginKpc.y(),
			_fOriginKpc.z());
	_fSampledMagneticField.reset(
			new SampledMagneticField(originKpc, _fSizeKpc));
	_fSampledMagneticField->init(_fStepSizeKpc);

	std::string s = _fFileName.substr(_fFileName.find_last_of('/') + 1);
	s.replace(s.rfind(".raw"), 4, "");
	stringstream ss;
	ss << "dump/" << s << "-" << _fStepSizeKpc << ".dump";
	std::cout << ss << std::endl;
	std::string dumpfilename = ss.str();

	std::cout << "[TSphMagField] checking previous field dump: " << dumpfilename
			<< std::endl;
	if (_fSampledMagneticField->restore(dumpfilename) == true)
		std::cout << "[TSphMagField] Successful" << std::endl;
	else {
		std::cout << "[TSphMagField] .. unsuccessful" << std::endl;
		std::cout << "[TSphMagField] loading field" << std::endl;
		_fSampledMagneticField->load(_fFileName);

		std::cout << "[TSphMagField] dumping field" << std::endl;
		std::stringstream tmpfilename;
		tmpfilename << dumpfilename << getpid();
		_fSampledMagneticField->dump(tmpfilename.str());
		::rename(tmpfilename.str().c_str(), dumpfilename.c_str());
	}
}

void TSphMagField::createDirectField() {
	std::cout << "[TSphMagField] initializing SPH field" << std::endl;
	Vector3f originKpc = Vector3f(_fOriginKpc.x(), _fOriginKpc.y(),
			_fOriginKpc.z());
	_fDirectMagneticField.reset(new DirectMagneticField(originKpc, _fSizeKpc));
	_fDirectMagneticField->init(_fStepSizeKpc);
	SmoothParticleHelper::read(_fFileName, _fSmoothParticles);
	_fDirectMagneticField->load(_fSmoothParticles);
}

TVector3D TSphMagField::getField(TVector3D position) const {
	Vector3f r = Vector3f(position.x(), position.y(), position.z());
	Vector3f b;
	if (_fFieldType == "Direct") {
		b = _fDirectMagneticField->getField(r/kpc);
	} else {
		b = _fSampledMagneticField->getField(r/kpc);
	}
	// SPH field strength is in Gauss
	TVector3D bField = TVector3D(b.x, b.y, b.z) * gauss;
//	std::cout << "[TSphMagField] field strength at " << position/Mpc << " Mpc, " << bField.mag() / gauss << " Gauss" << std::endl;
	return bField;
}




