#ifndef CRPROPA_PHOTONDINT1D_H
#define CRPROPA_PHOTONDINT1D_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"

#include <fstream>

namespace crpropa {

class PhotonDINT1DImpl;

class PhotonDINT1D: public Module {
private:
	std::string filename, dataPath;
	int IRFlag, RadioFlag;
	double Zmax, Cutcascade_Magfield;

	mutable PhotonDINT1DImpl *impl;
public:
	PhotonDINT1D(const std::string &filename);
	~PhotonDINT1D();
	void process(Candidate *candidate) const;
	std::string getDescription() const;

	/* Possible values:
	 * 0 -> IR_UNIFORM_HIGH (default)
	 * 1 -> IR_UNIFORM_LOW
	 * 2 -> IR_UNIFORM_PRIMACK
	 */
	void setIRFlag(int ir);

	/* Possible values:
	 * 0 -> High (default)
	 * 1 -> Med
	 * 2 -> Obs
	 * 3 -> Null
	 */
	void setRadioFlag(int radio);

	void setZmax(double zmax);
};

} // namespace crpropa

#endif // CRPROPA_PHOTONDINT1D_H
