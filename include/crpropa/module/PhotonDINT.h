#ifndef CRPROPA_PHOTONDINT_H
#define CRPROPA_PHOTONDINT_H

#include "crpropa/Module.h"
#include "crpropa/magneticField/MagneticField.h"

#include <fstream>

namespace crpropa {

class PhotonDINT: public Module {
private:
	std::string filename, dataPath;
	mutable std::ofstream fout;
	ref_ptr<MagneticField> field;

	int IRFlag, RadioFlag;
	double Zmax, Cutcascade_Magfield;

public:
	PhotonDINT(const std::string &filename, ref_ptr<MagneticField> field);
	~PhotonDINT();
	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_PHOTONDINT_H
