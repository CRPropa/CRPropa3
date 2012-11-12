#ifndef MPC_XML_EXECUTE_H
#define MPC_XML_EXECUTE_H

#include "mpc/ModuleList.h"
#include "mpc/Source.h"
#include "mpc/magneticField/MagneticField.h"

#include <vector>

namespace pugi {
class xml_node;
}

namespace mpc {

class XmlExecute {
	void loadUniformMagneticField(pugi::xml_node &node);
	void loadGridMagneticField(pugi::xml_node &node);
	void loadDeflectionCK(pugi::xml_node &node);
	void loadPeriodicBoundaries();
	void loadDiscreteSources(pugi::xml_node &node);
	void loadContinuousSources(pugi::xml_node &node);
	void loadSpectrumComposition(pugi::xml_node &node);
	void loadSourceNuclei(pugi::xml_node &node);
	void loadSophia(pugi::xml_node &node);
	void loadSpheresAroundObserver(pugi::xml_node &node);
	void loadSpheresAroundSource(pugi::xml_node &node);
	void loadOutput(pugi::xml_node &node);

	ModuleList modules;
	ref_ptr<MagneticField> magnetic_field;
	Source source;
	bool is1D;
	size_t nTrajectories;
	double Emin;
	double maxStep;
	Vector3d origin;
	Vector3d size;

public:
	XmlExecute();
	bool load(const std::string &filename);
	void run();
};

} // namespace mpc

#endif /* MPC_XML_EXECUTE_H */
