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
	ModuleList modules;
	ref_ptr<MagneticField> magnetic_field;
	Source source;

	void loadUniformMagneticField(pugi::xml_node &node);
	void loadGridMagneticField(pugi::xml_node &node);

	void loadDeflectionCK(pugi::xml_node &node);

	void loadDiscreteSources(pugi::xml_node &node);

	void loadSophia(pugi::xml_node &node);

	void loadSpheresAroundObserver(pugi::xml_node &node);
	void loadSpheresAroundSource(pugi::xml_node &node);

	void loadOutput(pugi::xml_node &node);

	bool is1D;
	size_t trajectories;
	size_t randomSeed;
	double Emin;
	Vector3d origin;
	Vector3d size;

public:
	bool load(const std::string &filename);
	void run();
};

} // namespace mpc

#endif /* MPC_XML_EXECUTE_H */
