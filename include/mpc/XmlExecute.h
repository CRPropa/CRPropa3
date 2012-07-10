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
	std::vector<Source> sources;

	void loadCashKarp(pugi::xml_node &);
	void loadLssGrid(pugi::xml_node &);
	void loadSophia(pugi::xml_node &);
	void loadSpheresAroundObserver(pugi::xml_node &);
	void loadDiscreteSources(pugi::xml_node &);
	void loadFullTrajectoryOutput(pugi::xml_node &);

	size_t trajectories;
	double minEnergyEeV;
	double maxTimeMpc;
	size_t randomSeed;

public:
	bool load(const std::string &filename);
	void run();
};

} // namespace mpc

#endif /* MPC_XML_EXECUTE_H */
