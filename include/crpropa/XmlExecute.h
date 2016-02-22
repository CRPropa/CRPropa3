#ifndef CRPROPA_XMLEXECUTE_H
#define CRPROPA_XMLEXECUTE_H

#include "crpropa/ModuleList.h"
#include "crpropa/Source.h"
#include "crpropa/module/Observer.h"
#include "crpropa/magneticField/MagneticField.h"

namespace pugi {
class xml_node;
}

namespace crpropa {

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

	ref_ptr<MagneticField> magnetic_field;
	ModuleList modules;
	Source source;
	Observer observer;
	bool is1D;
	bool hasRedshift;
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

} // namespace crpropa

#endif // CRPROPA_XMLEXECUTE_H
