#include "mpc/XmlExecute.h"
#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/Output.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/ModuleList.h"

#include "pugixml.hpp"

#include <fstream>

using namespace pugi;
using namespace std;

namespace mpc {

bool XmlExecute::load(const string &filename) {
	xml_document doc;
	xml_parse_result result = doc.load_file(filename.c_str());

	if (!result) {
		cout << "Error description: " << result.description() << "\n";
		cout << "Error offset: " << result.offset << "\n";
		return false;
	}

	xml_node root = doc.child("CRPropa");

	trajectories = root.child("TrajNumber").attribute("value").as_int();
	cout << "Trajectories: " << trajectories << endl;

	randomSeed = root.child("RandomSeed").attribute("value").as_int();
	cout << "RandomSeed: " << randomSeed << endl;

	// setup MagneticField
	xml_node field = root.child("MagneticField");
	if (field) {
		string type = field.attribute("type").as_string();
		cout << "MagenticField: " << type << endl;
		if (type == "LSS-Grid") {
			loadLssGrid(field);
		} else {
			cout << " --> unknown, set zero field" << endl;
			magnetic_field = new UniformMagneticField(Vector3f(0, 0, 0));
		}
	}

	// setup Integrator
	xml_node integrator = root.child("Integrator");
	if (integrator) {
		string type = integrator.attribute("type").as_string();
		cout << "Integrator: " << type << endl;
		if (type == "Cash-Karp RK")
			loadCashKarp(integrator);
		else
			cout << "  -> unknown" << endl;
	}

	// setup Interactions
	xml_node interactions = root.child("Interactions");
	if (interactions) {
		string type = interactions.attribute("type").as_string();
		cout << "Interactions: " << type << endl;
		if (type == "Sophia")
			loadSophia(interactions);
		else
			cout << "  -> unknown" << endl;
	}

	// observers
	xml_node observers = root.child("Observers");
	if (observers) {
		string type = observers.attribute("type").as_string();
		cout << "Observers: " << type << endl;
		if (type == "Spheres around Observers")
			loadSpheresAroundObserver(observers);
		else
			cout << "  -> unknown" << endl;
	}

	// output
	xml_node output = root.child("Output");
	if (output) {
		string type = output.attribute("type").as_string();
		cout << "Output: " << type << endl;
		if (type == "Full Trajectories")
			loadFullTrajectoryOutput(output);
		else
			cout << "  -> unknown" << endl;
	}

	// sources
	xml_node sources = root.child("Sources");
	if (sources) {
		string type = sources.attribute("type").as_string();
		cout << "Sources: " << type << endl;
		if (type == "Discrete")
			loadDiscreteSources(sources);
		else
			cout << "  -> unknown" << endl;
	}

	minEnergyEeV = root.child("MinEnergy_EeV").attribute("value").as_double()
			* EeV;
	cout << "Minimum Energy: " << minEnergyEeV / EeV << " EeV" << endl;
	modules.add(new MinimumEnergy(minEnergyEeV));

	maxTimeMpc = root.child("MaxTime_Mpc").attribute("value").as_double() * Mpc;
	cout << "Maximum Time: " << maxTimeMpc / Mpc << " Mpc" << endl;
	modules.add(new MaximumTrajectoryLength(maxTimeMpc));

	return true;
}

void XmlExecute::loadCashKarp(xml_node &node) {
	double epsilon = node.child("Epsilon").attribute("value").as_double();
	cout << "  - Epsilon: " << epsilon << endl;

	double minstep = node.child("MinStep_Mpc").attribute("value").as_double()
			* Mpc;
	cout << "  - Minimum Step: " << minstep / Mpc << " Mpc " << endl;

	ref_ptr<DeflectionCK> deflection = new DeflectionCK(magnetic_field, epsilon,
			minstep);
	modules.add(deflection);
}

void XmlExecute::loadLssGrid(xml_node &node) {
	int nx = node.child("Nx").attribute("value").as_int();
	int ny = node.child("Ny").attribute("value").as_int();
	int nz = node.child("Nz").attribute("value").as_int();
	if (nx != ny || nx != nz)
		throw runtime_error("Invalid grid size!");
	cout << "  - Samples: " << nx << endl;

	double step = node.child("Step_Mpc").attribute("value").as_double() * Mpc;
	double size = nx * step;
	cout << "  - Size: " << size / Mpc << " Mpc " << endl;

	xml_node origin_node = node.child("Origin");
	Vector3d origin;
	origin.x = origin_node.child("X_Mpc").attribute("value").as_double();
	origin.y = origin_node.child("Y_Mpc").attribute("value").as_double();
	origin.z = origin_node.child("Z_Mpc").attribute("value").as_double();
	cout << "  - Origin: " << origin << endl;

	magnetic_field = new MagneticFieldGrid(origin, size, nx);
}

void XmlExecute::loadSophia(xml_node &) {

}

void XmlExecute::loadSpheresAroundObserver(xml_node &) {

}

void XmlExecute::loadDiscreteSources(xml_node &) {

}

void XmlExecute::loadFullTrajectoryOutput(xml_node &node) {
	string filename = node.child("File").child_value();
	cout << "  - Filename: " << filename << endl;

	string option = node.child("File").attribute("option").as_string();
	if (option != "force") {
		ifstream ifile(filename.c_str());
		if (ifile) {
			throw runtime_error("Outputfile already exists!");
		}
	}

	ref_ptr<TrajectoryOutput> output = new TrajectoryOutput(filename);
	modules.add(output);
}

void XmlExecute::run() {
	//operator <<(cout, modules);
	ParticleState initial;
	initial.setId(getNucleusId(1, 1));
	initial.setEnergy(100 * EeV);
	initial.setPosition(Vector3d(0, 0, 0));
	initial.setDirection(Vector3d(1, 0, 0));

	for (size_t i = 0; i < trajectories; i++) {
		ref_ptr<Candidate> candidate = new Candidate(initial);
		if (i % (trajectories / 10) == 0)
			std::cout << i << std::endl;
		modules.run(candidate);
	}

}

} // namespace mpc
