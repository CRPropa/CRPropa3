#include "mpc/XmlExecute.h"
#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/Output.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/NuclearDecay.h"
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
	}

	xml_node root = doc.child("CRPropa");
	xml_node node;
	std::string type;

	trajectories = root.child("TrajNumber").attribute("value").as_int();
	cout << "Trajectories: " << trajectories << endl;

	double maxTime = root.child("MaxTime_Mpc").attribute("value").as_double() * Mpc;
	cout << "Maximum Time: " << maxTime / Mpc << " Mpc" << endl;

	Emin = root.child("MinEnergy_EeV").attribute("value").as_double() * EeV;
	cout << "Minimum Energy: " << Emin / EeV << " EeV" << endl;

	randomSeed = root.child("RandomSeed").attribute("value").as_int();
	cout << "RandomSeed: " << randomSeed << endl;

	// environment
	node = root.child("Environment");
	if (!node)
		throw runtime_error("Environment not specified");
	type = node.attribute("type").as_string();
	cout << "Environment: " << type << endl;
	if (type == "One Dimension")
		throw runtime_error("One Dimension not supported");
	else if (type == "LSS") {
		// will be overwritten if (re)defined by magnetic field
		origin.x = node.child("Xmin_Mpc").attribute("value").as_double() * Mpc;
		origin.y = node.child("Ymin_Mpc").attribute("value").as_double() * Mpc;
		origin.z = node.child("Zmin_Mpc").attribute("value").as_double() * Mpc;
		size.x = node.child("Xmax_Mpc").attribute("value").as_double() * Mpc;
		size.y = node.child("Ymax_Mpc").attribute("value").as_double() * Mpc;
		size.z = node.child("Zmax_Mpc").attribute("value").as_double() * Mpc;
		size -= origin;
	} else
		throw runtime_error("Unknown environment");

	// magnetic field
	node = root.child("MagneticField");
	if (!node)
		throw runtime_error("Magnetic field not specified");
	type = node.attribute("type").as_string();
	cout << "MagenticField: " << type << endl;
	if (type == "Uniform")
		loadUniformMagneticField(node);
	else if ((type == "LSS") or (type == "Kolmogoroff"))
		loadGridMagneticField(node);
	else {
		cout << " --> unknown, set zero field" << endl;
		magnetic_field = new UniformMagneticField(Vector3d(0, 0, 0));
	}

	// propagator
	node = root.child("Integrator");
	if (!node)
		throw runtime_error("Integrator not specified");
	type = node.attribute("type").as_string();
	cout << "Integrator: " << type << endl;
	if (type == "Cash-Karp RK")
		loadDeflectionCK(node);
	else
		throw runtime_error("Unknown integrator");

	// interactions
	node = root.child("Interactions");
	if (!node)
		throw runtime_error("Interactions not specified");
	type = node.attribute("type").as_string();
	cout << "Interactions: " << type << endl;
	if (type == "None")
		;
	else if (type == "F77-proton")
		cout << "  -> not supported" << endl;
	else if (type == "Sophia")
		loadSophia(node);
	else if (type == "Photon")
		cout << "  -> not supported" << endl;
	else
		cout << "  -> unknown" << endl;

	// minimum energy
	modules.add(new MinimumEnergy(Emin));

	// maximum trajectory length
	modules.add(new MaximumTrajectoryLength(maxTime));

	// periodic boundaries
	cout << "Periodic boundaries" << endl;
	cout << "  - Lower bounds: " << origin / Mpc << " Mpc" << endl;
	cout << "  - Upper bounds: " << (origin + size) / Mpc << " Mpc" << endl;
	modules.add(new PeriodicBox(origin, size));

	// sources
	node = root.child("Sources");
	if (!node)
		throw runtime_error("Source(s) not specified");
	type = node.attribute("type").as_string();
	cout << "Sources: " << type << endl;
	if (type == "Discrete")
		loadDiscreteSources(node);
	else if (type == "Continuous")
		throw runtime_error("Continuous sources not implemented");
	else
		throw runtime_error("Unknown source");

	// observers
	node = root.child("Observers");
	if (!node)
		cout << "Observer(s) not specified" << endl;
	else {
		string type = node.attribute("type").as_string();
		cout << "Observers: " << type << endl;
		if (type == "Spheres around Observers")
			loadSpheresAroundObserver(node);
		else if (type == "Spheres around Source")
			loadSpheresAroundSource(node);
		else
			cout << "  -> unknown observer" << endl;
	}

	// output
	node = root.child("Output");
	if (!node)
		cout << "No output" << endl;
	else
		loadOutput(node);

	return true;
}

void XmlExecute::loadDeflectionCK(xml_node &node) {
	double epsilon = node.child("Epsilon").attribute("value").as_double();
	cout << "  - Epsilon: " << epsilon << endl;
	double minstep = node.child("MinStep_Mpc").attribute("value").as_double()
			* Mpc;
	cout << "  - Minimum Step: " << minstep / Mpc << " Mpc " << endl;
	modules.add(new DeflectionCK(magnetic_field, epsilon, minstep));
}

void XmlExecute::loadUniformMagneticField(xml_node &node) {
	Vector3d bField;
	bField.x = node.child("Bx_nG").attribute("value").as_int() * nG;
	bField.y = node.child("By_nG").attribute("value").as_int() * nG;
	bField.z = node.child("Bz_nG").attribute("value").as_int() * nG;
	cout << "  - Bx, By, Bz: " << bField / nG << " nG" << endl;
	magnetic_field = new UniformMagneticField(bField);
}

void XmlExecute::loadGridMagneticField(xml_node &node) {
	xml_node origin_node = node.child("Origin");
	origin.x = origin_node.child("X_Mpc").attribute("value").as_double() * Mpc;
	origin.y = origin_node.child("Y_Mpc").attribute("value").as_double() * Mpc;
	origin.z = origin_node.child("Z_Mpc").attribute("value").as_double() * Mpc;
	cout << "  - Origin: " << origin / Mpc << endl;

	int Nx = node.child("Nx").attribute("value").as_int();
	int Ny = node.child("Ny").attribute("value").as_int();
	int Nz = node.child("Nz").attribute("value").as_int();
	cout << "  - Samples: " << Nx << " * " << Ny << " * " << Nz << endl;

	double spacing = node.child("Step_Mpc").attribute("value").as_double() * Mpc;
	cout << "  - Spacing: " << spacing / Mpc << " Mpc " << endl;

	size.x = Nx * spacing;
	size.y = Ny * spacing;
	size.z = Nz * spacing;

	ref_ptr<VectorGrid> field = new VectorGrid(origin, Nx, Ny, Nz, spacing);

	std::string type = node.attribute("type").as_string();
	if (type == "LSS-Grid") {
		string filetype = node.child("File").attribute("type").as_string();
		cout << "  - File type: " << filetype << endl;
		if (filetype != "ASCII")
			throw runtime_error("Only ASCII files supported");

		string fname = node.child_value("File");
		cout << "  - Loading file (values in [G]): " << fname << endl;
		loadTxt(field, fname, gauss);
	} else {
		double brms = node.child("RMS_muG").attribute("value").as_double() * 1e-6 * gauss;
		cout << "  - Brms : " << brms / nG << " nG" << endl;

		double kMin = node.child("Kmin").attribute("value").as_double();
		double kMax = node.child("Kmax").attribute("value").as_double();
		double lMin = spacing / kMax;
		double lMax = spacing / kMin;
		cout << "  - Turbulent range: " << lMin / Mpc << " - " << lMax / Mpc << " Mpc" << endl;

		double alpha = node.child("SpectralIndex").attribute("value").as_double();
		cout << "  - Turbulence spectral index: " << alpha << endl;
		cout << "  - Random seed: " << randomSeed << endl;

		initTurbulence(field, brms, lMin, lMax, alpha, randomSeed);
	}

	magnetic_field = new MagneticFieldGrid(field);
}

void XmlExecute::loadSophia(xml_node &node) {
	if (node.child("NoRedshift"))
		cout << "  - No redshift" << endl;
//	else
//		redshift

	if (node.child("NoPairProd"))
		cout << "  - No pair production" << endl;
	else
		modules.add(new ElectronPairProduction(CMB_IRB));

	if (node.child("NoPionProd"))
		cout << "  - No pion production" << endl;
	else {
		modules.add(new ElectronPairProduction(CMB));
		if (!node.child("NoIRPionProd"))
			modules.add(new ElectronPairProduction(IRB));
	}

	if (node.child("NoPhotoDisintegration"))
		cout << "  - No photo disintegration" << endl;
	else {
		modules.add(new PhotoDisintegration(CMB));
		modules.add(new PhotoDisintegration(IRB));
	}

	if (node.child("NoDecay"))
		cout << "  - No decay" << endl;
	else
		modules.add(new NuclearDecay);
}

void XmlExecute::loadSpheresAroundObserver(xml_node &node) {
	double r = node.child("Radius_Mpc").attribute("value").as_double() * Mpc;
	cout << "  - Radius: " << r / Mpc << " Mpc" << endl;

	int nObs = 0;
	for (xml_node n = node.child("SphereObserver"); n; n = n.next_sibling("SphereObserver")) {
		nObs += 1;
		Vector3d pos;
		pos.x = n.child("CoordX_Mpc").attribute("value").as_double() * Mpc;
		pos.y = n.child("CoordY_Mpc").attribute("value").as_double() * Mpc;
		pos.z = n.child("CoordZ_Mpc").attribute("value").as_double() * Mpc;
		cout << "  - Postion: " << pos / Mpc << " Mpc" << endl;
		modules.add(new SmallObserverSphere(pos, r));
	}
	if (nObs > 1)
		cout << "  -> Warning for multiple observers: propagation will stop after first detection." << endl;
}

void XmlExecute::loadSpheresAroundSource(pugi::xml_node &node) {
	int nObs = 0;
	for (xml_node n = node.child("Sphere"); n; n = n.next_sibling("Sphere")) {
		nObs += 1;
		Vector3d pos;
		pos.x = n.child("CoordX_Mpc").attribute("value").as_double() * Mpc;
		pos.y = n.child("CoordY_Mpc").attribute("value").as_double() * Mpc;
		pos.z = n.child("CoordZ_Mpc").attribute("value").as_double() * Mpc;
		double r = n.child("Radius_Mpc").attribute("value").as_double() * Mpc;
		cout << "  - Postion: " << pos / Mpc << " Mpc" << endl;
		cout << "  - Radius: " << r / Mpc << " Mpc" << endl;
		SphericalBoundary *sphere = new SphericalBoundary(pos, r, "Detected");
		modules.add(sphere);
	}
	if (nObs > 1)
		cout << "  -> Warning for multiple observers: propagation will stop after first detection." << endl;
}

void XmlExecute::loadDiscreteSources(xml_node &node) {
	// source positions
	SourceMultiplePositions *positions = new SourceMultiplePositions();
	for (xml_node n = node.child("PointSource"); n; n = n.next_sibling("PointSource")) {
		Vector3d pos;
		pos.x = n.child("CoordX_Mpc").attribute("value").as_double() * Mpc;
		pos.y = n.child("CoordY_Mpc").attribute("value").as_double() * Mpc;
		pos.z = n.child("CoordZ_Mpc").attribute("value").as_double() * Mpc;
		cout << "  - Position " << pos / Mpc << " Mpc" << endl;
		positions->add(pos);
	}
	source.addProperty(positions);

	// source spectrum
	xml_node spec = node.child("Spectrum");
	string spectrumType = spec.attribute("type").as_string();
	cout << "  - Spectrum: " << spectrumType << endl;

	if (spectrumType == "Monochromatic") {
		double E = spec.child("Energy_EeV").attribute("value").as_double() * EeV;
		source.addProperty(new SourceEnergy(E));
		cout << "  - Energy: " << E / EeV << " EeV" << endl;

		// source composition
		SourceNuclei *composition = new SourceNuclei;
		xml_node p = node.child("Particles");
		for (xml_node n = p.child("Species"); n; n = n.next_sibling("Species")) {
			int A = n.child("MassNumber").attribute("value").as_int();
			int Z = n.child("ChargeNumber").attribute("value").as_int();
			double ab = n.child("Abundance").attribute("value").as_double();
			composition->add(getNucleusId(A, Z), ab);
			cout << "  - Species: Z = " << Z << ", A = " << A << ", abundance = " << ab <<  endl;
		}
		source.addProperty(composition);

	} else if (spectrumType == "Power Law") {
		double alpha = spec.child("Alpha").attribute("value").as_double();
		double Rmax = spec.child("Rigidity_EeV").attribute("value").as_double() * EeV;
		cout << "  - Minimum energy: " << Emin / EeV << " EeV" << endl;
		cout << "  - Maximum rigidity: " << Rmax / EeV << " EeV" << endl;
		cout << "  - Power law index: " << alpha << endl;

		// source composition
		SourceComposition *composition = new SourceComposition(Emin, Rmax, alpha);
		xml_node p = node.child("Particles");
		for (xml_node n = p.child("Species"); n; n = n.next_sibling("Species")) {
			int A = n.attribute("MassNumber").as_int();
			int Z = n.attribute("ChargeNumber").as_int();
			double ab = n.attribute("Abundance").as_double();
			cout << "  - Species: Z = " << Z << ", A = " << A << ", abundance = " << ab <<  endl;
			composition->add(getNucleusId(A, Z), ab);
		}
		source.addProperty(composition);

	} else {
		throw runtime_error(" --> unknown source");
	}
}

void XmlExecute::loadOutput(xml_node &node) {
	string type = node.attribute("type").as_string();
	cout << "Output: " << type << endl;

	string filename = node.child("File").child_value();
	cout << "  - Filename: " << filename << endl;

	string option = node.child("File").attribute("option").as_string();
	if (option != "force") {
		ifstream ifile(filename.c_str());
		if (ifile) {
			throw runtime_error("Outputfile already exists!");
		}
	}

	if (type == "Full Trajectories")
		modules.add(new CRPropa2TrajectoryOutput(filename));
	else if (type == "Events")
		modules.add(new CRPropa2EventOutput(filename));
	else if (type == "None")
		return;
	else
		cout << "  -> unknown output" << endl;
}

void XmlExecute::run() {
	//operator <<(cout, modules);
	modules.setShowProgress(true);
	modules.run(&source, trajectories, true);
}

} // namespace mpc
