#include "mpc/XmlExecute.h"
#include "mpc/magneticField/UniformMagneticField.h"
#include "mpc/magneticField/MagneticFieldGrid.h"
#include "mpc/PeriodicGridTools.h"
#include "mpc/module/SimplePropagation.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoPionProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/Boundary.h"
#include "mpc/module/Observer.h"
#include "mpc/module/Output.h"
#include "mpc/ModuleList.h"

#include "pugixml.hpp"

#include <fstream>
#include <sstream>

using namespace pugi;
using namespace std;

namespace mpc {

double childValue(xml_node parent, string name, bool throwIfEmpty = true) {
	xml_node node = parent.child(name.c_str());
	if (!node and throwIfEmpty) {
		stringstream ss;
		ss << "Error reading XML card: " << name << " not found";
		throw runtime_error(ss.str());
	}
	return node.attribute("value").as_double();
}

string childType(xml_node parent, string name, bool throwIfEmpty = true) {
	xml_node node = parent.child(name.c_str());
	if (!node and throwIfEmpty) {
		stringstream ss;
		ss << "Error reading XML card: " << name << " not found";
		throw runtime_error(ss.str());
	}
	return node.attribute("type").as_string();
}

bool XmlExecute::load(const string &filename) {
	xml_document doc;
	xml_parse_result result = doc.load_file(filename.c_str());

	if (!result) {
		cout << "Error reading XML card: " << result.description() << "\n";
		cout << "Error position: " << result.offset << "\n";
		return false;
	}

	xml_node root = doc.child("CRPropa");
	if (!root)
		throw runtime_error("Error reading XML card: Root element not found");

	trajectories = (int)childValue(root, "TrajNumber");
	cout << "Trajectories: " << trajectories << endl;

	double maxTime = childValue(root, "MaxTime_Mpc") * Mpc;
	cout << "Maximum Time: " << maxTime / Mpc << " Mpc" << endl;

	Emin = childValue(root, "MinEnergy_EeV") * EeV;
	cout << "Minimum Energy: " << Emin / EeV << " EeV" << endl;

	randomSeed = (int)childValue(root, "RandomSeed");
	cout << "RandomSeed: " << randomSeed << endl;

	xml_node node;
	std::string type;

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
		origin.x = childValue(node, "Xmin_Mpc", false) * Mpc;
		origin.y = childValue(node, "Ymin_Mpc", false) * Mpc;
		origin.z = childValue(node, "Zmin_Mpc", false) * Mpc;
		size.x = childValue(node, "Xmax_Mpc", false) * Mpc;
		size.y = childValue(node, "Ymax_Mpc", false) * Mpc;
		size.z = childValue(node, "Zmax_Mpc", false) * Mpc;
		size -= origin;
	} else
		throw runtime_error("Unknown environment");

	// magnetic field
	node = root.child("MagneticField");
	if (!node)
		throw runtime_error("Magnetic field not specified");
	type = node.attribute("type").as_string();
	cout << "MagenticField: " << type << endl;
	if (type == "None")
		magnetic_field = new UniformMagneticField(Vector3d(0, 0, 0));
	else if (type == "Uniform")
		loadUniformMagneticField(node);
	else if ((type == "LSS-Grid") or (type == "Kolmogoroff"))
		loadGridMagneticField(node);
	else {
		cout << " --> unknown, set zero field" << endl;
		magnetic_field = new UniformMagneticField(Vector3d(0, 0, 0));
	}

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

	// minimum energy
	modules.add(new MinimumEnergy(Emin));

	// maximum trajectory length
	modules.add(new MaximumTrajectoryLength(maxTime));

	// periodic boundaries
	loadPeriodicBoundaries();

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
	double epsilon = childValue(node, "MinStep_Mpc");
	cout << "  - Epsilon: " << epsilon << endl;
	double minStep = childValue(node, "MinStep_Mpc") * Mpc;
	cout << "  - Minimum Step: " << minStep / Mpc << " Mpc" << endl;
	cout << "  - Maximum Step: " << maxStep / Mpc << " Mpc" << endl;
	if (minStep >= maxStep)
		throw runtime_error("MaxStep must be larger than MinStep");
	modules.add(new DeflectionCK(magnetic_field, epsilon, minStep, maxStep));
}

void XmlExecute::loadUniformMagneticField(xml_node &node) {
	Vector3d bField;
	bField.x = childValue(node, "Bx_nG") * nG;
	bField.y = childValue(node, "By_nG") * nG;
	bField.z = childValue(node, "Bz_nG") * nG;
	cout << "  - Bx, By, Bz: " << bField / nG << " nG" << endl;
	magnetic_field = new UniformMagneticField(bField);
}

void XmlExecute::loadGridMagneticField(xml_node &node) {
	xml_node origin_node = node.child("Origin");
	origin.x = childValue(origin_node, "X_Mpc") * Mpc;
	origin.y = childValue(origin_node, "Y_Mpc") * Mpc;
	origin.z = childValue(origin_node, "Z_Mpc") * Mpc;
	cout << "  - Origin: " << origin / Mpc << endl;

	int Nx = childValue(node, "Nx");
	int Ny = childValue(node, "Ny");
	int Nz = childValue(node, "Nz");
	cout << "  - Samples: " << Nx << " * " << Ny << " * " << Nz << endl;

	double spacing = childValue(node, "Step_Mpc") * Mpc;
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
		double brms = childValue(node, "RMS_muG") * 1e-6 * gauss;
		cout << "  - Brms : " << brms / nG << " nG" << endl;

		double kMin = childValue(node, "Kmin");
		double kMax = childValue(node, "Kmax");
		double lMin = spacing / kMax;
		double lMax = spacing / kMin;
		cout << "  - Turbulent range: " << lMin / Mpc << " - " << lMax / Mpc << " Mpc" << endl;

		double alpha = childValue(node, "SpectralIndex");
		cout << "  - Turbulence spectral index: " << alpha << endl;
		cout << "  - Random seed: " << randomSeed << endl;

#ifdef MPC_HAVE_FFTW3F
		initTurbulence(field, brms, lMin, lMax, alpha, randomSeed);
#endif // MPC_HAVE_FFTW3F
#ifndef MPC_HAVE_FFTW3F
		throw runtime_error("Turbulent field grid not available. Compile with FFTW3F.");
#endif // MPC_HAVE_FFTW3F
	}

	magnetic_field = new MagneticFieldGrid(field);
}

void XmlExecute::loadPeriodicBoundaries() {
	if ((size.x == 0) or (size.y == 0) or (size.z == 0))
		throw runtime_error("Environment boundaries not set.");
	cout << "Periodic boundaries" << endl;
	cout << "  - Lower bounds: " << origin / Mpc << " Mpc" << endl;
	cout << "  - Upper bounds: " << (origin + size) / Mpc << " Mpc" << endl;
	modules.add(new PeriodicBox(origin, size));
}

void XmlExecute::loadSophia(xml_node &node) {
	maxStep = childValue(node, "MaxStep_Mpc") * Mpc;

	if (node.child("NoRedshift"))
		cout << "  - No redshift" << endl;
	else
		cout << "  - Redshift not implemented" << endl;

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

	if (node.child("NoPhotodisintegration"))
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
	double r = childValue(node, "Radius_Mpc") * Mpc;
	cout << "  - Radius: " << r / Mpc << " Mpc" << endl;

	int nObs = 0;
	for (xml_node n = node.child("SphereObserver"); n; n = n.next_sibling("SphereObserver")) {
		nObs += 1;
		Vector3d pos;
		pos.x = childValue(n, "CoordX_Mpc") * Mpc;
		pos.y = childValue(n, "CoordY_Mpc") * Mpc;
		pos.z = childValue(n, "CoordZ_Mpc") * Mpc;
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
		pos.x = childValue(n, "CoordX_Mpc") * Mpc;
		pos.y = childValue(n, "CoordY_Mpc") * Mpc;
		pos.z = childValue(n, "CoordZ_Mpc") * Mpc;
		cout << "  - Postion: " << pos / Mpc << " Mpc" << endl;
		double r = childValue(n, "Radius_Mpc") * Mpc;
		cout << "  - Radius: " << r / Mpc << " Mpc" << endl;
		SphericalBoundary *sphere = new SphericalBoundary(pos, r, "Detected");
		modules.add(sphere);
	}
	if (nObs > 1)
		cout << "  -> Warning for multiple observers: propagation will stop after first detection." << endl;
}

void XmlExecute::loadDiscreteSources(xml_node &node) {
	// source isotropic emission
	source.addProperty(new SourceIsotropicEmission());

	// source positions
	SourceMultiplePositions *positions = new SourceMultiplePositions();
	for (xml_node n = node.child("PointSource"); n; n = n.next_sibling("PointSource")) {
		Vector3d pos;
		pos.x = childValue(n, "CoordX_Mpc") * Mpc;
		pos.y = childValue(n, "CoordY_Mpc") * Mpc;
		pos.z = childValue(n, "CoordZ_Mpc") * Mpc;
		cout << "  - Position " << pos / Mpc << " Mpc" << endl;
		positions->add(pos);
	}
	source.addProperty(positions);

	// source spectrum
	xml_node spec = node.child("Spectrum");
	string spectrumType = spec.attribute("type").as_string();
	cout << "  - Spectrum: " << spectrumType << endl;

	if (spectrumType == "Monochromatic") {
		if (spec.child("Energy_EeV")) {
			double E = childValue(spec, "Energy_EeV") * EeV;
			cout << "  - Energy: " << E / EeV << " EeV" << endl;
			source.addProperty(new SourceEnergy(E));
		} else if (spec.child("Rigidity_EeV"))
			throw runtime_error("Fixed rigidity not implemented");
		else
			throw runtime_error("Source energy missing");

		// source composition
		SourceNuclei *composition = new SourceNuclei;
		xml_node p = node.child("Particles");
		for (xml_node n = p.child("Species"); n; n = n.next_sibling("Species")) {
			int A = (int)childValue(n, "MassNumber");
			int Z = (int)childValue(n, "ChargeNumber");
			double ab = childValue(n, "Abundance");
			composition->add(getNucleusId(A, Z), ab);
			cout << "  - Species: Z = " << Z << ", A = " << A << ", abundance = " << ab <<  endl;
		}
		source.addProperty(composition);

	} else if (spectrumType == "Power Law") {
		double alpha = childValue(spec, "Alpha");
		cout << "  - Power law index: " << alpha << endl;
		cout << "  - Minimum energy: " << Emin / EeV << " EeV" << endl;

		// if the source is accelerated to a maximum rigidity
		if (spec.child("Rigidity_EeV")) {
			double Rmax = childValue(spec, "Rigidity_EeV") * EeV;
			cout << "  - Maximum rigidity: " << Rmax / EeV << " EeV" << endl;

			// combined source spectrum / composition
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

		} else if (spec.child("Ecut_EeV")) {
			double Emax = childValue(spec, "Ecut_EeV") * EeV;
			cout << "  - Maximum energy: " << Emax / EeV << " EeV" << endl;

			// source spectrum
			source.addProperty(new SourcePowerLawSpectrum(Emin, Emax, alpha));

			// source composition
			SourceNuclei *composition = new SourceNuclei();
			xml_node p = node.child("Particles");
			for (xml_node n = p.child("Species"); n; n = n.next_sibling("Species")) {
				int A = n.attribute("MassNumber").as_int();
				int Z = n.attribute("ChargeNumber").as_int();
				double ab = n.attribute("Abundance").as_double();
				cout << "  - Species: Z = " << Z << ", A = " << A << ", abundance = " << ab <<  endl;
				composition->add(getNucleusId(A, Z), ab);
			}
			source.addProperty(composition);

		} else
			throw runtime_error(" --> maximum source energy / rigidity missing");

	} else
		throw runtime_error(" --> unknown source spectrum");
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
