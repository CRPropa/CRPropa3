#include "crpropa/XmlExecute.h"
#include "crpropa/magneticField/MagneticFieldGrid.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"
#include "crpropa/PhotonBackground.h"
#include "crpropa/Cosmology.h"
#include "crpropa/ParticleID.h"
#include "crpropa/module/SimplePropagation.h"
#include "crpropa/module/PropagationCK.h"
#include "crpropa/module/Redshift.h"
#include "crpropa/module/ElectronPairProduction.h"
#include "crpropa/module/PhotoPionProduction.h"
#include "crpropa/module/PhotoDisintegration.h"
#include "crpropa/module/NuclearDecay.h"
#include "crpropa/module/BreakCondition.h"
#include "crpropa/module/Boundary.h"
#include "crpropa/module/TextOutput.h"
#include "crpropa/module/OutputROOT.h"
#include "crpropa/module/OutputCRPropa2.h"
#include "crpropa/ModuleList.h"

#include "pugixml.hpp"
#include "kiss/string.h"

#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>

using namespace pugi;
using namespace std;

namespace crpropa {

double childValue(xml_node parent, string childName, bool throwIfEmpty = true) {
	xml_node node = parent.child(childName.c_str());
	if (!node and throwIfEmpty) {
		stringstream ss;
		ss << "Error reading XML card: " << childName << " not specified";
		throw runtime_error(ss.str());
	}
	return node.attribute("value").as_double();
}

xml_node childNode(xml_node parent, string childName,
		bool throwIfEmpty = true) {
	xml_node node = parent.child(childName.c_str());
	if (!node and throwIfEmpty) {
		stringstream ss;
		ss << "Error reading XML card: " << childName << " not specified";
		throw runtime_error(ss.str());
	}
	return node;
}

ref_ptr<SourceUniformBox> loadSourceHomogeneousBox(pugi::xml_node &node) {
	Vector3d origin;
	origin.x = childValue(node, "Xmin_Mpc") * Mpc;
	origin.y = childValue(node, "Ymin_Mpc") * Mpc;
	origin.z = childValue(node, "Zmin_Mpc") * Mpc;
	cout << "  - Origin: " << origin / Mpc << endl;

	Vector3d size;
	size.x = childValue(node, "Xmax_Mpc") * Mpc;
	size.y = childValue(node, "Ymax_Mpc") * Mpc;
	size.z = childValue(node, "Zmax_Mpc") * Mpc;
	size -= origin;
	cout << "  - Size: " << size / Mpc << endl;

	return (new SourceUniformBox(origin, size));
}

ref_ptr<SourceDensityGrid> loadSourceDensityGrid(pugi::xml_node &node) {
	int nx = childValue(node, "Nx");
	int ny = childValue(node, "Ny");
	int nz = childValue(node, "Nz");
	cout << "  - Nx = " << nx << ", Ny = " << ny << ", Nz = " << nz << endl;

	double spacing = childValue(node, "Step_Mpc") * Mpc;
	cout << "  - Spacing = " << spacing / Mpc << " Mpc" << endl;

	xml_node origin_node = childNode(node, "Origin");
	Vector3d origin;
	origin.x = childValue(origin_node, "X_Mpc") * Mpc;
	origin.y = childValue(origin_node, "Y_Mpc") * Mpc;
	origin.z = childValue(origin_node, "Z_Mpc") * Mpc;
	cout << "  - Origin = " << origin / Mpc << " Mpc" << endl;

	ref_ptr<ScalarGrid> grid = new ScalarGrid(origin, nx, ny, nz, spacing);

	xml_node file_node = childNode(node, "File");
	string file_type = file_node.attribute("type").as_string();
	string file_name = file_node.child_value();
	cout << "  - File = " << file_name << endl;
	if (file_type == "ASCII")
		loadGridFromTxt(grid, file_name);
	else if (file_type == "FITS")
		throw runtime_error(" --> FITS files not supported");
	else
		throw runtime_error(" --> unknown file type");

	return (new SourceDensityGrid(grid));
}

ref_ptr<SourceDensityGrid1D> loadSourceDensityGrid1D(pugi::xml_node &node) {
	int nx = childValue(node, "Nx");
	cout << "  - Nx = " << nx << endl;

	double spacing = childValue(node, "Step_Mpc") * Mpc;
	cout << "  - Spacing = " << spacing / Mpc << " Mpc" << endl;

	// convert to comoving spacing
	spacing = lightTravel2ComovingDistance(spacing);

	ref_ptr<ScalarGrid> grid = new ScalarGrid(Vector3d(0, 0, 0), nx, 1, 1,
			spacing);

	xml_node file_node = childNode(node, "File");
	string file_type = file_node.attribute("type").as_string();
	string file_name = file_node.child_value();
	cout << "  - File = " << file_name << endl;
	if (file_type == "ASCII")
		loadGridFromTxt(grid, file_name);
	else if (file_type == "FITS")
		throw runtime_error(" --> FITS files not supported");
	else
		throw runtime_error(" --> unknown file type");

	return (new SourceDensityGrid1D(grid));
}

XmlExecute::XmlExecute() :
		is1D(false), hasRedshift(false), nTrajectories(0), Emin(0), maxStep(0) {
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
		throw runtime_error(
				"Error reading XML card: Root element CRPropa not found");

	// ----- general settings -----
	nTrajectories = (int) childValue(root, "TrajNumber");
	cout << "Number of particles: " << nTrajectories << endl;

	xml_node seed_node = root.child("RandomSeed");
	if (seed_node) {
		int seed = seed_node.attribute("value").as_int();
		Random::seedThreads(seed);
		cout << "Random seed: " << seed << endl;
	} else
		cout << "No random seed given. Using random random seed." << endl;

	// ----- environment -----
	xml_node node;
	node = childNode(root, "Environment");
	string type = node.attribute("type").as_string();
	cout << "Environment: " << type << endl;

	if (type == "One Dimension")
		is1D = true;
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
		throw runtime_error(" --> unknown environment");

	// ----- magnetic field -----
	node = childNode(root, "MagneticField", false);
	if (node) {
		type = node.attribute("type").as_string();
		cout << "MagenticField: " << type << endl;
		if ((type == "Null") or (type == "None"))
			magnetic_field = new UniformMagneticField(Vector3d(0, 0, 0));
		else if (type == "Uniform")
			loadUniformMagneticField(node);
		else if ((type == "LSS-Grid") or (type == "Kolmogoroff"))
			loadGridMagneticField(node);
		else if (type == "1D")
			cout << " --> not implemented" << endl;
		else {
			cout << " --> unknown, set zero field" << endl;
			magnetic_field = new UniformMagneticField(Vector3d(0, 0, 0));
		}
	} else {
		if (!is1D)
			throw runtime_error(" --> magnetic field not specified");
	}

	// ----- propagator -----
	xml_node interaction_node = childNode(root, "Interactions");
	maxStep = childValue(interaction_node, "MaxStep_Mpc", false) * Mpc;

	if (is1D) {
		cout << "Propagator: 1D" << endl;
		modules.add(new SimplePropagation());

		bool noRedshift = interaction_node.child("NoRedshift");
		hasRedshift = !(noRedshift);

		if (hasRedshift) {
			modules.add(new Redshift());

			double omegaM = 0.3;
			if (root.child("OmegaM"))
				omegaM = childValue(root, "OmegaM");

			double omegaL = 0.7;
			if (root.child("OmegaLambda"))
				omegaL = childValue(root, "OmegaLambda");

			double H0 = 70.;
			if (root.child("H0_km_s_Mpc"))
				H0 = childValue(root, "H0_km_s_Mpc");

			cout << "Cosmology: OmegaM = " << omegaM << ", OmegaLambda = "
					<< 1 - omegaM << ", H0 = " << H0 << " km/s/Mpc" << endl;
			setCosmologyParameters(H0 / 100, omegaM);
		} else {
			cout << "  - No redshift" << endl;
		}
	} else {
		node = childNode(root, "Integrator");
		type = node.attribute("type").as_string();
		cout << "Propagator: " << type << endl;
		if (type == "Cash-Karp RK")
			loadDeflectionCK(node);
		else
			throw runtime_error(" --> unknown integrator");
	}

	// ----- interactions -----
	type = interaction_node.attribute("type").as_string();
	cout << "Interactions: " << type << endl;
	if (type == "Sophia")
		loadSophia(interaction_node);

	// ----- minimum energy -----
	Emin = childValue(root, "MinEnergy_EeV") * EeV;
	cout << "Minimum energy: " << Emin / EeV << " EeV" << endl;
	modules.add(new MinimumEnergy(Emin));

	// ----- maximum trajectory length -----
	double maxTime = childValue(root, "MaxTime_Mpc") * Mpc;
	cout << "Maximum time: " << maxTime / Mpc << " Mpc" << endl;
	if (is1D)
		maxTime = lightTravel2ComovingDistance(maxTime); // convert to comoving distance
	modules.add(new MaximumTrajectoryLength(maxTime));

	// ----- periodic boundaries -----
	if (!is1D)
		loadPeriodicBoundaries();

	// ----- sources -----
	node = childNode(root, "Sources");

	// emission direction
	if (!is1D)
		source.add(new SourceIsotropicEmission());

	// position
	type = node.attribute("type").as_string();
	cout << "Source(s): " << type << endl;
	if (type == "Discrete")
		loadDiscreteSources(node);
	else if (type == "Continuous")
		loadContinuousSources(node);
	else
		throw runtime_error(" --> unknown source type");

	if (is1D and hasRedshift) {
		cout << "  - Redshift according to source distance" << endl;
		source.add(new SourceRedshift1D());
	}

	// spectrum + composition
	loadSpectrumComposition(node);

	// ----- observers -----
	ref_ptr<Observer> observer = new Observer();
	modules.add(observer);
	if (is1D) {
		observer->add(new ObserverPoint());
	} else {
		node = root.child("Observers");
		if (!node) {
			cout << "Observer(s) not specified" << endl;
		} else {
			string type = node.attribute("type").as_string();
			cout << "Observers: " << type << endl;
			if (type == "Spheres around Observers")
				loadSpheresAroundObserver(node);
			else if (type == "Spheres around Source")
				loadSpheresAroundSource(node);
			else
				cout << " --> unknown observer ('Spheres around Source' or "
						<< "'Spheres around Observer')" << endl;
		}
	}

	// ----- output -----
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
	cout << "  - Minimum step: " << minStep / Mpc << " Mpc" << endl;

	if (maxStep == 0)
		maxStep = numeric_limits<double>::max();
	cout << "  - Maximum step: " << maxStep / Mpc << " Mpc" << endl;

	if (minStep >= maxStep)
		throw runtime_error(
				" --> Maximum step must be larger than minimum step");

	modules.add(new PropagationCK(magnetic_field, epsilon, minStep, maxStep));
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

	string type = node.attribute("type").as_string();
	if (type == "LSS-Grid") {
		cout << "  - Loading field from file" << endl;

		string filetype = node.child("File").attribute("type").as_string();
		cout << "  - File type: " << filetype << endl;

		string filename = kiss::trim(node.child_value("File"));
		cout << "  - File name: " << filename << endl;

		if (filetype == "ASCII")
			loadGridFromTxt(field, filename, gauss);
		else if (filetype == "RAW")
			loadGrid(field, filename, gauss);
		else
			throw runtime_error("Unsupported file type");

	} else if (type == "Kolmogoroff") {
		cout << "  - Creating random turbulent field" << endl;

		double brms = childValue(node, "RMS_muG") * 1e-6 * gauss;
		cout << "  - Brms: " << brms / nG << " nG" << endl;

		double kMin = childValue(node, "Kmin");
		double kMax = childValue(node, "Kmax");
		double lMin = spacing / kMax;
		double lMax = spacing / kMin;
		cout << "  - Turbulent range: " << lMin / Mpc << " - " << lMax / Mpc
				<< " Mpc" << endl;

		double alpha = childValue(node, "SpectralIndex");
		cout << "  - Spectral index, <B^2(k)> ~ k^n, n:  " << alpha << endl;

#ifdef CRPROPA_HAVE_FFTW3F
		initTurbulence(field, brms, lMin, lMax, alpha);
#endif
#ifndef CRPROPA_HAVE_FFTW3F
		throw runtime_error(
				"Turbulent field not available: Compile with FFTW3F.");
#endif
	} else {
		throw runtime_error("Unknown field type");
	}

	magnetic_field = new MagneticFieldGrid(field);
}

void XmlExecute::loadPeriodicBoundaries() {
	if ((size.x == 0) or (size.y == 0) or (size.z == 0))
		throw runtime_error(" --> Environment boundaries not set");
	cout << "Periodic boundaries" << endl;
	cout << "  - Lower bounds: " << origin / Mpc << " Mpc" << endl;
	cout << "  - Upper bounds: " << (origin + size) / Mpc << " Mpc" << endl;
	modules.add(new PeriodicBox(origin, size));
}

void XmlExecute::loadSophia(xml_node &node) {
	bool pairprodCMB = true;
	bool pairprodIRB = true;
	bool pionprodCMB = true;
	bool pionprodIRB = true;
	bool photodisCMB = true;
	bool photodisIRB = true;
	bool decay = true;

	if (node.child("NoPairProd")) {
		cout << "  - No pair production" << endl;
		pairprodCMB = false;
		pairprodIRB = false;
	}
	if (node.child("NoPionProd")) {
		cout << "  - No pion production" << endl;
		pionprodCMB = false;
		pionprodIRB = false;
	}
	if (node.child("NoPhotodisintegration")) {
		cout << "  - No photo disintegration" << endl;
		photodisCMB = false;
		photodisIRB = false;
	}
	if (node.child("NoIRO")) {
		cout << "  - No interactions on IRB" << endl;
		pairprodIRB = false;
		pionprodIRB = false;
		photodisIRB = false;
	}
	if (node.child("NoIRPionProd")) {
		cout << "  - No pion production on IRB" << endl;
		pionprodIRB = false;
	}
	if (node.child("NoDecay")) {
		cout << "  - No decay" << endl;
		decay = false;
	}

	if (pairprodCMB)
		modules.add(new ElectronPairProduction(CMB));
	if (pairprodIRB)
		modules.add(new ElectronPairProduction(IRB));

	if (pionprodCMB)
		modules.add(new PhotoPionProduction(CMB));
	if (pionprodIRB)
		modules.add(new PhotoPionProduction(IRB));

	if (photodisCMB)
		modules.add(new PhotoDisintegration(CMB));
	if (photodisIRB)
		modules.add(new PhotoDisintegration(IRB));

	if (decay)
		modules.add(new NuclearDecay);
}

void XmlExecute::loadSpheresAroundObserver(xml_node &node) {
	double r = childValue(node, "Radius_Mpc") * Mpc;
	cout << "  - Radius: " << r / Mpc << " Mpc" << endl;

	for (xml_node n = node.child("SphereObserver"); n;
			n = n.next_sibling("SphereObserver")) {
		Vector3d pos;
		pos.x = childValue(n, "CoordX_Mpc") * Mpc;
		pos.y = childValue(n, "CoordY_Mpc") * Mpc;
		pos.z = childValue(n, "CoordZ_Mpc") * Mpc;
		cout << "  - Postion: " << pos / Mpc << " Mpc";
		cout << ", Detection does not stop propagation" << endl;
		observer.add(new ObserverSmallSphere(pos, r));
	}
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
		cout << "  - Detection does not stop propagation" << endl;
		observer.add(new ObserverLargeSphere(pos, r));
	}
}

void XmlExecute::loadDiscreteSources(pugi::xml_node &node) {
	ref_ptr<SourceMultiplePositions> sourcePositions =
			new SourceMultiplePositions();

	xml_node density_node = node.child("Density");
	if (density_node) {
		// draw positions from density distribution
		string type = density_node.attribute("type").as_string();
		cout << "  - Density: " << type << endl;
		int nSources = childValue(node, "Number");
		cout << "  - Number of sources = " << nSources << endl;

		ref_ptr<SourceFeature> sourceDistribution = new SourceFeature();
		if (type == "Uniform") {
			if (is1D) {
				double xmin = childValue(density_node, "Xmin_Mpc") * Mpc;
				double xmax = childValue(density_node, "Xmax_Mpc") * Mpc;
				xmin = lightTravel2ComovingDistance(xmin);
				xmax = lightTravel2ComovingDistance(xmax);
				sourceDistribution = new SourceUniform1D(xmin, xmax);
			} else {
				sourceDistribution = loadSourceHomogeneousBox(density_node);
			}
		} else if (type == "Grid") {
			if (is1D)
				sourceDistribution = loadSourceDensityGrid1D(density_node);
			else
				sourceDistribution = loadSourceDensityGrid(density_node);
		} else {
			throw runtime_error(" --> unknown source density type");
		}
		ParticleState p;
		for (int i = 0; i < nSources; i++) {
			sourceDistribution->prepareParticle(p);
			sourcePositions->add(p.getPosition());
			cout << "  - Position = " << p.getPosition() / Mpc << " Mpc"
					<< endl;
		}
	} else {
		// loop over point sources
		for (xml_node n = node.child("PointSource"); n;
				n = n.next_sibling("PointSource")) {
			Vector3d pos(0.);
			if (is1D) {
				// 1D
				double dlt = childValue(n, "CoordX_Mpc") * Mpc;
				pos.x = lightTravel2ComovingDistance(dlt);
				cout << "  - Light travel distance = " << dlt / Mpc << " Mpc"
						<< endl;
			} else {
				// 3D
				pos.x = childValue(n, "CoordX_Mpc") * Mpc;
				pos.y = childValue(n, "CoordY_Mpc") * Mpc;
				pos.z = childValue(n, "CoordZ_Mpc") * Mpc;
				cout << "  - Position = " << pos / Mpc << " Mpc" << endl;
			}
			sourcePositions->add(pos);
		}
	}
	source.add(sourcePositions);
}

void XmlExecute::loadContinuousSources(pugi::xml_node &node) {
	xml_node density_node = node.child("Density");
	string type = density_node.attribute("type").as_string();
	cout << "  - Density: " << type << endl;
	if (type == "Uniform")
		if (is1D) {
			double minD = childValue(density_node, "Xmin_Mpc") * Mpc;
			double maxD = childValue(density_node, "Xmax_Mpc") * Mpc;
			cout << "  - Minimum light travel distance: " << minD / Mpc
					<< " Mpc" << endl;
			cout << "  - Maximum light travel distance: " << maxD / Mpc
					<< " Mpc" << endl;
			minD = lightTravel2ComovingDistance(minD);
			maxD = lightTravel2ComovingDistance(maxD);
			source.add(new SourceUniform1D(minD, maxD));
		} else {
			source.add(loadSourceHomogeneousBox(density_node));
		}
	else if (type == "Grid") {
		if (is1D) {
			source.add(loadSourceDensityGrid1D(density_node));
		} else
			source.add(loadSourceDensityGrid(density_node));
	} else {
		throw runtime_error(" --> unknown source density type");
	}
}

void XmlExecute::loadSpectrumComposition(pugi::xml_node &node) {
	xml_node spectrum_node = node.child("Spectrum");
	string spectrumType = spectrum_node.attribute("type").as_string();
	cout << "  - Spectrum: " << spectrumType << endl;

	if (spectrumType == "Monochromatic") {
		if (spectrum_node.child("Energy_EeV")) {
			double E = childValue(spectrum_node, "Energy_EeV") * EeV;
			cout << "  - Energy: " << E / EeV << " EeV" << endl;
			source.add(new SourceEnergy(E));
		} else if (spectrum_node.child("Rigidity_EeV"))
			throw runtime_error(" --> Fixed rigidity not implemented");
		else {
			throw runtime_error(" --> Source energy missing");
		}
		loadSourceNuclei(node);
	} else if (spectrumType == "Power Law") {
		double alpha = -childValue(spectrum_node, "Alpha");
		cout << "  - Power law index E^a, a:" << alpha << endl;
		cout << "  - Minimum energy: " << Emin / EeV << " EeV" << endl;

		// if the source is accelerated to a maximum rigidity
		if (spectrum_node.child("Rigidity_EeV")) {
			double Rmax = childValue(spectrum_node, "Rigidity_EeV") * EeV;
			cout << "  - Maximum rigidity: " << Rmax / EeV << " EeV" << endl;

			// combined source spectrum / composition
			ref_ptr<SourceComposition> comp = new SourceComposition(Emin, Rmax,
					alpha);
			xml_node p = node.child("Particles");
			for (xml_node n = p.child("Species"); n;
					n = n.next_sibling("Species")) {
				int A = n.attribute("MassNumber").as_int();
				int Z = n.attribute("ChargeNumber").as_int();
				double ab = n.attribute("Abundance").as_double();
				cout << "  - Species: Z = " << Z << ", A = " << A
						<< ", abundance = " << ab << endl;
				comp->add(nucleusId(A, Z), ab);
			}
			source.add(comp);
		} else if (spectrum_node.child("Ecut_EeV")) {
			double Emax = childValue(spectrum_node, "Ecut_EeV") * EeV;
			cout << "  - Maximum energy: " << Emax / EeV << " EeV" << endl;

			// source spectrum
			source.add(new SourcePowerLawSpectrum(Emin, Emax, alpha));

			// source composition
			loadSourceNuclei(node);
		} else {
			throw runtime_error(
					" --> maximum source energy / rigidity missing");
		}
	} else {
		throw runtime_error(" --> unknown source spectrum");
	}
}

void XmlExecute::loadSourceNuclei(pugi::xml_node &node) {
	ref_ptr<SourceMultipleParticleTypes> composition =
			new SourceMultipleParticleTypes();
	xml_node p = node.child("Particles");
	for (xml_node n = p.child("Species"); n; n = n.next_sibling("Species")) {
		int A = n.attribute("MassNumber").as_int();
		int Z = n.attribute("ChargeNumber").as_int();
		double ab = n.attribute("Abundance").as_double();
		cout << "  - Species: Z = " << Z << ", A = " << A << ", abundance = "
				<< ab << endl;
		composition->add(nucleusId(A, Z), ab);
	}
	source.add(composition);
}

void XmlExecute::loadOutput(xml_node &node) {
	string type = node.attribute("type").as_string();
	cout << "Output: " << type << endl;

	xml_node file_node = node.child("File");
	string format = file_node.attribute("type").as_string();
	cout << "  - Filetype: " << format << endl;

	string filename = kiss::trim(node.child("File").child_value());
	cout << "  - Filename: " << filename << endl;

	string option = node.child("File").attribute("option").as_string();
	if (option != "force") {
		ifstream ifile(filename.c_str());
		if (ifile)
			throw runtime_error("Output file already exists!");
	}

	if (format == "ASCII") {
		if (type == "Full Trajectories") {
			if (is1D)
				modules.add(new CRPropa2TrajectoryOutput1D(filename));
			else
				modules.add(new CRPropa2TrajectoryOutput3D(filename));
		} else if (type == "Events") {
			if (is1D)
				observer.onDetection(new CRPropa2EventOutput1D(filename));
			else
				observer.onDetection(new CRPropa2EventOutput3D(filename));
		} else if (type == "None") {
			return;
		} else {
			cout << "  --> unknown output type "
					<< "('Events', 'Full Trajectories' or 'None')" << endl;
		}
	}
#ifdef CRPROPA_HAVE_ROOT
	else if (format == "ROOT") {
		if (type == "Full Trajectories")
			if (is1D)
				modules.add(new CRPropa2ROOTTrajectoryOutput1D(filename));
			else
				modules.add(new CRPropa2ROOTTrajectoryOutput3D(filename));
		else if (type == "Events")
			if (is1D)
				observer.onDetection(new CRPropa2ROOTEventOutput1D(filename));
			else
				observer.onDetection(new CRPropa2ROOTEventOutput3D(filename));
		else if (type == "None")
			return;
		else
			cout << "  --> unknown output type "
					<< "('Events', 'Full Trajectories' or 'None')" << endl;
	}
#endif // CRPROPA_HAVE_ROOT
	else {
		cout << "  --> unknown output format. "
				<< " Use 'ASCII' or 'ROOT' (if ROOT is set)" << endl;
	}
}

void XmlExecute::run() {
	cout << endl << "Active modules:" << endl;
	cout << modules.getDescription();
	modules.setShowProgress(true);
	modules.run(&source, nTrajectories, true);
}

} // namespace crpropa
