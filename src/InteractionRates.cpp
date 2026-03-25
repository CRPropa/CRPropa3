#include "crpropa/InteractionRates.h"
#include "crpropa/Units.h"
#include "crpropa/Random.h"

#include <nanoflann.hpp>

#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <filesystem>

#if defined(__APPLE__) && defined(_LIBCPP_VERSION)
	namespace fs = std::__fs::filesystem;
#else
	namespace fs = std::filesystem;
#endif

namespace crpropa {

InteractionRatesHomogeneous::InteractionRatesHomogeneous(std::string RateFile, std::string CumulativeRateFile) {
	
	this->ratesName = "interactionRatesHomogeneous";
	this->isPositionDependent = false;

	if (RateFile!="")
		initRate(RateFile);
	if (CumulativeRateFile!="")
		initCumulativeRate(CumulativeRateFile);
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedEnergy() const {
	return tabEnergy;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedRate() const {
	return tabRate;
}

std::vector<double> InteractionRatesHomogeneous::getTabulatedE() const {
	return tabE;
}

std::vector<double> InteractionRatesHomogeneous::getTabulateds() const {
	return tabs;
}

std::vector<std::vector<double>> InteractionRatesHomogeneous::getTabulatedCDF() const {
	return tabCDF;
}

double InteractionRatesHomogeneous::getProcessRate(const double E, const Vector3d &position) const {
	if (!this->isPositionDependent) {
		
		// compute the interaction rate for the given candidate energy, E
		double rate = interpolate(E, this->tabEnergy, this->tabRate);
		return rate;
		
	} else {
	
	  throw std::runtime_error("Error in boolean isPositionDependent!");
	  
	}
}

void InteractionRatesHomogeneous::loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const {
	tabE = this->getTabulatedE();
	tabs = this->getTabulateds();
	tabCDF = this->getTabulatedCDF();
}

void InteractionRatesHomogeneous::setTabulatedEnergy(std::vector<double>& tabEnergy) {
	this->tabEnergy = tabEnergy;
}

void InteractionRatesHomogeneous::setTabulatedRate(std::vector<double>& tabRate) {
	this->tabRate = tabRate;
}

void InteractionRatesHomogeneous::setTabulatedE(std::vector<double>& tabE) {
	this->tabE = tabE;
}

void InteractionRatesHomogeneous::setTabulateds(std::vector<double>& tabs) {
	this->tabs = tabs;
}

void InteractionRatesHomogeneous::setTabulatedCDF(std::vector<std::vector<double>>& tabCDF) {
	this->tabCDF = tabCDF;
}

void InteractionRatesHomogeneous::initRate(std::string filename){
	if(!fs::is_regular_file(filename))
		throw std::runtime_error("InteractionRatesHomogeneous: The given filename " + filename + " is a directory!\
			If you wanted to use position dependent photon fields instead, set the correct InteractionRates class first!");
	
	std::ifstream infile(filename.c_str());

	std::vector<double> tabEnergy;
	std::vector<double> tabRate;
		
	if (!infile.good())
		throw std::runtime_error("InteractionRatesHomogeneous: could not open file " + filename);

	while (infile.good()) {
		
		if (infile.peek() != '#') {	
			double a, b;
			infile >> a >> b;
			if (infile) {
				tabEnergy.push_back(pow(10, a) * eV);
				tabRate.push_back(b / Mpc);
			}
		}
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
		
	}
	
	infile.close();
		
	setTabulatedEnergy(tabEnergy);
	setTabulatedRate(tabRate);
}

void InteractionRatesHomogeneous::initCumulativeRate(std::string filename){
	if(!fs::is_regular_file(filename))
		throw std::runtime_error("InteractionRatesHomogeneous: The given filename " + filename + " is a directory!\
			If you wanted to use position dependent photon fields instead, set the correct InteractionRates class first!");

	std::ifstream infile(filename.c_str());
	
	std::vector<double> tabE;
	std::vector<double> tabs;
	std::vector<std::vector<double>> tabCDF;
	
	if (!infile.good())
		throw std::runtime_error("InteractionRatesHomogeneous: could not open file " + filename);
	
	// skip header
	while (infile.peek() == '#')
		infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
	
	// read s values in first line
	double a;
	infile >> a; // skip first value
	while (infile.good() and (infile.peek() != '\n')) {
		infile >> a;
		tabs.push_back(pow(10, a) * eV * eV);
	}
	
	// read all following lines: E, cdf values
	while (infile.good()) {
		infile >> a;
		if (!infile)
			break;  // end of file
		tabE.push_back(pow(10, a) * eV);
		std::vector<double> cdf;
		for (int i = 0; i < tabs.size(); i++) {
			infile >> a;
			cdf.push_back(a / Mpc);
		}
		tabCDF.push_back(cdf);
	}
	infile.close();
	
	setTabulatedE(tabE);
	setTabulateds(tabs);
	setTabulatedCDF(tabCDF);
}


InteractionRatesPositionDependent::InteractionRatesPositionDependent(
	std::string RateFilePath, std::string CumulativeRateFilePath, ref_ptr<Surface> surface) {
	
	this->ratesName = "interactionRatesPositionDependent";
	this->isPositionDependent = true;
	this->surface = surface;

	if (RateFilePath!="")
		initRate(RateFilePath);
	if (CumulativeRateFilePath!="")
		initCumulativeRate(CumulativeRateFilePath);
}

int InteractionRatesPositionDependent::findClosestGridPoint(const Vector3d &position) const {
  
  if (!tree) {
	throw std::runtime_error("KD-Tree not initialized!");
  }
  
  unsigned int closestIndex;
  double closestDistSquared;
  double queryPoint[3] = { position.x, position.y, position.z };
  
  this->tree->knnSearch(queryPoint, 1, &closestIndex, &closestDistSquared);
  return this->cloud.ids[closestIndex];
  
}

std::vector<double> InteractionRatesPositionDependent::getTabulatedEnergy() const {
	return tabEnergy;
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulatedRate() const {
	return tabRate;
}

std::vector<double> InteractionRatesPositionDependent::getTabulatedE() const {
	return tabE;
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getTabulateds() const {
	return tabs;
}

std::vector<std::vector<std::vector<double>>> InteractionRatesPositionDependent::getTabulatedCDF() const {
	return tabCDF;
}

std::unordered_map<int, Vector3d> InteractionRatesPositionDependent::getPhotonDict() const {
	return photonDict;
}

std::vector<double> InteractionRatesPositionDependent::getClosestRate(const Vector3d &position) const {
	int iMin = findClosestGridPoint(position);
	return tabRate[iMin];
}

std::vector<double> InteractionRatesPositionDependent::getClosests(const Vector3d &position) const {
	int iMin = findClosestGridPoint(position);
	return tabs[iMin];
}

std::vector<std::vector<double>> InteractionRatesPositionDependent::getClosestCDF(const Vector3d &position) const {
	int iMin = findClosestGridPoint(position);
	return tabCDF[iMin];
}

double InteractionRatesPositionDependent::getProcessRate(const double E, const Vector3d &position) const {
	if (!this->isPositionDependent) {
	
		throw std::runtime_error("Error in boolean isPositionDependent!");
		
	} else {
	  
	  std::vector<double> tabRate = this->getClosestRate(position);
	  
	  // compute the interaction rate for the given candidate energy, E
	  double rate = interpolate(E, this->tabEnergy, tabRate);
	  return rate;
	  
	}
}

void InteractionRatesPositionDependent::loadPerformInteractionTabs(const Vector3d &position, std::vector<double> &tabE, std::vector<double> &tabs, std::vector<std::vector<double>> &tabCDF) const {
	
	std::vector<double> E = this->getTabulatedE();
	std::vector<double> s = this->getClosests(position);
	std::vector<std::vector<double>> CDF = this->getClosestCDF(position);
	
	tabE = E;
	tabs = s;
	tabCDF = CDF;
}

void InteractionRatesPositionDependent::setTabulatedEnergy(std::vector<double>& tabEnergy) {
	this->tabEnergy = tabEnergy;
}

void InteractionRatesPositionDependent::setTabulatedRate(std::vector<std::vector<double>>& tabRate) {
	this->tabRate = tabRate;
}

void InteractionRatesPositionDependent::setTabulatedE(std::vector<double>& tabE) {
	this->tabE = tabE;
}

void InteractionRatesPositionDependent::setTabulateds(std::vector<std::vector<double>>& tabs) {
	this->tabs = tabs;
}

void InteractionRatesPositionDependent::setTabulatedCDF(std::vector<std::vector<std::vector<double>>>& tabCDF) {
	this->tabCDF = tabCDF;
}

void InteractionRatesPositionDependent::setPhotonDict(std::unordered_map<int, Vector3d>& photonDict) {
	this->photonDict = photonDict;
	
	// delete old clouds
	this->cloud.points.clear();
	this->cloud.ids.clear();
	
	for (const auto& el : this->photonDict) {
		this->cloud.ids.push_back(el.first);
		this->cloud.points.push_back(el.second);
	}
	
	// delete old tree
	if (this->tree) {
		delete this->tree;
	}
	
	int maxLeafTree = 20;
	int nThreads = 4;
	nanoflann::KDTreeSingleIndexAdaptorFlags flag;

	this->tree = new KDTree(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(maxLeafTree, flag, nThreads));
	this->tree->buildIndex();

}

void InteractionRatesPositionDependent::setSurface(ref_ptr<Surface> surface) {
		this->surface = surface;
}

ref_ptr<Surface> InteractionRatesPositionDependent::getSurface() const {
		return this->surface;
}

void InteractionRatesPositionDependent::initRate(std::string filepath){
	if(!fs::is_directory(filepath))
		throw std::runtime_error("InteractionRatesPositionDependent: The given path " + filepath + " is not a directory!\
			If you wanted to use homogeneous photon fields instead, set the correct InteractionRates class first!");

	std::vector<std::vector<double>> tabRate;
		
	fs::path dir = filepath;
	std::unordered_map<int, Vector3d> photonDict;
	int iFile = 0;
	
	if (!fs::exists(dir)) {
			std::cout << "Photon tables not found in " << dir << std::endl;
			return;
	}
	
	for (auto const& dir_entry : fs::directory_iterator{dir}) {
		
		// the input filename here should be a string
		//check if it is correct, i.e. a proper filename string
		std::string filename = dir_entry.path().string();
		std::ifstream infile(filename.c_str());
		
		std::vector<double> vecEnergy;
		std::vector<double> vecRate;
		
		if (!infile.good())
			throw
			std::runtime_error("InteractionRatesPositionDependent: could not open file " + filename);
		
		double x, y, z;
		std::string str;
		std::stringstream ss;
		
		std::string filename_split = splitFilename(dir_entry.path().string());
		ss << filename_split;
		
		int iLine = 0;
		
		std::locale::global(std::locale("C"));
		
		while (getline(ss, str, '_')) {
			if (iLine == 3) {
				x = -std::stod(str) * kpc;
			}
			if (iLine == 4) {
				y = std::stod(str) * kpc;
			}
			if (iLine == 5) {
				z = std::stod(str) * kpc;
			}
			iLine = iLine + 1;
		}
		
		Vector3d vPos(x, y, z);
		
		if (getSurface() and !getSurface()->isInside(vPos))
			continue;
		
		photonDict[iFile] = vPos;
		
		while (infile.good()) {
			if (infile.peek() != '#') {
				double a, b;
				infile >> a >> b;
				if (infile) {
					if (iFile == 0) {
						vecEnergy.push_back(pow(10, a) * eV);
						setTabulatedEnergy(vecEnergy);
					}
					vecRate.push_back(b / Mpc);
				}
			}
			infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
		}
		
		tabRate.push_back(vecRate);
		
		iFile = iFile + 1;
		infile.close();
		
	}
		
	if (tabRate.empty())
		throw std::runtime_error("Rate's table empty! Check if the surface is properly set.");
		
	setTabulatedRate(tabRate);
	setPhotonDict(photonDict);
}

void InteractionRatesPositionDependent::initCumulativeRate(std::string filepath){
	if(!fs::is_directory(filepath))
		throw std::runtime_error("InteractionRatesPositionDependent: The given path " + filepath + " is not a directory!\
			If you wanted to use homogeneous photon fields instead, set the correct InteractionRates class first!");

	std::vector<std::vector<double>> tabs;
	std::vector<std::vector<std::vector<double>>> tabCDF;
	
	fs::path dir = filepath;
	int iFile = 0;
	
	if (!fs::exists(dir)) {
			std::cout << "Photon tables not found in " << dir << std::endl;
			return;
	}
 
	for (auto const& dir_entry : fs::directory_iterator{dir}) {
		
		std::vector<double> vecE;
		std::vector<double> vecs;
		std::vector<std::vector<double>> vecCDF;
		
		std::string filename = dir_entry.path().string();
		std::ifstream infile(filename.c_str());
		
		if (!infile.good())
			throw std::runtime_error("InteractionRatesPositionDependent: could not open file " + filename);
		
		double x, y, z;
		std::string str;
		std::stringstream ss;
		
		std::string filename_split = splitFilename(dir_entry.path().string());
		ss << filename_split;
		
		int iLine = 0;
		
		std::locale::global(std::locale("C"));
		
		while (getline(ss, str, '_')) {
			if (iLine == 3) {
				x = -std::stod(str) * kpc;
			}
			if (iLine == 4) {
				y = std::stod(str) * kpc;
			}
			if (iLine == 5) {
				z = std::stod(str) * kpc;
			}
			iLine = iLine + 1;
		}
		
		Vector3d vPos(x, y, z);
		
		if (getSurface() and !getSurface()->isInside(vPos))
			continue;
		
		// skip header
		while (infile.peek() == '#')
			infile.ignore(std::numeric_limits < std::streamsize > ::max(), '\n');
		
		// read s values in first line
		double a;
		infile >> a; // skip first value
		while (infile.good() and (infile.peek() != '\n')) {
			infile >> a;
			vecs.push_back(pow(10, a) * eV * eV);
		}
		
		// read all following lines: E, cdf values
		while (infile.good()) {
			infile >> a;
			if (!infile)
				break;  // end of file
			if (iFile == 0) {
				vecE.push_back(pow(10, a) * eV);
				setTabulatedE(vecE);
			}
			std::vector<double> cdf;
			for (int i = 0; i < tabs.size(); i++) {
				infile >> a;
				cdf.push_back(a / Mpc);
			}
			vecCDF.push_back(cdf);
		}
		
		iFile = iFile + 1;
		
		tabs.push_back(vecs);
		tabCDF.push_back(vecCDF);
		infile.close();
	}
	
	setTabulateds(tabs);
	setTabulatedCDF(tabCDF);
}

} //namespace crpropa
