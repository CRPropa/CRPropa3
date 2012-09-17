#include "mpc/Common.h"

#include "kiss/path.h"
#include "kiss/logger.h"

#include <stdlib.h>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>

namespace mpc {

// Overall redshift scaling of the Kneiske IRB (see data/PhotonField/KneiskeIRB.py and Kneiske et al. 2004, astro-ph/0309141)
double a[9] = { 0, 0.2, 0.4, 0.6, 1, 2, 3, 4, 5 };
static std::vector<double> zKneiske(a, a + sizeof(a) / sizeof(double));
double b[9] = { 1, 1.6937, 2.5885, 3.6178, 5.1980, 7.3871, 8.5471, 7.8605, 0 };
static std::vector<double> sKneiske(b, b + sizeof(b) / sizeof(double));

double photonFieldScaling(int photonField, double z) {
	if (photonField == IRB)
		return interpolate(z, zKneiske, sKneiske);
	return pow(z + 1, 3); // CMB-like scaling
}

std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (dataPath.size())
		return concat_path(dataPath, filename);

	const char *env_path = getenv("MPC_DATA_PATH");
	if (env_path) {
		if (is_directory(env_path)) {
			dataPath = env_path;
			KISS_LOG_INFO << "getDataPath: use environment variable, "
					<< dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}

#ifdef MPC_INSTALL_PREFIX
	{
		std::string _path = MPC_INSTALL_PREFIX "/share/mpc";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO
			<< "getDataPath: use install prefix, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}
#endif

	{
		std::string _path = executable_path() + "../data";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO << "getDataPath: use executable path, " << dataPath
					<< std::endl;
			return concat_path(dataPath, filename);
		}
	}

	dataPath = "data";
	KISS_LOG_INFO << "getDataPath: use default, " << dataPath << std::endl;
	return concat_path(dataPath, filename);
}

double interpolate(double x, const std::vector<double> &X,
		const std::vector<double> &Y) {
	std::vector<double>::const_iterator it = std::upper_bound(X.begin(),
			X.end(), x);
	if (it == X.begin())
		return Y.front();
	if (it == X.end())
		return Y.back();

	size_t i = it - X.begin() - 1;
	return Y[i] + (x - X[i]) * (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]);
}

double interpolateEquidistant(double x, double lo, double hi,
		const std::vector<double> &Y) {
	if (x <= lo)
		return Y.front();
	if (x >= hi)
		return Y.back();

	double dx = (hi - lo) / (Y.size() - 1);
	double p = (x - lo) / dx;
	size_t i = floor(p);
	return Y[i] + (p - i) * (Y[i + 1] - Y[i]);
}

} // namespace mpc

