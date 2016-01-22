#include "crpropa/Common.h"

#include "kiss/path.h"
#include "kiss/logger.h"

#include <stdlib.h>
#include <fstream>
#include <string>
#include <math.h>
#include <algorithm>

#define index(i,j) ((j)+(i)*Y.size())

namespace crpropa {

std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (dataPath.size())
		return concat_path(dataPath, filename);

	const char *env_path = getenv("CRPROPA_DATA_PATH");
	if (env_path) {
		if (is_directory(env_path)) {
			dataPath = env_path;
			KISS_LOG_INFO << "getDataPath: use environment variable, "
					<< dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}

#ifdef CRPROPA_INSTALL_PREFIX
	{
		std::string _path = CRPROPA_INSTALL_PREFIX "/share/crpropa";
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
  
  double interpolate2d(double x, double y, const std::vector<double> &X,
		       const std::vector<double> &Y, const std::vector<double> &Z) {
    
    std::vector<double>::const_iterator itx = std::upper_bound(X.begin(), X.end(), x); // Redshift
    std::vector<double>::const_iterator ity = std::upper_bound(Y.begin(), Y.end(), y); // Lorentz
    
    if (x > X.back() || x < X.front())
      return 0;
    if (y > Y.back() || y < Y.front())
      return 0;
    
    if (itx == X.begin() && ity == Y.begin())
      return Z.front();
    if (itx == X.end() && ity == Y.end())
      return Z.back();
    
    size_t i = itx - X.begin() - 1;
    size_t j = ity - Y.begin() - 1;
    
    double Q11 = Z[index(i,j)];
    double Q12 = Z[index(i,j+1)];
    double Q21 = Z[index(i+1,j)];
    double Q22 = Z[index(i+1,j+1)];
    
    double R1 = ((X[i+1]-x)/(X[i+1]-X[i]))*Q11+((x-X[i])/(X[i+1]-X[i]))*Q21;
    double R2 = ((X[i+1]-x)/(X[i+1]-X[i]))*Q12+((x-X[i])/(X[i+1]-X[i]))*Q22;

    return ((Y[j+1]-y)/(Y[j+1]-Y[j]))*R1+((y-Y[j])/(Y[j+1]-Y[j]))*R2;
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

} // namespace crpropa

