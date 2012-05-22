#include "mpc/Common.h"

#include "kiss/path.h"
#include "kiss/logger.h"

#include <stdlib.h>
#include <fstream>
#include <string>
#include <math.h>

namespace mpc {

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

double interpolate(const double x, const double *xD, const double *yD) {
	size_t i = 0;
	while (x > xD[i])
		i++;
	i--;
	return yD[i] + (x - xD[i]) * (yD[i + 1] - yD[i]) / (xD[i + 1] - xD[i]);
}

double interpolateEquidistant(const double x, const double xLo, const double dx,
		const double *yD) {
	double p = (x - xLo) / dx;
	size_t i = floor(p);
	return yD[i] + (p - i) * (yD[i + 1] - yD[i]);
}

} // namespace mpc

