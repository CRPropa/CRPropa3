#include "mpc/module/common.h"

#include "kiss/path.h"
#include "kiss/logger.h"

#include <stdlib.h>
#include <fstream>
#include <string>

namespace mpc {

std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (dataPath.size())
		return concat_path(dataPath, filename);

	const char *env_path = getenv("MPC_DATA_PATH");
	if (env_path) {
		if (is_directory(env_path)) {
			dataPath = env_path;
			KISS_LOG_INFO
				<< "getDataPath: use environment variable, " << dataPath
						<< std::endl;
			return concat_path(dataPath, filename);
		}
	}

#ifdef MPC_INSTALL_PREFIX
	{
		std::string _path = MPC_INSTALL_PREFIX "/share/mpc";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO
				<< "getDataPath: use define, " << dataPath << std::endl;
			return concat_path(dataPath, filename);
		}
	}
#endif

	{
		std::string _path = executable_path() + "../data";
		if (is_directory(_path)) {
			dataPath = _path;
			KISS_LOG_INFO
				<< "getDataPath: use executable path, " << dataPath
						<< std::endl;
			return concat_path(dataPath, filename);
		}
	}

	dataPath = "data";
	KISS_LOG_INFO
		<< "getDataPath: use default, " << dataPath << std::endl;
	return concat_path(dataPath, filename);
}

} // namespace mpc

