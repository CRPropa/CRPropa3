#include "mpc/module/common.h"

#include "kiss/path.h"

#include <stdlib.h>
#include <fstream>
#include <string>

namespace mpc {

std::string getDataPath(std::string filename) {
	static std::string dataPath;
	if (dataPath.size())
		return dataPath;

	const char *env_path = getenv("MPC_DATA_PATH");
	if (env_path) {
		if (is_directory(env_path)) {
			dataPath = env_path;
			return dataPath;
		}
	}

#ifdef MPC_INSTALL_PREFIX
	{
		std::string _path = MPC_INSTALL_PREFIX "/share/mpc";
		if (is_directory(_path)) {
			dataPath = _path;
			return dataPath;
		}
	}
#endif

	{
		std::string _path = executable_path() + "../data";
		if (is_directory(_path)) {
			dataPath = _path;
			return dataPath;
		}
	}

	return "data";
}

} // namespace mpc

