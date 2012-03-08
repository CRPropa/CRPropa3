#include "mpc/module/common.h"

#include <stdlib.h>

namespace mpc {

std::string getDataPath(std::string filename) {
	std::string dataPath = "data";
	if (getenv("MPC_DATA_PATH"))
		dataPath = getenv("MPC_DATA_PATH");
	dataPath += "/" + filename;
	return dataPath;
}

} // namespace mpc
