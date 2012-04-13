#ifndef COMMON_H_
#define COMMON_H_

#include <string>

namespace mpc {

std::string getDataPath(std::string filename);

inline int digit(const int &value, const int &d) {
	return (value % (d * 10)) / d;
}


} // namespace mpc

#endif /* COMMON_H_ */
