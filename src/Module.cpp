#include "mpc/Module.h"

#include <typeinfo>

namespace mpc {

std::string Module::getDescription() const {
	const std::type_info &info = typeid(*this);
	return info.name();
}

}
