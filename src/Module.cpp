#include "mpc/Module.h"

#include <typeinfo>

namespace mpc {

Module::Module() {
	const std::type_info &info = typeid(*this);
	setDescription(info.name());
}

std::string Module::getDescription() const {
	return description;
}

void Module::setDescription(const std::string &d) {
	description = d;
}

} // namespace mpc
