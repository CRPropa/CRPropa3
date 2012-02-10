#ifndef MODULECHAINIMPORT_H_
#define MODULECHAINIMPORT_H_

#include <map>
#include <string>

#include "mpc/pugixml/pugixml.hpp"
#include "mpc/ModuleChain.h"

namespace mpc {

class ModuleProducer;

class ModuleFactory {
	std::map<std::string, ModuleProducer *> producers;
public:
	static ModuleFactory &instance();
	Module *create(const std::string type, pugi::xml_node &modules);
	void registerProducer(const std::string &type, ModuleProducer *producer);
};


class ModuleProducer {
public:

	inline ModuleProducer(const std::string &name) {
		ModuleFactory::instance().registerProducer(name, this);
	}

	virtual ~ModuleProducer() {
	}
	virtual Module*create(pugi::xml_node &modules) = 0;
};

class XMLImport {
	ModuleChain *chain;
public:
	XMLImport(ModuleChain *chain);
	void import(const std::string &filename);
	void import(pugi::xml_node &node);
};

} // namespace mpc

#endif /* MODULECHAINIMPORT_H_ */
