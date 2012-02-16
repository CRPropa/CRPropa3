#include "mpc/XMLImport.h"
#include "mpc/ModuleChain.h"
#include "mpc/magneticField/magneticField.h"
#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/GlutDisplay.h"

#include "pugixml.hpp"

#include <stdexcept>
#include <sstream>

namespace mpc {

Module *ModuleFactory::create(const std::string type,
		pugi::xml_node &moduleNode) {
	std::map<std::string, ModuleProducer *>::iterator i;
	i = producers.find(type);
	if (i == producers.end())
		return 0;

	return i->second->create(moduleNode);
}

void ModuleFactory::registerProducer(const std::string &type,
		ModuleProducer *producer) {
	producers.insert(make_pair(type, producer));
}

XMLImport::XMLImport(ModuleChain *chain) :
		chain(chain) {

}

void XMLImport::import(const std::string &filename) {
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(filename.c_str());

	if (!result) {
		std::stringstream sstr;
		sstr << "XML [" << filename << "] parsed with errors, attr value: ["
				<< doc.child("node").attribute("attr").value() << "]\n";
		sstr << "Error description: " << result.description() << "\n";
		throw std::runtime_error(sstr.str());
	}

	pugi::xml_node node = doc.child("mpc");
	if (!node) {
		std::stringstream sstr;
		sstr << "XML [" << filename << "] no 'mpc' element found!\n";
		throw std::runtime_error(sstr.str());
	}

	import(node);
}

ModuleFactory &ModuleFactory::instance() {
	static ModuleFactory factory;
	return factory;
}

void XMLImport::import(pugi::xml_node &modules) {
	for (pugi::xml_node moduleNode = modules.child("module"); moduleNode;
			moduleNode = moduleNode.next_sibling("module")) {
		std::string type = moduleNode.attribute("type").value();
		unsigned int priority = moduleNode.attribute("priority").as_uint();

		Module *module = ModuleFactory::instance().create(type, moduleNode);
		if (module)
			chain->add(priority, module);
		else
			std::cout << "No module '" << type << "' found!" << std::endl;
	}
}

class DeflectionCKProducer: public ModuleProducer {
public:

	DeflectionCKProducer() :
			ModuleProducer("DeflectionCK") {

	}

	Module *create(pugi::xml_node &module) {
		UniformMagneticField *field = new UniformMagneticField(
				Vector3(0, 1, 0));
		return new DeflectionCK(field, DeflectionCK::RMS, 0.005);
	}
};

class GlutDislpayProducer: public ModuleProducer {
public:

	GlutDislpayProducer() :
			ModuleProducer("GluDisplay") {

	}

	Module *create(pugi::xml_node &module) {
		return new GlutDisplay();
	}
};

static GlutDislpayProducer _glut_display_producer;

static DeflectionCKProducer _deflection_ck_producer;

} // namespace mpc
