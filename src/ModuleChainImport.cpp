#include "pugixml/pugixml.hpp"

#include "mpc/ModuleChain.h"
namespace mpc {

void import(ModuleChain &chain, pugi::xml_node &modules) {
	for (pugi::xml_node module = modules.child("module"); module; module =
			module.next_sibling("module")) {
		std::string type = module.attribute("type").value();
		unsigned int priority = module.attribute("priority").as_uint();

		std::cout << "Module: " << type << ", Prio: " << priority << std::endl;
	}
}
void import(ModuleChain &chain, std::string filename) {
	pugi::xml_document doc;
	pugi::xml_parse_result result = doc.load_file(filename.c_str());

	if (result) {

	}else
	{
	    std::cout << "XML [" << filename << "] parsed with errors, attr value: [" << doc.child("node").attribute("attr").value() << "]\n";
	    std::cout << "Error description: " << result.description() << "\n";
//	    std::cout << "Error offset: " << result.offset << " (error at [..." << (filename + result.offset) << "]\n\n";
	}

	pugi::xml_node node = doc.child("mpc");
	if (node)
		import(chain, node);
	else
		std::cout << "No mpc node found!" << std::endl;
}

} // namespace mpc
