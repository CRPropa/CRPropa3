#include "crpropa/XmlExecute.h"

#include <iostream>
#include <exception>

int main(int argc, char **argv) {
	if (argc < 2) {
		std::cout << "crpropa-run <XML File>" << std::endl;
		return 1;
	}

	try {
		crpropa::XmlExecute exe;
		if (!exe.load(argv[1]))
			return 1;
		exe.run();
	} catch (std::exception &e) {
		std::cerr << e.what() << std::endl;
		return 1;
	}

	return 0;
}
