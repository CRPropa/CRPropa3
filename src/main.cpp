#include "crpropa/XmlExecute.h"

#include <iostream>
#include <exception>

int main(int argc, char **argv) {
	std::cerr <<  "!!! Deprecation Warning !!! Support for legacy XML steering will soon be removed from future CRPropa versions. Please switch to python based steering.\n";
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

	std::cerr <<  "!!! Deprecation Warning !!! Support for legacy XML steering will soon be removed from future CRPropa versions. Please switch to python based steering.\n";
	return 0;
}
