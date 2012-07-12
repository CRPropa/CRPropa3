#include "mpc/XmlExecute.h"

#include <iostream>
#include <exception>

using namespace mpc;
using namespace std;

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "mpc-run <XML File>" << endl;
		return 1;
	}

	try {
		XmlExecute exe;
		if (!exe.load(argv[1]))
			return 1;
		exe.run();
	} catch (std::exception &e) {
		cerr << e.what() << endl;
		return 1;
	}

	return 0;
}
