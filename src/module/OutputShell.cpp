#include "crpropa/module/OutputShell.h"
#include "crpropa/Units.h"

#include <iomanip>

namespace crpropa {

void ShellOutput::process(Candidate* c) const {
#pragma omp critical(ShellOutput)
	{
		std::cout << std::fixed << std::showpoint << std::setprecision(3)
				<< std::setw(6);
		std::cout << c->getTrajectoryLength() / Mpc << " Mpc,  ";
		std::cout << c->getRedshift() << ",  ";
		std::cout << c->current.getId() << ",  ";
		std::cout << c->current.getEnergy() / EeV << " EeV,  ";
		std::cout << c->current.getPosition() / Mpc << " Mpc,  ";
		std::cout << c->current.getDirection();
		std::cout << std::endl;
	}
}

std::string ShellOutput::getDescription() const {
	return "Shell output";
}

void ShellOutput1D::process(Candidate* c) const {
#pragma omp critical(ShellOutput)
	{
		std::cout << std::fixed << std::showpoint << std::setprecision(3)
				<< std::setw(6);
		std::cout << c->current.getPosition().x / Mpc << " Mpc,  ";
		std::cout << c->getRedshift() << ",  ";
		std::cout << c->current.getId() << ",  ";
		std::cout << c->current.getEnergy() / EeV << " EeV";
		std::cout << std::endl;
	}
}

std::string ShellOutput1D::getDescription() const {
	return "Shell output for 1D";
}

void ShellPropertyOutput::process(Candidate* c) const {
	Candidate::PropertyMap::const_iterator i = c->properties.begin();
#pragma omp critical(ShellOutput)
	{
		for ( ; i != c->properties.end(); i++) {
			std::cout << "  " << i->first << ", " << i->second << std::endl;
		}
	}
}

std::string ShellPropertyOutput::getDescription() const {
	return "Shell property output";
}

} // namespace crpropa
