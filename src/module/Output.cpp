#include "crpropa/module/Output.h"
#include "crpropa/Units.h"

//#include <stdio.h>
//#include <stdexcept>

namespace crpropa {

Output::Output() : lengthScale(Mpc), energyScale(EeV), oneDimensional(false), begun(false), ended(false) {
	enableAll();
}

Output::Output(OutputType outputtype) : lengthScale(Mpc), energyScale(EeV), oneDimensional(false), begun(false), ended(false) {
		setOutputType(outputtype);
}

void Output::modify() {
	if (begun || ended)
		throw std::runtime_error("Output: cannot change Output parameters after begin/endRun.");
}

void Output::setOutputType(OutputType outputtype) {
	modify();
	if (outputtype == Trajectory1D) {
		// X, ID, E
		set(CurrentPositionColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set1D(true);
	} else if (outputtype == Event1D) {
		// ID, E, D, ID0, E0
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(TrajectoryLengthColumn, true);
		set(SourceIdColumn, true);
		set(SourceEnergyColumn, true);
		set1D(true);
	} else if (outputtype == Trajectory3D) {
		// D, ID, E, X, Y, Z, Px, Py, Pz
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(CurrentPositionColumn, true);
		set(CurrentDirectionColumn, true);
		set1D(false);
	} else if (outputtype == Event3D) {
		// D, ID, ID0, E, E0, X, Y, Z, X0, Y0, Z0, Px, Py, Pz, P0x, P0y, P0z,z
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(SourceIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(SourceEnergyColumn, true);
		set(CurrentPositionColumn, true);
		set(SourcePositionColumn, true);
		set(CurrentDirectionColumn, true);
		set(SourcePositionColumn, true);
		set1D(false);
	} else if (outputtype == Everything) {
		enableAll();
		set1D(false);
	} else {
		throw std::runtime_error(
				"TextOutput: unknown output type");
	}
}

void Output::setEnergyScale(double scale) {
	modify();
	energyScale = scale;
}

void Output::setLengthScale(double scale) {
	modify();
	lengthScale = scale;
}

void Output::set1D(bool value) {
	modify();
	oneDimensional = value;
}

void Output::enable(OutputColumn field) {
	modify();
	fields.set(field, true);
}

void Output::disable(OutputColumn field) {
	modify();
	fields.set(field, false);
}

void Output::set(OutputColumn field, bool value) {
	modify();
	fields.set(field, value);
}

void Output::enableAll() {
	modify();
	fields.set();
}

void Output::disableAll() {
	modify();
	fields.reset();
}

void  Output::process(Candidate *candidate) const {
	if (!begun)
		throw std::runtime_error("Output: call beginRun before using Output");
	if (ended)
		throw std::runtime_error("Output: cannot process output after endRun was called.");
}

void  Output::beginRun() {
	begun = true;
}

void  Output::endRun() {
	ended = true;
}

} // namespace crpropa
