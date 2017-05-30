#include "crpropa/module/Output.h"
#include "crpropa/Units.h"

#include <stdexcept>

namespace crpropa {

Output::Output() : lengthScale(Mpc), energyScale(EeV), oneDimensional(false), weights(false), count(0) {
	enableAll();
}

Output::Output(OutputType outputtype) : lengthScale(Mpc), energyScale(EeV), oneDimensional(false), weights(false), count(0) {
		setOutputType(outputtype);
}

void Output::modify() {
	if (count > 0)
		throw std::runtime_error("Output: cannot change Output parameters after data has been written to file.");
}

void Output::process(Candidate *c) const {
	count++;
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
		// (W), D, ID, E, ID0, E0
		if (weights)
			set(WeightColumn, true);	
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
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
		// (W), D, ID, E, X, Y, Z, Px, Py, Pz, ID0, E0, X0, Y0, Z0
		if (weights)
			set(WeightColumn, true);
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(CurrentPositionColumn, true);
		set(CurrentDirectionColumn, true);
		set(SourceIdColumn, true);
		set(SourceEnergyColumn, true);
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

void Output::enableWeights(bool value) {
	modify();
	weights = value;
	fields.set(WeightColumn, value);
	std::cout << "weights = " << value << std::endl;
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

size_t Output::getCount() const {
	return count;
}

} // namespace crpropa
