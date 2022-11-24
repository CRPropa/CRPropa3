#include "crpropa/module/Output.h"
#include "crpropa/Units.h"

#include <stdexcept>

namespace crpropa {

Output::Output() : outputName(OutputTypeName(Everything)), lengthScale(Mpc), energyScale(EeV), oneDimensional(false), count(0) {
	enableAll();
}

Output::Output(OutputType outputType) : outputName(OutputTypeName(outputType)), lengthScale(Mpc), energyScale(EeV), oneDimensional(false), count(0) {
	setOutputType(outputType);
}

std::string Output::OutputTypeName(OutputType outputType) {
	if (outputType == Trajectory1D)
		return "Trajectory1D";
	if (outputType == Event1D)
		return "Event1D";
	if (outputType == Trajectory3D)
		return "Trajectory3D";
	if (outputType == Event3D)
		return "Event3D";
	return "Everything";
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
		// D, ID, E, ID0, E0, tag
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(SourceIdColumn, true);
		set(SourceEnergyColumn, true);
		set(CandidateTagColumn, true);
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
		// D, ID, E, X, Y, Z, Px, Py, Pz, ID0, E0, X0, Y0, Z0, P0x, P0y, P0z, tag
		set(TrajectoryLengthColumn, true);
		set(CurrentIdColumn, true);
		set(CurrentEnergyColumn, true);
		set(CurrentPositionColumn, true);
		set(CurrentDirectionColumn, true);
		set(SourceIdColumn, true);
		set(SourceEnergyColumn, true);
		set(SourcePositionColumn, true);
		set(SourceDirectionColumn, true);
		set(CandidateTagColumn, true);
		set1D(false);
	} else if (outputtype == Everything) {
		enableAll();
		set1D(false);
	} else {
		throw std::runtime_error("Output: unknown output type");
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

size_t Output::size() const {
	return count;
}

void Output::enableProperty(const std::string &property, const Variant &defaultValue, const std::string &comment) {
	modify();
	Property prop;
	prop.name = property;
	prop.comment = comment;
	prop.defaultValue = defaultValue;
	properties.push_back(prop);
};

} // namespace crpropa
