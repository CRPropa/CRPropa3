#ifndef CRPROPA_ABSTRACT_OUTPUT_H
#define CRPROPA_ABSTRACT_OUTPUT_H

#include "crpropa/Module.h"

#include <bitset>

namespace crpropa {

/**
 @class Output
 @brief Configurable output base class.
 */
class Output: public Module {
protected:
	double lengthScale, energyScale;
	std::bitset<64> fields;
	bool oneDimensional;
	bool weights;
	mutable size_t count;
	
	void modify();

public:
	enum OutputColumn {
		WeightColumn,
		TrajectoryLengthColumn,
		RedshiftColumn,
		CurrentIdColumn,
		CurrentEnergyColumn,
		CurrentPositionColumn,
		CurrentDirectionColumn,
		SourceIdColumn,
		SourceEnergyColumn,
		SourcePositionColumn,
		SourceDirectionColumn,
		CreatedIdColumn,
		CreatedEnergyColumn,
		CreatedPositionColumn,
		CreatedDirectionColumn,
		SerialNumberColumn
	};

	enum OutputType {
		Trajectory1D,
		Trajectory3D,
		Event1D,
 		Event3D,
		Everything
	};

	Output();
	Output(OutputType outputtype);
	
	void setEnergyScale(double scale);
	void setLengthScale(double scale);

	void enableWeights(bool value);

	void setOutputType(OutputType outputtype);
	void set(OutputColumn field, bool value);
	void enable(OutputColumn field);
	void disable(OutputColumn field);
	void enableAll();
	void disableAll();
	void set1D(bool value);
	size_t getCount() const;

	void process(Candidate *) const;
};

} // namespace crpropa

#endif // CRPROPA_ABSTRACT_OUTPUT_H
