#ifndef CRPROPA_OUTPUT_H
#define CRPROPA_OUTPUT_H

#include "crpropa/Module.h"
#include "crpropa/Variant.h"

#include <bitset>
#include <vector>
#include <string>

namespace crpropa {

/**
 * \addtogroup Output
 * @{
 */

/**
 @class Output
 @brief Configurable output base class.
 The names of each quantity are provided in the table below.
 The left column corresponds to the labels of each quantity. 
 They are printed in the output files either as comments (in a text file) 
 or as columns (in a HDF5 file). The right columns are the names of each
 column for internal access.
 . D			 TrajectoryLengthColumn
 . SN			 SerialNumberColumn
 . ID			 CurrentIdColumn
 . E			 CurrentEnergyColumn
 . X/Y/Z		 CurrentPositionColumn
 . Px/Py/Pz		 CurrentDirectionColumn
 . SN0			 SourceSerialNumberColumn
 . ID0			 SourceIdColumn
 . E0			 SourceEnergyColumn
 . X0/Y0/Z0		 SourcePositionColumn
 . P0x/P0y/P0z	 SourceDirectionColumn
 . SN1			 CreatedSerialNumberColumn
 . ID1			 CreatedIdColumn
 . E1			 CreatedEnergyColumn
 . X1/Y1/Z1		 CreatedPositionColumn
 . P1x/P1y/P1z	 CreatedDirectionColumn
 . z			 RedshiftColumn
 . tag			 CandidateTagColumn
 . weight		 WeightColumn

 Some output types are pre-defined: 
 . Trajectory1D
 . Trajectory3D
 . Event1D
 . Event3D
 . Everything
 They can be easily customised by enabling/disabling specific columns.
 */
class Output: public Module {
protected:
	double lengthScale, energyScale;
	std::bitset<64> fields;

	struct Property
	{
		std::string name;
		std::string comment;
		Variant defaultValue;
	};
	std::vector<Property> properties;

	bool oneDimensional;
	mutable size_t count;

	void modify();

public:
	enum OutputColumn {
		TrajectoryLengthColumn,
		ColumnDensityColumn,
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
		CandidateTagColumn,
		SerialNumberColumn,
		WeightColumn
	};
	enum OutputType {
		Trajectory1D,
		Trajectory3D,
		Event1D,
 		Event3D,
		Everything
	};

	std::string OutputTypeName(OutputType outputtype);
	const std::string outputName;
	
	/** Default constructor. Output contains all the information available.
	 Units of energy and length are, by default, EeV and Mpc.
	 This can be changed with setEnergyScale and setLengthScale.
	 */
	Output();
	/** General constructor.
	 Units of energy and length are, by default, EeV and Mpc.
	 This can be changed with setEnergyScale and setLengthScale.
	 @param outputType	type of output: Trajectory1D, Trajectory3D, Event1D, Event3D, Everything
	 */
	Output(OutputType outputType);

	/** Set energy scale.
	 @param scale	energy scale (scale = 1 corresponds to 1 Joule)
	 */
	void setEnergyScale(double scale);
	/** Set length scale.
	 @param scale	length scale (scale = 1 corresponds to 1 meter)
	 */
	void setLengthScale(double scale);
	/** Set type of output.
	 @param outputType	type of output: Trajectory1D, Trajectory3D, Event1D, Event3D, Everything
	 */
	void setOutputType(OutputType outputType);
	/** Determines whether a given column will be displayed in the output.
	 @param field	name of the field to be added/removed from output
	 @param values	boolean flag adding (true) or removing (false) the field
	 */
	void set(OutputColumn field, bool value);
	/** Add a property to output. 
	 Default value is required to assign a type in the output.
	 @param property	string containing name of property
	 @param default		default value of property
	 @param	comment		string with a comment
	 */
	void enableProperty(const std::string &property, const Variant& defaultValue, const std::string &comment = "");
	/** Enable specific column in the output.
	 @param field	name of the field to be enabled
	 */
	void enable(OutputColumn field);
	/** Disable specific column in the output.
	 @param field	name of the field to be disabled
	 */
	void disable(OutputColumn field);
	/** Enable all fields.
	 Essentially a wrapper for set(field, true).
	 */
	void enableAll();
	/** Disable all fields.
	 Essentially a wrapper for set(field, false). 
	 */
	void disableAll();
	/** If true, output is of 1D type. 
	 3D quantities such as vectors will be reduced to the relevant components (x, by default).
	 @param value	boolean flag
	 */
	void set1D(bool value);
	/** Returns the size of the output
	 */
	size_t size() const;

	void process(Candidate *) const;
};

/** @}*/

} // namespace crpropa

#endif // CRPROPA_OUTPUT_H
