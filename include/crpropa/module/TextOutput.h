#ifndef CRPROPA_TEXTOUTPUT_H
#define CRPROPA_TEXTOUTPUT_H

#include "crpropa/module/Output.h"
#include "crpropa/module/ParticleCollector.h"

#include <fstream>

namespace crpropa {
/**
 * \addtogroup Output
 * @{
 */

/**
 @class TextOutput
 @brief Configurable plain text output for particle information.
 This type of output can also be used to generate a .tar.gz file if
 the library zlib is available. For details see:
 	http://zlib.net/
 */
class TextOutput: public Output {
protected:
	std::ostream *out;
	std::ofstream outfile;
	std::string filename;
	bool storeRandomSeeds;
	
	void printHeader() const;

public:
	/** Default constructor
	 */
	TextOutput();
	/** Constructor
	 @param outputType	type of output: Trajectory1D, Trajectory3D, Event1D, Event3D, Everything
	 */
	TextOutput(OutputType outputType);
	/** Constructor
	 @param out			output stream
	 */
	TextOutput(std::ostream &out);
	/** Constructor
	 @param out			output stream
	 @param outputType	type of output: Trajectory1D, Trajectory3D, Event1D, Event3D, Everything
	 */
	TextOutput(std::ostream &out, OutputType outputType);
	/** Constructor with the default OutputType (everything).
	 @param filename	string containing name of output text file
	 */
	TextOutput(const std::string &filename);
	/** Constructor
	 @param filename	string containing name of output text file
	 @param outputType	type of output: Trajectory1D, Trajectory3D, Event1D, Event3D, Everything
	 */
	TextOutput(const std::string &filename, OutputType outputType);
	/** Destructor
	 */
	~TextOutput();
	/** Whether to store the random seeds used in the simulation.
	 This enables reproducibility of each realisation of the simulation.
	 */
	void enableRandomSeeds() {storeRandomSeeds = true;};
	void close();
	void gzip();
	void process(Candidate *candidate) const;
	/** Loads a file to a particle collector.
	 This is useful for analysis involving, e.g., magnetic lenses.
	 @param filename	string containing the name of the file to be loaded
	 @param collector	object of type ParticleCollector that will store the information
	 */
	static void load(const std::string &filename, ParticleCollector *collector);
	std::string getDescription() const;

	void dumpIndexList(std::vector<int> indicies);
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_TEXTOUTPUT_H
