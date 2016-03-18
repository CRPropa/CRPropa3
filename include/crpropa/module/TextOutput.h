#ifndef CRPROPA_TEXTOUTPUT_H
#define CRPROPA_TEXTOUTPUT_H

#include "crpropa/module/Output.h"

#include <fstream>

namespace crpropa {

/**
 @class TextOutput
 @brief Configurable plain text output for cosmic ray information.
 */
class TextOutput: public Output {
protected:
	std::ostream *out;
	std::ofstream outfile;
	std::string filename;

	void printHeader() const;

public:
	TextOutput();
	TextOutput(OutputType outputtype);
	TextOutput(std::ostream &out);
	TextOutput(std::ostream &out, OutputType outputtype);
	TextOutput(const std::string &filename);
	TextOutput(const std::string &filename, OutputType outputtype);
	~TextOutput();

	void close();
	void gzip();

	void process(Candidate *candidate) const;
	std::string getDescription() const;
};

} // namespace crpropa

#endif // CRPROPA_TEXTOUTPUT_H
