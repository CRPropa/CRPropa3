#ifndef CRPROPA_OUTPUTSHELL_H
#define CRPROPA_OUTPUTSHELL_H

#include "crpropa/Module.h"
#include "crpropa/Variant.h"

#include <unordered_map>

namespace crpropa {
/**
 * \addtogroup Output
 * @{
 */

/**
 @class ShellOutput
 @brief Show the trajectory in the shell.
 */
class ShellOutput: public Module {
public:
	void process(ref_ptr<Candidate> candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellOutput1D
 @brief Show the trajectory in the shell.
 */
class ShellOutput1D: public Module {
public:
	void process(ref_ptr<Candidate> candidate) const;
	std::string getDescription() const;
};

/**
 @class ShellPropertyOutput
 @brief Show the candidate properties in the shell.
 */
class ShellPropertyOutput: public Module {
public:
	typedef std::unordered_map<std::string, Variant> PropertyMap;
	void process(ref_ptr<Candidate> candidate) const;
	std::string getDescription() const;
};
/** @}*/

} // namespace cprpropa

#endif // CRPROPA_OUTPUTSHELL_H
