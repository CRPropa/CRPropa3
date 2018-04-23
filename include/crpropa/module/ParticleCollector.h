#ifndef CRPROPA_PARTICLECOLLECTOR_H
#define CRPROPA_PARTICLECOLLECTOR_H
#include <vector>
#include <string>

#include "crpropa/Module.h"
#include "crpropa/ModuleList.h"

namespace crpropa {
/**
 * \addtogroup Tools
 * \addtogroup Output
 * @{
 */

/**
 @class ParticleCollector
 @brief A helper ouput mechanism to keep candidates in-memory and directly transfer them to Python
 */
class ParticleCollector: public Module {
protected:
        typedef std::vector<ref_ptr<Candidate> > tContainer;
        mutable tContainer container;
        std::size_t nBuffer;
	bool clone;
	bool recursive;

public:
        ParticleCollector();
        ParticleCollector(const std::size_t nBuffer);
        ParticleCollector(const std::size_t nBuffer, const bool clone);
        ParticleCollector(const std::size_t nBuffer, const bool clone, const bool recursive);
        ~ParticleCollector();

        void process(Candidate *candidate) const;
	void process(ref_ptr<Candidate> c) const;
	void reprocess(Module *action) const;
	void dump(const std::string &filename) const;
	void load(const std::string &filename);

        std::size_t size() const;
	ref_ptr<Candidate> operator[](const std::size_t i) const;
        void clearContainer();

	std::string getDescription() const;
	std::vector<ref_ptr<Candidate> > getAll() const;
	void setClone(bool b);
	bool getClone() const;

	/** iterator goodies */
        typedef tContainer::iterator iterator;
        typedef tContainer::const_iterator const_iterator;
        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;

	/**
	 Retrieves the trajectory of a detected particle
	 Procedure: takes the initial state of the particle, re-runs the ModuleList for that particle and captures trajectory
	*/
	void getTrajectory(ModuleList *mlist, std::size_t i, Module *output) const;
	void getTrajectory(ref_ptr<ModuleList> mlist, std::size_t i, ref_ptr<Module> output) const;
};
/** @}*/

} // namespace crpropa

#endif // CRPROPA_PARTICLECOLLECTOR_H
