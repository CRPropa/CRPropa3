#include <vector>

#include "crpropa/Module.h"

namespace crpropa {

/**
 @class ParticleCollector
 @brief A helper ouput mechanism to directly transfer candidates to Python
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
	void reprocess(Module *action) const;
        std::size_t getCount() const;
	ref_ptr<Candidate> operator[](const std::size_t i) const;
        void clearContainer();
        std::string getDescription() const;
	std::vector<ref_ptr<Candidate> > getAll() const;

	// iterator goodies
        typedef tContainer::iterator iterator;
        typedef tContainer::const_iterator const_iterator;
        iterator begin();
        const_iterator begin() const;
        iterator end();
        const_iterator end() const;
};

} // namespace crpropa
