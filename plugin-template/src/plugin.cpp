#include "myPlugin/plugin.h"

#include <string>


// also use the correct namespace in your .cpp files
namespace myPlugin{

// a global variable
const std::string myPropertyName = "counter";

// The parent's constructor need to be called on initialization!
MyModule::MyModule() : crpropa::Module()
{
	setDescription("MyPlugin::MyModule");
}

#if CRPROPA_VERSION<3003001
void MyModule::process(crpropa::Candidate* candidate) const
#else
void MyModule::process(crpropa::ref_ptr<crpropa::Candidate> candidate) const
#endif
{
	// To enable parallelization, the modules have to be stateless - the
	// process method should thus not modify internal variables!
	std::cout << "MyPlugin::MyModule::process() called\n";
	if(candidate->hasProperty(myPropertyName)) {
		uint32_t v = candidate->getProperty(myPropertyName);
		std::cout << " My property: " << myPropertyName << " = " <<  v << std::endl;
		candidate->setProperty(myPropertyName, crpropa::Variant::fromUInt32(v + 1));
	}
}

}

namespace crpropa{

// ------------------------------------------------------------------
// The parent's constructor need to be called on initialization!
AddMyProperty::AddMyProperty() : crpropa::SourceFeature()
{
	description = "AddMyProperty: Adds my property to the candidate";
}
void AddMyProperty::prepareCandidate(crpropa::Candidate &candidate) const
{
	candidate.setProperty(myPlugin::myPropertyName, crpropa::Variant::fromUInt32(0));
}

}