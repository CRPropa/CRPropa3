/// Example plugin for CRPropa.
///
/// Please consider sharing the awesome plugin with you fellow researchers by
/// creating a eperate repository for your project. We maintain a list of
/// plugins to CRPropa on our webpage and are happy to add a link to your
/// project, just send us: (name of the plugin, short description, url)

// start with header guards to ensure this file is only ever included once:
#ifndef MYPLUGIN_PLUGIN_H
#define MYPLUGIN_PLUGIN_H

#include <crpropa/Module.h>
#include <crpropa/Source.h>


// You should wrap everything in namespaces if possible, this leads to
// better usable code. Here you would call you module over myPlugin::MyModule
namespace myPlugin{

/// A custom C++ module
class MyModule : public crpropa::Module
{
public:
	/// The parent's constructor need to be called on initialization!
	MyModule();
	#if CRPROPA_VERSION<3003001
		void process(crpropa::Candidate* candidate) const override;
	#else
		void process(crpropa::ref_ptr<crpropa::Candidate> candidate) const override;
	#endif
};

}  // end namespace myPlugin


// You can also expand existing namespaces, here your new SourceFeature would be 
// callable like crpropa::AddMyProperty
namespace crpropa{

/// A custom source feature
class AddMyProperty: public SourceFeature
{
public:
	/// The parent's constructor need to be called on initialization!
	AddMyProperty();
	void prepareCandidate(Candidate &candidate) const;
};

}

#endif