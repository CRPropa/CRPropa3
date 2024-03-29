/* 3. Pretty print for Python */
/*  __repr__ << getDescription */

__REPR__( crpropa::ParticleState );
__REPR__( crpropa::Candidate );
__REPR__( crpropa::Module );
__REPR__( crpropa::ModuleList );
__REPR__( crpropa::Source );
__REPR__( crpropa::SourceList );
__REPR__( crpropa::SourceFeature );
__REPR__( crpropa::Observer );
__REPR__( crpropa::ObserverFeature );

%pythoncode %{
  __version__ = g_GIT_DESC
%}

