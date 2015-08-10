/* 3. Pretty print for Python */
/*  __repr__ << getDescription */


__REPR__( crpropa::ParticleState );
__REPR__( crpropa::Candidate );
__REPR__( crpropa::Module );
__REPR__( crpropa::ModuleList );

__REPR__( crpropa::Source );

__REPR__( crpropa::SourceList );
__REPR__( crpropa::SourceParticleType );
__REPR__( crpropa::SourceMultipleParticleTypes );
__REPR__( crpropa::SourceEnergy );
__REPR__( crpropa::SourcePowerLawSpectrum );
__REPR__( crpropa::SourceMultiplePositions );
__REPR__( crpropa::SourceComposition );
__REPR__( crpropa::SourcePosition );
__REPR__( crpropa::SourceUniform1D );
__REPR__( crpropa::SourceUniformBox );
__REPR__( crpropa::SourceUniformShell );
__REPR__( crpropa::SourceUniformSphere );
__REPR__( crpropa::SourceDensityGrid );
__REPR__( crpropa::SourceDensityGrid1D );
__REPR__( crpropa::SourceDirection );
__REPR__( crpropa::SourceIsotropicEmission );
__REPR__( crpropa::SourceEmissionCone );
__REPR__( crpropa::SourceRedshift );
__REPR__( crpropa::SourceRedshift1D );
__REPR__( crpropa::SourceUniformRedshift );

__REPR__( crpropa::Observer );
__REPR__( crpropa::ObserverPoint );
__REPR__( crpropa::ObserverSmallSphere );
__REPR__( crpropa::ObserverLargeSphere );
__REPR__( crpropa::ObserverRedshiftWindow );
__REPR__( crpropa::ObserverNucleusVeto );
__REPR__( crpropa::ObserverNeutrinoVeto );
__REPR__( crpropa::ObserverPhotonVeto );

VECTOR3__REPR__( crpropa::Vector3 );

%pythoncode %{
    DeflectionCK = PropagationCK  # legacy name
%}

