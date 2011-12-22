%module mpc

%include stl.i
%include std_set.i
%include std_multiset.i
%include std_map.i
%include std_pair.i
%include std_multimap.i
%include std_vector.i
%include std_string.i
%include stdint.i
%include std_container.i

%{
#include "mpc/interaction/Decay.h"
#include "mpc/interaction/ElectronPairProduction.h"
#include "mpc/interaction/PhotoDisintegration.h"

#include "mpc/magneticfield/MagneticField.h"
#include "mpc/magneticfield/MagneticFieldRing.h"
#include "mpc/magneticfield/SphMagneticField.h"
#include "mpc/magneticfield/TurbulentMagneticField.h"

#include "mpc/BreakCondition.h"
#include "mpc/Candidate.h"
#include "mpc/DeflectionCK.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/GlutDisplay.h"
#include "mpc/HepPID.h"
#include "mpc/MersenneTwister.h"
#include "mpc/Module.h"
#include "mpc/ModuleChain.h"
#include "mpc/Output.h"
#include "mpc/ParticleState.h"
#include "mpc/PhasePoint.h"
#include "mpc/Units.h"
#include "mpc/Vector3.h"
%}
 
/* Parse the header file to generate wrappers */
%include "mpc/Units.h"
%include "mpc/Vector3.h"
%include "mpc/ParticleState.h"
%include "mpc/Candidate.h"

%include "mpc/Module.h"

%include "mpc/magneticfield/MagneticField.h"
%include "mpc/magneticfield/MagneticFieldRing.h"
%include "mpc/magneticfield/SphMagneticField.h"
%include "mpc/magneticfield/TurbulentMagneticField.h"

%include "mpc/ExplicitRungeKutta.h"
%include "mpc/PhasePoint.h"

%include "mpc/BreakCondition.h"
%include "mpc/DeflectionCK.h"
%include "mpc/GlutDisplay.h"
%include "mpc/Output.h"

%include "mpc/ModuleChain.h"


%include "mpc/interaction/Decay.h"
%include "mpc/interaction/ElectronPairProduction.h"
%include "mpc/interaction/PhotoDisintegration.h"



 