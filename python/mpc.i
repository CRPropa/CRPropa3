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
#include "mpc/module/NuclearDecay.h"
#include "mpc/module/ElectronPairProduction.h"
#include "mpc/module/PhotoDisintegration.h"
#include "mpc/module/BreakCondition.h"
#include "mpc/module/DeflectionCK.h"
#include "mpc/module/Output.h"

#include "mpc/magneticField/magneticField.h"
#include "mpc/magneticField/uniformMagneticField.h"
#include "mpc/magneticField/turbulentMagneticFieldGrid.h"

#include "mpc/Candidate.h"
#include "mpc/ExplicitRungeKutta.h"
#include "mpc/module/GlutDisplay.h"
#include "mpc/HepPID.h"
#include "mpc/Random.h"
#include "mpc/Module.h"
#include "mpc/ModuleChain.h"
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

%include "mpc/magneticField/magneticField.h"
%include "mpc/magneticField/uniformMagneticField.h"
%include "mpc/magneticField/turbulentMagneticFieldGrid.h"

%include "mpc/ExplicitRungeKutta.h"
%include "mpc/PhasePoint.h"

%include "mpc/module/BreakCondition.h"
%include "mpc/module/DeflectionCK.h"
%include "mpc/module/GlutDisplay.h"
%include "mpc/module/Output.h"
%include "mpc/module/NuclearDecay.h"
%include "mpc/module/ElectronPairProduction.h"
%include "mpc/module/PhotoDisintegration.h"

%include "mpc/ModuleChain.h"


 
