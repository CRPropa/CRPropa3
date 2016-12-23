#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#include <crpropa/Module.h>
#include <crpropa/magneticField/MagneticField.h>
#include <crpropa/Units.h>

namespace crpropa {

/**
 @class DiffusionSDE
 @brief Propagates candidates as pseudo(!)-particles.
 The time integration of SDEs is used to solve the transport equation. 
 * Here an Euler-Mayurama integration scheme is used. The diffusion tensor
 * can be anisotropic with respect to the magnetic field line coordinates.
 * The integration of field lines is done via the CK-algorithm.
 */


class DiffusionSDE : public Module{

private:
	    ref_ptr<MagneticField> field;
	    double minStep;
	    double maxStep;
	    double tolerance;
	    double epsilon;
	    double alpha;
	    double scale;
	    bool isTurbulent;
	    

public:
	    DiffusionSDE(ref_ptr<crpropa::MagneticField> field, double tolerance = 1e-4, 
	    		    double minStep=(10*pc), double maxStep=(1*kpc), double epsilon=0.1);

	    void process(crpropa::Candidate *candidate) const;
	   
	    void tryStep(const Vector3d &Pos, Vector3d &POut, Vector3d &PosErr, Vector3d &PosTest, Vector3d &TVec,Vector3d &NVec,Vector3d &BVec, double z, double propStep ) const;
	    
	    void calculateBTensor(double rig, double BTen[], Vector3d pos, Vector3d dir, double z) const;

	    void setMinimumStep(double minStep);
	    void setMaximumStep(double maxStep);
	    void setTolerance(double tolerance);
	    void setEpsilon(double kappa);
	    void setAlpha(double alpha);
	    void setScale(double Scale);
	    void setField(ref_ptr<crpropa::MagneticField> field);

	    double getMinimumStep() const;
	    double getMaximumStep() const;
	    double getTolerance() const;
	    double getEpsilon() const;
	    double getAlpha() const;
	    double getScale() const;
	    std::string getDescription() const;

}; 

} //namespace crpropa
