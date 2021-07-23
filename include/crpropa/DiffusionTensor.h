#ifndef CRPROPA_DIFFUSIONCOEFFICENT_H
#define CRPROPA_DIFFUSIONCOEFFICENT_H

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"

#include <string>
#include <cmath>

namespace crpropa{

/**
 @class DiffusionTensor
 @brief Abstract base class for diffusion tensors
*/
class DiffusionTensor: public Referenced {
    private:
        std::string description;
    public:
        virtual ~DiffusionTensor(){
        }
	//void setDescription(const std::string &description);
    virtual double getKappaParallel(Candidate *cand){
        return 0;
    };
    virtual double getKappaPerpendicular(Candidate *cand){
        return 0;
    };
    virtual double getKappaPerpendicular2(Candidate *cand){
        return 0;
    };
    virtual double getEpsilon() {};
    //std::string getDescription() const;
};

/**
 @class QLTDiffusion
 @brief Diffusion Tensor for the quasi-linear theory.  
 */
class QLTDiffusion: public DiffusionTensor {
    private:
        double epsilon; // ratio between perpendicular and parallel diffusion coefficent
        double kappa0;  // norm value for the diffusioncoefficent at rigidity 4 GV
        double alpha;   // spectral index of the diffusion coefficent

    public:
        QLTDiffusion(double epsilon = 0.1 , double kappa0 = 6.1e24, double alpha = (1./3.) );
        
        double getKappaParallel(Candidate *cand);
        double getKappaPerpendicular(Candidate *cand);
        double getKappaPerpendicular2(Candidate *cand);

        void setEpsilon(double epsilon);
        void setKappa0(double kappa0);
        void setAlpha(double alpha);
        //void setDescription();

        double getEpsilon() const;
        double getAlpha() const;
        double getKappa0() const;
	    std::string getDescription() const;
};

class QLTTurbulent: public DiffusionTensor{
    private:
	    ref_ptr<MagneticField> backgroundField;
        ref_ptr<TurbulentField> turbulentField;
        double kappa0;      // value to norm the diffusioncoefficent at a rigidity of 4 GV
        double alphaPara;   // spectral index for the parallel component
        double alphaPerp;   // spectral index for the perpendicular component
        double normTurbulence;  // value to norm the turbulence (probably at earth)

    public:
        QLTTurbulent(ref_ptr<MagneticField> background, ref_ptr<TurbulentField> turbulent, double kappa0 = 6.1e24, double alphaPara=(1./3.), double alphaPerp=(1./3.));

        double getKappaParallel(Candidate *cand) const;
        double getKappaPerpendicular(Candidate *cand) const;
        double getKappaPerpendicular2(Candidate *cand) const;

        double getKappa0() const;
        double getAlphaPara() const;
        double getAlphaPerp() const;
        double getNormTurbulence() const;

        void setKappa0(double kappa0);
        void setAlphaPara(double alpha);
        void setAlphaPerp(double alpha);
        void setAlpha(double alpha);
        void setNormTurbulence(double eta);
        void normToEarthPosition(Vector3d posEarth = Vector3d(-8.5*kpc, 0., 0.));

        std::string getDescription() const;
};

class QLTRigidity: public DiffusionTensor{
    private:
        ref_ptr<MagneticField> backgroundField;
        ref_ptr<TurbulentField> turbulentField;
        bool hasTurbulentField;
        double kappa0;
        double normEta;
        double normB; 
        double alphaPara;
        double alphaPerp;
        double correlationLength = 59.96077*pc;
        Vector3d normPos; // position where the diffusion coefficent is normed. default at earth Vector3d(-8.5*kpc, 0, 0)
        double calculateLamorRadius(ParticleState &state) const;

    public:
        QLTRigidity(ref_ptr<MagneticField> magField, ref_ptr<TurbulentField> turbField, double kappa0=6.1e24, double alphaPara=(1./3.), double alpaPerp=(1./3.));
        QLTRigidity(ref_ptr<MagneticField> field, double kappa0=6.1e24, double alphaPara=(1./3.), double alphaPerp=(1./3.));

        void setMagneticField(ref_ptr<MagneticField> field);
        void setTurbulentField(ref_ptr<TurbulentField> field);
        void setKappa0(double kappa0);
        void setAlphaPara(double aPara);
        void setAlphaPerp(double aPerp);
        void setAlpha(double alpha);
        void normToPosition(const Vector3d &pos= Vector3d(-8.5*kpc, 0., 0.));
        void setNormEta(double eta);
        void setNormB(double B);
        void setHasTurbulentField(bool use);
};

} // namespace


#endif // CRPROPA_DIFFUSIONCOEFFICENT_H