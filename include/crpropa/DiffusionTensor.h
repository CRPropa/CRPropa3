#ifndef CRPROPA_DIFFUSIONCOEFFICENT_H
#define CRPROPA_DIFFUSIONCOEFFICENT_H

#include "crpropa/Candidate.h"
#include "crpropa/Referenced.h"

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
	void setDescription(const std::string &description);
    virtual void calculateBTensor(double BTen[], Candidate *cand){};
    std::string getDescription() const;
};


class QLTDiffusion: public DiffusionTensor {
    private:
        double epsilon;
        double kappa0;
        double alpha;

    public:
        QLTDiffusion(double epsilon = 0.1 , double kappa0 =6.1e24, double alpha= (1./3.) );

        void calculateBTensor(double BTen[], Candidate *cand);

        void setEpsilon(double epsilon);
        void setKappa0(double kappa0);
        void setAlpha(double alpha);
        void setDescription();

        double getEpsilon() const;
        double getAlpha() const;
        double getKappa0() const;
	    std::string getDescription() const;


};



} // namespace


#endif // CRPROPA_DIFFUSIONCOEFFICENT_H