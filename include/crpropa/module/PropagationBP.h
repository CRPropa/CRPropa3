#ifndef CRPROPA_PROPAGATIONBP_H
#define CRPROPA_PROPAGATIONBP_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

    class PropagationBP: public Module {

        ref_ptr<MagneticField> field;
        double step;


    public:
        PropagationBP(ref_ptr<MagneticField> field = NULL,
                      double step = (0.1 * kpc));

        void process(Candidate *candidate) const;

        void setField(ref_ptr<MagneticField> field);
        void setStep(double step);
        double getStep() const;
        std::string getDescription() const;
    };

}// namespace crpropa

#endif // CRPROPA_PROPAGATIONBP_H
