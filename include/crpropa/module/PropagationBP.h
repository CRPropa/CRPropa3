#ifndef CRPROPA_PROPAGATIONBP_H
#define CRPROPA_PROPAGATIONBP_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

    class PropagationBP: public Module {

        ref_ptr<MagneticField> field;
        double propStep;


    public:
        PropagationBP(ref_ptr<MagneticField> field = NULL,
                      double propStep = (0.1 * Mpc));

        void process(Candidate *candidate) const;

        void setField(ref_ptr<MagneticField> field);
        void setStep(double propStep);
        double getStep() const;
        std::string getDescription() const;
    };

}// namespace crpropa

#endif // CRPROPA_PROPAGATIONBP_H
