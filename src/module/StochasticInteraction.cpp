#include "crpropa/module/StochasticInteraction.h"

namespace crpropa {

void StochasticInteraction::process(Candidate* candidate) const {
	double step = candidate->getCurrentStep();
	while (step >= 0) {
		// get the interaction state, if there is one
		InteractionState interaction;
		bool noState = !candidate->getInteractionState(getDescription(),
				interaction);

		// if no interaction state, set a new one
		if (noState) {
			bool noNewState = !setNextInteraction(candidate, interaction);
			if (noNewState)
				return; // no new interaction; return
		}

		// if interaction distance not reached, reduce it and return
		if (interaction.distance > step) {
			interaction.distance -= step;
			candidate->limitNextStep(interaction.distance);
			candidate->setInteractionState(getDescription(), interaction);
			return;
		}

		// else: interaction distance reached; interact and repeat with remaining step
		performInteraction(candidate);
		step -= interaction.distance;
	}
}

} // namespace crpropa
