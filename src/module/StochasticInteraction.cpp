#include "mpc/module/StochasticInteraction.h"

namespace mpc {

void StochasticInteraction::process(Candidate* candidate) const {
	double step = candidate->getCurrentStep();
	InteractionState interaction;

	while (true) {
		// check if an interaction is set
		bool noState = !candidate->getInteractionState(getDescription(), interaction);
		if (noState) {
			// try to set a new interaction
			bool successful = setNextInteraction(candidate, interaction);
			if (not (successful))
				return; // no new interaction for this particle
		}

		// if not over, reduce distance and return
		if (interaction.distance > step) {
			interaction.distance -= step;
			candidate->limitNextStep(interaction.distance);
			candidate->setInteractionState(getDescription(), interaction);
			return;
		}

		// counter over: interact
		step -= interaction.distance;
		performInteraction(candidate);
	}
}

} // namespace mpc
