#include "mpc/module/StochasticInteraction.h"

namespace mpc {

StochasticInteraction::StochasticInteraction(std::string name) :
		name(name) {
}

void StochasticInteraction::process(Candidate* candidate) const {
	double step = candidate->getCurrentStep();
	InteractionState interaction;

	while (true) {
		// check if an interaction is set
		bool noState = !candidate->getInteractionState(name, interaction);
		if (noState) {
			// try to set a new interaction
			bool successful = setNextInteraction(candidate);
			if (not (successful))
				return;
			// get the new interaction
			candidate->getInteractionState(name, interaction);
		}

		// if not over, reduce distance and return
		if (interaction.distance > step) {
			interaction.distance -= step;
			candidate->limitNextStep(interaction.distance);
			candidate->setInteractionState(name, interaction);
			return;
		}

		// counter over: interact
		step -= interaction.distance;
		performInteraction(candidate);
	}
}

} // namespace mpc
