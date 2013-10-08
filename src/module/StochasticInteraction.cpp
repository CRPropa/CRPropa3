#include "crpropa/module/StochasticInteraction.h"

namespace crpropa {

void StochasticInteraction::process(Candidate* c) const {
	double step = c->getCurrentStep();

	while (step >= 0) {
		// get the interaction, if there is one
		InteractionState interaction;
		bool hasOne = c->getInteractionState(getDescription(), interaction);

		// if no interaction, try to set a new one
		if (not hasOne) {
			bool hasOneNow = randomInteraction(c, interaction);
			if (not hasOneNow) {
				return; // no new interaction, nothing more to do
			}
		}

		// if interaction distance not reached, reduce it and return
		if (interaction.distance > step) {
			interaction.distance -= step;
			c->limitNextStep(interaction.distance);
			c->setInteractionState(getDescription(), interaction);
			return;
		}

		// if interaction distance reached, interact and repeat with remaining step
		performInteraction(c, interaction);
		c->clearInteractionStates();
		step -= interaction.distance;
	}
}

} // namespace crpropa
