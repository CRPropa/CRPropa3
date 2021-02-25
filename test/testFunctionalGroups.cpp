#include "CRPropa.h"
#include "gtest/gtest.h"

namespace crpropa {

/*
 * Functional test which calculates the particle's gyroradius in a uniform field
 * r_g = R / (B*c) = 1 EV / (1 nG * c) \approx 1.08*Mpc
 */
TEST(testFunctionalGroups, gyroradius) {
	double energy = 1*EeV;
	double field = 1*nG;

	ParticleState p;
	p.setId(nucleusId(1, 1));
	p.setEnergy(energy);
	p.setPosition(Vector3d(0, 0, 0));
	p.setDirection(Vector3d(0, 0, 1));

	ref_ptr<Candidate> c = new Candidate(p);
	ref_ptr<PropagationCK> propa = new PropagationCK(new UniformMagneticField(Vector3d(field, 0, 0)));
	ref_ptr<ParticleCollector> collector = new ParticleCollector();
	collector->setClone(true);
	ref_ptr<ModuleList> sim = new ModuleList();

	Vector3d pos;
	double max_y = 0;

	sim->add(propa);
	sim->add(new MaximumTrajectoryLength(10*Mpc));
	sim->add(collector);

	sim->run(c);

	for (ParticleCollector::iterator itr = collector->begin(); itr != collector->end(); ++itr){
		pos = (*(itr->get())).current.getPosition();
		if (max_y < pos.getY())
			max_y = pos.getY();
	}

	EXPECT_NEAR(max_y/2.0, energy/(field * c_light * eplus), 0.01*Mpc);

}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa
