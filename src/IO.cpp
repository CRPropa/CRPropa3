#include "mpc/IO.h"

namespace mpc {

void write(kiss::Output &out, const mpc::Vector3d &vec3) {
	write(out, vec3.x);
	write(out, vec3.y);
	write(out, vec3.z);
}

void write(kiss::Output &out, const mpc::ParticleState &state) {
	write(out, state.getId());
	write(out, state.getEnergy());
	write(out, state.getPosition());
	write(out, state.getDirection());
}

void write(kiss::Output &out, const mpc::InteractionState &s) {
	write(out, s.distance);
	write(out, s.channel);
}

void write(kiss::Output &out, const mpc::Candidate &candidate) {
	write(out, candidate.current);
	write(out, candidate.initial);
	write(out, candidate.getRedshift());
	write(out, candidate.getTrajectoryLength());
	write(out, candidate.getCurrentStep());
	write(out, candidate.getNextStep());
}

bool read(kiss::Input &in, mpc::Vector3d &vec3) {
	kiss::read(in, vec3.x);
	kiss::read(in, vec3.y);
	kiss::read(in, vec3.z);
	return (in);
}

bool read(kiss::Input &in, mpc::ParticleState &state) {
	state.setId(kiss::read<int>(in));
	state.setEnergy(kiss::read<double>(in));
	state.setPosition(kiss::read<mpc::Vector3d>(in));
	state.setDirection(kiss::read<mpc::Vector3d>(in));
	return (in);
}

bool read(kiss::Input &in, mpc::InteractionState &s) {
	kiss::read(in, s.distance);
	kiss::read(in, s.channel);
	return (in);
}

bool read(kiss::Input &in, mpc::Candidate &candidate) {
	read(in, candidate.current);
	read(in, candidate.initial);
	candidate.setRedshift(kiss::read<double>(in));
	candidate.setTrajectoryLength(kiss::read<double>(in));
	candidate.setCurrentStep(kiss::read<double>(in));
	candidate.setNextStep(kiss::read<double>(in));
//	candidate.setStatus((mpc::Candidate::Status) kiss::read<int>(in));
//  read(in, candidate.getInteractionStates());
	return (in);
}

} // namespace mpc
