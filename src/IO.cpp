#include "mpc/IO.h"

#include <kiss/io.h>

namespace mpc {

void write(kiss::Output &out, const mpc::Vector3 &vec3) {
	write(out, vec3.getX());
	write(out, vec3.getY());
	write(out, vec3.getZ());
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
	write(out, candidate.getStatus());
//write(out, candidate.getInteractionStates());
}

bool read(kiss::Input &in, mpc::Vector3 &vec3) {
	kiss::read(in, vec3[0]);
	kiss::read(in, vec3[1]);
	kiss::read(in, vec3[2]);
	return (in);
}

bool read(kiss::Input &in, mpc::ParticleState &state) {
	state.setId(kiss::read<int>(in));
	state.setEnergy(kiss::read<double>(in));
	state.setPosition(kiss::read<mpc::Vector3>(in));
	state.setDirection(kiss::read<mpc::Vector3>(in));
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
	candidate.setStatus((mpc::Candidate::Status) kiss::read<int>(in));
//read(in, candidate.getInteractionStates());
	return (in);
}

} // namespace mpc
