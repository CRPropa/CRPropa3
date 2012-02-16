#include "mpc/IO.h"

#include "kissIO.h"

namespace mpc {

using kiss::write;

void write(std::ostream &out, const mpc::Vector3 &vec3) {
	write(out, vec3.getX());
	write(out, vec3.getY());
	write(out, vec3.getZ());
}

void write(std::ostream &out, const mpc::ParticleState &state) {
	write(out, state.getId());
	write(out, state.getEnergy());
	write(out, state.getPosition());
	write(out, state.getDirection());
}

void write(std::ostream &out, const mpc::InteractionState &s) {
	write(out, s.distance);
	write(out, s.channel);
}

void write(std::ostream &out, const mpc::Candidate &candidate) {
	write(out, candidate.current);
	write(out, candidate.last);
	write(out, candidate.initial);
	write(out, candidate.getRedshift());
	write(out, candidate.getTrajectoryLength());
	write(out, candidate.getCurrentStep());
	write(out, candidate.getNextStep());
	write(out, candidate.getStatus());
//write(out, candidate.getInteractionStates());
}

using kiss::read;

bool read(std::istream &in, mpc::Vector3 &vec3) {
	read(in, vec3[0]);
	read(in, vec3[1]);
	read(in, vec3[2]);
	return (in);
}

bool read(std::istream &in, mpc::ParticleState &state) {
	state.setId(read<int>(in));
	state.setEnergy(read<double>(in));
	state.setPosition(read<mpc::Vector3>(in));
	state.setDirection(read<mpc::Vector3>(in));
	return (in);
}

bool read(std::istream &in, mpc::InteractionState &s) {
	read(in, s.distance);
	read(in, s.channel);
	return (in);
}

bool read(std::istream &in, mpc::Candidate &candidate) {
	read(in, candidate.current);
	read(in, candidate.last);
	read(in, candidate.initial);
	candidate.setRedshift(read<double>(in));
	candidate.setTrajectoryLength(read<double>(in));
	candidate.setCurrentStep(read<double>(in));
	candidate.setNextStep(read<double>(in));
	candidate.setStatus((mpc::Candidate::Status) read<int>(in));
//read(in, candidate.getInteractionStates());
	return (in);
}

} // namespace mpc
