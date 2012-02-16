#ifndef MPC_IO_H_
#define MPC_IO_H_

#include "mpc/Candidate.h"

#include <iostream>

namespace mpc {

void write(std::ostream &out, const mpc::Vector3 &vec3);

void write(std::ostream &out, const mpc::ParticleState &state);

void write(std::ostream &out, const mpc::InteractionState &s);

void write(std::ostream &out, const mpc::Candidate &candidate);

bool read(std::istream &in, mpc::Vector3 &vec3);

bool read(std::istream &in, mpc::ParticleState &state);

bool read(std::istream &in, mpc::InteractionState &s);

bool read(std::istream &in, mpc::Candidate &candidate);

} // namespace mpc

#endif /* MPC_IO_H_ */
