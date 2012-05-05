#ifndef MPC_IO_H_
#define MPC_IO_H_

#include "mpc/Candidate.h"

#include <kiss/io.h>

namespace mpc {

void write(kiss::Output &out, const mpc::Vector3d &vec3);

void write(kiss::Output &out, const mpc::ParticleState &state);

void write(kiss::Output &out, const mpc::InteractionState &s);

void write(kiss::Output &out, const mpc::Candidate &candidate);

bool read(kiss::Input &in, mpc::Vector3d &vec3);

bool read(kiss::Input &in, mpc::ParticleState &state);

bool read(kiss::Input &in, mpc::InteractionState &s);

bool read(kiss::Input &in, mpc::Candidate &candidate);

} // namespace mpc

#endif /* MPC_IO_H_ */
