#ifndef ELECA_COMMON_H_
#define ELECA_COMMON_H_

namespace eleca {

double z2Mpc(double z);
double Mpc2z(double D);
double Uniform(double min, double max);

void setUniformCallback(double (*Uniform)(double min, double max));

} // namespace eleca

#endif // ELECA_COMMON_H_
