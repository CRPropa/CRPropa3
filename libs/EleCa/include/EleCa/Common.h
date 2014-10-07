#ifndef ELECA_COMMON_H_
#define ELECA_COMMON_H_

namespace eleca {

double z2Mpc(double z);
double Mpc2z(double D);
double Uniform(double min, double max);

// set the seed for the random generator. If 0, current ime is used
void setSeed(long int seedval=0);

void setUniformCallback(double (*Uniform)(double min, double max));


// integer pow implementation as template that is evaluated at compile time
template <unsigned int exponent>
inline double pow_integer(double base)
{
  return pow_integer<exponent >> 1>(base*base) * (((exponent & 1) > 0) ? base : 1);
}

template <>
inline double pow_integer<0>(double base)
{
  return 1;
}


} // namespace eleca

#endif // ELECA_COMMON_H_
