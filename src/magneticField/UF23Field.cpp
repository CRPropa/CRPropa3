#include "crpropa/magneticField/UF23Field.h"

#include <exception>
#include <limits>
#include <string>
#include <cmath>

// local helper functions and constants
namespace uf23 {

  template<typename T>
  crpropa::Vector3d CylToCart(const T v, const double cosPhi, const double sinPhi)
  {
    return crpropa::Vector3d(v[0] * cosPhi - v[1] * sinPhi,
                             v[0] * sinPhi + v[1] * cosPhi,
                             v[2]);
  }

  template<typename T>
  crpropa::Vector3d CartToCyl(const T v, const double cosPhi, const double sinPhi)
  {
    return crpropa::Vector3d(v[0] * cosPhi + v[1] * sinPhi,
                             -v[0] * sinPhi + v[1] * cosPhi,
                             v[2]);
  }

  // logistic sigmoid function
  inline
  double
  Sigmoid(const double x, const double x0, const double w)
  {
    return 1 / (1 + exp(-(x-x0)/w));
  }

  // angle between v0 = (cos(phi0), sin(phi0)) and v1 = (cos(phi1), sin(phi1))
  inline
  double
  DeltaPhi(const double phi0, const double phi1)
  {
    return acos(cos(phi1)*cos(phi0) + sin(phi1)*sin(phi0));
  }

  // Interal units used in this code.
  // Convert to crpropa with e.g. uf23::kpc / crpropa::kpc.
  // The conversion is, however, only needed in the single non-private
  // getField() method, all other functions use uf23 units.
  const double kPi = 3.1415926535897932384626;
  const double kTwoPi = 2*kPi;
  const double degree = kPi/180.;
  const double kpc = 1;
  const double microgauss = 1;
  const double megayear = 1;
  const double Gpc = 1e6*kpc;
  const double pc = 1e-3*kpc;
  const double second  = megayear / (1e6*60*60*24*365.25);
  const double kilometer = kpc / 3.0856775807e+16;
}

namespace crpropa {
UF23Field::UF23Field(const ModelType mt) :
  fModelType(mt),
  fMaxRadiusSquared(pow(30*uf23::kpc, 2))
{

  // all but expX model have a-->\infty, Eq.(38)
  fPoloidalA    =  1 * Gpc;

  switch (fModelType) {
    // ---------------------------------------------
  case base: {
    fDiskB1        =  1.0878565e+00 * uf23::microgauss;
    fDiskB2        =  2.6605034e+00 * uf23::microgauss;
    fDiskB3        =  3.1166311e+00 * uf23::microgauss;
    fDiskH         =  7.9408965e-01 * uf23::kpc;
    fDiskPhase1    =  2.6316589e+02 * uf23::degree;
    fDiskPhase2    =  9.7782269e+01 * uf23::degree;
    fDiskPhase3    =  3.5112281e+01 * uf23::degree;
    fDiskPitch     =  1.0106900e+01 * uf23::degree;
    fDiskW         =  1.0720909e-01 * uf23::kpc;
    fPoloidalB     =  9.7775487e-01 * uf23::microgauss;
    fPoloidalP     =  1.4266186e+00 * uf23::kpc;
    fPoloidalR     =  7.2925417e+00 * uf23::kpc;
    fPoloidalW     =  1.1188158e-01 * uf23::kpc;
    fPoloidalZ     =  4.4597373e+00 * uf23::kpc;
    fStriation     =  3.4557571e-01;
    fToroidalBN    =  3.2556760e+00 * uf23::microgauss;
    fToroidalBS    = -3.0914569e+00 * uf23::microgauss;
    fToroidalR     =  1.0193815e+01 * uf23::kpc;
    fToroidalW     =  1.6936993e+00 * uf23::kpc;
    fToroidalZ     =  4.0242749e+00 * uf23::kpc;
    break;
  }
  case cre10: {
    // ---------------------------------------------
    fDiskB1        =  1.2035697e+00 * uf23::microgauss;
    fDiskB2        =  2.7478490e+00 * uf23::microgauss;
    fDiskB3        =  3.2104342e+00 * uf23::microgauss;
    fDiskH         =  8.0844932e-01 * uf23::kpc;
    fDiskPhase1    =  2.6515882e+02 * uf23::degree;
    fDiskPhase2    =  9.8211313e+01 * uf23::degree;
    fDiskPhase3    =  3.5944588e+01 * uf23::degree;
    fDiskPitch     =  1.0162759e+01 * uf23::degree;
    fDiskW         =  1.0824003e-01 * uf23::kpc;
    fPoloidalB     =  9.6938453e-01 * uf23::microgauss;
    fPoloidalP     =  1.4150957e+00 * uf23::kpc;
    fPoloidalR     =  7.2987296e+00 * uf23::kpc;
    fPoloidalW     =  1.0923051e-01 * uf23::kpc;
    fPoloidalZ     =  4.5748332e+00 * uf23::kpc;
    fStriation     =  2.4950386e-01;
    fToroidalBN    =  3.7308133e+00 * uf23::microgauss;
    fToroidalBS    = -3.5039958e+00 * uf23::microgauss;
    fToroidalR     =  1.0407507e+01 * uf23::kpc;
    fToroidalW     =  1.7398375e+00 * uf23::kpc;
    fToroidalZ     =  2.9272800e+00 * uf23::kpc;
    break;
  }
  case nebCor: {
    // ---------------------------------------------
    fDiskB1        =  1.4081935e+00 * uf23::microgauss;
    fDiskB2        =  3.5292400e+00 * uf23::microgauss;
    fDiskB3        =  4.1290147e+00 * uf23::microgauss;
    fDiskH         =  8.1151971e-01 * uf23::kpc;
    fDiskPhase1    =  2.6447529e+02 * uf23::degree;
    fDiskPhase2    =  9.7572660e+01 * uf23::degree;
    fDiskPhase3    =  3.6403798e+01 * uf23::degree;
    fDiskPitch     =  1.0151183e+01 * uf23::degree;
    fDiskW         =  1.1863734e-01 * uf23::kpc;
    fPoloidalB     =  1.3485916e+00 * uf23::microgauss;
    fPoloidalP     =  1.3414395e+00 * uf23::kpc;
    fPoloidalR     =  7.2473841e+00 * uf23::kpc;
    fPoloidalW     =  1.4318227e-01 * uf23::kpc;
    fPoloidalZ     =  4.8242603e+00 * uf23::kpc;
    fStriation     =  3.8610837e-10;
    fToroidalBN    =  4.6491142e+00 * uf23::microgauss;
    fToroidalBS    = -4.5006610e+00 * uf23::microgauss;
    fToroidalR     =  1.0205288e+01 * uf23::kpc;
    fToroidalW     =  1.7004868e+00 * uf23::kpc;
    fToroidalZ     =  3.5557767e+00 * uf23::kpc;
    break;
  }
  case neCL: {
    // ---------------------------------------------
    fDiskB1        =  1.4259645e+00 * uf23::microgauss;
    fDiskB2        =  1.3543223e+00 * uf23::microgauss;
    fDiskB3        =  3.4390669e+00 * uf23::microgauss;
    fDiskH         =  6.7405199e-01 * uf23::kpc;
    fDiskPhase1    =  1.9961898e+02 * uf23::degree;
    fDiskPhase2    =  1.3541461e+02 * uf23::degree;
    fDiskPhase3    =  6.4909767e+01 * uf23::degree;
    fDiskPitch     =  1.1867859e+01 * uf23::degree;
    fDiskW         =  6.1162799e-02 * uf23::kpc;
    fPoloidalB     =  9.8387831e-01 * uf23::microgauss;
    fPoloidalP     =  1.6773615e+00 * uf23::kpc;
    fPoloidalR     =  7.4084361e+00 * uf23::kpc;
    fPoloidalW     =  1.4168192e-01 * uf23::kpc;
    fPoloidalZ     =  3.6521188e+00 * uf23::kpc;
    fStriation     =  3.3600213e-01;
    fToroidalBN    =  2.6256593e+00 * uf23::microgauss;
    fToroidalBS    = -2.5699466e+00 * uf23::microgauss;
    fToroidalR     =  1.0134257e+01 * uf23::kpc;
    fToroidalW     =  1.1547728e+00 * uf23::kpc;
    fToroidalZ     =  4.5585463e+00 * uf23::kpc;
    break;
  }
  case spur: {
    // ---------------------------------------------
    fDiskB1        = -4.2993328e+00 * uf23::microgauss;
    fDiskH         =  7.5019749e-01 * uf23::kpc;
    fDiskPhase1    =  1.5589875e+02 * uf23::degree;
    fDiskPitch     =  1.2074432e+01 * uf23::degree;
    fDiskW         =  1.2263120e-01 * uf23::kpc;
    fPoloidalB     =  9.9302987e-01 * uf23::microgauss;
    fPoloidalP     =  1.3982374e+00 * uf23::kpc;
    fPoloidalR     =  7.1973387e+00 * uf23::kpc;
    fPoloidalW     =  1.2262244e-01 * uf23::kpc;
    fPoloidalZ     =  4.4853270e+00 * uf23::kpc;
    fSpurCenter    =  1.5718686e+02 * uf23::degree;
    fSpurLength    =  3.1839577e+01 * uf23::degree;
    fSpurWidth     =  1.0318114e+01 * uf23::degree;
    fStriation     =  3.3022369e-01;
    fToroidalBN    =  2.9286724e+00 * uf23::microgauss;
    fToroidalBS    = -2.5979895e+00 * uf23::microgauss;
    fToroidalR     =  9.7536425e+00 * uf23::kpc;
    fToroidalW     =  1.4210055e+00 * uf23::kpc;
    fToroidalZ     =  6.0941229e+00 * uf23::kpc;
    break;
  }
  case synCG: {
    // ---------------------------------------------
    fDiskB1        =  8.1386878e-01 * uf23::microgauss;
    fDiskB2        =  2.0586930e+00 * uf23::microgauss;
    fDiskB3        =  2.9437335e+00 * uf23::microgauss;
    fDiskH         =  6.2172353e-01 * uf23::kpc;
    fDiskPhase1    =  2.2988551e+02 * uf23::degree;
    fDiskPhase2    =  9.7388282e+01 * uf23::degree;
    fDiskPhase3    =  3.2927367e+01 * uf23::degree;
    fDiskPitch     =  9.9034844e+00 * uf23::degree;
    fDiskW         =  6.6517521e-02 * uf23::kpc;
    fPoloidalB     =  8.0883734e-01 * uf23::microgauss;
    fPoloidalP     =  1.5820957e+00 * uf23::kpc;
    fPoloidalR     =  7.4625235e+00 * uf23::kpc;
    fPoloidalW     =  1.5003765e-01 * uf23::kpc;
    fPoloidalZ     =  3.5338550e+00 * uf23::kpc;
    fStriation     =  6.3434763e-01;
    fToroidalBN    =  2.3991193e+00 * uf23::microgauss;
    fToroidalBS    = -2.0919944e+00 * uf23::microgauss;
    fToroidalR     =  9.4227834e+00 * uf23::kpc;
    fToroidalW     =  9.1608418e-01 * uf23::kpc;
    fToroidalZ     =  5.5844594e+00 * uf23::kpc;
    break;
  }
  case twistX: {
    // ---------------------------------------------
    fDiskB1        =  1.3741995e+00 * uf23::microgauss;
    fDiskB2        =  2.0089881e+00 * uf23::microgauss;
    fDiskB3        =  1.5212463e+00 * uf23::microgauss;
    fDiskH         =  9.3806180e-01 * uf23::kpc;
    fDiskPhase1    =  2.3560316e+02 * uf23::degree;
    fDiskPhase2    =  1.0189856e+02 * uf23::degree;
    fDiskPhase3    =  5.6187572e+01 * uf23::degree;
    fDiskPitch     =  1.2100979e+01 * uf23::degree;
    fDiskW         =  1.4933338e-01 * uf23::kpc;
    fPoloidalB     =  6.2793114e-01 * uf23::microgauss;
    fPoloidalP     =  2.3292519e+00 * uf23::kpc;
    fPoloidalR     =  7.9212358e+00 * uf23::kpc;
    fPoloidalW     =  2.9056201e-01 * uf23::kpc;
    fPoloidalZ     =  2.6274437e+00 * uf23::kpc;
    fStriation     =  7.7616317e-01;
    fTwistingTime  =  5.4733549e+01 * uf23::megayear;
    break;
  }
  case expX: {
    // ---------------------------------------------
    fDiskB1        =  9.9258148e-01 * uf23::microgauss;
    fDiskB2        =  2.1821124e+00 * uf23::microgauss;
    fDiskB3        =  3.1197345e+00 * uf23::microgauss;
    fDiskH         =  7.1508681e-01 * uf23::kpc;
    fDiskPhase1    =  2.4745741e+02 * uf23::degree;
    fDiskPhase2    =  9.8578879e+01 * uf23::degree;
    fDiskPhase3    =  3.4884485e+01 * uf23::degree;
    fDiskPitch     =  1.0027070e+01 * uf23::degree;
    fDiskW         =  9.8524736e-02 * uf23::kpc;
    fPoloidalA     =  6.1938701e+00 * uf23::kpc;
    fPoloidalB     =  5.8357990e+00 * uf23::microgauss;
    fPoloidalP     =  1.9510779e+00 * uf23::kpc;
    fPoloidalR     =  2.4994376e+00 * uf23::kpc;
    // internally, xi is fitted and z = tan(xi)*a
    fPoloidalXi    =  2.0926122e+01 * uf23::degree;
    fPoloidalZ     =  fPoloidalA*tan(fPoloidalXi);
    fStriation     =  5.1440500e-01;
    fToroidalBN    =  2.7077434e+00 * uf23::microgauss;
    fToroidalBS    = -2.5677104e+00 * uf23::microgauss;
    fToroidalR     =  1.0134022e+01 * uf23::kpc;
    fToroidalW     =  2.0956159e+00 * uf23::kpc;
    fToroidalZ     =  5.4564991e+00 * uf23::kpc;
    break;
  }
  default: {
    throw std::runtime_error("unknown field model");
    break;
  }
  }

  fSinPitch = sin(fDiskPitch);
  fCosPitch = cos(fDiskPitch);
  fTanPitch = tan(fDiskPitch);

}

Vector3d
UF23Field::getField(const Vector3d &position)
  const
{
  Vector3d posInKpc = position / kpc;
  return (this->operator()(posInKpc)) / uf23::microgauss * microgauss;
}


Vector3d
UF23Field::operator()(const Vector3d& posInKpc)
  const
{
  const auto pos = posInKpc * uf23::kpc;
  if (pos.getR2() > fMaxRadiusSquared)
    return Vector3d(0, 0, 0);
  else {
    const auto diskField = getDiskField(pos);
    const auto haloField = getHaloField(pos);
    return (diskField + haloField) / uf23::microgauss;
  }
}

Vector3d
UF23Field::getDiskField(const Vector3d& pos)
  const
{
  if (fModelType == spur)
    return getSpurField(pos.x, pos.y, pos.z);
  else
    return getSpiralField(pos.x, pos.y, pos.z);
}


Vector3d
UF23Field::getHaloField(const Vector3d& pos)
  const
{
  if (fModelType == twistX)
    return getTwistedHaloField(pos.x, pos.y, pos.z);
  else
    return
      getToroidalHaloField(pos.x, pos.y, pos.z) +
      getPoloidalHaloField(pos.x, pos.y, pos.z);
}


Vector3d
UF23Field::getTwistedHaloField(const double x, const double y, const double z)
  const
{
  const double r = sqrt(x*x + y*y);
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;

  const Vector3d bXCart = getPoloidalHaloField(x, y, z);
  const double bXCartTmp[3] = {bXCart.x, bXCart.y, bXCart.z};
  const Vector3d bXCyl = uf23::CartToCyl(bXCartTmp, cosPhi, sinPhi);

  const double bZ = bXCyl.z;
  const double bR = bXCyl.x;

  double bPhi = 0;

  if (fTwistingTime != 0 && r != 0) {
    // radial rotation curve parameters (fit to Reid et al 2014)
    const double v0 = -240 * uf23::kilometer/uf23::second;
    const double r0 = 1.6 * uf23::kpc;
    // vertical gradient (Levine+08)
    const double z0 = 10 * uf23::kpc;

    // Eq.(43)
    const double fr = 1 - exp(-r/r0);
    // Eq.(44)
    const double t0 = exp(2*std::abs(z)/z0);
    const double gz = 2 / (1 + t0);

    // Eq. (46)
    const double signZ = z < 0 ? -1 : 1;
    const double deltaZ =  -signZ * v0 * fr / z0  * t0 * pow(gz, 2);
    // Eq. (47)
    const double deltaR = v0 * ((1-fr)/r0 - fr/r) * gz;

    // Eq.(45)
    bPhi = (bZ * deltaZ + bR * deltaR) * fTwistingTime;

  }
  const double bCylX[3] = {bR, bPhi , bZ};
  return uf23::CylToCart(bCylX, cosPhi, sinPhi);
}

Vector3d
UF23Field::getToroidalHaloField(const double x, const double y, const double z)
  const
{
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  const double absZ = std::abs(z);

  const double b0 = z >= 0 ? fToroidalBN : fToroidalBS;
  const double rh = fToroidalR;
  const double z0 = fToroidalZ;
  const double fwh = fToroidalW;
  const double sigmoidR = uf23::Sigmoid(r, rh, fwh);
  const double sigmoidZ = uf23::Sigmoid(absZ, fDiskH, fDiskW);

  // Eq. (21)
  const double bPhi = b0 * (1. - sigmoidR) * sigmoidZ * exp(-absZ/z0);

  const double bCyl[3] = {0, bPhi, 0};
  const double cosPhi = r > std::numeric_limits<double>::min() ? x / r : 1;
  const double sinPhi = r > std::numeric_limits<double>::min() ? y / r : 0;
  return uf23::CylToCart(bCyl, cosPhi, sinPhi);
}

Vector3d
UF23Field::getPoloidalHaloField(const double x, const double y, const double z)
  const
{
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);

  const double c = pow(fPoloidalA/fPoloidalZ, fPoloidalP);
  const double a0p = pow(fPoloidalA, fPoloidalP);
  const double rp = pow(r, fPoloidalP);
  const double abszp = pow(std::abs(z), fPoloidalP);
  const double cabszp = c*abszp;

  /*
    since $\sqrt{a^2 + b} - a$ is numerical unstable for $b\ll a$,
    we use $(\sqrt{a^2 + b} - a) \frac{\sqrt{a^2 + b} + a}{\sqrt{a^2
    + b} + a} = \frac{b}{\sqrt{a^2 + b} + a}$}
  */

  const double t0 = a0p + cabszp - rp;
  const double t1 = sqrt(pow(t0, 2) + 4*a0p*rp);
  const double ap = 2*a0p*rp / (t1  + t0);

  double a = 0;
  if (ap < 0) {
    if (r > std::numeric_limits<double>::min()) {
      // this should never happen
      throw std::runtime_error("ap = " + std::to_string(ap));
    }
    else
      a = 0;
  }
  else
    a = pow(ap, 1/fPoloidalP);

  // Eq.(29) and Eq.(32)
  const double radialDependence =
    fModelType == expX ?
    exp(-a/fPoloidalR) :
    1 - uf23::Sigmoid(a, fPoloidalR, fPoloidalW);

  // Eq.(28)
  const double Bzz = fPoloidalB * radialDependence;

  // (r/a)
  const double rOverA =  1 / pow(2*a0p / (t1  + t0), 1/fPoloidalP);

  // Eq.(35) for p=n
  const double signZ = z < 0 ? -1 : 1;
  const double Br =
    Bzz * c * a / rOverA * signZ * pow(std::abs(z), fPoloidalP - 1) / t1;

  // Eq.(36) for p=n
  const double Bz = Bzz * pow(rOverA, fPoloidalP-2) * (ap + a0p) / t1;

  if (r < std::numeric_limits<double>::min())
    return Vector3d(0, 0, Bz);
  else {
    const double bCylX[3] = {Br, 0 , Bz};
    const double cosPhi =  x / r;
    const double sinPhi =  y / r;
    return uf23::CylToCart(bCylX, cosPhi, sinPhi);
  }
}

Vector3d
UF23Field::getSpurField(const double x, const double y, const double z)
  const
{
  // reference approximately at solar radius
  const double rRef = 8.2*uf23::kpc;

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  const double r = sqrt(r2);
  if (r < std::numeric_limits<double>::min())
    return Vector3d(0, 0, 0);

  double phi = atan2(y, x);
  if (phi < 0)
    phi += uf23::kTwoPi;

  const double phiRef = fDiskPhase1;
  int iBest = -2;
  double bestDist = -1;
  for (int i = -1; i <= 1; ++i) {
    const double pphi = phi - phiRef + i*uf23::kTwoPi;
    const double rr = rRef*exp(pphi * fTanPitch);
    if (bestDist < 0 || std::abs(r-rr) < bestDist) {
      bestDist =  std::abs(r-rr);
      iBest = i;
    }
  }
  if (iBest == 0) {
    const double phi0 = phi - log(r/rRef) / fTanPitch;

    // Eq. (16)
    const double deltaPhi0 = uf23::DeltaPhi(phiRef, phi0);
    const double delta = deltaPhi0 / fSpurWidth;
    const double B = fDiskB1 * exp(-0.5*pow(delta, 2));

    // Eq. (18)
    const double wS = 5*uf23::degree;
    const double phiC = fSpurCenter;
    const double deltaPhiC = uf23::DeltaPhi(phiC, phi);
    const double lC = fSpurLength;
    const double gS = 1 - uf23::Sigmoid(std::abs(deltaPhiC), lC, wS);

    // Eq. (13)
    const double hd = 1 - uf23::Sigmoid(std::abs(z), fDiskH, fDiskW);

    // Eq. (17)
    const double bS = rRef/r * B * hd * gS;
    const double bCyl[3] = {bS * fSinPitch, bS * fCosPitch, 0};
    const double cosPhi = x / r;
    const double sinPhi = y / r;
    return uf23::CylToCart(bCyl, cosPhi, sinPhi);
  }
  else
    return Vector3d(0, 0, 0);

}

Vector3d
UF23Field::getSpiralField(const double x, const double y, const double z)
  const
{
  // reference radius
  const double rRef = 5*uf23::kpc;
  // inner boundary of spiral field
  const double rInner = 5*uf23::kpc;
  const double wInner = 0.5*uf23::kpc;
  // outer boundary of spiral field
  const double rOuter = 20*uf23::kpc;
  const double wOuter = 0.5*uf23::kpc;

  // cylindrical coordinates
  const double r2 = x*x + y*y;
  if (r2 == 0)
    return Vector3d(0, 0, 0);
  const double r = sqrt(r2);
  const double phi = atan2(y, x);

  // Eq.(13)
  const double hdz = 1 - uf23::Sigmoid(std::abs(z), fDiskH, fDiskW);

  // Eq.(14) times rRef divided by r
  const double rFacI = uf23::Sigmoid(r, rInner, wInner);
  const double rFacO = 1 - uf23::Sigmoid(r, rOuter, wOuter);
  // (using lim r--> 0 (1-exp(-r^2))/r --> r - r^3/2 + ...)
  const double rFac =  r > 1e-5*uf23::pc ? (1-exp(-r*r)) / r : r * (1 - r2/2);
  const double gdrTimesRrefByR = rRef * rFac * rFacO * rFacI;

  // Eq. (12)
  const double phi0 = phi - log(r/rRef) / fTanPitch;

  // Eq. (10)
  const double b =
    fDiskB1 * cos(1 * (phi0 - fDiskPhase1)) +
    fDiskB2 * cos(2 * (phi0 - fDiskPhase2)) +
    fDiskB3 * cos(3 * (phi0 - fDiskPhase3));

  // Eq. (11)
  const double fac = hdz * gdrTimesRrefByR;
  const double bCyl[3] =
    { b * fac * fSinPitch,
      b * fac * fCosPitch,
      0};

  const double cosPhi = x / r;
  const double sinPhi = y / r;
  return uf23::CylToCart(bCyl, cosPhi, sinPhi);
}
}
