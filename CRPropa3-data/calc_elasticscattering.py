import numpy as np
import os
import photonField
import interactionRate as iR


eV = 1.60217657e-19

# output folder
folder = 'data/ElasticScattering'
if not os.path.exists(folder):
    os.makedirs(folder)

gamma = np.logspace(6, 14, 201)  # tabulated UHECR Lorentz-factors

fields = [
    photonField.CMB(),
    photonField.EBL_Gilmore12()]


# load cross section data from TALYS
ddir = 'tables/PD_Talys1.8_Khan/'
eps = np.genfromtxt(ddir + 'eps_elastic.txt') * eV * 1e6  # nuclear rest frame photon energies [J]
data = np.genfromtxt(ddir + 'xs_elastic.txt', dtype=[('Z', int), ('N', int), ('xs', '%if8' % len(eps))])

# only consider TALYS cross sections for A >= 12
idx = (data['Z'] + data['N']) >= 12
data = data[idx]

# factor out the principal scaling given by the TRK formula: sigma_int ~ Z*N/A
data['xs'] /= (data['Z'] * data['N'] / (data['Z'] + data['N']))[:, np.newaxis]

# pad cross sections to next larger 2^n + 1 tabulation points for Romberg integration
eps = iR.romb_pad_logspaced(eps, 513)
xs = np.array([iR.romb_pad_zero(x, 513) for x in data['xs']]) * 1e-31


for field in fields:
    print(field.name)

    # calculate the interaction rate, averaged over all isotopes
    rate = np.mean([iR.calc_rate_eps(eps, x, gamma, field) for x in xs], axis=0)
    fname = folder + '/rate_%s.txt' % field.name.split('_')[0]
    header = 'Average interaction rate for elastic scattering of %s photons off nuclei\nScale with Z*N/A for nuclei\n1/lambda [1/Mpc] for log10(gamma) = 6-14 in 201 steps' % field.info
    np.savetxt(fname, rate, fmt='%g', header=header)

    # calculate CDF for background photon energies, averaged over all isotopes
    CDF = np.zeros((len(gamma), len(eps)))
    for x in xs:
        C = iR.calc_rate_eps(eps, x, gamma, field, cdf=True)
        CDF += C / np.max(C, axis=1, keepdims=True)
    CDF /= len(data)
    CDF = np.nan_to_num(CDF)

    fname = folder + '/cdf_%s.txt' % field.name.split('_')[0]
    data = np.c_[np.log10(gamma), CDF]
    fmt = '%g' + '\t%g' * len(eps)
    header = '# Average CDF(background photon energy) for elastic scattering with the %s\n# log10(gamma), (1/lambda)_cumulative for eps = log10(2 keV) - log10(263 MeV) in 513 steps' % field.info
    np.savetxt(fname, data, fmt=fmt, header=header)
