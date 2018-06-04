from __future__ import division
import numpy as np
import interactionRate
import photonField
import os


eV = 1.60217657E-19  # [J]
me2 = (510.998918E3 * eV)**2  # squared electron mass [J^2/c^4]
sigmaThompson = 6.6524E-29  # Thompson cross section [m^2]
alpha = 1 / 137.035999074  # fine structure constant


def sigmaPP(s):
    """ Pair production cross section (Bethe-Heitler), see Lee 1996 """
    smin = 4 * me2
    if (s < smin):
        return 0
    else:
        b = np.sqrt(1 - smin / s)
        return sigmaThompson * 3 / 16 * (1 - b**2) * ((3 - b**4) * np.log((1 + b) / (1 - b)) - 2 * b * (2 - b**2))


def sigmaDPP(s):
    """ Double-pair production cross section, see R.W. Brown eq. (4.5) with k^2 = q^2 = 0 """
    smin = 16 * me2
    if (s < smin):
        return 0
    else:
        return 6.45E-34 * (1 - smin / s)**6


def sigmaICS(s):
    """ Inverse Compton scattering cross sections, see Lee 1996 """
    smin = me2
    if (s < smin):  # numerically unstable close to smin
        return 0
    else:
        # note: formula unstable for (s - smin) / smin < 1E-5
        b = (s - smin) / (s + smin)
        A = 2 / b / (1 + b) * (2 + 2 * b - b**2 - 2 * b**3)
        B = (2 - 3 * b**2 - b**3) / b**2 * np.log((1 + b) / (1 - b))
        return sigmaThompson * 3 / 8 * smin / s / b * (A - B)


def sigmaTPP(s):
    """ Triplet-pair production cross section, see Lee 1996 """
    beta = 28 / 9 * np.log(s / me2) - 218 / 27
    if beta < 0:
        return 0
    else:
        return sigmaThompson * 3 / 8 / np.pi * alpha * beta


def getTabulatedXS(sigma, skin):
    """ Get crosssection for tabulated s_kin """
    if sigma in (sigmaPP, sigmaDPP):  # photon interactions
        return np.array([sigma(s) for s in skin])
    if sigma in (sigmaTPP, sigmaICS):  # electron interactions
        return np.array([sigma(s) for s in skin + me2])
    return False


def getSmin(sigma):
    """ Return minimum required s_kin = s - (mc^2)^2 for interaction """
    return {sigmaPP: 4 * me2,
            sigmaDPP: 16 * me2,
            sigmaTPP: np.exp((218 / 27) / (28 / 9)) * me2 - me2,
            sigmaICS: 1E-9 * me2
            }[sigma]


def getEmin(sigma, field):
    """ Return minimum required cosmic ray energy for interaction *sigma* with *field* """
    return getSmin(sigma) / 4 / field.getEmax()


def process(sigma, field, name):
    # output folder
    folder = 'data/' + name
    if not os.path.exists(folder):
        os.makedirs(folder)

    # tabulated energies, limit to energies where the interaction is possible
    Emin = getEmin(sigma, field)
    E = np.logspace(10, 23, 261) * eV
    E = E[E > Emin]

    # -------------------------------------------
    # calculate interaction rates
    # -------------------------------------------
    # tabulated values of s_kin = s - mc^2
    # Note: integration method (Romberg) requires 2^n + 1 log-spaced tabulation points
    s_kin = np.logspace(6, 23, 2049) * eV**2
    xs = getTabulatedXS(sigma, s_kin)
    rate = interactionRate.calc_rate_s(s_kin, xs, E, field)

    # save
    fname = folder + '/rate_%s.txt' % field.name
    data = np.c_[np.log10(E / eV), rate]
    fmt = '%.2f\t%.6g'
    header = '%s interaction rates\nphoton field: %s\nlog10(E/eV), 1/lambda [1/Mpc]' % (name, field.info)
    np.savetxt(fname, data, fmt=fmt, header=header)

    # -------------------------------------------
    # calculate cumulative differential interaction rates for sampling s values
    # -------------------------------------------
    # find minimum value of s_kin
    skin1 = getSmin(sigma)  # s threshold for interaction
    skin2 = 4 * field.getEmin() * E[0]  # minimum achievable s in collision with background photon (at any tabulated E)
    skin_min = max(skin1, skin2)

    # tabulated values of s_kin = s - mc^2, limit to relevant range
    # Note: use higher resolution and then downsample
    skin = np.logspace(6.2, 23, 1680 + 1) * eV**2
    skin = skin[skin > skin_min]

    xs = getTabulatedXS(sigma, skin)
    rate = interactionRate.calc_rate_s(skin, xs, E, field, cdf=True)

    # downsample
    skin_save = np.logspace(6.2, 23, 168 + 1) * eV**2
    skin_save = skin_save[skin_save > skin_min]
    rate_save = np.array([np.interp(skin_save, skin, r) for r in rate])

    # save
    data = np.c_[np.log10(E / eV), rate_save]  # prepend log10(E/eV) as first column
    row0 = np.r_[0, np.log10(skin_save / eV**2)][np.newaxis]
    data = np.r_[row0, data]  # prepend log10(s_kin/eV^2) as first row

    fname = folder + '/cdf_%s.txt' % field.name
    fmt = '%.2f' + '\t%.6g' * np.shape(rate_save)[1]
    header = '%s cumulative differential rate\nphoton field: %s\nlog10(E/eV), d(1/lambda)/ds_kin [1/Mpc/eV^2] for log10(s_kin/eV^2) as given in first row' % (name, field.info)
    np.savetxt(fname, data, fmt=fmt, header=header)


fields = [
    photonField.CMB(),
    photonField.EBL_Kneiske04(),
    photonField.EBL_Stecker05(),
    photonField.EBL_Franceschini08(),
    photonField.EBL_Finke10(),
    photonField.EBL_Dominguez11(),
    photonField.EBL_Gilmore12(),
    photonField.EBL_Stecker16('lower'),
    photonField.EBL_Stecker16('upper'),
    photonField.URB_Protheroe96()]

for field in fields:
    print(field.name)
    process(sigmaPP, field, 'EMPairProduction')
    process(sigmaDPP, field, 'EMDoublePairProduction')
    process(sigmaTPP, field, 'EMTripletPairProduction')
    process(sigmaICS, field, 'EMInverseComptonScattering')
