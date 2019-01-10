"""
purpose: convert IRB photon field data from external sources to 
a unified format for CRPropa to read. The 7 photon fields that can be converted
by this script are:
    IRB_Kneiske04, IRB_Stecker05, IRB_Finke10, IRB_Dominguez11,
    IRB_Gilmore12, IRB_Stecker16_upper, IRB_Stecker16_lower
usage: this file should be placed in the CRPropa-data/ folder
output: .txt files named by their field name
Written by Mario Hörbe (mario.hoerbe@rub.de)
"""

import numpy as np
import pandas as pd 
import os

eV = 1.60217657e-19  # [J]
c0 = 299792458  # [m/s]
h = 6.62606957e-34  # [m^2 kg / s]


def IRB_Stecker05(fileDir):
    name = 'IRB_Stecker05'
    info = '# cosmic infrared and optical background radiation model of Stecker at al. 2005'
    redshift = np.linspace(0, 5, 26)
    filePath = fileDir + "EBL_Stecker_2005/data2.txt"
    d = np.genfromtxt(filePath, unpack=True)
    eps = 10**d[0] # [eV]
    n = 10**d[1:]  # [1/cm^3]
    n /= eps       # [1/eVcm^3]
    n = pd.DataFrame(n)
    photonField = []
    energy = []
    for col in n.columns:
        energy.append(eps[col])
        photonField.append(n[col].values)
    createField(name, info, energy, redshift, photonField)


def IRB_Gilmore12(fileDir):
    name = "IRB_Gilmore12"
    info = "# These tables contain the data for the background flux and associated optical depths of gamma rays for the WMAP5+Fixed ('fixed') and Evolving Dust ('fiducial') models presented in Gilmore, Somerville, Primack, and Dominguez (2012), ArXiv:1104.0671v2"
    filePath = fileDir + "EBL_Gilmore_2012/eblflux_fiducial.dat"
    redshift = [0.0,0.015,0.025,0.044,0.05,0.2,0.4,0.5,0.6,0.8,1.0,1.25,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0]
    d = np.genfromtxt(filePath, unpack=True)
    wavelength = d[0]  # angstrom
    eps = 12.39842e3 / wavelength  # [eV]
    photonField = []
    energy = []
    d = pd.DataFrame(d)
    for i in range(len(eps)):
        fieldSlice = np.array(list(d[i])[1:]) * 4*np.pi /(100*c0) *1e-10 *6.2415091e11  # eV/cm^3
        fieldSlice /= eps[i]**2  # 1/eVcm^3
        photonField.append(fieldSlice / eV)  # /eV?!
        energy.append(eps[i])
    # invert, because lambda is antiprop to energy
    photonField = [x for x in reversed(photonField)]
    energy = [e for e in reversed(energy)]
    createField(name, info, energy, redshift, photonField)


def IRB_Finke10(fileDir):
    name = "IRB_Finke10"
    redshift = np.round(np.linspace(0,4.99,500), 2)
    info = "# Extragalactic background light model from Finke et al. 2010, DOI:10.1088/0004-637X/712/1/238, Files obtained from http://www.phy.ohiou.edu/~finke/EBL/"
    fileDir = fileDir + "EBL_Finke_2010/"
    fileList = os.listdir(fileDir)
    d = pd.DataFrame()
    col = 0
    for file in sorted(fileList):
        if "README.txt" not in file:
            data = np.genfromtxt(fileDir + file, unpack=True)
            # eps = data[0]
            eps = data[0]  # [eV]
            data = pd.DataFrame(data[1],columns=[col/100])
            col += 1
            d = pd.concat([d,data], axis=1, join_axes=[data.index])
    photonField = []
    energy = []
    for i,e in enumerate(eps):
        dens = np.array(list(d.iloc[i])) * 6.2415091e11 / eps[i]**2  # [1/eVcm^3]
        photonField.append(dens)
        energy.append(eps[i])
    createField(name, info, energy, redshift, photonField)


def IRB_Dominguez11(fileDir):
    name = "IRB_Dominguez11"
    info = "# EBL intensities for the paper >Extragalactic background light inferred from AEGIS galaxy-SED-type fractions<, A. Dominguez et al., 2011, MNRAS, 410, 2556"
    redshift = np.array([0,0.01,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.25,1.5,2.0,2.5,3.0,3.9])
    filePath = fileDir + "EBL_Dominguez_2011/ebl_dominguez11.out"
    d = np.genfromtxt(filePath, unpack=False)
    eps = [1.239842/fieldSlice[0]*eV for fieldSlice in d]  # [eV] | 1.238842 = h*c/1µ eV
    #                          nW->W : J->eV : 1/m³->1/cm³ : 1/sm²sr->1/m³
    n = np.array([fieldSlice[1:] *1e-9 *eV /1e6 /eps[i]**2 * (4*np.pi/c0) for i,fieldSlice in enumerate(d)])  # [1/eVcm^3] 
    energy = []
    for i,x in enumerate(n):
        energy.append(eps[i]/eV)
    photonField = [x for x in reversed(n)]
    energy = [e for e in reversed(energy)]
    createField(name, info, energy, redshift, photonField)


def IRB_Stecker16_lower(fileDir):
    name = "IRB_Stecker16_lower"
    info = "# Extragalactic background light model from Stecker et al. 2016, DOI:10.3847/0004-637X/827/1/6 <An Empirical Determination of the Intergalactic Background Light from UV to FIR Wavelengths Using FIR Deep Galaxy Surveys and the Gamma-ray Opacity of the Universe>"
    redshift = np.linspace(0, 5, 501)
    filePath = fileDir + "EBL_Stecker_2016/comoving_enerdens_lo.csv"
    d = np.genfromtxt(filePath, delimiter=',')
    eps = 10**np.arange(-2.84, 1.14001, 0.01)  # [eV]
    nu = eps * eV / h  # [Hz]
    photonField = []
    energy = []
    for i,dens in enumerate(d):
        d[i] = d[i] * 6.2415091e11 * nu[i] / eps[i]**2  # 1/eVcm^3
        photonField.append(d[i])
        energy.append(eps[i])
    createField(name, info, energy, redshift, photonField)


def IRB_Stecker16_upper(fileDir):
    name = "IRB_Stecker16_upper"
    info = "# Extragalactic background light model from Stecker et al. 2016, DOI:10.3847/0004-637X/827/1/6 <An Empirical Determination of the Intergalactic Background Light from UV to FIR Wavelengths Using FIR Deep Galaxy Surveys and the Gamma-ray Opacity of the Universe>"
    redshift = np.linspace(0, 5, 501)
    filePath = fileDir + "EBL_Stecker_2016/comoving_enerdens_up.csv"
    d = np.genfromtxt(filePath, delimiter=',')
    eps = 10**np.arange(-2.84, 1.14001, 0.01)  # [eV]
    nu = eps * eV / h  # [Hz]
    photonField = []
    energy = []
    for i,dens in enumerate(d):
        d[i] = d[i] * 6.2415091e11 * nu[i] / eps[i]**2  # 1/eVcm^3
        photonField.append(d[i])
        energy.append(eps[i])
    createField(name, info, energy, redshift, photonField)


def IRB_Kneiske04(fileDir):
    name = "IRB_Kneiske04"
    info = "# IRO spectrum from Tanja Kneiske et al. obtained from O. E. Kalashev (http://hecr.inr.ac.ru/)"
    redshift = np.linspace(0,5,51)
    filePath = fileDir + "EBL_Kneiske_2004/all_z"
    d = np.genfromtxt(filePath)
    photonField = []
    energy = []
    for i,fieldSlice in enumerate(d):
        for j,entry in enumerate(fieldSlice):
            if j != 0:
                d[i][j] /= 1e6  # 1/eVm^3 -> 1/eVcm^3
        photonField.append(fieldSlice[1:])
        energy.append(fieldSlice[0])
    createField(name, info, energy, redshift, photonField)


def createField(name, info, energy, redshift, photonField):
    with open(name+".txt", 'w') as f:
        f.write(info+"\n")
        f.write("# energy / eV:\n")
        np.savetxt(f, [energy], fmt="%1.6e", delimiter=" ")
        f.write("# redshift:\n")
        np.savetxt(f, [redshift], fmt="%1.2f", delimiter=" ")
        f.write("# photon field density / 1/eVcm^3:\n")
        np.savetxt(f, photonField, fmt="%1.6e", delimiter=" ")
    print("done: " + name)


if __name__ == "__main__":

    fileDir = "/tables/"
    IRB_Kneiske04(fileDir)
    IRB_Stecker05(fileDir)
    IRB_Finke10(fileDir)
    IRB_Dominguez11(fileDir)
    IRB_Gilmore12(fileDir)
    IRB_Stecker16_upper(fileDir)
    IRB_Stecker16_lower(fileDir)
