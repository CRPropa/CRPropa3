Photodisintegration cross sections for isotopes with A<12 from various references.
The selection of cross sections is briefly described in [1] and in more detail in [2].


Overview:
- H-2, H-3, He-3, He-4, Be-9 from [3]
    The cross sections for H-3 and He-3 are scaled by 1.7 and 0.66 (cf. [2], page 84; page 125 contains a typo)
    For Be-9 the parametrization is refitted to data (cf. [2] page 85, figure 3.2)
- Li-8, Be-7, Be-10, Be-11, B-8, B-10, B-11, C-10, C-11 from [4]
    The loss of one proton (neutron) is assumed if the neutron number N < Z (N > Z).
    For N = Z the loss of one neutron or proton is modeled with equal probability.
- Li-7 from interpolation of experimental data [6,7]


Changes:
- He-6, Li-8, Li-9 and B-5 which are present in the CRPropa 2.0 release were removed, because the life-time is smaller than 2 seconds, rendering photonuclear interactions negligible
- Li-6 is now also taken from [4] instead of from TALYS, modeling the loss of protons and neutrons with equal probability


Contents:
- isotopes.txt contains the list of isotopes
- <Z>-<A> folders contain the cross section data in TALYS format
- crosssection.cc is the original script by Nils Nierstenhoefer
- collect.py collects and summarizes the cross section data
- eps.txt    tabulated incident photon energies in MeV
- xs.txt     tabulated exclusive cross sections
- xs_sum.txt tabulated total cross sections


References:
[1] Kampert et al. 2013, DOI:10.1016/j.astropartphys.2012.12.001, arXiv:1206.3132
[2] Nils Nierstenhoefer PhD thesis
[3] J. Rachen PhD thesis
[4] Kossov 2002, European Physical Journal A, vol. 14, p. 377-92
[6] V.V. Varlamov et al., Photonuclear data. Photodisintegration of lithium. Evaluated cross sections of channels and reactions, Fotojad. Dannye Photodisint. Li. Suppl., Moscow, 1986.
[7] L.A. Kulchitskii, Y.M. Volkov, V.P. Denisov, V.I. Ogurtsov, Energy Levels of Li-7 observed in its photoemission, Izv. Rossiiskoi Akad. Nauk. Ser. Fiz. 27 (1963) 1412.
