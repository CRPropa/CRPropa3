#### Nuclear Masses
Nuclear masses are stored in data/nuclear_mass.txt
The masses are calculated from the atomic mass table of NIST, http://www.nist.gov/pml/data/comp.cfm

UHECRs are fully ionized so electron masses are substracted from the atomic masses.
<math>
m_{nucleus} = m_{atom} - Z * m_e
</math>

Electron binding energies (~keV) are thus neglected compared to atomic masses (~GeV) electron masses (~MeV) and nuclea binding energies (~MeV)

For isotopes which are not present in atomic mass table the nucleon appromation is used
<math>
m(A,Z) = A * amu - Z * m_e
</math>
