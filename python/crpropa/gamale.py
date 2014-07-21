"""
Galactic Magnetic Lens class

A Python reimplementation of the magnetic field lens technique from PARSEC.
PARSEC: A Parametrized Simulation Engine for Ultra-High Energy Cosmic Ray Protons
arXiv:1302.3761
http://web.physik.rwth-aachen.de/Auger_MagneticFields/PARSEC/
"""

from struct import pack, unpack
from bisect import bisect_left
from os import path

import numpy as np
import healpy
from scipy import sparse


def maxColumnSum(M):
    """
    Return the 1-norm (maximum of absolute sums of columns) of the given matrix.
    The absolute value can be omitted, since the matrix elements are all positive.
    """
    return M.sum(axis=0).max()

def maxRowSum(M):
    """
    Return the infinity-norm (maximum of sums of absolute rows) of the given matrix.
    The absolute value can be omitted, since the matrix elements are all positive.
    """
    return M.sum(axis=1).max()

def normalizeRowSum(Mcsc):
    """
    Normalize each row of a CSC matrix to a row sum of 1.
    """
    rowSum = np.array(Mcsc.sum(axis=1).transpose())[0]
    Mcsc.data /= rowSum[Mcsc.indices]

def generateLensPart(fname, nside=64):
    """
    Generate a lens part from the given CRPropa 3 simulation.
    """
    f = np.genfromtxt(fname, names=True)
    row = healpy.vec2pix(nside, f['P0x'], f['P0y'], f['P0z']) # earth
    col = healpy.vec2pix(nside, f['Px'],  f['Py'],  f['Pz'] ) # galaxy
    npix = healpy.nside2npix(nside)
    data = np.ones(len(row))
    M = sparse.coo_matrix((data, (row, col)), shape=(npix, npix))
    M = M.tocsc()
    normalizeRowSum(M)
    return M

def saveLensPart(Mcsc, fname):
    """
    Save the lens part in coordinate type sparse format.
    This format is compatible with PARSEC.
    """
    M = Mcsc.tocoo()
    fout = open(fname, 'wb')
    fout.write(pack('i4', M.nnz))
    fout.write(pack('i4', M.shape[0]))
    fout.write(pack('i4', M.shape[1]))
    data = np.zeros((M.nnz,), dtype=np.dtype([('row', 'i4'), ('col','i4'),('data','f8')]))
    data['row'] = M.row
    data['col'] = M.col
    data['data'] = M.data
    data.tofile(fout)
    fout.close()

def loadLensPart(fname):
    """
    Load a lens part from the given PARSEC file.
    """
    fin = open(fname, 'rb')
    nnz = unpack('i', fin.read(4))[0]
    nrows = unpack('i', fin.read(4))[0]
    ncols = unpack('i', fin.read(4))[0]
    data = np.fromfile(fin, dtype=np.dtype([('row','i4'), ('col','i4'), ('data','f8')]))
    fin.close()
    M = sparse.coo_matrix((data['data'],(data['row'], data['col'])), shape=(nrows, ncols))
    return M.tocsc()

def randVecInPix(nside, ipix, nest=False):
    """
    Draw vectors from a uniform distribution within a HEALpixel.
    nside : healpix nside parameter
    ipix  : pixel number(s)
    """
    if not(nest):
        ipix = healpy.ring2nest(nside, ipix=ipix)

    norder = nside2norder(nside)
    nUp = 29 - norder
    iUp = ipix * 4**nUp

    if np.iterable(ipix):
        iUp += np.random.randint(0, 4**nUp, size=np.size(ipix))
    else:
        iUp += np.random.randint(0, 4**nUp)
    vec = healpy.pix2vec(nside=2**29, ipix=iUp, nest=True)
    return vec

class Lens:
    """
    Galactic magnetic field lens class with the following conventions:
     - the lens maps directions at the galactic border (pointing outwards back to the source) to observed directions on Earth (pointing outwards)
     - the Galactic coordinate system is used
     - spherical coordinates are avoided
     - for each logarithmic energy bin there is a lens part represented by a matrix
     - energies are given in EeV
     - the matrices (lensParts) are in compressed sparse column format (scipy.sparse.csc)
     - for each matrix M_ij
        - the row number i indexes the observed direction
        - the column number j the direction at the Galactic edge
     - indices are HEALpixel in ring scheme.
    """
    lensParts = [] # list of matrices in order of ascending energy
    lRmins = [] # lower rigidity bounds per lens (log10(E/Z/[eV]))
    lRmax = 0 # upper rigidity bound of last lens (log10(E/Z/[eV]))
    nside = None # HEALpix nside parameter
    neutralLensPart = None # matrix for neutral particles
    maxColumnSum = 1 # maximum of column sums of all matrices

    def __init__(self, cfname=None):
        """
        If a configuration file is given, the lens will be loaded and normalized.
        Otherwise an empty lens is created.
        """
        if cfname == None:
            pass
        else:
            self.load(cfname)
            self.updateMaxColumnSum()

    def load(self, cfname):
        """
        Load and configure the lens from a config file
        filename minR maxR ... in order of ascending rigidity
        """
        dirname = path.dirname(cfname)
        data = np.genfromtxt(cfname, dtype=[('fname','S1000'),('E0','f'),('E1','f')])
        for fname, lR0, lR1 in data:
            M = loadLensPart(path.join(dirname, fname))
            self.checkLensPart(M)
            self.lensParts.append(M)
            self.lRmins.append(lR0)
            self.lRmax = max(self.lRmax, lR1)
        self.neutralLensPart = sparse.identity(healpy.nside2npix(self.nside), format='csc')

    def checkLensPart(self, M):
        """
        Perform sanity checks and set HEALpix nside parameter.
        """
        nrows, ncols = M.get_shape()
        if nrows != ncols:
            raise Exception("Matrix not square %i x %i"%(nrows, ncols))
        nside = healpy.npix2nside(nrows)
        if self.nside == None:
            self.nside = nside
        elif self.nside != int(nside):
            raise Exception("Matrix have different HEALpix schemes")

    def updateMaxColumnSum(self):
        """
        Update the maximum column sum
        """
        m = 0
        for M in self.lensParts:
            m = max(m, maxColumnSum(M))
        self.maxColumnSum = m

    def getLensPart(self, E, Z=1):
        """
        Return the matrix corresponding to a given energy E [EeV] and charge number Z
        """
        if Z == 0:
            return self.neutralLensPart
        if len(self.lensParts) == 0:
            raise Exception("Lens empty")        
        lR = np.log10(E / Z) + 18
        if (lR < self.lRmins[0]) or (lR > self.lRmax):
            raise ValueError("Rigidity %f/%i EeV not covered"%(E, Z))
        i = bisect_left(self.lRmins, lR) - 1
        return self.lensParts[i]

    def transformPix(self, j, E, Z=1):
        """
        Attempt to transform a pixel (ring scheme), given an energy E [EeV] and charge number Z.
        Returns a pixel (ring scheme) if successful or None if not.
        """
        M = self.getLensPart(E, Z)
        cmp_val = np.random.rand() * self.maxColumnSum
        sum_val = 0
        for i in range(M.indptr[j], M.indptr[j+1]):
            sum_val += M.data[i]
            if cmp_val < sum_val:
                return M.indices[i]
        return None

    def transformVec(self, x, y, z, E, Z=1):
        """
        Attempt to transform a galactic direction, given an energy E [EeV] and charge number Z.
        Returns a triple (x,y,z) if successful or None if not.
        """
        j = healpy.vec2pix(self.nside, x, y, z)
        i = self.transformPix(j, E, Z)
        if i == None:
            return None
        v = healpytools.randVecInPix(self.nside, i)
        return v

    def extragalacticVector(self, i, E, Z=1):
        """
        Return the HEALpix vector of extragalactic directions 
        for a given observed pixel i, energy E [EeV] and charge number Z.
        """
        M = self.getLensPart(E, Z)
        row = M.getrow(i)
        return np.array( row.todense() )[0]

    def observedVector(self, j, E, Z=1):
        """
        Return the HEALpix vector of observed directions 
        for a given extragalactic pixel j, energy E [EeV] and charge number Z.
        """
        M = self.getLensPart(E, Z)
        col = M.getcol(j)
        return np.array( col.transpose().todense() )[0]