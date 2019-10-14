from subprocess import call
from sys import argv
import os

current_dir = os.path.dirname(os.path.abspath(__file__))


def convert_4to3(fname):
    """
    Convert ipython notebook from nbformat=4 to 3
    >>> jupyter-nbconvert --to=notebook --nbformat=3 xxx --output=yyy
    """
    ofname = fname.replace('.v4', '.v3').split('/')[-1]
    call(['jupyter-nbconvert', '--to=notebook', '--nbformat=3', fname, '--output=%s'%ofname])


def convert_topython(fname):
    """
    Convert ipython notebook to python script
    >>> jupyter-nbconvert --to=python xxx --output=yyy
    """
    ofname = fname.replace('.v4.ipynb', '.py').split('/')[-1]
    call(['jupyter-nbconvert', '--to=python', fname, '--output=%s' % ofname])


files = (
    'basics/basics.v4.ipynb',
    'galactic_backtracking/galactic_backtracking.v4.ipynb',
    'galactic_lensing/lensing_cr.v4.ipynb',
    'galactic_lensing/lensing_maps.v4.ipynb',
    'galactic_trajectories/galactic_trajectories.v4.ipynb',
    'secondaries/neutrinos.v4.ipynb',
    'secondaries/photons.v4.ipynb',
    'sim1D/sim1D.v4.ipynb',
    'trajectories/trajectories.v4.ipynb',
    'Diffusion/DiffusionValidationI.v4.ipynb'
    )

if len(argv) == 1:
    # convert all listed notebooks from v4 to v3
    for f in files:
        fullf = os.path.join(current_dir, f)
        convert_4to3(fullf)
        convert_topython(fullf)
else:
    # convert specific notebooks
    for f in argv[1:]:
        fullf = os.path.join(current_dir, f)
        convert_4to3(fullf)
        convert_topython(fullf)
