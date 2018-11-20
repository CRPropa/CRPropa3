#Returns TRUE if matplotlib is available
#Silently exits with -1 otherwise
import sys
from distutils.version import LooseVersion
try:
  import matplotlib
  if LooseVersion(matplotlib.__version__) < LooseVersion("2.0.0"):
    sys.stdout.write('Need matplotlib version >= 2.0.0 - found ' + matplotlib.__version__ + "\n")
    sys.exit(-1)

except ImportError:
  sys.exit(-1)
