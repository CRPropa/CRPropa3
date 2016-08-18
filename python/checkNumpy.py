#Returns the numpy include path if numpy version > 1.6 is available
#Silently exits with -1 otherwise
import sys
try:
  import numpy
  
  from pkg_resources import parse_version
  if parse_version(numpy.__version__) < parse_version('1.6.0'):
    sys.exit(-1)
  
  sys.stdout.write(numpy.get_include())

except ImportError:
  sys.exit(-1)
