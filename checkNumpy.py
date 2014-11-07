#Returns the nupy include path if nupy is available
#Silently exits with -1 otherwise
import sys
try:
  import numpy
  sys.stdout.write(numpy.get_include())
except ImportError:
  sys.exit(-1)
