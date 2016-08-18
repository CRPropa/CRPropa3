#Returns TRUE if matplotlib is available
#Silently exits with -1 otherwise
import sys
try:
  import matplotlib
  sys.stdout.write('TRUE')

except ImportError:
  sys.exit(-1)
