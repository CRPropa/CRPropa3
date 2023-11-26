# Returns the numpy include path if numpy version > 1.6 is available
# Silently exits with -1 otherwise
import sys

try:
	import numpy

	if sys.version_info.major == 3 and sys.version_info.minor < 7:
		from pkg_resources import parse_version 
		if parse_version(numpy.__version__) < parse_version('1.6.0'):
			sys.exit(-1)

	else:
		from importlib.metadata import distribution
		vi = distribution('numpy').version.split('.')
		if int(vi[0]) == 0:
			sys.exit(-1)
		elif int(vi[0]) == 1:
			if int(vi[1]) < 6:
				sys.exit(-1)

	sys.stdout.write(numpy.get_include())

except ImportError:
	sys.exit(-1)
