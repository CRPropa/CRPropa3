# Returns the path to the crpropa swig_interface, install_prefix, or -1 if
# crpropa is not found
import sys
try:
	import crpropa
except ImportError:
	sys.exit(-1)

if sys.argv[1] == 'swig_interface':
	sys.stdout.write(crpropa.getDataPath('swig_interface'))
elif sys.argv[1] == 'install_prefix':
	sys.stdout.write(crpropa.getInstallPrefix())
