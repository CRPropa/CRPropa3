# CRPropa test script
# Testing ObserverTracking
#
import matplotlib
matplotlib.use('Agg')

import pylab
import numpy
import sys

def simulate():
	from crpropa import *
	vgrid = VectorGrid(Vector3d(0), 128, 1*Mpc)
	initTurbulence(vgrid, 1*nG, 2*Mpc, 5*Mpc, -11./3.)
	BField = MagneticFieldGrid(vgrid)

	m = ModuleList()
	m.add(PropagationCK(BField, 1e-4, 0.1*Mpc, 5*Mpc))
	m.add(MaximumTrajectoryLength(25 * Mpc))
	out = TextOutput("tracking.txt")
	o = Observer()
	o.add(ObserverTracking(Vector3d(64,64,64) * Mpc, 1 * Mpc, 0.1 * Mpc))
	o.setDeactivateOnDetection(False)
	o.onDetection(out)
	#m.add(out)
	m.add(o)

	# source setup
	source = Source()
	source.add(SourcePosition(Vector3d(70, 64, 64) * Mpc))
	source.add(SourceIsotropicEmission())
	source.add(SourceParticleType(- nucleusId(1,1)))
	source.add(SourceEnergy(1 * EeV))

	m.setShowProgress(True)
	m.run(source, 1000000)


def unique_ids(*args):
	srt = numpy.lexsort(args)
	total_delta = None
	for d in args:
		delta = numpy.diff(d[srt])
		if total_delta is None:
			total_delta = delta
		else:
			total_delta += delta
	total_delta = total_delta != 0
	idx = numpy.cumsum(numpy.hstack(([0], total_delta)))
	a = numpy.zeros(len(srt), dtype=int)
	numpy.put(a, srt, idx)
	return a

def plot_uhecrs(phi, theta, values, cmap='jet', **kwargs):
    f = kwargs.get("figure", pylab.figure(figsize=(10, 10 * 0.75)))
    s = f.add_subplot(111, projection='hammer')
    pylab.grid(lw=1, color='0.5', alpha=0.5)
    # s.invert_xaxis()
    if values is None:
        pylab.scatter(-phi, theta, edgecolors=kwargs.get('edge', 'none'),
                      marker='.', label="a", s=kwargs.get('size', 10))
    else:
        srt = pylab.argsort(values)
        pylab.scatter(-phi[srt], theta[srt], c=values[srt], vmin=min(values),
                      vmax=max(values), edgecolors=kwargs.get('edge', '0.7'),
                      lw=0.1, marker='.', label="a", s=kwargs.get('size', 10), cmap=cmap)
        cb = pylab.colorbar(
            format='%g',
            orientation='horizontal',
            aspect=30,
            shrink=0.8,
            pad=0.1,
            label=kwargs.get('valuelabel', None))
        cb.solids.set_edgecolor("face")
    pylab.gca().set_xticklabels([])
    pylab.tight_layout(pad=0.1)
    return f, s


def plot_uhecrs_xyz(x, y, z, values=None, **kwargs):
    phi = numpy.arctan2(y, x)
    theta = numpy.arctan2(z, (x * x + y * y) ** .5)
    return plot_uhecrs(phi, theta, values=values, **kwargs)


def plot():
	print "Load tracking.txt"
	d = numpy.genfromtxt("tracking.txt", names=True)
	n = len(d)

	# one per particle
	print "Find unique inital particles"
	ids_source = unique_ids(d['ID0'], d['E0'], d['X0'], d['Y0'], d['Z0'], d['P0x'], d['P0y'], d['P0z'], d['ID1'], d['E1'], d['X1'], d['Y1'], d['Z1'], d['P1x'], d['P1y'], d['P1z'])
	nids_source = ids_source.max() + 1
	print "Unique Source: ", nids_source

	# one per initial
	print "Find unique final particles"
	ids_particle = unique_ids(d['ID0'], d['E0'], d['X0'], d['Y0'], d['Z0'], d['P0x'], d['P0y'], d['P0z'])
	nids_particle = ids_particle.max() + 1
	print "Unique Particle: ", nids_source

	# collect data
	print "Detect at different distances"
	select_1000kpc = numpy.ones(nids_source, dtype=int)	* -1
	select_500kpc = numpy.ones(nids_source, dtype=int)	* -1
	select_200kpc = numpy.ones(nids_source, dtype=int)	* -1
	select_100kpc = numpy.ones(nids_source, dtype=int)	* -1
	select_50kpc = numpy.ones(nids_source, dtype=int)	* -1
	select_closest = numpy.ones(nids_source, dtype=int)	* -1

	dist = ((d['X'] - 64)**2 + (d['Y'] - 64)**2 + (d['Z'] - 64)**2)**0.5

	for i in xrange(nids_source):
		if i % 1000 == 0:
			sys.stdout.write(" %5.2f%%\r" % (100. * float(i) / nids_source))
			sys.stdout.flush()
		# select i'th initial particle
		s = numpy.nonzero(ids_source == i)[0]
		s_d = d[s]
		s_dist = dist[s]

		# distances sorted by trajectory
		D_sort = numpy.argsort(s_d['D'])
		dist_sort = s_dist[D_sort]

		idx = numpy.nonzero(dist_sort < 1.000)[0]
		if len(idx) > 0:
			select_1000kpc[i] = s[D_sort[idx[0]]]
		idx = numpy.nonzero(dist_sort < 0.500)[0]
		if len(idx) > 0:
			select_500kpc[i] = s[D_sort[idx[0]]]
		idx = numpy.nonzero(dist_sort < 0.200)[0]
		if len(idx) > 0:
			select_200kpc[i] = s[D_sort[idx[0]]]
		idx = numpy.nonzero(dist_sort < 0.100)[0]
		if len(idx) > 0:
			select_100kpc[i] = s[D_sort[idx[0]]]
		idx = numpy.nonzero(dist_sort < 0.05)[0]
		if len(idx) > 0:
			select_50kpc[i] = s[D_sort[idx[0]]]

		select_closest[i] = s[numpy.argmin(s_dist)]

    # plot first 10 tracks
	first_ten = ids_source < 10
	pylab.figure()
	pylab.scatter(d['X'][first_ten], d['Y'][first_ten], c=ids_source[first_ten], s=50*(d['Z'][first_ten]-63))

	first_ten_1000kpc = select_1000kpc[:10]
	pylab.scatter(d['X'][first_ten_1000kpc], d['Y'][first_ten_1000kpc], marker='+', s=1000)

	first_ten_500kpc = select_500kpc[:10]
	pylab.scatter(d['X'][first_ten_500kpc], d['Y'][first_ten_500kpc], marker='+', s=500)

	first_ten_200kpc = select_200kpc[:10]
	pylab.scatter(d['X'][first_ten_200kpc], d['Y'][first_ten_200kpc], marker='+', s=200)

	first_ten_100kpc = select_100kpc[:10]
	pylab.scatter(d['X'][first_ten_100kpc], d['Y'][first_ten_100kpc], marker='+', s=100)

	first_ten_50kpc = select_50kpc[:10]
	pylab.scatter(d['X'][first_ten_50kpc], d['Y'][first_ten_50kpc], marker='+', s=50)

	first_ten_closest = select_closest[:10]
	pylab.scatter(d['X'][first_ten_closest], d['Y'][first_ten_closest], marker='x', s=50)

	pylab.xlim(63, 65)
	pylab.ylim(63, 65)
	pylab.savefig("scatter.png")
	pylab.show()
	pylab.close()

	# plot skplots
	x, y, z = -d['Px'], -d['Py'], -d['Pz']
	phi = numpy.arctan2(y, x)
	theta = numpy.arctan2(z, (x * x + y * y) ** .5)

	plot_uhecrs(phi, theta, dist, cmap='jet_r')
	pylab.savefig("all.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi, theta, 1./dist)
	pylab.savefig("all_w.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi, theta, 1./(dist**2))
	pylab.savefig("all_w2.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_1000kpc], theta[select_1000kpc], None)
	pylab.savefig("1000kpc.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_500kpc], theta[select_500kpc], None)
	pylab.savefig("500kpc.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_200kpc], theta[select_200kpc], None)
	pylab.savefig("200kpc.png")
	pylab.show()
	pylab.close()


	plot_uhecrs(phi[select_100kpc], theta[select_100kpc], None)
	pylab.savefig("100kpc.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_50kpc], theta[select_50kpc], None)
	pylab.savefig("50kpc.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_closest], theta[select_closest], None)
	pylab.savefig("closest.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_closest], theta[select_closest], 1./dist[select_closest])
	pylab.savefig("closest_w.png")
	pylab.show()
	pylab.close()

	plot_uhecrs(phi[select_closest], theta[select_closest], 1./(dist[select_closest]**2))
	pylab.savefig("closest_w2.png")
	pylab.show()
	pylab.close()

		
if len(sys.argv) > 1:
	if sys.argv[1] == "sim":
		simulate()
	elif sys.argv[1] == "plot":
		plot()
else:
	simulate()
	plot()
