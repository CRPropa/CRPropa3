#include "crpropa/GridTools.h"
#include "crpropa/Random.h"
#include "crpropa/magneticField/MagneticField.h"

#include <fstream>
#include <sstream>

namespace crpropa {

void scaleGrid(ref_ptr<ScalarGrid> grid, double a) {
	for (int ix = 0; ix < grid->getNx(); ix++)
		for (int iy = 0; iy < grid->getNy(); iy++)
			for (int iz = 0; iz < grid->getNz(); iz++)
				grid->get(ix, iy, iz) *= a;
}

void scaleGrid(ref_ptr<VectorGrid> grid, double a) {
	for (int ix = 0; ix < grid->getNx(); ix++)
		for (int iy = 0; iy < grid->getNy(); iy++)
			for (int iz = 0; iz < grid->getNz(); iz++)
				grid->get(ix, iy, iz) *= a;
}

Vector3f meanFieldVector(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	Vector3f mean(0.);
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double meanFieldStrength(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double mean = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz).getR();
	return mean / Nx / Ny / Nz;
}

double meanFieldStrength(ref_ptr<ScalarGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double mean = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				mean += grid->get(ix, iy, iz);
	return mean / Nx / Ny / Nz;
}

double rmsFieldStrength(ref_ptr<VectorGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double sumV2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumV2 += grid->get(ix, iy, iz).getR2();
	return std::sqrt(sumV2 / Nx / Ny / Nz);
}

double rmsFieldStrength(ref_ptr<ScalarGrid> grid) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	double sumV2 = 0;
	for (int ix = 0; ix < Nx; ix++)
		for (int iy = 0; iy < Ny; iy++)
			for (int iz = 0; iz < Nz; iz++)
				sumV2 += pow(grid->get(ix, iy, iz), 2);
	return std::sqrt(sumV2 / Nx / Ny / Nz);
}

#ifdef CRPROPA_HAVE_FFTW3F
#include "fftw3.h"

void initTurbulence(ref_ptr<VectorGrid> grid, double Brms, double lMin,
		double lMax, double alpha, int seed) {
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	if ((Nx != Ny) or (Ny != Nz))
		throw std::runtime_error("turbulentField: only cubic grid supported");

	double spacing = grid->getSpacing();
	if (lMin < 2 * spacing)
		throw std::runtime_error("turbulentField: lMin < 2 * spacing");
	if (lMin >= lMax)
		throw std::runtime_error("turbulentField: lMin >= lMax");
	if (lMax > Nx * spacing / 2)
		throw std::runtime_error("turbulentField: lMax > size / 2");

	size_t n = Nx; // size of array
	size_t n2 = (size_t) floor(n / 2) + 1; // size array in z-direction in configuration space

	// arrays to hold the complex vector components of the B(k)-field
	fftwf_complex *Bkx, *Bky, *Bkz;
	Bkx = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bky = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);
	Bkz = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n * n * n2);

	Random random;
	if (seed != 0)
		random.seed(seed); // use given seed

	// calculate the n possible discrete wave numbers
	double K[n];
	for (int i = 0; i < n; i++)
		K[i] = (double) i / n - i / (n / 2);

	// construct the field in configuration space
	int i;
	double k, theta, phase, cosPhase, sinPhase;
	double kMin = spacing / lMax;
	double kMax = spacing / lMin;
	Vector3f b; // real b-field vector
	Vector3f ek, e1, e2; // orthogonal base
	Vector3f n0(1, 1, 1); // arbitrary vector to construct orthogonal base

	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n2; iz++) {

				i = ix * n * n2 + iy * n2 + iz;
				ek.setXYZ(K[ix], K[iy], K[iz]);
				k = ek.getR();

				// wave outside of turbulent range -> B(k) = 0
				if ((k < kMin) || (k > kMax)) {
					Bkx[i][0] = 0;
					Bkx[i][1] = 0;
					Bky[i][0] = 0;
					Bky[i][1] = 0;
					Bkz[i][0] = 0;
					Bkz[i][1] = 0;
					continue;
				}

				// construct an orthogonal base ek, e1, e2
				if (ek.isParallelTo(n0, float(1e-3))) {
					// ek parallel to (1,1,1)
					e1.setXYZ(-1., 1., 0);
					e2.setXYZ(1., 1., -2.);
				} else {
					// ek not parallel to (1,1,1)
					e1 = n0.cross(ek);
					e2 = ek.cross(e1);
				}
				e1 /= e1.getR();
				e2 /= e2.getR();

				// random orientation perpendicular to k
				theta = 2 * M_PI * random.rand();
				b = e1 * cos(theta) + e2 * sin(theta);

				// normal distributed amplitude with mean = 0 and sigma = k^alpha/2
				b *= random.randNorm() * pow(k, alpha / 2);

				// uniform random phase
				phase = 2 * M_PI * random.rand();
				cosPhase = cos(phase); // real part
				sinPhase = sin(phase); // imaginary part

				Bkx[i][0] = b.x * cosPhase;
				Bkx[i][1] = b.x * sinPhase;
				Bky[i][0] = b.y * cosPhase;
				Bky[i][1] = b.y * sinPhase;
				Bkz[i][0] = b.z * cosPhase;
				Bkz[i][1] = b.z * sinPhase;
			}
		}
	}

	// in-place, complex to real, inverse Fourier transformation on each component
	// note that the last elements of B(x) are unused now
	float *Bx = (float*) Bkx;
	fftwf_plan plan_x = fftwf_plan_dft_c2r_3d(n, n, n, Bkx, Bx, FFTW_ESTIMATE);
	fftwf_execute(plan_x);
	fftwf_destroy_plan(plan_x);

	float *By = (float*) Bky;
	fftwf_plan plan_y = fftwf_plan_dft_c2r_3d(n, n, n, Bky, By, FFTW_ESTIMATE);
	fftwf_execute(plan_y);
	fftwf_destroy_plan(plan_y);

	float *Bz = (float*) Bkz;
	fftwf_plan plan_z = fftwf_plan_dft_c2r_3d(n, n, n, Bkz, Bz, FFTW_ESTIMATE);
	fftwf_execute(plan_z);
	fftwf_destroy_plan(plan_z);

	// save to grid
	for (size_t ix = 0; ix < n; ix++) {
		for (size_t iy = 0; iy < n; iy++) {
			for (size_t iz = 0; iz < n; iz++) {
				i = ix * n * 2 * n2 + iy * 2 * n2 + iz;
				Vector3f &b = grid->get(ix, iy, iz);
				b.x = Bx[i];
				b.y = By[i];
				b.z = Bz[i];
			}
		}
	}

	fftwf_free(Bkx);
	fftwf_free(Bky);
	fftwf_free(Bkz);

	scaleGrid(grid, Brms / rmsFieldStrength(grid)); // normalize to Brms
}
#endif // CRPROPA_HAVE_FFTW3F

void fromMagneticField(ref_ptr<VectorGrid> grid, ref_ptr<MagneticField> field) {
	Vector3d origin = grid->getOrigin();
	double spacing = grid->getSpacing();
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	for (size_t ix = 0; ix < Nx; ix++)
		for (size_t iy = 0; iy < Ny; iy++)
			for (size_t iz = 0; iz < Nz; iz++) {
				Vector3d pos = Vector3d(double(ix) + 0.5, double(iy) + 0.5, double(iz) + 0.5) * spacing + origin;
				Vector3d B = field->getField(pos);
				grid->get(ix, iy, iz) = B;
	}
}

void fromMagneticFieldStrength(ref_ptr<ScalarGrid> grid, ref_ptr<MagneticField> field) {
	Vector3d origin = grid->getOrigin();
	double spacing = grid->getSpacing();
	size_t Nx = grid->getNx();
	size_t Ny = grid->getNy();
	size_t Nz = grid->getNz();
	for (size_t ix = 0; ix < Nx; ix++)
		for (size_t iy = 0; iy < Ny; iy++)
			for (size_t iz = 0; iz < Nz; iz++) {
				Vector3d pos = Vector3d(double(ix) + 0.5, double(iy) + 0.5, double(iz) + 0.5) * spacing + origin;
				double s = field->getField(pos).getR();
				grid->get(ix, iy, iz) = s;
	}
}

double turbulentCorrelationLength(double lMin, double lMax, double alpha) {
	double r = lMin / lMax;
	double a = -alpha - 2;
	return lMax / 2 * (a - 1) / a * (1 - pow(r, a)) / (1 - pow(r, a - 1));
}

void loadGrid(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg (0, fin.beg);

	size_t nx = grid->getNx();
	size_t ny = grid->getNy();
	size_t nz = grid->getNz();

	if (length != (3 * nx * ny * nz))
		throw std::runtime_error("loadGrid: file and grid size do not match");

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f &b = grid->get(ix, iy, iz);
				fin.read((char*) &(b.x), sizeof(float));
				fin.read((char*) &(b.y), sizeof(float));
				fin.read((char*) &(b.z), sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void loadGrid(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str(), std::ios::binary);
	if (!fin) {
		std::stringstream ss;
		ss << "load ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}

	// get length of file and compare to size of grid
	fin.seekg(0, fin.end);
	size_t length = fin.tellg() / sizeof(float);
	fin.seekg (0, fin.beg);

	size_t nx = grid->getNx();
	size_t ny = grid->getNy();
	size_t nz = grid->getNz();

	if (length != (nx * ny * nz))
		throw std::runtime_error("loadGrid: file and grid size do not match");

	for (int ix = 0; ix < nx; ix++) {
		for (int iy = 0; iy < ny; iy++) {
			for (int iz = 0; iz < nz; iz++) {
				float &b = grid->get(ix, iy, iz);
				fin.read((char*) &b, sizeof(float));
				b *= c;
			}
		}
	}
	fin.close();
}

void dumpGrid(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout) {
		std::stringstream ss;
		ss << "dump VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f b = grid->get(ix, iy, iz) * c;
				fout.write((char*) &(b.x), sizeof(float));
				fout.write((char*) &(b.y), sizeof(float));
				fout.write((char*) &(b.z), sizeof(float));
			}
		}
	}
	fout.close();
}

void dumpGrid(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str(), std::ios::binary);
	if (!fout) {
		std::stringstream ss;
		ss << "dump ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float b = grid->get(ix, iy, iz) * c;
				fout.write((char*) &b, sizeof(float));
			}
		}
	}
	fout.close();
}

void loadGridFromTxt(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	// skip header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f &b = grid->get(ix, iy, iz);
				fin >> b.x >> b.y >> b.z;
				b *= c;
				if (fin.eof())
					throw std::runtime_error("load VectorGrid: file too short");
			}
		}
	}
	fin.close();
}

void loadGridFromTxt(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ifstream fin(filename.c_str());
	if (!fin) {
		std::stringstream ss;
		ss << "load ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	// skip header lines
	while (fin.peek() == '#')
		fin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float &b = grid->get(ix, iy, iz);
				fin >> b;
				b *= c;
				if (fin.eof())
					throw std::runtime_error("load ScalarGrid: file too short");
			}
		}
	}
	fin.close();
}

void dumpGridToTxt(ref_ptr<VectorGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str());
	if (!fout) {
		std::stringstream ss;
		ss << "dump VectorGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				Vector3f b = grid->get(ix, iy, iz) * c;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

void dumpGridToTxt(ref_ptr<ScalarGrid> grid, std::string filename, double c) {
	std::ofstream fout(filename.c_str());
	if (!fout) {
		std::stringstream ss;
		ss << "dump ScalarGrid: " << filename << " not found";
		throw std::runtime_error(ss.str());
	}
	for (int ix = 0; ix < grid->getNx(); ix++) {
		for (int iy = 0; iy < grid->getNy(); iy++) {
			for (int iz = 0; iz < grid->getNz(); iz++) {
				float b = grid->get(ix, iy, iz) * c;
				fout << b << "\n";
			}
		}
	}
	fout.close();
}

} // namespace crpropa
