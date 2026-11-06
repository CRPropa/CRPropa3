#ifndef CRPROPA_UF26FIELD_H
#define CRPROPA_UF26FIELD_H

#include <string>
#include "crpropa/magneticField/MagneticField.h"
#include "crpropa/Grid.h"

namespace crpropa {
/**
 * \addtogroup MagneticFields
 * @{
 */

/**
 @class UF26Field
 @brief UF26 model of the random Galactic magnetic field

 Implements the two random-field models of UF26, expDisk and ringDisk
 (Table 2 of M. Unger and G.R. Farrar, arXiv:2608.21293). The coherent
 field of the same authors is available as UF23Field.

 The models describe the rms strength b(x) of the isotropic random field,
 i.e. <|B|^2> = b^2, where the average is over a volume large compared to
 the coherence length of the field. Each Cartesian component therefore has
 an rms of b/sqrt(3). The rms strength is given by getTurbulentStrength().

 As for the turbulent component of JF12Field, a realization of the field is
 obtained by scaling a turbulent grid with B_rms = 1 to the local rms
 strength, B(x) = b(x) * B_grid(x). The grid is created with
 randomTurbulent() or set with setTurbulentGrid(). The paper does not
 determine the coherence length; randomTurbulent() uses Lc = 30 pc, the value
 used for the deflection estimates in the paper. The scaling breaks
 div B = 0 at the level of Lc divided by the scale of variation of b
 (>~ 1 kpc).

 Usage together with the coherent UF23 field, e.g. for backtracking or
 lens production (Python):

 \code{.py}
 B = MagneticFieldList()
 B.addField(UF23Field(UF23Field.base))
 B26 = UF26Field(UF26Field.expDisk)
 B26.randomTurbulent(seed)
 B.addField(B26)
 \endcode

 For a different spectrum, create the grid e.g. with SimpleGridTurbulence
 (Brms = 1) and pass it with setTurbulentGrid().

 Assumes a galactocentric coordinate system with the Galactic center
 at the origin, the x-axis pointing in the opposite direction of the
 Sun, and the z-axis pointing towards Galactic North. The Sun is at
 a distance of 8.178 kpc from the Galactic center. UF23 and UF26 extend
 beyond 20 kpc (UF23 to 30 kpc, UF26 has b > 0.05 muG out to ~25 kpc).

 Reference implementation: github.com/MichaelUnger/GalacticRandomFieldUF26
 */
class UF26Field : public MagneticField {
public:
	/// models of Table 2 of the UF26 paper
	enum ModelType {
		expDisk,
		ringDisk
	};

	/**
	 @brief constructor
	 @param mt model type (expDisk or ringDisk)
	 */
	UF26Field(const ModelType mt);

	/// no default constructor
	UF26Field() = delete;

#ifdef CRPROPA_HAVE_FFTW3F
	/**
	 @brief create and set a random realization of the turbulent grid

	 Kolmogorov spectrum between 4 and 136 pc (Lc = 30 pc) with B_rms = 1 on
	 a periodic grid of 256^3 points with 2 pc spacing (512 pc, 200 MB).
	 */
	void randomTurbulent(int seed = 0);
#endif

	/**
	 Set a turbulent grid
	 @param grid	vector grid containing a random field of Brms = 1
	 */
	void setTurbulentGrid(ref_ptr<Grid3f> grid);
	ref_ptr<Grid3f> getTurbulentGrid();
	/// true if a turbulent grid is set
	bool isUsingTurbulentField() const;

	/// rms strength b of the random field (Eqs. 12 and 15 in quadrature)
	double getTurbulentStrength(const Vector3d &pos) const;

	/**
	 @brief realization of the random field, b(pos) * B_grid(pos)
	 @throws std::runtime_error if no turbulent grid is set
	 */
	Vector3d getTurbulentField(const Vector3d &pos) const;

	/// same as getTurbulentField()
	Vector3d getField(const Vector3d &pos) const;

	ModelType getModelType() const;
	std::string getModelName() const;

private:
	// rms strengths in microgauss for positions in kpc
	double rmsDisk(const Vector3d &posInKpc) const;
	double rmsAnnulus(const Vector3d &posInKpc) const;

	/// model type given in constructor
	const ModelType fModelType;

	// disk, Eqs. (12)-(14) [kpc, microgauss]
	double fDiskB = 0;     ///< b_d at the solar circle
	double fDiskZ = 0;     ///< 1/e height z_d
	double fDiskL = 0;     ///< radial scale length l_r
	double fDiskRc = 0;    ///< cusp smoothing radius r_c
	double fDiskRmax = 0;  ///< outer cutoff radius R_max
	double fDiskWmax = 0;  ///< width of outer cutoff w_max
	// annulus, Eq. (15) [kpc, microgauss]
	double fAnnulusB = 0;  ///< peak amplitude b_a,0
	double fAnnulusR = 0;  ///< peak radius r_a
	double fAnnulusW = 0;  ///< radial Gaussian width w_a
	double fAnnulusZ = 0;  ///< 1/e height z_a
	bool fHasAnnulus = false;

	// pre-calculated: sqrt(R_sun^2 + r_c^2)
	double fSunRadiusSmoothed = 0;

	// turbulent grid with B_rms = 1
	ref_ptr<Grid3f> turbulentGrid;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_UF26FIELD_H
