#ifndef CRPROPA_UF26FIELD_H
#define CRPROPA_UF26FIELD_H

#include <string>
#include "crpropa/magneticField/MagneticField.h"

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
 an rms of b/sqrt(3).

 - expDisk: vertically thin disk, exponential in radius with a smoothed
   central cusp, sech vertical profile and a logistic cutoff at 18 kpc.
 - ringDisk: radially flat disk that is phased out inside 4.5 kpc, plus an
   annular enhancement with a Gaussian vertical profile, added in quadrature.

 The models specify only b(x), not the spectrum of the field. getField()
 returns a realization of the random field by scaling a homogeneous
 turbulent field with unit rms to the local rms strength,

   B(x) = b(x) / B_rms * B_turb(x),

 where B_turb is set with setTurbulentField(), e.g. a PlaneWaveTurbulence or
 GridTurbulence with a user-defined spectrum. The paper does not determine
 the coherence length (Lc = 30 pc was used for the deflection estimates).
 The modulation breaks div B = 0 at the level of Lc divided by the scale of
 variation of b (>~ 1 kpc), as for the turbulent component of JF12Field.

 The rms strength is available without a turbulent field via getRmsField().

 Usage together with the coherent UF23 field, e.g. for backtracking or
 lens production (Python):

 \code{.py}
 B = MagneticFieldList()
 B.addField(UF23Field(UF23Field.base))
 spectrum = SimpleTurbulenceSpectrum(1, 4 * parsec, 136 * parsec)  # Kolmogorov, Lc = 30 pc
 grid = GridProperties(Vector3d(0), 256, 2 * parsec)  # periodic 512 pc box, 200 MB
 B.addField(UF26Field(UF26Field.expDisk, SimpleGridTurbulence(spectrum, grid, seed)))
 \endcode

 This corresponds to JF12Field with randomTurbulent(). SimpleGridTurbulence
 needs FFTW. PlaneWaveTurbulence works without it but is much slower per
 field evaluation. Note that UF23 and UF26 use R_sun = 8.178 kpc and extend
 beyond 20 kpc (UF23 to 30 kpc, UF26 has b > 0.05 muG out to ~25 kpc).

 Assumes a galactocentric coordinate system with the Galactic center
 at the origin, the x-axis pointing in the opposite direction of the
 Sun, and the z-axis pointing towards Galactic North. The Sun is at
 a distance of 8.178 kpc from the Galactic center.

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
	 @brief constructor without turbulent field (only getRmsField() is usable)
	 @param mt model type (expDisk or ringDisk)
	 */
	UF26Field(const ModelType mt);

	/**
	 @brief constructor with turbulent field, see setTurbulentField()
	 @param mt         model type (expDisk or ringDisk)
	 @param turbulence homogeneous turbulent field (a TurbulentField)
	 */
	UF26Field(const ModelType mt, ref_ptr<MagneticField> turbulence);

	/// no default constructor
	UF26Field() = delete;

	/**
	 @brief realization of the random field, b(pos) / B_rms * B_turb(pos)
	 @throws std::runtime_error if no turbulent field is set
	 */
	Vector3d getField(const Vector3d &pos) const;

	/// rms strength b of the random field (Eqs. 12 and 15 in quadrature)
	double getRmsField(const Vector3d &pos) const;
	/// rms strength of the disk component (Eq. 12)
	double getDiskRmsField(const Vector3d &pos) const;
	/// rms strength of the annular component (Eq. 15, zero for expDisk)
	double getAnnularRmsField(const Vector3d &pos) const;

	/**
	 @brief set the homogeneous turbulent field to be scaled by b(pos)
	 @param turbulence a TurbulentField (e.g. PlaneWaveTurbulence or
	                   GridTurbulence); its rms is taken from its spectrum
	 @throws std::runtime_error if turbulence is not a TurbulentField
	 */
	void setTurbulentField(ref_ptr<MagneticField> turbulence);

	/**
	 @brief set a homogeneous random field with known rms, e.g. a
	        MagneticFieldGrid read from file
	 @param field random field
	 @param brms  rms strength sqrt(<|B|^2>) of field
	 */
	void setTurbulentField(ref_ptr<MagneticField> field, double brms);

	ref_ptr<MagneticField> getTurbulentField() const;
	bool hasTurbulentField() const;

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

	// homogeneous turbulent field and its rms
	ref_ptr<MagneticField> fTurbulence;
	double fTurbulenceBrms = 0;
};

/** @}*/
} // namespace crpropa

#endif // CRPROPA_UF26FIELD_H
