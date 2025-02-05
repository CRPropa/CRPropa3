#ifndef _KST24_GMF_H_
#define _KST24_GMF_H_

#include <vector>
#include "crpropa/magneticField/MagneticField.h"

namespace crpropa {

/* 
The C++ implementation of the GMF model KST24 (A.Korochkin, D.Semikoz, P.Tinyakov 2024)
The model was presented in arXiv:2407.02148 and published in A&A
If you use the model, please cite A&A, 693, A284 (2025)

In KST24 GMF model the position of the Solar System is at {-8.2 kpc, 0, 0}
The Galactic Center is at {0, 0, 0}
z-axis points to the North pole
*/


class KST24Field : public MagneticField
{
private:
	double north_tor_B_gauss;
	double north_tor_zmin_kpc;
	double north_tor_zmax_kpc;
	double north_tor_rmin_kpc;
	double north_tor_rmax_kpc;

	double south_tor_B_gauss;
	double south_tor_zmin_kpc;
	double south_tor_zmax_kpc;
	double south_tor_rmin_kpc;
	double south_tor_rmax_kpc;

	double Xfield_B_gauss;
	double Xfield_rmin_kpc;
	double Xfield_rmax_kpc;
	double Xfield_theta_deg;
	double Xfield_theta_rad;

	double LB_B_gauss;
	double LB_lB_deg;
	double LB_bB_deg;
	double LB_rmin_kpc; 
	double LB_dr_kpc;
	double LB_x0_kpc;
	double LB_y0_kpc;
	double LB_z0_kpc;
	double LB_Bdir_x;
	double LB_Bdir_y;
	double LB_Bdir_z;

	double ScutumArm_B_gauss;
	double ScutumArm_pitch_deg;
	double ScutumArm_phi0_deg;
	double ScutumArm_x_shift_kpc;
	double ScutumArm_y_shift_kpc;
	double ScutumArm_arc_radius1_kpc;
	double ScutumArm_arc_radius2_kpc;
	double ScutumArm_arc_eps;
	double ScutumArm_arc_div_deg; 
	double ScutumArm_rmin_kpc;
	double ScutumArm_rmax_kpc;
	double ScutumArm_zmin_kpc;
	double ScutumArm_zmax_kpc;

	double CarinaSagittariusArm_B_gauss;
	double CarinaSagittariusArm_pitch_deg;
	double CarinaSagittariusArm_phi0_deg;
	double CarinaSagittariusArm_x_shift_kpc;
	double CarinaSagittariusArm_y_shift_kpc;
	double CarinaSagittariusArm_arc_radius1_kpc;
	double CarinaSagittariusArm_arc_radius2_kpc;
	double CarinaSagittariusArm_arc_eps;
	double CarinaSagittariusArm_arc_div_deg; 
	double CarinaSagittariusArm_rmin_kpc;
	double CarinaSagittariusArm_rmax_kpc;
	double CarinaSagittariusArm_zmin_kpc;
	double CarinaSagittariusArm_zmax_kpc;

	double LocalArm_B_gauss;
	double LocalArm_pitch_deg;
	double LocalArm_phi0_deg;
	double LocalArm_x_shift_kpc;
	double LocalArm_y_shift_kpc;
	double LocalArm_arc_radius1_kpc;
	double LocalArm_arc_radius2_kpc;
	double LocalArm_arc_eps;
	double LocalArm_arc_div_deg; 
	double LocalArm_rmin_kpc;
	double LocalArm_rmax_kpc;
	double LocalArm_zmin_kpc;
	double LocalArm_zmax_kpc;

	double PerseusArm1_B_gauss;
	double PerseusArm1_pitch_deg;
	double PerseusArm1_phi0_deg;
	double PerseusArm1_x_shift_kpc;
	double PerseusArm1_y_shift_kpc;
	double PerseusArm1_arc_radius1_kpc;
	double PerseusArm1_arc_radius2_kpc;
	double PerseusArm1_arc_eps;
	double PerseusArm1_arc_div_deg; 
	double PerseusArm1_rmin_kpc;
	double PerseusArm1_rmax_kpc;
	double PerseusArm1_zmin_kpc;
	double PerseusArm1_zmax_kpc;

	double PerseusArm2_B_gauss;
	double PerseusArm2_pitch_deg;
	double PerseusArm2_phi0_deg;
	double PerseusArm2_x_shift_kpc;
	double PerseusArm2_y_shift_kpc;
	double PerseusArm2_arc_radius1_kpc;
	double PerseusArm2_arc_radius2_kpc;
	double PerseusArm2_arc_eps;
	double PerseusArm2_arc_div_deg; 
	double PerseusArm2_rmin_kpc;
	double PerseusArm2_rmax_kpc;
	double PerseusArm2_zmin_kpc;
	double PerseusArm2_zmax_kpc;

public:
	Vector3d getField(const Vector3d& pos) const;
	KST24Field();

	Vector3d get_toroidal(const Vector3d pos_kpc, const double tor_B_gauss, 
						  const double tor_zmin_kpc, const double tor_zmax_kpc, 
						  const double tor_rmin_kpc, const double tor_rmax_kpc) const;

	Vector3d get_Xfield(const Vector3d pos_kpc, const double Xfield_B_gauss, 
						const double Xfield_rmin_kpc, const double Xfield_rmax_kpc, 
						const double Xfield_theta_rad) const;

	bool is_LB(const Vector3d pos_kpc, const double LB_rmin_kpc, const double LB_dr_kpc,
			   const double LB_x0_kpc, const double LB_y0_kpc, const double LB_z0_kpc) const;

	Vector3d get_LB(const Vector3d pos_kpc, const double LB_B_gauss, 
					const double LB_lB_deg, const double LB_bB_deg, 
					const double LB_rmin_kpc, const double LB_dr_kpc,
					const double LB_x0_kpc, const double LB_y0_kpc, const double LB_z0_kpc) const;

	Vector3d get_logspiral(const Vector3d pos_kpc, const double B_gauss,  
						   const double pitch_deg, const double phi0_deg,
						   const double x_shift_kpc, const double y_shift_kpc,
						   const double arc_radius1_kpc, const double arc_radius2_kpc, 
						   const double arc_eps , const double arc_div_deg, 
						   const double rmin_kpc, const double rmax_kpc, 
						   const double zmin_kpc, const double zmax_kpc) const;
};

}// namespace crpropa 

#endif /* _KST24_GMF_H_ */