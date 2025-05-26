#include "crpropa/magneticField/KST24Field.h"

/* 
The C++ implementation of the GMF model KST24 (A.Korochkin, D.Semikoz, P.Tinyakov 2024)
The model was presented in arXiv:2407.02148 and published in A&A
If you use the model, please cite A&A, 693, A284 (2025)

In KST24 GMF model the position of the Solar System is at {-8.2 kpc, 0, 0}
The Galactic Center is at {0, 0, 0}
z-axis points to the North pole
*/


namespace crpropa {
KST24Field::KST24Field()
{
	north_tor_B_gauss  = -3.2e-6;
	north_tor_zmin_kpc = 1.185;
	north_tor_zmax_kpc = 2.1;
	north_tor_rmin_kpc = 1.0;
	north_tor_rmax_kpc = 9.1;

	south_tor_B_gauss  = 3.2e-6;
	south_tor_zmin_kpc = -2.5;
	south_tor_zmax_kpc = -1.22;
	south_tor_rmin_kpc = 1.0;
	south_tor_rmax_kpc = 14.0;

	Xfield_B_gauss   = 1.8e-6;
	Xfield_rmin_kpc  = 1.0;
	Xfield_rmax_kpc  = 6.2;
	Xfield_theta_deg = 28.0;
	Xfield_theta_rad = Xfield_theta_deg*M_PI/180.;

	LB_B_gauss  = -3.5e-6;
	LB_lB_deg   = 50.0;
	LB_bB_deg   = 2.0;
	LB_rmin_kpc = 0.2; 
	LB_dr_kpc   = 0.03;
	LB_x0_kpc   = -8.2;
	LB_y0_kpc   = 0.095;
	LB_z0_kpc   = -0.05;
	LB_Bdir_x   = cos(LB_bB_deg*M_PI/180)*cos(LB_lB_deg*M_PI/180);
	LB_Bdir_y   = cos(LB_bB_deg*M_PI/180)*sin(LB_lB_deg*M_PI/180);
	LB_Bdir_z   = sin(LB_bB_deg*M_PI/180);


	ScutumArm_B_gauss         = 4.9e-6;
	ScutumArm_pitch_deg       = 20;
	ScutumArm_phi0_deg        = -134.0;
	ScutumArm_x_shift_kpc     = 1.0;
	ScutumArm_y_shift_kpc     = 0.0;
	ScutumArm_arc_radius1_kpc = 0.8;
	ScutumArm_arc_radius2_kpc = 1.0;
	ScutumArm_arc_eps         = 2.0;
	ScutumArm_arc_div_deg     = 0.0;
	ScutumArm_rmin_kpc        = 3.0;
	ScutumArm_rmax_kpc        = 16.0;
	ScutumArm_zmin_kpc        = -1.0;
	ScutumArm_zmax_kpc        = 1.0;	

	CarinaSagittariusArm_B_gauss         = 1.3e-6;
	CarinaSagittariusArm_pitch_deg       = 20;
	CarinaSagittariusArm_phi0_deg        = -80.0;
	CarinaSagittariusArm_x_shift_kpc     = 1.37;
	CarinaSagittariusArm_y_shift_kpc     = 0.0;
	CarinaSagittariusArm_arc_radius1_kpc = 1.0;
	CarinaSagittariusArm_arc_radius2_kpc = 0.79;
	CarinaSagittariusArm_arc_eps         = 2.3;
	CarinaSagittariusArm_arc_div_deg     = 3.0;
	CarinaSagittariusArm_rmin_kpc        = 3.0;
	CarinaSagittariusArm_rmax_kpc        = 15.0;
	CarinaSagittariusArm_zmin_kpc        = -1.2;
	CarinaSagittariusArm_zmax_kpc        = 1.2;

	LocalArm_B_gauss         = -3.5e-6;
	LocalArm_pitch_deg       = 20;
	LocalArm_phi0_deg        = -2.2;
	LocalArm_x_shift_kpc     = -0.15;
	LocalArm_y_shift_kpc     = 0.0;
	LocalArm_arc_radius1_kpc = 0.73;
	LocalArm_arc_radius2_kpc = 1.0;
	LocalArm_arc_eps         = 1.45;
	LocalArm_arc_div_deg     = 0.0;
	LocalArm_rmin_kpc        = 3.3;
	LocalArm_rmax_kpc        = 14.0;
	LocalArm_zmin_kpc        = -1.0;
	LocalArm_zmax_kpc        = 1.0;

	PerseusArm1_B_gauss         = -2.0e-6;
	PerseusArm1_pitch_deg       = 20;
	PerseusArm1_phi0_deg        = 46.0;
	PerseusArm1_x_shift_kpc     = -1.0;
	PerseusArm1_y_shift_kpc     = 0.0;
	PerseusArm1_arc_radius1_kpc = 0.4;
	PerseusArm1_arc_radius2_kpc = 1.1;
	PerseusArm1_arc_eps         = 2.0;
	PerseusArm1_arc_div_deg     = 0.0;
	PerseusArm1_rmin_kpc        = 3.0;
	PerseusArm1_rmax_kpc        = 17.0;
	PerseusArm1_zmin_kpc        = -1.0;
	PerseusArm1_zmax_kpc        = 1.0;

	PerseusArm2_B_gauss         = -3.5e-6;
	PerseusArm2_pitch_deg       = 20;
	PerseusArm2_phi0_deg        = 46.0;
	PerseusArm2_x_shift_kpc     = -1.0;
	PerseusArm2_y_shift_kpc     = 0.0;
	PerseusArm2_arc_radius1_kpc = 1.2;
	PerseusArm2_arc_radius2_kpc = 1.1;
	PerseusArm2_arc_eps         = 2.0;
	PerseusArm2_arc_div_deg     = 0.0;
	PerseusArm2_rmin_kpc        = 3.0;
	PerseusArm2_rmax_kpc        = 17.0;
	PerseusArm2_zmin_kpc        = -2.0;
	PerseusArm2_zmax_kpc        = 2.0;
}



Vector3d KST24Field::getField(const Vector3d& pos) const
{
	Vector3d vals(0, 0, 0);
	Vector3d pos_kpc = pos/kpc;

	if (pos_kpc.z > north_tor_zmin_kpc) // northern halo height lower boundary
	{
		vals += get_toroidal(pos_kpc, north_tor_B_gauss, north_tor_zmin_kpc, north_tor_zmax_kpc, 
							 north_tor_rmin_kpc, north_tor_rmax_kpc);
		vals += get_Xfield(pos_kpc, Xfield_B_gauss, Xfield_rmin_kpc, Xfield_rmax_kpc, Xfield_theta_rad);
		return vals*gauss;
	}

	if (pos_kpc.z < south_tor_zmax_kpc) // southern halo height upper boundary
	{
		vals += get_toroidal(pos_kpc, south_tor_B_gauss, south_tor_zmin_kpc, south_tor_zmax_kpc, 
							 south_tor_rmin_kpc, south_tor_rmax_kpc);
		vals += get_Xfield(pos_kpc, Xfield_B_gauss, Xfield_rmin_kpc, Xfield_rmax_kpc, Xfield_theta_rad);
		return vals*gauss;
	}

	if (is_LB(pos_kpc, LB_rmin_kpc, LB_dr_kpc, LB_x0_kpc, LB_y0_kpc, LB_z0_kpc))
	{
		vals += get_LB(pos_kpc, LB_B_gauss, LB_lB_deg, LB_bB_deg, 
					   LB_rmin_kpc, LB_dr_kpc, LB_x0_kpc, LB_y0_kpc, LB_z0_kpc);
		return vals*gauss;
	}

	vals += get_logspiral(pos_kpc, ScutumArm_B_gauss, ScutumArm_pitch_deg, ScutumArm_phi0_deg,
						  ScutumArm_x_shift_kpc, ScutumArm_y_shift_kpc, ScutumArm_arc_radius1_kpc,
						  ScutumArm_arc_radius2_kpc, ScutumArm_arc_eps, ScutumArm_arc_div_deg,
						  ScutumArm_rmin_kpc, ScutumArm_rmax_kpc, ScutumArm_zmin_kpc, ScutumArm_zmax_kpc);

	vals += get_logspiral(pos_kpc, CarinaSagittariusArm_B_gauss, CarinaSagittariusArm_pitch_deg, CarinaSagittariusArm_phi0_deg,
						  CarinaSagittariusArm_x_shift_kpc, CarinaSagittariusArm_y_shift_kpc, CarinaSagittariusArm_arc_radius1_kpc,
						  CarinaSagittariusArm_arc_radius2_kpc, CarinaSagittariusArm_arc_eps, CarinaSagittariusArm_arc_div_deg,
						  CarinaSagittariusArm_rmin_kpc, CarinaSagittariusArm_rmax_kpc, CarinaSagittariusArm_zmin_kpc, CarinaSagittariusArm_zmax_kpc);

	vals += get_logspiral(pos_kpc, LocalArm_B_gauss, LocalArm_pitch_deg, LocalArm_phi0_deg,
						  LocalArm_x_shift_kpc, LocalArm_y_shift_kpc, LocalArm_arc_radius1_kpc,
						  LocalArm_arc_radius2_kpc, LocalArm_arc_eps, LocalArm_arc_div_deg,
						  LocalArm_rmin_kpc, LocalArm_rmax_kpc, LocalArm_zmin_kpc, LocalArm_zmax_kpc);

	vals += get_logspiral(pos_kpc, PerseusArm1_B_gauss, PerseusArm1_pitch_deg, PerseusArm1_phi0_deg,
						  PerseusArm1_x_shift_kpc, PerseusArm1_y_shift_kpc, PerseusArm1_arc_radius1_kpc,
						  PerseusArm1_arc_radius2_kpc, PerseusArm1_arc_eps, PerseusArm1_arc_div_deg,
						  PerseusArm1_rmin_kpc, PerseusArm1_rmax_kpc, PerseusArm1_zmin_kpc, PerseusArm1_zmax_kpc);

	vals += get_logspiral(pos_kpc, PerseusArm2_B_gauss, PerseusArm2_pitch_deg, PerseusArm2_phi0_deg,
						  PerseusArm2_x_shift_kpc, PerseusArm2_y_shift_kpc, PerseusArm2_arc_radius1_kpc,
						  PerseusArm2_arc_radius2_kpc, PerseusArm2_arc_eps, PerseusArm2_arc_div_deg,
						  PerseusArm2_rmin_kpc, PerseusArm2_rmax_kpc, PerseusArm2_zmin_kpc, PerseusArm2_zmax_kpc);

	vals += get_Xfield(pos_kpc, Xfield_B_gauss, Xfield_rmin_kpc, Xfield_rmax_kpc, Xfield_theta_rad);

	return vals*gauss;
}

bool KST24Field::is_LB(const Vector3d pos_kpc, const double LB_rmin_kpc, const double LB_dr_kpc,
					   const double LB_x0_kpc, const double LB_y0_kpc, const double LB_z0_kpc) const
{
	double eRx, eRy, eRz;
	double cR;

	eRx = pos_kpc.x - LB_x0_kpc;
	eRy = pos_kpc.y - LB_y0_kpc;
	eRz = pos_kpc.z - LB_z0_kpc;
	cR = sqrt(eRx*eRx + eRy*eRy + eRz*eRz);
	if (cR < LB_rmin_kpc + LB_dr_kpc)
		return true;
	return false;
}

Vector3d KST24Field::get_toroidal(const Vector3d pos_kpc, const double tor_B_gauss, 
								  const double tor_zmin_kpc, const double tor_zmax_kpc, 
								  const double tor_rmin_kpc, const double tor_rmax_kpc) const 
{
	Vector3d vals_gauss(0, 0, 0);
	if ((pos_kpc.z <= tor_zmin_kpc) or (tor_zmax_kpc <= pos_kpc.z))
		return vals_gauss;
	
	double cR, theta;
	cR = sqrt(pow(pos_kpc.x, 2) + pow(pos_kpc.y, 2));
	if ((tor_rmin_kpc <= cR) and (cR <= tor_rmax_kpc))
	{
		theta = atan2(pos_kpc.y, pos_kpc.x);
		vals_gauss.x =  tor_B_gauss*sin(theta);
		vals_gauss.y = -tor_B_gauss*cos(theta);
		vals_gauss.z =  0;
	}
	return vals_gauss;
}

Vector3d KST24Field::get_Xfield(const Vector3d pos_kpc, const double Xfield_B_gauss, 
								const double Xfield_rmin_kpc, const double Xfield_rmax_kpc, 
								const double Xfield_theta_rad) const
{
	Vector3d vals_gauss(0, 0, 0);

	int sign = 1;
	if (pos_kpc.z < 0)
		sign = -1;
	if (fabs(pos_kpc.z) > 10)
		return vals_gauss;

	double cR, cR_0;
	cR = sqrt(pos_kpc.x*pos_kpc.x + pos_kpc.y*pos_kpc.y);
	cR_0 = cR - pos_kpc.z*sign*tan(Xfield_theta_rad);

	if ((cR_0 < Xfield_rmin_kpc) or (Xfield_rmax_kpc < cR_0))
		return vals_gauss;

	double mgn_frc = cR_0/cR;
	double xy_theta = atan2(pos_kpc.y, pos_kpc.x);
	vals_gauss.x =  Xfield_B_gauss*mgn_frc*cos(xy_theta)*sign*sin(Xfield_theta_rad);
	vals_gauss.y =  Xfield_B_gauss*mgn_frc*sin(xy_theta)*sign*sin(Xfield_theta_rad);
	vals_gauss.z =  Xfield_B_gauss*mgn_frc*cos(Xfield_theta_rad);

	return vals_gauss;
}

Vector3d KST24Field::get_LB(const Vector3d pos_kpc, const double LB_B_gauss, 
							const double LB_lB_deg, const double LB_bB_deg, 
							const double LB_rmin_kpc, const double LB_dr_kpc,
							const double LB_x0_kpc, const double LB_y0_kpc, const double LB_z0_kpc) const
{
	Vector3d vals_gauss(0, 0, 0);

	double eRx, eRy, eRz, cR;
	eRx = pos_kpc.x - LB_x0_kpc;
	eRy = pos_kpc.y - LB_y0_kpc;
	eRz = pos_kpc.z - LB_z0_kpc;
	cR = sqrt(eRx*eRx + eRy*eRy + eRz*eRz);
 	if ((cR < LB_rmin_kpc) or (LB_rmin_kpc + LB_dr_kpc < cR))
		return vals_gauss;
	eRx /= cR;
	eRy /= cR;
	eRz /= cR;

	double fdir_x, fdir_y, fdir_z, cosTheta;
	cosTheta = eRx*LB_Bdir_x + eRy*LB_Bdir_y + eRz*LB_Bdir_z; // expanding cross product: b(ac) - c(ab)
	fdir_x = LB_Bdir_x - eRx*cosTheta;
	fdir_y = LB_Bdir_y - eRy*cosTheta;
	fdir_z = LB_Bdir_z - eRz*cosTheta;

	// magnetic field amplification factor
	double ampl, ampl_elec;
	ampl = (1 + pow(LB_rmin_kpc, 2)/(2*LB_rmin_kpc*LB_dr_kpc + LB_dr_kpc*LB_dr_kpc));

	vals_gauss.x = LB_B_gauss*ampl*fdir_x;
	vals_gauss.y = LB_B_gauss*ampl*fdir_y;
	vals_gauss.z = LB_B_gauss*ampl*fdir_z;

	return vals_gauss;
}

Vector3d KST24Field::get_logspiral(const Vector3d pos_kpc, const double B_gauss,  
								   const double pitch_deg, const double phi0_deg,
								   const double x_shift_kpc, const double y_shift_kpc,
								   const double arc_radius1_kpc, const double arc_radius2_kpc, 
								   const double arc_eps , const double arc_div_deg, 
								   const double rmin_kpc, const double rmax_kpc, 
								   const double zmin_kpc, const double zmax_kpc) const
{
	Vector3d vals_gauss(0, 0, 0);

	if ((pos_kpc.z < zmin_kpc) or (zmax_kpc < pos_kpc.z))
		return vals_gauss;

	std::vector<double> pos_v;
	pos_v.push_back(pos_kpc.x + x_shift_kpc);
	pos_v.push_back(pos_kpc.y + y_shift_kpc);
	pos_v.push_back(pos_kpc.z);

	double r = sqrt(pos_v[0]*pos_v[0] + pos_v[1]*pos_v[1]);
	if ((r < rmin_kpc) or (rmax_kpc < r))
		return vals_gauss;

	double a_kpc = 3;
	double k = tan(pitch_deg*M_PI/180.);
	double cos_pitch = cos(pitch_deg*M_PI/180.);
	double sin_pitch = sin(pitch_deg*M_PI/180.);
	double phi0 = phi0_deg*M_PI/180.;
	double arc_div_rad = arc_div_deg*M_PI/180.;

	int nn;
	double phi = atan2(pos_v[1], pos_v[0]);	
	nn = floor((log(r/a_kpc) - k*(phi + phi0))/(2*M_PI*k));

	double r1, r2;
	r1 = a_kpc*exp(k*(phi + phi0 + 2*M_PI*nn));
	r2 = r1*exp(k*2*M_PI);

	if (fabs(r - r1) > fabs(r - r2))
	{
		r1 = r2;
		nn++;  // we need outer arm
	}


	// perpendicular to the spiral arm axis
	double xi, yi, dir1, dir2, dirtmp;
	xi = r1*cos(phi);
	yi = r1*sin(phi);
	dir1 = k*xi - yi;
	dir2 = k*yi + xi;
	dirtmp = sqrt(dir1*dir1 + dir2*dir2);
	dir1 /= dirtmp;
	dir2 /= dirtmp;
	// std::cout << ((xi - pos[0])*dir1 + (yi - pos[1])*dir2)/sqrt(pow((xi - pos[0]), 2) + pow((yi - pos[1]), 2)) << '\n';


	// first order correction 
	double delta_phi;
	delta_phi = -((xi - pos_v[0])*dir1 + (yi - pos_v[1])*dir2)/(r1/cos_pitch);
	r1 = a_kpc*exp(k*(phi + phi0 + delta_phi + 2*M_PI*nn));
	xi = r1*cos(phi + delta_phi);
	yi = r1*sin(phi + delta_phi);
	dir1 = k*xi - yi;
	dir2 = k*yi + xi;
	dirtmp = sqrt(dir1*dir1 + dir2*dir2);
	dir1 /= dirtmp;
	dir2 /= dirtmp;


	double d1, L1, d_at_L1;
	d1 = sqrt(pow((xi - pos_v[0]), 2) + pow((yi - pos_v[1]), 2));
	L1 = r1/sin_pitch;

	double r_scale = 5;
	double L1_at_r_scale, arc_r_plane, arc_r_z;
	L1_at_r_scale = r_scale/sin_pitch;
	arc_r_plane = std::max(arc_radius2_kpc, arc_radius2_kpc*(1 + (L1 - L1_at_r_scale)*arc_div_rad));
	// arc_r_plane = arc_radius2_kpc;
	arc_r_z     = std::max(arc_radius1_kpc, arc_radius1_kpc*(1 + (L1 - L1_at_r_scale)*arc_div_rad));
	if (arc_r_plane > 1.2)
		arc_r_plane = 1.2;
	if (arc_r_z > 1.2)
		arc_r_z = 1.2;


	d_at_L1 = sqrt(pow(fabs(d1/arc_r_plane), arc_eps) + pow(fabs(pos_v[2]/arc_r_z), arc_eps));
	if ((d_at_L1 > 1) or (d_at_L1 < 0.))
		return vals_gauss;
	else
	{
		double scaling_factor = (arc_radius1_kpc*arc_radius2_kpc)/(arc_r_plane*arc_r_z);
		vals_gauss.x = B_gauss*dir1*scaling_factor;
		vals_gauss.y = B_gauss*dir2*scaling_factor;
		vals_gauss.z = 0;
	}

	return vals_gauss;
}
}