#include "crpropa/magneticField/CMZField.h"
#include "crpropa/Units.h"

namespace crpropa {
//                                     Magnetic Field Model in the Galactic Center from M. Guenduez and J. B. Tjus (2019)                                    
// Poloidal Model (Model C) is taken from Katia Ferriere and Philippe Terral 2014 "Analytical models of X-shape magnetic fields in galactic halos"
// Azimuthal Model is taken from M.Guenduez, J.B. Tjus, K.Ferriere, R.-J. Dettmar (2019)
//                                                                                    
CMZField::CMZField() {
    useMCField = false;
    useICField = true;
    useNTFField = false;
    useRadioArc = false;
}

double CMZField::scale(double x, double d) const{
    return x*1./(pow(log(2),(1./d)));    
}

double CMZField::Br(double r, double z, double B1, double a, double L) const{//poloidal field
    double r1=1/(1+a*pow(z,2.));
    return 2*a*r*pow(r1,3.) *z *B1*exp(-r/L *r1);
}

double CMZField::Bz(double r, double z, double B1, double a, double L) const{//poloidal field
    double r1=1/(1+a*pow(z,2.));
    return pow(r1,2.) *B1*exp(-r/L *r1);
}

double CMZField::BrAz(double r, double phi, double z, double m, double B1, double eta, double R, double rr) const{//azimuthal field
	double A;
	if (r<rr)
	{
		A=  R/rr *(3.*r/rr -2.* pow(r/rr,2.));
	}
	else
	{
		A= R/r;
	}
	return B1 * exp(-pow(z/scale(R,2.),2.)) * cos(m/eta *log((r+1.)/(R+1.)) +m*phi) * A;
}

double CMZField::BphiAz(double r, double phi, double z, double m, double B1, double eta, double R, double rr) const{//azimuthal field
	double A;
	if (r<rr)
	{
		A=(1.+6.*(r-rr)/(2.*r-3.*rr) * ( tan(m/eta *log((r+1.)/(R+1.)) +m*phi) - sin(m/eta *log((r+1.)/(R+1.)))/cos(m/eta *log((r+1.)/(R+1.)) +m*phi) )) ;
	}
	else
	{
    	A= 1.;
	}
	return -1/eta *r/(r+1.) * BrAz(r,phi,z,m,B1,eta,R,rr) * A;
}

// azimuthal field in cartesian  coordinates// for molecular cloud region
double CMZField::BxAz(double x, double y, double z, double m, double B1, double eta, double R, double rr) const{
    double r= sqrt(pow(x,2.)+pow(y,2.));
    double phi= std::atan2(y,x);
    return (cos(phi)-1./eta * sin(phi)) *BrAz(r,phi,z,m,B1,eta,R,rr);
}

double CMZField::ByAz(double x, double y, double z, double m, double B1, double eta, double R, double rr) const{
    double r= sqrt(pow(x,2.)+pow(y,2.));
    double phi= std::atan2(y,x);
    return (sin(phi)+1./eta*cos(phi))*BrAz(r,phi,z,m,B1,eta,R,rr);
}

double CMZField::BzAz(double x, double y, double z, double m, double B1, double eta, double R,double rr) const{
    double r= sqrt(pow(x,2.)+pow(y,2.));
    double phi= std::atan2(y,x);
    return 0;//z component of the magnetic field is in molecular clouds zero
}

// Polodial field in cartesian coordinates// for non-thermal filament and the intercloud region
double CMZField::ByPol(double x, double y, double z, double B1, double a1, double a2, double eta) const{
    double r= sqrt(pow(x,2.)+pow(y,2.));
    double phi= std::atan2(y,x);
    double B2=B1/eta;// B1 is the observed magnetic field stremgth and eta the normalization factor otained from the integrations due to the normalization
    double L=scale(a2/2.,1);
    double a= 1/pow(scale(a1/2.,2),2.);
    return sin(phi)*Br(r,z,B2,a,L);
}

double CMZField::BxPol(double x, double y, double z, double B1, double a1, double a2, double eta) const{
    double r= sqrt(pow(x,2.)+pow(y,2.));
    double phi= std::atan2(y,x);
    double B2=B1/eta;
    double L=scale(a2/2.,1);
    double a= 1/pow(scale(a1/2.,2),2.);
    return cos(phi)*Br(r,z,B2,a,L);
}

double CMZField::BzPol(double x, double y, double z, double B1, double a1, double a2, double eta) const{
    double r= sqrt(pow(x,2.)+pow(y,2.));
    double phi= std::atan2(y,x);
    double B2=B1/eta;
    double L=scale(a2/2.,1);
    double a= 1/pow(scale(a1/2.,2),2.);
    return Bz(r,z,B2,a,L);
}
/////end model C + azimuthal component

bool CMZField::getUseMCField() const {
    return useMCField;
}
bool CMZField::getUseICField() const {
    return useICField;
}
bool CMZField::getUseNTFField() const {
    return useNTFField;
}
bool CMZField::getUseRadioArc() const {
    return useRadioArc;
}

void CMZField::setUseMCField(bool use) {
    useMCField = use;
}
void CMZField::setUseICField(bool use) {
    useICField = use;
}
void CMZField::setUseNTFField(bool use) {
    useNTFField = use;
}
void CMZField::setUseRaidoArc(bool use) {
    useRadioArc = use;
}

// Parameter identification:
// r: radius in zylindrical coordinates
// z= z component in zylindircal coordintes
// phi: phi component in zylindircal coordinates
// B1: observed magnetic field strength
// B2: normalization factor according to B1
// m: azimuthal wavenumber m=0 for axissymmetric, m=1 for bisymmetric and m=2 quadrosymmetric, ... field lines
// a: stricly positive free parameters governing the opening of field lines away from the z-axis
// L: exponential scale length length of the cloud radius Longitudinal extent
// Lp:scale length of the winding function set to a fraction of the cloud radius  Longitudes of peaks mainly at high |l|
// H: exponential scale high length of the cloud radius Latitudinal extent
// Hp: exponential scale high  set to a fraction of the cloud vertical size Longitude of peaks mainly at hight |b|
// p: pitch angle at the origin, i.e. the angle between the horizontal projection of a field line and the local azimuthal direction  (i.e., 10° or 20°) if you want the field to be tightly wound up, and hence predominantly azimuthal
// phi_s: orientation angle of the azimuthal pattern if g_s(r1(r,z,a),0,a,Lp,Hp,p), azimuthal angle at infinity of the crest surface if g_s(r,z,a,Lp,Hp,p)                                                                                  

Vector3d CMZField::getMCField(const Vector3d& pos) const {//Field in molecular clouds
	double pi=M_PI;
    double m=1;
    Vector3d b(0.);
    double eta=0.01;
    double N=59;
	// azimuthal component in dense clouds
    //A=SgrC 
    double y_mid=-81.59*pc;//  y midpoint of A
    double z_mid=-16.32*pc;//  z midpoint of A
    double xx=pos.x;
    double yy=pos.y-y_mid;// shifted coordinates
    double zz=pos.z-z_mid;// shiftes coordinates
    double R = 1.7*pc;// Radius of A
    double B1=2.1e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    // A=G0.253+0.016 Dust Ridge A
    y_mid=37.53*pc;//  y midpoint of A
    z_mid=-2.37*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=2.4*pc;// Radius of A
    B1=2.5e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //A=Dust Ridge B
    y_mid=50.44*pc;//  y midpoint of A
    z_mid=8.16*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=1.9*pc;// Radius of A
    B1=0.9e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //A=Dust Ridge C
    y_mid=56.37*pc;//  y midpoint of A
    z_mid=7.71*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=1.9*pc;// Radius of A
    B1=1.2e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //A=Dust Ridge D
    y_mid=60.82*pc;//  y midpoint of A
    z_mid=7.42*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=3.3*pc;// Radius of A
    B1=1.7e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //A=Dust Ridge E
    y_mid=70.91*pc;// y midpoint of A
    z_mid=0.74*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=3.5*pc;// Radius of A
    B1=4.1e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //A=Dust Ridge F
    y_mid=73.58*pc;//  y midpoint of A
    z_mid=2.97*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=2.4*pc;// Radius of A
    B1=3.9e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //Sgr D
    y_mid=166.14*pc;//  y midpoint of A
    z_mid=-10.38*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=1.8*pc;//
    B1=0.8e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);
    //Sgr B2
    y_mid=97.91*pc;//  y midpoint of A
    z_mid=-5.93;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=14*pc;// Radius of A
    B1=1.0e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.);                                                                                  
    //A=Inner R=5pc 
    y_mid=-8.3*pc;//  y midpoint of A
    z_mid=-6.9*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=5*pc;//2.*10*pc;// Radius of A
    B1=3.0e-3/0.91;
    b.x +=BxAz(xx,yy,zz,m,B1,0.77,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,0.77,R,R/10.);
    //20 km s^-1
    y_mid=-19.29*pc;//  y midpoint of A
    z_mid=-11.87*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=9.4*pc;// Radius of A
    B1=2.7e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.); 
    //50 km s^-1
    y_mid=-2.97*pc;//  y midpoint of A
    z_mid=-10.38*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=9.4*pc;// Radius of A
    B1=3.7e-3/N;
    b.x +=BxAz(xx,yy,zz,m,B1,eta,R,R/10.);
    b.y +=ByAz(xx,yy,zz,m,B1,eta,R,R/10.); 
    //SgrA*
    y_mid=-8.3*pc;//  y midpoint of A
    z_mid=-6.9*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    R=1.2e12;// Radius of A
    B1=65./3.07;
    double rr= sqrt(pow(xx,2.)+pow(yy,2.));
    double pphi= std::atan2(yy,xx);
    b.x +=-sin(pphi)*(-BrAz(rr,pphi,zz,0,B1,1.,R,R/10.));
    b.y +=cos(pphi)*(-BrAz(rr,pphi,zz,0,B1,1.,R,R/10.)); 
    //b.z+=Bz2(rr,pphi,zz,B2,m,a2,L2,Lp,Hp,p,phi_s);
    b.x=b.x*1.e-4;
    b.y=b.y*1.e-4;
    b.z=b.z*1.e-4;
    return b;//1.e-4 due to the conversion from Gauss to Tesla
} 
Vector3d CMZField::getICField(const Vector3d& pos) const {//Field in intercloud medium--> poloidal field
    Vector3d b(0.);
    double pi= M_PI;
    //poloidal field in the intercloud region
    double x=pos.x;
    double y=pos.y -(-8.3*pc);
    double z=pos.z -(-6.9*pc);
    double a2=2.*158.*pc;// Radius of A
    double a1=2.*35.*pc;
    double eta=0.85;
    double B1=1e-5*gauss;
    double B2 = B1/eta;

    double r = sqrt(x*x+y*y);
    double phi = std::atan2(y,x);
    double sPhi = sin(phi);
    double cPhi = cos(phi);
    double a = 4*log(2)/a1/a1;
    double L = a2/log(2)/2.;
    double r1 = 1/(1+a*z*z);
    double Br =2*a*r*pow(r1,3)*z*B2*exp(-r/L*r1);

    b.y += sPhi*Br;
    b.x += cPhi*Br;
    b.z += r1*r1 * B2 * exp(-r/L *r1);
    return b;
}                                                         
    //  
Vector3d CMZField::getNTFField(const Vector3d& pos) const {//Field in the non-thermal filaments--> predominantly poloidal field
    Vector3d b(0.);
    double pi= M_PI;
    //poloidal field in the non-thermal filament region (except "pelical"-> azimuthal)
    double arcmin=1./60.;

    //A=SgrC
    double y_mid=-81.59*pc;//  y midpoint of A
    double z_mid=-1.48*pc;//  z midpoint of A
    double xx=pos.x;
    double yy=pos.y-y_mid;// shifted coordinates
    double zz=pos.z-z_mid;// shiftes coordinates
    double a1=27.44*pc;// arcmin-> deg->cm
    double a2=1.73*pc;// arcmin-> deg-> cm
    double eta=0.48;
    double B1=1.e-4;
    b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G359.15-0.2 The Snake
    y_mid=-126.1*pc;//  y midpoint of A
    z_mid=-25.22*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    a1=12.86*pc;// arcmin-> deg->cm
    a2=2.22*pc;// arcmin-> deg-> cm
    B1=88.e-6;
    b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G359.54+0.18 Nonthermal Filament
    y_mid=-68.24*pc;//  y midpoint of A
    z_mid=25.22*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    a1=15.08*pc;// arcmin-> deg->cm
    a2=2.72;// arcmin-> deg-> cm
    B1=1.e-3;
    b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G359.79 +17 Nonthermal Filament
    y_mid=-31.15*pc;//  y midpoint of A
    z_mid=23.74*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    a1=16.07*pc;// arcmin-> deg->cm
    a2=3.46*pc;// arcmin-> deg-> cm
    B1=1.e-3;
    b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G359.91 -1.03 Possible Nonthermal Filament 
    // to use this filament delete comments below
    //x_mid=0;//  x midpoint of A
    //y_mid=sin(359.91*pi/180)*8.5*kpc;//  y midpoint of A
    //z_mid=sin(-1.03*pi/180)*8.5*kpc;//  z midpoint of A
    //xx=x-x_mid;// shifted coordinates
    //yy=y-y_mid;// shifted coordinates
    //zz=z-z_mid;// shiftes coordinates
    //a1=sin(2.3*arcmin*pi/180)*8.5*kpc;// arcmin-> deg->cm
    //a2=sin(0.6*arcmin*pi/180)*8.5*kpc;// arcmin-> deg-> cm
    //B1=1.e-3;
    //b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    //b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    //b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G359.96 +0.09  Nonthermal Filament Southern Thread
    y_mid=-5.93*pc;//  y midpoint of A
    z_mid=16.32*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    a1=28.68*pc;// arcmin-> deg->cm
    a2=1.73*pc;// arcmin-> deg-> cm
    B1=1.e-4;
    b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G0.09 +0.17  Nonthermal Filament Northern thread
    y_mid=13.35*pc;//  y midpoint of A
    z_mid=25.22*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    a1=29.42*pc;// arcmin-> deg->cm
    a2=2.23*pc;// arcmin-> deg-> cm
    B1=140.e-6;
    b.y+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=BzPol(xx,yy,zz,B1,a1,a2,eta);

    //A=G359.85+0.47  Nonthermal Filament The Pelican  is not poloidal but azimuthal
    y_mid=-22.25*pc;//  y midpoint of A
    z_mid=69.73*pc;//  z midpoint of A
    yy=pos.y-y_mid;// shifted coordinates
    zz=pos.z-z_mid;// shiftes coordinates
    a1=11.37*pc;// arcmin-> deg->cm
    a2=2.23*pc;// arcmin-> deg-> cm
    B1=70.e-6 /eta;// /by bz switched because pelican is differently oriented
    b.y+=BzPol(xx,yy,zz,B1,a1,a2,eta);
    b.x+=BxPol(xx,yy,zz,B1,a1,a2,eta);
    b.z+=ByPol(xx,yy,zz,B1,a1,a2,eta);
    
    b.x=b.x*1.e-4;
    b.y=b.y*1.e-4;
    b.z=b.z*1.e-4;
	return b;//1.e-4 due to the conversion from Gauss to Tesla
}
  
Vector3d CMZField::getRadioArcField(const Vector3d& pos) const {//Field in the non-thermal filaments--> predominantly poloidal field
    Vector3d b(0.);
    //poloidal field in the non-thermal filament region A=RadioArc
    double pi= M_PI;
    double arcmin=1./60.;
    double eta=0.48;
    double y_mid=26.7*pc;//  y midpoint of A
    double z_mid=10.38*pc;//  z midpoint of A
    double xx=pos.x;// shifted coordinates
    double yy=pos.y-y_mid;// shifted coordinates
    double zz=pos.z-z_mid;// shiftes coordinates
    double a1=70.47*pc;// arcmin-> deg->cm
    double a2=9.89*pc;// arcmin-> deg-> cm
    double B1=1.e-3;
    b.y=ByPol(xx,yy,zz,B1,a1,a2,eta)*1.e-4;//1.e-4 due to the conversion from Gauss to Tesla
    b.x=BxPol(xx,yy,zz,B1,a1,a2,eta)*1.e-4;//1.e-4 due to the conversion from Gauss to Tesla
    b.z=BzPol(xx,yy,zz,B1,a1,a2,eta)*1.e-4;//1.e-4 due to the conversion from Gauss to Tesla
    return b;
}

Vector3d CMZField::getField(const Vector3d& pos) const{
    Vector3d b(0.);

    if(useMCField){
        b += getMCField(pos);
    }
    if(useICField){
        b += getICField(pos);
    }
    if(useNTFField){
        b += getNTFField(pos);
    }
    if(useRadioArc){
        b += getRadioArcField(pos);
    }

    return b;
}

} // namespace crpropa
