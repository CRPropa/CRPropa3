#include "crpropa/module/HadronicInteraction.h"
#include "crpropa/Massdistribution/Ferrie07.h"
#include "crpropa/Units.h"
#include "crpropa/ParticleID.h"
#include "crpropa/ParticleMass.h"
#include "crpropa/Random.h"


#include <fstream>
#include <limits>
#include <cmath>
#include <stdexcept>

namespace crpropa {
    
    HadronicInteraction::HadronicInteraction() {
        
        setDescription("HadronicInteraction");
    }
    
    Vector3d HadronicInteraction::Position(double height, double radius) const{
        Random &random = Random::instance();
        int i=0;
        Vector3d pos(0,0,0);
        double phi=random.rand()*2*M_PI;
        int j=0;
        do{
            double r = random.rand()* radius;
            double yr = random.rand();
            double Fr=exp(-r*r/(2*4.2*4.2*kpc*kpc));
            
            if (yr < Fr){
                pos=Vector3d(cos(phi)*r, sin(phi)*r, 0);
                j++;
            }
        }while(j==0);
        do{
            double z = random.rand()* height;
            double yz = random.rand();
            double Fz=exp(-z/(10*pc));
            if (yz < Fz){
                double a = random.rand();
                if (a <= 0.5){
                    z=-z;
                }
                pos= pos + Vector3d(0,0,z);
                j++;
                
            }
        }while(j==1);
        return pos;
    }
    
    //Energy distribution for pions
    
    double HadronicInteraction::distribution_pi(double energy, double x) const{
        double L=log(energy/TeV);
        double Bp=5.58+0.78*L+0.1*L*L;
        double r=3.1/pow(Bp, 1.5);
        double a=0.89/(pow(Bp, 0.5)*(1-exp(-0.33*Bp)));
        double mp=134.9766*MeV;
        double A=4*a*Bp*pow(x, a-1);
        double B=(1-pow(x,a))/(1+r*pow(x,a)*(1-pow(x,a)));
        double C=(1/(1-pow(x,a))+(r*(1-2*pow(x,a))/(1+r*pow(x,a)*(1-pow(x,a)))))*pow(1-mp/(x*energy), 0.5);
        double F=A*pow(B, 4.)*C;
        
        return F;
    }
    
    double HadronicInteraction::number_pi(double energy) const{
        double x=1/1000.;
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_pi(energy,x);
            x=x+i;
            j++;
        }while(x < 1);
        return y/j*(x-1/1000.);
    }
    
    //Energy distribution for electron, electron neutrino and second muon neutrino based on Kelner 2006
    
    double HadronicInteraction::distribution_e(double energy, double x) const{
        
        double L=log(energy / TeV);
        double Be= 1/(69.5+2.65*L+0.3*pow(L,2.));
        double betae=1/pow((0.201+0.062*L+0.00042*pow(L,2.)), 0.25);
        double ke=(0.279 + 0.141 *L + 0.0172* pow(L, 2.))/(0.3+ pow((2.3+L), 2.));
        double F=Be*pow((1+ke*pow(log(x),2.)), 3.) /(x*(1+0.3/pow(x, betae)))*(pow(-log(x), 5.));
        
        return F;
    }
    
    double HadronicInteraction::number_e(double energy) const{
        double x=1/100000.;
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_e(energy,x);
            x=x+i;
            j++;
        }while(x < 1);
        return y/j*(x-1/1000.);
    }
    
    
    
    //Energy distribution for first muon neutrino based on Kelner 2006
    double HadronicInteraction::distribution_my1(double energy, double x) const{
        double L=log(energy / TeV);
        double Bm= 1.75+0.204*L+0.01 * pow(L,2.);
        double betam=1/(1.67+0.111*L+0.0038*pow(L,2.));
        double km=1.07-0.086*L+0.002*pow(L,2.);
        x=x/0.427;
        double aa=(1-pow(x,betam))/(1+km*pow(x, betam)*(1-pow(x,betam)));
        double A=Bm*log(x)/x*pow(aa, 4.);
        double B=1/log(x)-4*betam*pow(x,betam)/(1- pow(x,betam))-4*km*betam*pow(x, betam)*(1-2*pow(x,betam))/(1+km*pow(x,betam)*(1-pow(x,betam)));
        double F=A*B;
        
        return F;
    }
    
    double HadronicInteraction::number_my1(double energy) const{
        double x=1/(100000.);
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_my1(energy,x);
            x=x+i;
            j++;
        }while(x < 0.427);
        return y/j*(x-1/1000.);
    }
    //Energy distribution for gamma photons based on Kelner 2006
    double HadronicInteraction::distribution_gamma(double energy, double x) const{
        double L=log(energy / TeV);
        double Bg=1.3+0.14*L+0.011*L*L;
        double betag=1/(1.79+0.11*L+0.008*L*L);
        double kg=1/(0.801+0.049 *L+0.014 *L*L);
        double A=Bg*log(x)/x;
        double B=(1-pow(x, betag))/(1+kg*pow(x, betag)*(1-pow(x, betag)));
        double C=1/log(x)-4*betag*pow(x, betag)/(1-pow(x, betag))-4*kg*betag*pow(x, betag)*(1-2*pow(x, betag))/(1+kg*pow(x, betag)*(1-pow(x, betag)));
        double F=A*pow(B, 4.)*C;
        
        return F;
    }
    
    double HadronicInteraction::number_gamma(double energy) const{
        double x=1/100000.;
        double i=1/100000.;
        double y=0;
        double j=0;
        do{
            y=y+distribution_gamma(energy,x);
            x=x+i;
            j++;
        }while(x < 1);
        return y/j*(x-1/1000.);
    }
    //Energy distribution for lepton secondaries of pp interactions based on Carceller 2017
    
    double HadronicInteraction::distribution_Carceller(double energy, double x, double jcap, double a0, double b0) const{
        double a = a0 * (1+ 0.073*log(energy/ PeV)+ 0.0070*log(energy/PeV)*log(energy/PeV));
        double b = b0 * (1+0.020*log(energy/PeV)+0.0018*log(energy/PeV)*log(energy/PeV));
        double A = a * pow((1-jcap*x), 3.)/x;
        double B = exp(-b*pow(jcap*x, 0.43))/pow(1+pow(0.1*GeV/(x*energy), 0.5), 2.);
        double F=A*B;
        
        return F;
    }
    
    
    //Energy distribution for gamma photons based on Carceller 2017
    
    double HadronicInteraction::distribution_Carceller_g(double energy, double x, double jcap, double a0, double b0) const{
        double a = a0 * (1+ 0.073*log(energy/ PeV)+ 0.0070*log(energy/PeV)*log(energy/PeV));
        double b = b0 * (1+0.02*log(energy/PeV)+0.0018*log(energy/PeV)*log(energy/PeV));
        double A = a * pow((1-jcap*x), 3.)/x;
        double B = exp(-b*pow(jcap*x, 0.43))/pow(1+pow(0.2*GeV/(x*energy), 0.5), 2.);
        double F=A*B;
        
        return F;
    }
    
    //Cross Section of inelastic pp interaction based on Tan & Ng 1983 (Used in Galprop)
    double HadronicInteraction::CrossSection_Galprop(double energy) const{
        double cs_inel;
        double U = log(energy/ GeV * 1/200);
        if (U >= 0 and energy >= 3 * GeV){
            cs_inel=(32.2 * (1+0.0273*U))*1e-31+32.2*0.01*pow(U,2.)*1e-31;
        }
        if (U < 0 and energy >= 3 * GeV){
            cs_inel=(32.2 * (1+0.0273*U))*1e-31;
        }
        if (energy <= 0.3 * GeV){
            cs_inel = 0;
        }
        return cs_inel;
    }
    
    //Cross Section of inelastic pp interaction based on Kelner 2006
    double HadronicInteraction::CrossSection_Kelner(double energy) const{
        double L=log(energy / TeV);
        double A=1-pow(1.22*1e-3*TeV/energy, 4.);
        double cs_inel=(34.3 + 1.88*L+0.25 *L*L)*A*A*1e-31;
        return cs_inel;
    }
    
    //Cross Section of inelastic pp interaction based on Carceller 2017
    double HadronicInteraction::CrossSection_Carceller(double energy) const{
        double cs_inel=17.7*pow(energy/GeV, 0.082)*1e-31;
        return cs_inel;
    }
    
    //Process Function:
    
    void HadronicInteraction::process(Candidate *candidate) const {
        
        double step = candidate->getCurrentStep();
        double energy = candidate->current.getEnergy();
        double id = candidate->current.getId();
        //~ if (energy < 10 * TeV){
        //~ return;
        //~ }
        
        if (id < 500){
            return;
        }
        
        
        double cs_inel=0;
        double eG=energy / GeV;
        double jcap=1;
        
        
//        if (id == 1000010010) {
//            cs_inel=CrossSection_Galprop(energy);
//        }
        
        if (id == 1000010010) {
            cs_inel=CrossSection_Kelner(energy);
        }
        
        
        //~ if (id == 1000010010){
        //~ cs_inel=CrossSection_Carceller(energy);
        //~ }
        
        
        if (id == 1000020040 and energy > 1* TeV) {
            double cs_Hep=60.5*1e-31;
            cs_inel= cs_Hep*pow(eG, 0.062);
            jcap=4;
        }
        
        
        if (id == 1000260056 and energy > 1* TeV) {
            double cs_Fep=551*1e-31;
            cs_inel= cs_Fep*pow(eG, 0.026);
            jcap=56;
        }
        
        //~ Cross section for different target material. Calculations based on R. Silberberg and C. H. Tsao
        //~ else {
        //~ cs_inel=45 * pow(id, .7)*(1+0.016 * sin(5.3-2.63*log(id)))*1e-31 ;
        //~ cs_inel=cs_inel*(1-0.62* pow(exp(1), -energy/200)*sin(10.9 * pow(10.9*energy, -0.28)));
        //~ }
        
        
        Random &random = Random::instance();
        Vector3d pos= candidate->current.getPosition();
        Ferrie obj;
        
        double density=obj.getDensity(pos)*1000000;
        double p_pp=cs_inel*density*step;
        //std::cout << p_pp << std::endl;
        double ra = random.rand();
        
        if (ra > p_pp or energy < 1*GeV){
            
            return;
            
        }
        
        
        double limit = 1 / p_pp*0.1;
        std::cout << "Int" << std::endl;
        if (step > limit) {
            // limit next step to mean free path
            candidate->limitNextStep(limit);
            
        }
        double Eout=0;
        
        
        double Emo=0;
        double Ee=0;
        double Ene=0;
        double Emt=0;
        double Eg=0;
        double gamma=0;
        double j=0;

        //~ goto label2;
        if (jcap == 1)
        {
            goto label;
        }
          label:
        // Establish number of secondaries
        //Gammas
        double r=random.rand();
        double FG=number_gamma(energy);
        //std::cout<<"Fg="<<FG<<std::endl;
        double NG=std::floor(FG);
        if (r < FG-NG){
            NG=NG+1;
        }
        
        //First myon neutrino
        r=random.rand();
        double Fmy1=number_my1(energy);
//      std::cout<<"Fmy1="<<Fmy1<<std::endl;
        double Nmy1=std::floor(Fmy1);
        if (r < Fmy1-Nmy1){
            Nmy1=Nmy1+1;
        }
        
        //Electron
        r=random.rand();
        double FE=number_e(energy);
//        std::cout<<"FE="<<FE<<std::endl;
        double NE=std::floor(FE);
        if (r < FE-NE){
            NE=NE+1;
        }
        
        //Electron Neutrino
        r=random.rand();
        double NEN=std::floor(FE);
        
        if (r < FE-NEN){
            NEN=NEN+1;
        }
        
        //Second Myon Neutrino
        r=random.rand();
        double Nmy2=std::floor(FE);
        if (r < FE-Nmy2){
            Nmy2=Nmy2+1;
        }
        
        double N_tot=NG+Nmy1+NE+NEN+Nmy2;
        //std::cout<<N_tot<<std::endl;
        
  
        
        
        
        double Econ=0;
        
        double i=1;
        double iG=1;
        double imy1=1;
        double ie=1;
        double ien=1;
        double imy2=1;
        double test;
        
        do{
            //~ Gamma rays
            test=iG;
            j=0;
            if (iG <= NG){
                do{
                double x=random.rand()*(-5);
                    
                    j++;
                Eout=pow(10, x)*energy;
                double E=distribution_gamma(energy, pow(10, x));
                double Emax=distribution_gamma(energy, 0.00001);
                double y=random.rand()*Emax;
//                    if(j == 10000){
//                        x=-5;
//                        Eout=pow(10, x)*energy;
//                        candidate->addSecondary(22, Eout, pos);
//                        Eg=Eg+Eout;
//                        i++;
//                        iG++;
//                        //std::cout<<(Ee+Ene+Emt+Emo+Eg)/energy<<std::endl;
//                        std::cout<<"GammaE"<<std::endl;
//                    }
                if (y < E and (Ee+Ene+Emt+Emo+Eg+Eout)<energy ){
                    candidate->addSecondary(22, Eout, pos);
                    Eg=Eg+Eout;
                    i++;
                    iG++;
                }
                }while(test == iG);
            }

            
            //~ First myon neutrino 14
            test = imy1;
            j=0;
            if (imy1 <= Nmy1){
                do{
                label3:
                double x=-5+random.rand()*(log(0.427)+5);
                double my=pow(10,x);
                    j++;
                Eout=my*energy;
                double E=distribution_my1(energy, my);
                double Emax=distribution_my1(energy, 0.00001);
                double y=random.rand()*Emax;
//                if(j == 10000){
//                    x=-5;
//                    Eout=pow(10, x)*energy;
//                    candidate->addSecondary(14, Eout, pos);
//                    Emo=Emo+Eout;
//                    i++;
//                    imy1++;
//                    //std::cout<<(Ee+Ene+Emt+Emo+Eg)/energy<<std::endl;
//                    std::cout<<"My1"<<std::endl;
//                    }
                if (y < E and (Ee+Ene+Emt+Emo+Eg+Eout)<energy){
                    candidate->addSecondary(12, Eout, pos);
                    Emo=Emo+Eout;
                    //std::cout << Econ/energy << std::endl;
                    i++;
                    imy1++;
                }
            }while(test == imy1);
            }
  
            //~ Electron 11

        test=ie;
        j=0;
            
            if (ie <= NE){
                do{
                double x=random.rand()*(-5);
                    j++;
                Eout=pow(10, x)*energy;
                double E=distribution_e(energy, pow(10, x));
                double Emax=distribution_e(energy, 0.00001);
                double y=random.rand()*Emax;
//                    if(j == 10000){
//                        x=-5;
//                        Eout=pow(10, x)*energy;
//                        candidate->addSecondary(11, Eout, pos);
//                        Ee=Ee+Eout;
//                        i++;
//                        ie++;
//                        //std::cout<<(Ee+Ene+Emt+Emo+Eg)/energy<<std::endl;
//                        std::cout<<"eE"<<std::endl;
//                                    }
                if (y < E and (Ee+Ene+Emt+Emo+Eg+Eout)<energy ){
                    candidate->addSecondary(11, Eout, pos);
                    Ee=Ee+Eout;
                    //std::cout << Econ/energy << std::endl;
                    i++;
                    ie++;
                }
                }while(test == ie);
            }
            
            //~ Electron neutrino 12
            test=ien;
        j=0;
            if (ien <= NEN){
                do{
            double x=random.rand()*(-5);
                    j++;
            Eout=pow(10, x)*energy;
            double E=distribution_e(energy, pow(10, x));
            double Emax=distribution_e(energy, 0.00001);
            double y=random.rand()*Emax;
//                    if(j == 10000){
//                        x=-5;
//                        Eout=pow(10, x)*energy;
//                        y=distribution_gamma(energy, pow(10, x));
////                        candidate->addSecondary(12, Eout, pos);
//                        Ene=Ene+Eout;
//                        i++;
//                        ien++;
//                        //std::cout<<(Ee+Ene+Emt+Emo+Eg)/energy<<std::endl;
//                        std::cout<<"neE"<<std::endl;
//                        }
            if (y < E and (Ee+Ene+Emt+Emo+Eg+Eout)<energy){
//                candidate->addSecondary(12, Eout, pos);
                Ene=Ene+Eout;
                //std::cout << Econ/energy << std::endl;
                i++;
                ien++;
            }
                }while(ien==test);
            }

            //~ Second myon neutrino 14
        j=0;
            test=imy2;
            if(imy2 <= Nmy2){
                do{
            double x=random.rand()*(-5);
                    j++;
            Eout=pow(10, x)*energy;
            double E=distribution_e(energy, pow(10, x));
            double Emax=distribution_e(energy, 0.00001);
            double y=random.rand()*Emax;
//                    if(j == 10000){
//                        x=-5;
//                        Eout=pow(10, x)*energy;
//                        candidate->addSecondary(14, Eout, pos);
//                        Emt=Emt+Eout;
//                        i++;
//                        imy2++;
//                        //std::cout<<(Ee+Ene+Emt+Emo+Eg)/energy<<std::endl;
//                        std::cout<<"my2E"<<std::endl;
//                        }
            if (y < E and (Ee+Ene+Emt+Emo+Eg+Eout)<energy){
                candidate->addSecondary(14, Eout, pos);
                Emt=Emt+Eout;
                //std::cout << 's' << Econ/energy << std::endl;
                i++;
                imy2++;
            }
                }while(imy2==test);
            }
        }while (i <= N_tot);
        //std::cout<<"N="<< N_tot << std::endl;
 
        if (Ee+Ene+Emt+Emo+Eg> energy)
        {
            std::cout << (Ee+Ene+Emt+Emo+Eg)/energy << std::endl;
        }
        //std::cout << (Ee+Ene+Emt+Emo+Eg)/energy << std::endl;
        candidate->current.setEnergy(energy-(Ee+Ene+Emt+Emo+Eg));
        
        return;
        
        
    label2:
        test = 1;
        
        //gamma
        do {
            double x=random.rand()*(-3);
            
            
            double F=distribution_Carceller_g(energy, pow(10, x), jcap, 5.8, 5.3);
            double Fmax=distribution_Carceller_g(energy, 0.001, jcap, 5.8, 5.3);
            double y=random.rand()*Fmax;
            double Eout=pow(10, x)*energy;
            if (y < F ){
                candidate->addSecondary(22, Eout, pos);
                test=2;
            }
        } while (test == 1);
        
        //electron neutrino
        
        do {
            double x=random.rand()*(-3);
            
            
            double F=distribution_Carceller(energy, pow(10, x), jcap, 2.7, 7.7);
            double Fmax=distribution_Carceller(energy, 0.001, jcap, 2.7, 7.7);
            double y=random.rand()*Fmax;
            double Eout=pow(10, x)*energy;
            if (y < F ){
                candidate->addSecondary(12, Eout, pos);
                test=3;
            }
        } while (test == 2);
        
        //antielectron neutrino
        
        do {
            double x=random.rand()*(-3);
            
            
            double F=distribution_Carceller(energy, pow(10, x), jcap, 2.7, 8.7);
            double Fmax=distribution_Carceller(energy, 0.001, jcap, 2.7, 8.7);
            double y=random.rand()*Fmax;
            double Eout=pow(10, x)*energy;
            if (y < F ){
                candidate->addSecondary(-12, Eout, pos);
                test=4;
            }
        } while (test == 3);
        
        //muon neutrino
        
        do {
            double x=random.rand()*(-3);
            
            
            double F=distribution_Carceller(energy, pow(10, x), jcap, 5.1, 8.3);
            double Fmax=distribution_Carceller(energy, 0.001, jcap, 5.1, 8.3);
            double y=random.rand()*Fmax;
            double Eout=pow(10, x)*energy;
            if (y < F ){
                candidate->addSecondary(14, Eout, pos);
                test=4;
            }
        } while (test == 3);
        
        //antimuon neutrino
        do {
            double x=random.rand()*(-3);
            
            
            double F=distribution_Carceller(energy, pow(10, x), jcap, 5.1, 8.3);
            double Fmax=distribution_Carceller(energy, 0.001, jcap, 5.1, 8.3);
            double y=random.rand()*Fmax;
            double Eout=pow(10, x)*energy;
            if (y < F ){
                candidate->addSecondary(-14, Eout, pos);
                test=5;
            }
        } while (test == 4);
    }
    
    
    
    
} //~ namespace CRPropa



