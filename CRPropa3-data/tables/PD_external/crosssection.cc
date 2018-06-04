// Script for calculating the photodisintegration cross sections from various sources from Nils Nierstenhoefer

#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <TLegend.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <cmath>
#include <TAxis.h>
#include <TF1.h>
#include <TGraphErrors.h>

double mb_To_m_2 = 1;//1.e-31;
double MeV_to_J  = 1.602e-13;

double Emin = 0.001;
double Emax = 250.;
double Bins = 500.;

double SigmaLHat(double epsPrime, double E_0, double Gamma);

//returns the Bethe-Peierls cross section/m^2
//for an photon energy of Epsilon in MeV!
//see PhD thesis of Joerg Rachen page 71
double SigmaBP(double Epsilon /*MeV*/, double B /*MeV*/, double power){

  //all in si units!
  Epsilon*=MeV_to_J;
  B*=MeV_to_J;
  double sigma_Tp = 0.197*1.e-3*mb_To_m_2; //Rachen PhD thesis page 10
  double alpha_e  = 7.297e-3;
  double m_p      = 1.673e-27;               //kg
  double x        = Epsilon/B;
  double c        = 2.998e8;                 // m/s
  double sigma_BP = sigma_Tp*m_p*c*c*pow((x-1),1.5)/pow(x,power)/alpha_e/B;

  if(Epsilon < B){
    return(0.);
  }

  return(sigma_BP);
}

double LinApproxValTGraphErr(TGraphErrors* TheGraph, double xVal){

  double buffer;
  int NPoints = TheGraph->GetN();
  double MinX, MaxX;
  TheGraph->GetPoint(0, MinX, buffer);
  TheGraph->GetPoint(NPoints-1,MaxX, buffer);

  //std::cout<<"number of points="<<NPoints<<std::endl;
  //std::cout<<"First Point x="<<MinX<<std::endl;
  //std::cout<<"Last Point x="<<MaxX<<std::endl;

  if(xVal<=MinX || xVal>=MaxX) return(0.);

  double x_low = -1.;
  double y_low = -1.;
  double x_up  = -1.;
  double y_up = -1.;
  double xscan = MinX;
  double yscan = MinX;

  int pointNumber=0;
  while(xscan<xVal){
    x_low=xscan;
    y_low=yscan;
    TheGraph->GetPoint(pointNumber, xscan, yscan);
    if(pointNumber>NPoints) return(0.);
    pointNumber++;
  }

  TheGraph->GetPoint(pointNumber-1, x_up, y_up);

  //std::cout<<"x_low="<<x_low<<" y_low="<<y_low<<std::endl;
  //std::cout<<"x_up="<<x_up<<" y_up="<<y_up<<std::endl;

  double LinY=y_low-(y_low-y_up)*(x_low-xVal)/(x_low-x_up);
  //std::cout<<"LinY="<<LinY<<std::endl;
  return(LinY);
}


TGraph* OurLi6(){
  TCanvas* Li6XSDataCan = new TCanvas("Li6XSData","Li6XSData",1);
  TGraphErrors* gr1 = new TGraphErrors("dataFromDB/Li_6_3/xs100000.csv","%lg %lg %lg %lg");
  TGraphErrors* gr2 = new TGraphErrors("dataFromDB/Li_6_3/xs10000.csv","%lg %lg %lg %lg");
  TGraphErrors* grD = new TGraphErrors("dataFromDB/Li_6_3/xs100.csv","%lg %lg %lg %lg");

  std::vector<double> Energies;
  std::vector<double> xs100000Vec;
  std::vector<double> xs10000Vec;
  std::vector<double> xs100Vec;

  //OUTPUT FILES
  ofstream Li6_xs100000;
  Li6_xs100000.open("Li_6_3/xs100000.tot");
  ofstream Li6_xs10000;
  Li6_xs10000.open("Li_6_3/xs10000.tot");
  ofstream Li6_xs100;
  Li6_xs100.open("Li_6_3/xs100.tot");

  //Now the linear approximations
  for(double i=0; i<Bins; i++){

    double E = Emin + (Emax - Emin)  * i / (double) Bins;
    Energies.push_back(E);

    xs100Vec.push_back(LinApproxValTGraphErr(grD,E)*1000);
    Li6_xs100<<LinApproxValTGraphErr(grD,E)*1000.<<std::endl;

    xs100000Vec.push_back(LinApproxValTGraphErr(gr1,E)*1000);
    Li6_xs100000<<LinApproxValTGraphErr(gr1,E)*1000.<<std::endl;

    xs10000Vec.push_back(LinApproxValTGraphErr(gr2,E)*1000);
    Li6_xs10000<<LinApproxValTGraphErr(gr2,E)*1000.<<std::endl;

  }

  TGraph* xs100000GraphLinApp = new TGraph( xs100000Vec.size(), &(Energies[0]), &(xs100000Vec[0]));
  TGraph* xs10000GraphLinApp = new TGraph ( xs10000Vec.size(),  &(Energies[0]), &(xs10000Vec[0]));
  TGraph* xs100GraphLinApp = new TGraph ( xs100Vec.size(),  &(Energies[0]), &(xs100Vec[0]));

  xs10000GraphLinApp->SetMarkerStyle(24);
  xs10000GraphLinApp->SetMarkerColor(2);
  xs10000GraphLinApp->SetLineColor(2);
  xs10000GraphLinApp->Draw("Apl");

  xs100000GraphLinApp->SetMarkerStyle(24);
  xs100000GraphLinApp->SetMarkerColor(1);
  xs100000GraphLinApp->SetLineColor(1);
  xs100000GraphLinApp->Draw("pl");

  xs100GraphLinApp->SetMarkerStyle(24);
  xs100GraphLinApp->SetMarkerColor(4);
  xs100GraphLinApp->SetLineColor(4);
  xs100GraphLinApp->Draw("pl");

  gr1->SetMarkerColor(1);
  gr1->SetLineColor(1);
  gr1->SetMarkerStyle(22);
  gr2->SetMarkerColor(2);
  gr2->SetLineColor(2);
  gr2->SetMarkerStyle(23);
  gr2->Draw("p");
  gr1->Draw("p");

  grD->SetMarkerColor(4);
  grD->SetMarkerStyle(23);
  grD->Draw("p");

    Li6XSDataCan->SaveAs("dataFromDB/Li_6_3/Li6XSData.root");

  Li6_xs100000.close();
  Li6_xs10000.close();
  Li6_xs100.close();

}

TGraph* OurLi7(){
  TCanvas* Li7XSDataCan = new TCanvas("Li7XSData","Li6XSData",1);
  TGraphErrors* gr1 = new TGraphErrors("dataFromDB/Li_7_3/xs10000.csv","%lg %lg %lg %lg");
  TGraphErrors* gr2 = new TGraphErrors("dataFromDB/Li_7_3/gammaAbs.csv","%lg %lg %lg %lg");
  TGraphErrors* grT = new TGraphErrors("dataFromDB/Li_7_3/xs100.csv","%lg %lg %lg %lg");
  TGraphErrors* grA = new TGraphErrors("dataFromDB/Li_7_3/xs1.csv","%lg %lg %lg %lg");
  TGraphErrors* grND = new TGraphErrors("dataFromDB/Li_7_3/xs101000.csv","%lg %lg %lg %lg");

  // TF1* LorentzXs100000 = new TF1("LorentzXs100000","[0]*[1]*[1]/(pow(x-[2],2)+[1]*[1])",0,50);
//   LorentzXs100000->SetParameter(0,1.e-3);
//   LorentzXs100000->SetParameter(1,18);
//   LorentzXs100000->SetParameter(2,5);
//   gr1->Fit("LorentzXs100000");

  //   TF1* GaussXs100000 = new TF1("GaussXs100000","[0] * exp( -1*pow( (x-[1]) / [2] ,2 ) )",0,50);
  //   GaussXs100000->SetParameter(0,1.e-3);
  //   GaussXs100000->SetParameter(1,18);
  //   GaussXs100000->SetParameter(2,5);
  //   gr1->Fit("GaussXs100000");

  LinApproxValTGraphErr(gr1,12.);

  std::vector<double> Energies;
  std::vector<double> xs101000Vec;
  std::vector<double> xs100000Vec;
  std::vector<double> xs10000Vec;
  std::vector<double> xs100Vec;
  std::vector<double> xs1Vec;
  std::vector<double> gammaAbsVec;

  ofstream Li7_xs100;
  Li7_xs100.open("Li_7_3/xs000100.tot");
  ofstream Li7_xs10000;
  Li7_xs10000.open("Li_7_3/xs010000.tot");
  ofstream Li7_xs100000;
  Li7_xs100000.open("Li_7_3/xs100000.tot");
  ofstream Li7_xs101000;
  Li7_xs101000.open("Li_7_3/xs101000.tot");
  ofstream Li7_rp002004;
  Li7_rp002004.open("Li_7_3/rp002004.tot");


  //Now the linear approximations
  for(double i=0; i<Bins; i++){

    double E = Emin + (Emax - Emin)  * i / (double) Bins;
    Energies.push_back(E);

    xs1Vec.push_back(LinApproxValTGraphErr(grA,E));


    xs100Vec.push_back(LinApproxValTGraphErr(grT,E));
    Li7_xs100<<" "<<E<<" "<<LinApproxValTGraphErr(grT,E)*1000.<<" 0 0"<<std::endl;

    xs10000Vec.push_back(LinApproxValTGraphErr(gr1,E));
    Li7_xs10000<<" "<<E<<" "<<LinApproxValTGraphErr(gr1,E)*1000.<<" 0 0"<<std::endl;

    xs101000Vec.push_back(LinApproxValTGraphErr(grND,E));
    Li7_xs101000<<" "<<E<<" "<<LinApproxValTGraphErr(grND,E)*1000.<<" 0 0"<<std::endl;

    gammaAbsVec.push_back(LinApproxValTGraphErr(gr2,E));

    xs100000Vec.push_back(fabs(LinApproxValTGraphErr(gr2,E)
                   -LinApproxValTGraphErr(gr1,E)
                   //-GaussXs100000->Eval(E)
                   -LinApproxValTGraphErr(grT,E)));
    Li7_xs100000<<" "<<E<<" "<<fabs(LinApproxValTGraphErr(gr2,E)
               -LinApproxValTGraphErr(gr1,E)
               //-GaussXs100000->Eval(E)
               -LinApproxValTGraphErr(grT,E)
               -LinApproxValTGraphErr(grND,E))*1000.<<" 0 0"<<std::endl;

    Li7_rp002004<<" "<<E<<" "<<
      fabs(//LinApproxValTGraphErr(gr2,E)
      //-LinApproxValTGraphErr(gr1,E)
      //-GaussXs100000->Eval(E)
               LinApproxValTGraphErr(grT,E)
               +LinApproxValTGraphErr(grND,E)
       )*1000.<<" 0 0"<<std::endl;
  }

  TGraph* xs100000GraphLinApp = new TGraph (xs100000Vec.size(), &(Energies[0]), &(xs100000Vec[0]));
  TGraph* xs101000GraphLinApp = new TGraph (xs101000Vec.size(), &(Energies[0]), &(xs101000Vec[0]));
  TGraph* xs10000GraphLinApp  = new TGraph (xs10000Vec.size(),  &(Energies[0]), &(xs10000Vec[0]));
  TGraph* xs1GraphLinApp      = new TGraph (xs1Vec.size(),  &(Energies[0]), &(xs1Vec[0]));
  TGraph* gammaAbsGraphLinApp = new TGraph (gammaAbsVec.size(),  &(Energies[0]), &(gammaAbsVec[0]));
  TGraph* xs100GraphLinApp    = new TGraph  (xs100Vec.size(),  &(Energies[0]), &(xs100Vec[0]));

  gr1->SetMarkerColor(2);
  gr1->SetLineColor(2);
  gr1->SetMarkerStyle(22);

  gr2->SetMarkerColor(3);
  gr2->SetLineColor(3);
  gr2->SetMarkerStyle(23);

  grT->SetMarkerColor(6);
  grT->SetLineColor(6);
  grT->SetMarkerStyle(23);

  grA->SetMarkerColor(7);
  grA->SetLineColor(7);
  grA->SetMarkerStyle(23);

  grND->SetMarkerColor(9);
  grND->SetLineColor(9);
  grND->SetMarkerStyle(25);

  gr2->Draw("Ap");
  gr1->Draw("p");
  grT->Draw("p");
  grND->Draw("p");
  //grA->Draw("p");

  xs1GraphLinApp->SetMarkerColor(7);
  xs1GraphLinApp->SetLineColor(7);
  xs1GraphLinApp->SetMarkerStyle(23);
  //xs1GraphLinApp->Draw("p");

  xs100GraphLinApp->SetMarkerColor(6);
  xs100GraphLinApp->SetLineColor(6);
  xs100GraphLinApp->SetMarkerStyle(23);
  xs100GraphLinApp->Draw("p");

  xs10000GraphLinApp->SetMarkerStyle(24);
  xs10000GraphLinApp->SetMarkerColor(2);
  xs10000GraphLinApp->SetLineColor(2);
  xs10000GraphLinApp->Draw("p");

  gammaAbsGraphLinApp->SetMarkerColor(3);
  gammaAbsGraphLinApp->SetLineColor(3);
  gammaAbsGraphLinApp->SetMarkerStyle(25);
  gammaAbsGraphLinApp->Draw("p");

  xs100000GraphLinApp->SetMarkerStyle(24);
  xs100000GraphLinApp->SetMarkerColor(4);
  xs10000GraphLinApp->SetLineColor(4);
  xs100000GraphLinApp->Draw("p");

  xs101000GraphLinApp->SetMarkerStyle(24);
  xs101000GraphLinApp->SetMarkerColor(4);
  xs101000GraphLinApp->SetLineColor(4);
  xs101000GraphLinApp->Draw("p");

  Li7XSDataCan->SaveAs("dataFromDB/Li_7_3/Li7XSData.root");
}

//returns the cross section/m^2 of the reaction
//D+\gamma -> n+p
//for an photon energy of Epsilon in MeV!
//see PhD thesis of Joerg Rachen page 81
double Deuteron(double Epsilon /*MeV*/, double B /*MeV*/){

  //all in si units
  Epsilon*=MeV_to_J;
  B*=MeV_to_J;

  double m_p      = 1.673e-27;//938.3;                //MeV //1.673e-27;                //kg
  double sigma_BP = SigmaBP(Epsilon/MeV_to_J, B/MeV_to_J,3);  //m^2
  double h_bar    = 6.582e-16*1.e-6*MeV_to_J;                        //MeV
  double a        = sqrt(m_p*B)/h_bar;
  double r_eff    = 1.79e-15;                 //m
  double sigma_d  = sigma_BP/(1.-a*r_eff);

  return(sigma_d);
}


//???????????????????????????????
double Pl(double x, double x_th, double x_max, double alpha){

  double A       = alpha*x_max/x_th;
  double Factor1 = pow((x-x_th)/(x_max-x_th),A-alpha);
  double Factor2  = pow(x/x_max,-1.*A);

  if(x<x_th) return(0);

  return(Factor1*Factor2);
}

//Lorentz curve needed to paramterize Be-9 according to Rachen
//compare with eq. 1.1.56a from Rachen's PhD thesis (page 30).
double SigmaLHat(double epsPrime, double E_0, double Gamma){
  double gamma_2=Gamma*Gamma;
  return( gamma_2  / ( pow(epsPrime - E_0, 2) + gamma_2) );
}



//nparam[0] : Amplitude
//nparam[1] : gamma
//nparam[2] : E_0
Double_t SigmaBe9AllPeaksRoot(Double_t* x, Double_t* par){

  double gamma_2= par[1] * par[1];
  double Peak1  = par[0] * gamma_2  / ( pow(x[0] - par[2], 2) + gamma_2);

  double  gamma_3 = par[4] * par[4];
  double Peak2  = par[3] * gamma_3  / ( pow(x[0] - par[5], 2) + gamma_3);

  double gamma_4 = par[7] * par[7];
  double Peak3  = par[6] * gamma_4  / ( pow(x[0] - par[8], 2) + gamma_4);

  double Peak4  = par[9] * SigmaBP(x[0], par[10], 4);

  return(Peak1+Peak2+Peak3+Peak4);
  //return( nparam[0]*epsPrime[0]);
}

TF1* OurBe9Fit(){

  TF1* SummedUpPeaks = new TF1("SummedUpPeaks",SigmaBe9AllPeaksRoot,0.,90.,11);

  SummedUpPeaks->FixParameter(0,1.2);
  SummedUpPeaks->SetParameter(1,0.1);
  SummedUpPeaks->SetParameter(2,2.9);
  SummedUpPeaks->SetParameter(3,1.6);
  SummedUpPeaks->SetParameter(4,3.0);
  SummedUpPeaks->SetParameter(5,8.5);
  SummedUpPeaks->SetParameter(6,5.1);
  SummedUpPeaks->SetParameter(7,6.);
  SummedUpPeaks->SetParameter(8,27.);
  SummedUpPeaks->SetParameter(9,9.);
  SummedUpPeaks->SetParameter(10,16.9);

  SummedUpPeaks->SetNpx(5000);

  TCanvas* Be9XSDataCan = new TCanvas("Be9XSData","Be9XSData",1);
  TGraphErrors* gr1 = new TGraphErrors("be9-g-abs.csv","%lg,%lg,%lg,%lg");
  gr1->Draw("Ap");
  gr1->Fit("SummedUpPeaks","R");
  SummedUpPeaks->Draw();
  gr1->Draw("Ap");
  SummedUpPeaks->Draw();
  gr1->Draw("p");
  Be9XSDataCan->SaveAs("Be9XSDataAbs.root");

  return(SummedUpPeaks);
}


int main(){

  TLegend *legend = new TLegend(.75,.80,.95,.95);

  //The files in the TALYS output format as needed by CRPropa.
  //Beryllium
  ofstream Berillium_xs100002;
  Berillium_xs100002.open("Be_9_4/xs100001.tot");

  //Deuterium
  ofstream Deuterium_xs110000;
  Deuterium_xs110000.open("D_2_1/xs100000.tot");

  //Tritium
  ofstream Tritium_xs101000;
  Tritium_xs101000.open("T_3_1/xs100000.tot");
  ofstream Tritium_xs210000;
  Tritium_xs210000.open("T_3_1/xs200000.tot");

  //Alpha3=Tritium due to isospin symmetry
  ofstream Alpha3_xs101000;
  Alpha3_xs101000.open("He_3_2/xs010000.tot");
  ofstream Alpha3_xs210000;
  Alpha3_xs210000.open("He_3_2/xs110000.tot");

  //Alpha
  ofstream Alpha_xs100000;
  Alpha_xs100000.open("He_4_2/xs100000.tot");
  ofstream Alpha_xs110000;
  Alpha_xs110000.open("He_4_2/xs110000.tot");

  std::vector<double> Energies;
  std::vector<double> Cross_Sec;

  //Tritium
  std::vector<double> CS_Tritium_xs101000Vec;
  std::vector<double> CS_Tritium_xs210000Vec;

  //Alpha 3
  std::vector<double> CS_Alpha3_xs101000Vec;
  std::vector<double> CS_Alpha3_xs210000Vec;

  //Alpha
  std::vector<double> CS_Alpha_xs100000Vec;
  std::vector<double> CS_Alpha_xs110000Vec;

  //Be-9
  std::vector<double> Berillium_xs100002Vec;
  TF1* SummedUpPeaks = OurBe9Fit();

  OurLi6();
  OurLi7();

  for(double i=0; i<Bins; i++){

    double E = Emin + (Emax - Emin)  * i / (double) Bins;
    Energies.push_back(E);

    //Beryllium
    //Rachen
    /*  double BeCS=1.5*SigmaLHat(E, 1.7,  0.3)
      +
      1.6*SigmaLHat(E, 10., 10.)
      +
      3.5*SigmaLHat(E, 25., 15.)
      +
      9*SigmaBP(E, 16.9)
      ;*/

    //Our Fit (Apply hard cut-off at 1.7 MeV)
    double BeCS=0.;
    if(E>1.7){BeCS=SummedUpPeaks->Eval(E);}

    Berillium_xs100002Vec.push_back(BeCS);
    Berillium_xs100002<<" "<<E<<" "<<BeCS/mb_To_m_2<<" 0 0"<<std::endl;
    std::cout<<E<<" "<<BeCS/mb_To_m_2<<" 0 0"<<std::endl;

    //Deuteron
    Cross_Sec.push_back(Deuteron(E, 2.227)/mb_To_m_2);
    Deuterium_xs110000<<" "<<E<<" "<<Deuteron(E, 2.227)/mb_To_m_2<<" 0 0"<<std::endl;

    //Trition
    CS_Tritium_xs101000Vec.push_back(1.7*1.4*SigmaBP(E, 5.8, 3)/mb_To_m_2);
    CS_Tritium_xs210000Vec.push_back(1.7*1.7*SigmaBP(E, 7.3, 3)/mb_To_m_2);
    Tritium_xs101000<<" "<<E<<" "<<1.7*1.4*SigmaBP(E, 5.8, 3)/mb_To_m_2<<" 0 0"<<std::endl;
    Tritium_xs210000<<" "<<E<<" "<<1.7*1.7*SigmaBP(E, 7.3, 3)/mb_To_m_2<<" 0 0"<<std::endl;

    //Alpha3
    CS_Alpha3_xs101000Vec.push_back((1/1.5)*1.4*SigmaBP(E, 5.8, 3)/mb_To_m_2);
    CS_Alpha3_xs210000Vec.push_back((1/1.5)*1.7*SigmaBP(E, 7.3, 3)/mb_To_m_2);
    Alpha3_xs101000<<" "<<E<<" "<<(1/1.5)*1.4*SigmaBP(E, 5.8, 3)/mb_To_m_2<<" 0 0"<<std::endl;
    Alpha3_xs210000<<" "<<E<<" "<<(1/1.5)*1.7*SigmaBP(E, 7.3, 3)/mb_To_m_2<<" 0 0"<<std::endl;

    //Alpha
    Alpha_xs100000<<" "<<E<<" "<<3.6*Pl(E,19.8,27,5)<<" 0 0"<<std::endl;
    Alpha_xs110000<<" "<<E<<" "<<1.4*SigmaBP(E, 26.1, 3)/mb_To_m_2<<" 0 0"<<std::endl;
    CS_Alpha_xs100000Vec.push_back(3.6*Pl(E,19.8,27.,5));
    CS_Alpha_xs110000Vec.push_back(1.4*SigmaBP(E, 26.1,3)/mb_To_m_2);
  }

  std::cout<<"Prob="<<SummedUpPeaks->GetProb()<<std::endl;
  std::cout<<"Chi_2="<<SummedUpPeaks->GetChisquare()<<std::endl;
  std::cout<<"NDoF="<<SummedUpPeaks->GetNDF()<<std::endl;


  TCanvas* RachenCrossCan = new TCanvas("RachenCrossCan","RachenCrossCan",1);

  //Plot the cross sections
  TGraph* Be_g_To_2h_n_Graph = new TGraph(Berillium_xs100002Vec.size(),
                    &(Energies)[0],
                    &(Berillium_xs100002Vec)[0]);
  Be_g_To_2h_n_Graph->Draw("AL");
  Be_g_To_2h_n_Graph->SetLineColor(2);
  Be_g_To_2h_n_Graph->SetLineStyle(2);
  Be_g_To_2h_n_Graph->SetLineWidth(2);
  Be_g_To_2h_n_Graph->GetXaxis()->SetTitle("Photon energy #epsilon' in cms / MeV");
  Be_g_To_2h_n_Graph->GetYaxis()->SetTitle("cross section #sigma / mbarn");
  legend->AddEntry(Be_g_To_2h_n_Graph,"Be-9+#gamma#rightarrow2#alpha+n");

  //Deuterium
  TGraph* D_g_To_p_n_Graph = new TGraph(Cross_Sec.size(),
                    &(Energies)[0],
                    &(Cross_Sec)[0]);
  D_g_To_p_n_Graph->Draw("SAMEL");
  D_g_To_p_n_Graph->SetLineColor(2);
  D_g_To_p_n_Graph->SetLineStyle(2);
  D_g_To_p_n_Graph->SetLineWidth(2);
  //D_g_To_p_n_Graph->GetXaxis()->SetTitle("Photon energy #epsilon' in cms / MeV");
  //D_g_To_p_n_Graph->GetYaxis()->SetTitle("cross section #sigma / mbarn");
  legend->AddEntry(D_g_To_p_n_Graph,"d+#gamma#rightarrown+p");
  //Tritium
    //101000
  TGraph* T_g_To_101000_Graph = new TGraph(Cross_Sec.size(),
                       &(Energies)[0],
                       &(CS_Tritium_xs101000Vec)[0]);
  T_g_To_101000_Graph->Draw("SAMEL");
  T_g_To_101000_Graph->SetLineColor(3);
  T_g_To_101000_Graph->SetLineStyle(3);
  T_g_To_101000_Graph->SetLineWidth(2);
  legend->AddEntry(T_g_To_101000_Graph,"trinucleon (T )+#gamma#rightarrowd+n/p ");
    //210000
  TGraph* T_g_To_210000_Graph = new TGraph(Cross_Sec.size(),
                       &(Energies)[0],
                       &(CS_Tritium_xs210000Vec)[0]);
  T_g_To_210000_Graph->Draw("SAMEL");
  T_g_To_210000_Graph->SetLineColor(3);
  T_g_To_210000_Graph->SetLineStyle(4);
  T_g_To_210000_Graph->SetLineWidth(2);
  legend->AddEntry(T_g_To_210000_Graph,"T+#gamma#rightarrow2n+p or n+2p");




  //Alpha 3
  //101000
  TGraph* Alpha3_g_To_101000_Graph = new TGraph(Cross_Sec.size(),
                       &(Energies)[0],
                       &(CS_Alpha3_xs101000Vec)[0]);
  Alpha3_g_To_101000_Graph->Draw("SAMEL");
  Alpha3_g_To_101000_Graph->SetLineColor(6);
  Alpha3_g_To_101000_Graph->SetLineStyle(7);
  Alpha3_g_To_101000_Graph->SetLineWidth(2);
  legend->AddEntry(Alpha3_g_To_101000_Graph,"He-3+#gamma#rightarrowd+n/p ");
    //210000
  TGraph* Alpha3_g_To_210000_Graph = new TGraph(Cross_Sec.size(),
                       &(Energies)[0],
                       &(CS_Alpha3_xs210000Vec)[0]);
  Alpha3_g_To_210000_Graph->Draw("SAMEL");
  Alpha3_g_To_210000_Graph->SetLineColor(6);
  Alpha3_g_To_210000_Graph->SetLineStyle(8);
  Alpha3_g_To_210000_Graph->SetLineWidth(2);
  legend->AddEntry(Alpha3_g_To_210000_Graph,"He-3+#gamma#rightarrow2n+p or n+2p");




  //Alpha
    //100000
   TGraph* A_g_To_100000_Graph = new TGraph(Cross_Sec.size(),
                       &(Energies)[0],
                       &(CS_Alpha_xs100000Vec)[0]);
   A_g_To_100000_Graph->Draw("SAMEL");
   A_g_To_100000_Graph->SetLineColor(4);
   A_g_To_100000_Graph->SetLineStyle(4);
   A_g_To_100000_Graph->SetLineWidth(2);
   legend->AddEntry(A_g_To_100000_Graph,"#alpha+#gamma#rightarrowHe-3+n");
   //110000
   TGraph* A_g_To_110000_Graph = new TGraph(Cross_Sec.size(),
                       &(Energies)[0],
                       &(CS_Alpha_xs110000Vec)[0]);
   A_g_To_110000_Graph->Draw("SAMEL");
   A_g_To_110000_Graph->SetLineColor(4);
   A_g_To_110000_Graph->SetLineStyle(5);
   A_g_To_110000_Graph->SetLineWidth(2);
   legend->AddEntry(A_g_To_110000_Graph,"#alpha+#gamma#rightarrowd+p+n ");

   legend->Draw("same");

   RachenCrossCan->SaveAs("RachenCrossCan.root");

  //close files
  Deuterium_xs110000.close();
  Tritium_xs101000.close();
  Tritium_xs210000.close();

  //Create the files with the nucleus transitions by copying
  //Be-9
  system("cp Be_9_4/xs100001.tot Be_9_4/rp002004.tot");
  //Deuterium
  system("cp D_2_1/xs100000.tot D_2_1/rp001001.tot");
  //Tritium
  system("cp T_3_1/xs100000.tot T_3_1/rp001002.tot");
  system("cp T_3_1/xs200000.tot T_3_1/rp001001.tot");
  //Alpha3
  system("cp He_3_2/xs010000.tot He_3_2/rp001002.tot");
  system("cp He_3_2/xs110000.tot He_3_2/rp001001.tot");
  //Alpha
  system("cp He_4_2/xs100000.tot He_4_2/rp002003.tot");
  system("cp He_4_2/xs110000.tot He_4_2/rp001002.tot");
  //Li-7
  system("cp Li_7_3/xs100000.tot Li_7_3/rp003006.tot");
  system("cp Li_7_3/xs010000.tot Li_7_3/rp002006.tot");
}
