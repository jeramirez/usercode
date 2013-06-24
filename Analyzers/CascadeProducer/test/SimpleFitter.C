#define SimpleFitter_cxx
#include "SimpleFitter.h"
#include <TROOT.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>
#include <vector>


// .L SimpleFitter.C++
// fitter = new SimpleFitter ()
// fitter->doubleGaussianFit(histogram,cout)


void SimpleFitter::SetBins(int mybins)
{
bins = mybins;
return;
}
void SimpleFitter::SetXmin(double myxmin)
{
x_min = myxmin;
return;
}
void SimpleFitter::SetXmax(double myxmax)
{
x_max = myxmax;
return;
}

void SimpleFitter::doubleGaussianFit(TH1F *plot, std::ostream &out)
{
  TF1 *myfit = new TF1("myfit","[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[4])/[5])^2)");
  myfit->SetParameter(0, 1);
  myfit->SetParameter(1, plot->GetMean());
  myfit->SetParameter(2, plot->GetRMS()/2);
  myfit->SetParameter(3, 1);
  myfit->SetParameter(4, plot->GetMean());
  myfit->SetParameter(5, plot->GetRMS()*1.5);

  plot->Fit(myfit,"QO");
  out << myfit->GetParameter(2)<<" "<<myfit->GetParameter(5)<<std::endl;
  return;
  
}

void SimpleFitter::argus_gaus(TH1F *plot, std::ostream &out)
{

  int x_range = (x_max-x_min)*1000;
  int xxmin   = x_min*1000;


  TF1 *myfit = new TF1("f1", Form("%s%i%s%i%s%i%s%i%s","[0]*((x-", xxmin,  
  "/1000.)**(0.5)) +[1]*((x-", xxmin, "/1000.)**(3/2)) + ((", x_range,
  "/1000.)/", bins, ")*gausn(2)"), x_min, x_max);
  for (unsigned int ipar=0;ipar<Pars.size();ipar++){
  myfit->SetParameter(ipar, Pars[ipar]);
 }
  myfit->SetParName(0,"bkg param 1");
  myfit->SetParName(1,"bkg param 2");
  myfit->SetParName(2,"yield");
  myfit->SetParName(3,"mean");
  myfit->SetParName(4,"sigma");
  
  if ( FitMode=="") FitMode="LEM"; 
  plot->Fit(myfit,FitMode.c_str());
  out << myfit->GetParameter(2)<<" "<<myfit->GetParameter(4)<<std::endl;
  
 return;

}

void SimpleFitter::line_gaus(TH1F *plot, std::ostream &out)
{

//  int x_range = (x_max-x_min)*1000;
  int nbins = plot->GetNbinsX(); 
  double bin_width = plot->GetBinWidth(1);
  double histo_max = plot->GetBinLowEdge( plot->GetNbinsX() ) + bin_width;
  double histo_min = plot->GetBinLowEdge(1);
  int x_range = ( histo_max - histo_min )*1000;

  TF1 *myfit = new TF1("f2", Form("%s%i%s%i%s","[0] +  [1]*(x-[3]) + ((", x_range,
  "/1000.)/", nbins, ")*gausn(2)"), x_min, x_max);
  for (unsigned int ipar=0;ipar<Pars.size();ipar++){
  myfit->SetParameter(ipar, Pars[ipar]);
 }
  myfit->SetParName(0,"const");
  myfit->SetParName(1,"slope_m");
  myfit->SetParName(2,"yield");
  myfit->SetParName(3,"mean");
  myfit->SetParName(4,"sigma");
  
  if ( FitMode=="") FitMode="LEM"; 
  plot->Fit(myfit,FitMode.c_str());
  out << myfit->GetParameter(2)<<" "<<myfit->GetParameter(4)<<std::endl;
  
 return;

}

void SimpleFitter::SetPars(vector<double> par)
{
if  (par.size()>0) Pars = par;
return;
}
void SimpleFitter::SetFitMode(string fitmode)
{

FitMode =fitmode; 
return;
}

   

