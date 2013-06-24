#ifndef SimpleFitter_h
#define SimpleFitter_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include "TString.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>

  /**
   * Fitter
   */

class SimpleFitter {
   public :

   SimpleFitter();
   ~SimpleFitter();
   void SetBins(int mybins);
   void SetXmin(double myxmin);
   void SetXmax(double myxmax);
   void SetPars(vector<double> par);
   void doubleGaussianFit(TH1F *plot, ostream &out = std::cout);
   void argus_gaus(TH1F *plot, ostream &out = std::cout);
   void line_gaus(TH1F *plot, ostream &out = std::cout);
   void SetFitMode(string fitmode);

private:

   int bins;
   double x_max;
   double x_min;
   string FitMode;
   vector<double> Pars;
   
 
 
};

#endif

#ifdef SimpleFitter_cxx
SimpleFitter::SimpleFitter()
{
bins=100;
x_max=2;
x_min=1;
}
SimpleFitter::~SimpleFitter()
{
}
#endif // #ifdef SimpleFitter_cxx


