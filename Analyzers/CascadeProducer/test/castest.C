//#include "SimpleFitter.h"
//in batch
//root -b -q castest.C\(\"dca2\"\)
//This macro  plot mass for R different cuts
void castest(string CasDir="dca1",string sufix="test")
{
   // make sure use plain style
   gROOT->SetStyle("Plain");         // PAW style
   gStyle->SetLabelFont(132, "XYZ"); //times new roman
   gSystem->Load("SimpleFitter_C.so");
   TString theName = "mycasanalyzerlooseskim2010";   // Prefix of root File
     string theNameFile = theName;
//   string xicut    = "prmrevxivee";                 //common cut
//   string xicut   = "prmreIP3D(#pi)";                     //common cut
//   string xicut   = "prmreIP3D";                     //common cut
//   string xicut   = "prmrevxiveeIP3D";              //common cut
   string xicut    = "prmrevxiveels3";              //common cut
//   string xicut    = "prmrevxiveels4";              //common cut

   string filename = theNameFile + sufix + ".root";                // The Name of Root File.
   TFile * file = new TFile( filename.c_str() );     // Open Root File.

   TString theTitle = "mass vs R";

   const int NHisto = 6;
   string label[NHisto];

   string num;
   for (int ilabel=0;ilabel<NHisto;ilabel++){
    stringstream flujo;
    flujo << ilabel;
    num = flujo.str();
    ~flujo;
    label[ilabel]  = CasDir + "/masshistos/" + xicut  + "_#Xi_"+num;
   }
//++fitter
   vector<double> parametros;
   parametros.push_back(0.55);  //const at mean
   parametros.push_back(0.000); //slope at mean
//   parametros.push_back(500.);  //event dca1fit2 dca2fit2 dca5fit2
   parametros.push_back(1600.);  //event dca02fit dca1 dca2
//   parametros.push_back(1.3223);//pico dca02fit dca1fit2 dca5fit2 dca10fit2 dca02 dca05
   parametros.push_back(1.32201);//pico dca1 dca2fit2
   parametros.push_back(0.0027);//ancho
   int bins = 75;
   double  x_min=1.295;
   double  x_max=1.37;

   SimpleFitter fitter;
   fitter.SetPars(parametros);
   fitter.SetBins(bins);
   fitter.SetXmin(x_min);
   fitter.SetXmax(x_max);
   fitter.SetFitMode("LEMR");
//--fitter
   TH1F* h1[NHisto];                         //Create histograms ids
   for ( int i = 0; i < NHisto  ; i++ ) {    //Load histograms from root file 
    h1[i]  = (TH1F*)file->Get(label[i].c_str());
    fitter.line_gaus(h1[i]);                 // do the fit
   }
   //Set printing style
   TCanvas * c1 = new TCanvas(theName+"01",theTitle, 600, 800);
   TCanvas * c2 = new TCanvas(theName+"02",theTitle, 600, 800);
   c1->Divide(2, 3);
   c2->Divide(1, 2);

   //Print Canvas 1
   for ( int i = 0; i < NHisto  ; i++ ) {
     c1->cd(i+1);
     h1[i]->SetMinimum(0);
     h1[i]->Draw();
   }
   gStyle->SetOptFit(1111);
   string myprint1 = CasDir + "_" + xicut + "_" + sufix + "1.ps";
   c1->cd(0);
   // Draw a global picture title
   title1 = new TPaveLabel(0.43,0.66,0.62,0.69,CasDir.c_str());
   title1->SetFillColor(16);
   title1->Draw();
   c1->Update();
   c1-> Print (myprint1.c_str());

   //Print Canvas 2
   TLegend *legend = new TLegend(0.78,0.63,0.98,0.83);
   c2->cd(1);
   h1[1]->Draw();
   c2->cd(2);
   h1[2]->Draw();
   legend -> AddEntry(h1[1],label[1].c_str(),"l");
   for ( int i = 3; i < NHisto  ; i++ ) {
     h1[i]->SetLineColor(i-1);
     h1[i]->Draw("Same");
     legend -> AddEntry(h1[i],label[i].c_str(),"l");
   }
   legend -> Draw();
   string myprint2 = CasDir + "_" + xicut + "_" + sufix + "2.ps";
   c2->cd(0);
   // Draw a global picture title
   title2 = new TPaveLabel(0.35,0.49,0.70,0.52,CasDir.c_str());
   title2->SetFillColor(16);
   title2->Draw();
   c2->Update();

   c2-> Print (myprint2.c_str());

}//end casdca1 macro
