
#include <TH2.h>
#include <TStyle.h>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void fitBS(TString filename = "EarlyCollision_hlt.root") {
	TFile *f = TFile::Open(filename);

    TTree*tree = (TTree*)gDirectory->Get("mytree");
    TTree*treeBS = (TTree*)gDirectory->Get("mytreebeam");
 
 //  Init(treeBS);
   Double_t        ftrack;
   Double_t        beam_x0;
   Double_t        beam_y0;
   Double_t        beam_z0;
   Double_t        beam_sigmaZ;
   Double_t        beam_dxdz;
   Double_t        beam_dydz;
   Double_t        beam_cov[7][7];

   treeBS->SetBranchAddress("ftrack", &ftrack);
   treeBS->SetBranchAddress("beam_x0", &beam_x0);
   treeBS->SetBranchAddress("beam_y0", &beam_y0);
   treeBS->SetBranchAddress("beam_z0", &beam_z0);
   treeBS->SetBranchAddress("beam_sigmaZ", &beam_sigmaZ);
   treeBS->SetBranchAddress("beam_dxdz", &beam_dxdz);
   treeBS->SetBranchAddress("beam_dydz", &beam_dydz);
   treeBS->SetBranchAddress("beam_cov", beam_cov);  
   
   Int_t nevent = treeBS->GetEntries();
   
   const Int_t NDIM = nevent; 
   Float_t ntrk[NDIM];
   Float_t etrk[NDIM];
   Float_t x[NDIM];
   Float_t y[NDIM];
   Float_t z[NDIM];
   Float_t dxdz[NDIM];
   Float_t dydz[NDIM];
   Float_t ex[NDIM];
   Float_t ey[NDIM];
   Float_t ez[NDIM];
   Float_t edxdz[NDIM];
   Float_t edydz[NDIM];
 
 cout<<  nevent <<" "<< endl;
  for (Int_t i=0;i<nevent;i++) {
     treeBS->GetEvent(i);                  //read branch "fNtrack" only
      ntrk[i] = ftrack;
      x[i] = beam_x0;
      etrk[i] = 0;
      Float_t ex2 = 0;
      if (beam_cov[0][0] > 0) ex2 = TMath::Sqrt(beam_cov[0][0]);
      ex[i] = ex2;
      cout<< i <<" "<< ntrk[i] <<" "<< x[i] <<" "<< ex[i] <<endl;
      
      y[i] = beam_y0;
      Float_t ey2 = 0;
      if (beam_cov[0][0] > 0) ey2 = TMath::Sqrt(beam_cov[1][1]);
      ey[i] = ey2;
      cout<< i <<" "<< ntrk[i] <<" "<< y[i] <<" "<< ey[i] <<endl;
      
      z[i] = beam_z0;
      Float_t ez2 = 0;
      if (beam_cov[0][0] > 0) ez2 = TMath::Sqrt(beam_cov[2][2]);
      ez[i] = ez2;
      cout<< i <<" "<< ntrk[i] <<" "<< z[i] <<" "<< ez[i] <<endl;
      
      dxdz[i] = beam_dxdz;
      Float_t edxdz2 = 0;
      if (beam_cov[0][0] > 0) edxdz2 = TMath::Sqrt(beam_cov[4][4]);
      edxdz[i] = edxdz2;
      cout<< i <<" "<< ntrk[i] <<" "<< dxdz[i] <<" "<< edxdz[i] <<endl;
      
      dydz[i] = beam_dydz;
      Float_t edydz2 = 0;
      if (beam_cov[0][0] > 0) edydz2 = TMath::Sqrt(beam_cov[5][5]);
      edydz[i] = ez2;
      cout<< i <<" "<< ntrk[i] <<" "<< dydz[i] <<" "<< edydz[i] <<endl;
      
    //  ++i;
   }
   
   TGraphErrors* gr_beam_x0 = new TGraphErrors(nevent, ntrk, x, etrk, ex);
   gr_beam_x0->SetMarkerStyle(20);
   gr_beam_x0->SetMarkerColor(2);
   gr_beam_x0->SetTitle("x0 vs Tracks");
    
   TGraphErrors* gr_beam_y0 = new TGraphErrors(nevent, ntrk, y, etrk, ey);
   gr_beam_y0->SetMarkerStyle(20);
   gr_beam_y0->SetMarkerColor(2);
   gr_beam_y0->SetTitle("y0 vs Tracks");
   
   TGraphErrors* gr_beam_z0 = new TGraphErrors(nevent, ntrk, z, etrk, ez);
   gr_beam_z0->SetMarkerStyle(20);
   gr_beam_z0->SetMarkerColor(2);
   gr_beam_z0->SetTitle("z0 vs Tracks");

   TGraphErrors* gr_beam_dxdz = new TGraphErrors(nevent, ntrk, dxdz, etrk, edxdz);
   gr_beam_dxdz->SetMarkerStyle(20);
   gr_beam_dxdz->SetMarkerColor(2);
   gr_beam_dxdz->SetTitle("dxdz vs Tracks");
   
   TGraphErrors* gr_beam_dydz = new TGraphErrors(nevent, ntrk, dydz, etrk, edydz);
   gr_beam_dydz->SetMarkerStyle(20);
   gr_beam_dydz->SetMarkerColor(2);
   gr_beam_dydz->SetTitle("dydz vs Tracks");
 
  TCanvas *cv = new TCanvas("cv","cv",900,1200);
  	cv->Divide(2,3);
	cv->cd(1);	
      tree->Draw("d0:phi0");
	  	title = new TPaveLabel(-4.5,0.9,-2.8,0.81,"d0 vs phi0");
        //title->SetFillColor(2);
        title->Draw();

    cv->cd(2);
	  gStyle->SetOptFit(0111); //Style Show fit
      gr_beam_z0->Draw("awp");
      gr_beam_z0->GetXaxis()->SetTitle("Tracks");
      gr_beam_z0->GetYaxis()->SetTitle("z0[cm]");
      gr_beam_z0->Fit("pol0","LM");
      pol0->SetLineColor(4);
    cv->cd(3);	
      gStyle->SetOptFit(0111); //Style Show fit
      gr_beam_x0->Fit("pol0","LM");
      pol0->SetLineColor(4);
      gr_beam_x0->Draw("AP");
      gr_beam_x0->GetXaxis()->SetTitle("Tracks");
      gr_beam_x0->GetYaxis()->SetTitle("x0[cm]");
    cv->cd(4);
	  gStyle->SetOptFit(0111); //Style Show fit
      gr_beam_y0->Draw("awp");
      gr_beam_y0->GetXaxis()->SetTitle("Tracks");
      gr_beam_y0->GetYaxis()->SetTitle("y0[cm]");
      gr_beam_y0->Fit("pol0","LM");
      pol0->SetLineColor(4);
    cv->cd(5);
      gStyle->SetOptFit(0111); //Style Show fit
      gr_beam_dxdz->Draw("awp");
      gr_beam_dxdz->GetXaxis()->SetTitle("Tracks");
      gr_beam_dxdz->GetYaxis()->SetTitle("dx/dz");
      gr_beam_dxdz->Fit("pol0","LM");
      pol0->SetLineColor(4);
    cv->cd(6);
	 gStyle->SetOptFit(0111); //Style Show fit
     gr_beam_dydz->Draw("awp");
     gr_beam_dydz->GetXaxis()->SetTitle("Tracks");
     gr_beam_dydz->GetYaxis()->SetTitle("dy/dz");
     gr_beam_dydz->Fit("pol0","LM");
     pol0->SetLineColor(4);
	cv->Update();
	cv->cd(0);
 
	cv-> Print ("hltbeamspot.png");
	cv-> Print ("hltbeamspot.eps");

   new TCanvas();
   gr_beam_y0->Draw("awp");
   new TCanvas();
   gr_beam_z0->Draw("awp");
   new TCanvas();
   gr_beam_dxdz->Draw("awp");
   new TCanvas();
   gr_beam_dydz->Draw("awp");
}
