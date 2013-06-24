//in batch 
//root -b -q tracks_in_jet.C\(\"btagpatanalyzer.root\"\)
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void tracks_in_jet(TString filename="btagpatanalyzerpy_ttbar.root") {

	
	TFile *f = TFile::Open(filename);
        TH1D* histo[16];
	TCanvas *cv = new TCanvas("cv","cv",700,1000);
	TLegend *legend[4];

	legend[0] = new TLegend(0.72,0.62,0.98,0.82);
	legend[1] = new TLegend(0.72,0.62,0.98,0.82);
	legend[2] = new TLegend(0.72,0.62,0.98,0.82);
	legend[3] = new TLegend(0.72,0.62,0.98,0.82);
	
        histo[0]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_no_flavor");
        histo[1]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_b");
        histo[2]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_c");
        histo[3]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_udsg");

	//pt <30
        histo[4]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_no_flavor030");
        histo[5]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_b030");
        histo[6]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_c030");
        histo[7]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_udsg030");

	// 30< pt <50	
        histo[8]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_no_flavor3050");
        histo[9]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_b3050");
        histo[10]   = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_c3050");
        histo[11]   = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_udsg3050");

	//  pt > 50	
        histo[12]   = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_no_flavor50");
        histo[13]   = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_b50");
        histo[14]   = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_c50");
        histo[15]   = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/tracks_in_jet_udsg50");


	cv->Divide(2,2);
	for ( size_t i = 0; i < 4; ++i ){
		cv->cd(1+i);
		histo[i]->GetXaxis()->SetTitle("# tracks/jet");
		legend[i] -> AddEntry(histo[i],"all jet pt","l");
		histo[i]->Draw("LE");
		histo[i+4]->SetLineColor(2);
		legend[i] -> AddEntry(histo[i+4],"pt < 30","l");
		histo[i+4]->Draw("Same");
		histo[i+8]->SetLineColor(4);
		legend[i] -> AddEntry(histo[i+8],"30 < pt < 50","l");
		histo[i+8]->Draw("Same");
		histo[i+12]->SetLineColor(6);
		legend[i] -> AddEntry(histo[i+12],"pt > 50","l");
		histo[i+12]->Draw("Same");
		legend[i] ->Draw();
	}
	cv->Update();
	cv->cd(0);
       // Draw a global picture title
	title = new TPaveLabel(0.35,0.49,0.70,0.52,"# tracks/jet distribution");
        title->SetFillColor(16);
        title->Draw();
        cv->Update();
	cv-> Print ("BTagPATtracks_in_jet.eps");
	
}
