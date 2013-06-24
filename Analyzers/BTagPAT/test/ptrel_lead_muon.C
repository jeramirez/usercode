//in batch 
//root -b -q ptrel_lead_muon.C\(\"btagpatclosurepy.root\"\)
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void ptrel_lead_muon(TString filename="btagpatclosurepy.root", TString tagger="trackCountingHighEffBJetTags") {

        Int_t pTrelRebin = 1;
	
	TFile *f = TFile::Open(filename);
        TH1D* histo[12];
	TCanvas *cv = new TCanvas("cv","cv",700,1000);
	TLegend *legend[4];

        TString labeltag = "TCHE";
        if (tagger=="trackCountingHighPurBJetTags")labeltag = "TCHP";
        if (tagger=="jetProbabilityBJetTags")labeltag = "JP";
        if (tagger=="jetBProbabilityBJetTags")labeltag = "JBP";
        if (tagger=="simpleSecondaryVertexBJetTags")labeltag = "SSV";
        if (tagger=="combinedSecondaryVertexBJetTags")labeltag = "CSV";

	legend[0] = new TLegend(0.72,0.62,0.98,0.82);
	legend[1] = new TLegend(0.72,0.62,0.98,0.82);
	legend[2] = new TLegend(0.72,0.62,0.98,0.82);
	legend[3] = new TLegend(0.72,0.62,0.98,0.82);
	//all muons associated to a jet
        histo[0]    = (TH1D*) gDirectory->Get("BTagPATClosure/muDir/jet_ptrel_no_flavor");
        histo[1]    = (TH1D*) gDirectory->Get("BTagPATClosure/muDir/jet_ptrel_b");
        histo[2]    = (TH1D*) gDirectory->Get("BTagPATClosure/muDir/jet_ptrel_c");
        histo[3]    = (TH1D*) gDirectory->Get("BTagPATClosure/muDir/jet_ptrel_udsg");

	//leading muon asociated to a jet
        histo[4]    = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_leadmu_no_flavor");
        histo[5]    = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_leadmu_b");
        histo[6]    = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_leadmu_c");
        histo[7]    = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_leadmu_udsg");

	//leading muon asociated to a jet
        histo[8]    = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_no_flavor_"+tagger);
        histo[9]    = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_b_"+tagger);
        histo[10]   = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_c_"+tagger);
        histo[11]   = (TH1D*) gDirectory->Get("BTagPATClosure/ptrel/jet_ptrel_udsg_"+tagger);

	cv->Divide(2,2);
	for ( size_t i = 0; i < 4; i++ ){
		cv->cd(1+i);
		histo[i]->GetXaxis()->SetTitle("p_{T}rel muon tracks in Jet");
		legend[i] -> AddEntry(histo[i],"p_{T}rel all #mu","l");
		//gPad->SetLogy();
		//gPad->SetLogx();
                histo[i]->Rebin(pTrelRebin);
		histo[i]->Draw("LE");
		histo[i+4]->SetLineColor(2);
		legend[i] -> AddEntry(histo[i+4],"p_{T}rel leading #mu","l");
                histo[i+4]->Rebin(pTrelRebin);
		histo[i+4]->Draw("Same");
                histo[i+8]->Rebin(pTrelRebin);
		histo[i+8]->SetLineColor(4);
		legend[i] -> AddEntry(histo[i+8],"p_{T}rel leading #mu "+labeltag,"l");
		histo[i+8]->Draw("Same");
		legend[i] ->Draw();
	}
	cv->Update();
	cv->cd(0);
       // Draw a global picture title
	title = new TPaveLabel(0.35,0.49,0.70,0.52,"muons associated to a Jet");
        title->SetFillColor(16);
        title->Draw();
        cv->Update();
	cv-> Print ("BTagPATptrel_lead_muon.eps");
	
}
