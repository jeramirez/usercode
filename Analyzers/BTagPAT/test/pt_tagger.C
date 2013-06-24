//in batch 
//root -b -q pt_tagger.C\(\"btagpatanalyzerpy.root\",\"ttbar\",2\)
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void pt_tagger(TString filename="btagpatanalyzerpy_ttbar.root", TString ttbar="ttbar",size_t taggerchosen=0) {

	
	TString cmssw;
	// 210pre6
	
	cmssw = "$2.1.0_pre6$";
	
	TFile *f = TFile::Open(filename);

	std::vector< TString > taggers;
	taggers.push_back( "TC2" );
	taggers.push_back( "TC3" );
	taggers.push_back( "TP" );
	taggers.push_back( "SSV" );
	taggers.push_back( "CSV" );
	taggers.push_back( "MSV" );
	taggers.push_back( "SET" );
	taggers.push_back( "SMT" );

        if (taggers.size() -1 < taggerchosen ) taggerchosen=0; //protect against big values
	const int dim=taggers.size();
        TH1D* histo[dim];
        TH1D* histotag[dim];
        TH1D* eff[dim];
        TH1D* histoun[dim];
        TH1D* histotagun[dim];
        TH1D* effun[dim];

	TCanvas *cv = new TCanvas("cv","cv",700,700);
	TMultiGraph *mg =new TMultiGraph();
	TLegend *legend0 = new TLegend(0.68,0.72,0.98,0.92);
	
	for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {

		TString tag    = "g"+taggers[itagger]+"_udsg";
		histo[itagger]    = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_b");
		histotag[itagger] = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_b_tagged");
		histotag[itagger]->GetXaxis()->SetTitle("Pt(GeV/c^2)");
		histotag[itagger]->GetYaxis()->SetTitle("# b");

		histoun[itagger]    = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_uncorr_b");
		histotagun[itagger] = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_uncorr_b_tagged");
		histotagun[itagger]->SetLineColor(2);
		
	}
	cv->Divide(1,1);
        size_t itagger = taggerchosen;
	THStack hs;
		cv->cd(1);
		legend0 -> AddEntry(histotag[itagger],"p_{T} corrected","l");
		hs.Add(histotag[itagger]);
//		histotag[itagger]->Draw("LE");
		legend0 -> AddEntry(histotagun[itagger],"p_{T} uncorrected","l");
		hs.Add(histotagun[itagger]);
//		histotagun[itagger]->Draw("Same");
		hs.Draw("nostack");
		legend0 ->Draw();
	cv->Update();
	cv->cd(0);
       // Draw a global picture title
	title = new TPaveLabel(0.35,0.97,0.55,1.00,"#b vs Pt "+ttbar);
        title->SetFillColor(16);
        title->Draw();
        cv->Update();
	cv-> Print ("BTagPATpt_"+taggers[taggerchosen]+"_"+ttbar+".eps");
	
}
