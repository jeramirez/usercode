//in batch 
//root -b -q eff_vs_pt_tagger.C\(\"btagpatanalyzerpy.root\",\"ttbar\",2\)
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void eff_vs_pt_tagger(TString filename="btagpatanalyzerpy_ttbar.root", TString ttbar="ttbar",size_t taggerchosen=0) {

	
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
//        TH1D* eff[dim];
        TGraphAsymmErrors* eff[dim];
        TH1D* histoun[dim];
        TH1D* histotagun[dim];
//        TH1D* effun[dim];
        TGraphAsymmErrors* effun[dim];

	TCanvas *cv = new TCanvas("cv","cv",700,700);
	TMultiGraph *mg =new TMultiGraph();
	TLegend *legend0 = new TLegend(0.58,0.12,0.88,0.32);
	
	for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {

		TString tag    = "g"+taggers[itagger]+"_udsg";
		histo[itagger]    = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_b");
		histotag[itagger] = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_b_tagged");
//		eff[itagger]->SetName("Eff b "+taggers[itagger]);
//		eff[itagger] = (TH1D*) histo[itagger]->Clone();
//		eff[itagger]->Divide(histotag[itagger],histo[itagger],1,1,"B");      
		eff[itagger] = new TGraphAsymmErrors(histotag[itagger],histo[itagger]);
		eff[itagger]->GetXaxis()->SetTitle("Pt(GeV/c)");
		eff[itagger]->GetYaxis()->SetTitle("b-Eff");
		eff[itagger]->SetMaximum(1.05);
		eff[itagger]->SetMinimum(0.);

		histoun[itagger]    = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_uncorr_b");
		histotagun[itagger] = (TH1D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_pt_uncorr_b_tagged");
//		effun[itagger] = (TH1D*) histoun[itagger]->Clone();
//		effun[itagger]->SetName("Eff b "+taggers[itagger]);
//		effun[itagger]->Divide(histotagun[itagger],histoun[itagger],1,1,"BayesDivide");      
		effun[itagger]= new TGraphAsymmErrors(histotagun[itagger],histoun[itagger]);		
		effun[itagger]->SetLineColor(2);
	}
	cv->Divide(1,1);
        size_t itagger = taggerchosen;
		cv->cd(1);
		legend0 -> AddEntry(eff[itagger],"p_{T} corrected","l");
		eff[itagger]->Draw("AP");
		legend0 -> AddEntry(effun[itagger],"p_{T} uncorrected","l");
		effun[itagger]->Draw("P Same");
		legend0 ->Draw();
	cv->Update();
	cv->cd(0);
       // Draw a global picture title
	title = new TPaveLabel(0.35,0.97,0.55,1.00,"bEff vs Pt "+ttbar);
        title->SetFillColor(16);
        title->Draw();
        cv->Update();
	cv-> Print ("BTagPATeff_pt_"+taggers[taggerchosen]+"_"+ttbar+".eps");
	
}
