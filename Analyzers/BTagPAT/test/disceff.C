//para batch 
//root -b -q disceff.C\(\"btagpatanalyzer.root\"\)
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void disceff(TString filename) {
	gROOT->SetStyle("Plain");
	
	TString cmssw;
	// 167
	
	cmssw = "$3.1.0_pre9$";
		
	TFile *f = TFile::Open(filename);

	std::vector< TString > taggers;
	taggers.push_back( "TC2" );
	taggers.push_back( "TC3" );
	taggers.push_back( "TP" );
	taggers.push_back( "BTP" );
	taggers.push_back( "SSV" );
	taggers.push_back( "CSV" );
	taggers.push_back( "MSV" );
	taggers.push_back( "SMT" );
	taggers.push_back( "SETbyIP3d" );
	taggers.push_back( "SETbyPt" );
	taggers.push_back( "SMTbyIP3d" );
	taggers.push_back( "SMTbyPt" );

	std::vector< TString > discriminators;
	for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {
		discriminators.push_back( "disc"+taggers[itagger]+"_udsg" );
	}
//	discriminators.push_back( "discTC3_udsg" );
//	discriminators.push_back( "discTP_udsg" );

	const int dim=taggers.size();
        TCanvas *cv[dim];
	TMultiGraph* mg[dim];
	for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {
		TString tag = "g"+taggers[itagger]+"_udsg";
		TString tagb = "g"+taggers[itagger]+"_b";
		TString tagc = "g"+taggers[itagger]+"_c";
		cv[itagger] = new TCanvas("cv_"+taggers[itagger],"cv_"+taggers[itagger],700,700);
		TLegend *legend0 = new TLegend(0.68,0.70,0.88,0.90);

		TGraphErrors *agraph = (TGraphErrors*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/"+taggers[itagger]+"/"+tag);
		TGraphErrors *bgraph = (TGraphErrors*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/"+taggers[itagger]+"/"+tagb);
		TGraphErrors *cgraph = (TGraphErrors*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/"+taggers[itagger]+"/"+tagc);

		TGraph *dgraph = (TGraph*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/"+taggers[itagger]+"/"+discriminators[itagger]);
		TGraph *bvsdgraph = new TGraph(dgraph->GetN(),dgraph->GetY(),bgraph->GetY());
		TGraph *cvsdgraph = new TGraph(dgraph->GetN(),dgraph->GetY(),cgraph->GetY());
		TGraph *lvsdgraph = new TGraph(dgraph->GetN(),dgraph->GetY(),agraph->GetY());
		TGraph *udsgvsdgraph = new TGraph(dgraph->GetN(),dgraph->GetY(),dgraph->GetX());
		dgraph->Sort();
//		udsgvsdgraph->Sort();
//		udsgvsdgraph->SetLineColor(1);
//              legend0 -> AddEntry(udsgvsdgraph,"control","l");
		lvsdgraph->Sort();
		lvsdgraph->SetLineColor(2);
                legend0 -> AddEntry(lvsdgraph,tag,"l");
		cvsdgraph->Sort();
		cvsdgraph->SetLineColor(3);
                legend0 -> AddEntry(cvsdgraph,tagc,"l");
		bvsdgraph->Sort();
		bvsdgraph->SetLineColor(4);
                legend0 -> AddEntry(bvsdgraph,tagb,"l");
		mg[itagger]= new TMultiGraph();
//		mg[itagger]->Add(udsgvsdgraph);
		mg[itagger]->Add(lvsdgraph);
		mg[itagger]->Add(cvsdgraph);
		mg[itagger]->Add(bvsdgraph);
//		mg[itagger]->Add(dgraph);

		cv[itagger]->cd(1);
		mg[itagger]->Draw("ALP");
		mg[itagger]->GetYaxis()->SetTitle("eff");
		mg[itagger]->GetXaxis()->SetTitle("discriminant");
		legend0 -> Draw();
		cv[itagger]->Update();
	        cv[itagger]->cd(0);
		cv[itagger]-> Print ("BTagPATeff_vs_disc"+taggers[itagger]+".eps");
		cv[itagger]-> Print ("BTagPATeff_vs_disc"+taggers[itagger]+".ps");

		TGraphErrors *g = new TGraphErrors(agraph->GetN(),agraph->GetY(),agraph->GetX(),agraph->GetEY(),agraph->GetEX());
		g->Sort();
		g->SetLineColor(itagger+1);

		std::cout << " Tagger: " << tag << std::endl;
		std::cout << " Loose(10%): " << " cut > " << std::setprecision(4) << dgraph->Eval(0.1) << " b-eff = " << g->Eval(0.1) << std::endl;
		std::cout << " Medium(1%):  " << " cut > " << std::setprecision(4) << dgraph->Eval(0.01) << " b-eff = " << g->Eval(0.01) << std::endl;
		std::cout << " Tight(0.1%) = " << " cut > " << std::setprecision(4) << dgraph->Eval(0.001) << " b-eff = " << g->Eval(0.001) << std::endl;

	}//end for
 
	
}
