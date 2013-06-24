//para batch 
//root -b -q pointsed.C\(\"btagpatanalyzer.root\"\)
//root -b -q pointsed.C\(\"btagpatanalyzerpy_ttbar.root\",\"\Corr30t2_\"\)
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void pointsed(TString filename = "btagpatanalyzerpy.root", TString check="") {
	TString cmssw;
	// 167
	
	cmssw = "$3.1.0_pre9$";
	
	TFile *f = TFile::Open(filename);
/////////////////////////////////////////////////////////////////////////////////////////
	std::vector< TString > checks;
	checks.push_back("");
	checks.push_back("Corr30_");
	checks.push_back("Corr30t0_");
	checks.push_back("Corr30t1_");
	checks.push_back("Corr30t2_");
	checks.push_back("Corr30t0_nConstituent_");
	checks.push_back("Corr30t1_nConstituent_");
	checks.push_back("Corr30t2_nConstituent_");
	checks.push_back("Uncor10_");
	checks.push_back("Uncor10t0_");
	checks.push_back("Uncor10t1_");
	checks.push_back("Uncor10t2_");
	checks.push_back("Uncor10t0_nConstituent_");
	checks.push_back("Uncor10t1_nConstituent_");
	checks.push_back("Uncor10t2_nConstituent_");
/////////////////////////////////////////////////////////////////////////////////////////
	std::vector< TString > taggers;
	std::vector< TString > labeltag;
        std::vector< int> colorlines;
	taggers.push_back( "TC2" );colorlines.push_back( 1 );labeltag.push_back( "TCHE");
	taggers.push_back( "TC3" );colorlines.push_back( 2 );labeltag.push_back( "TCHP");
	taggers.push_back( "TP"  );colorlines.push_back( 3 );labeltag.push_back( "JP");
	taggers.push_back( "SSV" );colorlines.push_back( 4 );labeltag.push_back( "SSV");
	taggers.push_back( "CSV" );colorlines.push_back( 5 );labeltag.push_back( "CSV");
	taggers.push_back( "MSV" );colorlines.push_back( 6 );labeltag.push_back( "CSVMVAB");
	taggers.push_back( "SMT" );colorlines.push_back( 8 );labeltag.push_back( "SMT");
	taggers.push_back( "BTP" );colorlines.push_back( 9 );labeltag.push_back( "JBP");
	taggers.push_back( "SMTbyIP3d" );colorlines.push_back( 11 );labeltag.push_back( "SMTByIP3d");
	taggers.push_back( "SMTbyPt" );colorlines.push_back( 12 );labeltag.push_back( "SMTByPt");
	taggers.push_back( "SETbyIP3d" );colorlines.push_back( 13 );labeltag.push_back( "SETByIP3d");
	taggers.push_back( "SETbyPt" );colorlines.push_back( 14 );labeltag.push_back( "SETByPt");
//	taggers.push_back( "IPM" );colorlines.push_back( 11 );labeltag.push_back( "IPMVAB");
//	taggers.push_back( "SMNIPT" );colorlines.push_back( 12 );labeltag.push_back( "SMTNoIP");

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//ed	for ( size_t ichks = 0; ichks < checks.size(); ++ichks ) {
//ed		TString check = checks[ichks];
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
        TString PerfTitle= check;
        if (check == "") PerfTitle="uncorr pt>30";
	if (check == "Corr30_") PerfTitle="pt>30";
	if (check == "Corr30t0_") PerfTitle="pt>30, #tracks > 0";
	if (check == "Corr30t1_") PerfTitle="pt>30, #tracks > 1";
	if (check == "Corr30t2_") PerfTitle="pt>30, #tracks > 2";
	if (check == "Corr30t0_nConstituent_") PerfTitle="pt>30, #tracks > 0 & nConstituents >1";
	if (check == "Corr30t1_nConstituent_") PerfTitle="pt>30, #tracks > 1 & nConstituents >1";
	if (check == "Corr30t2_nConstituent_") PerfTitle="pt>30, #tracks > 2 & nConstituents >1";
	if (check == "Uncor10_") PerfTitle="uncorr pt>10";
	if (check == "Uncor10t0_") PerfTitle="uncorr pt>10 #tracks > 0";
	if (check == "Uncor10t1_") PerfTitle="uncorr pt>10 #tracks > 1";
	if (check == "Uncor10t2_") PerfTitle="uncorr pt>10 #tracks > 2";
	if (check == "Uncor10t0_nConstituent_") PerfTitle="uncorr pt>10 #tracks > 0 & nConstituents >1";
	if (check == "Uncor10t1_nConstituent_") PerfTitle="uncorr pt>10 #tracks > 1 & nConstituents >1";
	if (check == "Uncor10t2_nConstituent_") PerfTitle="uncorr pt>10 #tracks > 2 & nConstituents >1";

	std::vector< TString > discriminators;
        for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {
  		discriminators.push_back( check+"disc"+taggers[itagger]+"_udsg" );
	}
//	discriminators.push_back( "discTC3_udsg" );
//	discriminators.push_back( "discTP_udsg" );

	TCanvas *cv = new TCanvas("cv","cv",900,900);
	TMultiGraph *mg =new TMultiGraph();
	TLegend *legend0 = new TLegend(0.68,0.12,0.88,0.32);
        // inverted axis
	TMultiGraph *mginv =new TMultiGraph();
	TLegend *legend1 = new TLegend(0.65,0.18,0.85,0.48);
	std::ofstream salida("BTagPATop3"+check+".txt");
	
	for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {

		TString tag    = check+"g"+taggers[itagger]+"_udsg";
		TGraphErrors *agraph = (TGraphErrors*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/"+taggers[itagger]+"/"+tag);

		
		TGraph *dgraph = (TGraph*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/"+taggers[itagger]+"/"+discriminators[itagger]);
//		TGraph *udsgvsdgraph = new TGraph(dgraph->GetN(),dgraph->GetY(),dgraph->GetX());
		dgraph->Sort();
//		udsgvsdgraph->Sort();

		TGraphErrors *g = new TGraphErrors(agraph->GetN(),agraph->GetY(),agraph->GetX(),agraph->GetEY(),agraph->GetEX());
		g->Sort();
		g->SetLineColor(itagger+1);
                legend0 -> AddEntry(g,taggers[itagger],"l");
		mg->Add(g);

                //inverted axis
		TGraphErrors *ginv = new TGraphErrors(agraph->GetN(),agraph->GetX(),agraph->GetY(),agraph->GetEX(),agraph->GetEY());
		ginv->Sort();
		ginv->SetLineColor(colorlines[itagger]);
		ginv->SetMarkerStyle(8);
		ginv->SetMarkerColor(colorlines[itagger]);
                legend1 -> AddEntry(ginv,labeltag[itagger],"l");
		mginv->Add(ginv);

		std::cout << " Tagger: " << tag << std::endl;
		std::cout << " Loose(10%): " << " cut > " << std::setprecision(4) << dgraph->Eval(0.1) << " b-eff = " << g->Eval(0.1) << std::endl;
		std::cout << " Medium(1%):  " << " cut > " << std::setprecision(4) << dgraph->Eval(0.01) << " b-eff = " << g->Eval(0.01) << std::endl;
		std::cout << " Tight(0.1%) = " << " cut > " << std::setprecision(4) << dgraph->Eval(0.001) << " b-eff = " << g->Eval(0.001) << std::endl;
                salida << " " << taggers[itagger] << endl;
	        salida << " #bin \t b-eff  \t err-beff \t  non-beff \t err-non-b"<< endl;
        	for ( size_t i=0;i< agraph->GetN();i++){ 
			salida	<< "  " << (i+1) << " \t "
				<< agraph->GetX()[i] <<" \t " 
				<< agraph->GetEX()[i]<<" \t "
				<< agraph->GetY()[i] <<" \t "
				<< agraph->GetEY()[i]
                       		<< " "<< endl;
		}

	}
//	cv->SetLogx();
	mg->Draw("ALP");
	mg->GetYaxis()->SetTitle("b-eff");
	mg->GetXaxis()->SetTitle("udsg-mistagging");
	legend0 -> Draw();
        cv-> SetGrid();
	cv-> Print ("BTagPATop"+check+".eps");
	cv-> Print ("BTagPATop"+check+".png");
	cv->cd(0);
	cv->SetLogx();
	mg->Draw("ALP");
	legend0 -> Draw();
	cv-> Print ("BTagPATop1"+check+".eps");
	cv-> Print ("BTagPATop1"+check+".png");
        //inverted axis
	TCanvas *cvinv = new TCanvas("cvinv","cvinv",700,700);
        cvinv->SetLogy(0);
	mginv->Draw("ALP");
	mginv->GetXaxis()->SetTitle("b-eff");
	mginv->GetYaxis()->SetTitle("udsg-mistagging");
	legend1 -> Draw();
        cvinv->SetGrid();
	cvinv-> Print ("BTagPATop2"+check+".eps");
	cvinv-> Print ("BTagPATop2"+check+".png");
	cvinv->cd(0);
//	cvinv->Update();
//	cvinv->SetLogx();
//	mginv->SetXmin(10^-4);
//      axis range using a histogram helper
    TH2D*hpx = new TH2D("hpx",PerfTitle,200,0.0,1.005,200,0.00002,1.05);
        hpx->SetStats(kFALSE);
		hpx->Draw();
	hpx->GetXaxis()->SetTitle("b-eff");
	hpx->GetYaxis()->SetTitle("udsg-mistagging");
	cvinv->SetLogy(1);
	mginv->Draw("AP");
	legend1 -> Draw();
	cvinv-> Print ("BTagPATop3"+check+".eps");
	cvinv-> Print ("BTagPATop3"+check+".png");
	//ed}
}
