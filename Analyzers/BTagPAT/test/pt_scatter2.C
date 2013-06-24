//in batch 
//root -b -q pt_scatter2.C\(\"btagpatanalyzerpy.root\"\)
#include "TString.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void pt_scatter2(TString filename="btagpatanalyzerpy_ttbar.root") {

	
	TString cmssw;
	// 210pre6
	
	cmssw = "$2.1.0_pre9$";
	
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

	const int dim=taggers.size();
        TH2D* histob[dim];
        TH2D* histoc[dim];
        TH2D* histol[dim];
        TH2D* histonf[dim];

	TCanvas *cv      = new TCanvas("cv","cv",800,800);
	TLegend *legend0 = new TLegend(0.68,0.12,0.88,0.32);
        double x1=15.0,y1=0.0,x2=220.0,y2=150.0;
        double m=(y2-y1)/(x2-x1);
        double b=y1-m*x1;
        TLine *line1     = new TLine(x1,y1,x2,y2);
	
	for ( size_t itagger = 0; itagger < taggers.size(); ++itagger ) {

		TString tag    = "g"+taggers[itagger]+"_udsg";
		histob[itagger] = (TH2D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_scatter_pt_b");
		histoc[itagger] = (TH2D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_scatter_pt_c");
		histol[itagger] = (TH2D*) gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_scatter_pt_udsg");
		histonf[itagger] = (TH2D*) 
gDirectory->Get("BTagPATAnalyzer"+taggers[itagger]+"/jet_scatter_pt_no_flavor");
                histob[itagger]->GetXaxis()->SetTitle("Corrected pt");
                histoc[itagger]->GetXaxis()->SetTitle("Corrected pt");
                histob[itagger]->GetYaxis()->SetTitle("Uncorrected pt");
                histol[itagger]->GetYaxis()->SetTitle("Uncorrected pt");

                
//                histoc[itagger]->SetHistFillColor(2);		
	}
        line1->SetLineColor(2);
//	TH1D *scatterb = histob[0]->ProjectionX(" ",1,histob[0]->GetNbinsX(),"e");
	TProfile *scatterb = histob[0]->ProfileX();
	TProfile *scatterc = histoc[0]->ProfileX();
	TProfile *scatterl = histol[0]->ProfileX();
	TProfile *scattern = histonf[0]->ProfileX();
	cv->Divide(2,2);
        cv->cd(1); 
//        histob[0]->Fit("pol1","LM");
        scatterb->Fit("pol1","ELNQ");
	pol1->SetLineColor(4);
        scatterb->Fit("pol1","EL");
	gStyle->SetOptFit(0111); //Style Show fit
        histob[0]->Draw();
        scatterb->Draw("Same");
        line1->Draw("Same");
        cv->cd(2); 
        scatterc->Fit("pol1","EL");
	gStyle->SetOptFit(0111); //Style Show fit
	histoc[0]->Draw();
        scatterc->Draw("Same");
        line1->Draw("Same");
        cv->cd(3); 
        scatterl->Fit("pol1","ELQ");
	gStyle->SetOptFit(0111); //Style Show fit
	histol[0]->Draw();
        scatterl->Draw("Same");
        line1->Draw("Same");
        cv->cd(4); 
        scattern->Fit("pol1","ELQ");
	gStyle->SetOptFit(0111); //Style Show fit
	histonf[0]->Draw();
        scattern->Draw("Same");
        line1->Draw("Same");
	cv->Update();
	cv->cd(0);
       // Draw a global picture title
//	title = new TPaveLabel(0.15,0.97,0.35,1.00,"Pt uncorr(y) vs Pt corr(x)");
//        title->SetFillColor(16);
//        title->Draw();
        cv->Update();
        //Draw line parameters
        std::ofstream salida("BTagPATpt_scatter2.txt");
        std::ostringstream salida2;
        string sign=" + "; if (b<0) sign=" - ";
        salida << "m=" << m <<" ,b=" << b << endl;
        salida2 << "by eye y=" << m <<" x "<< sign << fabs(b) << endl;
	title = new TPaveLabel(0.07,0.82,0.27,0.85,salida2.str().c_str());
        title->SetFillColor(2);
        title->Draw();
	cv-> Print ("BTagPATpt_scatter2.eps");
	cv-> Print ("BTagPATpt_scatter2.ps");
	cv-> Print ("BTagPATpt_scatter2.png");
	
}
