//in batch 
//root -b -q eff_vs_taggability.C\(\"btagpatanalyzerpy.root\"\)
//root -b -q eff_vs_taggability.C\(\"btagpatanalyzerpy_qcd.root\",\"qcd\"\)
#include "TString.h"
#include "TFile.h"
#include "TCanvas.h"

#include <vector>
#include <iostream>
#include <iomanip>

void eff_eta_taggability(TString filename="btagpatanalyzerpy_ttbar.root", TString data="ttbar") {

	
	TFile *f = TFile::Open(filename);
        const int dim=27; 
        TH1D* histo[dim];
        TH1D* histotag[dim];
	TH1D* eff[dim];
        std::vector< TString > quarks;
        quarks.push_back( "b" );
        quarks.push_back( "c" );
        quarks.push_back( "udsg" );

        std::vector< TString > tracks;
        tracks.push_back( "#tracks>0" );
        tracks.push_back( "#tracks>1" );
        tracks.push_back( "#tracks>2" );

        std::vector< TString > qtaggable;
        for ( size_t itrk=0; itrk < tracks.size(); ++itrk){
		for ( size_t iquark=0; iquark< quarks.size(); ++iquark){
			qtaggable.push_back( quarks[iquark]+" "+tracks[itrk]);
		}
	}

	TCanvas *cv = new TCanvas("cv","cv",700,1000);
	TLegend *legend = new TLegend(0.18,0.72,0.28,0.82);
	TPaveLabel *titlerow[3];
	titlerow[0] = new TPaveLabel(-2.50,0.60,-0.50,0.65,"p < 30");
	titlerow[1] = new TPaveLabel(-2.50,0.905,-0.5,0.92,"30< p < 50");
	titlerow[2] = new TPaveLabel(-2.50,0.938,-0.5,0.948,"p > 50");
	
        histotag[0] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable0_030");
        histo[0]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b030");
        histotag[1] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable0_030");
        histo[1]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c030");
        histotag[2] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable0_030");
        histo[2]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg030");

        histotag[3] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable1_030");
        histo[3]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b030");
        histotag[4] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable1_030");
        histo[4]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c030");
        histotag[5] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable1_030");
        histo[5]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg030");

        histotag[6] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable2_030");
        histo[6]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b030");
        histotag[7] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable2_030");
        histo[7]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c030");
        histotag[8] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable2_030");
        histo[8]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg030");

        histotag[9] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable0_3050");
        histo[9]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b3050");
        histotag[10] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable0_3050");
        histo[10]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c3050");
        histotag[11] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable0_3050");
        histo[11]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg3050");

        histotag[12] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable1_3050");
        histo[12]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b3050");
        histotag[13] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable1_3050");
        histo[13]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c3050");
        histotag[14] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable1_3050");
        histo[14]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg3050");

        histotag[15] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable1_3050");
        histo[15]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b3050");
        histotag[16] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable1_3050");
        histo[16]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c3050");
        histotag[17] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable1_3050");
        histo[17]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg3050");

        histotag[18] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable0_50");
        histo[18]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b50");
        histotag[19] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable0_50");
        histo[19]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c50");
        histotag[20] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable0_50");
        histo[20]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg50");

        histotag[21] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable1_50");
        histo[21]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b50");
        histotag[22] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable1_50");
        histo[22]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c50");
        histotag[23] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable1_50");
        histo[23]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg50");

        histotag[24] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b_taggable2_50");
        histo[24]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_b50");
        histotag[25] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c_taggable2_50");
        histo[25]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_c50");
        histotag[26] = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg_taggable2_50");
        histo[26]    = (TH1D*) gDirectory->Get("BTagPATAnalyzerTC2/jet_eta_udsg50");

// Compute efficiencies
	for ( size_t i = 0; i < dim; ++i ){
                eff[i] = (TH1D*) histo[i]->Clone();
                eff[i]->SetName("Eff "+qtaggable[i%9]+" vs #eta");
                eff[i]->Divide(histotag[i],histo[i],1,1,"B");
		eff[i]->GetXaxis()->SetTitle("#eta");
		eff[i]->GetYaxis()->SetTitle("Eff");
	}


//We are ready to plot any histogram
	cv->Divide(3,3);
	for ( size_t irow= 0; irow < 3; ++irow){
  		for ( size_t icol = 0; icol < 3 ; ++icol ){
			cv->cd(1 + icol + 3*irow);
			eff[1 + 3*icol + 9*irow]->SetLineColor(2); //charm eff
			eff[2 + 3*icol + 9*irow]->SetLineColor(3); //udsg eff
			eff[3*icol+ 9*irow]->SetMinimum(0.6 + 0.45*sqrt(irow) - 0.15*irow);
			eff[3*icol + 9*irow]->Draw("PE0");
			eff[1 +3*icol + 9*irow]->Draw("Same");
			eff[2 +3*icol + 9*irow]->Draw("Same");
//		eff[irow + 3*icol]->Draw("PE0");
//                histotag[i]->SetLineColor(3);
//		histo[i]->Draw("LE0");
//		histotag[i]->Draw("Same");
		}
 		legend -> AddEntry(eff[irow],quarks[irow],"l");
		cv->Update();
		cv->cd(1 + 3*irow);
		titlerow[irow]->Draw();
	}
	cv->Update();
	cv->cd(0);
       // Draw a global picture title
	title = new TPaveLabel(0.00,0.97,0.20,1.00,"eff taggability");
        title->SetFillColor(16);
        title->Draw();
        cv->Update();
	legend ->Draw();
	cv-> Print ("BTagPAT_eff_eta_taggability0_"+data+".eps");
	cv-> Print ("BTagPAT_eff_eta_taggability0_"+data+".ps");
	cv-> Print ("BTagPAT_eff_eta_taggability0_"+data+".png");
	
}
