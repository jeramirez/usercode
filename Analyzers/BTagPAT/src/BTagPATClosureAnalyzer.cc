// -*- C++ -*-
//
// Class:      BTagPATClosureAnalyzer
// 
/**\class BTagPATClosureAnalyzer Analyzers/BTagPAT/src/BTagPATClosureAnalyzer.cc

 Description: <A CMSSW analyzer for PAT objects>

 Implementation:
 
 this analyzer shows how to loop over PAT objects and play with it. 
*/
//
// Original Author:  J.E. Ramirez
//
// $Id: BTagPATClosureAnalyzer.cc,v 1.5 2008/07/29 15:32:41 jramirez Exp $


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"


#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/Common/interface/View.h"
#include <string>
//ed++
#include "Analyzers/BTagPAT/interface/BTagPATClosureAnalyzer.h"

#include "Math/GenVector/VectorUtil.h"
#include "TVector3.h"
//ed--

//
// constructors and destructor
//
BTagPATClosureAnalyzer::BTagPATClosureAnalyzer(const edm::ParameterSet& iConfig):
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag"))
  ,muonLabel_(iConfig.getUntrackedParameter<edm::InputTag>("muonTag"))
  ,PtrelHistos(iConfig)
  ,DeltaRHistos(100,0.0,1.0)
  ,Chi2muHistos(110,-0.5,10.5)
  ,BTagger(iConfig.getParameter< edm::ParameterSet >("BJetOperatingPoints")) //this gives a warning
  ,nmuonhistos(11,-0.5,10.5)
  ,nmuonhits(31,-0.5,30.5)
  ,scatterjetpt(11,50,0.,5.)
  ,scatterjeteta(10,50,0.,5.)
  ,taggerslist(iConfig.getParameter< edm::ParameterSet >("BJetOperatingPoints").getParameter<std::vector<std::string> >("bdiscriminators"))
{
   //now do what ever initialization is needed
   edm::ParameterSet PatBjet_(iConfig.getParameter< edm::ParameterSet >("BjetTag"));
   BTagverbose       = PatBjet_.getUntrackedParameter<bool>("verbose",false);
   BTagtagger_       = PatBjet_.getParameter<std::string>("tagger");                //for away jet
   BTagpurity_       = PatBjet_.getParameter<std::string>("purity");                //for away jet

   // jet cuts
   MinJetPt_  = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinPt");
   MaxJetEta_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MaxEta");
   MaxDeltaR_ = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MaxDeltaR");
   MinPtRel_  = iConfig.getParameter<edm::ParameterSet>("jetcuts").getParameter<double>("MinPtRel");
 
   // muon cuts
   MinMuHits_= iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<int>("MinNHits");
   MinMuPt_  = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MinMuonPt");
   MaxMuEta_ = iConfig.getParameter<edm::ParameterSet>("muoncuts").getParameter<double>("MaxMuonEta"); 

  //get a list of taggers and fill a map of operating points and taggers
   std::vector<double>
   opLoose(iConfig.getParameter< edm::ParameterSet >("BJetOperatingPoints").getParameter<std::vector<double> >("discCutLoose"));
   std::vector<double>
   opMedium(iConfig.getParameter< edm::ParameterSet >("BJetOperatingPoints").getParameter<std::vector<double> >("discCutMedium"));
   std::vector<double>
   opTight(iConfig.getParameter< edm::ParameterSet >("BJetOperatingPoints").getParameter<std::vector<double> >("discCutTight"));
   std::cout << "we have " << taggerslist.size() << " taggers:" << std::endl; 
   for (unsigned int i=0; i<taggerslist.size(); i++){
     std::cout << i << " Loose OP(" << taggerslist[i] << ")= " << opLoose[i] << std::endl;
     opCut["Loose" ][taggerslist[i]] = opLoose[i];
     opCut["Medium"][taggerslist[i]] = opMedium[i];
     opCut["Tight" ][taggerslist[i]] = opTight[i];
   }
}


BTagPATClosureAnalyzer::~BTagPATClosureAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BTagPATClosureAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   // first: get all objects from the event.
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

   edm::Handle<edm::View<pat::Jet> > jetHandle;
   iEvent.getByLabel(jetLabel_,jetHandle);
   edm::View<pat::Jet> jets = *jetHandle;

   edm::Handle<edm::View<pat::Muon> > muonHandle;
   iEvent.getByLabel(muonLabel_,muonHandle);
   edm::View<pat::Muon> muons = *muonHandle;

   //Loop over jets
   for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
      //Jet requirements cuts
      double absjeteta = fabs(jet_iter->eta()); 
      double jetpt     = jet_iter->pt();
      if(
             jetpt     > MinJetPt_ 
          && absjeteta < MaxJetEta_
      ){

         int flavor   = jet_iter->partonFlavour();
         std::map <std::string,std::map<std::string,bool> > istagged;
         //method 1 using helper funtion
         //for (unsigned int i=0; i<taggerslist.size(); i++){
         //  istagged["Loose" ][taggerslist[i]] = BTagger.IsbTag(*jet_iter,"Loose",taggerslist[i]);
         //  istagged["Medium"][taggerslist[i]] = BTagger.IsbTag(*jet_iter,"Medium",taggerslist[i]);
         //  istagged["Tight" ][taggerslist[i]] = BTagger.IsbTag(*jet_iter,"Tight",taggerslist[i]);
         //}
         // method 2 using PairDiscri (faster?)
         for(unsigned int i=0; i!=jet_iter->getPairDiscri().size(); i++){
            istagged["Loose" ][jet_iter->getPairDiscri()[i].first] = 
               (jet_iter->getPairDiscri()[i].second > opCut["Loose" ][jet_iter->getPairDiscri()[i].first]);
            istagged["Medium"][jet_iter->getPairDiscri()[i].first] = 
               (jet_iter->getPairDiscri()[i].second > opCut["Medium"][jet_iter->getPairDiscri()[i].first]);
            istagged["Tight" ][jet_iter->getPairDiscri()[i].first] = 
               (jet_iter->getPairDiscri()[i].second > opCut["Tight" ][jet_iter->getPairDiscri()[i].first]);           
         }
         if (BTagverbose){
	    std::cout << "-------  Start Loop over Available discriminants--------" << std::endl;
            //std::vector<std::pair<std::string,float> > pairDiscri=jet_iter->getPairDiscri();
            for(unsigned int i=0; i!=jet_iter->getPairDiscri().size(); i++){
               std::cout << "Tagger= " << jet_iter->getPairDiscri()[i].first
                         << "Discriminant= " << jet_iter->getPairDiscri()[i].second
                         << std::endl; 
            }
	    std::cout << "-------  End Loop over Available discriminants--------" << std::endl;
         }


         //Get jet momentum vector
         TVector3 pjet(jet_iter->p4().Vect().X(),jet_iter->p4().Vect().Y(),jet_iter->p4().Vect().Z());
         //Loop over muons 
         double ptrel=0.;  
         double mupt_highest=0; 
         double ptrelleadmuon=0.; double deltaRleadmuon=0.; double normChi2leadmuon=0.;
         int nmuons=0; int nmuonhitsleadmu=0;
         if (BTagverbose){
	    std::cout << "-------  Start Muon Loop --------" << std::endl;
         }
         for(edm::View<pat::Muon>::const_iterator muon_iter = muons.begin(); muon_iter!=muons.end(); ++muon_iter){
            // muon cuts
            //double Chi2 = muon_iter->track()->chi2();
            double normChi2 = muon_iter->track()->normalizedChi2();
            int nmuhits     = muon_iter->track()->numberOfValidHits();
            if (
                 muon_iter->pt()        > MinMuPt_ 
              && fabs(muon_iter->eta()) < MaxMuEta_
              && nmuhits                > MinMuHits_
            ){
	       if (BTagverbose){
	          std::cout << "Muon candidate:" 
                            << " Chi2/ndof = "   << normChi2
                            << " #valid hits = " << nmuhits
                            << std::endl;
	       }
               // /\R(jet,mu) < 0.7
	       double DeltaR = ROOT::Math::VectorUtil::DeltaR(jet_iter->p4().Vect(),muon_iter->p4().Vect());
               //Compute Ptrel by
               //Getting the muon momentum vector, then adding the jet momentum
               TVector3 pmu(muon_iter->p4().Vect().X(),muon_iter->p4().Vect().Y(),muon_iter->p4().Vect().Z());
               TVector3 pmujet = pjet + pmu;
               ptrel = pmu.Perp(pmujet); 
               if (DeltaR < MaxDeltaR_ && ptrel > MinPtRel_ ){
                  if (BTagverbose){
                    std::cout << "Delta R="<< DeltaR 
                              << " and "
                              << "Ptrel=" << ptrel
                              << std::endl;
                  }
                  nmuons++;
                  //.....all....muons....in....jet....
                  //fill ptrel, DeltaR, normalized Chi2, #hits histograms
                  PtrelHistos.Fill(ptrel,"all");
                  DeltaRHistos.Fill(DeltaR,"all");
                  Chi2muHistos.Fill(normChi2,"all");
                  nmuonhits.Fill(nmuhits,"all");
                  std::string flvstr="";
                  if (flavor ==  0  )flvstr="no_flavor";
                  if (flavor ==  5 || flavor ==  -5 )flvstr="b";
                  if (flavor ==  4 || flavor ==  -4 )flvstr="c";
                  if ((-4 < flavor && flavor < 4 && flavor != 0 )||(flavor == 21 || flavor == -21 ))flvstr="udsg";
                  if (flvstr !=  ""  ) {
                    PtrelHistos.Fill(ptrel,flvstr);
		    DeltaRHistos.Fill(DeltaR,flvstr);
                    Chi2muHistos.Fill(normChi2,flvstr);
                    nmuonhits.Fill(nmuhits,flvstr);
                  }

                  //save lead muon (highest pt) info
                  if (muon_iter->pt() > mupt_highest) {
                    mupt_highest  = muon_iter->pt();
                    ptrelleadmuon = ptrel;
                    deltaRleadmuon= DeltaR;
                    normChi2leadmuon=normChi2;
                    nmuonhitsleadmu =nmuhits;

                  }//end if lead
               }//end if DeltaR
            }//end if muon cuts
         }//end muon loop
         if (BTagverbose){
	    std::cout << "-------  End Muon Loop  --------" << std::endl;
         }
     

         //Loop to find away jet
         bool AwayJetIstagged=false;
         for(edm::View<pat::Jet>::const_iterator awayjet_iter = jets.begin(); awayjet_iter!=jets.end(); ++awayjet_iter){
            //Jet requirements cuts
            double absawayjeteta = fabs(awayjet_iter->eta());
            if ( 
                  awayjet_iter->pt() > 30      //  pt > 30
               && absawayjeteta  < 2.4         // abs (eta) < 2.4
               && awayjet_iter != jet_iter     // make sure is different from original jet
               && nmuons > 0                   // make sure original jet has an associated muon
               ){
                AwayJetIstagged = BTagger.IsbTag(*awayjet_iter,BTagpurity_,BTagtagger_); //AwayJet tag
            } 
         }//end loop away jet

         //Fill jet histograms
         nmuonhistos.Fill(nmuons,"all");
         if (nmuons > 0){
            PtrelHistos.Fill(ptrelleadmuon,"leadmu_all");            //ptrel distribution
	    DeltaRHistos.Fill(deltaRleadmuon,"leadmu_all");          //DeltaR distribution 
            Chi2muHistos.Fill(normChi2leadmuon,"leadmu_all");        //chi2/ndof distribution
            nmuonhits.Fill(nmuonhitsleadmu,"leadmu_all");            //nmuhits distribution

            scatterjetpt.SetPrefix("n_pT_");            
            scatterjetpt.Fill(jetpt,ptrelleadmuon,"all");
            scatterjetpt.SetPrefix("p_pT_");            
            if(AwayJetIstagged)scatterjetpt.Fill(jetpt,ptrelleadmuon,"all");

            scatterjeteta.SetPrefix("n_eta_");            
            scatterjeteta.Fill(absjeteta,ptrelleadmuon,"all");
            scatterjeteta.SetPrefix("p_eta_");            
            if(AwayJetIstagged)scatterjeteta.Fill(absjeteta,ptrelleadmuon,"all");


            for (unsigned int i=0; i<taggerslist.size(); i++){
               if (istagged["Loose"][taggerslist[i]]){
                  PtrelHistos.Fill(ptrelleadmuon,"all_"+taggerslist[i]);            
		  DeltaRHistos.Fill(deltaRleadmuon,"all_"+taggerslist[i]);

                  scatterjetpt.SetPrefix("ntag_pT_");            
                  scatterjetpt.Fill(jetpt,ptrelleadmuon,"all_"+taggerslist[i]);
                  scatterjetpt.SetPrefix("ptag_pT_");            
                  if(AwayJetIstagged)scatterjetpt.Fill(jetpt,ptrelleadmuon,"all_"+taggerslist[i]);

                  scatterjeteta.SetPrefix("ntag_eta_");            
                  scatterjeteta.Fill(absjeteta,ptrelleadmuon,"all_"+taggerslist[i]);
                  scatterjeteta.SetPrefix("ptag_eta_");            
                  if(AwayJetIstagged)scatterjeteta.Fill(absjeteta,ptrelleadmuon,"all_"+taggerslist[i]);
               }//end if is tagged?
            }//end loop taggers
         }
         //Fill jet histograms divided by flavor
         std::string flavorstr="";
         if (flavor ==  0 ) flavorstr="no_flavor";
         if (flavor ==  5 || flavor ==  -5 ) flavorstr="b";
         if (flavor ==  4 || flavor ==  -4 ) flavorstr="c";
         if ((-4 < flavor && flavor < 4 && flavor != 0 )||(flavor == 21 || flavor == -21 ))flavorstr="udsg";

         if (flavorstr !=  ""  ) {
            nmuonhistos.Fill(nmuons,flavorstr);
            if (nmuons > 0){
               PtrelHistos.Fill(ptrelleadmuon,"leadmu_"+flavorstr);            //ptrel distribution
	       DeltaRHistos.Fill(deltaRleadmuon,"leadmu_"+flavorstr);          //DeltaR distribution
               Chi2muHistos.Fill(normChi2leadmuon,"leadmu_"+flavorstr);        //chi2/ndof distribution
               nmuonhits.Fill(nmuonhitsleadmu,"leadmu_"+flavorstr);            //nmuhits distribution
               //pt vs ptrel
               scatterjetpt.SetPrefix("n_pT_");            
               scatterjetpt.Fill(jetpt,ptrelleadmuon,flavorstr);
               scatterjetpt.SetPrefix("p_pT_");            
               if(AwayJetIstagged)scatterjetpt.Fill(jetpt,ptrelleadmuon,flavorstr);
               //eta vs ptrel
               scatterjeteta.SetPrefix("n_eta_");            
               scatterjeteta.Fill(absjeteta,ptrelleadmuon,flavorstr);
               scatterjeteta.SetPrefix("p_eta_");            
               if(AwayJetIstagged)scatterjeteta.Fill(absjeteta,ptrelleadmuon,flavorstr);

               for (unsigned int i=0; i<taggerslist.size(); i++){
		 if (istagged["Loose"][taggerslist[i]]){                                                //Loose Tagged Jet
		     PtrelHistos.Fill(ptrelleadmuon,flavorstr+"_"+taggerslist[i]);                      //ptrel distribution
		     DeltaRHistos.Fill(deltaRleadmuon,flavorstr+"_"+taggerslist[i]);                    //DeltaR distribution
                     Chi2muHistos.Fill(normChi2leadmuon,flavorstr+"_"+taggerslist[i]);                  //chi2/ndof distribution
                     nmuonhits.Fill(nmuonhitsleadmu,"leadmu_"+flavorstr+"_"+taggerslist[i]);            //nmuhits distribution

                     scatterjetpt.SetPrefix("ntag_pT_");            
                     scatterjetpt.Fill(jetpt,ptrelleadmuon,flavorstr+"_"+taggerslist[i]);
                     scatterjetpt.SetPrefix("ptag_pT_");            
                     if(AwayJetIstagged)scatterjetpt.Fill(jetpt,ptrelleadmuon,flavorstr+"_"+taggerslist[i]);

                     scatterjeteta.SetPrefix("ntag_eta_");            
                     scatterjeteta.Fill(absjeteta,ptrelleadmuon,flavorstr+"_"+taggerslist[i]);
                     scatterjeteta.SetPrefix("ptag_eta_");            
                     if(AwayJetIstagged)scatterjeteta.Fill(absjeteta,ptrelleadmuon,flavorstr+"_"+taggerslist[i]);
                  }//end if is tagged?
               }//end loop taggers
            }
         }

      }//end if jet cuts
   }//end loop over jets


}
// ------------ method called once each job just before starting event loop  ------------
void 
BTagPATClosureAnalyzer::beginJob(const edm::EventSetup&)
{
  
  edm::Service<TFileService> fs;
  
  //Directories 
  TFileDirectory ptreldir = fs->mkdir("ptrel");
  TFileDirectory mudir = fs->mkdir("muDir");
  TFileDirectory muinjet = fs->mkdir("muon_in_jet");
  std::vector<TFileDirectory> muinjettags;

  //variable size bin
  double jetetabins[10]={0.0,0.25,0.50,0.75,1.0,1.25,1.50,1.75,2.0,2.5};
  double jetptbins[11]={30.,40.,50.,60.,70.,80.,90.,100.,120.,140.,230.};

  std::vector<std::string> flavor;
  flavor.push_back("all");
  flavor.push_back("no_flavor");
  flavor.push_back("b");
  flavor.push_back("c");
  flavor.push_back("udsg");

  for (unsigned int i=0; i<taggerslist.size(); i++){
    muinjettags.push_back(muinjet.mkdir(taggerslist[i]));
  }  

  for (size_t iflv = 0; iflv < flavor.size(); ++iflv){

     //number of associated trk muons per jet
     nmuonhistos.Set(flavor[iflv],mudir,"mutrkjet_","# mu trks/ jet "+flavor[iflv]);

     //#hits, chi2, Ptrel, DeltaR histograms (all associated muons)
     nmuonhits.Set(flavor[iflv],mudir,"muhitsjet_","# valid track hits/ mu "+flavor[iflv]);
     Chi2muHistos.Set(flavor[iflv],mudir,"muchi2jet_","#chi^2 #mu distribution "+flavor[iflv]);
     PtrelHistos.Set(flavor[iflv],mudir,"jet_ptrel_","p_{Trel} "+flavor[iflv]+" [GeV/c]");
     DeltaRHistos.Set(flavor[iflv],mudir,"mudeltaR_","#Delta R distribution "+flavor[iflv]);
 
     //#hits, chi2,Ptrel using leading muon only
     nmuonhits.Set("leadmu_"+flavor[iflv],ptreldir,"muhitsjet_","# valid track hits/ mu "+flavor[iflv]);
     Chi2muHistos.Set("leadmu_"+flavor[iflv],ptreldir,"muchi2jet_","#chi^2 #mu distribution "+flavor[iflv]);
     PtrelHistos.Set("leadmu_"+flavor[iflv],ptreldir,"jet_ptrel_","p_{Trel} "+flavor[iflv]+" [GeV/c]");
     DeltaRHistos.Set("leadmu_"+flavor[iflv],ptreldir,"mudeltaR_","#Delta R distribution "+flavor[iflv]);

     //pt(eta) histograms vs ptrel using leading muon (pt is corrected)
     //no info yet on away jet
     scatterjetpt.Set(flavor[iflv],muinjet,jetptbins,"n_pT_","muon-jet p_{T} vs ptrel");
     scatterjeteta.Set(flavor[iflv],muinjet,jetetabins,"n_eta_","muon-jet eta vs ptrel");
     // now include
     //tagged jet away info
     scatterjeteta.Set(flavor[iflv],muinjet,jetetabins,"p_eta_","muon-jet eta vs ptrel");
     scatterjetpt.Set(flavor[iflv],muinjet,jetptbins,"p_pT_","muon-jet p_{T} vs ptrel");


     //Loose tagged jet
     for (unsigned int i=0; i<taggerslist.size(); i++){    
        //#hits, chi2,Ptrel using leading muon only Tagged Jet
        nmuonhits.Set("leadmu_"+flavor[iflv]+"_"+taggerslist[i],ptreldir,"muhitsjet_","# valid track hits/ mu "+flavor[iflv]);
	Chi2muHistos.Set(flavor[iflv]+"_"+taggerslist[i],ptreldir,"muchi2jet_","#chi^2 #mu distribution "+flavor[iflv]);
        PtrelHistos.Set(flavor[iflv]+"_"+taggerslist[i],ptreldir,"jet_ptrel_","p_{Trel} "+flavor[iflv]+" [GeV/c]");
	DeltaRHistos.Set(flavor[iflv]+"_"+taggerslist[i],ptreldir,"mudeltaR_","#Delta R distribution "+flavor[iflv]);

        scatterjetpt.Set(flavor[iflv]+"_"+taggerslist[i],muinjettags[i],jetptbins,"ntag_pT_","muon-jet p_{T} vs ptrel");
        scatterjeteta.Set(flavor[iflv]+"_"+taggerslist[i],muinjettags[i],jetetabins,"ntag_eta_","muon-jet eta vs ptrel");

        scatterjeteta.Set(flavor[iflv]+"_"+taggerslist[i],muinjettags[i],jetetabins,"ptag_eta_","muon-jet eta vs ptrel");
        scatterjetpt.Set(flavor[iflv]+"_"+taggerslist[i],muinjettags[i],jetptbins,"ptag_pT_","muon-jet p_{T} vs ptrel");
     }
  }

  // Set to save jet histogram errors
  PtrelHistos.Sumw2(); 
  DeltaRHistos.Sumw2();
  Chi2muHistos.Sumw2();  
  nmuonhistos.Sumw2();  
  nmuonhits.Sumw2();  
  scatterjeteta.Sumw2();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagPATClosureAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE( BTagPATClosureAnalyzer );
