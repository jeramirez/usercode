// -*- C++ -*-
//
// Package:    PatAlgos
// Class:      BTagPATAnalyzer
// 
/**\class BTagPATAnalyzer BTagPATAnalyzer.cc PhysicsTools/PatAlgos/test/BTagPATAnalyzer.cc

 Description: <A very (very) simple CMSSW analyzer for PAT objects>

 Implementation:
 
 this analyzer shows how to loop over PAT output. 
*/
//
// Original Author:  Freya Blekman
//         Created:  Mon Apr 21 10:03:50 CEST 2008
// $Id: BTagPATAnalyzer.cc,v 1.6 2009/06/15 22:27:28 jramirez Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "PhysicsTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/InputTag.h"


#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "TH1D.h"
#include "TH2D.h"
#include <map>

#include "DataFormats/Common/interface/View.h"
#include <string>
//ed++
#include "Analyzers/BTagPAT/interface/S8bPerformance.h"
#include "PhysicsTools/PatUtils/interface/bJetSelector.h"
#include "Analyzers/BTagPAT/interface/BTagPATCommonHistos.h"

#include "TGraphErrors.h"
//ed--
//
// class declaration
//

class BTagPATAnalyzer : public edm::EDAnalyzer {
   public:
      explicit BTagPATAnalyzer(const edm::ParameterSet&);
      ~BTagPATAnalyzer();


   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      // ----------member data ---------------------------

  
  std::map<std::string,TH1D*> histocontainer_;   // simple map to contain all histograms. Histograms are booked in the beginJob() method
  
  edm::InputTag jetLabel_;
//ed++
  edm::ParameterSet PatBjet_;
  std::string  BTagmethod_;
  std::string  BTagpurity_;
  std::string  BTagdiscriminator_;
  double  BTagdisccut_;
  double  BTagdiscmax_;
  bool    BTagverbose;
  S8bPerformance BTagPerf[18];
  std::string    discname[18];   
  std::string    bname[18];   
  std::string    cname[18];   
  std::map<int,std::string> udsgname;   
  size_t nbjetstotal      ;
  size_t ncjetstotal      ;
  size_t nlightjetstotal  ;
  size_t ngluonjetstotal  ;
  size_t nbtaggedjetstotal;
  size_t njets30          ;
  size_t nbjets30         ;
  std::map<std::string,TH2D*> h2_; // simple map to contain 2D  histograms. Histograms are booked in the beginJob() method
  std::map<std::string,TGraph*> graphcontainer_; // simple map to contain all graphs. Graphs are booked in the beginJob() method
  std::map<std::string,TGraphErrors*> grapherrorscontainer_; // simple map to contain all graphs. Graphs are booked in the beginJob() method
  bJetSelector BTagger;
  BTagPATCommonHistos BTagHistograms;
//ed--
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BTagPATAnalyzer::BTagPATAnalyzer(const edm::ParameterSet& iConfig):
  histocontainer_(),
  jetLabel_(iConfig.getUntrackedParameter<edm::InputTag>("jetTag"))
//ed++
  ,graphcontainer_()
  ,BTagger(iConfig.getParameter< edm::ParameterSet >("BJetOperatingPoints"))
  ,PatBjet_(iConfig.getParameter< edm::ParameterSet >("BjetTag"))
  ,BTagmethod_(PatBjet_.getUntrackedParameter<std::string>("tagger","TC2")),
  BTagHistograms(iConfig),
  BTagpurity_(PatBjet_.getParameter<std::string>("purity")),
  BTagdiscriminator_(PatBjet_.getParameter<std::string>("discriminator")),
  BTagdisccut_(PatBjet_.getUntrackedParameter<double>("mindiscriminatorcut",5.0)),
  BTagdiscmax_(PatBjet_.getUntrackedParameter<double>("maxdiscriminatorcut",15.0)),
  BTagverbose(PatBjet_.getUntrackedParameter<bool>("verbose",false))
//ed--
{
   //now do what ever initialization is needed
}


BTagPATAnalyzer::~BTagPATAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
void
BTagPATAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;


   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   // first: get all objects from the event.
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

   edm::Handle<edm::View<pat::Jet> > jetHandle;
   iEvent.getByLabel(jetLabel_,jetHandle);
   edm::View<pat::Jet> jets = *jetHandle;


   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   // example of a loop over objects... this works identical for all vectors defined above
   //   once you have a jet object you can use all methods defined in the header file 
   // (DataFormats/PatCandidates/interface/Jet.h and equivalent for other objects)
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
//ed++
   size_t nbjetscounter      = 0;
   size_t ncjetscounter      = 0;
   size_t nlightjetscounter  = 0;
   size_t ngluonjetscounter  = 0;
   size_t nbtaggedjetscounter= 0;
   for(edm::View<pat::Jet>::const_iterator jet_iter = jets.begin(); jet_iter!=jets.end(); ++jet_iter){
//       BTagTool.Add(jet_iter->bDiscriminator(BTagdiscriminator_),abs(jet_iter->partonFlavour()));

       float isb    = jet_iter->bDiscriminator(BTagdiscriminator_);
       int flavor   = jet_iter->partonFlavour();
       //
       // Fill in for performance standard pt(uncorrected) >30 and abs(eta)<2.4 
       if (jet_iter->correctedJet("raw").pt() > 30 &&
	   //       if( jet_iter->noCorrJet().pt() > 30  &&
	   fabs(jet_iter->eta()) < 2.4
	&& fabs(jet_iter->eta()) > 0
	 ){
          BTagPerf[0].Add(isb,abs(flavor));
          njets30++;
          if (flavor ==  5 || flavor ==  -5 ) nbjets30++;

       }
       //
       // Fill in for performance pt(uncorrected) >10 and abs(eta)<2.4
       // 4 cases 
       //   1)just those cuts
       //   2)# tracks > 0
       //   3)# tracks > 1
       //   4)# tracks > 2
       if ( jet_iter->correctedJet("raw").pt() > 10 &&
       //       if( jet_iter->noCorrJet().pt() > 10  &&
           fabs(jet_iter->eta()) < 2.4
        && fabs(jet_iter->eta()) > 0
         ){
	  BTagPerf[1].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 0) BTagPerf[2].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 1) BTagPerf[3].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 2) BTagPerf[4].Add(isb,abs(flavor));

	  if(jet_iter->associatedTracks().size() > 0 && jet_iter->nConstituents() > 1) BTagPerf[10].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 1 && jet_iter->nConstituents() > 1) BTagPerf[11].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 2 && jet_iter->nConstituents() > 1) BTagPerf[12].Add(isb,abs(flavor));
       }
       //
       // Fill in for performance pt(corrected) >30 and abs(eta)<2.4
       // 4 cases 
       //   1)just those cuts
       //   2)# tracks > 0
       //   3)# tracks > 1
       //   4)# tracks > 2
 
      if( jet_iter->pt() > 30  &&
           fabs(jet_iter->eta()) < 2.4
        && fabs(jet_iter->eta()) > 0
         ){

	  BTagPerf[5].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 0) BTagPerf[6].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 1) BTagPerf[7].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 2) BTagPerf[8].Add(isb,abs(flavor));

	  if(jet_iter->associatedTracks().size() > 0 && jet_iter->nConstituents() > 1) BTagPerf[13].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 1 && jet_iter->nConstituents() > 1) BTagPerf[14].Add(isb,abs(flavor));
	  if(jet_iter->associatedTracks().size() > 2 && jet_iter->nConstituents() > 1) BTagPerf[15].Add(isb,abs(flavor));

       }
       //
       // check auxiliar function compared when filtering data
       //
       if(
	   fabs(jet_iter->eta()) < 2.4
	&& fabs(jet_iter->eta()) > 0
	&& BTagger.IsbTag(*jet_iter,BTagpurity_,BTagdiscriminator_)
	 ){
	      BTagPerf[9].Add(isb,abs(flavor));
          nbtaggedjetscounter++;
          if (-4 < flavor && flavor < 4 && flavor != 0 ) nlightjetscounter++;          
          if (flavor == 21 || flavor == -21 ) ngluonjetscounter++;
          if (flavor ==  4 || flavor ==  -4 ) ncjetscounter++;          
          if (flavor ==  5 || flavor ==  -5 ) nbjetscounter++;          
	}//end if for tagged cut

	//Fill histograms
	BTagHistograms.Fill(jet_iter,"all");

//	BTagHistograms.Fill(jet_iter->nConstituents(),"nConstituents");

	if (flavor ==  0  ) BTagHistograms.Fill(jet_iter,"no_flavor");
	if (flavor ==  5 || flavor ==  -5 ) BTagHistograms.Fill(jet_iter,"b");
	if (flavor ==  4 || flavor ==  -4 ) BTagHistograms.Fill(jet_iter,"c");
        if ((-4 < flavor && flavor < 4 && flavor != 0 )||(flavor == 21 || flavor == -21 ))
		BTagHistograms.Fill(jet_iter,"udsg");


   }//end loop over jets

   nbtaggedjetstotal += nbtaggedjetscounter;
   ngluonjetstotal   += ngluonjetscounter;
   nlightjetstotal   += nlightjetscounter;
   ncjetstotal       += ncjetscounter;
   nbjetstotal       += nbjetscounter;
  if (BTagverbose){
      std::cout << "BTag discriminator("<< BTagdiscriminator_<<")>"
                << BTagdisccut_
                << std::endl;
      std::cout << "   # BJets Tagged(b)      =" << nbjetscounter << std::endl; 
      std::cout << "   # BJets Tagged(c)      =" << ncjetscounter << std::endl; 
      std::cout << "   # BJets Tagged(ubs)    =" << nlightjetscounter << std::endl; 
      std::cout << "   # BJets Tagged(g)      =" << ngluonjetscounter << std::endl; 
      std::cout << "   # Event BJets Tagged   =" << nbtaggedjetscounter << std::endl; 
      std::cout << "   # Total BJets Tagged   =" << nbtaggedjetstotal << std::endl; 
      if (nbtaggedjetscounter > 0) {
         double bjeteff = double(nbjetscounter)/double(nbtaggedjetscounter);     
         double mistagc = double(ncjetscounter)/double(nbtaggedjetscounter);
         double mistagl = double(nlightjetscounter)/double(nbtaggedjetscounter);
         double mistagg = double(ngluonjetscounter)/double(nbtaggedjetscounter);
         std::cout << "   # BJetEff    =" << bjeteff << std::endl;
         std::cout << "   # MissTagC   =" << mistagc << std::endl;
         std::cout << "   # MissTagUDS =" << mistagl << std::endl;
         std::cout << "   # MissTagG   =" << mistagg << std::endl;
         std::cout << "BTag method "<< BTagmethod_<< std::endl;
         std::cout << "BTag purity "<< BTagpurity_<< std::endl;
      }    
   }
//ed--

   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   //histocontainer_ is of type std::map<edm::View, TH1D*>. This means you can use it with this syntax:
   // histocontainer_["histname"]->Fill(x); 
   // histocontainer_["histname"]->Draw(); 
   // etc, etc. Essentially you use the histname string to look up a pointer to a TH1D* 
   // which you can do everything to you would normally do in ROOT.
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
   // for the other objects just quickly book the multiplicity. Again, just use the same infrastructure as for jets if you want to loop over them.
   //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
}
// ------------ method called once each job just before starting event loop  ------------
void 
BTagPATAnalyzer::beginJob(const edm::EventSetup&)
{
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // define some histograms using the framework tfileservice. Define the output file name in your .cfg.
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  edm::Service<TFileService> fs;
  
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  //histocontainer_ is of type std::map<std::string, TH1D*>. This means you can use it with this syntax:
  // histocontainer_["histname"]->Fill(x); 
  // histocontainer_["histname"]->Draw(); 
  // etc, etc. Essentially you use the histname string to look up a pointer to a TH1D* 
  // which you can do everything to you would normally do in ROOT.
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}

  
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}
  // here we book new histograms:
  //{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}


//ed++
  TString suffix1="_test";

  njets30          = 0;
  nbjets30         = 0;
  nbjetstotal      = 0;
  ncjetstotal      = 0;
  nlightjetstotal  = 0;
  ngluonjetstotal  = 0;
  nbtaggedjetstotal= 0;
//set performance variables collector
  for (int i=0; i < 16; i++){
    BTagPerf[i].Set(BTagmethod_);
    BTagPerf[i].SetMinDiscriminator(BTagdisccut_);
    BTagPerf[i].SetMaxDiscriminator(BTagdiscmax_);
  }
//  BTagTool.SetNcuts(200); //just 1 bin
  histocontainer_["njets"]=fs->make<TH1D>("njets","jet multiplicity for jets with p_{T} > 50 GeV/c",10,0,10);
// Std. 30 pt uncorr cut for performance
  discname[0]="disc"+BTagmethod_+"_udsg";   
  bname[0]   ="g"+BTagmethod_+"_b";   
  cname[0]   ="g"+BTagmethod_+"_c";   
  udsgname[0]="g"+BTagmethod_+"_udsg";   

// 10 pt uncorr for performance + all,>0,>1,>2 tracks
  discname[1]="Uncor10_disc"+BTagmethod_+"_udsg";
  bname[1]   ="Uncor10_g"+BTagmethod_+"_b";   
  cname[1]   ="Uncor10_g"+BTagmethod_+"_c";   
  udsgname[1]="Uncor10_g"+BTagmethod_+"_udsg";
  discname[2]="Uncor10t0_disc"+BTagmethod_+"_udsg";
  bname[2]   ="Uncor10t0_g"+BTagmethod_+"_b";   
  cname[2]   ="Uncor10t0_g"+BTagmethod_+"_c";   
  udsgname[2]="Uncor10t0_g"+BTagmethod_+"_udsg";
  discname[3]="Uncor10t1_disc"+BTagmethod_+"_udsg";
  bname[3]   ="Uncor10t1_g"+BTagmethod_+"_b";   
  cname[3]   ="Uncor10t1_g"+BTagmethod_+"_c";   
  udsgname[3]="Uncor10t1_g"+BTagmethod_+"_udsg";
  discname[4]="Uncor10t2_disc"+BTagmethod_+"_udsg";
  bname[4]   ="Uncor10t2_g"+BTagmethod_+"_b";   
  cname[4]   ="Uncor10t2_g"+BTagmethod_+"_c";   
  udsgname[4]="Uncor10t2_g"+BTagmethod_+"_udsg";
// 10 pt uncorr for performance + all,>0,>1,>2 tracks & nConstituent > 1
  discname[10]="Uncor10t0_nConstituent_disc"+BTagmethod_+"_udsg";
  bname[10]   ="Uncor10t0_nConstituent_g"+BTagmethod_+"_b";   
  cname[10]   ="Uncor10t0_nConstituent_g"+BTagmethod_+"_c";   
  udsgname[10]="Uncor10t0_nConstituent_g"+BTagmethod_+"_udsg";
  discname[11]="Uncor10t1_nConstituent_disc"+BTagmethod_+"_udsg";
  bname[11]   ="Uncor10t1_nConstituent_g"+BTagmethod_+"_b";   
  cname[11]   ="Uncor10t1_nConstituent_g"+BTagmethod_+"_c";   
  udsgname[11]="Uncor10t1_nConstituent_g"+BTagmethod_+"_udsg";
  discname[12]="Uncor10t2_nConstituent_disc"+BTagmethod_+"_udsg";
  bname[12]   ="Uncor10t2_nConstituent_g"+BTagmethod_+"_b";   
  cname[12]   ="Uncor10t2_nConstituent_g"+BTagmethod_+"_c";   
  udsgname[12]="Uncor10t2_nConstituent_g"+BTagmethod_+"_udsg";


// 30 pt corr for performance + all,>0,>1,>2 tracks   
  discname[5]="Corr30_disc"+BTagmethod_+"_udsg";
  bname[5]   ="Corr30_g"+BTagmethod_+"_b";   
  cname[5]   ="Corr30_g"+BTagmethod_+"_c";   
  udsgname[5]="Corr30_g"+BTagmethod_+"_udsg";
  discname[6]="Corr30t0_disc"+BTagmethod_+"_udsg";
  bname[6]   ="Corr30t0_g"+BTagmethod_+"_b";   
  cname[6]   ="Corr30t0_g"+BTagmethod_+"_c";   
  udsgname[6]="Corr30t0_g"+BTagmethod_+"_udsg";
  discname[7]="Corr30t1_disc"+BTagmethod_+"_udsg";
  bname[7]   ="Corr30t1_g"+BTagmethod_+"_b";   
  cname[7]   ="Corr30t1_g"+BTagmethod_+"_c";   
  udsgname[7]="Corr30t1_g"+BTagmethod_+"_udsg";
  discname[8]="Corr30t2_disc"+BTagmethod_+"_udsg";
  bname[8]   ="Corr30t2_g"+BTagmethod_+"_b";   
  cname[8]   ="Corr30t2_g"+BTagmethod_+"_c";   
  udsgname[8]="Corr30t2_g"+BTagmethod_+"_udsg";
// 30 pt corr for performance + all,>0,>1,>2 tracks  & nConstituent > 1  
  discname[13]="Corr30t0_nConstituent_disc"+BTagmethod_+"_udsg";
  bname[13]   ="Corr30t0_nConstituent_g"+BTagmethod_+"_b";   
  cname[13]   ="Corr30t0_nConstituent_g"+BTagmethod_+"_c";   
  udsgname[13]="Corr30t0_nConstituent_g"+BTagmethod_+"_udsg";
  discname[14]="Corr30t1_nConstituent_disc"+BTagmethod_+"_udsg";
  bname[14]   ="Corr30t1_nConstituent_g"+BTagmethod_+"_b";   
  cname[14]   ="Corr30t1_nConstituent_g"+BTagmethod_+"_c";   
  udsgname[14]="Corr30t1_nConstituent_g"+BTagmethod_+"_udsg";
  discname[15]="Corr30t2_nConstituent_disc"+BTagmethod_+"_udsg";
  bname[15]   ="Corr30t2_nConstituent_g"+BTagmethod_+"_b";   
  cname[15]   ="Corr30t2_nConstituent_g"+BTagmethod_+"_c";   
  udsgname[15]="Corr30t2_nConstituent_g"+BTagmethod_+"_udsg";

// check filter
  discname[9]="check_disc"+BTagmethod_+"_udsg";   
  bname[9]   ="check_g"+BTagmethod_+"_b";   
  cname[9]   ="check_g"+BTagmethod_+"_c";   
  udsgname[9]="check_g"+BTagmethod_+"_udsg";   


  for(int i=1; i<16;i++){
    graphcontainer_[discname[i]]      =fs->make<TGraph>(BTagPerf[i].GetN());       graphcontainer_[discname[i]]->SetName(discname[i].c_str());
    grapherrorscontainer_[bname[i]]   =fs->make<TGraphErrors>(BTagPerf[i].GetN()); grapherrorscontainer_[bname[i]]   ->SetName(bname[i].c_str());
    grapherrorscontainer_[cname[i]]   =fs->make<TGraphErrors>(BTagPerf[i].GetN()); grapherrorscontainer_[cname[i]]   ->SetName(cname[i].c_str());
    grapherrorscontainer_[udsgname[i]]=fs->make<TGraphErrors>(BTagPerf[i].GetN()); grapherrorscontainer_[udsgname[i]]->SetName(udsgname[i].c_str());   
  }
  //Define histograms
  BTagHistograms.Set("all");  
  BTagHistograms.Set("no_flavor");  
  BTagHistograms.Set("b");  
  BTagHistograms.Set("c");  
  BTagHistograms.Set("udsg");  

//  BTagHistograms.Set("nConstituents");  

  // Set to save histogram errors
  BTagHistograms.Sumw2();  

//ed--

}

// ------------ method called once each job just after ending the event loop  ------------
void 
BTagPATAnalyzer::endJob() {
//ed++
   edm::Service<TFileService> fs;

// Save performance plots as Tgraphs


   for (int i=1;i<16;i++){
      BTagPerf[i].Eval();
      for (int n=0; n<BTagPerf[i].GetN();  n++ ){
         graphcontainer_[discname[i]]       ->SetPoint(n,BTagPerf[i].GetArray("udsg")[n],BTagPerf[i].GetArray("discriminator")[n]);
         grapherrorscontainer_[bname[i]]    ->SetPoint(n,BTagPerf[i].GetArray("b")[n],BTagPerf[i].GetArray("b")[n]);
         grapherrorscontainer_[cname[i]]    ->SetPoint(n,BTagPerf[i].GetArray("b")[n],BTagPerf[i].GetArray("c")[n]);
         grapherrorscontainer_[udsgname[i]] ->SetPoint(n,BTagPerf[i].GetArray("b")[n],BTagPerf[i].GetArray("udsg")[n]);
         grapherrorscontainer_[bname[i]]    ->SetPointError(n,BTagPerf[i].GetArray("bErr")[n],BTagPerf[i].GetArray("bErr")[n]);
         grapherrorscontainer_[cname[i]]    ->SetPointError(n,BTagPerf[i].GetArray("bErr")[n],BTagPerf[i].GetArray("cErr")[n]);
         grapherrorscontainer_[udsgname[i]] ->SetPointError(n,BTagPerf[i].GetArray("bErr")[n],BTagPerf[i].GetArray("udsgErr")[n]);
      }//end for over BTagPerf[i] elements
      graphcontainer_[discname[i]]     ->SetTitle("Jet udsg-mistagging");
      grapherrorscontainer_[bname[i]]   ->SetTitle("Jet b-efficiency");
      grapherrorscontainer_[cname[i]]   ->SetTitle("Jet c-mistagging");
      grapherrorscontainer_[udsgname[i]]->SetTitle("discriminator vs udsg-mistagging");
   }//end for over [i]


// Save default cut performance plot
   BTagPerf[0].Eval();

//  TFileDirectory TaggerDir = fs->mkdir(BTagmethod_);
//  TGraphErrors *BTagger_b    = new TGraphErrors(BTagTool.GetN(),
  TGraphErrors *BTagger_b    = fs->mkdir(BTagmethod_).make<TGraphErrors>(BTagPerf[0].GetN(),
                                            BTagPerf[0].GetArray("b").GetArray(),BTagPerf[0].GetArray("b").GetArray(),
                                            BTagPerf[0].GetArray("bErr").GetArray(),BTagPerf[0].GetArray("bErr").GetArray());
        
  TGraphErrors *BTagger_c    = new TGraphErrors(BTagPerf[0].GetN(),
                                            BTagPerf[0].GetArray("b").GetArray(),BTagPerf[0].GetArray("c").GetArray(),
                                            BTagPerf[0].GetArray("bErr").GetArray(),BTagPerf[0].GetArray("cErr").GetArray());
                
  TGraphErrors *BTagger_udsg = new TGraphErrors(BTagPerf[0].GetN(),
                                               BTagPerf[0].GetArray("b").GetArray(),BTagPerf[0].GetArray("udsg").GetArray(),
                                               BTagPerf[0].GetArray("bErr").GetArray(),BTagPerf[0].GetArray("udsgErr").GetArray());
  TGraph *discBTagger_udsg   = new TGraph(BTagPerf[0].GetN(),
                                   BTagPerf[0].GetArray("udsg").GetArray(),
                                   BTagPerf[0].GetArray("discriminator").GetArray());
 
  BTagger_b->SetName(bname[0].c_str());
  BTagger_c->SetName(cname[0].c_str());
  BTagger_udsg->SetName(udsgname[0].c_str());
  discBTagger_udsg->SetName(discname[0].c_str());

  BTagger_b->SetTitle("Jet b-efficiency");
  BTagger_c->SetTitle("Jet c-mistagging");
  BTagger_udsg->SetTitle("Jet udsg-mistagging");
  discBTagger_udsg->SetTitle("discriminator vs udsg-mistagging");


  for (int i=1;i<16;i++){
   graphcontainer_[discname[i]]      ->Write();
   grapherrorscontainer_[bname[i]]   ->Write();
   grapherrorscontainer_[cname[i]]   ->Write();
   grapherrorscontainer_[udsgname[i]]->Write();
  }
  
  BTagger_b->Write();
  BTagger_c->Write();
  BTagger_udsg->Write();
  discBTagger_udsg->Write();
    

   if (BTagverbose){
      std::cout << "BTag discriminator("<< BTagdiscriminator_<<")>"
                << BTagdisccut_
                << " and Pt > 30 Gev/c"
                << std::endl;
      std::cout << "   #     b        Jets passed cut  TOTAL= " << nbjetstotal << std::endl; 
      std::cout << "   #     c        Jets passed cut  TOTAL= " << ncjetstotal << std::endl; 
      std::cout << "   #    uds       Jets passed cut  TOTAL= " << nlightjetstotal << std::endl; 
      std::cout << "   #     g        Jets passed cut  TOTAL= " << ngluonjetstotal << std::endl; 
      std::cout << "   # not matched  Jets passed cut  TOTAL= " 
                      << nbtaggedjetstotal-nbjetstotal-ncjetstotal-nlightjetstotal-ngluonjetstotal 
                      << std::endl; 
      std::cout << "   #    all       Jets passed cut  TOTAL= " << nbtaggedjetstotal << std::endl; 
      std::cout << "   #    all   Jets (Pt>30) no cut  TOTAL= " << njets30 << std::endl; 
      std::cout << "   #     b    Jets (Pt>30) no cut  TOTAL= " << nbjets30 << std::endl; 
      if (nbtaggedjetstotal > 0) {
         double bjeteff = double(nbjetstotal)/double(nbjets30);     
         double mistagc = double(ncjetstotal)/double(nbjets30);
         double mistagl = double(nlightjetstotal)/double(nbjets30);
         double mistagg = double(ngluonjetstotal)/double(nbjets30);
         std::cout << "   # BJetEff    TOTAL =" << bjeteff << std::endl;
         std::cout << "   # MissTagC   TOTAL =" << mistagc << std::endl;
         std::cout << "   # MissTagUDS TOTAL =" << mistagl << std::endl;
         std::cout << "   # MissTagG   TOTAL =" << mistagg << std::endl;         
      } 
      for (int i=0;i<BTagPerf[0].GetN();i++){
        std::cout << "  #   BJetEff Fran TOTAL["<<i+1<<"] =" << BTagPerf[0].GetMap("b")[i] << std::endl;
        std::cout << "  #   CJetEff Fran TOTAL["<<i+1<<"] =" << BTagPerf[0].GetMap("c")[i] << std::endl;
        std::cout << "  #udsgJetEff Fran TOTAL["<<i+1<<"] =" << BTagPerf[0].GetMap("udsg")[i] << std::endl;
      }
   }//endif BTagverbose
        std::cout << "  #  BJetall  Fran TOTAL =" << BTagPerf[0].Getb_all() << " check="<< BTagPerf[7].Getb_all()<< std::endl;
        std::cout << "  #  CJetall  Fran TOTAL =" << BTagPerf[0].Getc_all() << " check="<< BTagPerf[7].Getc_all()<<std::endl;
        std::cout << "  #udsgJetall Fran TOTAL =" << BTagPerf[0].Getudsg_all() << " check="<< BTagPerf[7].Getudsg_all()<< std::endl;
//ed--
}

//define this as a plug-in
DEFINE_FWK_MODULE(BTagPATAnalyzer);
