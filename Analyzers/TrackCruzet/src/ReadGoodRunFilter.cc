// -*- C++ -*-
//
// Package:    ReadGoodRunFilter
// Class:      ReadGoodRunFilter
// 
/**\class ReadMuonFilter ReadMuonFilter.cc Analyzers/ReadMuonFilter/src/ReadMuonFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Miguel Bonnett"
//         Created:  Thu Oct 16 14:12:58 CDT 2008
// $Id: ReadGoodRunFilter.cc,v 1.1 2008/12/03 17:27:28 bonnett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

 /* Collaborating Class Header */
 #include "FWCore/Framework/interface/ESHandle.h"
 #include "FWCore/MessageLogger/interface/MessageLogger.h"
 
 
 /* C++ Headers */
#include <iostream>
#include<fstream>
#include <vector>
#include <map>

using namespace std;
using namespace edm;
//
// class declaration
//

class ReadGoodRunFilter : public edm::EDFilter {
   public:
      explicit ReadGoodRunFilter(const edm::ParameterSet&);
      ~ReadGoodRunFilter();

    /* Operations */
   private:
      virtual void beginJob(const edm::EventSetup&) ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      std::ofstream filemuon;
      std::ofstream filebadrun;
//      std::ofstream fasciiFile;
      std::string fasciiGoodRun;

     int nl;
     std::map<int,int> run_n;
//     std::map<int,int> ev_n;

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
ReadGoodRunFilter::ReadGoodRunFilter(const edm::ParameterSet& iConfig)
{
    filemuon.open("Goodrun.txt");
    filebadrun.open("Badrun.txt");
    fasciiGoodRun = iConfig.getParameter<string>("listGoodRun");

}


ReadGoodRunFilter::~ReadGoodRunFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
ReadGoodRunFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    int run_nr = iEvent.id().run();
//    int ev_nr  = iEvent.id().event();
    
    std::string ntracker;
    bool accept = false;
    
    for(int i=0; i<nl ; i++){ 
//        if( run_nr == run_n[i] && ev_nr == ev_n[i]) {
        if( run_nr == run_n[i]) {
            std::cout<<" run  "<< run_nr << std::endl;
            filemuon<<" run  "<< run_nr << std::endl;
                accept = true;
                return accept;
        }
//        if( run_nr != run_n[i]) {filebadrun<<" run  "<< run_nr << std::endl;}
//       else {filebadrun<<" run  "<< run_nr << std::endl;}
   }    
    return accept;


}
 
// ------------ method called once each job just before starting event loop  ------------
void 
ReadGoodRunFilter::beginJob(const edm::EventSetup&)
{

     FILE * file;
     file = fopen(fasciiGoodRun.c_str(),"r");
     ifstream readfile;
     readfile.open(fasciiGoodRun.c_str());
     if(!readfile)cout<<"Can't read the file";
     
     int c;
     int run;
//     int event;
     nl =0;

     while((c=fgetc(file))!=EOF){
         if(c == '\n') nl++;
     }
     cout<<"nl = "<<nl<<endl;

     run_n[nl];
 //    ev_n[nl];
     for(int j=0 ; j<nl; j++){
         readfile >> run; 
     //    readfile >> event;
         run_n[j]= run; 
      //   ev_n[j]= event;     
     } 
     readfile.close();
    
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ReadGoodRunFilter::endJob() {
}
    

//define this as a plug-in
DEFINE_FWK_MODULE(ReadGoodRunFilter);
