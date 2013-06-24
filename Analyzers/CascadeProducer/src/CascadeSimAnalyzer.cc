// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeSimAnalyzer
// 
/**\class CascadeSimAnalyzer CascadeSimAnalyzer.cc Analyzers/CascadeProducer/src/CascadeSimAnalyzer.cc

 Description: analyzer for simulated info on cascade output 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadeSimAnalyzer.cc,v 1.1 2009/11/13 17:21:12 jramirez Exp $
//
//


// system include files
#include <memory>

//references
#include "FWCore/Framework/interface/ESHandle.h"
//MC
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "Analyzers/CascadeProducer/interface/CascadeSimAnalyzer.h"
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CascadeSimAnalyzer::CascadeSimAnalyzer(
//CascadeSimAnalyzer::CascadeSimAnalyzer(const edm::ParameterSet& iConfig,
                                       const edm::Event& iEvent, 
				       const edm::EventSetup& iSetup
):Lambda(3122),Xi(3312),Omega(3334),SigmaPlus(3222),SigmaMinus(3112)
,numLambda(0),numALambda(0)
,numXi(0),numAXi(0)
,numOmega(0),numAOmega(0)
,numSigmaP(0),numASigmaP(0)
,numSigmaM(0),numASigmaM(0)
{
   //now do what ever initialization is needed
   std::string SimTkLabel("g4SimHits");
   edm::Handle<edm::SimTrackContainer> SimTk;
   iEvent.getByLabel(SimTkLabel, SimTk);
   if (SimTk.isValid()){
      theSimTracks.insert( theSimTracks.end(), SimTk->begin(), SimTk->end() );
   }
   std::string SimVtxLabel("g4SimHits");
   edm::Handle<edm::SimVertexContainer> SimVtx;
   iEvent.getByLabel(SimVtxLabel, SimVtx);
   if (SimVtx.isValid()){
     theSimVerts.insert( theSimVerts.end(), SimVtx->begin(), SimVtx->end() );
   }
   // Count how many simulated /\0s and Cascade- and Omega- we have
   for(unsigned int ndx1 = 0; ndx1 < theSimTracks.size(); ndx1++) {
      if(theSimTracks[ndx1].type() ==  Lambda)
         numLambda  += 1;
      if(theSimTracks[ndx1].type() == -Lambda)
         numALambda  += 1;
      if(theSimTracks[ndx1].type() ==  Xi)
         numXi  += 1;
      if(theSimTracks[ndx1].type() == -Xi)
         numAXi  += 1;
      if(theSimTracks[ndx1].type() ==  Omega)
         numOmega  += 1;
      if(theSimTracks[ndx1].type() == -Omega)
         numAOmega  += 1;
      if(theSimTracks[ndx1].type() ==  SigmaPlus)
         numSigmaP  += 1;
      if(theSimTracks[ndx1].type() == -SigmaPlus)
         numASigmaP  += 1;
      if(theSimTracks[ndx1].type() ==  SigmaMinus)
         numSigmaM  += 1;
      if(theSimTracks[ndx1].type() == -SigmaMinus)
         numASigmaM  += 1;

   }
   
}


CascadeSimAnalyzer::~CascadeSimAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

bool CascadeSimAnalyzer::SimTrackerIsValid() const{
   return theSimVerts.size()>0?true:false;
} 

bool CascadeSimAnalyzer::SimVertexIsValid() const{
   return theSimTracks.size()>0?true:false;
} 
//
// member functions
//

void
CascadeSimAnalyzer::print() const
{ 

  std::cout << " num Sim Vtx=" <<theSimVerts.size() 
	    << " num Sim Trk=" <<theSimTracks.size()
	    << " num Sim Lambda0=" <<numLambda
	    << " num Sim Xi-=" <<numXi
	    << " num Sim Om-=" <<numOmega
	    << " num Sim ALambda0=" <<numALambda
	    << " num Sim Xi+=" <<numAXi
	    << std::endl;

}//end print


