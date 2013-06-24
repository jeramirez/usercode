// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeProducerPixellessCandidates
// 
/**\class CascadeProducerPixellessCandidates CascadeProducerPixellessCandidates.cc Analyzers/CascadeProducer/src/CascadeProducerPixellessCandidates.cc

 Description: select Pixelless  candidates.

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeProducerPixellessCandidates.cc,v 1.1 2011/10/11 19:18:19 jramirez Exp $
//
//


// system include files
#include <memory>
#include "Analyzers/CascadeProducer/interface/CascadeProducerPixellessCandidates.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h" //cascade collection
#include "DataFormats/TrackReco/interface/TrackFwd.h"                 //track collection
#include "Analyzers/CascadeProducer/interface/CascadePixelTrackFinder.h" //helper class
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
CascadeProducerPixellessCandidates::CascadeProducerPixellessCandidates(const edm::ParameterSet& iConfig):
   ConfigParameters(iConfig)
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/

   // Get the track reco algorithm from the ParameterSet
   tracksAlgo = iConfig.getParameter<edm::InputTag>("trackingAlgo");

   // Get the cas reco algorithm from the ParameterSet
   CasAlgo = iConfig.getParameter<std::string>("CasAlgo");

   // Cas Collection Label
   CasDecayName = iConfig.getUntrackedParameter<std::string>("CasDecayName","Cascade");


   produces< reco::VertexCompositeCandidateCollection >(CasDecayName);

}


CascadeProducerPixellessCandidates::~CascadeProducerPixellessCandidates()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CascadeProducerPixellessCandidates::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
/* This is an event example
   //Read 'ExampleData' from the Event
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);

   //Use the ExampleData to create an ExampleData2 which 
   // is put into the Event
   std::auto_ptr<ExampleData2> pOut(new ExampleData2(*pIn));
   iEvent.put(pOut);
*/

/* this is an EventSetup example
   //Read SetupData from the SetupRecord in the EventSetup
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
*/

  // Handles for tracks
  edm::Handle<reco::TrackCollection> theTrackHandle;
  // Get the tracks from the event
  bool found = iEvent.getByLabel(tracksAlgo, theTrackHandle);
  if( !found ) {
    LogError("CascadeProducerPixellessCandidates")
                              << "No tracks in event!, "
                              << "missing collection: " 
                              << tracksAlgo
                              << std::endl;
  }

  //Handles for Cascades
  edm::Handle<reco::VertexCompositeCandidateCollection> theCasHandle;
  found = iEvent.getByLabel(CasAlgo, CasDecayName, theCasHandle);
  if( !found ) {
    LogError("CascadeProducerPixellessCandidates")
                              << "No input Xi/Omega in event!, "
                              << "missing collection: " 
                              << CasAlgo << " " << CasDecayName 
                              << std::endl;
  }

  //Run CascadePixelTrackFinder
  CascadePixelTrackFinder MyXiTrack(tracksAlgo,CasAlgo,CasDecayName,iEvent,iSetup);

  //Load  Cascades into a vector of composite candidates
  std::vector<reco::VertexCompositeCandidate> theCascades;
  theCascades.insert( theCascades.end(),
                     theCasHandle->begin(),
                     theCasHandle->end()
                     );

  //Cascades to be saved into vector of composite candidates
  std::vector<reco::VertexCompositeCandidate> theCascadesSkimmed;
  std::vector<unsigned int> theCascadesIndexes;
  //Loop over cascades
  for(unsigned int casindex = 0;
      casindex < theCascades.size();
      casindex++) {
    if (MyXiTrack.NumPixelXiTracks(casindex)>0) {
      theCascadesSkimmed.push_back(theCascades[casindex]);
      theCascadesIndexes.push_back(casindex);
    }
  }

  //Output collection
  std:: auto_ptr< reco::VertexCompositeCandidateCollection > CasOutputColl(new reco::VertexCompositeCandidateCollection());
  CasOutputColl->resize(theCascadesSkimmed.size());

  //Save skimmed
  for (unsigned int iCopy=0;iCopy!=theCascadesSkimmed.size();++iCopy){
    unsigned int casindex = theCascadesIndexes[iCopy];
    const reco::VertexCompositeCandidate &iCas = (*theCasHandle)[casindex];
    (*CasOutputColl)[iCopy] = iCas; 
  }

  // Write the collections to the Event
  iEvent.put( CasOutputColl, CasDecayName );

  
   // Create CascadeFitter object which reconstructs the vertices and creates
   //  (and contains) collections of Cascades, Omegas
   // CascadeFitter theCascades(ConfigParameters, iEvent, iSetup);
   

   // Create auto_ptr for each collection to be stored in the Event
//   std::auto_ptr< reco::VertexCompositeCandidateCollection > 
//     XiCandidates( new reco::VertexCompositeCandidateCollection );
//   XiCandidates->reserve( theCascades.getCascades().size() ); 

//   std::auto_ptr< reco::VertexCompositeCandidateCollection >
//     OmegaCandidates( new reco::VertexCompositeCandidateCollection );
//   OmegaCandidates->reserve( theCascades.getOmegas().size() );
 
//   std::copy( theCascades.getCascades().begin(),
//	      theCascades.getCascades().end(),
//	      std::back_inserter(*XiCandidates) );
//   std::copy( theCascades.getOmegas().begin(),
//	      theCascades.getOmegas().end(),
//	      std::back_inserter(*OmegaCandidates) );

   // Write the collections to the Event
//   iEvent.put( XiCandidates, std::string("Cascade") );
//   iEvent.put( OmegaCandidates, std::string("Omega") );
}

// ------------ method called once each job just before starting event loop  ------------
void 
CascadeProducerPixellessCandidates::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CascadeProducerPixellessCandidates::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CascadeProducerPixellessCandidates);
