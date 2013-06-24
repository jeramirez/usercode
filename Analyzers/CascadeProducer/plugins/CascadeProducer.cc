// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeProducer
// 
/**\class CascadeProducer CascadeProducer.cc Analyzers/CascadeProducer/src/CascadeProducer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeProducer.cc,v 1.1 2011/06/29 22:03:56 jramirez Exp $
//
//


// system include files
#include <memory>
#include "Analyzers/CascadeProducer/interface/CascadeProducer.h"
#include "Analyzers/CascadeProducer/interface/CascadeFitter.h"
//
// constants, enums and typedefs
//


//
// static data member definitions
//

//
// constructors and destructor
//
CascadeProducer::CascadeProducer(const edm::ParameterSet& iConfig):
   ConfigParameters(iConfig)
//   ,LambdaAlgoLabel(iConfig.getParameter<std::string>("v0Algo"))
//   ,TracksAlgoLabel(iConfig.getParameter<edm::InputTag>("trackingAlgo"))
{
   //register your products
/* Examples
   produces<ExampleData2>();

   //if do put with a label
   produces<ExampleData2>("label");
*/
   produces< reco::VertexCompositeCandidateCollection >("Cascade");
   produces< reco::VertexCompositeCandidateCollection >("Omega");

   //now do what ever other initialization is needed
  
}


CascadeProducer::~CascadeProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
CascadeProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
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
   // Create CascadeFitter object which reconstructs the vertices and creates
   //  (and contains) collections of Cascades, Omegas
   CascadeFitter theCascades(ConfigParameters, iEvent, iSetup);
   

   // Create auto_ptr for each collection to be stored in the Event
   std::auto_ptr< reco::VertexCompositeCandidateCollection > 
     XiCandidates( new reco::VertexCompositeCandidateCollection );
   XiCandidates->reserve( theCascades.getCascades().size() ); 

   std::auto_ptr< reco::VertexCompositeCandidateCollection >
     OmegaCandidates( new reco::VertexCompositeCandidateCollection );
   OmegaCandidates->reserve( theCascades.getOmegas().size() );
 
   std::copy( theCascades.getCascades().begin(),
	      theCascades.getCascades().end(),
	      std::back_inserter(*XiCandidates) );
   std::copy( theCascades.getOmegas().begin(),
	      theCascades.getOmegas().end(),
	      std::back_inserter(*OmegaCandidates) );

   // Write the collections to the Event
   iEvent.put( XiCandidates, std::string("Cascade") );
   iEvent.put( OmegaCandidates, std::string("Omega") );
}

// ------------ method called once each job just before starting event loop  ------------
void 
CascadeProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
CascadeProducer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(CascadeProducer);
