// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeCandidatesProducer
// 
/**\class CascadeCandidatesProducer CascadeCandidatesProducer.h Analyzers/CascadeProducer/src/CascadeCandidatesProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeCandidatesProducer.h,v 1.1 2011/10/12 19:06:40 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_PRODUCER_CANDIDATES_H
#define RECOVERTEX__CASCADE_PRODUCER_CANDIDATES_H

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//
// class declaration
//

class CascadeCandidatesProducer : public edm::EDProducer {
   public:
      explicit CascadeCandidatesProducer(const edm::ParameterSet&);
      ~CascadeCandidatesProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
  //Input parameters
  edm::ParameterSet ConfigParameters; //Config Parameters passed to constructor and accessor funtions
  std::string CasAlgo;           //Label of Cascade collection
  std::string CasDecayName;      //Name of Cascade Collection (Cascade/Omega)
  edm::InputTag tracksAlgo;      //Label of TrackCollection
  double Xiip3dcut;              //Impact parameter cut at primary vertex.
  double XiCLcut;                //Xi CL at decay for non-pixel tracks candidates.
  double Xidelsig;               //Xi Significance Separation of decay from primary.
  double withpixels;             //Include in selection Xi/Omega which have real pixel track.
  double Xiip3ddecaycut;         //Impact parameter cut at decay vertex.
};

#endif
