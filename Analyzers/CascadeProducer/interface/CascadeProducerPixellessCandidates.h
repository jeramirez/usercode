// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeProducerPixellessCandidates
// 
/**\class CascadeProducerPixellessCandidates CascadeProducerPixellessCandidates.h Analyzers/CascadeProducer/src/CascadeProducerPixellessCandidates.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeProducerPixellessCandidates.h,v 1.1 2011/10/11 19:18:19 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_PRODUCER_PIXELLESS_CANDIDATES_H
#define RECOVERTEX__CASCADE_PRODUCER_PIXELLESS_CANDIDATES_H

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

class CascadeProducerPixellessCandidates : public edm::EDProducer {
   public:
      explicit CascadeProducerPixellessCandidates(const edm::ParameterSet&);
      ~CascadeProducerPixellessCandidates();

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

};

#endif
