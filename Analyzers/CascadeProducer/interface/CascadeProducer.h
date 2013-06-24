// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeProducer
// 
/**\class CascadeProducer CascadeProducer.h Analyzers/CascadeProducer/src/CascadeProducer.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeProducer.h,v 1.2 2009/07/07 00:49:32 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_PRODUCER_H
#define RECOVERTEX__CASCADE_PRODUCER_H

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

class CascadeProducer : public edm::EDProducer {
   public:
      explicit CascadeProducer(const edm::ParameterSet&);
      ~CascadeProducer();

   private:
      virtual void beginJob() ;
      virtual void produce(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      // ----------member data ---------------------------
      edm::ParameterSet ConfigParameters; //Config Parameters passed to fitter
      std::string   LambdaAlgoLabel;  //Vees input Collection
      edm::InputTag TracksAlgoLabel;  //Tracks input Collection
};

#endif
