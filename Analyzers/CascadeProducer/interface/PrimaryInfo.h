#ifndef RECOVERTEX_PRIMARYINFO_H
#define RECOVERTEX_PRIMARYINFO_H
// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      PrimaryInfo
//
/**\class PrimaryInfo PrimaryInfo.h Analyzers/CascadeProducer/interface/PrimaryInfo.h

 Description: Tool for Primary Vertices Information and Handling

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: PrimaryInfo.h,v 1.4 2012/07/17 21:01:53 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"

//
// class declaration
//

class PrimaryInfo {
   public:
      explicit PrimaryInfo(const edm::Event& iEvent
		          ,const std::string& PrimariesColl);
      explicit PrimaryInfo(const edm::Event& iEvent
		          ,const edm::ParameterSet& iConfig);
      ~PrimaryInfo();
      double cl(){ return PrimaryProb;}
      double chi2(){return PrimaryChi2;}
      double ndof(){return PrimaryNdof;}
      int size(){ return primariesHandle->size(); }
      reco::VertexCollection::const_iterator BestVertex(){return bestVertex;}
      reco::VertexCollection::const_iterator HigherMultiplicity();
      reco::VertexCollection::const_iterator HigherCosAlpha(reco::Vertex &SecondaryVertex, GlobalVector &p);
      reco::VertexCollection::const_iterator HigherCosAlpha(GlobalPoint &SecondaryVertex, GlobalVector &p);
   private:
      // ----------member data ---------------------------		   		   
      bool ThereArePrimaries;      //flag to check validity of primaries info
      std::string PrimColl;        //Label of Primaries Collection
      edm::Handle<reco::VertexCollection> primariesHandle;     //Handles for Primary Vertex
      void SetInfo();              //Based on Highest Multiplicity Primary
      double PrimaryChi2;
      double PrimaryNdof;
      double PrimaryProb;
      double CosAlpha(reco::VertexCollection::const_iterator prm, GlobalPoint &SecondaryVertex, GlobalVector &p);
      reco::VertexCollection::const_iterator bestVertex;
};
#endif
