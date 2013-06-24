#ifndef RECOVERTEX_VERTEXREFIT_H
#define RECOVERTEX_VERTEXREFIT_H
// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      VertexRefit
//
/**\class VertexRefit VertexRefit.h Analyzers/CascadeProducer/interface/VertexRefit.h

 Description: Tool for Refit Vertices using exclusion list of tracks

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: VertexRefit.h,v 1.2 2011/06/16 19:11:19 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
//Magnetic Field
#include "MagneticField/Engine/interface/MagneticField.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

//
// class declaration
//

class VertexRefit {
   public:
      explicit VertexRefit(reco::VertexCollection::const_iterator InputVertex
		          ,std::vector<reco::TrackRef>& ExclusionList
                          ,edm::ESHandle<MagneticField>& bFieldHandle
                          ,reco::BeamSpot &BeamSpot);
      ~VertexRefit();
      bool isValid(){return IsValid; };
      reco::Vertex Refitted(){return VertexRefitted; };
      void Dump(){std::cout << "Excluding " << numexcl << " track(s)" << std::endl;};
   private:
      // ----------member data ---------------------------		   		   
      bool IsValid;                //flag to check validity of vertex
      reco::Vertex VertexRefitted; //Vertex Refitted
      int unsigned numexcl;        //Number of excluded tracks from vertex
		   
};
#endif
