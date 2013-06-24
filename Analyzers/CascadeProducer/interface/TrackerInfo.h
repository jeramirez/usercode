#ifndef RECOVERTEX_TRACKERINFO_H
#define RECOVERTEX_TRACKERINFO_H
// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      TrackerInfo
//
/**\class TrackerInfo TrackerInfo.h Analyzers/CascadeProducer/interface/TrackerInfo.h

 Description: Tool for Getting Tracker Geometry Info

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: TrackerInfo.h,v 1.1 2011/06/29 22:06:57 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"


//
// class declaration
//

class TrackerInfo {
   public:
      explicit TrackerInfo(const edm::Event& iEvent
		          ,const edm::EventSetup& iSetup);
      ~TrackerInfo();
      bool IsBarrel(uint32_t subdetid);
      bool IsPixel(uint32_t subdetid);
      double rangeRZmax(uint32_t,uint32_t);
      double rangeRZmax(std::pair<uint32_t,uint32_t >);
      double rangeRmax(std::pair<uint32_t,uint32_t >);
      double rangeZmax(std::pair<uint32_t,uint32_t >);
      double rangeRmin(std::pair<uint32_t,uint32_t >);
      double rangeZmin(std::pair<uint32_t,uint32_t >);

   private:
      // ----------member data ---------------------------		   		   
      std::map<DetId,uint32_t> SubDetId;
      std::map<DetId,uint32_t> SubSubDetId;
      std::map< std::pair<uint32_t,uint32_t >,double> rangeRmax_; //given a pair of detector parameters get max R
      std::map< std::pair<uint32_t,uint32_t >,double> rangeRmin_; //given a pair of detector parameters get max R
      std::map< std::pair<uint32_t,uint32_t >,double> rangeZmax_; //given a pair of detector parameters get max R
      std::map< std::pair<uint32_t,uint32_t >,double> rangeZmin_; //given a pair of detector parameters get max R
		   
};
#endif
