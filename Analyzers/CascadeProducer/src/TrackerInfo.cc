// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      TrackerInfo
// 
/**\class TrackerInfo TrackerInfo.cc Analyzers/CascadeProducer/src/TrackerInfo.cc

 Description: Tool for Getting Tracker Geometry Info 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: TrackerInfo.cc,v 1.1 2011/06/29 22:06:41 jramirez Exp $
//
//

// Tools to unpacking Tracker Geometry
#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

// To convert detId to subdet/layer number.
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/DetId/interface/DetId.h"

//exceptions
#include "Utilities/General/interface/CMSexception.h"

#include "Analyzers/CascadeProducer/interface/TrackerInfo.h"

//
// constructors
//
TrackerInfo::TrackerInfo(const edm::Event& iEvent,
			 const edm::EventSetup& iSetup
                        )
{


//Unpack Tracker Geometry
  edm::ESHandle<TrackerGeometry> TGHandle;
  iSetup.get<TrackerDigiGeometryRecord>().get(TGHandle);
  if (!TGHandle.isValid()) {
     edm::LogWarning("TrackerInfo")
       << "Unable to find TrackerDigiGeometry in event!";
     return;
  }

//loop over det units
  for (TrackerGeometry::DetContainer::const_iterator it = TGHandle->dets().begin();
       it != TGHandle->dets().end();
       ++it) {
       DetId detid = ((*it)->geographicalId());
       SubDetId[detid] = detid.subdetId();
       if (detid.subdetId() == StripSubdetector::TIB) {
         SubSubDetId[detid] = TIBDetId(detid).layer();
       } else if (detid.subdetId() == StripSubdetector::TOB) {
         SubSubDetId[detid] = TOBDetId(detid).layer();
       } else if (detid.subdetId() == StripSubdetector::TID) {
         SubSubDetId[detid] = TIDDetId(detid).wheel();
       } else if (detid.subdetId() == StripSubdetector::TEC) {
         SubSubDetId[detid] = TECDetId(detid).wheel();
       } else if (detid.subdetId() == PixelSubdetector::PixelBarrel) {
         SubSubDetId[detid] = PXBDetId(detid).layer();
       } else if (detid.subdetId() == PixelSubdetector::PixelEndcap) {
         SubSubDetId[detid] = PXFDetId(detid).disk();
       } else {
         //never should happen
         throw cms::Exception("TrackerInfo") << "Problem unpacking Tracker Geometry";
       }   

       std::pair<uint32_t,uint32_t > DetID(SubDetId[detid],SubSubDetId[detid]);
       
       double rho = (*it)->position().perp();
       double z   = std::abs((*it)->position().z());
       
      //is first time? set defaults
       if (rangeRmax_.find(DetID) == rangeRmax_.end()){
          rangeRmax_[DetID] = 0;
	  rangeRmin_[DetID] = 999;  //set too|
	  rangeZmax_[DetID] = 0.;   //set too >creating map default
	  rangeZmin_[DetID] = 999.; //set too|
       }
       
       //set new min and max and save for a given subdetector
       if (rangeRmax_[DetID] < rho) rangeRmax_[DetID]=rho;
       if (rangeRmin_[DetID] > rho) rangeRmin_[DetID]=rho;
       if (rangeZmax_[DetID] < z)   rangeZmax_[DetID]=z;
       if (rangeZmin_[DetID] > z)   rangeZmin_[DetID]=z;
       
  }

}

//
// destructors
//
TrackerInfo::~TrackerInfo()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
//
// methods
//
bool TrackerInfo::IsBarrel(uint32_t subdetid){ 
  return (subdetid == StripSubdetector::TIB||
          subdetid == StripSubdetector::TOB||
	  subdetid == PixelSubdetector::PixelBarrel);
}
bool TrackerInfo::IsPixel(uint32_t subdetid){ 
  return (subdetid == PixelSubdetector::PixelEndcap||
	  subdetid == PixelSubdetector::PixelBarrel);
}
double TrackerInfo::rangeRZmax(std::pair<uint32_t,uint32_t > DetID)
{
  return (IsBarrel(DetID.first)?
               rangeRmax(DetID):
	       rangeZmax(DetID));
  
}
//overload rangeRZmax
double TrackerInfo::rangeRZmax(uint32_t subdetid,uint32_t layer)
{ 
  std::pair<uint32_t,uint32_t > DetID(subdetid,layer);
  return rangeRZmax(DetID);
}
double TrackerInfo::rangeRmax(std::pair<uint32_t,uint32_t > DetID)
{
  return (rangeRmax_.find(DetID)!= rangeRmax_.end()?
                                  rangeRmax_[DetID]:
				                 0.);
}
double TrackerInfo::rangeRmin(std::pair<uint32_t,uint32_t > DetID)
{
  return (rangeRmin_.find(DetID)!= rangeRmin_.end()?
                                  rangeRmin_[DetID]:
				               999.);
}
double TrackerInfo::rangeZmax(std::pair<uint32_t,uint32_t > DetID)
{
  return (rangeZmax_.find(DetID)!= rangeZmax_.end()?
                                  rangeZmax_[DetID]:
				                 0.);

}
double TrackerInfo::rangeZmin(std::pair<uint32_t,uint32_t > DetID)
{
  return (rangeZmin_.find(DetID)!= rangeZmin_.end()?
                                  rangeZmin_[DetID]:
				               999.);
}
