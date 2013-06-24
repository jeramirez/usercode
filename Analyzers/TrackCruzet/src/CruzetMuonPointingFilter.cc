// -*- C++ -*-CruzetMuonPointingFilter
//
// Package:    CruzetMuonPointingFilter
// Class:      CruzetMuonPointingFilter
// 
/**\class CruzetMuonPointingFilter CruzetMuonPointingFilter.cc Analyzers/CruzetMuonPointingFilter/src/CruzetMuonPointingFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Miguel Bonnett"
//         Created:  Thu Oct 16 18:28:03 CDT 2008
// $Id: CruzetMuonPointingFilter.cc,v 1.2 2008/12/03 17:27:28 bonnett Exp $
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

 /* Collaborating Class Header */
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GeometrySurface/interface/Cylinder.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4D.h"
#include "DataFormats/TrackingRecHit/interface/RecSegment.h"

#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"
#include "DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h"

#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"



 /* Collaborating Class Declarations */
 class Propagator;
 
 /* C++ Headers */
 #include <string>
 #include <iostream>
 #include<fstream>

 using namespace std;
 using namespace edm;
//
// class declaration
//

class CruzetMuonPointingFilter : public edm::EDFilter {
   public:
      explicit CruzetMuonPointingFilter(const edm::ParameterSet&);
      ~CruzetMuonPointingFilter();

   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      
      // ----------member data ---------------------------
      std::string theSTAMuonLabel; // label of muons 
      std::string thePropagatorName; // name of propagator to be used
    
      double theRadius;// radius of cylinder   
      double theMaxZ;// half lenght of cylinder

      std::ofstream fasciiFile;
      std::string fasciiFileName;
      
      mutable Propagator* thePropagator;
};


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
CruzetMuonPointingFilter::CruzetMuonPointingFilter(const edm::ParameterSet& iConfig)
{
   theSTAMuonLabel = iConfig.getParameter<string>("SALabel");

   thePropagatorName = iConfig.getParameter<std::string>("PropagatorName");
   thePropagator = 0;

   theRadius = iConfig.getParameter<double>("radius"); // cyl's radius (cm)
   theMaxZ = iConfig.getParameter<double>("maxZ"); // cyl's half lenght (cm)

   LogDebug("CruzetMuonFilter") << " SALabel : " << theSTAMuonLabel 
    << " Radius : " << theRadius
    << " Half lenght : " << theMaxZ;
        
   fasciiFileName = iConfig.getParameter<string>("PointigFilter");
   fasciiFile.open(fasciiFileName.c_str());
//    fasciiFile.open("filtromuon.txt");

}


CruzetMuonPointingFilter::~CruzetMuonPointingFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
CruzetMuonPointingFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//     double radius;
//     double max_z;
    int nrun = iEvent.id().run();
    int nevent  = iEvent.id().event();
    std::string ntracker;
    bool accept = false;

    if (!thePropagator)
    {
        ESHandle<Propagator> prop;
        iSetup.get<TrackingComponentsRecord>().get(thePropagatorName, prop);
        thePropagator = prop->clone();
        thePropagator->setPropagationDirection(anyDirection);
    }

    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

    ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

    // Get the RecTrack collection from the event
    Handle<reco::TrackCollection> staTracks;
    iEvent.getByLabel(theSTAMuonLabel, staTracks);

    if (staTracks->size() == 0)
    {
//   	std::cout<<"    NO theSTAMuon "<<std::endl;
        return accept;
    }

    reco::TrackCollection::const_iterator staTrack;

    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack)
    {
        reco::TransientTrack track(*staTrack, &*theMGField, theTrackingGeometry);

        TrajectoryStateOnSurface innerTSOS = track.innermostMeasurementState();

        LogDebug("CruzetMuonFilter") << " InnerTSOS " << innerTSOS;

        // Get a surface (here a cylinder of radius 1290mm) ECAL
        Cylinder::PositionType pos0;
        Cylinder::RotationType rot0;
        const Cylinder::CylinderPointer cyl = Cylinder::build(pos0, rot0, theRadius);

        TrajectoryStateOnSurface tsosAtCyl =
            thePropagator->propagate(*innerTSOS.freeState(), *cyl);

        if ( tsosAtCyl.isValid() )
        {
            LogDebug("CruzetMuonFilter") << " extrap TSOS " << tsosAtCyl;
            if (fabs(tsosAtCyl.globalPosition().z()) < theMaxZ)
            {
                 if (theRadius == 110) ntracker = "SiTrack"; 
                 if (theRadius == 15) ntracker = "Pixel  ";
                 if (theRadius == 4)    ntracker = "IP     ";
             //         std::cout<<ntracker<<"    radio = "<<radius<<" Z= "<<max_z<<"   RUN = "<<nrun<<"   Event = "<<nevent<<std::endl;
                fasciiFile<<ntracker <<  "    radio = " << theRadius << " Z= " << theMaxZ << "   RUN = " << nrun << "   Event = " << nevent << std::endl;
                accept = true;
                return accept;
            }
            else
            {
                LogDebug("CruzetMuonFilter") << " extrap TSOS z too big " << tsosAtCyl.globalPosition().z();
            }
        }
        else
        {
            LogDebug("CruzetMuonFilter") << " extrap to cyl failed ";
        }
    }
    return accept;





//    using namespace edm;
// #ifdef THIS_IS_AN_EVENT_EXAMPLE
//    Handle<ExampleData> pIn;
//    iEvent.getByLabel("example",pIn);
// #endif
// 
// #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
//    ESHandle<SetupData> pSetup;
//    iSetup.get<SetupRecord>().get(pSetup);
// #endif
//    return true;
}


//define this as a plug-in
DEFINE_FWK_MODULE(CruzetMuonPointingFilter);
