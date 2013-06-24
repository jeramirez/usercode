// -*- C++ -*-
//
// Package:    CruzetMuonFilter
// Class:      CruzetMuonFilter
// 
/**\class CruzetMuonFilter CruzetMuonFilter.cc Analyzers/CruzetMuonFilter/src/CruzetMuonFilter.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  "Miguel Bonnett"
//         Created:  Thu Oct 16 14:12:58 CDT 2008
// $Id: CruzetMuonFilter.cc,v 1.2 2008/12/03 17:27:28 bonnett Exp $
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
//  class Propagator;
 
 /* C++ Headers */
 #include <string>
 #include <iostream>
 #include<fstream>

 using namespace std;
 using namespace edm;
//
// class declaration
//

class CruzetMuonFilter : public edm::EDFilter {
   public:
      explicit CruzetMuonFilter(const edm::ParameterSet&);
      ~CruzetMuonFilter();

    /* Operations */
//      bool filterCilin(const edm::Event&, const edm::EventSetup& , double radius, double max_z);
   private:
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      
      // ----------member data ---------------------------
      std::string theSTAMuonLabel; // label of muons 

      double posX_;
      double posY_;
      double posZ_;
      int nmuon_;
      
      std::ofstream fasciiFile;
      std::string fasciiFileName;
      std::ofstream filefussy;

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
CruzetMuonFilter::CruzetMuonFilter(const edm::ParameterSet& iConfig)
{
   // the name of the STA rec hits collection
   theSTAMuonLabel = iConfig.getParameter<string>("SALabel");

   posX_ = iConfig.getParameter<double>("position_x"); // Region X (cm)
   posY_ = iConfig.getParameter<double>("position_y"); // Region Y (cm)
   posZ_ = iConfig.getParameter<double>("position_z"); // Region Z (cm)
   nmuon_ = iConfig.getParameter<unsigned int>("nmuon"); 
    
   fasciiFileName = iConfig.getParameter<string>("MuonONEFilter");
   fasciiFile.open(fasciiFileName.c_str());
   filefussy.open("fussymuon.txt");

}


CruzetMuonFilter::~CruzetMuonFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
CruzetMuonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    int run_nr = iEvent.id().run();
    int ev_nr  = iEvent.id().event();
    
    int Nmuonloose = 0;
    int NrecHitYP = 0;
    int NrecHitYN = 0;
    int NrecHitZ = 0;
    int NrecHitXYP = 0;
    int NrecHitXYN = 0;
    
    double recHitsGlobX = 0; 
    double recHitsGlobY = 0;
    double recHitsGlobZ = 0;
         
//    std::string ntracker;
    bool accept = false;
    
   // Get the RecTrack collection from the event
    ESHandle<MagneticField> theMGField;
    iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

    ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
    iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);
 
    edm::Handle<reco::TrackCollection> staTracks;
    iEvent.getByLabel(theSTAMuonLabel, staTracks);

    reco::TrackCollection::const_iterator staTrack;

    
    // loop on the Muons
    for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack)
    {
        reco::TransientTrack track(*staTrack, &*theMGField, theTrackingGeometry);

        int nrechits = 0;
        double deltaTheta = staTrack->outerMomentum().Theta() - staTrack->innerMomentum().Theta();
        // loop on the RecHits
        for (trackingRecHit_iterator it = track.recHitsBegin ();  it != track.recHitsEnd (); it++)
        {
            if ((*it)->isValid ())
            {               
                nrechits++;

                const GeomDet* geomDet = theTrackingGeometry->idToDet((*it)->geographicalId());
                recHitsGlobX = geomDet->toGlobal((*it)->localPosition()).x();
                recHitsGlobY = geomDet->toGlobal((*it)->localPosition()).y();
                recHitsGlobZ = geomDet->toGlobal((*it)->localPosition()).z();
//                cout<<" RecHitsGlobalPosition " <<" x: "<<recHitsGlobX <<" y: "<<recHitsGlobY <<" z: "<<recHitsGlobZ <<endl;
                
                if (recHitsGlobY > posY_) 
                {
                    NrecHitYP++;
//                    cout<<" RecHitsGlobY  "<<recHitsGlobY<<" RecHitsY+  "<<NrecHitYP <<endl;
                }
                if (recHitsGlobY < -posY_)
                {
                    NrecHitYN++;
                }
                if (fabs(recHitsGlobZ) < posZ_)
                {
                    NrecHitZ++;
                }
                if (recHitsGlobX > posX_ && fabs(recHitsGlobY) < posY_)
                {
                    NrecHitXYP++;
                }
                if (recHitsGlobX < -posX_ && fabs(recHitsGlobY) < posY_)
                {
                    NrecHitXYN++;
                }
               
                if((recHitsGlobX >= -200) && (recHitsGlobX <= 0) && (recHitsGlobY >= -200) && (recHitsGlobY <= 0))
                {
                    filefussy<<" Fussy  "<< run_nr <<"  Event = " << ev_nr << std::endl;
                }
            }                
        }// end loop on the RecHits
        
                
//        if (nrechits > 20 && (NrecHitXYP == 0 && NrecHitXYN == 0) && (NrecHitYP > 0 && NrecHitYN > 0) && NrecHitZ > 0)
        if (nrechits > 20 && (abs(deltaTheta) < 0.25) && (NrecHitXYP == 0 && NrecHitXYN == 0) && (NrecHitYP > 0 && NrecHitYN > 0) && NrecHitZ > 0)
        {
            Nmuonloose++;           
        }
       
        
    }// end loop on the Muons

    if (Nmuonloose > nmuon_)
    {
        cout<<" Nmuonloose =  "<< Nmuonloose <<endl;
        fasciiFile<<" Nmuonloose =  "<< Nmuonloose <<endl;
        fasciiFile<<" run  "<< run_nr <<"  Event = " << ev_nr << std::endl;
        accept = true;
        return accept;
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
DEFINE_FWK_MODULE(CruzetMuonFilter);
