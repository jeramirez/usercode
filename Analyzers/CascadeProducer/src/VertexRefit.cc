// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      VertexRefit
// 
/**\class VertexRefit VertexRefit.cc Analyzers/CascadeProducer/src/VertexRefit.cc

 Description: Tool for Refiting Vertices given a exclusion list  

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: VertexRefit.cc,v 1.2 2011/06/16 19:11:19 jramirez Exp $
//
//

#include "Analyzers/CascadeProducer/interface/VertexRefit.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

//
// constructors
//
VertexRefit::VertexRefit(reco::VertexCollection::const_iterator InputVertex
                        ,std::vector<reco::TrackRef>& ExclusionList
                        ,edm::ESHandle<MagneticField>& bFieldHandle
                        ,reco::BeamSpot& BeamSpot
                        ):
  IsValid(false)
  ,VertexRefitted(*InputVertex)
  ,numexcl(0)
{
  //check if we have a valid input vertex
  if (InputVertex->tracksSize()==0) return; //Quick exit 

  if (ExclusionList.size()==0){//if no exclusion list copy over and return
    IsValid = VertexRefitted.isValid();
    return; //nothing to do
  }

  //At this point copy the input vertex valid flag
  IsValid        = VertexRefitted.isValid();

  //Loop to exclude tracks from Vertex
  std::vector<reco::TransientTrack> newPrimVertexTracks;
  for ( std::vector<reco::TrackBaseRef >::const_iterator iprmTrack = InputVertex->tracks_begin();
        iprmTrack != InputVertex->tracks_end(); 
	++iprmTrack) {
        reco::TrackRef trackRef = iprmTrack->castTo<reco::TrackRef>();

        //look if track is in exclusion list
        bool exclude=false;
        for (std::vector<reco::TrackRef>::const_iterator iExcludedTrack = ExclusionList.begin();
             iExcludedTrack != ExclusionList.end();
             ++iExcludedTrack
            ){
            if (trackRef == *iExcludedTrack ){
              exclude = true;
              break;
            }
        }//loop exclusion list
        if (exclude) {
           numexcl++;
           //std::cout << "Excluding " << numexcl << " track(s)" << std::endl;
           continue;
        }

	reco::TransientTrack ttracktoprm(trackRef, &(*bFieldHandle) );
	newPrimVertexTracks.push_back(ttracktoprm);	    
  }//end loop over primary tracks

  // refit vertex if necessary
  if (  newPrimVertexTracks.size() > 0 
    && (InputVertex->tracksSize()!=newPrimVertexTracks.size() )
     ) {
       AdaptiveVertexFitter thePrmFitter;
       TransientVertex vtx = thePrmFitter.vertex(newPrimVertexTracks, BeamSpot);
       IsValid = vtx.isValid();
       if (vtx.isValid() ) {
	  VertexRefitted = vtx;
       }
		  
   }//if  of vertex refiting

}

//
// destructors
//
VertexRefit::~VertexRefit()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
//
// methods
//

