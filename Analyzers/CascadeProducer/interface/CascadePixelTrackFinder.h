#ifndef CASCADEPRODUCER_PIXELTRACKFINDER_H
#define CASCADEPRODUCER_PIXELTRACKFINDER_H
// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadePixelTrackFinder
//
/**\class CascadePixelTrackFinder CascadePixelTrackFinder.h Analyzers/CascadeProducer/interface/CascadePixelTrackFinder.h

 Description: Gives a list of tracks with pixel info which are identified as Xi/Omega particles

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadePixelTrackFinder.h,v 1.3 2011/10/25 18:38:34 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade collection
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track

//bfield
//#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"


//
// class declaration
//

class CascadePixelTrackFinder {
   public:
      explicit CascadePixelTrackFinder();
      explicit CascadePixelTrackFinder(const edm::ParameterSet& iConfig
                                      ,const edm::Event& iEvent
		                      ,const edm::EventSetup& iSetup);
      explicit CascadePixelTrackFinder(edm::InputTag
                                      ,std::string
                                      ,std::string
                                      ,const edm::Event& iEvent
		                      ,const edm::EventSetup& iSetup);
      ~CascadePixelTrackFinder();
      bool IsValid(){return isValid;};
      reco::TrackRef Piontrack(unsigned int casindex){
        return (dynamic_cast<reco::RecoChargedCandidate *>
                                ( theCascades[casindex].daughter(0) )
	       )->track();
      }
      reco::TrackRef ProtonVeetrack(unsigned int casindex){
        return (dynamic_cast<reco::RecoChargedCandidate *>
                                ( (dynamic_cast<reco::VertexCompositeCandidate*>
					( theCascades[casindex].daughter(1) )
				  )->daughter(0) 
				)
	       )->track();

      }
      reco::TrackRef PionVeetrack(unsigned int casindex){
        return (dynamic_cast<reco::RecoChargedCandidate *>
                                ( (dynamic_cast<reco::VertexCompositeCandidate*>
					( theCascades[casindex].daughter(1) )
				  )->daughter(1) 
				)
	       )->track();

      }
      const GlobalPoint XiVertexDecay(unsigned int casindex){
        return (GlobalPoint(theCascades[casindex].vx(),
                            theCascades[casindex].vy(),
                            theCascades[casindex].vz()
                           )
               ); 
      }
      int NumPixelXiTracks(unsigned int casindex){
        return NumberXiTracks.find(casindex)!=NumberXiTracks.end()?
	                                 NumberXiTracks[casindex]:0;};
      double BestCL(unsigned int casindex){
        return XiBestCL.find(casindex)!=XiBestCL.end()?
	                                 XiBestCL[casindex]:0;};
      double BestSignificanceImpactParameterAtPixelTrack(unsigned int casindex){//closest distance of xitrack to vertex decay
        return XiBestImpactParameter.find(casindex)!=XiBestImpactParameter.end()?
	                                 XiBestImpactParameter[casindex]:1999.;};
      double FirstBestImpactParameter(unsigned int casindex){
        std::vector<double> SortedList=XiImpactParameterList[casindex];
        std::sort(SortedList.begin(),SortedList.end());
        return XiImpactParameterList.find(casindex)!=XiImpactParameterList.end()
	     &&XiImpactParameterList[casindex].size()>0?
	                             SortedList[0]:1999.;};
      double SecondBestImpactParameter(unsigned int casindex){
        std::vector<double> SortedList=XiImpactParameterList[casindex];
        std::sort(SortedList.begin(),SortedList.end());
        return XiImpactParameterList.find(casindex)!=XiImpactParameterList.end()
	     &&XiImpactParameterList[casindex].size()>1?
//	                             XiImpactParameterList[casindex][1]:1999.;};
	                             SortedList[1]:1999.;};
      double LastBestImpactParameter(unsigned int casindex){
        std::vector<double> SortedList=XiImpactParameterList[casindex];
        int dim = SortedList.size()-1;
        std::sort(SortedList.begin(),SortedList.end());
        return XiImpactParameterList.find(casindex)!=XiImpactParameterList.end()
	     &&XiImpactParameterList[casindex].size()>0?
	                             SortedList[dim]:1999.;};
      unsigned int size_i(unsigned int casindex){
        return XiImpactParameterList.find(casindex)!=XiImpactParameterList.end()
	     &&XiImpactParameterList[casindex].size()>0?
	                             XiImpactParameterList[casindex].size():0;};
      double ImpactParameter_i(unsigned int casindex,unsigned int i){
        return XiImpactParameterList.find(casindex)!=XiImpactParameterList.end()
	     &&XiImpactParameterList[casindex].size()>i?
	                             XiImpactParameterList[casindex][i]:1999.;};
      double CL_i(unsigned int casindex,unsigned int i){
        return XiCLList.find(casindex)!=XiCLList.end()
	     &&XiCLList[casindex].size()>i?
	                             XiCLList[casindex][i]:1999.;};
      void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
      void LoadTracks(const edm::Event&,edm::InputTag);
      void LoadCascades(const edm::Event&,std::string,std::string);

   private:
      // ----------member data ---------------------------
      bool isValid;                                            //flag of validity
      bool TracksAreLoaded;                                    //flag of tracks
      bool CascadesAreLoaded ;                                 //flag of cascades
      int  hasPixelTrack;                                      //Number of Pixel tracks if non cero
      std::map<unsigned int,int>  pionpixhits;                 //Number of Pixel hits in Xi/Omega daughter track
      std::vector<reco::Track> theTracks;                      //tracks
      std::map<unsigned int,reco::TrackRef> TrkRefs;           //track reference
      std::vector<reco::VertexCompositeCandidate> theCascades; //cascades
      std::map<unsigned int,std::vector<reco::TrackRef> >XiTracks;//List of tracks linked to tracks
      std::map<unsigned int,reco::TrackRef >BestXiTracks;//List of tracks linked to tracks
      std::map<unsigned int,int >NumberXiTracks;//Number of pixel tracks linked to Xitracks
      std::map<unsigned int,double >XiBestCL;   //Confidence Level (prob)
      std::map<unsigned int,double >XiBestImpactParameter;   //Impact parameter to Xi decay vertex
      std::map<unsigned int,std::vector<double> >XiImpactParameterList;   //Impact parameter to Xi decay vertex
      std::map<unsigned int,std::vector<double> >XiCLList;   //CL of Xi decay vertex
};
#endif
