#ifndef CASCADEPRODUCER_CANDIDATES_H
#define CASCADEPRODUCER_CANDIDATES_H
// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeCandidates
//
/**\class CascadeCandidates CascadeCandidates.h Analyzers/CascadeProducer/interface/CascadeCandidates.h

 Description: Select Xi/Omega candidates based on configurable cuts

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: CascadeCandidates.h,v 1.5 2012/07/18 18:46:10 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "FWCore/Framework/interface/ESHandle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
// user include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"       //track   collection
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"    //cascade collection
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"    //cascade daughter track

//tools
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//bfield
#include "MagneticField/Engine/interface/MagneticField.h"

//primary
#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"

//
// class declaration
//

class CascadeCandidates {
   public:
      explicit CascadeCandidates(const edm::ParameterSet& iConfig
                                ,const edm::Event& iEvent
		                ,const edm::EventSetup& iSetup);
      explicit CascadeCandidates(edm::InputTag
                                ,std::string
                                ,std::string
                                ,const edm::Event& iEvent
		                ,const edm::EventSetup& iSetup);
      explicit CascadeCandidates(const edm::ParameterSet& iConfig
                                ,const edm::Event& iEvent);
      explicit CascadeCandidates(const edm::Event& iEvent);
      ~CascadeCandidates();      
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
      GlobalVector MomentumAtDecay(unsigned int casindex){
        return (GlobalVector(theCascades[casindex].momentum().x(),
                             theCascades[casindex].momentum().y(),
                             theCascades[casindex].momentum().z()
                            ) 
               );
      }
      double RhoTrkDaughter(unsigned int casindex){
           return  (Piontrack(casindex))->innerPosition().rho();
      }
      double RhoVertexDecay(unsigned int casindex){
           return  XiVertexDecay(casindex).perp();
      }
      reco::Vertex PrimaryRefitted(unsigned int casindex);
      reco::TransientTrack GetXiTransientTrack(unsigned int casindex);
      double delsig(unsigned int casindex);
      double delsig(unsigned int casindex, reco::Vertex &refitVertexPrim);
      double SignificanceImpactParameter3DAtPrimary(unsigned int casindex);
      double SignificanceImpactParameter3DAtVertex(unsigned int casindex,reco::Vertex &PrmVtx);
      double SignificanceAbsoluteImpactParameter3D(
                               reco::TransientTrack &TTrack,
                               reco::Vertex &Vertex) const;
      double XiCL(unsigned int casindex){
        double  chi2 =  theCascades[casindex].vertexChi2();
        double  ndof =  theCascades[casindex].vertexNdof();
        return ndof>1?ChiSquaredProbability(chi2,ndof-1)
                     :ChiSquaredProbability(chi2,ndof);

      }
      int NumPixelXiTracks(unsigned int casindex){
        return NumberXiTracks.find(casindex)!=NumberXiTracks.end()?
	                                 NumberXiTracks[casindex]:0;};
      double BestCL(unsigned int casindex){
        return XiBestCL.find(casindex)!=XiBestCL.end()?
	                                 XiBestCL[casindex]:0;};
      double BestSignificanceImpactParameterAtPixelTrack(unsigned int casindex){
        return XiBestImpactParameter.find(casindex)!=XiBestImpactParameter.end()?
	                                 XiBestImpactParameter[casindex]:1999.;};
      double BestImpactParameter(unsigned int casindex){//Duplicate of BestSignificanceImpactParameterAtPixelTrack
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
      double ImpactParameterPixelTrack_i(unsigned int casindex,unsigned int i){
        return XiImpactParameterList.find(casindex)!=XiImpactParameterList.end()
	     &&XiImpactParameterList[casindex].size()>i?
	                             XiImpactParameterList[casindex][i]:1999.;};
      double CLpixelTrack_i(unsigned int casindex,unsigned int i){
        return XiCLList.find(casindex)!=XiCLList.end()
	     &&XiCLList[casindex].size()>i?
	                             XiCLList[casindex][i]:1999.;};
      reco::VertexCollection::const_iterator ChoosePrimaryByCosAlpha(unsigned int casindex);
      void init(const edm::Event& iEvent,const edm::EventSetup& iSetup);
      void init(edm::InputTag
               ,std::string
               ,std::string
               ,const edm::Event& iEvent
	       ,const edm::EventSetup& iSetup);
      void analyze(const edm::Event& iEvent,const edm::EventSetup& iSetup);
      void LoadTracks(const edm::Event&,edm::InputTag);
      void LoadCascades(const edm::Event&,std::string,std::string);

   private:
      // ----------member data ---------------------------
      bool isValid;                                            //flag of validity
      bool TracksAreLoaded;                                    //flag of tracks
      bool CascadesAreLoaded ;                                 //flag of cascades
      int  hasPixelTrack;                                      //Number of Pixel tracks if non cero
      PrimaryInfo primary;                                     //Primary Vertex
      reco::BeamSpot beamSpot;                                 //BeamSpot
      edm::ESHandle<MagneticField> bFieldHandle;               //BField 
      edm::ESHandle<GlobalTrackingGeometry> globTkGeomHandle;  //Global tracker geometry
      edm::ESHandle<TransientTrackBuilder> theTTBHandle;       //Transient Track Builder
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
      std::map<unsigned int,std::vector<double> >XiCLList;   //CL of Xi decay vertex for pixeltrack
};
#endif
