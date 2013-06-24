// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      PrimaryInfo
// 
/**\class PrimaryInfo PrimaryInfo.cc Analyzers/CascadeProducer/src/PrimaryInfo.cc

 Description: Tool for Primary Vertices Handling 

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: PrimaryInfo.cc,v 1.3 2012/03/14 19:52:51 jramirez Exp $
//
//

//tools
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

#include "Analyzers/CascadeProducer/interface/PrimaryInfo.h"

//
// constructors
//
PrimaryInfo::PrimaryInfo(const edm::Event& iEvent
                        ,const std::string& PrimariesColl):
   ThereArePrimaries(false)
  ,PrimColl(PrimariesColl)
  ,PrimaryChi2(-1)
  ,PrimaryNdof(0)
  ,PrimaryProb(-1)
//  ,bestVertex(0) 
{
  iEvent.getByLabel(PrimColl, primariesHandle);
  if( primariesHandle->size() ) {
    ThereArePrimaries =true;
  }
  SetInfo();
}
PrimaryInfo::PrimaryInfo(const edm::Event& iEvent
                        ,const edm::ParameterSet& iConfig):
   ThereArePrimaries(false)
  ,PrimColl(
     iConfig.getUntrackedParameter<std::string>(
       "PrimaryCollection","offlinePrimaryVertices"))
  ,PrimaryChi2(-1)
  ,PrimaryNdof(0)
  ,PrimaryProb(-1)
//  ,bestVertex(0) 
{
  iEvent.getByLabel(PrimColl, primariesHandle);
  if( primariesHandle->size() ) {
    ThereArePrimaries =true;
  }
  SetInfo();

}

//
// destructors
//
PrimaryInfo::~PrimaryInfo()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
//
// methods
//
//default criteria for best primary vertex 
//is higher multiplicity
void PrimaryInfo::SetInfo()
{
  reco::VertexCollection::const_iterator prm;
  prm = HigherMultiplicity();
  return;
}
//criterias for Primary Selection
//select higher multiplicity
reco::VertexCollection::const_iterator  PrimaryInfo::HigherMultiplicity()
{
  if (ThereArePrimaries){
     bestVertex = primariesHandle->begin();
     unsigned int nprimTracks = 0;
     for (reco::VertexCollection::const_iterator prm = primariesHandle->begin();
          prm != primariesHandle->end(); 
	  ++prm){
		  if (  prm->isValid() 
		   && !(prm->isFake())
		   && nprimTracks < prm->tracksSize()  ){
			 nprimTracks = prm->tracksSize();
			 bestVertex  = prm; 
		  }

     }

     PrimaryChi2 = bestVertex->chi2();
     PrimaryNdof = bestVertex->ndof();
     PrimaryProb = ChiSquaredProbability(PrimaryChi2,PrimaryNdof);

  }
  return bestVertex;
}
//criteria for primary vertex 
//select higher cosalpha (most colinear with secondary)
reco::VertexCollection::const_iterator PrimaryInfo::HigherCosAlpha(reco::Vertex &SecondaryVertex, GlobalVector &p)
{
 GlobalPoint Secondary(SecondaryVertex.x(),SecondaryVertex.y(),SecondaryVertex.z());
 return HigherCosAlpha(Secondary,p);
}
reco::VertexCollection::const_iterator PrimaryInfo::HigherCosAlpha(GlobalPoint &SecondaryVertex, GlobalVector &p)
{
  if (ThereArePrimaries){
     bestVertex = primariesHandle->begin();
     double cosalpha = -1.;
     for (reco::VertexCollection::const_iterator prm = primariesHandle->begin();
          prm != primariesHandle->end(); 
	  ++prm){
                  double newcosalpha = CosAlpha(prm, SecondaryVertex, p);
		  if (  prm->isValid() 
		   && !(prm->isFake())
		   && cosalpha < newcosalpha  ){
			 cosalpha = newcosalpha;
			 bestVertex  = prm; 
		  }

     }

     PrimaryChi2 = bestVertex->chi2();
     PrimaryNdof = bestVertex->ndof();
     PrimaryProb = ChiSquaredProbability(PrimaryChi2,PrimaryNdof);

  }
  return bestVertex;
}
//compute cos "alpha"
double PrimaryInfo::CosAlpha(reco::VertexCollection::const_iterator prm, GlobalPoint &SecondaryVertex, GlobalVector &p)
{
   if (!ThereArePrimaries) return -1.; //skip if not primaries
   GlobalVector LineOfFlight (SecondaryVertex.x() - prm->x(),
                              SecondaryVertex.y() - prm->y(),
                              SecondaryVertex.z() - prm->z()
                             );
   return LineOfFlight.mag()>0?
           LineOfFlight.dot(p)/LineOfFlight.mag()/p.mag():0;
}

