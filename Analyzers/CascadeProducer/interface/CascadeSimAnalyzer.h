// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      CascadeSimAnalyzer
//
/**\class CascadeSimAnalyzer CascadeSimAnalyzer.h Analyzers/CascadeProducer/interface/CascadeSimAnalyzer.h

 Description: Find in simulation Cascades/Omegas

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Tue Jun 23 17:52:47 CDT 2009
// $Id: CascadeSimAnalyzer.h,v 1.1 2009/11/13 17:21:33 jramirez Exp $
//
//

#ifndef RECOVERTEX__CASCADE_SIMANALYZER_H
#define RECOVERTEX__CASCADE_SIMANALYZER_H

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
//MC
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"


//
// class declaration
//

class CascadeSimAnalyzer {
  public:
     CascadeSimAnalyzer( //const edm::ParameterSet& iConfig,
	           const edm::Event& iEvent, 
		   const edm::EventSetup& iSetup);
    ~CascadeSimAnalyzer();

  bool SimTrackerIsValid() const;
  bool SimVertexIsValid() const;
  void print() const;

 private:
  std::vector<SimTrack> theSimTracks;   //Collection of simulated tracks
  std::vector<SimVertex> theSimVerts;   //Collection of simulated vertices
  const int Lambda;                     // *
  const int Xi;                         //  *
  const int Omega;                      //   > pdg ids
  const int SigmaPlus;                  //  *
  const int SigmaMinus;                 // *
  int numLambda,numALambda;             // *
  int numXi,numAXi;                     //  *
  int numOmega,numAOmega;               //   > number of particles in events
  int numSigmaP,numASigmaP;             //  *
  int numSigmaM,numASigmaM;             // *

};
 
#endif
