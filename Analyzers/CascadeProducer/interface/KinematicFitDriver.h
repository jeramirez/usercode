#ifndef RECOVERTEX_KINEMATICFITFRIVER_H
#define RECOVERTEX_KINEMATICFITDRIVER_H
// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      KinematicFitDriver
//
/**\class KinematicFitDriver KinematicFitDriver.h Analyzers/CascadeProducer/interface/KinematicFitDriver.h

 Description: Tool for Wrapping the Kinematic fitting process

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: KinematicFitDriver.h,v 1.8 2013/05/17 17:21:13 jramirez Exp $
//
//

// system include files
#include <memory>
//references
#include "DataFormats/Candidate/interface/Candidate.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"

//
// class declaration
//

class KinematicFitDriver {
   public:
      explicit KinematicFitDriver(std::vector<RefCountedKinematicParticle>& Particle
                                 ,const std::string Name
                                 );
      ~KinematicFitDriver();
//    access functions
      bool isValid(){return IsValid; }
      //invariant mass of Particle
      double mass(){return (IsValid?Mother->currentState().mass():-1.); }
      //degrees of Freedom of Particle Decay Vertex
      double ndof(){return (IsValid?MotherDecayVertex->degreesOfFreedom():-1.); }
      //Chi Squared of Particle Decay Vertex
      double chi2(){return (IsValid?MotherDecayVertex->chiSquared():-1.); }
      //Mother Particle
      RefCountedKinematicParticle RKParent(){return Mother;}
      //Mother Vertex Decay
      RefCountedKinematicVertex RKVertex(){return MotherDecayVertex;}
      //Mother Vertex Decay Error (Covariance another fortmat)
      GlobalError error(){return MotherDecayVertex->error();}
      //4-Momentum Vector of Particle and Children
      reco::Particle::LorentzVector P4();
      reco::Particle::LorentzVector P4FirstChild();
      reco::Particle::LorentzVector P4NextChild();
      //Particle Decay Vertex
      reco::Particle::Point VertexDecay();
      //Particle Vertex
      reco::Vertex Vertex();
      //Particle Covariance Matrix Decay Vertex
      reco::Vertex::CovarianceMatrix VertexCov();
      //Vertex Tree
      RefCountedKinematicTree Tree(){return VertexFitTree;};
      //Mass constraint
      void AddMassConstraint(const ParticleMass, const float);
      //Multi Track Kinematic Constraint
      void AddMultiTrackKinematicConstraint(std::vector<RefCountedKinematicParticle>& Particle, const ParticleMass mass);
      //compute Cos "alpha" 
      double CosAlpha(reco::Vertex &PrimaryVertex);
      //compute significance separation from primary vertex
      double SignificanceSeparation(reco::Vertex &PrimaryVertex);
      //compute Separatio from primary
      double Separation(reco::Vertex &PrimaryVertex);
      //compute Lxy/sigma_xy (significant separation 2d)
      double SignificanceSeparation2D(reco::Vertex &PrimaryVertex);
      //compute Lxy (separation 2d)
      double Separation2D(reco::Vertex &PrimaryVertex);
      //cos "alpha" in plane x-y 
      double CosAlpha2D(reco::Vertex &PrimaryVertex);
   private:
      // ----------member data ---------------------------		   		   
      bool IsValid;                                   //flag to check validity of fit
      const std::string NameID;                       //decay id
      RefCountedKinematicTree      VertexFitTree;     //Tree formed by list of particles
      RefCountedKinematicParticle  Mother;            //Particle formed by list of particles
      RefCountedKinematicVertex    MotherDecayVertex; //Particle decay vertex
		   
};
#endif
