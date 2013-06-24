// -*- C++ -*-
//
// Package:    CascadeProducer
// Class:      KinematicFitDriver
// 
/**\class KinematicFitDriver KinematicFitDriver.cc Analyzers/CascadeProducer/src/KinematicFitDriver.cc

 Description: Tool for Wrapping process of Kinematic Fitting  

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Juan Eduardo Ramirez Vargas
//         Created:  Mon Jun 29 17:24:25 CDT 2009
// $Id: KinematicFitDriver.cc,v 1.10 2013/06/11 15:09:47 jramirez Exp $
//
//

#include "Analyzers/CascadeProducer/interface/KinematicFitDriver.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"

//
// constructors
//
KinematicFitDriver::KinematicFitDriver(std::vector<RefCountedKinematicParticle>& Particle
                                      ,const std::string Name
                                      ):
  IsValid(false)
 ,NameID(Name)
{
  //check if we have a valid input vector
  if (Particle.size()==0){
    return; //nothing to do
  }

  KinematicParticleVertexFitter fitter;
  try {
     VertexFitTree = fitter.fit(Particle);
  } catch(const std::exception& e) {
     std::cout <<  "KinematicFitDriver: fit failed and not catched by Caching Vertex" << std::endl;
     std::cout <<  "KinematicFitDriver: " << e.what() << std::endl; 
     LogDebug("KinematicFitDriver") << "caught an exception in the "
               << NameID
               << " vertex fit (uncached) see below:\n";
     LogDebug("KinematicFitDriver") << e.what() << std::endl; 
     return;
  }  
  if (!VertexFitTree->isValid()){ 
     LogDebug("KinematicFitDriver") << "caught an exception in the "
               << NameID
               << " vertex fit\n"; 
     return;
  }

 IsValid = VertexFitTree->isValid();
 VertexFitTree->movePointerToTheTop();
 Mother            = VertexFitTree->currentParticle();
 MotherDecayVertex = VertexFitTree->currentDecayVertex(); 

}

//
// destructors
//
KinematicFitDriver::~KinematicFitDriver()
{
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}
//
// methods
//
reco::Particle::LorentzVector KinematicFitDriver::P4(){
   return (IsValid?reco::Particle::LorentzVector(Mother->currentState().globalMomentum().x(),
                                                 Mother->currentState().globalMomentum().y(),
                                                 Mother->currentState().globalMomentum().z(),
                                                 Mother->currentState().kinematicParameters().energy()
                                                )
                  :reco::Particle::LorentzVector(0,0,0,0)
          );
}
reco::Particle::LorentzVector KinematicFitDriver::P4FirstChild(){
   return (VertexFitTree->movePointerToTheFirstChild()?
                   reco::Particle::LorentzVector(VertexFitTree->currentParticle()->currentState().globalMomentum().x(),
                                                 VertexFitTree->currentParticle()->currentState().globalMomentum().y(),
                                                 VertexFitTree->currentParticle()->currentState().globalMomentum().z(),
                                                 VertexFitTree->currentParticle()->currentState().kinematicParameters().energy()
                                                )
                  :reco::Particle::LorentzVector(0,0,0,0)
          );
}
reco::Particle::LorentzVector KinematicFitDriver::P4NextChild(){
   return (VertexFitTree->movePointerToTheNextChild()?
                   reco::Particle::LorentzVector(VertexFitTree->currentParticle()->currentState().globalMomentum().x(),
                                                 VertexFitTree->currentParticle()->currentState().globalMomentum().y(),
                                                 VertexFitTree->currentParticle()->currentState().globalMomentum().z(),
                                                 VertexFitTree->currentParticle()->currentState().kinematicParameters().energy()
                                                )
                  :reco::Particle::LorentzVector(0,0,0,0)
          );
}
reco::Particle::Point KinematicFitDriver::VertexDecay(){
   return (IsValid?reco::Particle::Point((*MotherDecayVertex).position().x(),
                                         (*MotherDecayVertex).position().y(),
                                         (*MotherDecayVertex).position().z()
                                        )
                  :reco::Particle::Point(0,0,0)
          );
}
//Vertex
reco::Vertex KinematicFitDriver::Vertex(){  
  return (IsValid?reco::Vertex( (*MotherDecayVertex) )
                 :reco::Vertex()
         );
}
//Vertex Covariance
reco::Vertex::CovarianceMatrix KinematicFitDriver::VertexCov(){  
  return (IsValid?reco::Vertex::CovarianceMatrix( MotherDecayVertex->error().matrix_new() )
                 :reco::Vertex::CovarianceMatrix()
         );
}
//include mass constraint
void KinematicFitDriver::AddMassConstraint(const ParticleMass mass, const float sigma_mass){
   if (!IsValid) return; //skip if not existing valid fit
   VertexFitTree->movePointerToTheTop();//constraint to Mother Particle

   KinematicConstraint * mass_constraint = new MassKinematicConstraint(mass,sigma_mass);
   KinematicParticleFitter Fitter;
   VertexFitTree = Fitter.fit(mass_constraint,VertexFitTree);
   IsValid = VertexFitTree->isValid(); //update validity  of fit output
   if (!VertexFitTree->isValid()){ 
     LogDebug("KinematicFitDriver") << "caught an exception in the "
               << NameID
               << " vertex constraint fit\n"; 
     return; //exit without updating vertex info
   }
   //Update vertex info
   VertexFitTree->movePointerToTheTop();
   Mother            = VertexFitTree->currentParticle();
   MotherDecayVertex = VertexFitTree->currentDecayVertex();   
   return;

}
//compute cos "alpha"
double KinematicFitDriver::CosAlpha(reco::Vertex &PrimaryVertex){
   if (!IsValid) return -1.; //skip if not existing valid fit
   GlobalVector LineOfFlight ((*MotherDecayVertex).position().x() - PrimaryVertex.x(),
                              (*MotherDecayVertex).position().y() - PrimaryVertex.y(),
                              (*MotherDecayVertex).position().z() - PrimaryVertex.z()
                             );
   GlobalVector MotherMomentum (Mother->currentState().globalMomentum().x(),
                                Mother->currentState().globalMomentum().y(),
                                Mother->currentState().globalMomentum().z()
                               );
   return LineOfFlight.mag()>0?
           LineOfFlight.dot(MotherMomentum)/LineOfFlight.mag()/MotherMomentum.mag():0;
}
//compute L/s
double KinematicFitDriver::SignificanceSeparation(reco::Vertex &PrimaryVertex){
   if (!IsValid) return -1.; //skip if not existing valid fit
   typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
   typedef ROOT::Math::SVector<double, 3> SVector3;

   SMatrixSym3D totalCov = PrimaryVertex.covariance() + VertexCov();

   SVector3 PrimarySep3D( (*MotherDecayVertex).position().x() - PrimaryVertex.x(),
                          (*MotherDecayVertex).position().y() - PrimaryVertex.y(),
                          (*MotherDecayVertex).position().z() - PrimaryVertex.z()
                        );

   return ROOT::Math::Similarity(totalCov, PrimarySep3D)>0?
            ROOT::Math::Dot(PrimarySep3D,PrimarySep3D)/sqrt(ROOT::Math::Similarity(totalCov, PrimarySep3D))
                                                            :-999;
}
//compute L (separation)
double KinematicFitDriver::Separation(reco::Vertex &PrimaryVertex){
   if (!IsValid) return -1.; //skip if not existing valid fit
   typedef ROOT::Math::SVector<double, 3> SVector3;

   SVector3 PrimarySep3D( (*MotherDecayVertex).position().x() - PrimaryVertex.x(),
                          (*MotherDecayVertex).position().y() - PrimaryVertex.y(),
                          (*MotherDecayVertex).position().z() - PrimaryVertex.z()
                        );

   return ROOT::Math::Mag(PrimarySep3D);
}
//compute L/s (2D x-y)
double KinematicFitDriver::SignificanceSeparation2D(reco::Vertex &PrimaryVertex){
   if (!IsValid) return -1.; //skip if not existing valid fit
   typedef ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > SMatrixSym3D;
   typedef ROOT::Math::SVector<double, 3> SVector3;

   SMatrixSym3D totalCov = PrimaryVertex.covariance() + VertexCov();

   SVector3 PrimarySep2D( (*MotherDecayVertex).position().x() - PrimaryVertex.x(),
                          (*MotherDecayVertex).position().y() - PrimaryVertex.y(),
                          0
                        );

   return ROOT::Math::Similarity(totalCov, PrimarySep2D)>0?
            ROOT::Math::Dot(PrimarySep2D,PrimarySep2D)/sqrt(ROOT::Math::Similarity(totalCov, PrimarySep2D))
                                                            :-999;
}
//compute Lxy (separation in 2 d)
double KinematicFitDriver::Separation2D(reco::Vertex &PrimaryVertex){
   if (!IsValid) return -1.; //skip if not existing valid fit
   typedef ROOT::Math::SVector<double, 3> SVector3;

   SVector3 PrimarySep2D( (*MotherDecayVertex).position().x() - PrimaryVertex.x(),
                          (*MotherDecayVertex).position().y() - PrimaryVertex.y(),
                          0
                        );

   return ROOT::Math::Mag(PrimarySep2D);
}
//compute cos "alpha" in x-y (2 dim)
double KinematicFitDriver::CosAlpha2D(reco::Vertex &PrimaryVertex){
   if (!IsValid) return -1.; //skip if not existing valid fit
   GlobalVector LineOfFlight ((*MotherDecayVertex).position().x() - PrimaryVertex.x(),
                              (*MotherDecayVertex).position().y() - PrimaryVertex.y(),
                              0
                             );
   GlobalVector MotherMomentum (Mother->currentState().globalMomentum().x(),
                                Mother->currentState().globalMomentum().y(),
                                0
                               );
   return LineOfFlight.mag()>0?
           LineOfFlight.dot(MotherMomentum)/LineOfFlight.mag()/MotherMomentum.mag():0;
}
//include MultiTrack Kinematic Constraint
void KinematicFitDriver::AddMultiTrackKinematicConstraint(std::vector<RefCountedKinematicParticle>& Particle, const ParticleMass mass){
   if (!IsValid) return; //skip if not existing valid fit
   VertexFitTree->movePointerToTheTop();//constraint to Mother Particle

   ParticleMass massconst = mass;
   MultiTrackKinematicConstraint * mass_constraint = new  TwoTrackMassKinematicConstraint(massconst);
   KinematicConstrainedVertexFitter Fitter;
   VertexFitTree = Fitter.fit(Particle,mass_constraint);
   IsValid = VertexFitTree->isValid(); //update validity  of fit output
   if (!VertexFitTree->isValid()){ 
     LogDebug("KinematicFitDriver") << "caught an exception in the "
               << NameID
               << " vertex constraint fit\n"; 
     return; //exit without updating vertex info
   }
   //Update vertex info
   VertexFitTree->movePointerToTheTop();
   Mother            = VertexFitTree->currentParticle();
   MotherDecayVertex = VertexFitTree->currentDecayVertex();   
   return;

}

