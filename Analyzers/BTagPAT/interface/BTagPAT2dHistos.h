#ifndef PatAlgos_BTagPAT2dHistos_H_
#define PatAlgos_BTagPAT2dHistos_H_

// -*- C++ -*-
//
// Package:    BTagPAT
// Class:      BTagPAT2dHistos
// 
/**\class BTagPAT2dHistos BTagPAT2dHistos.h

 Description: <Define and Fill a set of 2d histograms depending on flavor>

 Implementation:
 
 Create a container of 2d histograms. 
*/
//
// Original Author:  J. E. Ramirez
//
// $Id: BTagPAT2dHistos.h,v 1.1 2008/07/25 23:23:32 jramirez Exp $


// system include files
#include <memory>

// user include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "PhysicsTools/UtilAlgos/interface/TFileDirectory.h"
#include <string>
#include <map>
#include "TH2D.h"

//
// class declaration
//

class BTagPAT2dHistos {
   public:
      explicit BTagPAT2dHistos(const edm::ParameterSet& iConfig);
      BTagPAT2dHistos(const int pnxbin, const double pxmin, const double pxmax,
                      const int pnybin, const double pymin, const double pymax);
      BTagPAT2dHistos(const int pdim,
                      const int pnybin, const double pymin, const double pymax);

      ~BTagPAT2dHistos();

      void Set(const std::string& flavor, TFileDirectory& ptreldir, 
               const std::string& inprefix="n_", 
               const std::string& histotitle= "jet pt vs p_{Trel} [GeV/c]");
      void Set(const std::string& flavor, TFileDirectory& ptreldir,
               double*  , 
               const std::string& inprefix="n_", 
               const std::string& histotitle= "jet pt vs p_{Trel} [GeV/c]"
               );
      void Sumw2();
      void Fill(double&, double&, const std::string& );
      void SetPrefix(const std::string& inprefix="n_"){ prefix= inprefix;};
   private:

      // ----------member data ---------------------------

  const int nxbin,nybin, dim;
  const double xmin,xmax,ymin,ymax;
  std::map<std::string,TH2D*> histocontainer_;   // simple map to contain all histograms. Histograms are booked in the beginJob() method
  std::string prefix;

};

#endif
