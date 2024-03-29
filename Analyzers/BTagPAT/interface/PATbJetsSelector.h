#ifndef PatAlgos_PATbJetsSelector_H_
#define PatAlgos_PATbJetsSelector_H_

/** 
  \class    pat::PATbJetsSelector PATbJetsSelector.h "Analysis/BTagPAT/interface/PATbJetsSelector.h"

   The PATbJetsSelector is used to clean and sort a collection of jets.
   It allows a selection based on the discriminant of a given tagger using a 
   look up table of operanting points.

*/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

class PATbJetsSelector {
public:

  PATbJetsSelector(const edm::ParameterSet& cfg);
  bool IsbTag(const pat::Jet& JetCand,
              const std::string& operpoint,
              const std::string& tagger);
  bool IsbTag(const pat::Jet& JetCand,
              const std::string& operpoint);
  bool IsbTag(const pat::Jet& JetCand);
//
private:
  std::map <std::string,std::map<std::string,double> > discCut;    //map to associate cuts and taggers
  std::vector<double> discriminantCutsLoose_;                      //list of discriminant cut per tagger
  std::vector<double> discriminantCutsMedium_;                     //list of discriminant cut per tagger
  std::vector<double> discriminantCutsTight_;                      //list of discriminant cut per tagger
  std::vector<std::string> BTagdiscriminator_;                     //list of taggers
  std::string DefaultOp_;                                          //default operating point
  std::string DefaultTg_;                                          //default taggers
};

#endif

