#ifndef RecoVertex_CascadeProducer_VeeObjectCountEventSelector_h
#define RecoVertex_CascadeProducer_VeeObjectCountEventSelector_h

/** \class VeeObjectCountEventSelector
 *
 * Selects an event if a vee collection has at least N entries
 * 
 * \author J. Eduardo Ramirez, University of Puerto Rico
 *
 * \version $Revision: 1.1 $
 *
 * $Id: VeeObjectCountEventSelector.h,v 1.1 2010/03/03 20:56:17 jramirez Exp $
 *
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/UtilAlgos/interface/CollectionFilterTrait.h"
#include "CommonTools/UtilAlgos/interface/EventSelectorBase.h"

template<typename C, 
         typename S = AnySelector,
         typename N = MinNumberSelector,
         typename CS = typename helper::CollectionFilterTrait<C, S, N>::type>
class VeeObjectCountEventSelector : public EventSelectorBase
{
 public:
  /// constructor 

  explicit VeeObjectCountEventSelector( const edm::ParameterSet & cfg ) :
    algo_( cfg.template getParameter<std::string>( "algo" ) ),
    name_( cfg.template getParameter<std::string>( "name" ) ),
    select_( reco::modules::make<S>( cfg ) ),
    sizeSelect_( reco::modules::make<N>( cfg ) ) {
  }
 
  bool operator()(edm::Event& evt, const edm::EventSetup&) {
    edm::Handle<C> source;
    evt.getByLabel( algo_, name_, source);
    return CS::filter( * source, select_, sizeSelect_ );
  }
 
 private:
  /// algorithm collection label 
  std::string algo_;
  /// name collection label 
  std::string name_;
 
  /// object filter

  S select_;
 
  /// minimum number of entries in a collection

  N sizeSelect_;
};

#endif


