#ifndef RecoVertex_CascadeProducer_VeeObjectCountFilter_h
#define RecoVertex_CascadeProducer_VeeObjectCountFilter_h

/** \class VeeObjectCountFilter
 *
 * Filters an event if a vee collection has at least N entries
 * 
 * \author J. Eduardo Ramirez, INFN
 *
 * \version $Revision: 1.2 $
 *
 * $Id: VeeObjectCountFilter.h,v 1.2 2011/05/17 19:03:47 jramirez Exp $
 *
 */

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"
#include "CommonTools/UtilAlgos/interface/CollectionFilterTrait.h"
#include "CommonTools/UtilAlgos/interface/EventSelectorAdapter.h"
#include "Analyzers/CascadeProducer/interface/VeeObjectCountEventSelector.h"

template<typename C, 
         typename S = AnySelector,
         typename N = MinNumberSelector,
         typename CS = typename helper::CollectionFilterTrait<C, S, N>::type>
struct VeeObjectCountFilter {
  typedef EventSelectorAdapter< VeeObjectCountEventSelector<C, S, N, CS> > type;
};

#endif
