/* \class VeeCountFilter
 *
 * Filters events if at least N vees (Lambdas)
 *
 * \author: J. Eduardo Ramirez, Univ. of Puerto Rico
 *
 */
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "Analyzers/CascadeProducer/interface/VeeObjectCountFilter.h"
 
 typedef VeeObjectCountFilter<
           reco::VertexCompositeCandidateCollection
         >::type VeeCountFilter;
 
DEFINE_FWK_MODULE( VeeCountFilter );

