#ifndef RecoAlgos_CandidateProducer_h
#define RecoAlgos_CandidateProducer_h
/** \class CandidateProducer
 *
 * Framework module that produces a collection
 * of candidates from generic compoment
 *
 * \author Luca Lista, INFN
 *
 * \version $Revision: 1.4 $
 *
 * $Id: CandidateProducer.h,v 1.4 2010/02/11 00:10:53 wmtan Exp $
 *
 */
#include "CommonTools/UtilAlgos/interface/AnySelector.h"
#include "CommonTools/UtilAlgos/interface/EventSetupInitTrait.h"
#include "CommonTools/UtilAlgos/interface/MasterCollectionHelper.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace converter {
  namespace helper {
    template <typename T>
    struct CandConverter {};

    struct ConcreteCreator {
      template <typename CColl, typename Comp, typename Conv>
      static void create(size_t idx, CColl& cands, const Comp& components, Conv& converter) {
        typename Conv::Candidate c;
        typedef edm::Ref<std::vector<typename Conv::value_type> > ref_type;
        ref_type ref = components.template getConcreteRef<ref_type>(idx);
        converter.convert(ref, c);
        cands.push_back(c);
      }
    };

    struct PolymorphicCreator {
      template <typename CColl, typename Comp, typename Conv>
      static void create(size_t idx, CColl& cands, const Comp& components, Conv& converter) {
        typename Conv::Candidate* c = new typename Conv::Candidate;
        typedef edm::Ref<std::vector<typename Conv::value_type> > ref_type;
        ref_type ref = components.template getConcreteRef<ref_type>(idx);
        converter.convert(ref, *c);
        cands.push_back(c);
      }
    };

    template <typename CColl>
    struct CandCreator {
      typedef ConcreteCreator type;
    };

    template <>
    struct CandCreator<reco::CandidateCollection> {
      typedef PolymorphicCreator type;
    };
  }  // namespace helper
}  // namespace converter

template <typename TColl,
          typename CColl,
          typename Selector = AnySelector,
          typename Conv = typename converter::helper::CandConverter<typename TColl::value_type>::type,
          typename Creator = typename converter::helper::CandCreator<CColl>::type,
          typename Init = typename ::reco::modules::EventSetupInit<Selector>::type>
class CandidateProducer : public edm::stream::EDProducer<> {
public:
  /// constructor from parameter set
  CandidateProducer(const edm::ParameterSet& cfg)
      : srcToken_(consumes<TColl>(cfg.template getParameter<edm::InputTag>("src"))),
        converter_(cfg, consumesCollector()),
        selectorInit_(consumesCollector()),
        selector_(reco::modules::make<Selector>(cfg, consumesCollector())),
        initialized_(false) {
    produces<CColl>();
  }
  /// destructor
  ~CandidateProducer() override {}

  /// fillDescriptions
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("src", edm::InputTag(""));
    Conv::fillPSetDescription(desc);
    descriptions.addWithDefaultLabel(desc);
  }

private:
  /// begin job (first run)
  void beginRun(const edm::Run&, const edm::EventSetup& es) override {
    if (!initialized_) {
      converter_.beginFirstRun(es);
      initialized_ = true;
    }
  }
  /// process one event
  void produce(edm::Event& evt, const edm::EventSetup& es) override {
    edm::Handle<TColl> src;
    evt.getByToken(srcToken_, src);
    selectorInit_.init(selector_, evt, es);
    ::helper::MasterCollection<TColl> master(src, evt);
    std::unique_ptr<CColl> cands(new CColl);
    if (!src->empty()) {
      size_t size = src->size();
      cands->reserve(size);
      for (size_t idx = 0; idx != size; ++idx) {
        if (selector_((*src)[idx]))
          Creator::create(master.index(idx), *cands, master, converter_);
      }
    }
    evt.put(std::move(cands));
  }
  /// label of source collection and tag
  edm::EDGetTokenT<TColl> srcToken_;
  /// converter helper
  Conv converter_;
  /// selector
  Init selectorInit_;
  Selector selector_;
  /// particles initialized?
  bool initialized_;
};

#endif
