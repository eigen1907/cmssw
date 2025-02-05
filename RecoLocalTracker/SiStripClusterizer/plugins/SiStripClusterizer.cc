#include "RecoLocalTracker/SiStripClusterizer/plugins/SiStripClusterizer.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "RecoLocalTracker/SiStripClusterizer/interface/StripClusterizerAlgorithmFactory.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/transform.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

SiStripClusterizer::SiStripClusterizer(const edm::ParameterSet& conf)
    : inputTags(conf.getParameter<std::vector<edm::InputTag>>("DigiProducersList")),
      algorithm(StripClusterizerAlgorithmFactory::create(consumesCollector(),
                                                         conf.getParameter<edm::ParameterSet>("Clusterizer"))) {
  produces<edmNew::DetSetVector<SiStripCluster>>();
  inputTokens = edm::vector_transform(
      inputTags, [this](edm::InputTag const& tag) { return consumes<edm::DetSetVector<SiStripDigi>>(tag); });
}

void SiStripClusterizer::produce(edm::Event& event, const edm::EventSetup& es) {
  auto output = std::make_unique<edmNew::DetSetVector<SiStripCluster>>();
  output->reserve(10000, 4 * 10000);

  edm::Handle<edm::DetSetVector<SiStripDigi>> inputOld;
  //   edm::Handle< edmNew::DetSetVector<SiStripDigi> >  inputNew;

  algorithm->initialize(es);

  for (auto const& token : inputTokens) {
    if (findInput(token, inputOld, event))
      algorithm->clusterize(*inputOld, *output);
    //     else if( findInput( tag, inputNew, event) ) algorithm->clusterize(*inputNew, *output);
    else
      edm::LogError("Input Not Found") << "[SiStripClusterizer::produce] ";  // << tag;
  }

  LogDebug("Output") << output->dataSize() << " clusters from " << output->size() << " modules";
  output->shrink_to_fit();
  event.put(std::move(output));
}

void SiStripClusterizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<std::vector<edm::InputTag>>("DigiProducersList",
                                       {edm::InputTag("siStripDigis", "ZeroSuppressed"),
                                        edm::InputTag("siStripZeroSuppression", "VirginRaw"),
                                        edm::InputTag("siStripZeroSuppression", "ProcessedRaw"),
                                        edm::InputTag("siStripZeroSuppression", "ScopeMode")});

  {
    edm::ParameterSetDescription ClusterizerPSet;
    StripClusterizerAlgorithmFactory::fillDescriptions(ClusterizerPSet);
    desc.add<edm::ParameterSetDescription>("Clusterizer", ClusterizerPSet);
  }

  descriptions.addWithDefaultLabel(desc);
}

template <class T>
inline bool SiStripClusterizer::findInput(const edm::EDGetTokenT<T>& tag, edm::Handle<T>& handle, const edm::Event& e) {
  e.getByToken(tag, handle);
  return handle.isValid();
}
