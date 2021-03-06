#ifndef RemapRawDataCollection_H
#define RemapRawDataCollection_H


/** \class RemapRawDataCollection
 *
 */

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "DataFormats/FEDRawData/interface/FEDRawDataCollection.h" 
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RemapRawDataCollection: public edm::stream::EDProducer<> {
public:
    
    ///Constructor
    RemapRawDataCollection(const edm::ParameterSet& pset);
    
    ///Destructor
    ~RemapRawDataCollection() override;
    
    void produce(edm::Event & e, const edm::EventSetup& c) override; 
    
    void CopyFEDRawData(int, const FEDRawData&, int, FEDRawData&);

private:

    typedef std::vector<edm::InputTag>::const_iterator tag_iterator_t;
    typedef std::vector<edm::EDGetTokenT<FEDRawDataCollection> >::const_iterator tok_iterator_t;
    typedef std::vector<edm::ParameterSet> Parameters;

    std::vector<edm::InputTag> inputTags_ ;
    std::vector<edm::EDGetTokenT<FEDRawDataCollection> > inputTokens_;
    Parameters mappingList_;
    std::map<int,std::vector<int> > mapping_;
    std::set<int> dests_;
    int  verbose_ ;

};

#endif
