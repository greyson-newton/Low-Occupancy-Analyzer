#include <assert.h> 
#include <memory>
#include <cmath>
#include <iostream>
#include <sstream>
#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH


// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "RecoMuon/TrackingTools/interface/MuonSegmentMatcher.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/CSCRecHit/interface/CSCRecHit2D.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigiCollection.h>
#include <DataFormats/CSCDigi/interface/CSCCorrelatedLCTDigi.h>
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"

#include "DataFormats/GEMRecHit/interface/GEMRecHitCollection.h"
#include "Geometry/GEMGeometry/interface/GEMGeometry.h"
#include "Geometry/GEMGeometry/interface/GEMEtaPartitionSpecs.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"

#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TH1D.h"

using namespace std;
using namespace edm;

class Low_Occupancy_Analyzer : public edm::EDAnalyzer {
public:
  explicit Low_Occupancy_Analyzer(const edm::ParameterSet&);
  ~Low_Occupancy_Analyzer();

private:
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void beginJob() ;
  virtual void endJob() ;

  // ----------member data ---------------------------
  TH1D *th1d_tight_muon;

  edm::Service<TFileService> fs;
  MuonServiceProxy* theService_;

  edm::EDGetTokenT<edm::View<reco::Muon> > muons_;
  edm::ESHandle<Propagator> propagator_;
  edm::ESHandle<TransientTrackBuilder> ttrackBuilder_;

  edm::ESHandle<CSCGeometry> CSCGeometry_;
  edm::EDGetTokenT<CSCRecHit2DCollection> cscRecHits_;

  double maxMuonEta_, minMuonEta_;
  bool matchMuonwithCSCRechit_;

  unsigned int minTracks;
  nAllEvents = 0.0f;
  nInChamber = 0.0f;
  nGlobal = 0.0f;
  nStandAlone = 0.0f;
  nPT = 0.0f;
  nMatchedStations = 0.0f;
  nChi2 = 0.0f;
  nDB = 0.0f;
  nTightMuon = 0.0f;
};
Low_Occupancy_Analyzer::Low_Occupancy_Analyzer(const edm::ParameterSet& iConfig)
{
  th1d_tight_muon = fs->make<TH1D>("TightMuons", "Tight Muons", 7  , 0 , 7);
  th1d_tight_muon->GetXaxis()->SetBinLabel(1, "Global");
  th1d_tight_muon->GetXaxis()->SetBinLabel(2, "StandAlone");
  th1d_tight_muon->GetXaxis()->SetBinLabel(3, "PT");
  th1d_tight_muon->GetXaxis()->SetBinLabel(4, "Stations");
  th1d_tight_muon->GetXaxis()->SetBinLabel(5, "ValidDB");
  th1d_tight_muon->GetXaxis()->SetBinLabel(6, "Chi2");  

  edm::ParameterSet serviceParameters = iConfig.getParameter<edm::ParameterSet>("ServiceParameters");
  theService_ = new MuonServiceProxy(serviceParameters);

  muons_ = consumes<View<reco::Muon> >(iConfig.getParameter<InputTag>("muons"));

  cscRecHits_ = consumes<CSCRecHit2DCollection>(iConfig.getParameter<edm::InputTag>("cscRecHits"));;
  
  minMuonEta_ =  iConfig.getUntrackedParameter<double>("minMuonEta", 1.4);
  maxMuonEta_ =  iConfig.getUntrackedParameter<double>("maxMuonEta", 2.5);
  matchMuonwithCSCRechit_ =  iConfig.getUntrackedParameter<bool>("matchMuonwithCSCRechit", false);
}

void
Low_Occupancy_Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  iSetup.get<MuonGeometryRecord>().get(CSCGeometry_);
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",ttrackBuilder_);

  theService_->update(iSetup);
  auto propagator = theService_->propagator("SteppingHelixPropagatorAny");

  edm::Handle<View<reco::Muon> > muons;
  iEvent.getByToken(muons_, muons);

      
  bool hasCSCRechitcollection = false;
  edm::Handle<CSCRecHit2DCollection> cscRecHits;
  if (matchMuonwithCSCRechit_){
      try{
          iEvent.getByToken(cscRecHits_, cscRecHits);
          hasCSCRechitcollection = true;
      }catch (cms::Exception){
        std::cout<< "Error! Can't get CSC Rechit by label. " << std::endl;
        hasCSCRechitcollection = false;
      }
  }

  for (size_t i = 0; i < muons->size(); ++i) 
  {
    nAllEvents++;
    const reco::Muon* mu = muRef.get();
    const reco::Track* muonTrack = 0;
    if ( mu->globalTrack().isNonnull() ) muonTrack = mu->globalTrack().get();
    else if ( mu->outerTrack().isNonnull()  ) muonTrack = mu->outerTrack().get();
    else 
      continue;

    std::set<float> detLists;
    reco::TransientTrack ttTrack_gt = ttrackBuilder_->build(muonTrack);
    reco::TransientTrack ttTrack = ttrackBuilder_->build(standaloneMuon);
    reco::TransientTrack ttTrack_inner = ttrackBuilder_->build(innerTrack);

    for(const auto& ch : CSCGeometry_->layers())
    {
      bool isME11 = (ch->id().station() == 1 and (ch->id().ring() == 1 or ch->id().ring() == 4));
      TrajectoryStateOnSurface tsos = propagator->propagate(ttTrack.innermostMeasurementState(),ch->surface());
      TrajectoryStateOnSurface tsos_gt = propagator->propagate(ttTrack_gt.outermostMeasurementState(),ch->surface());
      TrajectoryStateOnSurface tsos_inner = propagator->propagate(ttTrack_inner.outermostMeasurementState(),ch->surface());
      if (!tsos.isValid()) continue;
      if (!tsos_gt.isValid()) continue;
      if (!tsos_inner.isValid()) continue;

      GlobalPoint tsosGP = tsos.globalPosition();
      GlobalPoint tsosGP_gt = tsos_gt.globalPosition();
      GlobalPoint tsosGP_inner = tsos_inner.globalPosition();

      if (tsosGP.eta() * mu->eta() < 0.0) continue;
      
      if (matchMuonwithCSCRechit_ and hasCSCRechitcollection) for (auto hit = cscRecHits->begin(); hit != cscRecHits->end(); hit++) 
      {
        if ((hit)->geographicalId().det() == DetId::Detector::Muon && (hit)->geographicalId().subdetId() == MuonSubdetId::CSC)
        {
          if ((hit)->rawId() == ch->id().rawId()  and isME11) 
          {
            CSCDetId cscid((hit)->geographicalId());
            for (auto muonhit = muonTrack->recHitsBegin(); muonhit != muonTrack->recHitsEnd(); muonhit++) 
            {
              if ( (*muonhit)->rawId() == ch->id().rawId() ) 
              {
                nInChamber++;
                cout <<"muonhit CSCid "<< CSCDetId((*muonhit)->geographicalId()) << endl;
                if(mu.isGlobalMuon())
                {
                  continue;
                }   
                  th1d_tight_muon->Fill(0.5);
                  nGlobal++;

                if(mu.isStandAloneMuon())
                {
                  continue;
                }
                  th1d_tight_muon->Fill(1.5);
                  nStandAlone++;

                if((mu.pt()>20) && (mu.pt()<200))
                {
                  continue;}
                  th1d_tight_muon->Fill(2.5);
                  nPT++;

                if(mu.numberOfMatchedStations() >= 2)
                {
                  continue;
                }
                  th1d_tight_muon->Fill(3.5);
                  nMatchedStations++;

                if(mu.dB() < 0.2)
                {
                        continue;
                }
                  th1d_tight_muon->Fill(4.5);
                  nDB++;
              }
            }
          }
        } 
      }
    }
  }        
}

void Low_Occupancy_Analyzer::beginJob(){}

void Low_Occupancy_Analyzer::endJob(){}

DEFINE_FWK_MODULE(Low_Occupancy_Analyzer);
