#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TTree.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>
#include <TVector3.h>
#include "DataFormats/Math/interface/LorentzVector.h"
#include <cmath>
/// Lorentz vector
typedef math::XYZTLorentzVector LorentzVector;
/// Lorentz vector
typedef math::PtEtaPhiMLorentzVector PolarLorentzVector;


#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "PhysicsTools/FWLite/interface/CommandLineParser.h"

// DataFormats
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/CorrMETData.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

// #include "MetTesting/FWLAnalyzer/interface/TreeInfo.h"

#include <map>

using std::vector;

// -------------------------------------------------//
// Class to hold the ntuple information 
// -------------------------------------------------//

class TreeInfo
{
	public:

		// constructors and destructor
		TreeInfo();
		~TreeInfo() {};

		void Reset();
		void SetBranches(TTree* tree);

		// event level info
		int run;
		int lumiblock;
		int event;
		// std::string sample;
		// std::string dataset;
		// std::string filename;
		// bool is_real_data;
		// double scale1fb;
		// double scale1fb_cms2;
		// double lumi;
		// double xsec;
		// double nevts_aod;
		// double nevts_cms2;
		// double nevts_file;
		// double kfactor;
		// double filt_eff;
		//
		// // Gen hypothesis specific info 
		// bool is_gen_z;
		// bool is_gen_ee;
		// bool is_gen_mm;
		// LorentzVector gen_p4;
		// LorentzVector gen_lep1_p4;
		// int gen_lep1_id;
		// int gen_lep1_charge;
		// LorentzVector gen_lep2_p4;
		// int gen_lep2_id;
		// int gen_lep2_charge;
		//
		// Met info
		double pfMetT1Pt;
		double pfMetT1Phi;
		double pfMetT1SumEt;
		double pfMetT1CHSPt;
		double pfMetT1CHSPhi;
		double pfMetT1CHSSumEt;

		vector<double> muonPts;
		vector<double> muonEtas;
		vector<double> muonPhis;

		bool validZ;
		double zMass;
		double zPt;
		double zEta;
		double zPhi;

		double uPar;
		double uPerp;
		double uParChs;
		double uPerpChs;

		vector<double> jetPts;
		vector<double> jetEtas;
		vector<double> jetPhis;
		vector<double> chsJetPts;
		vector<double> chsJetEtas;
		vector<double> chsJetPhis;

		vector<bool> vtxIsFakes;
		vector<bool> vtxIsValids;
		vector<double> vtxXs;
		vector<double> vtxYs;
		vector<double> vtxZs;
		vector<double> vtxRhos;
		vector<double> vtxNdfs;
		int nVtx;
		int nGoodVtx;
};

TreeInfo::TreeInfo()
{
	Reset();
}

void TreeInfo::Reset()
{
    run       = -999;
    lumiblock = -999;
    event     = -999;
    // sample               = ""; 
    // dataset              = ""; 
    // filename             = ""; 
    // is_real_data         = false; 
    // scale1fb             = 1.0; 
    // scale1fb_cms2        = 1.0; 
    // lumi                 = 1.0; 
    // xsec                 = -999; 
    // kfactor              = -999; 
    // filt_eff             = -999; 
    // is_gen_z             = false; 
    // is_gen_ee            = false; 
    // is_gen_mm            = false; 
    // gen_p4               = LorentzVector(0, 0, 0, 0); 
    // gen_lep1_p4          = LorentzVector(0, 0, 0, 0); 
    // gen_lep1_id          = -999; 
    // gen_lep1_charge      = -999; 
    // gen_lep2_p4          = LorentzVector(0, 0, 0, 0); 
    // gen_lep2_id          = -999; 
    // gen_lep2_charge      = -999; 
		pfMetT1Pt       = -999;
		pfMetT1Phi      = -999;
		pfMetT1SumEt    = -999;
		pfMetT1CHSPt    = -999;
		pfMetT1CHSPhi   = -999;
		pfMetT1CHSSumEt = -999;

		muonPts.     clear();
		muonEtas.    clear();
		muonPhis.    clear();

		validZ = false;
		zMass = -999;
		zPt   = -999;
		zEta  = -999;
		zPhi  = -999;

		uPar  = -999;
		uPerp = -999;
		uParChs  = -999;
		uPerpChs = -999;

		jetPts.      clear();
		jetEtas.     clear();
		jetPhis.     clear();
		chsJetPts.   clear();
		chsJetEtas.  clear();
		chsJetPhis.  clear();

		nVtx = -999;
		nGoodVtx = -999;
		vtxIsFakes.   clear();
		vtxIsValids.  clear();
		vtxXs.        clear();
		vtxYs.        clear();
		vtxZs.        clear();
		vtxRhos.      clear();
		vtxNdfs.      clear();
}
    
void TreeInfo::SetBranches(TTree* tree)
{
    tree->Branch("run"       , &run       ) ;
    tree->Branch("lumiblock" , &lumiblock ) ;
    tree->Branch("event"     , &event     ) ;
    // tree.Branch("sample"               , &sample               );
    // tree.Branch("dataset"              , &dataset              );
    // tree.Branch("filename"             , &filename             );
    // tree.Branch("is_real_data"         , &is_real_data         );
    // tree.Branch("scale1fb"             , &scale1fb             );
    // tree.Branch("scale1fb_cms2"        , &scale1fb_cms2        );
    // tree.Branch("lumi"                 , &lumi                 );
    // tree.Branch("xsec"                 , &xsec                 );
    // tree.Branch("filt_eff"             , &filt_eff             );
    // tree.Branch("is_gen_z"             , &is_gen_z             );
    // tree.Branch("is_gen_ee"            , &is_gen_ee            );
    // tree.Branch("is_gen_mm"            , &is_gen_mm            );
    // tree.Branch("gen_lep1_id"          , &gen_lep1_id          );
    // tree.Branch("gen_lep1_charge"      , &gen_lep1_charge      );
    // tree.Branch("gen_lep2_id"          , &gen_lep2_id          );
    // tree.Branch("gen_lep2_charge"      , &gen_lep2_charge      );
    //
    // tree.Branch("gen_p4"      , "LorentzVector" , &gen_p4     );
    // tree.Branch("gen_lep1_p4" , "LorentzVector" , &gen_lep1_p4);
    // tree.Branch("gen_lep2_p4" , "LorentzVector" , &gen_lep2_p4);

    tree->Branch("pfMetT1Pt"       , &pfMetT1Pt       ) ;
    tree->Branch("pfMetT1Phi"      , &pfMetT1Phi      ) ;
    tree->Branch("pfMetT1SumEt"    , &pfMetT1SumEt    ) ;
    tree->Branch("pfMetT1CHSPt"    , &pfMetT1CHSPt    ) ;
    tree->Branch("pfMetT1CHSPhi"   , &pfMetT1CHSPhi   ) ;
    tree->Branch("pfMetT1CHSSumEt" , &pfMetT1CHSSumEt ) ;

    tree->Branch("muonPts"  , &muonPts  ) ;
    tree->Branch("muonEtas" , &muonEtas ) ;
    tree->Branch("muonPhis" , &muonPhis ) ;

		tree->Branch("validZ" , &validZ ) ;
		tree->Branch("zMass"  , &zMass  ) ;
		tree->Branch("zPt"    , &zPt    ) ;
		tree->Branch("zEta"   , &zEta   ) ;
		tree->Branch("zPhi"   , &zPhi   ) ;

		tree->Branch("uPar"  , &uPar  ) ;
		tree->Branch("uPerp" , &uPerp ) ;
		tree->Branch("uParChs"  , &uParChs  ) ;
		tree->Branch("uPerpChs" , &uPerpChs ) ;

    tree->Branch("jetPts"     , &jetPts     ) ;
    tree->Branch("jetEtas"    , &jetEtas    ) ;
    tree->Branch("jetPhis"    , &jetPhis    ) ;
    tree->Branch("chsJetPts"  , &chsJetPts  ) ;
    tree->Branch("chsJetEtas" , &chsJetEtas ) ;
    tree->Branch("chsJetPhis" , &chsJetPhis ) ;

		tree->Branch("nVtx", &nVtx);
		tree->Branch("nGoodVtx", &nGoodVtx);

		tree->Branch("vtxIsFakes", &vtxIsFakes);
		tree->Branch("vtxIsValids", &vtxIsValids);
		tree->Branch("vtxXs", &vtxXs);
		tree->Branch("vtxYs", &vtxYs);
		tree->Branch("vtxZs", &vtxZs);
		tree->Branch("vtxRhos", &vtxRhos);
		tree->Branch("vtxNdfs", &vtxNdfs);
}


int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  using pat::Muon;
  using pat::MuonCollection;
  using pat::JetCollection;
  using reco::PFJet;
  using reco::PFJetCollection;

  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load( "libFWCoreFWLite" );
  AutoLibraryLoader::enable();

  // initialize command line parser
  optutl::CommandLineParser parser ("Analyze FWLite Histograms");

  // set defaults
  parser.integerValue ("maxEvents"  ) = 1000;
  parser.integerValue ("outputEvery") =  100;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";

  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  // TFileDirectory dir = fs.mkdir("analyzeBasicPat");
  // TH1D* pfMetT1Pt_  = dir.make<TH1D>("pfMetT1Pt"  , "pt"  ,   200,   0., 200.);
  // TH1D* pfMetT1CHSPt_  = dir.make<TH1D>("pfMetT1CHSPt"  , "pt"  ,   200,   0., 200.);
  // TH1D* pfMetT1PtDiff_  = dir.make<TH1D>("pfMetT1PtDiff"  , "pt"  ,   200,   0., 20.);
  // TH2D* pfMetT1PtComparison_  = dir.make<TH2D>("pfMetT1PtComparison"  , "pt"  ,   200,   0., 200., 200, 0., 200.);
	// TTree* tree = dir.make<TTree>("tree", "FWLAnalyzer tree");
	TTree* tree = fs.make<TTree>("FWLtree", "FWLAnalyzer tree");

	TreeInfo info;
	info.SetBranches(tree);

  // loop the events
  int ievt=0;  
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile){
    // open input file (can be located on castor)
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if( inFile ){
      // ----------------------------------------------------------------------
      // Second Part: 
      //
      //  * loop the events in the input file 
      //  * receive the collections of interest via fwlite::Handle
      //  * fill the histograms
      //  * after the loop close the input file
      // ----------------------------------------------------------------------      
      fwlite::Event ev(inFile);
      for(ev.toBegin(); !ev.atEnd(); ++ev, ++ievt){
				info.Reset();
				edm::EventBase const & event = ev;
				// break loop if maximal number of events is reached 
				if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
				// simple event counter
				if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
					std::cout << "  processing event: " << ievt << std::endl;

				// Handle to objects
				edm::Handle<reco::PFMETCollection> pfMetT1Handle;
				event.getByLabel(std::string("pfMetT1"), pfMetT1Handle);
				edm::Handle<reco::PFMETCollection> pfMetT1CHSHandle;
				event.getByLabel(std::string("pfMetT1CHS"), pfMetT1CHSHandle);

				edm::Handle<MuonCollection> muons;
				event.getByLabel(std::string("cleanPatMuons"), muons);

				// edm::Handle<PFJetCollection> jets;
				// edm::Handle<PFJetCollection> chsJets;
				// event.getByLabel(std::string("ak4PFJets"), jets);
				// event.getByLabel(std::string("ak4PFJetsCHS"), chsJets);
				edm::Handle<pat::JetCollection> jets;
				edm::Handle<pat::JetCollection> chsJets;
				event.getByLabel(std::string("cleanPatJetsPFlowNoChs"), jets);
				event.getByLabel(std::string("cleanPatJetsPFlowChs"), chsJets);

				edm::Handle<reco::VertexCollection> vertices;
				event.getByLabel(std::string("offlinePrimaryVertices"), vertices);
				
				// Fill Tree Info
				info.run       = event.id().run();
				info.event     = event.id().event();
				info.lumiblock = event.id().luminosityBlock();

				info.pfMetT1Pt       = pfMetT1Handle.product()->front().pt();
				info.pfMetT1Phi      = pfMetT1Handle.product()->front().phi();
				info.pfMetT1SumEt    = pfMetT1Handle.product()->front().sumEt();
				info.pfMetT1CHSPt    = pfMetT1CHSHandle.product()->front().pt();
				info.pfMetT1CHSPhi   = pfMetT1CHSHandle.product()->front().phi();
				info.pfMetT1CHSSumEt = pfMetT1CHSHandle.product()->front().sumEt();

				int nMuon = muons.product()->size();
				if (nMuon >=2){
					info.validZ = true;
					// PolarLorentzVector zVector = (muons.product())[0].polarP4() + muons.product()[1].polarP4();
					// auto it = muons->begin();
					PolarLorentzVector zVector = muons.product()->at(0).polarP4() + muons.product()->at(1).polarP4();
					// it++;
					// zVector += (*it).polarP4();
					info.zMass = zVector.M2();
					info.zPt = zVector.Pt();
					info.zEta = zVector.Eta();
					info.zPhi = zVector.Phi();
					// Dummy lorentz vector for met
					auto met = pfMetT1Handle.product()->front();
					PolarLorentzVector metVector(met.pt(), 0.0, met.phi(), met.pt());
					PolarLorentzVector uVector = - zVector - metVector;
					TVector3 uT(uVector.X(), uVector.Y() ,0.0);
					TVector3 qT(zVector.X(), zVector.Y(), 0.0);
					TVector3 qTunit = qT.Unit();
					// Assign sign based on uT.phi - qT phi. This is just convention. Not physics.
					info.uPerp = std::copysign(uT.Cross(qTunit).Mag(), uT.Phi()-qT.Phi());
					info.uPar = uT.Dot(qTunit);
					auto metChs = pfMetT1CHSHandle.product()->front();
					PolarLorentzVector metChsVector(metChs.pt(), 0.0, metChs.phi(), metChs.pt());
					PolarLorentzVector uChsVector = - zVector - metChsVector;
					TVector3 uTChs(uChsVector.X(), uChsVector.Y() ,0.0);
					// Assign sign based on uT.phi - qT phi. This is just convention. Not physics.
					info.uPerpChs = std::copysign(uTChs.Cross(qTunit).Mag(), uT.Phi()-qT.Phi());
					info.uParChs = uTChs.Dot(qTunit);
					// info.uPar = uVector
				}
				info.muonPts.reserve(nMuon);
				info.muonEtas.reserve(nMuon);
				info.muonPhis.reserve(nMuon);
				for (MuonCollection::const_iterator it = muons->begin(); it != muons->end(); it++) {
					Muon currentObj = *it;
					double pt  = currentObj.   pt();
					double eta = currentObj.  eta();
					double phi = currentObj.  phi();
					info.muonPts.push_back(pt);
					info.muonEtas.push_back(eta);
					info.muonPhis.push_back(phi);
				}

				info.jetPts.reserve(jets->size());
				info.jetEtas.reserve(jets->size());
				info.jetPhis.reserve(jets->size());
				for (auto it = jets->begin(); it != jets->end(); it++) {
					auto currentObj = *it;
					double pt  = currentObj.   pt();
					double eta = currentObj.  eta();
					double phi = currentObj.  phi();
					info.jetPts.push_back(pt);
					info.jetEtas.push_back(eta);
					info.jetPhis.push_back(phi);
				}

				info.chsJetPts.reserve(chsJets->size());
				info.chsJetEtas.reserve(chsJets->size());
				info.chsJetPhis.reserve(chsJets->size());
				for (auto it = chsJets->begin(); it != chsJets->end(); it++) {
					auto currentObj = *it;
					double pt  = currentObj.   pt();
					double eta = currentObj.  eta();
					double phi = currentObj.  phi();
					info.chsJetPts.push_back(pt);
					info.chsJetEtas.push_back(eta);
					info.chsJetPhis.push_back(phi);
				}

				// edm::Handle<int> nVtx;
				// void Event::getByLabel(std::string const& moduleLabel,
				//                        std::string const& productInstanceLabel, 
				//                        edm::Handle<T>&    result) 
				// edm::InputTag nVtxTag("metFilter:nVtx");
				// event.getByLabel(nVtxTag, nVtx);
				// info.nVtx = *nVtx;
				// info.nVtx = 0;
				int nVtx = 0;
				int nGoodVtx = 0;

				info.vtxIsFakes.   reserve(vertices->size());
				info.vtxIsValids.  reserve(vertices->size());
				info.vtxXs.        reserve(vertices->size());
				info.vtxYs.        reserve(vertices->size());
				info.vtxZs.        reserve(vertices->size());
				info.vtxRhos.      reserve(vertices->size());
				info.vtxNdfs.      reserve(vertices->size());
				for (reco::VertexCollection::const_iterator it = vertices->begin(); it != vertices->end(); it++) {
					reco::Vertex currentObj = *it;
					bool isFake  = currentObj.   isFake();
					bool isValid  = currentObj.   isValid();
					double X  = currentObj.   x();
					double Y  = currentObj.   y();
					double Z  = currentObj.   z();
					double Rho  = currentObj.position().Rho();
					double Ndf  = currentObj.ndof();
					info.vtxIsFakes.   push_back(isFake);
					info.vtxIsValids.  push_back(isValid);
					info.vtxXs.        push_back(X);
					info.vtxYs.        push_back(Y);
					info.vtxZs.        push_back(Z);
					info.vtxRhos.      push_back(Rho);
					info.vtxNdfs.      push_back(Ndf);
					nVtx++;

					bool vertexIsGood = (!isFake) && (Ndf>4) && (fabs(Z<=24.0)) && (Rho<=2.0);
					if (vertexIsGood) nGoodVtx++;

				}
				info.nVtx = nVtx;
				info.nGoodVtx = nGoodVtx;

				tree->Fill();
      }  
      // close input file
      inFile->Close();
			// tree->Write();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
