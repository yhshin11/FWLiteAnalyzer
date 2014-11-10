#include <memory>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TFile.h>
#include <TSystem.h>

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


int main(int argc, char* argv[]) 
{
  // define what muon you are using; this is necessary as FWLite is not 
  // capable of reading edm::Views
  using reco::Muon;

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
  parser.integerValue ("outputEvery") =   10;
  parser.stringValue  ("outputFile" ) = "analyzeFWLiteHistograms.root";

  // parse arguments
  parser.parseArguments (argc, argv);
  int maxEvents_ = parser.integerValue("maxEvents");
  unsigned int outputEvery_ = parser.integerValue("outputEvery");
  std::string outputFile_ = parser.stringValue("outputFile");
  std::vector<std::string> inputFiles_ = parser.stringVector("inputFiles");

  // book a set of histograms
  fwlite::TFileService fs = fwlite::TFileService(outputFile_.c_str());
  TFileDirectory dir = fs.mkdir("analyzeBasicPat");
  TH1D* pfMetT1Pt_  = dir.make<TH1D>("pfMetT1Pt"  , "pt"  ,   200,   0., 200.);
  TH1D* pfMetT1CHSPt_  = dir.make<TH1D>("pfMetT1CHSPt"  , "pt"  ,   200,   0., 200.);
  TH1D* pfMetT1PtDiff_  = dir.make<TH1D>("pfMetT1PtDiff"  , "pt"  ,   200,   0., 20.);
  TH2D* pfMetT1PtComparison_  = dir.make<TH2D>("pfMetT1PtComparison"  , "pt"  ,   200,   0., 200., 200, 0., 200.);

	std::auto_ptr< double > pfMetT1Pt           (new double);
	std::auto_ptr< double > pfMetT1CHSPt           (new double);

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
	edm::EventBase const & event = ev;
	// break loop if maximal number of events is reached 
	if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
	// simple event counter
	if(outputEvery_!=0 ? (ievt>0 && ievt%outputEvery_==0) : false) 
	  std::cout << "  processing event: " << ievt << std::endl;

	// Handle to the muon collection
	edm::Handle<std::vector<Muon> > muons;
	event.getByLabel(std::string("muons"), muons);
	edm::Handle<reco::PFMETCollection> pfMetT1Handle;
  event.getByLabel(std::string("pfMetT1"), pfMetT1Handle);
	edm::Handle<reco::PFMETCollection> pfMetT1CHSHandle;
  event.getByLabel(std::string("pfMetT1CHS"), pfMetT1CHSHandle);
	
	// loop muon collection and fill histograms
	*pfMetT1Pt = pfMetT1Handle.product()->front().pt(); 
	*pfMetT1CHSPt = pfMetT1CHSHandle.product()->front().pt(); 
  pfMetT1Pt_ ->Fill(*pfMetT1Pt);
  pfMetT1CHSPt_->Fill(*pfMetT1CHSPt);
  pfMetT1PtDiff_->Fill(*pfMetT1CHSPt-*pfMetT1Pt);
	pfMetT1PtComparison_->Fill(*pfMetT1Pt, *pfMetT1CHSPt);

      }  
      // close input file
      inFile->Close();
    }
    // break loop if maximal number of events is reached:
    // this has to be done twice to stop the file loop as well
    if(maxEvents_>0 ? ievt+1>maxEvents_ : false) break;
  }
  return 0;
}
