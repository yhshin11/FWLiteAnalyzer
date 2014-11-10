#include <string>
#include <map>
#include "TTree.h"
#include "MetTesting/FWLAnalyzer/interface/TreeInfo.h"

// C++
#include <iostream>
#include <vector>

// ROOT
// #include "TBenchmark.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TROOT.h"
// #include "TTreeCache.h"
// #include "Math/VectorUtil.h"


// -------------------------------------------------//
// Class to hold the ntuple information 
// -------------------------------------------------//

TreeInfo::TreeInfo():
		pfMetT1Pt  (-999999),
		pfMetT1Phi (-999999),
		pfMetT1SumEt (-999999),
		pfMetT1CHSPt (-999999),
		pfMetT1CHSPhi (-999999),
		pfMetT1CHSSumEt (-999999)
{
	// Reset();
}

void TreeInfo::Reset()
{
    // run                  = -999999; 
    // ls                   = -999999; 
    // evt                  = -999999; 
    // sample               = ""; 
    // dataset              = ""; 
    // filename             = ""; 
    // is_real_data         = false; 
    // scale1fb             = 1.0; 
    // scale1fb_cms2        = 1.0; 
    // lumi                 = 1.0; 
    // xsec                 = -999999; 
    // kfactor              = -999999; 
    // filt_eff             = -999999; 
    // is_gen_z             = false; 
    // is_gen_ee            = false; 
    // is_gen_mm            = false; 
    // gen_p4               = LorentzVector(0, 0, 0, 0); 
    // gen_lep1_p4          = LorentzVector(0, 0, 0, 0); 
    // gen_lep1_id          = -999999; 
    // gen_lep1_charge      = -999999; 
    // gen_lep2_p4          = LorentzVector(0, 0, 0, 0); 
    // gen_lep2_id          = -999999; 
    // gen_lep2_charge      = -999999; 
		pfMetT1Pt = -999999;
		pfMetT1Phi = -999999;
		pfMetT1SumEt = -999999;
		pfMetT1CHSPt = -999999;
		pfMetT1CHSPhi = -999999;
		pfMetT1CHSSumEt = -999999;
}
    
void TreeInfo::SetBranches(TTree* tree)
{
    // tree.Branch("run"                  , &run                  );
    // tree.Branch("ls"                   , &ls                   );
    // tree.Branch("evt"                  , &evt                  );
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

    tree->Branch("pfMetT1Pt"          , &pfMetT1Pt          );
    tree->Branch("pfMetT1Phi"          , &pfMetT1Phi          );
    tree->Branch("pfMetT1SumEt"          , &pfMetT1SumEt          );
    tree->Branch("pfMetT1CHSPt"          , &pfMetT1CHSPt          );
    tree->Branch("pfMetT1CHSPhi"          , &pfMetT1CHSPhi          );
    tree->Branch("pfMetT1CHSSumEt"          , &pfMetT1CHSSumEt          );
}


