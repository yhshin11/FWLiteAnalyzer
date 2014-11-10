#include <string>
#include <map>
#include "TTree.h"

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
        // int run;
        // int ls;
        // int evt;
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
};

