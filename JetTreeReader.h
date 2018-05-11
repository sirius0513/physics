#ifndef __JETTREETREEREADER_H__
#define __JETTREETREEREADER_H__

#include <TChain.h>
#include <TFile.h>
#include <TLorentzVector.h>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <string>
#include <map>

// Fixed size dimensions of array or collections stored in the TTree if any.

class JetTreeReader {

public :
  typedef std::map<unsigned int,std::map<unsigned int,long int> > map_t;
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  JetTreeReader(TTree *tree=0,bool hastruth=true, bool hasextra=true);
  virtual ~JetTreeReader();
  virtual void     Init(TTree *tree);
  static void BuildEventMap(JetTreeReader& reader, map_t& index);
  static void MatchEvents(JetTreeReader& overlay, JetTreeReader& bkgr, const map_t& index, TFile* fout, Long64_t nmax=-1);
  static void RunLoop(JetTreeReader& overlay, JetTreeReader& bkgr, TFile* fout, Long64_t nmax=-1);
  TTree* GetTree() {return fChain;};
  unsigned int GetNJets(unsigned int iColl) const;
  void GetJet(unsigned int iColl, unsigned int ijet, TLorentzVector& vec) const;
//  float PassedFR(unsigned int iColl, unsigned int ijet) const;
  static std::string GetCollectionName(unsigned int iColl);

  // Declaration of leaf types
  UInt_t          RunNumber;
  UInt_t          EventNumber;
  UInt_t          lbn;
  Float_t         mbtime_timeA;
  Float_t         mbtime_timeC;
  std::vector<float>   *cccEt_Et_Eh_FCal;
  std::vector<float>   *cccEt_p_Eh_FCal;
  std::vector<float>   *cccEt_n_Eh_FCal;
  std::vector<float>   *antikt4Track_pt;
  std::vector<float>   *antikt4Track_m;
  std::vector<float>   *antikt4Track_eta;
  std::vector<float>   *antikt4Track_phi;
  std::vector<float>   *EMFRClusterJets_pt;
  std::vector<float>   *EMFRClusterJets_m;
  std::vector<float>   *EMFRClusterJets_eta;
  std::vector<float>   *EMFRClusterJets_phi;

  std::vector<float>   *antikt2HIItrEMFR_pt;
  std::vector<float>   *antikt2HIItrEMFR_m;
  std::vector<float>   *antikt2HIItrEMFR_eta;
  std::vector<float>   *antikt2HIItrEMFR_phi;
  std::vector<float>   *antikt2HIItrEMFR_UncalibEt;
  std::vector<float>   *antikt2HIItrEMFR_PassedFR;

  std::vector<float>   *antikt4HIItrEMFR_pt;
  std::vector<float>   *antikt4HIItrEMFR_m;
  std::vector<float>   *antikt4HIItrEMFR_eta;
  std::vector<float>   *antikt4HIItrEMFR_phi;
  std::vector<float>   *antikt4HIItrEMFR_UncalibEt;
  std::vector<float>   *antikt4HIItrEMFR_PassedFR;
  std::vector<float>   *antikt2Truth_pt;
  std::vector<float>   *antikt2Truth_m;
  std::vector<float>   *antikt2Truth_eta;
  std::vector<float>   *antikt2Truth_phi;
  UInt_t          trig_DB_SMK;
  UInt_t          trig_DB_L1PSK;
  UInt_t          trig_DB_HLTPSK;
  std::vector<unsigned int> *trig_L1_TAV;
  std::vector<short>   *trig_L2_passedPhysics;
  std::vector<short>   *trig_EF_passedPhysics;
  std::vector<unsigned int> *trig_L1_TBP;
  std::vector<unsigned int> *trig_L1_TAP;
  std::vector<short>   *trig_L2_passedRaw;
  std::vector<short>   *trig_EF_passedRaw;
  std::vector<short>   *trig_L2_resurrected;
  std::vector<short>   *trig_EF_resurrected;
  std::vector<short>   *trig_L2_passedThrough;
  std::vector<short>   *trig_EF_passedThrough;
  std::vector<float>   *trig_EF_jet_pt;
  std::vector<float>   *trig_EF_jet_m;
  std::vector<float>   *trig_EF_jet_eta;
  std::vector<float>   *trig_EF_jet_phi;
  Int_t           vx_n;
  std::vector<int>     *vx_nTracks;
  std::vector<float>   *vx_sumPt;
  int                     mu_muid_n;
  std::vector<float>      *mu_muid_px;
  std::vector<float>      *mu_muid_py;
  std::vector<float>      *mu_muid_pz;
  std::vector<float>      *mu_muid_pt;
  std::vector<float>      *mu_muid_eta;
  std::vector<float>      *mu_muid_phi;
  std::vector<float>      *mu_muid_E;
  std::vector<float>      *mc_gen_E;
  std::vector<int>        *mc_gen_type;
  std::vector<int>        *mc_gen_mothertype;
  std::vector<int>        *mc_gen_motherbarcode;
  std::vector<float>      *mc_gen_pt;
  std::vector<float>      *mc_gen_eta;
  std::vector<float>      *mc_gen_phi;
  std::vector<int>        *mc_gen_barcode;
  int                     mc_n;

  // List of branches
  TBranch        *b_RunNumber;   //!
  TBranch        *b_EventNumber;   //!
  TBranch        *b_lbn;   //!
  TBranch        *b_mbtime_timeA;   //!
  TBranch        *b_mbtime_timeC;   //!
  TBranch        *b_cccEt_Et_Eh_FCal;   //!
  TBranch        *b_cccEt_p_Eh_FCal;   //!
  TBranch        *b_cccEt_n_Eh_FCal;   //!
  TBranch        *b_antikt4Track_pt;   //!
  TBranch        *b_antikt4Track_m;   //!
  TBranch        *b_antikt4Track_eta;   //!
  TBranch        *b_antikt4Track_phi;   //!
  TBranch        *b_EMFRClusterJets_pt;   //!
  TBranch        *b_EMFRClusterJets_m;   //!
  TBranch        *b_EMFRClusterJets_eta;   //!
  TBranch        *b_EMFRClusterJets_phi;   //!

  TBranch        *b_antikt2HIItrEMFR_pt;   //!
  TBranch        *b_antikt2HIItrEMFR_m;   //!
  TBranch        *b_antikt2HIItrEMFR_eta;   //!
  TBranch        *b_antikt2HIItrEMFR_phi;   //!
  TBranch        *b_antikt2HIItrEMFR_UncalibEt;   //!
  TBranch        *b_antikt2HIItrEMFR_PassedFR;   //!

  TBranch        *b_antikt4HIItrEMFR_pt;   //!
  TBranch        *b_antikt4HIItrEMFR_m;   //!
  TBranch        *b_antikt4HIItrEMFR_eta;   //!
  TBranch        *b_antikt4HIItrEMFR_phi;   //!
  TBranch        *b_antikt4HIItrEMFR_UncalibEt;   //!
  TBranch        *b_antikt4HIItrEMFR_PassedFR;   //!
  TBranch        *b_antikt2Truth_pt;   //!
  TBranch        *b_antikt2Truth_m;   //!
  TBranch        *b_antikt2Truth_eta;   //!
  TBranch        *b_antikt2Truth_phi;   //!
  TBranch        *b_trig_DB_SMK;   //!
  TBranch        *b_trig_DB_L1PSK;   //!
  TBranch        *b_trig_DB_HLTPSK;   //!
  TBranch        *b_trig_L1_TAV;   //!
  TBranch        *b_trig_L2_passedPhysics;   //!
  TBranch        *b_trig_EF_passedPhysics;   //!
  TBranch        *b_trig_L1_TBP;   //!
  TBranch        *b_trig_L1_TAP;   //!
  TBranch        *b_trig_L2_passedRaw;   //!
  TBranch        *b_trig_EF_passedRaw;   //!
  TBranch        *b_trig_L2_resurrected;   //!
  TBranch        *b_trig_EF_resurrected;   //!
  TBranch        *b_trig_L2_passedThrough;   //!
  TBranch        *b_trig_EF_passedThrough;   //!
  TBranch        *b_trig_EF_jet_pt;   //!
  TBranch        *b_trig_EF_jet_m;   //!
  TBranch        *b_trig_EF_jet_eta;   //!
  TBranch        *b_trig_EF_jet_phi;   //!
  TBranch        *b_vx_n;   //!
  TBranch        *b_vx_nTracks;   //!
  TBranch        *b_vx_sumPt;   //!
  TBranch        *b_mu_muid_n;
  TBranch        *b_mu_muid_E;
  TBranch        *b_mu_muid_pt;
  TBranch        *b_mu_muid_px;
  TBranch        *b_mu_muid_py;
  TBranch        *b_mu_muid_pz;
  TBranch        *b_mu_muid_eta;
  TBranch        *b_mu_muid_phi;
  TBranch        *b_mc_n;
  TBranch        *b_mc_gen_E;
  TBranch        *b_mc_gen_type;
  TBranch        *b_mc_gen_pt;
  TBranch        *b_mc_gen_eta;
  TBranch        *b_mc_gen_phi;
  TBranch        *b_mc_gen_mothertype;
  TBranch        *b_mc_gen_motherbarcode;
  TBranch        *b_mc_gen_barcode;
   
private:
  bool m_has_truth;
  bool m_has_extra_collections;
};

#endif
