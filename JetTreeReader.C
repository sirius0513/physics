
#include "JetTreeReader.h"
#include <TAxis.h>
#include <sstream>
#include <cmath>
#include <iostream>
#include <iomanip>

JetTreeReader::JetTreeReader(TTree *tree, bool hastruth, bool hasextra) : fChain(0), m_has_truth(hastruth), m_has_extra_collections(hasextra)
{
  Init(tree);
}

JetTreeReader::~JetTreeReader()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

void JetTreeReader::Init(TTree *tree)
{
  // Set object pointer
  cccEt_Et_Eh_FCal = 0;
  cccEt_p_Eh_FCal = 0;
  cccEt_n_Eh_FCal = 0;

  antikt4Track_pt=0;
  antikt4Track_m=0;
  antikt4Track_eta=0;
  antikt4Track_phi=0;
  EMFRClusterJets_pt=0;
  EMFRClusterJets_m=0;
  EMFRClusterJets_eta=0;
  EMFRClusterJets_phi=0;

  antikt2HIItrEMFR_pt=0;
  antikt2HIItrEMFR_m=0;
  antikt2HIItrEMFR_eta=0;
  antikt2HIItrEMFR_phi=0;
  antikt2HIItrEMFR_UncalibEt=0;
  antikt2HIItrEMFR_PassedFR=0;

  antikt4HIItrEMFR_pt = 0;
  antikt4HIItrEMFR_m = 0;
  antikt4HIItrEMFR_eta = 0;
  antikt4HIItrEMFR_phi = 0;
  antikt4HIItrEMFR_UncalibEt = 0;
  antikt4HIItrEMFR_PassedFR = 0;
  antikt2Truth_pt = 0;
  antikt2Truth_m = 0;
  antikt2Truth_eta = 0;
  antikt2Truth_phi = 0;

  mu_muid_E= 0;
  mu_muid_px = 0;
  mu_muid_py = 0;
  mu_muid_pz = 0;
  mu_muid_pt = 0;
  mu_muid_eta = 0;
  mu_muid_phi = 0;

  mc_gen_E= 0; 
  mc_gen_type = 0; 
  mc_gen_mothertype = 0; 
  mc_gen_motherbarcode= 0; 
  mc_gen_eta= 0;
  mc_gen_phi= 0;
  mc_gen_pt= 0;
  mc_gen_barcode= 0;
/*
  trig_L1_TAV = 0;
  trig_L2_passedPhysics = 0;
  trig_EF_passedPhysics = 0;
  trig_L1_TBP = 0;
  trig_L1_TAP = 0;
  trig_L2_passedRaw = 0;
  trig_EF_passedRaw = 0;
  trig_L2_resurrected = 0;
  trig_EF_resurrected = 0;
  trig_L2_passedThrough = 0;
  trig_EF_passedThrough = 0;
  trig_EF_jet_pt = 0;
  trig_EF_jet_m = 0;
  trig_EF_jet_eta = 0;
  trig_EF_jet_phi = 0;
  vx_nTracks = 0;
  vx_sumPt = 0;
*/
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
  fChain->SetBranchAddress("EventNumber", &EventNumber, &b_EventNumber);
  fChain->SetBranchAddress("lbn", &lbn, &b_lbn);
  fChain->SetBranchAddress("mbtime_timeA", &mbtime_timeA, &b_mbtime_timeA);
  fChain->SetBranchAddress("mbtime_timeC", &mbtime_timeC, &b_mbtime_timeC);
  fChain->SetBranchAddress("cccEt_Et_Eh_FCal", &cccEt_Et_Eh_FCal, &b_cccEt_Et_Eh_FCal);
  fChain->SetBranchAddress("cccEt_p_Eh_FCal", &cccEt_p_Eh_FCal, &b_cccEt_p_Eh_FCal);
  fChain->SetBranchAddress("cccEt_n_Eh_FCal", &cccEt_n_Eh_FCal, &b_cccEt_n_Eh_FCal);
  fChain->SetBranchAddress("antikt2HIItrEMFR_pt", &antikt2HIItrEMFR_pt, &b_antikt2HIItrEMFR_pt);
  fChain->SetBranchAddress("antikt2HIItrEMFR_m", &antikt2HIItrEMFR_m, &b_antikt2HIItrEMFR_m);
  fChain->SetBranchAddress("antikt2HIItrEMFR_eta", &antikt2HIItrEMFR_eta, &b_antikt2HIItrEMFR_eta);
  fChain->SetBranchAddress("antikt2HIItrEMFR_phi", &antikt2HIItrEMFR_phi, &b_antikt2HIItrEMFR_phi);
  fChain->SetBranchAddress("antikt2HIItrEMFR_UncalibEt", &antikt2HIItrEMFR_UncalibEt, &b_antikt2HIItrEMFR_UncalibEt);
  fChain->SetBranchAddress("antikt2HIItrEMFR_PassedFR", &antikt2HIItrEMFR_PassedFR, &b_antikt2HIItrEMFR_PassedFR);

  fChain->SetBranchAddress("mu_muid_n",&mu_muid_n,&b_mu_muid_n);
  fChain->SetBranchAddress("mu_muid_E",&mu_muid_E,&b_mu_muid_E);
  fChain->SetBranchAddress("mu_muid_px",&mu_muid_px,&b_mu_muid_px);
  fChain->SetBranchAddress("mu_muid_py",&mu_muid_py,&b_mu_muid_py);
  fChain->SetBranchAddress("mu_muid_pz",&mu_muid_pz,&b_mu_muid_pz);
  fChain->SetBranchAddress("mu_muid_pt",&mu_muid_pt,&b_mu_muid_pt);
  fChain->SetBranchAddress("mu_muid_eta",&mu_muid_eta,&b_mu_muid_eta);
  fChain->SetBranchAddress("mu_muid_phi",&mu_muid_phi,&b_mu_muid_phi);
    if(m_has_truth)
  {
    fChain->SetBranchAddress("antikt2Truth_pt", &antikt2Truth_pt, &b_antikt2Truth_pt);
    fChain->SetBranchAddress("antikt2Truth_m", &antikt2Truth_m, &b_antikt2Truth_m);
    fChain->SetBranchAddress("antikt2Truth_eta", &antikt2Truth_eta, &b_antikt2Truth_eta);
    fChain->SetBranchAddress("antikt2Truth_phi", &antikt2Truth_phi, &b_antikt2Truth_phi);
    fChain->SetBranchAddress("mc_n",&mc_n,&b_mc_n);
    fChain->SetBranchAddress("mc_gen_E",&mc_gen_E,&b_mc_gen_E);
    fChain->SetBranchAddress("mc_gen_type",&mc_gen_type,&b_mc_gen_type);
    fChain->SetBranchAddress("mc_gen_pt",&mc_gen_pt,&b_mc_gen_pt);
    fChain->SetBranchAddress("mc_gen_eta",&mc_gen_eta,&b_mc_gen_eta);
    fChain->SetBranchAddress("mc_gen_phi",&mc_gen_phi,&b_mc_gen_phi);
    fChain->SetBranchAddress("mc_gen_mothertype",&mc_gen_mothertype,&b_mc_gen_mothertype);
    fChain->SetBranchAddress("mc_gen_motherbarcode",&mc_gen_motherbarcode,&b_mc_gen_motherbarcode);
    fChain->SetBranchAddress("mc_gen_barcode",&mc_gen_barcode,&b_mc_gen_barcode); 


  }
 // if(m_has_extra_collections)
 // {

/*
    fChain->SetBranchAddress("antikt4HIItrEMFR_pt", &antikt4HIItrEMFR_pt, &b_antikt4HIItrEMFR_pt);
    fChain->SetBranchAddress("antikt4HIItrEMFR_m", &antikt4HIItrEMFR_m, &b_antikt4HIItrEMFR_m);
    fChain->SetBranchAddress("antikt4HIItrEMFR_eta", &antikt4HIItrEMFR_eta, &b_antikt4HIItrEMFR_eta);
    fChain->SetBranchAddress("antikt4HIItrEMFR_phi", &antikt4HIItrEMFR_phi, &b_antikt4HIItrEMFR_phi);
    fChain->SetBranchAddress("antikt4HIItrEMFR_UncalibEt", &antikt4HIItrEMFR_UncalibEt, &b_antikt4HIItrEMFR_UncalibEt);
    fChain->SetBranchAddress("antikt4HIItrEMFR_PassedFR", &antikt4HIItrEMFR_PassedFR, &b_antikt4HIItrEMFR_PassedFR);
*/

    fChain->SetBranchAddress("antikt4Track_pt", &antikt4Track_pt, &b_antikt4Track_pt);
    fChain->SetBranchAddress("antikt4Track_m", &antikt4Track_m, &b_antikt4Track_m);
    fChain->SetBranchAddress("antikt4Track_eta", &antikt4Track_eta, &b_antikt4Track_eta);
    fChain->SetBranchAddress("antikt4Track_phi", &antikt4Track_phi, &b_antikt4Track_phi);

  //}
/*
  fChain->SetBranchAddress("trig_DB_SMK", &trig_DB_SMK, &b_trig_DB_SMK);
  fChain->SetBranchAddress("trig_DB_L1PSK", &trig_DB_L1PSK, &b_trig_DB_L1PSK);
  fChain->SetBranchAddress("trig_DB_HLTPSK", &trig_DB_HLTPSK, &b_trig_DB_HLTPSK);
  fChain->SetBranchAddress("trig_L1_TAV", &trig_L1_TAV, &b_trig_L1_TAV);
  fChain->SetBranchAddress("trig_L2_passedPhysics", &trig_L2_passedPhysics, &b_trig_L2_passedPhysics);
  fChain->SetBranchAddress("trig_EF_passedPhysics", &trig_EF_passedPhysics, &b_trig_EF_passedPhysics);
  fChain->SetBranchAddress("trig_L1_TBP", &trig_L1_TBP, &b_trig_L1_TBP);
  fChain->SetBranchAddress("trig_L1_TAP", &trig_L1_TAP, &b_trig_L1_TAP);
  fChain->SetBranchAddress("trig_L2_passedRaw", &trig_L2_passedRaw, &b_trig_L2_passedRaw);
  fChain->SetBranchAddress("trig_EF_passedRaw", &trig_EF_passedRaw, &b_trig_EF_passedRaw);
  fChain->SetBranchAddress("trig_L2_resurrected", &trig_L2_resurrected, &b_trig_L2_resurrected);
  fChain->SetBranchAddress("trig_EF_resurrected", &trig_EF_resurrected, &b_trig_EF_resurrected);
  fChain->SetBranchAddress("trig_L2_passedThrough", &trig_L2_passedThrough, &b_trig_L2_passedThrough);
  fChain->SetBranchAddress("trig_EF_passedThrough", &trig_EF_passedThrough, &b_trig_EF_passedThrough);
  fChain->SetBranchAddress("trig_EF_jet_pt", &trig_EF_jet_pt, &b_trig_EF_jet_pt);
  fChain->SetBranchAddress("trig_EF_jet_m", &trig_EF_jet_m, &b_trig_EF_jet_m);
  fChain->SetBranchAddress("trig_EF_jet_eta", &trig_EF_jet_eta, &b_trig_EF_jet_eta);
  fChain->SetBranchAddress("trig_EF_jet_phi", &trig_EF_jet_phi, &b_trig_EF_jet_phi);
  fChain->SetBranchAddress("vx_nTracks", &vx_nTracks, &b_vx_nTracks);
  fChain->SetBranchAddress("vx_sumPt", &vx_sumPt, &b_vx_sumPt);
*/
  fChain->SetBranchAddress("vx_n", &vx_n, &b_vx_n);
}


void JetTreeReader::BuildEventMap(JetTreeReader& reader, map_t& index)
{
  TTree* tr=reader.GetTree();
  //index.clear();
  if (tr == 0) return;
  Long64_t nentries = tr->GetEntriesFast();
  for (long int jentry=0; jentry<nentries;jentry++)
  {
    tr->GetEntry(jentry);

    int statSize=1;
    if(jentry!=0)
    {
      double power=floor(log10(jentry));
      statSize=(int)pow(10.,power);
    }
    if(jentry%statSize==0) std::cout << "Event " << jentry << std::endl;

    index[reader.RunNumber][reader.EventNumber]=jentry;
  }
}




void JetTreeReader::MatchEvents(JetTreeReader& overlay, JetTreeReader& bkgr, const map_t& index, TFile* fout, Long64_t nmax)
{
  TTree* tree_overlay=overlay.GetTree();
  TTree* tree_bkgr=bkgr.GetTree();


  if(fout) fout->cd();
  TTree* out_tree=tree_overlay->CloneTree(0);
  TTree* missing_tree=new TTree("missing_tree","tree of missing events");
  TTree* overlap_tree=new TTree("overlap_tree","tree of overlaping between truth and overlay events");
  overlap_tree = (TTree*)tree_overlay->CloneTree(0);
  missing_tree->Branch("RunNumber",&overlay.RunNumber,"RunNumber/i");
  missing_tree->Branch("EventNumber",&overlay.EventNumber,"EventNumber/i");

  std::vector<bool> Truth_skip;
  std::vector<bool> Truth_Muon_skip;
  out_tree->Branch("Truth_skip",&Truth_skip);
  out_tree->Branch("Truth_Muon_skip",&Truth_Muon_skip);

  std::vector<std::vector<std::vector<float>* > > DR_moments(4,std::vector<std::vector<float>* > ());
  std::vector<std::vector<std::vector<float>* > > DR_moments_FR(4,std::vector<std::vector<float>* > ());

  double DR_bins[7]={0,0.05,0.1,0.2,0.4,0.6,0.8};
  TAxis DR_axis(6,DR_bins);
  
/*
  for(unsigned int iColl=0; iColl<=3; iColl++)
  {
    std::string coll_name=GetCollectionName(iColl);
    for(int i=1; i<=DR_axis.GetNbins(); i++)
    {
      std::vector<float>* temp=new std::vector<float>();
      DR_moments[iColl].push_back(temp);
      
      std::vector<float>* temp_FR=new std::vector<float>();
      DR_moments_FR[iColl].push_back(temp_FR);
      
      std::stringstream ss;
      ss << "antikt4Truth_" << coll_name << "_UEIso" << 100.*DR_axis.GetBinUpEdge(i);
      out_tree->Branch(ss.str().c_str(),"std::vector<float>",&(DR_moments[iColl][i-1]));
      ss << "_FR";
      out_tree->Branch(ss.str().c_str(),"std::vector<float>",&(DR_moments_FR[iColl][i-1]));
    }
  }
*/

  if (tree_overlay == 0) return;
  Long64_t nentries = tree_overlay->GetEntriesFast();
  std::cout << "Tree has " << nentries << " events." << std::endl;
  if(nmax > 0 && nmax < nentries) 
  {
    std::cout << "Processing " << nmax << " events." << std::endl;
    nentries=nmax;
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    tree_overlay->GetEntry(jentry);

    int statSize=1;
    if(jentry!=0)
    {
      double power=floor(log10(jentry));
      statSize=(int)pow(10.,power);
    }
    if(jentry%statSize==0) std::cout << "Event " << jentry << std::endl;

    
    //search for event in map
    map_t::const_iterator run_itr=index.find(overlay.RunNumber);
    if(run_itr==index.end())
    {
      missing_tree->Fill();
      continue;
    }
    std::map<unsigned int, long int>::const_iterator event_itr=run_itr->second.find(overlay.EventNumber);
    if(event_itr==run_itr->second.end())
    {
      missing_tree->Fill();
      continue;
    }


    //now we have the same event in both trees
    tree_bkgr->GetEntry(event_itr->second);    
    
    //initialize moments
/*
    for(unsigned int iColl=0; iColl<=3; iColl++)
    {
      for(int i=1; i<=DR_axis.GetNbins(); i++)
      {
	DR_moments[iColl][i-1]->assign(overlay.antikt4Truth_pt->size(),0);
	DR_moments_FR[iColl][i-1]->assign(overlay.antikt4Truth_pt->size(),0);
      }
    }
*/
    Truth_skip.clear();
    Truth_Muon_skip.clear();
    for(unsigned int iTruth=0; iTruth< overlay.antikt2Truth_pt->size(); iTruth++)
    {
      TLorentzVector truth_jet;
      truth_jet.SetPtEtaPhiM(overlay.antikt2Truth_pt->at(iTruth),
			     overlay.antikt2Truth_eta->at(iTruth),
			     overlay.antikt2Truth_phi->at(iTruth),
			     overlay.antikt2Truth_m->at(iTruth));

      TVector3 truth_jet_momentum(truth_jet.Vect());
      bool truth_skip = false;
/*
      std::cout << std::setw(15) << iTruth
		<< std::setw(15) << truth_jet.Pt()
		<< std::setw(15) << truth_jet.Eta()
		<< std::setw(15) << truth_jet.Phi()
		<< std::endl;
*/

      for(unsigned int iUEJet=0; iUEJet< bkgr.GetNJets(1); iUEJet++)
      {
        TLorentzVector ue_jet;
        bkgr.GetJet(1,iUEJet,ue_jet);
        double DR=truth_jet.DeltaR(ue_jet);
	float ue_jet_pt_cut = 0;
	if(truth_jet.Pt() < 40)
	  ue_jet_pt_cut = 0.5 * truth_jet.Pt();
	else if (truth_jet.Pt() >= 40)
	  ue_jet_pt_cut = 20e3;
	if(DR < 0.8 && ue_jet.Pt() > ue_jet_pt_cut)
	{
	  truth_skip = true;
	  break;
	}
      }//ue jet loop
      Truth_skip.push_back(truth_skip);
    }//truth jet loop

    for(int iMuonTruth = 0; iMuonTruth < overlay.mc_n; iMuonTruth++)
    {
      bool truth_muon_skip = false;
      if(std::abs(overlay.mc_gen_type->at(iMuonTruth)) == 13)
      {
  	TLorentzVector truth_muon;
        truth_muon.SetPtEtaPhiE(overlay.mc_gen_pt->at(iMuonTruth),
			       overlay.mc_gen_eta->at(iMuonTruth),
			       overlay.mc_gen_phi->at(iMuonTruth),
			       overlay.mc_gen_E->at(iMuonTruth));
        TVector3 truth_muon_momentum(truth_muon.Vect());
	for(int iUEMuon = 0; iUEMuon < bkgr.mu_muid_n; iUEMuon++)
	{
	  TLorentzVector ue_muon;
          ue_muon.SetPtEtaPhiE(bkgr.mu_muid_pt->at(iUEMuon),
			       bkgr.mu_muid_eta->at(iUEMuon),
			       bkgr.mu_muid_phi->at(iUEMuon),
			       bkgr.mu_muid_E->at(iUEMuon));
          double DR_muon=truth_muon.DeltaR(ue_muon);
	  if(DR_muon < 0.8 && ue_muon.Pt() > 0.5 * truth_muon.Pt())
	  {
	    truth_muon_skip = true;
	    break;
	  }
	}
      }
      
      Truth_Muon_skip.push_back(truth_muon_skip);
    }
    out_tree->Fill();

  }//event loop
  //save output trees
  if(fout) fout->cd();
  out_tree->AutoSave();
  missing_tree->Write();

}


void JetTreeReader::RunLoop(JetTreeReader& overlay, JetTreeReader& bkgr, TFile* fout, Long64_t nmax)
{

  map_t index;
  std::cout << "Indexing background events.." << std::endl;
  BuildEventMap(bkgr,index);
  std::cout << "Index built, now matching up events..." << std::endl;
  MatchEvents(overlay,bkgr,index,fout,nmax);
}


unsigned int JetTreeReader::GetNJets(unsigned int iColl) const
{
  if(iColl==0) return antikt4HIItrEMFR_m->size();
  else if(iColl==1)return antikt2HIItrEMFR_m->size();
  else if(iColl==2)return  antikt4Track_m->size();
  else if(iColl==3) return EMFRClusterJets_m->size();
  return 0;
}

void JetTreeReader::GetJet(unsigned int iColl, unsigned int ijet, TLorentzVector& vec) const
{
  if(ijet>=GetNJets(iColl)) return;
  if(iColl==0)
  {
    vec.SetPtEtaPhiM(antikt4HIItrEMFR_pt->at(ijet),
		     antikt4HIItrEMFR_eta->at(ijet),
		     antikt4HIItrEMFR_phi->at(ijet),
		     antikt4HIItrEMFR_m->at(ijet));
  }
  else if(iColl==1)
  {
    vec.SetPtEtaPhiM(antikt2HIItrEMFR_pt->at(ijet),
		     antikt2HIItrEMFR_eta->at(ijet),
		     antikt2HIItrEMFR_phi->at(ijet),
		     antikt2HIItrEMFR_m->at(ijet));
  }
  else if(iColl==2)
  {
    vec.SetPtEtaPhiM(antikt4Track_pt->at(ijet),
		     antikt4Track_eta->at(ijet),
		     antikt4Track_phi->at(ijet),
		     antikt4Track_m->at(ijet));
  }

}

/*
float JetTreeReader::PassedFR(unsigned int iColl, unsigned int ijet) const
{
  if(ijet>=GetNJets(iColl)) return 0.;
  if(iColl==0) return antikt4HIItrEMFR_PassedFR->at(ijet);
  else if(iColl==1) return antikt2HIItrEMFR_PassedFR->at(ijet);
  return 1;
}
*/
std::string JetTreeReader::GetCollectionName(unsigned int iColl)
{
  if(iColl==0) return "R4";
  else if(iColl==1)return "R2";
  else if(iColl==2)return  "TJ";
  else if(iColl==3) return "CL";
  return "";
}


