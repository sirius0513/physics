/* vim: set sw=4 ts=4 expandtab: */
#include <iomanip>
#include <iostream>
#include <vector>
#include <TTree.h>
#include <TFile.h>
#include <cmath>
#include <math.h>
#include <TMath.h>
#include <fstream>
#include <sstream>
#include "TRandom.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include <TVector3.h>
#include "AA_Data11_analysis.h"

AA_Data11_analysis::AA_Data11_analysis():m_copy_all(true),m_GRL_reader(0),m_centrality_getter(0),m_time_cut(5.)

{
}

AA_Data11_analysis::AA_Data11_analysis(TTree* inputTree, TTree* inputTree2, bool copy_all):m_copy_all(copy_all),m_GRL_reader(0),m_centrality_getter(0),m_time_cut(5.)
{
    fChain = inputTree;
    fChain2 = inputTree2;
    Init();
}

void AA_Data11_analysis::Init()
{
    if(!m_copy_all)
    {
        fChain->SetBranchStatus("*",0);

        fChain->SetBranchStatus("vx_*",1);

        fChain->SetBranchStatus("trig_EF_jet_n",1);
        fChain->SetBranchStatus("trig_EF_jet_pt",1);
        fChain->SetBranchStatus("trig_EF_jet_eta",1);
        fChain->SetBranchStatus("trig_EF_jet_phi",1);
        fChain->SetBranchStatus("trig_EF_jet_a2hi",1);

        fChain->SetBranchStatus("trig_EF_passedPhysics",1);
        fChain->SetBranchStatus("lbn",1);
        fChain->SetBranchStatus("RunNumber",1);
        fChain->SetBranchStatus("EventNumber",1);
        fChain->SetBranchStatus("vx_n",1);
        fChain->SetBranchStatus("mbtime_timeA",1);
        fChain->SetBranchStatus("mbtime_timeC",1);
        fChain->SetBranchStatus("cccEt_p_Eh_FCal",1);
        fChain->SetBranchStatus("cccEt_n_Eh_FCal",1);
        fChain->SetBranchStatus("cccEt_Et_Eh_FCal",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_n",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_m",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_pt",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_eta",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_phi",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_UncalibEt",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_PassedFR",1);
        fChain->SetBranchStatus("antikt2HIItrEM_n",1);
        fChain->SetBranchStatus("antikt2HIItrEM_m",1);
        fChain->SetBranchStatus("antikt2HIItrEM_pt",1);
        fChain->SetBranchStatus("antikt2HIItrEM_eta",1);
        fChain->SetBranchStatus("antikt2HIItrEM_phi",1);
        fChain->SetBranchStatus("antikt2HIItrEM_UncalibEt",1);
        fChain->SetBranchStatus("mu_muid_n",1);
        fChain->SetBranchStatus("mu_muid_pt",1);
        fChain->SetBranchStatus("mu_muid_px",1);
        fChain->SetBranchStatus("mu_muid_py",1);
        fChain->SetBranchStatus("mu_muid_pz",1);
        fChain->SetBranchStatus("mu_muid_eta",1);
        fChain->SetBranchStatus("mu_muid_phi",1);
        fChain->SetBranchStatus("mu_muid_id_qoverp",1);
        fChain->SetBranchStatus("mu_muid_id_theta",1);
        fChain->SetBranchStatus("mu_muid_me_qoverp",1);
        fChain->SetBranchStatus("mu_muid_me_theta",1);
        fChain->SetBranchStatus("mu_muid_matchchi2",1);
        fChain->SetBranchStatus("mu_muid_matchndof",1);
        fChain->SetBranchStatus("mu_muid_nBLHits",1);
        fChain->SetBranchStatus("mu_muid_nSCTHits",1);
        fChain->SetBranchStatus("mu_muid_nPixHoles",1);
        fChain->SetBranchStatus("mu_muid_nSCTHoles",1);
        fChain->SetBranchStatus("mu_muid_id_d0_exPV",1);
        fChain->SetBranchStatus("mu_muid_id_cov_d0_exPV",1);
        fChain->SetBranchStatus("mu_muid_id_z0_exPV",1);
        fChain->SetBranchStatus("mu_muid_nPixHits",1);
        fChain->SetBranchStatus("muid_scatSignificance",1);
        fChain->SetBranchStatus("mu_muid_isCombinedMuon",1);
/*
        fChain->SetBranchStatus("trk_d0_wrtPV",1);
        fChain->SetBranchStatus("trk_cov_d0_wrtPV",1);
        fChain->SetBranchStatus("trk_eta",1);
        fChain->SetBranchStatus("trk_phi_wrtPV",1);
        fChain->SetBranchStatus("trk_pt",1);
        fChain->SetBranchStatus("trk_chi2",1);
        fChain->SetBranchStatus("trk_ndof",1);
        fChain->SetBranchStatus("trk_patternReco1",1);
        fChain->SetBranchStatus("trk_patternReco2",1);
        fChain->SetBranchStatus("trk_nPixHits",1);
        fChain->SetBranchStatus("trk_nBLHits",1);
        fChain->SetBranchStatus("trk_expectBLayerHit",1);
*/

        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_assoctrk_signOfIP",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_assoctrk_signOfZIP",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_assoctrk_ud0wrtPriVtx",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_assoctrk_ud0ErrwrtPriVtx",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_assoctrk_uz0wrtPriVtx",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_assoctrk_uz0ErrwrtPriVtx",1);


        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_isValid",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_weight_SV0",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_weight_IP3D",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_mass",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_efrc",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_chi2",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_ndof",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_ntrk",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_x",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_y",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_z",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_err_x",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_err_y",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_err_z",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_cov_xy",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_cov_yz",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_sv0p_cov_xz",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_isValid",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_x",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_y",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_z",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_err_x",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_err_y",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_err_z",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_cov_xy",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_cov_yz",1);
        fChain->SetBranchStatus("antikt2HIItrEMFR_flavor_component_svp_cov_xz",1);
        fChain->SetBranchStatus("antikt4Track_n",1);
        fChain->SetBranchStatus("antikt4Track_pt",1);
        fChain->SetBranchStatus("antikt4Track_eta",1);
        fChain->SetBranchStatus("antikt4Track_phi",1);
        fChain->SetBranchStatus("antikt4Track_const_n",1);
        fChain->SetBranchStatus("antikt4Track_const_pt",1);
        fChain->SetBranchStatus("antikt4Track_const_eta",1);
        fChain->SetBranchStatus("antikt4Track_const_phi",1);
    }
    
    trig_EF_jet_pt = 0;
    trig_EF_jet_phi= 0;
    trig_EF_jet_eta= 0;
    trig_EF_jet_a2hi= 0;
    trig_EF_passedPhysics = 0;
    antikt2HIItrEMFR_pt = 0;
    antikt2HIItrEMFR_eta = 0;
    antikt2HIItrEMFR_phi = 0;
    antikt2HIItrEMFR_m= 0;
    antikt2HIItrEMFR_UncalibEt= 0;
    antikt2HIItrEMFR_PassedFR= 0;
    antikt2HIItrEM_pt = 0;
    antikt2HIItrEM_m= 0;
    antikt2HIItrEM_eta = 0;
    antikt2HIItrEM_phi = 0;
    antikt2HIItrEM_UncalibEt= 0;

    antikt4Track_pt = 0;;
    antikt4Track_eta= 0;;
    antikt4Track_phi= 0;;
    antikt4Track_const_pt= 0;;
    antikt4Track_const_eta= 0;;
    antikt4Track_const_phi= 0;;
    mu_muid_pt = 0;
    mu_muid_px = 0;
    mu_muid_py = 0;
    mu_muid_pz = 0;
    mu_muid_eta = 0;
    mu_muid_phi = 0;
 
    mu_muid_nPixHits = 0;
    mu_muid_nBLHits = 0;
    mu_muid_nSCTHits = 0;
    mu_muid_nPixHoles = 0;
    mu_muid_nSCTHoles = 0;
    mu_muid_id_d0_exPV = 0;
    mu_muid_id_cov_d0_exPV= 0;
    mu_muid_id_z0_exPV = 0;
    mu_muid_id_qoverp = 0;
    mu_muid_id_theta= 0;
    mu_muid_me_qoverp = 0;
    mu_muid_me_theta= 0;
    mu_muid_matchchi2 = 0;
    mu_muid_isCombinedMuon= 0;
    mu_muid_matchndof = 0;
    muid_scatSignificance = 0;

    cccEt_p_Eh_FCal = 0;
    cccEt_n_Eh_FCal = 0;
    cccEt_Et_Eh_FCal = 0;
    trk_d0_wrtPV= 0;
    trk_cov_d0_wrtPV= 0;
    trk_eta= 0;
    trk_phi_wrtPV= 0;
    trk_pt= 0;
    trk_chi2= 0;
    trk_ndof= 0;
    trk_patternReco1= 0;
    trk_patternReco2= 0;
    trk_nPixHits= 0;
    trk_nBLHits= 0;
    trk_expectBLayerHit= 0;

	vx_x = 0;
	vx_y = 0;
	vx_z = 0;
	vx_cov_x = 0;
	vx_cov_y = 0;
	vx_cov_z = 0;
	vx_cov_xy = 0;
	vx_cov_yz = 0;
	vx_cov_xz = 0;
    antikt2HIItrEMFR_flavor_component_assoctrk_signOfIP= 0;
    antikt2HIItrEMFR_flavor_component_assoctrk_signOfZIP= 0;
    antikt2HIItrEMFR_flavor_component_assoctrk_ud0wrtPriVtx= 0;
    antikt2HIItrEMFR_flavor_component_assoctrk_ud0ErrwrtPriVtx= 0;
    antikt2HIItrEMFR_flavor_component_assoctrk_uz0wrtPriVtx= 0;
    antikt2HIItrEMFR_flavor_component_assoctrk_uz0ErrwrtPriVtx= 0;
    antikt2HIItrEMFR_flavor_weight_SV0= 0;
    antikt2HIItrEMFR_flavor_weight_IP3D= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_ntrk= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_mass= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_efrc= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_isValid= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_chi2= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_ndof= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_x= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_y= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_z= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_err_x= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_err_y= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_err_z= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_cov_xy= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_cov_yz= 0;
    antikt2HIItrEMFR_flavor_component_sv0p_cov_xz= 0;
    antikt2HIItrEMFR_flavor_component_svp_isValid= 0;
    antikt2HIItrEMFR_flavor_component_svp_x= 0;
    antikt2HIItrEMFR_flavor_component_svp_y= 0;
    antikt2HIItrEMFR_flavor_component_svp_z= 0;
    antikt2HIItrEMFR_flavor_component_svp_err_x= 0;
    antikt2HIItrEMFR_flavor_component_svp_err_y= 0;
    antikt2HIItrEMFR_flavor_component_svp_err_z= 0;
    antikt2HIItrEMFR_flavor_component_svp_cov_xy= 0;
    antikt2HIItrEMFR_flavor_component_svp_cov_yz= 0;
    antikt2HIItrEMFR_flavor_component_svp_cov_xz= 0;

    Muon_eta = 0;
    Muon_phi= 0;
    Muon_pt= 0;
    Muon_composite = 0;

    fChain->SetBranchAddress("trig_EF_jet_n",&trig_EF_jet_n,&b_trig_EF_jet_n); 
    fChain->SetBranchAddress("trig_EF_jet_pt",&trig_EF_jet_pt,&b_trig_EF_jet_pt); 
    fChain->SetBranchAddress("trig_EF_jet_eta",&trig_EF_jet_eta,&b_trig_EF_jet_eta); 
    fChain->SetBranchAddress("trig_EF_jet_phi",&trig_EF_jet_phi,&b_trig_EF_jet_phi); 
    fChain->SetBranchAddress("trig_EF_jet_a2hi",&trig_EF_jet_a2hi,&b_trig_EF_jet_a2hi); 

    fChain->SetBranchAddress("trig_EF_passedPhysics",&trig_EF_passedPhysics,&b_trig_EF_passedPhysics); 
    fChain->SetBranchAddress("lbn",&lbn,&b_lbn); 
    fChain->SetBranchAddress("RunNumber",&RunNumber,&b_RunNumber); 
    fChain->SetBranchAddress("EventNumber",&EventNumber,&b_EventNumber); 
    fChain->SetBranchAddress("vx_n",&vx_n,&b_vx_n); 
    fChain->SetBranchAddress("mbtime_timeA",&mbtime_timeA,&b_mbtime_timeA); 
    fChain->SetBranchAddress("mbtime_timeC",&mbtime_timeC,&b_mbtime_timeC); 
    fChain->SetBranchAddress("cccEt_p_Eh_FCal",&cccEt_p_Eh_FCal,&b_cccEt_p_Eh_FCal); 
    fChain->SetBranchAddress("cccEt_n_Eh_FCal",&cccEt_n_Eh_FCal,&b_cccEt_n_Eh_FCal); 
    fChain->SetBranchAddress("cccEt_Et_Eh_FCal",&cccEt_Et_Eh_FCal,&b_cccEt_Et_Eh_FCal); 
    fChain->SetBranchAddress("antikt2HIItrEMFR_n",&antikt2HIItrEMFR_n,&b_antikt2HIItrEMFR_n);
    fChain->SetBranchAddress("antikt2HIItrEMFR_pt",&antikt2HIItrEMFR_pt,&b_antikt2HIItrEMFR_pt);
    fChain->SetBranchAddress("antikt2HIItrEMFR_eta",&antikt2HIItrEMFR_eta,&b_antikt2HIItrEMFR_eta);
    fChain->SetBranchAddress("antikt2HIItrEMFR_phi",&antikt2HIItrEMFR_phi,&b_antikt2HIItrEMFR_phi);
    fChain->SetBranchAddress("antikt2HIItrEMFR_m",&antikt2HIItrEMFR_m,&b_antikt2HIItrEMFR_m);
    fChain->SetBranchAddress("antikt2HIItrEMFR_UncalibEt",&antikt2HIItrEMFR_UncalibEt,&b_antikt2HIItrEMFR_UncalibEt);
    fChain->SetBranchAddress("antikt2HIItrEMFR_PassedFR",&antikt2HIItrEMFR_PassedFR,&b_antikt2HIItrEMFR_PassedFR);
    fChain->SetBranchAddress("antikt2HIItrEM_n",&antikt2HIItrEM_n,&b_antikt2HIItrEM_n);
    fChain->SetBranchAddress("antikt2HIItrEM_m",&antikt2HIItrEM_m,&b_antikt2HIItrEM_m);
    fChain->SetBranchAddress("antikt2HIItrEM_pt",&antikt2HIItrEM_pt,&b_antikt2HIItrEM_pt);
    fChain->SetBranchAddress("antikt2HIItrEM_eta",&antikt2HIItrEM_eta,&b_antikt2HIItrEM_eta);
    fChain->SetBranchAddress("antikt2HIItrEM_phi",&antikt2HIItrEM_phi,&b_antikt2HIItrEM_phi);
    fChain->SetBranchAddress("antikt2HIItrEM_UncalibEt",&antikt2HIItrEM_UncalibEt,&b_antikt2HIItrEM_UncalibEt);
    fChain->SetBranchAddress("mu_muid_n",&mu_muid_n,&b_mu_muid_n);
    fChain->SetBranchAddress("mu_muid_pt",&mu_muid_pt,&b_mu_muid_pt);
    fChain->SetBranchAddress("mu_muid_px",&mu_muid_px,&b_mu_muid_px);
    fChain->SetBranchAddress("mu_muid_py",&mu_muid_py,&b_mu_muid_py);
    fChain->SetBranchAddress("mu_muid_pz",&mu_muid_pz,&b_mu_muid_pz);
    fChain->SetBranchAddress("mu_muid_eta",&mu_muid_eta,&b_mu_muid_eta);
    fChain->SetBranchAddress("mu_muid_phi",&mu_muid_phi,&b_mu_muid_phi);
    fChain->SetBranchAddress("mu_muid_id_qoverp",&mu_muid_id_qoverp,&b_mu_muid_id_qoverp);
    fChain->SetBranchAddress("mu_muid_me_qoverp",&mu_muid_me_qoverp,&b_mu_muid_me_qoverp);
    fChain->SetBranchAddress("mu_muid_id_theta",&mu_muid_id_theta,&b_mu_muid_id_theta);
    fChain->SetBranchAddress("mu_muid_me_theta",&mu_muid_me_theta,&b_mu_muid_me_theta);
    fChain->SetBranchAddress("mu_muid_matchchi2",&mu_muid_matchchi2,&b_mu_muid_matchchi2);
    fChain->SetBranchAddress("mu_muid_matchndof",&mu_muid_matchndof,&b_mu_muid_matchndof);
    fChain->SetBranchAddress("mu_muid_nBLHits",&mu_muid_nBLHits,&b_mu_muid_nBLHits);
    fChain->SetBranchAddress("mu_muid_nSCTHits",&mu_muid_nSCTHits,&b_mu_muid_nSCTHits);
    fChain->SetBranchAddress("mu_muid_nPixHoles",&mu_muid_nPixHoles,&b_mu_muid_nPixHoles);
    fChain->SetBranchAddress("mu_muid_nSCTHoles",&mu_muid_nSCTHoles,&b_mu_muid_nSCTHoles);
    fChain->SetBranchAddress("mu_muid_id_d0_exPV",&mu_muid_id_d0_exPV,&b_mu_muid_id_d0_exPV);
    fChain->SetBranchAddress("mu_muid_id_cov_d0_exPV",&mu_muid_id_cov_d0_exPV,&b_mu_muid_id_cov_d0_exPV);
    fChain->SetBranchAddress("mu_muid_id_z0_exPV",&mu_muid_id_z0_exPV,&b_mu_muid_id_z0_exPV);
    fChain->SetBranchAddress("mu_muid_nPixHits",&mu_muid_nPixHits,&b_mu_muid_nPixHits);
    fChain->SetBranchAddress("muid_scatSignificance",&muid_scatSignificance,&b_muid_scatSignificance);
    fChain->SetBranchAddress("mu_muid_isCombinedMuon",&mu_muid_isCombinedMuon,&b_mu_muid_isCombinedMuon);
/*
    fChain->SetBranchAddress("trk_d0_wrtPV",&trk_d0_wrtPV,&b_trk_d0_wrtPV);
    fChain->SetBranchAddress("trk_cov_d0_wrtPV",&trk_cov_d0_wrtPV,&b_trk_cov_d0_wrtPV);
    fChain->SetBranchAddress("trk_eta",&trk_eta,&b_trk_eta);
    fChain->SetBranchAddress("trk_phi_wrtPV",&trk_phi_wrtPV,&b_trk_phi);
    fChain->SetBranchAddress("trk_pt",&trk_pt,&b_trk_pt);
    fChain->SetBranchAddress("trk_chi2",&trk_chi2,&b_trk_chi2);
    fChain->SetBranchAddress("trk_ndof",&trk_ndof,&b_trk_ndof);
    fChain->SetBranchAddress("trk_patternReco1",&trk_patternReco1,&b_trk_patternReco1);
    fChain->SetBranchAddress("trk_patternReco2",&trk_patternReco2,&b_trk_patternReco2);
    fChain->SetBranchAddress("trk_nPixHits",&trk_nPixHits,&b_trk_nPixHits);
    fChain->SetBranchAddress("trk_nBLHits",&trk_nBLHits,&b_trk_nBLHits);
    fChain->SetBranchAddress("trk_expectBLayerHit",&trk_expectBLayerHit,&b_trk_expectBLayerHit);
*/

	fChain->SetBranchAddress("vx_x",&vx_x,&b_vx_x);
	fChain->SetBranchAddress("vx_y",&vx_y,&b_vx_y);
	fChain->SetBranchAddress("vx_z",&vx_z,&b_vx_z);
	fChain->SetBranchAddress("vx_cov_x",&vx_cov_x,&b_vx_cov_x);
	fChain->SetBranchAddress("vx_cov_y",&vx_cov_y,&b_vx_cov_y);
	fChain->SetBranchAddress("vx_cov_z",&vx_cov_z,&b_vx_cov_z);
	fChain->SetBranchAddress("vx_cov_xy",&vx_cov_xy,&b_vx_cov_xy);
	fChain->SetBranchAddress("vx_cov_yz",&vx_cov_yz,&b_vx_cov_yz);
	fChain->SetBranchAddress("vx_cov_xz",&vx_cov_xz,&b_vx_cov_xz);

    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_assoctrk_signOfIP",&antikt2HIItrEMFR_flavor_component_assoctrk_signOfIP,&b_assoctrk_d0_sign);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_assoctrk_signOfZIP",&antikt2HIItrEMFR_flavor_component_assoctrk_signOfZIP,&b_assoctrk_z0_sign);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_assoctrk_ud0wrtPriVtx",&antikt2HIItrEMFR_flavor_component_assoctrk_ud0wrtPriVtx,&b_assoctrk_d0);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_assoctrk_ud0ErrwrtPriVtx",&antikt2HIItrEMFR_flavor_component_assoctrk_ud0ErrwrtPriVtx,&b_assoctrk_d0_err);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_assoctrk_uz0wrtPriVtx",&antikt2HIItrEMFR_flavor_component_assoctrk_uz0wrtPriVtx,&b_assoctrk_z0);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_assoctrk_uz0ErrwrtPriVtx",&antikt2HIItrEMFR_flavor_component_assoctrk_uz0ErrwrtPriVtx,&b_assoctrk_z0_err);


    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_isValid", &antikt2HIItrEMFR_flavor_component_sv0p_isValid, &b_antikt2HIItrEMFR_flavor_component_sv0p_isValid);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_weight_SV0", &antikt2HIItrEMFR_flavor_weight_SV0, &b_antikt2HIItrEMFR_flavor_weight_SV0);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_weight_IP3D", &antikt2HIItrEMFR_flavor_weight_IP3D, &b_antikt2HIItrEMFR_flavor_weight_IP3D);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_mass", &antikt2HIItrEMFR_flavor_component_sv0p_mass, &b_antikt2HIItrEMFR_flavor_component_sv0p_mass);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_efrc", &antikt2HIItrEMFR_flavor_component_sv0p_efrc, &b_antikt2HIItrEMFR_flavor_component_sv0p_efrc);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_chi2", &antikt2HIItrEMFR_flavor_component_sv0p_chi2, &b_antikt2HIItrEMFR_flavor_component_sv0p_chi2);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_ndof", &antikt2HIItrEMFR_flavor_component_sv0p_ndof, &b_antikt2HIItrEMFR_flavor_component_sv0p_ndof);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_ntrk", &antikt2HIItrEMFR_flavor_component_sv0p_ntrk, &b_antikt2HIItrEMFR_flavor_component_sv0p_ntrk);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_x", &antikt2HIItrEMFR_flavor_component_sv0p_x, &b_antikt2HIItrEMFR_flavor_component_sv0p_x);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_y", &antikt2HIItrEMFR_flavor_component_sv0p_y, &b_antikt2HIItrEMFR_flavor_component_sv0p_y);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_z", &antikt2HIItrEMFR_flavor_component_sv0p_z, &b_antikt2HIItrEMFR_flavor_component_sv0p_z);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_err_x", &antikt2HIItrEMFR_flavor_component_sv0p_err_x, &b_antikt2HIItrEMFR_flavor_component_sv0p_err_x);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_err_y", &antikt2HIItrEMFR_flavor_component_sv0p_err_y, &b_antikt2HIItrEMFR_flavor_component_sv0p_err_y);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_err_z", &antikt2HIItrEMFR_flavor_component_sv0p_err_z, &b_antikt2HIItrEMFR_flavor_component_sv0p_err_z);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_cov_xy", &antikt2HIItrEMFR_flavor_component_sv0p_cov_xy, &b_antikt2HIItrEMFR_flavor_component_sv0p_cov_xy);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_cov_yz", &antikt2HIItrEMFR_flavor_component_sv0p_cov_yz, &b_antikt2HIItrEMFR_flavor_component_sv0p_cov_yz);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_sv0p_cov_xz", &antikt2HIItrEMFR_flavor_component_sv0p_cov_xz, &b_antikt2HIItrEMFR_flavor_component_sv0p_cov_xz);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_isValid", &antikt2HIItrEMFR_flavor_component_svp_isValid, &b_antikt2HIItrEMFR_flavor_component_svp_isValid);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_x", &antikt2HIItrEMFR_flavor_component_svp_x, &b_antikt2HIItrEMFR_flavor_component_svp_x);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_y", &antikt2HIItrEMFR_flavor_component_svp_y, &b_antikt2HIItrEMFR_flavor_component_svp_y);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_z", &antikt2HIItrEMFR_flavor_component_svp_z, &b_antikt2HIItrEMFR_flavor_component_svp_z);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_err_x", &antikt2HIItrEMFR_flavor_component_svp_err_x, &b_antikt2HIItrEMFR_flavor_component_svp_err_x);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_err_y", &antikt2HIItrEMFR_flavor_component_svp_err_y, &b_antikt2HIItrEMFR_flavor_component_svp_err_y);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_err_z", &antikt2HIItrEMFR_flavor_component_svp_err_z, &b_antikt2HIItrEMFR_flavor_component_svp_err_z);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_cov_xy", &antikt2HIItrEMFR_flavor_component_svp_cov_xy, &b_antikt2HIItrEMFR_flavor_component_svp_cov_xy);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_cov_yz", &antikt2HIItrEMFR_flavor_component_svp_cov_yz, &b_antikt2HIItrEMFR_flavor_component_svp_cov_yz);
    fChain->SetBranchAddress("antikt2HIItrEMFR_flavor_component_svp_cov_xz", &antikt2HIItrEMFR_flavor_component_svp_cov_xz, &b_antikt2HIItrEMFR_flavor_component_svp_cov_xz);

    fChain->SetBranchAddress("antikt4Track_n",&antikt4Track_n,&b_antikt4Track_n);
    fChain->SetBranchAddress("antikt4Track_pt",&antikt4Track_pt,&b_antikt4Track_pt);
    fChain->SetBranchAddress("antikt4Track_eta",&antikt4Track_eta,&b_antikt4Track_eta);
    fChain->SetBranchAddress("antikt4Track_phi",&antikt4Track_phi,&b_antikt4Track_phi);
    fChain->SetBranchAddress("antikt4Track_const_n",&antikt4Track_const_n,&b_antikt4Track_const_n);
    fChain->SetBranchAddress("antikt4Track_const_pt",&antikt4Track_const_pt,&b_antikt4Track_const_pt);
    fChain->SetBranchAddress("antikt4Track_const_eta",&antikt4Track_const_eta,&b_antikt4Track_const_eta);
    fChain->SetBranchAddress("antikt4Track_const_phi",&antikt4Track_const_phi,&b_antikt4Track_const_phi);

    fChain2->SetBranchAddress("Muon_pt",&Muon_pt,&b_Muon_pt);
    fChain2->SetBranchAddress("Muon_eta",&Muon_eta,&b_Muon_eta);
    fChain2->SetBranchAddress("Muon_phi",&Muon_phi,&b_Muon_phi);
    fChain2->SetBranchAddress("Muon_composite",&Muon_composite,&b_Muon_composite);
    fChain2->SetBranchAddress("Muon_IP",&Muon_IP,&b_Muon_IP);
    fChain2->SetBranchAddress("Muon_IP_err",&Muon_IP_err,&b_Muon_IP_err);
    fChain2->SetBranchAddress("Muon_number",&Muon_number, &b_Muon_number);
    fChain2->SetBranchAddress("Centrality",&Centrality, &b_Centrality);
    fChain2->SetBranchAddress("Fcal",&Fcal, &b_Fcal);


}

void AA_Data11_analysis::SetJetTree(TTree* DataTree)
{

    DataTree->Branch("Flavour_tag",&Flavour_tag);
    DataTree->Branch("Centrality_hp",&Centrality_hp);
    DataTree->Branch("Jet_pt",&Jet_pt);
    DataTree->Branch("Jet_eta",&Jet_eta);
    DataTree->Branch("Jet_phi",&Jet_phi);
    DataTree->Branch("Muon_asso_eta",&Muon_asso_eta);
    DataTree->Branch("Muon_asso_pt",&Muon_asso_pt);
    DataTree->Branch("Muon_truth_dR",&Muon_truth_dR);
    DataTree->Branch("Muon_Jet_dR",&Muon_Jet_dR);
    DataTree->Branch("Jet_Trackjet_dR",&Jet_Trackjet_dR);
    DataTree->Branch("Muon_asso_phi",&Muon_asso_phi);
    DataTree->Branch("Muon_asso_composite",&Muon_asso_composite);
    DataTree->Branch("Muon_asso_sig",&Muon_asso_sig);
    DataTree->Branch("Ptrel",&Ptrel);
    DataTree->Branch("Number_Muon_asso",&Number_Muon_asso);
    DataTree->Branch("Jet_weight",&Jet_weight);
    DataTree->Branch("Jet_ip3d",&Jet_ip3d);
    DataTree->Branch("Jet_dl",&Jet_dl);
    DataTree->Branch("Jet_dl_err",&Jet_dl_err);
    DataTree->Branch("Jet_svmass",&Jet_svmass);

}

void AA_Data11_analysis::SetGraph()
{
    std::stringstream ss;

    test_jeta = new TH1F("tjeta","tjeta",50,-5,5);
    bjet_pt = new TH1F("bjpt","bjpt",20,0,100);
    bjet_eta = new TH1F("bjeta","bjeta",50,-5,5);
    hmuon_pt = new TH1F("mpt","mpt",50,0,50);
    hmuon_pt_bjet = new TH1F("mpt_bjet","mpt_bjet",50,0,50);
    hmuon_eta = new TH1F("meta","meta",50,-5,5);
    hmuon_eta_bjet = new TH1F("meta_bjet","meta_bjet",50,-5,5);
    
    hdR_bjet_mpt = new TH2F("dR_bjet_mpt","dR_bjet_mpt",50,0,1,50,0,50);
//    hdR_bjet_mpt_weight = new TH2F("dR_bjet_mpt_weight","dR_bjet_mpt_weight",50,0,1,50,0,50);
    hdR_bjet_pt = new TH2F("dR_bjet_pt","dR_bjet_pt",50,0,1,20,0,100);
    hdR_bjet_ptrel= new TH2F("dR_bjet_ptrel","dR_bjet_ptrel",50,0,1,50,0,1);
    hptrel_bjet_cent = new TH2F("ptrel_bjet_cent","ptrel_bjet_cent",50,0,5,14,-0.5,13.5);
    hjpt_bjet_cent = new TH2F("jpt_bjet_cent","jpt_bjet_cent",20,0,100,10,-0.5,9.5);
    hmpt_bjet_cent = new TH2F("mpt_bjet_cent","mpt_bjet_cent",50,0,50,10,-0.5,9.5);
    hptrel_bjet_pt_cent = new TH3F("ptrel_bjet_pt_cent","ptrel_bjet_pt_cent",50,0,5,20,0,100,10,-0.5,9.5);

    hweight_bjet= new TH1F("weight_bjet","weight_bjet",50,-30,30);
    hmass_bjet_7 = new TH1F("mass_bjet_7","mass_bjet_7",50,0,10);
    hmass_bjet_0 = new TH1F("mass_bjet_0","mass_bjet_0",50,0,10);
    d0 = new TH1F("d0","d0",200,-2,2);
    d0_err = new TH1F("d0_err","d0_err",100,-10,10);
    d0_sig = new TH1F("d0_sig","d0_sig",200,-50,50);
    d0_comp = new TH1F("d0_comp","d0_comp",200,-2,2);
    d0_err_comp = new TH1F("d0_err_comp","d0_err_comp",100,-10,10);
    d0_sig_comp = new TH1F("d0_sig_comp","d0_sig_comp",100,-10,10);

    Jet_jpt_jeta_cent = new TH3F("Jet_jpt_jeta_cent","Jet_jpt_jeta_cent",52,0,260,50,-2.5,2.5,14,-0.5,13.5);
    Jet_jpt_jeta_cent_check = new TH3F("Jet_jpt_jeta_cent_check","Jet_jpt_jeta_cent_check",52,0,260,50,-2.5,2.5,10,-0.5,9.5);
    Fcal_hp = new TH1F("Fcal_hp","Fcal_hp",360,-0.1,3.5);
    Fcal_mb = new TH1F("Fcal_mb","Fcal_mb",360,-0.1,3.5);
    Fcal_mb_weight = new TH1F("Fcal_mb_weight","Fcal_mb_weight",360,-0.1,3.5);

    for(int i =0; i<8; i++)
    {
        ss.str("");
        ss<<"Number_Event_"<<i;
        Number_Event[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),2,0,2,14,-0.5,13.5);

        ss.str("");
        ss<<"Number_Muon_"<<i;
        Number_Muon[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),2,0,2,14,-0.5,13.5);

        ss.str("");
        ss<<"Number_Muon_bkg_"<<i;
        Number_Muon_bkg[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),2,0,2,14,-0.5,13.5);

        ss.str("");
        ss<<"ptrel_bjet_"<<i;
        hptrel_bjet[i] = new TH1F(ss.str().c_str(),ss.str().c_str(),50,0,5);

        ss.str("");
        ss<<"muon_chi2_mpt_meta_"<<i;
        hmuon_chi2_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),500,0,50,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_chi2_mpt_meta_bkg_"<<i;
        hmuon_chi2_mpt_meta_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),500,0,50,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"bjeta_cent_"<<i;
        bjet_eta_cent[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,-5,5,14,-0.5,13.5);

        ss.str("");
        ss<<"bjeta_pt_"<<i;
        bjet_eta_pt[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,-5,5,52,0,260);

        ss.str("");
        ss<<"bjeta_pt_bkg_"<<i;
        bjet_eta_pt_bkg[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,-5,5,52,0,260);

        ss.str("");
        ss<<"meta_bjet_cent_"<<i;
        hmuon_eta_bjet_cent[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,-5,5,14,-0.5,13.5);

        ss.str("");
        ss<<"meta_bjet_cent_bkg_"<<i;
        hmuon_eta_bjet_cent_bkg[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,-5,5,14,-0.5,13.5);

        ss.str("");
        ss<<"ptrel_bjet_pt_"<<i;
        hptrel_bjet_pt[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,0,5,52,0,260);

        ss.str("");
        ss<<"ptrel_bjet_pt_bkg_"<<i;
        hptrel_bjet_pt_bkg[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),50,0,5,52,0,260);

        ss.str("");
        ss<<"muon_composite_bjet_ptrel_"<<i;
        hmuon_composite_bjet_ptrel[i] = new TH2F(ss.str().c_str(),ss.str().c_str(),300,0,3,50,0,5);

        ss.str("");
        ss<<"muon_composite_mpt_meta_"<<i;
        hmuon_composite_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),300,0,3,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_meqoverp_mpt_meta_"<<i;
        hmuon_meqoverp_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),500,-5,5,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_mept_mpt_meta_"<<i;
        hmuon_mept_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),100,0,50,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_mept_mpt_meta_bkg_"<<i;
        hmuon_mept_mpt_meta_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),100,0,50,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_meqoverp_idqoverp_meta_"<<i;
        hmuon_meqoverp_idqoverp_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),100,-50,50,100,-50,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_meqoverp_idqoverp_meta_bkg_"<<i;
        hmuon_meqoverp_idqoverp_meta_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),100,-50,50,100,-50,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_IP_mpt_meta_"<<i;
        hmuon_IP_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),200,-2,2,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_sig_mpt_meta_"<<i;
        hmuon_sig_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),100,-10,10,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"muon_eloss_mpt_meta_"<<i;
        hmuon_eloss_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),200,-5,5,50,0,50,50,-5,5);

        ss.str("");
        ss<<"muon_scat_mpt_meta_"<<i;
        hmuon_scat_mpt_meta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),200,-5,5,50,0,50,50,-5,5);

        ss.str("");
        ss<<"bjeta_phi_cent_"<<i;
        hbjeta_phi_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,-5,5,50,-TMath::Pi(),TMath::Pi(),14,-0.5,13.5);

        ss.str("");
        ss<<"bjeta_phi_pt_"<<i;
        hbjeta_phi_pt[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,-5,5,50,-TMath::Pi(),TMath::Pi(),52,0,260);

        ss.str("");
        ss<<"dR_bjet_ptrel_cent_"<<i;
        hdR_bjet_ptrel_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,5,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_ptrel_cent_bkg_"<<i;
        hdR_bjet_ptrel_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,5,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_mpt_weight_cent_"<<i;
        hdR_bjet_mpt_weight_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_mpt_weight_cent_bkg_"<<i;
        hdR_bjet_mpt_weight_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_mpt_jpt_"<<i;
        hdR_bjet_mpt_jpt[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,52,0,260);

        ss.str("");
        ss<<"dR_bjet_mpt_jpt_bkg_"<<i;
        hdR_bjet_mpt_jpt_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,52,0,260);

        ss.str("");
        ss<<"dR_bjet_mpt_jeta_"<<i;
        hdR_bjet_mpt_jeta[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"dR_bjet_mpt_jeta_bkg_"<<i;
        hdR_bjet_mpt_jeta_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,50,-2.5,2.5);

        ss.str("");
        ss<<"dR_bjet_mpt_cent_"<<i;
        hdR_bjet_mpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_mpt_cent_bkg_"<<i;
        hdR_bjet_mpt_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_pt_cent_"<<i;
        hdR_bjet_pt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,260,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_bjet_pt_cent_bkg_"<<i;
        hdR_bjet_pt_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,52,0,260,14,-0.5,13.5);

        ss.str("");
        ss<<"muon_IP_mpt_cent_"<<i;
        hmuon_IP_mpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),200,-2,2,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"muon_sig_mpt_cent_"<<i;
        hmuon_sig_mpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),100,-10,10,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"muon_err_mpt_cent_"<<i;
        hmuon_err_mpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),200,-50,50,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_mpt_cent_"<<i;
        dR_randm_mpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_mpt_cent_bkg_"<<i;
        dR_randm_mpt_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,0,50,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_jpt_cent_"<<i;
        dR_randm_jpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,52,0,260,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_jpt_cent_bkg_"<<i;
        dR_randm_jpt_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,52,0,260,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_jeta_cent_"<<i;
        dR_randm_jeta_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,-2.5,2.5,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_jeta_cent_bkg_"<<i;
        dR_randm_jeta_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,50,-2.5,2.5,14,-0.5,13.5);

        ss.str("");
        ss<<"ptrel_randm_jpt_cent_"<<i;
        ptrel_randm_jpt_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,5,52,0,260,14,-0.5,13.5);

        ss.str("");
        ss<<"ptrel_randm_jpt_cent_bkg_"<<i;
        ptrel_randm_jpt_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,5,52,0,260,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_num_jet_cent_"<<i;
        dR_randm_num_jet_cent[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,2,0,2,14,-0.5,13.5);

        ss.str("");
        ss<<"dR_randm_num_jet_cent_bkg_"<<i;
        dR_randm_num_jet_cent_bkg[i] = new TH3F(ss.str().c_str(),ss.str().c_str(),50,0,1,2,0,2,14,-0.5,13.5);
    }

}

void AA_Data11_analysis::Analysis(TH1F* h1_weight)
{
    int entries = fChain->GetEntries();
    std::cout<<entries<<std::endl;

    std::vector<Event*> pool[14];

    int  event_check = 0;

    DataJetTree = new TTree("DataJetTree","Tree with jet and associated muon information");
    DataJetTree_em = new TTree("DataJetTree_em","Tree with jet and associated muon information from eventmix");
    SetJetTree(DataJetTree);
    SetJetTree(DataJetTree_em);
    for(int l = 0; l< entries; l++)
    {
        fChain->GetEntry(l);
        if(!PassesEventCuts()) 
            continue;
        float jeta;
        float jphi;
        float jpt;
        float jeta_comp;
        float jphi_comp;
        float jpt_comp;
        float jpx;
        float jpy;
        float jpz;
        float jp;
        float weight;
        float ip3d;
        float svmass;
        float valid;
        float muon_eta;
        float muon_phi;
        float muon_pt;
        float muon_px;
        float muon_py;
        float muon_pz;
        float muon_proj;
        float muon_comp;
        float muon_scat;
        float muon_eloss;
        float ptrel_b;
        bool jetCut;
        bool muonCut;
        bool jet_muon_asso;
        bool jet_trig_asso;
        bool pass_trig;
        float dR_muon ;
        float dR_muon_inclu ;
        float dR_track_jet;
        float dR_jet;
        int particle_ID_muon;
        int particle_ID_track_jet;
        float d_0_unsign;
        float d_0;
        float d_0_err;
        float d_0_sig;
        int number_muon_asso;
        int centrality;
        int Centrality_comb;
        float fcal;
        bool skip_trash;
        bool skip_no_trash;
        bool isoCut;
        bool track_jet_cut;
        bool trig_EF[4];
        int trig_ID;
        float calib_factor_original;
        TLorentzVector Caljet;
        std::vector<int> *muon_check;
        bool muon_skip;

/*
        int ntrk_pv;
        int ntrk_sv;
        int ntrk;
        float sv0p_weight;
        float sv0p_mass;
        float sv0p_efrc;
        float sv0p_chi2;
        int sv0p_ndof;
        int sv0p_valid;
        float sign_d0;
        float sign_z0;
        float dR_trk;
        float trketa;
        float trkphi;
        float trkpt;
        float trkchi2;
        int trkndof;
        int trkpatternReco1;
        int trkpatternReco2;
        int trknPixHits;
        int trknBLHits;
        int trkexpectBLayerHit;;
        float d_0_comp;
        float d_0_err_comp;
        float d_0_sig_comp;

*/
        for(int i =0; i<4; i++)
            trig_EF[i] = false;
        trig_ID = -1;
        skip_trash = false;
        skip_no_trash = false;
        centrality = -1;
        fcal = 0;

        //std::cout<<GetCentrality(GetCentralityVariable())<<"GetCent"<<std::endl;
        //std::cout<<"Event :"<<l<<std::endl; 

        for(int c = 0; c<3; c++)
            fcal+=cccEt_Et_Eh_FCal->at(c)/1000000.0;
/*

        if (fcal> 2.3102) centrality =0; 
        else if (2.3102>fcal&&fcal>1.5707) centrality =1; 
        else if (1.5707>fcal&&fcal>1.0448) centrality =2; 
        else if (1.0448>fcal&&fcal>0.6624) centrality =3; 
        else if (0.6624>fcal&&fcal>0.3909) centrality =4; 
        else if (0.3909>fcal&&fcal>0.2118) centrality =5; 
        else if (0.2118>fcal&&fcal>0.1024) centrality =6; 
        else if (0.1024>fcal&&fcal>0.0438) centrality =7; 
        else if (0.0438>fcal&&fcal>0.0146) centrality =8; 
        else  centrality =9; 
 */      
        centrality = SetCentrality(GetCentrality(GetCentralityVariable()));

        for(unsigned int p = 0; p<trig_EF_passedPhysics->size(); p++)
        {
            if(trig_EF_passedPhysics->at(p)==837|| trig_EF_passedPhysics->at(p)==838)
                skip_trash = true;
            if(trig_EF_passedPhysics->at(p) ==807|| trig_EF_passedPhysics->at(p)== 937||trig_EF_passedPhysics->at(p)==945|| trig_EF_passedPhysics->at(p) == 1003)
                skip_no_trash = true;
            if(trig_EF_passedPhysics->at(p) == 807)
                trig_EF[0]= true;
            if(trig_EF_passedPhysics->at(p) == 945)
                trig_EF[1]= true;
            if(trig_EF_passedPhysics->at(p) == 937)
                trig_EF[2]= true;
            if(trig_EF_passedPhysics->at(p) == 1003)
                trig_EF[3]= true;
                
        }


        if(skip_trash && !skip_no_trash)
            continue;
        else
            for(int i = 0; i<2; i++)
                if(trig_EF[i])
                    Number_Event[i]->Fill(1,centrality);
////////////////////////////////////////////muon loop //////////////////////


        Muon this_muon;
        Event* this_event = new Event();
        this_event->RandMuon.clear() ;
        int depth[14]; 
        for(int d = 0; d< 14; d++)
            depth[d] = 10;

        this_event->cent = centrality; 

        for(int m = 0; m < mu_muid_n; m++)
        {

            muon_eta = mu_muid_eta->at(m);
            muon_phi= mu_muid_phi->at(m);
            muon_pt= mu_muid_pt->at(m)/1000.0;
            muon_scat = muid_scatSignificance->at(m);
            muon_eloss = 1.0 - std::abs(mu_muid_id_qoverp->at(m)/mu_muid_me_qoverp->at(m));

            muonCut = false;
           
            if(std::abs(mu_muid_id_qoverp->at(m)) > pow(10.0,-6) && mu_muid_matchndof->at(m) >=1)
                muonCut = SelectMuon(m);
            else if (std::abs(mu_muid_id_qoverp->at(m)) < pow(10.0,-6))
                std::cout<<mu_muid_id_qoverp->at(m)<<"qoverp"<<std::endl;

            if(muonCut)
            {
                this_muon.eta = mu_muid_eta->at(m);
                this_muon.phi= mu_muid_phi->at(m);
                this_muon.pt= mu_muid_pt->at(m);
                this_muon.composite = std::abs(muon_eloss) + 0.07*std::abs(muon_scat);

                if(trig_EF[0])
                {
                    if(this_muon.composite < 0.25)
                        Number_Muon[0]->Fill(1,centrality);
                    else
                        Number_Muon_bkg[0]->Fill(1,centrality);

                    this_event->RandMuon.push_back(this_muon);
                }
            }
        }


        for(int e = 0; e< antikt2HIItrEM_n; e++)
        {
            Fcal_hp->Fill(fcal);
            Caljet.SetPtEtaPhiM(antikt2HIItrEM_pt->at(e),antikt2HIItrEM_eta->at(e),antikt2HIItrEM_phi->at(e),antikt2HIItrEM_m->at(e));
            calib_factor_original=Caljet.Et()/antikt2HIItrEM_UncalibEt->at(e);
            Caljet*=1.0/calib_factor_original;
            Caljet*=1e-3; // convert from MeV to GeV
            m_calibrator->Calibrate(Caljet);

            jpt = Caljet.Pt();
            jeta = antikt2HIItrEM_eta->at(e);
            jphi= antikt2HIItrEM_phi->at(e);

            dR_muon = 100.0;
            particle_ID_muon= -1.0;
            jetCut = false;
            isoCut = false;
            jet_trig_asso = false;
            jetCut = std::abs(jeta) < 2.5 && jpt > 20.0 ;
            number_muon_asso = 0;

/////////////////////////////////////////////////////////jet isolation cut//////////////////////////////////////////////////

            for(int o = 0; o<antikt2HIItrEM_n; o++)
            {
                //std::cout<<"Event: "<<l<<"\t Jet :"<<e<<"\t Jet comp :"<<o<<std::endl;
                Caljet.SetPtEtaPhiM(antikt2HIItrEM_pt->at(o),antikt2HIItrEM_eta->at(o),antikt2HIItrEM_phi->at(o),antikt2HIItrEM_m->at(o));
                calib_factor_original=Caljet.Et()/antikt2HIItrEM_UncalibEt->at(o);
                Caljet*=1.0/calib_factor_original;
                Caljet*=1e-3; // convert from MeV to GeV
                m_calibrator->Calibrate(Caljet);

                jpt_comp = Caljet.Pt();
                jeta_comp = antikt2HIItrEM_eta->at(o);
                jphi_comp = antikt2HIItrEM_phi->at(o);
//                if( o != e && Association(jphi, jphi_comp, jeta, jeta_comp) < 0.8 && antikt2HIItrEM_PassedFR->at(o) > 0.99)
                if(Association(jphi, jphi_comp, jeta, jeta_comp) < 0.8 && jpt_comp > jpt)
                {
                    isoCut = true;
                    break;
                }
            
            }
            if(isoCut)
                continue;

///////////////////////////////////////////////////////Track jet correction//////////////////////////////////////////////
/*
            for(int c = 0; c < antikt4Track_n; c++)
            {
                bool track_jet_cut = antikt4Track_const_n->at(c) > 1 ;
                if(dR_track_jet > Association(jphi, antikt4Track_phi->at(c), jeta, antikt4Track_eta->at(c)) && track_jet_cut)
                {
                    dR_track_jet = Association(jphi, antikt4Track_phi->at(c), jeta, antikt4Track_eta->at(c));
                    particle_ID_track_jet = c;
                }
            }
*/

/////////////////////////////////////////////////////////jet fire jet trigger//////////////////////////////////////////////////
            for(int t = 0; t< trig_EF_jet_n; t++)
            {
                if(trig_EF_jet_pt->at(t)/1000.0 > 20.0 && trig_EF_jet_a2hi->at(t) == 1)
                {
                    if(Association(antikt2HIItrEM_phi->at(e), trig_EF_jet_phi->at(t), antikt2HIItrEM_eta->at(e), trig_EF_jet_eta->at(t)) < 0.4)
                    {
                        jet_trig_asso = true;
                        break;
                    }
                }
            }


            if(trig_EF[0] && jet_trig_asso && jetCut)
            {
                Jet_jpt_jeta_cent->Fill(jpt,jeta,centrality);

                for(int m = 0; m < mu_muid_n; m++)
                {
                    muon_eta = mu_muid_eta->at(m);
                    muon_phi= mu_muid_phi->at(m);
                    muon_pt= mu_muid_pt->at(m)/1000.0;
                    muon_scat = muid_scatSignificance->at(m);
                    muon_eloss = 1.0 - std::abs(mu_muid_id_qoverp->at(m)/mu_muid_me_qoverp->at(m));
                    muon_comp = std::abs(muon_eloss)+0.07*std::abs(muon_scat);
                    muonCut = false;
                    muon_skip = false;
                    dR_track_jet= 100.0;
                    particle_ID_track_jet= -1.0;
                    d_0_unsign = mu_muid_id_d0_exPV->at(m);

                    if((d_0_unsign > 0.0 && DeltaPhi(antikt2HIItrEM_phi->at(e), muon_phi) > 0.0) || (d_0_unsign < 0.0 && DeltaPhi(antikt2HIItrEM_phi->at(e), muon_phi) < 0.0))
                        d_0 = std::abs(d_0_unsign);
                    
                    else
                        d_0 = -std::abs(d_0_unsign);

                    d_0_err = sqrt(mu_muid_id_cov_d0_exPV->at(m));
                    d_0_sig = d_0/d_0_err;

                  //  dR_muon_inclu = Association(jphi, muon_phi, jeta, muon_eta);

                    if(std::abs(mu_muid_id_qoverp->at(m)) > pow(10.0,-6) && mu_muid_matchndof->at(m) >=1) 
                        muonCut = SelectMuon(m);

                    if(muonCut)
                    {
                        bool muonDynamiCut = (jpt >60 && jpt < 90 && muon_pt >6) || (jpt >90 && muon_pt >8) || (jpt < 60 && jpt >20);
                        muon_px = muon_pt*cos(muon_phi);
                        muon_py = muon_pt*sin(muon_phi);
                        muon_pz = muon_pt*sinh(muon_eta);
                        TVector3 Muon_reco(muon_px,muon_py,muon_pz);

                        float track_jet_near_eta,track_jet_near_phi;

                        for(int c = 0; c < antikt4Track_n; c++) 
                        {
                            bool have_muon_intrack = false;
                            float dR_track_muon = 100.0;
                            int particle_ID_track_muon= -1.0;
                            bool track_pt_cut;
                            float track_jet_eta, track_jet_phi;
                            for(int t = 0; t < antikt4Track_const_n->at(c); t++)
                            {
                                track_pt_cut = (antikt4Track_const_pt->at(c)[t]/1000.0-muon_pt)/muon_pt< 0.5;
                                if(dR_track_muon> Association(muon_phi, antikt4Track_const_phi->at(c)[t], muon_eta, antikt4Track_const_eta->at(c)[t]) && track_pt_cut)
                                {
                                    dR_track_muon= Association(muon_phi, antikt4Track_const_phi->at(c)[t], muon_eta, antikt4Track_const_eta->at(c)[t]);
                                    particle_ID_track_muon= t;

                                }
                            }
                            if(dR_track_muon < 0.02)
                            {

//                                std::cout<<"Event: "<<l<<"\t Jet :"<<e<<"\t Muon:"<<m<<std::endl;
//                                std::cout<<"muon phi:"<<muon_phi<<"\t Trackmuon phi:"<<antikt4Track_const_phi->at(c)[particle_ID_track_muon]<<std::endl;
 //                               std::cout<<"muon eta:"<<muon_eta<<"\t Trackmuon eta:"<<antikt4Track_const_eta->at(c)[particle_ID_track_muon]<<std::endl;
  //                              std::cout<<"muon pt:"<<muon_pt<<"\t Trackmuon eta:"<<antikt4Track_const_pt->at(c)[particle_ID_track_muon]/1000.0<<std::endl;
                
                                have_muon_intrack = true;
                                float trackjet_px = antikt4Track_pt->at(c)/1000.0*cos(antikt4Track_phi->at(c));
                                float trackjet_py = antikt4Track_pt->at(c)/1000.0*sin(antikt4Track_phi->at(c));
                                float trackjet_pz = antikt4Track_pt->at(c)/1000.0*sinh(antikt4Track_eta->at(c));
                                TVector3 Trackjet(trackjet_px,trackjet_py,trackjet_pz);
                                Trackjet -= Muon_reco;
                                track_jet_eta = Trackjet.Eta();
                                track_jet_phi= Trackjet.Phi();
                            }
                            else
                            {
                                track_jet_eta = antikt4Track_eta->at(c);
                                track_jet_phi= antikt4Track_phi->at(c);
                            }

                            track_jet_cut = antikt4Track_const_n->at(c) > 1 ;

                            jphi = antikt2HIItrEM_phi->at(e);
                            jeta= antikt2HIItrEM_eta->at(e);

                            if(dR_track_jet > Association(jphi, track_jet_phi, jeta, track_jet_eta) && track_jet_cut )
                            {
                                dR_track_jet = Association(jphi, track_jet_phi, jeta, track_jet_eta);
                                particle_ID_track_jet = c;
                                track_jet_near_eta = track_jet_eta;
                                track_jet_near_phi= track_jet_phi;
                            }
                        }

                        if(dR_track_jet < 0.2 && particle_ID_track_jet > -1.0)
                        {
                           // jphi= antikt4Track_phi->at(particle_ID_track_jet);
                           // jeta= antikt4Track_eta->at(particle_ID_track_jet);
/*
                            std::cout<<"Jet phi:"<<jphi<<"\t Track phi:"<<track_jet_near_phi<<"\t Track orig phi:"<<antikt4Track_phi->at(particle_ID_track_jet)<<std::endl;
                            std::cout<<"Jet eta:"<<jeta<<"\t Track eta:"<<track_jet_near_eta<<"\t Track orig eta:"<<antikt4Track_eta->at(particle_ID_track_jet)<<std::endl;
*/
                            jphi = track_jet_near_phi;
                            jeta= track_jet_near_eta;
                        }

                        jpx = jpt*cos(jphi) + muon_px;
                        jpy = jpt*sin(jphi) + muon_py;
                        jpz = jpt*sinh(jeta) + muon_pz ;
                        jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                        muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;

                        TVector3 Jet_reco(jpx,jpy,jpz);
                        TVector3 Cross_Jet_Muon_reco;
                        Cross_Jet_Muon_reco = Muon_reco.Cross(Jet_reco);

                        dR_muon_inclu= Association(muon_phi,Jet_reco.Phi(),muon_eta,Jet_reco.PseudoRapidity());
                        ptrel_b= muon_pt * sin(dR_muon_inclu);
                        //ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

                        if(dR_muon_inclu < 0.6)
                        {
                            number_muon_asso ++;
                            if(muon_comp < 0.25)
                            {
                                bjet_eta_pt[1]->Fill(jeta,jpt);
                                hmuon_eta_bjet_cent[1]->Fill(muon_eta,centrality);
                                hptrel_bjet_pt[1]->Fill(ptrel_b,jpt);
                                hdR_bjet_mpt_cent[1]->Fill(dR_muon_inclu,muon_pt,centrality);
                                hdR_bjet_mpt_jpt[1]->Fill(dR_muon_inclu,muon_pt,jpt);
                                hdR_bjet_mpt_jeta[1]->Fill(dR_muon_inclu,muon_pt,jeta);
                                hdR_bjet_ptrel_cent[1]->Fill(dR_muon_inclu,ptrel_b,centrality);

                                if(std::abs(muon_eta)< 1.05)
                                {
                                    bjet_eta_pt[7]->Fill(jeta,jpt);
                                    hmuon_eta_bjet_cent[7]->Fill(muon_eta,centrality);
                                    hdR_bjet_mpt_cent[7]->Fill(dR_muon_inclu,muon_pt,centrality);
                                    hdR_bjet_mpt_jpt[7]->Fill(dR_muon_inclu,muon_pt,jpt);
                                    hptrel_bjet_pt[7]->Fill(ptrel_b,jpt);
                                    hdR_bjet_ptrel_cent[7]->Fill(dR_muon_inclu,ptrel_b,centrality);
                                }
                                
                            }
                            else
                            {
                                bjet_eta_pt_bkg[1]->Fill(jeta,jpt);
                                hmuon_eta_bjet_cent_bkg[1]->Fill(muon_eta,centrality);
                                hptrel_bjet_pt_bkg[1]->Fill(ptrel_b,jpt);
                                hdR_bjet_mpt_cent_bkg[1]->Fill(dR_muon_inclu,muon_pt,centrality);
                                hdR_bjet_mpt_jpt_bkg[1]->Fill(dR_muon_inclu,muon_pt,jpt);
                                hdR_bjet_mpt_jeta_bkg[1]->Fill(dR_muon_inclu,muon_pt,jeta);
                                hdR_bjet_ptrel_cent_bkg[1]->Fill(dR_muon_inclu,ptrel_b,centrality);

                                if(std::abs(muon_eta)<1.05) 
                                {
                                    bjet_eta_pt_bkg[7]->Fill(jeta,jpt);
                                    hmuon_eta_bjet_cent_bkg[7]->Fill(muon_eta,centrality);
                                    hdR_bjet_mpt_cent_bkg[7]->Fill(dR_muon_inclu,muon_pt,centrality);
                                    hdR_bjet_mpt_jpt_bkg[7]->Fill(dR_muon_inclu,muon_pt,jpt);
                                    hptrel_bjet_pt_bkg[7]->Fill(ptrel_b,jpt);
                                    hdR_bjet_ptrel_cent_bkg[7]->Fill(dR_muon_inclu,ptrel_b,centrality);
                                }
                            }
                            if(d_0_err<5.0&&std::abs(d_0)<4.0)
                            {

                                if(d_0_sig >= 2.0)
                                {
                                    bjet_eta_pt[3]->Fill(jeta,jpt);
                                    hmuon_eta_bjet_cent[3]->Fill(muon_eta,centrality);
                                    hptrel_bjet_pt[3]->Fill(ptrel_b,jpt);
                                    hdR_bjet_mpt_cent[3]->Fill(dR_muon_inclu,muon_pt,centrality);
                                    hdR_bjet_mpt_jpt[3]->Fill(dR_muon_inclu,muon_pt,jpt);
                                    hdR_bjet_pt_cent[3]->Fill(dR_muon_inclu,jpt,centrality);
                                    hdR_bjet_ptrel_cent[3]->Fill(dR_muon_inclu,ptrel_b,centrality);
                                }
                                else
                                {
                                    bjet_eta_pt_bkg[3]->Fill(jeta,jpt);
                                    hmuon_eta_bjet_cent_bkg[3]->Fill(muon_eta,centrality);
                                    hptrel_bjet_pt_bkg[3]->Fill(ptrel_b,jpt);
                                    hdR_bjet_ptrel_cent_bkg[3]->Fill(dR_muon_inclu,ptrel_b,centrality);
                                    hdR_bjet_mpt_cent_bkg[3]->Fill(dR_muon_inclu,muon_pt,centrality);
                                    hdR_bjet_mpt_jpt_bkg[3]->Fill(dR_muon_inclu,muon_pt,jpt);
                                    hdR_bjet_pt_cent_bkg[3]->Fill(dR_muon_inclu,jpt,centrality);
                                }
                                if(muon_comp < 0.25 && d_0_sig >=2.0)
                                {
                                    bjet_eta_pt[5]->Fill(jeta,jpt);
                                    hmuon_eta_bjet_cent[5]->Fill(muon_eta,centrality);
                                    hptrel_bjet_pt[5]->Fill(ptrel_b,jpt);
                                    hdR_bjet_mpt_cent[5]->Fill(dR_muon_inclu,muon_pt,centrality);
                                    hdR_bjet_mpt_jpt[5]->Fill(dR_muon_inclu,muon_pt,jpt);
                                    hdR_bjet_pt_cent[5]->Fill(dR_muon_inclu,jpt,centrality);
                                    hdR_bjet_ptrel_cent[5]->Fill(dR_muon_inclu,ptrel_b,centrality);
                                }
                                else
                                {
                                    bjet_eta_pt_bkg[5]->Fill(jeta,jpt);
                                    hmuon_eta_bjet_cent_bkg[5]->Fill(muon_eta,centrality);
                                    hptrel_bjet_pt_bkg[5]->Fill(ptrel_b,jpt);
                                    hdR_bjet_mpt_cent_bkg[5]->Fill(dR_muon_inclu,muon_pt,centrality);
                                    hdR_bjet_mpt_jpt_bkg[5]->Fill(dR_muon_inclu,muon_pt,jpt);
                                    hdR_bjet_pt_cent_bkg[5]->Fill(dR_muon_inclu,jpt,centrality);
                                    hdR_bjet_ptrel_cent_bkg[5]->Fill(dR_muon_inclu,ptrel_b,centrality);
                                }
                            }

/////////////////////////////////////////Fill JetTree///////////////////////////////////////////////////////////////

                            Muon_asso_pt = muon_pt;
                            Muon_asso_eta = muon_eta;
                            Muon_asso_phi= muon_phi;
                            Muon_asso_composite= muon_comp;
                            Muon_asso_sig= d_0_sig;
                            Ptrel= ptrel_b;
                            Flavour_tag = 0;
                            Centrality_hp= centrality;
                            Jet_pt = jpt;
                            Jet_eta= jeta;
                            Jet_phi= jphi;
                            Muon_Jet_dR = dR_muon_inclu;
                            Jet_Trackjet_dR = dR_track_jet;
                            Number_Muon_asso = number_muon_asso;

                            float dR_EM_FR = 100.0;
                            int   particle_ID_EM_FR = -1.0;
                            for(int fr = 0; fr < antikt2HIItrEMFR_n; fr++)
                                if(Association(antikt2HIItrEMFR_phi->at(fr),jphi, antikt2HIItrEMFR_eta->at(fr), jeta) < dR_EM_FR)
                                {
                                    dR_EM_FR = Association(antikt2HIItrEMFR_phi->at(fr),jphi, antikt2HIItrEMFR_eta->at(fr), jeta);
                                    particle_ID_EM_FR = fr;
                                }
                                
                            if(particle_ID_EM_FR > -1.0 && dR_EM_FR < 0.1)
                            {
                                weight = antikt2HIItrEMFR_flavor_weight_SV0->at(particle_ID_EM_FR);
                                ip3d= antikt2HIItrEMFR_flavor_weight_IP3D->at(particle_ID_EM_FR);
                                valid = antikt2HIItrEMFR_flavor_component_sv0p_isValid->at(particle_ID_EM_FR);
                                svmass= antikt2HIItrEMFR_flavor_component_sv0p_mass->at(particle_ID_EM_FR);
                            }
                            else
                            {
                                weight = 0.0;
                                ip3d= 0.0;
                                svmass = 0.0;
                                valid = 0.0;
                            }
                            Jet_weight = weight;
                            Jet_ip3d= ip3d;
                            Jet_svmass = svmass;
                            if(valid !=0.0)
                            {
                                Jet_dl_err = Calerror(particle_ID_EM_FR);
                                Jet_dl = weight * Jet_dl_err;
                            }

                            DataJetTree->Fill();

                        }
                        if(dR_muon > Association(jphi, muon_phi, jeta, muon_eta))
                        {   
                            dR_muon = Association(jphi, muon_phi, jeta, muon_eta);
                            particle_ID_muon= m;
                        }   
                    }
                }

                /////////////////////////////////////////Nearest Muon result//////////////////////////////////////////////////////

                if(dR_muon < 0.6 && particle_ID_muon >=0)
                {
                    //std::cout<<"muon dR : "<<dR_muon<<"\t"<<"Event number: "<<l<<"\t"<<"Jet number: "<<e<<std::endl;
                    muon_eta = mu_muid_eta->at(particle_ID_muon);
                    muon_phi= mu_muid_phi->at(particle_ID_muon);
                    muon_pt= mu_muid_pt->at(particle_ID_muon)/1000.0;
                    muon_scat = muid_scatSignificance->at(particle_ID_muon);
                    muon_eloss = 1.0 - std::abs(mu_muid_id_qoverp->at(particle_ID_muon)/mu_muid_me_qoverp->at(particle_ID_muon));
                    muon_comp = std::abs(muon_eloss)+0.07*std::abs(muon_scat);

                    d_0_unsign = mu_muid_id_d0_exPV->at(particle_ID_muon);

                    if((d_0_unsign > 0.0 && DeltaPhi(antikt2HIItrEM_phi->at(e), muon_phi) > 0.0) || (d_0_unsign < 0.0 && DeltaPhi(antikt2HIItrEM_phi->at(e), muon_phi) < 0.0))
                        d_0 = std::abs(d_0_unsign);
                    
                    else
                        d_0 = -std::abs(d_0_unsign);

                    d_0_err = sqrt(mu_muid_id_cov_d0_exPV->at(particle_ID_muon));
                    d_0_sig = d_0/d_0_err;


                    jpx = jpt*cos(jphi);
                    jpy = jpt*sin(jphi);
                    jpz = jpt*sinh(jeta);
                    muon_px = muon_pt*cos(muon_phi);
                    muon_py = muon_pt*sin(muon_phi);
                    muon_pz = muon_pt*sinh(muon_eta);
                    jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                    muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;
                    ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

                
/////////////////////////////////////////////track jet//////////////////////////////////////////

                    for(int c = 0; c < antikt4Track_n; c++)
                    {
                        bool track_jet_cut = antikt4Track_pt->at(c)/1000.0 > 20.0 && antikt4Track_const_n->at(c) > 1 && std::abs(antikt4Track_eta->at(c) < 2.5);
                        if(dR_track_jet > Association(muon_phi, antikt4Track_phi->at(c), muon_eta, antikt4Track_eta->at(c)) && track_jet_cut)
                        {
                            dR_track_jet = Association(muon_phi, antikt4Track_phi->at(c), muon_eta, antikt4Track_eta->at(c));
                            particle_ID_track_jet = c;
                        }
                    }
                    if(dR_track_jet < 0.6 && particle_ID_track_jet > 0.0)
                    {
                            std::cout<<"muon_eta: "<<muon_eta<<"\t"<<"muon_phi: "<<muon_phi<<"\t"<<"muon_pt: "<<muon_pt<<std::endl;
                            std::cout<<"jet_eta: "<<jeta<<"\t"<<"jet_phi: "<<jphi<<"\t"<<"jet_pt: "<<jpt<<std::endl;
                            for(int t = 0; t< antikt4Track_const_n->at(particle_ID_track_jet); t++)
                            {
                                std::cout<<"particle_Id: "<<particle_ID_track_jet<<std::endl;
                                std::cout<<"const_n: "<<antikt4Track_const_n->at(particle_ID_track_jet)<<std::endl;
                                std::cout<<"Number: "<<t<<std::endl;
                                std::cout<<"track_jet_pt: "<<antikt4Track_const_pt->at(particle_ID_track_jet)[t]/1000.0<<"\t"<<"track_jet_phi: "<<antikt4Track_const_phi->at(particle_ID_track_jet)[t]<<"\t"<<"track_jet_phi: "<<antikt4Track_const_phi->at(particle_ID_track_jet)[t]<<"\t dR: "<<Association(muon_phi,antikt4Track_const_phi->at(particle_ID_track_jet)[t],muon_eta,antikt4Track_const_eta->at(particle_ID_track_jet)[t])<<std::endl;
                            }
                    }



                    hptrel_bjet[0]->Fill(ptrel_b);
                    bjet_eta_cent[0]->Fill(jeta,centrality);

                    hmuon_eloss_mpt_meta[0]->Fill(muon_eloss,muon_pt,muon_eta);
                    hmuon_scat_mpt_meta[0]->Fill(muon_scat,muon_pt,muon_eta);
                    hbjeta_phi_cent[0]->Fill(jeta,jphi,centrality);
                    hbjeta_phi_pt[0]->Fill(jeta,jphi,jpt);

                    hmuon_composite_bjet_ptrel[0]->Fill(muon_comp,ptrel_b);
                    hmuon_composite_mpt_meta[0]->Fill(muon_comp,muon_pt,muon_eta);
                    hmuon_meqoverp_mpt_meta[0]->Fill(mu_muid_me_qoverp->at(particle_ID_muon)*1000,muon_pt,muon_eta);
        
                    if(muon_comp < 0.25)
                    {
                        bjet_eta_pt[0]->Fill(jeta,jpt);
                        hmuon_eta_bjet_cent[0]->Fill(muon_eta,centrality);
                        hptrel_bjet_pt[0]->Fill(ptrel_b,jpt);
                        hdR_bjet_mpt_cent[0]->Fill(dR_muon,muon_pt,centrality);
                        hdR_bjet_mpt_jpt[0]->Fill(dR_muon,muon_pt,jpt);
                        hdR_bjet_mpt_jeta[0]->Fill(dR_muon,muon_pt,jeta);
                        hdR_bjet_pt_cent[0]->Fill(dR_muon,jpt,centrality);
                        hdR_bjet_ptrel_cent[0]->Fill(dR_muon,ptrel_b,centrality);
                        hmuon_meqoverp_idqoverp_meta[0]->Fill(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000,sin(mu_muid_id_theta->at(particle_ID_muon))/mu_muid_id_qoverp->at(particle_ID_muon)/1000,muon_eta);
                        hmuon_mept_mpt_meta[0]->Fill(std::abs(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000),muon_pt,muon_eta);
                        hmuon_chi2_mpt_meta[0]->Fill(mu_muid_matchchi2->at(particle_ID_muon)/mu_muid_matchndof->at(particle_ID_muon),muon_pt,muon_eta);

                        if(dR_muon > pow(10.0,-6))
                            hdR_bjet_mpt_weight_cent[0]->Fill(dR_muon, muon_pt,centrality,1/dR_muon); 
                        if(std::abs(muon_eta)< 1.05)
                        {
                            bjet_eta_pt[6]->Fill(jeta,jpt);
                            hmuon_eta_bjet_cent[6]->Fill(muon_eta,centrality);
                            hdR_bjet_mpt_cent[6]->Fill(dR_muon,muon_pt,centrality);
                            hdR_bjet_mpt_jpt[6]->Fill(dR_muon,muon_pt,jpt);
                            hdR_bjet_pt_cent[6]->Fill(dR_muon,jpt,centrality);
                            hptrel_bjet_pt[6]->Fill(ptrel_b,jpt);
                            hdR_bjet_ptrel_cent[6]->Fill(dR_muon,ptrel_b,centrality);
                            hmuon_meqoverp_idqoverp_meta[6]->Fill(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000,sin(mu_muid_id_theta->at(particle_ID_muon))/mu_muid_id_qoverp->at(particle_ID_muon)/1000,muon_eta);
                            hmuon_mept_mpt_meta[6]->Fill(std::abs(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000),muon_pt,muon_eta);
                            hmuon_chi2_mpt_meta[6]->Fill(mu_muid_matchchi2->at(particle_ID_muon)/mu_muid_matchndof->at(particle_ID_muon),muon_pt,muon_eta);
                        }
                        
                    }
                    else
                    {
                        bjet_eta_pt_bkg[0]->Fill(jeta,jpt);
                        hmuon_eta_bjet_cent_bkg[0]->Fill(muon_eta,centrality);
                        hptrel_bjet_pt_bkg[0]->Fill(ptrel_b,jpt);
                        hdR_bjet_mpt_cent_bkg[0]->Fill(dR_muon,muon_pt,centrality);
                        hdR_bjet_mpt_jpt_bkg[0]->Fill(dR_muon,muon_pt,jpt);
                        hdR_bjet_mpt_jeta_bkg[0]->Fill(dR_muon,muon_pt,jeta);
                        hdR_bjet_pt_cent_bkg[0]->Fill(dR_muon,jpt,centrality);
                        hdR_bjet_ptrel_cent_bkg[0]->Fill(dR_muon,ptrel_b,centrality);
                        hmuon_meqoverp_idqoverp_meta_bkg[0]->Fill(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000,sin(mu_muid_id_theta->at(particle_ID_muon))/mu_muid_id_qoverp->at(particle_ID_muon)/1000,muon_eta);
                        hmuon_mept_mpt_meta_bkg[0]->Fill(std::abs(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000),muon_pt,muon_eta);
                        hmuon_chi2_mpt_meta_bkg[0]->Fill(mu_muid_matchchi2->at(particle_ID_muon)/mu_muid_matchndof->at(particle_ID_muon),muon_pt,muon_eta);

                        if(dR_muon > pow(10.0,-6))
                            hdR_bjet_mpt_weight_cent_bkg[0]->Fill(dR_muon, muon_pt,centrality,1/dR_muon);
                        if(std::abs(muon_eta)<1.05) 
                        {
                            bjet_eta_pt_bkg[6]->Fill(jeta,jpt);
                            hmuon_eta_bjet_cent_bkg[6]->Fill(muon_eta,centrality);
                            hdR_bjet_mpt_cent_bkg[6]->Fill(dR_muon,muon_pt,centrality);
                            hdR_bjet_mpt_jpt_bkg[6]->Fill(dR_muon,muon_pt,jpt);
                            hdR_bjet_mpt_jeta_bkg[6]->Fill(dR_muon,muon_pt,jeta);
                            hdR_bjet_pt_cent_bkg[6]->Fill(dR_muon,jpt,centrality);
                            hptrel_bjet_pt_bkg[6]->Fill(ptrel_b,jpt);
                            hdR_bjet_ptrel_cent_bkg[6]->Fill(dR_muon,ptrel_b,centrality);
                            hmuon_meqoverp_idqoverp_meta_bkg[6]->Fill(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000,sin(mu_muid_id_theta->at(particle_ID_muon))/mu_muid_id_qoverp->at(particle_ID_muon)/1000,muon_eta);
                            hmuon_mept_mpt_meta_bkg[6]->Fill(std::abs(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000),muon_pt,muon_eta);
                            hmuon_chi2_mpt_meta_bkg[6]->Fill(mu_muid_matchchi2->at(particle_ID_muon)/mu_muid_matchndof->at(particle_ID_muon),muon_pt,muon_eta); 
                        }
                    }
                    if(d_0_err<5.0&&std::abs(d_0)<4.0)
                    {
                        hmuon_IP_mpt_cent[0]->Fill(d_0,muon_pt,centrality);
                        hmuon_err_mpt_cent[0]->Fill(d_0_err,muon_pt,centrality);
                        hmuon_sig_mpt_cent[0]->Fill(d_0_sig,muon_pt,centrality);
                        hmuon_IP_mpt_meta[0]->Fill(d_0,muon_pt,muon_eta);
                        hmuon_sig_mpt_meta[0]->Fill(d_0_sig,muon_pt,muon_eta);

                        if(d_0_sig >= 2.0)
                        {
                            hptrel_bjet_pt[2]->Fill(ptrel_b,jpt);
                            hdR_bjet_ptrel_cent[2]->Fill(dR_muon,ptrel_b,centrality);
                        }
                        else
                        {
                            hptrel_bjet_pt_bkg[2]->Fill(ptrel_b,jpt);
                            hdR_bjet_ptrel_cent_bkg[2]->Fill(dR_muon,ptrel_b,centrality);
                        }
                        if(muon_comp < 0.25 && d_0_sig >=2.0)
                        {
                            bjet_eta_pt[4]->Fill(jeta,jpt);
                            hmuon_eta_bjet_cent[4]->Fill(muon_eta,centrality);
                            hptrel_bjet_pt[4]->Fill(ptrel_b,jpt);
                            hdR_bjet_mpt_cent[4]->Fill(dR_muon,muon_pt,centrality);
                            hdR_bjet_mpt_jpt[4]->Fill(dR_muon,muon_pt,jpt);
                            hdR_bjet_pt_cent[4]->Fill(dR_muon,jpt,centrality);
                            hdR_bjet_ptrel_cent[4]->Fill(dR_muon,ptrel_b,centrality);
                            hmuon_meqoverp_idqoverp_meta[4]->Fill(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000,sin(mu_muid_id_theta->at(particle_ID_muon))/mu_muid_id_qoverp->at(particle_ID_muon)/1000,muon_eta);
                            hmuon_mept_mpt_meta[4]->Fill(std::abs(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000),muon_pt,muon_eta);
                            hmuon_chi2_mpt_meta[4]->Fill(mu_muid_matchchi2->at(particle_ID_muon)/mu_muid_matchndof->at(particle_ID_muon),muon_pt,muon_eta);
                        }
                        else
                        {
                            bjet_eta_pt_bkg[4]->Fill(jeta,jpt);
                            hmuon_eta_bjet_cent_bkg[4]->Fill(muon_eta,centrality);
                            hptrel_bjet_pt_bkg[4]->Fill(ptrel_b,jpt);
                            hdR_bjet_mpt_cent_bkg[4]->Fill(dR_muon,muon_pt,centrality);
                            hdR_bjet_mpt_jpt_bkg[4]->Fill(dR_muon,muon_pt,jpt);
                            hdR_bjet_pt_cent_bkg[4]->Fill(dR_muon,jpt,centrality);
                            hdR_bjet_ptrel_cent_bkg[4]->Fill(dR_muon,ptrel_b,centrality);
                            hmuon_meqoverp_idqoverp_meta_bkg[4]->Fill(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000,sin(mu_muid_id_theta->at(particle_ID_muon))/mu_muid_id_qoverp->at(particle_ID_muon)/1000,muon_eta);
                            hmuon_mept_mpt_meta_bkg[4]->Fill(std::abs(sin(mu_muid_me_theta->at(particle_ID_muon))/mu_muid_me_qoverp->at(particle_ID_muon)/1000),muon_pt,muon_eta);
                            hmuon_chi2_mpt_meta_bkg[4]->Fill(mu_muid_matchchi2->at(particle_ID_muon)/mu_muid_matchndof->at(particle_ID_muon),muon_pt,muon_eta);
                        }
                    }
                }

/////////////////////////////////////////////////minimum bias event mixing//////////////////////////////////////////////

                int b = 0;
                int move = 0;

                float fcal_bins[15] = {6.0,2.9188,2.405,1.9856,1.6351,1.3383,1.0877,0.8721,0.6896,0.5345,0.4070,0.3035,0.2205,0.0456,-0.04};
                float fcal_weight_cor[14];
                float h1_weight_bins[15];
                int h1_find_fcal_bins[15];

                for(int f = 0; f<14; f++)
                    fcal_weight_cor[f] = (float)(GetWeight(h1_weight ,fcal_bins[f])+GetWeight(h1_weight,fcal_bins[f+1]))/2;

                for(int f = 0; f< 15; f++)
                {
                    h1_find_fcal_bins[14-f] = h1_weight->FindBin(fcal_bins[f]);
                    h1_weight_bins[14-f] = fcal_bins[f];
                    if(f >= 1&& f <= 13)
                        h1_weight_bins[f] = h1_weight->GetXaxis()->GetBinLowEdge(h1_find_fcal_bins[f]);
                }

          //      TH1F* h1_weight_new = (TH1F*)h1_weight->Rebin(14,"h1_weight_new",h1_weight_bins);
           //     for(int r = 1; r<15; r++)
            //        std::cout<<h1_weight_new->GetBinContent(r)<<"Ha"<<std::endl;
                
                while(b < 5)
                {
                    int rand = b + move + event_check;
                    rand = rand%(fChain2->GetEntries());
                    fChain2->GetEntry(rand);
                   // std::cout<<" Event : "<<l<<"\t jet : "<<e<<std::endl;
                //    std::cout<<"Rand : "<<rand<<"\t b : "<<b<<"\t move : "<<move<<"\t Fcal: "<<Fcal<<std::endl;
                   // std::cout<<" Centrality: "<<Centrality_comb<<" centrality: "<<centrality<<" fcal_weight: "<<fcal_weight<<std::endl;
                    float fcal_weight = GetWeight(h1_weight , Fcal);
                    Centrality_comb = SetCentrality(GetCentrality(Fcal*1000.0));
                    if(Centrality_comb == centrality)
                    {

                        Fcal_mb->Fill(Fcal);
                        Fcal_mb_weight->Fill(Fcal, fcal_weight);

                        dR_muon = 100.0;
                        particle_ID_muon= -1.0;
                        for(int n = 0; n< Muon_number; n++)
                        {
                            muon_eta = Muon_eta->at(n);
                            muon_phi= Muon_phi->at(n);
                            muon_pt= Muon_pt->at(n)/1000.0;
                            muon_comp= Muon_composite->at(n);
                            d_0_sig = Muon_IP->at(n)/Muon_IP_err->at(n);

                            muon_px = muon_pt*cos(muon_phi);
                            muon_py = muon_pt*sin(muon_phi);
                            muon_pz = muon_pt*sinh(muon_eta);
                            TVector3 Muon_reco(muon_px,muon_py,muon_pz);

                            jpx = jpt*cos(jphi) + muon_px;
                            jpy = jpt*sin(jphi) + muon_py;
                            jpz = jpt*sinh(jeta) + muon_pz;
                            jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                            muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;

                            TVector3 Jet_reco(jpx,jpy,jpz);
                            TVector3 Cross_Jet_Muon_reco;
                            Cross_Jet_Muon_reco = Muon_reco.Cross(Jet_reco);

                            dR_muon_inclu= Association(muon_phi,Jet_reco.Phi(),muon_eta,Jet_reco.PseudoRapidity());
                            ptrel_b= muon_pt * sin(dR_muon_inclu);


                            bool muonDynamiCut = (jpt >60 && jpt < 90 && muon_pt >6) || (jpt >90 && muon_pt >8) || (jpt < 60 && jpt >20);
                            if(dR_muon_inclu < 0.6 )
                            {
                                if(muon_comp < 0.25)
                                {
                                    dR_randm_mpt_cent[1]->Fill(dR_muon_inclu, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    dR_randm_jpt_cent[1]->Fill(dR_muon_inclu, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    dR_randm_jeta_cent[1]->Fill(dR_muon_inclu, muon_eta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    ptrel_randm_jpt_cent[1]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                }
                                else
                                {
                                    dR_randm_mpt_cent_bkg[1]->Fill(dR_muon_inclu, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    dR_randm_jpt_cent_bkg[1]->Fill(dR_muon_inclu, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    dR_randm_jeta_cent_bkg[1]->Fill(dR_muon_inclu, muon_eta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    ptrel_randm_jpt_cent_bkg[1]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                } 
           

                                if(Muon_IP_err->at(n)<5.0&&std::abs(Muon_IP->at(n))<4.0)
                                {
                                    if(d_0_sig >= 2.0)
                                    {
                                        dR_randm_mpt_cent[3]->Fill(dR_muon_inclu, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jpt_cent[3]->Fill(dR_muon_inclu, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jeta_cent[3]->Fill(dR_muon_inclu, muon_eta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        ptrel_randm_jpt_cent[3]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    }
                                    else
                                    {
                                        dR_randm_mpt_cent_bkg[3]->Fill(dR_muon_inclu, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jpt_cent_bkg[3]->Fill(dR_muon_inclu, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jeta_cent_bkg[3]->Fill(dR_muon_inclu, muon_eta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        ptrel_randm_jpt_cent_bkg[3]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    } 
                                    if(d_0_sig >= 2.0 && muon_comp < 0.25)
                                    {
                                        dR_randm_mpt_cent[5]->Fill(dR_muon_inclu, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jpt_cent[5]->Fill(dR_muon_inclu, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jeta_cent[5]->Fill(dR_muon_inclu, muon_eta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        ptrel_randm_jpt_cent[5]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    }
                                    else
                                    {
                                        dR_randm_mpt_cent_bkg[5]->Fill(dR_muon_inclu, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jpt_cent_bkg[5]->Fill(dR_muon_inclu, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        dR_randm_jeta_cent_bkg[5]->Fill(dR_muon_inclu, muon_eta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                        ptrel_randm_jpt_cent_bkg[5]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                    } 
                        

                                }

                                Muon_asso_pt = muon_pt;
                                Muon_asso_eta = muon_eta;
                                Muon_asso_phi= muon_phi;
                                Muon_asso_composite= muon_comp;
                                Muon_asso_sig= d_0_sig;
                                Ptrel= ptrel_b;
                                Flavour_tag = 0;
                                Centrality_hp= centrality;
                                Jet_pt = jpt;
                                Jet_eta= jeta;
                                Jet_phi= jphi;
                                Muon_Jet_dR = dR_muon_inclu;
                                Jet_Trackjet_dR = dR_track_jet;
                                Jet_weight = weight;
                                Jet_ip3d= ip3d;
                                DataJetTree_em->Fill();
                            }

                            if(dR_muon > Association(jphi, muon_phi, jeta, muon_eta))
                            {
                                dR_muon = Association(jphi, muon_phi, jeta, muon_eta);
                                particle_ID_muon = n;
                            }
                        }
            
                        if(particle_ID_muon >= 0 && dR_muon < 0.6)
                        {
                            muon_eta = Muon_eta->at(particle_ID_muon);
                            muon_phi= Muon_phi->at(particle_ID_muon);
                            muon_pt= Muon_pt->at(particle_ID_muon)/1000.0;
                            muon_comp = Muon_composite->at(particle_ID_muon);

                            jpx = jpt*cos(jphi);
                            jpy = jpt*sin(jphi);
                            jpz = jpt*sinh(jeta);
                            muon_px = muon_pt*cos(muon_phi);
                            muon_py = muon_pt*sin(muon_phi);
                            muon_pz = muon_pt*sinh(muon_eta);
                            jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                            muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;
                            ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

      
                            if(muon_comp < 0.25 )
                            {
                                dR_randm_mpt_cent[0]->Fill(dR_muon, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                dR_randm_jpt_cent[0]->Fill(dR_muon, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                dR_randm_jeta_cent[0]->Fill(dR_muon, jeta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                ptrel_randm_jpt_cent[0]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                            }
                            else 
                            {
                                dR_randm_mpt_cent_bkg[0]->Fill(dR_muon, muon_pt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                dR_randm_jpt_cent_bkg[0]->Fill(dR_muon, jpt, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                dR_randm_jeta_cent_bkg[0]->Fill(dR_muon, jeta, centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                                ptrel_randm_jpt_cent_bkg[0]->Fill(ptrel_b, jpt,centrality,(float)fcal_weight/fcal_weight_cor[centrality]);
                            }

                        }

                        b++;
                    }
                    else
                        move ++;
                } 


                int Pool_size = pool[centrality].size();
                if(Pool_size > 0 && Pool_size <= depth[centrality])
            //    if( Pool_size == depth[centrality])
                {

                    for(int r = 0; r < Pool_size; r++)
                    {

                        Jet_jpt_jeta_cent_check->Fill(jpt, jeta, centrality);

                        Event* prev_event;
                        prev_event = pool[centrality].at(r);

                        dR_muon= 100.0;
                        particle_ID_muon = -1.0;

                        for(unsigned int n = 0; n< prev_event->RandMuon.size(); n++)
                        {
                            muon_eta = prev_event->RandMuon.at(n).eta;
                            muon_phi= prev_event->RandMuon.at(n).phi;
                            muon_pt= prev_event->RandMuon.at(n).pt/1000.0;
                            muon_comp = prev_event->RandMuon.at(n).composite;
                            dR_muon_inclu = Association(jphi, muon_phi, jeta, muon_eta);

                            jpx = jpt*cos(jphi);
                            jpy = jpt*sin(jphi);
                            jpz = jpt*sinh(jeta);
                            muon_px = muon_pt*cos(muon_phi);
                            muon_py = muon_pt*sin(muon_phi);
                            muon_pz = muon_pt*sinh(muon_eta);
                            jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                            muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;
                            ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

        
                            if(dR_muon_inclu < 0.6)
                            {
                                if(muon_comp < 0.25)
                                {
                                    dR_randm_mpt_cent[1]->Fill(dR_muon_inclu, muon_pt, centrality);
                                    dR_randm_jpt_cent[1]->Fill(dR_muon_inclu, jpt, centrality);
                                    ptrel_randm_jpt_cent[1]->Fill(ptrel_b, jpt,centrality);
                                }
                                else
                                {
                                    dR_randm_mpt_cent_bkg[1]->Fill(dR_muon_inclu, muon_pt, centrality);
                                    dR_randm_jpt_cent_bkg[1]->Fill(dR_muon_inclu, jpt, centrality);
                                    ptrel_randm_jpt_cent_bkg[1]->Fill(ptrel_b, jpt,centrality);
                                } 
                            }


                            if(dR_muon > Association(jphi, muon_phi, jeta, muon_eta))
                            {
                                dR_muon = Association(jphi, muon_phi, jeta, muon_eta);
                                particle_ID_muon = n;
                            }
                        }
            
                        if(particle_ID_muon >= 0 && dR_muon < 0.6)
                        {
                            muon_eta = prev_event->RandMuon.at(particle_ID_muon).eta;
                            muon_phi= prev_event->RandMuon.at(particle_ID_muon).phi;
                            muon_pt= prev_event->RandMuon.at(particle_ID_muon).pt/1000.0;
                            muon_comp = prev_event->RandMuon.at(particle_ID_muon).composite;

                            jpx = jpt*cos(jphi);
                            jpy = jpt*sin(jphi);
                            jpz = jpt*sinh(jeta);
                            muon_px = muon_pt*cos(muon_phi);
                            muon_py = muon_pt*sin(muon_phi);
                            muon_pz = muon_pt*sinh(muon_eta);
                            jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                            muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;
                            ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

      
                        if(muon_comp < 0.25 && muon_comp >= 0.0)
                            dR_randm_num_jet_cent[0]->Fill(dR_muon,1,centrality);
                        else if (muon_comp >=0.25) 
                            dR_randm_num_jet_cent_bkg[0]->Fill(dR_muon,1,centrality);


        
                            if(muon_comp < 0.25 )
                            {
                                dR_randm_mpt_cent[0]->Fill(dR_muon, muon_pt, centrality);
                                dR_randm_jpt_cent[0]->Fill(dR_muon, jpt, centrality);
                                ptrel_randm_jpt_cent[0]->Fill(ptrel_b, jpt,centrality);
                            }
                            else 
                            {
                                dR_randm_mpt_cent_bkg[0]->Fill(dR_muon, muon_pt, centrality);
                                dR_randm_jpt_cent_bkg[0]->Fill(dR_muon, jpt, centrality);
                                ptrel_randm_jpt_cent_bkg[0]->Fill(ptrel_b, jpt,centrality);
                            }
                        }
                    }
                }

            }
        }

        event_check += 5.0;
                
///////////////////////////////////////////////event mixing///////////////////////////////////////////////////
       
        for(int k =0; k< antikt2HIItrEM_n; k++)
        {
            Caljet.SetPtEtaPhiM(antikt2HIItrEM_pt->at(k),antikt2HIItrEM_eta->at(k),antikt2HIItrEM_phi->at(k),antikt2HIItrEM_m->at(k));
            calib_factor_original=Caljet.Et()/antikt2HIItrEM_UncalibEt->at(k);
            Caljet*=1.0/calib_factor_original;
            Caljet*=1e-3; // convert from MeV to GeV
            m_calibrator->Calibrate(Caljet);
            jetCut = false;
            jet_trig_asso = false;

            jeta = antikt2HIItrEM_eta->at(k);
            jphi = antikt2HIItrEM_phi->at(k);
            jpt = Caljet.Pt();
            jetCut = std::abs(jeta) < 2.5 && jpt > 20.0;

            for(int t = 0; t< trig_EF_jet_n; t++)
            {
                if(trig_EF_jet_pt->at(t)/1000.0 > 20.0 && trig_EF_jet_a2hi->at(t) == 1)
                {
                    if(Association(jphi, trig_EF_jet_phi->at(t), jeta, trig_EF_jet_eta->at(t)) < 0.4)
                    {
                        jet_trig_asso = true;
                        break;
                    }
                }
            }


            if(trig_EF[0] && jetCut && jet_trig_asso)
            {
                Jet_jpt_jeta_cent_check->Fill(jpt, jeta, centrality);

                int Pool_size = pool[centrality].size();
                if(Pool_size > 0 && Pool_size <= depth[centrality])
                {

                    for(int r = 0; r < Pool_size; r++)
                    {

                        Event* prev_event;
                        prev_event = pool[centrality].at(r);

                        dR_muon= 100.0;
                        particle_ID_muon = -1.0;

                        for(unsigned int n = 0; n< prev_event->RandMuon.size(); n++)
                        {
                            muon_eta = prev_event->RandMuon.at(n).eta;
                            muon_phi= prev_event->RandMuon.at(n).phi;
                            muon_pt= prev_event->RandMuon.at(n).pt/1000.0;
                            muon_comp = prev_event->RandMuon.at(n).composite;
                            dR_muon_inclu = Association(jphi, muon_phi, jeta, muon_eta);

                            jpx = jpt*cos(jphi);
                            jpy = jpt*sin(jphi);
                            jpz = jpt*sinh(jeta);
                            muon_px = muon_pt*cos(muon_phi);
                            muon_py = muon_pt*sin(muon_phi);
                            muon_pz = muon_pt*sinh(muon_eta);
                            jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                            muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;
                            ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

        
                            if(dR_muon_inclu < 0.6)
                            {
                                if(muon_comp < 0.25)
                                {
                                    dR_randm_mpt_cent[1]->Fill(dR_muon_inclu, muon_pt, centrality);
                                    dR_randm_jpt_cent[1]->Fill(dR_muon_inclu, jpt, centrality);
                                    ptrel_randm_jpt_cent[1]->Fill(ptrel_b, jpt,centrality);
                                }
                                else
                                {
                                    dR_randm_mpt_cent_bkg[1]->Fill(dR_muon_inclu, muon_pt, centrality);
                                    dR_randm_jpt_cent_bkg[1]->Fill(dR_muon_inclu, jpt, centrality);
                                    ptrel_randm_jpt_cent_bkg[1]->Fill(ptrel_b, jpt,centrality);
                                } 
                            }


                            if(dR_muon > Association(jphi, muon_phi, jeta, muon_eta))
                            {
                                dR_muon = Association(jphi, muon_phi, jeta, muon_eta);
                                particle_ID_muon = n;
                            }
                        }
            
                        if(particle_ID_muon >= 0 && dR_muon < 0.6)
                        {
                            muon_eta = prev_event->RandMuon.at(particle_ID_muon).eta;
                            muon_phi= prev_event->RandMuon.at(particle_ID_muon).phi;
                            muon_pt= prev_event->RandMuon.at(particle_ID_muon).pt/1000.0;
                            muon_comp = prev_event->RandMuon.at(particle_ID_muon).composite;

                            jpx = jpt*cos(jphi);
                            jpy = jpt*sin(jphi);
                            jpz = jpt*sinh(jeta);
                            muon_px = muon_pt*cos(muon_phi);
                            muon_py = muon_pt*sin(muon_phi);
                            muon_pz = muon_pt*sinh(muon_eta);
                            jp= sqrt(jpx*jpx+jpy*jpy+jpz*jpz);
                            muon_proj = (muon_px*jpx+muon_py*jpy+muon_pz*jpz)/jp;
                            ptrel_b = sqrt(muon_px*muon_px+muon_py*muon_py+muon_pz*muon_pz-muon_proj*muon_proj);

      
                        if(muon_comp < 0.25 && muon_comp >= 0.0)
                            dR_randm_num_jet_cent[0]->Fill(dR_muon,1,centrality);
                        else if (muon_comp >=0.25) 
                            dR_randm_num_jet_cent_bkg[0]->Fill(dR_muon,1,centrality);


        
                            if(muon_comp < 0.25 )
                            {
                                dR_randm_mpt_cent[0]->Fill(dR_muon, muon_pt, centrality);
                                dR_randm_jpt_cent[0]->Fill(dR_muon, jpt, centrality);
                                ptrel_randm_jpt_cent[0]->Fill(ptrel_b, jpt,centrality);
                            }
                            else 
                            {
                                dR_randm_mpt_cent_bkg[0]->Fill(dR_muon, muon_pt, centrality);
                                dR_randm_jpt_cent_bkg[0]->Fill(dR_muon, jpt, centrality);
                                ptrel_randm_jpt_cent_bkg[0]->Fill(ptrel_b, jpt,centrality);
                            }
                        }
                    }
                }
            }
        }




       
        int Pool_size2 = pool[centrality].size();
        if(Pool_size2 >= depth[centrality])
        {
            delete pool[centrality].at(0);
            pool[centrality].erase( pool[centrality].begin());
        }

   //     if (this_event->RandMuon.size() == 0) delete this_event;
    //    else 
        pool[centrality].push_back(this_event);

    }

}

void AA_Data11_analysis::WriteTree(int infile)
{
    std::stringstream output_file;
    output_file.str("");
    output_file<<"DataJetTree_"<<infile<<".root";
    TFile *fout=new TFile(output_file.str().c_str(),"RECREATE");
    
    DataJetTree->Write();

    fout->Close();

}
void AA_Data11_analysis::WriteGraph(int infile)
{
    std::stringstream output_file;
    output_file.str("");
    output_file<<"histo_"<<infile<<".root";
    TFile *fout=new TFile(output_file.str().c_str(),"RECREATE");

    DataJetTree->Write();
    DataJetTree_em->Write();
    //TFile *fout=new TFile("histogram.root","RECREATE");
    for(int w = 0; w<8; w++)
    {
        Number_Event[w]->Write();
        Number_Muon[w]->Write();
        Number_Muon_bkg[w]->Write();
        hptrel_bjet[w]->Write();
        bjet_eta_cent[w]->Write();
        hptrel_bjet_pt[w]->Write();
        bjet_eta_pt[w]->Write();
        hmuon_eta_bjet_cent[w]->Write();
        hdR_bjet_mpt_cent[w]->Write();
        hdR_bjet_mpt_jpt[w]->Write();
        hdR_bjet_mpt_jeta[w]->Write();
        hdR_bjet_pt_cent[w]->Write();
        hdR_bjet_ptrel_cent[w]->Write();
        dR_randm_mpt_cent[w]->Write();
        dR_randm_jpt_cent[w]->Write();
        dR_randm_jeta_cent[w]->Write();
        ptrel_randm_jpt_cent[w]->Write();
        dR_randm_num_jet_cent[w]->Write();
        hdR_bjet_mpt_weight_cent[w]->Write();
        hmuon_meqoverp_idqoverp_meta[w]->Write();
        hmuon_mept_mpt_meta[w]->Write();
        hmuon_chi2_mpt_meta[w]->Write();


        hptrel_bjet_pt_bkg[w]->Write();
        bjet_eta_pt_bkg[w]->Write();
        hmuon_eta_bjet_cent_bkg[w]->Write();
        hdR_bjet_mpt_cent_bkg[w]->Write();
        hdR_bjet_mpt_jpt_bkg[w]->Write();
        hdR_bjet_mpt_jeta_bkg[w]->Write();
        hdR_bjet_pt_cent_bkg[w]->Write();
        hdR_bjet_ptrel_cent_bkg[w]->Write();
        hdR_bjet_mpt_weight_cent_bkg[w]->Write();
        dR_randm_mpt_cent_bkg[w]->Write();
        dR_randm_jpt_cent_bkg[w]->Write();
        dR_randm_jeta_cent_bkg[w]->Write();
        ptrel_randm_jpt_cent_bkg[w]->Write();
        dR_randm_num_jet_cent_bkg[w]->Write();
        hmuon_meqoverp_idqoverp_meta_bkg[w]->Write();
        hmuon_mept_mpt_meta_bkg[w]->Write();
        hmuon_chi2_mpt_meta_bkg[w]->Write();

        hmuon_composite_bjet_ptrel[w]->Write();
        hmuon_composite_mpt_meta[w]->Write();
        hmuon_meqoverp_mpt_meta[w]->Write();
        hmuon_eloss_mpt_meta[w]->Write();
        hmuon_scat_mpt_meta[w]->Write();
        hbjeta_phi_cent[w]->Write();
        hbjeta_phi_pt[w]->Write();
        hmuon_IP_mpt_cent[w]->Write();
        hmuon_sig_mpt_cent[w]->Write();
        hmuon_err_mpt_cent[w]->Write();
        hmuon_IP_mpt_meta[w]->Write();
        hmuon_sig_mpt_meta[w]->Write();
    }
  
    Jet_jpt_jeta_cent->Write();   
    Jet_jpt_jeta_cent_check->Write();   
    Fcal_hp->Write();
    Fcal_mb->Write();
    Fcal_mb_weight->Write();
    test_jeta ->Write(); 
    hmuon_pt->Write();    
    hmuon_eta->Write();    

    bjet_pt->Write();
    bjet_eta->Write();
    hmuon_eta_bjet->Write();
    hmuon_pt_bjet->Write();
    hmass_bjet_0->Write();
    hmass_bjet_7->Write();
    hweight_bjet->Write();
    d0->Write();
    d0_err->Write();
    d0_sig->Write();
    d0_comp->Write();
    d0_err_comp->Write();
    d0_sig_comp->Write();

    hdR_bjet_pt->Write();
    hdR_bjet_mpt->Write();
    //hdR_bjet_mpt_weight->Write();
    hptrel_bjet_cent->Write();
    hjpt_bjet_cent->Write();
    hmpt_bjet_cent->Write();
    hptrel_bjet_pt_cent->Write();


    fout->Close();

}

float AA_Data11_analysis::Association(float phi1, float phi2,float eta1, float eta2)
{
    float dPhi = phi1 - phi2;
    float dEta = eta1 - eta2;
    float Delta_R;
    float pi = TMath::Pi();

    if(std::abs(dPhi)>= pi)  
        Delta_R=sqrt((dEta)*(dEta)+(2*pi-std::abs(dPhi))*(2*pi-std::abs(dPhi)));

    else 
        Delta_R=sqrt((dEta)*(dEta)+(dPhi)*(dPhi));
    
    return Delta_R;
}


bool AA_Data11_analysis::SelectMuon(int index)
{
    if(mu_muid_isCombinedMuon->at(index) != 1) return false;
    if( mu_muid_pt->at(index)/1000.<4 ) return false;
    if( std::abs(mu_muid_eta->at(index))>2.5 ) return false;
    if( mu_muid_nPixHits->at(index)<1 ) return false;
    if( mu_muid_nBLHits->at(index)<1 ) return false;
    if( mu_muid_nSCTHits->at(index)<7 ) return false;
    if( (mu_muid_nPixHoles->at(index)+mu_muid_nSCTHoles->at(index))>1 ) return false;
    if( mu_muid_nSCTHoles->at(index)>0 ) return false;
    if( std::abs(mu_muid_id_d0_exPV->at(index))>5. ) return false;
    if( std::abs(mu_muid_id_z0_exPV->at(index))>5. ) return false;
    if( std::abs(1./mu_muid_id_qoverp->at(index)/1000.)<3. ) return false;
    if( std::abs(1./mu_muid_me_qoverp->at(index)/1000.)<0.1 ) return false;
    if( mu_muid_matchchi2->at(index) < 0.0 ) return false;
    if( muid_scatSignificance->at(index) < -900.0 ) return false;
    if( (mu_muid_matchchi2->at(index)/mu_muid_matchndof->at(index))>10 ) return false;
    return true;

}


float AA_Data11_analysis::DeltaPhi(Float_t phi1, Float_t phi2)
{
    Float_t dPhi = phi1 - phi2;
    Float_t pi = TMath::Pi();
    while (dPhi < -pi) dPhi += 2*pi;
    while (dPhi >  pi) dPhi -= 2*pi;
    return dPhi;
}

//void AA_Data11_analysis::SetJESCalibrator(int R, std::string tag, std::string file_name)
void AA_Data11_analysis::SetJESCalibrator(std::string tag, std::string file_name)
{
    m_calibrator=new HIJetCalibTool(tag,file_name);
}

bool AA_Data11_analysis::PassesEventCuts() const
{

  //check if LumiBlock is good
  //if you didn't give the AA_Data11_analysis a GRL, don't apply LumiBlock selection
  if (m_GRL_reader)
  {
    if (! (m_GRL_reader->GetMergedGoodRunsList().HasRunLumiBlock(RunNumber, lbn)) ) 
    {
        std::cout<<"error GRL"<<std::endl;
        return false;
    }
  }

//   if(m_trigger_id!=0) //ie you chose to set it youself
//   {
//     bool has_trigger=false;
//     for(unsigned int itrig=0; itrig<trig_EF_passedPhysics->size(); itrig++)
//     {
//       if(trig_EF_passedPhysics->at(itrig)==m_trigger_id) //your trigger is in the list
//       {
    //      has_trigger=true;
    //      break;
//       }
//     }
//     if(!has_trigger) return false;
//   }

  //if there is no time measurement in one of the modules the time variable may be set to 0 or 75
  //make sure there is a good time in both sides and that the time difference is less than some timing cut.
  if ( (mbtime_timeA  == 0. ) ||    (mbtime_timeC  == 0. ) || 
       (std::abs(mbtime_timeA) == 75.) || (std::abs(mbtime_timeC) == 75.) || 
       (std::abs(mbtime_timeA - mbtime_timeC) > m_time_cut) ) 
    {
        std::cout<<"error mbtime"<<std::endl;
        std::cout<<m_time_cut<<"time_cut"<<std::endl;
        std::cout<<mbtime_timeA<<"timeA"<<std::endl;
        std::cout<<mbtime_timeC<<"timeC"<<std::endl;
        return false;

    }

  //require at least one reconstructed primary vertex
  if(vx_n <  2) 
  {
    std::cout<<"error vx_n"<<std::endl;
    return false;
  }
/*
  if(m_reject_pileup)
  {
    //reject pileup
    //count number of "strong" vertices- vertices built from tracks
    //whose sum pT > 6 GeV or have more than 10 tracks
    int vertexCount = 0;
    for (int ivx = 0; ivx < vx_n - 1; ivx++)
    {
      if (vx_sumPt->at(ivx)/1000. > 6. ||
          vx_nTracks->at(ivx) > 10)
      {
        vertexCount++;
      }
    }
    if(vertexCount > 1) return false; //expect 1 (primary) more than one is probably pileup
  }
*/

  //if event survives, its passed the cuts and return true
  return true;

}

void AA_Data11_analysis::FillFcal(int infile)
{
    int entries = fChain->GetEntries();
    std::cout<<entries<<std::endl;

    h1_fcal_hp = new TH1F("h1_fcal_hp","h1_fcal_hp",80,-0.1,6.1);

    for(int l = 0; l< entries; l++)
    {
        fChain->GetEntry(l);
        if(!PassesEventCuts()) 
            continue;
        h1_fcal_hp->Fill(GetCentralityVariable()/1000.0);
    }
    
    std::stringstream output_file;
    output_file.str("");
    output_file<<"Fcalhp_"<<infile<<".root";
    TFile *fout=new TFile(output_file.str().c_str(),"RECREATE");

    h1_fcal_hp->Write();

    fout->Close();

}

float AA_Data11_analysis::GetCentralityVariable() const
{
  //sum fcal layers and convert GeV to MeV
  float fcal_sum_positive_side=(cccEt_p_Eh_FCal->at(0)+cccEt_p_Eh_FCal->at(1)+cccEt_p_Eh_FCal->at(2))/1e3;
  float fcal_sum_negative_side=(cccEt_n_Eh_FCal->at(0)+cccEt_n_Eh_FCal->at(1)+cccEt_n_Eh_FCal->at(2))/1e3;
  if(RunNumber <=196000) //Pb+Pb collisions
  {
    float fcal_sum_total=fcal_sum_positive_side+fcal_sum_positive_side;
    if(RunNumber >= 192000) //2011 Pb+Pb
    {   
      fcal_sum_total*=1.04; //4% scale difference between 2010 and 2011
    }   
    return fcal_sum_total;
  }

  if(RunNumber >= 218677 && RunNumber <=219144) //Pb in beam 1, p in beam 2
  {
    return 0.983*fcal_sum_negative_side; //account for scale difference between two sides;
  }
  return fcal_sum_positive_side;
}

float AA_Data11_analysis::GetCentrality(float fcal_centrality) const
{
//  if(m_centrality_getter) return m_centrality_getter->GetCentralityBinCenter(GetCentralityVariable());
  if(m_centrality_getter) return m_centrality_getter->GetCentralityBinCenter(fcal_centrality);
  std::cerr << "Centrality getter not set" << std::endl;
  return -1; 
}

void AA_Data11_analysis::SetCentralityGetter(std::string cent_calib_file)    
{
  m_centrality_getter=new CentralityGetter(cent_calib_file);
  std::cout << "Using centrality getter with bins " << std::endl;
  m_centrality_getter->Print();
}

void AA_Data11_analysis::SetGRLReader(std::string grl_file_name)
{
  m_GRL_reader=new Root::TGoodRunsListReader();
  m_GRL_reader->SetXMLFile(grl_file_name.c_str());
  m_GRL_reader->Interpret();
}

float AA_Data11_analysis::GetWeight(TH1F* h1_weight, float fcal)
{
    return h1_weight->Interpolate(fcal);
}

int AA_Data11_analysis::SetCentrality(float centrality_getter)
{
    int centrality_trans = -1; 

    if (centrality_getter<5.0 && centrality_getter > 0.0 ) centrality_trans =0; 
    else if (10.0>centrality_getter&&centrality_getter>5.0) centrality_trans =1; 
    else if (15.0>centrality_getter&&centrality_getter>10.0) centrality_trans =2; 
    else if (20.0>centrality_getter&&centrality_getter>15.0) centrality_trans =3; 
    else if (25.0>centrality_getter&&centrality_getter>20.0) centrality_trans =4; 
    else if (30.0>centrality_getter&&centrality_getter>25.0) centrality_trans =5; 
    else if (35.0>centrality_getter&&centrality_getter>30.0) centrality_trans =6; 
    else if (40.0>centrality_getter&&centrality_getter>35.0) centrality_trans =7; 
    else if (45.0>centrality_getter&&centrality_getter>40.0) centrality_trans =8; 
    else if (50.0>centrality_getter&&centrality_getter>45.0) centrality_trans =9; 
    else if (55.0>centrality_getter&&centrality_getter>50.0) centrality_trans =10; 
    else if (60.0>centrality_getter&&centrality_getter>55.0) centrality_trans =11; 
    else if (80.0>centrality_getter&&centrality_getter>60.0) centrality_trans =12; 
    else  centrality_trans =13; 

    return centrality_trans;

}

float AA_Data11_analysis::Calerror(int index)
{
    float sv0p_x = antikt2HIItrEMFR_flavor_component_sv0p_x->at(index);
    float sv0p_y = antikt2HIItrEMFR_flavor_component_sv0p_y->at(index);
    float sv0p_z = antikt2HIItrEMFR_flavor_component_sv0p_z->at(index);
    float sv0p_err_x = antikt2HIItrEMFR_flavor_component_sv0p_err_x->at(index);
    float sv0p_err_y = antikt2HIItrEMFR_flavor_component_sv0p_err_y->at(index);
    float sv0p_err_z = antikt2HIItrEMFR_flavor_component_sv0p_err_z->at(index);
    float sv0p_cov_xy = antikt2HIItrEMFR_flavor_component_sv0p_cov_xy->at(index);
    float sv0p_cov_yz = antikt2HIItrEMFR_flavor_component_sv0p_cov_yz->at(index);
    float sv0p_cov_xz = antikt2HIItrEMFR_flavor_component_sv0p_cov_xz->at(index);
    float pvx_x = vx_x->at(0);
    float pvx_y = vx_y->at(0);
    float pvx_z = vx_z->at(0);
    float pvx_cov_x = sqrt(vx_cov_x->at(0));
    float pvx_cov_y = sqrt(vx_cov_y->at(0));
    float pvx_cov_z = sqrt(vx_cov_z->at(0));
    float pvx_cov_xy = vx_cov_xy->at(0);
    float pvx_cov_yz = vx_cov_yz->at(0);
    float pvx_cov_xz = vx_cov_xz->at(0);

	float length;
	float error;

	length = sqrt((pvx_x-sv0p_x)*(pvx_x-sv0p_x) + (pvx_y-sv0p_y)*(pvx_y-sv0p_y) + (pvx_z-sv0p_z)*(pvx_z-sv0p_z));

	error = sqrt((pvx_x-sv0p_x)*(pvx_x-sv0p_x)*sv0p_err_x*sv0p_err_x + (pvx_y-sv0p_y)*(pvx_y-sv0p_y)*sv0p_err_y*sv0p_err_y + (pvx_z-sv0p_z)*(pvx_z-sv0p_z)*sv0p_err_z*sv0p_err_z + 2*(pvx_x-sv0p_x)*(pvx_y-sv0p_y)*sv0p_cov_xy+ 2*(pvx_y-sv0p_y)*(pvx_z-sv0p_z)*sv0p_cov_yz+ 2*(pvx_x-sv0p_x)*(pvx_z-sv0p_z)*sv0p_cov_xz+ (pvx_x-sv0p_x)*(pvx_x-sv0p_x)*pvx_cov_x*pvx_cov_x+ (pvx_y-sv0p_y)*(pvx_y-sv0p_y)*pvx_cov_y*pvx_cov_y+ (pvx_z-sv0p_z)*(pvx_z-sv0p_z)*pvx_cov_z*pvx_cov_z+ 2*(pvx_x-sv0p_x)*(pvx_y-sv0p_y)*pvx_cov_xy+ 2*(pvx_y-sv0p_y)*(pvx_z-sv0p_z)*pvx_cov_yz+ 2*(pvx_x-sv0p_x)*(pvx_z-sv0p_z)*pvx_cov_xz)/length;

	
	return error;
}
