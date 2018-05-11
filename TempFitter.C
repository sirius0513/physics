/* vim: set sw=4 ts=4 expandtab: */

#include "TempFitter.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <TChain.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <RooWorkspace.h>
#include <RooKeysPdf.h>
#include <RooAbsPdf.h>
#include <RooAbsArg.h>
#include <RooDataSet.h>
#include <RooRealVar.h>
#include <RooNDKeysPdf.h>
#include <Roo2DKeysPdf.h>
#include <RooAddPdf.h>
#include <RooGaussian.h>
#include <RooFFTConvPdf.h>
#include <RooFitResult.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <RooExtendPdf.h>
#include <RooFormulaVar.h>
#include <RooCategory.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TH2D.h>
#include <RooStats/LikelihoodInterval.h>
#include <TLegendEntry.h>
#include "RooShiftedKeysPdf.h"
#include <RooChi2Var.h>
#include <RooMinuit.h>
#include <RooNLLVar.h>
#include <RooStats/ProfileLikelihoodCalculator.h>
#include "/atlas/u/tiwang/atlasstyle-00-03-05/AtlasStyle.C"

TempFitter::TempFitter(bool doBuildTemplate, bool doComposite, bool doIPSignificance, int doMuonPt, bool doPP):m_doBuildTemplate(doBuildTemplate),m_doComposite(doComposite),m_doIPSignificance(doIPSignificance),m_doMuonPt(doMuonPt),m_doPP(doPP)
{
    m_ptBins.push_back(20.0);

    m_ptBins.push_back(30.0);
    //m_ptBins.push_back(45.0);
    m_ptBins.push_back(50.0);
//   m_ptBins.push_back(60.0);
   
    m_ptBins.push_back(75.0);
    m_ptBins.push_back(105.);

    m_ptBins.push_back(150.);

    if(m_doPP)
    
    {
        m_centralityBins.push_back(-1.5);
       
        m_centralityBins.push_back(13.5);
    }

    else
    {
        m_centralityBins.push_back(-1.5);
        m_centralityBins.push_back(1.5);
        
        m_centralityBins.push_back(3.5);
        m_centralityBins.push_back(7.5);
        m_centralityBins.push_back(11.5);
        m_centralityBins.push_back(13.5);
        
    }

    m_muonptBins.push_back(4.0);
    m_muonptBins.push_back(5.0);
    m_muonptBins.push_back(6.0);
    m_muonptBins.push_back(8.0);
    m_muonptBins.push_back(15.0);
    m_muonptBins.push_back(50.0);

    file_FitPlot= new TFile("TemplateFittingPlot.root","RECREATE");
}

void TempFitter::DoTempFitting()
{
// get different JetTree from direct sample
    TFile* fin  = TFile::Open("/atlas/data4/tingting/data6/ReconCode/histogram_Data11_AA_minbia_randm_aveweight_addmuon_sig_ip3d.root");
    TFile* fin1  = TFile::Open("/atlas/data4/tingting/data6/ReconCode/histogram_track_inclusive_recopt_J1.root");
    TFile* fin2  = TFile::Open("/atlas/data4/tingting/data6/ReconCode/histogram_track_inclusive_recopt_J2.root");
    TFile* fin3  = TFile::Open("/atlas/data4/tingting/data6/ReconCode/histogram_track_inclusive_recopt_J3.root");

    TTree* JetTree = (TTree*)fin->Get("DataJetTree");
    TTree* MCTree_inclusive_1 = (TTree*)fin1->Get("MCJetTree");
    TTree* MCTree_inclusive_2 = (TTree*)fin2->Get("MCJetTree");
    TTree* MCTree_inclusive_3 = (TTree*)fin3->Get("MCJetTree");


// define dimensions for variable
    RooRealVar  Ptrel("Ptrel","Ptrel",0.0,4.0);
    RooRealVar  Jet_pt("Jet_pt","Jet_pt",19.5,1e9); 
    RooRealVar  Jet_eta("Jet_eta","Jet_eta",-3.0,+3.0);
    RooRealVar  Jet_weight("Jet_weight","Jet_weight",-100,+100);
    RooRealVar  Muon_asso_composite("Muon_asso_composite","Muon_asso_composite",0.0,10.0);
    RooRealVar  Muon_asso_sig("Muon_asso_sig","Muon_asso_sig",-10.0,10.0);
    RooRealVar  Muon_asso_pt("Muon_asso_pt","Muon_asso_pt",3.5,100);
    RooRealVar  Centrality_hp("Centrality_hp","Centrality_hp",-0.5,14);
    RooRealVar  Flavour_tag("Flavour_tag","Flavour_tag",-0.5,10);
    RooRealVar  Jet_Trackjet_dR("Jet_Trackjet_dR","Jet_Trackjet_dR",0.0,100);
    RooRealVar  Muon_Jet_dR("Muon_Jet_dR","Muon_Jet_dR",0.0,1);
    RooRealVar  weight("weight","weight",-1.0,10);


    RooRealVar  var = Ptrel;
    RooRealVar  var_2D_comp = Muon_asso_pt;
    RooRealVar  var_2D_sig= Muon_asso_sig;
    RooArgSet JetArgSet(var,Jet_pt,Muon_asso_composite,Jet_weight,Muon_asso_pt,Centrality_hp,Flavour_tag,Muon_Jet_dR,weight);
    RooArgSet JetArgSet_comb(var,Jet_pt,Muon_asso_composite,Muon_asso_sig,Muon_asso_pt,Centrality_hp,Flavour_tag,Jet_Trackjet_dR);
    RooRealVar sigRatio("sigRatio","sigRatio",0.3,0,1.0);
    RooRealVar bkgCRatio("bkgCRatio","bkgCRatio",0.3,0,1.0);
    RooRealVar bkgPRatio("bkgPRatio","bkgPRatio",0.3,0,1.0);
    RooRealVar bkgKRatio("bkgKRatio","bkgKRatio",0.3,0,1.0);

    RooDataSet *mcSet_inclusive_J1= new RooDataSet("mcSet_inclusive_J1","mcSet_inclusive_J1",MCTree_inclusive_1,JetArgSet);
    RooDataSet *mcSet_inclusive_J2= new RooDataSet("mcSet_inclusive_J2","mcSet_inclusive_J2",MCTree_inclusive_2,JetArgSet);
    RooDataSet *mcSet_inclusive_J3= new RooDataSet("mcSet_inclusive_J3","mcSet_inclusive_J3",MCTree_inclusive_3,JetArgSet);

    TFile *fCal_weight = TFile::Open("/atlas/data3/angerami/analysis_2011/reprocessed/skimmed_NTUP_HI/mc11_2TeV/JetSkim.013/histos_2011_bins_t1/calib/fcal.weights.root");
    TH1F* h1_weight = (TH1F*)fCal_weight->Get("h1_weights");

    TFile *fTrack_fraction= TFile::Open("/atlas/data4/tingting/data6/ReconCode/Trackjet_fraction_dataoverMC.root");
    TH1F* h1_track_ratio = (TH1F*)fTrack_fraction->Get("h1_jpt_ratio");

    TFile *fNoTrack_fraction= TFile::Open("/atlas/data4/tingting/data6/ReconCode/NoTrackjet_fraction_dataoverMC.root");
    TH1F* h1_notrack_ratio = (TH1F*)fNoTrack_fraction->Get("h1_jpt_ratio");

    TFile *fTrack_fraction_muon= TFile::Open("/atlas/data4/tingting/data6/ReconCode/Trackjet_fraction_dataoverMC_muon.root");
    TH1F* h1_track_ratio_muon = (TH1F*)fTrack_fraction_muon->Get("h1_jpt_ratio");

    TFile *fNoTrack_fraction_muon= TFile::Open("/atlas/data4/tingting/data6/ReconCode/NoTrackjet_fraction_dataoverMC_muon.root");
    TH1F* h1_notrack_ratio_muon = (TH1F*)fNoTrack_fraction_muon->Get("h1_jpt_ratio");
   // Fill data set and mc set 
    RooDataSet *mcSet_muon= new RooDataSet("mcSet_muon","mcSet_muon",JetArgSet,"weight");
    FillDataSet(0,0,0,JetArgSet, mcSet_muon, 0.0,1.0,1, h1_weight, h1_track_ratio_muon, h1_notrack_ratio_muon);
    FillDataSet(0,0,0,JetArgSet, mcSet_muon, 0.0,1.0,2, h1_weight, h1_track_ratio_muon, h1_notrack_ratio_muon);
    FillDataSet(0,0,0,JetArgSet, mcSet_muon, 0.0,1.0,3, h1_weight, h1_track_ratio_muon, h1_notrack_ratio_muon);
//    mcSet_muon->Print();
   
    RooDataSet *mcSet_inclusive= new RooDataSet("mcSet_inclusive","mcSet_inclusive",JetArgSet,"weight");
    FillDataSet(0,0,1,JetArgSet, mcSet_inclusive, 0.0,1.0,1, h1_weight, h1_track_ratio, h1_notrack_ratio);
    FillDataSet(0,0,1,JetArgSet, mcSet_inclusive, 0.0,1.0,2, h1_weight, h1_track_ratio, h1_notrack_ratio);
    if(!m_doPP)
        FillDataSet(0,0,1,JetArgSet, mcSet_inclusive, 0.0,1.0,3, h1_weight, h1_track_ratio, h1_notrack_ratio);
    RooDataSet *mcSet_inclusive_ref= new RooDataSet("mcSet_inclusive_ref","mcSet_inclusive_ref",JetArgSet,"weight");
    FillDataSet(0,0,1,JetArgSet, mcSet_inclusive_ref, 0.5,1.0,1, h1_weight, h1_track_ratio, h1_notrack_ratio);
    FillDataSet(0,0,1,JetArgSet, mcSet_inclusive_ref, 0.5,1.0,2, h1_weight, h1_track_ratio, h1_notrack_ratio);
    if(!m_doPP)
        FillDataSet(0,0,1,JetArgSet, mcSet_inclusive_ref, 0.5,1.0,3, h1_weight, h1_track_ratio, h1_notrack_ratio);
 //   mcSet_inclusive->Print();

    RooDataSet *dataSet = new RooDataSet("dataSet","dataSet",JetArgSet,"weight");
    RooDataSet *dataSet_EM = new RooDataSet("dataSet_EM","dataSet_EM",JetArgSet,"weight");
    RooDataSet *dataSet_Minus_EM = new RooDataSet("dataSet_Minus_EM","dataSet_Minus_EM",JetArgSet,"weight");
    FillDataSet(1,0,1, JetArgSet, dataSet, 0.0,1.0,1, h1_weight, h1_track_ratio, h1_notrack_ratio);
    FillDataSet(1,1,1, JetArgSet, dataSet_EM, 0.0,1.0,1, h1_weight, h1_track_ratio, h1_notrack_ratio);
    FillDataSet(1,0,1, JetArgSet, dataSet_Minus_EM, 0.0,1.0,3, h1_weight, h1_track_ratio, h1_notrack_ratio);
    FillDataSet(1,1,1, JetArgSet, dataSet_Minus_EM, 0.0,1.0,3, h1_weight, h1_track_ratio, h1_notrack_ratio);
///////////////////////////////////Check negative value///////////////////////////////////////////////
    //CheckNegative(dataSet_Minus_EM);
    //TH1F* hptrel= (TH1F*)dataSet_Minus_EM->createHistogram("Ptrel",Ptrel,Binning(25,0,4.0));

    //RooDataSet *dataSet = new RooDataSet("dataSet","dataSet",JetTree,JetArgSet_data);
  //  RooDataSet * dataSet = (RooDataSet*)CombineMCJ(1,JetArgSet_data,1);
   // RooDataSet *mcSet_muon= new RooDataSet("mcSet_muon","mcSet_muon",MCTree_muon,JetArgSet_data);
   // RooDataSet *mcSet_inclusive= (RooDataSet*)CombineMCJ(1,JetArgSet_comb,0,0);

    stringstream data_cut_name;
    data_cut_name.str();
    if(m_doComposite && !m_doIPSignificance)
        data_cut_name<<"Muon_asso_composite < 0.25 ";

    else if(!m_doComposite && m_doIPSignificance)
        data_cut_name<<"Muon_asso_sig > 2.0";
        //data_cut_name<<"Jet_weight> 3.0";
    
    else if(m_doComposite && m_doIPSignificance)
        data_cut_name<<"Muon_asso_sig > 2.0 && Muon_asso_composite < 0.25";
    else
        data_cut_name<<"Jet_pt >= 20.0 && Muon_asso_composite < 0.5  ";
    //    data_cut_name<<"Jet_pt >= 20.0 && Muon_asso_sig > 2.0 &&(Jet_pt < 60 || (Jet_pt > 60 &&Jet_pt < 75 && Muon_asso_pt > 6)|| (Jet_pt > 75 && Jet_pt < 90 && Muon_asso_pt > 6 )||( Jet_pt > 90 && Muon_asso_pt >8))";
    //    data_cut_name<<"Jet_pt >= 20.0 && (Jet_pt < 60 || (Jet_pt > 60 &&Jet_pt < 75 && Muon_asso_pt > 6)|| (Jet_pt > 75 && Jet_pt < 90 && Muon_asso_pt > 6 )||( Jet_pt > 90 && Muon_asso_pt >8))";
       //data_cut_name<<"Jet_pt >= 20.0 && Muon_asso_sig > 2.0";

    RooDataSet* dataSet_comp_sig = (RooDataSet*)dataSet->reduce(data_cut_name.str().c_str());
    RooDataSet* dataSet_EMcomp_sig = (RooDataSet*)dataSet_EM->reduce(data_cut_name.str().c_str());
    RooDataSet* dataSet_Minus_EM_comp_sig = (RooDataSet*)dataSet_Minus_EM->reduce(data_cut_name.str().c_str());
    RooDataSet* mcSet_inclusive_comp_sig = (RooDataSet*)mcSet_inclusive->reduce(data_cut_name.str().c_str());
    RooDataSet* mcSet_inclusive_J1_comp_sig = (RooDataSet*)mcSet_inclusive_J1->reduce(data_cut_name.str().c_str());
    RooDataSet* mcSet_inclusive_J2_comp_sig = (RooDataSet*)mcSet_inclusive_J2->reduce(data_cut_name.str().c_str());
    RooDataSet* mcSet_inclusive_J3_comp_sig = (RooDataSet*)mcSet_inclusive_J3->reduce(data_cut_name.str().c_str());
    RooDataSet* mcSet_inclusive_ref_comp_sig = (RooDataSet*)mcSet_inclusive_ref->reduce(data_cut_name.str().c_str());
    RooDataSet* mcSet_muon_comp_sig = (RooDataSet*)mcSet_muon->reduce(data_cut_name.str().c_str());

//    RooDataSet* dataSet_bkg= (RooDataSet*)dataSet->reduce("Muon_asso_composite >= 0.25");

    nPtBins = m_ptBins.size()-1;
    nMuonPtBins = m_muonptBins.size()-1;
    nCentBins = m_centralityBins.size()-1;
    float EMRatio[nPtBins][nCentBins];

    std::stringstream pp;

    RooDataSet* dataSubSet[nPtBins][nCentBins];
    RooDataSet* dataSubSet_comp_sig[nPtBins][nCentBins];
    RooDataSet* dataSubSet_Minus_EM_comp_sig[nPtBins][nCentBins];
    RooDataSet* dataSubSet_EM[nPtBins][nCentBins];
    RooDataSet* dataSubSet_EMcomp_sig[nPtBins][nCentBins];


    RooDataSet* mcSubSet_muon[nPtBins][nCentBins];
    RooDataSet* mcSubSet_muon_comp_sig[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_comp_sig[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_J1_comp_sig[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_J2_comp_sig[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_J3_comp_sig[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_ref_comp_sig[nPtBins][nCentBins];

    RooDataSet* mcSubSet_muon_comp_sig_B[nPtBins][nCentBins];
    RooDataSet* mcSubSet_muon_comp_sig_C[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_comp_sig_L[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_comp_sig_B[nPtBins][nCentBins];
    RooDataSet* mcSubSet_inclusive_comp_sig_C[nPtBins][nCentBins];
    RooDataSet* mcSubSet_muon_comp_sig_allcent[nPtBins];
    RooDataSet* mcSubSet_inclusive_comp_sig_allcent[nPtBins];
    RooDataSet* dataSubSet_comp_sig_allcent[nPtBins];

    //    ofstream out("fraction.txt");
    ofstream out;
    if(m_doMuonPt < -0.5)
    {
        if(!m_doPP)
            out.open("fraction.txt");
        else
            out.open("fraction_pp.txt");
    }
    else
    {
        if(!m_doPP)
            out.open("fraction_mpt_centralR4.txt");
        else 
            out.open("fraction_mpt_central2_pp.txt");
    }
    //out.open("fraction_error.txt", ios::out | ios::app);

    for(int i = 0; i< nPtBins; i++)
    {
        pp.str(""); 
        pp<<"Jet_pt < "<<m_ptBins[i+1] << "&& Jet_pt >"<<m_ptBins[i];

        mcSubSet_inclusive_comp_sig_allcent[i] = (RooDataSet*)mcSet_inclusive_comp_sig->reduce(pp.str().c_str());
        mcSubSet_muon_comp_sig_allcent[i] = (RooDataSet*)mcSet_muon_comp_sig->reduce(pp.str().c_str());
        dataSubSet_comp_sig_allcent[i] = (RooDataSet*)dataSet_comp_sig->reduce(pp.str().c_str());

        for(int k = 0; k < nCentBins; k++)
        {
            dataSubSet[i][k] = SelectPtEtaCent(dataSet,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            dataSubSet_comp_sig[i][k] = SelectPtEtaCent(dataSet_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            dataSubSet_Minus_EM_comp_sig[i][k] = SelectPtEtaCent(dataSet_Minus_EM_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            dataSubSet_EM[i][k] = SelectPtEtaCent(dataSet_EM,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            dataSubSet_EMcomp_sig[i][k] = SelectPtEtaCent(dataSet_EMcomp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);


            mcSubSet_muon[i][k] = SelectPtEtaCent(mcSet_muon,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_muon_comp_sig[i][k] = SelectPtEtaCent(mcSet_muon_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_inclusive[i][k] = SelectPtEtaCent(mcSet_inclusive,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_inclusive_comp_sig[i][k] = SelectPtEtaCent(mcSet_inclusive_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_inclusive_J1_comp_sig[i][k] = SelectPtEtaCent(mcSet_inclusive_J1_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_inclusive_J2_comp_sig[i][k] = SelectPtEtaCent(mcSet_inclusive_J2_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_inclusive_J3_comp_sig[i][k] = SelectPtEtaCent(mcSet_inclusive_J3_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);
            mcSubSet_inclusive_ref_comp_sig[i][k] = SelectPtEtaCent(mcSet_inclusive_ref_comp_sig,m_ptBins[i], m_ptBins[i+1],m_centralityBins[k],m_centralityBins[k+1]);

            mcSubSet_inclusive_comp_sig_B[i][k] = (RooDataSet*)mcSubSet_inclusive_comp_sig[i][k]->reduce("Flavour_tag ==1 || Flavour_tag == 7 ");
            mcSubSet_inclusive_comp_sig_C[i][k] = (RooDataSet*)mcSubSet_inclusive_comp_sig[i][k]->reduce("Flavour_tag ==2");
            //mcSubSet_inclusive_comp_sig_C[i][k] = (RooDataSet*)mcSubSet_inclusive_comp_sig[i][k]->reduce("Flavour_tag ==2");
            mcSubSet_inclusive_comp_sig_L[i][k] = (RooDataSet*)mcSubSet_inclusive_comp_sig[i][k]->reduce("Flavour_tag ==3|| Flavour_tag == 4 || Flavour_tag == 5 || Flavour_tag ==6 ");

            if(m_doMuonPt == i)
            {
                RooDataSet* dataSubSet_diff_mpt[nMuonPtBins];
                RooDataSet* dataSubSet_EM_diff_mpt[nMuonPtBins];

                RooDataSet* mcSubSet_muon_diff_mpt[nMuonPtBins];
                RooDataSet* mcSubSet_inclusive_diff_mpt[nMuonPtBins];
                RooDataSet* mcSubSet_inclusive_signal_diff_mpt[nMuonPtBins];
                RooDataSet* mcSubSet_muon_signal_diff_mpt[nMuonPtBins];

                for(int j = 0; j < nMuonPtBins; j++)
                {
                    dataSubSet_diff_mpt[j] = SelectMuonPtEtaCent(dataSubSet_comp_sig[i][k],m_muonptBins[j],m_muonptBins[j+1]);
                    dataSubSet_EM_diff_mpt[j] = SelectMuonPtEtaCent(dataSubSet_EMcomp_sig[i][k],m_muonptBins[j],m_muonptBins[j+1]);
                    mcSubSet_muon_diff_mpt[j] = SelectMuonPtEtaCent(mcSubSet_muon[i][k],m_muonptBins[j], m_muonptBins[j+1]);
                    mcSubSet_inclusive_diff_mpt[j] = SelectMuonPtEtaCent(mcSubSet_inclusive[i][k],m_muonptBins[j], m_muonptBins[j+1]);

                    mcSubSet_inclusive_signal_diff_mpt[j] =(RooDataSet*) mcSubSet_inclusive_diff_mpt[j]->reduce("Flavour_tag == 1 && Ptrel < 5.0");

                    EMRatio[j][k] = (float)dataSubSet_EM_diff_mpt[j]->numEntries()/5/dataSubSet_diff_mpt[j]->numEntries();
                    if(m_doBuildTemplate)
                        BuildTemplate(mcSubSet_muon_diff_mpt[j], mcSubSet_inclusive_diff_mpt[j], dataSubSet_EM_diff_mpt[j],var, j,k);

                    FitDataSet(dataSubSet_diff_mpt[j], j, k, file_FitPlot,EMRatio[j][k],out);
                    //FitDataSet(mcSubSet_inclusive_diff_mpt[j], j, k, file_FitPlot, 0, out);
                }
            }

    //////////////////////////////////////////Fit data in different jet pt bins//////////////////////////////////////////////

            if(m_doMuonPt < -0.5)
            {
                EMRatio[i][k] = (float)dataSubSet_EMcomp_sig[i][k]->numEntries()/5/dataSubSet_comp_sig[i][k]->numEntries();
                std::cout<<"muon jet number in different jpt: "<<i<<std::endl;
                std::cout<<dataSubSet_comp_sig[i][k]->numEntries()<<"number of Data jet"<<std::endl; 

                if(m_doBuildTemplate)
                {
                //    BuildTemplate(mcSubSet_muon_comp_sig[i][k], mcSubSet_inclusive_comp_sig[i][k], dataSubSet_EMcomp_sig[i][k],var, i,k);
            //        BuildTemplate(mcSubSet_muon_comp_sig[i][k], dataSubSet_Minus_EM_comp_sig[i][k], dataSubSet_EMcomp_sig[i][k],var, i,k);
                    //BuildTemplate_MC(mcSubSet_inclusive_J1_comp_sig[i][k], mcSubSet_inclusive_J2_comp_sig[i][k], mcSubSet_inclusive_J3_comp_sig[i][k],var, i,k);
                    BuildTemplate(mcSubSet_muon_comp_sig_allcent[i], mcSubSet_inclusive_comp_sig_allcent[i], dataSubSet_EMcomp_sig[i][k],var, i,k);
                }
                FitDataSet(dataSubSet_comp_sig[i][k], i,k, file_FitPlot,EMRatio[i][k], out);
                //FitDataSet(mcSubSet_inclusive_comp_sig[i][k], i,k, file_FitPlot, 0, out);
            }
    /*
            Roo2DKeysPdf tmpPdf("tmpPdf","tmpPdf",var,var_2D_comp,*mcSubSet_inclusive_comp_sig[i]);
            TH1 * hh_pdf = tmpPdf.createHistogram("hh_model",var,Binning(50),YVar(var_2D_comp,Binning(50))); 
            hh_pdf->Draw();

            RooFitResult* tmpResult = tmpPdf.fitTo(*dataSubSet[i], Minos(true), SumW2Error(true), Save(true), PrintLevel(-1),Warnings(false),PrintEvalErrors(0));
            tmpResult->Print("V"); 

    */

        }
    }
    out.close();

}
  
  
void TempFitter::FillDataSet(bool do_data, bool do_eventmix, bool do_inclusive, RooArgSet& JetArgSetMC , RooDataSet * dataSet,double startF, double endF, int J_sample, TH1F* h1_weight, TH1F* h1_track_ratio, TH1F* h1_notrack_ratio)
{
    double cross_sections[3];
    double filter_efficiency_inclusive[3];
    int event_Number_inclusive ;
    float jpt_cut[3];
    double cross_sections_AA[3]={1.8760E+05,8.2787E+03,2.9419E+02};
    double cross_sections_PP[3]={4.13E+07,3.2262E+03,0.0};
    double filter_efficiency_inclusive_AA[3]={1.0,1.0,1.0};
    double filter_efficiency_inclusive_PP[3]={1.5641E-03,3.6450E-02,0.0};
    float jpt_cut_AA[3] = {45.5,91,182};
    float jpt_cut_PP[3] = {65,200,300};
    if(!m_doPP)
    {
        event_Number_inclusive = 4500000;
        for(int J = 0; J < 3; J++)
        {
            cross_sections[J] = cross_sections_AA[J];
            filter_efficiency_inclusive[J] =filter_efficiency_inclusive_AA[J]; 
            jpt_cut[J] = jpt_cut_AA[J];
        }
    }
    else
    {
        event_Number_inclusive = 5000000;
        for(int J = 0; J < 3; J++)
        {
            cross_sections[J] = cross_sections_PP[J];
            filter_efficiency_inclusive[J] =filter_efficiency_inclusive_PP[J]; 
            jpt_cut[J] = jpt_cut_PP[J];
        }
    }
    //double filter_efficiency_inclusive_cor[3]={0.01500,0.03480,0.06732};
    //int event_Number_inclusive = 1274000;
    double filter_efficiency_muon[3]={6.6523E-03,1.4941E-02,2.4284E-02};
    //double filter_efficiency_muon_cor[3]={1.7745E-03,8.4775E-03,1.7822E-02};
    int event_Number_muon= 650000;
    double scale_factor[3];
 
    float Jet_eta;
    float Jet_pt;
    float Ptrel;
    float Ptrel_truth;
    int   Centrality_hp;
    int   Flavour_tag;
    float Muon_asso_composite;
    float Muon_asso_sig;
    float Muon_asso_pt;
    float Fcal;
    float Jet_Trackjet_dR;
    float Muon_Jet_dR;
    float Jet_weight;

    std::stringstream ss; 
    std::stringstream s; 

    TFile* files[3];
    TTree* JetTree[3];
    RooDataSet* mcSet[3];
    RooDataSet *mcSetComb ;

    RooRealVar* w[3] ;
    RooDataSet *mcSetWeight [3];

    for(int J=1; J<=3; J++)
    {   
        scale_factor[J-1]=cross_sections[J-1];
        if(do_inclusive) scale_factor[J-1]*=filter_efficiency_inclusive[J-1]/event_Number_inclusive;
        else scale_factor[J-1]*=filter_efficiency_muon[J-1]/event_Number_muon;
        //scale_factor[J-1] = 1.0;
    }    

//    scale_factor[J_sample-1] = (float)scale_factor[J_sample-1]/scale_factor[2];
    //std::cout<<scale_factor[J_sample-1]<<std::endl;
    ss.str("");

    if(do_inclusive && !do_data)
    {
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_test_inclusive_notrigoverlay_J" << J_sample << ".root";
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_track_inclusive_recopt_J" << J_sample << ".root";
        if(!m_doPP)
            ss << "/atlas/data4/tingting/data6/ReconCode/histogram_inclusive_recopt_notrkcor_J" << J_sample << ".root";
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_inclusive_recopt_notrkcor_smear_J" << J_sample << ".root";
        else
            ss << "/atlas/data4/tingting/data6/ReconCode/histogram_PP_inclusive_new_JZ" << J_sample << ".root";
    }
    else if(!do_inclusive && !do_data)
//        ss << "/atlas/data4/tingting/data6/ReconCode/histogram_test_muon_cal_trkcor_J" << J_sample << ".root";
      //  ss << "/atlas/data4/tingting/data6/ReconCode/histogram_track_muon_D3PD_loosetruth_J" << J_sample << ".root";
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_track_muon_recopt_J" << J_sample << ".root";
        ss << "/atlas/data4/tingting/data6/ReconCode/histogram_muon_recopt_notrkcor_J" << J_sample << ".root";
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_muon_recopt_notrkcor_smear_J" << J_sample << ".root";
    else if(do_data)
    {
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_Data11_AA_minbia_randm_aveweight_addmuon_sig_nodyn_nummuon.root"; 
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_Data11_AA_minbia_randm_aveweight_addmuon_sig_trackmuon.root"; 
        //ss << "/atlas/data4/tingting/data6/ReconCode/histogram_Data11_AA_minbia_randm_aveweight_addmuon_sig_weight_cor.root"; 
        if(!m_doPP)
            ss << "/atlas/data4/tingting/data6/ReconCode/histogram_Data11_AA_minbia_randm_aveweight_addmuon_sig_notrkcor.root"; 
        else
            ss << "/atlas/data4/tingting/data6/ReconCode/histogram_Data13_pp_cor.root"; 
    }

    files[J_sample-1]=TFile::Open(ss.str().c_str());

    if(!do_data)
        JetTree[J_sample-1] = (TTree*)files[J_sample-1]->Get("MCJetTree");
    else if(!do_eventmix)
        JetTree[J_sample-1] = (TTree*)files[J_sample-1]->Get("DataJetTree");
    else 
        JetTree[J_sample-1] = (TTree*)files[J_sample-1]->Get("DataJetTree_em");

    JetTree[J_sample-1]->SetBranchAddress("Jet_eta",&Jet_eta);
    JetTree[J_sample-1]->SetBranchAddress("Jet_pt",&Jet_pt);
    JetTree[J_sample-1]->SetBranchAddress("Centrality_hp",&Centrality_hp);
    JetTree[J_sample-1]->SetBranchAddress("Flavour_tag",&Flavour_tag);
    JetTree[J_sample-1]->SetBranchAddress("Ptrel",&Ptrel);
    JetTree[J_sample-1]->SetBranchAddress("Muon_asso_composite",&Muon_asso_composite);
    JetTree[J_sample-1]->SetBranchAddress("Muon_asso_sig",&Muon_asso_sig);
    JetTree[J_sample-1]->SetBranchAddress("Muon_asso_pt",&Muon_asso_pt);
    JetTree[J_sample-1]->SetBranchAddress("Fcal",&Fcal);
    JetTree[J_sample-1]->SetBranchAddress("Jet_Trackjet_dR",&Jet_Trackjet_dR);
    JetTree[J_sample-1]->SetBranchAddress("Muon_Jet_dR",&Muon_Jet_dR);
    JetTree[J_sample-1]->SetBranchAddress("Jet_weight",&Jet_weight);

    int start = (int)JetTree[J_sample-1]->GetEntries()*startF;
    int end = (int)JetTree[J_sample-1]->GetEntries()*endF;
    //int start = 18000;
    //int end = 230000;

    if(endF <= 1.0)
    {
        for ( int i = start; i < end; i++ ) 
        {
            JetTree[J_sample-1]->GetEntry(i);
          //  if(J_sample == 3 &&(i == 23529 || i == 9879 ||i == 72761|| i == 221910))
            if(isnan(Ptrel))
            {
                std::cout<<"Event :"<<i<<"Ptrel: "<<Ptrel<<std::endl;
                continue;
            }

        //    JetArgSetMC.setRealValue("Jet_eta",Jet_eta);
            JetArgSetMC.setRealValue("Jet_pt",Jet_pt);
            JetArgSetMC.setRealValue("Flavour_tag",Flavour_tag);
            JetArgSetMC.setRealValue("Centrality_hp",Centrality_hp);
            JetArgSetMC.setRealValue("Muon_asso_composite",Muon_asso_composite);
            JetArgSetMC.setRealValue("Muon_asso_sig",Muon_asso_sig);
            JetArgSetMC.setRealValue("Muon_asso_pt",Muon_asso_pt);
            //JetArgSetMC.setRealValue("Jet_Trackjet_dR",Jet_Trackjet_dR);
            JetArgSetMC.setRealValue("Muon_Jet_dR",Muon_Jet_dR);
            JetArgSetMC.setRealValue("Jet_weight",Jet_weight);
            JetArgSetMC.setRealValue("Ptrel",Ptrel);
            //else
            //JetArgSetMC.setRealValue("Ptrel",Ptrel_truth);

            float fcal_weight = h1_weight->Interpolate(Fcal);
            float ftrack_ratio = h1_track_ratio->Interpolate(Jet_pt);
            float fnotrack_ratio = h1_notrack_ratio->Interpolate(Jet_pt);

            if(do_data&& do_eventmix && J_sample == 3)
            {
                dataSet->add(JetArgSetMC, -0.20);
            }

            else if(do_data)
                dataSet->add(JetArgSetMC, 1.0);
            else if(Jet_pt < jpt_cut[J_sample-1]) 
            {
                dataSet->add(JetArgSetMC, scale_factor[J_sample-1]);
         //       dataSet->add(JetArgSetMC, 1.0);
    /*
                if(Jet_Trackjet_dR < 0.2)
                    dataSet->add(JetArgSetMC, scale_factor[J_sample-1]* fcal_weight * ftrack_ratio);
                else
                    dataSet->add(JetArgSetMC, scale_factor[J_sample-1]* fcal_weight * fnotrack_ratio);
    */
            }
        }
    }
    else
    {
        start = (int)JetTree[J_sample-1]->GetEntries()*startF;
        end = (int)JetTree[J_sample-1]->GetEntries()*1.0;
        for ( int i = start; i < end; i++ ) 
        {
            JetTree[J_sample-1]->GetEntry(i);
          //  if(J_sample == 3 &&(i == 23529 || i == 9879 ||i == 72761|| i == 221910))
            if(isnan(Ptrel))
            {
                std::cout<<"Event :"<<i<<"Ptrel: "<<Ptrel<<std::endl;
                continue;
            }

        //    JetArgSetMC.setRealValue("Jet_eta",Jet_eta);
            JetArgSetMC.setRealValue("Jet_pt",Jet_pt);
            JetArgSetMC.setRealValue("Flavour_tag",Flavour_tag);
            JetArgSetMC.setRealValue("Centrality_hp",Centrality_hp);
            JetArgSetMC.setRealValue("Ptrel",Ptrel);
            JetArgSetMC.setRealValue("Muon_asso_composite",Muon_asso_composite);
            JetArgSetMC.setRealValue("Muon_asso_sig",Muon_asso_sig);
            JetArgSetMC.setRealValue("Muon_asso_pt",Muon_asso_pt);
            JetArgSetMC.setRealValue("Jet_Trackjet_dR",Jet_Trackjet_dR);
//            JetArgSetMC.setRealValue("Muon_Jet_dR",Muon_Jet_dR);


            float fcal_weight = h1_weight->Interpolate(Fcal);
            float ftrack_ratio = h1_track_ratio->Interpolate(Jet_pt);
            float fnotrack_ratio = h1_notrack_ratio->Interpolate(Jet_pt);

            if(do_data)
                dataSet->add(JetArgSetMC, 1.0);
            else if(Jet_pt < jpt_cut[J_sample-1]) 
    //        else
            {
                dataSet->add(JetArgSetMC, scale_factor[J_sample-1]);
         //       dataSet->add(JetArgSetMC, 1.0);
    /*
                if(Jet_Trackjet_dR < 0.2)
                    dataSet->add(JetArgSetMC, scale_factor[J_sample-1]* fcal_weight * ftrack_ratio);
                else
                    dataSet->add(JetArgSetMC, scale_factor[J_sample-1]* fcal_weight * fnotrack_ratio);
    */
            }
        }

        start = (int)JetTree[J_sample-1]->GetEntries()*0.0;
        end = (int)JetTree[J_sample-1]->GetEntries()*(endF-1.0);
        for ( int i = start; i < end; i++ ) 
        {
            JetTree[J_sample-1]->GetEntry(i);
          //  if(J_sample == 3 &&(i == 23529 || i == 9879 ||i == 72761|| i == 221910))
            if(isnan(Ptrel))
            {
                std::cout<<"Event :"<<i<<"Ptrel: "<<Ptrel<<std::endl;
                continue;
            }

        //    JetArgSetMC.setRealValue("Jet_eta",Jet_eta);
            JetArgSetMC.setRealValue("Jet_pt",Jet_pt);
            JetArgSetMC.setRealValue("Flavour_tag",Flavour_tag);
            JetArgSetMC.setRealValue("Centrality_hp",Centrality_hp);
            JetArgSetMC.setRealValue("Ptrel",Ptrel);
            JetArgSetMC.setRealValue("Muon_asso_composite",Muon_asso_composite);
            JetArgSetMC.setRealValue("Muon_asso_sig",Muon_asso_sig);
            JetArgSetMC.setRealValue("Muon_asso_pt",Muon_asso_pt);
            JetArgSetMC.setRealValue("Jet_Trackjet_dR",Jet_Trackjet_dR);
           // JetArgSetMC.setRealValue("Muon_Jet_dR",Muon_Jet_dR);


            float fcal_weight = h1_weight->Interpolate(Fcal);
            float ftrack_ratio = h1_track_ratio->Interpolate(Jet_pt);
            float fnotrack_ratio = h1_notrack_ratio->Interpolate(Jet_pt);

            if(do_data)
                dataSet->add(JetArgSetMC, 1.0);
            else if(Jet_pt < jpt_cut[J_sample-1]) 
                dataSet->add(JetArgSetMC, scale_factor[J_sample-1]);
        }
    }
    if(do_data)
        dataSet->Print();

}

void TempFitter::BuildTemplate(RooDataSet* mcSubSet_muon, RooDataSet* mcSubSet_inclusive, RooDataSet* dataSubSet_EM,RooRealVar& variable ,int PtBins, int CentBins)
{

    RooDataSet *mcSubSetSigBUnbinned = (RooDataSet*)mcSubSet_muon->reduce("Flavour_tag == 1 || Flavour_tag ==7");
    RooDataSet *mcSubSetBkgCUnbinned= (RooDataSet*)mcSubSet_muon->reduce("Flavour_tag == 2");
    RooDataSet *mcSubSetBkgLUnbinned= (RooDataSet*)mcSubSet_inclusive->reduce(" Flavour_tag == 3 || Flavour_tag == 4 || Flavour_tag == 5 || Flavour_tag == 6");
  //  RooDataSet *mcSubSetBkgLUnbinned= (RooDataSet*)mcSubSet_inclusive->reduce(" Muon_asso_composite > 0.25");
    RooDataSet *mcSubSetBkgK= (RooDataSet*)mcSubSet_inclusive->reduce("Flavour_tag == 3 ");
    RooDataSet *mcSubSetBkgP= (RooDataSet*)mcSubSet_inclusive->reduce("Flavour_tag == 4 ");
    RooDataSet *mcSubSetBkgOther= (RooDataSet*)mcSubSet_inclusive->reduce("Flavour_tag == 5 || Flavour_tag ==6");

    variable.setBins(40); 
    RooDataHist* mcSubSetSigB= new RooDataHist("mcSubSetSigB","binned data",RooArgSet(variable),*mcSubSetSigBUnbinned);
    RooDataHist* mcSubSetBkgC= new RooDataHist("mcSubSetBkgC","binned data",RooArgSet(variable),*mcSubSetBkgCUnbinned);
    RooDataHist* mcSubSetBkgL= new RooDataHist("mcSubSetBkgL","binned data",RooArgSet(variable),*mcSubSetBkgLUnbinned);
    RooDataHist* dataSubSet_EMBinned = new RooDataHist("dataSubSet_EMBinned","binned data",RooArgSet(variable),*dataSubSet_EM);
  // --- Number of prompts, pions, kaons, others ---
    RooRealVar nSignal("nSignal","nSignal", mcSubSetSigB->numEntries() );
    RooRealVar nCQuarks("nCQuarks","nCQuarks", mcSubSetBkgC->numEntries() );
    RooRealVar nLQuarks("nLQuarks","nLQuarks", mcSubSetBkgL->numEntries() );
    RooRealVar nPions("nPions","nPions",    mcSubSetBkgP->numEntries() );
    RooRealVar nKaons("nKaons","nKaons",    mcSubSetBkgK->numEntries() );
    RooRealVar nOthers("nOthers","nOthers", mcSubSetBkgOther->numEntries() );
 
// build template    

    RooHistPdf mcSigBPdf("mcSigBPdf","mcSigBPdf",variable,*mcSubSetSigB);
    RooHistPdf mcBkgCPdf("mcBkgCPdf","mcBkgCPdf",variable,*mcSubSetBkgC);
    RooHistPdf mcBkgLPdf("mcBkgLPdf","mcBkgLPdf",variable,*mcSubSetBkgL);
    RooHistPdf dataEMPdf("dataEMPdf","dataEMPdf",variable,*dataSubSet_EMBinned);
   // RooKeysPdf mcSigBPdf("mcSigBPdf","mcSigBPdf",variable,*mcSubSetSigB);
   // RooKeysPdf mcBkgCPdf("mcBkgCPdf","mcBkgCPdf",variable,*mcSubSetBkgC);
   // RooKeysPdf mcBkgLPdf("mcBkgLPdf","mcBkgLPdf",variable,*mcSubSetBkgL);
    RooKeysPdf mcBkgKPdf("mcBkgKPdf","mcBkgKPdf",variable,*mcSubSetBkgK);
    RooKeysPdf mcBkgPPdf("mcBkgPPdf","mcBkgPPdf",variable,*mcSubSetBkgP);
    RooKeysPdf mcBkgOtherPdf("mcBkgOtherPdf","mcBkgOtherPdf",variable,*mcSubSetBkgOther);

    std::stringstream filename;
    filename.str("");
    if(m_doMuonPt > -0.5)
        filename<<"MCTemplateResult_Pt_diff_muonpt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
    else 
        filename<<"MCTemplateResult_Pt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
    RooWorkspace *w = new RooWorkspace("w","workspace");
    w->import(variable);
 
    // --- And write to file ---
    w->import(mcSigBPdf);
    w->import(mcBkgCPdf);
    w->import(mcBkgLPdf);
    w->import(mcBkgKPdf);
    w->import(mcBkgPPdf);
    w->import(mcBkgOtherPdf);
    w->import(dataEMPdf);
    w->import(nSignal);
    w->import(nPions);
    w->import(nKaons);
    w->import(nCQuarks);
    w->import(nOthers);
    w->writeToFile(filename.str().c_str());
    std::cout << "Done!" << std::endl;
  
    dataEMPdf.Print();
    // --- Clean up ---
    delete mcSubSetSigB;
    delete mcSubSetBkgC;
    delete mcSubSetBkgL;
    delete mcSubSetBkgK;
    delete mcSubSetBkgP;
    delete mcSubSetBkgOther;
    delete w; 


}

void TempFitter::BuildTemplate_MC(RooDataSet* mcSubSet_J1, RooDataSet* mcSubSet_J2, RooDataSet* mcSubSet_J3,RooRealVar& variable ,int PtBins, int CentBins)
{

    RooDataSet *mcSubSetSigB_J1= (RooDataSet*)mcSubSet_J1->reduce("Flavour_tag == 1 || Flavour_tag ==7");
    RooDataSet *mcSubSetBkgC_J1= (RooDataSet*)mcSubSet_J1->reduce("Flavour_tag == 2");
    RooDataSet *mcSubSetBkgL_J1= (RooDataSet*)mcSubSet_J1->reduce(" Flavour_tag == 3 || Flavour_tag == 4 || Flavour_tag == 5 || Flavour_tag == 6");

    RooDataSet *mcSubSetSigB_J2= (RooDataSet*)mcSubSet_J2->reduce("Flavour_tag == 1 || Flavour_tag ==7");
    RooDataSet *mcSubSetBkgC_J2= (RooDataSet*)mcSubSet_J2->reduce("Flavour_tag == 2");
    RooDataSet *mcSubSetBkgL_J2= (RooDataSet*)mcSubSet_J2->reduce(" Flavour_tag == 3 || Flavour_tag == 4 || Flavour_tag == 5 || Flavour_tag == 6");

    RooDataSet *mcSubSetSigB_J3= (RooDataSet*)mcSubSet_J3->reduce("Flavour_tag == 1 || Flavour_tag ==7");
    RooDataSet *mcSubSetBkgC_J3= (RooDataSet*)mcSubSet_J3->reduce("Flavour_tag == 2");
    RooDataSet *mcSubSetBkgL_J3= (RooDataSet*)mcSubSet_J3->reduce(" Flavour_tag == 3 || Flavour_tag == 4 || Flavour_tag == 5 || Flavour_tag == 6");

// build template    

    RooKeysPdf mcSigBJ1Pdf("mcSigBJ1Pdf","mcSigBJ1Pdf",variable,*mcSubSetSigB_J1);
    RooKeysPdf mcBkgCJ1Pdf("mcBkgCJ1Pdf","mcBkgCJ1Pdf",variable,*mcSubSetBkgC_J1);
    RooKeysPdf mcBkgLJ1Pdf("mcBkgLJ1Pdf","mcBkgLJ1Pdf",variable,*mcSubSetBkgL_J1);
    RooKeysPdf mcSigBJ2Pdf("mcSigBJ2Pdf","mcSigBJ2Pdf",variable,*mcSubSetSigB_J2);
    RooKeysPdf mcBkgCJ2Pdf("mcBkgCJ2Pdf","mcBkgCJ2Pdf",variable,*mcSubSetBkgC_J2);
    RooKeysPdf mcBkgLJ2Pdf("mcBkgLJ2Pdf","mcBkgLJ2Pdf",variable,*mcSubSetBkgL_J2);
    RooKeysPdf mcSigBJ3Pdf("mcSigBJ3Pdf","mcSigBJ3Pdf",variable,*mcSubSetSigB_J3);
    RooKeysPdf mcBkgCJ3Pdf("mcBkgCJ3Pdf","mcBkgCJ3Pdf",variable,*mcSubSetBkgC_J3);
    RooKeysPdf mcBkgLJ3Pdf("mcBkgLJ3Pdf","mcBkgLJ3Pdf",variable,*mcSubSetBkgL_J3);

    std::stringstream filename;
    filename.str("");
    if(m_doMuonPt > -0.5)
        filename<<"MCTemplateResult_J_Pt_diff_muonpt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
    else 
        filename<<"MCTemplateResult_J_Pt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
    RooWorkspace *w = new RooWorkspace("w","workspace");
    w->import(variable);
 
    // --- And write to file ---
    w->import(mcSigBJ1Pdf);
    w->import(mcBkgCJ1Pdf);
    w->import(mcBkgLJ1Pdf);
    w->import(mcSigBJ2Pdf);
    w->import(mcBkgCJ2Pdf);
    w->import(mcBkgLJ2Pdf);
    w->import(mcSigBJ3Pdf);
    w->import(mcBkgCJ3Pdf);
    w->import(mcBkgLJ3Pdf);
    w->writeToFile(filename.str().c_str());
    std::cout << "Done!" << std::endl;
  
    // --- Clean up ---
    delete mcSubSetSigB_J1;
    delete mcSubSetBkgC_J1;
    delete mcSubSetBkgL_J1;
    delete mcSubSetSigB_J2;
    delete mcSubSetBkgC_J2;
    delete mcSubSetBkgL_J2;
    delete mcSubSetSigB_J3;
    delete mcSubSetBkgC_J3;
    delete mcSubSetBkgL_J3;
    delete w; 


}
void TempFitter::FitDataSet(RooDataSet* dataSubSet, int PtBins, int CentBins,TFile* fitResult, float ratio, ofstream &out)
{
    std::stringstream filename;
    filename.str("");
    if(m_doMuonPt > -0.5)
        filename<<"MCTemplateResult_Pt_diff_muonpt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
    else 
        //filename<<"MCTemplateResult_J_Pt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
        filename<<"MCTemplateResult_Pt_"<<PtBins<<"_Cent_"<<CentBins<<".root";

    TFile *fMCJPdf =TFile::Open(filename.str().c_str()); 
    RooWorkspace* w_new = (RooWorkspace*) fMCJPdf->Get("w");

    filename.str("");
    if(m_doMuonPt > -0.5)
        filename<<"MCTemplateResult_Pt_diff_muonpt_"<<PtBins<<"_Cent_"<<CentBins<<".root";
    else 
        filename<<"MCTemplateResult_Pt_"<<PtBins<<"_Cent_"<<CentBins<<".root";

    TFile *fMCPdf =TFile::Open(filename.str().c_str()); 
    RooWorkspace* w= (RooWorkspace*) fMCPdf->Get("w");

    RooRealVar* var = (RooRealVar*)w->var("Ptrel");
    RooRealVar& variable = *var;
////////////////////////////////////////////try weight different J sample pdf////////////////////////////////////
/*
    RooKeysPdf* mcSigBJ1Pdf = (RooKeysPdf*)w->pdf("mcSigBJ1Pdf");
    RooKeysPdf* mcBkgCJ1Pdf= (RooKeysPdf*)w->pdf("mcBkgCJ1Pdf");
    RooKeysPdf* mcBkgLJ1Pdf= (RooKeysPdf*)w->pdf("mcBkgLJ1Pdf");
    RooKeysPdf* mcSigBJ2Pdf = (RooKeysPdf*)w->pdf("mcSigBJ2Pdf");
    RooKeysPdf* mcBkgCJ2Pdf= (RooKeysPdf*)w->pdf("mcBkgCJ2Pdf");
    RooKeysPdf* mcBkgLJ2Pdf= (RooKeysPdf*)w->pdf("mcBkgLJ2Pdf");
    RooKeysPdf* mcSigBJ3Pdf = (RooKeysPdf*)w->pdf("mcSigBJ3Pdf");
    RooKeysPdf* mcBkgCJ3Pdf= (RooKeysPdf*)w->pdf("mcBkgCJ3Pdf");
    RooKeysPdf* mcBkgLJ3Pdf= (RooKeysPdf*)w->pdf("mcBkgLJ3Pdf");
*/
////////////////////////////////////////////try not weight different J sample pdf////////////////////////////////////
    RooHistPdf* mcSigBPdf = (RooHistPdf*)w->pdf("mcSigBPdf");
    RooHistPdf* mcBkgCPdf= (RooHistPdf*)w->pdf("mcBkgCPdf");
    RooHistPdf* mcBkgLPdf= (RooHistPdf*)w->pdf("mcBkgLPdf");
   // RooKeysPdf* mcSigBPdf = (RooKeysPdf*)w->pdf("mcSigBPdf");
   // RooKeysPdf* mcBkgCPdf= (RooKeysPdf*)w->pdf("mcBkgCPdf");
   // RooKeysPdf* mcBkgLPdf= (RooKeysPdf*)w->pdf("mcBkgLPdf");
    RooKeysPdf* dataEMPdf = (RooKeysPdf*)w->pdf("dataEMPdf");


//////////////////////////////////////Shift and Stretch////////////////////////////////////////
    RooRealVar shift("shift","shift",0,-0.3,0.3);
    RooRealVar stretch("stretch","stretch",1.0,+0.8,+1.3);
    RooShiftedKeysPdf *sigBShiftPdfp, *bkgCShiftPdfp,*bkgLShiftPdfp;

    double cross_sections[3]={1.8760E+05,8.2787E+03,2.9419E+02};
    double filter_efficiency_inclusive[3]={1.0,1.0,1.0};
    int event_Number_inclusive = 4500000;
    int event_Number_muon= 650000;
    double filter_efficiency_muon_gen[3]={6.6523E-03,1.4941E-02,2.4284E-02};
    double scale_factor[3];
    double weight_value[2];
    for(int J=1; J<=3; J++)
    {   
        scale_factor[J-1]=cross_sections[J-1];
        scale_factor[J-1]*=filter_efficiency_inclusive[J-1]/event_Number_inclusive;
    }
/*
    if(PtBins < 2)
    {
        weight_value[0]= scale_factor[0]/(scale_factor[0]+scale_factor[1]+scale_factor[2]);
        weight_value[1]= scale_factor[1]/(scale_factor[0]+scale_factor[1]+scale_factor[2]);
        //RooRealVar weight1("weight1","weight1",scale_factor[0]/(scale_factor[0]+scale_factor[1]+scale_factor[2]));
        //RooRealVar weight2("weight2","weight2",scale_factor[1]/(scale_factor[0]+scale_factor[1]+scale_factor[2]));
        //RooRealVar weight3("weight3","weight3",scale_factor[2]);

    }
    else if (PtBins < 3)
    {
        weight_value[0]= 0.0;
        weight_value[1]=scale_factor[1]/(scale_factor[1]+scale_factor[2]); 
    }
    else 
    {
        weight_value[0]= 0.0;
        weight_value[1]= 0.0;
    }

    RooRealVar weight1("weight1","weight1",weight_value[0]);
    RooRealVar weight2("weight2","weight2",weight_value[1]);

    RooAddPdf *mcSigBPdf= new RooAddPdf("mcSigBPdf","sig J1+ J2+J3", RooArgList(*mcSigBJ1Pdf,*mcSigBJ2Pdf,*mcSigBJ3Pdf), RooArgList(weight1, weight2)); 
    RooAddPdf *mcBkgCPdf= new RooAddPdf("mcBkgCPdf","bkgC J1+ J2+J3", RooArgList(*mcBkgCJ1Pdf,*mcBkgCJ2Pdf,*mcBkgCJ3Pdf), RooArgList(weight1, weight2)); 
    RooAddPdf *mcBkgLPdf= new RooAddPdf("mcBkgLPdf","bkgL J1+ J2+J3", RooArgList(*mcBkgLJ1Pdf,*mcBkgLJ2Pdf,*mcBkgLJ3Pdf), RooArgList(weight1, weight2)); 
*/
    RooAbsPdf& sigBPdf = *mcSigBPdf;
    RooAbsPdf& bkgCPdf = *mcBkgCPdf;
    RooAbsPdf& bkgLPdf = *mcBkgLPdf;

    sigBShiftPdfp = new RooShiftedKeysPdf("sigBShiftPdf","sigBShiftPdf",variable,sigBPdf,shift,stretch);
    bkgCShiftPdfp = new RooShiftedKeysPdf("bkgCShiftPdf","bkgCShiftPdf",variable,bkgCPdf,shift,stretch);
    bkgLShiftPdfp = new RooShiftedKeysPdf("bkgLShiftPdf","bkgLShiftPdf",variable,bkgLPdf,shift,stretch);

//////////////////////////////////////Smear/////////////////////////////////////////////////////
    RooRealVar gaussMean("gaussMean","gaussMean",0.);
    RooRealVar gaussSigma("gaussSigma","gaussSigma",0.01,0.0001,1.);
  // gaussSigma.setConstant(true);
    RooGaussian gauss("gauss","gauss",variable,gaussMean,gaussSigma);

  
  //RooFFTConvPdf *sigConvPdf, *bkgConvPdf;
    // --- Smear PDF's ---
    RooAbsPdf *sigConvPdfp = new RooFFTConvPdf("sigConvPdfp","sigConvPdfp",variable,*sigBShiftPdfp,gauss); //with shift, stretch, and gaussian convolution
    RooAbsPdf *bkgCConvPdfp = new RooFFTConvPdf("bkgCConvPdfp","bkgCConvPdfp",variable,*bkgCShiftPdfp,gauss);
    RooAbsPdf *bkgLConvPdfp = new RooFFTConvPdf("bkgLConvPdfp","bkgLConvPdfp",variable,*bkgLShiftPdfp,gauss);

/*
    RooAbsPdf & sigConvPdf=*sigConvPdfp;
    RooAbsPdf & bkgCConvPdf=*bkgCConvPdfp;
    RooAbsPdf & bkgLConvPdf=*bkgLConvPdfp;

   
*/
    double OriginalEMRatio[4][6] = {{0.312587,0.174912,0.114999,0.0290559,0.0116935,0.00298507},
                                    {0.328693,0.19649,0.134754,0.0393399,0.01875,0.00526316},
                                    {0.163333,0.110857,0.0733747,0.0204152,0.0181818,0},
                                    {0.169378,0.114949,0.0770766,0.0236786,0.0209302,0}};
    int column;           //0: without cut; 1: composite cut; 2:significance cut; 3: both comp and sig cut.

    if(m_doComposite && !m_doIPSignificance)
        column = 1;

    else if(!m_doComposite && m_doIPSignificance)
        column = 2;
    
    else if(m_doComposite && m_doIPSignificance)
        column = 3;
    else
        column = 0;

    mcBkgCPdf->Print("V");

    FILE*  CLratio;
    if(m_doMuonPt < -0.5)
    {
        if(!m_doPP)
            CLratio= fopen("weighted_MC_CLratio.txt","r");
        else
            CLratio= fopen("weighted_MC_CLratio_pp.txt","r");
    }
    else
    {
        if(!m_doPP)
            CLratio= fopen("weighted_MC_CLratio_mpt_centralR4.txt","r");
        else
            CLratio= fopen("weighted_MC_CLratio_mpt_pp.txt","r");
    }
    double tmp;
    int jpt;
    int cent;
    double bkgCRatioFix[nPtBins][nCentBins];
    int loop;
    if(m_doMuonPt < -0.5)
        loop  = nPtBins*nCentBins;
    else
        loop = nMuonPtBins*nCentBins;
    

    for(int i = 0; i < loop; i++)
    {
        fscanf(CLratio, "pt = %d \t cent = %d \t ratio = %lf\n", &jpt, &cent, &tmp);
        bkgCRatioFix[jpt][cent] = tmp;
    }

    double sigRatioFix[5] = {0.374605,0.290132,0.0971048,0.0,0.400427};

    RooRealVar dataEMRatio("dataEMRatio","dataEMRatio",ratio);
    RooRealVar sigRatio("sigRatio","sigRatio",0.140,0.00,1.0);
    RooRealVar bkgLRatio("bkgLRatio","bkgLRatio",0.6,0,1.0);
    RooRealVar bkgCRatio("bkgCRatio","bkgCRatio",bkgCRatioFix[PtBins][CentBins]);
//    RooRealVar bkgCRatio("bkgCRatio","bkgCRatio",0.15,0,1.0);
  //  RooFormulaVar bkgCRatio("bkgCRatio","bkgCRatio","1.0-dataEMRatio-sigRatio-bkgLRatio",RooArgList(dataEMRatio, sigRatio, bkgLRatio));

   // RooAddPdf  *mcBkgKerPdf = new RooAddPdf("mcBkgKerPdf","combine C and Light bkg",RooArgList(*mcBkgCPdf,*mcBkgKPdf,*mcBkgPPdf,*mcBkgOtherPdf),RooArgList(*nCQuarks,*nKaons,*nPions,*nOthers));

//////////////////////////////////////////////try recursive fit////////////////////////////////////////////////

    RooAddPdf *bkgCombinePdf = new RooAddPdf("bkgCombinePdf","bkgL + bkgC", RooArgList(*mcBkgCPdf,*mcBkgLPdf), bkgCRatio); 
    RooAddPdf *sigCombinePdf = new RooAddPdf("sigCombinePdf","sigB +bkgL + bkgC", RooArgList(*mcSigBPdf,*bkgCombinePdf), sigRatio); 
  //  RooAddPdf *sigCombinePdf = new RooAddPdf("sigCombinePdf","sigB +bkgL + bkgC", RooArgList(*mcSigBPdf,*mcBkgLPdf), sigRatio); 
   // RooAddPdf *bkgCombinePdf = new RooAddPdf("bkgCombinePdf","bkgL + bkgC", RooArgList(*bkgCConvPdfp,*bkgLConvPdfp), bkgCRatio); 
    //RooAddPdf *sigCombinePdf = new RooAddPdf("sigCombinePdf","sigB +bkgL + bkgC", RooArgList(*sigConvPdfp,*bkgCombinePdf), sigRatio); 
    //RooAddPdf modelPdf("modelPdf","EM+sigB +bkgL + bkgC", RooArgList(*dataEMPdf,*sigCombinePdf), dataEMRatio); 
    RooAddPdf modelPdf("modelPdf","EM+sigB +bkgL + bkgC", RooArgList(*dataEMPdf,*sigCombinePdf), dataEMRatio); 
    //RooAddPdf modelPdf("modelPdf","EM+sigB +bkgL + bkgC", RooArgList(*mcSigBPdf,*bkgCombinePdf), sigRatio); 

//////////////////////////////////////////////try normal fit//////////////////////////////////////////////////
    //RooAddPdf modelPdf("modelPdf","bkgL + bkgC + sigpeak - EventMixing",RooArgList(*dataEMPdf,*mcSigBPdf,*mcBkgCPdf, *mcBkgLPdf),RooArgList(dataEMRatio, sigRatio, bkgCRatio)) ;     
    //RooAddPdf modelPdf("modelPdf","bkgL + bkgC + sigpeak - EventMixing",RooArgList(*dataEMPdf,*sigConvPdfp,*bkgLConvPdfp, *bkgCConvPdfp),RooArgList(dataEMRatio, sigRatio, bkgLRatio,bkgCRatio)) ;     

/////////////////////////////////////////////try chi2 fit////////////////////////////////////////////////////////


    variable.setBins(40); 
    RooDataHist* dataSubSetBinned = new RooDataHist("dataSubSetBinned","binned data",RooArgSet(variable),*dataSubSet);
   // RooAbsReal *nll_t = modelPdf.createNLL(*dataSubSetBinned, DataError(RooAbsData::SumW2),NumCPU(1));
    RooAbsReal *nll_t = modelPdf.createNLL(*dataSubSetBinned,Offset(1), NumCPU(1),Optimize(2));
   // RooChi2Var chi2("chi","chi2",modelPdf,*dataSubSetBinned,DataError(RooAbsData::SumW2),Offset(1));
    //RooMinuit minuit(chi2);
    RooMinuit minuit(*nll_t);
    minuit.migrad();
    minuit.minos();
    minuit.migrad();
    minuit.minos();
    RooFitResult* roofitResult = minuit.save() ;

//////////////////////////////////////////////Fit result/////////////////////////////////////////////////////////
/*
    variable.setBins(150); 
    RooDataHist* dataSubSetBinned = new RooDataHist("dataSubSetBinned","binned data",RooArgSet(variable),*dataSubSet);
    RooFitResult* roofitResult = modelPdf.fitTo(*dataSubSetBinned, Hesse(true), SumW2Error(true), Save(true), PrintLevel(-1),Warnings(false),PrintEvalErrors(0) );
*/
    roofitResult->Print("V");

    double fitdataEMRatio = dataEMRatio.getVal();
    double fitsigRatio = sigRatio.getVal() * (1-fitdataEMRatio);
    double fitsigRatioErr = sigRatio.getError();
    double fitbkgCRatio= bkgCRatio.getVal() * (1-sigRatio.getVal()) * (1-fitdataEMRatio);
    double fitbkgLRatio= (1-bkgCRatio.getVal()) * (1-sigRatio.getVal()) * (1-fitdataEMRatio);
    
    std::cout<<fitbkgLRatio<<"Light Ratio!"<<std::endl;
    std::cout<<fitsigRatioErr<<"Sig Ratio Error!"<<std::endl;
    
    out<<"pt = "<<PtBins<<"\t cent = "<<CentBins<<"\t sigRatio = "<<fitsigRatio<<"\t sigRatioErr = "<<fitsigRatioErr<<std::endl;
    
    TCanvas *c = PlotDataMCFitTemplate(variable, dataSubSet, modelPdf, dataEMPdf,mcSigBPdf,mcBkgCPdf, mcBkgLPdf,PtBins,CentBins,fitdataEMRatio, fitsigRatio, fitbkgCRatio);
    //TCanvas *c = PlotDataMCFitTemplate(variable, dataSubSet, modelPdf, dataEMPdf,sigConvPdfp,bkgCConvPdfp, bkgLConvPdfp,PtBins,fitdataEMRatio, fitsigRatio, fitbkgCRatio);
    std::stringstream plotname;
    plotname.str("");
    plotname<<c->GetName()<<".pdf";
    c->Print(plotname.str().c_str());
    fitResult->cd(); 
   // c->SetName((TString(fitResult.GetName())+="_")+=(c->GetName()));
    c->Write();

    RooNLLVar nll("nll","nll",modelPdf,*dataSubSet) ;

      // Plot likelihood in m0 in range that includes problematic values
      // In this configuration the number of errors per likelihood point 
      // evaluated for the curve is shown. A positive number in PrintEvalErrors(N)
      // will show details for up to N events. By default the values for likelihood
      // evaluations with errors are shown normally (unlike fitting), but the shape
      // of the curve can be erratic in these regions.

    RooPlot* frame2 = bkgCRatio.frame(Range(0.0,1.1),Title("-log(L) scan vs CRatio")) ;
    //nll.plotOn(frame2,PrintEvalErrors(0),ShiftToZero(),LineColor(kRed),Precision(1e-4)) ;
    //nll.plotOn(frame2,PrintEvalErrors(-1),EvalErrorValue(nll.getVal()+10),LineColor(kRed),Precision(1e-4)) ; 
 //   frame2->SetMaximum(15) ;
    frame2->SetMinimum(0) ;

    plotname.str("");
    plotname<<"nll_distribution_"<<PtBins<<".pdf";
   // TCanvas* c1 = new TCanvas(plotname.str().c_str(),plotname.str().c_str()) ; 
    //frame2->Draw();
    //c1->Print(plotname.str().c_str());


}

RooDataSet* TempFitter::SelectPtEtaCent(RooDataSet* dataSet, double ptLow, double ptUpp, double centLow, double centUpp) 
{
    std::stringstream cut;
    cut << "Jet_pt>" << ptLow << "&&Jet_pt<" << ptUpp << "&& Centrality_hp > "<< centLow << "&& Centrality_hp < "<< centUpp;
  
    return (RooDataSet*)dataSet->reduce(cut.str().c_str());
  
}

RooDataSet* TempFitter::SelectMuonPtEtaCent(RooDataSet* dataSet, double ptLow, double ptUpp) 
{
    std::stringstream cut;
    cut << "Muon_asso_pt>" << ptLow << "&&Muon_asso_pt<" << ptUpp;
  
    return (RooDataSet*)dataSet->reduce(cut.str().c_str());
  
}

TCanvas* TempFitter::PlotDataMCFitTemplate(RooRealVar& var, RooDataSet* dataSubSet,  RooAbsPdf& modelPdf, RooKeysPdf* dataEMPdf, RooAbsPdf* mcSigBPdf, RooAbsPdf* mcBkgCPdf, RooAbsPdf* mcBkgLPdf, int PtBins,int CentBins,double dataEMRatio,double sigRatio, double bkgCRatio) 
{
      
    SetAtlasStyle();
    TCanvas *c = new TCanvas("c","c");
    std::stringstream name;
    name.str("");
    if(m_doMuonPt > -0.5)
    {
        if(!m_doPP)
        name<<"TemplateFit_pt_diff_muonpt_"<<PtBins;
        else
        name<<"TemplateFit_pt_diff_muonpt_"<<PtBins<<"_pp";
    }
    else
    {
        if(!m_doPP)
        name<<"TemplateFit_pt_"<<PtBins<<"_cent_"<<CentBins;
        else
        name<<"TemplateFit_pt_"<<PtBins<<"_cent_"<<CentBins<<"_pp";
    }
    c->SetName(name.str().c_str());
    c->cd();
    RooPlot* fitFrame = var.frame();
    name<<"_rooPlot";
    fitFrame->SetName(name.str().c_str());
    dataSubSet->plotOn(fitFrame);
    modelPdf.plotOn(fitFrame,LineColor(kBlack));
    modelPdf.plotOn(fitFrame,Components(*mcSigBPdf),LineColor(2));
    modelPdf.plotOn(fitFrame,Components(*mcBkgCPdf),LineColor(kBlue+1));
    modelPdf.plotOn(fitFrame,Components(*mcBkgLPdf),LineColor(kGreen-8));
    modelPdf.plotOn(fitFrame,Components(*dataEMPdf),LineStyle(kDashed),LineColor(6));
    fitFrame->GetXaxis()->SetTitle(var.getTitle());
    fitFrame->Draw();

    //TLegend *leg = new TLegend(0.3,0.7,0.6,0.925);
    TLegend *leg = new TLegend(0.6,0.7,0.9,0.925);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    TLegendEntry* entry = leg->AddEntry("","Pb+Pb (#sqrt{#font[52]{s_{NN}}} = 2.76 TeV)","p");
    entry->SetMarkerStyle(20);
    entry->SetMarkerSize(1.2);
    entry = leg->AddEntry("","Best-fit template","l");
    entry->SetLineColor(kBlack);
    entry->SetLineWidth(2);
    entry = leg->AddEntry("","Data event mix","l");
    entry->SetLineColor(6);
    entry->SetLineWidth(2);
    entry = leg->AddEntry("","Signal","l");
    entry->SetLineColor(2);
    entry->SetLineWidth(2);
    entry = leg->AddEntry("","Bkg C Quark","l");
    entry->SetLineColor(kBlue+1);
    entry->SetLineWidth(2);
    entry = leg->AddEntry("","Bkg Light Quark","l");
    entry->SetLineColor(kGreen-8);
    entry->SetLineWidth(2);
    leg->Draw();
    TLatex la;
    la.SetNDC();
    std::stringstream namea;
    namea.str("");
    if(m_doMuonPt > -0.5)
        namea << m_muonptBins[PtBins] << " < p_{T} < " << m_muonptBins[PtBins+1] << " GeV";
    else
        namea << m_ptBins[PtBins] << " < p_{T} < " << m_ptBins[PtBins+1] << " GeV";
    la.DrawLatex(0.48,0.65,namea.str().c_str());
    namea.str("");
    namea << "Signal B Quark ratio: "<< sigRatio;
    la.SetTextSize(0.035);
    la.DrawLatex(0.48,0.57,namea.str().c_str());
    namea.str("");
    namea << "Bkg C Quark ratio: "<< bkgCRatio;
    la.SetTextSize(0.035);
    la.DrawLatex(0.48,0.54,namea.str().c_str());
    namea.str("");
    namea << "Event mix ratio: "<< dataEMRatio;
    la.SetTextSize(0.035);
    la.DrawLatex(0.48,0.51,namea.str().c_str());



    namea.str("");
    la.SetTextSize(0.045);
    namea<<m_centralityBins[CentBins]<<"-"<<m_centralityBins[CentBins+1]<<"%";
    //la.DrawLatex(0.48,0.45,namea.str().c_str());
    
    if(CentBins==0 && nCentBins > 1) la.DrawLatex(0.48,0.45,"0-10%");
    if(CentBins==0 && nCentBins == 1) la.DrawLatex(0.48,0.45,"0-100%");
    if(CentBins==1) la.DrawLatex(0.48,0.45,"10-20%");
    if(CentBins==2) la.DrawLatex(0.48,0.45,"20-40%");
    if(CentBins==3) la.DrawLatex(0.48,0.45,"40-60%");
    if(CentBins==4) la.DrawLatex(0.48,0.45,"60-100%");


/*
    namea.str("");
    namea << m_etaBins[iEta] << " < #eta < " << m_etaBins[iEta+1];
    la.DrawLatex(0.48,0.6,namea.str().c_str());

    //fitResult.addCanvas(c);
*/  
    return c;
  

}

RooDataSet* TempFitter::CombineMCJ(bool do_inclusive, RooArgSet& JetArgSetMC, bool do_data , bool do_eventmix)
{
    double cross_sections[3]={1.8760E+05,8.2787E+03,2.9419E+02};
    double filter_efficiency_inclusive[3]={1.0,1.0,1.0};
    int event_Number_inclusive = 4500000;
    int event_Number_muon= 650000;
    double filter_efficiency_muon_gen[3]={6.6523E-03,1.4941E-02,2.4284E-02};
//    double filter_efficiency_muon[3]={1.7745E-03,8.4775E-03,1.7822E-02};
    double scale_factor[3];
 
    std::stringstream ss; 
    std::stringstream s; 

    TFile* files[3];
    TTree* MCJetTree[3];
    TTree* JetTree;
    RooDataSet* mcSet[3];
    RooDataSet* dataSet[3];
    RooDataSet *mcSetComb ;
    RooDataSet *dataSetComb ;
    RooRealVar  Ptrel("Ptrel","Ptrel",0.0,0,5);
    RooRealVar var = Ptrel;

    RooRealVar* w[3] ;
    RooDataSet *mcSetWeight [3];
    RooDataSet *dataSetWeight [3];
   // RooRealVar roo_scale[3];

    if(do_data)
    {
/*
        ss.str("");
        ss << "/atlas/data4/tingting/data6/ReconCode/histogram_Data11_AA_minbia_randm_aveweight_addmuon_sig_trackmuon.root";

        files[0]=TFile::Open(ss.str().c_str());
        RooCategory datasample("datasample","datasample") ;
        if(!do_eventmix)
        {
            JetTree = (TTree*)files[0]->Get("DataJetTree");
            dataSet[0] = new RooDataSet("dataSet","dataSet",JetTree,JetArgSetMC);
            datasample.defineType("data") ;
        }
        else 
        {
            JetTree = (TTree*)files[0]->Get("DataJetTree_em");
            dataSet[1] = new RooDataSet("dataSet_eventmix","dataSet Event mix",JetTree,JetArgSetMC);
            datasample.defineType("data_eventmix") ;
        }

        RooRealVar  weight("weight","weight",-1.0,1.0);
        JetArgSetMC.add(weight);

        weight = 1.0;
        w[0] = (RooRealVar*) dataSet->addColumn(weight);

        dataSetWeight[0]= new RooDataSet(dataSet[0]->GetName(),dataSet[0]->GetTitle(),dataSet[0],*dataSet[0]->get(),0,w[0]->GetName()) ;

        dataSetComb = new RooDataSet("dataSetComb","combined data",JetArgSetMC,Index(datasample),Import("data",*dataSetWeight[0]),Import("data_eventmix",*dataSetWeight[1]),WeightVar("weight"));
     
        dataSetComb->Print();

        return dataSetComb;
*/
    }

    else
    {
        for(int J=1; J<=3; J++)
        {   
            ss.str("");
            if(do_inclusive)
                ss << "/atlas/data4/tingting/data6/ReconCode/histogram_track_inclusive_D3PD_loosetruth_J" << J << ".root";
            else
                ss << "/atlas/data4/tingting/data6/ReconCode/histogram_track_muon_D3PD_loosetruth_J" << J << ".root";

            files[J-1]=TFile::Open(ss.str().c_str());
            MCJetTree[J-1] = (TTree*)files[J-1]->Get("MCJetTree");
            scale_factor[J-1]=cross_sections[J-1];
            if(do_inclusive) scale_factor[J-1]*=filter_efficiency_inclusive[J-1]/event_Number_inclusive;
            else scale_factor[J-1]*=filter_efficiency_muon_gen[J-1]/event_Number_muon;


            ss.str("");
            ss<<"mcSet_different_J"<<J;
            mcSet[J-1] = new RooDataSet(ss.str().c_str(),ss.str().c_str(),MCJetTree[J-1],JetArgSetMC);
        }    
         //   wFunc[J-1] = new RooFormulaVar("w","event weight","(0.5)",RooArgList(var)) ;
           // w[J-1] = (RooRealVar*) mcSet[J-1]->addColumn(wFunc[J-1]);
/*
        for(int J=1; J <=3; J++)
            scale_factor[J-1] = (float)scale_factor[J-1]/scale_factor[2];
*/
        RooRealVar roo_scale_1("roo_scale_1","roo_scale_1",scale_factor[0]);
        RooRealVar roo_scale_2("roo_scale_2","roo_scale_2",scale_factor[1]);
        RooRealVar roo_scale_3("roo_scale_3","roo_scale_3",scale_factor[2]);

        RooRealVar  weight("weight","weight",0.0,1.0);
        JetArgSetMC.add(weight);

        weight = scale_factor[0];
        w[0] = (RooRealVar*) mcSet[0]->addColumn(weight);
        weight = scale_factor[1];
        w[1] = (RooRealVar*) mcSet[1]->addColumn(weight);
        weight = scale_factor[2];
        w[2] = (RooRealVar*) mcSet[2]->addColumn(weight);
        for(int J = 1; J<=3; J++)
        {
            mcSetWeight[J-1]= new RooDataSet(mcSet[J-1]->GetName(),mcSet[J-1]->GetTitle(),mcSet[J-1],*mcSet[J-1]->get(),0,w[J-1]->GetName()) ;
            mcSetWeight[J-1] ->Print();
        }   

        RooCategory mcsample("mcsample","mcsample") ;
        mcsample.defineType("J1") ;
        mcsample.defineType("J2") ;
        mcsample.defineType("J3") ;

        mcSetComb = new RooDataSet("mcSetComb","combined MC",JetArgSetMC,Index(mcsample),Import("J1",*mcSetWeight[0]),Import("J2",*mcSetWeight[1]),Import("J3",*mcSetWeight[2]),WeightVar("weight"));
     
        mcSetComb->Print();
        return mcSetComb;
    }

}
/*
void TempFitter::CheckNegative(RooDataSet* d0)
{
    RooArgSet* set;
    double val;
    RooRealVar
    int entries = d0->numEntries();
    for(int i = 0; i < entries; i++)
    {
        set = (RooArgSet*)d0->get(i); 
        xi = (RooRealVar*)set->find("Ptrel");
        val = xi->getVal();
        xi->setVal(val);
        prob = pdf->getVal(RooArgSet(*xi));
    }
}

*/

