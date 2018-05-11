/* vim: set sw=4 ts=4 expandtab: */

#ifndef TEMPFITTER_H_
#define TEMPFITTER_H_
#include <string>
#include <vector>
#include <algorithm>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <RooDataSet.h>
#include <RooArgSet.h>
#include <RooKeysPdf.h>
#include <RooRealVar.h>
#include <RooShiftedKeysPdf.h>

using namespace RooFit;
using namespace std;

class TempFitter
{
    public:
    TempFitter(bool doBuildTemplate = true, bool doComposite = false, bool doIPSignificance = false, int doMuonPt = -1.0, bool doPP = false);
    void            DoTempFitting();
    RooDataSet*     CombineMCJ(bool , RooArgSet&, bool, bool);
    void            FillDataSet(bool, bool, bool , RooArgSet& , RooDataSet* ,double , double , int, TH1F*, TH1F*, TH1F*);
    RooDataSet*     SelectPtEtaCent(RooDataSet*, double, double,double,double);
    RooDataSet*     SelectMuonPtEtaCent(RooDataSet*, double, double);
    void            FitDataSet(RooDataSet*,  int,int, TFile*, float ,ofstream&);
    void            BuildTemplate(RooDataSet*, RooDataSet*, RooDataSet*,RooRealVar& ,int,int);
    void            BuildTemplate_MC(RooDataSet*, RooDataSet*, RooDataSet*,RooRealVar& ,int,int);
    TCanvas*        PlotDataMCFitTemplate(RooRealVar& , RooDataSet* , RooAbsPdf& ,RooKeysPdf*, RooAbsPdf*, RooAbsPdf*, RooAbsPdf*,int, int, double, double, double);
//    void*           FillMCJSet(bool, RooArgSet&, RooDataSet*, int, TH1F*);
    //void            CheckNegative(RooDataSet*);

    std::vector<double> m_ptBins;
    std::vector<double> m_centralityBins;
    std::vector<double> m_muonptBins;

    int nPtBins;
    int nMuonPtBins;
    int nCentBins;

    TFile* file_FitPlot;
    TFile* file_FitResult;
// Flavour_tag :0->Data;1->MC bjet;2->MC cjet;3->MC light kjet;4->MC light pjet;
//5->MC light no usual jet; 6->MC fake jet(muon not asso truth muon);7-> Event mixing;
    
    private:
    bool            m_doBuildTemplate;
    bool            m_doComposite;
    bool            m_doIPSignificance;
    int             m_doMuonPt;
    bool            m_doPP;
    
};
#endif

