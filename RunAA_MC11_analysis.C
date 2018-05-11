
/* vim: set sw=4 ts=4 expandtab: */

#include <iostream>
#include "TChain.h"
#include "TString.h"
#include <sstream>
#include <vector>
void RunAA_MC11_analysis(int infile)
{
    gROOT->ProcessLine(".x $ROOTCOREBIN/scripts/load_packages.C");
    gROOT->ProcessLine(".L AA_MC11_analysis.C++");
//    gROOT->ProcessLine(".L AA_MC11_analysis_C.so");

    TChain* chain=new TChain("HeavyIonD3PD");
    TString file_name;
    ifstream fin;
    std::stringstream condor_infile;
    condor_infile.str("");
    condor_infile<<infile<<".infile";
    fin.open(condor_infile.str().c_str());

//    fin.open("MC11_check.txt");
    char ntuple_name[300];
    unsigned int filecount=0;
    std::vector<std::string> input_files;
    while(fin.good())
    {
        if(filecount >= 1000000) break;
        fin.getline(ntuple_name,300);
        std::string fname(ntuple_name);
        if(fname.find("skimmed_tree")==std::string::npos)
        //if(fname.find("NTUP")==std::string::npos)
            continue;
        chain->AddFile(ntuple_name);
        std::cout << "Adding "<<fname<< std::endl;
        filecount++;
    }


    AA_MC11_analysis filter(chain,false); //keep all branches
    filter.Init(chain);
    filter.SetGraph();
    //filter.SetJESCalibrator(2,"MC11","/atlas/data3/angerami/common_install/RootCore/data/HIJetCalibration/HIJetCalibConstants.root");
    filter.SetJESCalibrator("antikt2HIItrEMFR_MC11","/atlas/data3/angerami/common_install/RootCore/data/HIJetCalibration/HIJetCalibConstants_PbPb2p76TeV.root");
    filter.SetCentralityGetter("data11_hi.centrality.csv");

    filter.Loop();
    //filter.LoopTruth();
    filter.WriteGraph(infile);
 //   filter.FillFcal(infile);


}
