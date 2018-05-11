/* vim: set sw=4 ts=4 expandtab: */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <TFile.h>
#include <RooCategory.h>

void RunTempFitter()
{
    gROOT->ProcessLine(".L RooShiftedKeysPdf.cpp++");
    gROOT->ProcessLine(".L TempFitter.C+");

    TempFitter fitter(1,0,0,-1.0,0.0); 
    fitter.DoTempFitting();
}
