#ifndef _Rebin_v2_Data_
#define _Rebin_v2_Data_

using namespace std;

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <TChain.h>
#include "TLeaf.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TMath.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TKey.h"
#include "TObject.h"
#include "TObjArray.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include <string>
#include <stdlib.h>

TH1D *newHistFlow;
TH1D *newHistFlowBinCenter;
TH1D *newHistYield;

Char_t *inFile;
Float_t xbins[1000], xbins_p[1000];
float xx[1000], xx_err[1000];
float v2[1000], v2_err[1000], v2_rec[1000], v2_rec_err[1000];
// char fname[200];

TFile *lam_pt;

//TFile *outputfile_forRebinned;

static TGraphAsymmErrors *Rebin4Eta(TH1D *histFlow, TH1D *histYield, const Int_t newbins, const Float_t *xbins);
//int init_conditions(char *inputfile);
int init_conditions(int cint, int cend, string lam, char *proton_option);
void Mainfunction_rebin(int cen, int newbins, const char* prof_name1, const char* prof_name2, int ep_option);
void Print_rebin(int newbins);

#endif
