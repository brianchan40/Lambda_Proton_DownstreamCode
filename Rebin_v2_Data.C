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

#include <Rebin_v2_Data.h>

vector<float> Rebin_v2_Data(const char* prof_name1, const char* prof_name2, int cen, int ep_option)
{
    Int_t newbins = 1;
    xbins[0] = 0;
    xbins[1] = 4.0;

    Mainfunction_rebin(cen, newbins, prof_name1, prof_name2, ep_option);

    vector<float> final_v2;
    final_v2.clear();
    final_v2.push_back(newbins);
    for (int i = 0; i < newbins; i++)
    {
        final_v2.push_back(v2_rec[i]);
        // cout << "push_back: v2_rec[i] = " << v2_rec[i] << endl;
    }
    for (int i = 0; i < newbins; i++)
    {
        final_v2.push_back(v2_rec_err[i]);
        // cout << "i = " << i << endl;
    }

    return final_v2;
}

void Mainfunction_rebin(int cen, int newbins, const char* prof_name1, const char* prof_name2, int ep_option)
{
    bool y_err_zero = false;
    int y_err_index = -1;

    file->cd();
    TProfile *v2_vs_pt = (TProfile *)file->Get(prof_name1);
    // cout << "reso[ep_option] = " << reso[ep_option] << endl;
    v2_vs_pt->Scale(1.0/reso[ep_option]);
    output_File->cd();
    v2_vs_pt->Write(prof_name1);
    TProfile *yield_hist = (TProfile *)file->Get(prof_name2);

    TH1D *lam_pt_yield = new TH1D("yield", "Yield", 50, 0, 5);
    
    for (int j = 0; j < 50; j++)
    {
        lam_pt_yield->Fill((j * 0.1 + 0.05), yield_hist->GetBinEntries(j + 1));
    }

    TGraphAsymmErrors *gpt = Rebin4Eta(v2_vs_pt, lam_pt_yield, newbins, xbins);
    for (int j = 0; j < newbins; j++)
    {
        if (j == y_err_index)
            continue;
        // cout << "j = " << j << endl;
        double x, y, y_err;
        gpt->GetPoint(j, x, y);
        xx[j] = x;
        // cout << "x = " << x << endl;
        // cout << "y = " << y << endl;
        y_err = gpt->GetErrorY(j);
        // cout << "y_err = " << y_err << endl;
        // v2_rec[j] += y / y_err / y_err * reso[j];
        if (y_err == 0)
        {
            y_err_index = j;
            // cout << "y_err_index = " << y_err_index << endl;
            continue;
        }
        v2_rec[j] += y / y_err / y_err;
        // cout << "v2_rec[j] = " << v2_rec[j] << endl;
        // v2_rec_err[j] += 1 / y_err / y_err * reso[j] * reso[j];
        v2_rec_err[j] += 1 / y_err / y_err;
    }

    delete lam_pt_yield;
    lam_pt_yield = NULL;
    delete gpt;
    gpt = NULL;

    for (int j = 0; j < newbins; j++)
    {
        // cout << "j = " << j << endl;
        // res[j] /= res_err[j];
        // res_err[j] = 1 / sqrt(res_err[j]);
        if (j == y_err_index)
        {
            v2_rec[j] = 0;
            v2_rec_err[j] = 0;
            continue;
        }
        // cout << "2 v2_rec[j] = " << v2_rec[j] << endl;
        v2_rec[j] /= v2_rec_err[j];
        // cout << "3 v2_rec[j] = " << v2_rec[j] << endl;
        // cout << "v2_rec_err[j] = " << v2_rec_err[j] << endl;
        v2_rec_err[j] = 1 / sqrt(v2_rec_err[j]);
        // cout << "2 v2_rec_err[j] = " << v2_rec_err[j] << endl;
    }
}

void Print_rebin(int newbins)
{
    cout << "pt position:" << endl;
    for (int j = 0; j < newbins; j++)
        cout << xx[j] << ",";
    cout << endl;
    cout << endl;
    /*cout << "EP resolution:" << endl;
    for(int j = 0; j < 25; j++) cout << res[j] << ",";
    cout << endl;
    for(int j = 0; j < 25; j++) cout << res_err[j] << ",";
    cout << endl;
    cout << endl;*/
    cout << "v2 reconstructed:" << endl;
    for (int j = 0; j < newbins; j++)
        cout << v2_rec[j] << ",";
    cout << endl;
    for (int j = 0; j < newbins; j++)
        cout << v2_rec_err[j] << ",";
    cout << endl;
    cout << endl;
}

//-------------------------------------------------------------
TGraphAsymmErrors *Rebin4Eta(TH1D *histFlow, TH1D *histYield, const Int_t newbins, const Float_t *xbins)
{

    newHistFlow = new TH1D("Flow_vx", "plot", newbins, xbins);
    newHistFlowBinCenter = new TH1D("Flow_vx_center", "plot", newbins, xbins);
    newHistYield = new TH1D("new Yield x", "dN/dx", newbins, xbins);

    // cout << "Here?" << endl;
    Int_t nBins = histFlow->GetNbinsX();
    // cout << "Here? 1" << endl;

    double v, verr;
    double vSum;
    double content;
    double error;
    double error2sum;
    double vSqrSum;

    double yield;
    double yieldSum;
    double y;
    double pt;
    double x, xnew;

    double meanEtaSum;
    double meanEta;

    // cout << "Here?" << endl;

    for (int newbin = 1; newbin < newbins + 1; newbin++)
    {

        vSum = 0.;
        yieldSum = 0.;
        error2sum = 0.;
        vSqrSum = 0.;
        meanEtaSum = 0.;

        for (int bin = 1; bin < nBins + 1; bin++)
        {
            v = 0.;
            content = 0.;
            error = 0.;
            yield = 0.;
            meanEta = 0.;

            x = histYield->GetBinCenter(bin);
            xnew = newHistYield->GetBinLowEdge(newbin);

            if (x < xnew)
                continue;
            if (x > newHistYield->GetBinLowEdge(newbin + 1))
                break;

            yield = histYield->GetBinContent(bin);
            v = histFlow->GetBinContent(bin);
            verr = histFlow->GetBinError(bin);

            if (v != 0)
            {
                yieldSum += yield;
                // if(newbin == 1) cout << "yield = " << yield << endl;
                vSum += yield * v;
                error2sum += pow(yield * histFlow->GetBinError(bin), 2.);
                vSqrSum += yield * (yield * pow(histFlow->GetBinError(bin), 2.) + pow(v, 2.));
                meanEtaSum += yield * histFlow->GetBinCenter(bin);
            }
        }

        if (yieldSum)
        {
            content = vSum / yieldSum;
            error = sqrt((vSqrSum / yieldSum) - pow((vSum / yieldSum), 2.)) / sqrt(yieldSum);
            meanEta = meanEtaSum / yieldSum;

            if(newbin == 1)
            {
                // cout << "content = " << content << endl;
                // cout << "meanEta = " << meanEta << endl;
                // cout << "error = " << error << endl;
                // cout << "(vSqrSum / yieldSum) = " << (vSqrSum / yieldSum) << endl;
                // cout << "pow ( ( vSum / yieldSum), 2.) = " << pow ( ( vSum / yieldSum), 2.) << endl;
            }
        }

        newHistYield->SetBinContent(newbin, yieldSum * histYield->GetBinWidth(1) / newHistYield->GetBinWidth(newbin));
        newHistYield->SetBinError(newbin, sqrt(yieldSum * histYield->GetBinWidth(1) / newHistYield->GetBinWidth(newbin)));

        newHistFlow->SetBinContent(newbin, content);
        newHistFlow->SetBinError(newbin, error);
        newHistFlowBinCenter->SetBinContent(newbin, meanEta);
    }

    TGraphAsymmErrors *gr = new TGraphAsymmErrors();

    for (int i = 1; i < newHistFlow->GetNbinsX() + 1; i++)
    {
        // cout << "newHistFlowBinCenter->GetBinContent(i) = " << newHistFlowBinCenter->GetBinContent(i) << ", newHistFlow->GetBinContent(i) = " << newHistFlow->GetBinContent(i) << endl;

        gr->SetPoint(i - 1, newHistFlowBinCenter->GetBinContent(i), newHistFlow->GetBinContent(i));
        gr->SetPointError(i - 1,
                          (newHistFlowBinCenter->GetBinContent(i) - (newHistFlow->GetBinLowEdge(i))),
                          ((newHistFlow->GetBinLowEdge(i) + newHistFlow->GetBinWidth(i)) - newHistFlowBinCenter->GetBinContent(i)),
                          newHistFlow->GetBinError(i),
                          newHistFlow->GetBinError(i));
    }

    newHistFlow->Clear();
    newHistFlowBinCenter->Clear();
    newHistYield->Clear();

    return gr;
}