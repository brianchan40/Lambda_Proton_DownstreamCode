using namespace std;

#include "stdio.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TString.h"
#include <fstream>
#include <iostream>
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
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TMultiGraph.h"

double lam_purity[9][18];
double antilam_purity[9][18];
double lam_yield[9][18];
double antilam_yield[9][18];

double pt_bins[18] = {0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 1.95, 2.05, 2.15, 2.5};
double pt_bin_errors[18] = {0.};

void read_sig()
{
    std::fstream myfile("./purity.txt", std::ios_base::in);

    double a;
    int count = 0;
    while (myfile >> a)
    {
        if(count < 306)
        {
            //cout << "a = " << a << endl;
            if(count % 2 == 0)
            {
                lam_purity[count / 34][(count - (count / 34) * 34) / 2] = a;
                //cout << "(count-(count / 34)*34) = " << (count-(count / 34)*34) << endl;
                //cout << "count / 34 = " << count / 34 << endl;
                //cout << "lam_purity[count / 34][(count-(count / 34)*34)/2] = " << lam_purity[count / 34][(count-(count / 34)*34)/2];
            }
            else if(count % 2 == 1)
            {
                lam_yield[count / 34][(count - (count / 34) * 34) / 2] = a;
                //cout << "(count-(count / 34)*34) = " << endl;
                //cout << ", lam_yield[count / 34][(count-(count / 34)*34)/2] = " << lam_yield[count / 34][(count-(count / 34)*34)/2] << endl;
            }
        }
        else if(count < 612)
        {
            if(count % 2 == 0)
            {
                antilam_purity[(count - 306) / 34][(count - (count / 34) * 34) / 2] = a;
                //cout << "antilam_purity[(count-306) / 34][(count-(count / 34)*34)/2] = " << antilam_purity[(count-306) / 34][(count-(count / 34)*34)/2];
            }
            else if(count % 2 == 1)
            {
                antilam_yield[(count - 306) / 34][(count - (count / 34) * 34) / 2] = a;
                //cout << ", antilam_yield[(count-306) / 34][(count-(count / 34)*34)/2] = " << antilam_yield[(count-306) / 34][(count-(count / 34)*34)/2] << endl;
            }
        }
        else if(count == 612) break;

        count++;
    }

    std::fstream secondfile("./Without_ptsplit/purity.txt", std::ios_base::in);

    count = 0;
    while (secondfile >> a)
    {
        if(count < 18)
        {
            //cout << "a = " << a << endl;
            if(count % 2 == 0)
            {
                lam_purity[count / 2][17] = a;
                // if(debug_resolution == true) cout << "lam_purity[count / 2][17] = " << lam_purity[count / 2][17] << endl;
                // cout << "lam_purity[count / 2][17] = " << lam_purity[count / 2][17] << endl;
            }
            else if(count % 2 == 1)
            {
                lam_yield[count / 2][17] = a;
                // if(debug_resolution == true) cout << "lam_yield[count / 2][17] = " << lam_yield[count / 2][17] << endl;
            }
        }
        else if(count < 36)
        {
            if(count % 2 == 0)
            {
                antilam_purity[(count - 18) / 2][17] = a;
                // if(debug_resolution == true) cout << "antilam_purity[(count-18) / 2][17] = " << antilam_purity[(count - 18) / 2][17] << endl;
            }
            else if(count % 2 == 1)
            {
                antilam_yield[(count - 18) / 2][17] = a;
                // if(debug_resolution == true) cout << "antilam_yield[(count-18) / 2][17] = " << antilam_yield[(count - 18) / 2][17] << endl;
            }
        }
        else if(count == 36) break;

        count++;
    }
}

void plot_purity(){
    
    read_sig();

    TCanvas *c1 = new TCanvas("c1", "Purity by Lambda/AntiLambda pT", 1500, 800);
    c1->Divide(3,3);

    for(int i = 0; i < 9; i++){

        double lam_purity_tmp[18] = {0.};
        double antilam_purity_tmp[18] = {0.};

        for(int j = 0; j < 18; j++){
            lam_purity_tmp[j] = lam_purity[i][j];
            antilam_purity_tmp[j] = antilam_purity[i][j];
        }

        TMultiGraph *mg = new TMultiGraph();

        TGraph *graph1 = new TGraph(18, pt_bins, lam_purity_tmp);
        TGraph *graph2 = new TGraph(18, pt_bins, antilam_purity_tmp);
        graph1->SetMarkerSize(1);
        graph1->SetMarkerStyle(kFullSquare);
        graph2->SetMarkerSize(1);
        graph2->SetMarkerStyle(kFullCircle);
        graph2->SetMarkerColor(kRed);
        graph2->SetLineColor(kRed);

        // graph1->GetXaxis()->SetTitle("p_{T} GeV/c^{2}");
        // graph1->GetYaxis()->SetTitle("N_{sig}/N_{total}");

        mg->Add(graph1);
        mg->Add(graph2);

        c1->cd(i+1);
        mg->Draw("AP");
        mg->SetTitle(TString::Format("Centrality %d", (i+1)));
        mg->GetXaxis()->SetTitle("p_{T} GeV/c^{2}");
        mg->GetXaxis()->SetTitleSize(0.06);
        mg->GetXaxis()->SetLabelSize(0.05);
        mg->GetYaxis()->SetTitle("N_{sig}/N_{total}");
        mg->GetYaxis()->SetTitleSize(0.06);
        mg->GetYaxis()->SetLabelSize(0.05);
    }

    c1->Draw();
}




