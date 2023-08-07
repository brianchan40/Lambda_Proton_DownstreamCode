#include <iostream>
#include "TGraph.h"
#include "TFile.h"

TF1 *fit_eff_pT_final[9], *fit_tof_eff_pT_final[9];
TGraph *efficiency_pt[9];

void extract_proton_eff();
Double_t function_sum(Double_t *x, Double_t *par); 

void plot_proton_eff()
{
    TFile f("proton_eff_plots_27.root", "recreate");
    extract_proton_eff();

    TCanvas *all_graphs = new TCanvas("all_graphs", "Proton Efficiency Graphs", 1500, 800);

    for(int lmn = 0 ; lmn < 9 ; lmn++){

        // TF1 *temp_f1 = new TF1("Multiply",function_sum, 0.2, 2.5, 1);

        float xx[46] = {0.}, yy[46] = {0.};

        for(int x_c = 0; x_c < 46; x_c++){
            xx[x_c] = 0.2 + x_c * 0.05;
            yy[x_c] = fit_eff_pT_final[lmn]->Eval(xx[x_c])*fit_tof_eff_pT_final[lmn]->Eval(xx[x_c]);
        }

        efficiency_pt[lmn] = new TGraph(46, xx, yy);


        cout << lmn << endl;

        // temp_f1->SetParameters(lmn);

        cout << "past here?" << endl;

        // efficiency_pt[lmn] = temp_f1;

        efficiency_pt[lmn]->SetTitle("p/#bar{p} Efficiency(p_{T})");

        efficiency_pt[lmn]->GetXaxis()->SetTitle("p_{T} (GeV/c)");
        efficiency_pt[lmn]->GetYaxis()->SetTitle("Efficiency %");
        if((lmn+1)<5){
            efficiency_pt[lmn]->SetMarkerColor(lmn+1);
            efficiency_pt[lmn]->SetLineColor(lmn+1);
        }
        else if((lmn+2)<10){
            efficiency_pt[lmn]->SetMarkerColor(lmn+2);
            efficiency_pt[lmn]->SetLineColor(lmn+2);
        }
        else{
            efficiency_pt[lmn]->SetMarkerColor(lmn+3);
            efficiency_pt[lmn]->SetLineColor(lmn+3);
        }
        if(lmn == 0) efficiency_pt[lmn]->Draw();
        else efficiency_pt[lmn]->Draw("SAMES L");
    }

    TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->SetHeader("Centrality");
    leg->AddEntry(efficiency_pt[0],"70-80%","l");
    leg->AddEntry(efficiency_pt[1],"60-70%","l");
    leg->AddEntry(efficiency_pt[2],"50-60%","l");
    leg->AddEntry(efficiency_pt[3],"40-50%","l");
    leg->AddEntry(efficiency_pt[4],"30-40%","l");
    leg->AddEntry(efficiency_pt[5],"20-30%","l");
    leg->AddEntry(efficiency_pt[6],"10-20%","l");
    leg->AddEntry(efficiency_pt[7],"5-10%","l");
    leg->AddEntry(efficiency_pt[8],"0-5%","l");
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    

    leg->Draw();

    f.cd();
    all_graphs->Write();

    f.Close();
}

void extract_proton_eff()
{
    std::fstream eff_proton_file("./27GeV/proton_efficiency_coefficients.txt", std::ios_base::in);

    float a = 0;
    int count = 0;
    float par_sig[5] = {0.};

    while (eff_proton_file >> a)
    {
        par_sig[count % 5] = a;

        if (count % 5 == 4)
        {
            // if (debug_1)
                // cout << "count / 5 = " << count / 5 << endl;
            fit_eff_pT_final[count / 5] = new TF1("fit_eff_pT_final", "([0]+[3]*x+[4]*x*x)*exp(-pow([1]/x,[2]))", 0.2, 2.5);
            fit_eff_pT_final[count / 5]->SetParameters(par_sig[0], par_sig[1], par_sig[2], par_sig[3], par_sig[4]);
        }

        count++;
    }

    std::fstream tofeff_proton_file("./27GeV/proton_tofefficiency_coefficients.txt", std::ios_base::in);

    a = 0;
    count = 0;

    float par_sig_tof[4] = {0.};

    while (tofeff_proton_file >> a)
    {
        par_sig_tof[count % 4] = a;

        if (count % 4 == 3)
        {
            // if (debug_1)
                // cout << "count / 4 = " << count / 4 << endl;
            fit_tof_eff_pT_final[count / 4] = new TF1("fit_tof_eff_pT_final", "[0] + [1]*sqrt(x) + [2]*x + [3]*pow(x,2)", 0.2, 2.5);
            fit_tof_eff_pT_final[count / 4]->SetParameters(par_sig_tof[0], par_sig_tof[1], par_sig_tof[2], par_sig_tof[3]);

            // if (debug_1)
                // cout << "parameters are: " << par_sig_tof[0] << ", " << par_sig_tof[1] << ", " << par_sig_tof[2] << ", " << par_sig_tof[3] << endl;
        }

        count++;
    }
}

Double_t function_sum(Double_t *x, Double_t *par)
{
    const Double_t xx =x[0];
    //~ f1->SetParameters(par[0],par[1],par[2],par[3],par[4]);
    //~ f2->SetParameters(par[5],par[6],par[7],par[8],par[9]);
    cout << "par[0] = " << par[0] << endl;
    return fit_eff_pT_final[par[0]]->Eval(xx)*fit_tof_eff_pT_final[par[0]]->Eval(xx);
}





