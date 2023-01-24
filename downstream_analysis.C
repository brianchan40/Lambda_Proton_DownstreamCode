#include <TProfile>
// #include <TCanvas>
// #include <TString>
#include <fstream>
#include <iostream>
#include <TFile>
#include <TF1>
#include "Rebin_v2_Data.C"
#include <TGraphErrors>
#include <error_calc.C>

using namespace std;

TFile *file, *output_File, *output_File2;
ofstream result_file;
float reso[3] = {0.}, reso1[3] = {0.};
float cen9_eff[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
float cen9_err_eff[9] = {0.};
double npart[9];
double npart_err[9];

const char *method_name = "Traditional";

double lam_purity[9][18];
double antilam_purity[9][18];
double lam_yield[9][18];
double antilam_yield[9][18];

float lam_purity_Q2[3][9][100], antilam_purity_Q2[3][9][100];

std::vector<float> parent_v2_vect[3], alltrks_v2_vect[3], pionpion_v2_vect[3];
std::vector<TString> name_options3;
float reso_TPC[9] = {0.}, reso_EPD[9] = {0.}, reso_EPD1[9] = {0.};
float reso_TPC_err[9] = {0.}, reso_EPD_err[9] = {0.}, reso_EPD1_err[9] = {0.};
float gamma[2][9][3], gamma_err[2][9][3];
float gamma_ensemble[2][9][3], gamma_ensemble_err[2][9][3];

void extract_reso();
void read_sig();
void read_npart();
double calc_reso(double res);
double chi(double res);
double resEventPlane(double chi);

void downstream_analysis_bycen(int cen);

// void profile_divide_by_reso_fit_by_const(const char *profile_name, int cen);
void profile_divide_by_reso_corr_by_bkg(const char *profile_name, int cen);
// void profile_divide_by_reso(const char *profile_name, int cen, bool eff_corr);

void Q2_parent_results(int cen, int ep_option, int option112);

void plotting_TGraph_Helper(const char *title, const char *x_title, const char *y_title, int npts, float *x, float *x_err, float *y1, float *y1_err, float *y2 = NULL, float *y2_err = NULL, float *y3 = NULL, float *y3_err = NULL);

int min_cen = 3;
int max_cen = 5;

void downstream_analysis()
{
    read_sig();
    read_npart();
    read_purity_Q2();

    name_options3.push_back("TPC");
    name_options3.push_back("EPD");
    name_options3.push_back("EPD1");

    result_file.open(TString::Format("./%s_Results/output_v2.txt", method_name).Data());
    output_File2 = new TFile(TString::Format("./%s_Results/output_cen%d-%d.root", method_name, min_cen, max_cen).Data(), "recreate");

    for (int i = min_cen; i <= max_cen; i++)
    {
        cout << "Centrality " << i << endl;
        downstream_analysis_bycen(i);
    }

    float v2_averaged_TPC[9] = {0.}, v2_averaged_EPD[9] = {0.}, v2_averaged_EPD1[9] = {0.}, v2_averaged_TPC_err[9] = {0.}, v2_averaged_EPD_err[9] = {0.}, v2_averaged_EPD1_err[9] = {0.};

    for (int w = 0; w < 3; w++)
    {
        result_file << "const float v2_averaged_" << name_options3.at(w) << "[9] = {";

        for (int q = 0; q <= (max_cen - min_cen); q++)
        {
            // if(q != max_cen) result_file << alltrks_v2_vect[w].at(q) << ", ";
            // else result_file << alltrks_v2_vect[w].at(q) << "}; \n";

            // if(q != max_cen) result_file << alltrks_v2_vect[w].at(2*q) << " +/- " << alltrks_v2_vect[w].at(2*q + 1) << ", ";
            // else result_file << alltrks_v2_vect[w].at(2*q) << " +/- " << alltrks_v2_vect[w].at(2*q + 1) << "}; \n";

            if (q != max_cen)
                result_file << alltrks_v2_vect[w].at(2 * q) << ", ";
            else
                result_file << alltrks_v2_vect[w].at(2 * q) << "}; \n";

            if (w == 0)
            {
                v2_averaged_TPC[q] = alltrks_v2_vect[w].at(2 * q);
                v2_averaged_TPC_err[q] = alltrks_v2_vect[w].at(2 * q + 1);
            }
            else if (w == 1)
            {
                v2_averaged_EPD[q] = alltrks_v2_vect[w].at(2 * q);
                v2_averaged_EPD_err[q] = alltrks_v2_vect[w].at(2 * q + 1);
            }
            else if (w == 2)
            {
                v2_averaged_EPD1[q] = alltrks_v2_vect[w].at(2 * q);
                v2_averaged_EPD1_err[q] = alltrks_v2_vect[w].at(2 * q + 1);
            }
        }

        result_file << "const float v2_parent_averaged_" << name_options3.at(w) << "[9] = {";

        for (int q = 0; q <= (max_cen - min_cen); q++)
        {
            // if(q != max_cen) result_file << parent_v2_vect[w].at(q) << ", ";
            // else result_file << parent_v2_vect[w].at(q) << "}; \n";
            if (q != max_cen)
                result_file << parent_v2_vect[w].at(2 * q) << ", ";
            else
                result_file << parent_v2_vect[w].at(2 * q) << "}; \n";
        }

        result_file << "const float v2_averaged_pairpion_" << name_options3.at(w) << "[9] = {";

        for (int q = 0; q <= (max_cen - min_cen); q++)
        {
            // if(q != max_cen) result_file << parent_v2_vect[w].at(q) << ", ";
            // else result_file << parent_v2_vect[w].at(q) << "}; \n";
            if (q != max_cen)
                result_file << pionpion_v2_vect[w].at(2 * q) << ", ";
            else
                result_file << pionpion_v2_vect[w].at(2 * q) << "}; \n";
        }
    }

    TGraphErrors *reso_TPC_graph = new TGraphErrors(9, cen9_eff, reso_TPC, cen9_err_eff, reso_TPC_err);
    TGraphErrors *reso_EPD_graph = new TGraphErrors(9, cen9_eff, reso_EPD, cen9_err_eff, reso_EPD_err);
    TGraphErrors *reso_EPD1_graph = new TGraphErrors(9, cen9_eff, reso_EPD1, cen9_err_eff, reso_EPD1_err);

    TGraphErrors *v2_averaged_TPC_graph = new TGraphErrors(9, cen9_eff, v2_averaged_TPC, cen9_err_eff, v2_averaged_TPC_err);
    TGraphErrors *v2_averaged_EPD_graph = new TGraphErrors(9, cen9_eff, v2_averaged_EPD, cen9_err_eff, v2_averaged_EPD_err);
    TGraphErrors *v2_averaged_EPD1_graph = new TGraphErrors(9, cen9_eff, v2_averaged_EPD1, cen9_err_eff, v2_averaged_EPD1_err);

    output_File2->cd();
    TCanvas c1("resolution_graphs", "Resolution", 1500, 800);
    reso_TPC_graph->Draw("ALP");
    reso_EPD_graph->SetLineColor(kRed);
    reso_EPD_graph->Draw("SAMES LP");
    reso_EPD1_graph->SetLineColor(kBlue);
    reso_EPD1_graph->Draw("SAMES LP");
    c1.Write();

    TCanvas c2("v2_graphs", "v2", 1500, 800);
    v2_averaged_TPC_graph->Draw("ALP");
    v2_averaged_EPD_graph->SetLineColor(kRed);
    v2_averaged_EPD_graph->Draw("SAMES LP");
    v2_averaged_EPD1_graph->SetLineColor(kBlue);
    v2_averaged_EPD1_graph->Draw("SAMES LP");
    c2.Write();

    float gamma112_TPC[9] = {0.}, gamma112_EPD[9] = {0.}, gamma112_EPD1[9] = {0.};
    float gamma112_TPC_err[9] = {0.}, gamma112_EPD_err[9] = {0.}, gamma112_EPD1_err[9] = {0.};
    float gamma132_TPC[9] = {0.}, gamma132_EPD[9] = {0.}, gamma132_EPD1[9] = {0.};
    float gamma132_TPC_err[9] = {0.}, gamma132_EPD_err[9] = {0.}, gamma132_EPD1_err[9] = {0.};

    float gamma112_ensemble_TPC[9] = {0.}, gamma112_ensemble_EPD[9] = {0.}, gamma112_ensemble_EPD1[9] = {0.};
    float gamma112_ensemble_TPC_err[9] = {0.}, gamma112_ensemble_EPD_err[9] = {0.}, gamma112_ensemble_EPD1_err[9] = {0.};
    float gamma132_ensemble_TPC[9] = {0.}, gamma132_ensemble_EPD[9] = {0.}, gamma132_ensemble_EPD1[9] = {0.};
    float gamma132_ensemble_TPC_err[9] = {0.}, gamma132_ensemble_EPD_err[9] = {0.}, gamma132_ensemble_EPD1_err[9] = {0.};

    for (int j = min_cen; j <= max_cen; j++)
    {
        gamma112_TPC[j] = gamma[0][j][0] * npart[j];
        gamma112_EPD[j] = gamma[0][j][1] * npart[j];
        gamma112_EPD1[j] = gamma[0][j][2] * npart[j];
        gamma132_TPC[j] = gamma[1][j][0] * npart[j];
        gamma132_EPD[j] = gamma[1][j][1] * npart[j];
        gamma132_EPD1[j] = gamma[1][j][2] * npart[j];

        gamma112_TPC_err[j] = gamma_err[0][j][0] * npart[j];
        gamma112_EPD_err[j] = gamma_err[0][j][1] * npart[j];
        gamma112_EPD1_err[j] = gamma_err[0][j][2] * npart[j];
        gamma132_TPC_err[j] = gamma_err[1][j][0] * npart[j];
        gamma132_EPD_err[j] = gamma_err[1][j][1] * npart[j];
        gamma132_EPD1_err[j] = gamma_err[1][j][2] * npart[j];

        gamma112_ensemble_TPC[j] = gamma_ensemble[0][j][0] * npart[j];
        gamma112_ensemble_EPD[j] = gamma_ensemble[0][j][1] * npart[j];
        gamma112_ensemble_EPD1[j] = gamma_ensemble[0][j][2] * npart[j];
        gamma132_ensemble_TPC[j] = gamma_ensemble[1][j][0] * npart[j];
        gamma132_ensemble_EPD[j] = gamma_ensemble[1][j][1] * npart[j];
        gamma132_ensemble_EPD1[j] = gamma_ensemble[1][j][2] * npart[j];

        gamma112_ensemble_TPC_err[j] = gamma_ensemble_err[0][j][0] * npart[j];
        gamma112_ensemble_EPD_err[j] = gamma_ensemble_err[0][j][1] * npart[j];
        gamma112_ensemble_EPD1_err[j] = gamma_ensemble_err[0][j][2] * npart[j];
        gamma132_ensemble_TPC_err[j] = gamma_ensemble_err[1][j][0] * npart[j];
        gamma132_ensemble_EPD_err[j] = gamma_ensemble_err[1][j][1] * npart[j];
        gamma132_ensemble_EPD1_err[j] = gamma_ensemble_err[1][j][2] * npart[j];
    }

    // cout << "gamma_ensemble_err[0][0][0] = " << gamma_ensemble_err[0][0][0] << endl;

    plotting_TGraph_Helper("Compare ESE vs. Ensemble TPC #Delta#gamma_{112} * N_{part}", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_TPC, gamma112_TPC_err, gamma112_ensemble_TPC, gamma112_ensemble_TPC_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD #Delta#gamma_{112} * N_{part}", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD, gamma112_EPD_err, gamma112_ensemble_EPD, gamma112_ensemble_EPD_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD1 #Delta#gamma_{112} * N_{part}", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD1, gamma112_EPD1_err, gamma112_ensemble_EPD1, gamma112_ensemble_EPD1_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble TPC #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_TPC, gamma132_TPC_err, gamma132_ensemble_TPC, gamma132_ensemble_TPC_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_EPD, gamma132_EPD_err, gamma132_ensemble_EPD, gamma132_ensemble_EPD_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD1 #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_EPD1, gamma132_EPD1_err, gamma132_ensemble_EPD1, gamma132_ensemble_EPD1_err);

    TGraphErrors *gamma112_TPC_graph = new TGraphErrors(9, cen9_eff, gamma112_TPC, cen9_err_eff, gamma112_TPC_err);
    TGraphErrors *gamma112_EPD_graph = new TGraphErrors(9, cen9_eff, gamma112_EPD, cen9_err_eff, gamma112_EPD_err);
    TGraphErrors *gamma112_EPD1_graph = new TGraphErrors(9, cen9_eff, gamma112_EPD1, cen9_err_eff, gamma112_EPD1_err);
    TGraphErrors *gamma132_TPC_graph = new TGraphErrors(9, cen9_eff, gamma132_TPC, cen9_err_eff, gamma132_TPC_err);
    TGraphErrors *gamma132_EPD_graph = new TGraphErrors(9, cen9_eff, gamma132_EPD, cen9_err_eff, gamma132_EPD_err);
    TGraphErrors *gamma132_EPD1_graph = new TGraphErrors(9, cen9_eff, gamma132_EPD1, cen9_err_eff, gamma132_EPD1_err);

    TGraphErrors *gamma112_ensemble_TPC_graph = new TGraphErrors(9, cen9_eff, gamma112_ensemble_TPC, cen9_err_eff, gamma112_ensemble_TPC_err);
    TGraphErrors *gamma112_ensemble_EPD_graph = new TGraphErrors(9, cen9_eff, gamma112_ensemble_EPD, cen9_err_eff, gamma112_ensemble_EPD_err);
    TGraphErrors *gamma112_ensemble_EPD1_graph = new TGraphErrors(9, cen9_eff, gamma112_ensemble_EPD1, cen9_err_eff, gamma112_ensemble_EPD1_err);
    TGraphErrors *gamma132_ensemble_TPC_graph = new TGraphErrors(9, cen9_eff, gamma132_ensemble_TPC, cen9_err_eff, gamma132_ensemble_TPC_err);
    TGraphErrors *gamma132_ensemble_EPD_graph = new TGraphErrors(9, cen9_eff, gamma132_ensemble_EPD, cen9_err_eff, gamma132_ensemble_EPD_err);
    TGraphErrors *gamma132_ensemble_EPD1_graph = new TGraphErrors(9, cen9_eff, gamma132_ensemble_EPD1, cen9_err_eff, gamma132_ensemble_EPD1_err);

    TLine l1(0, 0, 80, 0);

    TCanvas c3("TPC", "TPC", 1500, 800);
    gamma112_TPC_graph->SetMarkerStyle(kFullCircle);
    gamma112_TPC_graph->GetXaxis()->SetTitle("Centrality %");
    gamma112_TPC_graph->GetYaxis()->SetTitle("#Delta#gamma * N_{part}");
    gamma112_TPC_graph->Draw("AP");
    gamma132_TPC_graph->SetMarkerStyle(kFullCircle);
    gamma132_TPC_graph->SetMarkerColor(kRed);
    gamma132_TPC_graph->SetLineColor(kRed);
    gamma132_TPC_graph->Draw("SAMES P");
    l1.Draw();
    c3.Write();

    TCanvas c4("EPD", "EPD", 1500, 800);
    gamma112_EPD_graph->SetMarkerStyle(kFullCircle);
    gamma112_EPD_graph->GetXaxis()->SetTitle("Centrality %");
    gamma112_EPD_graph->GetYaxis()->SetTitle("#Delta#gamma * N_{part}");
    gamma112_EPD_graph->Draw("AP");
    gamma132_EPD_graph->SetMarkerStyle(kFullCircle);
    gamma132_EPD_graph->SetMarkerColor(kRed);
    gamma132_EPD_graph->SetLineColor(kRed);
    gamma132_EPD_graph->Draw("SAMES P");
    l1.Draw();
    c4.Write();

    TCanvas c5("EPD1", "EPD1", 1500, 800);
    gamma112_EPD1_graph->SetMarkerStyle(kFullCircle);
    gamma112_EPD1_graph->GetXaxis()->SetTitle("Centrality %");
    gamma112_EPD1_graph->GetYaxis()->SetTitle("#Delta#gamma * N_{part}");
    gamma112_EPD1_graph->Draw("AP");
    gamma132_EPD1_graph->SetMarkerStyle(kFullCircle);
    gamma132_EPD1_graph->SetMarkerColor(kRed);
    gamma132_EPD1_graph->SetLineColor(kRed);
    gamma132_EPD1_graph->Draw("SAMES P");
    l1.Draw();
    c5.Write();

    TCanvas c6("Gamma112", "Gamma112", 1500, 800);
    gamma112_TPC_graph->Draw("AP");
    gamma112_TPC_graph->GetYaxis()->SetTitle("#Delta#gamma_{112} * N_{part}");
    gamma112_EPD_graph->SetMarkerColor(kRed);
    gamma112_EPD_graph->SetLineColor(kRed);
    gamma112_EPD_graph->Draw("SAMES P");
    gamma112_EPD1_graph->SetMarkerColor(kBlue);
    gamma112_EPD1_graph->SetLineColor(kBlue);
    gamma112_EPD1_graph->Draw("SAMES P");
    l1.Draw();
    c6.Write();

    TCanvas c7("Gamma132", "Gamma132", 1500, 800);
    gamma132_TPC_graph->SetMarkerColor(kBlack);
    gamma132_TPC_graph->SetLineColor(kBlack);
    gamma132_TPC_graph->GetXaxis()->SetTitle("Centrality %");
    gamma132_TPC_graph->GetYaxis()->SetTitle("#Delta#gamma_{132} * N_{part}");
    gamma132_TPC_graph->Draw("AP");
    gamma132_EPD_graph->SetMarkerColor(kRed);
    gamma132_EPD_graph->SetLineColor(kRed);
    gamma132_EPD_graph->Draw("SAMES P");
    gamma132_EPD1_graph->SetMarkerColor(kBlue);
    gamma132_EPD1_graph->SetLineColor(kBlue);
    gamma132_EPD1_graph->Draw("SAMES P");
    l1.Draw();
    c7.Write();

    output_File2->Close();
    result_file.close();
}

void plotting_TGraph_Helper(const char *title, const char *x_title, const char *y_title, int npts, float *x, float *x_err, float *y1, float *y1_err, float *y2, float *y2_err, float *y3, float *y3_err)
{
    TGraphErrors *tmp_graph1 = new TGraphErrors(npts, x, y1, x_err, y1_err);
    // cout << "y1_err[0] = " << y1_err[0] << endl;
    // cout << "y2_err[0] = " << y2_err[0] << endl;
    if (y2 != NULL)
        TGraphErrors *tmp_graph2 = new TGraphErrors(npts, x, y2, x_err, y2_err);
    if (y3 != NULL)
        TGraphErrors *tmp_graph3 = new TGraphErrors(npts, x, y3, x_err, y3_err);

    output_File2->cd();
    TCanvas c1(title, title, 1500, 800);
    tmp_graph1->SetMarkerStyle(kFullCircle);
    tmp_graph1->SetMarkerColor(kBlack);
    tmp_graph1->SetLineColor(kBlack);
    tmp_graph1->SetTitle(title);
    tmp_graph1->GetXaxis()->SetTitle(x_title);
    tmp_graph1->GetYaxis()->SetTitle(y_title);
    tmp_graph1->Draw("AP");
    if (y2 != NULL)
    {
        // cout << "plotting second" << endl;
        tmp_graph2->SetMarkerStyle(kFullCircle);
        tmp_graph2->SetMarkerColor(kRed);
        tmp_graph2->SetLineColor(kRed);
        tmp_graph2->Draw("SAMES P");
    }
    if (y3 != NULL)
    {
        tmp_graph3->SetMarkerStyle(kFullCircle);
        tmp_graph3->SetMarkerColor(kBlue);
        tmp_graph3->SetLineColor(kBlue);
        tmp_graph3->Draw("SAMES P");
    }
    TLine l1(0, 0, 80, 0);
    l1.Draw();
    c1.Write();
}

void downstream_analysis_bycen(int cen)
{
    file = new TFile(TString::Format("./%s/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", method_name, cen).Data());
    output_File = new TFile(TString::Format("./%s_Results/output_cen%d.root", method_name, cen).Data(), "recreate");

    extract_reso();

    reso_TPC[cen] = reso[0];
    reso_EPD[cen] = reso[1];
    reso_EPD1[cen] = reso[2];

    //     // profile_divide_by_reso_fit_by_const("Hist_v2parent_eta_obs5", cen);
    //     // profile_divide_by_reso_fit_by_const("Hist_v2parent_eta_obs5_rot", cen);

    profile_divide_by_reso("Parity_int_obs", cen, false);
    profile_divide_by_reso("Parity_int_obs", cen, true);
    profile_divide_by_reso("Parity_int_ss_obs", cen, false);
    profile_divide_by_reso("Parity_int_ss_obs", cen, true);
    profile_divide_by_reso("Delta_int_ss_obs", cen, false);
    profile_divide_by_reso("Delta_int_ss_obs", cen, true);

    char fname[200];
    for (int jk = 0; jk < 3; jk++)
    {
        sprintf(fname, "Hist_v2parent_pt_%s_obs5", name_options3.at(jk).Data());
        vector<float> temp_vec = Rebin_v2_Data(fname, fname, cen, jk);
        // parent_v2_vect[jk].push_back(temp_vec.at(0));
        parent_v2_vect[jk].push_back(temp_vec.at(1) / 100.0);
        parent_v2_vect[jk].push_back(temp_vec.at(2) / 100.0);
        temp_vec.clear();

        sprintf(fname, "Hist_v2_pt_obs2_%s_alltrks", name_options3.at(jk).Data());
        temp_vec = Rebin_v2_Data(fname, fname, cen, jk);
        // alltrks_v2_vect[jk].push_back(temp_vec.at(0));
        alltrks_v2_vect[jk].push_back(temp_vec.at(1) / 100.0);
        alltrks_v2_vect[jk].push_back(temp_vec.at(2) / 100.0);
        temp_vec.clear();

        sprintf(fname, "Hist_v2pion_pt_%s_obs5", name_options3.at(jk).Data());
        temp_vec = Rebin_v2_Data(fname, fname, cen, jk);
        // alltrks_v2_vect[jk].push_back(temp_vec.at(0));
        pionpion_v2_vect[jk].push_back(temp_vec.at(1) / 100.0);
        pionpion_v2_vect[jk].push_back(temp_vec.at(2) / 100.0);
        temp_vec.clear();
    }
    //     // vector<float> temp_vec = Rebin_v2_Data("Hist_v2parent_pt_obs5", "Hist_v2parent_pt_obs5", cen);
    //     // result_file << "Hist_v2parent_pt_obs5 cen " << cen << ": " << temp_vec.at(1) << ", " << temp_vec.at(2) << "\n";
    //     // temp_vec.clear();
    //     temp_vec = Rebin_v2_Data("Hist_v2parent_pt_obs5_rot", "Hist_v2parent_pt_obs5_rot", cen);
    //     result_file << "Hist_v2parent_pt_obs5_rot cen " << cen << ": " << temp_vec.at(1) << ", " << temp_vec.at(2) << "\n";
    //     temp_vec.clear();

    //     vector<float> temp_vec = Rebin_v2_Data("Hist_v2_pt_obs2_alltrks", "Hist_v2_pt_obs1_alltrks", cen);
    //     result_file << "Hist_v2_pt_obs2_alltrks cen " << cen << ": " << temp_vec.at(1) << ", " << temp_vec.at(2) << "\n";

    // profile_divide_by_reso_corr_by_bkg("Parity_int_obs", cen, false);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_ss_obs", cen, false);
    // profile_divide_by_reso_corr_by_bkg("Delta_int_ss_obs", cen, false);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_obs", cen, true);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_ss_obs", cen, true);
    // profile_divide_by_reso_corr_by_bkg("Delta_int_ss_obs", cen, true);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_obs3_splitpt_QQcut", cen);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_ss_obs3_splitpt_QQcut", cen);
    // profile_divide_by_reso_corr_by_bkg("Delta_int_ss_obs3_splitpt_QQcut", cen);

    for (int j = 0; j < 3; j++)
    {
        Q2_parent_results(cen, j, 0);
        Q2_parent_results(cen, j, 1);
        cout << "done with " << j << endl;
    }

    file->Close();
    output_File->Close();
}

// void profile_divide_by_reso_fit_by_const(const char *profile_name, int cen)
// {
//     cout << profile_name << endl;

//     TProfile *temp_profile = (TProfile *)file->Get(profile_name);

//     temp_profile->Scale(1.0 / reso);
//     // cout << "reso = " << reso << endl;
//     TF1 *temp_func = new TF1("temp_func", "[0]", -5, 5);
//     temp_profile->Fit(temp_func, "QR");
//     result_file << profile_name << " Fit Parameter = " << temp_func->GetParameter(0) << "\n";

//     output_File->cd();
//     temp_profile->Write(profile_name);
// }

void profile_divide_by_reso_corr_by_bkg(const char *profile_name, int cen, bool qq)
{
    cout << profile_name << endl;

    TString qqornot;

    if (qq)
        qqornot = "splitpt_QQcut";
    else
        qqornot = "splitpt";

    TProfile *result_profile = (TProfile *)file->Get(TString::Format("%s3_%s_0", profile_name, qqornot.Data()).Data());
    TProfile *first_bkg = (TProfile *)file->Get(TString::Format("%s3_%s_rot_0", profile_name, qqornot.Data()).Data());
    result_profile->Sumw2();
    first_bkg->Sumw2();
    result_profile->Add(result_profile, first_bkg, 1.0 / lam_purity[cen][0], -(1.0 - lam_purity[cen][0]) / lam_purity[cen][0]);

    TProfile *result_profile_err = (TProfile *)file->Get(TString::Format("%s1_%s_0", profile_name, qqornot.Data()).Data());
    TProfile *first_bkg_err = (TProfile *)file->Get(TString::Format("%s1_%s_rot_0", profile_name, qqornot.Data()).Data());
    result_profile_err->Sumw2();
    first_bkg_err->Sumw2();
    result_profile_err->Add(result_profile_err, first_bkg_err, 1.0 / lam_purity[cen][0], -(1.0 - lam_purity[cen][0]) / lam_purity[cen][0]);

    for (int i = 1; i < 15; i++)
    {
        TProfile *temp_profile = (TProfile *)file->Get(TString::Format("%s3_%s_%d", profile_name, qqornot.Data(), i).Data());
        TProfile *temp_profile_bkg = (TProfile *)file->Get(TString::Format("%s3_%s_rot_%d", profile_name, qqornot.Data(), i).Data());
        temp_profile->Sumw2();
        temp_profile_bkg->Sumw2();
        temp_profile->Add(temp_profile, temp_profile_bkg, 1.0 / lam_purity[cen][i], -(1.0 - lam_purity[cen][i]) / lam_purity[cen][i]);

        result_profile->Add((TProfile *)temp_profile->Clone());

        TProfile *temp_profile_err = (TProfile *)file->Get(TString::Format("%s1_%s_%d", profile_name, qqornot.Data(), i).Data());
        TProfile *temp_profile_bkg_err = (TProfile *)file->Get(TString::Format("%s1_%s_rot_%d", profile_name, qqornot.Data(), i).Data());
        temp_profile_err->Sumw2();
        temp_profile_bkg_err->Sumw2();
        temp_profile_err->Add(temp_profile_err, temp_profile_bkg_err, 1.0 / lam_purity[cen][i], -(1.0 - lam_purity[cen][i]) / lam_purity[cen][i]);

        result_profile_err->Add((TProfile *)temp_profile_err->Clone());
    }

    result_profile->Scale(1.0 / reso1);
    result_profile_err->Scale(1.0 / reso1);
    int n_bins = result_profile->GetNbinsX();
    double x[10] = {0.}, x_err[10] = {0.}, y[10] = {0.}, y_err[10] = {0.};
    // cout << "reso1 = " << reso1 << endl;

    for (int k = 1; k <= n_bins; k++)
    {
        x[k - 1] = k;
        y[k - 1] = result_profile->GetBinContent(k);
        y_err[k - 1] = result_profile_err->GetBinError(k);
    }

    TGraphErrors *result_graph = new TGraphErrors(n_bins, x, y, x_err, y_err);

    TProfile *result_profile_anti = (TProfile *)file->Get(TString::Format("%s3_%s_anti_0", profile_name, qqornot.Data()).Data());
    TProfile *first_bkg_anti = (TProfile *)file->Get(TString::Format("%s3_%s_anti_rot_0", profile_name, qqornot.Data()).Data());
    result_profile_anti->Sumw2();
    first_bkg_anti->Sumw2();
    result_profile_anti->Add(result_profile_anti, first_bkg_anti, 1.0 / antilam_purity[cen][0], -(1.0 - antilam_purity[cen][0]) / antilam_purity[cen][0]);

    TProfile *result_profile_anti_err = (TProfile *)file->Get(TString::Format("%s1_%s_anti_0", profile_name, qqornot.Data()).Data());
    TProfile *first_bkg_anti_err = (TProfile *)file->Get(TString::Format("%s1_%s_anti_rot_0", profile_name, qqornot.Data()).Data());
    result_profile_anti_err->Sumw2();
    first_bkg_anti_err->Sumw2();
    result_profile_anti_err->Add(result_profile_anti_err, first_bkg_anti_err, 1.0 / antilam_purity[cen][0], -(1.0 - antilam_purity[cen][0]) / antilam_purity[cen][0]);

    for (int i = 1; i < 15; i++)
    {
        TProfile *temp_profile = (TProfile *)file->Get(TString::Format("%s3_%s_anti_%d", profile_name, qqornot.Data(), i).Data());
        TProfile *temp_profile_bkg = (TProfile *)file->Get(TString::Format("%s3_%s_anti_rot_%d", profile_name, qqornot.Data(), i).Data());
        temp_profile->Sumw2();
        temp_profile_bkg->Sumw2();
        temp_profile->Add(temp_profile, temp_profile_bkg, 1.0 / antilam_purity[cen][i], -(1.0 - antilam_purity[cen][i]) / antilam_purity[cen][i]);

        result_profile_anti->Add((TProfile *)temp_profile->Clone());

        TProfile *temp_profile_err = (TProfile *)file->Get(TString::Format("%s1_%s_anti_%d", profile_name, qqornot.Data(), i).Data());
        TProfile *temp_profile_bkg_err = (TProfile *)file->Get(TString::Format("%s1_%s_anti_rot_%d", profile_name, qqornot.Data(), i).Data());
        temp_profile_err->Sumw2();
        temp_profile_bkg_err->Sumw2();
        temp_profile_err->Add(temp_profile_err, temp_profile_bkg_err, 1.0 / antilam_purity[cen][i], -(1.0 - antilam_purity[cen][i]) / antilam_purity[cen][i]);

        result_profile_anti_err->Add((TProfile *)temp_profile->Clone());
    }

    result_profile_anti->Scale(1.0 / reso1);
    result_profile_anti_err->Scale(1.0 / reso1);
    int n_bins_anti = result_profile_anti->GetNbinsX();
    double x_anti[10] = {0.}, x_anti_err[10] = {0.}, y_anti[10] = {0.}, y_anti_err[10] = {0.};
    // cout << "reso1 = " << reso1 << endl;

    for (int j = 1; j <= n_bins_anti; j++)
    {
        x_anti[j - 1] = j;
        y_anti[j - 1] = result_profile_anti->GetBinContent(j);
        y_anti_err[j - 1] = result_profile_anti_err->GetBinError(j);
    }

    TGraphErrors *result_graph_anti = new TGraphErrors(n_bins_anti, x_anti, y_anti, x_anti_err, y_anti_err);

    output_File->cd();
    result_graph->Write(TString::Format("%s3_%s_bkgcorr", profile_name, qqornot.Data()).Data());
    result_graph_anti->Write(TString::Format("%s3_%s_anti_bkgcorr", profile_name, qqornot.Data()).Data());

    result_profile->Add(result_profile_anti);
    result_profile_err->Add(result_profile_anti_err);
    double x_combined[10] = {0.}, x_combined_err[10] = {0.}, y_combined[10] = {0.}, y_combined_err[10] = {0.};
    // cout << "reso1 = " << reso1 << endl;

    for (int k = 1; k <= n_bins; k++)
    {
        x_combined[k - 1] = k;
        y_combined[k - 1] = result_profile->GetBinContent(k);
        y_combined_err[k - 1] = result_profile_err->GetBinError(k);
    }

    TGraphErrors *result_graph_combined = new TGraphErrors(n_bins, x_combined, y_combined, x_combined_err, y_combined_err);
    result_graph_combined->Write(TString::Format("%s3_%s_combined_bkgcorr", profile_name, qqornot.Data()).Data());
}

void profile_divide_by_reso(const char *profile_name, int cen, bool eff_corr)
{
    cout << profile_name << endl;
    const int n_entries = 25;

    int num = 0;
    if (eff_corr)
        num = 3;
    else
        num = 1;

    TProfile *temp_profile = (TProfile *)file->Get(TString::Format("%s%d", profile_name, num));
    TProfile *temp_profile_rot = (TProfile *)file->Get(TString::Format("%s%d_rot", profile_name, num));
    int n_bins = temp_profile->GetNbinsX();
    double x[n_entries] = {0.}, x_err[n_entries] = {0.}, y[n_entries] = {0.}, y_err[n_entries] = {0.};
    // cout << "Sumw2() checks 1" << endl;
    temp_profile->Sumw2();
    temp_profile_rot->Sumw2();

    // temp_profile->Add(temp_profile, temp_profile_rot, 1.0 / lam_purity[cen][17], -(1.0 - lam_purity[cen][17]) / lam_purity[cen][17]);

    // cout << "lam_purity[cen][17] = " << lam_purity[cen][17] << endl;

    // temp_profile->Scale(1.0 / reso1);

    TProfile *temp_profile_err = (TProfile *)file->Get(TString::Format("%s%d", profile_name, 1));
    TProfile *temp_profile_err_rot = (TProfile *)file->Get(TString::Format("%s%d_rot", profile_name, 1));
    // cout << "Sumw2() checks 2" << endl;
    // temp_profile_err->Sumw2();
    // temp_profile_err_rot->Sumw2();

    if (eff_corr)
    {
        // temp_profile_err->Add(temp_profile_err, temp_profile_err_rot, 1.0 / lam_purity[cen][17], -(1.0 - lam_purity[cen][17]) / lam_purity[cen][17]);

        // temp_profile_err->Scale(1.0 / reso1);

        for (int i = 1; i <= temp_profile->GetNbinsX(); i++)
        {
            x[i - 1] = i;
            y[i - 1] = temp_profile->GetBinContent(i);
            y_err[i - 1] = temp_profile_err->GetBinError(i);

            // cout << "y[" << i - 1 << "] = " << y[i - 1] << endl;
        }
    }
    else
    {
        for (int i = 1; i <= temp_profile->GetNbinsX(); i++)
        {
            x[i - 1] = i;
            y[i - 1] = temp_profile->GetBinContent(i);
            y_err[i - 1] = temp_profile->GetBinError(i);

            // cout << "y[" << i - 1 << "] = " << y[i - 1] << endl;
        }
    }

    TGraphErrors *temp_graph = new TGraphErrors(n_bins, x, y, x_err, y_err);

    // TProfile *temp_profile_QQcut = (TProfile *)file->Get(TString::Format("%s%d_QQcut", profile_name, num));
    // TProfile *temp_profile_QQcut_rot = (TProfile *)file->Get(TString::Format("%s%d_QQcut_rot", profile_name, num));
    // int n_bins_QQcut = temp_profile_QQcut->GetNbinsX();
    // double x_QQcut[10] = {0.}, x_err_QQcut[10] = {0.}, y_QQcut[10] = {0.}, y_err_QQcut[10] = {0.};

    // temp_profile_QQcut->Sumw2();
    // temp_profile_QQcut_rot->Sumw2();

    // temp_profile_QQcut->Add(temp_profile_QQcut, temp_profile_QQcut_rot, 1.0 / lam_purity[cen][17], -(1.0 - lam_purity[cen][17]) / lam_purity[cen][17]);

    // temp_profile_QQcut->Scale(1.0 / reso1);

    // TProfile *temp_profile_err_QQcut = (TProfile *)file->Get(TString::Format("%s%d_QQcut", profile_name, 1));
    // TProfile *temp_profile_err_QQcut_rot = (TProfile *)file->Get(TString::Format("%s%d_QQcut_rot", profile_name, 1));
    // temp_profile_err_QQcut->Sumw2();
    // temp_profile_err_QQcut_rot->Scale(1.0 / reso1);

    // if (eff_corr)
    // {
    //     temp_profile_err_QQcut->Add(temp_profile_err_QQcut, temp_profile_err_QQcut_rot, 1.0 / lam_purity[cen][17], -(1.0 - lam_purity[cen][17]) / lam_purity[cen][17]);

    //     temp_profile_err_QQcut->Scale(1.0 / reso1);

    //     for (int i = 1; i <= temp_profile_QQcut->GetNbinsX(); i++)
    //     {
    //         x_QQcut[i - 1] = i;
    //         y_QQcut[i - 1] = temp_profile_QQcut->GetBinContent(i);
    //         y_err_QQcut[i - 1] = temp_profile_err_QQcut->GetBinError(i);
    //     }
    // }
    // else
    // {
    //     for (int i = 1; i <= temp_profile_QQcut->GetNbinsX(); i++)
    //     {
    //         x_QQcut[i - 1] = i;
    //         y_QQcut[i - 1] = temp_profile_QQcut->GetBinContent(i);
    //         y_err_QQcut[i - 1] = temp_profile_QQcut->GetBinError(i);
    //     }
    // }

    // TGraphErrors *temp_graph_QQcut = new TGraphErrors(n_bins_QQcut, x_QQcut, y_QQcut, x_err_QQcut, y_err_QQcut);

    // antilambda
    TProfile *temp_profile_anti = (TProfile *)file->Get(TString::Format("%s%d_anti", profile_name, num));
    TProfile *temp_profile_anti_rot = (TProfile *)file->Get(TString::Format("%s%d_anti_rot", profile_name, num));
    int n_bins_anti = temp_profile_anti->GetNbinsX();
    double x_anti[n_entries] = {0.}, x_anti_err[n_entries] = {0.}, y_anti[n_entries] = {0.}, y_anti_err[n_entries] = {0.};
    // cout << "Sumw2() checks 3" << endl;
    temp_profile_anti->Sumw2();
    temp_profile_anti_rot->Sumw2();

    // temp_profile_anti->Add(temp_profile_anti, temp_profile_anti_rot, 1.0 / antilam_purity[cen][17], -(1.0 - antilam_purity[cen][17]) / antilam_purity[cen][17]);

    // temp_profile_anti->Scale(1.0 / reso1);

    TProfile *temp_profile_anti_err = (TProfile *)file->Get(TString::Format("%s%d_anti", profile_name, 1));
    TProfile *temp_profile_err_anti_rot = (TProfile *)file->Get(TString::Format("%s%d_anti_rot", profile_name, 1));
    // cout << "Sumw2() checks 4" << endl;
    // temp_profile_anti_err->Sumw2();
    // temp_profile_err_anti_rot->Sumw2();

    if (eff_corr)
    {
        // temp_profile_anti_err->Add(temp_profile_anti_err, temp_profile_err_anti_rot, 1.0 / antilam_purity[cen][17], -(1.0 - antilam_purity[cen][17]) / antilam_purity[cen][17]);

        // temp_profile_anti_err->Scale(1.0 / reso1);

        for (int i = 1; i <= temp_profile_anti->GetNbinsX(); i++)
        {
            x_anti[i - 1] = i;
            y_anti[i - 1] = temp_profile_anti->GetBinContent(i);
            y_anti_err[i - 1] = temp_profile_anti_err->GetBinError(i);
        }
    }
    else
    {
        for (int i = 1; i <= temp_profile_anti->GetNbinsX(); i++)
        {
            x_anti[i - 1] = i;
            y_anti[i - 1] = temp_profile_anti->GetBinContent(i);
            y_anti_err[i - 1] = temp_profile_anti->GetBinError(i);
        }
    }

    TGraphErrors *temp_graph_anti = new TGraphErrors(n_bins_anti, x_anti, y_anti, x_anti_err, y_anti_err);

    // TProfile *temp_profile_QQcut_anti = (TProfile *)file->Get(TString::Format("%s%d_QQcut_anti", profile_name, num));
    // TProfile *temp_profile_QQcut_anti_rot = (TProfile *)file->Get(TString::Format("%s%d_QQcut_anti_rot", profile_name, num));
    // int n_bins_QQcut_anti = temp_profile_QQcut_anti->GetNbinsX();
    // double x_QQcut_anti[10] = {0.}, x_err_QQcut_anti[10] = {0.}, y_QQcut_anti[10] = {0.}, y_err_QQcut_anti[10] = {0.};

    // temp_profile_QQcut_anti->Sumw2();
    // temp_profile_QQcut_anti_rot->Sumw2();

    // temp_profile_QQcut_anti->Add(temp_profile_QQcut_anti, temp_profile_QQcut_anti_rot, 1.0 / antilam_purity[cen][17], -(1.0 - antilam_purity[cen][17]) / antilam_purity[cen][17]);

    // temp_profile_QQcut_anti->Scale(1.0 / reso1);

    // TProfile *temp_profile_err_QQcut_anti = (TProfile *)file->Get(TString::Format("%s%d_QQcut_anti", profile_name, 1));
    // TProfile *temp_profile_err_QQcut_anti_rot = (TProfile *)file->Get(TString::Format("%s%d_QQcut_anti_rot", profile_name, 1));
    // temp_profile_err_QQcut_anti->Sumw2();
    // temp_profile_err_QQcut_anti_rot->Scale(1.0 / reso1);

    // if (eff_corr)
    // {
    //     temp_profile_err_QQcut_anti->Add(temp_profile_err_QQcut_anti, temp_profile_err_QQcut_anti_rot, 1.0 / antilam_purity[cen][17], -(1.0 - antilam_purity[cen][17]) / antilam_purity[cen][17]);

    //     temp_profile_err_QQcut_anti->Scale(1.0 / reso1);

    //     for (int i = 1; i <= temp_profile_QQcut_anti->GetNbinsX(); i++)
    //     {
    //         x_QQcut_anti[i - 1] = i;
    //         y_QQcut_anti[i - 1] = temp_profile_QQcut_anti->GetBinContent(i);
    //         y_err_QQcut_anti[i - 1] = temp_profile_err_QQcut_anti->GetBinError(i);
    //     }
    // }
    // else
    // {
    //     for (int i = 1; i <= temp_profile_QQcut_anti->GetNbinsX(); i++)
    //     {
    //         x_QQcut_anti[i - 1] = i;
    //         y_QQcut_anti[i - 1] = temp_profile_QQcut_anti->GetBinContent(i);
    //         y_err_QQcut_anti[i - 1] = temp_profile_QQcut_anti->GetBinError(i);
    //     }
    // }

    // TGraphErrors *temp_graph_QQcut_anti = new TGraphErrors(n_bins_QQcut_anti, x_QQcut_anti, y_QQcut_anti, x_err_QQcut_anti, y_err_QQcut_anti);

    temp_profile->Add(temp_profile_anti);
    double x_combined[n_entries] = {0.}, x_err_combined[n_entries] = {0.}, y_combined[n_entries] = {0.}, y_err_combined[n_entries] = {0.};
    if (eff_corr)
    {
        temp_profile_err->Add(temp_profile_anti_err);
        for (int i = 1; i <= temp_profile->GetNbinsX(); i++)
        {
            x_combined[i - 1] = i;
            y_combined[i - 1] = temp_profile->GetBinContent(i);
            y_err_combined[i - 1] = temp_profile_err->GetBinError(i);

            // cout << "y_combined[" << i - 1 << "] = " << y_combined[i - 1] << endl;
        }
    }
    else
    {
        for (int i = 1; i <= temp_profile->GetNbinsX(); i++)
        {
            x_combined[i - 1] = i;
            y_combined[i - 1] = temp_profile->GetBinContent(i);
            y_err_combined[i - 1] = temp_profile->GetBinError(i);
        }
    }

    TGraphErrors *temp_graph_combined = new TGraphErrors(n_bins, x_combined, y_combined, x_err_combined, y_err_combined);

    // cout << "n_bins = " << n_bins << endl;

    for (int tp1 = 0; tp1 < 3; tp1++)
    {
        if (profile_name == "Parity_int_obs")
        {
            gamma_ensemble[1][cen][tp1] = (float)((y_combined[4 + (8 * tp1)] - y_combined[3 + (8 * tp1)]) + (y_combined[8 + (8 * tp1)] - y_combined[7 + (8 * tp1)])) / 2.0;
            gamma_ensemble[1][cen][tp1] = gamma_ensemble[1][cen][tp1] / reso[tp1] / 100.0;

            float tmp_error1 = error_add(y_err_combined[4 + (8 * tp1)], y_err_combined[3 + (8 * tp1)]);
            float tmp_erorr2 = error_add(y_err_combined[8 + (8 * tp1)], y_err_combined[7 + (8 * tp1)]);
            gamma_ensemble_err[1][cen][tp1] = (float)error_add(tmp_error1, tmp_erorr2) / 2.0;
            gamma_ensemble_err[1][cen][tp1] = gamma_ensemble_err[1][cen][tp1] / reso[tp1] / 100.0;

            // cout << "gamma_ensemble[1][" << cen << "][" << tp1 << "] = " << gamma_ensemble[1][cen][tp1] << endl;
        }
        if (profile_name == "Parity_int_ss_obs")
        {
            // cout << "y_combined[" << 4 + (4 * tp1) << "] = " << y_combined[4 + (4 * tp1)] << endl;
            gamma_ensemble[0][cen][tp1] = y_combined[4 + (4 * tp1)] - y_combined[3 + (4 * tp1)];
            // cout << "gamma_ensemble[0][" << cen << "][" << tp1 << "] = " << gamma_ensemble[0][cen][tp1] << endl;
            gamma_ensemble[0][cen][tp1] = gamma_ensemble[0][cen][tp1] / reso[tp1] / 100.0;
            // cout << "reso1[" << tp1 << "] = " << reso1[tp1] << endl;
            // cout << "gamma_ensemble[0][" << cen << "][" << tp1 << "] = " << gamma_ensemble[0][cen][tp1] << endl;
            gamma_ensemble_err[0][cen][tp1] = error_add(y_err_combined[4 + (4 * tp1)], y_err_combined[3 + (4 * tp1)]);
            gamma_ensemble_err[0][cen][tp1] = gamma_ensemble_err[0][cen][tp1] / reso[tp1] / 100.0;
        }
    }

    // temp_profile_QQcut->Add(temp_profile_QQcut_anti);
    // double x_combined_QQcut[10] = {0.}, x_err_combined_QQcut[10] = {0.}, y_combined_QQcut[10] = {0.}, y_err_combined_QQcut[10] = {0.};
    // if (eff_corr)
    // {
    //     temp_profile_err_QQcut->Add(temp_profile_err_QQcut_anti);
    //     for (int i = 1; i <= temp_profile_QQcut->GetNbinsX(); i++)
    //     {
    //         x_combined_QQcut[i - 1] = i;
    //         y_combined_QQcut[i - 1] = temp_profile_QQcut->GetBinContent(i);
    //         y_err_combined_QQcut[i - 1] = temp_profile_err_QQcut->GetBinError(i);
    //     }
    // }
    // else{
    //     for (int i = 1; i <= temp_profile_QQcut->GetNbinsX(); i++)
    //     {
    //         x_combined_QQcut[i - 1] = i;
    //         y_combined_QQcut[i - 1] = temp_profile_QQcut->GetBinContent(i);
    //         y_err_combined_QQcut[i - 1] = temp_profile_QQcut->GetBinError(i);
    //     }
    // }

    // TGraphErrors *temp_graph_combined_QQcut = new TGraphErrors(n_bins_QQcut, x_combined_QQcut, y_combined_QQcut, x_err_combined_QQcut, y_err_combined_QQcut);

    output_File->cd();
    temp_graph->Write(TString::Format("%s%d", profile_name, num));
    // temp_graph_QQcut->Write(TString::Format("%s%d_QQcut", profile_name, num));
    temp_graph_anti->Write(TString::Format("%s%d_anti", profile_name, num));
    // temp_graph_QQcut_anti->Write(TString::Format("%s%d_QQcut_anti", profile_name, num));
    temp_graph_combined->Write(TString::Format("%s%d_combined", profile_name, num));
    // temp_graph_combined_QQcut->Write(TString::Format("%s%d_QQcut_combined", profile_name, num));
}

void extract_reso()
{
    file->cd();
    TProfile *Hist_cos = (TProfile *)file->Get("Hist_cos");
    TProfile *Hist_cos_EPD = (TProfile *)file->Get("Hist_cos_EPD");

    reso1[0] = calc_reso(Hist_cos->GetBinContent(2)); // for full event plane
    reso[0] = sqrt(Hist_cos->GetBinContent(1));       // for eta sub TPC
    if (Hist_cos_EPD->GetBinContent(1) >= 0)
        reso[1] = sqrt(Hist_cos_EPD->GetBinContent(1)); // for eta sub EPD
    else
        reso[1] = 1;
    reso[2] = Hist_cos_EPD->GetBinContent(4); // for eta sub EPD1

    for (int i = 0; i < 3; i++)
        cout << "reso[" << i << "] = " << reso[i] << endl;

    // reso_TPC.push_back(reso[0]);
    // reso_EPD.push_back(reso[1]);
    // reso_EPD1.push_back(reso[2]);

    delete Hist_cos;
    Hist_cos = NULL;
    delete Hist_cos_EPD;
    Hist_cos_EPD = NULL;
}

double calc_reso(double res)
{
    double chis = chi(sqrt(res));
    double final_res = resEventPlane(sqrt(2.) * chis);

    return final_res;
}

double chi(double res)
{
    double chi = 2.0;
    double delta = 1.0;

    for (int i = 0; i < 15; i++)
    {
        while (resEventPlane(chi) < res)
        {
            chi += delta;
        }
        delta = delta / 2.;
        while (resEventPlane(chi) > res)
        {
            chi -= delta;
        }
        delta = delta / 2.;
    }

    return chi;
}

//-----------------------------------------------------------------------
double resEventPlane(double chi)
{
    double con = 0.626657;
    double arg = chi * chi / 4.;

    Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

    return res;
}

void read_sig()
{
    std::fstream myfile(TString::Format("./Fit_and_subtract/%s/purity.txt", method_name).Data(), std::ios_base::in);

    float a;
    int count = 0;
    while (myfile >> a)
    {
        if (count < 306)
        {
            // cout << "a = " << a << endl;
            if (count % 2 == 0)
            {
                lam_purity[count / 34][(count - (count / 34) * 34) / 2] = a;
                // cout << "(count-(count / 34)*34) = " << (count-(count / 34)*34) << endl;
                // cout << "count / 34 = " << count / 34 << endl;
                // cout << "lam_purity[count / 34][(count-(count / 34)*34)/2] = " << lam_purity[count / 34][(count-(count / 34)*34)/2];
            }
            else if (count % 2 == 1)
            {
                lam_yield[count / 34][(count - (count / 34) * 34) / 2] = a;
                // cout << "(count-(count / 34)*34) = " << endl;
                // cout << ", lam_yield[count / 34][(count-(count / 34)*34)/2] = " << lam_yield[count / 34][(count-(count / 34)*34)/2] << endl;
            }
        }
        else if (count < 612)
        {
            if (count % 2 == 0)
            {
                antilam_purity[(count - 306) / 34][(count - (count / 34) * 34) / 2] = a;
                // cout << "antilam_purity[(count-306) / 34][(count-(count / 34)*34)/2] = " << antilam_purity[(count-306) / 34][(count-(count / 34)*34)/2];
            }
            else if (count % 2 == 1)
            {
                antilam_yield[(count - 306) / 34][(count - (count / 34) * 34) / 2] = a;
                // cout << ", antilam_yield[(count-306) / 34][(count-(count / 34)*34)/2] = " << antilam_yield[(count-306) / 34][(count-(count / 34)*34)/2] << endl;
            }
        }
        else if (count == 612)
            break;

        count++;
    }

    std::fstream secondfile(TString::Format("./Fit_and_subtract/%s/Without_ptsplit/purity.txt", method_name).Data(), std::ios_base::in);

    count = 0;
    while (secondfile >> a)
    {
        if (count < 18)
        {
            // cout << "a = " << a << endl;
            if (count % 2 == 0)
            {
                lam_purity[count / 2][17] = a;
                // if(debug_resolution == true) cout << "lam_purity[count / 2][17] = " << lam_purity[count / 2][17] << endl;
                // cout << "lam_purity[count / 2][17] = " << lam_purity[count / 2][17] << endl;
            }
            else if (count % 2 == 1)
            {
                lam_yield[count / 2][17] = a;
                // if(debug_resolution == true) cout << "lam_yield[count / 2][17] = " << lam_yield[count / 2][17] << endl;
            }
        }
        else if (count < 36)
        {
            if (count % 2 == 0)
            {
                antilam_purity[(count - 18) / 2][17] = a;
                // if(debug_resolution == true) cout << "antilam_purity[(count-18) / 2][17] = " << antilam_purity[(count - 18) / 2][17] << endl;
            }
            else if (count % 2 == 1)
            {
                antilam_yield[(count - 18) / 2][17] = a;
                // if(debug_resolution == true) cout << "antilam_yield[(count-18) / 2][17] = " << antilam_yield[(count - 18) / 2][17] << endl;
            }
        }
        else if (count == 36)
            break;

        count++;
    }
}

void Q2_parent_results(int cen, int ep_option, int option112)
{
    float Q2_range[3][9] = {{3.5, 5, 5, 5, 5, 5, 3.5, 1.5, 0.6}, {3.5, 3.5, 5, 5, 5, 5, 5, 5, 5}, {7.5, 5, 5, 5, 5, 5, 5, 5, 5}};
    float Q2_range_min[3][9] = {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}};
    // cout << "Q2_range_min[" << ep_option << "][" << cen << "] = " << Q2_range_min[ep_option][cen] << endl;

    int rebin[3][9] = {{5, 5, 5, 5, 5, 5, 5, 2, 1}, {5, 5, 5, 5, 5, 5, 5, 5, 5}, {5, 5, 5, 5, 5, 5, 5, 5, 5}};
    TString ep_option_name[3] = {"", "_EPD", "_EPD1"};

    // resolution
    char printname1[200];
    sprintf(printname1, "cen%d.Resolution.parent.%s.%d.pdf", cen, ep_option_name[ep_option].Data(), option112);

    // q2-g112-ssos
    double YRangeMax_ew_Q2[9] = {2, 1, 0.8, 0.4, 0.2, 0.1, 0.2, 0.2, 0.2};
    char printname4[200];
    sprintf(printname4, "cen%d.osss.Q2.parent.%s.%d.pdf", cen, ep_option_name[ep_option].Data(), option112);

    // DG112-v2
    double Yrange_delGamma[9] = {1, 0.4, 0.1, 0.1, 0.04, 0.02, 0.02, 0.02, 0.02};
    // double Xrange_delGamma[9] = {2, 0.25, 0.3, 0.22, 0.18, 0.15, 0.05, 0.03, 0.02};
    double Xrange_delGamma[9] = {1.2, 1.2, 0.6, 1.0, 0.5, 0.4, 0.4, 0.4, 0.3};
    // double Xfitting_delGamma[9] = {0.1, 0.1, 0.3, 0.22, 0.18, 0.15, 0.05, 0.03, 0.02};
    double Xfitting_delGamma[3][9] = {{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}, {0.4, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}, {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}};
    // double Xfitting_delGamma_min[2][3][9] = {{{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0.035, 0, 0, 0, 0, 0, 0, 0}, {0, 0.02, 0, 0, 0, 0, 0, 0, 0}},
    //                                          {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.01, 0, 0, 0, 0, 0, 0, 0, 0}}};
    double Xfitting_delGamma_min[2][3][9] = {{{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}},
                                             {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}}};
    char printname5[200];
    sprintf(printname5, "cen%d.ese.parent.Q2range%d.rebin%d.%s.%d.pdf", cen, Q2_range[ep_option][cen], rebin[ep_option][cen], ep_option_name[ep_option].Data(), option112);

    ////////////////Q2 parent results
    // TProfile *p_res = (TProfile *)file->Get("p_cos_Q2_parent");
    // p_cos_Q2, p_cos_Q2_EPD, p_cos_Q2_EPD1
    TProfile *p_res = (TProfile *)file->Get(TString::Format("p_cos_Q2%s", ep_option_name[ep_option].Data()).Data());
    TF1 *f_res = new TF1("f_res", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", Q2_range_min[ep_option][cen], Q2_range[ep_option][cen]);
    p_res->Fit("f_res", "0Q", "", Q2_range_min[ep_option][cen], Q2_range[ep_option][cen]);

    // p_v2parente_Q2parent_obs1, p_v2parente_Q2parent_obs2, p_v2parentw_Q2parent_obs1, p_v2parentw_Q2parent_obs2
    // p_v2parente_Q2parent_EPD_obs1, p_v2parente_Q2parent_EPD_obs2, p_v2parentw_Q2parent_EPD_obs1, p_v2parentw_Q2parent_EPD_obs2
    // p_v2parente_Q2parent_EPD1_obs1, p_v2parente_Q2parent_EPD1_obs2
    // TProfile *p_v2_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2parente_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2parente_Q2parent%s_obs2", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2parentw_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2parentw_Q2parent%s_obs2", ep_option_name[ep_option].Data()));

    TProfile *p_v2_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2e_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    TProfile *p_v2_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2e_Q2parent%s_obs2", ep_option_name[ep_option].Data()));
    TProfile *p_v2w_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2w_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    TProfile *p_v2w_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2w_Q2parent%s_obs2", ep_option_name[ep_option].Data()));

    // cout << "p_v2_Q_obs1 entries = " << p_v2_Q_obs1->GetEntries() << endl;
    // cout << "p_v2w_Q_obs1 entries = " << p_v2w_Q_obs1->GetEntries() << endl;

    p_v2_Q_obs1->Scale(1);
    if (ep_option <= 1)
        p_v2_Q_obs1->Add(p_v2w_Q_obs1);
    p_v2_Q_obs2->Scale(1);
    if (ep_option <= 1)
        p_v2_Q_obs2->Add(p_v2w_Q_obs2);

    // cout << "after addition p_v2_Q_obs1 entries = " << p_v2_Q_obs1->GetEntries() << endl;

    TH1 *p_v2_Q_obs_new = (TH1 *)p_v2_Q_obs1->ProjectionX(TString::Format("p_v2_Q_obs_new%s%d", ep_option_name[ep_option].Data(), option112));
    // cout << "p_v2_Q_obs_new (PROJECTION) entries = " << p_v2_Q_obs_new->GetEntries() << endl;
    int Nbin = p_v2_Q_obs1->GetNbinsX();
    int Nbin_min = p_v2_Q_obs1->GetXaxis()->FindBin(Q2_range_min[ep_option][cen]);

    for (int i = (Nbin_min - 1); i < Nbin; i++)
    {
        float cont = p_v2_Q_obs2->GetBinContent(i + 1);
        float err = p_v2_Q_obs1->GetBinError(i + 1);
        float res_q = 1;

        // if (p_v2_Q_obs2->GetBinCenter(i + 1) < Q2_range[cen])
        //     res_q = f_res->Eval(p_v2_Q_obs2->GetBinCenter(i + 1));
        // else
        //     res_q = f_res->Eval(Q2_range[cen]);

        if (p_v2_Q_obs2->GetBinCenter(i + 1) < Q2_range[ep_option][cen])
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(p_v2_Q_obs2->GetBinCenter(i + 1)));
        else
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]));

        if (ep_option <= 1)
            res_q = sqrt(res_q);

        p_v2_Q_obs_new->SetBinContent(i + 1, cont / res_q / 100.);
        p_v2_Q_obs_new->SetBinError(i + 1, err / res_q / 100.);

        if (cen == 8 && ep_option == 0 && option112 == 0 && (cont / res_q / 100.) != 0)
        {
            // cout << "cont = " << cont << endl;
            // cout << "err = " << err << endl;
            // cout << "res_q = " << res_q << endl;
            // cout << "cont / res_q / 100. = " << cont / res_q / 100. << endl;
            // cout << "err / res_q / 100. = " << err / res_q / 100. << endl;
        }
    }

    TProfile2D *Parity_Q_obs1 = (TProfile2D *)file->Get("pParity_e_Q2parent_obs1");
    TProfile2D *Parity_Q_obs2 = (TProfile2D *)file->Get("pParity_e_Q2parent_obs2");
    TProfile2D *Parity_w_Q_obs1 = (TProfile2D *)file->Get("pParity_w_Q2parent_obs1");
    TProfile2D *Parity_w_Q_obs2 = (TProfile2D *)file->Get("pParity_w_Q2parent_obs2");
    Parity_Q_obs1->Scale(1);
    Parity_Q_obs1->Add(Parity_w_Q_obs1);
    Parity_Q_obs2->Scale(1);
    Parity_Q_obs2->Add(Parity_w_Q_obs2);

    int index_for_data[12] = {3, 4, 7, 8, 15, 16, 19, 20, 27, 28, 31, 32};
    TH1 *Parity_Q_ss1 = (TH1 *)Parity_Q_obs1->ProjectionY("Parity_Q_ss1", index_for_data[0 + 2 * option112 + 4 * ep_option], index_for_data[0 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_os1 = (TH1 *)Parity_Q_obs1->ProjectionY("Parity_Q_os1", index_for_data[1 + 2 * option112 + 4 * ep_option], index_for_data[1 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_ss2 = (TH1 *)Parity_Q_obs2->ProjectionY("Parity_Q_ss2", index_for_data[0 + 2 * option112 + 4 * ep_option], index_for_data[0 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_os2 = (TH1 *)Parity_Q_obs2->ProjectionY("Parity_Q_os2", index_for_data[1 + 2 * option112 + 4 * ep_option], index_for_data[1 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_ss = (TH1 *)Parity_Q_ss1->Clone();
    TH1 *Parity_Q_os = (TH1 *)Parity_Q_os1->Clone();

    for (int i = (Nbin_min - 1); i < Nbin; i++)
    {
        float res_q = 1;

        // if (Parity_Q_ss2->GetBinCenter(i + 1) < Q2_range[cen])
        //     res_q = f_res->Eval(Parity_Q_ss2->GetBinCenter(i + 1));
        // else
        //     res_q = f_res->Eval(Q2_range[cen]);

        // cout << "Parity_Q_ss2->GetBinCenter(i + 1) = " << Parity_Q_ss2->GetBinCenter(i + 1) << endl;
        // cout << "p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i + 1)) = " << p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i + 1)) << endl;

        // cout << "Parity_Q_ss2->GetBinCenter(i + 1) = " << Parity_Q_ss2->GetBinCenter(i + 1) << endl;
        // cout << "Q2_range[ep_option][cen] = " << Q2_range[ep_option][cen] << endl;

        if (Parity_Q_ss2->GetBinCenter(i + 1) < Q2_range[ep_option][cen])
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i + 1)));
        else
        {
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]));
            // if(ep_option == 1) cout << "p_res->GetBinContent(41) = " << p_res->GetBinContent(41) << endl;
        }

        if (ep_option <= 1)
            res_q = sqrt(res_q);

        // if(ep_option == 1) cout << "Q2_range[ep_option][cen] = " << Q2_range[ep_option][cen] << ", bin center = " << Parity_Q_ss2->GetBinCenter(i + 1) << endl;

        float cont = Parity_Q_ss2->GetBinContent(i + 1);
        float err = Parity_Q_ss1->GetBinError(i + 1);
        Parity_Q_ss->SetBinContent(i + 1, cont / res_q / 100.);
        Parity_Q_ss->SetBinError(i + 1, err / res_q / 100.);

        cont = Parity_Q_os2->GetBinContent(i + 1);
        err = Parity_Q_os1->GetBinError(i + 1);
        Parity_Q_os->SetBinContent(i + 1, cont / res_q / 100.);
        Parity_Q_os->SetBinError(i + 1, err / res_q / 100.);

        if (cen == 8 && ep_option == 0 && (cont / res_q / 100.) != 0)
        {
            // cout << "cont = " << cont << endl;
            // cout << "err = " << err << endl;
            // cout << "res_q = " << res_q << endl;
            // cout << "cont / res_q / 100. = " << cont / res_q / 100. << endl;
            // cout << "err / res_q / 100. = " << err / res_q / 100. << endl;
        }
    }

    // rebin the data, three choices...
    const float diff = Q2_range[ep_option][cen] - Q2_range_min[ep_option][cen];
    const int N_bins = diff * 10 / rebin[ep_option][cen];
    if (cen == 8 && ep_option == 0 && option112 == 0)
    {
        // cout << "Q2_range[" << ep_option << "][" << cen << "] = " << Q2_range[ep_option][cen] << endl;
        // cout << "Q2_range_min[" << ep_option << "][" << cen << "] = " << Q2_range_min[ep_option][cen] << endl;
        // cout << "diff = " << diff << endl;
        // cout << "rebin = " << rebin[ep_option][cen] << endl;
        // cout << " N_bins = " << N_bins << endl;
    }
    // cout << "p_v2_Q_obs_new->GetBinCenter(" << N_bins << ") = " << p_v2_Q_obs_new->GetXaxis()->GetBinCenter(N_bins) << endl;
    // cout << endl;

    float v2q[N_bins], v2q_err[N_bins], d_gq[N_bins], d_gq_err[N_bins];

    if (rebin[ep_option][cen] == 1)
    {
        for (int i = 0; i < N_bins; i++)
        {
            v2q[i] = p_v2_Q_obs_new->GetBinContent(i + 1);
            v2q_err[i] = p_v2_Q_obs_new->GetBinError(i + 1);

            d_gq[i] = Parity_Q_os->GetBinContent(i + 1) - Parity_Q_ss->GetBinContent(i + 1);
            d_gq_err[i] = sqrt(pow(Parity_Q_ss->GetBinError(i + 1), 2) + pow(Parity_Q_ss->GetBinError(i + 1), 2));

            // cout << "v2q[" << i << "] = " << v2q[i] << endl;
            // cout << "v2q_err[" << i << "] = " << v2q_err[i] << endl;
            // cout << "d_gq[" << i << "] = " << d_gq[i] << endl;
            // cout << "d_gq_err[" << i << "] = " << d_gq_err[i] << endl;
        }
    }

    else if (rebin[ep_option][cen] == 2)
    {
        for (int i = 0; i < N_bins; i++)
        {
            double temp1 = p_v2_Q_obs_new->GetBinContent(2 * i + 1);
            double temp1_err = p_v2_Q_obs_new->GetBinError(2 * i + 1);
            double temp2 = p_v2_Q_obs_new->GetBinContent(2 * i + 2);
            double temp2_err = p_v2_Q_obs_new->GetBinError(2 * i + 2);
            // cout << "temp1 = " << temp1 << endl;
            // cout << "temp1_err = " << temp1_err << endl;
            // cout << "temp2 = " << temp2 << endl;
            // cout << "temp2_err = " << temp2_err << endl;
            v2q[i] = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err);
            v2q_err[i] = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err));
            // cout << "v2q[" << i << "] = " << v2q[i] << endl;
            // cout << "v2q_err[" << i << "] = " << v2q_err[i] << endl;

            temp1 = Parity_Q_os->GetBinContent(2 * i + 1);
            temp1_err = Parity_Q_os->GetBinError(2 * i + 1);
            temp2 = Parity_Q_os->GetBinContent(2 * i + 2);
            temp2_err = Parity_Q_os->GetBinError(2 * i + 2);
            // cout << "temp1 = " << temp1 << endl;
            // cout << "temp1_err = " << temp1_err << endl;
            // cout << "temp2 = " << temp2 << endl;
            // cout << "temp2_err = " << temp2_err << endl;
            double os = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err);
            double os_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err));
            // cout << "os = " << os << endl;
            // cout << "os_err = " << os_err << endl;

            temp1 = Parity_Q_ss->GetBinContent(2 * i + 1);
            temp1_err = Parity_Q_ss->GetBinError(2 * i + 1);
            temp2 = Parity_Q_ss->GetBinContent(2 * i + 2);
            temp2_err = Parity_Q_ss->GetBinError(2 * i + 2);
            // cout << "temp1 = " << temp1 << endl;
            // cout << "temp1_err = " << temp1_err << endl;
            // cout << "temp2 = " << temp2 << endl;
            // cout << "temp2_err = " << temp2_err << endl;
            double ss = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err);
            double ss_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err));
            // cout << "ss = " << ss << endl;
            // cout << "ss_err = " << ss_err << endl;

            d_gq[i] = os - ss;
            d_gq_err[i] = sqrt(pow(os_err, 2) + pow(ss_err, 2));
        }
    }

    else if (rebin[ep_option][cen] == 5)
    {
        double temp1, temp2, temp3, temp4, temp5, temp1_err, temp2_err, temp3_err, temp4_err, temp5_err;
        for (int i = 0; i < N_bins; i++)
        {
            temp1 = p_v2_Q_obs_new->GetBinContent(5 * i + 1);
            temp1_err = p_v2_Q_obs_new->GetBinError(5 * i + 1);
            temp2 = p_v2_Q_obs_new->GetBinContent(5 * i + 2);
            temp2_err = p_v2_Q_obs_new->GetBinError(5 * i + 2);
            temp3 = p_v2_Q_obs_new->GetBinContent(5 * i + 3);
            temp3_err = p_v2_Q_obs_new->GetBinError(5 * i + 3);
            temp4 = p_v2_Q_obs_new->GetBinContent(5 * i + 4);
            temp4_err = p_v2_Q_obs_new->GetBinError(5 * i + 4);
            temp5 = p_v2_Q_obs_new->GetBinContent(5 * i + 5);
            temp5_err = p_v2_Q_obs_new->GetBinError(5 * i + 5);

            v2q[i] = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err + temp3 / temp3_err / temp3_err + temp4 / temp4_err / temp4_err + temp5 / temp5_err / temp5_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err);
            v2q_err[i] = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err));

            temp1 = Parity_Q_os->GetBinContent(5 * i + 1);
            temp1_err = Parity_Q_os->GetBinError(5 * i + 1);
            temp2 = Parity_Q_os->GetBinContent(5 * i + 2);
            temp2_err = Parity_Q_os->GetBinError(5 * i + 2);
            temp3 = Parity_Q_os->GetBinContent(5 * i + 3);
            temp3_err = Parity_Q_os->GetBinError(5 * i + 3);
            temp4 = Parity_Q_os->GetBinContent(5 * i + 4);
            temp4_err = Parity_Q_os->GetBinError(5 * i + 4);
            temp5 = Parity_Q_os->GetBinContent(5 * i + 5);
            temp5_err = Parity_Q_os->GetBinError(5 * i + 5);
            double os = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err + temp3 / temp3_err / temp3_err + temp4 / temp4_err / temp4_err + temp5 / temp5_err / temp5_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err);
            double os_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err));
            temp1 = Parity_Q_ss->GetBinContent(5 * i + 1);
            temp1_err = Parity_Q_ss->GetBinError(5 * i + 1);
            temp2 = Parity_Q_ss->GetBinContent(5 * i + 2);
            temp2_err = Parity_Q_ss->GetBinError(5 * i + 2);
            temp3 = Parity_Q_ss->GetBinContent(5 * i + 3);
            temp3_err = Parity_Q_ss->GetBinError(5 * i + 3);
            temp4 = Parity_Q_ss->GetBinContent(5 * i + 4);
            temp4_err = Parity_Q_ss->GetBinError(5 * i + 4);
            temp5 = Parity_Q_ss->GetBinContent(5 * i + 5);
            temp5_err = Parity_Q_ss->GetBinError(5 * i + 5);
            double ss = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err + temp3 / temp3_err / temp3_err + temp4 / temp4_err / temp4_err + temp5 / temp5_err / temp5_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err);
            double ss_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err));

            d_gq[i] = os - ss;
            d_gq_err[i] = sqrt(pow(os_err, 2) + pow(ss_err, 2));
        }
    }

    int left[9] = {70, 60, 50, 40, 30, 20, 10, 5, 0};
    int right[9] = {80, 70, 60, 50, 40, 30, 20, 10, 5};
    char centrality[200];
    sprintf(centrality, "%d - %d %% Au+Au, 27 GeV", left[cen], right[cen]);

    // Q2-res
    TCanvas *c1 = new TCanvas(TString::Format("Resolution_%s_cen%d_%d", ep_option_name[ep_option].Data(), cen, option112).Data());
    c1->SetTopMargin(0.06);
    c1->SetRightMargin(0.1);
    c1->SetLeftMargin(0.15);
    c1->SetBottomMargin(0.17);
    c1->Draw();
    gPad->SetGridx(0);
    gPad->SetGridy(0);
    gPad->SetTickx();
    gPad->SetTicky();

    p_res->SetLineColor(4);
    p_res->SetMarkerColor(4);
    p_res->GetXaxis()->SetRangeUser(0, 50);
    p_res->GetXaxis()->SetTitle("Q_{2}^{2}");
    p_res->GetXaxis()->CenterTitle();
    p_res->GetYaxis()->SetTitle("Resolution");
    p_res->Draw();
    f_res->SetLineColor(2);
    f_res->Draw("same");

    /// Q2-G112 ssos
    TCanvas *c4 = new TCanvas(TString::Format("Parity Q2 %s cen %d %d", ep_option_name[ep_option].Data(), cen, option112).Data());
    c4->SetTopMargin(0.06);
    c4->SetRightMargin(0.1);
    c4->SetLeftMargin(0.15);
    c4->SetBottomMargin(0.17);
    c4->Draw();
    Parity_Q_ss->SetLineColor(2);
    Parity_Q_ss->SetMarkerColor(2);
    Parity_Q_ss->GetXaxis()->SetRangeUser(0, 50);
    Parity_Q_ss->SetMaximum(YRangeMax_ew_Q2[cen]);
    Parity_Q_ss->SetMinimum(-YRangeMax_ew_Q2[cen]);
    Parity_Q_ss->GetXaxis()->SetTitle("Q_{2}^{2}");
    Parity_Q_ss->GetXaxis()->CenterTitle();
    if (option112 == 0)
        Parity_Q_ss->GetYaxis()->SetTitle("#gamma^{112}");
    if (option112 == 1)
        Parity_Q_ss->GetYaxis()->SetTitle("#gamma^{132}");
    Parity_Q_ss->Draw();
    Parity_Q_os->SetLineColor(4);
    Parity_Q_os->SetMarkerColor(4);
    Parity_Q_os->Draw("same");

    TF1 *fu = new TF1("fu", "[0]+[1]*x+[2]*x*x", 0, Q2_range[ep_option][cen]);
    Parity_Q_ss->Fit("fu", "0QR", "", 0, Q2_range[ep_option][cen]);
    // Parity_Q_ss->Fit("fu", "0", "", 0, 50);
    float sig_ss = fu->GetParameter(0);
    float sig_ss_err = fu->GetParError(0);
    fu->SetLineColor(2);
    fu->Draw("sames");
    Parity_Q_os->Fit("fu", "0QR", "", 0, Q2_range[ep_option][cen]);
    // Parity_Q_os->Fit("fu", "0", "", 0, 50);
    float sig_os = fu->GetParameter(0);
    float sig_os_err = fu->GetParError(0);
    fu->SetLineColor(4);
    fu->Draw("sames");
    cout << endl;
    cout << "#gamma_{ss} at zero q = " << sig_ss << " +/- " << sig_ss_err << endl;
    cout << "#gamma_{os} at zero q = " << sig_os << " +/- " << sig_os_err << endl;
    cout << "difference = " << sig_os - sig_ss << " +/- " << sqrt(sig_ss_err * sig_ss_err + sig_os_err * sig_os_err) << endl;

    TLatex *tex1 = new TLatex(1, -YRangeMax_ew_Q2[cen] * 0.8, "same sign");
    tex1->SetTextSize(0.06);
    tex1->SetTextColor(2);
    tex1->Draw();
    TLatex *tex2 = new TLatex(1, -YRangeMax_ew_Q2[cen] * 0.6, "opposite sign");
    tex2->SetTextSize(0.06);
    tex2->SetTextColor(4);
    tex2->Draw();

    tex1->Draw();
    tex2->Draw();
    // c4->Print(printname4);

    // Q2-DG112-v2
    TCanvas *can5 = new TCanvas(TString::Format("parent_v2(Q2)_DG(Q2)_%s_cen%d_%d", ep_option_name[ep_option].Data(), cen, option112).Data(), TString::Format("parent_v2(Q2)_DG(Q2)_%s_cen%d", ep_option_name[ep_option].Data(), cen).Data(), 40, 40, 740, 500);
    can5->SetTopMargin(0.06);
    can5->SetRightMargin(0.1);
    can5->SetLeftMargin(0.15);
    can5->SetBottomMargin(0.17);
    can5->Draw();

    gPad->SetGridx(0);
    gPad->SetGridy(0);
    gPad->SetTickx();
    gPad->SetTicky();
    TH1F *histGraph2 = new TH1F(TString::Format("parent_v2(Q2)_DG(Q2)_%s_%d", ep_option_name[ep_option].Data(), option112).Data(), "", 100, 0, Xrange_delGamma[cen]);
    histGraph2->SetMaximum(Yrange_delGamma[cen]);
    histGraph2->SetMinimum(-Yrange_delGamma[cen]);
    histGraph2->SetLineColor(kBlack);
    histGraph2->GetYaxis()->SetTitleOffset(0.67);
    histGraph2->GetYaxis()->SetTitleSize(0.065);
    histGraph2->GetXaxis()->SetTitleSize(0.08);
    histGraph2->GetXaxis()->SetTitleOffset(0.90);
    histGraph2->GetYaxis()->CenterTitle();
    histGraph2->GetXaxis()->SetTitle("v_{2, parent}");
    histGraph2->GetXaxis()->CenterTitle();
    if (option112 == 0)
        histGraph2->GetYaxis()->SetTitle("#Delta#gamma^{112}");
    if (option112 == 1)
        histGraph2->GetYaxis()->SetTitle("#Delta#gamma^{132}");
    histGraph2->GetXaxis()->SetNdivisions(6);
    histGraph2->GetYaxis()->SetNdivisions(605);
    // histGraph2->GetYaxis()->SetLabelSize(lsize * 1.0);
    // histGraph2->GetXaxis()->SetLabelSize(lsize * 1.0);
    histGraph2->Draw();

    TGraphErrors *g112_v2_2 = new TGraphErrors(N_bins, v2q, d_gq, v2q_err, d_gq_err);
    g112_v2_2->SetMarkerStyle(kOpenStar);
    g112_v2_2->SetMarkerSize(1.5);
    g112_v2_2->SetMarkerColor(2);
    g112_v2_2->SetLineColor(2);
    g112_v2_2->SetFillColor(2);
    g112_v2_2->SetLineStyle(1);
    g112_v2_2->SetLineWidth(2);
    g112_v2_2->Draw("pe1");

    // if(ep_option == 1){
    //     cout << "N_bins = " << N_bins << endl;
    //     for(int lmn = 0 ; lmn < N_bins; lmn++){
    //         cout << "v2q[" << lmn << "] = " << v2q[lmn] << endl;
    //         cout << "v2q_err[" << lmn << "] = " << v2q_err[lmn] << endl;
    //         cout << "d_gq[" << lmn << "] = " << d_gq[lmn] << endl;
    //         cout << "d_gq_err[" << lmn << "] = " << d_gq_err[lmn] << endl;
    //     }
    //     cout << "Xfitting_delGamma[ep_option][cen] = " << Xfitting_delGamma[ep_option][cen] << endl;
    // }

    TF1 *fit_g112_v2_2 = new TF1("fit_g112_v2_2", "[0]+[1]*x", Xfitting_delGamma_min[option112][ep_option][cen], Xfitting_delGamma[ep_option][cen]);
    fit_g112_v2_2->SetLineStyle(4);
    fit_g112_v2_2->SetLineColor(2);
    fit_g112_v2_2->SetLineWidth(4);
    // g112_v2_2->Fit("fit_g112_v2_2", "EQR", "", 0, Xfitting_delGamma[ep_option][cen]);
    g112_v2_2->Fit("fit_g112_v2_2", "QR");

    float fit_sig2 = fit_g112_v2_2->GetParameter(0);
    float fit_sig2_err = fit_g112_v2_2->GetParError(0);
    cout << fit_sig2 << endl;
    cout << "ESE signal(Q2) = " << fit_sig2 << " +/- " << fit_sig2_err << endl
         << endl;

    TLatex *tex_new = new TLatex(Xrange_delGamma[cen] * 0.1, Yrange_delGamma[cen] * 0.8, centrality);
    tex_new->SetTextSize(0.08);
    tex_new->SetTextColor(1);
    tex_new->Draw();

    // TLatex *tex_new1 = new TLatex(Xrange_delGamma[cen] * 0.1, Yrange_delGamma[cen] * 0.8, TString::Format("ESE signal(Q2) = %e +/- %e", fit_sig2, fit_sig2_err).Data());
    // tex_new->SetTextSize(0.04);
    // tex_new->SetTextColor(1);
    // tex_new1->Draw();

    TLegend *legend2 = new TLegend(0.27, 0.22, 0.33, 0.37);
    legend2->SetFillColor(0);
    legend2->SetTextSize(0.035);
    legend2->SetLineColor(0);
    legend2->SetBorderSize(0);
    legend2->SetLineStyle(3);
    if (option112 == 0)
        legend2->AddEntry(g112_v2_2, "#Delta#gamma_{112}", "p");
    if (option112 == 1)
        legend2->AddEntry(g112_v2_2, "#Delta#gamma_{132}", "p");
    legend2->AddEntry((TObject *)0, TString::Format("ESE signal(Q2) = %e +/- %e", fit_sig2, fit_sig2_err).Data(), "");
    legend2->Draw();

    gamma[option112][cen][ep_option] = fit_sig2;
    gamma_err[option112][cen][ep_option] = fit_sig2_err;

    output_File->cd();
    c1->Write(printname1);
    c4->Write(printname4);
    can5->Write(printname5);
    p_v2_Q_obs_new->GetXaxis()->SetTitle("Q_{2}^{2}");
    p_v2_Q_obs_new->GetYaxis()->SetTitle("v_{2}");
    p_v2_Q_obs_new->Write();

    delete p_v2_Q_obs1;
    p_v2_Q_obs1 = NULL;
    delete p_v2_Q_obs2;
    p_v2_Q_obs2 = NULL;
    delete p_v2w_Q_obs1;
    p_v2w_Q_obs1 = NULL;
    delete p_v2w_Q_obs2;
    p_v2w_Q_obs2 = NULL;
    delete p_v2_Q_obs_new;
    p_v2_Q_obs_new = NULL;

    delete Parity_Q_obs1;
    Parity_Q_obs1 = NULL;
    delete Parity_Q_obs2;
    Parity_Q_obs2 = NULL;
    delete Parity_w_Q_obs1;
    Parity_w_Q_obs1 = NULL;
    delete Parity_w_Q_obs2;
    Parity_w_Q_obs2 = NULL;

    delete Parity_Q_ss1;
    Parity_Q_ss1 = NULL;
    delete Parity_Q_os1;
    Parity_Q_os1 = NULL;
    delete Parity_Q_ss2;
    Parity_Q_ss2 = NULL;
    delete Parity_Q_os2;
    Parity_Q_os2 = NULL;
    delete Parity_Q_ss;
    Parity_Q_ss = NULL;
    delete Parity_Q_os;
    Parity_Q_os = NULL;

    delete p_res;
    p_res = NULL;
    delete f_res;
    f_res = NULL;
}

void read_npart()
{
    std::fstream myfile("./npart.txt", std::ios_base::in);

    float a;
    int count = 0;
    string line_text;
    while (myfile >> a)
    {
        if (count % 2 == 0)
        {
            npart[count / 2] = a;
        }
        else if (count % 2 == 1)
        {
            npart_err[count / 2] = a;
        }

        count++;
    }
}

void read_purity_Q2()
{
    for (int curr_c = 0; curr_c < 9; curr_c++)
    {
        std::fstream myfile(TString::Format("./Fit_and_subtract/Q2_purity/output_cen%d.txt", curr_c), std::ios_base::in);

        // cout << "Centrality = " << curr_c << endl;
        
        float a = 0;
        int count = 0, count2 = 0;
        string line_text;
        
        while (getline(myfile, line_text))
        {
            
            istringstream ss(line_text);
            while (ss >> a)
            {
                if (a < 0)
                    break;

                // cout << a << ", ";

                lam_purity_Q2[count][curr_c][count2] = a;
            }
            // cout << endl;
            count++;
            count2 = 0;
        }

        std::fstream myfile2(TString::Format("./Fit_and_subtract/Q2_purity/output_cen%d_anti.txt", curr_c), std::ios_base::in);

        // cout << "Centrality = " << curr_c << endl;
        
        a = 0;
        count = 0;
        count2 = 0;
        line_text = "";
        
        while (getline(myfile2, line_text))
        {
            
            istringstream ss2(line_text);
            while (ss2 >> a)
            {
                if (a < 0)
                    break;

                // cout << a << ", ";

                antilam_purity_Q2[count][curr_c][count2] = a;
            }
            // cout << endl;
            count++;
            count2 = 0;
        }
    }

    myfile.close();
    myfile2.close();
}