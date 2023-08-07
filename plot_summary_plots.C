#include <fstream>
#include <iostream>
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TCanvas.h"

using namespace std;

float values[4][16][9], stats_errs[4][16][9], sys_errs[16][9];
float cen9_eff[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
float cen9_err_eff[9] = {0.};
float cen9_eff2[9] = {74.5, 64.5, 54.5, 44.5, 34.5, 24.5, 14.5, 7, 2};
float cen9_eff3[9] = {75.5, 65.5, 55.5, 45.5, 35.5, 25.5, 15.5, 8, 3};
float cen9_eff4[9] = {74, 64, 54, 44, 34, 24, 14, 6.5, 1.5};
float cen9_eff5[9] = {76, 66, 56, 46, 36, 26, 16, 8.5, 3.5};
double npart[9] = {0.}, npart_err[9] = {0.};
TFile *output_File;
std::fstream systematic_file, sys_0_file_recipe1, sys_0_file_recipe2, sys_0_file_recipe3, sys_0_file_recipe4;

TString plot_names[16] = {"#Delta#gamma^{112}_{EPD}", 
                            "#Delta#gamma^{132}_{EPD}", 
                            "#Delta#gamma^{112}_{EPD1}", 
                            "#Delta#gamma^{132}_{EPD1}", 
                            "#Delta#gamma^{112}_{EPD, ESS}", 
                            "#Delta#gamma^{132}_{EPD, ESS}", 
                            "#Delta#gamma^{112}_{EPD1, ESS}", 
                            "#Delta#gamma^{132}_{EPD1, ESS}", 
                            "#kappa^{112}_{EPD}", 
                            "#kappa^{132}_{EPD}", 
                            "#kappa^{112}_{EPD1}", 
                            "#kappa^{132}_{EPD1}", 
                            "#kappa^{112}_{EPD, ESS}", 
                            "#kappa^{132}_{EPD, ESS}", 
                            "#kappa^{112}_{EPD1, ESS}", 
                            "#kappa^{132}_{EPD1, ESS}"};

TString y_axis_options[2] = {"#Delta#gamma * N_{part}", "#kappa"};

int first_index[12]  = {0, 0, 1, 2, 0, 1, 4, 0, 1, 4, 5, 6};
int second_index[12] = {1, 2, 3, 3, 4, 5, 5, 6, 7, 6, 7, 7};

void read_npart();
void read_sys_0_results(std::fstream sys_0_file, int recipe);
void read_sys_errs();
void plotting_TGraph_Helper(bool first, bool last, 
                            const char *title, const char *x_title, const char *y_title, int npts, 
                            float *x, float *x_err, 
                            float *y1, float *y1_err, float *y1_2_err, const char *leg_1 = NULL, 
                            float *y2 = NULL, float *y2_err = NULL, float *y2_2_err = NULL, const char *leg_2 = NULL, 
                            float *y3 = NULL, float *y3_err = NULL, float *y3_2_err = NULL, const char *leg_3 = NULL,
                            float *y4 = NULL, float *y4_err = NULL, float *y4_2_err = NULL, const char *leg_4 = NULL,
                            float *y5 = NULL, float *y5_err = NULL, float *y5_2_err = NULL, const char *leg_5 = NULL);

void plot_summary_plots(){
    
    output_File = new TFile("./KFParticle_Results/summary_plots.root", "recreate");
    // systematic_file.open("./KFParticle_Results/sys_errs.csv", std::ios_base::in);
    sys_0_file_recipe1.open("./KFParticle_Results/4_recipes/19GeV/sys_0/gamma_results_recipe1.txt", std::ios_base::in);
    sys_0_file_recipe2.open("./KFParticle_Results/4_recipes/19GeV/sys_0/gamma_results_recipe2.txt", std::ios_base::in);
    sys_0_file_recipe3.open("./KFParticle_Results/4_recipes/19GeV/sys_0/gamma_results_recipe3.txt", std::ios_base::in);
    sys_0_file_recipe4.open("./KFParticle_Results/4_recipes/19GeV/sys_0/gamma_results_recipe4.txt", std::ios_base::in);

    read_npart();
    read_sys_0_results(sys_0_file_recipe1, 1);
    read_sys_0_results(sys_0_file_recipe2, 2);
    read_sys_0_results(sys_0_file_recipe3, 3);
    read_sys_0_results(sys_0_file_recipe4, 4);
    // read_sys_errs();

    for (int r = 0; r < 4; r++){
        for (int i = 0; i < 12; i++){
            int first = first_index[i];
            int second = second_index[i];
            bool first_bool = false;
            if ((i == 0) && (r == 0)) first_bool = true;
            plotting_TGraph_Helper(first_bool, false, TString::Format("%s vs. %s Recipe %d", plot_names[first].Data(), plot_names[second].Data(), r), "Centrality %", y_axis_options[0], 9, cen9_eff, cen9_err_eff, values[r][first], stats_errs[r][first], sys_errs[first], plot_names[first].Data(), values[r][second], stats_errs[r][second], sys_errs[second], plot_names[second].Data());
        }

        for (int i = 0; i < 12; i++){
            int first = first_index[i] + 8;
            int second = second_index[i] + 8;
            bool second_bool = false;
            // if ((i == 11) && (r == 3)) second_bool = true;
            plotting_TGraph_Helper(false, second_bool, TString::Format("%s vs. %s Recipe %d", plot_names[first].Data(), plot_names[second].Data(), r), "Centrality %", y_axis_options[1], 9, cen9_eff, cen9_err_eff, values[r][first], stats_errs[r][first], sys_errs[first], plot_names[first].Data(), values[r][second], stats_errs[r][second], sys_errs[second], plot_names[second].Data());
        }
    }

    for(int i = 4; i < 8; i++){
        plotting_TGraph_Helper(false, false, plot_names[i], "Centrality %", y_axis_options[0], 9, cen9_eff, cen9_err_eff, 
                                values[0][i], stats_errs[0][i], sys_errs[i], "ESE: v_{2,pair} q^{2}_{2,pair}", 
                                values[1][i], stats_errs[1][i], sys_errs[i], "ESE: v_{2} q^{2}_{2,pair}",
                                values[2][i], stats_errs[2][i], sys_errs[i], "ESE: v_{2,pair} q^{2}_{2}",
                                values[3][i], stats_errs[3][i], sys_errs[i], "ESE: v_{2} q^{2}_{2}",
                                values[0][i-4], stats_errs[0][i-4], sys_errs[i-4], "Ensemble Average");
    }

    for(int i = 12; i < 16; i++){
        bool second_bool = false;
        if(i == 15) second_bool = true;

        plotting_TGraph_Helper(false, second_bool, plot_names[i], "Centrality %", y_axis_options[1], 9, cen9_eff, cen9_err_eff, 
                                values[0][i], stats_errs[0][i], sys_errs[i], "ESE: v_{2,pair} q^{2}_{2,pair}", 
                                values[1][i], stats_errs[1][i], sys_errs[i], "ESE: v_{2} q^{2}_{2,pair}",
                                values[2][i], stats_errs[2][i], sys_errs[i], "ESE: v_{2,pair} q^{2}_{2}",
                                values[3][i], stats_errs[3][i], sys_errs[i], "ESE: v_{2} q^{2}_{2}",
                                values[0][i-4], stats_errs[0][i-4], sys_errs[i-4], "Ensemble Average");
    }

    
    
    sys_0_file_recipe1.close();
    sys_0_file_recipe2.close();
    sys_0_file_recipe3.close();
    sys_0_file_recipe4.close();
    // systematic_file.close();   
    output_File->Close();
}

void read_sys_0_results(std::fstream& sys_0_file, int recipe)
{
    int row_c = 0;
    string line, word, temp;
  
    while (sys_0_file >> temp) {
  
        // used for breaking words
        stringstream s(temp);
  
        row_c++;
        if (row_c <= 33) continue;
        // read every column data of a row 
        while (getline(s, word, ', ')) {

            // 34 - 66 is one row --> corresponds to 33 values
            // 67 (33) - 99 second row
            // use row_c to keep track, if (row_c - 34) % 33 == 0 --> Centrality, skip
            // cent is (row_c - 34) / 33
            // graph is [(row_c - 34) / 33] - 1
            // But need to split into values and errors
            // So we check if {[(row_c - 34) / 33] - 1} % 2 == 0 then it is a value
            // if {[(row_c - 34) / 33] - 1} % 2 == 1 then it is an error

            float value = ::atof(word.c_str());
            
            if ((row_c - 34) % 33 == 0) continue;

            if( ((((row_c - 34) % 33) - 1) / 2) < 8){
                if ( (((row_c - 34) % 33) - 1) % 2 == 0 ) values[recipe-1][ ((((row_c - 34) % 33) - 1) / 2) ][ ((row_c - 34) / 33) ] = value * npart[((row_c - 34) / 33)];
                else stats_errs[recipe-1][ ((((row_c - 34) % 33) - 1) / 2) ][ ((row_c - 34) / 33) ] = value * npart[((row_c - 34) / 33)];
            }
            else{
                if ( (((row_c - 34) % 33) - 1) % 2 == 0 ) values[recipe-1][ ((((row_c - 34) % 33) - 1) / 2) ][ ((row_c - 34) / 33) ] = value;
                else stats_errs[recipe-1][ ((((row_c - 34) % 33) - 1) / 2) ][ ((row_c - 34) / 33) ] = value;
            }
            

            // cout << "Centrality = " << ((row_c - 34) / 33) << endl;
            // cout << "npart = " << npart[((row_c - 34) / 33)] << endl;

        }
  
    }
}

void read_sys_errs()
{
    int graph_c = 0, row_c = 0;
    
    string line, word, temp;
  
    while (systematic_file >> temp) {
  
        // used for breaking words
        stringstream s(temp);
  
        row_c++;
        // cout << "row_c = " << row_c << endl;
        if (row_c <= 17) continue;
        // read every column data of a row 
        while (getline(s, word, ', ')) {
  
            float value = ::atof(word.c_str());

            // row_c - 18 = Current Centrality
            // graph_c - 2 = Current Graph
            graph_c++;
            if (graph_c == 1) continue;

            if ((graph_c - 2) < 8) sys_errs[graph_c - 2][row_c - 18] = value * npart[row_c - 18];
            else sys_errs[graph_c - 2][row_c - 18] = value;

        }

        graph_c = 0;
  
    }
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

void plotting_TGraph_Helper(bool first, bool last, 
                            const char *title, const char *x_title, const char *y_title, int npts, 
                            float *x, float *x_err, 
                            float *y1, float *y1_err, float *y1_2_err, const char *leg_1, 
                            float *y2, float *y2_err, float *y2_2_err, const char *leg_2, 
                            float *y3, float *y3_err, float *y3_2_err, const char *leg_3,
                            float *y4, float *y4_err, float *y4_2_err, const char *leg_4,
                            float *y5, float *y5_err, float *y5_2_err, const char *leg_5)
{
    TGraphErrors *tmp_graph1 = new TGraphErrors(npts, x, y1, x_err, y1_err);
    if (y1_2_err != NULL){}
        TGraphErrors *tmp_graph1_2 = new TGraphErrors(npts, x, y1, x_err, y1_2_err);
    // cout << "y1_err[0] = " << y1_err[0] << endl;
    // cout << "y2_err[0] = " << y2_err[0] << endl;
    if (y2 != NULL)
        TGraphErrors *tmp_graph2 = new TGraphErrors(npts, cen9_eff2, y2, x_err, y2_err);
    if (y2_2_err != NULL)
        TGraphErrors *tmp_graph2_2 = new TGraphErrors(npts, cen9_eff2, y2, x_err, y2_2_err);
    if (y3 != NULL)
        TGraphErrors *tmp_graph3 = new TGraphErrors(npts, cen9_eff3, y3, x_err, y3_err);
    if (y3_2_err != NULL)
        TGraphErrors *tmp_graph3_2 = new TGraphErrors(npts, cen9_eff3, y3, x_err, y3_2_err);
    if (y4 != NULL)
        TGraphErrors *tmp_graph4 = new TGraphErrors(npts, cen9_eff4, y4, x_err, y4_err);
    if (y4_2_err != NULL)
        TGraphErrors *tmp_graph4_2 = new TGraphErrors(npts, cen9_eff4, y4, x_err, y4_2_err);
    if (y5 != NULL)
        TGraphErrors *tmp_graph5 = new TGraphErrors(npts, cen9_eff5, y5, x_err, y5_err);
    if (y5_2_err != NULL)
        TGraphErrors *tmp_graph5_2 = new TGraphErrors(npts, cen9_eff5, y5, x_err, y5_2_err);

    // for(int k = 0; k < npts; k++){
    //     cout << "x[" << k << "] = " << x[k] << endl;
    //     cout << "y1[" << k << "] = " << y1[k] << endl;
    //     cout << "y1_err[" << k << "] = " << y1_err[k] << endl;
    //     cout << "y1_2_err[" << k << "] = " << y1_2_err[k] << endl;

    //     cout << "y2[" << k << "] = " << y2[k] << endl;
    //     cout << "y2_err[" << k << "] = " << y2_err[k] << endl;
    //     cout << "y2_2_err[" << k << "] = " << y2_2_err[k] << endl;
    // }

    output_File->cd();
    TCanvas c1(title, title, 1500, 800);
    tmp_graph1_2->SetMarkerStyle(kFullCircle);
    tmp_graph1_2->SetMarkerColor(kBlack);
    tmp_graph1_2->SetLineWidth(4);
    tmp_graph1_2->SetLineColor(kBlack);
    tmp_graph1_2->SetTitle(title);
    tmp_graph1_2->GetXaxis()->SetTitle(x_title);
    tmp_graph1_2->GetYaxis()->SetTitle(y_title);
    tmp_graph1_2->Draw("[] AP");
    tmp_graph1->SetMarkerStyle(kFullCircle);
    tmp_graph1->SetMarkerColor(kBlack);
    tmp_graph1->SetLineColor(kBlack);
    tmp_graph1->Draw("SAMES P");
    if (y2 != NULL)
    {
        // cout << "plotting second" << endl;
        tmp_graph2_2->SetMarkerStyle(kFullCircle);
        tmp_graph2_2->SetMarkerColor(kRed);
        tmp_graph2_2->SetLineWidth(4);
        tmp_graph2_2->SetLineColor(kRed);
        tmp_graph2_2->Draw("SAMES [] P");
        tmp_graph2->SetMarkerStyle(kFullCircle);
        tmp_graph2->SetMarkerColor(kRed);
        tmp_graph2->SetLineColor(kRed);
        tmp_graph2->Draw("SAMES P");
    }
    if (y3 != NULL)
    {
        tmp_graph3_2->SetMarkerStyle(kOpenSquare);
        tmp_graph3_2->SetMarkerColor(kBlack);
        tmp_graph3_2->SetLineWidth(4);
        tmp_graph3_2->SetLineColor(kBlack);
        tmp_graph3_2->Draw("SAMES [] P");
        tmp_graph3->SetMarkerStyle(kOpenSquare);
        tmp_graph3->SetMarkerColor(kBlack);
        tmp_graph3->SetLineColor(kBlack);
        tmp_graph3->Draw("SAMES P");
    }
    if (y4 != NULL)
    {
        tmp_graph4_2->SetMarkerStyle(kOpenSquare);
        tmp_graph4_2->SetMarkerColor(kRed);
        tmp_graph4_2->SetLineWidth(4);
        tmp_graph4_2->SetLineColor(kRed);
        tmp_graph4_2->Draw("SAMES [] P");
        tmp_graph4->SetMarkerStyle(kOpenSquare);
        tmp_graph4->SetMarkerColor(kRed);
        tmp_graph4->SetLineColor(kRed);
        tmp_graph4->Draw("SAMES P");
    }
    if (y5 != NULL)
    {
        tmp_graph5_2->SetMarkerStyle(kFullTriangleUp);
        tmp_graph5_2->SetMarkerColor(kBlue);
        tmp_graph5_2->SetLineWidth(4);
        tmp_graph5_2->SetLineColor(kBlue);
        tmp_graph5_2->Draw("SAMES [] P");
        tmp_graph5->SetMarkerStyle(kFullTriangleUp);
        tmp_graph5->SetMarkerColor(kBlue);
        tmp_graph5->SetLineColor(kBlue);
        tmp_graph5->Draw("SAMES P");
    }
    TLine l1(0, 0, 70, 0);
    l1.Draw();

    TLegend *legend = new TLegend(0.1,0.7,0.25,0.9);
    legend->SetBorderSize(0);
    legend->SetFillColorAlpha(0,0);
    // legend->SetHeader("The Legend Title","C"); // option "C" allows to center the header
    if(leg_1 != NULL) {legend->AddEntry((TObject*)0, "", ""); legend->AddEntry(tmp_graph1,leg_1,"lep"); legend->AddEntry((TObject*)0, "", "");}
    if(leg_2 != NULL) {legend->AddEntry(tmp_graph2,leg_2,"lep"); legend->AddEntry((TObject*)0, "", "");}
    if(leg_3 != NULL) {legend->AddEntry(tmp_graph3,leg_3,"lep"); legend->AddEntry((TObject*)0, "", "");}
    if(leg_4 != NULL) {legend->AddEntry(tmp_graph4,leg_4,"lep"); legend->AddEntry((TObject*)0, "", "");}
    if(leg_5 != NULL) {legend->AddEntry(tmp_graph5,leg_5,"lep"); legend->AddEntry((TObject*)0, "", "");}
    legend->Draw();

    c1.Write();

    if(first) c1.Print("./KFParticle_Results/summary.pdf(","pdf");
    if(last) c1.Print("./KFParticle_Results/summary.pdf)","pdf");
    if(!first && !last) c1.Print("./KFParticle_Results/summary.pdf","pdf");
}