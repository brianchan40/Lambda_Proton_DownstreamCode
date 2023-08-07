#include <TGraph.h>
#include <TCanvas.h>
#include <TFile.h>
#include <iostream>
#include <fstream>

void plot_purity()
{
    TCanvas *c1 = new TCanvas("purity_by_cen", "purity_by_cen", 1000, 800);
    c1->Divide(2, 3);
    TCanvas *c2 = new TCanvas("purity_by_cen2", "purity_by_cen2", 1200, 640);
    c2->Divide(2, 2);

    for (int i = 0; i < 9; i++)
    {
        std::fstream myfile(TString::Format("output_cen%d_single_pion.txt", i), std::ios_base::in);

        float a = 0;
        int count = 0, count2 = 0, TPC_count = 0, EPD_count = 0, EPD1_count = 0;
        string line_text;
        // std::vector<float> TPC_tmp, EPD_tmp, EPD1_tmp, Q2_value;

        float TPC_tmp[1000] = {0.}, EPD_tmp[1000] = {0.}, EPD1_tmp[1000] = {0.}, Q2_value[1000] = {0.};

        while (getline(myfile, line_text))
        {
            // cout << line_text << endl;
            istringstream ss(line_text);
            while (ss >> a)
            {
                if (a < 0) break;
                // cout << "a = " << a << endl;

                // cout << a << endl;
                if (count == 0){
                    TPC_tmp[TPC_count] = a;
                    TPC_count++;
                }
                else if (count == 1){
                    EPD_tmp[EPD_count] = a;
                    EPD_count++;
                }
                else if (count == 2){
                    EPD1_tmp[EPD1_count] = a;
                    EPD1_count++;
                }

                Q2_value[count2] = 0.05 + 0.1 * count2;

                count2++;
            }
            // cout << endl;

            count++;
        }

        std::fstream myfile_anti(TString::Format("output_cen%d_single_pion_anti.txt", i), std::ios_base::in);

        a = 0;
        int count_anti = 0, count2_anti = 0, TPC_count_anti = 0, EPD_count_anti = 0, EPD1_count_anti = 0;
        string line_text_anti;
        float TPC_tmp_anti[1000] = {0.}, EPD_tmp_anti[1000] = {0.}, EPD1_tmp_anti[1000] = {0.}, Q2_value_anti[1000] = {0.};

        while (getline(myfile_anti, line_text_anti))
        {
            // cout << line_text << endl;
            istringstream ss(line_text_anti);
            while (ss >> a)
            {
                if (a < 0) break;
                // cout << "anti a = " << a << endl;

                // cout << a << endl;
                if (count_anti == 0){
                    TPC_tmp_anti[TPC_count_anti] = a;
                    TPC_count_anti++;
                }
                else if (count_anti == 1){
                    EPD_tmp_anti[EPD_count_anti] = a;
                    EPD_count_anti++;
                }
                else if (count_anti == 2){
                    EPD1_tmp_anti[EPD1_count_anti] = a;
                    EPD1_count_anti++;
                }

                Q2_value_anti[count2_anti] = 0.05 + 0.1 * count2_anti;

                count2_anti++;
            }
            // cout << endl;

            count_anti++;
        }

        TGraph *tpc_graph = new TGraph(TPC_count, Q2_value, TPC_tmp);
        TGraph *epd_graph = new TGraph(EPD_count, Q2_value, EPD_tmp);
        TGraph *epd_graph_anti = new TGraph(EPD_count_anti, Q2_value_anti, EPD_tmp_anti);
        TGraph *epd1_graph = new TGraph(EPD1_count, Q2_value, EPD1_tmp);

        if (i < 6) c1->cd(i + 1);
        else c2->cd(i-5);
        
        epd_graph->SetMarkerStyle(kCircle);
        epd_graph->SetMarkerSize(0.5);
        epd_graph->SetTitle(TString::Format("Centrality Bin %d", i+1).Data());
        epd_graph->GetXaxis()->SetTitle("q^{2}_{pion-pair}");
        // epd_graph->GetXaxis()->SetTitleSize(0.55);
        // epd_graph->GetXaxis()->SetTitleOffset(0.7);
        epd_graph->GetYaxis()->SetRangeUser(0.9, 1.02);
        epd_graph->GetYaxis()->SetTitle("purity % (# Signal / # Total)");
        // epd_graph->GetYaxis()->SetTitleSize(0.55);
        epd_graph->Draw("AP");
        epd_graph_anti->SetMarkerStyle(kCircle);
        epd_graph_anti->SetMarkerSize(0.5);
        epd_graph_anti->SetMarkerColor(kRed);
        epd_graph_anti->SetLineColor(kRed);
        epd_graph_anti->Draw("SAMES P");
        epd1_graph->SetMarkerStyle(kCircle);
        epd1_graph->SetMarkerColor(kBlue);
        epd1_graph->SetLineColor(kBlue);
        // epd1_graph->Draw("SAMES P");
    }

    TFile f("purity_plot.root", "recreate");
    c1->Write();
    c2->Write();
    f.Close();
}