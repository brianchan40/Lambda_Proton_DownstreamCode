#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TTree.h>
#include <TLine.h>

void find_scale();
double find_finalscale(TH1F *hist_norm, TH1F *hist_rot, double min_sc[], double max_sc[], int min_b, int max_b);
double efficiency_calc(int mom_ind, int cent_ind);
void fit_hist();

void fit_hist()
{

    //variables to decide what to draw
    bool compare_norm_rot = true;
    bool show_scales = false;

    TString type = "lambda";

    TFile lam_nor("plam.root");
    TFile lam_rot("plamrot.root");
    TFile antilam_nor("pantilam.root");
    TFile antilam_rot("pantilamrot.root");

    TH1F *V0Mass_lam_nor;
    TH1F *V0Mass_lam_rot;
    TH1F *V0Mass_antilam_nor;
    TH1F *V0Mass_antilam_rot;

    TFile *output = new TFile("output.root", "recreate");

    TH1F RefMult[2];

    RefMult[0] = (TH1F)lam_nor.Get("RefMult");
    RefMult[1] = (TH1F)antilam_nor.Get("RefMult");
    TH1I *cent_check_9_temp = (TH1I *)lam_nor.Get("cent_check_9");
    Double_t normalize_const[7];
    normalize_const[0] = cent_check_9_temp->GetBinContent(1) + cent_check_9_temp->GetBinContent(2);
    normalize_const[1] = cent_check_9_temp->GetBinContent(3) + cent_check_9_temp->GetBinContent(4);
    for(int i = 2; i < 7; i++)
    {
        normalize_const[i] = cent_check_9_temp->GetBinContent(i + 3);
    }

    V0Mass_lam_nor = (TH1F *)lam_nor.Get("V0Mass")->Clone();
    V0Mass_lam_rot = (TH1F *)lam_rot.Get("V0Mass")->Clone();
    V0Mass_antilam_nor = (TH1F *)antilam_nor.Get("V0Mass")->Clone();
    V0Mass_antilam_rot = (TH1F *)antilam_rot.Get("V0Mass")->Clone();

    double bmin_value = 0;
    double bmax_value = 0;

    if(type == "lambda")
    {
        bmin_value = 1.112;
        bmax_value = 1.119;
    }
    else
    {
        cout << "Something went wrong with defining the bounds!" << endl;
    }

    TAxis *axis = V0Mass_lam_nor->GetXaxis();
    int bmin = axis->FindBin(bmin_value);
    int bmax = axis->FindBin(bmax_value);
    cout << bmin << ", " << bmax << endl;
    int bmin_peak = axis->FindBin(bmin_value);
    int bmax_peak = axis->FindBin(bmax_value);

    ////////////////////////////////////////// STARTING TO FIND SCALE TO SCALE //////////////////////////////////////////
    TH1F *V0Mass_wo_back[2];
    int backgroundcount[2] = {0};
    double min_scale[2] = {0.}; // scale of the lower range
    double max_scale[2] = {0.}; // scale of the higher range

    /// minimum bound - to find min_scale
    /// first we look at the nonzero bins of the rotated graph and record them
    int curr_index_min_rot = 0;
    int x_bins_rot_min[100] = {0.};

    //cout << "Recording nonzero bins in rotated" << endl;

    for(int j = 1; j <= bmin; j++)
    {
        if(V0Mass_lam_rot->GetBinContent(j) != 0)
        {
            x_bins_rot_min[curr_index_min_rot] = j;
            //cout << "min = " << x_bins_rot_min[curr_index_min_rot] << endl;
            curr_index_min_rot++;
        }
    }

    double curr_scale_min = 0;
    int num_nonzero_min = 0;
    int curr_index_min_nor = 0;

    /// Finished finding the scale
    //min_scale = curr_scale_min/(double)num_nonzero_min;
    find_scale(V0Mass_lam_nor, V0Mass_lam_rot, curr_index_min_rot, x_bins_rot_min, min_scale);

    /// maximum bound - to find max_scale
    /// first we look at the nonzero bins of the rotated graph and record them
    int curr_index_max_rot = 0;
    int x_bins_rot_max[100] = {0.};

    //cout << "Recording nonzero bins in rotated" << endl;
    for(int j = bmax; j <= 200; j++)
    {
        if(V0Mass_lam_rot->GetBinContent(j) != 0)
        {
            //cout << "curr_index_max_rot = " << curr_index_max_rot << endl;
            x_bins_rot_max[curr_index_max_rot] = j;
            //cout << x_bins_rot_max[curr_index_max_rot] << endl;
            curr_index_max_rot++;
        }
    }

    double curr_scale_max = 0;
    int num_nonzero_max = 0;
    int curr_index_max_nor = 0;

    /// Finished finding the scale
    //max_scale = curr_scale_max/(double)num_nonzero_max;
    find_scale(V0Mass_lam_nor, V0Mass_lam_rot, curr_index_max_rot, x_bins_rot_max, max_scale);

    //double final_scale = -(min_scale+max_scale)/2.0;
    //double final_scale = -TMath::Max(max_scale[0], min_scale[0]);
    double final_scale = -find_finalscale(V0Mass_lam_nor, V0Mass_lam_rot, min_scale, max_scale, bmin, bmax);

    if(show_scales)
    {
        cout << "histogram Lambda" << endl;
        cout << "min_scale = " << min_scale[0] << endl;
        cout << "max_scale = " << max_scale[0] << endl;
        cout << "final_scale = " << final_scale << endl;
    }

    TH1F *temp_nor = (TH1F*)V0Mass_lam_nor->Clone();
    TH1F *temp_rot = (TH1F*)V0Mass_lam_rot->Clone();

    temp_nor->Add(temp_rot, final_scale);

    V0Mass_wo_back[0] = temp_nor;

    temp_rot->Scale(-final_scale);
    // for(int j = 0; j < 200; j++)
    // {
    //     backgroundcount[0] += (double)temp_rot->GetBinContent(j);
    // }

    backgroundcount[0] = temp_rot->Integral(bmin_peak, bmax_peak);

    TCanvas *lambda_normal = new TCanvas("lambda_normal", "Lambda", 1500, 800);
    V0Mass_lam_nor->GetXaxis()->SetTitle("Mass (GeV/c2)");
    V0Mass_lam_nor->GetYaxis()->SetTitle("Count");
    V0Mass_lam_nor->Draw();
    temp_rot->SetLineColor(kRed);
    temp_rot->Draw("SAMES");
    TLine *min_line = new TLine(bmin_value, 0, bmin_value, 29000000);
    TLine *max_line = new TLine(bmax_value, 0, bmax_value, 29000000);
    min_line->SetLineColor(kBlue);
    max_line->SetLineColor(kBlue);
    min_line->Draw();
    max_line->Draw();

    output->cd();
    lambda_normal->Write();

    temp_rot->Clear();
    temp_nor->Clear();
    min_line->Clear();
    max_line->Clear();

    ////////////////////////////////////////// ENDING TO FIND SCALE TO SCALE //////////////////////////////////////////


    ////////////////////////////////////////// STARTING TO FIND SCALE TO SCALE //////////////////////////////////////////
    memset(min_scale, 0, sizeof(min_scale));
    memset(min_scale, 0, sizeof(max_scale));

    /// minimum bound - to find min_scale
    /// first we look at the nonzero bins of the rotated graph and record them
    curr_index_min_rot = 0;
    memset(x_bins_rot_min, 0, sizeof(x_bins_rot_min));

    //cout << "Recording nonzero bins in rotated" << endl;

    for(int j = 1; j <= bmin; j++)
    {
        if(V0Mass_antilam_rot->GetBinContent(j) != 0)
        {
            x_bins_rot_min[curr_index_min_rot] = j;
            //cout << "min = " << x_bins_rot_min[curr_index_min_rot] << endl;
            curr_index_min_rot++;
        }
    }

    curr_scale_min = 0;
    num_nonzero_min = 0;
    curr_index_min_nor = 0;

    /// Finished finding the scale
    //min_scale = curr_scale_min/(double)num_nonzero_min;
    find_scale(V0Mass_antilam_nor, V0Mass_antilam_rot, curr_index_min_rot, x_bins_rot_min, min_scale);

    /// maximum bound - to find max_scale
    /// first we look at the nonzero bins of the rotated graph and record them
    curr_index_max_rot = 0;
    memset(x_bins_rot_max, 0, sizeof(x_bins_rot_max));

    //cout << "Recording nonzero bins in rotated" << endl;
    for(int j = bmax; j <= 200; j++)
    {
        if(V0Mass_antilam_rot->GetBinContent(j) != 0)
        {
            //cout << "curr_index_max_rot = " << curr_index_max_rot << endl;
            x_bins_rot_max[curr_index_max_rot] = j;
            //cout << x_bins_rot_max[curr_index_max_rot] << endl;
            curr_index_max_rot++;
        }
    }

    curr_scale_max = 0;
    num_nonzero_max = 0;
    curr_index_max_nor = 0;

    /// Finished finding the scale
    //max_scale = curr_scale_max/(double)num_nonzero_max;
    find_scale(V0Mass_antilam_nor, V0Mass_antilam_rot, curr_index_max_rot, x_bins_rot_max, max_scale);

    //double final_scale = -(min_scale+max_scale)/2.0;
    //double final_scale = -TMath::Max(max_scale[0], min_scale[0]);
    final_scale = -find_finalscale(V0Mass_antilam_nor, V0Mass_antilam_rot, min_scale, max_scale, bmin, bmax);

    if(show_scales)
    {
        cout << "histogram AntiLambda" << endl;
        cout << "min_scale = " << min_scale[0] << endl;
        cout << "max_scale = " << max_scale[0] << endl;
        cout << "final_scale = " << final_scale << endl;
    }

    temp_nor = (TH1F*)V0Mass_antilam_nor->Clone();
    temp_rot = (TH1F*)V0Mass_antilam_rot->Clone();

    temp_nor->Add(temp_rot, final_scale);

    V0Mass_wo_back[1] = temp_nor;

    temp_rot->Scale(-final_scale);
    // for(int j = 0; j < 200; j++)
    // {
    //     backgroundcount[1] += (double)temp_rot->GetBinContent(j);
    // }
    backgroundcount[1] = temp_rot->Integral(bmin_peak, bmax_peak);

    TCanvas *antilambda_normal = new TCanvas("antilambda_normal", "AntiLambda", 1500, 800);
    V0Mass_antilam_nor->GetXaxis()->SetTitle("Mass (GeV/c2)");
    V0Mass_antilam_nor->GetYaxis()->SetTitle("Count");
    V0Mass_antilam_nor->Draw();
    temp_rot->SetLineColor(kRed);
    temp_rot->Draw("SAMES");
    min_line = new TLine(bmin_value, 0, bmin_value, 8800000);
    max_line = new TLine(bmax_value, 0, bmax_value, 8800000);
    min_line->SetLineColor(kBlue);
    max_line->SetLineColor(kBlue);
    min_line->Draw();
    max_line->Draw();

    output->cd();
    antilambda_normal->Write();

    temp_rot->Clear();
    temp_nor->Clear();
    min_line->Clear();
    max_line->Clear();

    ////////////////////////////////////////// ENDING TO FIND SCALE TO SCALE //////////////////////////////////////////

    double yield_fit[2] = {0.};

    /// Fitting with 2 Gaussian Functions and a third order polynomial
    /*for(int i = 0; i < 2; i++)
    {
        TF1 *gaus1 = new TF1("gaus1", "gaus", 1.105, 1.125);
        TF1 *gaus2 = new TF1("gaus2", "gaus", 1.105, 1.125);
        TF1 *third_poly = new TF1("third_poly", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.07, 1.18);

        //TF1 *final = new TF1("final", "gaus(0) + gaus(3) + [6]*pow(x, 3) + [7]*pow(x, 2) + [8]*x + [9]", (1.115684-0.07), (1.115684+0.07));
        TF1 *final = new TF1("final", "gaus1 + gaus2 + third_poly", 1.07, 1.18);
        Double_t par[10];
        gaus1->SetLineColor(1);
        gaus2->SetLineColor(1);
        third_poly->SetLineColor(1);
        final->SetLineColor(2);

        V0Mass_wo_back[i]->Fit(gaus1, "QR");
        V0Mass_wo_back[i]->Fit(gaus2, "QR");
        V0Mass_wo_back[i]->Fit(third_poly, "QR");
        gaus1->GetParameters(&par[0]);
        gaus2->GetParameters(&par[3]);
        third_poly->GetParameters(&par[6]);

        final->SetParameters(par);

        V0Mass_wo_back[i]->Fit(final, "QR");

        //Double_t par_temp[10];
        //final->GetParameters(&par_temp[0]);
        //yield_fit[i] = par_temp[0];

        yield_fit[i] = final->Integral(1.1125, 1.1185);

        gaus1->Clear();
        gaus2->Clear();
        third_poly->Clear();
    }*/

    output->cd();
    for(int i = 0; i < 2; i++)
    {
        V0Mass_wo_back[i]->GetXaxis()->SetTitle("V0 Mass (GeV/c2)");
        V0Mass_wo_back[i]->GetYaxis()->SetTitle("Count");
        if(i == 0) V0Mass_wo_back[i]->SetTitle("Lambda");
        else V0Mass_wo_back[i]->SetTitle("AntiLambda");
        //if(cent >= 4){
        V0Mass_wo_back[i]->Write();
        //}
    }

    cout << "here?" << endl;

    double yield_count[2] = {0.};
    int numbins = bmax_peak - bmin_peak;


    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < numbins; j++)
        {
            yield_count[i] += (double)V0Mass_wo_back[i]->GetBinContent(j + bmin_peak);
        }

        cout << "yield_count = " << yield_count[i] << ", backgroundcount = " << backgroundcount[i] << endl;
        cout << "significance = " << yield_count[i]/(yield_count[i] + backgroundcount[i]) << endl;
        cout << "sqrt = " << yield_count[i]/(sqrt(yield_count[i] + backgroundcount[i]) * sqrt(RefMult[i].GetEntries())) << endl;
        cout << "RefMult = " << RefMult[i].GetEntries() << endl;
    }

    //cout << "number of events = " << num_events << endl;

    // Double_t spectrum_x[14] = {0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.15, 2.45, 2.8, 3.2, 3.6, 4.0};
    // Double_t spectrum_y[14][7];
    // Double_t error_x[14];
    // Double_t error_y[14][7];

    // for(int cent = 0; cent < 7; cent++)
    // {
    //     for(int i = 0; i < 8; i++)
    //     {
    //         spectrum_y[i][cent] = yield_count[i][cent] / (2 * TMath::Pi() * normalize_const[cent] * spectrum_x[i] * 0.1 * 1) * pow(10, (cent - 6));
    //         error_x[i] = 0.1;
    //         error_y[i][cent] = sqrt(yield_count[i][cent]) / yield_count[i][cent] * spectrum_y[i][cent];
    //     }
    //     for(int i = 8; i < 10; i++)
    //     {
    //         spectrum_y[i][cent] = yield_count[i][cent] / (2 * TMath::Pi() * normalize_const[cent] * spectrum_x[i] * 0.15 * 1) * pow(10, (cent - 6));
    //         error_x[i] = 0.15;
    //         error_y[i][cent] = sqrt(yield_count[i][cent]) / yield_count[i][cent] * spectrum_y[i][cent];
    //     }
    //     for(int i = 10; i < 14; i++)
    //     {
    //         spectrum_y[i][cent] = yield_count[i][cent] / (2 * TMath::Pi() * normalize_const[cent] * spectrum_x[i] * 0.2 * 1) * pow(10, (cent - 6));
    //         error_x[i] = 0.2;
    //         error_y[i][cent] = sqrt(yield_count[i][cent]) / yield_count[i][cent] * spectrum_y[i][cent];
    //     }

    //     //cout << "normalize_const " << cent << " = " << normalize_const[cent] << endl;
    // }

    // TGraphErrors *spectrum[7];
    // TMultiGraph *all_graphs = new TMultiGraph();
    // TLegend *legend = new TLegend(0.1, 0.7, 0.48, 0.9);

    // output->cd();

    // for(int cent = 0; cent < 7; cent++)
    // {
    //     Double_t temp_spectrumy[14];
    //     Double_t temp_errory[14];

    //     for(int m = 0; m < 14; m++)
    //     {
    //         temp_spectrumy[m] = spectrum_y[m][cent];
    //         temp_errory[m] = error_y[m][cent];
    //     }
    //     spectrum[cent] = new TGraphErrors(14, spectrum_x, temp_spectrumy, error_x, temp_errory);
    //     spectrum[cent]->GetXaxis()->SetTitle("Pt (GeV)");
    //     spectrum[cent]->GetYaxis()->SetTitle("Yield");
    //     spectrum[cent]->SetTitle(TString::Format("Raw Spectrum of Cent - %d", cent));
    //     spectrum[cent]->SetLineColor(cent + 1);
    //     spectrum[cent]->SetMarkerColor(cent + 1);
    //     spectrum[cent]->Write();

    //     all_graphs->Add(spectrum[cent]);
    //     legend->AddEntry(spectrum[cent], TString::Format("Centrality %d", cent));
    // }

    // output->cd();
    // all_graphs->SetTitle("Raw Spectrum of Lambda; Pt (GeV); Yield");
    // all_graphs->Write();

    // //Getting Efficiencies
    // Double_t efficencies[14][7];
    // Double_t spectrum_y_eff[14][7];
    // Double_t error_y_eff[14][7];

    // for(int i = 0; i < 14; i++)
    // {
    //     for(int j = 0; j < 7; j++)
    //     {
    //         double efficiency_temp = efficiency_calc(i, j);
    //         efficencies[i][j] = efficiency_temp;

    //         if(efficencies[i][j] != 0)
    //         {
    //             spectrum_y_eff[i][j] = spectrum_y[i][j] * efficencies[i][j];
    //             error_y_eff[i][j] = error_y[i][j] * efficencies[i][j];

    //             cout << "i = " << i << " j = " << j << " is not 0!" << endl;
    //         }
    //         else
    //         {
    //             spectrum_y_eff[i][j] = spectrum_y[i][j];
    //             error_y_eff[i][j] = error_y[i][j];


    //         }
    //     }
    // }

    // TGraphErrors *spectrum_eff[7];
    // TMultiGraph *all_graphs_eff = new TMultiGraph();
    // TLegend *legend_eff = new TLegend(0.1, 0.7, 0.48, 0.9);

    // output->cd();

    // for(int cent = 0; cent < 7; cent++)
    // {
    //     Double_t temp_spectrumy[14];
    //     Double_t temp_errory[14];

    //     for(int m = 0; m < 14; m++)
    //     {
    //         temp_spectrumy[m] = spectrum_y_eff[m][cent];
    //         temp_errory[m] = error_y_eff[m][cent];
    //     }
    //     spectrum_eff[cent] = new TGraphErrors(14, spectrum_x, temp_spectrumy, error_x, temp_errory);
    //     spectrum_eff[cent]->GetXaxis()->SetTitle("Pt (GeV)");
    //     spectrum_eff[cent]->GetYaxis()->SetTitle("Yield");
    //     spectrum_eff[cent]->SetTitle(TString::Format("Efficiency Corrected Spectrum of Cent - %d", cent));
    //     spectrum_eff[cent]->SetLineColor(cent + 1);
    //     spectrum_eff[cent]->SetMarkerColor(cent + 1);
    //     spectrum_eff[cent]->Write();

    //     all_graphs_eff->Add(spectrum_eff[cent]);
    //     legend_eff->AddEntry(spectrum_eff[cent], TString::Format("Centrality %d", cent));
    // }

    // output->cd();
    // all_graphs_eff->SetTitle("Efficiency Corrected Spectrum of Lambda; Pt (GeV); Yield");
    // all_graphs_eff->Write();
}

void find_scale(TH1F *hist_norm, TH1F *hist_rot, int num_nonzero_bins, int x_bins_rot[], double scale[])
{
    TH1F *temp_norm = hist_norm->Clone();
    TH1F *temp_rotate = hist_rot->Clone();

    double curr_scale_1 = 0;
    int num_nonzero_1 = 0;
    double curr_scale_2 = 0;
    int curr_index_nor = 0;

    for(int k = 0; k < num_nonzero_bins; k++)
    {
        int curr_index_checked = x_bins_rot[k];
        double curr_scale_checked = 0;

        if(temp_norm->GetBinContent(curr_index_checked) != 0)
        {
            curr_scale_checked = ((double)temp_norm->GetBinContent(curr_index_checked)) / ((double)temp_rotate->GetBinContent(curr_index_checked));
            //cout << "min = " << curr_index_checked << endl;
            curr_scale_1 += curr_scale_checked;
            num_nonzero_1++;
        }
        else
        {
            continue;
        }

        if((curr_scale_checked < curr_scale_2) || (curr_scale_2 == 0))
        {
            curr_scale_2 = curr_scale_checked;
            curr_index_nor = curr_index_checked;
        }
    }



    if((double)num_nonzero_1 == 0)
    {
        //cout << "ERRRORRRRR!!!!" << endl;
        scale[0] = 0;
    }
    else
    {
        scale[0] = curr_scale_1 / (double)num_nonzero_1;
    }

    scale[1] = curr_scale_2;
}

double find_finalscale(TH1F *hist_norm, TH1F *hist_rot, double min_sc[], double max_sc[], int min_b, int max_b)
{
    TH1F *temp_norm = hist_norm->Clone();
    TH1F *temp_rotate = hist_rot->Clone();
    double all_scales[4] = {min_sc[0], min_sc[1], max_sc[0], max_sc[1]};

    double bin_content_neg = 0;
    double neg_sc = 0;
    double bin_content_pos = 0;
    double pos_sc = 0;

    double our_scale = 0;

    for(int i = 0; i < 4; i++)
    {
        TH1F *double_temp_norm = temp_norm->Clone();
        TH1F *double_temp_rotate = temp_rotate->Clone();

        double_temp_norm->Add(double_temp_rotate, -all_scales[i]);

        cout << "all_scales = " << all_scales[i] << endl;


        double bin_content_curr = 0;
        for(int j = 1; j <= min_b; j++)
        {
            bin_content_curr += double_temp_norm->GetBinContent(j);
        }
        for(int k = max_b; k <= 200; k++)
        {
            bin_content_curr += double_temp_norm->GetBinContent(k);
        }

        if(bin_content_curr < 0)
        {
            if((bin_content_curr >= bin_content_neg) || bin_content_neg == 0)
            {
                bin_content_neg = bin_content_curr;
                neg_sc = all_scales[i];
                //cout << "Do I even get in here? " << endl;
            }
        }
        else if(bin_content_curr > 0)
        {
            if((bin_content_curr <= bin_content_pos) || bin_content_pos == 0)
            {
                bin_content_pos = bin_content_curr;
                pos_sc = all_scales[i];
            }
        }
        else
        {
            our_scale = all_scales[i];
        }
    }

    cout << "STARTING NOW!" << endl;
    //cout << "neg_sc = " << neg_sc << endl;
    //cout << "pos_sc = " << pos_sc << endl;

    if(neg_sc == 0)
    {
        neg_sc = pos_sc + pos_sc * 0.5;
    }

    if(our_scale != 0)
    {
        continue;
    }
    else
    {
        double step_size = fabs(neg_sc - pos_sc) / 100.0;
        bool keepgoing = true;

        while(keepgoing)
        {
            double mid = (neg_sc + pos_sc) / 2.0;

            cout << "neg_sc = " << neg_sc << endl;
            cout << "pos_sc = " << pos_sc << endl;
            cout << "mid = " << mid << endl;

            TH1F *double_temp_norm = (TH1F *)temp_norm->Clone();
            TH1F *double_temp_rotate = (TH1F *)temp_rotate->Clone();

            double_temp_norm->Add(double_temp_rotate, -mid);

            double bin_content_curr = 0;
            for(int j = 1; j <= min_b; j++)
            {
                bin_content_curr += double_temp_norm->GetBinContent(j);
            }
            for(int k = max_b; k <= 200; k++)
            {
                bin_content_curr += double_temp_norm->GetBinContent(k);
            }

            bool pos_sc_choice = true;

            cout << "bin_content_curr = " << bin_content_curr << endl;

            if(bin_content_curr < 0)
            {
                neg_sc = mid;
                pos_sc_choice = false;

                //cout << "Replace neg" << endl;
            }
            else if(bin_content_curr > 0)
            {
                pos_sc = mid;

                //cout << "Replace pos" << endl;
            }
            else
            {
                keepgoing = false;
                our_scale = mid;

                cout << "Done early!" << endl;
                continue;
            }

            cout << "fabs(neg_sc - pos_sc) = " << fabs(neg_sc - pos_sc) << endl;
            cout << "step_size = " << step_size << endl;

            if(fabs(neg_sc - pos_sc) <= step_size)
            {
                keepgoing = false;

                //cout << "Last step in finding scale" << endl;

                double check_now_sc = 0;
                TH1F *triple_temp_norm = temp_norm->Clone();
                TH1F *triple_temp_rotate = temp_rotate->Clone();

                if(pos_sc_choice == true)
                {

                    //cout << "Pos was switched" << endl;

                    triple_temp_norm->Add(triple_temp_rotate, -neg_sc);

                    //cout << "Was I okay adding?" << endl;

                    double bin_content_now = 0;
                    //cout << "min_b = " << min_b << endl;
                    for(int l = 1; l <= min_b; l++)
                    {
                        bin_content_now += triple_temp_norm->GetBinContent(l);
                    }
                    //cout << "max_b = " << max_b << endl;
                    for(int m = max_b; m <= 200; m++)
                    {
                        bin_content_now += triple_temp_norm->GetBinContent(m);
                        //cout << "m = " << m << endl;
                    }

                    //cout << "couldn't get past here?" << endl;

                    if(fabs(bin_content_curr) >= fabs(bin_content_now))
                    {
                        our_scale = neg_sc;
                    }
                    else
                    {
                        our_scale = pos_sc;
                    }

                    //cout << "How about here?" << endl;
                }
                else
                {
                    //cout << "Neg was switched" << endl;

                    triple_temp_norm->Add(triple_temp_rotate, -pos_sc);

                    double bin_content_now = 0;
                    for(int l = 1; l <= min_b; l++)
                    {
                        bin_content_now += triple_temp_norm->GetBinContent(l);
                    }
                    for(int m = max_b; m <= 200; m++)
                    {
                        bin_content_now += triple_temp_norm->GetBinContent(m);
                    }

                    if(fabs(bin_content_curr) <= fabs(bin_content_now))
                    {
                        our_scale = neg_sc;
                    }
                    else
                    {
                        our_scale = pos_sc;
                    }
                }
            }
        }
    }

    return our_scale;
}

double efficiency_calc(int mom_ind, int cent_ind)
{
    bool debug_efficiency_calc = false;

    TFile tree_file("auau200GeV_run11_la_fp.0.la.picodst.root");
    TTree *data_tree = (TTree *)tree_file.Get("McV0PicoDst");

    int centBin9 = 0;
    float mcv0pt[100], v0pt[100];
    int nmcv0 = 0;
    int nv0 = 0;

    data_tree->SetBranchAddress("centBin9", &centBin9);
    data_tree->SetBranchAddress("nmcv0", &nmcv0);
    data_tree->SetBranchAddress("nv0", &nv0);
    data_tree->SetBranchAddress("mcv0pt", mcv0pt);
    data_tree->SetBranchAddress("v0pt", v0pt);

    int nentries = data_tree->GetEntries();
    int mc_count = 0;
    int rec_count = 0;

    for(int n = 0; n < nentries; n++)
    {
        data_tree->GetEntry(n);


        int cent_interest[2] = {0.};
        int cent_interest_count = 0;

        if(debug_efficiency_calc == true) cout << "Step 1" << endl;

        if(cent_ind <= 1)
        {
            cent_interest_count = 2;
            cent_interest[0] = cent_ind * 2;

            //cout << "before cent_interest = " << cent_interest[0] << endl;

            cent_interest[1] = cent_ind * 2 + 1;
        }
        else
        {
            cent_interest_count = 1;
            cent_interest[0] = cent_ind + 2;
            //cout << "before cent_interest = " << cent_interest[0] << endl;
        }


        if(debug_efficiency_calc == true) cout << "Step 2" << endl;


        float mom_min = 0;
        float mom_max = 0;

        if(mom_ind == 0)
        {
            mom_min = 0.4;
            mom_max = 0.6;
        }
        else if(mom_ind == 1)
        {
            mom_min = 0.6;
            mom_max = 0.8;
        }
        else if(mom_ind == 2)
        {
            mom_min = 0.8;
            mom_max = 1.0;
        }
        else if(mom_ind == 3)
        {
            mom_min = 1.0;
            mom_max = 1.2;
        }
        else if(mom_ind == 4)
        {
            mom_min = 1.2;
            mom_max = 1.4;
        }
        else if(mom_ind == 5)
        {
            mom_min = 1.4;
            mom_max = 1.6;
        }
        else if(mom_ind == 6)
        {
            mom_min = 1.6;
            mom_max = 1.8;
        }
        else if(mom_ind == 7)
        {
            mom_min = 1.8;
            mom_max = 2.0;
        }
        else if(mom_ind == 8)
        {
            mom_min = 2.0;
            mom_max = 2.3;
        }
        else if(mom_ind == 9)
        {
            mom_min = 2.3;
            mom_max = 2.6;
        }
        else if(mom_ind == 10)
        {
            mom_min = 2.6;
            mom_max = 3.0;
        }
        else if(mom_ind == 11)
        {
            mom_min = 3.0;
            mom_max = 3.4;
        }
        else if(mom_ind == 12)
        {
            mom_min = 3.4;
            mom_max = 3.8;
        }
        else if(mom_ind == 13)
        {
            mom_min = 3.8;
            mom_max = 4.2;
        }

        //cout << "cent_interest_count = " << cent_interest_count << endl;

        for(int i = 0; i < cent_interest_count; i++)
        {

            //cout << "centBin9 = " << centBin9 << endl;
            //cout << "cent_interest = " << cent_interest[i] << endl;

            if(centBin9 == cent_interest[i])
            {

                cout << "nmcv0 = " << nmcv0 << endl;

                for(int j = 0; j < nmcv0; j++)
                {

                    cout << "centBin9 = " << centBin9 << ", mcv0pt = " << mcv0pt[j] << endl;

                    if((mcv0pt[j] <= mom_max) && (mcv0pt[j] >= mom_min))
                    {
                        mc_count++;
                    }
                }
                for(int k = 0; k < nv0; k++)
                {
                    if((v0pt[k] <= mom_max) && (v0pt[k] >= mom_min))
                    {
                        rec_count++;
                    }
                }
            }
        }

    }

    if(debug_efficiency_calc == true) cout << "Step 3" << endl;

    if(rec_count != 0)
    {
        double calc_eff = (double)mc_count / (double)rec_count;
    }
    else
    {
        double calc_eff = 0;
    }

    if(debug_efficiency_calc == true) cout << "Step 4" << endl;

    //return inverted efficiency
    return calc_eff;

}