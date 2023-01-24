#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TTree.h>
#include <TLine.h>
#include <iostream>
#include <fstream>

// variables to decide what to draw
bool compare_norm_rot = true;
bool show_scales = false;
int min_cen = 0;
int max_cen = 0;
TFile *output;
ofstream my_file;

double histfunc(Double_t *v, Double_t *par);
double histfunc2(Double_t *v, Double_t *par);
void compute_purity_stats(int cen, double bmin_value, double bmax_value, TString option);

void my_func_nosplit(int cen, const char* option_lam = "", double bmin_value = 1.113, double bmax_value = 1.119)
{
    // TString type = "lambda";

    min_cen = cen;
    max_cen = cen;

    if(option_lam == "lam"){
        option_lam = "";
    }
    else if(option_lam == "antilam"){
        option_lam = "_anti";
    }

    for (int i = min_cen; i <= max_cen; i++)
    {
        output = new TFile(TString::Format("output_cen%d%s.root", i, option_lam).Data(), "recreate");
        my_file.open(TString::Format("output_cen%d%s.txt", i, option_lam).Data());
        compute_purity_stats(i, bmin_value, bmax_value, TString::Format("0%s", option_lam));
        compute_purity_stats(i, bmin_value, bmax_value, TString::Format("1%s", option_lam));
        compute_purity_stats(i, bmin_value, bmax_value, TString::Format("2%s", option_lam));
        output->Close();
        my_file.close();
    }
}

void compute_purity_stats(int cen, double bmin_value, double bmax_value, TString option)
{
    TFile lam_nor(TString::Format("../../Traditional/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", cen).Data());

    TH2F *V0Mass = (TH2F *)lam_nor.Get(TString::Format("V0Mass_%s", option.Data()));

    TAxis *axis = V0Mass->ProjectionY("py", 1, 1)->GetXaxis();
    int bmin = axis->FindBin(bmin_value);
    int bmax = axis->FindBin(bmax_value);
    cout << bmin << ", " << bmax << endl;

    double yield_fit[2][9];
    double background_fit[2][9];
    TCanvas *canvasss[2][9], *full_fit[2][9], *subtracted[2][9];

    double sl_beg = 1.108;
    double sl_end = 1.136;

    int max_bins_to_fit = 80;

    for (int b = 1; b <= max_bins_to_fit; b++)
    {
        TH1F *tmp_V0Mass = V0Mass->ProjectionY("py", b, b)->Clone();
        if (tmp_V0Mass->GetEntries() == 0)
            break;

        /// Fitting with 2 Gaussian Functions and a third order polynomial
        TF1 *gaus1 = new TF1("gaus1", "gaus", 1.105, 1.125);
        TF1 *gaus2 = new TF1("gaus2", "gaus", 1.105, 1.125);
        TF1 *third_poly = new TF1("third_poly", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.09, 1.14);
        TF1 *third_poly1 = new TF1("third_poly1", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.09, 1.105);
        TF1 *third_poly2 = new TF1("third_poly2", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.125, 1.14);
        // TF1 *myfunc = new TF1("myfunc", histfunc, 1.07, 1.185, 2);
        // myfunc->SetParameters(i, 0.25);

        cout << "X max = " << tmp_V0Mass->GetXaxis()->GetXmax() << ", X min = " << tmp_V0Mass->GetXaxis()->GetXmin() << endl;

        TF1 *straight_line = new TF1("straight_line", "[0]*x + [1]", sl_beg, sl_end);
        double par4[2];
        par4[0] = ((tmp_V0Mass->GetBinContent(axis->FindBin(sl_end)) - tmp_V0Mass->GetBinContent(axis->FindBin(sl_beg))) / (sl_end - sl_beg));
        par4[1] = tmp_V0Mass->GetBinContent(axis->FindBin(sl_end)) - ((tmp_V0Mass->GetBinContent(axis->FindBin(sl_end)) - tmp_V0Mass->GetBinContent(axis->FindBin(sl_beg))) / (sl_end - sl_beg)) * sl_end;
        straight_line->SetParameters(par4);
        TH1F *V0Mass_nor_temp = (TH1F *)tmp_V0Mass->Clone();
        V0Mass_nor_temp->Add(straight_line, -1);

        TH1F *GausDist = new TH1F("GausDist", "Gaussian Peak", 200, tmp_V0Mass->GetXaxis()->GetXmin(), tmp_V0Mass->GetXaxis()->GetXmax());
        for (int k = 1; k < 40; k++)
        {
            GausDist->Fill((sl_beg + k * 0.0007), V0Mass_nor_temp->GetBinContent((axis->FindBin(sl_beg + k * 0.0007))));
        }
        V0Mass_nor_temp->Clear();
        V0Mass_nor_temp = (TH1F *)tmp_V0Mass->Clone();
        V0Mass_nor_temp->Add(GausDist, -1);

        // cout << "Let's save now!" << endl;
        output->cd();

        TCanvas *subtracted_temp = new TCanvas(TString::Format("subtracted_%d_%d_%s", cen, b, option.Data()).Data(), TString::Format("subtracted_%d_%d", cen, b).Data(), 1500, 800);
        V0Mass_nor_temp->Fit(third_poly, "QR");
        V0Mass_nor_temp->Draw();
        straight_line->SetLineColor(kBlue);
        straight_line->Draw("LSAME");

        TLegend *lam_leg = new TLegend(0.6, 0.6, 0.8, 0.8);
        Double_t line_height = tmp_V0Mass->GetBinContent(tmp_V0Mass->GetMaximumBin()) + 100;
        tmp_V0Mass->SetTitle("Lambda Mass Distribution");
        tmp_V0Mass->GetXaxis()->SetTitle("p #pi^{-} Invariant Mass (GeV/c^{2})");
        TLine *maxline = new TLine(bmax_value, 0, bmax_value, line_height);
        TLine *minline = new TLine(bmin_value, 0, bmin_value, line_height);
        lam_leg->SetHeader("(a)");

        TCanvas *canvas_temp = new TCanvas(TString::Format("canvas_Lambda_%d_%d_%s", cen, b, option.Data()).Data(), TString::Format("canvas_Lambda_%d_%d", cen, b).Data(), 1500, 800);
        tmp_V0Mass->GetYaxis()->SetTitle("Count");
        tmp_V0Mass->Draw();

        maxline->SetLineColor(kBlue);
        minline->SetLineColor(kBlue);
        maxline->Draw();
        minline->Draw();
        third_poly->Draw("LSAME");

        cout << "Done with finishing up Histogram" << endl;

        TLegend *leg1 = new TLegend(0.6, 0.6, 0.8, 0.8);
        leg1->AddEntry(minline, "Mass Cuts for Signal", "l");
        leg1->AddEntry(third_poly, "Background Fit", "l");
        leg1->SetBorderSize(0);
        leg1->SetFillStyle(0);
        leg1->Draw();
        lam_leg->SetBorderSize(0);
        lam_leg->SetFillStyle(0);
        lam_leg->Draw();
        TLegend *lam_leg2 = new TLegend(0.6, 0.6, 0.7, 0.7);
        lam_leg2->SetHeader("27 GeV Au+Au");
        lam_leg2->SetBorderSize(0);
        lam_leg2->SetFillStyle(0);
        lam_leg2->Draw();
        TLegend *lam_leg3 = new TLegend(0.6, 0.6, 0.7, 0.7);
        lam_leg3->SetHeader("20-30%");
        lam_leg3->SetBorderSize(0);
        lam_leg3->SetFillStyle(0);
        lam_leg3->Draw();
        cout << "Lines done" << endl;
        canvas_temp->Write();

        TCanvas *full_fit_temp = new TCanvas(TString::Format("canvas_Lambda_%d_%d_fullfit_%s", cen, b, option.Data()).Data(), TString::Format("canvas_Lambda_%d_%d_fullfit", cen, b).Data(), 1500, 800);
        tmp_V0Mass->Draw();
        third_poly->Draw("LSAME");

        TLine *maxline_low = new TLine(1.105, 0, 1.105, line_height);
        TLine *minline_low = new TLine(1.09, 0, 1.09, line_height);
        TLine *maxline_high = new TLine(1.14, 0, 1.14, line_height);
        TLine *minline_high = new TLine(1.125, 0, 1.125, line_height);
        maxline_low->SetLineColor(kGreen);
        minline_low->SetLineColor(kGreen);
        maxline_low->Draw();
        minline_low->Draw();
        maxline_high->SetLineColor(kGreen);
        minline_high->SetLineColor(kGreen);
        maxline_high->Draw();
        minline_high->Draw();

        TLegend *leg2 = new TLegend(0.6, 0.6, 0.8, 0.8);
        leg2->AddEntry(maxline_low, "Mass Cuts for Background Bands", "l");
        leg2->AddEntry(third_poly, "Background Fit", "l");
        leg2->Draw();

        full_fit_temp->Write();

        float yield_fit_tmp = 0;
        float background_fit_tmp = 0;
        for (int ddd = 0; ddd < (bmax - bmin); ddd++)
        {
            yield_fit_tmp += tmp_V0Mass->GetBinContent(bmin + ddd);
            background_fit_tmp += third_poly->Eval(bmin_value + ddd * 0.0007);
        }
        my_file << (yield_fit_tmp - background_fit_tmp)/yield_fit_tmp << " ";

        maxline->Clear();
        minline->Clear();
        maxline_low->Clear();
        minline_low->Clear();
        maxline_high->Clear();
        minline_high->Clear();

        gaus1->Clear();
        gaus2->Clear();
        // myfunc->Clear();
        third_poly->Clear();
        third_poly1->Clear();
        third_poly2->Clear();
        V0Mass_nor_temp->Clear();
        GausDist->Clear();
        straight_line->Clear();
        leg1->Clear();
        leg2->Clear();
        lam_leg->Clear();
        lam_leg2->Clear();
        lam_leg3->Clear();
    }

    my_file << "\n";

    // yield_fit[i][j] = 0;
    // background_fit[i][j] = 0;

    //

    // ofstream myfile;
    // myfile.open ("./Without_ptsplit/purity.txt");

    // for(int i = 0; i < 2; i++)
    // {
    //     for(int j = 0; j < 9; j++)
    //     {
    //         cout << "yield_count = " << yield_fit[i][j] << ", backgroundcount = " << background_fit[i][j] << endl;
    //         cout << "significance = " << (yield_fit[i][j] - background_fit[i][j]) / yield_fit[i][j] << endl;
    //         myfile << (yield_fit[i][j] - background_fit[i][j]) / yield_fit[i][j];
    //         myfile << " ";
    //         myfile << yield_fit[i][j];
    //         myfile << " ";
    //         //cout << "sqrt = " << yield_count[i]/(sqrt(yield_count[i] + backgroundcount[i]) * sqrt(RefMult[i].GetEntries())) << endl;

    //         myfile << "\n";
    //     }
    //     myfile << "\n";
    //     cout << "RefMult = " << RefMult[i].GetEntries() << endl;
    // }

    // myfile.close();
}

double histfunc(Double_t *v, Double_t *par)
{
    double i = par[0];
    double m = v[0];
    if (i == 0)
    {
        TFile lam_rot("plamrot.root");
        TH1F *V0Mass_rot = (TH1F *)lam_rot.Get("V0Mass")->Clone();
        TAxis *axis = V0Mass_rot->GetXaxis();
        int value_bin = axis->FindBin(m);
        double value = V0Mass_rot->GetBinContent(value_bin);
        V0Mass_rot->Clear();
        axis->Clear();
        return par[1] * value;
    }
    else if (i == 1)
    {
        TFile antilam_rot("pantilamrot.root");
        TH1F *V0Mass_rot = (TH1F *)antilam_rot.Get("V0Mass")->Clone();
        TAxis *axis = V0Mass_rot->GetXaxis();
        int value_bin = axis->FindBin(m);
        double value = V0Mass_rot->GetBinContent(value_bin);
        V0Mass_rot->Clear();
        axis->Clear();
        return par[1] * value;
    }
}

double histfunc2(Double_t *v, Double_t *par)
{
    double i = par[0];
    double m = v[0];

    double exponential_value = par[2] * exp(-0.5 * pow(((v[0] - par[3]) / par[4]), 2));
    // double exponential_value2 = par[5]*exp(-0.5*pow(((v[0]-par[6])/par[7]), 2));

    if (i == 0)
    {
        TFile lam_rot("plamrot.root");
        TH1F *V0Mass_rot = (TH1F *)lam_rot.Get("V0Mass")->Clone();
        TAxis *axis = V0Mass_rot->GetXaxis();
        int value_bin = axis->FindBin(m);
        double value = V0Mass_rot->GetBinContent(value_bin);
        V0Mass_rot->Clear();
        axis->Clear();
        // return (par[1] * value + exponential_value + exponential_value2);
        return (par[1] * value + exponential_value);
    }
    else if (i == 1)
    {
        TFile antilam_rot("pantilamrot.root");
        TH1F *V0Mass_rot = (TH1F *)antilam_rot.Get("V0Mass")->Clone();
        TAxis *axis = V0Mass_rot->GetXaxis();
        int value_bin = axis->FindBin(m);
        double value = V0Mass_rot->GetBinContent(value_bin);
        V0Mass_rot->Clear();
        axis->Clear();
        // return (par[1] * value + exponential_value + exponential_value2);
        return (par[1] * value + exponential_value);
    }
}