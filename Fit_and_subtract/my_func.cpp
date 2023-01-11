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

TFile *output;
ofstream myfile, yieldfile;
double my_func_main(double bmin_value, double bmax_value);

double histfunc(Double_t *v, Double_t *par)
{
    double i = par[0];
    double m = v[0];
    if(i == 0)
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
    else if(i == 1)
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
    //double exponential_value2 = par[5]*exp(-0.5*pow(((v[0]-par[6])/par[7]), 2));

    if(i == 0)
    {
        TFile lam_rot("plamrot.root");
        TH1F *V0Mass_rot = (TH1F *)lam_rot.Get("V0Mass")->Clone();
        TAxis *axis = V0Mass_rot->GetXaxis();
        int value_bin = axis->FindBin(m);
        double value = V0Mass_rot->GetBinContent(value_bin);
        V0Mass_rot->Clear();
        axis->Clear();
        //return (par[1] * value + exponential_value + exponential_value2);
        return (par[1] * value + exponential_value);
    }
    else if(i == 1)
    {
        TFile antilam_rot("pantilamrot.root");
        TH1F *V0Mass_rot = (TH1F *)antilam_rot.Get("V0Mass")->Clone();
        TAxis *axis = V0Mass_rot->GetXaxis();
        int value_bin = axis->FindBin(m);
        double value = V0Mass_rot->GetBinContent(value_bin);
        V0Mass_rot->Clear();
        axis->Clear();
        //return (par[1] * value + exponential_value + exponential_value2);
        return (par[1] * value + exponential_value);
    }
}

void my_func(double bmin_value, double bmax_value){
    output = new TFile("output.root", "recreate");
    myfile.open ("purity.txt");
    yieldfile.open("yield.txt");

    for(int cen = 0; cen < 9; cen++){
        cout << "centrality = " << cen << endl;
        my_func_main(bmin_value, bmax_value, cen);
    }

    output->Close();
    myfile.close();
    yieldfile.close();

}

double my_func_main(double bmin_value, double bmax_value, int c)
{

    //variables to decide what to draw
    // bool compare_norm_rot = true;
    // bool show_scales = false;

    // TString type = "lambda";

    TFile lam_nor("plam.root");
    // TFile antilam_nor("pantilam.root");

    TH1F *V0Mass_nor[2][17];

    // TH1F *V0Mass_antilam_nor;
    // TH1F *V0Mass_antilam_rot;

    TH1F RefMult[2];

    RefMult[0] = (TH1F)lam_nor.Get("RefMult");
    // RefMult[1] = (TH1F)antilam_nor.Get("RefMult");
    RefMult[1] = (TH1F)lam_nor.Get("RefMult");
    // TH1I *cent_check_9_temp = (TH1I *)lam_nor.Get("cent_check_9");
    // Double_t normalize_const[7];
    // normalize_const[0] = cent_check_9_temp->GetBinContent(1) + cent_check_9_temp->GetBinContent(2);
    // normalize_const[1] = cent_check_9_temp->GetBinContent(3) + cent_check_9_temp->GetBinContent(4);
    // for(int i = 2; i < 7; i++)
    // {
    //     normalize_const[i] = cent_check_9_temp->GetBinContent(i + 3);
    // }

    // for(int j = 0; j < 9; j++)
    // {
        for(int l = 0; l < 17; l++)
        {
            V0Mass_nor[0][l] = (TH1F *)lam_nor.Get(TString::Format("V0Mass_%d_%d", c, l))->Clone();
            // V0Mass_nor[1][j][l] = (TH1F *)antilam_nor.Get(TString::Format("V0Mass_%d_%d", j, l))->Clone();
            V0Mass_nor[1][l] = (TH1F *)lam_nor.Get(TString::Format("V0Mass_anti_%d_%d", c, l))->Clone();
        }
    // }

    //double bmin_value = 0;
    //double bmax_value = 0;

    // if(type == "lambda")
    // {
    //     bmin_value = 1.112;
    //     bmax_value = 1.119;
    // }
    // else
    // {
    //     cout << "Something went wrong with defining the bounds!" << endl;
    // }

    TAxis *axis = V0Mass_nor[0][0]->GetXaxis();
    int bmin = axis->FindBin(bmin_value);
    int bmax = axis->FindBin(bmax_value);
    cout << bmin << ", " << bmax << endl;

    double yield_fit[2][17];
    double background_fit[2][17];
    TCanvas *canvasss[2][17], *full_fit[2][17], *subtracted[2][17];

    double sl_beg = 1.108;
    double sl_end = 1.136;

    /// Fitting with 2 Gaussian Functions and a third order polynomial
    for(int i = 0; i < 2; i++)
    {
        cout << "i = " << i << endl;
        // for(int j = 0; j < 9; j++)
        // {
            for(int m = 0 ; m < 17 ; m++)
            {
                cout << "m = " << m << endl;

                TF1 *gaus1 = new TF1("gaus1", "gaus", 1.105, 1.125);
                TF1 *gaus2 = new TF1("gaus2", "gaus", 1.105, 1.125);
                TF1 *third_poly = new TF1("third_poly", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.09, 1.14);
                TF1 *third_poly1 = new TF1("third_poly1", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.09, 1.105);
                TF1 *third_poly2 = new TF1("third_poly2", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.125, 1.14);
                TF1 *myfunc = new TF1("myfunc", histfunc, 1.07, 1.185, 2);
                myfunc->SetParameters(i, 0.25);

                //cout << "Setting of function?" << endl;

                //TF1 *final = new TF1("final", "gaus(0) + gaus(3) + [6]*pow(x, 3) + [7]*pow(x, 2) + [8]*x + [9]", (1.115684-0.07), (1.115684+0.07));
                // TF1 *final = new TF1("final", histfunc2, 1.07, 1.18, 5);
                // //Double_t par[8];
                // Double_t par[5];
                // gaus1->SetLineColor(1);
                // gaus2->SetLineColor(1);
                // myfunc->SetLineColor(1);
                // final->SetLineColor(2);

                // //cout << "Before Fitting" << endl;

                // V0Mass_nor[i]->Fit(gaus1, "QR");
                // //V0Mass_nor[i]->Fit(gaus2, "QR");
                // V0Mass_nor[i]->Fit(myfunc, "QR");
                // gaus1->GetParameters(&par[2]);
                // //gaus2->GetParameters(&par[5]);
                // myfunc->GetParameters(&par[0]);

                // final->SetParameters(par);

                // //cout << "Done with Round 1" << endl;

                // V0Mass_nor[i]->Fit(final, "QR");
                // final->GetParameters(&par[0]);
                // //myfunc->SetParameters(par[0], par[1]);
                // gaus1->GetParameters(&par[2]);
                //gaus2->GetParameters(&par[5]);

                // cout << "X max = " << V0Mass_nor[i][m]->GetXaxis()->GetXmax() << ", X min = " << V0Mass_nor[i][m]->GetXaxis()->GetXmin() << endl;

                TF1 *straight_line = new TF1("straight_line", "[0]*x + [1]", sl_beg, sl_end);
                double par4[2];
                par4[0] = ((V0Mass_nor[i][m]->GetBinContent(axis->FindBin(sl_end)) - V0Mass_nor[i][m]->GetBinContent(axis->FindBin(sl_beg))) / (sl_end - sl_beg));
                par4[1] = V0Mass_nor[i][m]->GetBinContent(axis->FindBin(sl_end)) - ((V0Mass_nor[i][m]->GetBinContent(axis->FindBin(sl_end)) - V0Mass_nor[i][m]->GetBinContent(axis->FindBin(sl_beg))) / (sl_end - sl_beg)) * sl_end;
                straight_line->SetParameters(par4);
                TH1F *V0Mass_nor_temp = (TH1F *)V0Mass_nor[i][m]->Clone();
                V0Mass_nor_temp->Add(straight_line, -1);

                TH1F *GausDist = new TH1F("GausDist", "Gaussian Peak", 200, V0Mass_nor[i][m]->GetXaxis()->GetXmin(), V0Mass_nor[i][m]->GetXaxis()->GetXmax()) ;
                for (int k = 1; k < 40; k++)
                {
                    GausDist->Fill((sl_beg + k * 0.0007), V0Mass_nor_temp->GetBinContent((axis->FindBin(sl_beg + k * 0.0007))));
                }
                V0Mass_nor_temp->Clear();
                V0Mass_nor_temp = (TH1F *)V0Mass_nor[i][m]->Clone();
                V0Mass_nor_temp->Add(GausDist, -1);

                //Double_t par_temp[10];
                //final->GetParameters(&par_temp[0]);
                //yield_fit[i] = par_temp[0];

                //yield_fit[i] = final->Integral(bmin_value, bmax_value);
                //background_fit[i] = myfunc->Integral(bmin_value, bmax_value);


                //cout << "Let's save now!" << endl;
                output->cd();

                char sub_name[200];
                sprintf(sub_name, "subtracted_%d_%d_%d", i, c, m);
                subtracted[i][m] = new TCanvas(sub_name, sub_name, 1500, 800);
                V0Mass_nor_temp->Fit(third_poly, "QR");
                V0Mass_nor_temp->Draw();
                straight_line->SetLineColor(kBlue);
                straight_line->Draw("LSAME");
                //subtracted[i]->Write();


                string c_name = "canvas_";
                char cc_name[200];
                TLegend *lam_leg = new TLegend(0.6, 0.6, 0.8, 0.8);
                Double_t line_height = V0Mass_nor[i][m]->GetBinContent(V0Mass_nor[i][m]->GetMaximumBin()) + 100;
                if(i == 0)
                {
                    V0Mass_nor[i][m]->SetTitle("Lambda Mass Distribution");
                    V0Mass_nor[i][m]->GetXaxis()->SetTitle("p #pi^{-} Invariant Mass (GeV/c^{2})");
                    c_name.append(TString::Format("Lambda_%d_%d", c, m));
                    TLine *maxline = new TLine(bmax_value, 0, bmax_value, line_height);
                    TLine *minline = new TLine(bmin_value, 0, bmin_value, line_height);

                    lam_leg->SetHeader("(a)");
                }
                else
                {
                    V0Mass_nor[i][m]->SetTitle("AntiLambda Mass Distribution");
                    V0Mass_nor[i][m]->GetXaxis()->SetTitle("#bar{p} #pi^{+} Invariant Mass (GeV/c^{2})");
                    c_name.append(TString::Format("AntiLambda_%d_%d", c, m));
                    TLine *maxline = new TLine(bmax_value, 0, bmax_value, line_height);
                    TLine *minline = new TLine(bmin_value, 0, bmin_value, line_height);

                    lam_leg->SetHeader("(b)");

                }
                strcpy(cc_name, c_name.c_str());
                canvasss[i][m] = new TCanvas(cc_name, cc_name, 1500, 800);
                //V0Mass_nor[i][j]->GetXaxis()->SetTitle("V0 Mass (GeV/c2)");
                V0Mass_nor[i][m]->GetYaxis()->SetTitle("Count");
                V0Mass_nor[i][m]->Draw();

                cout << "Done with finishing up Histogram" << endl;

                maxline->SetLineColor(kBlue);
                minline->SetLineColor(kBlue);
                maxline->Draw();
                minline->Draw();
                third_poly->Draw("LSAME");

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
                //straight_line->Draw("LSAME");
                cout << "Lines done" << endl;
                //myfunc->Draw("SAME");
                //cout << "Now to write it!" << endl;
                canvasss[i][m]->Write();

                // cout << "1" << endl;
                c_name.append("_fullfit");
                // cout << "2" << endl;
                strcpy(cc_name, c_name.c_str());
                // cout << "3" << endl;
                full_fit[i][m] = new TCanvas(cc_name, cc_name, 1500, 800);
                // cout << "4" << endl;
                V0Mass_nor[i][m]->Draw();
                // cout << "5" << endl;
                third_poly->Draw("LSAME");

                cout << "Done Drawing" << endl;

                if(i == 0)
                {
                    TLine *maxline_low = new TLine(1.105, 0, 1.105, 29000000);
                    TLine *minline_low = new TLine(1.09, 0, 1.09, 29000000);
                    TLine *maxline_high = new TLine(1.14, 0, 1.14, 29000000);
                    TLine *minline_high = new TLine(1.125, 0, 1.125, 29000000);
                }
                else
                {
                    TLine *maxline_low = new TLine(1.105, 0, 1.105, 8800000);
                    TLine *minline_low = new TLine(1.09, 0, 1.09, 8800000);
                    TLine *maxline_high = new TLine(1.14, 0, 1.14, 8800000);
                    TLine *minline_high = new TLine(1.125, 0, 1.125, 8800000);
                }
                maxline_low->SetLineColor(kGreen);
                minline_low->SetLineColor(kGreen);
                maxline_low->Draw();
                minline_low->Draw();
                maxline_high->SetLineColor(kGreen);
                minline_high->SetLineColor(kGreen);
                maxline_high->Draw();
                minline_high->Draw();

                cout << "done bkg lines" << endl;

                TLegend *leg2 = new TLegend(0.6, 0.6, 0.8, 0.8);
                leg2->AddEntry(maxline_low, "Mass Cuts for Background Bands", "l");
                leg2->AddEntry(third_poly, "Background Fit", "l");
                leg2->Draw();

                full_fit[i][m]->Write();
                //third_poly->Write();

                yield_fit[i][m] = 0;
                background_fit[i][m] = 0;

                for(int ddd = 0; ddd < (bmax - bmin) ; ddd++)
                {
                    yield_fit[i][m] += V0Mass_nor[i][m]->GetBinContent(bmin + ddd);
                    background_fit[i][m] += third_poly->Eval(bmin_value + ddd * 0.0007);
                }
                //background_fit[i] = third_poly->Integral(bmin_value, bmax_value);

                cout << "before clearing" << endl;

                maxline->Clear();
                minline->Clear();
                maxline_low->Clear();
                minline_low->Clear();
                maxline_high->Clear();
                minline_high->Clear();

                gaus1->Clear();
                gaus2->Clear();
                myfunc->Clear();
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

                cout << "done clearing" << endl;
            }
        // }
    }

    // output->Close();

    //cout << "here?" << endl;

    //double yield_count[2] = {0.};
    //int numbins = bmax - bmin;

    yieldfile << "Centrality " << c << "\n";

    for(int i = 0; i < 2; i++)
    {
        if(i == 0) yieldfile << "Lambda" << "\n";
        if(i == 1) yieldfile << "AntiLambda" << "\n";

        double total_yield = 0;

        // for(int j = 0; j < 9; j++)
        // {
            for(int m = 0; m < 17; m++)
            {
                cout << "yield_count = " << yield_fit[i][m] << ", backgroundcount = " << background_fit[i][m] << endl;
                cout << "significance = " << (yield_fit[i][m] - background_fit[i][m]) / yield_fit[i][m] << endl;
                myfile << (yield_fit[i][m] - background_fit[i][m]) / yield_fit[i][m];
                myfile << " ";
                myfile << yield_fit[i][m];
                myfile << " ";

                total_yield += (yield_fit[i][m] - background_fit[i][m]);
                //cout << "sqrt = " << yield_count[i]/(sqrt(yield_count[i] + backgroundcount[i]) * sqrt(RefMult[i].GetEntries())) << endl;
            }
            myfile << "\n";
        // }
        myfile << "\n";
        cout << "RefMult = " << RefMult[i].GetEntries() << endl;

        yieldfile << "total yield = " << total_yield << "\n";
    }
    //return (((yield_fit[0] - background_fit[0])/yield_fit[0]) + ((yield_fit[1] - background_fit[1])/yield_fit[1]))/2;
    // return (((yield_fit[0] - background_fit[0])) + ((yield_fit[1] - background_fit[1])))/2;
    return 2.0;
}
