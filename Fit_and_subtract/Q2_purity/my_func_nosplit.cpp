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

double my_func_nosplit(double bmin_value, double bmax_value)
{

    //variables to decide what to draw
    bool compare_norm_rot = true;
    bool show_scales = false;

    TString type = "lambda";

    TFile lam_nor("plam.root");
    // TFile antilam_nor("pantilam.root");

    TH1F *V0Mass_nor[2][9];

    // TH1F *V0Mass_antilam_nor;
    // TH1F *V0Mass_antilam_rot;

    TFile *output = new TFile("./Without_ptsplit/output.root", "recreate");

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

    for(int j = 0; j < 9; j++)
    {
        TH1F *V0Mass_nor_lam_temp, *V0Mass_nor_antilam_temp;

        for(int l = 0; l < 17; l++)
        {
            if(l == 0)
            {
                V0Mass_nor_lam_temp = (TH1F *)lam_nor.Get(TString::Format("V0Mass_%d_%d", j, l))->Clone();
                // V0Mass_nor_antilam_temp = (TH1F *)antilam_nor.Get(TString::Format("V0Mass_%d_%d", j, l))->Clone();
                V0Mass_nor_antilam_temp = (TH1F *)lam_nor.Get(TString::Format("V0Mass_anti_%d_%d", j, l))->Clone();
            }
            else
            {
                V0Mass_nor_lam_temp->Add((TH1F *)lam_nor.Get(TString::Format("V0Mass_%d_%d", j, l))->Clone());
                // V0Mass_nor_antilam_temp->Add((TH1F *)antilam_nor.Get(TString::Format("V0Mass_%d_%d", j, l))->Clone());
                V0Mass_nor_antilam_temp->Add((TH1F *)lam_nor.Get(TString::Format("V0Mass_anti_%d_%d", j, l))->Clone());
            }

        }

        V0Mass_nor[0][j] = (TH1F *)V0Mass_nor_lam_temp->Clone();
        V0Mass_nor[1][j] = (TH1F *)V0Mass_nor_antilam_temp->Clone();

        delete V0Mass_nor_lam_temp;
        V0Mass_nor_lam_temp = NULL;
        delete V0Mass_nor_antilam_temp;
        V0Mass_nor_antilam_temp = NULL;
    }

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

    double yield_fit[2][9];
    double background_fit[2][9];
    TCanvas *canvasss[2][9], *full_fit[2][9], *subtracted[2][9];

    double sl_beg = 1.108;
    double sl_end = 1.136;

    /// Fitting with 2 Gaussian Functions and a third order polynomial
    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 9; j++)
        {
            cout << "i = " << i << endl;

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

            cout << "X max = " << V0Mass_nor[i][j]->GetXaxis()->GetXmax() << ", X min = " << V0Mass_nor[i][j]->GetXaxis()->GetXmin() << endl;

            TF1 *straight_line = new TF1("straight_line", "[0]*x + [1]", sl_beg, sl_end);
            double par4[2];
            par4[0] = ((V0Mass_nor[i][j]->GetBinContent(axis->FindBin(sl_end)) - V0Mass_nor[i][j]->GetBinContent(axis->FindBin(sl_beg))) / (sl_end - sl_beg));
            par4[1] = V0Mass_nor[i][j]->GetBinContent(axis->FindBin(sl_end)) - ((V0Mass_nor[i][j]->GetBinContent(axis->FindBin(sl_end)) - V0Mass_nor[i][j]->GetBinContent(axis->FindBin(sl_beg))) / (sl_end - sl_beg)) * sl_end;
            straight_line->SetParameters(par4);
            TH1F *V0Mass_nor_temp = (TH1F *)V0Mass_nor[i][j]->Clone();
            V0Mass_nor_temp->Add(straight_line, -1);

            TH1F *GausDist = new TH1F("GausDist", "Gaussian Peak", 200, V0Mass_nor[i][j]->GetXaxis()->GetXmin(), V0Mass_nor[i][j]->GetXaxis()->GetXmax()) ;
            for (int k = 1; k < 40; k++)
            {
                GausDist->Fill((sl_beg + k * 0.0007), V0Mass_nor_temp->GetBinContent((axis->FindBin(sl_beg + k * 0.0007))));
            }
            V0Mass_nor_temp->Clear();
            V0Mass_nor_temp = (TH1F *)V0Mass_nor[i][j]->Clone();
            V0Mass_nor_temp->Add(GausDist, -1);

            //Double_t par_temp[10];
            //final->GetParameters(&par_temp[0]);
            //yield_fit[i] = par_temp[0];

            //yield_fit[i] = final->Integral(bmin_value, bmax_value);
            //background_fit[i] = myfunc->Integral(bmin_value, bmax_value);


            //cout << "Let's save now!" << endl;
            output->cd();

            char sub_name[200];
            sprintf(sub_name, "subtracted_%d_%d", i, j);
            subtracted[i][j] = new TCanvas(sub_name, sub_name, 1500, 800);
            V0Mass_nor_temp->Fit(third_poly, "R");
            V0Mass_nor_temp->Draw();
            straight_line->SetLineColor(kBlue);
            straight_line->Draw("LSAME");
            //subtracted[i]->Write();


            string c_name = "canvas_";
            char cc_name[200];
            TLegend *lam_leg = new TLegend(0.6, 0.6, 0.8, 0.8);
            Double_t line_height = V0Mass_nor[i][j]->GetBinContent(V0Mass_nor[i][j]->GetMaximumBin()) + 100;
            if(i == 0)
            {
                V0Mass_nor[i][j]->SetTitle("Lambda Mass Distribution");
                V0Mass_nor[i][j]->GetXaxis()->SetTitle("p #pi^{-} Invariant Mass (GeV/c^{2})");
                c_name.append(TString::Format("Lambda_%d", j));
                TLine *maxline = new TLine(bmax_value, 0, bmax_value, line_height);
                TLine *minline = new TLine(bmin_value, 0, bmin_value, line_height);

                lam_leg->SetHeader("(a)");
            }
            else
            {
                V0Mass_nor[i][j]->SetTitle("AntiLambda Mass Distribution");
                V0Mass_nor[i][j]->GetXaxis()->SetTitle("#bar{p} #pi^{+} Invariant Mass (GeV/c^{2})");
                c_name.append(TString::Format("AntiLambda_%d", j));
                TLine *maxline = new TLine(bmax_value, 0, bmax_value, line_height);
                TLine *minline = new TLine(bmin_value, 0, bmin_value, line_height);

                lam_leg->SetHeader("(b)");

            }
            strcpy(cc_name, c_name.c_str());
            canvasss[i][j] = new TCanvas(cc_name, cc_name, 1500, 800);
            //V0Mass_nor[i][j]->GetXaxis()->SetTitle("V0 Mass (GeV/c2)");
            V0Mass_nor[i][j]->GetYaxis()->SetTitle("Count");
            V0Mass_nor[i][j]->Draw();

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
            canvasss[i][j]->Write();

            //straight_line->Write();


            //GausDist->Write();

            // V0Mass_nor[i]->Fit(third_poly1, "R");
            // Double_t par1[4];
            // third_poly1->GetParameters(&par1[0]);
            // V0Mass_nor[i]->Fit(third_poly2, "R");
            // Double_t par2[4];
            // third_poly2->GetParameters(&par2[0]);
            // Double_t par3[4];
            // for(int l = 0; l < 4; l++) {
            //     par3[l] = (par1[l] + par2[l])/2.0;
            //     cout << "par3[l] = " <<  par3[l] << "par1[l] = " <<  par1[l] << "par2[l] = " <<  par2[l] << endl;
            // }
            // third_poly->SetParameters(par3);
            // third_poly->SetLineColor(kBlue);

            c_name.append("_fullfit");
            strcpy(cc_name, c_name.c_str());
            full_fit[i][j] = new TCanvas(cc_name, cc_name, 1500, 800);
            V0Mass_nor[i][j]->Draw();
            third_poly->Draw("LSAME");

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

            TLegend *leg2 = new TLegend(0.6, 0.6, 0.8, 0.8);
            leg2->AddEntry(maxline_low, "Mass Cuts for Background Bands", "l");
            leg2->AddEntry(third_poly, "Background Fit", "l");
            leg2->Draw();

            full_fit[i][j]->Write();
            //third_poly->Write();

            yield_fit[i][j] = 0;
            background_fit[i][j] = 0;

            for(int ddd = 0; ddd < (bmax - bmin) ; ddd++)
            {
                yield_fit[i][j] += V0Mass_nor[i][j]->GetBinContent(bmin + ddd);
                background_fit[i][j] += third_poly->Eval(bmin_value + ddd * 0.0007);
            }
            //background_fit[i] = third_poly->Integral(bmin_value, bmax_value);


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
        }
    }

    output->Close();

    //cout << "here?" << endl;

    //double yield_count[2] = {0.};
    //int numbins = bmax - bmin;

    ofstream myfile;
    myfile.open ("./Without_ptsplit/purity.txt");

    for(int i = 0; i < 2; i++)
    {
        for(int j = 0; j < 9; j++)
        {
            cout << "yield_count = " << yield_fit[i][j] << ", backgroundcount = " << background_fit[i][j] << endl;
            cout << "significance = " << (yield_fit[i][j] - background_fit[i][j]) / yield_fit[i][j] << endl;
            myfile << (yield_fit[i][j] - background_fit[i][j]) / yield_fit[i][j];
            myfile << " ";
            myfile << yield_fit[i][j];
            myfile << " ";
            //cout << "sqrt = " << yield_count[i]/(sqrt(yield_count[i] + backgroundcount[i]) * sqrt(RefMult[i].GetEntries())) << endl;

            myfile << "\n";
        }
        myfile << "\n";
        cout << "RefMult = " << RefMult[i].GetEntries() << endl;
    }

    myfile.close();
    //return (((yield_fit[0] - background_fit[0])/yield_fit[0]) + ((yield_fit[1] - background_fit[1])/yield_fit[1]))/2;
    // return (((yield_fit[0] - background_fit[0])) + ((yield_fit[1] - background_fit[1])))/2;
    return 2.0;
}
