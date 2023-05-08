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
#include <downstream_analysis.h>

using namespace std;

void downstream_analysis()
{
    // 1. Read in purity statistics for Lambdas and Protons 
    // 2. Create Files to store results across centralities
    initialize_macro(); 

    for (int i = min_cen; i <= max_cen; i++)
    {
        cout << "Centrality " << i << endl;
        downstream_analysis_bycen(i);
    }

    cout << "after loop delta[0] = " << delta[0] << endl;

    plotting_TGraph_Helper("TPC vs. EPD vs. EPD1 Resolution", "Centrality %", "Resolution", 9, cen9_eff, cen9_err_eff, reso_TPC, reso_TPC_err, reso_EPD, reso_EPD_err, reso_EPD1, reso_EPD1_err);
    organize_and_plot_v2_results();
    organize_and_plot_gamma_results();

    output_File2->Close();
    result_file.close();
    gamma_results.close();
}

void initialize_macro()
{
    read_sig();       // getting purity statistics by pT
    read_npart();     // getting npart statistics
    read_purity_Q2(); // getting purity statistics by Q2

    name_options3.push_back("TPC");
    name_options3.push_back("EPD");
    name_options3.push_back("EPD1");

    result_file.open(TString::Format("./%s_Results/output_v2.txt", method_name).Data());
    gamma_results.open(TString::Format("./%s_Results/gamma_results.txt", method_name).Data());
    output_File2 = new TFile(TString::Format("./%s_Results/output_cen%d-%d.root", method_name, min_cen, max_cen).Data(), "recreate");
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
        std::fstream myfile(TString::Format("./Fit_and_subtract/Q2_purity/%s/pion/output_cen%d_pion.txt", method_name, curr_c), std::ios_base::in);

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

        std::fstream myfile2(TString::Format("./Fit_and_subtract/Q2_purity/%s/pion/output_cen%d_pion_anti.txt", method_name, curr_c), std::ios_base::in);

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

void organize_and_plot_v2_results()
{

    float v2_averaged_TPC[9] = {0.}, v2_averaged_EPD[9] = {0.}, v2_averaged_EPD1[9] = {0.}, v2_averaged_TPC_err[9] = {0.}, v2_averaged_EPD_err[9] = {0.}, v2_averaged_EPD1_err[9] = {0.};

    for (int w = 0; w < 3; w++)
    {
        result_file << "const float v2_averaged_" << name_options3.at(w) << "[9] = {";

        for (int q = 0; q <= (max_cen - min_cen); q++)
        {
            if (q != (max_cen - min_cen))
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
            if (q != (max_cen - min_cen))
                result_file << parent_v2_vect[w].at(2 * q) << ", ";
            else
                result_file << parent_v2_vect[w].at(2 * q) << "}; \n";
        }

        result_file << "const float v2_averaged_pairpion_" << name_options3.at(w) << "[9] = {";

        for (int q = 0; q <= (max_cen - min_cen); q++)
        {
            if (q != (max_cen - min_cen))
                result_file << pionpion_v2_vect[w].at(2 * q) << ", ";
            else
                result_file << pionpion_v2_vect[w].at(2 * q) << "}; \n";
        }
    }

    plotting_TGraph_Helper("TPC vs. EPD vs. EPD1 Single Charged Particles v_{2}", "Centrality %", "Resolution", 9, cen9_eff, cen9_err_eff, v2_averaged_TPC, v2_averaged_TPC_err, v2_averaged_EPD, v2_averaged_EPD_err, v2_averaged_EPD1, v2_averaged_EPD1_err);
}

void organize_and_plot_gamma_results()
{
    gamma_results << "Centrality, Gamma112_EPD, Gamma112_EPD_err, Gamma132_EPD, Gamma132_EPD_err, Gamma112_EPD1, Gamma112_EPD1_err, Gamma132_EPD1, Gamma132_EPD1_err, ";
    gamma_results << "Gamma112_EPD_ESE, Gamma112_EPD_ESE_err, Gamma132_EPD_ESE, Gamma132_EPD_ESE_err, Gamma112_EPD1_ESE, Gamma112_EPD1_ESE_err, Gamma132_EPD1_ESE, Gamma132_EPD1_ESE_err, ";
    gamma_results << "kappa112_EPD, kappa112_EPD_err, kappa132_EPD, kappa132_EPD_err, kappa112_EPD1, kappa112_EPD1_err, kappa132_EPD1, kappa132_EPD1_err, ";
    gamma_results << "kappa112_EPD_ESE, kappa112_EPD_ESE_err, kappa132_EPD_ESE, kappa132_EPD_ESE_err, kappa112_EPD1_ESE, kappa112_EPD1_ESE_err, kappa132_EPD1_ESE, kappa132_EPD1_ESE_err\n";

    float gamma112_TPC[9] = {0.}, gamma112_EPD[9] = {0.}, gamma112_EPD1[9] = {0.};
    float gamma112_TPC_err[9] = {0.}, gamma112_EPD_err[9] = {0.}, gamma112_EPD1_err[9] = {0.};
    float gamma132_TPC[9] = {0.}, gamma132_EPD[9] = {0.}, gamma132_EPD1[9] = {0.};
    float gamma132_TPC_err[9] = {0.}, gamma132_EPD_err[9] = {0.}, gamma132_EPD1_err[9] = {0.};

    float gamma112_ensemble_TPC[9] = {0.}, gamma112_ensemble_EPD[9] = {0.}, gamma112_ensemble_EPD1[9] = {0.};
    float gamma112_ensemble_TPC_err[9] = {0.}, gamma112_ensemble_EPD_err[9] = {0.}, gamma112_ensemble_EPD1_err[9] = {0.};
    float gamma132_ensemble_TPC[9] = {0.}, gamma132_ensemble_EPD[9] = {0.}, gamma132_ensemble_EPD1[9] = {0.};
    float gamma132_ensemble_TPC_err[9] = {0.}, gamma132_ensemble_EPD_err[9] = {0.}, gamma132_ensemble_EPD1_err[9] = {0.};

    float v2_lam_TPC[9] = {0.}, v2_p_TPC[9] = {0.}, v2_lam_TPC_err[9] = {0.}, v2_p_TPC_err[9] = {0.};
    float v2_lam_EPD[9] = {0.}, v2_p_EPD[9] = {0.}, v2_lam_EPD_err[9] = {0.}, v2_p_EPD_err[9] = {0.};
    float v2_lam_EPD1[9] = {0.}, v2_p_EPD1[9] = {0.}, v2_lam_EPD1_err[9] = {0.}, v2_p_EPD1_err[9] = {0.};

    float kappa112_TPC[9] = {0.}, kappa112_TPC_err[9] = {0.}, kappa112_EPD[9] = {0.}, kappa112_EPD_err[9] = {0.}, kappa112_EPD1[9] = {0.}, kappa112_EPD1_err[9] = {0.};
    float kappa132_TPC[9] = {0.}, kappa132_TPC_err[9] = {0.}, kappa132_EPD[9] = {0.}, kappa132_EPD_err[9] = {0.}, kappa132_EPD1[9] = {0.}, kappa132_EPD1_err[9] = {0.};
    float kappa112_TPC_ESE[9] = {0.}, kappa112_TPC_ESE_err[9] = {0.}, kappa112_EPD_ESE[9] = {0.}, kappa112_EPD_ESE_err[9] = {0.}, kappa112_EPD1_ESE[9] = {0.}, kappa112_EPD1_ESE_err[9] = {0.};
    float kappa132_TPC_ESE[9] = {0.}, kappa132_TPC_ESE_err[9] = {0.}, kappa132_EPD_ESE[9] = {0.}, kappa132_EPD_ESE_err[9] = {0.}, kappa132_EPD1_ESE[9] = {0.}, kappa132_EPD1_ESE_err[9] = {0.};

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

        gamma_results << j << ", " << gamma_ensemble[0][j][1] << ", " << gamma_ensemble_err[0][j][1] << ", " << gamma_ensemble[1][j][1] << ", " << gamma_ensemble_err[1][j][1] << ", " << gamma_ensemble[0][j][2] << ", " << gamma_ensemble_err[0][j][2] << ", " << gamma_ensemble[1][j][2] << ", " << gamma_ensemble_err[1][j][2] << ", ";
        gamma_results << gamma[0][j][1] << ", " << gamma_err[0][j][1] << ", " << gamma[1][j][1] << ", " << gamma_err[1][j][1] << ", " << gamma[0][j][2] << ", " << gamma_err[0][j][2] << ", " << gamma[1][j][2] << ", " << gamma_err[1][j][2] << ", ";

        v2_lam_TPC[j] = v2_lam[j][0];
        v2_lam_EPD[j] = v2_lam[j][1];
        v2_lam_EPD1[j] = v2_lam[j][2];
        v2_lam_TPC_err[j] = v2_lam_err[j][0];
        v2_lam_EPD_err[j] = v2_lam_err[j][1];
        v2_lam_EPD1_err[j] = v2_lam_err[j][2];

        v2_p_TPC[j] = v2_p[j][0];
        v2_p_EPD[j] = v2_p[j][1];
        v2_p_EPD1[j] = v2_p[j][2];
        v2_p_TPC_err[j] = v2_p_err[j][0];
        v2_p_EPD_err[j] = v2_p_err[j][1];
        v2_p_EPD1_err[j] = v2_p_err[j][2];

        // cout << "delta[" << j << "] = " << delta[j] << endl;
        // cout << "sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)) = " << sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)) << endl;

        kappa112_EPD[j] = gamma_ensemble[0][j][1] / (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)));
        float tmp_error1 = error_mult(v2_lam[j][1], v2_lam[j][1], v2_lam_err[j][1], v2_lam_err[j][1]);
        float tmp_error2 = error_mult(v2_p[j][1], v2_p[j][1], v2_p_err[j][1], v2_p_err[j][1]);
        float tmp_error3 = error_add(tmp_error1, tmp_error2);
        float tmp_error4 = error_mult(delta[j], sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)), delta_err[j], tmp_error3);
        kappa112_EPD_err[j] = error_divide(gamma_ensemble[0][j][1], (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2))), gamma_ensemble_err[0][j][1], tmp_error4);

        kappa132_EPD[j] = gamma_ensemble[1][j][1] / (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)));
        kappa132_EPD_err[j] = error_divide(gamma_ensemble[1][j][1], (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2))), gamma_ensemble_err[1][j][1], tmp_error4);

        kappa112_EPD_ESE[j] = gamma[0][j][1] / (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)));
        kappa112_EPD_ESE_err[j] = error_divide(gamma[0][j][1], (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2))), gamma_err[0][j][1], tmp_error4);

        kappa132_EPD_ESE[j] = gamma[1][j][1] / (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2)));
        kappa132_EPD_ESE_err[j] = error_divide(gamma[1][j][1], (delta[j] * sqrt(pow(v2_lam[j][1], 2) + pow(v2_p[j][1], 2))), gamma_err[1][j][1], tmp_error4);

        kappa112_EPD1[j] = gamma_ensemble[0][j][2] / (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2)));
        tmp_error1 = error_mult(v2_lam[j][2], v2_lam[j][2], v2_lam_err[j][2], v2_lam_err[j][2]);
        tmp_error2 = error_mult(v2_p[j][2], v2_p[j][2], v2_p_err[j][2], v2_p_err[j][2]);
        tmp_error3 = error_add(tmp_error1, tmp_error2);
        tmp_error4 = error_mult(delta[j], sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2)), delta_err[j], tmp_error3);
        kappa112_EPD1_err[j] = error_divide(gamma_ensemble[0][j][2], (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2))), gamma_ensemble_err[0][j][2], tmp_error4);

        kappa132_EPD1[j] = gamma_ensemble[1][j][2] / (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2)));
        kappa132_EPD1_err[j] = error_divide(gamma_ensemble[1][j][2], (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2))), gamma_ensemble_err[1][j][2], tmp_error4);

        kappa112_EPD1_ESE[j] = gamma[0][j][2] / (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2)));
        kappa112_EPD1_ESE_err[j] = error_divide(gamma[0][j][2], (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2))), gamma_err[0][j][2], tmp_error4);

        kappa132_EPD1_ESE[j] = gamma[1][j][2] / (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2)));
        kappa132_EPD1_ESE_err[j] = error_divide(gamma[1][j][2], (delta[j] * sqrt(pow(v2_lam[j][2], 2) + pow(v2_p[j][2], 2))), gamma_err[1][j][2], tmp_error4);

        gamma_results << kappa112_EPD[j] << ", " << kappa112_EPD_err[j] << ", " << kappa132_EPD[j] << ", " << kappa132_EPD_err[j] << ", " << kappa112_EPD1[j] << ", " << kappa112_EPD1_err[j] << ", " << kappa132_EPD1[j] << ", " << kappa132_EPD1_err[j] << ", ";
        gamma_results << kappa112_EPD_ESE[j] << ", " << kappa112_EPD_ESE_err[j] << ", " << kappa132_EPD_ESE[j] << ", " << kappa132_EPD_ESE_err[j] << ", " << kappa112_EPD1_ESE[j] << ", " << kappa112_EPD1_ESE_err[j] << ", " << kappa132_EPD1_ESE[j] << ", " << kappa132_EPD1_ESE_err[j] << "\n";
    }

    // plotting_TGraph_Helper("Compare ESE vs. Ensemble TPC #Delta#gamma_{112} * N_{part}", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_TPC, gamma112_TPC_err, gamma112_ensemble_TPC, gamma112_ensemble_TPC_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD #Delta#gamma_{112} * N_{part}", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD, gamma112_EPD_err, gamma112_ensemble_EPD, gamma112_ensemble_EPD_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD1 #Delta#gamma_{112} * N_{part}", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD1, gamma112_EPD1_err, gamma112_ensemble_EPD1, gamma112_ensemble_EPD1_err);
    // plotting_TGraph_Helper("Compare ESE vs. Ensemble TPC #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_TPC, gamma132_TPC_err, gamma132_ensemble_TPC, gamma132_ensemble_TPC_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_EPD, gamma132_EPD_err, gamma132_ensemble_EPD, gamma132_ensemble_EPD_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD1 #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_EPD1, gamma132_EPD1_err, gamma132_ensemble_EPD1, gamma132_ensemble_EPD1_err);

    // plotting_TGraph_Helper("ESE: TPC #Delta#gamma_{112} * N_{part} vs. #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_TPC, gamma112_TPC_err, gamma132_TPC, gamma132_TPC_err);
    plotting_TGraph_Helper("ESE: EPD #Delta#gamma_{112} * N_{part} vs. #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD, gamma112_EPD_err, gamma132_EPD, gamma132_EPD_err);
    plotting_TGraph_Helper("ESE: EPD1 #Delta#gamma_{112} * N_{part} vs. #Delta#gamma_{132} * N_{part}", "Centrality %", "#Delta#gamma * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD1, gamma112_EPD1_err, gamma132_EPD1, gamma132_EPD1_err);
    // plotting_TGraph_Helper("ESE: #Delta#gamma_{112} * N_{part} TPC vs. EPD vs. EPD1", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_TPC, gamma112_TPC_err, gamma112_EPD, gamma112_EPD_err, gamma112_EPD1, gamma112_EPD1_err);
    // plotting_TGraph_Helper("ESE: #Delta#gamma_{132} * N_{part} TPC vs. EPD vs. EPD1", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_TPC, gamma132_TPC_err, gamma132_EPD, gamma132_EPD_err, gamma132_EPD1, gamma132_EPD1_err);
    plotting_TGraph_Helper("ESE: #Delta#gamma_{112} * N_{part} EPD vs. EPD1", "Centrality %", "#Delta#gamma_{112} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma112_EPD, gamma112_EPD_err, gamma112_EPD1, gamma112_EPD1_err);
    plotting_TGraph_Helper("ESE: #Delta#gamma_{132} * N_{part} EPD vs. EPD1", "Centrality %", "#Delta#gamma_{132} * N_{part}", 9, cen9_eff, cen9_err_eff, gamma132_EPD, gamma132_EPD_err, gamma132_EPD1, gamma132_EPD1_err);

    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD #kappa_{112}", "Centrality %", "#kappa_{112}", 9, cen9_eff, cen9_err_eff, kappa112_EPD_ESE, kappa112_EPD_ESE_err, kappa112_EPD, kappa112_EPD_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD1 #kappa_{112}", "Centrality %", "#kappa_{112}", 9, cen9_eff, cen9_err_eff, kappa112_EPD1_ESE, kappa112_EPD1_ESE_err, kappa112_EPD1, kappa112_EPD1_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD #kappa_{132}", "Centrality %", "#kappa_{132}", 9, cen9_eff, cen9_err_eff, kappa132_EPD_ESE, kappa132_EPD_ESE_err, kappa132_EPD, kappa132_EPD_err);
    plotting_TGraph_Helper("Compare ESE vs. Ensemble EPD1 #kappa_{132}", "Centrality %", "#kappa_{132}", 9, cen9_eff, cen9_err_eff, kappa132_EPD1_ESE, kappa132_EPD1_ESE_err, kappa132_EPD1, kappa132_EPD1_err);

    plotting_TGraph_Helper("Ensemble: EPD #kappa_{112} vs. #kappa_{132}", "Centrality %", "#kappa", 9, cen9_eff, cen9_err_eff, kappa112_EPD, kappa112_EPD_err, kappa132_EPD, kappa132_EPD_err);
    plotting_TGraph_Helper("Ensemble: EPD1 #kappa_{112} vs. #kappa_{132}", "Centrality %", "#kappa", 9, cen9_eff, cen9_err_eff, kappa112_EPD1, kappa112_EPD1_err, kappa132_EPD1, kappa132_EPD1_err);
    plotting_TGraph_Helper("Ensemble: #kappa_{112} EPD vs. EPD1", "Centrality %", "#kappa_{112}", 9, cen9_eff, cen9_err_eff, kappa112_EPD, kappa112_EPD_err, kappa112_EPD1, kappa112_EPD1_err);
    plotting_TGraph_Helper("Ensemble: #kappa_{132} EPD vs. EPD1", "Centrality %", "#kappa_{132}", 9, cen9_eff, cen9_err_eff, kappa132_EPD, kappa132_EPD_err, kappa132_EPD1, kappa132_EPD1_err);

    plotting_TGraph_Helper("ESE: EPD #kappa_{112} vs. #kappa_{132}", "Centrality %", "#kappa", 9, cen9_eff, cen9_err_eff, kappa112_EPD_ESE, kappa112_EPD_ESE_err, kappa132_EPD_ESE, kappa132_EPD_ESE_err);
    plotting_TGraph_Helper("ESE: EPD1 #kappa_{112} vs. #kappa_{132}", "Centrality %", "#kappa", 9, cen9_eff, cen9_err_eff, kappa112_EPD1_ESE, kappa112_EPD1_ESE_err, kappa132_EPD1_ESE, kappa132_EPD1_ESE_err);
    plotting_TGraph_Helper("ESE: #kappa_{112} EPD vs. EPD1", "Centrality %", "#kappa_{112}", 9, cen9_eff, cen9_err_eff, kappa112_EPD_ESE, kappa112_EPD_ESE_err, kappa112_EPD1_ESE, kappa112_EPD1_ESE_err);
    plotting_TGraph_Helper("ESE: #kappa_{132} EPD vs. EPD1", "Centrality %", "#kappa_{132}", 9, cen9_eff, cen9_err_eff, kappa132_EPD_ESE, kappa132_EPD_ESE_err, kappa132_EPD1_ESE, kappa132_EPD1_ESE_err);

    plotting_TGraph_Helper("Delta", "Centrality %", "delta", 9, cen9_eff, cen9_err_eff, delta, delta_err);
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
    // 1. Create Files to store Centrality Specific Results
    // 2. Store Resolution info
    initialize_each_cen(cen);

    organize_ensemble_average_results(cen);
    obtain_v2_average(cen);
    // old_ptsplit_method(cen);

    for (int j = 1; j < 3; j++)
    {
        Q2_parent_results(cen, j, 0);
        Q2_parent_results(cen, j, 1);
        cout << "done with " << j << endl;
    }

    file->Close();
    output_File->Close();
}

void initialize_each_cen(int cen)
{
    file = new TFile(TString::Format("./%s/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", method_name, cen).Data());
    output_File = new TFile(TString::Format("./%s_Results/output_cen%d.root", method_name, cen).Data(), "recreate");

    extract_reso();

    reso_TPC[cen] = reso[0];
    reso_EPD[cen] = reso[1];
    reso_EPD1[cen] = reso[2];
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

    reso[3] = 1;

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

double resEventPlane(double chi)
{
    double con = 0.626657;
    double arg = chi * chi / 4.;

    Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

    return res;
}

void organize_ensemble_average_results(int cen)
{
    // Without Efficiency Corrections
    // profile_divide_by_reso("Parity_int_obs", cen, false);
    // profile_divide_by_reso("Parity_int_ss_obs", cen, false);
    // profile_divide_by_reso("Delta_int_ss_obs", cen, false);
    
    //  With Efficiency Corrections
    profile_divide_by_reso("Parity_int_obs", cen, true);
    profile_divide_by_reso("Parity_int_ss_obs", cen, true);
    profile_divide_by_reso("Delta_int_ss_obs", cen, true);

    // With Efficiency Corrections and QQ > 0.8 Cut
    profile_divide_by_reso("Parity_int_obs", cen, true, true);
    profile_divide_by_reso("Parity_int_ss_obs", cen, true, true);
    profile_divide_by_reso("Delta_int_ss_obs", cen, true, true);
}

void profile_divide_by_reso(const char *profile_name, int cen, bool eff_corr, bool QQ_cut)
{
    // grabs lam and antilam results, correct by background, combined, and store resolution corrected results in array
    
    const int n_entries = 25;
    char *qq_cut = (QQ_cut) ? "_QQcut" : "";
    int num = (eff_corr) ? 3 : 1;
    cout << "profile_divide_by_reso Profile Name = " << TString::Format("%s%d%s", profile_name, num, qq_cut) << endl;

    //////////////////////////// Lambda Results ////////////////////////////
    TProfile *temp_profile = (TProfile *)file->Get(TString::Format("%s%d%s", profile_name, num, qq_cut));
    TProfile *temp_profile_rot = (TProfile *)file->Get(TString::Format("%s%d%s_rot", profile_name, num, qq_cut));
    int n_bins = temp_profile->GetNbinsX();
    double x[n_entries] = {0.}, x_err[n_entries] = {0.}, y[n_entries] = {0.}, y_err[n_entries] = {0.};
    temp_profile->Sumw2();
    temp_profile_rot->Sumw2();

    temp_profile->Add(temp_profile, temp_profile_rot, 1.0 / lam_purity[cen][17], -(1.0 - lam_purity[cen][17]) / lam_purity[cen][17]);
    
    // To get errors that are not affected by efficiency corrections
    TProfile *temp_profile_err = (TProfile *)file->Get(TString::Format("%s%d%s", profile_name, 1, qq_cut));
    TProfile *temp_profile_err_rot = (TProfile *)file->Get(TString::Format("%s%d%s_rot", profile_name, 1, qq_cut));
    temp_profile_err->Sumw2();
    temp_profile_err_rot->Sumw2();

    if (eff_corr) temp_profile_err->Add(temp_profile_err, temp_profile_err_rot, 1.0 / lam_purity[cen][17], -(1.0 - lam_purity[cen][17]) / lam_purity[cen][17]);
    for (int i = 1; i <= temp_profile->GetNbinsX(); i++)
    {
        x[i - 1] = i;
        y[i - 1] = temp_profile->GetBinContent(i);
        y_err[i - 1] = (eff_corr)? temp_profile_err->GetBinError(i) : temp_profile->GetBinError(i);
    }

    TGraphErrors *temp_graph = new TGraphErrors(n_bins, x, y, x_err, y_err);

    //////////////////////////// End of Lambda Results ////////////////////////////    

    //////////////////////////// AntiLambda Results ////////////////////////////
    TProfile *temp_profile_anti = (TProfile *)file->Get(TString::Format("%s%d%s_anti", profile_name, num, qq_cut));
    TProfile *temp_profile_anti_rot = (TProfile *)file->Get(TString::Format("%s%d%s_anti_rot", profile_name, num, qq_cut));
    int n_bins_anti = temp_profile_anti->GetNbinsX();
    double x_anti[n_entries] = {0.}, x_anti_err[n_entries] = {0.}, y_anti[n_entries] = {0.}, y_anti_err[n_entries] = {0.};
    temp_profile_anti->Sumw2();
    temp_profile_anti_rot->Sumw2();

    temp_profile_anti->Add(temp_profile_anti, temp_profile_anti_rot, 1.0 / antilam_purity[cen][17], -(1.0 - antilam_purity[cen][17]) / antilam_purity[cen][17]);
    
    // To get errors that are not affected by efficiency corrections
    TProfile *temp_profile_anti_err = (TProfile *)file->Get(TString::Format("%s%d%s_anti", profile_name, 1, qq_cut));
    TProfile *temp_profile_err_anti_rot = (TProfile *)file->Get(TString::Format("%s%d%s_anti_rot", profile_name, 1, qq_cut));
    temp_profile_anti_err->Sumw2();
    temp_profile_err_anti_rot->Sumw2();

    if (eff_corr) temp_profile_anti_err->Add(temp_profile_anti_err, temp_profile_err_anti_rot, 1.0 / antilam_purity[cen][17], -(1.0 - antilam_purity[cen][17]) / antilam_purity[cen][17]);
    for (int i = 1; i <= temp_profile_anti->GetNbinsX(); i++)
    {
        x_anti[i - 1] = i;
        y_anti[i - 1] = temp_profile_anti->GetBinContent(i);
        y_anti_err[i - 1] = (eff_corr)? temp_profile_anti_err->GetBinError(i) : temp_profile_anti->GetBinError(i);
    }

    TGraphErrors *temp_graph_anti = new TGraphErrors(n_bins_anti, x_anti, y_anti, x_anti_err, y_anti_err);

    //////////////////////////// End of AntiLambda Results ////////////////////////////

    //////////////////////////// Combining Lambda and AntiLambda Results ////////////////////////////
    temp_profile->Add(temp_profile_anti);
    double x_combined[n_entries] = {0.}, x_err_combined[n_entries] = {0.}, y_combined[n_entries] = {0.}, y_err_combined[n_entries] = {0.};
    if (eff_corr) temp_profile_err->Add(temp_profile_anti_err);
    for (int i = 1; i <= temp_profile->GetNbinsX(); i++)
    {
        x_combined[i - 1] = i;
        y_combined[i - 1] = temp_profile->GetBinContent(i);
        y_err_combined[i - 1] = (eff_corr)? temp_profile_err->GetBinError(i) : temp_profile->GetBinError(i);
    }

    TGraphErrors *temp_graph_combined = new TGraphErrors(n_bins, x_combined, y_combined, x_err_combined, y_err_combined);

    //////////////////////////// End of Combining Lambda and AntiLambda Results ////////////////////////////

    //////////////////////////// Recording for Cross-Centrality Analysis ////////////////////////////
    for (int tp1 = 0; tp1 < 3; tp1++)
    {
        if (QQ_cut) continue;
        
        // Gamma132 Results
        if (profile_name == "Parity_int_obs")
        {
            // Value
            gamma_ensemble[1][cen][tp1] = (float)((y_combined[4 + (8 * tp1)] - y_combined[3 + (8 * tp1)]) + (y_combined[8 + (8 * tp1)] - y_combined[7 + (8 * tp1)])) / 2.0;
            if (tp1 == 0) gamma_ensemble[1][cen][tp1] = gamma_ensemble[1][cen][tp1] / reso1[tp1] / 100.0;
            else gamma_ensemble[1][cen][tp1] = gamma_ensemble[1][cen][tp1] / reso[tp1] / 100.0;

            // Error
            float tmp_error1 = error_add(y_err_combined[4 + (8 * tp1)], y_err_combined[3 + (8 * tp1)]);
            float tmp_erorr2 = error_add(y_err_combined[8 + (8 * tp1)], y_err_combined[7 + (8 * tp1)]);
            gamma_ensemble_err[1][cen][tp1] = (float)error_add(tmp_error1, tmp_erorr2) / 2.0;
            if (tp1 == 0) gamma_ensemble_err[1][cen][tp1] = gamma_ensemble_err[1][cen][tp1] / reso1[tp1] / 100.0;
            else gamma_ensemble_err[1][cen][tp1] = gamma_ensemble_err[1][cen][tp1] / reso[tp1] / 100.0;
        }

        // Gamma112 Results
        if (profile_name == "Parity_int_ss_obs")
        {
            // Value
            gamma_ensemble[0][cen][tp1] = y_combined[4 + (4 * tp1)] - y_combined[3 + (4 * tp1)];
            if (tp1 == 0) gamma_ensemble[0][cen][tp1] = gamma_ensemble[0][cen][tp1] / reso1[tp1] / 100.0;
            else gamma_ensemble[0][cen][tp1] = gamma_ensemble[0][cen][tp1] / reso[tp1] / 100.0;
            
            // Error
            gamma_ensemble_err[0][cen][tp1] = error_add(y_err_combined[4 + (4 * tp1)], y_err_combined[3 + (4 * tp1)]);
            if (tp1 == 0) gamma_ensemble_err[0][cen][tp1] = gamma_ensemble_err[0][cen][tp1] / reso1[tp1] / 100.0;
            else gamma_ensemble_err[0][cen][tp1] = gamma_ensemble_err[0][cen][tp1] / reso[tp1] / 100.0;
        }
    }

    output_File->cd();
    temp_graph->Write(TString::Format("%s%d", profile_name, num));
    temp_graph_anti->Write(TString::Format("%s%d_anti", profile_name, num));
    temp_graph_combined->Write(TString::Format("%s%d_combined", profile_name, num));

    if ((profile_name == "Delta_int_ss_obs"))
    {
        if (!QQ_cut)
        {
            delta[cen] = y_combined[4] - y_combined[3];
            delta[cen] = delta[cen] / 100.0;
            delta_err[cen] = error_add(y_err_combined[4], y_err_combined[3]);
            delta_err[cen] = delta_err[cen] / 100.0;
        }
    }
}

void obtain_v2_average(int cen)
{
    TProfile *v2_single_lamda = (TProfile *)file->Get("v2_single_lamda");
    TProfile *v2_single_proton = (TProfile *)file->Get("v2_single_proton");

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

        // sprintf(fname, "Hist_v2pion_pt_%s_obs5", name_options3.at(jk).Data());
        // temp_vec = Rebin_v2_Data(fname, fname, cen, jk);
        // pionpion_v2_vect[jk].push_back(temp_vec.at(1) / 100.0);
        // pionpion_v2_vect[jk].push_back(temp_vec.at(2) / 100.0);
        // temp_vec.clear();

        sprintf(fname, "Hist_v2pion_pt_%s_obs5", name_options3.at(jk).Data());
        cout << fname << endl;
        temp_vec = Rebin_v2_Data(fname, fname, cen, 3);
        cout << "sqrt(temp_vec.at(1)) = " << sqrt(fabs(temp_vec.at(1))) << endl;
        pionpion_v2_vect[jk].push_back(sqrt(fabs(temp_vec.at(1))));
        pionpion_v2_vect[jk].push_back(temp_vec.at(2));
        // cout << temp_vec.at(1) / 100.0 << endl;
        temp_vec.clear();

        v2_lam[cen][jk] = v2_single_lamda->GetBinContent(jk);
        v2_lam_err[cen][jk] = v2_single_lamda->GetBinError(jk);
        v2_p[cen][jk] = v2_single_proton->GetBinContent(jk);
        v2_p_err[cen][jk] = v2_single_proton->GetBinError(jk);
    }

    // Lambda: Hist_v2_pt_obs2, Hist_v2_pt_EPD1_obs, Hist_v2_pt_EPD_obs
    temp_vec = Rebin_v2_Data("Hist_v2_pt_obs2", "Hist_v2_pt_obs2", cen, 0);
    lam_v2_vect.push_back(temp_vec.at(1) / 100.0);
    lam_v2_vect.push_back(temp_vec.at(2) / 100.0);
    temp_vec.clear();
    // Proton: Hist_v2_pt_obs2_p
    temp_vec = Rebin_v2_Data("Hist_v2_pt_obs2_p", "Hist_v2_pt_obs2_p", cen, 0);
    proton_v2_vect.push_back(temp_vec.at(1) / 100.0);
    proton_v2_vect.push_back(temp_vec.at(2) / 100.0);
    temp_vec.clear();

    // profile_divide_by_reso_fit_by_const("Hist_v2parent_eta_obs5", cen);
    // profile_divide_by_reso_fit_by_const("Hist_v2parent_eta_obs5_rot", cen);
}

void profile_divide_by_reso_fit_by_const(const char *profile_name, int cen)
{
    cout << profile_name << endl;

    TProfile *temp_profile = (TProfile *)file->Get(profile_name);

    temp_profile->Scale(1.0 / reso);
    // cout << "reso = " << reso << endl;
    TF1 *temp_func = new TF1("temp_func", "[0]", -5, 5);
    temp_profile->Fit(temp_func, "QR");
    result_file << profile_name << " Fit Parameter = " << temp_func->GetParameter(0) << "\n";

    output_File->cd();
    temp_profile->Write(profile_name);
}

void Q2_parent_results(int cen, int ep_option, int option112)
{

    bool pion_results = true;

    //////////////// Setting Up for Function ////////////////
    // float Q2_range[3][9] = {{3.5, 5, 5, 5, 5, 5, 3.5, 1.5, 0.6}, {2, 3.5, 5, 5, 5, 5, 5, 5, 5}, {5, 5, 5, 5, 5, 5, 5, 5, 5}};
    float Q2_range[3][9] = {{3.5, 5, 5, 5, 5, 5, 3.5, 1.5, 0.6}, {5, 5, 5, 5, 10, 10, 10, 10, 10}, {5, 5, 5, 5, 10, 10, 10, 10, 10}};
    // int rebin[3][9] = {{5, 5, 5, 5, 5, 5, 5, 2, 1}, {5, 5, 5, 5, 5, 5, 5, 5, 5}, {5, 5, 5, 5, 5, 5, 5, 5, 5}};
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
    double Xrange_delGamma[9] = {1.2, 1.2, 0.6, 1.0, 0.5, 0.4, 0.4, 0.4, 0.3};
    double Xfitting_delGamma[3][9] = {{0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}, {0.4, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}, {0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2}};
    // double Xfitting_delGamma_min[2][3][9] = {{{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0.035, 0, 0, 0, 0, 0, 0, 0}, {0, 0.02, 0, 0, 0, 0, 0, 0, 0}},
    //                                          {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.01, 0, 0, 0, 0, 0, 0, 0, 0}}};
    double Xfitting_delGamma_min[2][3][9] = {{{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}},
                                             {{0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0}}};
    char printname5[200];
    sprintf(printname5, "cen%d.ese.parent.Q2range%d.rebin%d.%s.%d.pdf", cen, Q2_range[ep_option][cen], rebin[ep_option][cen], ep_option_name[ep_option].Data(), option112);
    ////////////////////////////////////////////////////////////////////////////////

    //////////////////////////////// Analysis ////////////////////////////////

    //////// Resolution: p_cos_Q2, p_cos_Q2_EPD, p_cos_Q2_EPD1
    if (!pion_results)
        TProfile *p_res = (TProfile *)file->Get(TString::Format("p_cos_Q2%s", ep_option_name[ep_option].Data()).Data());
    if (pion_results)
        TProfile *p_res = (TProfile *)file->Get(TString::Format("p_cos_Q2%s_pion", ep_option_name[ep_option].Data()).Data());
    TF1 *f_res = new TF1("f_res", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, Q2_range[ep_option][cen]);
    p_res->Fit("f_res", "0Q", "", 0, Q2_range[ep_option][cen]);

    cout << "resolution loaded" << endl;

    //////// Obtaining v2 plots: 1st Set - Lambda-Proton pair v2; 2nd set - charged-particles-pair v2; 3rd set - single charged particle v2; 4th set - pion-pair v2

    // TProfile *p_v2_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2parente_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2parente_Q2parent%s_obs2", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2parentw_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2parentw_Q2parent%s_obs2", ep_option_name[ep_option].Data()));

    // TProfile *p_v2_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2e_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2e_Q2parent%s_obs2", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2w_Q2parent%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2w_Q2parent%s_obs2", ep_option_name[ep_option].Data()));

    // TProfile *p_v2_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2e_Q2%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2e_Q2%s_obs2", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs1 = (TProfile *)file->Get(TString::Format("p_v2w_Q2%s_obs1", ep_option_name[ep_option].Data()));
    // TProfile *p_v2w_Q_obs2 = (TProfile *)file->Get(TString::Format("p_v2w_Q2%s_obs2", ep_option_name[ep_option].Data()));

    TProfile *p_v2_Q_obs1 = (TProfile *)file->Get(TString::Format("p_pionv2e_Q2%s_obs1", ep_option_name[ep_option].Data()));
    TProfile *p_v2_Q_obs2 = (TProfile *)file->Get(TString::Format("p_pionv2e_Q2%s_obs2", ep_option_name[ep_option].Data()));
    TProfile *p_v2w_Q_obs1 = (TProfile *)file->Get(TString::Format("p_pionv2w_Q2%s_obs1", ep_option_name[ep_option].Data()));
    TProfile *p_v2w_Q_obs2 = (TProfile *)file->Get(TString::Format("p_pionv2w_Q2%s_obs2", ep_option_name[ep_option].Data()));

    p_v2_Q_obs1->Scale(1);
    if (ep_option <= 1)
        p_v2_Q_obs1->Add(p_v2w_Q_obs1);
    p_v2_Q_obs2->Scale(1);
    if (ep_option <= 1)
        p_v2_Q_obs2->Add(p_v2w_Q_obs2);
    cout << "after addition p_v2_Q_obs1 entries = " << p_v2_Q_obs1->GetEntries() << endl;

    TH1 *p_v2_Q_obs_new = (TH1 *)p_v2_Q_obs1->ProjectionX(TString::Format("p_v2_Q_obs_new%s%d", ep_option_name[ep_option].Data(), option112));

    int Nbin = p_v2_Q_obs1->GetNbinsX();
    // cout << "Nbin = " << Nbin << endl;
    for (int i = 0; i < Nbin; i++)
    {
        // cout << "i = " << i << endl;
        float cont = p_v2_Q_obs2->GetBinContent(i + 1);
        float err = p_v2_Q_obs1->GetBinError(i + 1);
        float res_q = 1;

        if (p_v2_Q_obs2->GetBinCenter(i + 1) < Q2_range[ep_option][cen])
        {
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(p_v2_Q_obs2->GetBinCenter(i + 1)));

            int i_tmp = i + 2;
            int i_tmp2 = i;

            while (res_q <= 0)
            {
                res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(p_v2_Q_obs2->GetBinCenter(i_tmp)));
                i_tmp++;

                if ((res_q <= 0) && (i_tmp2 > 0)){
                    res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(p_v2_Q_obs2->GetBinCenter(i_tmp2)));
                    i_tmp2--;
                }

                if (i_tmp > (i+7)) break;
            }
        }
        else{
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]));

            int i_tmp = p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]) + 1;
            int i_tmp2 = p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]) - 1;

            while (res_q <= 0)
            {
                res_q = p_res->GetBinContent(i_tmp);
                i_tmp++;

                if ((res_q <= 0) && (i_tmp2 > 0)){
                    res_q = p_res->GetBinContent(i_tmp2);
                    i_tmp2--;
                }

                if (i_tmp > (i+7)) break;
            }
        }

        if (ep_option <= 1)
            res_q = sqrt(res_q);

        p_v2_Q_obs_new->SetBinContent(i + 1, cont / res_q / 100.);
        p_v2_Q_obs_new->SetBinError(i + 1, err / res_q / 100.);

        // if (cen == 8 && ep_option == 0 && option112 == 0 && (cont / res_q / 100.) != 0)
        // {
        // cout << "cont = " << cont << endl;
        // cout << "err = " << err << endl;
        // cout << "res_q = " << res_q << endl;
        // cout << "cont / res_q / 100. = " << cont / res_q / 100. << endl;
        // cout << "err / res_q / 100. = " << err / res_q / 100. << endl;
        // }
    }

    //////// Obtaining delta gamma plots:
    TProfile2D *Parity_Q_obs1 = (TProfile2D *)file->Get("pParity_e_Q2parent_obs1");
    TProfile2D *Parity_Q_obs2 = (TProfile2D *)file->Get("pParity_e_Q2parent_obs2");
    TProfile2D *Parity_w_Q_obs1 = (TProfile2D *)file->Get("pParity_w_Q2parent_obs1");
    TProfile2D *Parity_w_Q_obs2 = (TProfile2D *)file->Get("pParity_w_Q2parent_obs2");
    Parity_Q_obs1->Scale(1);
    Parity_Q_obs1->Add(Parity_w_Q_obs1);
    Parity_Q_obs2->Scale(1);
    Parity_Q_obs2->Add(Parity_w_Q_obs2);

    // cout << "sig loaded" << endl;

    if (!pion_results)
    {
        int index_for_data[12] = {3, 4, 7, 8, 15, 16, 19, 20, 27, 28, 31, 32};
    }
    else if (pion_results)
    {
        int index_for_data[12] = {39, 40, 43, 44, 51, 52, 55, 56, 63, 64, 67, 68};
    }
    TH1 *Parity_Q_ss1 = (TH1 *)Parity_Q_obs1->ProjectionY("Parity_Q_ss1", index_for_data[0 + 2 * option112 + 4 * ep_option], index_for_data[0 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_os1 = (TH1 *)Parity_Q_obs1->ProjectionY("Parity_Q_os1", index_for_data[1 + 2 * option112 + 4 * ep_option], index_for_data[1 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_ss2 = (TH1 *)Parity_Q_obs2->ProjectionY("Parity_Q_ss2", index_for_data[0 + 2 * option112 + 4 * ep_option], index_for_data[0 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_os2 = (TH1 *)Parity_Q_obs2->ProjectionY("Parity_Q_os2", index_for_data[1 + 2 * option112 + 4 * ep_option], index_for_data[1 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_ss = (TH1 *)Parity_Q_ss1->Clone();
    TH1 *Parity_Q_os = (TH1 *)Parity_Q_os1->Clone();

    // background plots
    TProfile2D *Parity_Q_obs1_rot = (TProfile2D *)file->Get("pParity_e_Q2parent_obs1_rot");
    TProfile2D *Parity_Q_obs2_rot = (TProfile2D *)file->Get("pParity_e_Q2parent_obs2_rot");
    TProfile2D *Parity_w_Q_obs1_rot = (TProfile2D *)file->Get("pParity_w_Q2parent_obs1_rot");
    TProfile2D *Parity_w_Q_obs2_rot = (TProfile2D *)file->Get("pParity_w_Q2parent_obs2_rot");
    Parity_Q_obs1_rot->Scale(1);
    Parity_Q_obs1_rot->Add(Parity_w_Q_obs1_rot);
    Parity_Q_obs2_rot->Scale(1);
    Parity_Q_obs2_rot->Add(Parity_w_Q_obs2_rot);

    cout << "bkg loaded" << endl;

    TH1 *Parity_Q_ss1_rot = (TH1 *)Parity_Q_obs1_rot->ProjectionY("Parity_Q_ss1_rot", index_for_data[0 + 2 * option112 + 4 * ep_option], index_for_data[0 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_os1_rot = (TH1 *)Parity_Q_obs1_rot->ProjectionY("Parity_Q_os1_rot", index_for_data[1 + 2 * option112 + 4 * ep_option], index_for_data[1 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_ss2_rot = (TH1 *)Parity_Q_obs2_rot->ProjectionY("Parity_Q_ss2_rot", index_for_data[0 + 2 * option112 + 4 * ep_option], index_for_data[0 + 2 * option112 + 4 * ep_option]);
    TH1 *Parity_Q_os2_rot = (TH1 *)Parity_Q_obs2_rot->ProjectionY("Parity_Q_os2_rot", index_for_data[1 + 2 * option112 + 4 * ep_option], index_for_data[1 + 2 * option112 + 4 * ep_option]);
    // end of background plots

    float last_nonzero_purity = 0;

    for (int i = 0; i < Nbin; i++)
    {
        float res_q = 1;

        // cout << "Parity_Q_ss2->GetBinCenter(i + 1) = " << Parity_Q_ss2->GetBinCenter(i + 1) << endl;
        // cout << "p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i + 1)) = " << p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i + 1)) << endl;
        // cout << "Q2_range[ep_option][cen] = " << Q2_range[ep_option][cen] << endl;

        if (Parity_Q_ss2->GetBinCenter(i + 1) < Q2_range[ep_option][cen])
        {
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i + 1)));

            int i_tmp = i + 2;
            int i_tmp2 = i;

            while (res_q <= 0)
            {
                res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i_tmp)));
                i_tmp++;

                if ((res_q <= 0) && (i_tmp2 > 0)){
                    res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Parity_Q_ss2->GetBinCenter(i_tmp2)));
                    i_tmp2--;
                }

                if (i_tmp > (i+7)) break;
            }
        }
        else
        {
            res_q = p_res->GetBinContent(p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]));

            int i_tmp = p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]) + 1;
            int i_tmp2 = p_res->GetXaxis()->FindBin(Q2_range[ep_option][cen]) - 1;

            while (res_q <= 0)
            {
                res_q = p_res->GetBinContent(i_tmp);
                i_tmp++;

                if ((res_q <= 0) && (i_tmp2 > 0)){
                    res_q = p_res->GetBinContent(i_tmp2);
                    i_tmp2--;
                }

                if (i_tmp > (i+7)) break;
            }

            // if(ep_option == 1) cout << "p_res->GetBinContent(41) = " << p_res->GetBinContent(41) << endl;
        }

        if (ep_option <= 1)
            res_q = sqrt(res_q);

        float purity = 0;

        if (i < 100)
            purity = (lam_purity_Q2[ep_option][cen][i] + antilam_purity_Q2[ep_option][cen][i]) / 2.0;

        if (purity != 0)
        {
            last_nonzero_purity = purity;
        }
        else
        {
            purity = last_nonzero_purity;
        }

        // purity = 1;

        float cont = (Parity_Q_ss2->GetBinContent(i + 1) - (1.0 - purity) * Parity_Q_ss2_rot->GetBinContent(i + 1)) / purity;
        float err = error_add(Parity_Q_ss1->GetBinError(i + 1), (1.0 - purity) * Parity_Q_ss1_rot->GetBinError(i + 1)) / purity;
        Parity_Q_ss->SetBinContent(i + 1, cont / res_q / 100.);
        Parity_Q_ss->SetBinError(i + 1, err / res_q / 100.);

        cont = (Parity_Q_os2->GetBinContent(i + 1) - (1.0 - purity) * Parity_Q_os2_rot->GetBinContent(i + 1)) / purity;
        err = error_add(Parity_Q_os1->GetBinError(i + 1), (1.0 - purity) * Parity_Q_os1_rot->GetBinError(i + 1)) / purity;
        Parity_Q_os->SetBinContent(i + 1, cont / res_q / 100.);
        Parity_Q_os->SetBinError(i + 1, err / res_q / 100.);

        // if (cen == 8 && ep_option == 0 && (cont / res_q / 100.) != 0)
        // {
        // cout << "cont = " << cont << endl;
        // cout << "err = " << err << endl;
        // cout << "res_q = " << res_q << endl;
        // cout << "cont / res_q / 100. = " << cont / res_q / 100. << endl;
        // cout << "err / res_q / 100. = " << err / res_q / 100. << endl;
        // }
    }

    //////// Rebinning Data:
    const int N_bins = Q2_range[ep_option][cen] * 10 / rebin[ep_option][cen];

    float v2q[N_bins], v2q_err[N_bins], d_gq[N_bins], d_gq_err[N_bins];
    float new_x_numerator[N_bins] = {0.}, new_x_denomenator[N_bins] = {0.}, new_y_numerator[N_bins] = {0.}, new_y_denomenator[N_bins] = {0.};

    for (int i = 0; i < N_bins; i++)
    {

        for (int j = 0; j < rebin[ep_option][cen]; j++)
        {
            if (p_v2_Q_obs_new->GetBinError(rebin[ep_option][cen] * i + 1 + j) == 0 || Parity_Q_os->GetBinError(rebin[ep_option][cen] * i + 1 + j) == 0)
                continue;
            new_x_numerator[i] += p_v2_Q_obs_new->GetBinContent(rebin[ep_option][cen] * i + 1 + j) / pow(p_v2_Q_obs_new->GetBinError(rebin[ep_option][cen] * i + 1 + j), 2);
            new_x_denomenator[i] += 1.0 / pow(p_v2_Q_obs_new->GetBinError(rebin[ep_option][cen] * i + 1 + j), 2);
            new_y_numerator[i] += Parity_Q_os->GetBinContent(rebin[ep_option][cen] * i + 1 + j) / pow(Parity_Q_os->GetBinError(rebin[ep_option][cen] * i + 1 + j), 2);
            new_y_denomenator[i] += 1.0 / pow(Parity_Q_os->GetBinError(rebin[ep_option][cen] * i + 1 + j), 2);
        }

        v2q[i] = new_x_numerator[i] / new_x_denomenator[i];
        v2q_err[i] = sqrt(1.0 / new_x_denomenator[i]);
        d_gq[i] = new_y_numerator[i] / new_y_denomenator[i];
        d_gq_err[i] = sqrt(1.0 / new_y_denomenator[i]);
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
    p_res->SetMarkerStyle(1);
    p_res->GetXaxis()->SetRangeUser(0, 50);
    p_res->GetXaxis()->SetTitle("Q_{2}^{2}");
    p_res->GetXaxis()->CenterTitle();
    p_res->GetYaxis()->SetTitle("Resolution");
    p_res->GetYaxis()->SetTitleOffset(1.1);
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

void old_ptsplit_method(int cen)
{
    // profile_divide_by_reso_corr_by_bkg("Parity_int_obs", cen, false);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_ss_obs", cen, false);
    // profile_divide_by_reso_corr_by_bkg("Delta_int_ss_obs", cen, false);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_obs", cen, true);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_ss_obs", cen, true);
    // profile_divide_by_reso_corr_by_bkg("Delta_int_ss_obs", cen, true);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_obs3_splitpt_QQcut", cen);
    // profile_divide_by_reso_corr_by_bkg("Parity_int_ss_obs3_splitpt_QQcut", cen);
    // profile_divide_by_reso_corr_by_bkg("Delta_int_ss_obs3_splitpt_QQcut", cen);
}

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