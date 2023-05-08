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

// Settings for Macro
const char *method_name = "KFParticle";
int min_cen = 0;
int max_cen = 8;

// Global Files
TFile *file, *output_File; // Files by centrality; file --> input; output_File --> output
TFile *output_File2;       // output File for all centralities in range
ofstream result_file, gamma_results;      // text File output for all centralities in range

// Global Variables for Functional Purposes
float reso[4] = {0.}, reso1[3] = {0.}; // reso = Sub EP Reso; reso1 = Full EP Reso
float cen9_eff[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
float cen9_err_eff[9] = {0.};
double npart[9] = {0.}, npart_err[9] = {0.};
double lam_purity[9][18], antilam_purity[9][18], lam_yield[9][18], antilam_yield[9][18];
float lam_purity_Q2[3][9][100], antilam_purity_Q2[3][9][100];
std::vector<TString> name_options3;

// Global Variables for Plotting all Centralities
std::vector<float> parent_v2_vect[3], alltrks_v2_vect[3], pionpion_v2_vect[3]; // Getting average v2
std::vector<float> proton_v2_vect, lam_v2_vect;
float delta[9], delta_err[9];
float reso_TPC[9] = {0.}, reso_EPD[9] = {0.}, reso_EPD1[9] = {0.};
float reso_TPC_err[9] = {0.}, reso_EPD_err[9] = {0.}, reso_EPD1_err[9] = {0.};
float gamma[2][9][3], gamma_err[2][9][3];
float gamma_ensemble[2][9][3], gamma_ensemble_err[2][9][3];
float v2_lam[9][3], v2_p[9][3], v2_lam_err[9][3], v2_p_err[9][3];

// Beginning of Macro, Establishing Analysis
void initialize_macro();
void read_sig();
void read_npart();
void read_purity_Q2();

// Centrality by Centrality Analysis
void downstream_analysis_bycen(int cen);
//
void initialize_each_cen(int cen);
void extract_reso();
double calc_reso(double res);
double chi(double res);
double resEventPlane(double chi);
//
void organize_ensemble_average_results(int cen);
void profile_divide_by_reso(const char *profile_name, int cen, bool eff_corr, bool QQ_cut = false);
//
void obtain_v2_average(int cen);
void profile_divide_by_reso_fit_by_const(const char *profile_name, int cen); // used for getting average v2(eta)
//
void Q2_parent_results(int cen, int ep_option, int option112);
float float finding_res_from_resq(TH1 *profile_need_corr, TProfile *res_prof, int b_num, float q2range);

//
// Old method that is not used:
void old_ptsplit_method(int cen);
void profile_divide_by_reso_corr_by_bkg(const char *profile_name, int cen);

// Wrapping up Analysis to put results from all centralities together
void organize_and_plot_v2_results();
void organize_and_plot_gamma_results();
void plotting_TGraph_Helper(const char *title, const char *x_title, const char *y_title, int npts, float *x, float *x_err, float *y1, float *y1_err, float *y2 = NULL, float *y2_err = NULL, float *y3 = NULL, float *y3_err = NULL);
