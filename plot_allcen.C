#include <TProfile>
#include <TCanvas>
// #include <TString>
#include <fstream>
#include <iostream>
#include <TFile>
#include <TF1>

using namespace std;

TFile *output_File;

double npart[9];
double npart_err[9];
float cen9_eff[9] = {75, 65, 55, 45, 35, 25, 15, 7.5, 2.5};
float cen9_err_eff[9] = {0.};

const char* method_name = "KFParticle";

void read_npart();

void plot_TH1D(const char* name);
void plot_TH2D(const char* name);
void plot_TProfile(const char* name);

void combine_processed_results(const char* name);
void combine_proton_tof_eff(const char* name);
void combine_proton_tof_eff_2D(const char* name);
TGraphErrors combine_cen_into_one_graph(const char* name, bool save, int bin);
void plot_two_TGraphErrors(TGraphErrors graph1, TGraphErrors graph2, const char* name);
TGraphErrors subtract_two_graphs(TGraphErrors graph1, TGraphErrors graph2, bool save, const char* name);
void plot_TH1D_plam(const char* name);
void plot_TH2D_plam(const char* name);
void plot_TH1D_systematics(const char* name, float cut_value, const char* x_name);

int max_cen = 8;
int min_cen = 0;

void plot_allcen(){
    output_File = new TFile(TString::Format("./%s_Results/output_allcen.root", method_name).Data(), "recreate");

    read_npart();

    // plot_TH1D("proton_overlap_ratio");
    plot_TH1D("num_lam_tree"); // num of lam in analysis
    plot_TH1D("num_antilam_tree"); // num of antilam in analysis
    plot_TH1D("num_lam_rot_tree"); // num of bkg lam in analysis
    plot_TH1D("num_antilam_rot_tree"); // num of bkg antilam in analysis
    plot_TH1D("num_proton_tree");
    plot_TH1D("num_antiproton_tree");
    // plot_TH1D("num_lam_final");
    // plot_TH1D("num_proton_final");
    // plot_TH1D("num_proton_used");

    /************** Systematic Error Cuts **************/
    plot_TH1D_systematics("nHitsFitQA", 20, "nHitsFit");
    plot_TH1D_systematics("nSigma_dauQA", 2, "nSigma - Daughters of Reconstructed Baryons");
    plot_TH1D_systematics("dca_protonQA", 1, "DCA of Daughter Protons");
    plot_TH1D_systematics("dca_pionQA", 2, "DCA of Daughter Pions");
    // plot_TH1D_systematics("dca_LambdaQA");
    plot_TH1D_systematics("Hist_Pt2", 0.5, "Transverse Momentum of Primary Protons"); // Pt of Protons + Antiprotons (multiplied by number of lam+antilam)
    plot_TH1D_systematics("nSigma_prim_protonQA", 1, "nSigma - Primary Protons");
    /**************************************************/

    plot_TH1D("Hist_Pt"); // Pt of Lambdas + AntiLambdas
    plot_TH1D("Hist_Pt_rot"); // Pt of Bkg Lambdas + AntiLambdas
    // plot_TH1D("Hist_Pt2_rot");
    plot_TH1D("hDpt"); // Abs difference between Pt of (lam/antilam) and (p/antip)
    plot_TH1D("hDpt_rot"); // Abs difference between Pt of (bkg lam/antilam) and (p/antip)

    plot_TH1D("Hist_Q2_TPC");
    plot_TH1D("Hist_Q2_EPD");
    plot_TH1D("Hist_Q2_EPD1");
    plot_TH1D("Hist_Q2_TPC_pion");
    plot_TH1D("Hist_Q2_EPD_pion");
    plot_TH1D("Hist_Q2_EPD1_pion");
    plot_TH1D("Hist_lam_phi_before");
    plot_TH1D("Hist_lam_phi_after");
    plot_TH1D("Hist_p_phi_before");
    plot_TH1D("Hist_p_phi_after");
    // plot_TH1D("Hist_Q2_parent");
    // plot_TH1D("Hist_Q2_parent_remove1");
    // plot_TH1D("Hist_Q2_parent_QQcut");
    // plot_TH1D("num_gamma_final");
    // plot_TH1D_plam("SelectRefMultM0");
    // plot_TH1D("Hist_parent_phi_low_pT");
    // plot_TH1D("Hist_parent_phi_high_pT");

    plot_TProfile("v2_2_pion");

    // plot_TH1D("Hist_Q2_parent_tp");
    // plot_TH1D("Hist_Q2_parent_p");
    // plot_TH1D("Hist_Q2_parent_ap");
    // plot_TH1D("Hist_Q2_parent_lam");
    // plot_TH1D("Hist_ngamma_wrange");
    // plot_TH1D("Hist_nlam_wrange");
    // plot_TH1D("Hist_np_wrange");
    // plot_TH1D("Hist_nap_wrange");
    // plot_TH1D("Hist_ntp_wrange");
    // plot_TH2D("Hist_check_trksplitting_flowvector_ss");
    // plot_TH2D("Hist_check_trksplitting_flowvector_os");
    // plot_TH1D("Hist_check_trksplitting_phi_ss");
    // plot_TH1D("Hist_check_trksplitting_phi_os");
    // plot_TH1D("Hist_check_trksplitting_phi_all_ss");
    // plot_TH1D("Hist_check_trksplitting_phi_all_os");
    // plot_TH1D("Hist_Q2_parent_diff");
    // plot_TH1D("Hist_check_trksplitting_phi_p");
    // plot_TH1D("Hist_check_trksplitting_phi_ap");
    // plot_TH1D("Hist_check_trksplitting_mag_ss");
    // plot_TH1D("Hist_check_trksplitting_mag_os");
    // plot_TH1D("Hist_check_trksplitting_mag_peak_ss");
    // plot_TH1D("Hist_check_trksplitting_mag_peak_os");
    
    // plot_TH2D("Hist_check_trksplitting_pmom_ss");
    // plot_TH2D("Hist_check_trksplitting_pmom_os");
    // plot_TH1D("Hist_check_trksplitting_ptrkid_ss");
    // plot_TH1D("Hist_check_trksplitting_ptrkid_os");
    // plot_TH1D("Hist_check_trksplitting_ptrkid2_ss");
    // plot_TH1D("Hist_check_trksplitting_ptrkid2_os");

    // plot_TH1D("Hist_phi_below1_ss");
    // plot_TH1D("Hist_phi_after1_ss");
    // plot_TH1D("Hist_phi_both1_ss");
    // plot_TH1D("Hist_phi_below1_os");
    // plot_TH1D("Hist_phi_after1_os");
    // plot_TH1D("Hist_phi_both1_os");
    // plot_TH1D("Hist_Q2_parent_pair_ss");
    // plot_TH1D("Hist_Q2_parent_pair_os");
    // plot_TH1D("nHits_ratio_under_peak");

    plot_TH2D("Ref_TOF");
    plot_TH2D("V0Mass_0");
    plot_TH2D("V0Mass_1");
    plot_TH2D("V0Mass_2");
    plot_TH2D("V0Mass_0_anti");
    plot_TH2D("V0Mass_1_anti");
    plot_TH2D("V0Mass_2_anti");
    plot_TH2D("V0Mass_0_pion");
    plot_TH2D("V0Mass_1_pion");
    plot_TH2D("V0Mass_2_pion");
    plot_TH2D("V0Mass_0_pion_anti");
    plot_TH2D("V0Mass_1_pion_anti");
    plot_TH2D("V0Mass_2_pion_anti");
    plot_TH2D("Hist_EPD_EP_east_flat");
    plot_TH2D("Hist_EPD_EP_west_flat");
    plot_TH2D("Hist_EPD_EP1_east_flat");
    plot_TH2D("Hist_EPD_EP1_west_flat");
    // plot_TH2D("EtaPtDist");
    // plot_TH2D("EtaPtDist_anti");
    // plot_TH2D("Hist_RefMult_Q2_parent");
    // plot_TH2D("Hist_RefMult_Q2");

    // plot_TProfile("Hist_Q2_parent_vs_Qcount_parent");

    // plot_TH1D_plam("proton_mom_diff");
    // plot_TH1D_plam("percentLeft");

    /* plot_TH1D_plam("RefMult"); // number events after all event cuts
    plot_TH1D_plam("RefMultA"); // number events in particular centrality
    plot_TH2D_plam("VertexXY_kf");
    plot_TH1D_plam("centrality");
    plot_TH1D_plam("lam_dist");
    plot_TH1D_plam("antilam_dist");*/
    plot_TH1D_plam("VertexZ_kf");
    /*plot_TH2D_plam("Ref_TOF_before");
    plot_TH2D_plam("Ref_TOF_after");*/

    // combine_processed_results("Hist_v2parent_eta_obs5");
    // combine_processed_results("Hist_v2parent_eta_obs5_rot");
    // combine_processed_results("Hist_v2parent_pt_obs5");
    // combine_processed_results("Hist_v2parent_pt_obs5_rot");

    // combine_proton_tof_eff("TOF_eff");
    // combine_proton_tof_eff_2D("m2_proton_vs_pT");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1", false, 4), combine_cen_into_one_graph("Parity_int_obs1", false, 3), "Gamm132 w/o Efficiency Corrections - os vs. ss [Lambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3", false, 4), combine_cen_into_one_graph("Parity_int_obs3", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss [Lambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1", false, 8), combine_cen_into_one_graph("Parity_int_obs1", false, 7), "Gamm132 flipped w/o Efficiency Corrections - os vs. ss [Lambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3", false, 8), combine_cen_into_one_graph("Parity_int_obs3", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss [Lambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs1", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1", false, 3), "Gamm112 w/o Efficiency Corrections - os vs. ss [Lambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss [Lambda]");
    // // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs1", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1", false, 3), "Delta w/o Efficiency Corrections - os vs. ss [Lambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3", false, 3), "Delta w Efficiency Corrections - os vs. ss [Lambda]");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs1_anti", false, 3), "Gamm132 w/o Efficiency Corrections - os vs. ss [AntiLambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs3_anti", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss [AntiLambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs1_anti", false, 7), "Gamm132 flipped w/o Efficiency Corrections - os vs. ss [AntiLambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs3_anti", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss [AntiLambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs1_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_anti", false, 3), "Gamm112 w/o Efficiency Corrections - os vs. ss [AntiLambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_anti", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss [AntiLambda]");
    // // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs1_anti", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1_anti", false, 3), "Delta w/o Efficiency Corrections - os vs. ss [AntiLambda]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_anti", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_anti", false, 3), "Delta w Efficiency Corrections - os vs. ss [AntiLambda]");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs1_combined", false, 3), "Gamm132 w/o Efficiency Corrections - os vs. ss [Combined]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs3_combined", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss [Combined]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs1_combined", false, 7), "Gamm132 flipped w/o Efficiency Corrections - os vs. ss [Combined]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs3_combined", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss [Combined]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs1_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_combined", false, 3), "Gamm112 w/o Efficiency Corrections - os vs. ss [Combined]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_combined", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss [Combined]");
    // // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs1_combined", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1_combined", false, 3), "Delta w/o Efficiency Corrections - os vs. ss [Combined]");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_combined", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_combined", false, 3), "Delta w Efficiency Corrections - os vs. ss [Combined]");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 3), "Gamm132 w/o Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 8), combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 7), "Gamm132 flipped w/o Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 8), combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut", false, 3), "Gamm112 w/o Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs1_QQcut", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1_QQcut", false, 3), "Delta w/o Efficiency Corrections - os vs. ss [Lambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_QQcut", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_QQcut", false, 3), "Delta w Efficiency Corrections - os vs. ss [Lambda] - QQcut");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 3), "Gamm132 w/o Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 7), "Gamm132 flipped w/o Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_anti", false, 3), "Gamm112 w/o Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_anti", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs1_QQcut_anti", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1_QQcut_anti", false, 3), "Delta w/o Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_QQcut_anti", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_QQcut_anti", false, 3), "Delta w Efficiency Corrections - os vs. ss [AntiLambda] - QQcut");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 3), "Gamm132 w/o Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 7), "Gamm132 flipped w/o Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_combined", false, 3), "Gamm112 w/o Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_combined", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs1_QQcut_combined", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1_QQcut_combined", false, 3), "Delta w/o Efficiency Corrections - os vs. ss [Combined] - QQcut");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_QQcut_combined", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_QQcut_combined", false, 3), "Delta w Efficiency Corrections - os vs. ss [Combined] - QQcut");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_bkgcorr", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss Lambda (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_bkgcorr", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss Lambda (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_bkgcorr", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss Lambda (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_bkgcorr", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_bkgcorr", false, 3), "Delta w Efficiency Corrections - os vs. ss Lambda (BKG Corrected)");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_anti_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_anti_bkgcorr", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss AntiLambda (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_anti_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_anti_bkgcorr", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss AntiLambda (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_anti_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_anti_bkgcorr", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss AntiLambda (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_anti_bkgcorr", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_anti_bkgcorr", false, 3), "Delta w Efficiency Corrections - os vs. ss AntiLambda (BKG Corrected)");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss Combined (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss Combined (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_combined_bkgcorr", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss Combined (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_combined_bkgcorr", false, 3), "Delta w Efficiency Corrections - os vs. ss Combined (BKG Corrected)");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_bkgcorr", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss Lambda [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_bkgcorr", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss Lambda [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_bkgcorr", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss Lambda [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_QQcut_bkgcorr", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_QQcut_bkgcorr", false, 3), "Delta w Efficiency Corrections - os vs. ss Lambda [QQ Cut Included] (BKG Corrected)");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_anti_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_anti_bkgcorr", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss AntiLambda [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_anti_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_anti_bkgcorr", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss AntiLambda [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_anti_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_anti_bkgcorr", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss AntiLambda [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_QQcut_anti_bkgcorr", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_QQcut_anti_bkgcorr", false, 3), "Delta w Efficiency Corrections - os vs. ss AntiLambda [QQ Cut Included] (BKG Corrected)");

    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 3), "Gamm132 w Efficiency Corrections - os vs. ss Combined [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 7), "Gamm132 flipped w Efficiency Corrections - os vs. ss Combined [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 3), "Gamm112 w Efficiency Corrections - os vs. ss Combined [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 3), "Delta w Efficiency Corrections - os vs. ss Combined [QQ Cut Included] (BKG Corrected)");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1", false, 4), combine_cen_into_one_graph("Parity_int_obs1", false, 3), false, "Gamma132 w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w/o Efficiency Corrections (Lambda)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3", false, 4), combine_cen_into_one_graph("Parity_int_obs3", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (Lambda)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1", false, 8), combine_cen_into_one_graph("Parity_int_obs1", false, 7), false, "Gamma132 flipped w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w/o Efficiency Corrections (Lambda)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3", false, 8), combine_cen_into_one_graph("Parity_int_obs3", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (Lambda)");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 3), false, "Gamma132 w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w/o Efficiency Corrections (Lambda) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (Lambda) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 8), combine_cen_into_one_graph("Parity_int_obs1_QQcut", false, 7), false, "Gamma132 flipped w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w/o Efficiency Corrections (Lambda) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 8), combine_cen_into_one_graph("Parity_int_obs3_QQcut", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (Lambda) [QQ cut]");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_anti", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs1_anti", false, 3), false, "Gamma132 w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w/o Efficiency Corrections (AntiLambda)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_anti", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs3_anti", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (AntiLambda)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_anti", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs1_anti", false, 7), false, "Gamma132 flipped w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w/o Efficiency Corrections (AntiLambda)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_anti", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs3_anti", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (AntiLambda)");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_anti", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 3), false, "Gamma132 w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w/o Efficiency Corrections (AntiLambda) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_anti", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (AntiLambda) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_anti", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs1_QQcut_anti", false, 7), false, "Gamma132 flipped w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w/o Efficiency Corrections (AntiLambda) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_anti", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_anti", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 8), combine_cen_into_one_graph("Parity_int_obs3_QQcut_anti", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (AntiLambda) [QQ cut]");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_combined", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs1_combined", false, 3), false, "Gamma132 w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w/o Efficiency Corrections (Combined)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_combined", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs3_combined", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (Combined)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_combined", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs1_combined", false, 7), false, "Gamma132 flipped w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w/o Efficiency Corrections (Combined)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_combined", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs3_combined", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (Combined)");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_combined", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 3), false, "Gamma132 w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w/o Efficiency Corrections (Combined) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_combined", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (Combined) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs1_QQcut_combined", false, 3), false, "Gamma112 w/o Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs1_QQcut_combined", false, 7), false, "Gamma132 flipped w/o Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w/o Efficiency Corrections (Combined) [QQ cut]");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_combined", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_QQcut_combined", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 8), combine_cen_into_one_graph("Parity_int_obs3_QQcut_combined", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (Combined) [QQ cut]");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_combined_bkgcorr", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections (BKG Corrected)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_combined_bkgcorr", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_combined_bkgcorr", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections (BKG Corrected)");

    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 3), false, "Gamma132 w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 w Efficiency Corrections [QQ Cut Included] (BKG Corrected)");
    // plot_two_TGraphErrors(subtract_two_graphs(combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 4), combine_cen_into_one_graph("Parity_int_ss_obs3_splitpt_QQcut_combined_bkgcorr", false, 3), false, "Gamma112 w Efficiency Correstions os-ss"), subtract_two_graphs(combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 8), combine_cen_into_one_graph("Parity_int_obs3_splitpt_QQcut_combined_bkgcorr", false, 7), false, "Gamma132 flipped w Efficiency Correstions os-ss"), "Delta Gamma 112 vs. Delta Gamma132 flipped w Efficiency Corrections [QQ Cut Included] (BKG Corrected)");

    // subtract_two_graphs(combine_cen_into_one_graph("Delta_int_ss_obs1", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs1", false, 3), true, "Delta Delta w/o Efficiency Corrections");
    // subtract_two_graphs(combine_cen_into_one_graph("Delta_int_ss_obs3", false, 4), combine_cen_into_one_graph("Delta_int_ss_obs3", false, 3), true, "Delta Delta w Efficiency Corrections");

    output_File->Close();
}

void plot_TH1D(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    for(int i = min_cen; i <= max_cen; i++){
        cout << "i = " << i << endl;
        TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", method_name, i).Data());
        TH1D *temp_hist = (TH1D*) file->Get(name);
        c1.cd(i+1);
        temp_hist->GetYaxis()->SetTitle("Number of Events");
        temp_hist->GetYaxis()->SetTitleOffset(1.3);
        temp_hist->Draw();
    }

    output_File->cd();
    c1.Write(name);
}

void plot_TH1D_systematics(const char* name, float cut_value, const char* x_name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    for(int i = min_cen; i <= max_cen; i++){
        cout << "i = " << i << endl;
        TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", method_name, i).Data());
        TH1D *temp_hist = (TH1D*) file->Get(name);
        c1.cd(i+1);
        temp_hist->SetTitle(TString::Format("Centrality Bin %d", i+1));
        temp_hist->GetYaxis()->SetTitle("Number of Events");
        temp_hist->GetYaxis()->SetTitleOffset(1.3);
        temp_hist->GetXaxis()->SetTitle(x_name);
        temp_hist->Draw();
        TLine *l1 = new TLine(cut_value, 0, cut_value, temp_hist->GetMaximum());
        l1->SetLineColor(kRed);
        l1->Draw();
    }

    output_File->cd();
    c1.Write(name);
}

void plot_TProfile(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    for(int i = min_cen; i <= max_cen; i++){
        cout << "i = " << i << endl;
        TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", method_name, i).Data());
        TProfile *temp_hist = (TProfile*) file->Get(name);
        c1.cd(i+1);
        temp_hist->GetYaxis()->SetTitle("Qcount_Parent");
        temp_hist->GetYaxis()->SetTitleOffset(1.3);
        temp_hist->GetXaxis()->SetTitle("Q2_Parent");
        temp_hist->Draw();
    }

    output_File->cd();
    c1.Write(name);
}

void plot_TH1D_plam(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    for(int i = min_cen; i <= max_cen; i++){
        cout << "i = " << i << endl;
        // TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_plam.root", method_name, i).Data());
        // TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_output.root", method_name, i).Data());
        TFile *file = new TFile(TString::Format("./%s/plam_27.root", method_name, i).Data());
        TH1D *temp_hist = (TH1D*) file->Get(name);
        c1.cd(i+1);
        // temp_hist->GetYaxis()->SetTitle("Number of p/#bar{p} Pairs");
        // temp_hist->GetYaxis()->SetTitleOffset(1.3);
        // temp_hist->GetXaxis()->SetTitle("p_{T} Difference (GeV/c)");
        temp_hist->SetTitle(TString::Format("Centrality Bin %d", i+1));
        temp_hist->GetYaxis()->SetTitle("Number of Events");
        temp_hist->GetYaxis()->SetTitleOffset(1.3);
        temp_hist->GetXaxis()->SetTitle("z-component of Primary Vertex (cm)");
        temp_hist->Draw();
    }

    output_File->cd();
    c1.Write(name);
}

// void plot_TH2D_plam(const char* name){

//     cout << name << endl;

//     TCanvas c1("c1", name, 1500, 800);
//     c1.Divide(3, 3);

//     for(int i = min_cen; i <= max_cen; i++){
//         cout << "i = " << i << endl;
//         // TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_plam.root", method_name, i).Data());
//         TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_output.root", method_name, i).Data());
//         TH2D *temp_hist = (TH2D*) file->Get(name);
//         c1.cd(i+1);
//         // temp_hist->GetYaxis()->SetTitle("Number of p/#bar{p} Pairs");
//         // temp_hist->GetYaxis()->SetTitleOffset(1.3);
//         // temp_hist->GetXaxis()->SetTitle("p_{T} Difference (GeV/c)");
//         temp_hist->Draw();
//     }

//     output_File->cd();
//     c1.Write(name);
// }

void plot_TH2D(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    TCanvas c2("c2", TString::Format("%s_px", name).Data(), 1500, 800);
    c2.Divide(3, 3);

    TCanvas c3("c3", TString::Format("%s_py", name).Data(), 1500, 800);
    c3.Divide(3, 3);

    TH2D *temp_hist;
    TH1D *temp_hist_projection;

    for(int i = min_cen; i <= max_cen; i++){
        TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d.gamma112_fullEP_eff_pT02_module.root", method_name, i).Data());
        temp_hist = (TH2D*) file->Get(name);
        c1.cd(i+1);
        gPad->SetLogz();
        (TH2D *)(temp_hist->Clone())->Draw("COLZ");

        temp_hist_projection = (TH1D*)((TH2D *)temp_hist->Clone())->ProjectionX()->Clone();
        c2.cd(i+1);
        temp_hist_projection->GetYaxis()->SetRangeUser(0, temp_hist_projection->GetBinContent(temp_hist_projection->GetMaximumBin()) * 1.1);
        (TH1D*)(temp_hist_projection->Clone())->Draw();

        temp_hist_projection_y = (TH1D*)((TH2D *)temp_hist->Clone())->ProjectionY()->Clone();
        c3.cd(i+1);
        (TH1D*)(temp_hist_projection_y->Clone())->Draw();

        delete temp_hist;
        temp_hist = NULL;
        delete temp_hist_projection;
        temp_hist_projection = NULL;
        delete temp_hist_projection_y;
        temp_hist_projection_y = NULL;

        // file.Close();
    }

    output_File->cd();
    c1.Write(name);
    c2.Write(TString::Format("%s_px", name));
    c3.Write(TString::Format("%s_py", name));
}

void plot_TH2D_plam(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    TCanvas c2("c2", TString::Format("%s_px", name).Data(), 1500, 800);
    c2.Divide(3, 3);

    TH2D *temp_hist;
    TH1D *temp_hist_projection;

    for(int i = min_cen; i <= max_cen; i++){
        // TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_plam.root", method_name, i).Data());
        TFile *file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_output.root", method_name, i).Data());
        temp_hist = (TH2D*) file->Get(name);
        c1.cd(i+1);
        gPad->SetLogz();
        (TH2D *)(temp_hist->Clone())->Draw("COLZ");

        temp_hist_projection = (TH1D*)((TH2D *)temp_hist->Clone())->ProjectionY()->Clone();
        c2.cd(i+1);
        (TH1D*)(temp_hist_projection->Clone())->Draw();

        delete temp_hist;
        temp_hist = NULL;
        delete temp_hist_projection;
        temp_hist_projection = NULL;

        // file.Close();
    }

    output_File->cd();
    c1.Write(name);
    c2.Write(TString::Format("%s_px", name));
}

void combine_processed_results(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    for(int i = min_cen; i <= max_cen; i++){
        TFile *file = new TFile(TString::Format("./%s_Results/output_cen%d.root", method_name, i).Data());
        TProfile *temp_hist = (TProfile*) file->Get(name);
        c1.cd(i+1);
        temp_hist->SetMarkerStyle(kFullCircle);
        temp_hist->Draw();
    }

    output_File->cd();
    c1.Write(name);
}

TGraphErrors combine_cen_into_one_graph(const char* name, bool save, int bin){

    cout << name << ", " << bin << endl;

    double temp_x[9] = {0.}, temp_y[9] = {0.}, temp_x_err[9] = {0.}, temp_y_err[9] = {0.};
    
    for(int i = min_cen; i <= max_cen; i++){
        // cout << TString::Format("./%s_Results/output_cen%d.root", method_name, i).Data() << endl;
        TFile *file = new TFile(TString::Format("./%s_Results/output_cen%d.root", method_name, i).Data());
        // cout << name << endl;
        TGraphErrors *temp_hist = (TGraphErrors*) file->Get(name);

        temp_x[i] = i;
        temp_x_err[i] = 0;
        temp_y[i] = temp_hist->GetY()[bin-1] / 100.0 * npart[i];
        temp_y_err[i] = temp_hist->GetErrorY(bin-1) / 100.0 * npart[i];

        // cout << "temp_x[" << i << "] = " << temp_x[i] << endl;
        // cout << "temp_y[" << i << "] = " << temp_y[i]/npart[i] * 100.0 << endl;
    }

    TGraphErrors tempgraph((max_cen-min_cen+1), temp_x, temp_y, temp_x_err, temp_y_err);

    if(save){
        tempgraph.SetMarkerStyle(kFullCircle);
        output_File->cd();
        tempgraph.Write(name);
    }

    return tempgraph;
}

void plot_two_TGraphErrors(TGraphErrors graph1, TGraphErrors graph2, const char* name){
    TCanvas c1("c1", name, 1500, 800);
    c1.cd();
    graph1.SetMarkerStyle(kFullCircle);
    graph1.SetTitle(name);
    graph1.Draw("AP");
    graph2.SetMarkerStyle(kFullCircle);
    graph2.SetMarkerColor(kRed);
    graph2.SetLineColor(kRed);
    graph2.Draw("SAMES P");

    TLine line(0,0,9,0);
    line.Draw();

    output_File->cd();
    c1.Write(name);
}

TGraphErrors subtract_two_graphs(TGraphErrors graph1, TGraphErrors graph2, bool save, const char* name){
    
    double temp_x[9] = {0.}, temp_y[9] = {0.}, temp_x_err[9] = {0.}, temp_y_err[9] = {0.};
    
    for(int i = 0; i < graph1.GetN(); i++){
        temp_x[i] = graph1.GetX()[i];
        temp_x_err[i] = graph1.GetErrorX(i);
        temp_y[i] = graph1.GetY()[i] - graph2.GetY()[i];
        temp_y_err[i] = sqrt(pow(graph1.GetErrorY(i),2) + pow(graph2.GetErrorY(i),2));
    }

    TGraphErrors tempgraph(graph1.GetN(), temp_x, temp_y, temp_x_err, temp_y_err);

    if(save){
        tempgraph.SetMarkerStyle(kFullCircle);
        output_File->cd();
        tempgraph.Write(name);
    }

    return tempgraph;
}

void combine_proton_tof_eff(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    TFile *file = new TFile(TString::Format("./%s_Results/tof_eff.root", method_name).Data());

    for(int i = min_cen; i <= max_cen; i++){
        TH1F *temp_hist = (TH1F*) file->Get(TString::Format("%s_cen%d", name, i));
        c1.cd(i+1);
        temp_hist->Draw();
    }

    output_File->cd();
    c1.Write(name);
}

void combine_proton_tof_eff_2D(const char* name){

    cout << name << endl;

    TCanvas c1("c1", name, 1500, 800);
    c1.Divide(3, 3);

    TFile *file = new TFile(TString::Format("./%s_Results/tof_eff.root", method_name).Data());

    for(int i = min_cen; i <= max_cen; i++){
        TH2D *temp_hist = (TH2D*) file->Get(TString::Format("%s_cen%d", name, i));
        c1.cd(i+1);
        temp_hist->Draw("COLZ");
    }

    output_File->cd();
    c1.Write(name);
}

void read_npart()
{
    std::fstream myfile("./npart.txt", std::ios_base::in);

    float a;
    int count = 0;
    string line_text;
    while (myfile >> a)
    {
        if(count % 2 == 0)
        {
            npart[count / 2] = a;
        }
        else if(count % 2 == 1)
        {
            npart_err[count / 2] = a;
        }

        count++;
    }
}