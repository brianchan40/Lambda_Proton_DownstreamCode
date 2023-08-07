#include <TProfile>
#include <TCanvas>
// #include <TString>
#include <fstream>
#include <iostream>
#include <TF1>

using namespace std;

TFile *file, *output_File;
ofstream result_file;

const char* method_name = "Traditional";

void extract_tof_eff(int cen);
void extract_m2vspT(int cen);
void proton_tof_eff_bycen(int cen);

void proton_tof_eff(){
    output_File = new TFile(TString::Format("./%s_Results/tof_eff.root", method_name).Data(), "recreate");
    result_file.open(TString::Format("./%s_Results/output.txt", method_name).Data());

    for(int i = 0; i < 9; i++){
        cout << "Centrality " << i << endl;
        proton_tof_eff_bycen(i);
    }

    output_File->Close();
    result_file.close();
}

void proton_tof_eff_bycen(int cen){
    file = new TFile(TString::Format("./%s/Results_lam_18/cen%d_plam.root", method_name, cen).Data());
    
    extract_tof_eff(cen);
    extract_m2vspT(cen);

    file->Close();
}

void extract_tof_eff(int cen){
    TH1F *pTDist_b4_TOF = (TH1F*)file->Get("pTDist_b4_TOF");
    TH1F *pTDist_after_TOF = (TH1F*)file->Get("pTDist_after_TOF");

    pTDist_after_TOF->Divide(pTDist_b4_TOF);
    TF1 *temp_func = new TF1("temp_func", "[0] + [1]*sqrt(x) + [2]*x + [3]*pow(x,2)", 0.2, 4);
    pTDist_after_TOF->Fit(temp_func, "R");

    result_file << temp_func->GetParameter(0) << " " << temp_func->GetParameter(1) << " " << temp_func->GetParameter(2) << " " << temp_func->GetParameter(3) << "\n";

    output_File->cd();
    pTDist_after_TOF->Write(TString::Format("TOF_eff_cen%d", cen));
}

void extract_m2vspT(int cen){
    TH2D *m2_proton_vs_pT = (TH2D*)file->Get("m2_proton_vs_pT");

    output_File->cd();
    m2_proton_vs_pT->Write(TString::Format("m2_proton_vs_pT_cen%d", cen));
}