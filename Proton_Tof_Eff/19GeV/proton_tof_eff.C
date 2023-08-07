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
    output_File = new TFile(TString::Format("./tof_eff.root", method_name).Data(), "recreate");
    result_file.open(TString::Format("./output.txt", method_name).Data());

    for(int i = 1; i < 10; i++){
        cout << "Centrality " << i << endl;
        proton_tof_eff_bycen(i);
    }

    output_File->Close();
    result_file.close();
}

void proton_tof_eff_bycen(int cen){
    file = new TFile(TString::Format("./TOF_matchEfficeny_cen%d.root", cen).Data());
    
    extract_tof_eff(cen);

    file->Close();
}

void extract_tof_eff(int cen){

    TCanvas *c1 = (TCanvas *) file->Get("c1");
    TProfile *Hist_Pt_TOF = (TProfile *) c1->FindObject("Hist_Pt_TOF");

    TF1 *temp_func = new TF1("temp_func", "[0] + [1]*sqrt(x) + [2]*x + [3]*pow(x,2)", 0.2, 4);
    Hist_Pt_TOF->Fit(temp_func, "R");

    result_file << temp_func->GetParameter(0) << " " << temp_func->GetParameter(1) << " " << temp_func->GetParameter(2) << " " << temp_func->GetParameter(3) << "\n";

    output_File->cd();
    Hist_Pt_TOF->Write(TString::Format("TOF_eff_cen%d", cen-1));
}