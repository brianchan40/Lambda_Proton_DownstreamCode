void test(){

    cout << "no Cut, with Cut, full" << endl;

    for(int i = 0; i < 9; i++){
        TFile f1(TString::Format("./KFParticle/Results_lam_18/test/15_noCut/cen%d.gamma112_fullEP_eff_pT02_module.root", i).Data());
        TFile f2(TString::Format("./KFParticle/Results_lam_18/test/15_wCut/cen%d.gamma112_fullEP_eff_pT02_module.root", i).Data());
        TFile f3(TString::Format("./KFParticle/Results_lam_18/sys_0/cen%d.gamma112_fullEP_eff_pT02_module.root", i).Data());

        cout << ((TProfile*) f1.Get("Hist_cos_EPD"))->GetBinEntries(1) << ", " << ((TProfile*) f2.Get("Hist_cos_EPD"))->GetBinEntries(1) << ", " << ((TProfile*) f3.Get("Hist_cos_EPD"))->GetBinEntries(1) << endl;

        f1.Close();
        f2.Close();
        f3.Close();
    }
}