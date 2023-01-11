#include <TH1>
#include <TRandom3>

float q2parent(int num);
int count = 0;

void simpleToyMC(){
    std::vector<float> proton1, proton2;
    TFile f("simpleToyMC.root", "recreate");

    TH1D *proton1_phi_dist = new TH1D("proton1_phi_dist", "proton1_phi_dist", 800, -4, 4);
    TH1D *proton2_phi_dist = new TH1D("proton2_phi_dist", "proton2_phi_dist", 800, -4, 4);
    TH1D *proton_phidiff_dist = new TH1D("proton_phidiff_dist", "proton_phidiff_dist", 800, -4, 4);
    TH1D *Q2_parent = new TH1D("Q2_parent", "Q2_parent", 5000, 0, 10);
    TH1D *Q2_parent_simp = new TH1D("Q2_parent_simp", "Q2_parent_simp", 5000, 0, 10);

    for(int i = 0; i < 700000 ; i++){
        // proton1.push_back(rand_gen->Uniform(-TMath::Pi(), TMath::Pi()));
        // proton1_phi_dist->Fill(proton1.at(i));
        // proton2.push_back(rand_gen->Uniform(-TMath::Pi(), TMath::Pi()));
        // proton2_phi_dist->Fill(proton2.at(i));

        // float phi_diff = proton1.at(i) - proton2.at(i);
        // if(phi_diff < -TMath::Pi()) phi_diff += 2*TMath::Pi();
        // else if(phi_diff > TMath::Pi()) phi_diff -= 2*TMath::Pi();
        // proton_phidiff_dist->Fill(phi_diff);

        // float Q2 = (double)(pow((TMath::Cos(proton1.at(i)) + TMath::Cos(proton2.at(i))), 2) + pow((TMath::Sin(proton1.at(i)) + TMath::Sin(proton2.at(i))), 2)) / 2.0;
        // Q2_parent->Fill(Q2);
        // Q2_parent_simp->Fill(1+TMath::Cos(phi_diff));

        Q2_parent->Fill(q2parent(2));
    }
    for(int i = 0; i < 700000 ; i++) Q2_parent->Fill(q2parent(3));
    for(int i = 0; i < 500000 ; i++) Q2_parent->Fill(q2parent(4));
    for(int i = 0; i < 300000 ; i++) Q2_parent->Fill(q2parent(5));
    for(int i = 0; i < 150000 ; i++) Q2_parent->Fill(q2parent(6));
    for(int i = 0; i < 50000 ; i++) Q2_parent->Fill(q2parent(7));

    f.cd();
    proton1_phi_dist->Write();
    proton2_phi_dist->Write();
    proton_phidiff_dist->Write();
    Q2_parent->Write();
    Q2_parent_simp->Write();
    f.Close();
}

float q2parent(int num){

    float mQQx = 0, mQQy = 0;
    TRandom3 *rand_gen;

    for(int i = 0 ; i < num ; i++){
        rand_gen = new TRandom3(count);
        float proton1_phi = (rand_gen->Rndm() - 0.5) * 2 * TMath::Pi();
        proton1_phi_dist->Fill(proton1_phi);

        mQQx += TMath::Cos(proton1_phi);
        mQQy += TMath::Sin(proton1_phi);
        rand_gen->Clear();
        count++;
    }

    return (double)(mQQx * mQQx + mQQy * mQQy)/2.0;
}