#include<map>
#include<string>

void TOF_eff_calculator(){
    // map<int, string> cen_to_centrality;
    string cen_to_centrality[10];
    cen_to_centrality[0]="80-100";
    cen_to_centrality[1]="70-80";
    cen_to_centrality[2]="60-70";
    cen_to_centrality[3]="50-60";
    cen_to_centrality[4]="40-50";
    cen_to_centrality[5]="30-40";
    cen_to_centrality[6]="20-30";
    cen_to_centrality[7]="10-20";
    cen_to_centrality[8]="5-10";
    cen_to_centrality[9]="0-5";
    char efffile[200];
    char name[200];
    TFile *f;
    TH1* mc;
    TH1* rc;
    TCanvas* canvas=new TCanvas();
    for(int cent=1;cent<10;cent++){
        //cout<<cent<<endl;
        snprintf(efffile,200,"cen%d.v2_pion.root",cent);
        f = new TFile(efffile,"READ");
        mc = (TH1*)f->Get("Hist_Pt");
        rc = (TH1*)f->Get("Hist_Pt_TOF");
        rc->Divide(mc);
        int Nbins=rc->GetNbinsX();
        for(int bin=1;bin<Nbins+1;bin++){
           cout<<"cen="<<cent<<" centrality="<<cen_to_centrality[cent]<< "% pt="<<rc->GetBinCenter(bin)<<" TOF_match_eff="<<rc->GetBinContent(bin)<<" TOF_match_eff_error="<<rc->GetBinError(bin)<<endl;
        }
        snprintf(name,200,"TOF_matchEfficeny_cen%d.root",cent);
        rc->SetTitle(name);
        rc->GetXaxis()->SetTitle("p_{T}");
        rc->GetYaxis()->SetTitle("");
        rc->Draw();
        canvas->SaveAs(name);
    }
}
