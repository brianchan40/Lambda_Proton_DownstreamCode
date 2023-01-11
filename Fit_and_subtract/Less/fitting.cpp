#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TTree.h>
#include <TLine.h>

void fitting(){
	TFile lam_rot("plamrot.root");
	TH1F *V0Mass_lam_rot;
	TFile *output = new TFile("output_fitting.root", "recreate");
	V0Mass_lam_rot = (TH1F *)lam_rot.Get("V0Mass")->Clone();

	//TF1 *final = new TF1("final", "[0]*TMath::DiLog([1]*x + [2])+[3]", 1.075, 1.185);
	TF1 *final = new TF1("final", "[0]*pow(x, 3) + [1]*pow(x, 2) + [2]*x + [3]", 1.075, 1.185);
	//TF1 *final = new TF1("final", "[0]/pow(([1] + [2]*x), [3]) + [4]", 1.075, 1.185);
	Double_t par[5];
	// par[0] = 1000000;
	// par[1] = 10;
	// par[2] = -10;
	// par[3] = -1.25;
	// final->SetParameters(par);
	V0Mass_lam_rot->Fit(final, "R");
	final->GetParameters(&par[0]);
	final->SetParameters(par);
	V0Mass_lam_rot->Fit(final, "R");

	//TCanvas *trial = new TCanvas("trial", "Try", 800, 800);
	//TF1 *fa2 = new TF1("fa2","TMath::DiLog(x)",0,10);
   	//fa2->Draw();
}