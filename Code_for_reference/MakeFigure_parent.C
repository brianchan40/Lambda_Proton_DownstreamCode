static Double_t chi(double res);
static Double_t resEventPlane(double chi);

void MakeFigure_parent(int cen = 1)
{ // 0 means TPC, 1 means EPD outer, 11 means EPD inner

  int Q2_range = 10; // range of fitting in resolution, also range of plotting
  int rebin = 2;     // choose between 1,2,5
  cout << " rebin = " << rebin << endl;

  // resolution
  char printname1[200];
  sprintf(printname1, "cen%d.Resolution.parent.pdf", cen);

  // q2-v2
  double Yrange_v2[9] = {1, 1, 1, 1, 0.1, 1, 1, 1, 1};
  char printname2[200];
  sprintf(printname2, "cen%d.Q2.v2.parent.pdf", cen);

  // raw plot
  double Xrange_v2_raw = 25;
  double Yrange_delGamma0[9] = {1, 0.6, 0.4, 0.3, 0.4, 0.2, 0.2, 0.2, 0.2};
  char printname3[200];
  sprintf(printname3, "cen%d.raw.parent.fitting.pdf", cen);

  // q2-g112-ssos
  double YRangeMax_ew_Q2[9] = {2, 1, 0.8, 0.4, 0.2, 0.1, 0.2, 0.2, 0.2};
  char printname4[200];
  sprintf(printname4, "cen%d.osss.Q2.parent.pdf", cen);

  // DG112-v2
  double Yrange_delGamma[9] = {0.1, 0.05, 0.02, 0.008, 0.004, 0.002, 0.0005, 0.001, 0.001};
  double Xrange_delGamma[9] = {0.3, 0.25, 0.3, 0.22, 0.18, 0.15, 0.05, 0.03, 0.02};
  double Xfitting_delGamma[9] = {0.1, 0.1, 0.3, 0.22, 0.18, 0.15, 0.05, 0.03, 0.02};
  char printname5[200];
  sprintf(printname5, "cen%d.ese.parent.Q2range%d.rebin%d.pdf", cen, Q2_range, rebin);

  int left[9] = {70, 60, 50, 40, 30, 20, 10, 5, 0};
  int right[9] = {80, 70, 60, 50, 40, 30, 20, 10, 5};
  char centrality[200];
  sprintf(centrality, "%d - %d %% Au+Au, 27 GeV", left[cen - 1], right[cen - 1]);

  /////////////////Read in Files
  char fname[200];
  sprintf(fname, "result-parent-TPC-Q2-N+v22N2/cen%d.gamma112_fullEP_eff_pT02_module.root", cen);

  TFile *f = new TFile(fname);

  ////////////////overall Resolution
  TProfile *Hist_cos = (TProfile *)f->Get("Hist_cos_EPD");
  cout << "our" << endl;
  float res = Hist_cos->GetBinContent(1);
  double res_v2 = sqrt(res);
  // if(opt_EP==0) {
  // 	Hist_cos = (TProfile*)f->Get("Hist_cos");
  // 	res = Hist_cos->GetBinContent(1);
  // 	res_v2 = Hist_cos->GetBinContent(2);
  // 	res_v2 = resEventPlane(sqrt(2)*chi(sqrt(res_v2)));
  // }
  res = sqrt(res);
  // if(opt_EP==11) {
  res = Hist_cos->GetBinContent(4);
  res_v2 = res;
  // }
  cout << "ESE sub EP res = " << res << endl;
  cout << "v2 sub EP res = " << res_v2 << endl;

  ////////////////inclusive v2
  TF1 *fun = new TF1("fun", "[0]", 0, 2);
  TProfile *Hist_v2parent_eta_obs5 = (TProfile *)f->Get("Hist_v2parent_eta_obs5");
  Hist_v2parent_eta_obs5->Fit("fun", "0", "", -1, 1);
  float v2obs = fun->GetParameter(0) / 100 / res_v2;
  float v2obs_err = fun->GetParError(0) / 100 / res_v2;

  TProfile *Hist_v2parent_eta_obs6 = (TProfile *)f->Get("Hist_v2parent_eta_obs6");
  Hist_v2parent_eta_obs6->Fit("fun", "0", "", -1, 1);
  float v2obs_ESE = fun->GetParameter(0) / 100 / res;
  float v2obs_ESE_err = fun->GetParError(0) / 100 / res;

  cout << endl;

  cout << "w.r.t. the participant plane" << endl;
  cout << "v2 parent = " << v2obs << " +/- " << v2obs_err << endl;
  cout << "v2 parent ESE = " << v2obs_ESE << " +/- " << v2obs_ESE_err << endl;

  ////////////////Direct plots
  TProfile *p_Parity_v2e_parent_obs1 = (TProfile *)f->Get("p_Parity_v2e_parent_obs1");
  TProfile *p_Parity_v2e_parent_obs2 = (TProfile *)f->Get("p_Parity_v2e_parent_obs2");
  TProfile *p_Parity_v2w_parent_obs1 = (TProfile *)f->Get("p_Parity_v2w_parent_obs1");
  TProfile *p_Parity_v2w_parent_obs2 = (TProfile *)f->Get("p_Parity_v2w_parent_obs2");
  p_Parity_v2e_parent_obs1->Scale(1);
  p_Parity_v2e_parent_obs1->Add(p_Parity_v2w_parent_obs1);
  p_Parity_v2e_parent_obs2->Scale(1);
  p_Parity_v2e_parent_obs2->Add(p_Parity_v2w_parent_obs2);

  TProfile *p_Parity_v2e_parent_obs3 = (TProfile *)f->Get("p_Parity_v2e_parent_obs3");
  TProfile *p_Parity_v2e_parent_obs4 = (TProfile *)f->Get("p_Parity_v2e_parent_obs4");
  TProfile *p_Parity_v2w_parent_obs3 = (TProfile *)f->Get("p_Parity_v2w_parent_obs3");
  TProfile *p_Parity_v2w_parent_obs4 = (TProfile *)f->Get("p_Parity_v2w_parent_obs4");
  p_Parity_v2e_parent_obs3->Scale(1);
  p_Parity_v2e_parent_obs3->Add(p_Parity_v2w_parent_obs3);
  p_Parity_v2e_parent_obs4->Scale(1);
  p_Parity_v2e_parent_obs4->Add(p_Parity_v2w_parent_obs4);

  int Nbin = p_Parity_v2e_parent_obs1->GetNbinsX();
  float v2_parent_raw[Nbin], v2_parent_raw_err[Nbin], DG112_raw[Nbin], DG112_raw_err[Nbin];
  for (int i = 0; i < Nbin; i++)
  {
    v2_parent_raw[i] = p_Parity_v2e_parent_obs2->GetBinCenter(i + 1);
    v2_parent_raw_err[i] = p_Parity_v2e_parent_obs2->GetBinError(i + 1);

    float cont_ss = p_Parity_v2e_parent_obs2->GetBinContent(i + 1);
    float err_ss = p_Parity_v2e_parent_obs1->GetBinError(i + 1);
    float cont_os = p_Parity_v2e_parent_obs4->GetBinContent(i + 1);
    float err_os = p_Parity_v2e_parent_obs3->GetBinError(i + 1);

    DG112_raw[i] = cont_os - cont_ss;
    DG112_raw_err[i] = sqrt(err_ss * err_ss + err_os * err_os);
  }

  ////////////////Q2 parent results
  TProfile *p_res = (TProfile *)f->Get("p_cos_Q2_parent");
  TF1 *f_res = new TF1("f_res", "[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x", 0, Q2_range);
  p_res->Fit("f_res", "0", "", 0, Q2_range);

  TProfile *p_v2_Q_obs1 = (TProfile *)f->Get("p_v2parente_Q2parent_obs1");
  TProfile *p_v2_Q_obs2 = (TProfile *)f->Get("p_v2parente_Q2parent_obs2");
  TProfile *p_v2w_Q_obs1 = (TProfile *)f->Get("p_v2parentw_Q2parent_obs1");
  TProfile *p_v2w_Q_obs2 = (TProfile *)f->Get("p_v2parentw_Q2parent_obs2");
  p_v2_Q_obs1->Scale(1);
  p_v2_Q_obs1->Add(p_v2w_Q_obs1);
  p_v2_Q_obs2->Scale(1);
  p_v2_Q_obs2->Add(p_v2w_Q_obs2);

  TH1 *p_v2_Q_obs_new = (TH1 *)p_v2_Q_obs1->ProjectionX("p_v2_Q_obs_new");
  Nbin = p_v2_Q_obs1->GetNbinsX();
  for (int i = 0; i < Nbin; i++)
  {
    float cont = p_v2_Q_obs2->GetBinContent(i + 1);
    float err = p_v2_Q_obs1->GetBinError(i + 1);
    float res_q = 1;
    if (p_v2_Q_obs2->GetBinCenter(i + 1) < Q2_range)
      res_q = f_res->Eval(p_v2_Q_obs2->GetBinCenter(i + 1));
    else
      res_q = f_res->Eval(Q2_range);
    // if(opt_EP != 11) res_q = sqrt(res_q);
    p_v2_Q_obs_new->SetBinContent(i + 1, cont / res_q / 100.);
    p_v2_Q_obs_new->SetBinError(i + 1, err / res_q / 100.);
  }
  TProfile2D *Parity_Q_obs1 = (TProfile2D *)f->Get("pParity_e_Q2parent_obs1");
  TProfile2D *Parity_Q_obs2 = (TProfile2D *)f->Get("pParity_e_Q2parent_obs2");
  TProfile2D *Parity_w_Q_obs1 = (TProfile2D *)f->Get("pParity_w_Q2parent_obs1");
  TProfile2D *Parity_w_Q_obs2 = (TProfile2D *)f->Get("pParity_w_Q2parent_obs2");
  Parity_Q_obs1->Scale(1);
  Parity_Q_obs1->Add(Parity_w_Q_obs1);
  Parity_Q_obs2->Scale(1);
  Parity_Q_obs2->Add(Parity_w_Q_obs2);
  TH1 *Parity_Q_ss1 = (TH1 *)Parity_Q_obs1->ProjectionY("Parity_Q_ss1", 3, 3);
  TH1 *Parity_Q_os1 = (TH1 *)Parity_Q_obs1->ProjectionY("Parity_Q_os1", 4, 4);
  TH1 *Parity_Q_ss2 = (TH1 *)Parity_Q_obs2->ProjectionY("Parity_Q_ss2", 3, 3);
  TH1 *Parity_Q_os2 = (TH1 *)Parity_Q_obs2->ProjectionY("Parity_Q_os2", 4, 4);
  TH1 *Parity_Q_ss = (TH1 *)Parity_Q_ss1->Clone();
  TH1 *Parity_Q_os = (TH1 *)Parity_Q_os1->Clone();
  for (int i = 0; i < Nbin; i++)
  {
    float res_q = 1;
    if (Parity_Q_ss2->GetBinCenter(i + 1) < Q2_range)
      res_q = f_res->Eval(Parity_Q_ss2->GetBinCenter(i + 1));
    else
      res_q = f_res->Eval(Q2_range);
    // if(opt_EP != 11) res_q = sqrt(res_q);

    float cont = 0.;
    float err = 0.;
    cont = Parity_Q_ss2->GetBinContent(i + 1);
    err = Parity_Q_ss1->GetBinError(i + 1);
    Parity_Q_ss->SetBinContent(i + 1, cont / res_q / 100.);
    Parity_Q_ss->SetBinError(i + 1, err / res_q / 100.);

    cont = Parity_Q_os2->GetBinContent(i + 1);
    err = Parity_Q_os1->GetBinError(i + 1);
    Parity_Q_os->SetBinContent(i + 1, cont / res_q / 100.);
    Parity_Q_os->SetBinError(i + 1, err / res_q / 100.);
  }

  // rebin the data, three choices...
  const int N_bins = Q2_range * 10 / rebin;
  cout << endl;
  cout << " N bins = " << N_bins << endl;

  float v2q[N_bins], v2q_err[N_bins], d_gq[N_bins], d_gq_err[N_bins];

  if (rebin == 1)
  {
    for (int i = 0; i < N_bins; i++)
    {
      v2q[i] = p_v2_Q_obs_new->GetBinContent(i + 1);
      v2q_err[i] = p_v2_Q_obs_new->GetBinError(i + 1);

      d_gq[i] = Parity_Q_os->GetBinContent(i + 1) - Parity_Q_ss->GetBinContent(i + 1);
      d_gq_err[i] = sqrt(pow(Parity_Q_ss->GetBinError(i + 1), 2) + pow(Parity_Q_ss->GetBinError(i + 1), 2));
    }
  }

  else if (rebin == 2)
  {
    for (int i = 0; i < N_bins; i++)
    {
      double temp1 = p_v2_Q_obs_new->GetBinContent(2 * i + 1);
      double temp1_err = p_v2_Q_obs_new->GetBinError(2 * i + 1);
      double temp2 = p_v2_Q_obs_new->GetBinContent(2 * i + 2);
      double temp2_err = p_v2_Q_obs_new->GetBinError(2 * i + 2);
      v2q[i] = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err);
      v2q_err[i] = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err));

      temp1 = Parity_Q_os->GetBinContent(2 * i + 1);
      temp1_err = Parity_Q_os->GetBinError(2 * i + 1);
      temp2 = Parity_Q_os->GetBinContent(2 * i + 2);
      temp2_err = Parity_Q_os->GetBinError(2 * i + 2);
      double os = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err);
      double os_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err));

      temp1 = Parity_Q_ss->GetBinContent(2 * i + 1);
      temp1_err = Parity_Q_ss->GetBinError(2 * i + 1);
      temp2 = Parity_Q_ss->GetBinContent(2 * i + 2);
      temp2_err = Parity_Q_ss->GetBinError(2 * i + 2);
      double ss = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err);
      double ss_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err));

      d_gq[i] = os - ss;
      d_gq_err[i] = sqrt(pow(os_err, 2) + pow(ss_err, 2));
    }
  }

  else if (rebin == 5)
  {
    double temp1, temp2, temp3, temp4, temp5, temp1_err, temp2_err, temp3_err, temp4_err, temp5_err;
    for (int i = 0; i < N_bins; i++)
    {
      temp1 = p_v2_Q_obs_new->GetBinContent(5 * i + 1);
      temp1_err = p_v2_Q_obs_new->GetBinError(5 * i + 1);
      temp2 = p_v2_Q_obs_new->GetBinContent(5 * i + 2);
      temp2_err = p_v2_Q_obs_new->GetBinError(5 * i + 2);
      temp3 = p_v2_Q_obs_new->GetBinContent(5 * i + 3);
      temp3_err = p_v2_Q_obs_new->GetBinError(5 * i + 3);
      temp4 = p_v2_Q_obs_new->GetBinContent(5 * i + 4);
      temp4_err = p_v2_Q_obs_new->GetBinError(5 * i + 4);
      temp5 = p_v2_Q_obs_new->GetBinContent(5 * i + 5);
      temp5_err = p_v2_Q_obs_new->GetBinError(5 * i + 5);

      v2q[i] = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err + temp3 / temp3_err / temp3_err + temp4 / temp4_err / temp4_err + temp5 / temp5_err / temp5_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err);
      v2q_err[i] = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err));

      temp1 = Parity_Q_os->GetBinContent(5 * i + 1);
      temp1_err = Parity_Q_os->GetBinError(5 * i + 1);
      temp2 = Parity_Q_os->GetBinContent(5 * i + 2);
      temp2_err = Parity_Q_os->GetBinError(5 * i + 2);
      temp3 = Parity_Q_os->GetBinContent(5 * i + 3);
      temp3_err = Parity_Q_os->GetBinError(5 * i + 3);
      temp4 = Parity_Q_os->GetBinContent(5 * i + 4);
      temp4_err = Parity_Q_os->GetBinError(5 * i + 4);
      temp5 = Parity_Q_os->GetBinContent(5 * i + 5);
      temp5_err = Parity_Q_os->GetBinError(5 * i + 5);
      double os = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err + temp3 / temp3_err / temp3_err + temp4 / temp4_err / temp4_err + temp5 / temp5_err / temp5_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err);
      double os_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err));
      temp1 = Parity_Q_ss->GetBinContent(5 * i + 1);
      temp1_err = Parity_Q_ss->GetBinError(5 * i + 1);
      temp2 = Parity_Q_ss->GetBinContent(5 * i + 2);
      temp2_err = Parity_Q_ss->GetBinError(5 * i + 2);
      temp3 = Parity_Q_ss->GetBinContent(5 * i + 3);
      temp3_err = Parity_Q_ss->GetBinError(5 * i + 3);
      temp4 = Parity_Q_ss->GetBinContent(5 * i + 4);
      temp4_err = Parity_Q_ss->GetBinError(5 * i + 4);
      temp5 = Parity_Q_ss->GetBinContent(5 * i + 5);
      temp5_err = Parity_Q_ss->GetBinError(5 * i + 5);
      double ss = (temp1 / temp1_err / temp1_err + temp2 / temp2_err / temp2_err + temp3 / temp3_err / temp3_err + temp4 / temp4_err / temp4_err + temp5 / temp5_err / temp5_err) / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err);
      double ss_err = sqrt(1 / (1 / temp1_err / temp1_err + 1 / temp2_err / temp2_err + 1 / temp3_err / temp3_err + 1 / temp4_err / temp4_err + 1 / temp5_err / temp5_err));

      d_gq[i] = os - ss;
      d_gq_err[i] = sqrt(pow(os_err, 2) + pow(ss_err, 2));
    }
  }

  ///////////////////
  /**     Draw    **/
  ///////////////////

  gStyle->SetPalette(1);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetOptTitle(0);
  gStyle->SetOptDate(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetFrameFillColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetPadBorderSize(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetLabelSize(.05, "X");
  gStyle->SetLabelSize(.05, "Y");
  gStyle->SetTitleSize(.06, "X");
  gStyle->SetTitleSize(.06, "Y");

  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(kWhite);
  gStyle->SetCanvasColor(kWhite);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetTitleBorderSize(0);

  // Q2-res
  TCanvas *c1 = new TCanvas("Resolution");
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
  p_res->GetXaxis()->SetRangeUser(0, 50);
  p_res->GetXaxis()->SetTitle("Q_{2, parent}^{2}");
  p_res->GetXaxis()->CenterTitle();
  p_res->GetYaxis()->SetTitle("Resolution");
  p_res->Draw();
  f_res->SetLineColor(2);
  f_res->Draw("same");

  c1->Print(printname1);

  /// Q2-v2
  TCanvas *c2 = new TCanvas("v2 parent - Q2 parent");
  c2->SetTopMargin(0.06);
  c2->SetRightMargin(0.1);
  c2->SetLeftMargin(0.15);
  c2->SetBottomMargin(0.17);
  c2->Draw();
  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  p_v2_Q_obs_new->GetXaxis()->SetRangeUser(0, 50);
  // p_v2_Q_obs_new->SetMaximum(Yrange_v2[cen-1]);
  // p_v2_Q_obs_new->SetMinimum(-Yrange_v2[cen-1]);
  p_v2_Q_obs_new->GetXaxis()->SetTitle("Q_{2, parent}^{2}");
  p_v2_Q_obs_new->GetXaxis()->CenterTitle();
  p_v2_Q_obs_new->GetYaxis()->SetTitle("v_{2,parent}");
  p_v2_Q_obs_new->Draw();

  c2->Print(printname2);

  // raw plots
  TCanvas *can3 = new TCanvas("Flow_v", "Flow_v", 40, 40, 740, 500);
  can3->SetTopMargin(0.06);
  can3->SetRightMargin(0.1);
  can3->SetLeftMargin(0.15);
  can3->SetBottomMargin(0.17);
  can3->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  TString *histGraphName = new TString("Fl");
  TH1F *histGraph = new TH1F(histGraphName->Data(), "", 2 * Xrange_v2_raw, -Xrange_v2_raw, Xrange_v2_raw);
  histGraph->SetMaximum(Yrange_delGamma0[cen - 1]);
  histGraph->SetMinimum(-Yrange_delGamma0[cen - 1]);
  histGraph->SetLineColor(kBlack);
  histGraph->GetYaxis()->SetTitleOffset(0.9);
  histGraph->GetYaxis()->SetTitleSize(0.065);
  histGraph->GetXaxis()->SetTitleSize(0.08);
  histGraph->GetXaxis()->SetTitleOffset(0.90);
  histGraph->GetYaxis()->CenterTitle();
  histGraph->GetXaxis()->SetTitle("v_{2, parent} raw");
  histGraph->GetXaxis()->CenterTitle();
  histGraph->GetYaxis()->SetTitle("#Delta#gamma^{112} raw");
  histGraph->GetXaxis()->SetNdivisions(6);
  histGraph->GetYaxis()->SetNdivisions(605);
  double lsize = histGraph->GetLabelSize();
  histGraph->GetYaxis()->SetLabelSize(lsize * 1.0);
  histGraph->GetXaxis()->SetLabelSize(lsize * 1.0);
  histGraph->Draw();

  TGraphErrors *g112_v2 = new TGraphErrors(Nbin, v2_parent_raw, DG112_raw, 0, DG112_raw_err); // turn off the v2 error
  g112_v2->SetMarkerStyle(kOpenStar);
  g112_v2->SetMarkerSize(1.5);
  g112_v2->SetMarkerColor(2);
  g112_v2->SetLineColor(2);
  g112_v2->SetFillColor(2);
  g112_v2->SetLineStyle(1);
  g112_v2->SetLineWidth(2);
  g112_v2->Draw("pe1");

  TF1 *fit_g112_v2 = new TF1("fit_g112_v2", "[0]+[1]*x", -Xrange_v2_raw, Xrange_v2_raw);
  fit_g112_v2->SetLineStyle(4);
  fit_g112_v2->SetLineColor(2);
  fit_g112_v2->SetLineWidth(4);
  g112_v2->Fit("fit_g112_v2", "E", "", -Xrange_v2_raw, Xrange_v2_raw);

  float fit_sig = fit_g112_v2->GetParameter(0);
  float fit_sig_err = fit_g112_v2->GetParError(0);
  cout << fit_sig << endl;
  cout << "raw ESE signal = " << fit_sig << " +/- " << fit_sig_err << endl
       << endl;

  TLatex *tex = new TLatex(-Xrange_v2_raw * 0.8, Yrange_delGamma0[cen - 1] * 0.8, centrality);
  tex->SetTextSize(0.08);
  tex->SetTextColor(1);
  tex->Draw();

  TLegend *legend = new TLegend(0.67, 0.2, 0.77, 0.26);
  legend->SetFillColor(0);
  legend->SetTextSize(0.065);
  legend->SetLineColor(0);
  legend->SetBorderSize(0);
  legend->SetLineStyle(3);
  legend->AddEntry(g112_v2, "#Delta#gamma_{112} #pi-TPC", "p");
  legend->Draw();

  can3->Print(printname3);

  /// Q2-G112 ssos
  TCanvas *c4 = new TCanvas("Parity parent Q2");
  c4->SetTopMargin(0.06);
  c4->SetRightMargin(0.1);
  c4->SetLeftMargin(0.15);
  c4->SetBottomMargin(0.17);
  c4->Draw();
  Parity_Q_ss->SetLineColor(2);
  Parity_Q_ss->SetMarkerColor(2);
  Parity_Q_ss->GetXaxis()->SetRangeUser(0, 50);
  Parity_Q_ss->SetMaximum(YRangeMax_ew_Q2[cen - 1]);
  Parity_Q_ss->SetMinimum(-YRangeMax_ew_Q2[cen - 1]);
  Parity_Q_ss->GetXaxis()->SetTitle("Q_{2,parent}^{2}");
  Parity_Q_ss->GetXaxis()->CenterTitle();
  Parity_Q_ss->GetYaxis()->SetTitle("#gamma^{112}");
  Parity_Q_ss->Draw();
  Parity_Q_os->SetLineColor(4);
  Parity_Q_os->SetMarkerColor(4);
  Parity_Q_os->Draw("same");

  TF1 *fu = new TF1("fu", "[0]+[1]*x+[2]*x*x", 0, 50);
  Parity_Q_ss->Fit("fu", "0", "", 0, 50);
  float sig_ss = fu->GetParameter(0);
  float sig_ss_err = fu->GetParError(0);
  fu->SetLineColor(2);
  fu->Draw("same");
  Parity_Q_os->Fit("fu", "0", "", 0, 50);
  float sig_os = fu->GetParameter(0);
  float sig_os_err = fu->GetParError(0);
  fu->SetLineColor(4);
  fu->Draw("same");
  cout << endl;
  cout << "#gamma_{ss} at zero q = " << sig_ss << " +/- " << sig_ss_err << endl;
  cout << "#gamma_{os} at zero q = " << sig_os << " +/- " << sig_os_err << endl;
  cout << "difference = " << sig_os - sig_ss << " +/- " << sqrt(sig_ss_err * sig_ss_err + sig_os_err * sig_os_err) << endl;

  TLatex *tex1 = new TLatex(1, -YRangeMax_ew_Q2[cen - 1] * 0.8, "same sign");
  tex1->SetTextSize(0.06);
  tex1->SetTextColor(2);
  tex1->Draw();
  TLatex *tex2 = new TLatex(1, -YRangeMax_ew_Q2[cen - 1] * 0.6, "opposite sign");
  tex2->SetTextSize(0.06);
  tex2->SetTextColor(4);
  tex2->Draw();

  tex1->Draw();
  tex2->Draw();
  c4->Print(printname4);

  // Q2-DG112-v2
  TCanvas *can5 = new TCanvas("parent_v2(Q2)_DG(Q2)", "parent v2(Q2) DG(Q2)", 40, 40, 740, 500);
  can5->SetTopMargin(0.06);
  can5->SetRightMargin(0.1);
  can5->SetLeftMargin(0.15);
  can5->SetBottomMargin(0.17);
  can5->Draw();

  gPad->SetGridx(0);
  gPad->SetGridy(0);
  gPad->SetTickx();
  gPad->SetTicky();
  TString *histGraphName2 = new TString("parent_v2(Q2)_DG(Q2)");
  TH1F *histGraph2 = new TH1F(histGraphName2->Data(), "", 100, 0, Xrange_delGamma[cen - 1]);
  histGraph2->SetMaximum(Yrange_delGamma[cen - 1]);
  histGraph2->SetMinimum(-Yrange_delGamma[cen - 1]);
  histGraph2->SetLineColor(kBlack);
  histGraph2->GetYaxis()->SetTitleOffset(0.67);
  histGraph2->GetYaxis()->SetTitleSize(0.065);
  histGraph2->GetXaxis()->SetTitleSize(0.08);
  histGraph2->GetXaxis()->SetTitleOffset(0.90);
  histGraph2->GetYaxis()->CenterTitle();
  histGraph2->GetXaxis()->SetTitle("v_{2, parent}");
  histGraph2->GetXaxis()->CenterTitle();
  histGraph2->GetYaxis()->SetTitle("#Delta#gamma^{112}");
  histGraph2->GetXaxis()->SetNdivisions(6);
  histGraph2->GetYaxis()->SetNdivisions(605);
  histGraph2->GetYaxis()->SetLabelSize(lsize * 1.0);
  histGraph2->GetXaxis()->SetLabelSize(lsize * 1.0);
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

  TF1 *fit_g112_v2_2 = new TF1("fit_g112_v2_2", "[0]+[1]*x", 0, Xfitting_delGamma[cen - 1]);
  fit_g112_v2_2->SetLineStyle(4);
  fit_g112_v2_2->SetLineColor(2);
  fit_g112_v2_2->SetLineWidth(4);
  g112_v2_2->Fit("fit_g112_v2_2", "E", "", 0, Xfitting_delGamma[cen - 1]);

  float fit_sig2 = fit_g112_v2_2->GetParameter(0);
  float fit_sig2_err = fit_g112_v2_2->GetParError(0);
  cout << fit_sig2 << endl;
  cout << "ESE signal(Q2) = " << fit_sig2 << " +/- " << fit_sig2_err << endl
       << endl;

  TLatex *tex_new = new TLatex(Xrange_delGamma[cen - 1] * 0.1, Yrange_delGamma[cen - 1] * 0.8, centrality);
  tex_new->SetTextSize(0.08);
  tex_new->SetTextColor(1);
  tex_new->Draw();
  tex_new->Draw();

  TLegend *legend2 = new TLegend(0.27, 0.22, 0.37, 0.42);
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.065);
  legend2->SetLineColor(0);
  legend2->SetBorderSize(0);
  legend2->SetLineStyle(3);
  legend2->AddEntry(g112_v2, "#Delta#gamma_{112}", "p");
  legend2->Draw();

  can5->Print(printname5);
}
/////////////////////////////////////////////
static Double_t chi(double res)
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
//-------------------------
static Double_t resEventPlane(double chi)
{
  double con = 0.626657; // sqrt(pi/2)/2
  double arg = chi * chi / 4.;

  Double_t res = con * chi * exp(-arg) * (TMath::BesselI0(arg) + TMath::BesselI1(arg));

  return res;
}
