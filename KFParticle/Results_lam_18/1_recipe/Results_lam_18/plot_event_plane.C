void plot_event_plane(){
	TFile f("event_plane_plots.root", "recreate");

	// TFile g("./19GeV_results/sys_0/cen0.gamma112_fullEP_eff_pT02_module.root");
	// TFile g("../27GeV/sys_1/cen0.gamma112_fullEP_eff_pT02_module.root");
	TFile g("../../0cen4.gamma112_fullEP_eff_pT02_module.root");

	TH1D *Hist_TPC_EP_for_before = ((TH2F *) g.Get("Hist_TPC_EP_for"))->ProjectionX();
	Hist_TPC_EP_for_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_TPC_EP_for_before->GetYaxis()->SetTitle("Counts");
	Hist_TPC_EP_for_before->GetYaxis()->SetRangeUser(0, Hist_TPC_EP_for_before->GetMaximum() * 1.02);
	Hist_TPC_EP_for_before->SetTitle("TPC Event Plane Distribution (#eta > 0.1, 19.6GeV, Before Flattening)");

	TH1D *Hist_TPC_EP_bac_before = ((TH2F *) g.Get("Hist_TPC_EP_bac"))->ProjectionX();
	Hist_TPC_EP_bac_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_TPC_EP_bac_before->GetYaxis()->SetTitle("Counts");
	Hist_TPC_EP_bac_before->GetYaxis()->SetRangeUser(0, Hist_TPC_EP_bac_before->GetMaximum() * 1.02);
	Hist_TPC_EP_bac_before->SetTitle("TPC Event Plane Distribution (#eta < 0.1, 19.6GeV, Before Flattening)");

	TH1D *Hist_TPC_EP_full_before = ((TH2F *) g.Get("Hist_TPC_EP_full"))->ProjectionX();
	Hist_TPC_EP_full_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_TPC_EP_full_before->GetYaxis()->SetTitle("Counts");
	Hist_TPC_EP_full_before->GetYaxis()->SetRangeUser(0, Hist_TPC_EP_full_before->GetMaximum() * 1.02);
	Hist_TPC_EP_full_before->SetTitle("Full TPC Event Plane Distribution (19.6GeV, Before Flattening)");

	TH1D *Hist_EPD_EP1_east_before = ((TH2F *) g.Get("Hist_EPD_EP1_east"))->ProjectionX();
	Hist_EPD_EP1_east_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP1_east_before->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP1_east_before->GetYaxis()->SetRangeUser(0, Hist_EPD_EP1_east_before->GetMaximum() * 1.02);
	Hist_EPD_EP1_east_before->SetTitle("1st Order EPD Event Plane Distribution (#eta > 0.1, 19.6GeV, Before Flattening)");

	TH1D *Hist_EPD_EP1_west_before = ((TH2F *) g.Get("Hist_EPD_EP1_west"))->ProjectionX();
	Hist_EPD_EP1_west_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP1_west_before->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP1_west_before->GetYaxis()->SetRangeUser(0, Hist_EPD_EP1_west_before->GetMaximum() * 1.02);
	Hist_EPD_EP1_west_before->SetTitle("1st Order EPD Event Plane Distribution (#eta < 0.1, 19.6GeV, Before Flattening)");

	TH1D *Hist_EPD_EP_east_before = ((TH2F *) g.Get("Hist_EPD_EP_east"))->ProjectionX();
	Hist_EPD_EP_east_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP_east_before->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP_east_before->GetYaxis()->SetRangeUser(0, Hist_EPD_EP_east_before->GetMaximum() * 1.02);
	Hist_EPD_EP_east_before->SetTitle("2nd Order EPD Event Plane Distribution (#eta > 0.1, 19.6GeV, Before Flattening)");

	TH1D *Hist_EPD_EP_west_before = ((TH2F *) g.Get("Hist_EPD_EP_west"))->ProjectionX();
	Hist_EPD_EP_west_before->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP_west_before->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP_west_before->GetYaxis()->SetRangeUser(0, Hist_EPD_EP_west_before->GetMaximum() * 1.02);
	Hist_EPD_EP_west_before->SetTitle("2nd Order EPD Event Plane Distribution (#eta < 0.1, 19.6GeV, Before Flattening)");


	// TH1D *Hist_TPC_EP_for_after = ((TH2F *) g.Get("Hist_TPC_EP_for_flat"))->ProjectionX();
	// Hist_TPC_EP_for_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_TPC_EP_for_after->GetYaxis()->SetTitle("Counts");
	// Hist_TPC_EP_for_after->SetTitle("TPC Event Plane Distribution (#eta > 0.1, 19.6GeV, After Flattening)");

	// TH1D *Hist_TPC_EP_bac_after = ((TH2F *) g.Get("Hist_TPC_EP_bac_flat"))->ProjectionX();
	// Hist_TPC_EP_bac_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_TPC_EP_bac_after->GetYaxis()->SetTitle("Counts");
	// Hist_TPC_EP_bac_after->SetTitle("TPC Event Plane Distribution (#eta < 0.1, 19.6GeV, After Flattening)");

	// TH1D *Hist_TPC_EP_full_after = ((TH2F *) g.Get("Hist_TPC_EP_full_flat"))->ProjectionX();
	// // TH1D *Hist_TPC_EP_full_after = ((TH1D *) g.Get("Hist_TPC_EP_full_m2_flat"));
	// Hist_TPC_EP_full_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_TPC_EP_full_after->GetYaxis()->SetTitle("Counts");
	// Hist_TPC_EP_full_after->SetTitle("Full TPC Event Plane Distribution (19.6GeV, After Flattening)");

	// TH1D *Hist_EPD_EP1_east_after = ((TH2F *) g.Get("Hist_EPD_EP1_east_flat"))->ProjectionX();
	// Hist_EPD_EP1_east_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_EPD_EP1_east_after->GetYaxis()->SetTitle("Counts");
	// Hist_EPD_EP1_east_after->SetTitle("1st Order EPD Event Plane Distribution (#eta > 0.1, 19.6GeV, After Flattening)");

	// TH1D *Hist_EPD_EP1_west_after = ((TH2F *) g.Get("Hist_EPD_EP1_west_flat"))->ProjectionX();
	// Hist_EPD_EP1_west_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_EPD_EP1_west_after->GetYaxis()->SetTitle("Counts");
	// Hist_EPD_EP1_west_after->SetTitle("1st Order EPD Event Plane Distribution (#eta < 0.1, 19.6GeV, After Flattening)");

	// TH1D *Hist_EPD_EP_east_after = ((TH2F *) g.Get("Hist_EPD_EP_east_flat"))->ProjectionX();
	// Hist_EPD_EP_east_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_EPD_EP_east_after->GetYaxis()->SetTitle("Counts");
	// Hist_EPD_EP_east_after->SetTitle("2nd Order EPD Event Plane Distribution (#eta > 0.1, 19.6GeV, After Flattening)");

	// TH1D *Hist_EPD_EP_west_after = ((TH2F *) g.Get("Hist_EPD_EP_west_flat"))->ProjectionX();
	// Hist_EPD_EP_west_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	// Hist_EPD_EP_west_after->GetYaxis()->SetTitle("Counts");
	// Hist_EPD_EP_west_after->SetTitle("2nd Order EPD Event Plane Distribution (#eta < 0.1, 19.6GeV, After Flattening)");

	TH1D *Hist_TPC_EP_for_after = ((TH2F *) g.Get("Hist_TPC_EP_for_flat"))->ProjectionX();
	Hist_TPC_EP_for_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_TPC_EP_for_after->GetYaxis()->SetTitle("Counts");
	Hist_TPC_EP_for_after->GetYaxis()->SetRangeUser(0, Hist_TPC_EP_for_after->GetMaximum() * 1.02);
	Hist_TPC_EP_for_after->SetTitle("TPC Event Plane Distribution (#eta > 0.1, 27GeV, After Flattening)");

	TH1D *Hist_TPC_EP_bac_after = ((TH2F *) g.Get("Hist_TPC_EP_bac_flat"))->ProjectionX();
	Hist_TPC_EP_bac_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_TPC_EP_bac_after->GetYaxis()->SetTitle("Counts");
	Hist_TPC_EP_bac_after->GetYaxis()->SetRangeUser(0, Hist_TPC_EP_bac_after->GetMaximum() * 1.02);
	Hist_TPC_EP_bac_after->SetTitle("TPC Event Plane Distribution (#eta < 0.1, 27GeV, After Flattening)");

	TH1D *Hist_TPC_EP_full_after = ((TH2F *) g.Get("Hist_TPC_EP_full_flat"))->ProjectionX();
	// TH1D *Hist_TPC_EP_full_after = ((TH1D *) g.Get("Hist_TPC_EP_full_m2_flat"));
	Hist_TPC_EP_full_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_TPC_EP_full_after->GetYaxis()->SetTitle("Counts");
	Hist_TPC_EP_full_after->GetYaxis()->SetRangeUser(0, Hist_TPC_EP_full_after->GetMaximum() * 1.02);
	Hist_TPC_EP_full_after->SetTitle("Full TPC Event Plane Distribution (27GeV, After Flattening)");

	TH1D *Hist_EPD_EP1_east_after = ((TH2F *) g.Get("Hist_EPD_EP1_east_flat"))->ProjectionX();
	Hist_EPD_EP1_east_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP1_east_after->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP1_east_after->GetYaxis()->SetRangeUser(0, Hist_EPD_EP1_east_after->GetMaximum() * 1.02);
	Hist_EPD_EP1_east_after->SetTitle("1st Order EPD Event Plane Distribution (#eta > 0.1, 27GeV, After Flattening)");

	TH1D *Hist_EPD_EP1_west_after = ((TH2F *) g.Get("Hist_EPD_EP1_west_flat"))->ProjectionX();
	Hist_EPD_EP1_west_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP1_west_after->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP1_west_after->GetYaxis()->SetRangeUser(0, Hist_EPD_EP1_west_after->GetMaximum() * 1.02);
	Hist_EPD_EP1_west_after->SetTitle("1st Order EPD Event Plane Distribution (#eta < 0.1, 27GeV, After Flattening)");

	TH1D *Hist_EPD_EP_east_after = ((TH2F *) g.Get("Hist_EPD_EP_east_flat"))->ProjectionX();
	Hist_EPD_EP_east_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP_east_after->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP_east_after->GetYaxis()->SetRangeUser(0, Hist_EPD_EP_east_after->GetMaximum() * 1.02);
	Hist_EPD_EP_east_after->SetTitle("2nd Order EPD Event Plane Distribution (#eta > 0.1, 27GeV, After Flattening)");

	TH1D *Hist_EPD_EP_west_after = ((TH2F *) g.Get("Hist_EPD_EP_west_flat"))->ProjectionX();
	Hist_EPD_EP_west_after->GetXaxis()->SetTitle("Event Plane Orientation/ Angle");
	Hist_EPD_EP_west_after->GetYaxis()->SetTitle("Counts");
	Hist_EPD_EP_west_after->GetYaxis()->SetRangeUser(0, Hist_EPD_EP_west_after->GetMaximum() * 1.02);
	Hist_EPD_EP_west_after->SetTitle("2nd Order EPD Event Plane Distribution (#eta < 0.1, 27GeV, After Flattening)");
	

	f.cd();
	Hist_TPC_EP_for_before->Write("Hist_TPC_EP_for_before");
	Hist_TPC_EP_bac_before->Write("Hist_TPC_EP_bac_before");
	Hist_TPC_EP_full_before->Write("Hist_TPC_EP_full_before");
	Hist_EPD_EP1_west_before->Write("Hist_EPD_EP1_west_before");
	Hist_EPD_EP1_east_before->Write("Hist_EPD_EP1_east_before");
	Hist_EPD_EP_west_before->Write("Hist_EPD_EP_west_before");
	Hist_EPD_EP_east_before->Write("Hist_EPD_EP_east_before");

	Hist_TPC_EP_for_after->Write("Hist_TPC_EP_for_after");
	Hist_TPC_EP_bac_after->Write("Hist_TPC_EP_bac_after");
	Hist_TPC_EP_full_after->Write("Hist_TPC_EP_full_after");
	Hist_EPD_EP1_west_after->Write("Hist_EPD_EP1_west_after");
	Hist_EPD_EP1_east_after->Write("Hist_EPD_EP1_east_after");
	Hist_EPD_EP_west_after->Write("Hist_EPD_EP_west_after");
	Hist_EPD_EP_east_after->Write("Hist_EPD_EP_east_after");

	f.Close();
	g.Close();
	// h.Close();
}