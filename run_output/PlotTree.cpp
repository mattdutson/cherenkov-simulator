void PlotTree(const char* csv_file)
{
    // Read the CSV file into a TTree and open a file for writing
    TTree tree;
    const char* branch_desc = "seed:id:energy:psi:im:trig:mono_psi:mono_im:chkv:chkv_psi:chkv_im";
    tree.ReadFile(csv_file, branch_desc, ',');
    TFile file("PlotTreeOut.root", "RECREATE");

    // Make two 1d histograms of percent error in impact
    tree.Draw("(mono_im - im) / im >> mono_im(40, -1, 2)", "chkv");
    tree.Draw("(chkv_im - im) / im >> chkv_im(40, -1, 2)", "chkv");
    TH1* mono_im = (TH1*) gDirectory->Get("mono_im");
    TH1* chkv_im = (TH1*) gDirectory->Get("chkv_im");
    mono_im->Write();
    chkv_im->Write();
    TCanvas c_hist("c_hist", "Histogram Canvas", 432, 500);
    c_hist.Divide(1, 2);
    mono_im->SetTitle("Fractional Monocular Impact Errors");
    chkv_im->SetTitle("Fractional Cherenkov Impact Errors");
    c_hist.cd(1);
    mono_im->Draw();
    c_hist.cd(2);
    chkv_im->Draw();
    c_hist.Write();

    // Make two TProfile histograms of absolute error in psi vs psi
    tree.Draw("mono_psi - psi : psi >> mono_psi_psi(50, 0, 3.14)", "chkv", "prof");
    tree.Draw("chkv_psi - psi : psi >> chkv_psi_psi(50, 0, 3.14)", "chkv", "prof");
    TProfile* mono_psi_psi = (TProfile*) gDirectory->Get("mono_psi_psi");
    TProfile* chkv_psi_psi = (TProfile*) gDirectory->Get("chkv_psi_psi");
    mono_psi_psi->Write();
    chkv_psi_psi->Write();
    TCanvas c_prof("c_prof", "Profile Canvas", 432, 500);
    c_prof.Divide(1, 2);
    mono_psi_psi->SetTitle("Monocular Angular Errors");
    mono_psi_psi->SetXTitle("Actual Angle (rad)");
    mono_psi_psi->SetYTitle("Angle Error (rad)");
    chkv_psi_psi->SetTitle("Cherenkov Angular Errors");
    chkv_psi_psi->SetXTitle("Actual Angle (rad)");
    chkv_psi_psi->SetYTitle("Angle Error (rad)");
    c_prof.cd(1);
    mono_psi_psi->Draw();
    c_prof.cd(2);
    chkv_psi_psi->Draw();
    c_prof.Write();

    // Make two 1d histograms of absolute error in psi
    // tree.Draw("mono_psi - psi >> mono_psi", "chkv");
    // tree.Draw("chkv_psi - psi >> chkv_psi", "chkv");
    // gDirectory->Get("mono_psi")->Write();
    // gDirectory->Get("chkv_psi")->Write();

    // Make two TProfile histograms of absolute error in psi vs impact
    // tree.Draw("mono_psi - psi : im >> mono_psi_im", "chkv", "prof");
    // tree.Draw("chkv_psi - psi : im >> chkv_psi_im", "chkv", "prof");
    // gDirectory->Get("mono_psi_im")->Write();
    // gDirectory->Get("chkv_psi_im")->Write();

    // Make two TProfile histograms of percent error in impact vs psi
    // tree.Draw("(mono_im - im) / im : psi >> mono_im_psi", "chkv", "prof");
    // tree.Draw("(chkv_im - im) / im : psi >> chkv_im_psi", "chkv", "prof");
    // gDirectory->Get("mono_im_psi")->Write();
    // gDirectory->Get("chkv_im_psi")->Write();

    // Make two TProfile histograms of percent error in impact vs impact
    // tree.Draw("(mono_im - im) / im : im >> mono_im_im", "chkv", "prof");
    // tree.Draw("(chkv_im - im) / im : im >> chkv_im_im", "chkv", "prof");
    // gDirectory->Get("mono_im_im")->Write();
    // gDirectory->Get("chkv_im_im")->Write();
}
