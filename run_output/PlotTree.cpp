void MakeDoubleProfile(const char* filter1, const char* filter2, const char* plot1, const char* plot2,
    const char* name1, const char* name2,
    const char* title1, const char* title2, const char* xlab1, const char* xlab2, const char* ylab1, const char* ylab2)
{
    tree.Draw(filter1, plot1, "prof");
    tree.Draw(filter2, plot2, "prof");
    TProfile* prof_1 = (TProfile*) gDirectory->Get(name1);
    TProfile* prof_2 = (TProfile*) gDirectory->Get(name2);
    prof_1->Write();
    prof_2->Write();
    TCanvas c_prof("c_prof", "Profile Canvas", 432, 500);
    c_prof.Divide(1, 2);
    prof_1->SetTitle(title1);
    prof_1->SetXTitle(xlab1);
    prof_1->SetYTitle(ylab1);
    prof_2->SetTitle(title2);
    prof_2->SetXTitle(xlab2);
    prof_2->SetYTitle(ylab2);
    c_prof.cd(1);
    prof_1->Draw();
    c_prof.cd(2);
    prof_2->Draw();
    c_prof.Write();
}

void MakeDoubleHisto(const char* filter1, const char* filter2, const char* plot1, const char* plot2,
    const char* name1, const char* name2,
    const char* title1, const char* title2)
{
    tree.Draw(filter1, plot1);
    tree.Draw(filter2, plot2);
    TH1* mono_im = (TH1*) gDirectory->Get(name1);
    TH1* chkv_im = (TH1*) gDirectory->Get(name2);
    mono_im->Write();
    chkv_im->Write();
    TCanvas c_hist("c_hist", "Histogram Canvas", 432, 500);
    c_hist.Divide(1, 2);
    mono_im->SetTitle(title1);
    chkv_im->SetTitle(title2);
    c_im_hist.cd(1);
    mono_im->Draw();
    c_im_hist.cd(2);
    chkv_im->Draw();
    c_im_hist.Write();
}

void PlotTree(const char* csv_file)
{
    // Read the CSV file into a TTree and open a file for writing
    TTree tree;
    const char* branch_desc = "seed:id:energy:psi:im:trig:mono_psi:mono_im:chkv:chkv_psi:chkv_im";
    tree.ReadFile(csv_file, branch_desc, ',');
    TFile file("PlotTreeOut.root", "RECREATE");

    // Make two 1d histograms of percent error in impact
    MakeDoubleHisto("chkv", "chkv", "(mono_im - im) / im >> mono_im(40, -1, 2)", "(chkv_im - im) / im >> chkv_im(40, -1, 2)",
        "mono_im", "chkv_im", "Fractional Monocular Impact Errors", "Fractional Cherenkov Impact Errors");

    // Make two 1d histograms of absolute error in psi
    MakeDoubleHisto("chkv", "chkv", "mono_psi - psi >> mono_psi", "chkv_psi - psi >> chkv_psi",
        "mono_psi", "chkv_psi", "Monocular Angular Errors", "Cherenkov Angular Errors");

    // Make two TProfile histograms of absolute error in psi vs psi
    MakeDoubleProfile("chkv", "chkv", "mono_psi - psi : psi >> mono_psi_psi(50, 0, 3.14)", "chkv_psi - psi : psi >> chkv_psi_psi(50, 0, 3.14)",
        "mono_psi_psi", "chkv_psi_psi",
        "Monocular Angular Errors", "Cherenkov Angular Errors", "Actual Angle (rad)", "Actual Angle (rad)", "Angle Error (rad)", "Angle Error (rad)");

    // Make two TProfile histograms of percent error in impact vs impact
    MakeDoubleProfile("chkv", "chkv", "(mono_im - im) / im : im >> mono_im_im", "(chkv_im - im) / im : im >> chkv_im_im",
        "mono_im_im", "chkv_im_im",
        "Fractional Monocular Impact Errors", "Fractional Cherenkov Impact Errors", "Impact Parameter (km)", "Impact Parameter (km)", "", "");

    // TODO: Make two TProfile histograms of absolute error in psi vs log(E)

    // TODO: Make two TProfile histograms of percent error in impact vs log(E)

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
}
