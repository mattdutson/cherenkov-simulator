string CreationString(const char* strn, const char* name, size_t n_bins, double min, double max)
{
    stringstream stream;
    stream << strn << " >> " << name << "(" << n_bins << ", " << min << ", " << max << ")";
    return stream.str().c_str();
}

struct Params
{
    const char* mono_name;
    const char* chkv_name;
    const char* mono_strn;
    const char* chkv_strn;
    const char* mono_titl;
    const char* chkv_titl;
    const char* canv_name;
    const char* filter = "chkv > 0";
    int n_bins;
    double min;
    double max;
};

void MakeDoubleHisto(TTree& tree, Params par)
{
    tree.Draw(CreationString(par.mono_strn, par.mono_name, par.n_bins, par.min, par.max).c_str(), par.filter);
    tree.Draw(CreationString(par.chkv_strn, par.chkv_name, par.n_bins, par.min, par.max).c_str(), par.filter);
    TH1* mono_im = (TH1*) gDirectory->Get(par.mono_name);
    TH1* chkv_im = (TH1*) gDirectory->Get(par.chkv_name);
    mono_im->Write();
    chkv_im->Write();
    TCanvas c_hist(par.canv_name, "Histogram Canvas", 432, 500);
    c_hist.Divide(1, 2);
    mono_im->SetTitle(par.mono_titl);
    chkv_im->SetTitle(par.chkv_titl);
    c_hist.cd(1);
    mono_im->Draw();
    c_hist.cd(2);
    chkv_im->Draw();
    c_hist.Write();
}

void MakeDoubleProfile(TTree& tree, Params par, const char* x_lab, const char* y_lab, double min_y, double max_y)
{
    tree.Draw(CreationString(par.mono_strn, par.mono_name, par.n_bins, par.min, par.max).c_str(), par.filter, "prof");
    tree.Draw(CreationString(par.chkv_strn, par.chkv_name, par.n_bins, par.min, par.max).c_str(), par.filter, "prof");
    TProfile* prof_1 = (TProfile*) gDirectory->Get(par.mono_name);
    TProfile* prof_2 = (TProfile*) gDirectory->Get(par.chkv_name);
    prof_1->Write();
    prof_2->Write();
    TCanvas c_prof(par.canv_name, "Profile Canvas", 432, 500);
    c_prof.Divide(1, 2);
    prof_1->SetTitle(par.mono_titl);
    prof_1->SetXTitle(x_lab);
    prof_1->SetYTitle(y_lab);
    prof_1->SetMinimum(min_y);
    prof_1->SetMaximum(max_y);
    prof_2->SetTitle(par.chkv_titl);
    prof_2->SetXTitle(x_lab);
    prof_2->SetYTitle(y_lab);
    prof_2->SetMinimum(min_y);
    prof_2->SetMaximum(max_y);
    c_prof.cd(1);
    prof_1->Draw();
    c_prof.cd(2);
    prof_2->Draw();
    c_prof.Write();

}

void PlotResults(const char* csv_file)
{
    TTree tree;
    const char* branch_desc = "seed:id:energy:psi:im:gnd:trig:mono_psi:mono_im:mono_gnd:chkv:chkv_psi:chkv_im:chkv_gnd";
    tree.ReadFile(csv_file, branch_desc, ',');
    TFile file("Results.root", "RECREATE");
    Params par;

    par.mono_name = "mono_im_err";
    par.chkv_name = "chkv_im_err";
    par.mono_strn = "(mono_im - im) / im";
    par.chkv_strn = "(chkv_im - im) / im";
    par.mono_titl = "Monocular Impact Errors";
    par.chkv_titl = "Cherenkov Impact Errors";
    par.canv_name = "impact_err";
    par.n_bins = 80;
    par.min = -1;
    par.max =  1;
    MakeDoubleHisto(tree, par);

    par.mono_name = "mono_psi_err";
    par.chkv_name = "chkv_psi_err";
    par.mono_strn = "mono_psi - psi";
    par.chkv_strn = "chkv_psi - psi";
    par.mono_titl = "Monocular Angular Errors";
    par.chkv_titl = "Cherenkov Angular Errors";
    par.canv_name = "psi_err";
    par.n_bins = 80;
    par.min = -90;
    par.max =  90;
    MakeDoubleHisto(tree, par);

    par.mono_name = "mono_im_psi";
    par.chkv_name = "chkv_im_psi";
    par.mono_strn = "(mono_im - im) / im : psi";
    par.chkv_strn = "(chkv_im - im) / im : psi";
    par.mono_titl = "Monocular Impact Error vs Angle";
    par.chkv_titl = "Cherenkov Impact Error vs Angle";
    par.canv_name = "im_err_psi";
    par.n_bins = 40;
    par.min = 0;
    par.max = 180;
    MakeDoubleProfile(tree, par, "Shower Angle", "Impact Error", -1, 1);

    par.mono_name = "mon_im_im";
    par.chkv_name = "ckv_im_im";
    par.mono_strn = "(mono_im - im) / im : im";
    par.chkv_strn = "(chkv_im - im) / im : im";
    par.mono_titl = "Monocular Impact Error vs Impact";
    par.chkv_titl = "Cherenkov Impact Error vs Impact";
    par.canv_name = "im_err_im";
    par.n_bins = 40;
    par.min = 0;
    par.max = 40;
    MakeDoubleProfile(tree, par, "Impact Parameter", "Impact Error", -1, 1);

    par.mono_name = "mon_im_gnd";
    par.chkv_name = "ckv_im_gnd";
    par.mono_strn = "(mono_im - im) / im : gnd";
    par.chkv_strn = "(chkv_im - im) / im : gnd";
    par.mono_titl = "Monocular Impact Error vs Ground Distance";
    par.chkv_titl = "Cherenkov Impact Error vs Ground Distance";
    par.canv_name = "im_err_gnd";
    par.n_bins = 40;
    par.min = 0;
    par.max = 60;
    MakeDoubleProfile(tree, par, "Ground Distance", "Impact Error", -1, 1);

    par.mono_name = "mon_im_en";
    par.chkv_name = "ckv_im_en";
    par.mono_strn = "(mono_im - im) / im : log(energy) / log(10)";
    par.chkv_strn = "(chkv_im - im) / im : log(energy) / log(10)";
    par.mono_titl = "Monocular Impact Error vs Log(Energy)";
    par.chkv_titl = "Cherenkov Impact Error vs Log(Energy)";
    par.canv_name = "im_err_en";
    par.n_bins = 40;
    par.min = 17;
    par.max = 21;
    MakeDoubleProfile(tree, par, "Log(Energy)", "Impact Error", -1, 1);

    // We would have to make plots against angle, but as this profile shows, but impact and angle errors are highly
    // correlated.
    par.mono_name = "mon_imerr_psierr";
    par.chkv_name = "ckv_imerr_psierr";
    par.mono_strn = "(mono_im - im) / im : mono_psi - psi";
    par.chkv_strn = "(chkv_im - im) / im : chkv_psi - psi";
    par.mono_titl = "Monocular Impact Error vs Angular Error";
    par.chkv_titl = "Cherenkov Impact Error vs Angular Error";
    par.canv_name = "im_err_psi_err";
    par.n_bins = 40;
    par.min = -90;
    par.max =  90;
    MakeDoubleProfile(tree, par, "Angular Error", "Impact Error", -1, 1);
}
