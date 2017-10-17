string CreationString(const char* strn, const char* name, size_t n_bins, double min, double max)
{
    stringstream stream;
    stream << strn << " >> " << name << "(" << n_bins << ", " << min << ", " << max << ")";
    return stream.str().c_str();
}

struct Params
{
    const char* mon_name;
    const char* ckv_name;
    const char* mon_strn;
    const char* ckv_strn;
    const char* mon_titl;
    const char* ckv_titl;
    const char* can_name;
    const char* filter = "chkv";
    int n_bins;
    double min;
    double max;
};

void MakeDoubleHisto(TTree& tree, Params par)
{
    tree.Draw(CreationString(par.mon_strn, par.mon_name, par.n_bins, par.min, par.max).c_str(), par.filter);
    tree.Draw(CreationString(par.ckv_strn, par.ckv_name, par.n_bins, par.min, par.max).c_str(), par.filter);
    TH1* mono_im = (TH1*) gDirectory->Get(par.mon_name);
    TH1* chkv_im = (TH1*) gDirectory->Get(par.ckv_name);
    mono_im->Write();
    chkv_im->Write();
    TCanvas c_hist(par.can_name, "Histogram Canvas", 432, 500);
    c_hist.Divide(1, 2);
    mono_im->SetTitle(par.mon_titl);
    chkv_im->SetTitle(par.ckv_titl);
    c_hist.cd(1);
    mono_im->Draw();
    c_hist.cd(2);
    chkv_im->Draw();
    c_hist.Write();
}

void MakeDoubleProfile(TTree& tree, Params par, const char* x_lab, const char* y_lab)
{
    tree.Draw(CreationString(par.mon_strn, par.mon_name, par.n_bins, par.min, par.max).c_str(), par.filter, "prof");
    tree.Draw(CreationString(par.ckv_strn, par.ckv_name, par.n_bins, par.min, par.max).c_str(), par.filter, "prof");
    TProfile* prof_1 = (TProfile*) gDirectory->Get(par.mon_name);
    TProfile* prof_2 = (TProfile*) gDirectory->Get(par.ckv_name);
    prof_1->Write();
    prof_2->Write();
    TCanvas c_prof(par.can_name, "Profile Canvas", 432, 500);
    c_prof.Divide(1, 2);
    prof_1->SetTitle(par.mon_titl);
    prof_1->SetXTitle(x_lab);
    prof_1->SetYTitle(y_lab);
    prof_2->SetTitle(par.ckv_titl);
    prof_2->SetXTitle(x_lab);
    prof_2->SetYTitle(y_lab);
    c_prof.cd(1);
    prof_1->Draw();
    c_prof.cd(2);
    prof_2->Draw();
    c_prof.Write();
}

void PlotResults(const char* csv_file)
{
    TTree tree;
    const char* branch_desc = "seed:id:energy:psi:im:trig:mono_psi:mono_im:chkv:chkv_psi:chkv_im";
    tree.ReadFile(csv_file, branch_desc, ',');
    TFile file("Results.root", "RECREATE");
    Params par;

    par.mon_name = "mono_im_err";
    par.ckv_name = "chkv_im_err";
    par.mon_strn = "(mono_im - im) / im";
    par.ckv_strn = "(chkv_im - im) / im";
    par.mon_titl = "Monocular Impact Errors";
    par.ckv_titl = "Cherenkov Impact Errors";
    par.can_name = "impact_err";
    par.n_bins = 80;
    par.min = -1;
    par.max =  1;
    MakeDoubleHisto(tree, par);

    par.mon_name = "mono_psi_err";
    par.ckv_name = "chkv_psi_err";
    par.mon_strn = "mono_psi - psi";
    par.ckv_strn = "chkv_psi - psi";
    par.mon_titl = "Monocular Angular Errors";
    par.ckv_titl = "Cherenkov Angular Errors";
    par.can_name = "psi_err";
    par.n_bins = 80;
    par.min = -1;
    par.max =  1;
    MakeDoubleHisto(tree, par);

    par.mon_name = "mono_im_psi";
    par.ckv_name = "chkv_im_psi";
    par.mon_strn = "(mono_im - im) / im : psi";
    par.ckv_strn = "(chkv_im - im) / im : psi";
    par.mon_titl = "Monocular Impact Error vs Angle";
    par.ckv_titl = "Cherenkov Impact Error vs Angle";
    par.can_name = "im_psi_err";
    par.n_bins = 40;
    par.min = 0;
    par.max = 3.14;
    MakeDoubleProfile(tree, par, "Shower Angle", "Impact Error");

    par.mon_name = "mon_im_im";
    par.ckv_name = "ckv_im_im";
    par.mon_strn = "(mono_im - im) / im : im";
    par.ckv_strn = "(chkv_im - im) / im : im";
    par.mon_titl = "Monocular Impact Error vs Impact";
    par.ckv_titl = "Cherenkov Impact Error vs Impact";
    par.can_name = "im_im_err";
    par.n_bins = 40;
    par.min = 0;
    par.max = 40;
    MakeDoubleProfile(tree, par, "Impact Parameter", "Impact Error");

    par.mon_name = "mon_im_en";
    par.ckv_name = "ckv_im_en";
    par.mon_strn = "(mono_im - im) / im : log(energy) / log(10)";
    par.ckv_strn = "(chkv_im - im) / im : log(energy) / log(10)";
    par.mon_titl = "Monocular Impact Error vs Log(Energy)";
    par.ckv_titl = "Cherenkov Impact Error vs Log(Energy)";
    par.can_name = "im_en_err";
    par.n_bins = 40;
    par.min = 17;
    par.max = 21;
    MakeDoubleProfile(tree, par, "Log(Energy)", "Impact Error");

    // We would have to make plots against angle, but as this profile shows, but impact and angle errors are highly
    // correlated.
    par.mon_name = "mon_imerr_psierr";
    par.ckv_name = "ckv_imerr_psierr";
    par.mon_strn = "(mono_im - im) / im : mono_psi - psi";
    par.ckv_strn = "(chkv_im - im) / im : chkv_psi - psi";
    par.mon_titl = "Monocular Impact Error vs Angular Error";
    par.ckv_titl = "Cherenkov Impact Error vs Angular Error";
    par.can_name = "imerr_psierr";
    par.n_bins = 40;
    par.min = -3.14;
    par.max = 3.14;
    MakeDoubleProfile(tree, par, "Angular Error", "Impact Error");
}
