{
    TFile * ofile = new TFile("xt.test.root","RECREATE");
    TF1 * fl = new TF1("fl_5_0","(x<216)*x*7.75/216+(x>=216)*(7.75+(x-216)*0.55/(400-216.))",0,800);
    TF1 * fr = new TF1("fr_5_0","-(x<216)*x*7.75/216-(x>=216)*(7.75+(x-216)*0.55/(400-216.))",0,800);
    fl->Write();
    fr->Write();
    ofile->Close();
}
