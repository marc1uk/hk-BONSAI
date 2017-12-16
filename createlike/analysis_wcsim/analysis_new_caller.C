void analysis_new_caller(const char* inpath){
gROOT->ProcessLine(".L $ROOTSYS/tutorials/fit/langaus.C");
std::string bonsaidir = gSystem->Getenv("BONSAIDIR");
std::string theline = ".x " + bonsaidir + "/../../createlike_nickwp/analysis_wcsim/analysis_new.C(\"" + inpath + "\")";
gROOT->ProcessLine(theline.c_str());
}
