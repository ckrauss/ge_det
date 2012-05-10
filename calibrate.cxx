
#include "iostream"
#include "TString.h"
#include "TFile.h"
#include "MCA_data.h"
#include "TApplication.h"
#include "TCanvas.h"


using namespace std;

int main(int argc, char** argv){
  vector<TString> names;
  for (int i = 1; i<argc;i++){
    names.push_back(TString(argv[i]));
  }
  TApplication t("app",0,0);
  GeCalibrate cal;
  cal.Init(names);
  cal.Calibrate();
  TString name = TString("Calibration") + cal.GetCalibrationDate() + TString(".root");
  TFile y(name,"recreate");
  cal.Write();
  y.Close();
  t.Run();

}
