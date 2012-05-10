
#include "iostream"
#include "TString.h"
#include "TFile.h"
#include "MCA_data.h"
#include "TApplication.h"
#include "TCanvas.h"


using namespace std;

int main(int argc, char** argv){

  if (argc<2) {
    cerr << " give .root or .IEC file name as argument!" << endl;
    return -1;
  }
  TApplication t("app",0,0);
  TFile yc("Calibration.root");
  cout << " using Calibration.root file for calibration." << endl;
  GeCalibrate *cal = (GeCalibrate*) yc.Get("GeCalibrate");
  if (cal !=0) cout << " got calibration" << endl;
  CalMCA_data *cd = cal->ApplyCalibration(TString(argv[1]));
  cd->Analyse();
  t.Run();

}
