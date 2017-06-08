
#include "iostream"
#include "TString.h"
#include "TKey.h"
#include "TFile.h"
#include "MCA_data.h"
#include "TApplication.h"
#include "TCanvas.h"


using namespace std;

int main(int argc, char** argv){

  if (argc<3) {
    cerr << " give calibration.root file and a .root, .IEC or .CNF file name as argument!" << endl;
    return -1;
  }
  //  TApplication t("app",0,0);
  TFile yc(argv[1]);
  cout << " using "<<argv[1]<<" file for calibration." << endl;
  //  yc.ls();
  GeCalibrate *cal = (GeCalibrate*) yc.Get("GeCalibrate;1");
  //cout << " DEBUG "<< endl;
  if (cal != 0) cout << " got calibration" << endl;
  cal->ApplyCalibration(TString(argv[2]));
  //  t.Run();

}
