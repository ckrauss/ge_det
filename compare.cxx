
#include "iostream"
#include "TString.h"
#include "TFile.h"
#include "MCA_data.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TKey.h"


using namespace std;

int main(int argc, char** argv){

  if (argc<3) {
    cerr << " give name of two _cal.root files as argument! first is the sample, the second should be the background file to be subtracted. " << endl;
    return -1;
  }
  TApplication t("app",0,0);
  TFile y1(argv[1]);
  TFile y2(argv[2]);
  if (y1.IsZombie() || y2.IsZombie()) return (-1);
  GeAnalyse a;
  
  CalMCA_data* data1 = (CalMCA_data*)((TKey*)y1.GetListOfKeys()->First())->ReadObj();
  CalMCA_data* data2 = (CalMCA_data*)((TKey*)y2.GetListOfKeys()->First())->ReadObj();
  GeCalibrate *cal = 0;
  if (argc==4){
    TFile yc(argv[3]);
    cout << " using "<<argv[3]<<" file for calibration." << endl;
    yc.ls();
    cal = (GeCalibrate*) yc.Get("GeCalibrate;1");
  }


  a.Compare(data1,data2,cal);

  t.Run();

}
