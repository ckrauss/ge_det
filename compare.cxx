
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
    cerr << " give name of two _cal.root files as argument! first is sample, second should be background to be subtracted. " << endl;
    return -1;
  }
  TApplication t("app",0,0);
  TFile y1(argv[1]);
  TFile y2(argv[2]);
  if (y1.IsZombie() || y2.IsZombie()) return (-1);
  GeAnalyse a;
  
  CalMCA_data* data1 = (CalMCA_data*)((TKey*)y1.GetListOfKeys()->First())->ReadObj();
  CalMCA_data* data2 = (CalMCA_data*)((TKey*)y2.GetListOfKeys()->First())->ReadObj();
  
  a.Compare(data1,data2);

  t.Run();

}
