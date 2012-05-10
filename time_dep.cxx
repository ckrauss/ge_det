
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
    cerr << " give name of several _cal.root files as argument!" << endl;
    return -1;
  }
  TApplication t("app",0,0);
  GeAnalyse a;

  CalMCA_data *list[argc-1];
  
  for (int i = 1 ; i < argc ; i++){
    TFile y(argv[i]);
    if (y.IsZombie()) {
      cout << " File " << argv[i] << " could not be opened. " << endl;
      return (-1);
    }
    
    list[i-1] = (CalMCA_data*)((TKey*)y.GetListOfKeys()->First())->ReadObj();
  }
  //  TCanvas *c = new TCanvas("c","Time Dependence", 1000,750);
  a.TimeDependence(argc-1,list);

  t.Run();

}
