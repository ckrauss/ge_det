// Test for isotope data base
//
// Carsten Krauss April 2012
//


#include "iostream"
#include "fstream"

#include "time.h"
#include "TH1I.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TFile.h"
#include "TKey.h"


#include "MCA_data.h"

using namespace std;


int main(int argc, char **argv){

  //  if (argc<2) {
  //    cerr << " give argument!" << endl;
  //    return -1;
  //  }
  if (argc<3) {
    cerr << " give name of two _cal.root files as argument! first is the sample, the second should be the background file to be subtracted. The name of the calibration file, if used, may be supplied as an optional third argument." << endl;
    return -1;
  }
 
  IsotopeDB base;
  base.InitFromFile("isotopedb.txt");
  cout << " Found " << base.GetNumberOfIsotopes() << " isotopes in database"  << endl;

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
    //    yc.ls();
    cal = (GeCalibrate*) yc.Get("GeCalibrate;1");
  }


  a.ActivityDetermination(data1,data2,base,cal);

  t.Run();
  return 0;
}
