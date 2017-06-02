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
#include "TFile.h"

#include "MCA_data.h"

using namespace std;


int main(int argc, char **name){

  //  if (argc<2) {
  //    cerr << " give argument!" << endl;
  //    return -1;
  //  }
  IsotopeDB base;
  base.InitFromFile("isotopedb.txt");
  cout << " Found " << base.GetNumberOfIsotopes() << " isotopes in database"  << endl;
  //for (int i = 11 ; i < 12 ; i++){
      for (int i = 0 ; i < base.GetNumberOfIsotopes() ; i++){
    Isotope Is = base.GetIsotopeByNumber(i);
    cout << Is.GetName() << endl;
    double E, I;
    for (int y = 0; y < Is.GetNumberOfLines() ; y++){
      Is.GetEnergyAndIntensityByIntensity(y,E,I);
      cout << " ---> " << I << " " << E << " keV" << endl;
    }

  }

  return 0;
}
