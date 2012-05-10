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
  cout << "Found " << base.GetNumberofIsotopes() << endl;
  return 0;
}
