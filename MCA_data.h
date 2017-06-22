
#ifndef _MCA_data_
#define _MCA_data_

#include "vector"
#include "TH1I.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

using namespace std;

class MCA_data : public TNamed {
 private:
  double dreal_time; // seconds
  double dlive_time; // seconds
  time_t istart_time; // run start time (UNIX time)
  int ichannels; // MCA channels
  TH1I *hdata;
  char *sshortname;//! internal basename placeholder
  char *soutfilename;//! internal name placeholder
 public:
  MCA_data(const char *name);
  MCA_data(){hdata=0;};
  ~MCA_data();
  double GetRealTime(){return dreal_time;};
  double GetLiveTime(){return dlive_time;};
  time_t GetStartTime(){return istart_time;};
  int GetChannels(){return ichannels;};
  TH1I* GetData(){return hdata;};
  TH1D* GetNormalizedData();
  void SetRealTime(double v){dreal_time= v;};
  void SetLiveTime(double v){dlive_time = v;};
  void SetStartTime(time_t t){istart_time = t;};
  void SetChannels(int c){ichannels = c;};
  void SetData(TH1I* h){hdata = h;};
  int Read(char* filename);
  int Read_IEC(char* filename);
  int Read_CNF(char* filename);
  char* GetOutputFilename(){return soutfilename;};
  char* GetBasename(){return sshortname;};
  ClassDef(MCA_data,2);
};

class CalMCA_data : public MCA_data {
 private:
  TH1D* hCalData;
 public:
  TH1D* GetCalData(){return hCalData;};
  void SetCalData(TH1D* cald){hCalData = cald;};
  CalMCA_data(){hCalData = 0;};
  void Analyse();
  CalMCA_data(MCA_data &dat);
  ~CalMCA_data(){};
  ClassDef(CalMCA_data,1);
};

class SourceData : public TObject {
 private:
  vector<double> vdGammaEnergy;
  vector<double> vdLineIntensity;
  vector<double> vdChannelNumber;
  double dSourceStrength;// activity in MBq
  double dHalfLife; // in seconds
  time_t tReferenceDate;
  TString        sIsotopeName;
  TString        sFileName;
  MCA_data       data;
  bool           GetFile();
 public:
  TString GetFilename(){return sFileName;};
  TString GetIsotopeName(){return sIsotopeName;};
  double GetEnergy(int line){return vdGammaEnergy[line];};
  double GetIntensity(int line){return vdLineIntensity[line];};
  int GetNumLines(){return vdGammaEnergy.size();};
  time_t GetReferenceDate(){return tReferenceDate;};
  double GetSourceStrength(){return dSourceStrength;};
  double GetHalfLife(){return dHalfLife;};
  void InitCo60(TString filename, bool old = false);
  void InitCo57(TString filename);
  void InitMn54(TString filename);
  void InitBa133(TString filename, bool old = false);
  void InitCd109(TString filename);
  void InitCs137(TString filename, bool old = false);
  void InitFe55(TString filename);
  void InitSr90(TString filename);
  void InitTl204(TString filename);
  void InitZn65(TString filename);
  void InitAm241(TString filename, bool old = false);
  void InitHg203(TString filename, bool old = false);
  void SetChannelNumber(int line,double channel){vdChannelNumber[line] = channel;};
  double GetChannelNumber(int line){return vdChannelNumber[line];};
  MCA_data *GetData(){return &data;};
  ClassDef(SourceData,4);
};

class GeCalibrate: public TObject {
  vector<SourceData*> vpData;
  TGraphErrors *calplot;
  double dSlope; //calibration constant
  double dOffset;// calibration constant
  bool valid_calibration;
  TString cal_date;//! name tag for output file.
  TGraphErrors *linplot; // Plot to show deviation from linearity in the energy calibration
  TGraphErrors *effplot; // Plot to show  efficiency as function of energy.
  TF1 *fEff;// Efficiency function as determined by calibration with known sources.
  double EffFit(double *x, double *par);
  TCanvas *pefffitplot;
 public:
  GeCalibrate(){valid_calibration=false; calplot=0;linplot=0;effplot=0;fEff=0;pefffitplot=0;};
  ~GeCalibrate(){if (calplot) delete calplot; if (linplot) delete linplot; if (effplot) delete effplot; if (fEff) delete fEff; if (pefffitplot) delete pefffitplot;};
  void Init(vector <TString> filenames);
  void Calibrate();
  CalMCA_data* ApplyCalibration(TString filename);
  double GetEnergyfromChannel(double channel){if (valid_calibration) return dOffset+dSlope*channel; else return -1;};
  double GetChannelfromEnergy(double energy){if(valid_calibration) return (energy-dOffset)/dSlope; else return -1;};
  TString GetCalibrationDate(){return cal_date;};
  double GetEfficiency(double energy){return fEff->Eval(energy);};// energy in keV
  double FindFWHM(TH1I* data, double X, double Y);
  ClassDef(GeCalibrate,6);
};

class IsotopeDB;

class GeAnalyse {
 private:

 public:
  void Compare(CalMCA_data *d1, CalMCA_data* d2, GeCalibrate *cal);
  void BackgroundRemove(CalMCA_data *data, CalMCA_data* background);
  void TimeDependence(int ld_count,CalMCA_data** ld);
  void ActivityDetermination(CalMCA_data *d1, CalMCA_data* d2, IsotopeDB &db, GeCalibrate *cal);
};

// Class for analysis of meterials, similar to source class, but used for different purpose
class Isotope : public TObject {
 private:
  TString tsname;// isotope name
  vector<double> vdenergy;// gamma lines (in keV)
  vector<double> vdintensity;// intensity for each line (in %)
  double dhalflife; // for the isotope in seconds;
 public:
  Isotope(){};
  ~Isotope(){};
  int GetNumberOfLines(){return vdenergy.size();};
  double GetEnergy(int line){return vdenergy[line];};
  double GetIntensity(int line){return vdintensity[line];};
  void GetEnergyAndIntensityByIntensity(int line, double &energy,double &intensity);
  void SetEnergyAndIntensity(double energy, double intensity){vdenergy.push_back(energy); vdintensity.push_back(intensity);};
  void SetHalfLife(double halflife){dhalflife = halflife;};
  double GetHalfLife(){return dhalflife;};
  TString GetName(){return tsname;};
  TString GetLatexName();
  double GetAtomicMass();
  void SetName(TString name){tsname = name;};
  ClassDef(Isotope,1);
};

// class database of all known isotopes to be analysed...
class IsotopeDB : public TObject {
 private:
  vector <Isotope> vIdb;// database (list) of all known and implemented isotopes.
 public:
  IsotopeDB();
  ~IsotopeDB(){};
  void Init();
  void InitFromFile(TString fn);
  Isotope &GetIsotopeByName(TString name);
  Isotope &GetIsotopeByNumber(int n);
  // Isotope &GetIsotopeByEnergy(double e);
  int GetNumberOfIsotopes(){return vIdb.size();};
  ClassDef(IsotopeDB,1);
};

#endif
