// started by C. Krauss in 2008
// updated to include calibration and analysis in 2012
//
// Compare spectra
// database of background lines

#include "sstream"
#include "math.h"
#include "TMath.h"
#include "MCA_data.h"
#include "TCanvas.h"
#include "TRegexp.h"
#include "TFile.h"
#include "TLegend.h"
#include "TF1.h"
#include "iostream"
#include "fstream"
#include "iomanip"
#include "TSpectrum.h"
#include "TString.h"
#include "TPolyMarker.h"
#include "TGraphErrors.h"

using namespace std;

// define reference times for source strength
#define JAN2011 1293904800 // (Jan 1 2011, 12:00 MT)
#define NOV2010 1288634400 // (Nov 1 2010, 12:00 MT)
#define APR1974  134049600 // (Apr 1 1974, 12:00 GMT)
#define DAYSINYEAR 365.2425
#define SECONDSINDAY (3600.0*24.0)
#define SECONDSINYEAR (DAYSINYEAR*SECONDSINDAY)

ClassImp(MCA_data);

ClassImp(CalMCA_data);

ClassImp(GeCalibrate);

ClassImp(SourceData);

MCA_data::MCA_data(const char *name):TNamed(name, "MCA data"){ 
  hdata = 0;
};

MCA_data::~MCA_data(){
  //  if (hdata!=0)
  //    delete hdata;
};


int MCA_data::Read(char *filename){
  SetName(filename);
  char* ext = strrchr(filename,'.');
  if (ext !=0){
    if (!strcmp(ext,".CNF")||!strcmp(ext,".cnf"))
      Read_CNF(filename);
    else
      if (!strcmp(ext,".IEC")||!strcmp(ext,".iec"))
	Read_IEC(filename);
      else {
	cout << "File format unknown as input type." << endl;
	return -1;
      }
  }
  else{
    cout << " filename needs to have an extension." << endl;
    return -1;
  }
  return 0;
}

int MCA_data::Read_CNF(char *filename){
  static int instance = 0;
  instance++;
  //  char basename[300];
  //  char outfilename[300];
  soutfilename = new char[300];
  sshortname = new char[300];
  strncpy(sshortname,filename,300);
  char *ext = strrchr(sshortname,'.');
  if (!strcmp(ext,".CNF")||!strcmp(ext,".cnf")){
    ext[0] = 0;
    cout<< " opening CNF file: "<<sshortname<<endl;
    strcpy(soutfilename,sshortname);
    strncpy(soutfilename+strlen(soutfilename),".root",10);
  }
  if (strrchr(sshortname,'/')){
    strncpy(sshortname,strrchr(sshortname,'/')+1,strlen(sshortname)-(sshortname-strrchr(sshortname,'/')));
  }
  SetName(sshortname);
  ifstream inf(filename,ios::binary);
  unsigned short header[8];
  inf.read((char*)header,16);
  //  header[3] = (header[3] & 0xF) >> 4 |(header[3]&0xF0)>>4| (header[3]&0xF00)<<4 | (header[3]&0xF000)>>4;
  cout << " 1: " << header[0] << " 2: "  << header[1] << " 3: " << header[2] << endl;
  cout << " 4: " << header[3] << " 5: " << header[4] << " 6: " << header[5] << endl;
  //  int channels = header[3];
  int channels = 0;
  if (header[5] == 38400){
    cout << " 32768 channels in file !?" << endl;
    channels = 32768;
  }
  if (header[5] == 37888){
    cout << " 32768 channels in file !?" << endl;
    channels = 32768;
  }
  if (header[5] == 39936){
    cout << " 32768 channels in file !?" << endl;
    channels = 32768;
  }
  if (header[5] == 5632){
    cout << " 8192 channels in file !?" << endl;
    channels = 32768/4.;// 8k file
    header[5] = 38400;
  }    
  if (channels == 0){
    cerr << " unkonwn file type! aborting." << endl;
    exit(0);
  }
  // if (header[1] == 4 && header[2] == 1028) {
  //   channels = 32768/4.;
  //   header[5] = 38400;
  // }
  // if (header[1] == 0 && header[2] == 1024) {
  //   channels = 32768;
  // }
  cout << " number of channels: " << channels << " " << "  data offset: " << header[5]<< endl;
  //  int channel = 32768;
  int date_offs = 0xb07;  // Current version of CNF file has date and time info start here, time is consecutive
  int realt_offs = 0xb0f;
  int livet_offs = 0xb17;
  inf.seekg(date_offs);
  unsigned long long date, livet, realt;
  inf.read((char*)&date,8);
  double ddate = double(date)/1e7;
  double realdate = ddate - 3506691600; // date in CNF file is in seconds since 1858, Nov 16, 7:00 GMT, the offset to unix time is given here.  (time in Mountain time)
  //  realdate -= 2* 3600; // convert to eastern time
  time_t real_date_t = (time_t(realdate));
  SetStartTime(real_date_t);
  char *time_c = ctime(&real_date_t);
  cout << " converted, extracted time: " <<   time_c;
  char *tmp_c = strchr(time_c,'\n');
  tmp_c[0]= '\0';
  inf.read((char*)&realt,8);
  inf.read((char*)&livet,8);
  double drealt = double(~(realt)+1)/1e7; // Calculate 2's complement and convert to seconds...
  SetRealTime(drealt);
  double dlivet = double(~(livet)+1)/1e7;
  SetLiveTime(dlivet);
  SetChannels(channels);
  cout << " data taking time (real): " << setprecision(10)<< drealt << "s, (live): " << dlivet << "s" << endl;
  inf.seekg(header[5]);
  unsigned int mydata[channels];
  inf.read((char*)mydata,channels*4);
  char title[400];
  sprintf(title,"MCA Data, %s, %7.2lf s live time",time_c, GetLiveTime());
  char hn[400];
  sprintf(hn,"data_%d",instance-1);
  SetData(new TH1I(hn,title, GetChannels(),0,GetChannels()+1));
  int entries = 0;
  for (int i = 0 ; i<channels ; i++){
    GetData()->SetBinContent(1+i,mydata[i]);
    entries += mydata[i];
  }
  GetData()->SetEntries(entries);
  inf.close();
  //  TFile yo(outfilename,"recreate");
  //  h->Write();
  //  yo.Close();
  cout << " Done, success." << endl;
  return 1;
}


int MCA_data::Read_IEC(char *filename){
  static int instance = 0;
  instance++;
  char buffer[1000];
  ifstream inf(filename);
  if (!inf.is_open()) {
    cout << " file " << filename << " not found " << endl;
    return 0;
  }
  //  char outname[300],shortname[300];
  soutfilename = new char[300];
  sshortname = new char[300];
  if (strrchr(filename,'.')!=0){
    strncpy(sshortname,filename,strrchr(filename,'.')-filename);
    sshortname[strrchr(filename,'.')-filename] = '\0';
    //    cout << "DEBUG: " << sshortname << " " << filename<< " " << strchr(filename,'.')-filename<< endl;
  }
  else
    strncpy(sshortname,filename,strlen(filename));
  //  cout << " DEBUG: " << sshortname << endl;
  sprintf(soutfilename,"%s.root",sshortname);
  //  TFile *yo = new TFile(soutfilename,"recreate");
  // remove leading path
  if (strrchr(sshortname,'/')){
    strncpy(sshortname,strrchr(sshortname,'/')+1,strlen(sshortname)-(sshortname-strrchr(sshortname,'/')));
  }
  SetName(sshortname);
  //  MCA_data mydata(sshortname);
  // first line contains system information, 
  inf.getline(buffer,999);
  // second line contains live time, real time and channel number
  inf.getline(buffer,999);
  char tmp[50], *tmp2;
  strncpy(tmp,buffer+4,14);
  tmp[14] = '\0';
  SetLiveTime(strtod(tmp,&tmp2));
  strncpy(tmp,buffer+18,14);
  tmp[14] = '\0';
  SetRealTime(strtod(tmp,&tmp2));
  strncpy(tmp,buffer+32,strlen(buffer)-32);
  SetChannels(atoi(tmp));
  // third line contains start time and date
  inf.getline(buffer,999);
  struct tm tdate;
  tdate.tm_mday = atoi(buffer+4); // day
  tdate.tm_mon = atoi(buffer+7)-1; // month
  tdate.tm_year = atoi(buffer+10)+100; // year
  tdate.tm_hour = atoi(buffer+13); //hour
  tdate.tm_min = atoi(buffer+16); // minutes
  tdate.tm_sec = atoi(buffer+19); //seconds
  SetStartTime(mktime(&tdate));
  time_t tmp_t = GetStartTime();
  char *time_c = ctime(&tmp_t);
  char *tmp_c = strchr(time_c,'\n');
  tmp_c[0]= '\0';
  cout << " got time: " << time_c << endl;
  // fourth line contains energy calibration (unused here)
  inf.getline(buffer,999);

  cout << " got channels: " << GetChannels() << endl;
  sprintf(buffer,"MCA Data, %s, %7.2lf s live time",time_c, GetLiveTime());
  char hn[400];
  sprintf(hn,"data_%d",instance-1);
  SetData(new TH1I(hn,buffer, GetChannels(),0,GetChannels()+1));

  bool found = false;
  while (!inf.eof()){
    inf.getline(buffer,999);
    //    cout << buffer << endl;
    if (!strncmp(buffer,"A004USERDEFINED",15)) {
      found = true;
      break;
    }
  }
  if (found){
    // data entries start here
    int data[5];
    int channel = 0;
    int entries = 0;
    while (!inf.eof() && channel+5<GetChannels()){
      inf.getline(buffer,999);
      //      cout << buffer << endl;
      sscanf(buffer+4,"%d %d %d %d %d %d",&channel,&data[0],&data[1],&data[2],&data[3],&data[4]);
      GetData()->SetBinContent(channel+1,data[0]);
      GetData()->SetBinContent(channel+2,data[1]);
      GetData()->SetBinContent(channel+3,data[2]);
      GetData()->SetBinContent(channel+4,data[3]);
      GetData()->SetBinContent(channel+5,data[4]);
      entries += data[0] + data[1] + data[2] + data[3] + data[4];
      if (inf.eof()) break;
    }
    GetData()->SetEntries(entries);
    //    mydata.Write();
    cout << " done, success." << endl;
  }
  else{
    cout << " data elements not found, no output." << endl;
    return 0;
  }
  //  yo->Close();


  return 1;
}

TH1D* MCA_data::GetNormalizedData(){
  hdata->Sumw2();
  TH1D* rv = new TH1D("data_norm",hdata->GetTitle(),hdata->GetNbinsX(),hdata->GetXaxis()->GetXmin(),hdata->GetXaxis()->GetXmax());
  for (int i = 0 ; i < hdata->GetNbinsX()+2 ; i++){ // include under- and overflow bins.
    rv->SetBinContent(i,hdata->GetBinContent(i)/dreal_time);
    rv->SetBinError(i,hdata->GetBinError(i)/dreal_time);
  }
  rv->SetEntries(double(hdata->GetEntries()));
  return rv;
}

void SourceData::InitCo60(TString filename, bool old){
  // Data from http://nucleardata.nuclear.lu.se/nucleardata/toi/nuclide.asp?iZA=270060
  //Eg (keV)    	Ig (%)    	Decay mode
  //346.93 7 	0.0076 5 	b- 
  //826.06 3 	0.0076 8 	b- 
  //1173.237 4 	99.9736 7 	b- 
  //1332.501 5 	99.9856 4 	b- 
  //2158.57 10 	0.00111 18 	b- 
  //2505 	2.0E-6 4 	b- 
  vdGammaEnergy.push_back(346.93);
  vdLineIntensity.push_back(0.0076);
  vdGammaEnergy.push_back(826.06);
  vdLineIntensity.push_back(0.0076);
  vdGammaEnergy.push_back(1173.237);
  vdLineIntensity.push_back(99.9736);
  vdGammaEnergy.push_back(1332.501);
  vdLineIntensity.push_back(99.9856);
  vdGammaEnergy.push_back(2158.57);
  vdLineIntensity.push_back(0.00111);
  vdGammaEnergy.push_back(3505);
  vdLineIntensity.push_back(2.0E-6);
  if (old){
    dSourceStrength = 0.39701; // source stength in MBq (10.73uCi)
    tReferenceDate = APR1974;
  }
  else{
    dSourceStrength = 0.037; // source stength in MBq
    tReferenceDate = JAN2011;
  }
  dHalfLife = 1925.28*SECONDSINDAY; // in seconds, 1925.28d (NNDC)
  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Co60");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitAm241(TString filename, bool old){
  // Data from http://nucleardata.nuclear.lu.se/nucleardata/toi/nuclide.asp?iZA=95024
  // list truncated, only 60keV line implemented...
  //Eg (keV)    	Ig (%)    	Decay mode
  //13.81 2 	 	a 
  //26.3448 2 	2.40 2 	a 
  //27.03 	 	a 
  //31.4 	 	a 
  //32.183 	0.0174 4 	a 
  //33.1964 3 	0.126 3 	a 
  //38.54 3 	 	a 
  //42.73 5 	0.0055 11 	a 
  //43.423 10 	0.073 8 	a 
  //51.01 3 	0.000026 12 	a 
  //54.0 	 	a 
  //55.56 2 	0.0181 18 	a 
  //56.8 	 	a 
  //57.85 5 	0.0052 15 	a 
  //59.5412 2 	35.9 4 	a 
  //61.46 	 	a 
  //64.83 2 	0.000145 18 	a 
  //67.45 5 	0.00042 10 	a 
  //69.76 3 	0.0029 4 	a 
  //75.8 2 	~0.0006 	a 
  //78.1 	 	a  

  vdGammaEnergy.push_back(59.5412);
  vdLineIntensity.push_back(39.9);
  if (old){
    dSourceStrength = 0.43216; // source stength in MBq (11.68uCi)
    tReferenceDate = APR1974;
  }
  else {
    cerr << " new Am241 Source does not exist in our inventory." << endl;
    exit(0);
  }
  dHalfLife = 432.6*SECONDSINYEAR; // in seconds, 432.6y (NNDC)
  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Am241");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitCo57(TString filename){
  //Eg (keV)    	Ig (%)    	Decay mode
  //14.41300 15 	9.16 15 	e 
  //122.0614 4 	85.60 17 	e 
  //136.4743 5 	10.68 8 	e 
  //230.29 2 	0.0004 4 	e 
  //339.54 18 	0.0139 3 	e 
  //352.36 1 	0.0132 3 	e 
  //366.75 1 	0.0013 3 	e 
  //569.92 4 	0.017 1 	e 
  //692.03 2 	0.157 9 	e 
  //706.40 20 	0.0253 5 	e 
  vdGammaEnergy.push_back(14.413);
  vdLineIntensity.push_back(9.16);
  vdGammaEnergy.push_back(122.0614);
  vdLineIntensity.push_back(85.6);
  vdGammaEnergy.push_back(136.4743);
  vdLineIntensity.push_back(10.68);
  vdGammaEnergy.push_back(230.29);
  vdLineIntensity.push_back(0.0004);
  vdGammaEnergy.push_back(339.54);
  vdLineIntensity.push_back(0.0139);
  vdGammaEnergy.push_back(352.36);
  vdLineIntensity.push_back(0.0132);
  vdGammaEnergy.push_back(366.75);
  vdLineIntensity.push_back(0.0013);
  vdGammaEnergy.push_back(569.92);
  vdLineIntensity.push_back(0.017);
  vdGammaEnergy.push_back(692.03);
  vdLineIntensity.push_back(0.157);
  vdGammaEnergy.push_back(706.40);
  vdLineIntensity.push_back(0.0253);
  dSourceStrength = 0.037; // source stength in MBq
  tReferenceDate = JAN2011;
  dHalfLife = 271.74*SECONDSINDAY; // in seconds, 271.74d (NNDC)


  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Co57");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitBa133(TString filename, bool old){
  // Data from http://nucleardata.nuclear.lu.se/nucleardata/toi/nuclide.asp?iZA=560133
// Eg (keV)    	Ig (%)    	Decay mode
// 53.161 1 	2.199 22 	e 
// 79.6139 26 	2.62 6 	e 
// 80.9971 14 	34.06 27 	e 
// 160.613 8 	0.645 8 	e 
// 223.234 12 	0.450 4 	e 
// 276.398 2 	7.164 22 	e 
// 302.853 1 	18.33 6 	e 
// 356.017 2 	62.05 19 	e 
// 383.851 3 	8.94 3 	e   
  vdGammaEnergy.push_back(53.161);
  vdLineIntensity.push_back(2.199);
  vdGammaEnergy.push_back(79.6139);
  vdLineIntensity.push_back(2.62);
  vdGammaEnergy.push_back(80.9971);
  vdLineIntensity.push_back(34.06);
  vdGammaEnergy.push_back(160.613);
  vdLineIntensity.push_back(0.645);
  vdGammaEnergy.push_back(223.234);
  vdLineIntensity.push_back(0.450);
  vdGammaEnergy.push_back(276.394);
  vdLineIntensity.push_back(7.164);
  vdGammaEnergy.push_back(302.853);
  vdLineIntensity.push_back(18.33);
  vdGammaEnergy.push_back(356.017);
  vdLineIntensity.push_back(62.05);
  vdGammaEnergy.push_back(383.851);
  vdLineIntensity.push_back(8.94);
  if (old){
    tReferenceDate = APR1974;
    dSourceStrength = 0.39627; // source stength in MBq (10.71uCi)
  }
  else{
    tReferenceDate = JAN2011;
    dSourceStrength = 0.037; // source stength in MBq
  }
  dHalfLife = 10.551*SECONDSINYEAR; // in seconds, 10.551y (NNDC)
  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Ba133");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

bool SourceData::GetFile(){
  if (sFileName.Contains(".root")){
    TFile y(sFileName);
    if (y.IsZombie()) cerr<< "File is zombie!" << endl;
    sFileName.Resize(sFileName.Last('.'));
    if (sFileName.Last('/')!=0)
      sFileName.Replace(0,1+sFileName.Last('/'),0,0);
    //  cout << sFileName << endl;
    MCA_data *pdata = (MCA_data*) y.Get(sFileName);
    if (pdata != 0){
      data = *pdata;
      return true;
    }
  }
  else{
    if (sFileName.Contains(".IEC") || sFileName.Contains(".iec")){
      MCA_data *pdata = new MCA_data();
      char *fn = new char[300];
      strncpy(fn,sFileName.Data(),300);
      if (pdata->Read_IEC(fn)!=0){
	data = *pdata;//small memory leak here. ignore.
	return true;
      }
    }
    if (sFileName.Contains(".CNF") || sFileName.Contains(".cnf")){
      MCA_data *pdata = new MCA_data();
      char *fn = new char[300];
      strncpy(fn,sFileName.Data(),300);
      if (pdata->Read_CNF(fn)!=0){
	data = *pdata;//small memory leak here. ignore.
	return true;
      }
    }
  }
  return false;
}

void SourceData::InitMn54(TString filename){
  // Data from http://nucleardata.nuclear.lu.se/nucleardata/toi/nuclide.asp?iZA=250054
  //Eg (keV)    	Ig (%)    	Decay mode
  //834.848 3 	99.976 1 	e+b+ 
  vdGammaEnergy.push_back(834.848);
  vdLineIntensity.push_back(99.976);
  dSourceStrength = 0.037; // source stength in MBq
  tReferenceDate = JAN2011;
  dHalfLife = 312.12*SECONDSINDAY; // in seconds, 312.12d (NNDC)

  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Mn54");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitCs137(TString filename, bool old){
  //Eg (keV)    	Ig (%)    	Decay mode
  //283.53 4 	0.00058 8 	b- 
  //661.657 3 	85.1 2 	b- 
  vdGammaEnergy.push_back(283.53);
  vdLineIntensity.push_back(0.00058);
  vdGammaEnergy.push_back(661.657);
  vdLineIntensity.push_back(85.1);
  if (old){
    dSourceStrength = 0.41033; // source stength in MBq (11.09uCi)
    tReferenceDate = APR1974;
  }
  else{
    dSourceStrength = 0.037; // source stength in MBq
    tReferenceDate = JAN2011;
  }
  dHalfLife = 30.8*SECONDSINYEAR; // in seconds, 30.8y (NNDC)

  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Cs137");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitCd109(TString filename){
  //Eg (keV)    	Ig (%)    	Decay mode
  //88.04 5 	3.61 10 	e 
  vdGammaEnergy.push_back(88.04);
  vdLineIntensity.push_back(3.61);
  dSourceStrength = 0.037; // source stength in MBq
  tReferenceDate = JAN2011;
  dHalfLife = 461.4*SECONDSINDAY; // in seconds, 461.4 days (NNDC)

  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Cd109");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}



void SourceData::InitSr90(TString filename){
  cerr << " Not implemented. " << endl;

}

void SourceData::InitTl204(TString filename){
  //XR l	     9.99	      0.81 % 3 	  8.1E-5 3 
  //XR kα2	    68.894	      0.469 % 15 	  3.23E-4 10 
  //XR kα1	    70.818	      0.789 % 24 	  5.59E-4 17 
  //XR kβ3	    79.824	      0.095 % 3 	  7.60E-5 24 
  //XR kβ1	    80.225	      0.182 % 6 	  1.46E-4 5 
  //XR kβ2	    82.473	      0.0659 % 21 	  5.44E-5 17 
  // from NNDC database....
  vdGammaEnergy.push_back(9.99);
  vdLineIntensity.push_back(0.810);
  vdGammaEnergy.push_back(68.894);
  vdLineIntensity.push_back(0.469);
  vdGammaEnergy.push_back(70.818);
  vdLineIntensity.push_back(0.789);
  vdGammaEnergy.push_back(79.824);
  vdLineIntensity.push_back(0.095);
  vdGammaEnergy.push_back(80.225);
  vdLineIntensity.push_back(0.182);
  vdGammaEnergy.push_back(82.473);
  vdLineIntensity.push_back(0.0659);
  dSourceStrength = 0.037; // source stength in MBq
  tReferenceDate = NOV2010;
  dHalfLife = 3.783*365.2425*SECONDSINDAY; // in seconds, 3.783y (NNDC)

  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Tl204");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitZn65(TString filename){
  //344.95 20 	0.0030 3 	e+b+ 
  //770.6 2 	0.0030 3 	e+b+ 
  //1115.546 4 	50.60 24 	e+b+ 
  vdGammaEnergy.push_back(344.95);
  vdLineIntensity.push_back(0.0030);
  vdGammaEnergy.push_back(770.6);
  vdLineIntensity.push_back(0.0030);
  vdGammaEnergy.push_back(1115.546);
  vdLineIntensity.push_back(50.60);
  dSourceStrength = 0.037; // source stength in MBq
  tReferenceDate = JAN2011;
  dHalfLife = 243.93*SECONDSINDAY; // in seconds, 243.93d (NNDC)

  for (int i = 0 ; i<vdGammaEnergy.size(); i++)
    vdChannelNumber.push_back(0);
  sFileName = filename;
  sIsotopeName = TString("Zn65");
  if (!GetFile() ) cerr << " File not opened. Error" <<endl;
}

void SourceData::InitFe55(TString filename){
  cerr << " Not implemented. " << endl;

}

void GeCalibrate::Init(vector<TString> filenames){
  for (vector<TString>::iterator i=filenames.begin() ; i!=filenames.end();i++){
    cout << "File: " << *i << endl;
    
    // file must start with isotope name:
    if ((*i).Contains("Co-60")){
      cout << "found Co60 file" << endl;
      SourceData *d = new SourceData();
      if ((*i).Contains("old"))
	d->InitCo60(*i, true);
      else
	d->InitCo60(*i);
      vpData.push_back(d);
      char tmp[60];
      strcpy(tmp,(*i).Data()+6);
      char* t2 = strrchr(tmp,'.');
      t2[0] = 0;
      cal_date = TString(tmp);
    }
    if ((*i).Contains("Am-241")){
      cout << "found Am241 file" << endl;
      SourceData *d = new SourceData();
      if ((*i).Contains("old"))
	d->InitAm241(*i, true);
      else
	d->InitAm241(*i);
      vpData.push_back(d);
      char tmp[60];
      strcpy(tmp,(*i).Data()+6);
      char* t2 = strrchr(tmp,'.');
      t2[0] = 0;
      cal_date = TString(tmp);
    }
    if ((*i).Contains("Mn-54")){
      cout << "found Mn54 file" << endl;
      SourceData *d = new SourceData();
      d->InitMn54(*i);
      vpData.push_back(d);
      char tmp[60];
      strcpy(tmp,(*i).Data()+6);
      char* t2 = strrchr(tmp,'.');
      t2[0] = 0;
      cal_date = TString(tmp);
    }
    if ((*i).Contains("Ba-133")){
      cout << "found Ba133 file" << endl;
      SourceData *d = new SourceData();
      if ((*i).Contains("old"))
	d->InitBa133(*i, true);
      else
	d->InitBa133(*i);
      vpData.push_back(d);
      char tmp[60];
      strcpy(tmp,(*i).Data()+7);
      char* t2 = strrchr(tmp,'.');
      t2[0] = 0;
      cal_date = TString(tmp);
    }
    if ((*i).Contains("Zn-65")){
      cout << "found Zn65 file" << endl;
      SourceData *d = new SourceData();
      d->InitZn65(*i);
      vpData.push_back(d);
      char tmp[60];
      strcpy(tmp,(*i).Data()+6);
      char* t2 = strrchr(tmp,'.');
      t2[0] = 0;
      cal_date = TString(tmp);
    }
    if ((*i).Contains("Cd-109")){
      cout << "found Cd109 file" << endl;
      SourceData *d = new SourceData();
      d->InitCd109(*i);
      vpData.push_back(d);
    }
    if ((*i).Contains("Co-57")){
      cout << "found Co-57 file" << endl;
      SourceData *d = new SourceData();
      d->InitCo57(*i);
      vpData.push_back(d);
    }
    if ((*i).Contains("Cs-137")){
      cout << "found Cs-137 file" << endl;
      SourceData *d = new SourceData();
      if ((*i).Contains("old"))
	d->InitCs137(*i, true);
      else
	d->InitCs137(*i);
      vpData.push_back(d);
    }
    if ((*i).Contains("Tl-204")){
      cout << "found Tl-204 file" << endl;
      SourceData *d = new SourceData();
      d->InitTl204(*i);
      vpData.push_back(d);
    }
    if ((*i).Contains("Fe-55")){
      cout << "found Fe-55 file" << endl;
      SourceData *d = new SourceData();
      d->InitFe55(*i);
      vpData.push_back(d);
    }
  }
  if (vpData.size()==0){
    cerr<< " Error: no source run data found. Aborting." << endl;
    exit(0);
  }
  if (vpData.size()<filenames.size()){
    cerr<< " Error: some source run data not found. Aborting. Check file name format (has to start with NN-xxx, e.g. Mn-54 or Ba-133)! " << endl;
    exit(0);
  }
}

CalMCA_data::CalMCA_data(MCA_data &dat) : MCA_data(dat) {
  hCalData=0;
}

void CalMCA_data::Analyse(){
  // ...
  TCanvas *cc = new TCanvas("cc","Spectra");
  cc->Divide(1,2);
  cc->cd(1);
  hCalData->Draw();
  TSpectrum *s = new TSpectrum(160);
  TH1D *b = (TH1D*) s->Background(hCalData,20,"Compton");
  TH1D* d = (TH1D*)hCalData->Clone("diff");
  d->Add(b,-1);
  d->Rebin(16);
  s->Search(d,10);
  s->Search(hCalData,10);
  hCalData->Draw();
  b->Draw("same");
  gPad->SetLogy(1);
  cc->cd(2);
  d->Draw();
  TList *functions = hCalData->GetListOfFunctions();
  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
  double *X;
  double *Y;
  X= pm->GetX();
  Y= pm->GetY();
  for (int i = 0 ; i<pm->GetN() ; i++){
    cout << i << " " << X[i] << " " << Y[i] << endl;
  }
  
}

CalMCA_data* GeCalibrate::ApplyCalibration(TString filename){
  // Calculate livetime corrected and calibrated histogram from raw (root) datafile.

  if (!valid_calibration) {
    cerr << "calibration in file is not valid!" << endl;
    exit(0);
  }
  MCA_data *pdata = 0;
  if (filename.Contains(".root")){
    TFile y(filename);
    if (y.IsZombie()) {
      cerr<< "File is zombie!" << endl;
      exit(0);
    }
    filename.Resize(filename.Last('.'));
    if (filename.Last('/')!=0)
      filename.Replace(0,1+filename.Last('/'),0,0);
    //  cout << sFileName << endl;
    pdata = (MCA_data*) y.Get(filename);
  }
  else{
    if (filename.Contains(".IEC") || filename.Contains(".iec")){
      pdata = new MCA_data(filename.Data());
      char *fn = new char[200];
      strncpy(fn,filename.Data(),filename.Length());
      if (pdata->Read_IEC(fn)==0){
	cerr << " File " << filename << " could not be read."<< endl;
	exit(0);
      }
    }
    if (filename.Contains(".CNF") || filename.Contains(".cnf")){
      pdata = new MCA_data(filename.Data());
      char *fn = new char[200];
      strncpy(fn,filename.Data(),filename.Length());
      if (pdata->Read_CNF(fn)==0){
	cerr << " File " << filename << " could not be read."<< endl;
	exit(0);
      }
    }
  }
  CalMCA_data *caldat = 0;
  if (pdata != 0){
    // got the raw data...
    TH1D *hist = pdata->GetNormalizedData();
    // time corrected data
    caldat = new CalMCA_data(*pdata);
    int bins = hist->GetNbinsX();
    double bmin = GetEnergyfromChannel(hist->GetXaxis()->GetXmin());
    double bmax = GetEnergyfromChannel(hist->GetXaxis()->GetXmax());
    if (bmin==bmax) return 0;
    caldat->SetCalData(new TH1D("cal_data",hist->GetTitle(),bins,bmin,bmax));
    for (int i = 0; i < bins+2 ; i++){
      caldat->GetCalData()->SetBinContent(i,hist->GetBinContent(i));
      caldat->GetCalData()->SetBinError(i,hist->GetBinError(i));
    }
    caldat->GetCalData()->SetEntries(hist->GetEntries());
    TString outfilename(filename+TString("_cal.root"));
    cout << " Trying to write file to : " << outfilename << endl;
    TFile yo(outfilename,"recreate");
    caldat->Write();
    yo.Close();
  }
  
  return caldat;
}

void GeCalibrate::Calibrate(){ 	
  int num_peaks[vpData.size()];
  int alllines = 0;
  TCanvas *c[vpData.size()];
  for (int source = 0 ; source < vpData.size() ; source++){
    char idn[20],name[300];
    sprintf(idn,"can_%d",source);
    sprintf(name,"Spectrum %s",vpData[source]->GetIsotopeName().Data());
    c[source] = new TCanvas(idn,name);
    TSpectrum ts(vpData[source]->GetNumLines(),2);
    ts.SetResolution(6);
    int peaks = ts.Search(vpData[source]->GetData()->GetData(),6,"",0.01);
    vpData[source]->GetData()->GetData()->Draw();
    cout << peaks << " peaks found for source " << vpData[source]->GetIsotopeName()<< endl;
    num_peaks[source] = peaks;
    alllines += num_peaks[source];
  }
  cout << " using " << alllines << " lines for calibration " << endl;
  int min_peaks = 1000;
  int i_min = -1;
  for (int i = 0 ; i<vpData.size();i++){
    if (num_peaks[i]<min_peaks) {min_peaks = num_peaks[i];i_min = i;}
  }
  SourceData *firstI = vpData[i_min];// get isotope data
  cout << "Starting calibration with " << firstI->GetIsotopeName() << endl;
  int max_intensity=0;
  int max_line = -1;
  for (int ii = 0; ii < firstI->GetNumLines() ;ii++){//loop over all lines to find larges intensity
    if (firstI->GetIntensity(ii)>max_intensity){
      max_intensity = firstI->GetIntensity(ii);
      max_line = ii;
    }
  }
  MCA_data* firstD = firstI->GetData();// get measured spectrum
  TH1I *firstH = firstD->GetData(); // get histogram
  firstH->Draw();
  TList *functions = firstH->GetListOfFunctions();
  TPolyMarker *pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
  double *X;
  double *Y;
  X= pm->GetX();
  Y= pm->GetY();
  int dmax_intensity = 0;
  int max_marker = -1;
  for (int i = 0 ; i < pm->GetN(); i++){
    cout << i << " " << X[i] << " " << Y[i] << endl;
    if (Y[i]>dmax_intensity) {
      dmax_intensity = Y[i];
      max_marker = i;
    }
  }
  cout << " marker " << max_marker << " " << max_line << endl;
  if (max_marker>=0 && max_line >= 0){// found match..
    if (firstI->GetChannelNumber(max_line)>0) cout << "DOUBLE ASSIGNMENT OF LINE. Fix me."<< endl;
    else
      firstI->SetChannelNumber(max_line,X[max_marker]);
  }
  // now calculate two point estimate for calibration trough 0,0
  double initialoffset = 0;
  double initialslope = firstI->GetEnergy(max_line)/firstI->GetChannelNumber(max_line);
  cout << "initial guess: E = 0 + " << initialslope << " * ch " << endl;
  int line = 0;
  double energies[alllines], intensities[alllines];
  double channels[alllines], peakh[alllines];
  double ch_e[alllines], rel_int[alllines];
  double eff[alllines];
  //  energies[line] = firstI->GetEnergy(max_line);
  //  channels[line] = firstI->GetChannelNumber(max_line);
  //  ch_e[line] = sqrt(channels[line]);
  //  line++;
  // now process all other lines, but only when the points are close to the initial calibration line...
  for (int isotope = 0; isotope<vpData.size();isotope++){
    SourceData* source = vpData[isotope];
    cout << " *****************************************************************"<< endl;
    cout << source->GetIsotopeName() << endl;
    MCA_data* dat = source ->GetData();// get measured spectrum
    TH1I *hist = dat->GetData(); // get histogram
    hist->Draw();
    functions = hist->GetListOfFunctions();
    pm = (TPolyMarker*)functions->FindObject("TPolyMarker");
    X= pm->GetX();
    Y= pm->GetY();
    bool used[pm->GetN()];
    for (int i = 0; i < pm->GetN() ; i++){ used[i] = false;}// initialize array
    for (int iline = 0 ; iline < source->GetNumLines() ; iline++){
      if (source->GetEnergy(iline) < 20) continue;
      cout << " processing line " << iline <<" " << source->GetEnergy(iline)<<"keV Intensity: "<< source->GetIntensity(iline)<< endl;
      // calculate the channel range for this line:
      double ch = source->GetEnergy(iline) * 1/initialslope;
      double min_ch = ch - (12 + 18*ch/1000);
      double max_ch = ch + (12 + 18*ch/1000);
      cout << " Search window: " << min_ch << " " << max_ch << endl;
      if (iline==0){
	char newtitle[300];
	sprintf (newtitle,"%s(%d) %s", source->GetIsotopeName().Data(),source->GetNumLines(),hist->GetTitle());
	hist->SetTitle(newtitle);
      }
      //   if (isotope == i_min && iline == max_line) continue;//already used this line for initial guess
      double intensity = 0;
      bool match = false;
      for (int i = 0 ; i < pm->GetN(); i++){
	if (X[i]<max_ch && X[i]>min_ch) {
	  intensity = Y[i];
	  // match
	  match = true;
	  if (used[i] == true){
	    cout << "DOUBLE ASSIGNMENT OF MARKER. FIX ME." << endl;
	    continue;
	  }
	  used[i] = true;
	  energies[line] = source->GetEnergy(iline);
	  intensities[line] = source->GetIntensity(iline);
	  double fit_min = X[i] - (12 + 18*ch/10000);
	  double fit_max = X[i] + (12 + 18*ch/10000);
	  double fwhm = FindFWHM(hist, X[i], Y[i]);
	  double area = fwhm/2.35482 * Y[i] * sqrt(2.0*TMath::Pi());
	  double delta_t = dat->GetStartTime() - source->GetReferenceDate();
	    double exp_rate = dat->GetLiveTime()*source->GetSourceStrength()*1e6*exp(-delta_t*log(2)/source->GetHalfLife())*source->GetIntensity(iline);
	  double efficiency = area/exp_rate;
	  eff[line] = efficiency;
	  cout << "&&& *** $$$ fwhm: " << fwhm << " area: " << area<<" exp: " <<exp_rate<< " " << efficiency<< endl;
	  channels[line] = X[i];
	  if (source->GetChannelNumber(iline)>0)
	    cout << " DOUBLE ASSIGNMENT OF LINE " << source->GetIsotopeName() << " " 
		 << source->GetEnergy(iline) << " keV (#"<< iline << "). Fix me." << endl; 
	  else
	    source->SetChannelNumber(iline,X[i]);
	  peakh[line] = Y[i];
	  rel_int[line] = source->GetIntensity(iline)/double(Y[i]);
	  hist->Fit("gaus","Q+","",fit_min,fit_max);//X[i]-15,X[i]+15);// fit gaussian in narrow window around peak
	  TF1 *ff = hist->GetFunction("gaus");
	  ch_e[line] = ff->GetParError(2);//(X[i]-min_ch)+(max_ch-X[i])/2.;
	  //	  channels[line] = ff->GetParameter(1);
	  line++;
	  cout << "Match! " << i << " ch: " << X[i] << " cts: " << Y[i] << " min ch: " << min_ch << " max ch: " << max_ch<< endl;
	  if ((fabs(X[i]-ff->GetParameter(1))/X[i])>.1)
	    cout << " LARGE difference between gaussian and tspectrum peak pos: " << X[i] << " " << ff->GetParameter(1) << endl;
	  break;
	}
      }
      if (match = false) {
	cout << " No match for line " << source->GetEnergy(iline) << " keV (" << source->GetIntensity(iline) << "%) found"<< endl;
      }
    }
  }
  cout <<  " found " << line << " matches" << endl;
  // Add double match removal here:
  
  //
  calplot = new TGraphErrors (line,channels,energies,ch_e,0);
  TCanvas *cc = new TCanvas("Cal","Calibration");
  cc->Divide(1,3);
  cc->cd(1);
  calplot->Draw("a*");
  calplot->GetXaxis()->SetTitle("Channel");
  calplot->GetYaxis()->SetTitle("Energy [keV]");
  calplot->Fit("pol1");
  TF1* f = (TF1*) calplot->GetFunction("pol1");
  if (f!=0){
    dSlope = f->GetParameter(1);
    dOffset = f->GetParameter(0);
    cout << " final function: " << dOffset << " + " << dSlope << " * ch" << endl;
    valid_calibration = true;
  } 
  // Add efficiency curve here!
  TGraph *gi = new TGraph(line,intensities,peakh);
  cc->cd(2);
  gi->Draw("a*");
  TGraph *gi2 = new TGraph(line,intensities,rel_int);
  cc->cd(3);
  gi2->Draw("a*");
  // energy linearity plot
  double diff[line],d_e[line];
  for (int i = 0 ; i < line ; i++){
    cout << " data: " << i << " " << channels[i] << " " << f->Eval(channels[i]) <<" " << energies[i]<<  endl;
    diff[i] = 100*(f->Eval(channels[i])-energies[i])/energies[i];
    d_e[i] = 100*(dSlope*ch_e[i])/energies[i];
  }
  linplot = new TGraphErrors(line, energies, diff, 0, d_e);
  linplot->SetTitle("Linearity");
  TCanvas *ccc = new TCanvas("lin","Linearity");
  linplot->GetXaxis()->SetTitle("Energy [keV]");
  linplot->GetYaxis()->SetTitle("Deviation from linear calibration [%]");
  linplot->Draw("a*");

  effplot = new TGraphErrors(line, energies, eff, 0, 0);
  effplot->SetTitle("Efficiency");
  TCanvas *cc4 = new TCanvas("eff","Efficiency");
  effplot->GetXaxis()->SetTitle("Energy [keV]");
  effplot->GetYaxis()->SetTitle("Efficiency");
  effplot->Draw("a*");
}

double GeCalibrate::FindFWHM(TH1I* data, double X, double Y){
  double fw = 0;
  for (int i = int(X); i>0; i--){
    if (data->GetBinContent(1+i)*2<Y) {
      fw = double(X-i); 
      break;
    }
  }
  for (int i = int(X); i<data->GetNbinsX(); i++){
    if (data->GetBinContent(1+i)*2<Y) {
      fw += double(i-X); 
      break;
    }
  }
  return fw;

}


void GeAnalyse::Compare(CalMCA_data* d1,CalMCA_data* d2){
  d1->Analyse();
  d2->Analyse();
  TH1D* data1 = d1->GetCalData();
  TH1D* data2 = d2->GetCalData();
  
  TH1D* diff = (TH1D*) data1->Clone("diff");
  diff->Add(data2,-1);
  TCanvas *c = new TCanvas("c","Comparison",1000,750);
  c->Divide(2,2);
  c->cd(1);
  data1->Draw("");
  c->cd(2);
  data2->Draw("");
  c->cd(3);
  diff->Draw();
  c->cd(4);
  data1->Draw("hist");
  data1->SetLineColor(kRed);
  data2->Draw("same hist");
  data1->GetXaxis()->SetTitle("Energy [keV]");
  data1->GetYaxis()->SetTitle("Counts/second");
  data1->Draw("same hist");
  TCanvas *x = new TCanvas("x","Comparison",900,720);  
  x->cd();
  data1->Draw();
  data2->Draw("same hist");
  data1->Draw("same hist"); 
  TLegend *leg = new TLegend(.2,.65,.8,.8);
  leg->AddEntry(data1,d1->GetName(),"L");
  leg->AddEntry(data2,d2->GetName(),"L");
  leg->Draw();
  gPad->SetLogy(1);     	    
  x->Print("compare.pdf");
  TCanvas *xx = new TCanvas("xx","Rate Difference",900,720);  
  xx->cd();
  diff->Draw();
  xx->Print("difference.pdf");
  // potassium 40 peak
  double K40_peak_err;
  double K40_peak = diff->IntegralAndError(diff->FindBin(1452),diff->FindBin(1459),K40_peak_err);
  cout << " K-40 " << K40_peak << "+/-" << K40_peak_err << endl;
  // Pb-212 238.6 keV
  double Pb212_peak_err;
  double Pb212_peak = diff->IntegralAndError(diff->FindBin(236),diff->FindBin(239.5),Pb212_peak_err);
  cout << " Pb-212(238.63) " << Pb212_peak << "+/-" << Pb212_peak_err << endl;

  double Pb214_peak_err;
  double Pb214_peak = diff->IntegralAndError(diff->FindBin(293),diff->FindBin(296),Pb214_peak_err);
  cout << " Pb-214(295.21) " << Pb214_peak << "+/-" << Pb214_peak_err << endl;

  double Pb214_peak2_err;
  double Pb214_peak2 = diff->IntegralAndError(diff->FindBin(349),diff->FindBin(353),Pb214_peak2_err);
  cout << " Pb-214(351.92) " << Pb214_peak2 << "+/-" << Pb214_peak2_err << endl;

  double Bi214_peak_err;
  double Bi214_peak = diff->IntegralAndError(diff->FindBin(605),diff->FindBin(609.5),Bi214_peak_err);
  cout << " Bi-214(609.31) " << Bi214_peak << "+/-" << Bi214_peak_err << endl;

  double Bi214_peak2_err;
  double Bi214_peak2 = diff->IntegralAndError(diff->FindBin(1114.5),diff->FindBin(1118.5),Bi214_peak2_err);
  cout << " Bi-214(1120.29) " << Bi214_peak2 << "+/-" << Bi214_peak2_err << endl;
}

void GeAnalyse::ActivityDetermination(CalMCA_data* d1,CalMCA_data* d2, IsotopeDB &db){
  d1->Analyse();
  d2->Analyse();
  TH1D* data1 = d1->GetCalData();
  TH1D* data2 = d2->GetCalData();
  
  TH1D* diff = (TH1D*) data1->Clone("diff");
  diff->Add(data2,-1);
  TCanvas *c = new TCanvas("c","Activity",1000,750);
  diff->Draw();
  TSpectrum  *ts = new TSpectrum();
  
 
}


void GeAnalyse::TimeDependence(int ld_count, CalMCA_data** ld){
  // takes list of data files and makes time dependent rate plot.

  double time[ld_count];
  double total_rate[ld_count];
  double t_rate_err[ld_count];
  double rate_208[ld_count];
  double rate_208_e[ld_count];
  double peak208_pos[ld_count];
  double peak208_height[ld_count];
  double peak208_pos_e[ld_count];
  double peak208_height_e[ld_count];
  double peak1460_pos[ld_count];
  double peak1460_height[ld_count];
  double peak1460_pos_e[ld_count];
  double peak1460_height_e[ld_count];
  double peak186_pos[ld_count];
  double peak186_height[ld_count];
  double peak186_pos_e[ld_count];
  double peak186_height_e[ld_count];

  for (int i = 0 ; i < ld_count ; i++){ //loop over all files
    time[i] = double(ld[i]->GetStartTime())+ld[i]->GetRealTime();
    if (i==0){
      TH1D* d = ld[i]->GetCalData();
      cout << "events " << d->Integral()<< " ev/s"<< endl;;
      total_rate[i] = d->Integral();
      t_rate_err[i] = sqrt(d->Integral()*ld[1]->GetLiveTime())/ld[i]->GetLiveTime();
      rate_208[i] = d->Integral(d->FindBin(205.),d->FindBin(211.));
      rate_208_e[i] = sqrt(ld[1]->GetLiveTime()*d->Integral(d->FindBin(205.),d->FindBin(211.)))/ld[i]->GetLiveTime();
      cout << " rate208: " << rate_208[i] <<" " << d->Integral(d->FindBin(205.),d->FindBin(211.))<<"delta t: " <<ld[i]->GetLiveTime() << endl;
    }
    else{
      TH1D* d = ld[i]->GetCalData();
      TH1D* dold = ld[i-1]->GetCalData();
      cout << "events" << d->Integral()<< endl;;
      total_rate[i] = (d->Integral()*ld[i]->GetLiveTime()-dold->Integral()*ld[i-1]->GetLiveTime())/(ld[i]->GetLiveTime()-ld[i-1]->GetLiveTime());
      t_rate_err[i] = sqrt(d->Integral()*ld[i]->GetLiveTime()-dold->Integral()*ld[i-1]->GetLiveTime())/(ld[i]->GetLiveTime()-ld[i-1]->GetLiveTime());
      rate_208[i] = ((ld[i]->GetLiveTime()*d->Integral(d->FindBin(205.),d->FindBin(211.)))-(ld[i-1]->GetLiveTime()*dold->Integral(dold->FindBin(205.),dold->FindBin(211.))))/(ld[i]->GetLiveTime()-ld[i-1]->GetLiveTime());
      rate_208_e[i] = sqrt(ld[i]->GetLiveTime()*d->Integral(d->FindBin(205.),d->FindBin(211.))-ld[i-1]->GetLiveTime()*dold->Integral(dold->FindBin(205.),dold->FindBin(211.)))/(ld[i]->GetLiveTime()-ld[i-1]->GetLiveTime());
      cout << " rate208: " << rate_208[i] <<" " << ld[i]->GetLiveTime()*d->Integral(d->FindBin(205.),d->FindBin(211.))<<" " <<ld[i-1]->GetLiveTime()*dold->Integral(d->FindBin(205.),d->FindBin(211.))<<" delta t: " << (ld[i]->GetLiveTime()-ld[i-1]->GetLiveTime()) << endl;
    }
  }
  double time_0 = ld[0]->GetStartTime();
  for (int i = 0 ; i < ld_count ; i++){ //remove time offset in plot.
    time[i] = (time[i] - time_0)/(3600.*24.);
  }
  TCanvas *c = new TCanvas("c","Time Dependence",1000,750);
  c->Divide(1,3);
  c->cd(1);
  TGraphErrors *g = new TGraphErrors(ld_count,time,total_rate,0,t_rate_err);
  g->Draw("a*");
  g->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g->GetYaxis()->SetTitle("Total rate [1/s]");
  c->cd(2);
  TGraphErrors *g208 = new TGraphErrors(ld_count,time,rate_208,0,rate_208_e);
  // g208->Draw("a*");
  g208->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g208->GetYaxis()->SetTitle("rate  under 208keV peak [1/s]");
  //  c->Update();
  TH1D* diff[ld_count];
  TCanvas *ct = new TCanvas("ct","spectra",1000,750);
  ct->Divide(4,(ld_count+1)/2);
  for (int i = 0 ; i < ld_count ; i++){
    // plot 208 keV peak from all files (difference only).
    char name[200];
    sprintf(name,"diff_%d",i);
    double lt = 0;
    if (i ==0){
      diff[i] = (TH1D*)ld[i]->GetCalData()->Clone(name);
      diff[i]->Scale(ld[i]->GetLiveTime());
      lt = ld[i]->GetLiveTime();
    }
    else{
      diff[i] = (TH1D*)ld[i]->GetCalData()->Clone(name);
      diff[i]->Scale(ld[i]->GetLiveTime());
      diff[i]->Add(ld[i-1]->GetCalData(),-1*ld[i-1]->GetLiveTime());
      lt = ld[i]->GetLiveTime()-ld[i-1]->GetLiveTime();
    }
    ct->cd(1+2*i);
    diff[i]->Draw();
    diff[i]->GetXaxis()->SetRangeUser(205.,211.);
    diff[i]->Fit("gaus","","",205,211);
    TF1 *f = (TF1*) diff[i]->GetFunction("gaus");
    peak208_pos[i] = f->GetParameter(1);
    peak208_pos_e[i] = f->GetParError(1);
    peak208_height[i] = sqrt(2*TMath::Pi())*f->GetParameter(2)*f->GetParameter(0)/(lt*diff[i]->GetBinWidth(10));
    peak208_height_e[i] = sqrt(pow(f->GetParError(2)/f->GetParameter(2)*peak208_height[i],2)+pow(f->GetParError(0)/f->GetParameter(0)*peak208_height[i],2));

    diff[i]->GetXaxis()->SetRangeUser(1455.,1464.);
    diff[i]->Fit("gaus","","",1455,1464);
    TF1 *f2 = (TF1*) diff[i]->GetFunction("gaus");
    peak1460_pos[i] = f2->GetParameter(1);
    peak1460_pos_e[i] = f2->GetParError(1);
    peak1460_height[i] = sqrt(2*TMath::Pi())*f2->GetParameter(2)*f2->GetParameter(0)/(lt*diff[i]->GetBinWidth(10));
    peak1460_height_e[i] = sqrt(pow(f2->GetParError(2)/f2->GetParameter(2)*peak1460_height[i],2)+pow(f2->GetParError(0)/f2->GetParameter(0)*peak1460_height[i],2));

    diff[i]->GetXaxis()->SetRangeUser(182.5,185.);
    diff[i]->Fit("gaus","","",182.5,185);
    TF1 *f3 = (TF1*) diff[i]->GetFunction("gaus");
    peak186_pos[i] = f3->GetParameter(1);
    peak186_pos_e[i] = f3->GetParError(1);
    peak186_height[i] = sqrt(2*TMath::Pi())*f3->GetParameter(2)*f3->GetParameter(0)/(lt*diff[i]->GetBinWidth(10));
    peak186_height_e[i] = sqrt(pow(f3->GetParError(2)/f3->GetParameter(2)*peak186_height[i],2)+pow(f3->GetParError(0)/f3->GetParameter(0)*peak186_height[i],2));

    cout << "DEBUG: " << i << " " << peak208_height[i] << " " << sqrt(2*TMath::Pi())<< endl;
    ct->cd(2+2*i);
    ld[i]->GetCalData()->Draw();
    ld[i]->GetCalData()->GetXaxis()->SetRangeUser(205.,211.);
    //    ld[i]->GetCalData()->Fit("gaus","","",205,211);

  }
  ct->Update();
  c->cd(2);
  TGraphErrors *g208_rate = new TGraphErrors(ld_count,time,peak208_height,0,peak208_height_e);
  g208->Draw("a*");
  g208_rate->SetMarkerStyle(24);
  g208_rate->SetMarkerColor(kRed);
  g208_rate->Draw("P");
  g208_rate->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g208_rate->GetYaxis()->SetTitle("gaussian height under 208keV peak [1/s]");
  c->cd(3);
  TGraphErrors *g208_pos = new TGraphErrors(ld_count,time,peak208_pos,0,peak208_pos_e);
  g208_pos->Draw("a*");
  g208_pos->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g208_pos->GetYaxis()->SetTitle("Peak position of 208 keV peak [keV]");
  c->Update();
  TCanvas *c3 = new TCanvas("c3","K40 peak plots", 1000,750);
  c3->Divide(1,2);
  c3->cd(1);
  TGraphErrors *g1460_rate = new TGraphErrors(ld_count,time,peak1460_height,0,peak1460_height_e);
  g1460_rate->SetMarkerStyle(24);
  g1460_rate->SetMarkerColor(kRed);
  g1460_rate->Draw("AP");
  g1460_rate->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g1460_rate->GetYaxis()->SetTitle("gaussian height under 208keV peak [1/s]");
  c3->cd(2);
  TGraphErrors *g1460_pos = new TGraphErrors(ld_count,time,peak1460_pos,0,peak1460_pos_e);
  g1460_pos->Draw("a*");
  g1460_pos->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g1460_pos->GetYaxis()->SetTitle("Peak position of 1460 keV peak [keV]");
  c3->Update();

  TCanvas *c4 = new TCanvas("c4","Ra226 peak plots", 1000,750);
  c4->Divide(1,2);
  c4->cd(1);
  TGraphErrors *g186_rate = new TGraphErrors(ld_count,time,peak186_height,0,peak186_height_e);
  g186_rate->SetMarkerStyle(24);
  g186_rate->SetMarkerColor(kRed);
  g186_rate->Draw("AP");
  g186_rate->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g186_rate->GetYaxis()->SetTitle("gaussian height under 208keV peak [1/s]");
  c4->cd(2);
  TGraphErrors *g186_pos = new TGraphErrors(ld_count,time,peak186_pos,0,peak186_pos_e);
  g186_pos->Draw("a*");
  g186_pos->GetXaxis()->SetTitle("Day since Jan 26 2012");
  g186_pos->GetYaxis()->SetTitle("Peak position of 186 226-Ra keV peak [keV]");
  c4->Update();

  TCanvas *csp = new TCanvas("csp","Final Spectrum",1000,750);
  
  ld[ld_count-1]->GetCalData()->Draw();
  csp->Update();
  
  return;
}

void GeAnalyse::BackgroundRemove(CalMCA_data* d1,CalMCA_data* d2){



}

IsotopeDB::IsotopeDB(){
}

void IsotopeDB::Init(){
}

void IsotopeDB::InitFromFile(TString fn){
  ifstream indat(fn.Data());
  char buffer[500];
  if (!indat.is_open()) {
    cerr<< "Isotope database file does not open!" << endl;
    exit (0);
  }
  while (!indat.eof()){
    indat.getline(buffer,499);
    TString line = TString(buffer);
    if (line.BeginsWith("#")) continue;
    if (line.BeginsWith("[")) {// first line of isotope entry
	Isotope i;
	
	i.SetName(line(1,line.Length()-2));
	double dat;
	indat.getline(buffer,499);
	stringstream ss(buffer);
	ss >> dat;
	//	cout << dat << endl;
	i.SetHalfLife(dat);
	int num;
	indat.getline(buffer,499);
	ss.str(buffer);
	ss  >> num;
	cout << i.GetName()<< " " << num << " lines."<< endl;
	double en,inten;
	for (int line = 0 ; line < num ; line++){
	  indat >> en >> inten ;
	  //	  cout << line << " " << en << " " << inten << endl;
	  i.SetEnergy(en);
	  i.SetIntensity(inten);
	}
	indat.getline(buffer,499);
	//	cout << buffer << endl;
	vIdb.push_back(i);
	
    }
  }
}

Isotope &IsotopeDB::GetIsotopeByName(TString name){

}

Isotope &IsotopeDB::GetIsotopeByNumber(int n){

}

Isotope &IsotopeDB::GetIsotopeByEnergy(double e){

}
