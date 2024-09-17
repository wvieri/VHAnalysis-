#include <string>
#include <iostream>
#include <fstream>
#include <TH1.h>
#include <TFile.h>

class readHist{
    
 public:
  
  readHist(string);
  TH1D* getHist(string);

 private:

  TFile* thisFile;
  string thisFileName;
  float crossSection(string);
   
};

readHist::readHist(string rootFile){

  string tmpStr = rootFile.erase(rootFile.find_last_not_of("/")+1);
  thisFileName = tmpStr.substr(tmpStr.find_last_of("/")+1);
  thisFile = TFile::Open(rootFile.data());

}

float readHist::crossSection(string thisPath){

  ifstream textFile("/afs/cern.ch/work/y/yuchang/ZH_analysis_with_2016_data/CMSSW_8_0_8/src/DataMC_comparison/interface/xSec.txt");
//  ifstream textFile("/afs/cern.ch/work/h/htong/ZpZHllbb_13TeV/xSec.txt");
  string token;
  float crosssection = 0., thisNum = 0.;

  while( textFile >> token >> thisNum ){

    if( thisPath.find(token) != string::npos )
      crosssection = thisNum;

  }
  
  return crosssection;

}

TH1D* readHist::getHist(string hname){

  TH1D* thisHist = (TH1D*)(thisFile->Get(Form("%s", hname.c_str())));  

//  cout<<"here is readHist::getHist"<< endl;
//  cout<<"thisFileName: "<< thisFileName << endl;

//  cout<<"histo open name: "<<Form( "%s", hname.c_str()  ) ;
  

//  cout<<"crossSection(thisFileName.data()): "<< crossSection(thisFileName.data()) << endl;
//  cout<<"(TH1D*)(thisFile->Get(totalEvents)))->Integral(): "<< ((TH1D*)(thisFile->Get("totalEvents")))->Integral()<< endl;

//  thisHist->Scale(thisFileName.find("Run2016") != string::npos ? 1 : 2512.*crossSection(thisFileName.data())/((TH1D*)(thisFile->Get("totalEvents")))->Integral());

//double lumi = 4327.13;
//double lumi = 1/178*4327.13;// approx for 1M events
double lumi = 432.713;

  thisHist->Scale(thisFileName.find("Run2016") != string::npos ? 1 : lumi*crossSection(thisFileName.data())/((TH1D*)(thisFile->Get("totalEvents")))->Integral());

  return thisHist;

}
