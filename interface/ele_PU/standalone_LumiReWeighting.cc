/* information for how I get pileup weight
 
  take MC's pu_nTrueInt distribution from 
    https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_X/SimGeneral/MixingModule/python/mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi.py 

  data's pu_nTrueInt distribution are computed with
    My Json file: /afs/cern.ch/work/y/yuchang/produce_samples/CMSSW_8_0_8/src/DelPanj/CrabUtilities/crab_20160419/crab_SingleElectron-Run2016B-PromptReco-v2/results/processedLumis.json

    INPUTLUMIJSON file : /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt
    minBiasXsec : 71300

    see more in https://twiki.cern.ch/twiki/bin/viewauth/CMS/ExoDiBosonResonancesRun2#PU_weights
   
  finaly get PU weight in standalone_LumiReWeighting.cc   
    see more in https://github.com/syuvivida/xtozh_common/tree/76X_analysis/macro_examples/pileup_reweight

*/

/**
   \class    standalone_LumiReWeighting standalone_LumiReWeighting.h "PhysicsTools/Utilities/interface/standalone_LumiReWeighting.h"
   \brief    Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data

   This class will trivially take two histograms:
   1. The generated "flat-to-N" distributions from a given processing (or any other generated input)
   2. A histogram generated from the "estimatePileup" macro here:

   https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc#How_to_use_script_estimatePileup

   and produce weights to convert the input distribution (1) to the latter (2).

   \author Shin-Shan Eiko Yu, Salvatore Rappoccio, modified by Mike Hildreth
  
*/
#ifndef standalone_LumiReWeighting_cxx
#define standalone_LumiReWeighting_cxx
#include "TH1.h"
#include "TFile.h"
#include <string>
#include "standalone_LumiReWeighting.h"

//=======================================================
// For 2012 Data and MC
//=======================================================


Double_t Summer2012_S10[60] = {
		0.000829312873542,
 		0.00124276120498,
 		0.00339329181587,
 		0.00408224735376,
 		0.00383036590008,
		0.00659159288946,
 		0.00816022734493,
 		0.00943640833116,
 		0.0137777376066,
 		0.017059392038,
 		0.0213193035468,
 		0.0247343174676,
 		0.0280848773878,
 		0.0323308476564,
 		0.0370394341409,
 		0.0456917721191,
 		0.0558762890594,
 		0.0576956187107,
 		0.0625325287017,
 		0.0591603758776,
 		0.0656650815128,
 		0.0678329011676,
 		0.0625142146389,
 		0.0548068448797,
 		0.0503893295063,
 		0.040209818868,
 		0.0374446988111,
 		0.0299661572042,
 		0.0272024759921,
 		0.0219328403791,
 		0.0179586571619,
 		0.0142926728247,
 		0.00839941654725,
 		0.00522366397213,
 		0.00224457976761,
 		0.000779274977993,
 		0.000197066585944,
 		7.16031761328e-05,
     };

double Data2012[60]={
1363.19,
24112.4,
160204,
381537,
714502,
956099,
1.18106e+06,
3.536e+06,
1.48677e+07,
2.83833e+07,
4.19268e+07,
6.41916e+07,
9.44924e+07,
1.36196e+08,
1.94289e+08,
2.62877e+08,
3.2434e+08,
3.64754e+08,
3.84137e+08,
3.86116e+08,
3.72314e+08,
3.44915e+08,
3.07339e+08,
2.63746e+08,
2.18194e+08,
1.73797e+08,
1.32715e+08,
9.65661e+07,
6.65305e+07,
4.31757e+07,
2.6302e+07,
1.50233e+07,
8.05577e+06,
4.06978e+06,
1.9485e+06,
891383,
394079,
171222,
75165.2,
34932.4,
18357.1,
11523.6,
8626.26,
7313.52,
6655.46,
6289.87,
6073.54,
5944.48,
5869.7,
5826.94,
5798.89,
5771.83,

};


double Data2012Up[60]={
820.261,
20621.1,
135995,
297808,
645034,
834386,
1.04507e+06,
1.90678e+06,
8.64011e+06,
2.17854e+07,
3.26788e+07,
4.84338e+07,
7.18593e+07,
1.02528e+08,
1.45086e+08,
2.01703e+08,
2.64158e+08,
3.16932e+08,
3.50622e+08,
3.66567e+08,
3.67644e+08,
3.5511e+08,
3.30759e+08,
2.9737e+08,
2.58341e+08,
2.1709e+08,
1.76352e+08,
1.38051e+08,
1.03621e+08,
7.41648e+07,
5.03637e+07,
3.23302e+07,
1.958e+07,
1.11873e+07,
6.04232e+06,
3.09687e+06,
1.51486e+06,
712727,
325989,
147260,
67417,
32628,
17684.7,
11261,
8432.19,
7113.22,
6440.98,
6063.63,
5836.9,
5697.82,
5614.14,
5565.17,

};

double Data2012Down[60]={
2144.18,
28780.9,
187588,
487719,
794113,
1.08132e+06,
1.53711e+06,
7.18794e+06,
2.25803e+07,
3.61842e+07,
5.59784e+07,
8.58628e+07,
1.26612e+08,
1.85178e+08,
2.59488e+08,
3.30831e+08,
3.79617e+08,
4.0336e+08,
4.06527e+08,
3.91256e+08,
3.60234e+08,
3.17715e+08,
2.68831e+08,
2.18411e+08,
1.70002e+08,
1.2605e+08,
8.83804e+07,
5.81902e+07,
3.57878e+07,
2.05019e+07,
1.09415e+07,
5.45622e+06,
2.55726e+06,
1.13627e+06,
484466,
201789,
84588.1,
37579.2,
19083.5,
11801.9,
8839.96,
7537.79,
6895.75,
6543.22,
6338.84,
6221.11,
6155.62,
6117.94,
6088.98,
6053.17,
5994.32,
5909.29,

};



standalone_LumiReWeighting::standalone_LumiReWeighting(int mode) {

//  std::cout << "=======================================================================" << std::endl;
  
  std::vector<double> MC_distr;
  std::vector<double> Lumi_distr;

  MC_distr.clear();
  Lumi_distr.clear();
  switch (mode)
    {
    case 0:
//      std::cout << "Using central value " << std::endl;
      break;
    case 1:
      std::cout << "Using +1 sigma 5% value " << std::endl;
      break;
    case -1:
      std::cout << "Using -1 sigma 5% value " << std::endl;
      break;
    default:
      std::cout << "Using central value " << std::endl;
      break;
    } // end of switch

  Int_t NBins = 60;
  
  for( int i=0; i< NBins; ++i) {
    switch (mode){
    case 0:
      Lumi_distr.push_back(Data2012[i]);
      break;
    case 1:
      Lumi_distr.push_back(Data2012Up[i]);
      break;
    case -1:
      Lumi_distr.push_back(Data2012Down[i]);
      break;
    default:
      Lumi_distr.push_back(Data2012[i]);
      break;
    } // end of switch

    MC_distr.push_back(Summer2012_S10[i]);
  } // end of loop over bins

  // no histograms for input: use vectors
  
  // now, make histograms out of them:

  // first, check they are the same size...

  if( MC_distr.size() != Lumi_distr.size() ){   
    std::cout << "MC_distr.size() = " << MC_distr.size() << std::endl;
    std::cout << "Lumi_distr.size() = " << Lumi_distr.size() << std::endl;
    std::cerr <<"ERROR: standalone_LumiReWeighting: input vectors have different sizes. Quitting... \n";

  }


  weights_ = new TH1D(Form("luminumer_%d",mode),
 		      Form("luminumer_%d",mode),
 		      NBins,0.0, double(NBins));

  TH1D* den = new TH1D(Form("lumidenom_%d",mode),
 		       Form("lumidenom_%d",mode),
 		       NBins,0.0, double(NBins));


  
  for(int ibin = 1; ibin<NBins+1; ++ibin ) {
    weights_->SetBinContent(ibin, Lumi_distr[ibin-1]);
    den->SetBinContent(ibin,MC_distr[ibin-1]);
  }

/*
  std::cout << "Data Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }
  std::cout << "MC Input " << std::endl;
  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << den->GetBinContent(ibin) << std::endl;
  }
*/

  // check integrals, make sure things are normalized

  double deltaH = weights_->Integral();
  if(fabs(1.0 - deltaH) > 0.02 ) { //*OOPS*...
    weights_->Scale( 1.0/ weights_->Integral() );
  }
  double deltaMC = den->Integral();
  if(fabs(1.0 - deltaMC) > 0.02 ) {
    den->Scale(1.0/ den->Integral());
  }

  weights_->Divide( den );  // so now the average weight should be 1.0    

/*
  std::cout << "Reweighting: Computed Weights per In-Time Nint " << std::endl;


  for(int ibin = 1; ibin<NBins+1; ++ibin){
    std::cout << "   " << ibin-1 << " " << weights_->GetBinContent(ibin) << std::endl;
  }

  std::cout << "=======================================================================" << std::endl;
*/

  delete den; 
}

standalone_LumiReWeighting::~standalone_LumiReWeighting()
{
}



double standalone_LumiReWeighting::weight( double npv ) {
  int bin = weights_->GetXaxis()->FindBin( npv );
  double weight_value = weights_->GetBinContent( bin );
  delete weights_;
  return weight_value;

//  return weights_->GetBinContent( bin );
}


#endif
