#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TChain.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>

#include "interface/untuplizer.h"
#include "interface/isPassZee.h"
#include "interface/readSample.h"
//#include "interface/ele_PU_71300/standalone_LumiReWeighting.cc"
#include "interface/ele_PU_69200/standalone_LumiReWeighting.cc"
#include "interface/dataFilter.h"


void xAna_ele(std::string inputFile, std::string outputFolder, std::string outputFile){

  cout<<"finish compiling and start to run macro"<< endl;
  // --------  open  ntuples ------------ 

  // get paths of every root files
  std::vector<string> infiles; 
  readSample(inputFile, infiles);

  cout<<"infiles.size():"<< infiles.size() << endl;

  // combine TTree
  TChain* tree = new TChain("tree/treeMaker");

  for(unsigned int i=0 ; i<infiles.size() ; i++ ){
    std::string open_root="root://" + infiles[i] ;
    tree->Add( open_root.c_str()  );
  }

  cout<<"tree->GetEntries(): "<< tree->GetEntries() << endl;

  // read the TTree

  TreeReader data( tree );
  cout<<"finish TreeReader data( tree );"<< endl; 

  // ------- Declare the histogram ----------------

  TH1D* h_totalEvents      = new TH1D("h_totalEv",          "totalEvents",       10,     0,   10);

  TH1D* h_pu_nTrueInt      = new TH1D("h_pu_nTrueInt",      "h_pu_nTrueInt",     60,     0,   60);
  TH1D* h_nVtx             = new TH1D("h_nVtx",             "nVtx",              30,  -0.5, 29.5); 

  TH1D* h_leadElePt        = new TH1D("h_leadElePt",        "leadElePt",         16,     0,  800);
  TH1D* h_leadEleEta       = new TH1D("h_leadEleEta",       "leadEleEta",        40,    -4,    4);
  TH1D* h_subleadElePt     = new TH1D("h_subleadElePt",     "subleadElePt",      25,     0,  500);
  TH1D* h_subleadEleEta    = new TH1D("h_subleadEleEta",    "subleadEleEta",     40,    -4,    4);

  TH1D* h_Zmass            = new TH1D("h_Zmass",            "Zmass",             30,    60,  120);
  TH1D* h_Zpt              = new TH1D("h_Zpt",              "Zpt",               50,     0, 1000);
  TH1D* h_Zeta             = new TH1D("h_Zeta",             "Zeta",              40,    -4,    4);
  TH1D* h_ZRapidity        = new TH1D("h_ZRapidity",        "ZRapidity",         40,    -4,    4);

  TH1D* h_FATjetPt         = new TH1D("h_FATjetPt",         "FATjetPt",           20,  100, 1000);
  TH1D* h_FATjetEta        = new TH1D("h_FATjetEta",        "FATjetEta",          20,   -4,    4);
  TH1D* h_FATjetCISVV2     = new TH1D("h_FATjetCISVV2",     "FATjetCISVV2",       20,    0,  1.2);
  TH1D* h_FATjetSDmass     = new TH1D("h_FATjetSDmass",     "FATjetSDmass",       20,    0,  200);
  TH1D* h_FATjetPRmass     = new TH1D("h_FATjetPRmass",     "FATjetPRmass",       20,    0,  200);
  TH1D* h_FATjetPRmassCorr = new TH1D("h_FATjetPRmassCorr", "FATjetPRmassCorr",   20,    0,  200);
  TH1D* h_FATjetTau1       = new TH1D("h_FATjetTau1",       "FATjetTau1",         20,    0,    1);
  TH1D* h_FATjetTau2       = new TH1D("h_FATjetTau2",       "FATjetTau2",         20,    0,    1);
  TH1D* h_FATjetTau2dvTau1 = new TH1D("h_FATjetTau2dvTau1", "FATjetTau2dvTau1",   20,    0,    1);

  TH1D* h_FATsubjetPt      = new TH1D("h_FATsubjetPt",      "FATsubjetPt",        20,    0,  800);
  TH1D* h_FATsubjetEta     = new TH1D("h_FATsubjetEta",     "FATsubjetEta",       20,   -4,    4);  
  TH1D* h_FATsubjetSDCSV   = new TH1D("h_FATsubjetSDCSV",   "FATsubjetSDCSV",     20,    0,  1.2);

  TH1D* h_ZHmass           = new TH1D("h_ZHmass",           "ZHmass",             50,    0, 5000);

  h_pu_nTrueInt     ->Sumw2();
  h_nVtx            ->Sumw2();

  h_leadElePt       ->Sumw2();
  h_leadEleEta      ->Sumw2();
  h_subleadElePt    ->Sumw2();
  h_subleadEleEta   ->Sumw2();

  h_Zmass           ->Sumw2();
  h_Zpt             ->Sumw2();
  h_Zeta            ->Sumw2();
  h_ZRapidity       ->Sumw2();

  h_FATjetPt        ->Sumw2();   
  h_FATjetEta       ->Sumw2();
  h_FATjetCISVV2    ->Sumw2();
  h_FATjetSDmass    ->Sumw2();
  h_FATjetPRmass    ->Sumw2();
  h_FATjetPRmassCorr->Sumw2();
  h_FATjetTau1      ->Sumw2();
  h_FATjetTau2      ->Sumw2();
  h_FATjetTau2dvTau1->Sumw2();
  h_FATsubjetPt     ->Sumw2();
  h_FATsubjetEta    ->Sumw2();
  h_FATsubjetSDCSV  ->Sumw2();

  h_ZHmass          ->Sumw2();

  int counter_0=0;
  int counter_1=0;
  int counter_2=0;    
  int counter_3=0;

  cout<<"finish define Histogram"<< endl;

  // ------ begin of event loop ---------------

//  int Number_to_print    =  100000;
  int Number_to_print    =  5000;

//  int Max_number_to_read = 1000000;
  int Max_number_to_read = 30000;
//  int Max_number_to_read =  10;

  bool PU_weight_flag = true;
//  PU_weight_flag = false; // turn off the PU weight

  if( PU_weight_flag ) cout<<"Will use Pile Up weight"<<endl;
  else cout<<"Will not use Pile Up weight"<<endl;

  cout<<"staring loop events"<< endl;

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev %   Number_to_print == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    if( ev > Max_number_to_read) break;

    data.GetEntry(ev);

    h_totalEvents->Fill(1,1);

    // get Branch

    Bool_t   isData                  = data.GetBool("isData");
    Float_t  ntrue                   = data.GetFloat("pu_nTrueInt");
    Int_t    nVtx                    = data.GetInt("nVtx");

    const TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");

    Int_t          FATnJet           = data.GetInt("FATnJet");    
    Int_t*         FATnSubSDJet      = data.GetPtrInt("FATnSubSDJet");
    Float_t*       FATjetCISVV2      = data.GetPtrFloat("FATjetCISVV2");
    Float_t*       FATjetSDmass      = data.GetPtrFloat("FATjetSDmass");
    Float_t*       FATjetPRmass      = data.GetPtrFloat("FATjetPRmass");
    Float_t*       FATjetPRmassCorr  = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Float_t*       FATjetTau1        = data.GetPtrFloat("FATjetTau1");
    Float_t*       FATjetTau2        = data.GetPtrFloat("FATjetTau2");
    TClonesArray*  FATjetP4          = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    vector<bool>&  FATjetPassIDLoose = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    vector<float>* FATsubjetSDPx     = data.GetPtrVectorFloat("FATsubjetSDPx", FATnJet);
    vector<float>* FATsubjetSDPy     = data.GetPtrVectorFloat("FATsubjetSDPy", FATnJet);
    vector<float>* FATsubjetSDPz     = data.GetPtrVectorFloat("FATsubjetSDPz", FATnJet);
    vector<float>* FATsubjetSDE      = data.GetPtrVectorFloat("FATsubjetSDE", FATnJet);
    vector<float>* FATsubjetSDCSV    = data.GetPtrVectorFloat("FATsubjetSDCSV", FATnJet);

    // get Pile-Up weight 


    double PU_weight_central =1;// weight=1 for data

    if(!isData && PU_weight_flag){// for MC

      standalone_LumiReWeighting LumiWeights_central(0);
      PU_weight_central = LumiWeights_central.weight(ntrue);

    }
   
//    cout<<"PU_weight_central: "<< PU_weight_central << endl;
    double eventWeight = PU_weight_central;   

    h_pu_nTrueInt   ->Fill(ntrue , eventWeight);
 
    // data filter 

    bool CSCT       = FilterStatus(data, "Flag_CSCTightHaloFilter");
    bool eeBadSc    = FilterStatus(data, "Flag_eeBadScFilter");
    bool Noise      = FilterStatus(data, "Flag_HBHENoiseFilter");
    bool NoiseIso   = FilterStatus(data, "Flag_HBHENoiseIsoFilter");

    bool filter_pass = false; 
    if( isData && CSCT && eeBadSc && Noise && NoiseIso ) filter_pass = true;// Data 
    else if (!isData) filter_pass = true; // MC, doesn't apply data filter

    if( !filter_pass ) continue;



    //  apply HLT 

    bool eleTrigger1 = TriggerStatus(data, "HLT_Ele105_CaloIdVT_GsfTrkIdT_v");

    bool HLT_pass = false;
    if( isData && eleTrigger1 ) HLT_pass= true;// Data
    else if (!isData) HLT_pass = true; // MC, doesn't apply HLT

    if( !HLT_pass ) continue;

    counter_0++;



    //  nVtx>=1 

    if( nVtx < 1 ) continue;

    counter_1++;



    // select good electrons and Z boson 


    vector<Int_t> goodEleID;
    if( !isPassZee(data, goodEleID) ) continue;

    counter_2++;

    TLorentzVector* thisEle = (TLorentzVector*)eleP4->At(goodEleID[0]);
    TLorentzVector* thatEle = (TLorentzVector*)eleP4->At(goodEleID[1]);


 
    // select Fat jet for Higgs 

    Int_t goodFATJetID = -1;
    TLorentzVector* thisJet = NULL;

    for(Int_t ij = 0; ij < FATnJet; ++ij){

      thisJet = (TLorentzVector*)FATjetP4->At(ij);

      if( thisJet->Pt() < 200 ) continue;
      if( fabs(thisJet->Eta()) > 2.4 ) continue;
      if( !FATjetPassIDLoose[ij] ) continue;
      if( thisJet->DeltaR(*thisEle) < 0.8 || thisJet->DeltaR(*thatEle) < 0.8 ) continue;
      if( !( FATjetPRmassCorr[ij] < 65 || FATjetPRmassCorr[ij] > 135) ) continue;// side band region

      goodFATJetID = ij;
      break;

    }

    if( goodFATJetID < 0 ) continue;

    counter_3++;



    // fill histogram in physical object level 

    // nVtx

    h_nVtx            ->Fill(nVtx,eventWeight);

    // Leading and SubLeading Electron 

    if( thisEle->Pt() > thatEle->Pt() ){

      h_leadElePt    ->Fill(thisEle->Pt(),eventWeight);
      h_leadEleEta   ->Fill(thisEle->Eta(),eventWeight);
      h_subleadElePt ->Fill(thatEle->Pt(),eventWeight);
      h_subleadEleEta->Fill(thatEle->Eta(),eventWeight);

    }else{

      h_leadElePt    ->Fill(thatEle->Pt(),eventWeight);
      h_leadEleEta   ->Fill(thatEle->Eta(),eventWeight);
      h_subleadElePt ->Fill(thisEle->Pt(),eventWeight);
      h_subleadEleEta->Fill(thisEle->Eta(),eventWeight);

    }

    // Z boson

    TLorentzVector l4_Z = (*thisEle+*thatEle);

    h_Zmass    ->Fill( l4_Z.M()       ,eventWeight);
    h_Zpt      ->Fill( l4_Z.Pt()      ,eventWeight);
    h_Zeta     ->Fill( l4_Z.Eta()     ,eventWeight);
    h_ZRapidity->Fill( l4_Z.Rapidity(),eventWeight);

    // Fat Jet for Higgs

    h_FATjetPt        ->Fill(thisJet->Pt(),eventWeight);
    h_FATjetEta       ->Fill(thisJet->Eta(),eventWeight);
    h_FATjetCISVV2    ->Fill(FATjetCISVV2[goodFATJetID],eventWeight);
    h_FATjetSDmass    ->Fill(FATjetSDmass[goodFATJetID],eventWeight);
    h_FATjetPRmass    ->Fill(FATjetPRmass[goodFATJetID],eventWeight);
    h_FATjetPRmassCorr->Fill(FATjetPRmassCorr[goodFATJetID],eventWeight);
    h_FATjetTau1      ->Fill(FATjetTau1[goodFATJetID],eventWeight);
    h_FATjetTau2      ->Fill(FATjetTau2[goodFATJetID],eventWeight);
    h_FATjetTau2dvTau1->Fill(FATjetTau2[goodFATJetID]/FATjetTau1[goodFATJetID],eventWeight);

    // SubJet 

    for(Int_t is = 0; is < FATnSubSDJet[goodFATJetID]; ++is){

      TLorentzVector l4_subjet(0,0,0,0);

      l4_subjet.SetPxPyPzE(FATsubjetSDPx[goodFATJetID][is],
			   FATsubjetSDPy[goodFATJetID][is],
			   FATsubjetSDPz[goodFATJetID][is],
			   FATsubjetSDE [goodFATJetID][is]);

      h_FATsubjetPt   ->Fill(l4_subjet.Pt(),eventWeight);
      h_FATsubjetEta  ->Fill(l4_subjet.Eta(),eventWeight);
      h_FATsubjetSDCSV->Fill(FATsubjetSDCSV[goodFATJetID][is],eventWeight);

    }

    // ZH invariant mass
    
    TLorentzVector l4_ZH = l4_Z + ( *thisJet );     

    h_ZHmass ->Fill( l4_ZH.M()   ,eventWeight);

  } // ------------------ end of event loop ------------------

  cout<<"counter_0: "<< counter_0 << endl;
  cout<<"counter_1: "<< counter_1 << endl;
  cout<<"counter_2: "<< counter_2 << endl;
  cout<<"counter_3: "<< counter_3 << endl;


  fprintf(stderr, "Processed all events\n");

  std::string root_name = Form("%s.root",outputFile.c_str() ) ;
  TString save_path = outputFolder +"/"+  root_name;

  TFile* outFile = new TFile( save_path , "recreate");

  h_pu_nTrueInt     ->Write("pu_nTrueInt");
  h_nVtx            ->Write("nVtx");

  h_leadElePt       ->Write("leadElePt");
  h_leadEleEta      ->Write("leadEleEta");
  h_subleadElePt    ->Write("subleadElePt");
  h_subleadEleEta   ->Write("subleadEleEta");

  h_Zmass           ->Write("Zmass");
  h_Zpt             ->Write("Zpt");
  h_Zeta            ->Write("Zeta");
  h_ZRapidity       ->Write("ZRapidity");
 
  h_FATjetPt        ->Write("FATjetPt");
  h_FATjetEta       ->Write("FATjetEta");
  h_FATjetCISVV2    ->Write("FATjetCISVV2");
  h_FATjetSDmass    ->Write("FATjetSDmass");
  h_FATjetPRmass    ->Write("FATjetPRmass");
  h_FATjetPRmassCorr->Write("FATjetPRmassCorr");
  h_FATjetTau1      ->Write("FATjetTau1");
  h_FATjetTau2      ->Write("FATjetTau2");
  h_FATjetTau2dvTau1->Write("FATjetTau2dvTau1");

  h_FATsubjetPt     ->Write("FATsubjetPt");
  h_FATsubjetEta    ->Write("FATsubjetEta");
  h_FATsubjetSDCSV  ->Write("FATsubjetSDCSV");

  h_ZHmass          ->Write("ZHmass");

  h_totalEvents     ->Write("totalEvents");

  outFile->Write();  
  delete outFile;


}




