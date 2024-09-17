#include "XSecsLumi.h"
#include <string>
#include <iostream>
#include <TPad.h>
#include <TH1D.h>
#include <TFile.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <THStack.h>
#include <TLegend.h>
#include <TTree.h>
#include <TKey.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TObject.h"
#include "TColor.h"
#include "TFile.h"
#include "TLine.h"
#include "TMarker.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "THStack.h"
#include "TMath.h"
#include "TVectorD.h"
#include "TMatrix.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TVirtualFitter.h"
#include "TString.h"

#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <TSystemDirectory.h>


//string path = "/afs/cern.ch/work/y/yuchang/public/show_vieri_DataMC_comparison/output_with_PUweight/";

string path = "/afs/cern.ch/work/y/yuchang/ZH_analysis_with_2016_data/CMSSW_8_0_8/src/DataMC_comparison/test_ele_with_PUweight/output_ele/";

//string path = "/afs/cern.ch/user/v/vieri/work/ZH_Analysis/CMSSW_8_0_15/src/ZHAnalysis/DataMC_comparison/";

//void VHPlotter(string title="", int plot=0) {
void VHPlotter(string title="", int plot=0,   TCanvas* c1=0 ) {


string subdir="0";
string postfix="";
string dirbSel="";

double norm11 = ((Lumi2016_ele * Xsec_dy_amc1) / Ngen_dy_amc1);
double norm12 = ((Lumi2016_ele * Xsec_dy_amc2) / Ngen_dy_amc2);
double norm13 = ((Lumi2016_ele * Xsec_dy_amc3) / Ngen_dy_amc3);
double norm14 = ((Lumi2016_ele * Xsec_dy_amc4) / Ngen_dy_amc4);

double norm2 = ((Lumi2016_ele * Xsec_tt) / Ngen_tt);
double norm3 = ((Lumi2016_ele * Xsec_zz) / Ngen_zz);
double norm4 = ((Lumi2016_ele * Xsec_wz) / Ngen_wz);
double norm6 = ((Lumi2016_ele * Xsec_ww) / Ngen_ww);
double norm7 = ((Lumi2016_ele * Xsec_zhm) / Ngen_zhm);
double norm8 = ((Lumi2016_ele * Xsec_zhp) / Ngen_zhp);

//if (ilepton==1) Lumi2016 = Lumi2016_ele;

/*double enorm1 = ((Lumi2012 * eXsec_dy) / Ngen_dy);
double enorm1_amc = ((Lumi2012 * eXsec_dy_amc) / Ngen_dy_amc);
double enorm1_1 = ((Lumi2012 * eXsec_dy_1) / Ngen_dy_1);
double enorm1_2=0;
double enorm2 = ((Lumi2012 * eXsec_tt) / Ngen_tt);
double enorm4 = ((Lumi2012 * eXsec_wz) / Ngen_wz);
double enorm5 = ((Lumi2012 * eXsec_qcd) / Ngen_qcd);
double enorm6 = ((Lumi2012 * eXsec_ww) / Ngen_ww);
double enorm7 = ((Lumi2012 * eXsec_wj) / Ngen_wj);
double enorm8 = ((Lumi2012 * eXsec_tS) / Ngen_tS);
double enorm9 = ((Lumi2012 * eXsec_tT) / Ngen_tT);
double enorm10 = ((Lumi2012 * eXsec_tW) / Ngen_tW);
double enorm11 = ((Lumi2012 * eXsec_tSb) / Ngen_tSb);
double enorm12 = ((Lumi2012 * eXsec_tTb) / Ngen_tTb);
double enorm13 = ((Lumi2012 * eXsec_tWb) / Ngen_tWb);
*/

TFile* mc11 = TFile::Open((path + "DYJetsToLL_M-50_HT-100to200_13TeV.root").c_str());
TH1F* h_mc11 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc12 = TFile::Open((path + "DYJetsToLL_M-50_HT-200to400_13TeV.root").c_str());
TH1F* h_mc12 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc13 = TFile::Open((path + "DYJetsToLL_M-50_HT-400to600_13TeV.root").c_str());
TH1F* h_mc13 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc14 = TFile::Open((path + "DYJetsToLL_M-50_HT-600toInf_13TeV.root").c_str());
TH1F* h_mc14 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc2 = TFile::Open((path + "TT_TuneCUETP8M1_13TeV.root").c_str());
TH1F* h_mc2 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc3 = TFile::Open((path + "ZZ_TuneCUETP8M1_13TeV.root").c_str());
TH1F* h_mc3 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc4 = TFile::Open((path + "WZ_TuneCUETP8M1_13TeV.root").c_str());
TH1F* h_mc4 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc6 = TFile::Open((path + "WW_TuneCUETP8M1_13TeV.root").c_str());
TH1F* h_mc6 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc7 = TFile::Open((path + "ZH_HToBB_ZToLL_M125_13TeV_amcatnlo.root").c_str());
TH1F* h_mc7 = (TH1F*)gDirectory->Get(title.c_str());
TFile* mc8 = TFile::Open((path + "ZH_HToBB_ZToLL_M125_13TeV_powheg.root").c_str());
TH1F* h_mc8 = (TH1F*)gDirectory->Get(title.c_str());

TFile* data = TFile::Open((path + "SingleElectron-Run2016B-v2.root").c_str());;
TH1F* h_data = (TH1F*)gDirectory->Get(title.c_str());

h_mc2 -> SetLineColor(kBlack);
h_mc2 -> SetFillColor(kBlue);
//h_mc2 -> SetFillStyle(3004);

h_mc3 -> SetLineColor(kBlack);
h_mc3 -> SetFillColor(kGray+2);
//h_mc3 -> SetFillStyle(3004);

h_mc4 -> SetLineColor(kBlack);
h_mc4 -> SetFillColor(kGray+3);
//h_mc4 -> SetFillStyle(3004);

h_mc6 -> SetLineColor(kBlack);
h_mc6 -> SetFillColor(kRed+2);
//h_mc6 -> SetFillStyle(3004);

h_mc7 -> SetLineColor(kBlack);
h_mc7 -> SetFillColor(kGray);
//h_mc7 -> SetFillStyle(3004);

h_mc8 -> SetLineColor(kBlack);
h_mc8 -> SetFillColor(kPink);
//h_mc8 -> SetFillStyle(3004);

h_mc11->Scale(norm11);
h_mc12->Scale(norm12);
h_mc13->Scale(norm13);
h_mc14->Scale(norm14);
h_mc2->Scale(norm2);
h_mc3->Scale(norm3);
h_mc4->Scale(norm4);
h_mc6->Scale(norm6);
h_mc7->Scale(norm7);
h_mc8->Scale(norm8);

TH1F* h_mcDY = (TH1F*)h_mc11->Clone("h_mcDY");
h_mcDY->Reset();
h_mcDY->Add(h_mc14);
h_mcDY->Add(h_mc13);
h_mcDY->Add(h_mc12);
h_mcDY->Add(h_mc11);

h_mcDY -> SetLineColor(kBlack);
h_mcDY -> SetFillColor(kYellow-4);
//h_mcDY -> SetFillStyle(3004);

TH1F *ht = (TH1F*)h_mcDY->Clone("ht");
ht->Reset();
//ht->Add(h_mc8);
ht->Add(h_mc7);
ht->Add(h_mc6);
ht->Add(h_mc4);
ht->Add(h_mc3);
ht->Add(h_mc2);
ht->Add(h_mcDY);

THStack *hs = new THStack("hs","");
hs->Add(h_mc8);
hs->Add(h_mc7);
hs->Add(h_mc6);
hs->Add(h_mc4);
hs->Add(h_mc3);
hs->Add(h_mc2);
hs->Add(h_mcDY);

//TCanvas* c1 = 0;
//c1 = new TCanvas("c","c",10,10,800,600);
c1->cd();

TPad *pad1 = new TPad("pad1","pad1",0.0,0.3,1.0,1.0);
pad1->SetBottomMargin(0.001);
pad1->Draw();
pad1->cd();
//pad1->SetLogy();

hs->Draw("HIST");
hs->GetYaxis()->SetTitle("Events");
hs->GetYaxis()->SetTitleSize(0.05);
hs->GetYaxis()->SetLabelSize(0.045);
hs->GetYaxis()->SetTitleOffset(1.0);// 0.7
hs->SetMinimum(8);
hs->SetMaximum(1.2*hs->GetMaximum());
if (title=="ZHmass") { hs->GetXaxis()->SetRangeUser(0, 3000); }

h_data->Draw("EPX0SAMES");
h_data->SetMarkerColor(kBlack);
h_data->SetMarkerStyle(20);
h_data->SetMarkerSize (1.0);
h_data->SetStats(0);

TLegend *leg;

if (title=="FATjetTau2dvTau1") { leg = new TLegend(0.15, 0.547, 0.41, 0.88); }
else { leg = new TLegend(0.65, 0.547, 0.91, 0.88);}
//leg = new TLegend(0.65, 0.547, 0.91, 0.88);
leg->SetBorderSize(0);
leg->SetEntrySeparation(0.01);
leg->SetFillColor(0);
leg->SetFillStyle(0);

leg->AddEntry(h_data,"Data","p");
leg->AddEntry(h_mcDY,"DY","f");
leg->AddEntry(h_mc2,"t#bar{t}","f");
leg->AddEntry(h_mc3,"ZZ","f");
leg->AddEntry(h_mc4,"WZ","f");
leg->AddEntry(h_mc6,"WW","f");
leg->AddEntry(h_mc7,"ZH madgraph","f");
leg->AddEntry(h_mc8,"ZH powheg","f");

leg->Draw();

TLatex *lar = new TLatex();
lar->SetNDC(kTRUE);
lar->SetTextSize(0.04);
lar->SetLineWidth(5);
lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2016}}");
lar->DrawLatex(0.60, 0.94, "L = 4.327 fb^{-1} at #sqrt{s} = 13 TeV");

pad1->Update();
c1->Update();

c1->cd();

TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");

TPad *pad2 = new TPad("pad2","pad2",0,0,1,0.3);// 0.3
pad2->SetTopMargin(0);
pad2->SetBottomMargin(0.3);
pad2->Draw();
pad2->cd();
h_ratio->SetTitle("");
h_ratio->SetStats(0);

if (title=="Zpt") {
  h_ratio->GetXaxis ()->SetTitle("p_{T}(Z) [GeV/c]");
} else if (title=="Zmass") {
  h_ratio->GetXaxis ()->SetTitle("M(#ell#ell) [GeV/c^{2}]");
} else if (title=="ZRapidity") {
  h_ratio->GetXaxis ()->SetTitle("y(Z) [GeV/c^{2}]");
} else if (title=="Zeta") {
  h_ratio->GetXaxis ()->SetTitle("y(Z) [GeV/c^{2}]");
} else if (title=="nVtx") {
  h_ratio->GetXaxis ()->SetTitle("Number of Vertices");
} else if (title=="leadElePt") {
  h_ratio->GetXaxis ()->SetTitle("Leading Electron p_{T} [GeV/c]");
} else if (title=="subleadElePt") {
  h_ratio->GetXaxis ()->SetTitle("Sub-Leading Electron p_{T} [GeV/c]");
} else if (title=="leadEleEta") {
  h_ratio->GetXaxis ()->SetTitle("Leading Electron #Eta");
} else if (title=="subleadEleEta") {
  h_ratio->GetXaxis ()->SetTitle("Sub-Leading Electron #Eta");
} else if (title=="FATjetPt") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet p_{T} [GeV/c]");
} else if (title=="FATjetEta") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet #Eta [GeV/c]");
} else if (title=="FATjetCISVV2") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet CSVv2 [GeV/c]");
} else if (title=="FATjetSDmass") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet Soft Drop Mass [GeV/c^{2}]");
} else if (title=="FATjetPRmass") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet Pruned Mass [GeV/c^{2}]");
} else if (title=="FATjetPRmassCorr") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet Corrected Pruned Mass [GeV/c^{2}]");
} else if (title=="FATjetTau1") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet #tau_{1}");
} else if (title=="FATjetTau2") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet #tau_{2}");
} else if (title=="FATjetTau2dvTau1") {
  h_ratio->GetXaxis ()->SetTitle("FAT Jet #tau_{21}");
} else if (title=="FATsubjetPt") {
  h_ratio->GetXaxis ()->SetTitle("FAT SubJet p_{T} [GeV/c]");
} else if (title=="FATsubjetEta") {
  h_ratio->GetXaxis ()->SetTitle("FAT SubJet #Eta [GeV/c]");
} else if (title=="FATsubjetSDCSV") {
  h_ratio->GetXaxis ()->SetTitle("FAT SubJet Soft Drop CSV [GeV/c]");
} else if (title=="ZHmass") {
  h_ratio->GetXaxis ()->SetTitle("ZH invariant mass [GeV/c^{2}]");
}

if (title=="ZHmass") { h_ratio->GetXaxis()->SetRangeUser(0, 3000); }
h_ratio->GetXaxis()->SetTitleSize(0.11);
h_ratio->GetXaxis()->SetLabelFont(42);
h_ratio->GetXaxis()->SetLabelSize(0.10);
h_ratio->GetXaxis()->SetTitleFont(42);
h_ratio->GetYaxis()->SetTitle("Data/MC");
h_ratio->GetYaxis()->SetNdivisions(505);
h_ratio->GetYaxis()->SetTitleSize(0.11);
h_ratio->GetYaxis()->SetLabelSize(0.10);
h_ratio->GetYaxis()->SetRangeUser(0.2, 2);
h_ratio->GetYaxis()->SetTitleOffset(0.33);

h_ratio->Divide(ht);
h_ratio->SetMarkerStyle(20);
h_ratio->Draw("E0PX0");

TLine *OLine = new TLine(h_ratio->GetXaxis()->GetXmin(),1.,h_ratio->GetXaxis()->GetXmax(),1.);
OLine->SetLineColor(kRed);
OLine->SetLineWidth(2);
OLine->Draw();

/*
TLatex *lar = new TLatex();
lar->SetNDC(kTRUE);
lar->SetTextSize(0.04);
lar->SetLineWidth(5);
lar->DrawLatex(0.14, 0.94, "CMS #it{#bf{2016}}");
lar->DrawLatex(0.60, 0.94, "L = 4.327 fb^{-1} at #sqrt{s} = 13 TeV");
*/

c1->cd();
c1->SaveAs((path + title + ".pdf").c_str());

}
