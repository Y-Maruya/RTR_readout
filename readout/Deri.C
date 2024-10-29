#define Deri_cxx
#include "Deri.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TFile.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TCanvas.h>


void Deri::Loop(config_data config)
{
//   In a ROOT session, you can do:
//      root> .L Deri.C
//      root> Deri t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   gROOT->SetStyle("ATLAS");
   Long64_t nentries = fChain->GetEntriesFast();
   TFile* before_laser_val_file = new TFile(config.before_laser_val_filename.c_str(),"READ");
   TFile* after_laser_val_file = new TFile(config.after_laser_val_filename.c_str(),"READ");
   TFile* before_dark_val_file = new TFile(config.before_dark_val_filename.c_str(),"READ");
   TFile* after_dark_val_file = new TFile(config.after_dark_val_filename.c_str(),"READ");
   if(!before_laser_val_file->IsOpen() || !after_laser_val_file->IsOpen() || !before_dark_val_file->IsOpen() || !after_dark_val_file->IsOpen()){
      std::cerr << "Error: some validation files are not opened" << std::endl;
      return;
   }
      before_laser_val_file->cd();
      TGraphErrors* before_laser_val_bias = (TGraphErrors*)before_laser_val_file->Get("bias");
      TGraphErrors* before_laser_val_gain = (TGraphErrors*)before_laser_val_file->Get("gain");
      TGraphErrors* before_laser_val_sigma0 = (TGraphErrors*)before_laser_val_file->Get("sigma0");
      TGraphErrors* before_laser_val_sigma1 = (TGraphErrors*)before_laser_val_file->Get("sigma1");
      double before_cut_value[64];
      for(int i = 0; i<64; ++i){
         before_cut_value[i] = before_laser_val_bias->Eval(i)+before_laser_val_gain->Eval(i)*0.5;
      }
      after_laser_val_file->cd();
      TGraphErrors* after_laser_val_bias = (TGraphErrors*)after_laser_val_file->Get("bias");
      TGraphErrors* after_laser_val_gain = (TGraphErrors*)after_laser_val_file->Get("gain");
      TGraphErrors* after_laser_val_sigma0 = (TGraphErrors*)after_laser_val_file->Get("sigma0");
      TGraphErrors* after_laser_val_sigma1 = (TGraphErrors*)after_laser_val_file->Get("sigma1");
      double after_cut_value[64];
      for(int i = 0; i<64; ++i){
         after_cut_value[i] = after_laser_val_bias->Eval(i)+after_laser_val_gain->Eval(i)*0.5;
      }
      double x_ch[64];
      double x_ch_err[64] = {0};
      for(int i = 0; i<64; ++i){
         x_ch[i] = i;
      }
      TGraph* before_cut_value_graph = new TGraph(64,x_ch,before_cut_value);
      TGraph* after_cut_value_graph = new TGraph(64,x_ch,after_cut_value);
      TCanvas* c0 = new TCanvas("c0","c0",800,800);
      before_cut_value_graph->SetTitle("cut value before;channel;cut value");
      before_cut_value_graph->SetMarkerStyle(20);
      before_cut_value_graph->SetMarkerSize(1.0);
      before_cut_value_graph->SetMarkerColor(kRed);
      before_cut_value_graph->Draw("AP");
      after_cut_value_graph->SetTitle("cut value after;channel;cut value");
      after_cut_value_graph->SetMarkerStyle(20);
      after_cut_value_graph->SetMarkerSize(1.0);
      after_cut_value_graph->SetMarkerColor(kBlue);
      after_cut_value_graph->Draw("P");
      c0->BuildLegend();
      c0->Write();
      std::string filename0 = config.easiroc_data_filename+"_cut_value.png";
      c0->SaveAs(filename0.c_str());
      filename0 = config.easiroc_data_filename+"_cut_value.root";
      c0->SaveAs(filename0.c_str());
      before_dark_val_file->cd();
      double before_dark_rate_err[64] = {0};
      double before_dark_rate[64] = {0};
      for(int i = 0; i<64; ++i){
         TH1F* before_dark_val_hist_tmp = (TH1F*)before_dark_val_file->Get(Form("ch%d",i));
         if(before_dark_val_hist_tmp == NULL){
            std::cerr << "Error: before_dark_val_hist_tmp is NULL" << std::endl;
            return;
         }
         before_dark_rate[i] = before_dark_val_hist_tmp->IntegralAndError(before_dark_val_hist_tmp->FindBin(before_cut_value[i]),before_dark_val_hist_tmp->FindBin(4096),before_dark_rate_err[i])/before_dark_val_hist_tmp->Integral(0,4096);
         before_dark_rate_err[i] /= before_dark_val_hist_tmp->Integral(0,4096);
         before_dark_val_hist_tmp->Delete();
      }
      after_dark_val_file->cd();
      double after_dark_rate_err[64] = {0};
      double after_dark_rate[64] = {0};
      for(int i = 0; i<64; ++i){
         TH1F* after_dark_val_hist_tmp = (TH1F*)after_dark_val_file->Get(Form("ch%d",i));
         if(after_dark_val_hist_tmp == NULL){
            std::cerr << "Error: after_dark_val_hist_tmp is NULL" << std::endl;
            return;
         }
         after_dark_rate[i] = after_dark_val_hist_tmp->IntegralAndError(after_dark_val_hist_tmp->FindBin(after_cut_value[i]),after_dark_val_hist_tmp->FindBin(4096),after_dark_rate_err[i])/after_dark_val_hist_tmp->Integral(0,4096);
         after_dark_rate_err[i] /= after_dark_val_hist_tmp->Integral(0,4096);
         after_dark_val_hist_tmp->Delete();
      }
      TCanvas* c1 = new TCanvas("c1","c1",800,800);
      TGraphErrors* before_dark_rate_graph = new TGraphErrors(64,x_ch,before_dark_rate,x_ch_err,before_dark_rate_err);
      before_dark_rate_graph->SetTitle("dark rate before;channel;rate");
      before_dark_rate_graph->SetMarkerStyle(20);
      before_dark_rate_graph->SetMarkerSize(1.0);
      before_dark_rate_graph->SetMarkerColor(kRed);
      before_dark_rate_graph->Draw("AP");
      TGraphErrors* after_dark_rate_graph = new TGraphErrors(64,x_ch,after_dark_rate,x_ch_err,after_dark_rate_err);
      after_dark_rate_graph->SetTitle("dark rate after;channel;rate");
      after_dark_rate_graph->SetMarkerStyle(20);
      after_dark_rate_graph->SetMarkerSize(1.0);
      after_dark_rate_graph->SetMarkerColor(kBlue);
      after_dark_rate_graph->Draw("P");
      // c1->Write();
      c1->BuildLegend();
      std::string filename = config.easiroc_data_filename+"_dark_rate.png";
      c1->SaveAs(filename.c_str());
      filename = config.easiroc_data_filename+"_dark_rate.root";
      c1->SaveAs(filename.c_str());
   Long64_t nbytes = 0, nb = 0;
   int overcut_before[64] = {0};
   int overcut_after[64] = {0};
   double rate_before[64] = {0};
   double rate_after[64] = {0};
   double rate_before_err[64] = {0};
   double rate_after_err[64] = {0};
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%1000 == 0) std::cout << "Processing event " << jentry << std::endl;
      for(int i = 0; i<64; ++i){
         if(adc[i] > before_cut_value[i]){
            overcut_before[i]++;
         }
         if(adc[i] > after_cut_value[i]){
            overcut_after[i]++;
         }
      } 
   }
   for(int i = 0; i<64; ++i){
      rate_before[i] = (double)overcut_before[i]/nentries-before_dark_rate[i];
      rate_after[i] = (double)overcut_after[i]/nentries-after_dark_rate[i];
      rate_before_err[i] = std::hypot(sqrt((double)overcut_before[i])/nentries,before_dark_rate_err[i]);
      rate_after_err[i] = std::hypot(sqrt((double)overcut_after[i])/nentries,after_dark_rate_err[i]);
   }
   TCanvas* c2 = new TCanvas("c2","c2",800,800);
   TGraphErrors* rate_before_graph = new TGraphErrors(64,x_ch,rate_before,x_ch_err,rate_before_err);
   TGraphErrors* rate_after_graph = new TGraphErrors(64,x_ch,rate_after,x_ch_err,rate_after_err);
   rate_before_graph->SetTitle("rate before;channel;rate");
   rate_before_graph->SetMarkerStyle(20);
   rate_before_graph->SetMarkerSize(1.0);
   rate_before_graph->SetMarkerColor(kRed);
   rate_before_graph->Draw("AP");
   rate_after_graph->SetTitle("rate after;channel;rate");
   rate_after_graph->SetMarkerStyle(20);
   rate_after_graph->SetMarkerSize(1.0);
   rate_after_graph->SetMarkerColor(kBlue);
   rate_after_graph->Draw("P");
   c2->BuildLegend();
   std::string filename2 = config.easiroc_data_filename+"_rate.png";
   c2->SaveAs(filename2.c_str());
   before_laser_val_file->Close();
   after_laser_val_file->Close();
   before_dark_val_file->Close();
   after_dark_val_file->Close();
   std::string filename3 = config.easiroc_data_filename+"_rate.root";
   TFile* outfile = new TFile(filename3.c_str(),"RECREATE");
   outfile->cd();
   c1->Write();
   c2->Write();
   before_dark_rate_graph->Write();
   after_dark_rate_graph->Write();
   rate_before_graph->Write();
   rate_after_graph->Write();
   outfile->Close();
}

int main(int argc, char* argv[]){
   if(argc != 2){
      std::cout << "Usage: " << argv[0] << " <config file name>" << std::endl;
      return 1;
   }
   config_data config;
   config.read_config(std::string(argv[1]));
   std::string filename = config.easiroc_data_filename+".root";
   Deri *d = new Deri(filename, config.fCAEN);
   d->Loop(config);
   return 0;
}