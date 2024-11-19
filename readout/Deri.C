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
double cut = 6;
double x_ch_hist[64] ={
  -3.5,-3.5,-3.5,-3.5,-2.5,-2.5,-2.5,-2.5,-1.5,-1.5,-1.5,-1.5,-0.5,-0.5,-0.5,-0.5,-3.5,-3.5,-3.5,-3.5,-2.5,-2.5,-2.5,-2.5,-1.5,-1.5,-1.5,-1.5,-0.5,-0.5,-0.5,-0.5,
  3.5,3.5,3.5,3.5,2.5,2.5,2.5,2.5,1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,
  3.5,3.5,3.5,3.5,2.5,2.5,2.5,2.5,1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5};
double y_ch_hist[64] ={
  2.5,3.5,1.5,0.5,3.5,2.5,0.5,1.5,2.5,3.5,1.5,0.5,3.5,2.5,0.5,1.5,
  -1.5,-0.5,-2.5,-3.5,-0.5,-1.5,-3.5,-2.5,-1.5,-0.5,-2.5,-3.5,-0.5,-1.5,-3.5,-2.5,
  1.5,0.5,2.5,3.5,0.5,1.5,3.5,2.5,1.5,0.5,2.5,3.5,0.5,1.5,3.5,2.5,
  -2.5,-3.5,-1.5,-0.5,-3.5,-2.5,-0.5,-1.5,-2.5,-3.5,-1.5,-0.5,-3.5,-2.5,-0.5,-1.5};

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
   // gROOT->SetStyle("ATLAS");
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
      TH2D* before_dark_rate_hist = new TH2D("before_dark_rate_hist","before_dark_rate_hist",8,-4,4,8,-4,4);
      for(int i =0;i<64;++i){
         if(i==44){
            continue;
         }
         if(i == 29){
            continue;
         }
         before_dark_rate_hist->Fill(x_ch_hist[i],y_ch_hist[i],before_dark_rate[i]);
      }
      TGraphErrors* before_dark_rate_graph = new TGraphErrors(64,x_ch,before_dark_rate,x_ch_err,before_dark_rate_err);
      before_dark_rate_graph->SetTitle("dark rate before;channel;rate");
      before_dark_rate_graph->SetMarkerStyle(20);
      before_dark_rate_graph->SetMarkerSize(1.0);
      before_dark_rate_graph->SetMarkerColor(kRed);
      before_dark_rate_graph->Draw("AP");
      TH2D* after_dark_rate_hist = new TH2D("after_dark_rate_hist","after_dark_rate_hist",8,-4,4,8,-4,4);
      for(int i =0;i<64;++i){
         if(i==44){
            continue;
         }
         if(i == 29){
            continue;
         }
         after_dark_rate_hist->Fill(x_ch_hist[i],y_ch_hist[i],after_dark_rate[i]);
      }
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
      TCanvas* c2 = new TCanvas("c2","c2",800,800);
      before_dark_rate_hist->Draw("colz");
      filename = config.easiroc_data_filename+"_before_dark_rate_hist.png";
      c2->SaveAs(filename.c_str());
      filename = config.easiroc_data_filename+"_before_dark_rate_hist.root";
      c2->SaveAs(filename.c_str());
      TCanvas* c3 = new TCanvas("c3","c3",800,800);
      after_dark_rate_hist->Draw("colz");
      filename = config.easiroc_data_filename+"_after_dark_rate_hist.png";
      c3->SaveAs(filename.c_str());
      filename = config.easiroc_data_filename+"_after_dark_rate_hist.root";
      c3->SaveAs(filename.c_str());
   Long64_t nbytes = 0, nb = 0;
   int overcut_before[64] = {0};
   int overcut_after[64] = {0};
   int overcut_before_withcut[64] = {0};
   int overcut_after_withcut[64] = {0};
   double rate_before[64] = {0};
   double rate_after[64] = {0};
   double rate_before_withcut[64] = {0};
   double rate_after_withcut[64] = {0};
   double rate_before_err[64] = {0};
   double rate_after_err[64] = {0};
   double rate_before_withcut_err[64] = {0};
   double rate_after_withcut_err[64] = {0};
   TH1D*before_sumhist = new TH1D("before_sumhist","before_sumhist",300,0,300);
   TH1D*after_sumhist = new TH1D("after_sumhist","after_sumhist",300,0,300);
   double before_sum_ch[64] = {0};
   double after_sum_ch[64] = {0};
   double before_sum_ch_withcut[64] = {0};
   double after_sum_ch_withcut[64] = {0};
   double n_event_before=0;
   double n_event_after=0;
   double n_event_before_withcut=0;
   double n_event_after_withcut=0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      if(jentry%1000 == 0) std::cout << "Processing event " << jentry << std::endl;
      double tmp_before=0;
      double tmp_after=0;
      for(int i = 0; i<64; ++i){
         if(adc[i] > before_cut_value[i]){
            overcut_before[i]++;
            tmp_before+= (adc[i]-before_laser_val_bias->Eval(i))/before_laser_val_gain->Eval(i);
            // before_sumhist->Fill((adc[i]-before_laser_val_bias->Eval(i))/before_laser_val_gain->Eval(i));
         }
         if(adc[i] > after_cut_value[i]){
            overcut_after[i]++;
            tmp_after+= (adc[i]-after_laser_val_bias->Eval(i))/after_laser_val_gain->Eval(i);
            // after_sumhist->Fill((adc[i]-after_laser_val_bias->Eval(i))/after_laser_val_gain->Eval(i));
         }
      }
      if(tmp_before < cut){
         for(int i = 0; i<64; ++i){
            if(adc[i] > before_cut_value[i]){
               overcut_before_withcut[i]++;
            }
         }
         n_event_before_withcut+=1;
      }
      if(tmp_after < cut){
         for(int i = 0; i<64; ++i){
            if(adc[i] > after_cut_value[i]){
               overcut_after_withcut[i]++;
            }
         }
         n_event_after_withcut+=1;
      }
      before_sumhist->Fill(tmp_before);
      after_sumhist->Fill(tmp_after);
      for(int i =0; i<64; ++i){
         before_sum_ch[i] += (adc[i]-before_laser_val_bias->Eval(i))/before_laser_val_gain->Eval(i);
         if(tmp_before < cut ){
            before_sum_ch_withcut[i] += (adc[i]-before_laser_val_bias->Eval(i))/before_laser_val_gain->Eval(i);
         }
         after_sum_ch[i] += (adc[i]-after_laser_val_bias->Eval(i))/after_laser_val_gain->Eval(i);
         if(tmp_after<cut){
            after_sum_ch_withcut[i] += (adc[i]-after_laser_val_bias->Eval(i))/after_laser_val_gain->Eval(i);
         }
      }
   }
   for(int i=0; i<64; ++i){
      std::cout<< before_sum_ch[i] << " " << n_event_before << std::endl;
      before_sum_ch[i] /= n_event_before;
      std::cout<< before_sum_ch_withcut[i] << " " << n_event_before_withcut << std::endl;
      before_sum_ch_withcut[i] /= n_event_before_withcut;
      std::cout<< after_sum_ch[i] << " " << n_event_after << std::endl;
      after_sum_ch[i] /= n_event_after;
      std::cout<< after_sum_ch_withcut[i] << " " << n_event_after_withcut << std::endl;
      after_sum_ch_withcut[i] /= n_event_after_withcut;
   }
   TH2D* before_sum_ch_hist = new TH2D("before_sum_ch_hist","before_sum_ch_hist",8,-4,4,8,-4,4);
   TH2D* before_sum_ch_withcut_hist = new TH2D("before_sum_ch_withcut_hist","before_sum_ch_withcut_hist",8,-4,4,8,-4,4);
   TH2D* after_sum_ch_hist = new TH2D("after_sum_ch_hist","after_sum_ch_hist",8,-4,4,8,-4,4);
   TH2D* after_sum_ch_withcut_hist = new TH2D("after_sum_ch_withcut_hist","after_sum_ch_withcut_hist",8,-4,4,8,-4,4);
   for(int i = 0; i<64; ++i){
      if(i==44){
         continue;
      }
      if(i == 29){
         continue;
      }
      std::cout<<before_sum_ch[i]<<std::endl;
      std::cout<<after_sum_ch[i]<<std::endl;
      std::cout<<x_ch_hist[i]<<" "<<y_ch_hist[i]<<std::endl;
      before_sum_ch_hist->Fill(x_ch_hist[i],y_ch_hist[i],before_sum_ch[i]);
      before_sum_ch_withcut_hist->Fill(x_ch_hist[i],y_ch_hist[i],before_sum_ch_withcut[i]);
      after_sum_ch_hist->Fill(x_ch_hist[i],y_ch_hist[i],after_sum_ch[i]);
      after_sum_ch_withcut_hist->Fill(x_ch_hist[i],y_ch_hist[i],after_sum_ch_withcut[i]);
   }
   TCanvas* c31 = new TCanvas("c31","c31",800,800);
   before_sum_ch_hist->Draw("colz");
   std::string filename1 = config.easiroc_data_filename+"_before_sum_ch.png";
   c31->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_before_sum_ch.root";
   c31->SaveAs(filename1.c_str());
   TCanvas* c4 = new TCanvas("c4","c4",800,800);
   before_sum_ch_withcut_hist->Draw("colz");
   filename1 = config.easiroc_data_filename+"_before_sum_ch_withcut.png";
   c4->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_before_sum_ch_withcut.root";
   c4->SaveAs(filename1.c_str());
   TCanvas* c5 = new TCanvas("c5","c5",800,800);
   after_sum_ch_hist->Draw("colz");
   filename1 = config.easiroc_data_filename+"_after_sum_ch.png";
   c5->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_after_sum_ch.root";
   c5->SaveAs(filename1.c_str());
   TCanvas* c6 = new TCanvas("c6","c6",800,800);
   after_sum_ch_withcut_hist->Draw("colz");
   filename1 = config.easiroc_data_filename+"_after_sum_ch_withcut.png";
   c6->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_after_sum_ch_withcut.root";
   c6->SaveAs(filename1.c_str());
   TCanvas*c7 = new TCanvas("c7","c7",800,800);
   before_sumhist->Draw();
   filename1 = config.easiroc_data_filename+"_before_sum.png";
   c7->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_before_sum.root";
   c7->SaveAs(filename1.c_str());
   TCanvas*c8 = new TCanvas("c8","c8",800,800);
   after_sumhist->Draw();
   filename1 = config.easiroc_data_filename+"_after_sum.png";
   c8->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_after_sum.root";
   c8->SaveAs(filename1.c_str());
   filename1 = config.easiroc_data_filename+"_sumhist.root";
   TFile* outfile_forsum = new TFile(filename1.c_str(),"RECREATE");
   outfile_forsum->cd();
   before_sumhist->Write();
   after_sumhist->Write();
   before_sum_ch_hist->Write();
   before_sum_ch_withcut_hist->Write();
   after_sum_ch_hist->Write();
   after_sum_ch_withcut_hist->Write();
   outfile_forsum->Close();
   std::cout<<overcut_after_withcut[0]<<std::endl;
   std::cout<<n_event_after_withcut<<std::endl;
   for(int i = 0; i<64; ++i){
      rate_before[i] = (double)overcut_before[i]/nentries-before_dark_rate[i];
      rate_after[i] = (double)overcut_after[i]/nentries-after_dark_rate[i];
      rate_before_withcut[i] = (double)overcut_before_withcut[i]/n_event_before_withcut - before_dark_rate[i];
      rate_after_withcut[i] = (double)overcut_after_withcut[i]/n_event_after_withcut - after_dark_rate[i];
      rate_before_err[i] = std::hypot(sqrt((double)overcut_before[i])/nentries,before_dark_rate_err[i]);
      rate_after_err[i] = std::hypot(sqrt((double)overcut_after[i])/nentries,after_dark_rate_err[i]);
      rate_before_withcut_err[i] = std::hypot(sqrt((double)overcut_before_withcut[i])/n_event_before_withcut,before_dark_rate_err[i]);
      rate_after_withcut_err[i] =std::hypot(sqrt((double)overcut_after_withcut[i])/n_event_after_withcut,after_dark_rate_err[i]);
   }
   TH2D* rate_before_hist = new TH2D("rate_before_hist","rate_before_hist",8,-4,4,8,-4,4);
   TH2D* rate_after_hist = new TH2D("rate_after_hist","rate_after_hist",8,-4,4,8,-4,4);
   TH2D* rate_before_withcut_hist = new TH2D("rate_before_withcut_hist","rate_before_withcut_hist",8,-4,4,8,-4,4);
   TH2D* rate_after_withcut_hist = new TH2D("rate_after_withcut_hist","rate_after_withcut_hist",8,-4,4,8,-4,4);
   for(int i =0;i<64;++i){
      if(i==44){
         continue;
      }
      if(i == 29){
         continue;
      }
      rate_before_hist->Fill(x_ch_hist[i],y_ch_hist[i],rate_before[i]);
      rate_after_hist->Fill(x_ch_hist[i],y_ch_hist[i],rate_after[i]);
      rate_before_withcut_hist->Fill(x_ch_hist[i],y_ch_hist[i],rate_before_withcut[i]);
      rate_after_withcut_hist->Fill(x_ch_hist[i],y_ch_hist[i],rate_after_withcut[i]);
   }
   TCanvas* c10 = new TCanvas("c10","c10",800,800);
   rate_before_hist->Draw("colz");
   std::string filename10 = config.easiroc_data_filename+"_rate_before_hist.png";
   c10->SaveAs(filename10.c_str());
   filename10 = config.easiroc_data_filename+"_rate_before_hist.root";
   c10->SaveAs(filename10.c_str());
   TCanvas* c11 = new TCanvas("c11","c11",800,800);
   rate_after_hist->Draw("colz");
   std::string filename11 = config.easiroc_data_filename+"_rate_after_hist.png";
   c11->SaveAs(filename11.c_str());
   filename11 = config.easiroc_data_filename+"_rate_after_hist.root";
   c11->SaveAs(filename11.c_str());
   TCanvas* c12 = new TCanvas("c12","c12",800,800);
   rate_before_withcut_hist->Draw("colz");
   std::string filename12 = config.easiroc_data_filename+"_rate_before_withcut_hist.png";
   c12->SaveAs(filename12.c_str());
   filename12 = config.easiroc_data_filename+"_rate_before_withcut_hist.root";
   c12->SaveAs(filename12.c_str());
   TCanvas* c13 = new TCanvas("c13","c13",800,800);
   rate_after_withcut_hist->Draw("colz");
   std::string filename13 = config.easiroc_data_filename+"_rate_after_withcut_hist.png";
   c13->SaveAs(filename13.c_str());
   filename13 = config.easiroc_data_filename+"_rate_after_withcut_hist.root";
   c13->SaveAs(filename13.c_str());
   TCanvas* c20 = new TCanvas("c20","c20",800,800);
   TGraphErrors* rate_before_withcut_graph = new TGraphErrors(64,x_ch,rate_before_withcut,x_ch_err,rate_before_withcut_err);
   TGraphErrors* rate_after_withcut_graph = new TGraphErrors(64,x_ch,rate_after_withcut,x_ch_err,rate_after_withcut_err);
   rate_before_withcut_graph->SetTitle("rate before with cut;channel;rate");
   rate_before_withcut_graph->SetMarkerStyle(20);
   rate_before_withcut_graph->SetMarkerSize(1.0);
   rate_before_withcut_graph->SetMarkerColor(kRed);
   rate_before_withcut_graph->Draw("AP");
   rate_after_withcut_graph->SetTitle("rate after with cut;channel;rate");
   rate_after_withcut_graph->SetMarkerStyle(20);
   rate_after_withcut_graph->SetMarkerSize(1.0);
   rate_after_withcut_graph->SetMarkerColor(kBlue);
   rate_after_withcut_graph->Draw("P");
   c20->BuildLegend();
   std::string filename20 = config.easiroc_data_filename+"_rate_withcut.png";
   c20->SaveAs(filename20.c_str());
   filename20 = config.easiroc_data_filename+"_rate_withcut.root";
   c20->SaveAs(filename20.c_str());
   TCanvas* c21 = new TCanvas("c21","c21",800,800);
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
   c21->BuildLegend();
   std::string filename2 = config.easiroc_data_filename+"_rate.png";
   c21->SaveAs(filename2.c_str());
   before_laser_val_file->Close();
   after_laser_val_file->Close();
   before_dark_val_file->Close();
   after_dark_val_file->Close();
   std::string filename3 = config.easiroc_data_filename+"_rate.root";
   TFile* outfile = new TFile(filename3.c_str(),"RECREATE");
   outfile->cd();
   c1->Write();
   c21->Write();
   before_dark_rate_graph->Write();
   after_dark_rate_graph->Write();
   rate_before_graph->Write();
   rate_after_graph->Write();
   rate_before_hist->Write();
   rate_after_hist->Write();
   rate_before_withcut_graph->Write();
   rate_after_withcut_graph->Write();
   rate_before_withcut_hist->Write();
   rate_after_withcut_hist->Write();

   outfile->Close();
}
void val(string filename = "laser_5kHz30p8_HV56p22", string valfilename = "laser_HV56p24_3.root"){
  gROOT->SetStyle("Plain");
  
  double mean = 0;
  double RMS = 0;
  
  std::ifstream hoge(filename,std::ios::binary);
  if(!hoge.is_open()){
    cout << "no file" << endl;
    return;
  }
  string outfile = filename + "_sumhist.root";
  TFile*file = new TFile(outfile.c_str(),"RECREATE");
  std::cout<< "file open" << std::endl;
  TFile*val_file = new TFile(valfilename.c_str(),"READ");
  std::cout<< "val_file open" << std::endl;
  TGraphErrors*gains = (TGraphErrors*)val_file->Get("gain");
  TGraphErrors*biass = (TGraphErrors*)val_file->Get("bias");
  if(gains == NULL || biass == NULL){
    std::cout<< "val_TGraph Error" << std::endl;
    return;
  }
  std::cout<< "val_TGraph Got" << std::endl;
  file->cd();
  TH1F *hist[64];
  TH1D* sumhist = new TH1D("sumhist","sumhist",300,0,300);
  TH2F* chhist = new TH2F("chhist","",64,0-0.5,64-0.5,1500,500,2000);
  std::vector<TH2D*> chhist2;
  for(int i = 0; i<64; ++i){
    std::ostringstream _name;
    _name << "ch" << i;
    std::string hist_name = _name.str();
    hist[i] = new TH1F(hist_name.c_str(), hist_name.c_str(), 4096, 0, 4096);
    hist[i]->GetXaxis()->SetTitle("ADC count");
    hist[i]->GetYaxis()->SetTitle("Counts");
  }
  int event_num =0;
  while(!hoge.eof()&& event_num<10000){
    UInt_t val;
    hoge.read((char*)&val, sizeof(int));
    //    std::cout << std::hex << val << std::endl;
    
    if(val == 0xffffea0c){
      hoge.read((char*)&val, sizeof(int));
      std::cout<< Form("event %d",event_num) << std::endl;
      TH2D*chhist2_tmp = new TH2D(Form("event_%d",event_num),Form("event_%d",event_num),8,-4,4,8,-4,4);
      for(int i = 0; i<65; ++i){
       	hoge.read((char*)&val, sizeof(int));
      	if(i>0){
        	int buffer = val & 0xffff;
          if(i-1 == 44){
            continue;
          }
         if(i-1 == 29){
            continue;
          }
         	hist[i-1]->Fill(buffer);
        	chhist->Fill(i-1,buffer);
          if((buffer-biass->Eval(i-1))/gains->Eval(i-1) > 0.5){
            chhist2_tmp->Fill(x_ch_hist[i-1],y_ch_hist[i-1],(buffer-biass->Eval(i-1))/gains->Eval(i-1));
          }
      	}
      }
      chhist2.push_back(chhist2_tmp);
      std::cout<< chhist2_tmp->Integral() << std::endl;
      sumhist->Fill(chhist2_tmp->Integral());
      event_num++;
    }
  }
  int NMAX=chhist2.size();
  // if(NMAX == 0){
  //   std::cout<< "no event" << std::endl;
  //   return;
  // }else if(NMAX == 1){
  //   std::cout<< "only 1 event" << std::endl;
  //   chhist2[0]->Write();
  //   string histname;
  //   gStyle->SetOptStat(0);
  //   TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  //         TText *text = new TText(0.7,0.93,Form("sum=%g",chhist2[0]->Integral()));
  //     text->SetNDC();
  //     text->SetTextSize(0.04);
  //     chhist2[0]->Draw("colz");
  //     text->Draw("same");
  //   chhist2[0]->GetZaxis()->SetRangeUser(0,30);
  //   histname = filename+"_2D.pdf";
  //   c1->SaveAs(histname.c_str(),"pdf");
  // }else{
  //   for(int i = 0; i< NMAX; i++){
  //     chhist2[i]->Write();
  //     string histname;
  //     gStyle->SetOptStat(0);
  //     TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  //     TText *text = new TText(0.7,0.93,Form("sum=%g",chhist2[i]->Integral()));
  //     text->SetNDC();
  //     text->SetTextSize(0.04);
  //     chhist2[i]->Draw("colz");
  //     text->Draw("same");
  //     chhist2[i]->GetZaxis()->SetRangeUser(0,30);
  //     if(i==0){histname = filename+"_2D.pdf(";}
  //     else if(i==NMAX-1){histname = filename+"_2D.pdf)";}
  //     else{histname = filename+"_2D.pdf";}
  //     c1->SaveAs(histname.c_str(),"pdf");
  //   }
  // }
  TCanvas*c2 = new TCanvas("c2","c2",800,800);
  sumhist->Draw();
  std::string histnames = filename+"_sumhist.pdf";
  c2->SaveAs(histnames.c_str(),"pdf");
  // file->Write();
  sumhist->Write();
  file->Close();
  val_file->Close();
  return;
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
   // val(config.easiroc_data_filename,config.before_laser_val_filename);
   return 0;
}