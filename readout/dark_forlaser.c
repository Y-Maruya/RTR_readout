#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TMath.h>
#include <TLine.h>
#include <TString.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TLatex.h>
using namespace std;
double x_ch[64] ={
  -3.5,-3.5,-3.5,-3.5,-2.5,-2.5,-2.5,-2.5,-1.5,-1.5,-1.5,-1.5,-0.5,-0.5,-0.5,-0.5,-3.5,-3.5,-3.5,-3.5,-2.5,-2.5,-2.5,-2.5,-1.5,-1.5,-1.5,-1.5,-0.5,-0.5,-0.5,-0.5,
  3.5,3.5,3.5,3.5,2.5,2.5,2.5,2.5,1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5,
  3.5,3.5,3.5,3.5,2.5,2.5,2.5,2.5,1.5,1.5,1.5,1.5,0.5,0.5,0.5,0.5};
double y_ch[64] ={
  2.5,3.5,1.5,0.5,3.5,2.5,0.5,1.5,2.5,3.5,1.5,0.5,3.5,2.5,0.5,1.5,
  -1.5,-0.5,-2.5,-3.5,-0.5,-1.5,-3.5,-2.5,-1.5,-0.5,-2.5,-3.5,-0.5,-1.5,-3.5,-2.5,
  1.5,0.5,2.5,3.5,0.5,1.5,3.5,2.5,1.5,0.5,2.5,3.5,0.5,1.5,3.5,2.5,
  -2.5,-3.5,-1.5,-0.5,-3.5,-2.5,-0.5,-1.5,-2.5,-3.5,-1.5,-0.5,-3.5,-2.5,-0.5,-1.5};


void val(string filename = "laser_5kHz30p8_HV56p22", string val_filename = "laser_HV56p24_903_2.root",int threshold = 850){
  gROOT->SetStyle("Plain");
  double mean = 0;
  double RMS = 0;
  
  std::ifstream hoge(filename,std::ios::binary);
  if(!hoge.is_open()){
    cout << "no file" << endl;
    return;
  }
  string outfile = filename + ".root";
  TFile*file = new TFile(outfile.c_str(),"RECREATE");
  TFile*val_file = new TFile(val_filename.c_str(),"READ");
  std::cout<< "val_file open" << std::endl;
  TGraphErrors*gains = (TGraphErrors*)val_file->Get("gain");
  TGraphErrors*biass = (TGraphErrors*)val_file->Get("bias");
  // TGraphErrors*sigmas = (TGraphErrors*)val_file->Get("sigma");
  // TGraphErrors*sigmas1 = (TGraphErrors*)val_file->Get("sigma1");
  if(gains == NULL || biass == NULL){
    std::cout<< "val_TGraph Error" << std::endl;
    return;
  }
  std::cout<< "val_TGraph Got" << std::endl;
  
  val_file->Close();
  file->cd();
  std::vector<Int_t> notgood;
  for(int i = 0; i<64; i++){
    if(std::isnan(biass->Eval(i))){
      notgood.push_back(i);
      std::cout<<"not good : "<< i<<std::endl;
    }
  }
  const int N_notgood = notgood.size();
  TH1F *hist[64];
  TH2F* chhist = new TH2F("chhist","",64,0-0.5,64-0.5,1500,500,2000);
  for(int i = 0; i<64; ++i){
    std::ostringstream _name;
    _name << "ch" << i;
    std::string hist_name = _name.str();
    hist[i] = new TH1F(hist_name.c_str(), hist_name.c_str(), 4096, 0, 4096);
    hist[i]->GetXaxis()->SetTitle("ADC count");
    hist[i]->GetYaxis()->SetTitle("Counts");
  }
  TH1D*sumhist= new TH1D("sumhist","sumhist",1000,0,20);
  while(!hoge.eof()){
    UInt_t val;
    hoge.read((char*)&val, sizeof(int));
    //    std::cout << std::hex << val << std::endl;
    double sum=0;
    if(val == 0xffffea0c){
      hoge.read((char*)&val, sizeof(int));
      for(int i = 0; i<65; ++i){
       	hoge.read((char*)&val, sizeof(int));
      	if(i>0){
        	int buffer = val & 0xffff;
         	hist[i-1]->Fill(buffer);
        	chhist->Fill(i-1,buffer);
          if((buffer-biass->Eval(i-1))/gains->Eval(i-1) > 0.5){
            sum+=(buffer-biass->Eval(i-1))/gains->Eval(i-1);
          }
      	}
      }
    }
    sumhist->Fill(sum);
  }
  // sumhist->Write();
  file->Write();
  // file->Close();
  TFile* new_val_file = new TFile(val_filename.c_str(),"READ");
  new_val_file->cd();
  TH1F* hist2[64];
  double n_L_true[64];
  double n_L_true_err[64];
  for(int i = 0; i<64; ++i){
    hist2[i] = (TH1F*)new_val_file->Get(Form("ch%d",i));
    if(hist2[i] == NULL){
      std::cout << "no hist" << i << std::endl;
    }
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    gStyle->SetOptFit(0000);
    gStyle->SetOptStat(0000);
    c1->Divide(1,2);
    c1->cd(1);
    hist[i]->GetXaxis()->SetRangeUser(700,1000);
    // hist[i]->Scale(1./hist[i]->GetMaximum());
    hist[i]->SetLineColor(2);
    hist[i]->SetLineWidth(2);
    hist[i]->Draw("hist");
    // hist2[i]->Scale(1./hist2[i]->GetMaximum());
    hist2[i]->SetLineColor(4);
    hist2[i]->SetLineWidth(2);
    hist2[i]->Draw("histsame");
    c1->cd(2);

    TH1F *hist3 = (TH1F*)hist2[i]->Clone();
    hist3->Add(hist[i],-1);
    // hist3->Divide(hist[i]);
    hist3->GetXaxis()->SetRangeUser(700,1000);
    hist3->SetLineColor(1);
    hist3->SetLineWidth(2);
    hist3->Draw("hist");
    string histname;
    if(i==0){histname = filename+".pdf(";}
    else if(i==63){histname = filename+".pdf)";}
    else{histname = filename+".pdf";}
    std::cout << histname << std::endl;
    c1->SaveAs(histname.c_str(),"pdf");
    int n_D;
    double err;
    hist[i]->GetXaxis()->SetRangeUser(0,4096);
    n_D=hist[i]->IntegralAndError(850,4096,err);
    int n_L;
    double err_L;
    hist2[i]->GetXaxis()->SetRangeUser(0,4096);
    n_L=hist2[i]->IntegralAndError(850,4096,err_L);
    std::cout << "ch" << i << " : " << n_D << " : " << n_L << std::endl;
    n_L_true[i] = (n_L - n_D/hist[i]->Integral(0,4096)*hist2[i]->Integral(0,4096))/hist2[i]->Integral(0,4096);
    n_L_true_err[i] = sqrt(err_L*err_L + pow(err/hist[i]->Integral(0,4096)*hist2[i]->Integral(0,4096),2))/hist2[i]->Integral(0,4096);
  }
  double x[64];
  double xerr[64];
  for(int i = 0; i<64; ++i){
    x[i] = i;
    xerr[i] = 0;
  }
  TGraphErrors* n_L_true_graph = new TGraphErrors(64,x,n_L_true,xerr,n_L_true_err);
  n_L_true_graph->SetName("p_Laser");
  n_L_true_graph->SetTitle("p_Laser");
  n_L_true_graph->GetXaxis()->SetTitle("Channel");
  n_L_true_graph->GetYaxis()->SetTitle("n_L_true");
  n_L_true_graph->SetMarkerStyle(20);
  n_L_true_graph->SetMarkerSize(1);
  n_L_true_graph->SetMarkerColor(2);
  TH2F*frame = new TH2F("2Dmap","",8,-4,4,8,-4,4);
  frame->GetXaxis()->SetTitle("x");
  frame->GetYaxis()->SetTitle("y");
  for(int i = 0; i<64; ++i){
    frame->Fill(x_ch[i],y_ch[i],n_L_true[i]);
  }
  frame->Draw("colz");
  file->cd();
  frame->Write();
  n_L_true_graph->Write();
  new_val_file->Close();
  file->Close();
  return;
}
int main(int argc, char* argv[]){
  if(argc == 3){ val(string(argv[1]),string(argv[2]));}
  else{ return 1;}
  return 0;
}