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


void val(string filename = "laser_5kHz30p8_HV56p22", string val_filename = "laser_HV56p24_903_2.root"){
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
  sumhist->Write();
  file->Write();
  file->Close();
  val_file->cd();
  TH1F* hist2[64];
  for(int i = 0; i<64; ++i){
    hist2[i] = (TH1F*)val_file->Get(Form("ch%d",i));
    if(hist2[i] == NULL){
      std::cout << "no hist" << i << std::endl;
    }
    TCanvas *c1 = new TCanvas("c1","c1",800,800);
    hist[i]->GetXaxis()->SetRangeUser(700,1000);
    hist[i]->Scale(1/hist[i]->GetMaximum());
    hist[i]->SetLineColor(2);
    hist[i]->SetLineWidth(2);
    hist[i]->Draw("hist");
    hist2[i]->Scale(1/hist2[i]->GetMaximum());
    hist2[i]->SetLineColor(4);
    hist2[i]->SetLineWidth(2);
    hist2[i]->Draw("histsame");
    string histname;
    if(i==0){histname = filename+".pdf(";}
    else if(i==63){histname = filename+".pdf)";}
    else{histname = filename+".pdf";}
    c1->SaveAs(histname.c_str(),"pdf"); 
  }
  val_file->Close();
  return;
}
int main(int argc, char* argv[]){
  if(argc == 3){ val(string(argv[1]),string(argv[2]));}
  else{ return 1;}
  return 0;
}