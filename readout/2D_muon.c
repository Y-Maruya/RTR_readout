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

void val(string filename = "laser_5kHz30p8_HV56p22",int num = 1000000, int trigger = 100, int time = 47){
  gROOT->SetStyle("Plain");
  string outfile = filename + ".root";
  TFile*file = new TFile(outfile.c_str(),"RECREATE");
  std::cout<< "file open" << std::endl;
  TFile*val_file = new TFile("/rhome/maruya/disklink/MPPC_kousei/data/ak_laser_kousei_HV56p24_3.root","READ");
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
  for(int i = 0; i< num; i++){
    double mean = 0;
    double RMS = 0;
    string filename_string = filename + "/"+to_string(i)+".dat";
    string filename_string2 = filename + "/"+to_string(i)+".txt";
    std::ifstream hoge(filename_string,std::ios::binary);
    std::ifstream hoge2(filename_string2);
    if(!hoge.is_open() || !hoge2.is_open()){
      cout << "no file" << endl;
      continue;
    }
    std::string buff;
    std::string time_mu_string;
    std::getline(hoge2,buff);
    std::getline(hoge2,buff);
    std::getline(hoge2,buff);
    std::getline(hoge2,buff);
    std::getline(hoge2,buff);
    std::getline(hoge2,buff);
    std::getline(hoge2,time_mu_string);
    std::getline(hoge2,buff);
    hoge2.close();
    time_mu = std::stoi(time_mu_string);
    if(time_mu/trigger == time){
      continue;
    }
    while(!hoge.eof()){
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
         	  hist[i-1]->Fill(buffer);
        	  chhist->Fill(i-1,buffer);
            if((buffer-biass->Eval(i-1))/gains->Eval(i-1) > 0.5){
              chhist2_tmp->Fill(x_ch[i-1],y_ch[i-1],(buffer-biass->Eval(i-1))/gains->Eval(i-1));
            }
      	  }
        }
        chhist2.push_back(chhist2_tmp);
        event_num++;
      }
    }
  }
  int NMAX=chhist2.size();
  if(NMAX == 0){
    std::cout<< "no event" << std::endl;
    return;
  }else if(NMAX == 1){
    std::cout<< "only 1 event" << std::endl;
    chhist2[0]->Write();
    string histname;
    gStyle->SetOptStat(0);
    TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
    chhist2[0]->Draw("colz");
    chhist2[0]->GetZaxis()->SetRangeUser(0,30);
    histname = filename+"_2D.pdf";
    c1->SaveAs(histname.c_str(),"pdf");
  }else{
    for(int i = 0; i< NMAX; i++){
      chhist2[i]->Write();
      string histname;
      gStyle->SetOptStat(0);
      TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
      chhist2[i]->Draw("colz");
      chhist2[i]->GetZaxis()->SetRangeUser(0,30);
      if(i==0){histname = filename+"_2D.pdf(";}
      else if(i==NMAX-1){histname = filename+"_2D.pdf)";}
      else{histname = filename+"_2D.pdf";}
      c1->SaveAs(histname.c_str(),"pdf");
    }
  }
  file->Write();
  file->Close();
  val_file->Close();
  return;
}
int main(int argc, char* argv[]){
  if(argc == 5){ val(string(argv[1]),std::stoi(argv[2]),std::stoi(argv[3])),std::stoi(argv[4]);}
  else{ return 1;}
  return 0;
}