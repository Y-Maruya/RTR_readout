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


void val(string filename = "laser_5kHz30p8_HV56p22"){
  gROOT->SetStyle("Plain");
  std::vector<int> notgood;
  double mean = 0;
  double RMS = 0;
  
  std::ifstream hoge(filename,std::ios::binary);
  if(!hoge.is_open()){
    cout << "no file" << endl;
    return;
  }
  string outfile = filename + ".root";
  TFile*file = new TFile(outfile.c_str(),"RECREATE");
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

  while(!hoge.eof()){
    UInt_t val;
    hoge.read((char*)&val, sizeof(int));
    //    std::cout << std::hex << val << std::endl;
    
    if(val == 0xffffea0c){
      hoge.read((char*)&val, sizeof(int));
      for(int i = 0; i<65; ++i){
       	hoge.read((char*)&val, sizeof(int));
      	if(i>0){
        	int buffer = val & 0xffff;
         	hist[i-1]->Fill(buffer);
        	chhist->Fill(i-1,buffer);
      	}
      }
    }
  }
  double x[64];
  double xe[64];
  double bias[64];
  double biaserr[64];
  double gain[64];
  double gainerr[64];
  double sigma0[64];
  double sigma0err[64];
  double sigma1[64];
  double sigma1err[64];
  for(int i=0;i<64;i++){
    int Max_bin = hist[i]->GetMaximumBin();
    int second_Max_bin=0;
    for(int j=0;j<25;j++){
      hist[i]->GetXaxis()->SetRangeUser(Max_bin+j*13,4096);
      if(hist[i]->GetMaximumBin()-(Max_bin+j*13)>4){
        second_Max_bin = hist[i]->GetMaximumBin();
        break;
      }
      if(j==24){second_Max_bin = 0;}
    }
    int third_Max_bin=0;
    if(second_Max_bin!=0){
      for(int j=0;j<25;j++){
        hist[i]->GetXaxis()->SetRangeUser(second_Max_bin+j*13,4096);
        if(hist[i]->GetMaximumBin()-(second_Max_bin+j*13)>5){
          third_Max_bin =hist[i]->GetMaximumBin();
          break;
        }
        if(j==24){third_Max_bin = 0;}
      }
    }
    if(std::abs(second_Max_bin-third_Max_bin)<10){
      std::cout << "channel" << i << "is not good" << std::endl;
      std::cout << "second_Max_bin" << second_Max_bin << std::endl;
      std::cout << "third_Max_bin" << third_Max_bin << std::endl;
    }
    hist[i]->GetXaxis()->SetRangeUser(0,4096);
    TF1 *f1 = new TF1("f1","[0]*exp(-(x-[2])^2/(2*[1]^2))",Max_bin-10,Max_bin+10);
    // f1->SetParLimitss(1000,5,30,770,1000,5,30,830,1000,5,30,902,1000,5,30,1000);
    f1->SetParameter(0,hist[i]->GetBinContent(Max_bin));
    f1->SetParLimits(0,hist[i]->GetBinContent(Max_bin)*0.5,hist[i]->GetBinContent(Max_bin)*1.5);
    f1->SetParameter(1,3);
    f1->SetParLimits(1,3,40);
    f1->SetParameter(2,Max_bin);
    f1->SetParLimits(2,Max_bin-5,Max_bin+5);
    // f1->SetParameter(3,hist[i]->GetBinContent(second_Max_bin));
    // f1->SetParLimits(3,hist[i]->GetBinContent(second_Max_bin)*0.5,hist[i]->GetBinContent(second_Max_bin)*1.5);
    // f1->SetParameter(4,3);
    // f1->SetParLimits(4,3,15);
    // f1->SetParameter(5,second_Max_bin);
    // f1->SetParLimits(5,second_Max_bin-5,second_Max_bin+5);
    // f1->SetParameter(6,hist[i]->GetBinContent(third_Max_bin));
    // f1->SetParLimits(6,hist[i]->GetBinContent(third_Max_bin)*0.8,1000000);
    // f1->SetParameter(7,5);
    // f1->SetParLimits(7,5,25);
    // f1->SetParameter(8,third_Max_bin);
    // f1->SetParLimits(8,third_Max_bin-5,third_Max_bin+5);
    // f1->SetParameter(9,1000);
    // f1->SetParameter(10,5);
    // f1->SetParLimits(10,5,30);
    // f1->SetParameter(11,1000);
    // f1->SetParLimits(11,1000,1140);
    TF1*f2 = new TF1("f2","[0]*exp(-(x-[2])^2/(2*[1]^2))",second_Max_bin-10,second_Max_bin+10);
    f2->SetParameter(0,hist[i]->GetBinContent(second_Max_bin));
    f2->SetParLimits(0,hist[i]->GetBinContent(second_Max_bin)*0.5,hist[i]->GetBinContent(second_Max_bin)*1.5);
    f2->SetParameter(1,3);
    f2->SetParLimits(1,3,9);
    f2->SetParameter(2,second_Max_bin);
    f2->SetParLimits(2,second_Max_bin-5,second_Max_bin+5);
    gStyle->SetOptFit(1111);
    hist[i]->Fit("f1","0M");
    hist[i]->Fit("f2","0M");
    double n_photons[2] ={0,1};
    double n_photons_gosa[2] ={0,0};
    double mean[2] ={f1->GetParameter(2),f2->GetParameter(2)};
    double mean_gosa[2] ={f1->GetParError(2),f2->GetParError(2)};
    // double sigma[3] ={f1->GetParameter(1),f1->GetParameter(4),f1->GetParameter(7)};
    // double sigma_gosa[3] ={f1->GetParError(1),f1->GetParError(4),f1->GetParError(7)};
    f1->SetNpx(10000);
    f2->SetNpx(10000);
    f1->SetLineColor(2);
    f2->SetLineColor(3);
    TCanvas *c1 = new TCanvas("c1","c1",2000,1000);
    c1->Divide(2,1);
    c1->cd(1);
    hist[i]->GetXaxis()->SetRangeUser(700,1300);
    hist[i]->Draw();
    f1->Draw("same");
    f2->Draw("same");
    if(second_Max_bin==0 && third_Max_bin==0){
      std::cout << "channel" << i << "is not good" << std::endl;
      notgood.push_back(i);
      TLatex *tex = new TLatex(0.3,0.3,"not good");
      tex->SetNDC();
      tex->SetTextSize(0.1);
      tex->Draw();
    }
    c1->cd(2);
    TF1*pol = new TF1("pol","pol1",0,2);
    TGraphErrors *g1 = new TGraphErrors(2,n_photons,mean,n_photons_gosa,mean_gosa);
    pol->SetParameters(0,600);
    pol->SetLineWidth(1);
    g1->Fit("pol");
    x[i]=i;
    xe[i]=0;
    bias[i]=pol->GetParameter(0);
    biaserr[i]=pol->GetParError(0);
    gain[i]=pol->GetParameter(1);
    gainerr[i]=pol->GetParError(1);
    sigma0[i]=f1->GetParameter(1);
    sigma0err[i]=f1->GetParError(1);
    sigma1[i]=f2->GetParameter(1);
    sigma1err[i]=f2->GetParError(1);
    if(second_Max_bin==0 && third_Max_bin==0){
      bias[i]=std::numeric_limits<double>::quiet_NaN();
      biaserr[i]=std::numeric_limits<double>::quiet_NaN();
      gain[i]=std::numeric_limits<double>::quiet_NaN();
      gainerr[i]=std::numeric_limits<double>::quiet_NaN();
      sigma0[i]=std::numeric_limits<double>::quiet_NaN();
      sigma0err[i]=std::numeric_limits<double>::quiet_NaN();
      sigma1[i]=std::numeric_limits<double>::quiet_NaN();
      sigma1err[i]=std::numeric_limits<double>::quiet_NaN();
    }
    g1->SetTitle("gain");
    g1->GetXaxis()->SetTitle("n_photons");
    g1->GetYaxis()->SetTitle("mean");
    g1->SetMarkerStyle(20);
    g1->SetMarkerSize(1.5);
    g1->GetXaxis()->SetLimits(-0.5,1.5);
    g1->Draw("AP");
    string histname;
    if(i==0){histname = filename+".pdf(";}
    else if(i==63){histname = filename+".pdf)";}
    else{histname = filename+".pdf";}
    c1->SaveAs(histname.c_str(),"pdf"); 
  }
  TCanvas *c2 = new TCanvas("c2","c2",2000,1000);
  TGraphErrors *biass = new TGraphErrors(64,x,bias,xe,biaserr);
  TGraphErrors *gains = new TGraphErrors(64,x,gain,xe,gainerr);
  TGraphErrors *sigmas = new TGraphErrors(64,x,sigma0,xe,sigma0err);
  TGraphErrors *sigmas1 = new TGraphErrors(64,x,sigma1,xe,sigma1err);
  biass->SetTitle("bias");
  biass->SetName("bias");
  biass->GetXaxis()->SetTitle("channel");
  biass->GetYaxis()->SetTitle("bias");
  biass->SetMarkerStyle(20);
  biass->SetMarkerSize(1.5);
  biass->GetXaxis()->SetLimits(-0.5,63.5);
  gains->SetTitle("gain");
  gains->SetName("gain");
  gains->GetXaxis()->SetTitle("channel");
  gains->GetYaxis()->SetTitle("gain");
  gains->SetMarkerStyle(20);
  gains->SetMarkerSize(1.5);
  gains->GetXaxis()->SetLimits(-0.5,63.5);
  sigmas->SetTitle("sigma0");
  sigmas->SetName("sigma0");
  sigmas->GetXaxis()->SetTitle("channel");
  sigmas->GetYaxis()->SetTitle("sigma0");
  sigmas->SetMarkerStyle(20);
  sigmas->SetMarkerSize(1.5);
  sigmas1->SetTitle("sigma1");
  sigmas1->SetName("sigma1");
  sigmas1->GetXaxis()->SetTitle("channel");
  sigmas1->GetYaxis()->SetTitle("sigma1");
  sigmas1->SetMarkerStyle(20);
  sigmas1->SetMarkerSize(1.5);

  c2->Divide(2,2);
  c2->cd(1);
  for(int oo: notgood){
    TLine *line = new TLine(oo,770,oo,820);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw();
  }
  biass->Draw("APsame");
  biass->Write();
  c2->cd(2);
  for(int oo: notgood){
    TLine *line = new TLine(oo,60,oo,80);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw();
  }
  gains->Draw("APsame");
  gains->Write();
  c2->cd(3);
  for(int oo: notgood){
    TLine *line = new TLine(oo,3,oo,5);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw();
  }
  sigmas->Draw("APsame");
  sigmas->Write();
  c2->cd(4);
  for(int oo: notgood){
    TLine *line = new TLine(oo,3,oo,5);
    line->SetLineColor(2);
    line->SetLineWidth(2);
    line->Draw();
  }
  sigmas1->Draw("APsame");
  sigmas1->Write();
  c2->SaveAs((filename+"_gainbias.png").c_str(),"png");
  file->Write();
  file->Close();
  for(int oo: notgood){
    std::cout << oo << std::endl;
  }
  return;
}
int main(int argc, char* argv[]){
  if(argc == 2){ val(string(argv[1]));}
  else{ return 1;}
  return 0;
}