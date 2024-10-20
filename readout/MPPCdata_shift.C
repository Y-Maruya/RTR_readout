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
#include <iostream>
#include "RooWorkspace.h"
#include "RooFitResult.h"
#include "RooRealVar.h"
#include "RooRealSumPdf.h"
#include "RooAbsCollection.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooCrystalBall.h"
#include "RooAddPdf.h"
#include "RooPlot.h"
#include "RooHist.h"
#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooCFunction1Binding.h"
#include "RooCFunction3Binding.h"
#include "RooCFunction4Binding.h"
#include "RooPlot.h"
#include "RooTFnPdfBinding.h"
#include "TFile.h"
#include "TH1D.h"
#include "TCanvas.h"

using namespace RooFit;
using namespace RooStats;
TF1*f1;
double fmul(double *x, double *p){
    return p[4]*f1->EvalPar(x,p);
}
double erfc(double x,double mean, double sigma,double mean_0){
    if(mean<mean_0){
        return 0;
    }
    if(TMath::Erfc((x-mean)/(sigma))-TMath::Erfc((x-mean_0)/(4.5))<0){
        return 0;
    }
    return TMath::Erfc((x-mean)/(sigma))-TMath::Erfc((x-mean_0)/(4.5));
}
Double_t DSCB(double* x, double* p)
{
  double u   = (x[0]-p[0])/p[1];
  double A1  = TMath::Power(p[3]/TMath::Abs(p[2]),p[3])*TMath::Exp(-p[2]*p[2]/2);
  double B1  = p[3]/TMath::Abs(p[2]) - TMath::Abs(p[2]);

  double result(p[4]);
  if      (u<-p[2]) result *= A1*TMath::Power(B1-u,-p[3]);
  else result *= TMath::Exp(-u*u/2);
  return result;
}
void MPPCdata_DSCB(std::string filename = "laser_HV56p24_0913_1.root",std::string randomfilename = "random_HV56p24_0913_1.root",std::string valfilename = "laser_HV56p24_0913_1.root"){ 
    // ファイルを開く
    double cut = 0;
    TFile *file = TFile::Open(filename.c_str());
    if (!file) {
        std::cerr << "Error: File not found." << std::endl;
        return;
    }
    TFile *randomfile = TFile::Open(randomfilename.c_str());
    if (!randomfile) {
        std::cerr << "Error: random File not found." << std::endl;
        return;
    }
    TFile* valfile = TFile::Open(valfilename.c_str());
    if (!valfile) {
        std::cerr << "Error: val File not found." << std::endl;
        return;
    }
    TGraphErrors *g = (TGraphErrors*)valfile->Get("bias");
    if (!g) {
        std::cerr << "Error: Graph not found." << std::endl;
        return;
    }
    TGraphErrors *g2 = (TGraphErrors*)valfile->Get("gain");
    if (!g2) {
        std::cerr << "Error: Graph not found." << std::endl;
        return;
    }
    std::vector<Int_t> notgood;
    for(int i=0;i<64;i++){
        if(std::isnan(g->Eval(i))){
            std::cout << "channel" << i << "is not good" << std::endl;
            notgood.push_back(i);
        }
    }
    std::string HVstr = filename.substr(filename.find("HV")+2,5);
    TH1F *randomhist[64];
    TH1F *targethist[64];
    TH1F *valhist[64];
    TH1F *randomhist_pre[64];
    TH1F *targethist_pre[64];
    TH1F *valhist_pre[64];
    int start_bin = -100;
    int end_bin = 250;
    int nbin = 350;
    for(int i = 0; i<64; ++i){
        if(std::find(notgood.begin(),notgood.end(),i)!=notgood.end()){
            continue;
        }
        std::ostringstream _name;
        _name << "ch" << i;
        std::string hist_name = _name.str();
        randomhist_pre[i] = (TH1F*)randomfile->Get(hist_name.c_str());
        targethist_pre[i] = (TH1F*)file->Get(hist_name.c_str());
        valhist_pre[i] = (TH1F*)valfile->Get(hist_name.c_str());
        randomhist[i] = new TH1F(Form("random_ch%d",i),Form("random trigger data of ch%d",i),nbin,start_bin,end_bin);
        targethist[i] = new TH1F(Form("target_ch%d",i),Form("target data of ch%d",i),nbin,start_bin,end_bin);
        valhist[i] = new TH1F(Form("val_ch%d",i),Form("val data of ch%d",i),nbin,start_bin,end_bin);
        int target_max = targethist_pre[i]->GetMaximumBin();
        TF1 *fit_target = new TF1("fit_target","gaus",target_max-20,target_max+20);
        targethist_pre[i]->Fit(fit_target,"Q","",target_max-20,target_max+20);
        int target_mean = std::round(fit_target->GetParameter(1));
        
        int random_max = randomhist_pre[i]->GetMaximumBin();
        TF1 *fit_random = new TF1("fit_random","gaus",random_max-20,random_max+20);
        randomhist_pre[i]->Fit(fit_random,"Q","",random_max-20,random_max+20);
        int random_mean = std::round(fit_random->GetParameter(1));
        int val_max = valhist_pre[i]->GetMaximumBin();
        TF1 *fit_val = new TF1("fit_val","gaus",val_max-20,val_max+20);
        valhist_pre[i]->Fit(fit_val,"Q","",val_max-20,val_max+20);
        int val_mean =std::round( fit_val->GetParameter(1));
        for(int j=random_mean+start_bin;j<random_mean+end_bin;j++){
            randomhist[i]->SetBinContent(j-random_mean-start_bin+1,randomhist_pre[i]->GetBinContent(j));
        }
        for(int j=target_mean+start_bin;j<target_mean+end_bin;j++){
            targethist[i]->SetBinContent(j-target_mean-start_bin+1,targethist_pre[i]->GetBinContent(j));
        }
        for(int j=val_mean+start_bin;j<val_mean+end_bin;j++){
            valhist[i]->SetBinContent(j-val_mean-start_bin+1,valhist_pre[i]->GetBinContent(j));
        }
        int Max_bin = valhist[i]->GetMaximumBin();
        int second_Max_bin=0;
        for(int j=0;j<25;j++){
            if(Max_bin+j*13<nbin){
                valhist[i]->GetXaxis()->SetRange(Max_bin+j*13,nbin);
                if(valhist[i]->GetMaximumBin()-(Max_bin+j*13)>3){
                    second_Max_bin = valhist[i]->GetMaximumBin();
                    break;
                }
                if(j==24){second_Max_bin = 0;}
            }
        }
        int third_Max_bin=0;
        if(second_Max_bin!=0){
            for(int j=0;j<25;j++){
                if(second_Max_bin+j*13<nbin){
                    valhist[i]->GetXaxis()->SetRange(second_Max_bin+j*13,nbin);
                    if(valhist[i]->GetMaximumBin()-(second_Max_bin+j*13)>5){
                        third_Max_bin =valhist[i]->GetMaximumBin();
                        break;
                    }
                    if(j==24){third_Max_bin = 0;}
                }
            }
        }
        if(std::abs(second_Max_bin-third_Max_bin)<10){
            std::cout << "channel" << i << "is not good" << std::endl;
            std::cout << "second_Max_bin" << second_Max_bin << std::endl;
            std::cout << "third_Max_bin" << third_Max_bin << std::endl;
        }
        int nbin = randomhist[i]->GetNbinsX();
        if(nbin != targethist[i]->GetNbinsX()){
            std::cerr << "Error: The number of bins in the random histogram and the target histogram are different." << std::endl;
            return;
        }
        targethist[i]->Fit("gaus","Q","",second_Max_bin-20+start_bin,second_Max_bin+20+start_bin);
        TF1 *fit0 = targethist[i]->GetFunction("gaus");
        double n_gaus1_pre = fit0->GetParameter(0)*TMath::Sqrt(2*3.141592)*fit0->GetParameter(2);
        double mean_1_pre = fit0->GetParameter(1);
        double sigma_1_pre = fit0->GetParameter(2);
        targethist[i]->Fit("gaus","Q","",third_Max_bin-20+start_bin,third_Max_bin+20+start_bin);
        TF1 *fit = targethist[i]->GetFunction("gaus");
        double n_gaus2_pre = fit->GetParameter(0)*TMath::Sqrt(2*3.141592)*fit->GetParameter(2);
        double mean_2_pre = fit->GetParameter(1);
        double sigma_2_pre = fit->GetParameter(2);
        RooRealVar ADCvalue(Form("ADCvalue_ch%d",i),Form("ADC value of ch%d",i),start_bin,end_bin);
        RooDataHist random(Form("random_ch%d",i),Form("random trigger data of ch%d",i),RooArgList(ADCvalue),randomhist[i]);
        RooDataHist target(Form("target_ch%d",i),Form("target data of ch%d",i),RooArgList(ADCvalue),targethist[i]);
        RooHistPdf random_pdf(Form("random_pdf_ch%d",i),Form("random PDF of ch%d",i),RooArgSet(ADCvalue),random);
        RooRealVar mean_1(Form("mean_1_ch%d",i),Form("mean of 1 photon gaussian of ch%d",i),mean_1_pre-5,mean_1_pre+5);
        RooRealVar mean_2(Form("mean_2_ch%d",i),Form("mean of 2 photon gaussian of ch%d",i),mean_2_pre-5,mean_2_pre+5);
        // RooRealVar mean_3(Form("mean_3_ch%d",i),Form("mean of 3 photon gaussian of ch%d",i),third_Max_bin,4096);
        RooRealVar sigma_L_1(Form("sigma_L_1_ch%d",i),Form("sigma L of 1 photon gaussian ch%d",i),sigma_1_pre-0.1,sigma_1_pre+0.1);
        RooRealVar sigma_L_2(Form("sigma_L_2_ch%d",i),Form("sigma L of 2 photon gaussian ch%d",i),sigma_2_pre-0.2,sigma_2_pre+0.2);
        // RooRealVar sigma_3(Form("sigma_3_ch%d",i),Form("sigma of 3 photon gaussian ch%d",i),0,15);
        RooRealVar sigma_R_1(Form("sigma_R_1_ch%d",i),Form("sigma R of 1 photon gaussian ch%d",i),sigma_1_pre-0.1,sigma_1_pre+0.1);
        RooRealVar sigma_R_2(Form("sigma_R_2_ch%d",i),Form("sigma R of 2 photon gaussian ch%d",i),sigma_2_pre-0.2,sigma_2_pre+0.2);
        // RooRealVar sigma_3(Form("sigma_3_ch%d",i),Form("sigma of 3 photon gaussian ch%d",i),0,15);
        RooRealVar alpha_L_1(Form("alpha_L_1_ch%d",i),Form("alpha L of 1 photon gaussian ch%d",i),0.1,10);
        RooRealVar alpha_L_2(Form("alpha_L_2_ch%d",i),Form("alpha L of 2 photon gaussian ch%d",i),0.1,10);
        // RooRealVar alpha_L_3(Form("alpha_L_3_ch%d",i),Form("alpha L of 3 photon gaussian ch%d",i),1,10);
        RooRealVar alpha_R_1(Form("alpha_R_1_ch%d",i),Form("alpha R of 1 photon gaussian ch%d",i),1,10);
        RooRealVar alpha_R_2(Form("alpha_R_2_ch%d",i),Form("alpha R of 2 photon gaussian ch%d",i),1,10);
        // RooRealVar alpha_R_3(Form("alpha_R_3_ch%d",i),Form("alpha R of 3 photon gaussian ch%d",i),1,10);
        RooRealVar n_L_1(Form("n_L_1_ch%d",i),Form("n L of 1 photon gaussian ch%d",i),0,10);
        RooRealVar n_L_2(Form("n_L_2_ch%d",i),Form("n L of 2 photon gaussian ch%d",i),0,10);
        // RooRealVar n_L_3(Form("n_L_3_ch%d",i),Form("n L of 3 photon gaussian ch%d",i),0,1000000);
        RooRealVar n_R_1(Form("n_R_1_ch%d",i),Form("n R of 1 photon gaussian ch%d",i),0,10);
        RooRealVar n_R_2(Form("n_R_2_ch%d",i),Form("n R of 2 photon gaussian ch%d",i),0,10);
        // RooRealVar n_R_3(Form("n_R_3_ch%d",i),Form("n R of 3 photon gaussian ch%d",i),0,1000000);
        RooCrystalBall CB1_pdf(Form("CB1_pdf_ch%d",i),Form("1 photon Crystal Ball PDF of ch%d",i),ADCvalue,mean_1,sigma_L_1,alpha_L_1,n_L_1);
        RooCrystalBall CB2_pdf(Form("CB2_pdf_ch%d",i),Form("2 photon Crystal Ball PDF of ch%d",i),ADCvalue,mean_2,sigma_L_2,alpha_L_2,n_L_2);
        RooRealVar n_random(Form("n_random_ch%d",i),Form("number of event at pedestal PDF of ch%d",i),0,targethist[i]->Integral());
        RooRealVar n_gaus1(Form("n_gaus1_ch%d",i),Form("number of event at 1 photon gaussian PDF of ch%d",i),n_gaus1_pre*2/3,n_gaus1_pre*3/2);
        RooRealVar n_gaus2(Form("n_gaus2_ch%d",i),Form("number of event at 2 photon gaussian PDF of ch%d",i),n_gaus2_pre*2/3,n_gaus2_pre*3/2);
        RooAddPdf model(Form("model_ch%d",i),Form("model of ch%d",i),RooArgList(random_pdf,CB1_pdf,CB2_pdf),RooArgList(n_random,n_gaus1,n_gaus2));
        RooFitResult* result = model.fitTo(target, RooFit::Save(),Range(start_bin,end_bin));
        // result->Print();
        std::cout<<"====================="<<std::endl;
        std::cout<<"channel"<<i<<std::endl;
        std::cout<<mean_1.getVal()<<std::endl;
        std::cout<<mean_2.getVal()<<std::endl;
        // std::cout<<mean_3.getVal()<<std::endl;
        std::cout<<sigma_L_1.getVal()<<std::endl;
        std::cout<<sigma_L_2.getVal()<<std::endl;
        // std::cout<<sigma_R_1.getVal()<<std::endl;
        // std::cout<<sigma_R_2.getVal()<<std::endl;
        std::cout<<alpha_L_1.getVal()<<std::endl;
        std::cout<<alpha_L_2.getVal()<<std::endl;
        // std::cout<<alpha_R_1.getVal()<<std::endl;
        // std::cout<<alpha_R_2.getVal()<<std::endl;
        std::cout<<n_L_1.getVal()<<std::endl;
        std::cout<<n_L_2.getVal()<<std::endl;
        // std::cout<<n_R_1.getVal()<<std::endl;
        // std::cout<<n_R_2.getVal()<<std::endl;
        // std::cout<<sigma_3.getVal()<<std::endl;
        std::cout<<n_random.getVal()<<std::endl;
        std::cout<<n_gaus1.getVal()<<std::endl;
        std::cout<<n_gaus2.getVal()<<std::endl;
        // std::cout<<n_erf1.getVal()<<std::endl;
        // std::cout<<n_erf2.getVal()<<std::endl;
        // std::cout<<n_gaus3.getVal()<<std::endl;
        std::cout<<"====================="<<std::endl;
        double post_mean_1 = mean_1.getVal();
        double post_mean_2 = mean_2.getVal();
        double post_sigma_L_1 = sigma_L_1.getVal();
        double post_sigma_L_2 = sigma_L_2.getVal();
        double post_alpha_L_1 = alpha_L_1.getVal();
        double post_alpha_L_2 = alpha_L_2.getVal();
        double post_n_L_1 = n_L_1.getVal();
        double post_n_L_2 = n_L_2.getVal();
        double post_n_random = n_random.getVal();
        double post_n_gaus1 = n_gaus1.getVal();
        double post_n_gaus2 = n_gaus2.getVal();
        TCanvas *c = new TCanvas(Form("c_ch%d",i),Form("c_ch%d",i),800,1400);
        c->Divide(1,3);
        c->cd(1);
        RooPlot* frame = ADCvalue.frame(Title(Form("ch%d",i)),Range(start_bin,end_bin),Bins(nbin));
        target.plotOn(frame,Name("target"));
        model.plotOn(frame,LineColor(kMagenta),Name("model"));
        model.plotOn(frame,Components(random_pdf),LineStyle(kDashed),LineColor(kRed),Name("random"));
        model.plotOn(frame,Components(CB1_pdf),LineStyle(kDashed),LineColor(kGreen),Name("CB1"));
        model.plotOn(frame,Components(CB2_pdf),LineStyle(kDashed),LineColor(kBlue),Name("CB2"));
        frame->Draw();
        c->cd(2)->SetLogy();
        frame->Draw();
        c->cd(3);
        RooHist * h_resid = frame->residHist("target","model");
        RooPlot* frame2 = ADCvalue.frame(Title(Form("ch%d residual",i)),Range(start_bin,end_bin),Bins(nbin));
        frame2->addPlotable(h_resid,"L");
        frame2->Draw();
        if(i==0){
            c->Print("fit_result_test/fit_CB_shift.pdf(","pdf");
        }else if(i==63){
            c->Print("fit_result_test/fit_CB_shift.pdf)","pdf");
        }else{
            c->Print("fit_result_test/fit_CB_shift.pdf","pdf");
        }
    }
    file->Close();
    randomfile->Close();
    valfile->Close();
    return;
}

void MPPCdata(std::string filename = "laser_HV56p24_0913_1.root",std::string randomfilename = "random_HV56p24_0913_1.root",std::string valfilename = "laser_HV56p24_0913_1.root"){ 
    // ファイルを開く
    double cut = 0;
    TFile *file = TFile::Open(filename.c_str());
    if (!file) {
        std::cerr << "Error: File not found." << std::endl;
        return;
    }
    TFile *randomfile = TFile::Open(randomfilename.c_str());
    if (!randomfile) {
        std::cerr << "Error: random File not found." << std::endl;
        return;
    }
    TFile* valfile = TFile::Open(valfilename.c_str());
    if (!valfile) {
        std::cerr << "Error: val File not found." << std::endl;
        return;
    }
    TGraphErrors *g = (TGraphErrors*)valfile->Get("bias");
    if (!g) {
        std::cerr << "Error: Graph not found." << std::endl;
        return;
    }
    TGraphErrors *g2 = (TGraphErrors*)valfile->Get("gain");
    if (!g2) {
        std::cerr << "Error: Graph not found." << std::endl;
        return;
    }
    std::vector<Int_t> notgood;
    for(int i=0;i<64;i++){
        if(std::isnan(g->Eval(i))){
            std::cout << "channel" << i << "is not good" << std::endl;
            notgood.push_back(i);
        }
    }
    std::string HVstr = filename.substr(filename.find("HV")+2,5);
    TH1F *randomhist_pre[64];
    TH1F *targethist_pre[64];
    TH1F *valhist_pre[64];
    TH1F *randomhist[64];
    TH1F *targethist[64];
    TH1F *valhist[64];
    int start_bin = -100;
    int end_bin = 250;
    int nbin = 350;
    for(int i = 0; i<64; ++i){
        if(std::find(notgood.begin(),notgood.end(),i)!=notgood.end()){
            continue;
        }
        std::ostringstream _name;
        _name << "ch" << i;
        std::string hist_name = _name.str();
        randomhist_pre[i] = (TH1F*)randomfile->Get(hist_name.c_str());
        targethist_pre[i] = (TH1F*)file->Get(hist_name.c_str());
        valhist_pre[i] = (TH1F*)valfile->Get(hist_name.c_str());
        randomhist[i] = new TH1F(Form("random_ch%d",i),Form("random trigger data of ch%d",i),nbin,start_bin,end_bin);
        targethist[i] = new TH1F(Form("target_ch%d",i),Form("target data of ch%d",i),nbin,start_bin,end_bin);
        valhist[i] = new TH1F(Form("val_ch%d",i),Form("val data of ch%d",i),nbin,start_bin,end_bin);
        int target_max = targethist_pre[i]->GetMaximumBin();
        int random_max = randomhist_pre[i]->GetMaximumBin();
        int val_max = valhist_pre[i]->GetMaximumBin();
        for(int j=random_max+start_bin;j<random_max+end_bin;j++){
            randomhist[i]->SetBinContent(j-random_max-start_bin+1,randomhist_pre[i]->GetBinContent(j));
        }
        for(int j=target_max+start_bin;j<target_max+end_bin;j++){
            targethist[i]->SetBinContent(j-target_max-start_bin+1,targethist_pre[i]->GetBinContent(j));
        }
        for(int j=val_max+start_bin;j<val_max+end_bin;j++){
            valhist[i]->SetBinContent(j-val_max-start_bin+1,valhist_pre[i]->GetBinContent(j));
        }

        int Max_bin = valhist[i]->GetMaximumBin();
        int second_Max_bin=0;
        for(int j=0;j<25;j++){
            if(Max_bin+j*13<nbin){
                valhist[i]->GetXaxis()->SetRange(Max_bin+j*13,nbin);
                if(valhist[i]->GetMaximumBin()-(Max_bin+j*13)>3){
                    second_Max_bin = valhist[i]->GetMaximumBin();
                    break;
                }
                if(j==24){second_Max_bin = 0;}
            }
        }
        int third_Max_bin=0;
        if(second_Max_bin!=0){
            for(int j=0;j<25;j++){
                if(second_Max_bin+j*13<nbin){
                    valhist[i]->GetXaxis()->SetRange(second_Max_bin+j*13,nbin);
                    if(valhist[i]->GetMaximumBin()-(second_Max_bin+j*13)>5){
                        third_Max_bin =valhist[i]->GetMaximumBin();
                        break;
                    }
                    if(j==24){third_Max_bin = 0;}
                }
            }

        }
        std::cout << "channel" << i << std::endl;
        std::cout << "Max_bin" << Max_bin << std::endl;
        std::cout << "second_Max_bin" << second_Max_bin << std::endl;
        std::cout << "third_Max_bin" << third_Max_bin << std::endl;
        if(std::abs(second_Max_bin-third_Max_bin)<10){
            std::cout << "channel" << i << "is not good" << std::endl;
            std::cout << "second_Max_bin" << second_Max_bin << std::endl;
            std::cout << "third_Max_bin" << third_Max_bin << std::endl;
        }
        int nbin = randomhist[i]->GetNbinsX();
        if(nbin != targethist[i]->GetNbinsX()){
            std::cerr << "Error: The number of bins in the random histogram and the target histogram are different." << std::endl;
            return;
        }
        targethist[i]->Fit("gaus","Q","",second_Max_bin-20+start_bin,second_Max_bin+20+start_bin);
        TF1 *fit0 = targethist[i]->GetFunction("gaus");
        double n_gaus1_pre = fit0->GetParameter(0)*TMath::Sqrt(2*3.141592)*fit0->GetParameter(2);
        double mean_1_pre = fit0->GetParameter(1);
        std::cout << "mean_1_pre : "<< mean_1_pre << std::endl;
        double sigma_1_pre = fit0->GetParameter(2);
        targethist[i]->Fit("gaus","Q","",third_Max_bin-20+start_bin,third_Max_bin+20+start_bin);
        TF1 *fit = targethist[i]->GetFunction("gaus");
        double n_gaus2_pre = fit->GetParameter(0)*TMath::Sqrt(2*3.141592)*fit->GetParameter(2);
        double mean_2_pre = fit->GetParameter(1);
        std::cout << "mean_2_pre : "<< mean_2_pre << std::endl;
        double sigma_2_pre = fit->GetParameter(2);
        RooRealVar ADCvalue(Form("ADCvalue_ch%d",i),Form("ADC value of ch%d",i),start_bin,end_bin);
        RooDataHist random(Form("random_ch%d",i),Form("random trigger data of ch%d",i),RooArgList(ADCvalue),randomhist[i]);
        RooDataHist target(Form("target_ch%d",i),Form("target data of ch%d",i),RooArgList(ADCvalue),targethist[i]);
        RooHistPdf random_pdf(Form("random_pdf_ch%d",i),Form("random PDF of ch%d",i),RooArgSet(ADCvalue),random);
        RooRealVar mean_1(Form("mean_1_ch%d",i),Form("mean of 1 photon gaussian of ch%d",i),mean_1_pre-5,mean_1_pre+5);
        RooRealVar mean_2(Form("mean_2_ch%d",i),Form("mean of 2 photon gaussian of ch%d",i),mean_2_pre-5,mean_2_pre+5);
        // RooRealVar mean_3(Form("mean_3_ch%d",i),Form("mean of 3 photon gaussian of ch%d",i),third_Max_bin,4096);
        RooRealVar sigma_1(Form("sigma_1_ch%d",i),Form("sigma of 1 photon gaussian ch%d",i),sigma_1_pre-0.1,sigma_1_pre+0.1);
        RooRealVar sigma_2(Form("sigma_2_ch%d",i),Form("sigma of 2 photon gaussian ch%d",i),sigma_2_pre-0.1,sigma_2_pre+0.1);
        // RooRealVar sigma_3(Form("sigma_3_ch%d",i),Form("sigma of 3 photon gaussian ch%d",i),0,15);
        RooGaussian gaus1_pdf(Form("gaus1_pdf_ch%d",i),Form("1 photon gaussian PDF of ch%d",i),ADCvalue,mean_1,sigma_1);
        RooGaussian gaus2_pdf(Form("gaus2_pdf_ch%d",i),Form("2 photon gaussian PDF of ch%d",i),ADCvalue,mean_2,sigma_2);
        // RooGaussian gaus3_pdf(Form("gaus3_pdf_ch%d",i),Form("3 photon gaussian PDF of ch%d",i),ADCvalue,mean_3,sigma_3);
        
        RooRealVar n_random(Form("n_random_ch%d",i),Form("number of event at pedestal PDF of ch%d",i),0,targethist[i]->Integral());
        RooRealVar n_gaus1(Form("n_gaus1_ch%d",i),Form("number of event at 1 photon gaussian PDF of ch%d",i),n_gaus1_pre*2/3,n_gaus1_pre*3/2);
        RooRealVar n_gaus2(Form("n_gaus2_ch%d",i),Form("number of event at 2 photon gaussian PDF of ch%d",i),n_gaus2_pre*2/3,n_gaus2_pre*3/2);
        // RooRealVar n_gaus3(Form("n_gaus3_ch%d",i),Form("number of event at 3 photon gaussian PDF of ch%d",i),0,1000000);
        RooAddPdf model(Form("model_ch%d",i),Form("model of ch%d",i),RooArgList(random_pdf,gaus1_pdf,gaus2_pdf),RooArgList(n_random,n_gaus1,n_gaus2));
        RooFitResult* result = model.fitTo(target, RooFit::Save(),Range(start_bin,end_bin));
        // result->Print();
        std::cout<<"====================="<<std::endl;
        std::cout<<"channel"<<i<<std::endl;
        std::cout<<mean_1.getVal()<<std::endl;
        std::cout<<mean_2.getVal()<<std::endl;
        // std::cout<<mean_3.getVal()<<std::endl;
        std::cout<<sigma_1.getVal()<<std::endl;
        std::cout<<sigma_2.getVal()<<std::endl;
        // std::cout<<sigma_3.getVal()<<std::endl;
        std::cout<<n_random.getVal()<<std::endl;
        std::cout<<n_gaus1.getVal()<<std::endl;
        std::cout<<n_gaus2.getVal()<<std::endl;
        // std::cout<<n_erf1.getVal()<<std::endl;
        // std::cout<<n_erf2.getVal()<<std::endl;
        // std::cout<<n_gaus3.getVal()<<std::endl;
        std::cout<<"====================="<<std::endl;
        RooWorkspace w(Form("w_ch%d",i));
        w.import(target);
        w.import(model);
        TCanvas *c = new TCanvas(Form("c_ch%d",i),Form("c_ch%d",i),800,1400);
        c->Divide(1,3);
        c->cd(1);
        RooPlot* frame = ADCvalue.frame(Title(Form("ch%d",i)),Range(start_bin,end_bin),Bins(nbin));
        target.plotOn(frame,Name("target"));
        model.plotOn(frame,LineColor(kMagenta),Name("model"));
        model.plotOn(frame,Components(random_pdf),LineStyle(kDashed),LineColor(kRed),Name("random"));
        model.plotOn(frame,Components(gaus1_pdf),LineStyle(kDashed),LineColor(kGreen),Name("gaus1"));
        model.plotOn(frame,Components(gaus2_pdf),LineStyle(kDashed),LineColor(kBlue),Name("gaus2"));
        frame->Draw();
        c->cd(2)->SetLogy();
        frame->Draw();
        c->cd(3);
        RooHist * h_resid = frame->residHist("target","model");
        RooPlot* frame2 = ADCvalue.frame(Title(Form("ch%d residual",i)),Range(start_bin,end_bin),Bins(nbin));
        frame2->addPlotable(h_resid,"L");
        frame2->Draw();
        if(i==0){
            c->Print("fit_result_test/fit_shift.pdf(","pdf");
        }else if(i==63){
            c->Print("fit_result_test/fit_shift.pdf)","pdf");
        }else{
            c->Print("fit_result_test/fit_shift.pdf","pdf");
        }
        // TCanvas *c2 = new TCanvas(Form("c2_ch%d",i),Form("c2_ch%d",i),800,1000);
        // c2->Divide(1,3);
        // c2->cd(1)->SetLogy();
        // TF1 *fit1 = new TF1("fit1","gaus",700,1000);
        // fit1->SetParameters(n_gaus1.getVal()/TMath::Sqrt(2*3.141592)/sigma_1.getVal(),mean_1.getVal(),sigma_1.getVal());
        // fit1->SetLineColor(kGreen);
        // TF1 *fit2 = new TF1("fit2","gaus",700,1000);
        // fit2->SetParameters(n_gaus2.getVal()/TMath::Sqrt(2*3.141592)/sigma_2.getVal(),mean_2.getVal(),sigma_2.getVal());
        // fit2->SetLineColor(kBlue);
        // randomhist[i]->Scale(n_random.getVal()/randomhist[i]->Integral());
        // // randomhist[i]->Scale(targethist[i]->GetMaximum()/randomhist[i]->GetMaximum());
        // randomhist[i]->SetLineColor(kRed);
        // targethist[i]->SetLineColor(kBlack);
        // targethist[i]->Draw("hist");
        // randomhist[i]->Draw("hist same");
        // // fit1->Draw("same");
        // // fit2->Draw("same");
        // c2->cd(2)->SetLogy();
        // TH1F *residual = (TH1F*)targethist[i]->Clone("residual");
        // residual->Add(randomhist[i],-1);
        // // residual->Add(fit1,-1);
        // // residual->Add(fit2,-1);
        // residual->GetYaxis()->SetTitle("residual (data-pedestal)");
        // residual->Draw("hist");
        // c2->cd(3);
        // residual->Draw("hist");

        // c2->Print(Form("fit_result_test/residual_ch%d.png",i));

    }
    file->Close();
    randomfile->Close();
    valfile->Close();
    return;
}

int main() {
    MPPCdata_DSCB();
    MPPCdata();
    return 0;
}
