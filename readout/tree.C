#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<vector>
#include<TText.h>
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TCanvas.h"
#include"TGraph.h"
using namespace std;


const int MaxEventsLength = 5000; // 時間の長さの最大値
const int MaxEventsNum = 50000;
void tree(int flag,Int_t Triggerstart_value,Int_t Triggerend_value, std::string file_dir){
    TFile*fout=new TFile(Form("%stree_%d.root",file_dir.c_str(),flag),"recreate");
    string file_name = file_dir + Form("wave%d.txt",flag);
    if (fout == nullptr || !fout->IsOpen()) {
        cout << "Failed to open the output file." << endl;
        return;
    }
    TTree *tree =new TTree("tree","test_event_w_fit_data");
    int NUMBER=MaxEventsNum;
    /*ifstream ifs;
    ifs.open("souchi.txt");
    ifs>>file_name;
    ifs>>NUMBER;
    ifs.close();*/
    vector<int> adc_vector;
    vector<int> time_vector;
    //double discri[MaxEventsLength];
    double mean,hensa=0,Qs=0,Ql=0,Qsn=0,Qln=0,psd=-1,meanvol,MAX,MAXvol;
    int trigger=0;
    int trigger_times=0;
    int th =20;
    int eNum;
    ULong64_t TimeStump;
    vector<TH1D*> base_vector;
    //TH1D *basevol[MaxEventsNum];
    //TH1D *hist[MaxEventsNum];
    vector<TH1D*> hist_vector;
    vector<TCanvas*> c;
    string buf,str_data,str_conma_data,hanbetu_buf,ch;
    Int_t h,Leng;
    string file_d("/mnt/hkpmt-kraid1/data/maruya/");
    string slash("/");
    string ubar("_");
    string event("event.pdf");
    string ev =event ;
    string evi =ev+"(" ;
    string evf =ev + ")";
    tree->Branch(Form("Leng_%d",flag),&Leng,"Leng/I");
    tree->Branch(Form("TimeStump_%d",flag),&TimeStump,"TimeStump/I");
    tree->Branch(Form("mean_%d",flag),&mean,"mean/D");
    tree->Branch(Form("hensa_%d",flag),&hensa,"hensa/D");
    tree->Branch(Form("trigger_%d",flag),&trigger,"trigger/I");
    tree->Branch(Form("trigger_times_%d",flag),&trigger_times,"trigger_times/I");
    tree->Branch(Form("Qs_%d",flag),&Qs,"Qs/D");
    tree->Branch(Form("Ql_%d",flag),&Ql,"Ql/D");
    //tree->Branch("psd",&psd,"psd/D");
    tree->Branch(Form("MAX_%d",flag),&MAX,"MAX/D");
    //tree->Branch("MAXvol",&MAX,"MAXvol/D");
    //tree->Branch("maxadc",maxadc,"maxadc[Leng]/I");
    //tree->Branch("minadc",minadc,"minadc[Leng]/I");
    //tree->Branch("maxvol",maxvol,"maxvol[Leng]/I");
    //tree->Branch("minvol",minvol,"minvol[Leng]/I");
    tree->Branch(Form("adc_%d",flag),&adc_vector);
    //tree->Branch("minadcd",minadc,"minadcd[Leng]/I");
    //tree->Branch("maxvold",maxvol,"maxvold[Leng]/I");
    //tree->Branch("minvold",minvol,"minvold[Leng]/I");
    tree->Branch(Form("time_%d",flag),&time_vector);
    //tree->Branch("discri",discri,"discri[Leng]/D");
    //TH1D *q_pmt_dist_70000 =new TH1D("Ql_pmt_dist_70000",Form("Ql_pmt_psd_%s;Ql;psd",file_name.c_str()),100,0,70000);
    //TH1D *q_sci_dist_70000 =new TH1D("Ql_sci_dist_70000",Form("Ql_sci_psd_%s;Ql;psd",file_name.c_str()),100,0,70000);
    //TH1D *q_pmt_dist_20000 =new TH1D("Ql_pmt_dist_20000",Form("Ql_pmt_psd_%s;Ql;psd",file_name.c_str()),100,0,20000);
    //TH1D *q_sci_dist_20000 =new TH1D("Ql_sci_dist_20000",Form("Ql_sci_psd_%s;Ql;psd",file_name.c_str()),100,0,20000);
    TH1D *MAX_dist =new TH1D("MAX_dist",Form("MAX_%s;MAX;count",file_name.c_str()),2000,0,1000);
    TH1D *Ql_dist =new TH1D("Ql_dist",Form("Ql_%s;Ql;count",file_name.c_str()),2000,0,20000);
    //TH1D *MAXvol_dist =new TH1D("MAXvol_dist",Form("MAXvol_%s;MAXvol;count",file_name.c_str()),100,0,1000);
    TH1D *stuck1 =new TH1D("stuck1",Form("stuck_%s;time[ns];sum(mean-wavehigh)",file_name.c_str()),100000,0,400000);
    //TH2D *ql_psd_dist =new TH2D("Ql_psd_dist",Form("Ql_psd_%s;Ql;psd",file_name.c_str()),200,0,20000,50,0,1);
    //TH2D *ql_psd_dist_70000 =new TH2D("Ql_psd_dist_70000",Form("Ql_psd_%s;Ql;psd",file_name.c_str()),700,0,70000,50,0,1);
    TH2D *ql_MAX_dist_20000 =new TH2D("Ql_MAX_dist_20000",Form("Ql_MAX_%s;Ql;MAX",file_name.c_str()),200,0,20000,100,0,1000);
    TH2D *ql_MAX_dist_70000=new TH2D("Ql_MAX_dist_70000",Form("Ql_MAX_%s;Ql;MAX",file_name.c_str()),700,0,70000,100,0,1000);
    ifstream ifs1;
    string file_path = file_dir + Form("wave%d.txt",flag);
        
    ifs1.open(file_path);
    if (!ifs1) {
        cout << "Failed to open file: " << file_path << endl;
        return;
    }    
    ULong64_t before_TimeStump=0;
    int Event_id=0;
    for( Int_t i = 0; i <NUMBER; i++){
        if(i%1000==0){
            cout<<i<<endl;
        }
        TCanvas*c1 = new TCanvas();
        c.push_back(c1);
        c[i]->Divide(1,2);
        c[i]->cd(1);
        //string num(Form("%04d.txt",i));
        if(ifs1.eof()){
            break;
        }
        ifs1>>buf>>buf>>Leng;
        ifs1>>buf>>buf;
        ifs1>>buf>>ch;
        ifs1>>buf>>buf>>eNum;
        getline(ifs1, buf);
        ifs1>>buf>>buf>>buf>>buf>>buf>>TimeStump;
        // cout<<TimeStump<<endl;
        Long64_t TimeStump_res = TimeStump-before_TimeStump;
        // cout<<TimeStump_res*8e-9<<endl;
        before_TimeStump=TimeStump;
        if(abs(TimeStump_res*8e-9-14.1)<1.0){
            Event_id++;
        }
        getline(ifs1, buf);
        getline(ifs1, buf);
        Int_t j=0;
        Int_t hhhh=0;
        TH1D*base= new TH1D(Form("base%d",i+1),";ADC channel;count",2*th+1,2048-th,2048+th);
        base_vector.push_back(base); 
        //basevol[i] = new TH1D(Form("basevol%d",i+1),Form("event_%d;channel;count",i+1),100,-50,50);
        TH1D*hist=new TH1D(Form("Event_adc%d",i+1),";time[ns];count",MaxEventsLength,0,MaxEventsLength*4);
        hist_vector.push_back(hist);
        if(ifs1.eof()){
            break;
        }
        adc_vector.clear();
        time_vector.clear();
        while(j<Leng){
            int maxadcd;
            getline(ifs1,str_data);
            stringstream ss(str_data);
            ss>>maxadcd;
            adc_vector.push_back(maxadcd);
            //cout<<maxadcd[j]<<endl;
            time_vector.push_back(j*4);
            if(maxadcd<2048+th &&maxadcd>2048-th && j < 1500){
                base_vector[i]->Fill(maxadcd);
            }
            j++;
        }
        if(hhhh==1){
            continue;
        }else{
            //G[i] = new TF1(Form("G%d",i+1),"[0]*exp(-0.5*((x-[1])/[2])^2)",-1000,1000);
            //G[i]->SetParameter(0,6000.);
            //G[i]->SetParLimits(0,6000.,20000.);
            //G[i]->SetParameter(1,130.);
            //G[i]->SetParLimits(1,130.,350);
            //G[i]->SetParameter(2,160.);
            //G[i]->SetParLimits(2,160.,2000.);
            //base[i]->Fit(Form("G%d",i+1),"M","",-450,750);
            mean= base_vector[i]->GetMean();
            hensa=base_vector[i]->GetRMS();
            Double_t x[MaxEventsLength],y[MaxEventsLength];
            int a=0;
            trigger = -1;
            for(Int_t k =0; k<Leng;k++){
                if(adc_vector[k]<mean-6*hensa && a==0 ){
                    a=1;
                    trigger=k;
                }
            }
            trigger_times=0;
            for(Int_t k =Triggerstart_value; k<Triggerend_value;k++){
                if(adc_vector[k]<mean-20*hensa ){
                    trigger_times++;
                    if(trigger_times==1){trigger=k;}
                }
            }
            for(Int_t k =0; k<Leng;k++){
                hist_vector[i]->Fill(time_vector[k],mean-adc_vector[k]);
                if(trigger>49800){
                    stuck1->Fill(time_vector[k],mean-adc_vector[k]);
                }
                //histvol[i]->Fill(time[k],meanvol-(maxvol[k]+minvol[k])/2);
            }
            Qs=0;
            Ql=0;
            Qsn=0;
            Qln=0;
            Int_t l=trigger,o=trigger,m =0,ln=trigger,on=trigger;
            while( time_vector[l]- trigger*4<9){
                Qs += (time_vector[l+1]-time_vector[l])*(mean-(adc_vector[l]));
                l++;
            }
            while( time_vector[o]- trigger*4<100){ 
                Ql += (time_vector[o+1]-time_vector[o])*(mean-adc_vector[o]);
                o++;
            }
            /*
            while( time[ln]- trigger*4<20){
                if((maxadcd[ln]+minadcd[ln])/2<-18){
                    Qsn += (time[ln+1]-time[ln])*(mean-(maxadcd[ln]+minadcd[ln])/2);
                }
                ln++;
            }
            while( time[on]- trigger*4<400){ 
                if((maxadcd[on]+minadcd[on])/2<-18){
                    Qln += (time[on+1]-time[on])*(meanadcd-(maxadcd[on]+minadcd[on])/2);
                }
                on++;
            }
            if(Qln!=0){
                psd =1-Qsn/Qln;
            }else{
                psd=0;
            }*/
            /*
            if(trigger>480 && trigger<490){
                while( time[m]- trigger_time + 100<2000){ 
                    stuck1->Fill((time[m]-trigger_time +100),mean-(maxadcd[m]+minadcd[m])/2);
                    m++;
                }
            }*/
            MAX=hist_vector[i]->GetMaximum();
            hist_vector[i]->Write();
            //MAXvol=histvol[i]->GetMaximum();
            //ql_psd_dist->Fill(Ql,psd);
            //ql_psd_dist_70000->Fill(Ql,psd);
            MAX_dist->Fill(MAX);
            Ql_dist->Fill(Ql);
            //MAXvol_dist->Fill(MAXvol);
            ql_MAX_dist_20000->Fill(Ql,MAX);
            ql_MAX_dist_70000->Fill(Ql,MAX);
            tree->Fill();
            hist_vector[i]->GetXaxis()->SetTitle("time[ns]");
            hist_vector[i]->GetYaxis()->SetTitle("base-(waveheigh)[adc channel]");
            gStyle->SetOptStat(0);
            hist_vector[i]->Draw("hist");
            /*TText *t = new TText(0.6,0.4,Form("Ql=%g,trigger=%d",Ql,trigger));
            t->SetTextAlign(22);
            t->SetTextColor(kRed);
            t->SetTextFont(43);
            t->SetTextSize(15);
            t->SetNDC(1);
            t->Draw("SAME");*/
            c[i]->cd(2);
            base_vector[i]->Draw("hist");
            // if(trigger_times>0){
            //     c[i]->Print(Form("%d_event_%d.pdf",flag,i));
            //     TCanvas*c4 = new TCanvas();
            //     hist_vector[i]->GetXaxis()->SetTitle("time[ns]");
            //     hist_vector[i]->GetYaxis()->SetTitle("base-(waveheigh)[adc channel]");
            //     gStyle->SetOptStat(0);
            //     hist_vector[i]->GetXaxis()->SetRangeUser(Triggerstart_value*4,Triggerend_value*4);
            //     hist_vector[i]->Draw("hist");
            //     c4->Print(Form("%d_event_%d_%d_zoom_%g.pdf",flag,i,Event_id,Ql));
            // }
            /*
            if(i==0){
                c[i]->Print(Form("%d_%s",evi.c_str()));
            }else if(ifs1.eof() || i==1548){
                c[i]->Print(Form("%d_%s",evf.c_str()));
            }else{
                c[i]->Print(Form("%d_%s",ev.c_str()));
            }*/
        }
    }
    ifs1.close();
    TCanvas*stu = new TCanvas();
    gPad->SetLogy();
    stuck1->Draw("hist");
    stu->SaveAs(Form("%s%d_stuck.png",file_dir.c_str(),flag));
    /*TCanvas*qp = new TCanvas();
    q_pmt_dist_70000->Draw();
    qp->SaveAs("ql_pmt_dist_70000.png");
    TCanvas*qs = new TCanvas();
    q_sci_dist_70000->Draw();
    qs->SaveAs("ql_sci_dist_70000.png");
    TCanvas*qpq = new TCanvas();
    q_pmt_dist_20000->Draw();
    qpq->SaveAs("ql_pmt_dist_20000.png");
    TCanvas*qsq = new TCanvas();
    q_sci_dist_20000->Draw();
    qsq->SaveAs("ql_sci_dist_20000.png");*/
    TCanvas*max = new TCanvas();
    MAX_dist->Draw();
    max->SaveAs(Form("%s%d_MAX_dist.png",file_dir.c_str(),flag));
    TCanvas*qax = new TCanvas();
    Ql_dist->Draw();
    qax->SaveAs(Form("%s%d_Ql_dist.png",file_dir.c_str(),flag));
    /*TCanvas*vs = new TCanvas();
    MAXvol_dist->Draw();
    vs->SaveAs("MAXvol_dist.png");
    TCanvas*ql = new TCanvas();
    ql_psd_dist->Draw("colz");
    ql->SaveAs("ql_psd_2ddist_20000.png");
    TCanvas*ql_70000 = new TCanvas();
    ql_psd_dist_70000->Draw("colz");
    ql_70000->SaveAs("ql_psd_2ddist_70000.png");*/
    TCanvas*ql_20000 = new TCanvas();
    ql_MAX_dist_20000->Draw("colz");
    ql_20000->SaveAs(Form("%s%d_ql_MAXvol_2ddist_20000.png",file_dir.c_str(),flag));
    TCanvas*q111_20000 = new TCanvas();
    ql_MAX_dist_70000->Draw("colz");
    q111_20000->SaveAs(Form("%s%d_ql_MAXvol_2ddist_70000.png",file_dir.c_str(),flag));    
    tree->Write();
    ql_MAX_dist_20000->Write();
    ql_MAX_dist_70000->Write();
    Ql_dist->Write();
    stuck1->Write();
    MAX_dist->Write();
    fout->Close();
}
    //TF1*Exp =new TF1("Exp","[0]*exp(-x/[1])",0,50);
    //Exp->SetParameter(0,8000.);
    //Exp->SetParLimits(0,8000.,100000);
    //Exp->SetParameter(1,30.);
    //Exp->SetParLimits(1,30.,100.);
    //stuck->Fit("Exp","","",0,50);

int main(int argc, char *argv[]){
    std::string file_dir =argv[1];
    Int_t Triggerstart_value_1 = 2000;
    Int_t Triggerend_value_1 = 3000;
    Int_t Triggerstart_value_2 = Triggerstart_value_1;
    Int_t Triggerend_value_2 = Triggerend_value_1;
    if(argc==3){
        Triggerstart_value_1 = atoi(argv[2]);
        Triggerend_value_1 = atoi(argv[3]);
    }else if(argc==5){
        Triggerstart_value_1 = atoi(argv[2]);
        Triggerend_value_1 = atoi(argv[3]);
        Triggerstart_value_2 = atoi(argv[4]);
        Triggerend_value_2 = atoi(argv[5]);
    }
    tree(0,Triggerstart_value_1,Triggerend_value_1,file_dir);
    // tree(1,Triggerstart_value_2,Triggerend_value_2,file_dir);
    return 0;
}