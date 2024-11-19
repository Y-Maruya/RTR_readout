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
#include <vector>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1F.h>
#include <TCanvas.h>

std::vecotor<std::string> split(std::string& input, char deliminater){
    std::istringstream stream(input);
    std::string field;
    std::vector<std::string> result;
    while(getline(stream, field, deliminater)){
        reslut.push_back(field);
    }
    return result;
}

int main(){
    std::string path_run = "/media/kuzelab/HDD/run29/run.txt";
    std::ifstream ifs_run(path_run);
    std::string path_wave = "/media/kuzelab/HDD/run29/wave0.txt";
    std::ifstream ifs_wave(path_wave);
    const int N_event = 360000;
    const int bunch_interval = 1000;//need to optimize
    std::string tmp_run;
    int i_event_run = 0;
    int i_event_wave = 0;
    int timestamp_bunch_run[100] = {};
    int i_event_bunch = 0;
    bool now_in_bunch = false;
    int first_timestamp_bunch_run;
    int first_timestamp_bunch_wave;
    int i_bunch = 0;
    while(std::getline(ifs_run tmp_run)){
        std::vector<std::string> vec_tmp = split(tmp_run, ' ');
        int timestamp_run = std::stoi(vec_tmp[0])*1.E+6 +std::stoi(vec_tmp[1]);//unit: us
        if(now_in_bunch){
            if(timestamp_run - timestamp_bunch_run[i_event_bunch-1] < bunch_interval){
                timestamp_bunch_run[i_event_bunch] = timestamp_run;
                i_event_bunch++;
            }
            else{
                int size_bunch = i_event_bunch;//bunch 
                int timestamp_bunch_wave[size_bunch] = {};
                int i_event_bunch_wave = 1;
                while(std::getline(ifs_wave, tmp_wave)){
                    std::vector<std::string> vec_tmp = split(tmp_wave, ' ');
                    if(vec_temp[2] == "Stamp:"){
                        int timestamp_wave = std::stoi(vec_tmp[3]);//unit: 8 ns
                        if(i_event_bunch_wave == size_bunch){
                            if(i_event_wave == 0){
                                first_timestamp_bunch_wave = timestamp_wave;
                            }
                            else{
                                int delta_timestamp_bunch_run = 125 * (first_timestamp_bunch_run - timestamp_run); //unit: 8ns
                                int delta_timestamp_bunch_wave = first_timestamp_bunch_wave - timestamp_wave;
                            }
                        }
                        else timestamp_bunch_wave[i_event_bunch_wave] = timestamp_wave;
                        i_event_wave++;
                        i_event_bunch_wave++;
                    }
                }


                now_in_bunch = false;
                i_event_bunch = 0;
            }
        }
        else{
            if(i_event_run == 0) first_timestamp_bunch_run = timestamp_run;
            timestamp_bunch[i_event_bunch] = timestamp_run;
            now_in_bunch = true;
            i_event_bunch = 1;
        }
        i_event_run++;
    }
    return EXIT_SUCCESS;
}