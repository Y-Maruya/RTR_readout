#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <numeric>
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
#include <TText.h>
#include <TPad.h>
#include <TTree.h>
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
class config_data{
    public:
        config_data();
        ~config_data();
        void read_config(std::string filename);
        std::string daqstarttime;
        std::string before_laser_val_filename;
        std::string after_laser_val_filename;
        std::string before_dark_val_filename;
        std::string after_dark_val_filename;
        std::string easiroc_data_filename;
        double fVoltage;
        double after_Voltage;
        bool fCAEN;
        bool fCAENbinary;
        std::string fCAENfilename[2];
        int fCAENPMTchannel[2];
        bool fCAMACTDC;
        std::string fCAMACTDCfiledir;
        bool fCAMACADC;
        std::string fCAMACADCfiledir;
};
config_data::config_data(){
    std::string before_laser_val_filename = "laser_HV56p24_903_2.root";
    std::string after_laser_val_filename = "laser_HV56p24_903_2.root";
    std::string before_dark_val_filename = "dark_HV56p24_903_2.root";
    std::string after_dark_val_filename = "dark_HV56p24_903_2.root";
    std::string easiroc_data_filename = "ak-MANNAKA-HV56P24";
    double fVoltage = 56.24;
    double after_Voltage = 56.24;
    bool fCAEN = false;
    bool fCAENbinary = false;
    std::string fCAENfilename[2] = {"",""};
    int fCAENPMTchannel[2] = {0,1};
    bool fCAMACTDC = false;
    std::string fCAMACTDCfiledir = "";
    bool fCAMACADC = false;
    std::string fCAMACADCfiledir = "";
}
config_data::~config_data(){
}
void config_data::read_config(std::string filename){
    std::ifstream ifs(filename);
    if(!ifs){
        std::cerr << "Error: file not opened" << std::endl;
        return;
    }
    std::string str;
    while(getline(ifs,str)){
        std::istringstream iss(str);
        std::string key;
        if(!(iss >> key)){
            continue;
        }
        if(key == "before_laser_val_filename"){
            iss >> before_laser_val_filename;
        }else if(key == "after_laser_val_filename"){
            iss >> after_laser_val_filename;
        }else if(key == "before_dark_val_filename"){
            iss >> before_dark_val_filename;
        }else if(key == "after_dark_val_filename"){
            iss >> after_dark_val_filename;
        }else if(key == "easiroc_data_filename"){
            iss >> easiroc_data_filename;
        }else if(key == "Voltage"){
            iss >> fVoltage;
        }else if(key == "after_Voltage"){
            iss >> after_Voltage;
        }else if(key == "fCAEN"){
            iss >> std::boolalpha >> fCAEN;
        }else if(key == "fCAENbinary"){
            iss >> std::boolalpha >> fCAENbinary;
        }else if(key == "fCAENfilename"){
            iss >> fCAENfilename[0] >> fCAENfilename[1];
        }else if(key == "fCAENPMTchannel"){
            iss >> fCAENPMTchannel[0] >> fCAENPMTchannel[1];
        }else if(key == "fCAMACTDC"){
            iss >> std::boolalpha >> fCAMACTDC;
        }else if(key == "fCAMACADC"){
            iss >> std::boolalpha >> fCAMACADC;
        }
    }
}
void treemake(config_data config){
    std::ifstream hoge(config.easiroc_data_filename,std::ios::binary);
    if(!hoge.is_open()){
        cout << "no file" << endl;
        return;
    }
    std::ifstream ifs_caen[2];
    if(config.fCAEN && config.fCAENbinary){
        ifs_caen[0]=ifstream(config.fCAENfilename[0],std::ios::binary);
        ifs_caen[1]=ifstream(config.fCAENfilename[1],std::ios::binary);
    }else if( config.fCAEN){
        ifs_caen[0]=ifstream(config.fCAENfilename[0]);
        ifs_caen[1]=ifstream(config.fCAENfilename[1]);
    }
    if((!ifs_caen[0] || !ifs_caen[1]) && config.fCAEN){
        std::cerr << "Error: file not opened" << std::endl;
        return;
    }
    std::string outfile = config.easiroc_data_filename + ".root";
    std::cout << "making TTrees : " << outfile << std::endl; 
    TFile*file = new TFile(outfile.c_str(),"RECREATE");
    TTree *tree = new TTree("tree","tree");
    int adc[64];
    int event_number;
    int CAEN_wavelength;
    int CAEN_PMTch_0;
    int CAEN_PMTch_1;
    
    ULong64_t TimeStump_8ns;
    // std::vector<int> CAEN_wave0;
    // std::vector<int> CAEN_wave0_base;
    // std::vector<int> CAEN_wave1;
    // std::vector<int> CAEN_wave1_base;
    double CAEN_wave0_max;
    double CAEN_wave0_base_mean;
    double CAEN_wave1_max;
    double CAEN_wave1_base_mean;
    int CAEN_wave0_max_index;
    int CAEN_wave1_max_index;
    tree->Branch("event_number",&event_number,"event_number/I");
    tree->Branch("adc",adc,"adc[64]/I");
    if(config.fCAEN){
        tree->Branch("CAEN_wavelength",&CAEN_wavelength,"CAEN_wavelength/I");
        tree->Branch("CAEN_PMTch_0",&CAEN_PMTch_0,"CAEN_PMTch_0/I");
        tree->Branch("CAEN_PMTch_1",&CAEN_PMTch_1,"CAEN_PMTch_1/I");
        tree->Branch("TimeStump_8ns",&TimeStump_8ns,"TimeStump_8ns/l");
        // tree->Branch("CAEN_wave0",&CAEN_wave0);
        // tree->Branch("CAEN_wave1",&CAEN_wave1);
        tree->Branch("CAEN_wave0_base_mean",&CAEN_wave0_base_mean,"CAEN_wave0_base_mean/D");
        tree->Branch("CAEN_wave1_base_mean",&CAEN_wave1_base_mean,"CAEN_wave1_base_mean/D");
        tree->Branch("CAEN_wave0_max",&CAEN_wave0_max,"CAEN_wave0_max/D");
        tree->Branch("CAEN_wave1_max",&CAEN_wave1_max,"CAEN_wave1_max/D");
        tree->Branch("CAEN_wave0_max_index",&CAEN_wave0_max_index,"CAEN_wave0_max_index/I");
        tree->Branch("CAEN_wave1_max_index",&CAEN_wave1_max_index,"CAEN_wave1_max_index/I");
    }
    int CAMACTDC[8];
    int CAMACTDC_Time;
    int CAMACTDC_uTime;
    if(config.fCAMACTDC){
        tree->Branch("CAMACTDC",CAMACTDC,"CAMACTDC[8]/I");
        tree->Branch("CAMACTDC_Time",&CAMACTDC_Time,"CAMACTDC_Time/I");
        tree->Branch("CAMACTDC_uTime",&CAMACTDC_uTime,"CAMACTDC_uTime/I");
    }
    int CAMACADC[8];
    int CAMACADC_Time;
    int CAMACADC_uTime;
    if(config.fCAMACADC){
        tree->Branch("CAMACADC",CAMACADC,"CAMACADC[8]/I");
        tree->Branch("CAMACADC_Time",&CAMACADC_Time,"CAMACADC_Time/I");
        tree->Branch("CAMACADC_uTime",&CAMACADC_uTime,"CAMACADC_uTime/I");
    }
    event_number = 0;
    int j =0;
    while(!hoge.eof() && (((!ifs_caen[0].eof() || !ifs_caen[1].eof()) && config.fCAEN)|| !config.fCAEN) ){
            std::cout << "event_number : " << event_number << std::endl;
        
        UInt_t val;
        hoge.read((char*)&val, sizeof(int));
        if(config.fCAEN && config.fCAENbinary){
            std::cout << "not supported" << std::endl;
        //     CAEN_PMTch_0 = config.fCAENPMTchannel[0];
        //     CAEN_PMTch_1 = config.fCAENPMTchannel[1];
        //     char buf[2];
        //     ifs_caen[0].read(buf, 2);
        //     unsigned int x = (buf[0] << 24); x = x >> 24;
        //     unsigned int y = (buf[1] << 24); y = y >> 16;
        //     int num = x + y;

        //     while(j<CAEN_wavelength){

        //     }
        //     if(j<12){
        //         j+=1;
        //     }
        }
        if(config.fCAEN && !config.fCAENbinary){
            CAEN_PMTch_0 = config.fCAENPMTchannel[0];
            CAEN_PMTch_1 = config.fCAENPMTchannel[1];
            std::string buf;
            ifs_caen[0]>>buf>>buf>>CAEN_wavelength;
            ifs_caen[0]>>buf>>buf;
            ifs_caen[0]>>buf>>buf;
            int event_number_CAEN;
            ifs_caen[0]>>buf>>buf>>event_number_CAEN;
            // std::cout << event_number_CAEN << std::endl;
            if(event_number_CAEN != event_number){
                std::cerr << "Error: event number is not matched" << std::endl;
                return;
            }
            getline(ifs_caen[0], buf);
            ifs_caen[0]>>buf>>buf>>buf>>buf>>buf>>TimeStump_8ns;
            getline(ifs_caen[0], buf);
            getline(ifs_caen[0], buf);
            int j=0;
            std::string str_data;
            int th = 10;
            int wave0_min =4096;
            int wave0_min_index = 0;
            int sum0_base = 0;
            int sum0_base_num = 0;
            while(j<CAEN_wavelength){
                int maxadcd;
                getline(ifs_caen[0],str_data);
                stringstream ss(str_data);
                ss>>maxadcd;
                // CAEN_wave0.push_back(maxadcd);
                if(maxadcd<wave0_min){
                    wave0_min = maxadcd;
                    wave0_min_index = j;
                }
                j++;
                if(maxadcd<2048+th &&maxadcd>2048-th && j < CAEN_wavelength*3./10){
                    // CAEN_wave0_base.push_back(maxadcd);
                    sum0_base += maxadcd;
                    sum0_base_num++;
                }
            }
            ifs_caen[1]>>buf>>buf>>CAEN_wavelength;
            if(CAEN_wavelength != CAEN_wavelength){
                std::cerr << "Error: wavelength is not matched" << std::endl;
                return;
            }
            ifs_caen[1]>>buf>>buf;
            ifs_caen[1]>>buf>>buf;
            ifs_caen[1]>>buf>>buf>>event_number_CAEN;
            if(event_number_CAEN != event_number){
                std::cerr << "Error: event number is not matched" << std::endl;
                return;
            }
            getline(ifs_caen[1], buf);
            ifs_caen[1]>>buf>>buf>>buf>>buf>>buf>>TimeStump_8ns;
            if(TimeStump_8ns != TimeStump_8ns){
                std::cerr << "Error: TimeStump is not matched" << std::endl;
                return;
            }
            getline(ifs_caen[1], buf);
            getline(ifs_caen[1], buf);
            j=0;
            int wave1_min =4096;
            int wave1_min_index = 0;
            int sum1_base = 0;
            int sum1_base_num = 0;
            while(j<CAEN_wavelength){
                int maxadcd;
                getline(ifs_caen[1],str_data);
                stringstream ss(str_data);
                ss>>maxadcd;
                if(maxadcd<wave1_min){
                    wave1_min = maxadcd;
                    wave1_min_index = j;
                }
                // CAEN_wave1.push_back(maxadcd);
                j++;
                if(maxadcd<2048+th &&maxadcd>2048-th && j < CAEN_wavelength*3./10){
                    // CAEN_wave1_base.push_back(maxadcd);
                    sum1_base += maxadcd;
                    sum1_base_num++;
                }
            }
            CAEN_wave0_base_mean = sum0_base/sum0_base_num;
            CAEN_wave1_base_mean = sum1_base/sum1_base_num;
            CAEN_wave0_max_index = wave0_min_index;
            CAEN_wave1_max_index = wave1_min_index;
            CAEN_wave0_max = CAEN_wave0_base_mean - wave0_min;
            CAEN_wave1_max = CAEN_wave1_base_mean - wave1_min;
        }
        if(config.fCAMACTDC){
            std::string buf;
            std::ifstream ifs_camactdc(Form("%s/%d.txt",config.fCAMACADCfiledir,event_number));
            if(!ifs_camactdc){
                std::cerr << "Error: file not opened" << std::endl;
                return;
            }
            ifs_camactdc>>buf>>CAMACTDC_Time>>CAMACTDC_uTime>>CAMACTDC[0]>>CAMACTDC[1]>>CAMACTDC[2]>>CAMACTDC[3]>>CAMACTDC[4]>>CAMACTDC[5]>>CAMACTDC[6]>>CAMACTDC[7];
            ifs_camactdc.close();
        }
        if(config.fCAMACADC){
            std::string buf;
            std::ifstream ifs_camacadc(Form("%s/%d.txt",config.fCAMACADCfiledir,event_number));
            if(!ifs_camacadc){
                std::cerr << "Error: file not opened" << std::endl;
                return;
            }
            ifs_camacadc>>buf>>CAMACADC_Time>>CAMACADC_uTime>>CAMACADC[0]>>CAMACADC[1]>>CAMACADC[2]>>CAMACADC[3]>>CAMACADC[4]>>CAMACADC[5]>>CAMACADC[6]>>CAMACADC[7];
            ifs_camacadc.close();
        }
        if(val == 0xffffea0c){
            hoge.read((char*)&val, sizeof(int));
            event_number++;
            for(int i = 0; i<65; ++i){
                hoge.read((char*)&val, sizeof(int));
                if(i>0){
                    int buffer = val & 0xffff;
                    adc[i-1] = buffer;
                }
            }
        }
        tree->Fill();
    }
    tree->Write();
    file->Close();
}
int main(int argc, char** argv){
    config_data config;
    config.read_config(std::string(argv[1]));
    treemake(config);
    return 0;
}