//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 21 11:54:29 2024 by ROOT version 6.32.04
// from TTree tree/tree
// found on file: /rhome/maruya/disklink/MPPC_kousei/data/ak_HASHI_HV56p24_0611_1.root
//////////////////////////////////////////////////////////

#ifndef Deri_h
#define Deri_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
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

class Deri {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   Bool_t          fCAEN;

   // Declaration of leaf types
   Int_t           event_number;
   Int_t           adc[64];
   Int_t           CAEN_wavelength;
   Int_t           CAEN_PMTch_0;
   Int_t           CAEN_PMTch_1;
   ULong64_t       TimeStump_8ns;
   // vector<int>     *CAEN_wave0;
   // vector<int>     *CAEN_wave1;
   Double_t        CAEN_wave0_base_mean;
   Double_t        CAEN_wave1_base_mean;
   Double_t        CAEN_wave0_max;
   Double_t        CAEN_wave1_max;
   Int_t           CAEN_wave0_max_index;
   Int_t           CAEN_wave1_max_index;

   // List of branches
   TBranch        *b_event_number;   //!
   TBranch        *b_adc;   //!
   TBranch        *b_CAEN_wavelength;   //!
   TBranch        *b_CAEN_PMTch_0;   //!
   TBranch        *b_CAEN_PMTch_1;   //!
   TBranch        *b_TimeStump_8ns;   //!
   // TBranch        *b_CAEN_wave0;   //!
   // TBranch        *b_CAEN_wave1;   //!
   TBranch        *b_CAEN_wave0_base_mean;   //!
   TBranch        *b_CAEN_wave1_base_mean;   //!
   TBranch        *b_CAEN_wave0_max;   //!
   TBranch        *b_CAEN_wave1_max;   //!
   TBranch        *b_CAEN_wave0_max_index;   //!
   TBranch        *b_CAEN_wave1_max_index;   //!

   Deri( std::string tree_file, Bool_t fCAEN_ar);
   virtual ~Deri();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(config_data config);
   virtual bool     Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Deri_cxx
Deri::Deri( std::string tree_file, Bool_t fCAEN_ar) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TTree* tree = 0;
   fCAEN = fCAEN_ar;
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(tree_file.c_str());
   if (!f || !f->IsOpen()) {
      f = new TFile(tree_file.c_str());
   }
   f->GetObject("tree",tree);
   Init(tree);
}

Deri::~Deri()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Deri::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Deri::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void Deri::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   // CAEN_wave0 = 0;
   // CAEN_wave1 = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event_number", &event_number, &b_event_number);
   fChain->SetBranchAddress("adc", adc, &b_adc);
   if(fCAEN){
      fChain->SetBranchAddress("CAEN_wavelength", &CAEN_wavelength, &b_CAEN_wavelength);
      fChain->SetBranchAddress("CAEN_PMTch_0", &CAEN_PMTch_0, &b_CAEN_PMTch_0);
      fChain->SetBranchAddress("CAEN_PMTch_1", &CAEN_PMTch_1, &b_CAEN_PMTch_1);
      fChain->SetBranchAddress("TimeStump_8ns", &TimeStump_8ns, &b_TimeStump_8ns);
      // fChain->SetBranchAddress("CAEN_wave0", &CAEN_wave0, &b_CAEN_wave0);
      // fChain->SetBranchAddress("CAEN_wave1", &CAEN_wave1, &b_CAEN_wave1);
      fChain->SetBranchAddress("CAEN_wave0_base_mean", &CAEN_wave0_base_mean, &b_CAEN_wave0_base_mean);
      fChain->SetBranchAddress("CAEN_wave1_base_mean", &CAEN_wave1_base_mean, &b_CAEN_wave1_base_mean);
      fChain->SetBranchAddress("CAEN_wave0_max", &CAEN_wave0_max, &b_CAEN_wave0_max);
      fChain->SetBranchAddress("CAEN_wave1_max", &CAEN_wave1_max, &b_CAEN_wave1_max);
      fChain->SetBranchAddress("CAEN_wave0_max_index", &CAEN_wave0_max_index, &b_CAEN_wave0_max_index);
      fChain->SetBranchAddress("CAEN_wave1_max_index", &CAEN_wave1_max_index, &b_CAEN_wave1_max_index);
   }
   Notify();
}

bool Deri::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return true;
}

void Deri::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Deri::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Deri_cxx
