void ehist64(string filename = "ak-MANNAKA-HV56P24"){
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
  TH1F *hist[64];
  TH2F* chhist = new TH2F("chhist","",64,0-0.5,64-0.5,1500,500,2000);
  for(int i = 0; i<64; ++i){
    std::ostringstream _name;
    _name << "ch" << i;
    std::string hist_name = _name.str();
    hist[i] = new TH1F(hist_name.c_str(), hist_name.c_str(), 4096, 0, 4096);   
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
 
  for(int i=0;i<32;i++){
    mean = hist[i]->GetMean();
    //RMS = hist[i]->GetRMS();
    std::cout<<mean<<std::endl;
    //std::cout<<"ch"<<i<<" mean = "<<mean<<"  RMS = "<<RMS<<std::endl;
  }
 for(int i=0;i<32;i++){
   //mean = hist[i]->GetMean();
    RMS = hist[i]->GetRMS();
    std::cout<<RMS<<std::endl;
    //std::cout<<"ch"<<i<<" mean = "<<mean<<"  RMS = "<<RMS<<std::endl;
  }
 file->Write();
 file->Close();
  return;
}
