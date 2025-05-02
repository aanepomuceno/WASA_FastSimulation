#include <iostream>  
#include <fstream>
#include<string>  
#include "TVector3.h"


void muon_analysis(){

   const char *filename = "";
   TFile *f = NULL;
   const char *tName = "";
   TTree *t1 = NULL;
   TH1F *hem  = new TH1F("EMCAL","Muon Energy",60,0,240);
   TVector3 v1; 
   TVector3 v2; 
   double em = 0;
   double x=0; double y=0; double z = 0;
   int nfiles = 4;   
   for (int ifile =0; ifile < nfiles; ifile++) {
          std::string str1 = "WASAFastOutput_t"+to_string(ifile)+".root";
          filename = str1.c_str();
          f = new TFile(filename);
  
     for (int i=0;i<10000;i++) {
      std::string str2 = "Event_"+to_string(i);
      tName = str2.c_str();
      t1 = (TTree*)f->Get(tName);
      if (t1 == NULL) continue;
      t1->SetBranchAddress("emcal_E",&em);
      t1->SetBranchAddress("emcal_X",&x);
      t1->SetBranchAddress("emcal_Y",&y);
      t1->SetBranchAddress("emcal_Z",&z);
      if (t1->GetEntries() < 2) continue;
        for (int j=0; j< t1->GetEntries(); j++) {
        t1->GetEntry(j);
        if (em < 1.0 ) continue;
        hem->Fill(em);
        }
     }
  }

    TCanvas * c1 = new TCanvas("c1", "c1", 600, 500);
    hem->Draw();
   
 }
