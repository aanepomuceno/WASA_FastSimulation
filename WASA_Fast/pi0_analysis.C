#include <iostream>  
#include <fstream>
#include<string>  
#include "TVector3.h"

void pi0_analysis(){

   const char *filename = "";
   TFile *f = NULL;
   const char *tName = "";
   TTree *t1 = NULL;
   TH1F *hmass  = new TH1F("InvMass","InvMass",50,0,250);
   TH1F *hem  = new TH1F("EMCAL","Photon Energy",50,20,500);
   TH1F *hpz  = new TH1F("MC_pz","MC Pz",100,-1000,1000);
   TH1F *hmcE  = new TH1F("MC_E","MC E",100,0,1.);
   TH2  *hE2d = new TH2F("E2D", "E2D", 50, 0, 400, 50, 0, 400);
   TH2* hmc = new TH2F("h2", "MC Truth", 50, -3.00, 3.00, 50, -3.00, 3.00);
   TVector3 v1; 
   TVector3 v2; 
   int nfiles = 4;
   double em = 0;
   double emsum = 0;
   double mc_x = 0;
   double mc_y = 0;
   double mc_z = 0;
   double pdg = 0;
   double p = 0;
   double x=0; double y=0; double z = 0;
   double angle = 0;
   double mass = 0;
   vector<double> v;
   vector<double> energy;
   Int_t ev;

   for (int ifile =0; ifile < nfiles; ifile++) {
          std::string str1 = "WASAFastOutput_t"+to_string(ifile)+".root";
          filename = str1.c_str();
       f = new TFile(filename);
  
     for (int i=0;i<10000;i++) {
   //cout << i << endl;
     std::string str2 = "Event_"+to_string(i);
     tName = str2.c_str();
     t1 = (TTree*)f->Get(tName);
     if (t1 == NULL) continue;
     t1->SetBranchAddress("emcal_E",&em);
     t1->SetBranchAddress("emcal_X",&x);
     t1->SetBranchAddress("emcal_Y",&y);
     t1->SetBranchAddress("emcal_Z",&z);
     t1->SetBranchAddress("emcal_PDG",&pdg);
     t1->SetBranchAddress("MC_pX",&mc_x);
     t1->SetBranchAddress("MC_pY",&mc_y);
     t1->SetBranchAddress("MC_pZ",&mc_z);
     emsum = 0;
     if (t1->GetEntries() < 3) continue;
     v.clear();
     v1.Clear();
     v2.Clear();
     energy.clear();
     for (int j=0; j< t1->GetEntries(); j++) {
     t1->GetEntry(j);
     if (em < 1.0 ) continue;
     emsum += em;
     v.push_back(x);
     v.push_back(y);
     v.push_back(z);
     energy.push_back(em);
     }
      t1->DropBaskets();  
      delete t1;
     v1.SetXYZ(v[0],v[1],v[2]);
     v2.SetXYZ(v[3],v[4],v[5]);
     angle = v1.Angle(v2);
     mass = sqrt( 2*energy[0]*energy[1]*(1.-cos(angle) ) );
     p = sqrt(pow(mc_x,2) + pow(mc_y,2) + pow(mc_z,2) );
     hem->Fill(emsum);
     hmass->Fill(mass);
     hE2d->Fill(energy[0],energy[1]);
     hpz->Fill(mc_z);
     hmc->Fill(mc_x,mc_y);
     hmcE->Fill(p/1000.);
    
     }
     f->Close();
     delete f;
  }

   TCanvas * c1 = new TCanvas("c1", "c1", 600, 500);
   hE2d->Draw();
   TCanvas * c2 = new TCanvas("c2", "c2", 600, 500);
   hmass->Draw();


   
 }
