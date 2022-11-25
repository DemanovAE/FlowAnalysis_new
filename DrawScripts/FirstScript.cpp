// C++ headers
#include <string>
#include <vector>
#include <iostream>
#include <fstream> 
#include <map>

// ROOT headers
#include <TStopwatch.h>
#include "TProfile.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TF2.h"
#include "TF1.h"

#include "../macro/ConstantsNew.h"
#include "./ArticlePoint.C"

//#include "./Draw_Picture_new.h"
void AddHistoAllRun2(TFile *f1, TFile *out, std::vector<int> BadRunId);
Bool_t DuplicateRunId(std::vector<int> vec, int r);
void makeplotstyle();


void FirstScript(const Char_t *inFileName = "flow_27GeV_hadd.root",
                 const Char_t *outFileName = "Flow27GeV.root",
                 const Int_t energy = 28){
	
	std::map<int, string> Energy_scan_text={{7, "7.7"},{11, "11.5"},{14, "14.5"},{19, "19.6"},{27, "27"},{39, "39"},{62, "62.4"}};
	string particle[]={"#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","Hadrons","","PID"};
	string bukva[] = {"(a)","(b)","(c)","(d)","(i)", "(f)", "(g)","(h)", "(i)", "(j)","(k)","(l)","(m)","(n)","(o)","(p)", "(q)"};
	int cent[10]={80,70,60,50,40,30,20,10,5,0};
 
  TFile *fileOut1 = new TFile(outFileName,"recreate");
  TFile *fileIn1 = new TFile(inFileName,"read");

  AddHistoAllRun2(fileIn1, fileOut1, BadRuns.at(energy));
	
  WritePubl();

  fileOut1->Close();
  fileIn1->Close();
}

void AddHistoAllRun2(TFile *f1, TFile *out, std::vector<int> BadRunId){
  
  std::cout<<"Start Merging!"<<std::endl;
  TStopwatch timer1;
  timer1.Start();

  std::vector<int> *aa;
  std::vector<TString> *bb;
  std::vector<int> RunId;
  std::vector<TString> NameObject;
  
  ///////////////////// Find RunId ///////////////////////////////
  for(Int_t aaii=1; aaii<2000; aaii++){
    f1->GetObject(Form("VecRunId;%i",aaii),aa);
    if(aa!=nullptr){
      for(Int_t iv=0; iv<aa[0].size(); iv++){
        if(DuplicateRunId(RunId,aa[0][iv])==true) continue;
        if( std::find(BadRunId.begin(), BadRunId.end(), aa[0][iv] ) != BadRunId.end() ) continue;
        RunId.push_back(aa[0][iv]);
      }
    }else{
      aaii=2001;
    }
  }
  for(Int_t i=0; i<RunId.size(); i++){
    std::cout<<"Run: "<<RunId[i]<<"\tN= "<<i+1<<std::endl;
  }

  ///////////////////// Find NameInFile ////////////////////////////
  f1->GetObject(Form("VecNameInFile;1"),bb);
  if(bb==nullptr)std::cout<<"Not find: VecNameInFile"<<std::endl;
  for(Int_t iv=0; iv<bb[0].size(); iv++){
    NameObject.push_back(bb[0][iv]);
  }

  //////////////////// Hadd object /////////////////////////////////
  out->cd();
  for(Int_t in=0; in<NameObject.size(); in++ ){
    if( strncmp(NameObject[in], "tp_",3) == 0){
      auto *write = (TProfile*)(f1->FindObjectAny( Form("%s_Run%i",NameObject[in].Data(),RunId[0])) );
      for(Int_t ir=1; ir<RunId.size(); ir++){
        write->Add( (TProfile*)(f1->FindObjectAny( Form("%s_Run%i",NameObject[in].Data(),RunId[ir]))) );
      }
      write->SetName(NameObject[in]);
      write->Write();
      std::cout<<"Write:\t"<<write->GetName()<<std::endl;
      delete write;
    }
    
    if( strncmp(NameObject[in], "tp2_",4) == 0){
      auto *write = (TProfile2D*)(f1->FindObjectAny( Form("%s_Run%i",NameObject[in].Data(),RunId[0])) );
      for(Int_t ir=1; ir<RunId.size(); ir++){
        write->Add( (TProfile2D*)(f1->FindObjectAny( Form("%s_Run%i",NameObject[in].Data(),RunId[ir]))) );
      }
      write->SetName(NameObject[in]);
      write->Write();
      std::cout<<"Write:\t"<<write->GetName()<<std::endl;
      delete write;
    }

    if( strncmp(NameObject[in], "tp3_",4) == 0){
      auto *write = (TProfile3D*)(f1->FindObjectAny( Form("%s_Run%i",NameObject[in].Data(),RunId[0])) );
      for(Int_t ir=1; ir<RunId.size(); ir++){
        write->Add( (TProfile3D*)(f1->FindObjectAny( Form("%s_Run%i",NameObject[in].Data(),RunId[ir]))) );
      }
      write->SetName(NameObject[in]);
      write->Write();
      std::cout<<"Write:\t"<<write->GetName()<<std::endl;
      delete write;
    }

  }

  //std::cout<<"Merging completed! Numb RunId: "<<first<<std::endl;
  timer1.Stop();
  timer1.Print();

}

Bool_t DuplicateRunId(std::vector<int> vec, int r){
  for(Int_t i=0; i<vec.size(); i++){
    if(vec[i]==r)return true;
  }
  return false;
}

void makeplotstyle(){
  
    TStyle *mystyle = new TStyle("PlottingInStyle", "Style for Summary Plots");
    mystyle->SetLineWidth(2);
    mystyle->SetPalette(1);
    mystyle->SetCanvasColor(10);
    mystyle->SetHistFillColor(10);
    mystyle->SetHistFillStyle(0);
    mystyle->SetOptTitle(0);
    mystyle->SetOptStat(0);
    mystyle->SetCanvasBorderMode(0);//removes the yellow frame around the canvas
    mystyle->SetPadLeftMargin(0.16);
    mystyle->SetPadBottomMargin(0.15);
    mystyle->SetPadTickX(2);
    mystyle->SetPadTickY(2);
    mystyle->SetAxisColor(1, "X");
    mystyle->SetAxisColor(1, "Y");
    mystyle->SetLabelColor(1, "X");
    mystyle->SetLabelColor(1, "Y");
    mystyle->SetTickLength(0.03, "X");
    mystyle->SetTickLength(0.03, "Y");
    mystyle->SetTitleXSize(0.05);
    mystyle->SetTitleYSize(0.05);
    mystyle->SetNdivisions(508, "X");
    mystyle->SetNdivisions(505, "Y");
    mystyle->SetTitleXOffset(1.2);
    mystyle->SetTitleYOffset(1.4);
    mystyle->SetLabelOffset(0.02, "X");
    mystyle->SetLabelOffset(0.02, "Y");
    mystyle->SetLabelSize(0.05, "X");
    mystyle->SetLabelSize(0.05, "Y");
    mystyle->SetEndErrorSize(8);
    //mystyle->SetGridx();

    TFile f("style.root", "RECREATE");
    f.cd();
    mystyle->Write();
    f.Close();
}