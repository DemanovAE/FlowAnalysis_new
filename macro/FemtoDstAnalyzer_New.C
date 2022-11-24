/**
 * \brief Example of how to read a file (list of files) using StFemtoEvent classes
 *
 * RunFemtoDstAnalyzer.C is an example of reading FemtoDst format.
 * One can use either FemtoDst file or a list of femtoDst files (inFile.lis or
 * inFile.list) as an input, and preform physics analysis
 *
 * \author Grigory Nigmatkulov, Povarov Alexey, Demanov Alexandr
 * \date May 29, 2018
 */

// This is needed for calling standalone classes
#define _VANILLA_ROOT_

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

// FemtoDst headers

#include "StFemtoDstReader.h"
#include "StFemtoEpdHit.h"
#include "StFemtoDst.h"
#include "StFemtoEvent.h"
#include "StFemtoTrack.h"
#include "StFemtoV0.h"
#include "StFemtoXi.h"

#include "ConstantsNew.h"
#include "functions.C"
#include "FlowEtaSubEP.h"
#include "FlowScalarProduct.h"
#include "StEpdGeom.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStFemtoDst.so)
R__LOAD_LIBRARY(libFlowAnalysis.so)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

// Used function
Bool_t isGoodEvent(StFemtoEvent *const &event, const Int_t _energy, const Bool_t BadRunIdKeyFlowStage);
Bool_t isGoodTrackEP(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);
void GetVectorRunId(const Char_t *inFile, std::vector<int> &runIdVec);
int PID_TPC_TOF(StFemtoTrack *const &track, Double_t DCA, const Int_t _energy);
void WriteCheckHisto2DForAllRun(std::vector<std::vector<FlowEtaSubEP*>> ClassSubEP, Int_t EW);
void WriteCheckHisto1DForAllRun(std::vector<std::vector<FlowEtaSubEP*>> ClassSubEP, Int_t EW);
Double_t GetRapidity2(StFemtoTrack *const &track, Int_t iPID);

//_________________
void FemtoDstAnalyzer_New(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                          const Char_t *outFileName = "oTest.root",
                          const Char_t *FileNameAfterFirstRun = "iRec.root",
                          const Char_t *FileNameAfterSecondRun = "iFlat.root",
                          const Char_t *mode = "QAmode",
                          const Int_t energy = 39) {

  TStopwatch timer1;
  timer1.Start();
  
  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
    gSystem->Load("libFlowAnalysis.so");
  #endif

  StFemtoDstReader* femtoReader = new StFemtoDstReader(inFile);
  femtoReader->Init();

  // This is a way if you want to spead up IO
  std::cout << "Explicit read status for some branches" << std::endl;
  femtoReader->SetStatus("*",0);
  femtoReader->SetStatus("Event",1);
  femtoReader->SetStatus("Track",1);
  femtoReader->SetStatus("EpdHit",1);
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  std::vector<int> vRunId;
  std::map<int, int> mRunId;
  TFile *outFile = new TFile(outFileName, "RECREATE");
  TFile *FileRec, *FileFlow;
  outFile->cd();
  // Eta-sub Event Plane
  std::vector<std::vector<FlowEtaSubEP*>> flowEtaSub;
  std::vector<std::vector<FlowScalarProduct*>> flowSP;
  // Mode works
  Bool_t mode_raw = false;
  Bool_t mode_rec = false;
  Bool_t mode_flow = false;
  Bool_t CalFlow = false;
  Bool_t CheckHistoQvectorAndPsi = false;
  Bool_t CalEPFromTPC = true;
  Bool_t CalEPFromEPD = false;

  if( strncmp(mode, "raw",3) == 0){
    mode_raw = true;
  }
  if( strncmp(mode, "rec",3) == 0){
    mode_rec = true;
  }
  if( strncmp(mode, "flow",4) == 0){
    CalFlow = true;
  }

  GetVectorRunId(inFile,vRunId);
  for(Int_t i=0; i<(int)vRunId.size();i++){
    mRunId.insert(std::make_pair( vRunId[i], i));
    std::cout<<vRunId[i]<<"\t map"<<mRunId.at(vRunId[i])<<"\n";
  }

  flowEtaSub.resize(Nharm);
  for(Int_t ih=0; ih<Nharm; ih++){
    flowEtaSub[ih].resize((int)vRunId.size());
    for(Int_t ir=0; ir<(int)vRunId.size(); ir++){
      flowEtaSub[ih][mRunId.at(vRunId[ir])] = new FlowEtaSubEP();
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetFirstRun(mode_raw);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetSecondRun(mode_rec);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetCalFlow(CalFlow);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetChekQVectorAndPsi(CheckHistoQvectorAndPsi);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetHarmonic(vec_harmonic.at(ih));
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetEtaGap(EtaGapSubEP);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetRunId(vRunId[ir]);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetNCentBins(nBinCent);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetNBinsVtxZ(nBinVtxZ);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetVtxZ(CutVtxZ.at(energy));
      if(ih==3){
        flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetNameSys("v42");
        flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetHarmonicNK(4,2);
      }
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetInputFileFromFirstRun((TString)FileNameAfterFirstRun);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetMeanQVector();
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetNumbFlattening(4);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetInputFileFromSecondRun((TString)FileNameAfterSecondRun);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetMeanSinCosPsi();
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->SetMapPidCode(PIDmap);
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->Init();
    }
  }

  flowSP.resize(Nharm);
  for(Int_t ih=0; ih<Nharm; ih++){
    flowSP[ih].resize((int)vRunId.size());
    for(Int_t ir=0; ir<(int)vRunId.size(); ir++){
      flowSP[ih][mRunId.at(vRunId[ir])] = new FlowScalarProduct();
      flowSP[ih][mRunId.at(vRunId[ir])]->SetFirstRun(mode_raw);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetCalFlow(CalFlow);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetChekQVectorAndPsi(CheckHistoQvectorAndPsi);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetHarmonic(vec_harmonic.at(ih));
      flowSP[ih][mRunId.at(vRunId[ir])]->SetEtaGap(EtaGapSubEP);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetRunId(vRunId[ir]);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetNCentBins(nBinCent);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetNBinsVtxZ(nBinVtxZ);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetVtxZ(CutVtxZ.at(energy));
      flowSP[ih][mRunId.at(vRunId[ir])]->SetNameSys("SP");
      flowSP[ih][mRunId.at(vRunId[ir])]->SetInputFileFromFirstRun((TString)FileNameAfterFirstRun);
      flowSP[ih][mRunId.at(vRunId[ir])]->SetMeanQVector();
      flowSP[ih][mRunId.at(vRunId[ir])]->SetMapPidCode(PIDmap);
      flowSP[ih][mRunId.at(vRunId[ir])]->Init();
    }
  }

  StEpdGeom *EpdGeom = new StEpdGeom(); 

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    if ( iEvent % 10000 == 0){std::cout << "Working on event #[" << (iEvent+1) << "/" << events2read << "]" << std::endl;}
    if (iEvent == events2read-1) {std::cout << "Working on event #[" << (events2read) << "/" << events2read << "]" << std::endl;}

    Bool_t readEvent = femtoReader->readFemtoEvent(iEvent);
    if( !readEvent ) {
      std::cout << "Something went wrong, Master! Nothing to analyze..." << std::endl;
      break;
    }

    // Retrieve femtoDst
    StFemtoDst *dst = femtoReader->femtoDst();

    // Retrieve event information
    StFemtoEvent *event = dst->event();
    if( !event ) {
      std::cout << "Something went wrong, Master! Event is hiding from me..." << std::endl;
      break;
    }

    if( !isGoodEvent( event,  energy, false) ) continue;
    if( !TofMatchedCut(dst, event, 4, energy) ) continue;

    Int_t RunID = event -> runId();
    Int_t cent=0; 
    if( nBinCent == 9 ){
      cent = event -> cent9();
    }
    if( nBinCent == 16 ){
      cent = event -> cent16();
    }

    for(Int_t iharm=0; iharm<Nharm; iharm++){
      flowEtaSub[iharm][mRunId.at(RunID)]->Zero();
      flowSP[iharm][mRunId.at(RunID)]->Zero();
    }

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
      
    if(CalEPFromTPC){
      // Track loop
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
        // Retrieve i-th femto track
        StFemtoTrack *femtoTrack = dst->track(iTrk);

        if( !isGoodTrackEP(event, femtoTrack, energy) ) continue; 

        //if(iEvent==1001)std::cout<<femtoTrack -> eta()<<"\t"<<femtoTrack -> phi()<<"\t"<<femtoTrack -> pt()<<std::endl;

        for(Int_t iharm=0; iharm<Nharm; iharm++){
          //flowEtaSub[iharm][mRunId.at(RunID)]->ProcessFirstTrackLoop(femtoTrack->eta(),femtoTrack->phi(),GetWeight(femtoTrack));
          flowEtaSub[iharm][mRunId.at(RunID)]->ProcessFirstTrackLoop(femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt());
          flowSP[iharm][mRunId.at(RunID)]    ->ProcessFirstTrackLoop(femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt());
        }

      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    }// TPC EP determ

    if(CalEPFromEPD){

      Int_t nHitEpd = dst->numberOfEpdHits();
      TVector3 PrimaryVertex = event->primaryVertex();

      for(Int_t iHit=0; iHit<nHitEpd; iHit++) {
        StFemtoEpdHit *femtoEpdHit = dst->epdHit(iHit);
        if( !femtoEpdHit->isGood() ) continue; 

        TVector3 EpdPoint = EpdGeom->RandomPointOnTile(femtoEpdHit->id());
        TVector3 StraightLine = EpdPoint-PrimaryVertex;
        double nMip = femtoEpdHit->nMIP();
        double nMipEff = nMip; //put a max on nMip in calc.
        if(nMipEff > 2.) nMipEff = 2.;
        else if (nMipEff < 0.3) continue;
        double RandomEta = StraightLine.Eta();
        double RandomPhi = StraightLine.Phi(); while(RandomPhi < 0.) RandomPhi += 2.*pi;
        int ew = (femtoEpdHit->id()<0)?0:1;  //is EPD east or west
        if (!(TMath::Abs(RandomEta)>0)||(TMath::Abs(RandomEta)>5.0) || (TMath::Abs(RandomEta)<2.1)) continue;
        //std::cout<<ew<<"\t"<<RandomEta<<"\t"<<RandomPhi<<"\t"<<nMip<<std::endl;
        for(Int_t iharm=0; iharm<Nharm; iharm++){
          //flowEtaSub[iharm][mRunId.at(RunID)]->ProcessFirstTrackLoop(femtoTrack->eta(),femtoTrack->phi(),GetWeight(femtoTrack));
          flowEtaSub[iharm][mRunId.at(RunID)]->ProcessFirstTrackLoop(RandomEta,RandomPhi,nMip);
          flowSP[iharm][mRunId.at(RunID)]->ProcessFirstTrackLoop(RandomEta,RandomPhi,nMip);
        }
      
      }

    }// EPD EP Determ

    for(Int_t iharm=0; iharm<Nharm; iharm++){
      flowEtaSub[iharm][mRunId.at(RunID)]->ProcessEventAfterFirstTrackLoop((Double_t)cent,(event->primaryVertex()).Z());
      flowSP[iharm][mRunId.at(RunID)]->ProcessEventAfterFirstTrackLoop((Double_t)cent,(event->primaryVertex()).Z());
    }

    if(CalFlow==true){
      //Second Track loop
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {

        // Retrieve i-th femto track
        StFemtoTrack *femtoTrack = dst->track(iTrk);
        
        if( !isGoodTrackFlow(event, femtoTrack, energy) ) continue; 

        Int_t parTPCandTOF = PID_TPC_TOF(femtoTrack, femtoTrack -> gDCA(event->primaryVertex()).Mag(), energy);      
        Int_t parHadrons = -1;
        if(femtoTrack->charge()>0) parHadrons=0;
        if(femtoTrack->charge()<0) parHadrons=4;
        
        if(  parTPCandTOF < 0 || parHadrons < 0) continue;

        for(Int_t iharm=0; iharm<Nharm; iharm++){
          flowEtaSub[iharm][mRunId.at(RunID)]->ProcessSecondTrackLoop(femtoTrack->eta(),femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt(),parHadrons,(Double_t)cent);
          flowEtaSub[iharm][mRunId.at(RunID)]->ProcessSecondTrackLoop(femtoTrack->eta(),femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt(),parTPCandTOF,(Double_t)cent);
          //flowEtaSub[iharm][mRunId.at(RunID)]->ProcessSecondTrackLoop(GetRapidity2(femtoTrack, parTPCandTOF),femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt(),parTPCandTOF,(Double_t)cent);

          flowSP[iharm][mRunId.at(RunID)]->ProcessSecondTrackLoop(femtoTrack->eta(),femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt(),parHadrons,(Double_t)cent);
          flowSP[iharm][mRunId.at(RunID)]->ProcessSecondTrackLoop(femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt(),parTPCandTOF,(Double_t)cent);
          //flowSP[iharm][mRunId.at(RunID)]->ProcessSecondTrackLoop(GetRapidity2(femtoTrack, parTPCandTOF),femtoTrack->eta(),femtoTrack->phi(),femtoTrack->pt(),parTPCandTOF,(Double_t)cent);

        }

      } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
    }
  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  
  outFile->cd();
  for(Int_t ih=0; ih<Nharm; ih++){
    for(Int_t ir=0; ir<(int)vRunId.size(); ir++){
      flowEtaSub[ih][mRunId.at(vRunId[ir])]->WtiteHist();
      flowSP[ih][mRunId.at(vRunId[ir])]->WtiteHist();
    }
  }

  if(CheckHistoQvectorAndPsi){
    WriteCheckHisto2DForAllRun(flowEtaSub, -1);
    WriteCheckHisto2DForAllRun(flowEtaSub, +1);
    WriteCheckHisto1DForAllRun(flowEtaSub, -1);
    WriteCheckHisto1DForAllRun(flowEtaSub, +1);
  }  

  
  std::vector<TString> VecNameInFIle;
  for(Int_t ih=0; ih<Nharm; ih++){
    for(Int_t in=0; in<(flowSP[ih][0]->GetNameInFile()).size(); in++ ){
      VecNameInFIle.push_back((flowSP[ih][0]->GetNameInFile())[in]);
    }
    for(Int_t in=0; in<(flowEtaSub[ih][0]->GetNameInFile()).size(); in++ ){
      VecNameInFIle.push_back((flowEtaSub[ih][0]->GetNameInFile())[in]);
    }
  }
  outFile->WriteObjectAny(&VecNameInFIle,"std::vector<TString>","VecNameInFile");
  outFile->WriteObjectAny(&vRunId,"std::vector<int>","VecRunId");
  outFile->Close();

  femtoReader->Finish();

  for(Int_t ih=0; ih<Nharm; ih++){
    for(Int_t ir=0; ir<(int)vRunId.size(); ir++){
      delete flowEtaSub[ih][ir];
      delete flowSP[ih][ir];
    }
  }
  flowEtaSub.clear();
  flowSP.clear();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;

  timer1.Stop();
  timer1.Print();
}

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________DESCRIPTION OF FUNCTIONS_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

//********************CHECK EVENT ON GOOD********************//
Bool_t isGoodEvent(StFemtoEvent *const &event, const Int_t _energy, const Bool_t BadRunIdKeyFlowStage) {
  
  TVector3 pVtx = event->primaryVertex();
  // Reject vertices that are far from the central membrane along the beam
  if( TMath::Abs( pVtx.Z() ) > CutVtxZ.at(_energy) ) return false;
  //if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + 0.8847, 2) ) > 1 ) check = false; 14.5
  if( sqrt( pow( pVtx.X(), 2) + pow(pVtx.Y() + CutDeltaVtxY.at(_energy), 2) ) > CutVtxR.at(_energy) ) return false;
  return true;
}// isGoodEvent(){}

//********************CHECK TRACK ON GOOD********************//
Bool_t isGoodTrackEP(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy) {

  TVector3 pVtx = event->primaryVertex();

  if ( !track ) return false;
  // Must be a primary track
  if ( !track->isPrimary() ) return false;
  if ( ( track -> dEdx() ) == 0. ) return false;
  // Simple single-track cut
  if( track -> gMom().Mag() < 0.1) return false; 
  if( track -> gDCA(pVtx).Mag() > CutDCApidEP.at(_energy) ) return false;    
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( TMath::Abs( track -> eta() ) < 0.05 ) return false; 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtotEPMin_PID) return false;
  if( track -> pt() > CutPtotEPMax_PID) return false; 
  if( ( (Double_t)track -> nHits() )/( (Double_t)track -> nHitsPoss() )  < CutnHitsRatio ) return false;

  return true; 
}// isGoodTrack(){}

//********************CHECK TRACK FLOW ON GOOD********************//
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy) {

  TVector3 pVtx = event->primaryVertex();

  if ( !track ) return false;
  // Must be a primary track
  if ( !track->isPrimary() ) return false;
  if ( ( track -> dEdx() ) == 0. ) return false;
  // Simple single-track cut
  if( track -> gMom().Mag() < 0.1) return false; 
  if( track -> gDCA(pVtx).Mag() > CutDCApidEP.at(_energy) ) return false;
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( TMath::Abs( track -> eta() ) < 0.05 ) return false; 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtotEPMin_PID) return false;
  if( track -> pt() > CutPtotFlowMax_PID) return false; 
  if( track -> p() < CutPtotEPMin_PID) return false;
  if( track -> p() > CutPtotFlowMax_PID) return false; 
  if( ( (Double_t)track -> nHits() )/( (Double_t)track -> nHitsPoss() )  < CutnHitsRatio ) return false;

  return true; 
}// isGoodTrack(){}

void GetVectorRunId(const Char_t *inFile, std::vector<int> &runIdVec){
    
    std::string s; // сюда будем класть считанные строки
    std::string s_f = "st_physics_"; // часть строки, которую ищем
    std::ifstream file(inFile); 

    while(std::getline(file, s)){
      //std::cout << s << std::endl; // выводим на экран
      //std::cout<<s.substr((s.find(s_f))+s_f.size(),8)<<"\n";
      int runId = std::stoi(s.substr((s.find(s_f))+s_f.size(),8));
      bool k=false;
      for(Int_t i=0; i<(int)runIdVec.size(); i++){
        if(runIdVec[i]==runId)k=true;
      }
      if(k==false)runIdVec.push_back(runId);
    }

    file.close(); // обязательно закрываем файл что бы не повредить его
}

int PID_TPC_TOF(StFemtoTrack *const &track, const Double_t DCA ,const Int_t _energy){

  Int_t ch = -99;
  
  if( DCA > CutDCApidFlow.at(_energy) ) return -1;    

  if ( track->isTofTrack() ){
    if(track -> charge() > 0)ch=0;
    if(track -> charge() < 0)ch=4;
    if( TMath::Abs( track->nSigmaPion() ) < nSigmaTofTpc.at(_energy) && track->massSqr() > SqMdown[0] && track->massSqr() < SqMup[0]){
      return 1 + ch;
    }
    if( TMath::Abs( track->nSigmaKaon() ) < nSigmaTofTpc.at(_energy) && track->massSqr()> SqMdown[1] && track->massSqr() < SqMup[1]){
      return 2 + ch;
    }
    if( TMath::Abs( track->nSigmaProton() ) < nSigmaTofTpc.at(_energy) && track->massSqr()> SqMdown[2] && track->massSqr() < SqMup[2]){
      return 3 + ch;
    }
  }
  return -1;
}// PID

Double_t GetRapidity2(StFemtoTrack *const &track, Int_t iPID){
  Double_t rap = -999;
  double sqrt_mass_pion   = 0.0194797849;
  double sqrt_mass_kaon   = 0.244036;
  double sqrt_mass_proton = 0.8803505929;
  if( iPID==1 || iPID==5 )
  {
   rap = 0.5*TMath::Log( (TMath::Sqrt((track->pMom())*(track->pMom())+ sqrt_mass_pion)+(track->pMom().Z()))/ (TMath::Sqrt((track->pMom())*(track->pMom())+ sqrt_mass_pion)-(track->pMom().Z()))   );
  }
  else if( iPID==2 || iPID==6 )
  {
    rap = 0.5*TMath::Log( (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_kaon)+(track->pMom().Z()))/ (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_kaon)-(track->pMom().Z()))   );
  }
  else if( iPID==3 || iPID==7 )
  {
    rap = 0.5*TMath::Log( (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_proton)+(track->pMom().Z()))/ (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_proton)-(track->pMom().Z()))   );
  }
  return rap;
}

void WriteCheckHisto2DForAllRun(std::vector<std::vector<FlowEtaSubEP*>> ClassSubEP, Int_t EW){
  std::vector<TH2D*> WriteHisto2D;
  for(Int_t ih=0; ih<(int)ClassSubEP.size(); ih++){
    for(Int_t ir=0; ir<(int)ClassSubEP[ih].size(); ir++){
      if(ir==0){
        WriteHisto2D = ClassSubEP[ih][ir]->GetCheckHistoQVector(EW);
        for(Int_t ic=0; ic<nBinCent; ic++){
          std::string NameHisto = (std::string)WriteHisto2D[ic]->GetName();
          WriteHisto2D[ic]->SetName( NameHisto.substr(0,(NameHisto.find("_Run"))).data() );
        }
      }
      if(ir>0){
        for(Int_t ic=0; ic<nBinCent; ic++){
          WriteHisto2D[ic]->Add( (ClassSubEP[ih][ir]->GetCheckHistoQVector(EW))[ic], 1.0);
        }
      }
    }
    for(Int_t i=0; i<(int)WriteHisto2D.size(); i++){
      WriteHisto2D[i]->Write();
    }
    WriteHisto2D.clear();
  }
}

void WriteCheckHisto1DForAllRun(std::vector<std::vector<FlowEtaSubEP*>> ClassSubEP, Int_t EW){
  std::vector<TH1D*> WriteHisto1D;
  for(Int_t ih=0; ih<(int)ClassSubEP.size(); ih++){
    for(Int_t ir=0; ir<(int)ClassSubEP[ih].size(); ir++){
      if(ir==0){
        WriteHisto1D = ClassSubEP[ih][ir]->GetCheckHistoPsi(EW);
        for(Int_t ic=0; ic<nBinCent; ic++){
          std::string NameHisto = (std::string)WriteHisto1D[ic]->GetName();
          WriteHisto1D[ic]->SetName( NameHisto.substr(0,(NameHisto.find("_Run"))).data() );
        }
      }
      if(ir>0){
        for(Int_t ic=0; ic<nBinCent; ic++){
          WriteHisto1D[ic]->Add( (ClassSubEP[ih][ir]->GetCheckHistoPsi(EW))[ic], 1.0);
        }
      }
    }
    for(Int_t i=0; i<(int)WriteHisto1D.size(); i++){
      WriteHisto1D[i]->Write();
    }
    WriteHisto1D.clear();
  }
}