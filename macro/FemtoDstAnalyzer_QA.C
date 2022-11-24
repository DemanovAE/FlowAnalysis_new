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
#include "StEpdGeom.h"
#include "QAFemtoDst2.h"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStFemtoDst.so)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

// Used function
void GetVectorRunId(const Char_t *inFile, std::vector<int> &runIdVec);
Double_t GetRapidity(StFemtoTrack *const &track);

//_________________
void FemtoDstAnalyzer_QA(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                          const Char_t *outFileName = "oTest.root",
                          const Char_t *mode = "QAmode",
                          const Int_t energy = 39) {

  TStopwatch timer1;
  timer1.Start();
  
  std::cout << "Hi! Lets do some physics, Master!" << std::endl;

  #if ROOT_VERSION_CODE < ROOT_VERSION(6,0,0)
    gSystem->Load("libStFemtoDst.so");
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
  outFile->cd();
  // Mode works
  Bool_t CalEPFromTPC = true;
  Bool_t CalEPFromEPD = false;

  GetVectorRunId(inFile,vRunId);
  for(Int_t i=0; i<(int)vRunId.size();i++){
    mRunId.insert(std::make_pair( vRunId[i], i));
    std::cout<<vRunId[i]<<"\t map"<<mRunId.at(vRunId[i])<<"\n";
  }
  
  StEpdGeom *EpdGeom = new StEpdGeom(); 

  Int_t NumberCut = 6;
  std::vector<TString>  NameSuf = {"Raw","WithCut","AddVtxR","AddTofMatched","AddEta","AddNHits","AddNHitsRatio","AddDCA","AddPtot","AddPt","AddDCApid"};
  std::vector<QAFemtoDst2*> QAFemtoDstObj;
  for(Int_t i=0; i<NumberCut; i++){
    QAFemtoDstObj.push_back(new QAFemtoDst2());
    QAFemtoDstObj[i]->SetRunIdRange(RunIdMin.at(energy),RunIdMax.at(energy));
    if(i>0)QAFemtoDstObj[i]->SetEventCutVtxZ(CutVtxZ.at(energy));    
    if(i>0)QAFemtoDstObj[i]->SetEventCutVtxR(CutVtxR.at(energy));    
    if(i>1)QAFemtoDstObj[i]->SetTrackCutPseudorapidity(1.0);
    if(i>1)QAFemtoDstObj[i]->SetTrackCutNHits(15);
    if(i>1)QAFemtoDstObj[i]->SetTrackCutNHitsRatioNHitsPoss(0.52);
    if(i>1)QAFemtoDstObj[i]->SetTrackCutDCA(2.0);
    if(i>1)QAFemtoDstObj[i]->SetTrackCutPtotMin(0.1);
    if(i>1)QAFemtoDstObj[i]->SetTrackCutPtMin(0.2);
    if(i==0)QAFemtoDstObj[i]->SetNameSuf("Raw");
    if(i==1)QAFemtoDstObj[i]->SetNameSuf("RawTr");
    if(i==2)QAFemtoDstObj[i]->SetNameSuf("Charge_hadrons");
    if(i==3){
      QAFemtoDstObj[i]->SetNameSuf("Pion");   
      QAFemtoDstObj[i]->SetTrackCutPIDcode(1);
      QAFemtoDstObj[i]->SetTrackCutDCA(1.0);
    }
    if(i==4){
      QAFemtoDstObj[i]->SetNameSuf("Kaon");
      QAFemtoDstObj[i]->SetTrackCutPIDcode(2);
      QAFemtoDstObj[i]->SetTrackCutDCA(1.0);
    }
    if(i==5){
      QAFemtoDstObj[i]->SetNameSuf("Proton"); 
      QAFemtoDstObj[i]->SetTrackCutPIDcode(3);
      QAFemtoDstObj[i]->SetTrackCutDCA(1.0);
    }
    QAFemtoDstObj[i]->SetTrackCutNSigmaTPC({3.,3.,3.,3.});
    QAFemtoDstObj[i]->SetTrackCutM2({0.00000025, -0.15, 0.2, 0.74},{0.00000027, 0.1, 0.32, 1.20});
    QAFemtoDstObj[i]->Init();
  }
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

    // Track analysis
    Int_t nTracks = dst->numberOfTracks();
  
    /////////////// BBS ADC SUM ///////////////////////
    Float_t bbcAdcSum = 0;
    for( Int_t iTile=0; iTile<24; iTile++ ) {
      bbcAdcSum += event->bbcAdcEast(iTile);
      bbcAdcSum += event->bbcAdcWest(iTile);
    }
    ///////////// TOF MAtched ////////////////////////
    Int_t numberTofMatched = 0;
    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      
      if (!femtoTrack) continue;
      if (!femtoTrack->isPrimary() ) continue;
      if (!femtoTrack->isTofTrack() ) continue;
      numberTofMatched++;
    }
    /////////////////////////////////////////////////

    for(Int_t i=0; i<NumberCut; i++){
      QAFemtoDstObj[i]->FillEventHisto(
        event->runId(), 
        event->cent9(),
        event->cent16(),
        event->primaryVertex(),
        event->ranking(),
        numberTofMatched,
        event->refMult(),
        event->refMult2(),
        event->gRefMult(),
        event->numberOfPrimaryTracks(),
        event->numberOfGlobalTracks(),
        event->numberOfBTofHit(),
        event->numberOfTofMatched(),
        numberTofMatched,
        event->vpdVz(),
        event->zdcSumAdcEast(),
        event->zdcSumAdcWest(),
        event->transverseSphericity(),
        event->transverseSphericity2(),
        event->numberOfPrimaryVertices(),
        bbcAdcSum );
    }

    if(CalEPFromTPC){
      // Track loop
      for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
        // Retrieve i-th femto track
        StFemtoTrack *femtoTrack = dst->track(iTrk);

        for(Int_t i=0; i<NumberCut; i++){
          QAFemtoDstObj[i]->FillTrackHisto(
            event->runId(), 
            event->cent9(),
            event->primaryVertex(),
            numberTofMatched,
            femtoTrack->nHits(),
            femtoTrack->nHitsPoss(),
            GetRapidity(femtoTrack),
            femtoTrack->gDCA(event->primaryVertex()),
            femtoTrack->gMom(),
            femtoTrack->pMom(),
            femtoTrack->chi2(),
            femtoTrack->charge(),
            femtoTrack->dEdx(),
            femtoTrack->nSigmaElectron(),
            femtoTrack->nSigmaPion(),
            femtoTrack->nSigmaKaon(),
            femtoTrack->nSigmaProton(),
            femtoTrack->isTofTrack(),
            femtoTrack->beta(),
            femtoTrack->invBeta(),
            femtoTrack->massSqr() );
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
        double RandomEta = StraightLine.Eta();
        double RandomPhi = StraightLine.Phi(); while(RandomPhi < 0.) RandomPhi += 2.*pi;
        int ew = (femtoEpdHit->id()<0)?0:1;  //is EPD east or west
        if (!(TMath::Abs(RandomEta)>0)||(TMath::Abs(RandomEta)>1000)) continue;
        //std::cout<<"good\n";
        //std::cout<<RandomEta<<"\t"<<RandomPhi<<"\t"<<nMip<<std::endl;
      }

    }// EPD EP Determ

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  outFile->cd();
  outFile->Write();
  outFile->Close();

  femtoReader->Finish();

  for(Int_t ir=0; ir<NumberCut; ir++){
    delete QAFemtoDstObj[ir];
  }
  QAFemtoDstObj.clear();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;

  timer1.Stop();
  timer1.Print();
}

  /*////////////////////////////////////////////////////////////////////////////////////////*/
 /*___________________________DESCRIPTION OF FUNCTIONS_____________________________________*/
/*////////////////////////////////////////////////////////////////////////////////////////*/

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

Double_t GetRapidity(StFemtoTrack *const &track){
  Double_t rap = -999;
  double sqrt_mass_pion   = 0.0194797849;
  double sqrt_mass_kaon   = 0.244036;
  double sqrt_mass_proton = 0.8803505929;
  if( track->massSqr()>=-0.15 && track->massSqr()<=0.10 && abs(track->nSigmaPion())<=3.0 )
  {
   rap = 0.5*TMath::Log( (TMath::Sqrt((track->pMom())*(track->pMom())+ sqrt_mass_pion)+(track->pMom().Z()))/ (TMath::Sqrt((track->pMom())*(track->pMom())+ sqrt_mass_pion)-(track->pMom().Z()))   );
  }
  else if( (track->massSqr())>0.2 && (track->massSqr())<0.32 && abs((track->nSigmaKaon()))<3)
  {
    rap = 0.5*TMath::Log( (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_kaon)+(track->pMom().Z()))/ (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_kaon)-(track->pMom().Z()))   );
  }
  else if( (track->massSqr())>0.74 && (track->massSqr())<1.2 && abs((track->nSigmaProton()))<3)
  {
    rap = 0.5*TMath::Log( (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_proton)+(track->pMom().Z()))/ (TMath::Sqrt((track->pMom())*(track->pMom())+sqrt_mass_proton)-(track->pMom().Z()))   );
  }
  return rap;
}