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
#include "TProfile.h"
#include "TProfile2D.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"
#include "TF2.h"
#include "TF1.h"
#include "TGraphErrors.h"


// FemtoDst headers

#include "StFemtoDstReader.h"
#include "StFemtoDst.h"
#include "StFemtoEvent.h"
#include "StFemtoTrack.h"
#include "StFemtoV0.h"
#include "StFemtoXi.h"

// Constant
#include "Constants.h"
#include "functions.C"

// Load libraries (for ROOT_VERSTION_CODE >= 393215)
#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(libStFemtoDst.so)
#endif

// inFile - is a name of name.FemtoDst.root file or a name
//          of a name.lis(t) files, that contains a list of
//          name1.FemtoDst.root, name2.FemtoDst.root, ... files

// Used function
Bool_t isGoodEvent(StFemtoEvent *const &event, const Int_t _energy, const Bool_t BadRunIdKeyFlowStage);
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy);

TVector2 CalcNewXY(StFemtoTrack *const &track, Double_t width_nsigma, Double_t width_m2, 
                                               Double_t mean_nsigma_Kaon, Double_t mean_m2_Kaon,
                                               Double_t mean_nsigma_Pion, Double_t mean_m2_Pion);


//_________________
void PIDcomb(const Char_t *inFile = "st_physics_12150008_raw_4030001.femtoDst.root",
                          const Char_t *outFileName = "oTest.root",
                          const Char_t *mode = "QAmode",
                          const Int_t energy = 39) {

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
  std::cout << "Status has been set" << std::endl;

  std::cout << "Now I know what to read, Master!" << std::endl;

  Int_t RunIdRange = RunIdMax.at(energy) - RunIdMin.at(energy);
  
  Int_t cent;
  Int_t RunID;

  // Mode works
  Bool_t check_TofMatch = false;
  Bool_t mode_first = false;
  Bool_t mode_NewXY = false;
  
  if( strncmp(mode, "first",5) == 0){
    mode_first = true;
  } 
  if( strncmp(mode, "NewXY",6) == 0){
    mode_NewXY = true;
  }
  Double_t factor=2.;
  if(energy==27 || energy==62){
    factor=1.;
  }
  //create file 
  TFile *outFile = new TFile(outFileName, "RECREATE");
  TFile *FileRec, *FileFlow;
  TFile *FileFitM2;
  outFile->cd();

  TH2D *h2_m2VsPt_all[2];
  TH1D *h1_m2_ptBin[2][30];
  
  TH2D *h2_m2VsPt_all_new[2][3];
  TH2D *h2_m2VsPt_all_new_fit[2][3];

  //___________
  TGraphErrors *gr_ReadGraph;
  TGraphErrors *gr_MeanNSigma[2][3];
  TGraphErrors *gr_SigmaNSigma[2][3];
  TGraphErrors *gr_MeanM2[2][3];
  TGraphErrors *gr_SigmaM2[2][3];

  TH2D *h2_m2vsnSigmaPion_piKp[2][30];
  TH2D *h2_m2vsnSigmaPion_piKp_nSigmaM2[2][30];
  
  TH2D *h2_m2vsnSigmaPion_piK_new[2][30];
  TH2D *h2_m2vsnSigmaPion_piK_new_cut1[2][30];
  TH2D *h2_m2vsnSigmaPion_piK_new_cut2[2][30];
  TH2D *h2_m2vsnSigmaPion_piK_new_cut3[2][30];
  TH2D *h2_m2vsnSigmaPion_piK_new_cut4[2][30];

  TH2D *h2_m2vsnSigmaPion_piK[2][30];
  
  if(mode_first==true){
    
    h2_m2VsPt_all[0] = new TH2D("h2_m2VsPt_all_ch0","m^2 vs p_{T} (charge > 0); p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",1000,0,5,2000,-0.5,1.5);
    h2_m2VsPt_all[1] = new TH2D("h2_m2VsPt_all_ch1","m^2 vs p_{T} (charge < 0); p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",1000,0,5,2000,-0.5,1.5);

    for(Int_t ch=0; ch<2; ch++){
      for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
        h2_m2vsnSigmaPion_piKp[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piKp_%i_pt%i",ch,pti),Form("m^{2} vs n#sigma(%s) %.1f<p_T<%.1f GeV/c;n#sigma(%s);m^{2},(GeV/c^{2})^{2}", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch]), factor*300, -6*factor, 6*factor, 300, -1.0, 2.0 );
      }
    }

  }//mode_nSigma

  if(mode_NewXY==true){
    
    h2_m2VsPt_all[0] = new TH2D("h2_m2VsPt_all_ch0","m^2 vs p_{T} (charge > 0); p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",1000,0,5,400,-0.5,1.5);
    h2_m2VsPt_all[1] = new TH2D("h2_m2VsPt_all_ch1","m^2 vs p_{T} (charge < 0); p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}",1000,0,5,400,-0.5,1.5);

    for(Int_t ch=0; ch<2; ch++){
      for(Int_t pti=0; pti<(int)ptBinRange.size()-1; pti++){
        //h2_m2vsnSigmaPion_piKp[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piKp_%i_pt%i",ch,pti),Form("m^{2} vs n#sigma(%s) %.1f<p_T<%.1f GeV/c;n#sigma(%s);m^{2},(GeV/c^{2})^{2}", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch]), factor*600, -6*factor, 6*factor, factor*600, -1.0, 2.0 );
     
        h2_m2vsnSigmaPion_piK_new[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piK_%i_pt%i_new",ch,pti),Form("m^{2} vs n#sigma(%s) New coord %.1f<p_T<%.1f GeV/c;x(n#sigma(%s),m^{2});y(n#sigma(%s),m^{2})", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch],partLateX[ch]), 800, -2, 2, 800, -2.0, 2.0 );
        h2_m2vsnSigmaPion_piK_new_cut1[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piK_%i_pt%i_new_cut1",ch,pti),Form("m^{2} vs n#sigma(%s) New coord %.1f<p_T<%.1f GeV/c;x(n#sigma(%s),m^{2});y(n#sigma(%s),m^{2})", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch],partLateX[ch]), 800, -2, 2, 800, -2.0, 2.0 );
        //h2_m2vsnSigmaPion_piK_new_cut2[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piK_%i_pt%i_new_cut2",ch,pti),Form("m^{2} vs n#sigma(%s) New coord %.1f<p_T<%.1f GeV/c;x(n#sigma(%s),m^{2});y(n#sigma(%s),m^{2})", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch],partLateX[ch]), 800, -2, 2, 800, -2.0, 2.0 );
        //h2_m2vsnSigmaPion_piK_new_cut3[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piK_%i_pt%i_new_cut3",ch,pti),Form("m^{2} vs n#sigma(%s) New coord %.1f<p_T<%.1f GeV/c;x(n#sigma(%s),m^{2});y(n#sigma(%s),m^{2})", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch],partLateX[ch]), 800, -2, 2, 800, -2.0, 2.0 );
        h2_m2vsnSigmaPion_piK_new_cut4[ch][pti] = new TH2D(Form("h2_m2VsnSigma_piK_%i_pt%i_new_cut4",ch,pti),Form("m^{2} vs n#sigma(%s) New coord %.1f<p_T<%.1f GeV/c;x(n#sigma(%s),m^{2});y(n#sigma(%s),m^{2})", partLateX[ch], ptBinRange[pti], ptBinRange[pti+1] ,partLateX[ch],partLateX[ch]), 800, -2, 2, 800, -2.0, 2.0 );

      }
    }

    /*
    for(Int_t par=0; par<3; par++){
      h2_m2VsPt_all_new_fit[0][par] = new TH2D(Form("h2_m2VsPt_all_%s_ch0_new_fit",particles[par]),"n#sigma(m^{2}) distributionm vs p_{T} (charge > 0); p_{T} (GeV/c);n#sigma(m^{2})",640,0,3.2,1500,-15,15);
      h2_m2VsPt_all_new_fit[1][par] = new TH2D(Form("h2_m2VsPt_all_%s_ch1_new_fit",particles[par]),"n#sigma(m^{2}) distributionm vs p_{T} (charge < 0); p_{T} (GeV/c);n#sigma(m^{2})",640,0,3.2,1500,-15,15);
    }
    */

    FileFitM2 = new TFile(Form("%s/OUT/%iGeV/FitFun_%iGeVRun10.root",path,energy,energy),"READ");
    FileFitM2->cd();

    for(Int_t ch=0; ch<2; ch++){
      for(Int_t par=0; par<3; par++){
        
        gr_ReadGraph = (TGraphErrors*)FileFitM2->Get(Form("gr_%s_meanNSigma_%s_ch%i","old",particles[par],ch));
        gr_MeanNSigma[ch][par] = (TGraphErrors*)gr_ReadGraph->Clone();
        delete gr_ReadGraph;

        gr_ReadGraph = (TGraphErrors*)FileFitM2->Get(Form("gr_%s_sigmaNSigma_%s_ch%i","old",particles[par],ch));
        gr_SigmaNSigma[ch][par] = (TGraphErrors*)gr_ReadGraph->Clone();
        delete gr_ReadGraph;

        gr_ReadGraph = (TGraphErrors*)FileFitM2->Get(Form("gr_%s_meanM2_%s_ch%i","old",particles[par],ch));
        gr_MeanM2[ch][par] = (TGraphErrors*)gr_ReadGraph->Clone();
        delete gr_ReadGraph;

        gr_ReadGraph = (TGraphErrors*)FileFitM2->Get(Form("gr_%s_sigmaM2_%s_ch%i","old",particles[par],ch));
        gr_SigmaM2[ch][par] = (TGraphErrors*)gr_ReadGraph->Clone();
        delete gr_ReadGraph;

      }
    }
    
    FileFitM2->Close();

  }//mode_nSigma

  outFile->cd();

  if( !femtoReader->chain() ) {
    std::cout << "No chain has been found." << std::endl;
  }
  Long64_t eventsInTree = femtoReader->tree()->GetEntries();
  std::cout << "eventsInTree: "  << eventsInTree << std::endl;
  Long64_t events2read = femtoReader->chain()->GetEntries();

  std::cout << "Number of events to read: " << events2read << std::endl;

  // Loop over events
  for(Long64_t iEvent=0; iEvent<events2read; iEvent++) {

    if ( iEvent % 10000 == 0) {
      std::cout << "Working on event #[" << (iEvent+1)
                << "/" << events2read << "]" << std::endl;
    }

    if (iEvent == events2read-1) {
      std::cout << "Working on event #[" << (events2read)
                << "/" << events2read << "]" << std::endl;
    }

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
    
    Int_t nTracks = dst->numberOfTracks();

    //Work place

    //Event cut
    if( !isGoodEvent( event,  energy, true) ) continue;
    if( !TofMatchedCut(dst, event, 4, energy) ) continue;
    
    // Event cut

    for(Int_t iTrk=0; iTrk<nTracks; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);

      if( !femtoTrack ) continue;
      if( !isGoodTrackFlow(event, femtoTrack, energy) ) continue;
      if( !femtoTrack->isTofTrack() ) continue;

      Int_t charge = GetCharge(femtoTrack);
      Int_t ptBin = GetBinPtRange(femtoTrack);
      Double_t pt = femtoTrack->pt();

      if( charge == -1 || ptBin == -1 ) continue;

      //h2_m2vsnSigmaPion_piKp[charge][ptBin]->Fill(femtoTrack->nSigmaPion(), femtoTrack->massSqr());

      if(mode_first==true){
        h2_m2VsPt_all[charge]->Fill(femtoTrack -> pt(), femtoTrack->massSqr());
        h2_m2vsnSigmaPion_piKp[charge][ptBin]->Fill(femtoTrack->nSigmaPion(), femtoTrack->massSqr());

      }//mode_nSigma

      if(mode_NewXY==true){
        h2_m2VsPt_all[charge]->Fill(femtoTrack -> pt(), femtoTrack->massSqr());

        //h2_m2VsPt_all_new_fit[charge][0]->Fill(femtoTrack -> pt(), (femtoTrack->massSqr() - gr_MeanM2[charge][0]->Eval(pt,0,"S")) / gr_SigmaM2[charge][0]->Eval(pt,0,"S"));
        //h2_m2VsPt_all_new_fit[charge][1]->Fill(femtoTrack -> pt(), (femtoTrack->massSqr() - gr_MeanM2[charge][1]->Eval(pt,0,"S")) / gr_SigmaM2[charge][1]->Eval(pt,0,"S"));
        //h2_m2VsPt_all_new_fit[charge][2]->Fill(femtoTrack -> pt(), (femtoTrack->massSqr() - gr_MeanM2[charge][2]->Eval(pt,0,"S")) / gr_SigmaM2[charge][2]->Eval(pt,0,"S"));

        Double_t GetCenterPtBin = (ptBinRange[ptBin] + ptBinRange[ptBin+1])/2.0;

        Double_t NewX = CalcNewXY(femtoTrack, gr_SigmaNSigma[charge][0]->Eval( GetCenterPtBin ,0,"S"), 
                                                  gr_SigmaM2[charge][0]->Eval( GetCenterPtBin ,0,"S"), 
                                               gr_MeanNSigma[charge][1]->Eval( GetCenterPtBin ,0,"S"),
                                                   gr_MeanM2[charge][1]->Eval( GetCenterPtBin ,0,"S"), 
                                               gr_MeanNSigma[charge][0]->Eval( GetCenterPtBin ,0,"S"),
                                                   gr_MeanM2[charge][0]->Eval( GetCenterPtBin ,0,"S") ).X();
        Double_t NewY = CalcNewXY(femtoTrack, gr_SigmaNSigma[charge][0]->Eval( GetCenterPtBin ,0,"S"), 
                                                  gr_SigmaM2[charge][0]->Eval( GetCenterPtBin ,0,"S"), 
                                               gr_MeanNSigma[charge][1]->Eval( GetCenterPtBin ,0,"S"),
                                                   gr_MeanM2[charge][1]->Eval( GetCenterPtBin ,0,"S"), 
                                               gr_MeanNSigma[charge][0]->Eval( GetCenterPtBin ,0,"S"),
                                                   gr_MeanM2[charge][0]->Eval( GetCenterPtBin ,0,"S") ).Y();

        //h2_m2vsnSigmaPion_piKp[charge][ptBin]->Fill(femtoTrack->nSigmaPion(), femtoTrack->massSqr());
        h2_m2vsnSigmaPion_piK_new[charge][ptBin]->Fill(NewX, NewY);

        if( femtoTrack->massSqr() < 0.4 ){
          h2_m2vsnSigmaPion_piK_new_cut1[charge][ptBin]->Fill(NewX, NewY);
        }
        /*
        if( TMath::Abs( (femtoTrack->massSqr() - gr_MeanM2[charge][2]->Eval(pt,0,"S")) / gr_SigmaM2[charge][2]->Eval(pt,0,"S") ) > 3.0 ){
          h2_m2vsnSigmaPion_piK_new_cut2[charge][ptBin]->Fill(NewX, NewY);
        }

        if( TMath::Abs( (femtoTrack->massSqr() - gr_MeanM2[charge][2]->Eval(pt,0,"S")) / gr_SigmaM2[charge][2]->Eval(pt,0,"S") ) > 3.0 && femtoTrack->massSqr() < 0.4){
          h2_m2vsnSigmaPion_piK_new_cut3[charge][ptBin]->Fill(NewX, NewY);
        }
        */
        if( femtoTrack->massSqr() < ( gr_MeanM2[charge][2]->Eval(pt,0,"S") - 3.0 * gr_SigmaM2[charge][2]->Eval(pt,0,"S"))  ){
          h2_m2vsnSigmaPion_piK_new_cut4[charge][ptBin]->Fill(NewX, NewY);
        }

      }// mode_NewXY;

    } //for(Int_t iTrk=0; iTrk<nTracks; iTrk++)
      

  } //for(Long64_t iEvent=0; iEvent<events2read; iEvent++)

  outFile->Write();
  outFile->Close();

  femtoReader->Finish();

  std::cout << "I'm done with analysis. We'll have a Nobel Prize, Master!"
	          << std::endl;
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
  
  std::vector<Int_t> Bad = BadRuns.at(_energy);
  if(BadRunIdKeyFlowStage==true){
    if( std::find(Bad.begin(), Bad.end(), event -> runId()) != Bad.end() ) return false;
  }
  // if( event->numberOfTofMatched() <= 4) return false;
/*  
std::vector<Int_t> trig = Trigger.at(_energy);
  for(Int_t i = 0; i < trig.size(); i++){
    if( event->isTrigger(trig[i]) ) return true;
  }
  */
  return true;
}// isGoodEvent(){}


//********************CHECK TRACK FLOW ON GOOD********************//
Bool_t isGoodTrackFlow(StFemtoEvent *const &event, StFemtoTrack *const &track, const Int_t _energy) {

  TVector3 pVtx = event->primaryVertex();

  if ( !track ) return false;
  // Must be a primary track
  if ( !track->isPrimary() ) return false;
  if ( ( track -> dEdx() ) == 0. ) return false;
  // Simple single-track cut
  if( track -> gMom().Mag() < 0.1) return false; 
  if( track -> gDCA(pVtx).Mag() > CutDCApidFlow.at(_energy) ) return false;
  if( TMath::Abs( track -> eta() ) > 1.0 ) return false; 
  if( TMath::Abs( track -> eta() ) < 0.05 ) return false; // for pid 
  if( track -> nHits() < CutnHits) return false;
  if( track -> pt() < CutPtotFlowMin_PID) return false;
  if( track -> pt() > CutPtotFlowMax_PID) return false; 
  if( ( (Double_t)track -> nHits() )/( (Double_t)track -> nHitsPoss() )  < CutnHitsRatio ) return false;

  return true; 
}// isGoodTrack(){}

TVector2 CalcNewXY(StFemtoTrack *const &track, Double_t width_nsigma, Double_t width_m2, 
                                               Double_t mean_nsigma_Kaon, Double_t mean_m2_Kaon,
                                               Double_t mean_nsigma_Pion, Double_t mean_m2_Pion){

  TVector2 qv(0.,0.);

  Double_t f_scale = width_nsigma / width_m2;

  Double_t den = (mean_nsigma_Kaon - mean_nsigma_Pion)/f_scale;
  
  Double_t alpha = -1*TMath::ATan2(mean_m2_Kaon - mean_m2_Pion , den); 
  
  Double_t x1 = ( track->nSigmaPion() - mean_nsigma_Pion ) / f_scale;
  Double_t y1 = track->massSqr() - mean_m2_Pion;

  qv.Set(TMath::Cos(alpha)*x1 - TMath::Sin(alpha)*y1 , TMath::Sin(alpha)*x1 + TMath::Cos(alpha)*y1 );

  return qv;
}