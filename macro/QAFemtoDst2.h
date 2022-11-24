#ifndef QAFEMTODST2_H
#define QAFEMTODST2_H
#include <iostream>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TDirectoryFile.h>

class QAFemtoDst2
{
public:
  QAFemtoDst2();
  virtual ~QAFemtoDst2();

  void Init();
  void SetEnergy(Int_t i){ this->fEnergy = i; }
  void SetRunIdRange(Int_t imin, Int_t imax){ this->fRunIdMin=imin; this->fRunIdMax=imax; }
  void SetNameSuf(TString str) { this->fNameSuf = str; }
  void SetEventCutVtxZ(Double_t d){ this->fEventCutVtxZ = d; }
  void SetEventCutVtxR(Double_t d){ this->fEventCutVtxR = d; }
  void SetEventCutTofMatched(Double_t d){ this->fEventCutTofMatched = d; }
  void SetTrackCutNHits(Double_t d){ this->fTrackCutNHits = d; }
  void SetTrackCutNHitsRatioNHitsPoss(Double_t d){ this->fTrackCutNHitsRatioNHitsPoss = d; }
  void SetTrackCutDCA(Double_t d){ this->fTrackCutDCA = d; }
  void SetTrackCutPseudorapidity(Double_t d){ this->fTrackCutPseudorapidity = d; }
  void SetTrackCutPtotMin(Double_t d){ this->fTrackCutPtotMin = d; }
  void SetTrackCutPtMin(Double_t d){ this->fTrackCutPtMin = d; }
  void SetTrackCutIsPrimary(Bool_t bl){ this->fTrackCutIsPrimary = bl; }
  void SetTrackCutNSigmaTPC(std::vector<Double_t> vec){ this->fTrackCutNSigmaTPC = vec; }
  void SetTrackCutM2(std::vector<Double_t> vecMin, std::vector<Double_t> vecMax){ this->fTrackCutM2Min = vecMin; this->fTrackCutM2Max = vecMax; }
  void SetTrackCutPIDcode(Int_t i){this->fPIDcode = i;}
  
  void PidTpcAndTof(Double_t nSigma, Double_t m2); 

  void FillEventHisto(
    const Int_t run, 
    const Int_t cent9,
    const Int_t cent16,
    const TVector3 pVtx,
    const Float_t ranking,
    const Int_t TofMatched,
    const Int_t refMult,
    const Int_t refMult2,
    const Int_t gRefMult,
    const Int_t numberOfPrimaryTracks,
    const Int_t numberOfGlobalTracks,
    const Int_t numberOfBTofHit,
    const Int_t numberOfTofMatched,
    const Int_t numberOfTofMatchedOld,
    const Float_t vpdVz,
    const Float_t zdcSumAdcEast,
    const Float_t zdcSumAdcWest,
    const Float_t transverseSphericity,
    const Float_t transverseSphericity2,
    const UInt_t numberOfPrimaryVertices,
    const Float_t bbcAdcSum);

  void FillTrackHisto(
    const Int_t run, 
    const Int_t cent,
    const TVector3 pVtx,
    const Int_t TofMatched,
    const Int_t nHits,
    const Int_t nHitsPoss,
    const Double_t Rapidity,
    const TVector3 gDCA,
    const TVector3 gMom,
    const TVector3 pMom,
    const Float_t chi2,
    const Int_t charge,
    const Double_t dEdx,
    const Double_t nSigmaElectron,
    const Double_t nSigmaPion,
    const Double_t nSigmaKaon,
    const Double_t nSigmaProton,
    const Bool_t isTofTrack,
    const Double_t beta,
    const Double_t invBeta,
    const Double_t massSqr);

  Int_t GetEnergy(){ return this->fEnergy; }

private:

  const Float_t electron_mass_sqr = 0.000000301;
  const Float_t pion_mass_sqr = 0.019479955;
  const Float_t kaon_mass_sqr = 0.24371698;
  const Float_t proton_mass_sqr = 0.880354499;

  Int_t fEnergy;
  TString fNameSuf;
  Double_t fRunIdMax;
  Double_t fRunIdMin;
  Double_t fEventCutVtxZ;
  Double_t fEventCutVtxR;
  Double_t fEventCutTofMatched;
  Int_t fPIDcode; //0-pi,1-k,2-p;
  Bool_t fTrackCutIsPrimary;
  Double_t fTrackCutNHits;
  Double_t fTrackCutNHitsRatioNHitsPoss;
  Double_t fTrackCutDCA;
  Double_t fTrackCutPseudorapidity;
  Double_t fTrackCutPtotMin;
  Double_t fTrackCutPtMin;
  std::vector<Double_t> fTrackCutNSigmaTPC;
  std::vector<Double_t> fTrackCutM2Min;
  std::vector<Double_t> fTrackCutM2Max;

  // Events
  TH1D *QA_hRefMult;
  TH1D *QA_hVtxZ;
  TH1D *QA_hVtxX;
  TH1D *QA_hVtxY;
  TH1D *QA_hRefMult2;
  TH1D *QA_hGRefMult;
  TH1D *QA_hNumberOfPrimaries;
  TH1D *QA_hNumberOfGlobals;
  TH1D *QA_hCent9;
  TH1D *QA_hCent16;
  TH1D *QA_hBTofHit[2];
  TH1D *QA_hBTofMatched[2];
  TH1D *QA_hBemcMatched;
  TH1D *QA_hRanking;

  TH1D *QA_hTransSphericity;
  TH1D *QA_hTransSphericity2;
  TH1D *QA_hNumberOfVertices;

  TH2D *QA_hVtxXvsY;
  TH2D *QA_hVpdVzDiffVsVz;
  TH2D *QA_hZdcAdcWestEast;
  TH2D *QA_hZdcXVsnBTofMatch;
  TH2D *QA_hBTofTrayMultVsRefMult;
  TH2D *QA_hBEMCMatchVsRefMult;
  TH2D *QA_hBTofMatchedVsRefMult[2];

  TProfile *QA_hEventProfile[20];

  // Track
  TH1D *QA_hGlobalPtot;
  TH1D *QA_hPrimaryPtot;
  TH1D *QA_hGlobalPt;
  TH1D *QA_hPrimaryPt;
  TH1D *QA_hNHits;
  TH1D *QA_hNHitsRatio;
  TH1D *QA_hChi2;
  TH1D *QA_hDca;
  TH1D *QA_hDcaZ;
  TH1D *QA_hDcaX;
  TH1D *QA_hDcaY;
  TH2D *QA_hDcaVsPt;
  TH1D *QA_hPhi;
  TH1D *QA_hRapidity;
  TH2D *QA_hRapidityPt;
  TH1D *QA_hEta;
  TH1D *QA_hEtaG;
  TH2D *QA_hPtVsEta;
  TH2D *QA_hPhiVsEta;
  TH2D *QA_hPrimaryPhiVsPt[2];
  TH1D *QA_hDedx;
  TH2D *QA_hDedxVsPt;
  TH2D *QA_hNSigmaPionVsPt;
  TH2D *QA_hNSigmaElectronVsPt;
  TH2D *QA_hNSigmaKaonVsPt;
  TH2D *QA_hNSigmaProtonVsPt;

  // TofPidTrait
  TH1D *QA_hTofBeta;
  TH2D *QA_hInvBetaVsPt;
  TH1D *QA_hMassSqr;
  TH2D *QA_hMassSqrVsPt;
  TH2D *QA_hDedxVsMassSqr[2];

  TH2D *QA_hInvBetaDiffElectronVsPt;
  TH2D *QA_hInvBetaDiffPionVsPt;
  TH2D *QA_hInvBetaDiffKaonVsPt;
  TH2D *QA_hInvBetaDiffProtonVsPt;
  TH2D *QA_hDedxVsPtPID[4];

  TProfile *QA_hTrackProfile[20];
  ClassDef(QAFemtoDst2,0);
};

#endif