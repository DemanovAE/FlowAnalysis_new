#ifndef FlowScalarProduct_H
#define FlowScalarProduct_H
#include <iostream>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include <TString.h>
#include <TDirectoryFile.h>

#include "QVector.h"

class FlowScalarProduct
{
public:
  FlowScalarProduct();
  virtual ~FlowScalarProduct();
  
  // Set
  void SetMeanQVector();
  void SetMeanSinCosPsi();
  void SetFirstRun(Bool_t bl) { this->fFirstRun = bl; }
  void SetCalFlow(Bool_t bl) { this->fCalFlow = bl; }
  void SetChekQVectorAndPsi(Bool_t bl) { this->fChekHisto = bl; }
  void SetHarmonic(Int_t i) { this->fNHarmonic = i; }
  void SetHarmonicNK(Int_t n, Int_t k);
  void SetNumbFlattening(Int_t i) { this->fNumbFlattening = i; }
  void SetEtaGap(Double_t d) { this->fEtaGap = d; }
  void SetRunId(Int_t i) { this->fRunId = i; }
  void SetNCentBins(Int_t i) { this->fNCentBins = i; }
  void SetNBinsVtxZ(Int_t i) { this->fNBinsVtxZ = i; }
  void SetVtxZ(Double_t d) { this->fVtxZ = d; }
  void SetInputFileFromFirstRun(TString str) { this->fstrInputFileFromFirstRun = str; }
  void SetInputFileFromSecondRun(TString str) { this->fstrInputFileFromSecondRun = str; }
  void SetNameSys(TString str) { this->fNameSys = str;}
  void SetMapPidCode(std::map<Int_t,TString> mp){ this->fPidCode = mp; }

  // Analisys
  void Zero();
  void Init();
  void ProcessFirstTrackLoop(const Double_t &eta, const Double_t &phi, const Double_t &pt);
  void ProcessSecondTrackLoop(const Double_t &rapidity, const Double_t &eta, const Double_t &phi, const Double_t &pt, const Int_t &pid, const Double_t &dCent);
  void ProcessEventAfterFirstTrackLoop(const Double_t &dCent, const Double_t &dVtxZ);
  void WtiteHist();
  
  // Get
  int GetHarmonic() const { return this->fNHarmonic; }
  TProfile* GetTpSqRes(){ return this->tp_SqRes; };
  std::vector<TProfile2D*> GetTp2MeanPt() { return this->tp2_MeanPt; }
  std::vector<TProfile*> GetTpVnCent() { return this->tp_VnCent; }
  std::vector<TProfile2D*> GetTp2VnPtCent() { return this->tp2_VnPtCent; }
  std::vector<TString> GetNameInFile(){ return this->fNameInFile; }

private:
  std::vector<TString> fNameInFile;
  TString fstrInputFileFromFirstRun;
  TString fstrInputFileFromSecondRun;
  Bool_t fFirstRun;
  Bool_t fCalFlow;
  Bool_t fMeanQVectorForReсentrening;
  Bool_t fMeanSinCosPsiForFlattening;
  Bool_t fChekHisto;
  Bool_t fKHarmonicThroughNHarmonic;
  TString fNameSys;
  Int_t fRunId;
  Int_t fNHarmonic;
  Int_t fKHarmonic;
  Int_t fNBinsVtxZ;
  Int_t fNumbFlattening;
  Double_t fVtxZ;
  Int_t fNCentBins;
  Double_t fEtaGap;
  Double_t fPsi_E;
  Double_t fPsi_W;
  Double_t fDeltaPsi_E;
  Double_t fDeltaPsi_W;
  Double_t fResolution;
  Double_t fVn;
  QVector *fQvector_E;
  QVector *fQvector_W;
  TH1D *h_VtxZ;
  TProfile2D *tp2_QxE;
  TProfile2D *tp2_QxW;
  TProfile2D *tp2_QyE;
  TProfile2D *tp2_QyW;
  TProfile2D tp2_readQxE;
  TProfile2D tp2_readQxW;
  TProfile2D tp2_readQyE;
  TProfile2D tp2_readQyW;
  std::vector<TProfile2D*> tp2_SinPsiEast;
  std::vector<TProfile2D*> tp2_CosPsiEast;
  std::vector<TProfile2D*> tp2_SinPsiWest;
  std::vector<TProfile2D*> tp2_CosPsiWest;
  std::vector<TH2D*> h2_QEast;
  std::vector<TH2D*> h2_QWest;
  std::map<Int_t,TString> fPidCode;
  TProfile *tp_SqRes;
  TProfile *tp_SqResX;
  TProfile *tp_SqResY;
  std::vector<TProfile2D*> tp2_MeanPt;
  std::vector<TProfile*> tp_VnCent;
  std::vector<TProfile*> tp_VnCentX;
  std::vector<TProfile*> tp_VnCentY;
  std::vector<TProfile2D*> tp2_VnPtCent;
  std::vector<TProfile2D*> tp2_VnPtCentE;
  std::vector<TProfile2D*> tp2_VnPtCentW;
  std::vector<TProfile3D*> tp3_VnPtCentRapidity;
  std::vector<TProfile2D*> tp2_VnRapidityCent;
  std::vector<TProfile2D*> tp2_VnPtCentX;
  std::vector<TProfile2D*> tp2_VnPtCentY;
  std::vector<TProfile3D*> tp3_VnPtCent;
  std::vector<TProfile3D*> tp3_VnPtCentX;
  std::vector<TProfile3D*> tp3_VnPtCentY;

  ClassDef(FlowScalarProduct,0);
};

#endif