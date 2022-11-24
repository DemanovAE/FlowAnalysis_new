#include <FlowEtaSubEP.h>
ClassImp(FlowEtaSubEP);

FlowEtaSubEP::FlowEtaSubEP() :
  fFirstRun(kFALSE),
  fSecondRun(kFALSE),
  fCalFlow(kFALSE),
  fMeanQVectorForReсentrening(kFALSE),
  fMeanSinCosPsiForFlattening(kFALSE),
  fChekHisto(kTRUE),
  fKHarmonicThroughNHarmonic(kFALSE),
  fNameSys(""),
  fRunId(0),
  fNHarmonic(0),
  fKHarmonic(0),
  fNumbFlattening(0),
  fNBinsVtxZ(2),
  fVtxZ(40.),
  fNCentBins(9),
  fEtaGap(0.),
  fPsi_E(0.),
  fPsi_W(0.),
  fDeltaPsi_E(0.),
  fDeltaPsi_W(0.),
  fResolution(0.),
  fVn(0.),
  fQvector_E(nullptr),
  fQvector_W(nullptr),
  h_VtxZ(nullptr),
  tp2_QxE(nullptr),
  tp2_QxW(nullptr),
  tp2_QyE(nullptr),
  tp2_QyW(nullptr),
  tp2_readQxE(),
  tp2_readQxW(),
  tp2_readQyE(),
  tp2_readQyW(),
  tp2_SinPsiEast(),
  tp2_CosPsiEast(),
  tp2_SinPsiWest(),
  tp2_CosPsiWest(),
  tp2_ReadSinPsiEast(),
  tp2_ReadCosPsiEast(),
  tp2_ReadSinPsiWest(),
  tp2_ReadCosPsiWest(),
  h2_QEast(),
  h2_QWest(),
  h_PsiEast(),
  h_PsiWest(),
  fPidCode(),
  tp_SqRes(),
  tp2_MeanPt(),
  tp_VnCent(),
  tp2_VnPtCent(),
  tp3_VnPtCentRapidity(),
  tp_VnCent1(),
  tp2_VnPtCent1(),
  tp_VnCent2(),
  tp2_VnPtCent2()
{
}

FlowEtaSubEP::~FlowEtaSubEP()
{
}

void FlowEtaSubEP::SetHarmonicNK(Int_t n, Int_t k){
  fNHarmonic = k;
  fKHarmonic = n;
  fKHarmonicThroughNHarmonic = true;
}


void FlowEtaSubEP::Zero(){
  fPsi_E = 0.;
  fPsi_W = 0.;
  fDeltaPsi_E = 0.;
  fDeltaPsi_W = 0.;
  fQvector_E->Zero();
  fQvector_W->Zero();
}

void FlowEtaSubEP::Init(){
  
  fQvector_E = new QVector(fNHarmonic);
  fQvector_W = new QVector(fNHarmonic);

  h_VtxZ = new TH1D();//(Form("h_VtxZ%iEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),"", fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ);

  if(fFirstRun){
    tp2_QxE = new TProfile2D(Form("tp2_Qx%iEastEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{x}> for #psi_{%i} East, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
    tp2_QxW = new TProfile2D(Form("tp2_Qx%iWestEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{x}> for #psi_{%i} West, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
    tp2_QyE = new TProfile2D(Form("tp2_Qy%iEastEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{y}> for #psi_{%i} East, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
    tp2_QyW = new TProfile2D(Form("tp2_Qy%iWestEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{y}> for #psi_{%i} West, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
  }

  if(fSecondRun){
    for(Int_t ifl=0; ifl<fNumbFlattening;ifl++){
      tp2_SinPsiEast.push_back(new TProfile2D(Form("tp2_Sin%iPsi%iEastEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <sin(%i*%i#psi_{%i})> East, |#eta|>%.2f ;VtxZ;cent",fRunId,(ifl+1),fNHarmonic,fNHarmonic,fEtaGap),fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_CosPsiEast.push_back(new TProfile2D(Form("tp2_Cos%iPsi%iEastEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <cos(%i*%i#psi_{%i})> East, |#eta|>%.2f ;VtxZ;cent",fRunId,(ifl+1),fNHarmonic,fNHarmonic,fEtaGap),fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_SinPsiWest.push_back(new TProfile2D(Form("tp2_Sin%iPsi%iWestEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <sin(%i*%i#psi_{%i})> West, |#eta|>%.2f ;VtxZ;cent",fRunId,(ifl+1),fNHarmonic,fNHarmonic,fEtaGap),fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_CosPsiWest.push_back(new TProfile2D(Form("tp2_Cos%iPsi%iWestEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <cos(%i*%i#psi_{%i})> West, |#eta|>%.2f ;VtxZ;cent",fRunId,(ifl+1),fNHarmonic,fNHarmonic,fEtaGap),fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5));
    }
  }

  double binCent[9] = {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,8.5};

  if(fCalFlow){
    fNameInFile.push_back( Form("tp_SqRes%iTPC%i%s", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data()) );
    tp_SqRes = new TProfile(Form("tp_SqRes%iTPC%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("Resolution^{2} for v_{%i}, |#eta|>%.2f;cent bin;Res^{2}",fNHarmonic,fEtaGap),fNCentBins, -0.5, fNCentBins-0.5);
    for(Int_t i=0; i<(int)fPidCode.size(); i++){
      tp2_MeanPt.push_back(   new TProfile2D(Form("tp2_meanPtV%i%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s, <p_{T}> for bins, v_{%i}, |#eta|>%.2f;p_{T} [GeV/c]; cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp_VnCent.push_back(    new TProfile(Form("tp_v%iewTPC%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(cent), |#eta|>%.2f;cent bin; v_{%i}", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap, fNHarmonic), fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCent.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp3_VnPtCentRapidity.push_back( new TProfile3D(Form("tp3_v%iewTPC%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent, #eta), |#eta|>%.2f;p_{T},[GeV/c];cent,%%; #eta", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5, 100, -1.,1.0));

      tp_VnCent1.push_back(    new TProfile(Form("tp_v%iewTPC%sEta%iTPCandTOF%s1_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(cent), |#eta|>%.2f;cent bin; v_{%i}", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap, fNHarmonic), fNCentBins, -0.5, fNCentBins-0.5));
      tp_VnCent2.push_back(    new TProfile(Form("tp_v%iewTPC%sEta%iTPCandTOF%s2_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(cent), |#eta|>%.2f;cent bin; v_{%i}", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap, fNHarmonic), fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCent1.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%s1_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCent2.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%s2_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      
      fNameInFile.push_back( Form("tp2_meanPtV%i%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp_v%iewTPC%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp_v%iewTPC%sEta%iTPCandTOF%s1", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp_v%iewTPC%sEta%iTPCandTOF%s2", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp3_v%iewTPC%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%s1", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%s2", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
 
    }
  }

  if(fChekHisto==true){
    for(int ic=0; ic<fNCentBins; ic++){
      h2_QEast.push_back(new TH2D(Form("h2_Q%iEastEta%icent%i%s_Run%i",fNHarmonic,(int)(fEtaGap*100),ic, fNameSys.Data(), fRunId),Form("Q for #psi_{%i} East cent%i;Q_{x};Q_{y}",fNHarmonic,ic),300,-1.5,1.5,300,-1.5,1.5));
      h2_QWest.push_back(new TH2D(Form("h2_Q%iWestEta%icent%i%s_Run%i",fNHarmonic,(int)(fEtaGap*100),ic, fNameSys.Data(), fRunId),Form("Q for #psi_{%i} West cent%i;Q_{x};Q_{y}",fNHarmonic,ic),300,-1.5,1.5,300,-1.5,1.5));
      h_PsiEast.push_back(new TH1D(Form("h_Psi%iEastEta%icent%i%s_Run%i",fNHarmonic,(int)(fEtaGap*100),ic, fNameSys.Data(), fRunId),Form("#psi_{%i} East cent%i;#psi_{%i}; count",fNHarmonic,ic,fNHarmonic), (int)( 2.*(TMath::Pi()+0.1) / fNHarmonic  / 0.02), -0.1, 2.*TMath::Pi()/fNHarmonic + 0.1));
      h_PsiWest.push_back(new TH1D(Form("h_Psi%iWestEta%icent%i%s_Run%i",fNHarmonic,(int)(fEtaGap*100),ic, fNameSys.Data(), fRunId),Form("#psi_{%i} West cent%i;#psi_{%i}; count",fNHarmonic,ic,fNHarmonic), (int)( 2.*(TMath::Pi()+0.1) / fNHarmonic  / 0.02), -0.1, 2.*TMath::Pi()/fNHarmonic + 0.1));
    }
  }

}

void FlowEtaSubEP::SetMeanQVector(){
  
  if (fstrInputFileFromFirstRun == "") { 
    std::cerr << "Warning: in FlowAnalysisWithEtaSubEventPlane::GetMeanQVector() fstrInputFileFromFirstRun="" " << std::endl;
  }
  if(!fFirstRun){
    
    TFile *fi = new TFile(fstrInputFileFromFirstRun.Data(), "read");
    tp2_QxE = dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Qx%iEastEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId)));
    tp2_QxW = dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Qx%iWestEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId)));
    tp2_QyE = dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Qy%iEastEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId)));
    tp2_QyW = dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Qy%iWestEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId)));
    
    if (!tp2_QxE || !tp2_QxW || !tp2_QyE || !tp2_QyW){
      std::cerr << "Cannot find histograms from first run for eta-sub event plane method." << std::endl;
      return;
    }   
    tp2_QxE->Copy(tp2_readQxE);
    tp2_QxW->Copy(tp2_readQxW);
    tp2_QyE->Copy(tp2_readQyE);
    tp2_QyW->Copy(tp2_readQyW);
    fMeanQVectorForReсentrening = true;
    fi->Close();
  }
}

void FlowEtaSubEP::SetMeanSinCosPsi(){
  
  if (fstrInputFileFromSecondRun == "") { 
    std::cerr << "Warning: in FlowAnalysisWithEtaSubEventPlane::GetMeanQVector() fstrInputFileFromSecondRun="" " << std::endl;
  }

  if(!fFirstRun && !fSecondRun){
    
    TFile *fi = new TFile(fstrInputFileFromSecondRun.Data(), "read");

    for(Int_t ifl=0; ifl<fNumbFlattening;ifl++){
      tp2_SinPsiEast.push_back( dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Sin%iPsi%iEastEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId))));
      tp2_CosPsiEast.push_back( dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Cos%iPsi%iEastEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId))));
      tp2_SinPsiWest.push_back( dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Sin%iPsi%iWestEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId))));
      tp2_CosPsiWest.push_back( dynamic_cast<TProfile2D*> (fi->FindObjectAny(Form("tp2_Cos%iPsi%iWestEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId))));

      if (!tp2_SinPsiEast[ifl] || !tp2_CosPsiEast[ifl] || !tp2_SinPsiWest[ifl] || !tp2_CosPsiWest[ifl]){
        std::cerr << "Cannot find histograms from first run for eta-sub event plane method, Number iFlat="<<ifl<< std::endl;
        //std::cout<<Form("tp2_Sin%iPsi%iEastEta%i%s_Run%i", (ifl+1), fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId)<<"\n";
        return;
      }
      tp2_ReadSinPsiEast.push_back(*tp2_SinPsiEast[ifl]);
      tp2_ReadCosPsiEast.push_back(*tp2_CosPsiEast[ifl]);
      tp2_ReadSinPsiWest.push_back(*tp2_SinPsiWest[ifl]);
      tp2_ReadCosPsiWest.push_back(*tp2_CosPsiWest[ifl]);
      fMeanSinCosPsiForFlattening = true;
    }
    fi->Close();
  }
}

void FlowEtaSubEP::ProcessFirstTrackLoop(const Double_t &eta, const Double_t &phi, const Double_t &pt){
  if(eta < - fEtaGap){
    fQvector_E->CalQVector(phi, pt);
  }
  if(eta > fEtaGap){
    fQvector_W->CalQVector(phi, pt);
  }
}

void FlowEtaSubEP::ProcessSecondTrackLoop(const Double_t &rapidity, const Double_t &eta, const Double_t &phi, const Double_t &pt, const Int_t &pid, const Double_t &dCent){
  
  Double_t VnEtaSubEventPlane = -999.0;
  if(eta < - fEtaGap){
    VnEtaSubEventPlane = TMath::Cos( ((fKHarmonicThroughNHarmonic) ? fKHarmonic : fNHarmonic) * (phi - fPsi_W) );
  }
  if(eta > fEtaGap){
    VnEtaSubEventPlane = TMath::Cos( ((fKHarmonicThroughNHarmonic) ? fKHarmonic : fNHarmonic) * (phi - fPsi_E) );
  }
  if(VnEtaSubEventPlane != -999.0){
    tp2_MeanPt[pid]  ->Fill(pt,dCent,pt);
    tp_VnCent[pid]   ->Fill(dCent,VnEtaSubEventPlane);
    tp2_VnPtCent[pid]->Fill(pt,dCent,VnEtaSubEventPlane);
    tp3_VnPtCentRapidity[pid]->Fill(pt,dCent,rapidity,VnEtaSubEventPlane);

    tp_VnCent1[pid]   ->Fill(dCent,TMath::Cos( ((fKHarmonicThroughNHarmonic) ? fKHarmonic : fNHarmonic) * (phi - fPsi_W) ));
    tp2_VnPtCent1[pid]->Fill(pt,dCent,TMath::Cos( ((fKHarmonicThroughNHarmonic) ? fKHarmonic : fNHarmonic) * (phi - fPsi_W) ));

    tp_VnCent2[pid]   ->Fill(dCent,TMath::Cos( ((fKHarmonicThroughNHarmonic) ? fKHarmonic : fNHarmonic) * (phi - fPsi_E) ));
    tp2_VnPtCent2[pid]->Fill(pt,dCent,TMath::Cos( ((fKHarmonicThroughNHarmonic) ? fKHarmonic : fNHarmonic) * (phi - fPsi_E) ));

  }
}

void FlowEtaSubEP::ProcessEventAfterFirstTrackLoop(const Double_t &dCent, const Double_t &dVtxZ){
  if(fQvector_E->GetMult()!=0 && fQvector_W->GetMult()!=0 && fQvector_E->GetMod()!=0. && fQvector_W->GetMod()!=0.){
    if(fMeanQVectorForReсentrening){
      fQvector_E->SetMeanQVector(tp2_readQxE.GetBinContent(tp2_readQxE.FindBin(dVtxZ,dCent)), tp2_readQyE.GetBinContent(tp2_readQyE.FindBin(dVtxZ,dCent)));
      fQvector_W->SetMeanQVector(tp2_readQxW.GetBinContent(tp2_readQxW.FindBin(dVtxZ,dCent)), tp2_readQyW.GetBinContent(tp2_readQyW.FindBin(dVtxZ,dCent)));
    }
    fQvector_E->WeightQVector();
    fQvector_W->WeightQVector();
    fPsi_E = fQvector_E->GetPsi() / fNHarmonic;
    fPsi_W = fQvector_W->GetPsi() / fNHarmonic;

    if(fMeanSinCosPsiForFlattening){
      fDeltaPsi_E = 0.;
      fDeltaPsi_W = 0.;
      for(Int_t ifl=0; ifl<fNumbFlattening; ifl++){  
        fDeltaPsi_E += ( 2.0 / (Double_t)(ifl+1) )*(-1. * tp2_ReadSinPsiEast[ifl].GetBinContent(tp2_ReadSinPsiEast[ifl].FindBin(dVtxZ,dCent)) * TMath::Cos( (Double_t)(ifl+1)*fNHarmonic*fPsi_E ) + 
                                                          tp2_ReadCosPsiEast[ifl].GetBinContent(tp2_ReadCosPsiEast[ifl].FindBin(dVtxZ,dCent)) * TMath::Sin( (Double_t)(ifl+1)*fNHarmonic*fPsi_E ));
        fDeltaPsi_W += ( 2.0 / (Double_t)(ifl+1) )*(-1. * tp2_ReadSinPsiWest[ifl].GetBinContent(tp2_ReadSinPsiWest[ifl].FindBin(dVtxZ,dCent)) * TMath::Cos( (Double_t)(ifl+1)*fNHarmonic*fPsi_W ) + 
                                                          tp2_ReadCosPsiWest[ifl].GetBinContent(tp2_ReadCosPsiWest[ifl].FindBin(dVtxZ,dCent)) * TMath::Sin( (Double_t)(ifl+1)*fNHarmonic*fPsi_W ));
      }
    }

    fPsi_E += fDeltaPsi_E / fNHarmonic;
    fPsi_W += fDeltaPsi_W / fNHarmonic;

    if (fFirstRun){
      tp2_QxE->Fill(dVtxZ, dCent, fQvector_E->X());
      tp2_QxW->Fill(dVtxZ, dCent, fQvector_W->X());
      tp2_QyE->Fill(dVtxZ, dCent, fQvector_E->Y());
      tp2_QyW->Fill(dVtxZ, dCent, fQvector_W->Y());
    }
    if (fSecondRun){
      for(Int_t ifl=0; ifl<fNumbFlattening; ifl++){
        tp2_SinPsiEast[ifl]->Fill(dVtxZ, dCent, TMath::Sin( (Double_t)(ifl+1)*fNHarmonic*fPsi_E ));
        tp2_CosPsiEast[ifl]->Fill(dVtxZ, dCent, TMath::Cos( (Double_t)(ifl+1)*fNHarmonic*fPsi_E ));
        tp2_SinPsiWest[ifl]->Fill(dVtxZ, dCent, TMath::Sin( (Double_t)(ifl+1)*fNHarmonic*fPsi_W ));
        tp2_CosPsiWest[ifl]->Fill(dVtxZ, dCent, TMath::Cos( (Double_t)(ifl+1)*fNHarmonic*fPsi_W ));
      }
    }
    //TMath::Cos( 2*(Psi2TPC[1][i] - Psi2TPC[0][i]) )
    if (fCalFlow){
      if(!fKHarmonicThroughNHarmonic)tp_SqRes->Fill(dCent, TMath::Cos(fNHarmonic*(fPsi_E - fPsi_W)) );
      if( fKHarmonicThroughNHarmonic)tp_SqRes->Fill(dCent, TMath::Cos(fKHarmonic*(fPsi_E - fPsi_W)) );
    }
    if(fChekHisto && TMath::Abs(dVtxZ)<fVtxZ){
      h2_QEast[(int)dCent] ->Fill(fQvector_E->X(),fQvector_E->Y());
      h2_QWest[(int)dCent] ->Fill(fQvector_W->X(),fQvector_W->Y());
      h_PsiEast[(int)dCent]->Fill(fPsi_E);
      h_PsiWest[(int)dCent]->Fill(fPsi_W);
      //[h_VtxZ->FindBin(dVtxZ) - 1]
    }
  }
}

void FlowEtaSubEP::WtiteHist(){
  if(fFirstRun){
    tp2_QxE->Write();
    tp2_QxW->Write();
    tp2_QyE->Write();
    tp2_QyW->Write();
  }
  if(fSecondRun){
    for(Int_t ifl=0; ifl<fNumbFlattening; ifl++){
      tp2_SinPsiEast[ifl]->Write();
      tp2_CosPsiEast[ifl]->Write();
      tp2_SinPsiWest[ifl]->Write();
      tp2_CosPsiWest[ifl]->Write();
    }
  }
  if(fCalFlow){
    tp_SqRes->Write(); 
    for(Int_t i=0; i<(int)fPidCode.size(); i++){
      tp2_MeanPt[i]   ->Write();
      tp_VnCent[i]    ->Write();
      tp2_VnPtCent[i] ->Write();
      tp3_VnPtCentRapidity[i] ->Write();
      tp_VnCent1[i]   ->Write();
      tp2_VnPtCent1[i]->Write();
      tp_VnCent2[i]   ->Write();
      tp2_VnPtCent2[i]->Write();
    }
  }
  /*
  if(fChekHisto){
    for(int ic=0; ic<fNCentBins; ic++){
      h2_QEast[ic]->Write();
      h2_QWest[ic]->Write();
      h_PsiEast[ic]->Write();
      h_PsiWest[ic]->Write();
    }
  }
  */
}