#include <FlowScalarProduct.h>
ClassImp(FlowScalarProduct);

FlowScalarProduct::FlowScalarProduct() :
  fFirstRun(kFALSE),
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
  h2_QEast(),
  h2_QWest(),
  fPidCode(),
  tp_SqRes(),
  tp_SqResX(),
  tp_SqResY(),
  tp2_MeanPt(),
  tp_VnCent(),
  tp2_VnPtCent(),
  tp2_VnPtCentE(),
  tp2_VnPtCentW(),
  tp2_VnRapidityCent(),
  tp3_VnPtCentRapidity(),
  tp_VnCentX(),
  tp2_VnPtCentX(),
  tp_VnCentY(),
  tp2_VnPtCentY(),
  fNameInFile()
{
}

FlowScalarProduct::~FlowScalarProduct()
{
}

void FlowScalarProduct::SetHarmonicNK(Int_t n, Int_t k){
  fNHarmonic = k;
  fKHarmonic = n;
  fKHarmonicThroughNHarmonic = true;
}


void FlowScalarProduct::Zero(){
  fPsi_E = 0.;
  fPsi_W = 0.;
  fDeltaPsi_E = 0.;
  fDeltaPsi_W = 0.;
  fQvector_E->Zero();
  fQvector_W->Zero();
}

void FlowScalarProduct::Init(){
  
  fQvector_E = new QVector(fNHarmonic);
  fQvector_W = new QVector(fNHarmonic);

  h_VtxZ = new TH1D();//(Form("h_VtxZ%iEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),"", fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ);

  if(fFirstRun){
    tp2_QxE = new TProfile2D(Form("tp2_Qx%iEastEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{x}> for #psi_{%i} East, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
    tp2_QxW = new TProfile2D(Form("tp2_Qx%iWestEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{x}> for #psi_{%i} West, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
    tp2_QyE = new TProfile2D(Form("tp2_Qy%iEastEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{y}> for #psi_{%i} East, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
    tp2_QyW = new TProfile2D(Form("tp2_Qy%iWestEta%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("RunId=%i, <Q_{y}> for #psi_{%i} West, |#eta|>%.2f ;VtxZ;cent", fRunId,fNHarmonic,fEtaGap), fNBinsVtxZ, -1.*fVtxZ, 1.*fVtxZ, fNCentBins, -0.5, fNCentBins-0.5);
  }

  if(fCalFlow){
    
    fNameInFile.push_back( Form("tp_SqRes%iTPC%i%s", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data()) );
    fNameInFile.push_back( Form("tp_SqRes%iTPC%i%sX", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data()) );
    fNameInFile.push_back( Form("tp_SqRes%iTPC%i%sY", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data()) );
    
    tp_SqRes  = new TProfile(Form("tp_SqRes%iTPC%i%s_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("Resolution^{2} for v_{%i}, |#eta|>%.2f;cent bin;Res^{2}",fNHarmonic,fEtaGap),fNCentBins, -0.5, fNCentBins-0.5);
    tp_SqResX = new TProfile(Form("tp_SqRes%iTPC%i%sX_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("Resolution^{2} for v_{%i},X, |#eta|>%.2f;cent bin;Res^{2}",fNHarmonic,fEtaGap),fNCentBins, -0.5, fNCentBins-0.5);
    tp_SqResY = new TProfile(Form("tp_SqRes%iTPC%i%sY_Run%i", fNHarmonic, (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("Resolution^{2} for v_{%i},Y, |#eta|>%.2f;cent bin;Res^{2}",fNHarmonic,fEtaGap),fNCentBins, -0.5, fNCentBins-0.5);
    
    for(Int_t i=0; i<(int)fPidCode.size(); i++){
      tp2_MeanPt.push_back(   new TProfile2D(Form("tp2_meanPtV%i%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s, <p_{T}> for bins, v_{%i}, |#eta|>%.2f;p_{T} [GeV/c]; cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp_VnCent.push_back(    new TProfile(Form("tp_v%iewTPC%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(cent), |#eta|>%.2f;cent bin; v_{%i}", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap, fNHarmonic), fNCentBins, -0.5, fNCentBins-0.5));
      tp_VnCentX.push_back(   new TProfile(Form("tp_v%iewTPC%sEta%iTPCandTOF%sX_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(cent), |#eta|>%.2f;cent bin; v_{%i}", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap, fNHarmonic), fNCentBins, -0.5, fNCentBins-0.5));
      tp_VnCentY.push_back(   new TProfile(Form("tp_v%iewTPC%sEta%iTPCandTOF%sY_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(cent), |#eta|>%.2f;cent bin; v_{%i}", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap, fNHarmonic), fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCent.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCentW.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%sW_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T} East, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCentE.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%sE_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T} West, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnRapidityCent.push_back( new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%sRapidity_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;rapidity;cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, -1.0, 1.0, fNCentBins, -0.5, fNCentBins-0.5));
      tp3_VnPtCentRapidity.push_back( new TProfile3D(Form("tp3_v%iewTPC%sEta%iTPCandTOF%s_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin; y", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5, 100, -1., 1.));
      tp2_VnPtCentX.push_back(new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%sX_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      tp2_VnPtCentY.push_back(new TProfile2D(Form("tp2_v%iewTPC%sEta%iTPCandTOF%sY_Run%i", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data(), fRunId),Form("%s v_{%i}(p_{T}, cent), |#eta|>%.2f;p_{T} [GeV/c];cent bin", (fPidCode.at(i)).Data(), fNHarmonic, fEtaGap), 100, 0.15, 5.15, fNCentBins, -0.5, fNCentBins-0.5));
      
      fNameInFile.push_back( Form("tp2_meanPtV%i%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp_v%iewTPC%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp_v%iewTPC%sEta%iTPCandTOF%sX", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp_v%iewTPC%sEta%iTPCandTOF%sY", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%sE", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%sW", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%sRapidity", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp3_v%iewTPC%sEta%iTPCandTOF%s", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%sX", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
      fNameInFile.push_back( Form("tp2_v%iewTPC%sEta%iTPCandTOF%sY", fNHarmonic, (fPidCode.at(i)).Data(), (int)(fEtaGap*100), fNameSys.Data()) );  
 
    }
  }

  if(fChekHisto==true){
    for(int ic=0; ic<fNCentBins; ic++){
      h2_QEast.push_back(new TH2D(Form("h2_Q%iEastEta%icent%i%s_Run%i",fNHarmonic,(int)(fEtaGap*100),ic, fNameSys.Data(), fRunId),Form("Q for #psi_{%i} East cent%i;Q_{x};Q_{y}",fNHarmonic,ic),300,-1.5,1.5,300,-1.5,1.5));
      h2_QWest.push_back(new TH2D(Form("h2_Q%iWestEta%icent%i%s_Run%i",fNHarmonic,(int)(fEtaGap*100),ic, fNameSys.Data(), fRunId),Form("Q for #psi_{%i} West cent%i;Q_{x};Q_{y}",fNHarmonic,ic),300,-1.5,1.5,300,-1.5,1.5));
    }
  }

}

void FlowScalarProduct::SetMeanQVector(){
  
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

void FlowScalarProduct::ProcessFirstTrackLoop(const Double_t &eta, const Double_t &phi, const Double_t &pt){
  if(eta < - fEtaGap){
    fQvector_E->CalQVector(phi, pt);
  }
  if(eta > fEtaGap){
    fQvector_W->CalQVector(phi, pt);
  }
}

void FlowScalarProduct::ProcessSecondTrackLoop(const Double_t &rapidity, const Double_t &eta, const Double_t &phi, const Double_t &pt, const Int_t &pid, const Double_t &dCent){
   
  TVector2 un = TVector2( TMath::Cos(fNHarmonic * phi), TMath::Sin(fNHarmonic * phi) );
  TVector2 vn = TVector2( -999.0, -999.0);

  if(eta < - fEtaGap){
    vn.Set( un.X()*fQvector_W->X(), un.Y() * fQvector_W->Y());
  }
  if(eta > fEtaGap){
    vn.Set( un.X()*fQvector_E->X(), un.Y() * fQvector_E->Y());  
  }

  if(vn.X() != -999.0 && vn.Y() != -999.0  && TMath::Abs(rapidity)<=1.0){
    tp2_MeanPt[pid]   ->Fill(pt,dCent,pt);
    tp_VnCent[pid]    ->Fill(dCent, vn.X() + vn.Y());
    tp_VnCentX[pid]   ->Fill(dCent, vn.X() );
    tp_VnCentY[pid]   ->Fill(dCent, vn.Y() );
    tp2_VnPtCent[pid] ->Fill(pt,dCent,vn.X() + vn.Y());
    if(eta < - fEtaGap)tp2_VnPtCentE[pid] ->Fill(pt,dCent,vn.X() + vn.Y());
    if(eta > fEtaGap)tp2_VnPtCentW[pid] ->Fill(pt,dCent,vn.X() + vn.Y());
    if(pt>0.2 && pt<2.0){
      tp2_VnRapidityCent[pid] ->Fill(rapidity,dCent, vn.X() + vn.Y());
    }
    tp3_VnPtCentRapidity[pid] ->Fill(pt,dCent,rapidity,vn.X() + vn.Y());
    tp2_VnPtCentX[pid]->Fill(pt,dCent,vn.X() );
    tp2_VnPtCentY[pid]->Fill(pt,dCent,vn.Y() );
  }
}

void FlowScalarProduct::ProcessEventAfterFirstTrackLoop(const Double_t &dCent, const Double_t &dVtxZ){
  if(fQvector_E->GetMult()!=0 && fQvector_W->GetMult()!=0 && fQvector_E->GetMod()!=0. && fQvector_W->GetMod()!=0.){
    if(fMeanQVectorForReсentrening){
      fQvector_E->SetMeanQVector(tp2_readQxE.GetBinContent(tp2_readQxE.FindBin(dVtxZ,dCent)), tp2_readQyE.GetBinContent(tp2_readQyE.FindBin(dVtxZ,dCent)));
      fQvector_W->SetMeanQVector(tp2_readQxW.GetBinContent(tp2_readQxW.FindBin(dVtxZ,dCent)), tp2_readQyW.GetBinContent(tp2_readQyW.FindBin(dVtxZ,dCent)));
    }

    fQvector_E->WeightQVector();
    fQvector_W->WeightQVector();

    if (fFirstRun){
      tp2_QxE->Fill(dVtxZ, dCent, fQvector_E->X());
      tp2_QxW->Fill(dVtxZ, dCent, fQvector_W->X());
      tp2_QyE->Fill(dVtxZ, dCent, fQvector_E->Y());
      tp2_QyW->Fill(dVtxZ, dCent, fQvector_W->Y());
    }
    if (fCalFlow){
      tp_SqRes ->Fill(dCent, fQvector_E->GetQvector() * fQvector_W->GetQvector() );
      tp_SqResX->Fill(dCent, fQvector_E->X()          * fQvector_W->X());
      tp_SqResY->Fill(dCent, fQvector_E->Y()          * fQvector_W->Y());
    }
    if(fChekHisto && TMath::Abs(dVtxZ)<fVtxZ){
      h2_QEast[(int)dCent] ->Fill(fQvector_E->X(),fQvector_E->Y());
      h2_QWest[(int)dCent] ->Fill(fQvector_W->X(),fQvector_W->Y());
    }
  }
}

void FlowScalarProduct::WtiteHist(){
  if(fFirstRun){
    tp2_QxE->Write();
    tp2_QxW->Write();
    tp2_QyE->Write();
    tp2_QyW->Write();
  }
  if(fCalFlow){
    tp_SqRes->Write(); 
    tp_SqResX->Write(); 
    tp_SqResY->Write(); 
    for(Int_t i=0; i<(int)fPidCode.size(); i++){
      tp2_MeanPt[i]   ->Write();
      tp_VnCent[i]    ->Write();
      tp_VnCentX[i]   ->Write();
      tp_VnCentY[i]   ->Write();
      tp2_VnPtCent[i] ->Write();
      tp2_VnPtCentE[i] ->Write();
      tp2_VnPtCentW[i] ->Write();
      tp2_VnRapidityCent[i] ->Write();
      tp3_VnPtCentRapidity[i] ->Write();
      tp2_VnPtCentX[i]->Write();
      tp2_VnPtCentY[i]->Write();
    }
  }
}