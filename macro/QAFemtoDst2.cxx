#include <QAFemtoDst2.h>

ClassImp(QAFemtoDst2);

QAFemtoDst2::QAFemtoDst2() : 
  fEnergy(-6666),
  fNameSuf(-6666),
  fRunIdMax(-6666),
  fRunIdMin(-6666),
  fEventCutVtxZ(-6666),
  fEventCutVtxR(-6666),
  fEventCutTofMatched(-6666),
  fTrackCutIsPrimary(kFALSE),
  fPIDcode(-6666),
  fTrackCutNHits(-6666),
  fTrackCutNHitsRatioNHitsPoss(-6666),
  fTrackCutDCA(-6666),
  fTrackCutPseudorapidity(-6666),
  fTrackCutPtotMin(-6666),
  fTrackCutPtMin(-6666),
  fTrackCutNSigmaTPC(),
  fTrackCutM2Min(),
  fTrackCutM2Max(),
  QA_hRefMult(nullptr),
  QA_hVtxZ(nullptr),
  QA_hRefMult2(nullptr),
  QA_hGRefMult(nullptr),
  QA_hNumberOfPrimaries(nullptr),
  QA_hNumberOfGlobals(nullptr),
  QA_hCent9(nullptr),
  QA_hCent16(nullptr),
  QA_hBTofHit(),
  QA_hBTofMatched(),
  QA_hBemcMatched(nullptr),
  QA_hRanking(nullptr),
  QA_hTransSphericity(nullptr),
  QA_hTransSphericity2(nullptr),
  QA_hNumberOfVertices(nullptr),
  QA_hVtxXvsY(nullptr),
  QA_hVtxX(nullptr),
  QA_hVtxY(nullptr),
  QA_hVpdVzDiffVsVz(nullptr),
  QA_hZdcAdcWestEast(nullptr),
  QA_hZdcXVsnBTofMatch(nullptr),
  QA_hBTofTrayMultVsRefMult(nullptr),
  QA_hBEMCMatchVsRefMult(nullptr),
  QA_hBTofMatchedVsRefMult(),
  QA_hEventProfile(),
  QA_hGlobalPtot(nullptr),
  QA_hPrimaryPtot(nullptr),
  QA_hGlobalPt(nullptr),
  QA_hPrimaryPt(nullptr),
  QA_hNHits(nullptr),
  QA_hNHitsRatio(nullptr),
  QA_hChi2(nullptr),
  QA_hDca(nullptr),
  QA_hDcaZ(nullptr),
  QA_hDcaX(nullptr),
  QA_hDcaY(nullptr),
  QA_hDcaVsPt(nullptr),
  QA_hPhi(nullptr),
  QA_hRapidity(nullptr),
  QA_hRapidityPt(nullptr),
  QA_hEta(nullptr),
  QA_hEtaG(nullptr),
  QA_hPtVsEta(nullptr),
  QA_hPhiVsEta(nullptr),
  QA_hPrimaryPhiVsPt(),
  QA_hDedx(nullptr),
  QA_hDedxVsPt(nullptr),
  QA_hNSigmaPionVsPt(nullptr),
  QA_hNSigmaElectronVsPt(nullptr),
  QA_hNSigmaKaonVsPt(nullptr),
  QA_hNSigmaProtonVsPt(nullptr),
  QA_hTofBeta(nullptr),
  QA_hInvBetaVsPt(nullptr),
  QA_hMassSqr(nullptr),
  QA_hMassSqrVsPt(nullptr),
  QA_hDedxVsMassSqr(),
  QA_hInvBetaDiffElectronVsPt(nullptr),
  QA_hInvBetaDiffPionVsPt(nullptr),
  QA_hInvBetaDiffKaonVsPt(nullptr),
  QA_hInvBetaDiffProtonVsPt(nullptr),
  QA_hDedxVsPtPID(),
  QA_hTrackProfile()
{
}

QAFemtoDst2::~QAFemtoDst2(){
}

void QAFemtoDst2::Init(){
  
  QA_hRefMult                 = new TH1D(Form("Event.hRefMult_%s",fNameSuf.Data()), "Reference multiplicity;RefMult;Entries", 1000, -0.5, 999.5);
  QA_hVtxXvsY                 = new TH2D(Form("Event.hVtxXvsY_%s",fNameSuf.Data()), "Vertex X vs Y;x (cm);y (cm)", 1000,-5.,5.,1000,-5.,5.);
  QA_hVtxX                    = new TH1D(Form("Event.hVtxX_%s",fNameSuf.Data()), "Vertex X;x (cm)", 500,-3.,3.);
  QA_hVtxY                    = new TH1D(Form("Event.hVtxY_%s",fNameSuf.Data()), "Vertex Y;y (cm)", 500,-3.,3.);
  //QA_hVtxZ                    = new TH1D(Form("Event.hVtxZ_%s",fNameSuf.Data()),"Vertex Z;z (cm); Entries", 400, -100., 100.);
  QA_hVtxZ                    = new TH1D(Form("Event.hVtxZ_%s",fNameSuf.Data()),"Vertex Z;z (cm); Entries", 500, -80., 80.);
  QA_hRefMult2                = new TH1D(Form("Event.hRefMult2_%s",fNameSuf.Data()),"Reference multiplicity in |#eta|<1;RefMult2;Entries", 2000, -0.5, 1999.5);
  QA_hGRefMult                = new TH1D(Form("Event.hGRefMult_%s",fNameSuf.Data()),"Reference multiplicity of global tracks;gRefMult;Entries", 2000, -0.5, 1999.5);
  QA_hNumberOfPrimaries       = new TH1D(Form("Event.hNumberOfPrimaries_%s",fNameSuf.Data()),"Number of primary tracks;Number of primary tracks;Entries", 1000, -0.5, 999.5);
  QA_hNumberOfGlobals         = new TH1D(Form("Event.hNumberOfGlobals_%s",fNameSuf.Data()),"Number of global tracks;Number of global tracks;Entries", 1000, -0.5, 999.5);
  QA_hCent9                   = new TH1D(Form("Event.hCent9_%s",fNameSuf.Data()),"Centralitity;Cent9;Entries", 13, -1.5, 11.5);
  QA_hCent16                  = new TH1D(Form("Event.hCent16_%s",fNameSuf.Data()),"Centralitity;Cent16;Entries", 19, -1.5, 17.5);
  QA_hBTofHit[0]              = new TH1D(Form("Event.hBTofHit_%s",fNameSuf.Data()),"Number of hits in TOF;bTofTrayMult;Entries", 1500, -0.5, 1499.5);
  QA_hBTofHit[1]              = new TH1D(Form("Event.hBTofHit_old_%s",fNameSuf.Data()),"Number of hits in TOF;bTofTrayMult;Entries", 1500, -0.5, 1499.5);
  QA_hBTofMatched[0]          = new TH1D(Form("Event.hBTofMatched_%s",fNameSuf.Data()),"Number of TOF-matched tracks;bTofMatched;Entries", 700, -0.5, 699.5);
  QA_hBTofMatched[1]          = new TH1D(Form("Event.hBTofMatched_old_%s",fNameSuf.Data()),"Number of TOF-matched tracks;bTofMatched;Entries", 600, -0.5, 599.5);
  QA_hBemcMatched             = new TH1D(Form("Event.hBemcMatched_%s",fNameSuf.Data()),"Number of BEMC-matched tracks;bEmcMatched;Entries", 400, -0.5, 399.5);
  QA_hRanking                 = new TH1D(Form("Event.hRanking_%s",fNameSuf.Data()),"Primary vertex ranking;Primary vertex ranking;Entries", 21, -10.5, 10.5);
  QA_hVpdVzDiffVsVz           = new TH2D(Form("Event.hVpdVzDiffVsVz_%s",fNameSuf.Data()),"v_{z}(TPC) - v_{z}(VPD) vs. v_{z}(TPC);v_{z}(TPC);v_{z}(TPC) - v_{z}(VPD)",280, -70., 70., 80, -20., 20.);
  QA_hBTofTrayMultVsRefMult   = new TH2D(Form("Event.hBTofTrayMultVsRefMult_%s",fNameSuf.Data()),"TOF tray multiplicity vs. refMult;refMult;bTofTrayMult",1600, -0.5, 1599.5, 1500, -0.5, 1499.5);
  QA_hBTofMatchedVsRefMult[0] = new TH2D(Form("Event.hBTofMatchedVsRefMult_%s",fNameSuf.Data()),"TOF-matched tracks vs. refMult;refMult;TOF-matched",1000, -0.5, 999.5, 1000, -0.5, 999.5);
  QA_hBTofMatchedVsRefMult[1] = new TH2D(Form("Event.hBTofMatchedVsRefMult_old_%s",fNameSuf.Data()),"TOF-matched tracks vs. refMult;refMult;TOF-matched",1000, -0.5, 999.5, 1000, -0.5, 999.5);
  QA_hTransSphericity         = new TH1D(Form("Event.hTransSphericity_%s",fNameSuf.Data()),"Transverse sphericity;Sphericity;Entries", 10, 0., 1.);
  QA_hTransSphericity2        = new TH1D(Form("Event.hTransSphericity2_%s",fNameSuf.Data()),"Transverse sphericity in |#eta|<1;Sphericity;Entries", 10, 0., 1.);
  QA_hNumberOfVertices        = new TH1D(Form("Event.hNumberOfVertices_%s",fNameSuf.Data()),"Number of primary vertices;Number of primary vertices;Entries", 15, -0.5, 14.5);
  
  QA_hEventProfile[0]   = new TProfile(Form("Event.TProfile_0_%s",fNameSuf.Data()),"Profile of refMult;Run ID;<refMult>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[1]   = new TProfile(Form("Event.TProfile_1_%s",fNameSuf.Data()),"Profile of TOF tray multiplicity;Run ID;<bTofTrayMultiplicity>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[2]   = new TProfile(Form("Event.TProfile_2_%s",fNameSuf.Data()),"Profile of TOF-matched tracks;Run ID;<bTofMatched>",  fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[3]   = new TProfile(Form("Event.TProfile_3_%s",fNameSuf.Data()),"Profile of number of primary tracks;Run ID;<nPrimTracks>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[4]   = new TProfile(Form("Event.TProfile_4_%s",fNameSuf.Data()),"Profile of number of global tracks;Run ID;<nGlobTracks>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[5]   = new TProfile(Form("Event.TProfile_5_%s",fNameSuf.Data()),"Profile of ZDC ADC;Run ID; <ADC_{ZDC}>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[6]   = new TProfile(Form("Event.TProfile_6_%s",fNameSuf.Data()),"Profile of BBC ADC;Run ID; <ADC_{BBC}>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[7]   = new TProfile(Form("Event.TProfile_7_%s",fNameSuf.Data()),"Profile of primary vertex Z position;Run ID; <VtxZ> (cm)", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[8]   = new TProfile(Form("Event.TProfile_8_%s",fNameSuf.Data()),"Profile of primary vertex X position;Run ID; <VtxX> (cm)", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[9]   = new TProfile(Form("Event.TProfile_9_%s",fNameSuf.Data()),"Profile of primary vertex Y position;Run ID; <VtxY> (cm)", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hEventProfile[10]  = new TProfile(Form("Event.TProfile_2_old_%s",fNameSuf.Data()),"Profile of TOF-matched tracks;Run ID;<bTofMatched>",  fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  
  for(int iHist=0; iHist<11; iHist++) {
    Int_t mColor = 1;
    QA_hEventProfile[iHist]->SetMarkerStyle(20);    // filled circle
    QA_hEventProfile[iHist]->SetMarkerColor(mColor);
    QA_hEventProfile[iHist]->SetMarkerSize(1.1);
    QA_hEventProfile[iHist]->SetLineWidth(2);
    QA_hEventProfile[iHist]->SetLineColor(mColor);  // black
  }

  // Track
  QA_hGlobalPtot  = new TH1D(Form("TpcTrack.hGlobalPtot_%s",fNameSuf.Data()),"Global track momentum;p (GeV/c);Entries", 600, 0., 6. );
  QA_hPrimaryPtot = new TH1D(Form("TpcTrack.hPrimaryPtot_%s",fNameSuf.Data()),"Primary track momentum;p (GeV/c);Entries", 600, 0., 6. );
  QA_hGlobalPt    = new TH1D(Form("TpcTrack.hGlobalPt_%s",fNameSuf.Data()),"Global track transverse momentum;p_{T} (GeV/c)", 600, 0., 6.);
  QA_hPrimaryPt   = new TH1D(Form("TpcTrack.hPrimaryPt_%s",fNameSuf.Data()),"Primary track transverse momentum;p_{T} (GeV/c)", 600, 0., 6.);
  QA_hNHits       = new TH1D(Form("TpcTrack.hNHits_%s",fNameSuf.Data()),"Number of hits;nHits;Entries", 50, 0.0, 50.0);
  QA_hNHitsRatio  = new TH1D(Form("TpcTrack.hNHitsRatio_%s",fNameSuf.Data()),"nHitsFit to nHitsPoss ratio;nHitsFit/nHitsPoss;Entries", 100, 0., 1. );
  QA_hChi2        = new TH1D(Form("TpcTrack.hChi2_%s",fNameSuf.Data()),"#chi^{2} of the track;#chi^{2};Entries", 200, 0., 20.);
  QA_hDca         = new TH1D(Form("TpcTrack.hDca_%s",fNameSuf.Data()),"DCA to primary vertex;DCA (cm);Entries", 200, 0., 5.);
  QA_hDcaZ        = new TH1D(Form("TpcTrack.hDcaZ_%s",fNameSuf.Data()),"DCA to primary vertex;DCA (cm);Entries", 200, -5., 5.);
  QA_hDcaX        = new TH1D(Form("TpcTrack.hDcaX_%s",fNameSuf.Data()),"DCA to primary vertex;DCA (cm);Entries", 200, -5., 5.);
  QA_hDcaY        = new TH1D(Form("TpcTrack.hDcaY_%s",fNameSuf.Data()),"DCA to primary vertex;DCA (cm);Entries", 200, -5., 5.);
  QA_hDcaVsPt     = new TH2D(Form("TpcTrack.hDcaVsPt_%s",fNameSuf.Data()),"charge*p_{T} vs. DCA;charge * p_{T} (GeV/c);DCA (cm)", 840, -2.1, 2.1, 200, 0., 10.);
  QA_hPhi         = new TH1D(Form("TpcTrack.hPhi_%s",fNameSuf.Data()),"Azimuthal angle distribution;#phi;Entries", 1280, -3.2, 3.2 );
  QA_hEta         = new TH1D(Form("TpcTrack.hEta_%s",fNameSuf.Data()),"Track pseudorapidity;#eta;Entries", 400, -2., 2. );
  QA_hRapidity    = new TH1D(Form("TpcTrack.hRapidity_%s",fNameSuf.Data()),"Track rapidity;#eta;Entries", 400, -2., 2. );
  QA_hRapidityPt  = new TH2D(Form("TpcTrack.hRapidityPt_%s",fNameSuf.Data()),"Track rapidity vs p_{T};#rapidity;p_{T}", 400, -2., 2., 600, 0., 6.);
  QA_hEtaG        = new TH1D(Form("TpcTrack.hEtaG_%s",fNameSuf.Data()),"Track pseudorapidity of global track;#eta;Entires", 440, -1., 1. );
  QA_hPtVsEta     = new TH2D(Form("TpcTrack.hPtVsEta_%s",fNameSuf.Data()),"p_{T} vs. #eta of primary track;#eta;p_{T} (GeV/c)", 440, -1.1, 1.1, 160, 0.05, 2.05);
  QA_hPhiVsEta    = new TH2D(Form("TpcTrack.hPhiVsEta_%s",fNameSuf.Data()),"Azimuthal angle vs. #eta of primary track;#eta;#phi", 440, -1.1, 1.1, 1280,-3.2, 3.2);
  
  QA_hPrimaryPhiVsPt[2];
  for(int i=0; i<2; i++) {
    QA_hPrimaryPhiVsPt[i] = new TH2D(Form("TpcTrack.hPrimaryPhiVsPt_%d%s",i,fNameSuf.Data()),
         Form("#phi vs. p_{T} for charge: %d;p_{T} (GeV/c);#phi (rad)", (i==0) ? 1 : -1),
         300, 0., 3., 630, -3.15, 3.15 );
  }
  QA_hDedx                = new TH1D(Form("TpcTrack.hDedx_%s",fNameSuf.Data()),"dE/dx;dE/dx (keV/cm);Entries", 125, 0., 12.5);
  QA_hDedxVsPt            = new TH2D(Form("TpcTrack.hDedxVsPt_%s",fNameSuf.Data()),"dE/dx vs. charge*p_{T};charge * p_{T} (GeV/c);dE/dx (keV/cm)", 840, -2.1, 2.1, 600, 0., 12.);
  QA_hNSigmaPionVsPt      = new TH2D(Form("TpcTrack.hNSigmaPionVsPt_%s",fNameSuf.Data()),"n#sigma(#pi) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(#pi)", 840, -2.1, 2.1, 200, -10., 10.);
  QA_hNSigmaElectronVsPt  = new TH2D(Form("TpcTrack.hNSigmaElectronVsPt_%s",fNameSuf.Data()),"n#sigma(e) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(e)", 840, -2.1, 2.1, 200, -10., 10.);
  QA_hNSigmaKaonVsPt      = new TH2D(Form("TpcTrack.hNSigmaKaonVsPt_%s",fNameSuf.Data()),"n#sigma(K) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(K)", 840, -2.1, 2.1, 200, -10., 10.);
  QA_hNSigmaProtonVsPt    = new TH2D(Form("TpcTrack.hNSigmaProtonVsPt_%s",fNameSuf.Data()),"n#sigma(p) vs. charge*p_{T};charge * p_{T} (GeV/c);n#sigma(p)", 840, -2.1, 2.1, 200, -10., 10.);

  for ( int i=0; i<4; i++ ) {
    TString name = "TpcTrack.hDedxVsPtPID_";
    name += i;
    name += "_" + fNameSuf;
    TString title = "dE/dx vs. charge*p_{T} ";
    switch (i) {
      case 0: title += "|n#sigma(e)| #leq 2;"; break;
      case 1: title += "|n#sigma(#pi)| #leq 2;"; break;
      case 2: title += "|n#sigma(K)| #leq 2;"; break;
      case 3: title += "|n#sigma(p)| #leq 2;"; break;
      default: title += "unknown PID;";
    }
    title += "charge*p_{T} (GeV/c);dE/dx (keV/cm)";
    QA_hDedxVsPtPID[i] = new TH2D(name.Data(), title.Data(), 840, -2.1, 2.1, 600, 0., 12.);
  }

  // TofPidTrait
  QA_hTofBeta     = new TH1D(Form("TpcTrack.hTofBeta_%s",fNameSuf.Data()),"BTofPidTraits #beta;#beta", 2000, 0., 2.);
  QA_hInvBetaVsPt = new TH2D(Form("TpcTrack.hInvBetaVsPt_%s",fNameSuf.Data()),"1/#beta vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta", 840, -2.1, 2.1, 200, 0.8, 2.8);
  QA_hMassSqr     = new TH1D(Form("TpcTrack.hMassSqr_%s",fNameSuf.Data()),"m^{2};m^{2} (GeV/c^{2})^{2};dN/dm^{2} (entries)", 520, -0.1, 5.1 );
  QA_hMassSqrVsPt = new TH2D(Form("TpcTrack.hMassSqrVsPt_%s",fNameSuf.Data()),"m^{2} vs. charge*p_{T};charge * p_{T} (GeV/c);m^{2} (GeV/c^{2})^{2}", 840, -2.1, 2.1, 200, -0.2, 1.8);
  
  QA_hDedxVsMassSqr[0]        = new TH2D(Form("TpcTrack.hDedxVsMassSqr_0_%s",fNameSuf.Data()),"dE/dx vs. mass^{2} charge>0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)", 440, -0.4, 1.8, 250, 0., 12.5 );
  QA_hDedxVsMassSqr[1]        = new TH2D(Form("TpcTrack.hDedxVsMassSqr_1_%s",fNameSuf.Data()),"dE/dx vs. mass^{2} charge<0;m^{2} (GeV/c^{2})^{2};dE/dx (keV/cm)", 440, -0.4, 1.8, 250, 0., 12.5 );
  QA_hInvBetaDiffElectronVsPt = new TH2D(Form("TpcTrack.hInvBetaDiffElectronVsPt_%s",fNameSuf.Data()),"1/#beta - 1/#beta(electron) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(e)", 840, -2.1, 2.1, 200, -0.1, 0.1);
  QA_hInvBetaDiffPionVsPt     = new TH2D(Form("TpcTrack.hInvBetaDiffPionVsPt_%s",fNameSuf.Data()),"1/#beta - 1/#beta(pion) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(#pi)", 840, -2.1, 2.1, 200, -0.1, 0.1);
  QA_hInvBetaDiffKaonVsPt     = new TH2D(Form("TpcTrack.hInvBetaDiffKaonVsPt_%s",fNameSuf.Data()),"1/#beta - 1/#beta(kaon) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(K)", 840, -2.1, 2.1, 200, -0.1, 0.1);
  QA_hInvBetaDiffProtonVsPt   = new TH2D(Form("TpcTrack.hInvBetaDiffProtonVsPt_%s",fNameSuf.Data()),"1/#beta - 1/#beta(p) vs. charge*p_{T};charge * p_{T} (GeV/c);1/#beta - 1/#beta(p)", 840, -2.1, 2.1, 200, -0.1, 0.1);
  
  QA_hTrackProfile[0] = new TProfile(Form("TpcTrack.TProfile_0_%s",fNameSuf.Data()),"Profile of track #phi;Run ID;<#phi>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hTrackProfile[1] = new TProfile(Form("TpcTrack.TProfile_1_%s",fNameSuf.Data()),"Profile of track p_{T};Run ID;<p_{T}>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hTrackProfile[2] = new TProfile(Form("TpcTrack.TProfile_2_%s",fNameSuf.Data()),"Profile of track nHits;Run ID;<nHits>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hTrackProfile[3] = new TProfile(Form("TpcTrack.TProfile_3_%s",fNameSuf.Data()),"Profile of track DCA;Run ID;<DCA>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hTrackProfile[4] = new TProfile(Form("TpcTrack.TProfile_4_%s",fNameSuf.Data()),"Profile of track #beta;Run ID;<#beta>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hTrackProfile[5] = new TProfile(Form("TpcTrack.TProfile_5_%s",fNameSuf.Data()),"Profile of track dE/dx;Run ID;<dE/dx> (keV/cm)", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );
  QA_hTrackProfile[6] = new TProfile(Form("TpcTrack.TProfile_6_%s",fNameSuf.Data()),"Profile of track nHitsFit to nHitsPoss ratio;Run ID; <nHitsFit/nHitsPoss>", fRunIdMax - fRunIdMin, fRunIdMin, fRunIdMax );

  for(int iTrk=0; iTrk<7; iTrk++) {
    Int_t mColor = 1;
    QA_hTrackProfile[iTrk]->SetMarkerStyle(20);    // filled circle
    QA_hTrackProfile[iTrk]->SetMarkerColor(mColor);
    QA_hTrackProfile[iTrk]->SetMarkerSize(1.1);
    QA_hTrackProfile[iTrk]->SetLineWidth(2);
    QA_hTrackProfile[iTrk]->SetLineColor(mColor);  // black
  }

}

void QAFemtoDst2::FillEventHisto(
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
  const Float_t bbcAdcSum)
{  
 
  Bool_t fFillKey = true;
  
  if(fEventCutVtxZ != -6666 && TMath::Abs(pVtx.Z()) > fEventCutVtxZ) fFillKey = false;
  //if(fEventCutVtxR != -6666 && TMath::Sqrt(pVtx.X()*pVtx.X() + pVtx.Y()*pVtx.X()) > fEventCutVtxR) fFillKey = false;
  if(fEventCutVtxR != -6666 && TMath::Abs(pVtx.X()) > fEventCutVtxR) fFillKey = false;
  if(fEventCutVtxR != -6666 && TMath::Abs(pVtx.Y()) > fEventCutVtxR) fFillKey = false;
  if(fEventCutTofMatched != -6666 && TofMatched < fEventCutTofMatched) fFillKey = false;

  if(fFillKey==true){
    // Fill event histograms
    QA_hRefMult                 ->Fill( refMult );
    QA_hVtxXvsY                 ->Fill( pVtx.X(), pVtx.Y() );
    QA_hVtxX                    ->Fill( pVtx.X() );
    QA_hVtxY                    ->Fill( pVtx.Y() );
    QA_hVtxZ                    ->Fill( pVtx.Z() );
    QA_hRefMult2                ->Fill( refMult2 );
    QA_hGRefMult                ->Fill( gRefMult );
    QA_hNumberOfPrimaries       ->Fill( numberOfPrimaryTracks );
    QA_hNumberOfGlobals         ->Fill( numberOfGlobalTracks );
    QA_hCent9                   ->Fill( cent9 );
    QA_hCent16                  ->Fill( cent16 );
    QA_hBTofHit[0]              ->Fill( numberOfBTofHit );
    QA_hBTofMatched[0]          ->Fill( numberOfTofMatched );
    QA_hBTofMatched[1]          ->Fill( numberOfTofMatchedOld );
    //QA_hBemcMatched->Fill( event->numberOfBEMCMatched() );
    QA_hRanking                 ->Fill( ranking );
    QA_hVpdVzDiffVsVz           ->Fill( pVtx.Z(), pVtx.Z() - vpdVz );
    QA_hBTofTrayMultVsRefMult   ->Fill( refMult, numberOfBTofHit );
    QA_hBTofMatchedVsRefMult[0] ->Fill( refMult, numberOfTofMatched );
    QA_hBTofMatchedVsRefMult[1] ->Fill( refMult, numberOfTofMatchedOld );
    QA_hTransSphericity         ->Fill( transverseSphericity );
    QA_hTransSphericity2        ->Fill( transverseSphericity2 );
    QA_hNumberOfVertices        ->Fill( numberOfPrimaryVertices );
    QA_hEventProfile[0]         ->Fill( run, gRefMult );
    QA_hEventProfile[1]         ->Fill( run, numberOfBTofHit );
    QA_hEventProfile[2]         ->Fill( run, numberOfTofMatched );
    QA_hEventProfile[3]         ->Fill( run, numberOfPrimaryTracks );
    QA_hEventProfile[4]         ->Fill( run, numberOfGlobalTracks );
    QA_hEventProfile[5]         ->Fill( run, zdcSumAdcEast + zdcSumAdcWest );    
    QA_hEventProfile[6]         ->Fill( run, bbcAdcSum );
    QA_hEventProfile[7]         ->Fill( run, pVtx.Z() );
    QA_hEventProfile[8]         ->Fill( run, pVtx.X() );
    QA_hEventProfile[9]         ->Fill( run, pVtx.Y() );

  }

}

void QAFemtoDst2::FillTrackHisto(
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
  const Double_t massSqr)
{

  Bool_t fFillKey = true;

  if(fEventCutVtxZ != -6666 && TMath::Abs(pVtx.Z()) > fEventCutVtxZ) fFillKey = false;
  //if(fEventCutVtxR != -6666 && TMath::Sqrt(pVtx.X()*pVtx.X() + pVtx.Y()*pVtx.X()) > fEventCutVtxR) fFillKey = false;
  if(fEventCutVtxR != -6666 && TMath::Abs(pVtx.X()) > fEventCutVtxR) fFillKey = false;
  if(fEventCutVtxR != -6666 && TMath::Abs(pVtx.Y()) > fEventCutVtxR) fFillKey = false;
  if(fEventCutTofMatched != -6666 && TofMatched < fEventCutTofMatched) fFillKey = false;

  //if(TMath::Abs(pVtx.Z()) >= fEventCutVtxZ) fFillKey = false;
  //if(TMath::Sqrt( pVtx.X()*pVtx.X() + pVtx.Y()*pVtx.X() ) >= fEventCutVtxR) fFillKey = false;
  //if(TofMatched <= fEventCutTofMatched) fFillKey = false;

  if(fTrackCutNHits != -6666 && nHits < fTrackCutNHits) fFillKey = false;
  if(fTrackCutNHitsRatioNHitsPoss != -6666 && (Float_t)nHits/nHitsPoss < fTrackCutNHitsRatioNHitsPoss) fFillKey = false;
  if(fTrackCutDCA != -6666 && gDCA.Mag() > fTrackCutDCA) fFillKey = false;
  if(fTrackCutPseudorapidity != -6666 && TMath::Abs(pMom.Eta()) > fTrackCutPseudorapidity) fFillKey = false;
  if(fTrackCutPtotMin != -6666 && pMom.Mag() < fTrackCutPtotMin) fFillKey = false;
  if(fTrackCutPtMin != -6666 && pMom.Pt() < fTrackCutPtMin) fFillKey = false;
  
  if(fPIDcode != -6666){
    Double_t ArreySigma[4] = {nSigmaElectron,nSigmaPion,nSigmaKaon, nSigmaProton};
    if( !(TMath::Abs( ArreySigma[fPIDcode] ) <= fTrackCutNSigmaTPC[fPIDcode] && massSqr>=fTrackCutM2Min[fPIDcode] && massSqr<=fTrackCutM2Max[fPIDcode]) )fFillKey = false;
  }
  
  if(fFillKey == true){

    // Fill track histograms
    QA_hGlobalPtot->Fill( gMom.Mag() );
    QA_hPrimaryPtot->Fill( pMom.Mag() );
    QA_hGlobalPt->Fill( gMom.Pt() );
    QA_hPrimaryPt->Fill( pMom.Pt() );
    QA_hNHits->Fill( nHits );
    QA_hNHitsRatio->Fill( (Float_t)nHits/nHitsPoss );
    QA_hChi2->Fill( chi2 );
    QA_hDca ->Fill( gDCA.Mag());
    QA_hDcaZ->Fill( gDCA.Z() );
    QA_hDcaX->Fill( gDCA.X() );
    QA_hDcaY->Fill( gDCA.Y() );
    QA_hDcaVsPt->Fill( charge * pMom.Pt(), gDCA.Mag() );
    QA_hPhiVsEta->Fill( pMom.Eta(), pMom.Phi() );
    QA_hPtVsEta->Fill( pMom.Eta(), pMom.Pt() );
    QA_hPhi->Fill( pMom.Phi() );
    QA_hRapidity->Fill( Rapidity );
    QA_hRapidityPt->Fill( Rapidity, pMom.Pt() );
    QA_hEta->Fill( pMom.Eta() );
    QA_hEtaG->Fill( gMom.Eta() );
    QA_hDedx->Fill( dEdx * 1e6 );
    QA_hDedxVsPt->Fill( charge * pMom.Pt(), dEdx * 1e6 );

    // If electron has passed PID nsigma cut
    if ( TMath::Abs( nSigmaElectron ) <= fTrackCutNSigmaTPC[0] && massSqr>=fTrackCutM2Min[0] && massSqr<=fTrackCutM2Max[0]) {
      QA_hDedxVsPtPID[0]->Fill( charge * pMom.Pt(), dEdx * 1e6 );
    }
    // If pion has passed PID nsigma cut
    if ( TMath::Abs( nSigmaPion ) <= fTrackCutNSigmaTPC[1] && massSqr>=fTrackCutM2Min[1] && massSqr<=fTrackCutM2Max[1] ) {
      QA_hDedxVsPtPID[1]->Fill( charge * pMom.Pt(), dEdx * 1e6 );
    }
    // If kaon has passed PID nsigma cut
    if ( TMath::Abs( nSigmaKaon ) <= fTrackCutNSigmaTPC[2] && massSqr>=fTrackCutM2Min[2] && massSqr<=fTrackCutM2Max[2] ) {
      QA_hDedxVsPtPID[2]->Fill( charge * pMom.Pt(), dEdx * 1e6 );
    }
    // If proton has passed PID nsigma cut
    if ( TMath::Abs( nSigmaProton ) <= fTrackCutNSigmaTPC[3] && massSqr>=fTrackCutM2Min[3] && massSqr<=fTrackCutM2Max[3] ) {
      QA_hDedxVsPtPID[3]->Fill( charge * pMom.Pt(), dEdx * 1e6 );
    }

    if( charge > 0 ) {
      QA_hPrimaryPhiVsPt[0]->Fill( pMom.Pt(), pMom.Phi() );
    }
    else {
      QA_hPrimaryPhiVsPt[1]->Fill( pMom.Pt(), pMom.Phi() );
    }

    QA_hNSigmaElectronVsPt  ->Fill( charge * pMom.Pt(), nSigmaElectron );
    QA_hNSigmaPionVsPt      ->Fill( charge * pMom.Pt(), nSigmaPion );
    QA_hNSigmaKaonVsPt      ->Fill( charge * pMom.Pt(), nSigmaKaon );
    QA_hNSigmaProtonVsPt    ->Fill( charge * pMom.Pt(), nSigmaProton );
    QA_hTrackProfile[0]     ->Fill( run, pMom.Phi() );
    QA_hTrackProfile[1]     ->Fill( run, pMom.Pt() );
    QA_hTrackProfile[2]     ->Fill( run, nHits );
    QA_hTrackProfile[3]     ->Fill( run, gDCA.Mag() );
    QA_hTrackProfile[5]     ->Fill( run, dEdx * 1e6 );
    QA_hTrackProfile[6]     ->Fill( run, (Float_t)nHits/nHitsPoss );

    // Check if track has TOF signal
    if ( isTofTrack ) {

      QA_hTofBeta                   ->Fill( beta );
      QA_hInvBetaVsPt               ->Fill( charge * pMom.Pt(), invBeta );
      QA_hMassSqr                   ->Fill( massSqr );
      QA_hMassSqrVsPt               ->Fill( charge * pMom.Pt(), massSqr );
      QA_hInvBetaDiffElectronVsPt   ->Fill( charge * pMom.Pt(), invBeta - TMath::Sqrt(electron_mass_sqr + pMom.Mag2() )/ pMom.Mag() );
      QA_hInvBetaDiffPionVsPt       ->Fill( charge * pMom.Pt(), invBeta - TMath::Sqrt(pion_mass_sqr     + pMom.Mag2() )/ pMom.Mag() );
      QA_hInvBetaDiffKaonVsPt       ->Fill( charge * pMom.Pt(), invBeta - TMath::Sqrt(kaon_mass_sqr     + pMom.Mag2() )/ pMom.Mag() );
      QA_hInvBetaDiffProtonVsPt     ->Fill( charge * pMom.Pt(), invBeta - TMath::Sqrt(proton_mass_sqr   + pMom.Mag2() )/ pMom.Mag() );
      if ( charge > 0 ) {
        QA_hDedxVsMassSqr[0]->Fill( massSqr, dEdx * 1e6 );
      }
      else {
        QA_hDedxVsMassSqr[1]->Fill( massSqr, dEdx * 1e6 );
      }
      QA_hTrackProfile[4]->Fill( run, beta );
    } //if( isTofTrack() )

  }

}