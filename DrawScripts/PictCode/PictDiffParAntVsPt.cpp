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

//#include "../macro/ConstantsNew.h"

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,0,0)
R__LOAD_LIBRARY(../build/libDrawPicture.so)
#endif

void PictDiffParAntVsPt(const Char_t *inFileName,const Int_t energy = 28){
	
    TFile *fstyle = new TFile("style.root");
    TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
    tsty->cd();

    gROOT->ForceStyle();
    gStyle->SetOptStat(0);

	std::map<int, string> Energy_scan_text={{7, "7.7"},{11, "11.5"},{14, "14.5"},{19, "19.6"},{27, "27"},{39, "39"},{62, "62.4"}};
	TString particle[]={"h^{+}","h^{-}","#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","Hadrons","","PID"};
    const Int_t cent[10]={80,70,60,50,40,30,20,10,5,0};
  
    TFile *fileIn = new TFile(inFileName,"read");
    //TFile *fileIn = new TFile(Form("/home/aleksandr/STAR_BES/result/root_sys_27GeV/flow_28GeV_PID.root"),"READ");

    Int_t color[]={2, 1, 4, 6, 4, 2, 46};
    Int_t style[]={23, 26, 22, 32, 20, 24, 34, 28};
    Float_t size[]={2, 2, 1.6, 2};
  
    std::vector<Double_t> AxisUp   = {-0.01,0.21};
    std::vector<Double_t> AxisDown = {-0.01,0.01};

    std::vector<MyDrawObject *> drawOdj;
    std::vector<AnalysisFunction *> AnalysOdj;
    Int_t CentMin=2;
    Int_t CentMax=8;
    std::vector<int> fNHarmonic ={2,3,4};
    std::vector<TString> fPidCode = {"HadronsPos","HadronsNeg","PionPos","PionNeg","KoanPos","KaonNeg","ProtonPos","ProtonNeg"};
    //std::vector<TString> fPidCode = {"PionPos","PionNeg","KaonPos","KaonNeg","ProtonPos","ProtonNeg"};
    std::vector<float> fEtaGap = {0.075};
    std::vector<TString> fNameSys = {"","v42"};

    std::vector<std::vector<double>> ptBin         = {
        {0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.4,2.8,3.2},
        {0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.6,3.2},
        {0.2,0.6,1.0,1.4,2.0,2.6,3.2},
        {0.2,0.4,0.6,0.8,1.2,1.6,2.0,2.6,3.2},
    };

    for(Int_t pid=0; pid<4; pid++){
        for(Int_t ih=0; ih<fNHarmonic.size(); ih++){    
            for(Int_t ch=0; ch<2; ch++){
                AnalysOdj.push_back(new AnalysisFunction());
                drawOdj.push_back(new MyDrawObject());
                drawOdj.back()->SetMarkerColor(color[ih]);
                drawOdj.back()->SetMarkerStyle(style[ch+pid*2]);
                drawOdj.back()->SetMarkerSize(size[pid]);
                drawOdj.back()->SetNumberPadDraw(pid);
                if(pid==0 && ch==0)drawOdj.back()->SetLegendText(Form("#font[42]{#scale[1.0]{v_{%i}: full-X, open-#bar{X}}}",fNHarmonic[ih]));
                drawOdj.back()->SetDrawOdject( AnalysOdj.back()->FlowVsPt_EtaSub(fileIn,
                        Form("tp_SqRes%iTPC%i",fNHarmonic[ih],(int)(fEtaGap[0]*100)),
                        Form("tp2_v%iewTPC%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("tp2_meanPtV%i%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("gt_v%i_%s_pt",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        CentMin,CentMax, ptBin[ih]) );

                /*
                                drawOdj.back()->SetDrawOdject( AnalysOdj.back()->FlowVsPt_EtaSub(fileIn,
                        Form("tp_SqRes%iTPC%i",fNHarmonic[ih],(int)(fEtaGap[0]*100)),
                        Form("tp2_v%iewTPC%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("tp2_meanPtV%i%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("gt_v%i_%s_pt",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        CentMin,CentMax, ptBin[ih]) );

                        drawOdj.back()->SetDrawOdject( AnalysOdj.back()->FlowVsPt_EtaSub(fileIn,
                        Form("tp_SqRes%iTPCEta01",fNHarmonic[ih]),
                        Form("tp_v%iewTPC%sEta01TPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        Form("tp_meanPt%sEta01TPCandTOF",fPidCode[pid*2+ch].Data()),
                        Form("gt_v%i_%s_pt",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        CentMin,CentMax, ptBin[ih]) );
                */
            }
            AnalysOdj.push_back(new AnalysisFunction());
            drawOdj.push_back(new MyDrawObject());
            drawOdj.back()->SetMarkerColor(color[ih]);
            drawOdj.back()->SetMarkerStyle(style[pid*2]);
            drawOdj.back()->SetMarkerSize(size[pid]);
            drawOdj.back()->SetNumberPadDraw(pid+4);
            //if(pid==0 && ch==0)drawOdj.back()->SetLegendText(Form("#font[42]{#scale[1.0]{v_{%i}: full-X, open-#bar{X}}}",fNHarmonic[ih]));
            drawOdj.back()->SetDrawOdject( AnalysOdj.back()->DifferenceGraphPointToPoint(
                drawOdj[drawOdj.size()-3]->GetDrawObjectGraph(),
                drawOdj[drawOdj.size()-2]->GetDrawObjectGraph(),
                Form("gtDif_v%i_%s_pt",fNHarmonic[ih],fPidCode[pid*2].Data())));
        }
    }

    std::cout<<energy<<std::endl;

    Draw_Picture_new *canvas = new Draw_Picture_new();
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ v_{n} }}"),"Y", 1,0.3,0.65,0.5,0.0);
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ v_{n}(X)-v_{n}(#bar{X}) }}"),"Y", 1,0.5,0.15,0.5,90.0);
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ p_{T}, GeV/c }}"),"X", 1,0.5,0.5,0.5,0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{Au+Au, #sqrt{S_{NN}}=27 GeV Run18}}" ),"Pad", 0, 0.1, 0.9*AxisUp[1], 0.07, 0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{ Centrality %i-%i %%}}",cent[CentMax+1], cent[CentMin]),"Pad", 0, 0.1, 0.8*AxisUp[1], 0.07, 0.0);
    
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{h^{#pm} }}"  ),"Pad", 0, 1, 0.7*AxisUp[1], 0.07, 0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{#pi^{#pm} }}"),"Pad", 1, 1, 0.9*AxisUp[1], 0.07, 0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{K^{#pm} }}"  ),"Pad", 2, 1, 0.9*AxisUp[1], 0.07, 0.0);
    //canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{p#bar{p} }}" ),"Pad", 3, 1, 0.9*AxisUp[1], 0.07, 0.0);

  
    canvas->SetTLineOnPad(-0.01, 0., 3.2, 0., 1, 3, 7, -1);
    canvas->SetAxisToPad(-0.01, 3.2, 0.1, AxisUp[0], AxisUp[1], 0.1);
    canvas->SetAxisToPad2(-0.01, 3.2, 0.1, AxisDown[0], AxisDown[1], 0.15);
    canvas->SetLegendPar(0.01,0.6,0.5,0.86,1,1);

    TCanvas *result = canvas->CanvasWithRatio(520,360, 4, drawOdj, 0.07);
    result->SaveAs(Form("./PictDiffParAntVsPt_%i_%i_%i.png",energy,cent[CentMax+1], cent[CentMin]));

    delete result;
    delete canvas;

    fileIn->Close();
}