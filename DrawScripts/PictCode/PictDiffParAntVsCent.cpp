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
R__LOAD_LIBRARY(/home/aleksandr/STAR_BES/DrawScript/build/libDrawPicture.so)
#endif

void PictDiffParAntVsCent(const Char_t *inFileName, const Int_t energy = 28){
	
    TFile *fstyle = new TFile("style.root");
    TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
    tsty->cd();

    gROOT->ForceStyle();
    gStyle->SetOptStat(0);

	std::map<int, string> Energy_scan_text={{7, "7.7"},{11, "11.5"},{14, "14.5"},{19, "19.6"},{27, "27"},{39, "39"},{62, "62.4"}};
	TString particle[]={"h^{+}","h^{-}","#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","Hadrons","","PID"};
  
    TFile *fileIn = new TFile(inFileName,"read");
    //TFile *fileIn = new TFile(Form("/home/aleksandr/STAR_BES/result/root_sys_27GeV/flow_28GeV_PID.root"),"READ");

    Int_t color[]={2, 1, 4, 6, 4, 2, 46};
    Int_t style[]={23, 26, 22, 32, 20, 24, 34, 28};
    Float_t size[]={2, 2, 1.6, 2};
  

    std::vector<MyDrawObject *> drawOdj;
    std::vector<AnalysisFunction *> AnalysOdj;
    std::vector<double> minPt={0.2,0.2,0.2,0.5};
    std::vector<double> maxPt={3.2,2.4,2.4,3.2};
    std::vector<int> fNHarmonic ={2,3,4};
    std::vector<TString> fPidCode = {"HadronsPos","HadronsNeg","PionPos","PionNeg","KoanPos","KaonNeg","ProtonPos","ProtonNeg"};
    //std::vector<TString> fPidCode = {"PionPos","PionNeg","KaonPos","KaonNeg","ProtonPos","ProtonNeg"};
    std::vector<float> fEtaGap = {0.05};
    std::vector<TString> fNameSys = {"","v42"};

    std::vector<std::vector<double>> CentBinFemtoDst = {
        {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,1.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,3.5,4.5,5.5,6.5,7.5,8.5}
    };
    std::vector<std::vector<double>> CentBin         = {
        { 80.,70.,60.,50.,40.,30.,20.,10.,5.,0.},
        { 80.,70.,60.,50.,40.,30.,20.,10.,5.,0.},
        { 80.,70.,60.,50.,40.,30.,20.,10.,5.,0.},
        { 80.,70.,60.,50.,40.,30.,20.,10.,5.,0.},
        { 80.,70.,60.,50.,40.,30.,20.,10.,5.,0.},
        { 80.,60.,40.,30.,20.,10.,5.,0.},
        { 80.,40.,30.,20.,10.,5.,0.},
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
                if(ih<3)drawOdj.back()->SetDrawOdject( AnalysOdj.back()->FlowVsCent_EtaSub(fileIn,
                        Form("tp_SqRes%iTPC%i",fNHarmonic[ih],(int)(fEtaGap[0]*100)),
                        Form("tp2_v%iewTPC%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("tp2_meanPtV%i%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("gt_v%i_%s_cent",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        minPt[pid],maxPt[pid], CentBinFemtoDst[ih], CentBin[ih]) );
                if(ih==3){
                    AnalysOdj.back()->SetKeyAnaliticalRes(true,2,4);
                    drawOdj.back()->SetDrawOdject( AnalysOdj.back()->FlowVsCent_EtaSub(fileIn,
                    Form("tp_SqRes%iTPC%i",2,(int)(fEtaGap[0]*100)),
                    Form("tp2_v%iewTPC%sEta%iTPCandTOFv42",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                    Form("tp2_meanPtV%i%sEta%iTPCandTOFv42",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                    Form("gt_v%i_%s_centv42",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                    minPt[pid],maxPt[pid], CentBinFemtoDst[ih], CentBin[ih]) );
                }
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
                Form("gtDif_v%i_%s_cent",fNHarmonic[ih],fPidCode[pid*2].Data())));
        }
    }
                       /* Form("tp_SqRes%iTPC01",fNHarmonic[ih]),
                        Form("tp2_v%iewTPC%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("tp2_meanPtV%i%sEta%iTPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data(),(int)(fEtaGap[0]*100)),
                        Form("gt_v%i_%s_cent",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),


                         Form("tp_SqRes%iTPCEta01",fNHarmonic[ih]),
                        Form("tp_v%iewTPC%sEta01TPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        Form("tp_meanPtV%i%sEta01TPCandTOF",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
                        Form("gt_v%i_%s_cent",fNHarmonic[ih],fPidCode[pid*2+ch].Data()),
	   */
    std::cout<<energy<<std::endl;

    Draw_Picture_new *canvas = new Draw_Picture_new();
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ v_{n} }}"),"Y", 1,0.3,0.65,0.5,0.0);
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ v_{n}(X)-v_{n}(#bar{X}) }}"),"Y", 1,0.5,0.15,0.5,90.0);
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ centtrality}}"),"X", 1,0.5,0.5,0.5,0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{Au+Au, #sqrt{S_{NN}}=27 GeV Run18}}" ),"Pad", 0, 0.1, 0.11, 0.07, 0.0);
    
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{h^{#pm}, %.1f<p_{T}<%.1f GeV/c}}",   minPt[0],maxPt[0]),"Pad", 0, 1, 0.095, 0.07, 0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{#pi^{#pm}, %.1f<p_{T}<%.1f GeV/c}}", minPt[1],maxPt[1]),"Pad", 1, 1, 0.115, 0.07, 0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{K^{#pm}, %.1f<p_{T}<%.1f GeV/c}}",   minPt[2],maxPt[2]),"Pad", 2, 1, 0.115, 0.07, 0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{p#bar{p}, %.1f<p_{T}<%.1f GeV/c}}",  minPt[3],maxPt[3]),"Pad", 3, 1, 0.115, 0.07, 0.0);

  
    canvas->SetTLineOnPad(-1, 0., 85, 0., 1, 3, 7, -1);
    canvas->SetAxisToPad(-1, 85, 0.1, -0.01, 0.13, 0.1);
    canvas->SetAxisToPad2(-1, 85, 0.1, -0.004, 0.02, 0.15);
    canvas->SetLegendPar(0.01,0.6,0.5,0.86,1,1);

    TCanvas *result = canvas->CanvasWithRatio(520,360, 4, drawOdj, 0.07);
    result->SaveAs(Form("./pict/PictDiffParAntVsCent_%i.png",energy));

    delete result;
    delete canvas;

    fileIn->Close();
}