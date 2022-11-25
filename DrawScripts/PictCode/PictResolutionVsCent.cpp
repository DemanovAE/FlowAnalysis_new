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

void PictResolutionVsCent(const Char_t *inFileName,const Int_t energy = 28){
	
    TFile *fstyle = new TFile("style.root");
    TStyle *tsty = (TStyle *)fstyle->Get("PlottingInStyle");
    tsty->cd();

    gROOT->ForceStyle();
    gStyle->SetOptStat(0);

	std::map<int, string> Energy_scan_text={{7, "7.7"},{11, "11.5"},{14, "14.5"},{19, "19.6"},{27, "27"},{39, "39"},{62, "62.4"}};
	TString particle[]={"h^{+}","h^{-}","#pi^{+}","#pi^{-}","K^{+}","K^{-}","p","#bar{p}","#pi^{+}#pi^{-}","K^{+}K^{-}","p#bar{p}","Hadrons","","PID"};
  
    TFile *fileIn = new TFile(inFileName,"read");
    //TFile *fileIn = new TFile(Form("/home/aleksandr/STAR_BES/result/root_sys_27GeV/flow_28GeV_PID.root"),"READ");

    Int_t color[]={2, 1, 4, 6, 2, 46, 2, 1};
    Int_t style[]={23, 26, 22, 32, 20, 24, 34, 28, 29, 30};
    Float_t size[]={2, 2, 2, 2, 3, 2};
  

    std::vector<MyDrawObject *> drawOdj;
    std::vector<AnalysisFunction *> AnalysOdj;
    std::vector<double> minPt={0.2,0.2,0.2,0.5};
    std::vector<double> maxPt={3.2,2.4,2.4,3.2};
    std::vector<int> fNHarmonic ={2,3,4,2,4};
    std::vector<TString> fPidCode = {"HadronsPos","HadronsNeg","PionPos","PionNeg","KoanPos","KaonNeg","ProtonPos","ProtonNeg"};
    //std::vector<TString> fPidCode = {"PionPos","PionNeg","KaonPos","KaonNeg","ProtonPos","ProtonNeg"};
    std::vector<float> fEtaGap = {0.075};
    std::vector<TString> fNameSys = {"","","","v42",""};

    std::vector<std::vector<double>> CentBinFemtoDst = {
        {-0.5,0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,1.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,3.5,4.5,5.5,6.5,7.5,8.5},
        {-0.5,3.5,4.5,5.5,6.5,7.5,8.5}
    };
    std::vector<std::vector<double>> CentBin         = {
        { 80.,70.,60.,50.,40.,30.,20.,10.,5.,0.},
        { 80.,60.,40.,30.,20.,10.,5.,0.},
        { 80.,40.,30.,20.,10.,5.,0.},
        { 80.,40.,30.,20.,10.,5.,0.}
    };


    for(Int_t ih=0; ih<fNHarmonic.size(); ih++){    
        
        AnalysOdj.push_back(new AnalysisFunction());
        if(ih==4)AnalysOdj.back()->SetScale(sqrt(2.));
        //if(ih==3)AnalysOdj.back()->SetScale(1./2.);
        drawOdj.push_back(new MyDrawObject());
        drawOdj.back()->SetMarkerColor(color[ih]);
        drawOdj.back()->SetMarkerStyle(style[2*ih]);
        drawOdj.back()->SetMarkerSize(size[ih]);
        drawOdj.back()->SetNumberPadDraw(0);
        drawOdj.back()->SetLegendText(Form("#font[42]{#scale[1.0]{R(#Psi_{%i}): #sqrt{<cos(%i*(#Psi_{%i}^{E}-#Psi_{%i}^{W}))>} }}",fNHarmonic[ih],fNHarmonic[ih],fNHarmonic[ih],fNHarmonic[ih]));
        if(ih==3)drawOdj.back()->SetLegendText(Form("#font[42]{#scale[1.0]{R(#Psi_{4}): #sqrt{<cos(4*(#Psi_{2}^{E}-#Psi_{2}^{W}))>} }}"));
        if(ih==4)drawOdj.back()->SetLegendText(Form("#font[42]{#scale[1.0]{R(#Psi_{4}): estimate, R(#chi)}}"));
        if(ih<4)drawOdj.back()->SetDrawOdject( AnalysOdj.back()->ResolutionEP_EtaSub(fileIn,
                Form("tp_SqRes%iTPC%i%s",fNHarmonic[ih],(int)(fEtaGap[0]*100),fNameSys[ih].Data()),
                Form("gt_Res%i_cent%s",fNHarmonic[ih],fNameSys[ih].Data())
                ));
        if(ih==4)drawOdj.back()->SetDrawOdject( AnalysOdj.back()->ResolutionEPCent9_Analitic(fileIn,
                Form("tp_SqRes%iTPC%i%s",2,(int)(fEtaGap[0]*100),fNameSys[ih].Data()),2,
                Form("gt_Res%i_cent%s",fNHarmonic[ih],fNameSys[ih].Data()),4
                ));
    }

    std::cout<<energy<<std::endl;

    Draw_Picture_new *canvas = new Draw_Picture_new();
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ Resolution }}"),"Y", 1,0.9,0.65,0.3,90.0);
    canvas->SetTextInfo(Form("#font[42]{#scale[1.0]{ centtrality}}"),"X", 1,0.5,0.8,0.3,0.0);
    canvas->SetTextInfo(Form("#font[42]{ #scale[1.0]{Au+Au, #sqrt{S_{NN}}=27 GeV Run18}}" ),"Pad", 0, 1, 0.60, 0.07, 0.0);    

    canvas->SetTLineOnPad(-1, 0., 85, 0., 1, 3, 7, -1);
    
    canvas->SetAxisToPad(-1, 85, 0.05, -0.01, 0.66, 0.05);
    canvas->SetLegendPar(0.55,0.65,0.96,0.96,1,0);

    TCanvas *result = canvas->CanvasNxM(820,660, 1,1, drawOdj, 0.05);
    result->SaveAs(Form("./PictResolutionVsCent2_%i.png",energy));

    delete result;
    delete canvas;

    fileIn->Close();
}