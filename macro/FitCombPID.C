// C++ headers
#include <vector>
#include <iostream>
#include <fstream> 
#include <map>
#include <string>

// ROOT headers
#include <TStopwatch.h>
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
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
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TVector2.h"
#include "TF2.h"
#include "TF1.h"

#include "PIDflow.h"
#include "ConstantsNew.h"

void TGraphParametrsVsPt(TFile *outFile, const Char_t *NamePID, Int_t charge, const Char_t *prefix, const Char_t *pathToFile);
void DrawParametrsVsPt(TFile *inFile, const Char_t *NamePID, Int_t charge, const Char_t *prefix, Int_t energy);


void FitCombPID(const Char_t *inFileNameWithHisto = "iTest.root", const Char_t *inFileNameWithFitFun = "iTest1.root", const Char_t *outFileName = "oTest2.root", const Int_t energy = 39, Int_t step = 30, const Int_t pt=2){

	TStopwatch timer1;
  	timer1.Start();

  	TFile *inFileWithHisto;
  	TFile *inFileWithFitFun;
  	TFile *outFile;
  	TFile *inFileAfterFirstFit;

  	inFileWithHisto = new TFile(inFileNameWithHisto, "READ");
  	outFile = new TFile(outFileName, "RECREATE");

  	PIDflow *PIDInPtBin[2];

  	if(step==1){
		PIDInPtBin[0] = new PIDflow(energy,ptBinRange[pt],ptBinRange[pt+1],true,false,false);
		PIDInPtBin[0]->SetInitialValForFit_InPtBin(inFileWithHisto, outFile, Form("h2_m2VsnSigma_piKp_%i_pt%i",0,pt),0);
		PIDInPtBin[0]->Fit_nSigmaM2_InPtBin(inFileWithHisto, outFile, Form("h2_m2VsnSigma_piKp_%i_pt%i",0,pt),0);
		
		PIDInPtBin[1] = new PIDflow(energy,ptBinRange[pt],ptBinRange[pt+1],true,false,false);
		PIDInPtBin[1]->SetInitialValForFit_InPtBin(inFileWithHisto, outFile, Form("h2_m2VsnSigma_piKp_%i_pt%i",1,pt),1);
		PIDInPtBin[1]->Fit_nSigmaM2_InPtBin(inFileWithHisto, outFile, Form("h2_m2VsnSigma_piKp_%i_pt%i",1,pt),1);
	}



	if(step==12){
		TString pathToFitFile="/mnt/pool/rhic/1/demanov/cherenkov/NewSTAR/BES/OUT_new/39GeV/HistoM2_39GeV_1_361282/root";
		TGraphParametrsVsPt(outFile, "Pion", 0, "old", pathToFitFile.Data());
		TGraphParametrsVsPt(outFile, "Pion", 1, "old", pathToFitFile.Data());
		TGraphParametrsVsPt(outFile, "Kaon", 0, "old", pathToFitFile.Data());
		TGraphParametrsVsPt(outFile, "Kaon", 1, "old", pathToFitFile.Data());
		TGraphParametrsVsPt(outFile, "Proton", 0, "old", pathToFitFile.Data());
		TGraphParametrsVsPt(outFile, "Proton", 1, "old", pathToFitFile.Data());
		
	}

	if(step==13){
 	 	inFileWithFitFun = new TFile(inFileNameWithFitFun, "READ");
		DrawParametrsVsPt(inFileWithFitFun, "Proton", 1, "old",energy);
		DrawParametrsVsPt(inFileWithFitFun, "Pion", 1, "old",energy);
		DrawParametrsVsPt(inFileWithFitFun, "Kaon", 1, "old",energy);

		DrawParametrsVsPt(inFileWithFitFun, "Proton", 0, "old",energy);
		DrawParametrsVsPt(inFileWithFitFun, "Pion", 0, "old",energy);
		DrawParametrsVsPt(inFileWithFitFun, "Kaon", 0, "old",energy);

		inFileWithFitFun->Close();

	}

	std::cout<<"goos\n";

	inFileWithHisto->Close();
  	outFile->Close();

  	timer1.Stop();
  	timer1.Print();

}

void TGraphParametrsVsPt(TFile *outFile, const Char_t *NamePID, Int_t charge, const Char_t *prefix, const Char_t *pathToFile){

	TFile *inFile[100];
	TF2 *tf2_fit[100];

	Double_t xpt[100]={0.};
	Double_t xpt_err[100]={0.};
	Double_t con[2][100]={0.};
	Double_t mean[2][2][100]={0.};
	Double_t sigma[2][2][100]={0.}; // [X,Y][val, err][point]

	TString NameGraph[] = {	Form("gr_%s_const_%s_ch%i",prefix,NamePID,charge),
							Form("gr_%s_meanNSigma_%s_ch%i",prefix,NamePID,charge),
							Form("gr_%s_sigmaNSigma_%s_ch%i",prefix,NamePID,charge),
							Form("gr_%s_meanM2_%s_ch%i",prefix,NamePID,charge),
							Form("gr_%s_sigmaM2_%s_ch%i",prefix,NamePID,charge)};

	TString txt_val[]={"#mu","#sigma"};
	TString txt_PID[]={"#pi^{+}","K^{+}","p","#pi^{-}","K^{-}","#bar{p}"};
	TString txt_XYcoord[]={"n#sigma_{#pi^{+}}","n#sigma_{#pi^{-}}","m^{2}"};
	TString txt_Title;

	TGraphErrors *gr_point[5];

	PIDflow *PIDInPtBin = new PIDflow();

	Int_t PID=-1;
	if( strncmp(NamePID, "Pion",4) == 0){ PID=0;}
	if( strncmp(NamePID, "Kaon",4) == 0){ PID=1;}
	if( strncmp(NamePID, "Proton",6) == 0){ PID=2;}

	std::cout<<Form("Start write on file: %s ch %i \n",NamePID, charge)<<std::endl;

	for(Int_t i=0; i<(int)ptBinRange.size()-1; i++){
		inFile[i] = new TFile(Form("%s/JOB_pt%i.root",pathToFile,i),"read");
		
		xpt[i] = (ptBinRange[i]+ptBinRange[i+1])/2.0;
		con[0][i] = PIDInPtBin->GetParTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 0+5*PID);
		mean[0][0][i] = PIDInPtBin->GetParTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 1+5*PID);
		sigma[0][0][i] = PIDInPtBin->GetParTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 2+5*PID);
		mean[1][0][i] = PIDInPtBin->GetParTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 3+5*PID);
		sigma[1][0][i] = PIDInPtBin->GetParTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 4+5*PID);

		con[1][i] = PIDInPtBin->GetParErrTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 0+5*PID);
		mean[0][1][i] = PIDInPtBin->GetParErrTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 1+5*PID);
		sigma[0][1][i] = PIDInPtBin->GetParErrTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 2+5*PID);
		mean[1][1][i] = PIDInPtBin->GetParErrTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 3+5*PID);
		sigma[1][1][i] = PIDInPtBin->GetParErrTF2FromFile(inFile[i], Form("Three2x2DGaus_pt%i%i_ch%i",(int)(ptBinRange[i]*10),(int)(ptBinRange[i+1]*10),charge), 4+5*PID);

		std::cout<<"Work pt bin:\t ["<<ptBinRange[i]<<"\t, "<< ptBinRange[i+1]<<"]\n";

	}

	gr_point[0] = new TGraphErrors((int)ptBinRange.size()-1,xpt,con[0],xpt_err,con[1]);
	gr_point[1] = new TGraphErrors((int)ptBinRange.size()-1,xpt,mean[0][0],xpt_err,mean[0][1]);
	gr_point[2] = new TGraphErrors((int)ptBinRange.size()-1,xpt,sigma[0][0],xpt_err,sigma[0][1]);
	gr_point[3] = new TGraphErrors((int)ptBinRange.size()-1,xpt,mean[1][0],xpt_err,mean[1][1]);
	gr_point[4] = new TGraphErrors((int)ptBinRange.size()-1,xpt,sigma[1][0],xpt_err,sigma[1][1]);

	gr_point[1]->SetTitle(Form("%s_{%s}(%s)",txt_val[0].Data(),txt_PID[PID+charge*3].Data(),txt_XYcoord[0+charge].Data()));
	gr_point[2]->SetTitle(Form("%s_{%s}(%s)",txt_val[1].Data(),txt_PID[PID+charge*3].Data(),txt_XYcoord[0+charge].Data()));
	gr_point[3]->SetTitle(Form("%s_{%s}(%s)",txt_val[0].Data(),txt_PID[PID+charge*3].Data(),txt_XYcoord[2].Data()));
	gr_point[4]->SetTitle(Form("%s_{%s}(%s)",txt_val[1].Data(),txt_PID[PID+charge*3].Data(),txt_XYcoord[2].Data()));

	outFile->cd();
	for(Int_t i=0; i<5; i++){
		gr_point[i]->SetName(Form("%s",NameGraph[i].Data()));
		gr_point[i]->Write();
	}

	for(Int_t i=0; i<(int)ptBinRange.size()-1; i++){
		inFile[i] -> Close();
	}
}

void DrawParametrsVsPt(TFile *inFile, const Char_t *NamePID, Int_t charge, const Char_t *prefix, Int_t energy){

	gROOT->SetStyle("Pub");

	TF1 *fit[5];
	TF1 *fit_draw[5];

	TString NameGraph[] = {	Form("gr_%s_const_%s_ch%i",prefix,NamePID,charge),
						Form("gr_%s_meanNSigma_%s_ch%i",prefix,NamePID,charge),
						Form("gr_%s_sigmaNSigma_%s_ch%i",prefix,NamePID,charge),
						Form("gr_%s_meanM2_%s_ch%i",prefix,NamePID,charge),
						Form("gr_%s_sigmaM2_%s_ch%i",prefix,NamePID,charge)};

	TString fitFun[] = {"[0]*TMath::Power(1.0/x,[1]) +[2] + [3]*x + [4]*x*x + [5]*x*x*x + [6]*x*x*x*x",
						"[0]*TMath::Power(1.0/x,[1]) +[2] + [3]*x + [4]*x*x + [5]*x*x*x + [6]*x*x*x*x",
						"[0]*TMath::Power(1.0/x,[1]) +[2] + [3]*x + [4]*x*x + [5]*x*x*x + [6]*x*x*x*x",
						"[0]*TMath::Power(1.0/x,[1]) +[2] + [3]*x + [4]*x*x + [5]*x*x*x + [6]*x*x*x*x",
						"[0]*TMath::Power(1.0/x,[1]) +[2] + [3]*x + [4]*x*x + [5]*x*x*x + [6]*x*x*x*x"};
	//Int_t nPar[5]={6,6,4,6,4};
	Int_t nPar[5]={7,7,7,7,7};
//[0] + [1]*x + [2]*x*x + [3]*x*x*x
	TGraphErrors *gr_draw[5];
	for(Int_t i=0; i<5; i++){
		gr_draw[i] = (TGraphErrors*)inFile->Get(Form("%s",NameGraph[i].Data()));
		
		fit[i]=new TF1(Form("fit_%s",NameGraph[i].Data()),fitFun[i].Data(),0.2,2.8);
		fit[i]->SetLineColor(kRed);
		//gr_draw[i]->Fit(fit[i],"RM");

		fit_draw[i]=new TF1(Form("fit_dr_%s",NameGraph[i].Data()),fitFun[i].Data(),0.2,3.2);
		for(Int_t jp=0; jp<nPar[i];jp++){
			fit_draw[i]->SetParameter(jp,fit[i]->GetParameter(jp));
		}
		fit_draw[i]->SetLineColor(kRed);


	}

	TCanvas *canvas = new TCanvas(Form("can_%s_%i",NamePID,charge),"",1260,960);
	canvas->Divide(2,2);

	for(Int_t i=0; i<4;i++){
	    canvas->cd(i+1);

    	gr_draw[i+1]->SetMarkerColor(kBlue);
    	gr_draw[i+1]->SetMarkerStyle(20);
    	gr_draw[i+1]->SetMarkerSize(1);
    	gr_draw[i+1]->SetLineWidth(2);
    	gr_draw[i+1]->GetYaxis()->SetTitle(gr_draw[i+1]->GetTitle());
    	gr_draw[i+1]->GetXaxis()->SetTitle("p_{T} [GeV/c]");
	    gr_draw[i+1]->Draw("AP");
	    //fit_draw[i+1]->Draw("same");


	}

	canvas->SaveAs(Form("%ican_%s_ch%i.png",energy,NamePID,charge));

}