#ifndef PIDFLOW_H
#define PIDFLOW_H
// C++ headers
#include <vector>
#include <iostream>
#include <fstream> 
#include <map>
#include <string>

// ROOT headers
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TF2.h"
#include "TF1.h"

class PIDflow
{

	public:

		PIDflow();
		
		PIDflow(const Int_t &Energy, const Double_t pt_left, const Double_t pt_right, const Bool_t fixCuts, const Bool_t nSigmaM2, const Bool_t CombPID);
		TF1 *fitFun1DGaus(TH1D *histo, Double_t cons, Double_t mean, Double_t meanL, Double_t meanR, Double_t sigma, Double_t min, Double_t max, Int_t color, Int_t drawFlag = 0);
		TF1 *fitFunThree1DGaus(TH1D *histo,TF1 *FirstGaus, TF1 *SecondGaus, TF1 *ThirdGaus, Double_t min, Double_t max);
		TF1 *fitFunThree2x1DGaus(TH1D *histo,TF1 *Three1DGaus, Double_t min, Double_t max );
		void FitSquareOfMassInPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto = "h2_m2VsPt_all_ch0", Int_t charge = 0);
		void SetInitialValForFit_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto = "h2_m2VsPt_all_ch0", Int_t charge = 0);
		void SetInitialValForFit_NewXY_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto = "h2_m2VsPt_all_ch0", Int_t charge = 0);
		void Fit_nSigmaM2_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto = "h2_m2VsPt_all_ch0", Int_t charge = 0);
		void Fit_NewXY_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto = "h2_m2VsPt_all_ch0", Int_t charge = 0);
		void Fit_TH2(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto = "h2_m2VsPt_all_ch0", Int_t charge = 0, const Char_t *prefix ="", Double_t AxisXMin = 0., Double_t AxisXMax = 0., Double_t AxisYMin = 0., Double_t AxisYMax = 0.);


		Double_t GetParTF2FromFile(TFile *inFile, const Char_t *NameFitFun = "fitFunct", Int_t numberPar = 0);
		Double_t GetParErrTF2FromFile(TFile *inFile, const Char_t *NameFitFun = "fitFunct", Int_t numberPar = 0);


		virtual ~PIDflow(){ /* empty */};

		TFile *InputFile;
		TFile *OutputFile;

		Double_t ptRight;
		Double_t ptLeft;

		TCanvas *canvas;

		TF1 *tf1_gausParticle[3];
		TF1 *tf1_Three1DGaus;
		TF1 *tf1_Three2x1DGaus;

		TH1D *th1_projX;		
		TH1D *th1_projY;

		TH1D *th1_projX_dEdx_M2;		
		TH1D *th1_projY_dEdx_M2;

		Double_t tf2_init_vals[2][3][3];  // x - nSigma,y-m2 // const, mean, sigma // particle
    	Double_t tf2_init_vals_range[2][3][3][2]; // x,y // const, mean, sigma // particle // min,max

    	TF2 *Three2DGaus;
		TF2 *Three2x2DGaus;

		TCanvas *GetCanvas(){return canvas;};


	private:

		//Масса и квадрат массы частиц
		const Float_t electron_mass = 0.0005485799;
		const Float_t pion_mass = 0.13957061;
		const Float_t kaon_mass = 0.493677;
		const Float_t proton_mass = 0.9382720813;

		Int_t Energy_;
		Bool_t fixCuts_;
		Bool_t nSigmaM2_;
		Bool_t CombPID_;

		ClassDef(PIDflow,0)

};

#endif