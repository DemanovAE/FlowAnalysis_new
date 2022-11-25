#ifndef ANALYSISFUNCTION_H
#define ANALYSISFUNCTION_H
#include <string>
#include <iostream>
#include <vector>
#include <TGraph.h>
#include <Math/SpecFuncMathMore.h>
#include <TMath.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TProfile3D.h>
#include "TFile.h"
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>
#include <algorithm>
#include <TGraphErrors.h>

class AnalysisFunction
{
public:
	AnalysisFunction();
    // Деструктор
    ~AnalysisFunction();
	void ClearVecGraph();
    void SetScale(double d){this->fScale=d;};
    void SetKeyAnaliticalRes(bool bl, Int_t ii, Int_t io);
    void SetKeyСomponentSP(){ this->fKeyСomponentSP=true; }
	double GetRes(double _chi, double _harm);
	double GetChi(double _res, double _harm, int accuracy);
	double ErrorYRatio(double num, double den, double num_err, double den_err){	return TMath::Sqrt( TMath::Power( num_err/den, 2) + TMath::Power( (num*den_err) / (den*den) , 2) ); }

	TGraphErrors* ResolutionEPCent9_Analitic(TFile *InputFile, TString InputResolutionTPName, Int_t InputHarm, TString OutputResolutionTPName, Int_t OutputHarm);	
	TGraphErrors* ResolutionEP_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString OutputResolutionTPName);
	
	TGraphErrors* CalculateFlowVsPtEtaSub(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);
	TGraphErrors* CalculateFlowVsPtScalarProduct(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);
	TGraphErrors* CalculateFlowVsRapidity(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector);
	TGraphErrors* CalculateFlowVsRapidity2(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);


	TGraphErrors* CalculateFlowVsCentEtaSub(TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);
	TGraphErrors* CalculateFlowVsCentEtaSub2(TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);
	TGraphErrors* CalculateFlowVsCentScalarProduct(TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);

	TGraphErrors* FlowVsPt_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);
	TGraphErrors* FlowVsPt_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);
	TGraphErrors* FlowVsPt_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);
	TGraphErrors* FlowVsPt_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector);
	
	TGraphErrors* FlowVsRapidity_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, Double_t PtMin, Double_t PtMax, std::vector<double> rebin_vector);
	TGraphErrors* FlowVsRapidity_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, Double_t PtMin, Double_t PtMax, std::vector<double> rebin_vector);
	TGraphErrors* FlowVsRapidity_ScalarProduct2(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, std::vector<double> rebin_vector);
	TGraphErrors* FlowVsRapidity_ScalarProduct2(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, std::vector<double> rebin_vector);

	TGraphErrors* FlowVsCent_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);
	TGraphErrors* FlowVsCent_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);
	TGraphErrors* FlowVsCent_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);
	TGraphErrors* FlowVsCent_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector);

	TGraphErrors* RatioGraphPointToPoint(TGraphErrors *numerator, TGraphErrors *denominator, TString OutputGraphName);
	TGraphErrors* RatioGraphEval(TGraphErrors *numerator, TGraphErrors *denominator, TString OutputGraphName);
	TGraphErrors* DifferenceGraphPointToPoint(TGraphErrors *numerator, TGraphErrors *denominator, TString OutputGraphName);
	TGraphErrors* RatioGraphToOnePoint(TGraphErrors *numerator, TGraphErrors *denominator, Int_t point, TString OutputGraphName);


	//TGraphErrors* FlowVsPt_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector );

private:
	
	double fScale;
	Bool_t fKeyСomponentSP;
	Bool_t fKeyAnaliticalRes;
	Double_t fInputHarm;
	Double_t fOutputHarm;

	TFile *fInputFile;
	TFile *fOutputFile;

	TProfile *fResolutionTP;
	TProfile *fFlowVsCentTP;
	TProfile2D *fFlowVsPtVsCentTP2D;
	TProfile3D *fFlowVsPtVsCentTP3D;
	TProfile2D *fMeanPtVsPtVsCentTP2D;

	std::vector<Double_t> fGraphX;
	std::vector<Double_t> fGraphXErrors;
	std::vector<Double_t> fGraphXErrorsSys;
	std::vector<Double_t> fGraphY;
	std::vector<Double_t> fGraphYErrors;
	std::vector<Double_t> fGraphYErrorsSys;

	ClassDef(AnalysisFunction,0);
};

#endif