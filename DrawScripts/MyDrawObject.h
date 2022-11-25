#ifndef MYDRAWOBJECT_H
#define MYDRAWOBJECT_H

#include <string>
#include <iostream>
#include <vector>
#include <TFile.h>
#include <TMath.h>
#include <TVector3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TString.h>
#include <TDirectoryFile.h>
#include <TGraphErrors.h>
#include <TF1.h>

class MyDrawObject
{
public:
  MyDrawObject();
  // Деструктор
  virtual ~MyDrawObject();

  void SetNameGraph(TString str){ this->fNameGraph = str; }
  void SetLegendText(TString str){ this->fLegend = str; fKeyLegend=true; }
  void SetMarkerColor(int i){ this->fMarkerColor = i; }
  void SetMarkerStyle(int i){ this->fMarkerStyle = i; }
  void SetMarkerSize(float f){ this->fMarkerSize = f; }
  void SetLineWidth(float f){ this->fLineWidth = f; }
  void SetNumberPadDraw(int i){ this->fNumberPadDraw = i; }
  
  void SetDrawOdject(TF1 *tf1){ this->fFitFunction = tf1; fKeyTF1 = true; }
  void SetDrawOdject(TH1D *th1){ this->fHisto = th1; fKeyTH1D = true; }
  void SetDrawOdject(TGraphErrors *tge);

  void SetPointGraph(double x, double y, double ex, double ey);
  void SetPointGraph(int i, double x, double y, double ex, double ey);
  void SetVectorGraph(std::vector<double> x, std::vector<double> y, std::vector<double> ex, std::vector<double> ey);
  void ChangePointGraph(int i, double x, double y, double ex, double ey);
  void ClearGraph();
  void SetParametrs();


  int GetNumberPadDraw(){ return this->fNumberPadDraw; }
  int GetSizeVector(){ return this->fAxisGraphX.size(); }
  TF1* GetDrawObjectTF1(){ return this->fFitFunction; }
  TH1D* GetDrawObjectTH1D(){ return this->fHisto; }
  TGraphErrors* GetDrawObjectGraph(){ return this->fGraphErrors; }
  std::vector<double> GetPointGraph(int i){ return ( (int)fAxisGraphX.size()<i ) ? std::vector<double>{-999., -999., -999., -999.} : std::vector<double>{this->fAxisGraphX[i],this->fAxisGraphY[i],this->fAxisGraphXerror[i],this->fAxisGraphYerror[i]}; }
  void GetPointDrawObject(int i, double &x, double &y, double &ex, double &ey);

  TString GetLegendText(){ return this->fLegend; }
  TString GetNameGraph(){ return this->fNameGraph; }
  Bool_t GetKeyGraph(){ return this->fKeyTGraphErrors; }
  Bool_t GetKeyLegend(){ return this->fKeyLegend; }
  Bool_t GetKeyFun(){ return this->fKeyTF1; }
  Bool_t GetKeyTH(){ return this->fKeyTH1D; }

private:

  Bool_t fKeyTF1;
  Bool_t fKeyTH1D;
  Bool_t fKeyTGraphErrors;
  Bool_t fKeyLegend;

  std::vector<double> fAxisGraphX;
  std::vector<double> fAxisGraphY;
  std::vector<double> fAxisGraphXerror;
  std::vector<double> fAxisGraphYerror;
  std::vector<double> fAxisGraphSysYerror;

  TGraphErrors *fGraphErrors;
  TGraphErrors *fGraphErrorsSystematic;
  TF1 *fFitFunction;
  TH1D *fHisto;

  TString fNameGraph;
  TString fLegend;

  int fMarkerColor;
  int fMarkerStyle;
  float fLineWidth;
  float fMarkerSize;
  int fNumberPadDraw;

  ClassDef(MyDrawObject,0);
};

#endif