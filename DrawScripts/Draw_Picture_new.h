#ifndef CANVASPLOT_H
#define CANVASPLOT_H

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "TLegend.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "MyDrawObject.h"
#include "TLine.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TPad.h"
#include "TGraphErrors.h"
#include "TFrame.h"
#include <TF1.h>

class Draw_Picture_new
{
public:
	Draw_Picture_new();
    // Деструктор
    ~Draw_Picture_new();
    
    struct TextInfo
    {
    	TString text;
    	TString place; //X , Y , Pad, Info
    	int NumberPad;
    	float x;
    	float y;
    	float size;
    	float angle;
    };

   	void SetTLineOnPad(Float_t mXmin, Float_t mYmin, Float_t mXmax, Float_t mYmax, Int_t mColour, Int_t mLineWidth, Int_t mLineStyle, Int_t NumberPad);
    void SetAxisToPad(Double_t Xmin, Double_t Xmax, Double_t XLabelSize, Double_t YminUp, Double_t YmaxUp, Double_t YLabelSize);
	void SetAxisToPad2(Double_t Xmin, Double_t Xmax, Double_t XLabelSize, Double_t YminUp, Double_t YmaxUp, Double_t YLabelSize);
	void SetLegendPar(Double_t xMin, Double_t yMin, Double_t xMax, Double_t yMax, Int_t N, Int_t NumberPad);
	void SetTextInfo(TString text, TString place, Int_t NumberPad, Float_t x, Float_t y, Float_t size, Float_t angle);
	TCanvas *CanvasNxM(Int_t xCanvas, Int_t yCanvas, Int_t line, Int_t column, std::vector<MyDrawObject*> DrawResult, Float_t size_bukva);
    TCanvas *CanvasWithRatio(Int_t xCanvas, Int_t yCanvas, Int_t column, std::vector<MyDrawObject*> DrawResult, Float_t size_bukva);
private:

    std::vector<TextInfo> labels;
    bool Sys_mod=false;
    bool legend_flag=false;

    TFile *fileIn;
    TFile *fileOut;

    TH2F *hPadAxis;
    TH2F *hPadAxis2;
    TH2F *hPadUp;
    TH2F *hPadDown;
    TH2F *hPadOne;
    std::vector<TH2F *> hPadNM;

    TLegend *legend;
    int legend_nPad;
	TLatex _textInfo;
	int type_canvas;

	std::vector<TLine *> Line_vector;
	std::vector<int> Line_vector_DrawPad;

	std::vector<TGraphErrors *> graph_up;
	std::vector<TGraphErrors *> graph_up_sys;
	std::vector<TGraphErrors *> graph_down;

	ClassDef(Draw_Picture_new,0);

};
#endif