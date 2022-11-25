#include <MyDrawObject.h>

ClassImp(MyDrawObject);

MyDrawObject::~MyDrawObject(){}

MyDrawObject::MyDrawObject():
  fKeyTF1(kFALSE),
  fKeyTH1D(kFALSE),
  fKeyTGraphErrors(kFALSE),
  fKeyLegend(kFALSE),
  fAxisGraphX(),
  fAxisGraphY(),
  fAxisGraphXerror(),
  fAxisGraphYerror(),
  fAxisGraphSysYerror(),
  fGraphErrors(nullptr),
  fGraphErrorsSystematic(nullptr),
  fFitFunction(nullptr),
  fHisto(nullptr),
  fNameGraph(""),
  fLegend(""),
  fMarkerColor(-1),
  fMarkerStyle(-1),
  fLineWidth(-1),
  fMarkerSize(-1),
  fNumberPadDraw(-1)
{  
}

void MyDrawObject::ClearGraph(){
  fAxisGraphX.clear();
  fAxisGraphY.clear();
  fAxisGraphXerror.clear();
  fAxisGraphYerror.clear();
  fAxisGraphSysYerror.clear();
}

void MyDrawObject::SetDrawOdject(TGraphErrors *tge)
{
  fGraphErrors = tge;
  if(!fGraphErrors){
    std::cout<<"Errors!!! Not Graph"<<std::endl;
    return;
  }
  this->ClearGraph();
  for(int i=0; i<(int)tge->GetN(); i++){
    fAxisGraphX.push_back( tge->GetPointX(i) );
    fAxisGraphY.push_back( tge->GetPointY(i) );
    fAxisGraphXerror.push_back( tge->GetErrorX(i) );
    fAxisGraphYerror.push_back( tge->GetErrorY(i) );
  }
  fKeyTGraphErrors = true;
  MyDrawObject::SetParametrs();
}

void MyDrawObject::SetParametrs()
{
  if(fNameGraph!="")fGraphErrors->SetName(fNameGraph);
  if(fMarkerColor!=-1)fGraphErrors->SetMarkerColor(fMarkerColor);
  if(fMarkerColor!=-1)fGraphErrors->SetLineColor(fMarkerColor);
  if(fMarkerSize!=-1)fGraphErrors->SetMarkerSize(fMarkerSize);
  if(fMarkerStyle!=-1)fGraphErrors->SetMarkerStyle(fMarkerStyle);
  if(fLineWidth!=-1)fGraphErrors->SetLineWidth(fLineWidth);
}


void MyDrawObject::SetPointGraph(double x, double y, double ex, double ey)
{
  fAxisGraphX.push_back(x);
  fAxisGraphY.push_back(y);
  fAxisGraphXerror.push_back(ex);
  fAxisGraphYerror.push_back(ey);

}

void MyDrawObject::SetPointGraph(int i,double x, double y, double ex, double ey)
{
  fAxisGraphX[i]      = x;
  fAxisGraphY[i]      = y;
  fAxisGraphXerror[i] = ex;
  fAxisGraphYerror[i] = ey;
}

void MyDrawObject::SetVectorGraph(std::vector<double> x, std::vector<double> y, std::vector<double> ex, std::vector<double> ey)
{
  fAxisGraphX      = x;
  fAxisGraphY      = y;
  fAxisGraphXerror = ex;
  fAxisGraphYerror = ey;
  fGraphErrors = new TGraphErrors((int)fAxisGraphX.size(),&fAxisGraphX[0],&fAxisGraphY[0],&fAxisGraphXerror[0],&fAxisGraphYerror[0]);
  fKeyTGraphErrors = true;
  MyDrawObject::SetParametrs();
}

void MyDrawObject::GetPointDrawObject(int i, double &x, double &y, double &ex, double &ey)
{
  x  = fAxisGraphX[i];
  y  = fAxisGraphY[i];
  ex = fAxisGraphXerror[i];
  ey = fAxisGraphYerror[i];
}