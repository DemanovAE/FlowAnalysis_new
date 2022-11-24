#ifndef MYDRAWOBJECT_H
#define MYDRAWOBJECT_H

#include <string>
#include <iostream>
#include <vector>
#include <TGraph.h>
#include <Math/SpecFuncMathMore.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TPad.h>

class MyDrawObject
{
public:
  MyDrawObject();
  // Деструктор
  virtual ~MyDrawObject();
  void ClearObject();

private:

  std::vector<double> axisX;
  std::vector<double> axisY;
  std::vector<double> axisXerror;
  std::vector<double> axisYerror;
  std::vector<double> axisSysYerror;

  ClassDef(MyDrawObject,0);
};

#endif