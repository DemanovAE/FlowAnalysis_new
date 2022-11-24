#ifndef QVECTOR_H
#define QVECTOR_H
#include <iostream>
#include <TMath.h>
#include "TVector2.h"

class QVector
{
public:
  QVector();
  QVector(Double_t n);
  virtual ~QVector();
  void Zero();
  void CalQVector(const Double_t &phi, const Double_t &weight);
  void SetMeanQVector(const Double_t &MeanQx,const Double_t &MeanQy);
  void RecentreningQVector();
  void WeightQVector();
  void MultQVector();
  
  Double_t X() const { return this->fQvector.X(); }
  Double_t Y() const { return this->fQvector.Y(); }
  int GetMult() const { return this->fMult; }
  Double_t GetWeight() const { return this->fWeight; }

  Double_t GetMod() const { return this->fQvector.Mod(); }
  Double_t GetPsi() const { return this->fQvector.Phi(); }
  TVector2 GetQvector() const { return this->fQvector; }
  TVector2 GetMeanQvector() const { return this->fMeanQVector; }


private:
  Double_t fNHarmonic;
  Double_t fWeight;
  TVector2 fQvector;
  TVector2 fMeanQVector;
  int fMult;
  ClassDef(QVector,0);
};

#endif