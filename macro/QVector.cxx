#include <QVector.h>

ClassImp(QVector);

QVector::QVector() : fNHarmonic(2), fMult(0), fWeight(0), fQvector(0.,0.), fMeanQVector(0.,0.){}

QVector::QVector(Double_t n) : fNHarmonic(n), fMult(0), fWeight(0), fQvector(0.,0.), fMeanQVector(0.,0.){}

QVector::~QVector(){}

void QVector::Zero(){
  fQvector.Set(0.,0.);
  fMeanQVector.Set(0.,0.);
  fMult = 0;
  fWeight = 0.; 
}

void QVector::SetMeanQVector(const Double_t &MeanQx,const Double_t &MeanQy){
  fMeanQVector.Set(MeanQx,MeanQy);
}

void QVector::CalQVector(const Double_t &phi, const Double_t &weight){
  fQvector += weight * TVector2(TMath::Cos( fNHarmonic * phi ), TMath::Sin( fNHarmonic * phi ));
  fWeight += weight;
  fMult++;
}

void QVector::RecentreningQVector(){
  fQvector = fQvector - fMeanQVector;
}

void QVector::WeightQVector(){
  if (fMult == 0) { std::cout << "Warning! fMult==0!" << std::endl;}
  if (fMult != 0)
  {
    if (fWeight<0) {std::cout<<"Warning! fWeight<0!" << std::endl;} /*fWeight *= -1.;*/
    fQvector = fQvector / fWeight - fMeanQVector;
  }
}


void QVector::MultQVector(){
  fWeight = fMult;
  WeightQVector();
}