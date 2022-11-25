#include <AnalysisFunction.h>
ClassImp(AnalysisFunction);
AnalysisFunction::AnalysisFunction():
	fInputHarm(1.),
	fOutputHarm(1.),
	fKeyAnaliticalRes(kFALSE),
	fKeyСomponentSP(kFALSE),
	fScale(1.0),
	fInputFile(nullptr),
	fOutputFile(nullptr),
	fResolutionTP(nullptr),
	fFlowVsCentTP(nullptr),
	fFlowVsPtVsCentTP2D(nullptr),
	fFlowVsPtVsCentTP3D(nullptr),
	fMeanPtVsPtVsCentTP2D(nullptr),
	fGraphX(),
	fGraphXErrors(),
	fGraphXErrorsSys(),
	fGraphY(),
	fGraphYErrors(),
	fGraphYErrorsSys()
{
};
// Деструктор
AnalysisFunction::~AnalysisFunction(){
};

void AnalysisFunction::ClearVecGraph(){
	fGraphX.clear();
	fGraphXErrors.clear();
	fGraphXErrorsSys.clear();
	fGraphY.clear();
	fGraphYErrors.clear();
	fGraphYErrorsSys.clear();
};

void AnalysisFunction::SetKeyAnaliticalRes(bool bl, Int_t ii, Int_t io){
	fKeyAnaliticalRes=bl;
	fInputHarm=ii;
	fOutputHarm=io;
}
// Resolution calculation
//S----------------------------------------------------------------
double AnalysisFunction::GetRes(double _chi, double _harm)
{

    double con = TMath::Sqrt(TMath::Pi() / 2) / 2;
    double arg = _chi * _chi / 4.;
    double order1 = (_harm - 1) / 2.;
    double order2 = (_harm + 1) / 2.;
    double res = con * _chi * exp(-arg) * (ROOT::Math::cyl_bessel_i(order1, arg) + ROOT::Math::cyl_bessel_i(order2, arg));
    return res;
}
// Chi2 calculation
//S----------------------------------------------------------------
double AnalysisFunction::GetChi(double _res, double _harm, int accuracy)
{
    double chi = 2.0;
    double delta = 1;
    for (int i = 0; i < accuracy; i++)
    {
        if (GetRes(chi, _harm) < _res)
            chi = chi + delta;
        else
            chi = chi - delta;
        delta = delta / 2.;
    }
    return chi;
}
//E-----

TGraphErrors* AnalysisFunction::CalculateFlowVsPtEtaSub(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	AnalysisFunction::ClearVecGraph();
	TProfile *tp_vn_bin[9];

	TProfile *tp_pt2 = (TProfile*)fMeanPtVsPtVsCentTP2D->ProfileX(Form("ProfileX_%s",OutputGraphName.Data()),CentBinMin+1,CentBinMax+1);
	TH1 *tp_meanPt = tp_pt2->Rebin((int)rebin_vector.size()-1, Form("hRebin_%s",OutputGraphName.Data()), &rebin_vector[0]);

	//set v_n axis
	for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
		tp_vn_bin[j] = (TProfile*)fFlowVsPtVsCentTP2D -> ProfileX(Form("profX_%s_cent%i%i", OutputGraphName.Data(), j, j),j+1,j+1);
		tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( TMath::Abs(fResolutionTP->GetBinContent(fResolutionTP->FindBin((Double_t)j)))) );
	}
	
	if(CentBinMin != CentBinMax){
		for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
			tp_vn_bin[CentBinMin]->Add(tp_vn_bin[j]);
		}
	}

	TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin((int)rebin_vector.size()-1, Form("RebinTH1_%s",OutputGraphName.Data()), &rebin_vector[0]);

	for(Int_t i=0; i < (int)tp_meanPt->GetNbinsX()-1; i++){
		fGraphX.push_back(tp_meanPt->GetBinContent(i+1));
		//fGraphX.push_back(tp_meanPt->GetBinCenter(i+1));
		fGraphXErrors.push_back(0.);
		fGraphY.push_back(tp_flow->GetBinContent(tp_flow->FindBin(fGraphX[i])));
		fGraphXErrors.push_back(tp_flow->GetBinError(tp_flow->FindBin(fGraphX[i])));
	}

	TGraphErrors *gr = new TGraphErrors((int)fGraphX.size(), &fGraphX[0], &fGraphY[0], &fGraphXErrors[0], &fGraphXErrors[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::CalculateFlowVsPtScalarProduct(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	TProfile *tp_vn_bin[9];

	TProfile *tp_pt2 = (TProfile*)fMeanPtVsPtVsCentTP2D->ProfileX(Form("ProfileX_%s",OutputGraphName.Data()),CentBinMin+1,CentBinMax+1);
	TH1 *tp_meanPt = tp_pt2->Rebin((int)rebin_vector.size()-1, Form("hRebin_%s",OutputGraphName.Data()), &rebin_vector[0]);

	Double_t Resol;

	//set v_n axis
	for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
		tp_vn_bin[j] = (TProfile*)fFlowVsPtVsCentTP2D -> ProfileX(Form("profX_%s_cent%i%i", OutputGraphName.Data(), j, j),j+1,j+1);
		
		if(fKeyAnaliticalRes==false){
			Resol = TMath::Sqrt( TMath::Abs(fResolutionTP->GetBinContent(fResolutionTP->FindBin((Double_t)j)) ));
		}
		if(fKeyAnaliticalRes==true ){
			Resol = GetRes( fScale*GetChi( TMath::Sqrt(TMath::Abs(fResolutionTP -> GetBinContent(fResolutionTP->FindBin((Double_t)j) ))), fInputHarm, 100), fOutputHarm);
		}
		tp_vn_bin[j] -> Scale( 1.0 / (TMath::Sqrt( ((fKeyСomponentSP) ? 2 : 1) ) * Resol) );
		//tp_vn_bin[j] -> Scale( 1.0 /  Resol );
	}
	
	if(CentBinMin != CentBinMax){
		for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
			tp_vn_bin[CentBinMin]->Add(tp_vn_bin[j]);
		}
	}

	TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin((int)rebin_vector.size()-1, Form("RebinTH1_%s",OutputGraphName.Data()), &rebin_vector[0]);

	for(Int_t i=0; i < (int)tp_meanPt->GetNbinsX()-1; i++){
		fGraphX.push_back(tp_meanPt->GetBinContent(i+1));
		//fGraphX.push_back(tp_meanPt->GetBinCenter(i+1));
		fGraphXErrors.push_back(0.);
		fGraphY.push_back(tp_flow->GetBinContent(tp_flow->FindBin(fGraphX[i])));
		fGraphYErrors.push_back(tp_flow->GetBinError(tp_flow->FindBin(fGraphX[i])));
	}

	TGraphErrors *gr = new TGraphErrors((int)fGraphX.size(), &fGraphX[0], &fGraphY[0], &fGraphXErrors[0], &fGraphYErrors[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

/*
TGraphErrors* AnalysisFunction::CalculateFlowVsPtEtaSub2(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	std::vector<Double_t> pt;
	std::vector<Double_t> pt_err;
	std::vector<Double_t> flow;
	std::vector<Double_t> flow_err;

	TProfile *tp_vn_bin[9];

	Double_t centReb[3];
	if(CentBinMin==-0.5 && CentBinMax!=8.5){
		centReb[0]=CentBinMin;
		centReb[1]=CentBinMax;
		centReb[2]=8.5;
	}

	if(CentBinMin!=-0.5){
		centReb[0]=-0.5;
		centReb[1]=CentBinMin;
		centReb[2]=CentBinMax;
	}

	int vn_CentMinBin = fFlowVsPtVsCentTP2D->GetYaxis()->FindBin(CentBinMin);
	int vn_CentMaxBin = fFlowVsPtVsCentTP2D->GetYaxis()->FindBin(CentBinMax);
	int pt_CentMinBin = fMeanPtVsPtVsCentTP2D->GetYaxis()->FindBin(CentBinMin);
	int pt_CentMaxBin = fMeanPtVsPtVsCentTP2D->GetYaxis()->FindBin(CentBinMax);

	TH1 *th_Res = fResolutionTP->Rebin(2, Form("ResRebin_%s",OutputGraphName.Data()), centReb);
	
	TProfile *tp_vn  = (TProfile*)fFlowVsPtVsCentTP2D->ProfileX(Form("ProfileX_%s",OutputGraphName.Data()),vn_CentMinBin,vn_CentMaxBin);
	TH1 *tp_flow = tp_vn->Rebin((int)rebin_vector.size()-1, OutputGraphName.Data(), &rebin_vector[0]);
	tp_flow -> Scale( 1.0 / TMath::Sqrt( TMath::Abs(fResolutionTP->GetBinContent(fResolutionTP->FindBin( (CentBinMin+CentBinMax)/2. )))) );

	TProfile *tp_pt2 = (TProfile*)fMeanPtVsPtVsCentTP2D->ProfileX(Form("ProfileX_%s",OutputGraphName.Data()),pt_CentMinBin,pt_CentMaxBin);
	TH1 *tp_meanPt = tp_pt2->Rebin((int)rebin_vector.size()-1, Form("hRebin_%s",OutputGraphName.Data()), &rebin_vector[0]);

	//set v_n axis
	
	for(Int_t i=0; i < (int)tp_meanPt->GetNbinsX()-1; i++){
		pt.push_back(tp_meanPt->GetBinContent(i+1));
		//pt.push_back(tp_meanPt->GetBinCenter(i+1));
		pt_err.push_back(0.);
		flow.push_back(tp_flow->GetBinContent(tp_flow->FindBin(pt[i])));
		flow_err.push_back(tp_flow->GetBinError(tp_flow->FindBin(pt[i])));
	}

	TGraphErrors *gr = new TGraphErrors((int)pt.size(), &pt[0], &flow[0], &pt_err[0], &flow_err[0]);
	gr->SetName(OutputGraphName);
	return gr;

}
*/
TGraphErrors* AnalysisFunction::CalculateFlowVsCentEtaSub(TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector ){

	std::vector<Double_t> cent;
	std::vector<Double_t> cent_err;
	std::vector<Double_t> flow;
	std::vector<Double_t> flow_err;

	int PtMinBin = fFlowVsPtVsCentTP2D->GetXaxis()->FindBin(PtMin);
	int PtMaxBin = fFlowVsPtVsCentTP2D->GetXaxis()->FindBin(PtMax);

	TProfile *tp_flowNoRebin = (TProfile*)fFlowVsPtVsCentTP2D -> ProfileY(Form("profY_%s_cent_pt_%.2f_%.2f", OutputGraphName.Data(), PtMin, PtMax),PtMinBin,PtMaxBin);
	TH1 *tp_flow = tp_flowNoRebin->Rebin((int)rebin_vector.size()-1, Form("h2_%s",OutputGraphName.Data()), &rebin_vector[0]);
	TH1 *tp_res = 	fResolutionTP->Rebin((int)rebin_vector.size()-1, Form("h2_%s",OutputGraphName.Data()), &rebin_vector[0]);
	
	Double_t Resol = 0.;

	for(Int_t ic=1; ic<=tp_flow->GetNbinsX(); ic++){
		//cent.push_back(tp_flow->GetBinCenter(ic));
		cent.push_back( (cent_vector[ic-1]+cent_vector[ic])/2. );
		cent_err.push_back(0.);
		if(fKeyAnaliticalRes==false) Resol = TMath::Sqrt( TMath::Abs(tp_res->GetBinContent(ic)));
		if(fKeyAnaliticalRes==true ){
			Resol = GetRes( fScale*GetChi( TMath::Sqrt(TMath::Abs(fResolutionTP -> GetBinContent(ic))), fInputHarm, 100), fOutputHarm);
		}
		flow.push_back( tp_flow->GetBinContent(ic) / Resol );
		//flow.push_back(tp_flow->GetBinContent(ic)); // TMath::Sqrt( TMath::Abs(tp_res->GetBinContent(ic))));
		flow_err.push_back(AnalysisFunction::ErrorYRatio(tp_flow->GetBinContent(ic), 
														 Resol, 
														 tp_flow->GetBinError(ic), 
														 TMath::Abs(tp_res->GetBinError(ic)) / Resol ));

	}

	TGraphErrors *gr = new TGraphErrors((int)cent.size(), &cent[0], &flow[0], &cent_err[0], &flow_err[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::CalculateFlowVsCentScalarProduct(TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector ){

	std::vector<Double_t> cent;
	std::vector<Double_t> cent_err;
	std::vector<Double_t> flow;
	std::vector<Double_t> flow_err;

	int PtMinBin = fFlowVsPtVsCentTP2D->GetXaxis()->FindBin(PtMin);
	int PtMaxBin = fFlowVsPtVsCentTP2D->GetXaxis()->FindBin(PtMax);

	TProfile *tp_flowNoRebin = (TProfile*)fFlowVsPtVsCentTP2D -> ProfileY(Form("profY_%s_cent_pt_%.2f_%.2f", OutputGraphName.Data(), PtMin, PtMax),PtMinBin,PtMaxBin);
	TH1 *tp_flow = tp_flowNoRebin->Rebin((int)rebin_vector.size()-1, Form("h2_%s",OutputGraphName.Data()), &rebin_vector[0]);
	TH1 *tp_res = 	fResolutionTP->Rebin((int)rebin_vector.size()-1, Form("h2_%s",OutputGraphName.Data()), &rebin_vector[0]);
	
	Double_t Resol = 0.;

	for(Int_t ic=1; ic<=tp_flow->GetNbinsX(); ic++){
		//cent.push_back(tp_flow->GetBinCenter(ic));
		cent.push_back( (cent_vector[ic-1]+cent_vector[ic])/2. );
		cent_err.push_back(0.);
		if(fKeyAnaliticalRes==false) Resol = TMath::Sqrt(TMath::Abs(tp_res->GetBinContent(ic)));
		if(fKeyAnaliticalRes==true ){
			Resol = GetRes( fScale*GetChi( TMath::Sqrt(TMath::Abs(fResolutionTP -> GetBinContent(ic))), fInputHarm, 100), fOutputHarm);
		}
		Resol = TMath::Sqrt( ((fKeyСomponentSP) ? 0.5 : 1) ) * Resol;
		flow.push_back( tp_flow->GetBinContent(ic) / Resol );
		//flow.push_back(tp_flow->GetBinContent(ic)); // TMath::Sqrt( TMath::Abs(tp_res->GetBinContent(ic))));
		flow_err.push_back(AnalysisFunction::ErrorYRatio(tp_flow->GetBinContent(ic), 
														 Resol, 
														 tp_flow->GetBinError(ic), 
														 TMath::Abs(tp_res->GetBinError(ic)) / Resol ));

	}

	TGraphErrors *gr = new TGraphErrors((int)cent.size(), &cent[0], &flow[0], &cent_err[0], &flow_err[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::CalculateFlowVsRapidity(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector){

	Int_t nBinPtMin = fFlowVsPtVsCentTP3D->GetXaxis()->FindBin(PtMin);
	Int_t nBinPtMax = fFlowVsPtVsCentTP3D->GetXaxis()->FindBin(PtMax);
	fFlowVsPtVsCentTP3D->GetXaxis()->SetRange(nBinPtMin,nBinPtMax);
	
	fFlowVsPtVsCentTP2D = (TProfile2D*)fFlowVsPtVsCentTP3D->Project3DProfile("YZ");
	TProfile *tp_vn_bin[9];
	//set v_n axis
	for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
		tp_vn_bin[j] = (TProfile*)fFlowVsPtVsCentTP2D -> ProfileX(Form("profX_%s_cent%i%i", OutputGraphName.Data(), j, j),j+1,j+1);
		tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( TMath::Abs(fResolutionTP->GetBinContent(fResolutionTP->FindBin((Double_t)j)))) );
	}
	
	if(CentBinMin != CentBinMax){
		for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
			tp_vn_bin[CentBinMin]->Add(tp_vn_bin[j]);
		}
	}

	TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin((int)rebin_vector.size()-1, Form("RebinTH1_%s",OutputGraphName.Data()), &rebin_vector[0]);

	for(Int_t i=0; i < (int)tp_flow->GetNbinsX(); i++){
		fGraphX.push_back(tp_flow->GetBinCenter(i+1));
		//pt.push_back(tp_meanPt->GetBinCenter(i+1));
		fGraphXErrors.push_back(0.);
		fGraphY.push_back(tp_flow->GetBinContent(tp_flow->FindBin(fGraphX[i])));
		fGraphYErrors.push_back(tp_flow->GetBinError(tp_flow->FindBin(fGraphX[i])));
	}

	TGraphErrors *gr = new TGraphErrors((int)fGraphX.size(), &fGraphX[0], &fGraphY[0], &fGraphXErrors[0], &fGraphYErrors[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::CalculateFlowVsRapidity2(TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, std::vector<double> rebin_vector){

	TProfile *tp_vn_bin[9];
	//set v_n axis
	for(Int_t j=CentBinMin;  j<=CentBinMax; j++){
		tp_vn_bin[j] = (TProfile*)fFlowVsPtVsCentTP2D -> ProfileX(Form("profX_%s_cent%i%i", OutputGraphName.Data(), j, j),j+1,j+1);
		tp_vn_bin[j] -> Scale( 1.0 / TMath::Sqrt( TMath::Abs(fResolutionTP->GetBinContent(fResolutionTP->FindBin((Double_t)j)))) );
	}
	
	if(CentBinMin != CentBinMax){
		for(Int_t j=CentBinMin+1; j<=CentBinMax; j++){
			tp_vn_bin[CentBinMin]->Add(tp_vn_bin[j]);
		}
	}

	TH1 *tp_flow = tp_vn_bin[CentBinMin]->Rebin((int)rebin_vector.size()-1, Form("RebinTH1_%s",OutputGraphName.Data()), &rebin_vector[0]);

	for(Int_t i=0; i < (int)tp_flow->GetNbinsX(); i++){
		fGraphX.push_back(tp_flow->GetBinCenter(tp_flow->FindBin(rebin_vector[i])));
		//pt.push_back(tp_meanPt->GetBinCenter(i+1));
		fGraphXErrors.push_back(0.);
		fGraphY.push_back(tp_flow->GetBinContent(tp_flow->FindBin(fGraphX[i])));
		fGraphYErrors.push_back(tp_flow->GetBinError(tp_flow->FindBin(fGraphX[i])));
	}

	TGraphErrors *gr = new TGraphErrors((int)fGraphX.size(), &fGraphX[0], &fGraphY[0], &fGraphXErrors[0], &fGraphYErrors[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::ResolutionEPCent9_Analitic(TFile *InputFile, TString InputResolutionTPName, Int_t InputHarm, TString OutputGraphName, Int_t OutputHarm){

	std::vector<Double_t> x = {2.5,7.5,15,25,35,45,55,65,75};
	std::vector<Double_t> y;
	std::vector<Double_t> x_err;
	std::vector<Double_t> y_err;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName.Data());
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	for(Int_t i=8; i>=0; i--){
		Int_t nBin = fResolutionTP->FindBin(i);
		x_err.push_back(0.0);
		y.push_back(GetRes( fScale*GetChi( TMath::Sqrt(TMath::Abs(fResolutionTP -> GetBinContent(nBin))), InputHarm, 100), OutputHarm) );
		y_err.push_back(  fResolutionTP->GetBinError(nBin) / TMath::Sqrt(TMath::Abs(fResolutionTP->GetBinContent(nBin))) );
	}

	TGraphErrors *gr = new TGraphErrors(9, &x[0], &y[0], &x_err[0], &y_err[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::ResolutionEP_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString OutputGraphName){

	std::vector<Double_t> x = {2.5,7.5,15,25,35,45,55,65,75};
	std::vector<Double_t> y;
	std::vector<Double_t> x_err;
	std::vector<Double_t> y_err;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	for(Int_t i=8; i>=0; i--){
		Int_t nBin = fResolutionTP->FindBin(i);
		x_err.push_back(0.0);
		y.push_back( TMath::Sqrt(TMath::Abs(fScale * fResolutionTP->GetBinContent(nBin))));
		y_err.push_back(  fScale * fResolutionTP -> GetBinError(nBin) / TMath::Sqrt(TMath::Abs(fScale * fResolutionTP -> GetBinContent(nBin))) );
	}

	TGraphErrors *gr = new TGraphErrors(9, &x[0], &y[0], &x_err[0], &y_err[0]);
	gr->SetName(OutputGraphName);
	return gr;

}

TGraphErrors* AnalysisFunction::RatioGraphPointToPoint(TGraphErrors *numerator, TGraphErrors *denominator, TString OutputGraphName){

	std::vector<Double_t> x;
	std::vector<Double_t> x_err;
	std::vector<Double_t> y;
	std::vector<Double_t> y_err;

	if(!numerator){
		std::cout<<"Error, num"<<std::endl;
		return nullptr;
	}
	if(!denominator){
		std::cout<<"Error, den"<<std::endl;
		return nullptr;
	}

	for(Int_t i = 0; i < numerator->GetN(); i++){
		
		x.push_back(numerator->GetPointX(i));
		x_err.push_back(0.);
		y.push_back(numerator->GetPointY(i) / denominator->GetPointY(i));
		y_err.push_back(AnalysisFunction::ErrorYRatio(numerator->GetPointY(i),denominator->GetPointY(i),numerator->GetErrorY(i),denominator->GetErrorY(i)));

	}

	TGraphErrors *gr = new TGraphErrors((int)x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
	gr->SetName(OutputGraphName);
	return gr;
}

TGraphErrors* AnalysisFunction::RatioGraphEval(TGraphErrors *numerator, TGraphErrors *denominator, TString OutputGraphName){

	TGraphErrors *graph_errors_num = new TGraphErrors();
	TGraphErrors *graph_errors_den = new TGraphErrors();
	std::vector<Double_t> x;
	std::vector<Double_t> x_err;
	std::vector<Double_t> y;
	std::vector<Double_t> y_err;

	if(!numerator){
		std::cout<<"Error, num"<<std::endl;
		return nullptr;
	}
	if(!denominator){
		std::cout<<"Error, den"<<std::endl;
		return nullptr;
	}
	for(Int_t i = 0; i < numerator->GetN(); i++){
		graph_errors_num->SetPoint(i,numerator->GetPointX(i),numerator->GetErrorY(i));
	}
	for(Int_t i = 0; i < denominator->GetN(); i++){
		graph_errors_den->SetPoint(i,denominator->GetPointX(i),denominator->GetErrorY(i));
	}

	for(Int_t i = 0; i < numerator->GetN(); i++){
		x.push_back(numerator->GetPointX(i));
		x_err.push_back(0.);
		y.push_back(numerator->Eval(x.back()) / denominator->Eval(x.back()) );
		y_err.push_back(AnalysisFunction::ErrorYRatio(numerator->Eval(x.back()),denominator->Eval(x.back()),graph_errors_num->Eval(x.back()),graph_errors_den->Eval(x.back())));
		std::cout<<i<<"\t"<<numerator->Eval(x.back())<<"\t"<<denominator->Eval(x.back())<<"\t"<<y.back()<<"\n";

	}

	TGraphErrors *gr = new TGraphErrors((int)x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
	gr->SetName(OutputGraphName);
	return gr;
}

TGraphErrors* AnalysisFunction::RatioGraphToOnePoint(TGraphErrors *numerator, TGraphErrors *denominator, Int_t point, TString OutputGraphName){

	std::vector<Double_t> x;
	std::vector<Double_t> x_err;
	std::vector<Double_t> y;
	std::vector<Double_t> y_err;

	if(!numerator){
		std::cout<<"Error, num"<<std::endl;
		return nullptr;
	}
	if(!denominator){
		std::cout<<"Error, den"<<std::endl;
		return nullptr;
	}

	for(Int_t i = 0; i < numerator->GetN(); i++){
		
		x.push_back(numerator->GetPointX(i));
		x_err.push_back(0.);
		y.push_back(numerator->GetPointY(i) / denominator->GetPointY(point));
		y_err.push_back(AnalysisFunction::ErrorYRatio(numerator->GetPointY(i),denominator->GetPointY(point),numerator->GetErrorY(i),denominator->GetErrorY(point)));

	}

	TGraphErrors *gr = new TGraphErrors((int)x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
	gr->SetName(OutputGraphName);
	return gr;
}

TGraphErrors* AnalysisFunction::DifferenceGraphPointToPoint(TGraphErrors *numerator, TGraphErrors *denominator, TString OutputGraphName){
	
	std::vector<Double_t> x;
	std::vector<Double_t> x_err;
	std::vector<Double_t> y;
	std::vector<Double_t> y_err;

	if(!numerator){
		std::cout<<"Error, num"<<std::endl;
		return nullptr;
	}
	if(!denominator){
		std::cout<<"Error, den"<<std::endl;
		return nullptr;
	}

	for(Int_t i = 0; i < numerator->GetN(); i++){
		
		x.push_back(numerator->GetPointX(i));
		x_err.push_back(0.);
		y.push_back(numerator->GetPointY(i) - denominator->GetPointY(i));
		y_err.push_back(numerator->GetErrorY(i) + denominator->GetErrorY(i));

	}

	TGraphErrors *gr = new TGraphErrors((int)x.size(), &x[0], &y[0], &x_err[0], &y_err[0]);
	gr->SetName(OutputGraphName);
	return gr;
}











////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
TGraphErrors* AnalysisFunction::FlowVsRapidity_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, Double_t PtMin, Double_t PtMax, std::vector<double> rebin_vector)
{

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DName<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	fFlowVsPtVsCentTP3D = (TProfile3D*)InputFile->FindObjectAny(InputFlowTP2DName);
	if(!fFlowVsPtVsCentTP3D){
		std::cout<<"Error!!! Not find TProfile3D:\t"<<InputFlowTP2DName<<std::endl;
		return nullptr;
	}

	return AnalysisFunction::CalculateFlowVsRapidity(OutputGraphName, CentBinMin, CentBinMax, PtMin, PtMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsRapidity_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, Double_t PtMin, Double_t PtMax, std::vector<double> rebin_vector)
{

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DNamePos<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile3D *tp_vnpos = (TProfile3D*)InputFile->FindObjectAny(InputFlowTP2DNamePos);
	TProfile3D *tp_vnneg = (TProfile3D*)InputFile->FindObjectAny(InputFlowTP2DNameNeg);
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile3D:\t"<<InputFlowTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_vnneg){
		std::cout<<"Error!!! Not find TProfile3D:\t"<<InputFlowTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP3D = (TProfile3D*)tp_vnpos->Clone();
	fFlowVsPtVsCentTP3D->Add( (TProfile3D*)tp_vnneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsRapidity(OutputGraphName, CentBinMin, CentBinMax, PtMin, PtMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsRapidity_ScalarProduct2(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, std::vector<double> rebin_vector)
{

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DNamePos<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile2D *tp_vnpos = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNamePos);
	TProfile2D *tp_vnneg = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNameNeg);
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_vnneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP2D = (TProfile2D*)tp_vnpos->Clone();
	fFlowVsPtVsCentTP2D->Add( (TProfile2D*)tp_vnneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsRapidity2(OutputGraphName, CentBinMin, CentBinMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsRapidity_ScalarProduct2(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax, std::vector<double> rebin_vector)
{

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DName<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile2D *tp_vnpos = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DName);
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DName<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP2D = (TProfile2D*)tp_vnpos->Clone();
	//fFlowVsPtVsCentTP2D->Add( (TProfile2D*)tp_vnneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsRapidity2(OutputGraphName, CentBinMin, CentBinMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsPt_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	std::cout<<InputResolutionTPName<<"\n"<<InputFlowTP2DNamePos<<"\t"<<InputMeanPtTP2DNamePos<<std::endl;
	std::cout<<InputFlowTP2DNameNeg<<"\t"<<InputMeanPtTP2DNameNeg<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile2D *tp_vnpos = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNamePos);
	TProfile2D *tp_vnneg = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNameNeg);
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_vnneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP2D = (TProfile2D*)tp_vnpos->Clone();
	fFlowVsPtVsCentTP2D->Add( (TProfile2D*)tp_vnneg->Clone(), 1.0 );
	

	TProfile2D *tp_ptpos = (TProfile2D*)InputFile->FindObjectAny(InputMeanPtTP2DNamePos);
	TProfile2D *tp_ptneg = (TProfile2D*)InputFile->FindObjectAny(InputMeanPtTP2DNameNeg);
	if(!tp_ptpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputMeanPtTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_ptneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputMeanPtTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fMeanPtVsPtVsCentTP2D = (TProfile2D*)tp_ptpos->Clone();
	fMeanPtVsPtVsCentTP2D->Add( (TProfile2D*)tp_ptneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsPtEtaSub(OutputGraphName, CentBinMin, CentBinMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsPt_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DName<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	fFlowVsPtVsCentTP2D = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DName);
	if(!fFlowVsPtVsCentTP2D){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DName<<std::endl;
		return nullptr;
	}

	fMeanPtVsPtVsCentTP2D = (TProfile2D*)InputFile->FindObjectAny(InputMeanPtTP2DName);
	if(!fMeanPtVsPtVsCentTP2D){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputMeanPtTP2DName<<std::endl;
		return nullptr;
	}

	return AnalysisFunction::CalculateFlowVsPtEtaSub(OutputGraphName, CentBinMin, CentBinMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsPt_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DName<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	fFlowVsPtVsCentTP2D = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DName);
	if(!fFlowVsPtVsCentTP2D){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DName<<std::endl;
		return nullptr;
	}

	fMeanPtVsPtVsCentTP2D = (TProfile2D*)InputFile->FindObjectAny(InputMeanPtTP2DName);
	if(!fMeanPtVsPtVsCentTP2D){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputMeanPtTP2DName<<std::endl;
		return nullptr;
	}

	return AnalysisFunction::CalculateFlowVsPtScalarProduct(OutputGraphName, CentBinMin, CentBinMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsPt_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Int_t CentBinMin, Int_t CentBinMax,  std::vector<double> rebin_vector){

	std::cout<<InputResolutionTPName<<"\n"<<InputFlowTP2DNamePos<<"\t"<<InputMeanPtTP2DNamePos<<std::endl;
	std::cout<<InputFlowTP2DNameNeg<<"\t"<<InputMeanPtTP2DNameNeg<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile2D *tp_vnpos = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNamePos);
	TProfile2D *tp_vnneg = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNameNeg);
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_vnneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP2D = (TProfile2D*)tp_vnpos->Clone();
	fFlowVsPtVsCentTP2D->Add( (TProfile2D*)tp_vnneg->Clone(), 1.0 );
	

	TProfile2D *tp_ptpos = (TProfile2D*)InputFile->FindObjectAny(InputMeanPtTP2DNamePos);
	TProfile2D *tp_ptneg = (TProfile2D*)InputFile->FindObjectAny(InputMeanPtTP2DNameNeg);
	if(!tp_ptpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputMeanPtTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_ptneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputMeanPtTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fMeanPtVsPtVsCentTP2D = (TProfile2D*)tp_ptpos->Clone();
	fMeanPtVsPtVsCentTP2D->Add( (TProfile2D*)tp_ptneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsPtScalarProduct(OutputGraphName, CentBinMin, CentBinMax, rebin_vector);

}

TGraphErrors* AnalysisFunction::FlowVsCent_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector){

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DName<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	fFlowVsPtVsCentTP2D = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DName);
	if(!fFlowVsPtVsCentTP2D){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DName<<std::endl;
		return nullptr;
	}

	return AnalysisFunction::CalculateFlowVsCentEtaSub(OutputGraphName, PtMin, PtMax, rebin_vector, cent_vector);

}

TGraphErrors* AnalysisFunction::FlowVsCent_EtaSub(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector){

	std::cout<<InputResolutionTPName<<"\n"<<InputFlowTP2DNamePos<<"\t"<<InputMeanPtTP2DNamePos<<std::endl;
	std::cout<<InputFlowTP2DNameNeg<<"\t"<<InputMeanPtTP2DNameNeg<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile2D *tp_vnpos = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNamePos);
	TProfile2D *tp_vnneg = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNameNeg);
	
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_vnneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP2D = (TProfile2D*)tp_vnpos->Clone();
	fFlowVsPtVsCentTP2D->Add( (TProfile2D*)tp_vnneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsCentEtaSub(OutputGraphName, PtMin, PtMax, rebin_vector, cent_vector);

}

TGraphErrors* AnalysisFunction::FlowVsCent_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DName, TString InputMeanPtTP2DName, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector){

	std::cout<<InputResolutionTPName<<"\t"<<InputFlowTP2DName<<"\t"<<InputMeanPtTP2DName<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	fFlowVsPtVsCentTP2D = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DName);
	if(!fFlowVsPtVsCentTP2D){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DName<<std::endl;
		return nullptr;
	}

	return AnalysisFunction::CalculateFlowVsCentScalarProduct(OutputGraphName, PtMin, PtMax, rebin_vector, cent_vector);

}

TGraphErrors* AnalysisFunction::FlowVsCent_ScalarProduct(TFile *InputFile, TString InputResolutionTPName, TString InputFlowTP2DNamePos, TString InputFlowTP2DNameNeg, TString InputMeanPtTP2DNamePos, TString InputMeanPtTP2DNameNeg, TString OutputGraphName, Double_t PtMin, Double_t PtMax,  std::vector<double> rebin_vector, std::vector<double> cent_vector){

	std::cout<<InputResolutionTPName<<"\n"<<InputFlowTP2DNamePos<<"\t"<<InputMeanPtTP2DNamePos<<std::endl;
	std::cout<<InputFlowTP2DNameNeg<<"\t"<<InputMeanPtTP2DNameNeg<<"\n"<<OutputGraphName<<std::endl;

	fResolutionTP = (TProfile*)InputFile->FindObjectAny(InputResolutionTPName);
	if(!fResolutionTP){
		std::cout<<"Error!!! Not find TProfile:\t"<<InputResolutionTPName<<std::endl;
		return nullptr;
	}

	TProfile2D *tp_vnpos = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNamePos);
	TProfile2D *tp_vnneg = (TProfile2D*)InputFile->FindObjectAny(InputFlowTP2DNameNeg);
	
	if(!tp_vnpos){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNamePos<<std::endl;
		return nullptr;
	}
	if(!tp_vnneg){
		std::cout<<"Error!!! Not find TProfile2D:\t"<<InputFlowTP2DNameNeg<<std::endl;
		return nullptr;
	}
	fFlowVsPtVsCentTP2D = (TProfile2D*)tp_vnpos->Clone();
	fFlowVsPtVsCentTP2D->Add( (TProfile2D*)tp_vnneg->Clone(), 1.0 );

	return AnalysisFunction::CalculateFlowVsCentScalarProduct(OutputGraphName, PtMin, PtMax, rebin_vector, cent_vector);

}
