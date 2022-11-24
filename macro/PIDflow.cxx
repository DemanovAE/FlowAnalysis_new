#include "PIDflow.h"

ClassImp(PIDflow);

PIDflow::PIDflow() : 
	Energy_(0),
	fixCuts_(false),
	nSigmaM2_(false),
	CombPID_(false),
	th1_projX(nullptr),
	th1_projY(nullptr)
{
}

PIDflow::PIDflow(const Int_t &Energy, const Double_t pt_left, const Double_t pt_right, const Bool_t fixCuts, const Bool_t nSigmaM2, const Bool_t CombPID) : 
	Energy_(Energy),
	ptRight(pt_right),
	ptLeft(pt_left),
	fixCuts_(fixCuts),
	nSigmaM2_(nSigmaM2),
	CombPID_(CombPID),
	th1_projX(nullptr),
	th1_projY(nullptr)
{
}

TF1 *PIDflow::fitFun1DGaus(TH1D *histo, Double_t cons, Double_t mean, Double_t meanL, Double_t meanR, Double_t sigma, Double_t min, Double_t max, Int_t color, Int_t drawFlag){
	
	TF1 *gaus = new TF1("","gaus",min,max);
	gaus->SetParameter(0,cons);
	gaus->SetParameter(1,mean);
	gaus->SetParameter(2,sigma);

	gaus->SetParLimits(0, 0., (Double_t)histo->GetEntries());
	gaus->SetParLimits(1, meanL, meanR);
	gaus->SetParLimits(2, 0., 10.);

	gaus->SetLineColor(color);

	if(drawFlag==1)histo->Fit(gaus,"MR");
	return gaus;

}

TF1 *PIDflow::fitFunThree1DGaus(TH1D *histo,TF1 *FirstGaus, TF1 *SecondGaus, TF1 *ThirdGaus, Double_t min, Double_t max){

	TF1 *gaus = new TF1("","gaus(0)+gaus(3)+gaus(6)",min,max);
	gaus->SetParameter(0,FirstGaus->GetParameter(0));
	gaus->SetParameter(1,FirstGaus->GetParameter(1));
	gaus->SetParameter(2,FirstGaus->GetParameter(2));
	gaus->SetParameter(3,SecondGaus->GetParameter(0));
	gaus->SetParameter(4,SecondGaus->GetParameter(1));
	gaus->SetParameter(5,SecondGaus->GetParameter(2));
	gaus->SetParameter(6,ThirdGaus->GetParameter(0));
	gaus->SetParameter(7,ThirdGaus->GetParameter(1));
	gaus->SetParameter(8,ThirdGaus->GetParameter(2));

	gaus->SetParLimits(0, 0., 1.1*FirstGaus->GetParameter(0));
	gaus->SetParLimits(1, 0.9*FirstGaus->GetParameter(1), 1.1*FirstGaus->GetParameter(1));
	gaus->SetParLimits(2, 0.5*FirstGaus->GetParameter(2), 1.1*FirstGaus->GetParameter(2));

	gaus->SetParLimits(3, 0., 1.1*SecondGaus->GetParameter(0));
	gaus->SetParLimits(4, 0.9*SecondGaus->GetParameter(1), 1.1*SecondGaus->GetParameter(1));
	gaus->SetParLimits(5, 0.5*SecondGaus->GetParameter(2), 1.1*SecondGaus->GetParameter(2));
	
	gaus->SetParLimits(6, 0., 1.1*ThirdGaus->GetParameter(0));
	gaus->SetParLimits(7, 0.9*ThirdGaus->GetParameter(1), 1.1*ThirdGaus->GetParameter(1));
	gaus->SetParLimits(8, 0.5*ThirdGaus->GetParameter(2), 1.1*ThirdGaus->GetParameter(2));

	gaus->SetLineColor(kBlack);
	histo->Fit(gaus,"RM");
	return gaus;

}

TF1 *PIDflow::fitFunThree2x1DGaus(TH1D *histo,TF1 *Three1DGaus, Double_t min, Double_t max ){

	TF1 *gaus = new TF1("6gs","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",min,max);
	
	gaus->SetParameter(0,Three1DGaus->GetParameter(0));
	gaus->FixParameter(1,Three1DGaus->GetParameter(1));
	gaus->SetParameter(2,Three1DGaus->GetParameter(2));

	gaus->SetParameter(3,Three1DGaus->GetParameter(0)/50);
	gaus->FixParameter(4,Three1DGaus->GetParameter(1));
	gaus->SetParameter(5,Three1DGaus->GetParameter(2)*5.0);
	gaus->SetParLimits(3,1,Three1DGaus->GetParameter(0)/10);
	gaus->SetParLimits(5,Three1DGaus->GetParameter(2)*2.0,Three1DGaus->GetParameter(2)*10.0);

	gaus->SetParameter(6,Three1DGaus->GetParameter(3));
	gaus->FixParameter(7,Three1DGaus->GetParameter(4));
	gaus->SetParameter(8,Three1DGaus->GetParameter(5));
	
	gaus->SetParameter(9,Three1DGaus->GetParameter(3)/50);
	gaus->FixParameter(10,Three1DGaus->GetParameter(4));
	gaus->SetParameter(11,Three1DGaus->GetParameter(5)*5.0);
	gaus->SetParLimits(9,1,Three1DGaus->GetParameter(3)/10);
	gaus->SetParLimits(11,Three1DGaus->GetParameter(5)*2.0,Three1DGaus->GetParameter(5)*10.0);


	gaus->SetParameter(12,Three1DGaus->GetParameter(6));
	gaus->FixParameter(13,Three1DGaus->GetParameter(7));
	gaus->SetParameter(14,Three1DGaus->GetParameter(8));
	
	gaus->SetParameter(15,Three1DGaus->GetParameter(6)/50);
	gaus->FixParameter(16,Three1DGaus->GetParameter(7));
	gaus->SetParameter(17,Three1DGaus->GetParameter(8)*5.0);
	gaus->SetParLimits(15,1,Three1DGaus->GetParameter(6)/10);
	gaus->SetParLimits(17,Three1DGaus->GetParameter(7)*2.0,Three1DGaus->GetParameter(7)*10.0);

	gaus->SetLineColor(kBlack);

	histo->Fit(gaus,"RM");
	return gaus;

}

void PIDflow::FitSquareOfMassInPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto, Int_t charge){

	//gStyle->SetOptStat("n");
    //gStyle->SetOptFit(0);

    TH2D *h2_m2VsPt = (TH2D*)fileWithHisto->Get(NameHisto);
    th1_projX = (TH1D*)h2_m2VsPt->ProjectionX();
    th1_projY = (TH1D*)h2_m2VsPt->ProjectionY(Form("%s_py_%i_%i_ch%i", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projX->FindBin(ptLeft),th1_projX->FindBin(ptRight));
    th1_projY->SetTitle(Form("m^{2}, %.1f < p_{T} < %.1f [GeV/c]",ptLeft,ptRight));


    Double_t min[4]={-0.2, 0.2, 0.6};
    Double_t max[4]={ 0.1, 0.5, 1.4};
    Double_t massSqrt[3]={pion_mass*pion_mass, kaon_mass*kaon_mass, proton_mass*proton_mass};
    Double_t massSqrtLimitsLeft[3] ={-0.06, 0.18, 0.8};
    Double_t massSqrtLimitsRight[3]={ 0.02, 0.24, 0.9};
    Int_t colorFitLine[4]={2,3,4,1};

    fileOut->cd();

    if(ptRight < 1.00){
    	tf1_gausParticle[0] = fitFun1DGaus(th1_projY, 10., massSqrt[0], 0.010, 0.021, 0.01,  0.0, 0.15, colorFitLine[0],1);
    	tf1_gausParticle[1] = fitFun1DGaus(th1_projY, 10., massSqrt[1], 0.225, 0.255, 0.02,  0.2, 0.50, colorFitLine[1],1);
    	tf1_gausParticle[2] = fitFun1DGaus(th1_projY, 10., massSqrt[2], 0.800, 0.900, 0.01,  0.6, 1.20, colorFitLine[2],1);
    }

	if(ptRight < 2.21){
		tf1_gausParticle[0] = fitFun1DGaus(th1_projY, 10., massSqrt[0],-0.020, 0.020, 0.1, -0.1, 0.15, colorFitLine[0],1);
    	tf1_gausParticle[1] = fitFun1DGaus(th1_projY, 10., massSqrt[1], 0.215, 0.240, 0.2,  0.2, 0.50, colorFitLine[1],1);
    	tf1_gausParticle[2] = fitFun1DGaus(th1_projY, 10., massSqrt[2], 0.800, 0.900, 0.1,  0.6, 1.20, colorFitLine[2],1);
    }

    if(ptRight > 2.21){
    	tf1_gausParticle[0] = fitFun1DGaus(th1_projY, 10., massSqrt[0],-0.030,-0.010, 0.1, -0.1, 0.10, colorFitLine[0],1);
    	tf1_gausParticle[1] = fitFun1DGaus(th1_projY, 10., massSqrt[1], 0.205, 0.235, 0.2,  0.2, 0.40, colorFitLine[1],1);
    	tf1_gausParticle[2] = fitFun1DGaus(th1_projY, 10., massSqrt[2], 0.800, 0.900, 0.1,  0.7, 1.40, colorFitLine[2],1);
    }

	tf1_Three1DGaus = fitFunThree1DGaus(th1_projY, tf1_gausParticle[0], tf1_gausParticle[1], tf1_gausParticle[2], -0.5, 1.5);
    tf1_Three1DGaus->SetName(Form("tf1_FitM2_Three1DGaus_charge%i_pt%i_%i",charge, (int)(ptLeft*10),(int)(ptRight*10))); 
	tf1_Three1DGaus->SetTitle(Form("Three 1D Gaus, %.1f < p_{T} < %.1f [GeV/c]",ptLeft,ptRight));

	tf1_Three2x1DGaus = fitFunThree2x1DGaus(th1_projY, tf1_Three1DGaus, -0.5, 1.5);
    tf1_Three2x1DGaus->SetName(Form("tf1_FitM2_Three2x1DGaus_charge%i_pt%i_%i",charge, (int)(ptLeft*10),(int)(ptRight*10))); 
	tf1_Three2x1DGaus->SetTitle(Form("Three 2x1D Gaus, %.1f < p_{T} < %.1f [GeV/c]",ptLeft,ptRight));

	canvas = new TCanvas(Form("canvas_pt_%i_%i_ch%i",(int)(ptLeft*10),(int)(ptRight*10),charge),"",1024,1024);
	canvas->SetLogy();
	th1_projY->Draw();
	tf1_Three1DGaus->Draw("Same");
	tf1_gausParticle[0] = fitFun1DGaus(th1_projY, tf1_Three1DGaus->GetParameter(0), tf1_Three1DGaus->GetParameter(1),-0.050, 1.020, tf1_Three1DGaus->GetParameter(2), -0.5, 1.5, colorFitLine[0],0);
	tf1_gausParticle[1] = fitFun1DGaus(th1_projY, tf1_Three1DGaus->GetParameter(3), tf1_Three1DGaus->GetParameter(4),-0.050, 1.020, tf1_Three1DGaus->GetParameter(5), -0.5, 1.5, colorFitLine[1],0);
	tf1_gausParticle[2] = fitFun1DGaus(th1_projY, tf1_Three1DGaus->GetParameter(6), tf1_Three1DGaus->GetParameter(7),-0.050, 1.020, tf1_Three1DGaus->GetParameter(8), -0.5, 1.5, colorFitLine[2],0);
	
	tf1_gausParticle[0]->Draw("Same");
	tf1_gausParticle[1]->Draw("Same");
	tf1_gausParticle[2]->Draw("Same");

	//canvas->SaveAs(Form("./pict/Can_pt_%i_%i_ch%i.png",(int)(ptLeft*10),(int)(ptRight*10),charge));
	canvas->Write();
	th1_projY->Write();
	tf1_Three1DGaus->Write();
	tf1_Three2x1DGaus->Write();
			
}

void PIDflow::SetInitialValForFit_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto, Int_t charge){

    Double_t massSqrt[3]={pion_mass*pion_mass, kaon_mass*kaon_mass, proton_mass*proton_mass};
    Int_t colorFitLine[4]={2,3,4,1};
	
	Double_t factor=2.0; // 27GeV Run10, 62.4 GeV Run10
	if(Energy_==27 || Energy_==62){
		factor=1.0;
	}

	TF1 *tf1_gausParticle_nSigma[3];
	TF1 *tf1_gausParticle_M2[3];
	TF1 *tf1_Three1DGaus_M2;

	TH1D *th1_projY_M2[3];
	TH1D *th1_projX_nSigma[3];

	TH2D *h2_m2V_nSigma = (TH2D*)fileWithHisto->Get(NameHisto);
    th1_projY_dEdx_M2 = (TH1D*)h2_m2V_nSigma->ProjectionY();
    th1_projX_dEdx_M2 = (TH1D*)h2_m2V_nSigma->ProjectionX();
    
    th1_projX_nSigma[0] = (TH1D*)h2_m2V_nSigma->ProjectionX(Form("%s_px_%i_%i_ch%i_par0", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projY_dEdx_M2->FindBin(massSqrt[0]),th1_projY_dEdx_M2->FindBin(massSqrt[0]));
    th1_projX_nSigma[1] = (TH1D*)h2_m2V_nSigma->ProjectionX(Form("%s_px_%i_%i_ch%i_par1", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projY_dEdx_M2->FindBin(massSqrt[1]),th1_projY_dEdx_M2->FindBin(massSqrt[1]));
    th1_projX_nSigma[2] = (TH1D*)h2_m2V_nSigma->ProjectionX(Form("%s_px_%i_%i_ch%i_par2", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projY_dEdx_M2->FindBin(massSqrt[2]),th1_projY_dEdx_M2->FindBin(massSqrt[2]));
    
    th1_projX_nSigma[0]->SetTitle(Form("n#sigma(#pi), %.1f < p_{T} < %.1f [GeV/c], #pi, %i",ptLeft,ptRight,charge));
    th1_projX_nSigma[1]->SetTitle(Form("n#sigma(#pi), %.1f < p_{T} < %.1f [GeV/c], K, %i",ptLeft,ptRight, charge));
    th1_projX_nSigma[2]->SetTitle(Form("n#sigma(#pi), %.1f < p_{T} < %.1f [GeV/c], p, %i",ptLeft,ptRight, charge));

    Double_t mean_dEdx_init_vals_range[3][2] = {{-2.0,2.0},{0.0,8.0},{1.0,20.0}};
    Double_t mean_dEdx_init_vals[3] = {0.0, 4.0, 8.0};
    Double_t mean_M2_init_vals_range[3][2] = {{0.01,0.02},{0.22,0.245},{0.83,0.94}};
    Double_t mean_M2_init_vals[3] = {massSqrt[0], massSqrt[1], massSqrt[2]};

    if(ptRight>0.6 && ptRight<=2.0){
    	
    	mean_dEdx_init_vals_range[0][0]=-2.0 * factor;
    	mean_dEdx_init_vals_range[0][1]= 2.0 * factor;
    	mean_dEdx_init_vals_range[1][0]=-1.5 * factor;
    	mean_dEdx_init_vals_range[1][1]= 1.5 * factor;
		mean_dEdx_init_vals_range[2][0]=-2.0 * factor;
    	mean_dEdx_init_vals_range[2][1]= 10.;

     	mean_dEdx_init_vals[0]= 0.0;
     	mean_dEdx_init_vals[1]=-1.0;
     	mean_dEdx_init_vals[2]= 0.0;

     	mean_M2_init_vals_range[0][0] =-0.02;
     	mean_M2_init_vals_range[0][1] = 0.02;
    	mean_M2_init_vals_range[1][0] = 0.22;
     	mean_M2_init_vals_range[1][1] = 0.245;
     	mean_M2_init_vals_range[2][0] = 0.85;
     	mean_M2_init_vals_range[2][1] = 0.91;
    }

    if(ptRight>2.0){

    	mean_dEdx_init_vals_range[0][0]=-1.0 * factor;
    	mean_dEdx_init_vals_range[0][1]= 1.0 * factor;
    	mean_dEdx_init_vals_range[1][0]=-2.5 * factor;
    	mean_dEdx_init_vals_range[1][1]= 0.0 * factor;
		mean_dEdx_init_vals_range[2][0]=-2.0 * factor;
    	mean_dEdx_init_vals_range[2][1]= 2.0 * factor;

     	mean_dEdx_init_vals[0]= 0.0;
     	mean_dEdx_init_vals[1]=-1.0;
     	mean_dEdx_init_vals[2]= 0.0;

     	mean_M2_init_vals_range[0][0] =-0.03;
     	mean_M2_init_vals_range[0][1] = 0.01;
    	mean_M2_init_vals_range[1][0] = 0.185;
     	mean_M2_init_vals_range[1][1] = 0.24;
     	mean_M2_init_vals_range[2][0] = 0.85;
     	mean_M2_init_vals_range[2][1] = 0.91;
    }

	tf1_gausParticle_nSigma[0] = fitFun1DGaus(th1_projX_nSigma[0], th1_projX_nSigma[0]->GetMaximum(), mean_dEdx_init_vals[0], mean_dEdx_init_vals_range[0][0],mean_dEdx_init_vals_range[0][1], 0.1, -3., 3., colorFitLine[0],1);
    tf1_gausParticle_nSigma[1] = fitFun1DGaus(th1_projX_nSigma[1], th1_projX_nSigma[1]->GetMaximum(), mean_dEdx_init_vals[1], mean_dEdx_init_vals_range[1][0],mean_dEdx_init_vals_range[1][1], 0.1, -3., 10., colorFitLine[1],1);
    tf1_gausParticle_nSigma[2] = fitFun1DGaus(th1_projX_nSigma[2], th1_projX_nSigma[2]->GetMaximum(), mean_dEdx_init_vals[2], mean_dEdx_init_vals_range[2][0],mean_dEdx_init_vals_range[2][1], 0.1, -3., 20., colorFitLine[2],1);

    th1_projY_M2[0] = (TH1D*)h2_m2V_nSigma->ProjectionY(Form("%s_py_%i_%i_ch%i_par0_m2", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projX_dEdx_M2->FindBin(tf1_gausParticle_nSigma[0]->GetParameter(1)),th1_projX_dEdx_M2->FindBin(tf1_gausParticle_nSigma[0]->GetParameter(1)));
    th1_projY_M2[1] = (TH1D*)h2_m2V_nSigma->ProjectionY(Form("%s_py_%i_%i_ch%i_par1_m2", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projX_dEdx_M2->FindBin(tf1_gausParticle_nSigma[1]->GetParameter(1)),th1_projX_dEdx_M2->FindBin(tf1_gausParticle_nSigma[1]->GetParameter(1)));
    th1_projY_M2[2] = (TH1D*)h2_m2V_nSigma->ProjectionY(Form("%s_py_%i_%i_ch%i_par2_m2", NameHisto,(int)(ptLeft*10),(int)(ptRight*10),charge),th1_projX_dEdx_M2->FindBin(tf1_gausParticle_nSigma[2]->GetParameter(1)),th1_projX_dEdx_M2->FindBin(tf1_gausParticle_nSigma[2]->GetParameter(1)));
    

	tf1_gausParticle_M2[0] = fitFun1DGaus(th1_projY_M2[0], th1_projY_M2[0]->GetMaximum(), mean_M2_init_vals[0], mean_M2_init_vals_range[0][0],mean_M2_init_vals_range[0][1], 0.1,-0.3, 0.1, colorFitLine[0],1);
    tf1_gausParticle_M2[1] = fitFun1DGaus(th1_projY_M2[1], th1_projY_M2[1]->GetMaximum(), mean_M2_init_vals[1], mean_M2_init_vals_range[1][0],mean_M2_init_vals_range[1][1], 0.1, 0.15, 0.55, colorFitLine[1],1);
    tf1_gausParticle_M2[2] = fitFun1DGaus(th1_projY_M2[2], th1_projY_M2[2]->GetMaximum(), mean_M2_init_vals[2], mean_M2_init_vals_range[2][0],mean_M2_init_vals_range[2][1], 0.1, -0.6, 1.5, colorFitLine[2],1);

    Double_t constM2[3]={tf1_gausParticle_M2[0]->GetParameter(0),tf1_gausParticle_M2[1]->GetParameter(0),tf1_gausParticle_M2[2]->GetParameter(0)};
    Double_t meanM2[3] ={tf1_gausParticle_M2[0]->GetParameter(1),tf1_gausParticle_M2[1]->GetParameter(1),tf1_gausParticle_M2[2]->GetParameter(1)};
    Double_t sigmaM2[3]={tf1_gausParticle_M2[0]->GetParameter(2),tf1_gausParticle_M2[1]->GetParameter(2),tf1_gausParticle_M2[2]->GetParameter(2)};

    if(ptRight>0.8){
    	tf1_Three1DGaus_M2 = fitFunThree1DGaus(th1_projY_M2[1], tf1_gausParticle_M2[0], tf1_gausParticle_M2[1], tf1_gausParticle_M2[2], -0.5, 1.5);
    	tf1_Three1DGaus_M2->SetName(Form("tf1_FitM2_Three1DGaus_charge%i_pt%i_%i",charge, (int)(ptLeft*10),(int)(ptRight*10))); 
		tf1_Three1DGaus_M2->SetTitle(Form("Three 1D Gaus, %.1f < p_{T} < %.1f [GeV/c]",ptLeft,ptRight));
    
    	for(Int_t i=0; i<3; i++){
			constM2[i]= tf1_gausParticle_M2[i]->GetParameter(0);
     		meanM2[i] = tf1_gausParticle_M2[i]->GetParameter(1);
     		sigmaM2[i]= tf1_gausParticle_M2[i]->GetParameter(2);
     	}
    }

    tf2_init_vals[2][3][3];  // x - nSigma,y-m2 // const, mean, sigma // particle
    tf2_init_vals_range[2][3][3][2]; // x,y // const, mean, sigma // particle // min,max

    // const nSigma
    tf2_init_vals[0][0][0] = tf1_gausParticle_nSigma[0]->GetParameter(0);
    tf2_init_vals[0][0][1] = tf1_gausParticle_nSigma[1]->GetParameter(0);
    tf2_init_vals[0][0][2] = tf1_gausParticle_nSigma[2]->GetParameter(0);
    
    //const M2
    tf2_init_vals[1][0][0] = constM2[0]; 
    tf2_init_vals[1][0][1] = constM2[1];
    tf2_init_vals[1][0][2] = constM2[2];

    // mean nSigma
    tf2_init_vals[0][1][0] = tf1_gausParticle_nSigma[0]->GetParameter(1);
    tf2_init_vals[0][1][1] = tf1_gausParticle_nSigma[1]->GetParameter(1);
    tf2_init_vals[0][1][2] = tf1_gausParticle_nSigma[2]->GetParameter(1);

    //mean M2
    tf2_init_vals[1][1][0] = meanM2[0]; 
    tf2_init_vals[1][1][1] = meanM2[1];
    tf2_init_vals[1][1][2] = meanM2[2];

    // sigma nSigma
    tf2_init_vals[0][2][0] = tf1_gausParticle_nSigma[0]->GetParameter(2);
    tf2_init_vals[0][2][1] = tf1_gausParticle_nSigma[1]->GetParameter(2);
    tf2_init_vals[0][2][2] = tf1_gausParticle_nSigma[2]->GetParameter(2);

    //sigma M2
    tf2_init_vals[1][2][0] = sigmaM2[0]; 
    tf2_init_vals[1][2][1] = sigmaM2[1];
    tf2_init_vals[1][2][2] = sigmaM2[2];

    //const
    tf2_init_vals_range[0][0][0][0] = 0.0;
    tf2_init_vals_range[0][0][1][0] = 0.0;
    tf2_init_vals_range[0][0][2][0] = 0.0;
    tf2_init_vals_range[1][0][0][0] = 0.0;
    tf2_init_vals_range[1][0][1][0] = 0.0;
    tf2_init_vals_range[1][0][2][0] = 0.0;

    tf2_init_vals_range[0][0][0][1] = 10000 * tf2_init_vals[0][0][0];
    tf2_init_vals_range[0][0][1][1] = 10000 * tf2_init_vals[0][0][1];
    tf2_init_vals_range[0][0][2][1] = 10000 * tf2_init_vals[0][0][2];
    tf2_init_vals_range[1][0][0][1] = 10000 * tf2_init_vals[1][0][0];
    tf2_init_vals_range[1][0][1][1] = 10000 * tf2_init_vals[1][0][1];
    tf2_init_vals_range[1][0][2][1] = 10000 * tf2_init_vals[1][0][2];

    //mean
    for(Int_t par=0; par<3; par++){
    	for(Int_t j=0; j<2; j++){
    		tf2_init_vals_range[0][1][par][j] = mean_dEdx_init_vals_range[par][j];
    		tf2_init_vals_range[1][1][par][j] = mean_M2_init_vals_range[par][j];
    	}
    	tf2_init_vals_range[0][2][par][0] = 0.0001;
    	tf2_init_vals_range[1][2][par][0] = 0.0001;
    	tf2_init_vals_range[0][2][par][1] = 2.0 * tf2_init_vals[0][2][par];
    	tf2_init_vals_range[1][2][par][1] = 2.0 * tf2_init_vals[1][2][par];
    }

    fileOut->cd();
    //h2_m2V_nSigma->Write();
    for(Int_t i=0; i<3; i++){
    	th1_projX_nSigma[i]->Write();
		th1_projY_M2[i]->Write();

    	tf1_gausParticle_nSigma[i]->Write();
		tf1_gausParticle_M2[i]->Write();
	}

    if(ptRight>0.8){
    	tf1_Three1DGaus_M2->Write();
    }

}

void PIDflow::SetInitialValForFit_NewXY_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto, Int_t charge){
	Int_t par=0;
}

void PIDflow::Fit_nSigmaM2_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto, Int_t charge){

	Double_t x_min = -5.0;
	Double_t x_max =  6.0;
	Double_t y_min = -0.5;
	Double_t y_max =  1.5;

	if(ptRight<0.61){
		x_min = -3.0;
		x_max = 12.0;
	}

	if(ptRight<0.81){
		x_min = -3.0;
		x_max = 10.0;
	}

	PIDflow::Fit_TH2(fileWithHisto, fileOut, NameHisto, charge, "",x_min,x_max,y_min,y_max);

}

void PIDflow::Fit_NewXY_InPtBin(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto, Int_t charge){

	Double_t factor=2.0;
	if(Energy_==27 || Energy_==62){
		factor=1.0;
	}
	Double_t x_min = -5.0 * factor;
	Double_t x_max =  6.0 * factor;
	Double_t y_min = -0.5;
	Double_t y_max =  1.5;

	if(ptRight<0.61){
		x_min = -4.0 * factor;
		x_max = 12.0 * factor;
	}

	if(ptRight<0.81){
		x_min = -4.0 * factor;
		x_max = 10.0 * factor;
	}

	PIDflow::Fit_TH2(fileWithHisto, fileOut, NameHisto, charge, "",x_min,x_max,y_min,y_max);

}


void PIDflow::Fit_TH2(TFile *fileWithHisto, TFile *fileOut, const Char_t *NameHisto, Int_t charge, const Char_t *prefix, Double_t AxisXMin, Double_t AxisXMax, Double_t AxisYMin, Double_t AxisYMax){

	TH2D *h2_m2V_nSigma = (TH2D*)fileWithHisto->Get(NameHisto);

	Three2DGaus = new TF2(Form("%sThree2DGaus_pt%i%i_ch%i", prefix,(int)(ptLeft*10),(int)(ptRight*10),charge),"xygaus(0)+xygaus(5)+xygaus(10)", AxisXMin,AxisXMax,AxisYMin,AxisYMax);
	Three2x2DGaus = new TF2(Form("%sThree2x2DGaus_pt%i%i_ch%i", prefix,(int)(ptLeft*10),(int)(ptRight*10),charge),"xygaus(0)+xygaus(5)+xygaus(10)+xygaus(15)+xygaus(20)+xygaus(25)", AxisXMin,AxisXMax,AxisYMin,AxisYMax);

	for(Int_t par=0; par<3; par++){
		std::cout<<"Const \tpar "<<(int)(0+par*5)<<"\t\t"<<tf2_init_vals[0][0][par]<<"\n";
		std::cout<<"MeanS \tpar "<<(int)(1+par*5)<<"\t\t"<<tf2_init_vals[0][1][par]<<"\n";
		std::cout<<"Sigma \tpar "<<(int)(2+par*5)<<"\t\t"<<tf2_init_vals[0][2][par]<<"\n";
		std::cout<<"MeanM2\tpar "<<(int)(3+par*5)<<"\t\t"<<tf2_init_vals[1][1][par]<<"\n";
		std::cout<<"Sigma \tpar "<<(int)(4+par*5)<<"\t\t"<<tf2_init_vals[1][2][par]<<"\n";
	}

	for(Int_t par=0; par<3; par++){
		// x,y // const, mean, sigma // particle
		// x,y // const, mean, sigma // particle // min,max
		Three2DGaus->SetParameter(0 + 5*par, tf2_init_vals[0][0][par]);
		Three2DGaus->SetParameter(1 + 5*par, tf2_init_vals[0][1][par]);
		Three2DGaus->SetParameter(2 + 5*par, tf2_init_vals[0][2][par]);
		Three2DGaus->SetParameter(3 + 5*par, tf2_init_vals[1][1][par]);
		Three2DGaus->SetParameter(4 + 5*par, tf2_init_vals[1][2][par]);
		
		Three2DGaus->SetParLimits(0 + 5*par, tf2_init_vals_range[0][0][par][0], tf2_init_vals_range[0][0][par][1]);
		Three2DGaus->SetParLimits(1 + 5*par, tf2_init_vals_range[0][1][par][0], tf2_init_vals_range[0][1][par][1]);
		Three2DGaus->SetParLimits(2 + 5*par, tf2_init_vals_range[0][2][par][0], tf2_init_vals_range[0][2][par][1]);
		Three2DGaus->SetParLimits(3 + 5*par, tf2_init_vals_range[1][1][par][0], tf2_init_vals_range[1][1][par][1]);
		Three2DGaus->SetParLimits(4 + 5*par, tf2_init_vals_range[1][2][par][0], tf2_init_vals_range[1][2][par][1]);

	}

	std::cout<<"Fit step 1.1:\n";
	h2_m2V_nSigma->Fit(Three2DGaus,"MR");
	std::cout<<"Fit step 1.2:\n";
	h2_m2V_nSigma->Fit(Three2DGaus,"MR");

	Double_t parThree2x2DGaus[30]={0.0};
	Three2DGaus->GetParameters(&parThree2x2DGaus[0]);
	Three2DGaus->GetParameters(&parThree2x2DGaus[15]);
	
	Three2x2DGaus->SetParameters(parThree2x2DGaus);

	for(Int_t i=0; i<6; i++){
		Three2x2DGaus->FixParameter(1+5*i,parThree2x2DGaus[1+5*i]);
		Three2x2DGaus->FixParameter(3+5*i,parThree2x2DGaus[3+5*i]);
	}

	for(Int_t i=0; i<3; i++){
		Three2x2DGaus->SetParLimits(0+5*i, 0.9*parThree2x2DGaus[0+5*i], 1.1*parThree2x2DGaus[0+5*i]);
		Three2x2DGaus->SetParLimits(2+5*i, 0.9*parThree2x2DGaus[2+5*i], 1.1*parThree2x2DGaus[2+5*i]);		
		Three2x2DGaus->SetParLimits(4+5*i, 0.9*parThree2x2DGaus[4+5*i], 1.1*parThree2x2DGaus[4+5*i]);
	}

	for(Int_t i=3; i<6; i++){
		Three2x2DGaus->SetParameter(0+5*i, parThree2x2DGaus[0+5*i]/20.);
		Three2x2DGaus->SetParameter(2+5*i, 5.*parThree2x2DGaus[2+5*i]);
		Three2x2DGaus->SetParameter(4+5*i, 5.*parThree2x2DGaus[4+5*i]);

		Three2x2DGaus->SetParLimits(0+5*i, 0.0, parThree2x2DGaus[0+5*i]/10);
		Three2x2DGaus->SetParLimits(2+5*i, 1.1*parThree2x2DGaus[2+5*i], 20.0*parThree2x2DGaus[2+5*i]);		
		Three2x2DGaus->SetParLimits(4+5*i, 1.1*parThree2x2DGaus[4+5*i], 20.0*parThree2x2DGaus[4+5*i]);
	}
	TH2D *h_clone = (TH2D*)h2_m2V_nSigma->Clone("cl");
	h_clone->SetName(Form("clone_pt%i%i_ch%i",(int)(ptLeft*10),(int)(ptRight*10), charge));
	std::cout<<"\nStep 2.1\n";
	h_clone->Fit(Three2x2DGaus,"RM");
	std::cout<<"\nStep 2.2\n";
	h_clone->Fit(Three2x2DGaus,"RM");


	h2_m2V_nSigma->Write();
	Three2DGaus->Write();

	h_clone->Write();
	Three2x2DGaus->Write();

}

Double_t PIDflow::GetParTF2FromFile(TFile *inFile, const Char_t *NameFitFun, Int_t numberPar){
	Three2x2DGaus = (TF2*)inFile->Get(NameFitFun);
	return Three2x2DGaus->GetParameter(numberPar);
}

Double_t PIDflow::GetParErrTF2FromFile(TFile *inFile, const Char_t *NameFitFun, Int_t numberPar){
	Three2x2DGaus = (TF2*)inFile->Get(NameFitFun);
	return Three2x2DGaus->GetParError(numberPar);
}