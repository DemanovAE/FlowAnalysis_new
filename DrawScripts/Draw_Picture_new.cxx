#include <Draw_Picture_new.h>

Draw_Picture_new::Draw_Picture_new():
    labels(),
    Sys_mod(kFALSE),
    legend_flag(kFALSE),
    fileIn(nullptr),
    fileOut(nullptr),
    hPadAxis(nullptr),
    hPadAxis2(nullptr),
    hPadUp(nullptr),
    hPadDown(nullptr),
    hPadOne(nullptr),
    hPadNM(),
    legend(nullptr),
    legend_nPad(-1),
	_textInfo(),
	type_canvas(-1),
	Line_vector(),
	Line_vector_DrawPad(),
	graph_up(),
	graph_up_sys(),
	graph_down()
{
};
// Деструктор
Draw_Picture_new::~Draw_Picture_new(){
};

void Draw_Picture_new::SetTLineOnPad(Float_t mXmin, Float_t mYmin, Float_t mXmax, Float_t mYmax, Int_t mColour, Int_t mLineWidth, Int_t mLineStyle, Int_t NumberPad) {

	Line_vector.push_back(new TLine( mXmin, mYmin, mXmax, mYmax ));
	Line_vector_DrawPad.push_back(NumberPad);
	Line_vector.back() -> SetLineColor( mColour );
	Line_vector.back() -> SetLineWidth( mLineWidth );
	Line_vector.back() -> SetLineStyle( mLineStyle );

}// 

void Draw_Picture_new::SetAxisToPad(Double_t Xmin, Double_t Xmax, Double_t XLabelSize, Double_t YminUp, Double_t YmaxUp, Double_t YLabelSize){

	hPadAxis = new TH2F(Form("hAxis%.4f", YmaxUp),"", 2,Xmin,Xmax,2,YminUp,YmaxUp);
    hPadAxis->GetXaxis()->SetLabelSize(XLabelSize);
    hPadAxis->GetYaxis()->SetLabelSize(YLabelSize);

}

void Draw_Picture_new::SetAxisToPad2(Double_t Xmin, Double_t Xmax, Double_t XLabelSize, Double_t YminUp, Double_t YmaxUp, Double_t YLabelSize){

	hPadAxis2 = new TH2F(Form("hAxis%.4f", YmaxUp),"", 2,Xmin,Xmax,2,YminUp,YmaxUp);
    hPadAxis2->GetXaxis()->SetLabelSize(XLabelSize);
    hPadAxis2->GetYaxis()->SetLabelSize(YLabelSize);

}

void Draw_Picture_new::SetLegendPar(Double_t xMin, Double_t yMin, Double_t xMax, Double_t yMax, Int_t N, Int_t NumberPad){
	
	legend = new TLegend(xMin,yMin,xMax,yMax);
	legend->SetNColumns(N);
	legend_flag=true;
	legend_nPad=NumberPad;
	
}

void Draw_Picture_new::SetTextInfo(TString text, TString place, Int_t NumberPad, Float_t x, Float_t y, Float_t size, Float_t angle){
	labels.push_back({text,place,NumberPad,x,y,size, angle});
}


TCanvas* Draw_Picture_new::CanvasNxM(Int_t xCanvas, Int_t yCanvas, Int_t line, Int_t column, std::vector<MyDrawObject*> DrawResult, Float_t size_bukva){


	TString bukva[] = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)","(p)", "(q)","(r)","(s)","(t)","(u)","(v)"};
	//Float_t
    
    Double_t coeff=0.2;
	Float_t xx = xCanvas * (column + 0.4 + 0.05);
	Float_t yy = yCanvas * (line + 0.4 + 0.05);

	Float_t otstupUp = 0.05 * yCanvas / yy;
	Float_t otstupDown = 0.4 * yCanvas / yy;
	Float_t otstupLeft = 0.4 * xCanvas / xx;
	Float_t otstupRight = 0.05 * xCanvas / xx;
	Float_t h = (1.0 - otstupUp - otstupDown )/line;
	Float_t l = (1.0 - otstupRight - otstupLeft )/column;
	//TCanvas *canvas;

	//TCanvas *canvas = new TCanvas("canvas","plot",xCanvas*column,yCanvas*line);
	TCanvas *canvas = new TCanvas("canvas","plot", xx, yy);
    canvas->GetFrame()->SetFillColor(21);
    canvas->GetFrame()->SetBorderSize(115);
    canvas->cd();
    
    std::vector<TPad *> pad;
    std::vector<TPad *> pad_Axis;
    Int_t n=0;

    for(Int_t y=line; y>0; y--){
    	for(Int_t x=0; x<column; x++){

		  	//std::cout<<0.1+l*x<<"\t"<< 1.0-h*y<<"\t"<< 0.1+l*(x+1)<<"\t"<< 1.0-h*(y-1)<<"\n\n";
    		pad.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown + h*(y-1)), 
    														(Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
        	pad.back()->SetFrameBorderMode(1);
        	pad.back()->SetLeftMargin(0);
        	pad.back()->SetRightMargin(0);
        	pad.back()->SetBottomMargin(0);
        	pad.back()->SetTopMargin(0);

    		if(x==0){
				pad_Axis.push_back(new TPad(Form("PadAxisLeft_%i_%i",x,y),"",(Float_t)(otstupLeft - coeff*l + l*x),     (Float_t)( otstupDown + h*(y-1)), 
        													            	 (Float_t)(otstupLeft - coeff*l + l*(x+1)), (Float_t)( otstupDown + h*y), 0,0,0));
				pad_Axis.back()->SetFrameBorderMode(1);
        		pad_Axis.back()->SetLeftMargin(coeff);
        		pad_Axis.back()->SetRightMargin(0);
        		pad_Axis.back()->SetBottomMargin(0);
        		pad_Axis.back()->SetTopMargin(0);

            }

            if(y==1){
    			//std::cout<<(Float_t)(otstupLeft + l*x)<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*y)<<"\t"<< (Float_t)(otstupLeft + l*(x+1))<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*(y-1))<<"\n\n";
				pad_Axis.push_back(new TPad(Form("PadAxisDown_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown - coeff*h + h*(y-1)), 
        										  		    		    	 (Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown - coeff*h + h*y), 0,0,0));
				pad_Axis.back()->SetFrameBorderMode(1);
        		pad_Axis.back()->SetLeftMargin(0.0);
        		pad_Axis.back()->SetRightMargin(0);
        		pad_Axis.back()->SetBottomMargin(coeff);
        		pad_Axis.back()->SetTopMargin(0);

            }

    	}
    	
    }

    TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.0, otstupLeft - coeff*l, 1.0 - otstupUp, 0, 0, 0);
	TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, otstupDown - coeff*h, 0, 0, 0);

	pad_AxisY->SetFrameBorderMode(0);
    pad_AxisY->SetRightMargin(0.23);
    pad_AxisY->SetBottomMargin(0);
    
    pad_AxisX->SetTopMargin(0);
    pad_AxisX->SetFrameBorderMode(0);
    pad_AxisX->SetBottomMargin(0);

    pad_AxisY->Draw();
    pad_AxisX->Draw();
	
	for(Int_t i=0; i<pad_Axis.size(); i++){
    	canvas->cd();
    	pad_Axis[i]->Draw();
    	pad_Axis[i]->cd();
    	hPadAxis->Draw();
    }
    canvas->cd();
	for(Int_t i=0; i<pad.size(); i++){
    	pad[i]->Draw();
    }

    pad_AxisY->cd();
	for(Int_t i=0; i<labels.size(); i++){
		if(strncmp(labels[i].place, "Y",1)==0 ){
			_textInfo.SetTextSize(labels[i].size);
    		_textInfo.SetTextAngle(labels[i].angle);
			_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
		}
	}

	pad_AxisX->cd();
	for(Int_t i=0; i<labels.size(); i++){
		if(strncmp(labels[i].place, "X",1)==0 ){
			_textInfo.SetTextSize(labels[i].size);
    		_textInfo.SetTextAngle(labels[i].angle);
			_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
		}
	}

    Int_t p=0;
    std::cout<<pad.size()<<"\n";
    for(Int_t i=0; i<pad.size(); i++){
    	pad[i]->cd();
    	hPadAxis->Draw();
    	_textInfo.SetTextSize(size_bukva);
    	_textInfo.DrawLatex(hPadAxis->GetXaxis()->GetXmin() + 0.86*(hPadAxis->GetXaxis()->GetXmax() - hPadAxis->GetXaxis()->GetXmin()),
    						hPadAxis->GetYaxis()->GetXmin() + 0.85*(hPadAxis->GetYaxis()->GetXmax() - hPadAxis->GetYaxis()->GetXmin()), Form("#font[42]{#scale[1.0]{(%i)}}", i+1));
    	for(Int_t iline=0; iline<Line_vector.size(); iline++){
    		if(Line_vector_DrawPad[iline]==i || Line_vector_DrawPad[iline]==-1){
    			Line_vector[iline]->Draw("same");
    		}
    	}
    }

    std::cout<<DrawResult.size()<<"goog\n";
    for(Int_t i=0; i<(int)DrawResult.size(); i++){
    	
    	Int_t currentPad = DrawResult[i]->GetNumberPadDraw();
		pad[currentPad]->cd();
		
		if(DrawResult[i]->GetKeyGraph()==true){
			(DrawResult[i]->GetDrawObjectGraph())->DrawClone("Psame");
			//Set legend
			if(legend_flag==true && DrawResult[i]->GetKeyLegend()==true){
				legend->AddEntry( DrawResult[i]->GetDrawObjectGraph(), DrawResult[i]->GetLegendText(), "p");
			}
		}
		
		if(DrawResult[i]->GetKeyFun()==true){
			(DrawResult[i]->GetDrawObjectTF1())->Draw("same");
			//Set legend
    		if(legend_flag==true && DrawResult[i]->GetKeyLegend()==true){
				legend->AddEntry(DrawResult[i]->GetDrawObjectTF1(), DrawResult[i]->GetLegendText(), "l");
    		}
		}
		
		// draw labes
		for(Int_t il=0; il<labels.size(); il++){
				if(strncmp(labels[il].place, "Pad",3)==0 && labels[il].NumberPad == currentPad){
				_textInfo.SetTextSize(labels[il].size);
        		_textInfo.SetTextAngle(labels[il].angle);
				_textInfo.DrawLatex(labels[il].x, labels[il].y, (labels[il].text).Data() );
			}
		}	
    }
    
    if(legend_flag==true){
    	pad[legend_nPad]->cd();
        legend->Draw("same");
    }
    std::cout<<"goog\n";

    return canvas;

}

TCanvas* Draw_Picture_new::CanvasWithRatio(Int_t xCanvas, Int_t yCanvas, Int_t column, std::vector<MyDrawObject*> DrawResult, Float_t size_bukva){


	TString bukva[] = {"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)","(m)","(n)","(o)","(p)", "(q)","(r)","(s)","(t)","(u)","(v)"};
	//Float_t
    Int_t line=2;
    Double_t coeff=0.2;
	Float_t xx = xCanvas * (column + 0.4 + 0.05);
	Float_t yy = yCanvas * (line + 0.4 + 0.05);

	Float_t otstupUp = 0.05 * yCanvas / yy;
	Float_t otstupDown = 0.4 * yCanvas / yy;
	Float_t otstupLeft = 0.4 * xCanvas / xx;
	Float_t otstupRight = 0.05 * xCanvas / xx;
	Float_t h = (1.0 - otstupUp - otstupDown )/line;
	Float_t l = (1.0 - otstupRight - otstupLeft )/column;
	//TCanvas *canvas;

	//TCanvas *canvas = new TCanvas("canvas","plot",xCanvas*column,yCanvas*line);
	TCanvas *canvas = new TCanvas("canvas","plot", xx, yy);
    canvas->GetFrame()->SetFillColor(21);
    canvas->GetFrame()->SetBorderSize(115);
    canvas->cd();
    
    std::vector<TPad *> pad;
    std::vector<TPad *> pad_Axis;
    std::vector<int> UpDown;

    Int_t n=0;

    for(Int_t y=line; y>0; y--){
    	Float_t UU = (y==2)?1.0:0.7;
    	Float_t DD = (y==2)?0.7:2.0;
    	for(Int_t x=0; x<column; x++){

		  	//std::cout<<0.1+l*x<<"\t"<< 1.0-h*y<<"\t"<< 0.1+l*(x+1)<<"\t"<< 1.0-h*(y-1)<<"\n\n";
    		pad.push_back(new TPad(Form("Pad_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown + DD*h*(y-1)), 
    														(Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown + UU*h*y), 0,0,0));
        	pad.back()->SetFrameBorderMode(1);
        	pad.back()->SetLeftMargin(0);
        	pad.back()->SetRightMargin(0);
        	pad.back()->SetBottomMargin(0);
        	pad.back()->SetTopMargin(0);

    		if(x==0){
				pad_Axis.push_back(new TPad(Form("PadAxisLeft_%i_%i",x,y),"",(Float_t)(otstupLeft - coeff*l + l*x),     (Float_t)( otstupDown  + DD*h*(y-1)), 
        													            	 (Float_t)(otstupLeft - coeff*l + l*(x+1)), (Float_t)( otstupDown  + UU*h*y), 0,0,0));
				pad_Axis.back()->SetFrameBorderMode(1);
        		pad_Axis.back()->SetLeftMargin(coeff);
        		pad_Axis.back()->SetRightMargin(0);
        		pad_Axis.back()->SetBottomMargin(0);
        		pad_Axis.back()->SetTopMargin(0);

            }

            if(y==1){
    			//std::cout<<(Float_t)(otstupLeft + l*x)<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*y)<<"\t"<< (Float_t)(otstupLeft + l*(x+1))<<"\t"<< (Float_t)( 1.0 - otstupUp - otstupDown - h*(y-1))<<"\n\n";
				pad_Axis.push_back(new TPad(Form("PadAxisDown_%i_%i",x,y),"",(Float_t)(otstupLeft + l*x),     (Float_t)( otstupDown - coeff*h + h*(y-1)), 
        										  		    		    	 (Float_t)(otstupLeft + l*(x+1)), (Float_t)( otstupDown - coeff*h + h*y), 0,0,0));
				pad_Axis.back()->SetFrameBorderMode(1);
        		pad_Axis.back()->SetLeftMargin(0.0);
        		pad_Axis.back()->SetRightMargin(0);
        		pad_Axis.back()->SetBottomMargin(coeff);
        		pad_Axis.back()->SetTopMargin(0);

            }

    	}
    	
    }

    TPad *pad_AxisY = new TPad("pad_AxisY", "", 0.0, 0.0, otstupLeft - coeff*l, 1.0 - otstupUp, 0, 0, 0);
	TPad *pad_AxisX = new TPad("pad_AxisX", "", 0.0, 0.0, 1.0, otstupDown - coeff*h, 0, 0, 0);

	pad_AxisY->SetFrameBorderMode(0);
    pad_AxisY->SetRightMargin(0.23);
    pad_AxisY->SetBottomMargin(0);
    
    pad_AxisX->SetTopMargin(0);
    pad_AxisX->SetFrameBorderMode(0);
    pad_AxisX->SetBottomMargin(0);

    pad_AxisY->Draw();
    pad_AxisX->Draw();
	
	for(Int_t i=0; i<pad_Axis.size(); i++){
    	canvas->cd();
    	pad_Axis[i]->Draw();
    	pad_Axis[i]->cd();
    	if(i==0)hPadAxis->Draw();
    	if(i!=0)hPadAxis2->Draw();
    }
    canvas->cd();
	for(Int_t i=0; i<pad.size(); i++){
    	pad[i]->Draw();
    }

    pad_AxisY->cd();
	for(Int_t i=0; i<labels.size(); i++){
		if(strncmp(labels[i].place, "Y",1)==0 ){
			_textInfo.SetTextSize(labels[i].size);
    		_textInfo.SetTextAngle(labels[i].angle);
			_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
		}
	}

	pad_AxisX->cd();
	for(Int_t i=0; i<labels.size(); i++){
		if(strncmp(labels[i].place, "X",1)==0 ){
			_textInfo.SetTextSize(labels[i].size);
    		_textInfo.SetTextAngle(labels[i].angle);
			_textInfo.DrawLatex(labels[i].x, labels[i].y, (labels[i].text).Data() );
		}
	}

    Int_t p=0;
    std::cout<<pad.size()<<"\n";
    for(Int_t i=0; i<pad.size(); i++){
    	pad[i]->cd();
    	if(i<column)hPadAxis->Draw();
    	if(i>=column)hPadAxis2->Draw();
    	_textInfo.SetTextSize(size_bukva);
    	_textInfo.DrawLatex(hPadAxis->GetXaxis()->GetXmin() + 0.86*(hPadAxis->GetXaxis()->GetXmax() - hPadAxis->GetXaxis()->GetXmin()),
    						hPadAxis->GetYaxis()->GetXmin() + 0.85*(hPadAxis->GetYaxis()->GetXmax() - hPadAxis->GetYaxis()->GetXmin()), Form("#font[42]{#scale[1.0]{(%i)}}", i+1));
    	for(Int_t iline=0; iline<Line_vector.size(); iline++){
    		if(Line_vector_DrawPad[iline]==i || Line_vector_DrawPad[iline]==-1){
    			Line_vector[iline]->Draw("same");
    		}
    	}
    }

    std::cout<<DrawResult.size()<<"goog\n";
    for(Int_t i=0; i<(int)DrawResult.size(); i++){
    	
    	Int_t currentPad = DrawResult[i]->GetNumberPadDraw();
		pad[currentPad]->cd();
		
		if(DrawResult[i]->GetKeyGraph()==true){
			(DrawResult[i]->GetDrawObjectGraph())->DrawClone("Psame");
			//Set legend
			if(legend_flag==true && DrawResult[i]->GetKeyLegend()==true){
				legend->AddEntry( DrawResult[i]->GetDrawObjectGraph(), DrawResult[i]->GetLegendText(), "p");
			}
		}
		
		if(DrawResult[i]->GetKeyFun()==true){
			(DrawResult[i]->GetDrawObjectTF1())->Draw("same");
			//Set legend
    		if(legend_flag==true && DrawResult[i]->GetKeyLegend()==true){
				legend->AddEntry(DrawResult[i]->GetDrawObjectTF1(), DrawResult[i]->GetLegendText(), "l");
    		}
		}
		
		// draw labes
		for(Int_t il=0; il<labels.size(); il++){
				if(strncmp(labels[il].place, "Pad",3)==0 && labels[il].NumberPad == currentPad){
				_textInfo.SetTextSize(labels[il].size);
        		_textInfo.SetTextAngle(labels[il].angle);
				_textInfo.DrawLatex(labels[il].x, labels[il].y, (labels[il].text).Data() );
			}
		}	
    }
    
    if(legend_flag==true){
    	pad[legend_nPad]->cd();
        legend->Draw("same");
    }
    std::cout<<"goog\n";

    return canvas;

}