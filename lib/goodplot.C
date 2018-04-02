#include <iostream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TStyle.h"
#include "TF1.h"
#include "TPad.h"
#include "TProfile.h"

void goodplot(const std::string& inFileName, const std::string& ch1)
{
	//string ch1 = "15";
	string ch2 = "9";
	
	string time_ch1 = "LP2_10";
	string time_ch2 = "gaus_mean";

	float amp_low = 200.0;//ch16
	if(ch1=="15" || ch1=="14") amp_low = 100.0;//ch15
	float amp_high = 400.0;
	
	float x_low = 5.0;
	float x_high = 30.0;
	float y_low = 0.0;
	float y_high = 25.0;

	//for ch16
	
	float x_tile_low = 17.0;
	float x_tile_low_fit = 17.0;
	float x_tile_high = 26.0;
	float x_tile_high_fit = 26.0;
	float y_tile_low = 3.0;
	float y_tile_low_fit = 3.0;
	float y_tile_high = 13.0;
	float y_tile_high_fit = 13.0;
	
	//for ch15
	if(ch1=="15" || ch1=="14")
	{	
	x_tile_low = 15.0;
	x_tile_low_fit = 17.0;
	x_tile_high = 25.0;
	x_tile_high_fit = 25.0;
	y_tile_low = 3.0;
	y_tile_low_fit = 3.0;
	y_tile_high = 13.0;
	y_tile_high_fit = 13.0;
	}

	string inputDir = "/eos/uscms/store/group/cmstestbeam/BTL/March2018/OTSDAQ/CMSTiming/RECO/V3/";
	string plotDir = "plots/";

	mkdir(plotDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

	cout<<"plotting run:"<<inFileName<<endl;
	
	string cut_noPos = "amp["+ch1+"]>"+std::to_string(amp_low)+" && amp["+ch1+"] < "+std::to_string(amp_high)+" && amp["+ch2+"]>40.0 && amp["+ch2+"]<400.0 && gaus_mean["+ch2+"]>1 && gaus_mean["+ch2+"]<100";

        string cut = "amp["+ch1+"]>"+std::to_string(amp_low)+" && amp["+ch1+"] < "+std::to_string(amp_high)+" && amp["+ch2+"]>40.0 && amp["+ch2+"]<400.0 && x_dut[0] > "+std::to_string(x_tile_low)+" && x_dut[0]<"+std::to_string(x_tile_high)+" && y_dut[0]>"+std::to_string(y_tile_low)+" && y_dut[0]<"+std::to_string(y_tile_high)+"&& gaus_mean["+ch2+"]>1 && gaus_mean["+ch2+"]<100";


	//ch1 and ch2 very basic cut
	string cut1 = "amp["+ch1+"]>10";
	string cut2 = "amp["+ch2+"]>10";
 	
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(111);
	gStyle->SetPalette(1);

	TFile* inFile = TFile::Open((inputDir+"/DataCMSVMETiming_Run"+inFileName+".root").c_str(),"READ");
	TTree* tree = (TTree*)( inFile->Get("pulse") );

	Int_t NEntries = tree->GetEntries();
        Int_t NEntries_cut_noPos = tree->GetEntries(cut_noPos.c_str());
        Int_t NEntries_cut = tree->GetEntries(cut.c_str());
        Int_t NEntries_cut1 = tree->GetEntries(cut1.c_str());
        Int_t NEntries_cut2 = tree->GetEntries(cut2.c_str());

        cout<<"[tot]: "<<NEntries<<"   [cut_noPos]: "<<NEntries_cut_noPos<<"   [cut]: "<<NEntries_cut<<"   [cut"+ch1+"]: "<<NEntries_cut1<<"   [cut"+ch2+"]"<<NEntries_cut2<<endl;

	float axisTitleSizeX = 0.06;
	float axisTitleSizeY = 0.05;
	float axisTitleOffsetX = 0.9;
	float axisTitleOffsetY = 1.2;

	float axisTitleSizeRatioX   = 0.18;
	float axisLabelSizeRatioX   = 0.12;
	float axisTitleOffsetRatioX = 0.94;
	float axisTitleSizeRatioY   = 0.15;
	float axisLabelSizeRatioY   = 0.108;
	float axisTitleOffsetRatioY = 0.32;

	float leftMargin   = 0.12;
	float rightMargin  = 0.14;
	float topMargin    = 0.07;
	float bottomMargin = 0.12;


	TCanvas *myC = new TCanvas( "myC", "myC", 200, 10, 800, 600 );
	myC->SetHighLightColor(2);
	myC->SetFillColor(0);
	myC->SetBorderMode(0);
	myC->SetBorderSize(2);
	myC->SetLeftMargin( leftMargin );
	myC->SetRightMargin( rightMargin );
	myC->SetTopMargin( topMargin );
	myC->SetBottomMargin( bottomMargin );
	myC->SetFrameBorderMode(0);
	myC->SetFrameBorderMode(0);



	// beam spot selection
	myC->SetGridy(1);
	myC->SetGridx(1);
	// ch1
	TH2F * h2_beamSpotch1 = new TH2F("h2_beamSpotch1","h2_beamSpotch1",100,x_low, x_high,100,y_low, y_high);
	TH2F * h2_beamSpotch1Num = new TH2F("h2_beamSpotch1Num","h2_beamSpotch1Num",100,x_low, x_high,100,y_low,y_high);
	tree->Draw("y_dut[0]:x_dut[0]>>h2_beamSpotch1Num",cut1.c_str());
	tree->Draw("y_dut[0]:x_dut[0]>>h2_beamSpotch1",("amp["+ch1+"]*("+cut1+")").c_str());
	h2_beamSpotch1->Divide(h2_beamSpotch1Num);
	h2_beamSpotch1->Draw("colz");
	
	h2_beamSpotch1->GetXaxis()->SetTitle("beam position X [mm]");
	h2_beamSpotch1->GetYaxis()->SetTitle("beam position Y [mm]");
	h2_beamSpotch1->GetZaxis()->SetTitle("amplitude [mV]");
	h2_beamSpotch1->SetTitle("");
	h2_beamSpotch1->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h2_beamSpotch1->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h2_beamSpotch1->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h2_beamSpotch1->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
 
	myC->SaveAs((plotDir+"/Run"+inFileName+"_beamSpot_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_beamSpot_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_beamSpot_ch"+ch1+".C").c_str());



	// ch2
	TH2F * h2_beamSpotch2 = new TH2F("h2_beamSpotch2","h2_beamSpotch2",100,x_low, x_high,100,y_low,y_high);
	TH2F * h2_beamSpotch2Num = new TH2F("h2_beamSpotch2Num","h2_beamSpotch2Num",100,x_low, x_high,100,y_low,y_high);
	tree->Draw("y_dut[0]:x_dut[0]>>h2_beamSpotch2Num",cut2.c_str());
	tree->Draw("y_dut[0]:x_dut[0]>>h2_beamSpotch2",("amp["+ch2+"]*("+cut2+")").c_str());
	h2_beamSpotch2->Divide(h2_beamSpotch2Num);
	h2_beamSpotch2->Draw("colz");
	
	h2_beamSpotch2->GetXaxis()->SetTitle("beam position X [mm]");
	h2_beamSpotch2->GetYaxis()->SetTitle("beam position Y [mm]");
	h2_beamSpotch2->GetZaxis()->SetTitle("amplitude [mV]");
	h2_beamSpotch2->SetTitle("");
	h2_beamSpotch2->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h2_beamSpotch2->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h2_beamSpotch2->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h2_beamSpotch2->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
 
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_beamSpot_ch"+ch2+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_beamSpot_ch"+ch2+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_beamSpot_ch"+ch2+".C").c_str());


	//amplitude
	myC->SetGridy(0);
	myC->SetGridx(0);
	TH1F * h_ch1E = new TH1F("h_ch1E","h_ch1E",400,15, 2000.0);
	TH1F * h_ch2E = new TH1F("h_ch2E","h_ch2E",400,15,2000.0);

	tree->Draw(("amp["+ch1+"]>>h_ch1E").c_str(), cut1.c_str());	
	tree->Draw(("amp["+ch2+"]>>h_ch2E").c_str(), cut2.c_str());	
	
	h_ch1E->SetTitle("");
	float maxY_E = 1.2*std::max(h_ch1E->GetMaximum(), h_ch2E->GetMaximum());

	float maxX_ch1E = h_ch1E->GetBinCenter(h_ch1E->GetMaximumBin());
	float highch1E=h_ch1E->GetBinCenter(h_ch1E->FindLastBinAbove(int(0.1*h_ch1E->GetMaximum())));
	float maxX_ch1 = std::max(2.0*maxX_ch1E, highch1E+200.0);

	float maxX_ch2E = h_ch2E->GetBinCenter(h_ch2E->GetMaximumBin());
	float highch2E=h_ch2E->GetBinCenter(h_ch2E->FindLastBinAbove(int(0.1*h_ch2E->GetMaximum())));
	float maxX_ch2 = std::max(2.0*maxX_ch2E, highch2E+200.0);

	float maxX_E = std::max(maxX_ch1, maxX_ch2);

	h_ch1E->SetMarkerStyle( 20 );
	h_ch1E->SetMarkerColor( 2 );
	h_ch1E->SetLineColor( 2 );
	h_ch2E->SetMarkerStyle( 20 );
	h_ch2E->SetMarkerColor( 4 );
	h_ch2E->SetLineColor( 4 );
	h_ch1E->GetXaxis()->SetTitle("Amplitude [mV]");
	h_ch1E->GetYaxis()->SetTitle("Events");
	h_ch1E->SetTitle("");
	h_ch1E->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h_ch1E->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h_ch1E->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h_ch1E->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	h_ch1E->GetYaxis()->SetRangeUser(0.0,maxY_E  );
	h_ch1E->GetXaxis()->SetRangeUser(0.0,maxX_E );
        h_ch1E->Draw("histE");
	
	h_ch2E->Draw("samehistE");
	
	TLegend * leg_E  = new TLegend (0.65,0.7,0.85,0.9);
	leg_E->SetBorderSize(0);
	leg_E->SetTextSize(0.03);
	leg_E->SetLineColor(1);
	leg_E->SetLineStyle(1);
	leg_E->SetLineWidth(1);
	leg_E->SetFillColor(0);
	leg_E->SetFillStyle(1001);
		
	leg_E->AddEntry(h_ch1E, ("ch "+ch1).c_str(), "lp");
	leg_E->AddEntry(h_ch2E, ("ch "+ch2).c_str(), "lp");
	leg_E->Draw();

	myC->SaveAs((plotDir+"/Run"+inFileName+"_amplitude_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_amplitude_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_amplitude_ch"+ch1+".C").c_str());


	// time resolution
	TH1F * h_deltaT = new TH1F("h_deltaT","h_deltaT",4000,-50.0, 50.0);
	tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"]>>h_deltaT").c_str(), cut.c_str());	
	h_deltaT->SetTitle("");
	h_deltaT->SetMarkerStyle( 20 );
	h_deltaT->SetMarkerColor( 1 );
	h_deltaT->SetLineColor( 1 );
	h_deltaT->GetXaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	h_deltaT->GetYaxis()->SetTitle("Events");
	h_deltaT->SetTitle("");
	h_deltaT->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h_deltaT->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h_deltaT->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h_deltaT->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	float maxY_t = h_deltaT->GetMaximum();
	float maxX_t = h_deltaT->GetBinCenter(h_deltaT->GetMaximumBin());
       	h_deltaT->GetXaxis()->SetRangeUser(maxX_t-0.6, maxX_t+0.6);
        h_deltaT->Draw("E");

	float highDeltaT=h_deltaT->GetBinCenter(h_deltaT->FindLastBinAbove(int(0.1*maxY_t)));
	float lowDeltaT=h_deltaT->GetBinCenter(h_deltaT->FindFirstBinAbove(int(0.1*maxY_t)));
	if(highDeltaT - maxX_t > 2.0*(maxX_t-lowDeltaT)) highDeltaT = maxX_t + (maxX_t - lowDeltaT);
	TF1 * tf1_gaus = new TF1("tf1_gaus","gaus", maxX_t - 1.0, maxX_t + 1.0);
	tf1_gaus->SetParameter(1, h_deltaT->GetMean());
	h_deltaT->Fit("tf1_gaus","","",lowDeltaT, highDeltaT);

	gPad->Modified();
	gPad->Update();

	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_ch"+ch1+".C").c_str());

	//amplitude dependency of the deltaT
	myC->SetGridy(1);
	myC->SetGridx(1);
	TH2F * h2_deltaT_vs_amp = new TH2F("h2_deltaT_vs_amp","h2_deltaT_vs_amp", 100, amp_low, amp_high, 100, maxX_t-0.6, maxX_t+0.6);
	tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"]:amp["+ch1+"]>>h2_deltaT_vs_amp").c_str(),cut_noPos.c_str());
	h2_deltaT_vs_amp->GetXaxis()->SetTitle("amplitude [mV]");
	h2_deltaT_vs_amp->GetYaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	h2_deltaT_vs_amp->SetTitle("");
	h2_deltaT_vs_amp->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h2_deltaT_vs_amp->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h2_deltaT_vs_amp->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h2_deltaT_vs_amp->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	
	h2_deltaT_vs_amp->Draw();
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_vs_amp_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_vs_amp_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_vs_amp_ch"+ch1+".C").c_str());

	TProfile * p_deltaT_vs_amp = h2_deltaT_vs_amp->ProfileX(); //new ("p_deltaT_vs_amp","p_deltaT_vs_amp",100, x_low, x_high);
	p_deltaT_vs_amp->GetXaxis()->SetTitle("amplitude [mV]");
	p_deltaT_vs_amp->GetYaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	p_deltaT_vs_amp->SetTitle("");
	
	p_deltaT_vs_amp->GetXaxis()->SetTitleSize( axisTitleSizeX );
	p_deltaT_vs_amp->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	p_deltaT_vs_amp->GetYaxis()->SetTitleSize( axisTitleSizeY );
	p_deltaT_vs_amp->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
       	p_deltaT_vs_amp->GetYaxis()->SetRangeUser(maxX_t-0.2, maxX_t+0.3);
       	p_deltaT_vs_amp->GetXaxis()->SetRangeUser(amp_low, amp_high);
        p_deltaT_vs_amp->SetMarkerStyle( 20 );
        p_deltaT_vs_amp->SetMarkerColor( 1 );
        p_deltaT_vs_amp->SetLineColor( 1 );
	p_deltaT_vs_amp->Draw();


	TF1 * tf1_pol1_amp = new TF1("tf1_pol1_amp","pol1", amp_low, amp_high);
	p_deltaT_vs_amp->Fit("tf1_pol1_amp","","",amp_low, amp_high);
	float amp_cor_p0 = tf1_pol1_amp->GetParameter(0);
	float amp_cor_p1 = tf1_pol1_amp->GetParameter(1);

	cout<<"p0: "<<amp_cor_p0<<" ;  p1: "<<amp_cor_p1<<endl;

	gPad->Modified();
	gPad->Update();
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_vs_amp_profile_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_vs_amp_profile_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_vs_amp_profile_ch"+ch1+".C").c_str());

	//time walk correction
	myC->SetGridy(0);
	myC->SetGridx(0);
	
	TH1F * h_deltaT_TWcorr = new TH1F("h_deltaT_TWcorr","h_deltaT_TWcorr",4000,-50.0, 50.0);
	tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"]- ("+std::to_string(amp_cor_p1)+"*amp["+ch1+"])>>h_deltaT_TWcorr").c_str(), cut.c_str());	

	h_deltaT_TWcorr->SetTitle("");
	h_deltaT_TWcorr->SetMarkerStyle( 20 );
	h_deltaT_TWcorr->SetMarkerColor( 1 );
	h_deltaT_TWcorr->SetLineColor( 1 );
	h_deltaT_TWcorr->GetXaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	h_deltaT_TWcorr->GetYaxis()->SetTitle("Events");
	h_deltaT_TWcorr->SetTitle("");
	h_deltaT_TWcorr->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h_deltaT_TWcorr->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h_deltaT_TWcorr->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h_deltaT_TWcorr->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	float maxY_t_TWcorr = h_deltaT_TWcorr->GetMaximum();
	float maxX_t_TWcorr = h_deltaT_TWcorr->GetBinCenter(h_deltaT_TWcorr->GetMaximumBin());
       	h_deltaT_TWcorr->GetXaxis()->SetRangeUser(maxX_t_TWcorr-0.6, maxX_t_TWcorr+0.6);
        h_deltaT_TWcorr->Draw("E");

	float highDeltaT_TWcorr=h_deltaT_TWcorr->GetBinCenter(h_deltaT_TWcorr->FindLastBinAbove(int(0.1*maxY_t)));
	float lowDeltaT_TWcorr=h_deltaT_TWcorr->GetBinCenter(h_deltaT_TWcorr->FindFirstBinAbove(int(0.1*maxY_t)));
	if(highDeltaT_TWcorr - maxX_t_TWcorr > 2.0*(maxX_t_TWcorr-lowDeltaT_TWcorr)) highDeltaT_TWcorr = maxX_t_TWcorr + (maxX_t_TWcorr - lowDeltaT_TWcorr);
	TF1 * tf1_gaus_TWcorr = new TF1("tf1_gaus_TWcorr","gaus", maxX_t_TWcorr - 1.0, maxX_t_TWcorr + 1.0);
	tf1_gaus_TWcorr->SetParameter(1, h_deltaT_TWcorr->GetMean());
	h_deltaT_TWcorr->Fit("tf1_gaus_TWcorr","","",lowDeltaT_TWcorr, highDeltaT_TWcorr);

	gPad->Modified();
	gPad->Update();

	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_TWcorr_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_TWcorr_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_TWcorr_ch"+ch1+".C").c_str());



	// impact point dependency of the deltaT
	myC->SetGridy(1);
	myC->SetGridx(1);
	TH2F * h2_deltaT_vs_x = new TH2F("h2_deltaT_vs_x","h2_deltaT_vs_x", 100, x_low, x_high, 100, maxX_t-0.6, maxX_t+0.6);
	//tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"]:x_dut[0]>>h2_deltaT_vs_x").c_str(),cut_noPos.c_str());
	tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"] - ("+std::to_string(amp_cor_p1)+"*amp["+ch1+"]) : x_dut[0]>>h2_deltaT_vs_x").c_str(),cut_noPos.c_str());
	h2_deltaT_vs_x->GetXaxis()->SetTitle("beam position X [mm]");
	h2_deltaT_vs_x->GetYaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	h2_deltaT_vs_x->SetTitle("");
	h2_deltaT_vs_x->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h2_deltaT_vs_x->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h2_deltaT_vs_x->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h2_deltaT_vs_x->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	
	h2_deltaT_vs_x->Draw();
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamX_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamX_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamX_ch"+ch1+".C").c_str());

	TProfile * p_deltaT_vs_x = h2_deltaT_vs_x->ProfileX(); //new ("p_deltaT_vs_x","p_deltaT_vs_x",100, x_low, x_high);
	p_deltaT_vs_x->GetXaxis()->SetTitle("beam position X [mm]");
	p_deltaT_vs_x->GetYaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	p_deltaT_vs_x->SetTitle("");
	
	p_deltaT_vs_x->GetXaxis()->SetTitleSize( axisTitleSizeX );
	p_deltaT_vs_x->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	p_deltaT_vs_x->GetYaxis()->SetTitleSize( axisTitleSizeY );
	p_deltaT_vs_x->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
       	//p_deltaT_vs_x->GetYaxis()->SetRangeUser(maxX_t-0.2, maxX_t+0.3);
       	p_deltaT_vs_x->GetYaxis()->SetRangeUser(maxX_t_TWcorr-0.2, maxX_t_TWcorr+0.3);
       	p_deltaT_vs_x->GetXaxis()->SetRangeUser(x_tile_low, x_tile_high);
        p_deltaT_vs_x->SetMarkerStyle( 20 );
        p_deltaT_vs_x->SetMarkerColor( 1 );
        p_deltaT_vs_x->SetLineColor( 1 );
	p_deltaT_vs_x->Draw();

	TF1 * tf1_pol2_x = new TF1("tf1_pol2_x","pol2", x_tile_low-1, x_tile_high+1.0);
	p_deltaT_vs_x->Fit("tf1_pol2_x","","",x_tile_low_fit, x_tile_high_fit);
	float x_cor_p0 = tf1_pol2_x->GetParameter(0);
	float x_cor_p1 = tf1_pol2_x->GetParameter(1);
	float x_cor_p2 = tf1_pol2_x->GetParameter(2);

	cout<<"p0: "<<x_cor_p0<<" ;  p1: "<<x_cor_p1<<" ;  p2: "<<x_cor_p2<<endl;

	gPad->Modified();
	gPad->Update();
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamX_profile_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamX_profile_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamX_profile_ch"+ch1+".C").c_str());


	TH2F * h2_deltaT_vs_y = new TH2F("h2_deltaT_vs_y","h2_deltaT_vs_y", 100, y_low, y_high, 100, maxX_t-0.6, maxX_t+0.6);
	tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"] - ("+std::to_string(amp_cor_p1)+"*amp["+ch1+"]) : y_dut[0]>>h2_deltaT_vs_y").c_str(),cut_noPos.c_str());
	h2_deltaT_vs_y->GetXaxis()->SetTitle("beam position X [mm]");
	h2_deltaT_vs_y->GetYaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	h2_deltaT_vs_y->SetTitle("");
	h2_deltaT_vs_y->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h2_deltaT_vs_y->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h2_deltaT_vs_y->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h2_deltaT_vs_y->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	
	h2_deltaT_vs_y->Draw();
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamY_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamY_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamY_ch"+ch1+".C").c_str());

	TProfile * p_deltaT_vs_y = h2_deltaT_vs_y->ProfileX(); //new ("p_deltaT_vs_y","p_deltaT_vs_y",100, y_low, y_high);
	p_deltaT_vs_y->GetXaxis()->SetTitle("beam position Y [mm]");
	p_deltaT_vs_y->GetYaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	p_deltaT_vs_y->SetTitle("");
	
	p_deltaT_vs_y->GetXaxis()->SetTitleSize( axisTitleSizeX );
	p_deltaT_vs_y->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	p_deltaT_vs_y->GetYaxis()->SetTitleSize( axisTitleSizeY );
	p_deltaT_vs_y->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
       	p_deltaT_vs_y->GetYaxis()->SetRangeUser(maxX_t_TWcorr-0.2, maxX_t_TWcorr+0.3);
       	p_deltaT_vs_y->GetXaxis()->SetRangeUser(y_tile_low, y_tile_high);
        p_deltaT_vs_y->SetMarkerStyle( 20 );
        p_deltaT_vs_y->SetMarkerColor( 1 );
        p_deltaT_vs_y->SetLineColor( 1 );
	p_deltaT_vs_y->Draw();

	TF1 * tf1_pol2_y = new TF1("tf1_pol2_y","pol2", y_tile_low-1, y_tile_high+1.0);
	p_deltaT_vs_y->Fit("tf1_pol2_y","","",y_tile_low_fit, y_tile_high_fit);
	float y_cor_p0 = tf1_pol2_y->GetParameter(0);
	float y_cor_p1 = tf1_pol2_y->GetParameter(1);
	float y_cor_p2 = tf1_pol2_y->GetParameter(2);

	cout<<"p0: "<<y_cor_p0<<" ;  p1: "<<y_cor_p1<<" ;  p2: "<<y_cor_p2<<endl;

	gPad->Modified();
	gPad->Update();
	
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamY_profile_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamY_profile_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_afterTW_vs_beamY_profile_ch"+ch1+".C").c_str());

	//TW+impact point correction
	myC->SetGridy(0);
	myC->SetGridx(0);
	
	TH1F * h_deltaT_TWXYcorr = new TH1F("h_deltaT_TWXYcorr","h_deltaT_TWXYcorr",4000,-50.0, 50.0);
	tree->Draw((time_ch1+"["+ch1+"]-"+time_ch2+"["+ch2+"]- ("+std::to_string(amp_cor_p1)+"*amp["+ch1+"] + "+std::to_string(x_cor_p1)+"*x_dut[0] + "+std::to_string(x_cor_p2)+"*x_dut[0]*x_dut[0] + "+std::to_string(y_cor_p1)+"*y_dut[0] + "+std::to_string(y_cor_p2)+"*y_dut[0]*y_dut[0]"+")>>h_deltaT_TWXYcorr").c_str(), cut.c_str());	

	h_deltaT_TWXYcorr->SetTitle("");
	h_deltaT_TWXYcorr->SetMarkerStyle( 20 );
	h_deltaT_TWXYcorr->SetMarkerColor( 1 );
	h_deltaT_TWXYcorr->SetLineColor( 1 );
	h_deltaT_TWXYcorr->GetXaxis()->SetTitle(("#Delta T (ch"+ch1+", ch"+ch2+") / ns").c_str());
	h_deltaT_TWXYcorr->GetYaxis()->SetTitle("Events");
	h_deltaT_TWXYcorr->SetTitle("");
	h_deltaT_TWXYcorr->GetXaxis()->SetTitleSize( axisTitleSizeX );
	h_deltaT_TWXYcorr->GetXaxis()->SetTitleOffset( axisTitleOffsetX );
	h_deltaT_TWXYcorr->GetYaxis()->SetTitleSize( axisTitleSizeY );
	h_deltaT_TWXYcorr->GetYaxis()->SetTitleOffset( axisTitleOffsetY );
	float maxY_t_TWXYcorr = h_deltaT_TWXYcorr->GetMaximum();
	float maxX_t_TWXYcorr = h_deltaT_TWXYcorr->GetBinCenter(h_deltaT_TWXYcorr->GetMaximumBin());
       	h_deltaT_TWXYcorr->GetXaxis()->SetRangeUser(maxX_t_TWXYcorr-0.6, maxX_t_TWXYcorr+0.6);
        h_deltaT_TWXYcorr->Draw("E");

	float highDeltaT_TWXYcorr=h_deltaT_TWXYcorr->GetBinCenter(h_deltaT_TWXYcorr->FindLastBinAbove(int(0.1*maxY_t)));
	float lowDeltaT_TWXYcorr=h_deltaT_TWXYcorr->GetBinCenter(h_deltaT_TWXYcorr->FindFirstBinAbove(int(0.1*maxY_t)));
	if(highDeltaT_TWXYcorr - maxX_t_TWXYcorr > 2.0*(maxX_t_TWXYcorr-lowDeltaT_TWXYcorr)) highDeltaT_TWXYcorr = maxX_t_TWXYcorr + (maxX_t_TWXYcorr - lowDeltaT_TWXYcorr);
	TF1 * tf1_gaus_TWXYcorr = new TF1("tf1_gaus_TWXYcorr","gaus", maxX_t_TWXYcorr - 1.0, maxX_t_TWXYcorr + 1.0);
	tf1_gaus_TWXYcorr->SetParameter(1, h_deltaT_TWXYcorr->GetMean());
	h_deltaT_TWXYcorr->Fit("tf1_gaus_TWXYcorr","","",lowDeltaT_TWXYcorr, highDeltaT_TWXYcorr);

	gPad->Modified();
	gPad->Update();

	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_TW_then_XYcorr_ch"+ch1+".pdf").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_TW_then_XYcorr_ch"+ch1+".png").c_str());
	myC->SaveAs((plotDir+"/Run"+inFileName+"_deltaT_TW_then_XYcorr_ch"+ch1+".C").c_str());

}
