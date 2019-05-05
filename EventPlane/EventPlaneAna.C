#include "EventPlaneAna.h"

EventPlaneAna::EventPlaneAna():
	N(0),
	ptg(),
	etag(),
	phig()
	{
	TH1::SetDefaultSumw2();
	double ptarr[] = {0, 0.2, 0.5, 0.8, 1.5, 3.0, 6.0};
        ptbin = std::vector<double>(ptarr, ptarr + sizeof(ptarr)/sizeof(ptarr[0]));
        npt = ptbin.size()-1;
	hReso = new TH1F("hReso","hReso",100,-1.1,1.1);
	hpTvnRaw = new TProfile("hpTvnRaw","hpTvnRaw",npt,ptbin.data(),-1.1,1.1);
	hpTvnCorr = new TH1F("hpTvnCorr","hpTvnCorr",npt,ptbin.data());
	nharm = 2;
}

EventPlaneAna::~EventPlaneAna() {
	delete tree;
}

void EventPlaneAna::Init() {
	tree = new TChain("tree");
	for(unsigned int ifile = 0; ifile < filelist.size(); ifile++){
		tree -> Add(filelist[ifile].c_str());
	}
	tree -> SetBranchAddress("n", &N);
	tree -> SetBranchAddress("ptg", ptg);
	tree -> SetBranchAddress("etag", etag);
	tree -> SetBranchAddress("phig", phig);
}

void EventPlaneAna::ProcessEvents() {
	for (int ievent = 0; ievent < tree -> GetEntries(); ievent ++) {
		if(ievent % 1000 == 0) std::cout << "processed " << ievent << " events"<< std::endl;
		tree -> GetEntry(ievent);
		float Qx_A = 0;
		float Qy_A = 0; 		
		float Qx_B = 0;
		float Qy_B = 0; 
		for (int ipart = 0; ipart < N; ipart ++) {
			float phi = phig[ipart];
			float eta = etag[ipart];
			float pt = ptg[ipart];
			if(eta < -1.5){
				Qx_A += cos(nharm * phi);
				Qy_A += sin(nharm * phi);
			}			
			if(eta > 1.5){
				Qx_B += cos(nharm * phi);
				Qy_B += sin(nharm * phi);
			}
		} // first particle loop
		float Psi_A = atan2(Qy_A, Qx_A) / nharm;
		float Psi_B = atan2(Qy_B, Qx_B) / nharm;

		hReso -> Fill(cos(nharm * (Psi_A - Psi_B)));

		for (int ipart = 0; ipart < N; ipart ++) {
			float phi = phig[ipart];
			float eta = etag[ipart];
			float pt = ptg[ipart];
			if(eta > 0 && eta < 1.5){
				hpTvnRaw -> Fill (pt, cos(nharm * (phi - Psi_A)));
			}			
			if(eta > -1.5 && eta < 0){
				hpTvnRaw -> Fill (pt, cos(nharm * (phi - Psi_B)));
			}
		} // second particle loop
	} // event loop
}

void EventPlaneAna::End() {
	float reso = sqrt(hReso -> GetMean());
	std::cout << "Event plane " << nharm <<"th resolution = " << reso << std::endl;
	hpTvnCorr = (TH1F*)hpTvnRaw -> ProjectionX();
	hpTvnCorr -> Scale(1./reso);
}


void EventPlaneAna::Plot() {
	gStyle -> SetOptStat(0);
	TCanvas *c1 = new TCanvas();
	TH1F* hempty = new TH1F("","",12,0,6);
	hempty -> SetMaximum(0.20);
	hempty -> SetTitle("1k events event plane method STEG");
	hempty -> GetXaxis() -> SetTitle("p_{T}");
	hempty -> GetYaxis() -> SetTitle(Form("v_{%d}",nharm));
	hempty -> Draw();
	hpTvnCorr -> SetMarkerStyle(20);
	hpTvnCorr -> SetMarkerSize(0.8);
	hpTvnCorr -> SetMarkerColor(1);
	hpTvnCorr -> Draw("Psame");
        TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        TF1 *V3vsPt = new TF1("V3vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        TF1 *V4vsPt = new TF1("V4vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        V2vsPt -> SetParameters(3.31699,2.35142,3.49188,3.54429,.00005,1.50600);
        V3vsPt -> SetParameters(3.2,2.3,3.4,2.1,.00005,1.4);
        V4vsPt -> SetParameters(4.8,2.1,3.4,2.1,.00005,1.4);
	if (nharm == 2) {
		V2vsPt -> Draw("same");
	}	
	if (nharm == 3) {
		V3vsPt -> Draw("same");
	}
	if (nharm == 4) {
		V4vsPt -> Draw("same");
	}
	TLegend *leg = new TLegend(0.1,0.75,0.45,0.88);
	leg -> SetTextSize(0.05);
	leg -> SetBorderSize(0);
	leg -> SetFillStyle(0);
	leg -> AddEntry(V2vsPt, Form("input v_{%d}",nharm),"l");
	leg -> AddEntry(hpTvnCorr, Form("v_{%d} from EP method",nharm),"p");
	leg -> Draw("same");

	c1 -> Print(Form("V%dvsPT.png",nharm));
}
