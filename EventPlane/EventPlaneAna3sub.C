#include "EventPlaneAna3sub.h"

EventPlaneAna3sub::EventPlaneAna3sub():
	N(0),
	ptg(),
	etag(),
	phig()
	{
	TH1::SetDefaultSumw2();
	double ptarr[] = {0, 0.2, 0.5, 0.8, 1.5, 3.0, 6.0};
        ptbin = std::vector<double>(ptarr, ptarr + sizeof(ptarr)/sizeof(ptarr[0]));
        npt = ptbin.size()-1;
	hReso_AB = new TH1F("hReso_AB","hReso_AB",100,-1.1,1.1);
	hReso_BC = new TH1F("hReso_BC","hReso_BC",100,-1.1,1.1);
	hReso_AC = new TH1F("hReso_AC","hReso_AC",100,-1.1,1.1);
	hpTvnRaw_B = new TProfile("hpTvnRaw_B","hpTvnRaw_B",npt,ptbin.data(),-1.1,1.1);
	hpTvnRaw_C = new TProfile("hpTvnRaw_C","hpTvnRaw_C",npt,ptbin.data(),-1.1,1.1);
	hpTvnCorr_B = new TH1F("hpTvnCorr_B","hpTvnCorr_B",npt,ptbin.data());
	hpTvnCorr_C = new TH1F("hpTvnCorr_C","hpTvnCorr_C",npt,ptbin.data());
	nharm = 2;
}

EventPlaneAna3sub::~EventPlaneAna3sub() {
	delete tree;
}

void EventPlaneAna3sub::Init() {
	tree = new TChain("tree");
	for(unsigned int ifile = 0; ifile < filelist.size(); ifile++){
		tree -> Add(filelist[ifile].c_str());
	}
	tree -> SetBranchAddress("n", &N);
	tree -> SetBranchAddress("ptg", ptg);
	tree -> SetBranchAddress("etag", etag);
	tree -> SetBranchAddress("phig", phig);
}

void EventPlaneAna3sub::ProcessEvents() {
	for (int ievent = 0; ievent < tree -> GetEntries(); ievent ++) {
		if(ievent % 1000 == 0) std::cout << "processed " << ievent << " events"<< std::endl;
		tree -> GetEntry(ievent);
		float Qx_A = 0;
		float Qy_A = 0; 		
		float Qx_B = 0;
		float Qy_B = 0; 
		float Qx_C = 0;
		float Qy_C = 0; 
		for (int ipart = 0; ipart < N; ipart ++) {
			float phi = phig[ipart];
			float eta = etag[ipart];
			float pt = ptg[ipart];
			if(eta > -0.35 && eta < 0.35 && pt > 0.2){
				Qx_A += cos(nharm * phi);
				Qy_A += sin(nharm * phi);
			}			
			if(eta > -3.9 && eta < -3.1){
				Qx_B += cos(nharm * phi);
				Qy_B += sin(nharm * phi);
			}
			if(eta > -3 && eta < -1){
				Qx_C += cos(nharm * phi);
				Qy_C += sin(nharm * phi);
			}
		} // first particle loop
		float Psi_A = atan2(Qy_A, Qx_A) / nharm;
		float Psi_B = atan2(Qy_B, Qx_B) / nharm;
		float Psi_C = atan2(Qy_C, Qx_C) / nharm;

		hReso_AB -> Fill(cos(nharm * (Psi_A - Psi_B)));
		hReso_AC -> Fill(cos(nharm * (Psi_A - Psi_C)));
		hReso_BC -> Fill(cos(nharm * (Psi_B - Psi_C)));
                for (int ipart = 0; ipart < N; ipart ++) {
			float phi = phig[ipart];
			float eta = etag[ipart];
			float pt = ptg[ipart];
			if(eta > -0.35 && eta < 0.35){
				hpTvnRaw_B -> Fill (pt, cos(nharm * (phi - Psi_B)));
				hpTvnRaw_C -> Fill (pt, cos(nharm * (phi - Psi_C)));
			}
		} // second particle loop
	} // event loop
}

void EventPlaneAna3sub::End() {
	float reso_A = sqrt(hReso_AB -> GetMean() * hReso_AC -> GetMean() / hReso_BC -> GetMean());
	float reso_B = sqrt(hReso_AB -> GetMean() * hReso_BC -> GetMean() / hReso_AC -> GetMean());
	float reso_C = sqrt(hReso_AC -> GetMean() * hReso_BC -> GetMean() / hReso_AB -> GetMean());
	std::cout << "Event plane A " << nharm <<"th resolution = " << reso_A << std::endl;
	std::cout << "Event plane B " << nharm <<"th resolution = " << reso_B << std::endl;
	std::cout << "Event plane C " << nharm <<"th resolution = " << reso_C << std::endl;
	hpTvnCorr_B = (TH1F*)hpTvnRaw_B -> ProjectionX();
	hpTvnCorr_B -> Scale(1./reso_B);
	hpTvnCorr_C = (TH1F*)hpTvnRaw_C -> ProjectionX();
	hpTvnCorr_C -> Scale(1./reso_C);
}


void EventPlaneAna3sub::Plot() {
	gStyle -> SetOptStat(0);
	TCanvas *c1 = new TCanvas();
	TH1F* hempty = new TH1F("","",12,0,6);
	hempty -> SetMaximum(0.20);
	hempty -> SetTitle("1k events event plane method STEG");
	hempty -> GetXaxis() -> SetTitle("p_{T}");
	hempty -> GetYaxis() -> SetTitle(Form("v_{%d}",nharm));
	hempty -> Draw();
	hpTvnCorr_B -> SetMarkerStyle(20);
	hpTvnCorr_B -> SetMarkerSize(0.8);
	hpTvnCorr_B -> SetMarkerColor(1);
	hpTvnCorr_B -> Draw("Psame");
	hpTvnCorr_C -> SetMarkerStyle(24);
	hpTvnCorr_C -> SetMarkerSize(0.8);
	hpTvnCorr_C -> SetMarkerColor(4);
	hpTvnCorr_C -> Draw("Psame");
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
	leg -> AddEntry(hpTvnCorr_B, Form("v_{%d}{BBCS} from EP method",nharm),"p");
	leg -> AddEntry(hpTvnCorr_C, Form("v_{%d}{FVTXS} from EP method",nharm),"p");
	leg -> Draw("same");

	c1 -> Print(Form("V%dvsPT.png",nharm));
}
