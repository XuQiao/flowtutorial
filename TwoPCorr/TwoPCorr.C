#include "TwoPCorr.h"

TwoPCorr::TwoPCorr():
	N(0),
	ptg(),
	etag(),
	phig()
	{
	TH1::SetDefaultSumw2();
	pi = acos(-1);
	nDepth = 3;
	float ptarr[] = {0, 0.2, 0.5, 0.8, 1.5, 3.0, 6.0};
        ptbin = std::vector<float>(ptarr, ptarr + sizeof(ptarr)/sizeof(ptarr[0]));
        npt = ptbin.size()-1;
	for (int ipt = 0; ipt < npt; ipt ++){
		hSignal[ipt] = new TH2F(Form("hSignal_%d",ipt),"",40,-8,8,20,-0.5*pi,1.5*pi);
		hBackgd[ipt] = new TH2F(Form("hBackgd_%d",ipt),"",40,-8,8,20,-0.5*pi,1.5*pi);
		fcos[ipt] = new TF1(Form("fcos_%d",ipt),"[0]*(1+2*[1]*cos(1*x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x))",-0.5*pi,1.5*pi);
	}
	hv2 = new TH1F("hv2","",npt,ptbin.data());
	hv3 = new TH1F("hv3","",npt,ptbin.data());
}

TwoPCorr::~TwoPCorr() {
	delete tree;
}

void TwoPCorr::Init() {
	tree = new TChain("tree");
	for(unsigned int ifile = 0; ifile < filelist.size(); ifile++){
		tree -> Add(filelist[ifile].c_str());
	}	
	tree -> SetBranchAddress("n", &N);
	tree -> SetBranchAddress("ptg", ptg);
	tree -> SetBranchAddress("etag", etag);
	tree -> SetBranchAddress("phig", phig);
}

void TwoPCorr::ProcessEvents() {
	for (int ievent = 0; ievent < tree -> GetEntries(); ievent ++) {
		if(ievent % 1000 == 0) std::cout << "processed " << ievent << " events"<< std::endl;
		tree -> GetEntry(ievent);
		vector<float> phivec;
		vector<float> etavec;
		vector<float> ptvec;
		for (int ipart = 0; ipart < N; ipart ++) {
			float i_phi = phig[ipart];
			float i_eta = etag[ipart];
			float i_pt = ptg[ipart];
			int ipt = -1;
			for (ipt = 0; ipt < npt; ipt ++) {
				if (i_pt >= ptbin[ipt] && i_pt < ptbin[ipt+1])
					break;
			}
			phivec.push_back(i_phi);
			etavec.push_back(i_eta);
			ptvec.push_back(i_pt);
			for (int jpart = 0; jpart < N; jpart ++) {
				float j_phi = phig[jpart];
				float j_eta = etag[jpart];
				float j_pt = ptg[jpart];
				float deta = j_eta - i_eta;
				float dphi = j_phi - i_phi;
				if (fabs(deta) < 0.01) continue; // remove auto-correlation
				if(dphi < -0.5 * pi) dphi += 2 * pi;
				if(dphi > 1.5 * pi) dphi -= 2 * pi;
				if (j_pt >= ptbin[ipt] && j_pt < ptbin[ipt+1])
                                    hSignal[ipt] -> Fill(deta, dphi);
			} // first particle loop
			for(unsigned int jevent = 0; jevent < phivec_mix.size(); jevent ++) {
				for(unsigned int jpart = 0; jpart < phivec_mix[jevent].size(); jpart ++) {
					float jmix_phi = phivec_mix[jevent][jpart];
					float jmix_eta = etavec_mix[jevent][jpart];
					float jmix_pt = ptvec_mix[jevent][jpart];
					float deta = jmix_eta - i_eta;
					float dphi = jmix_phi - i_phi;
					if(dphi < -0.5 * pi) dphi += 2 * pi;
					if(dphi > 1.5 * pi) dphi -= 2 * pi;
				        if (jmix_pt >= ptbin[ipt] && jmix_pt < ptbin[ipt+1])
                                            hBackgd[ipt] -> Fill(deta, dphi);
				}
			} // second particle loop
		}
		phivec_mix.push_back(phivec);
		etavec_mix.push_back(etavec);
		ptvec_mix.push_back(ptvec);
		if (phivec_mix.size() > nDepth) {
			phivec_mix.erase (phivec_mix.begin());
			etavec_mix.erase (etavec_mix.begin());
			ptvec_mix.erase (ptvec_mix.begin());
		} 
	} // event loop
}

void TwoPCorr::End() {
	for (int ipt = 0; ipt < npt; ipt ++) {
		hCorr[ipt] = (TH2F*)hSignal[ipt] -> Clone(Form("hCorr_pt%d",ipt));
		hCorr[ipt] -> Scale(1./hCorr[ipt]->GetBinContent(hCorr[ipt]->FindBin(0,0)));
		if (hBackgd[ipt] -> GetEntries() > 0)
			hCorr[ipt] -> Divide(hBackgd[ipt]);
		hCorr[ipt] -> Scale(hBackgd[ipt]->GetBinContent(hBackgd[ipt]->FindBin(0,0)));
		//Integrated 
		hIntCorr[ipt] = (TH1F*)hSignal[ipt] -> ProjectionY();
		TH1F* hIntBackgd = (TH1F*)hBackgd[ipt] -> ProjectionY();
		hIntCorr[ipt] -> Scale(1./hIntCorr[ipt]-> GetBinContent(hIntCorr[ipt]->FindBin(0)));
		if (hIntBackgd -> GetEntries() > 0)
			hIntCorr[ipt] -> Divide(hIntBackgd);
		hIntCorr[ipt] -> Scale(hIntBackgd->GetBinContent(hIntBackgd->FindBin(0)));
		fcos[ipt] -> SetParameters(1,0,0,0,0);
		hIntCorr[ipt] -> Fit(fcos[ipt],"NQ0");
		if (fcos[ipt] -> GetParameter(2) > 0){
                    float v2 = sqrt(fcos[ipt] -> GetParameter(2));
                    float v2_err = 	1./2*fcos[ipt] -> GetParError(2)/v2;
                    hv2 -> SetBinContent(ipt+1, v2);
                    hv2 -> SetBinError(ipt+1, v2_err);
                }
                if(fcos[ipt] -> GetParameter(3) > 0){
                    float v3 = sqrt(fcos[ipt] -> GetParameter(3));
                    float v3_err = 1./2*fcos[ipt] -> GetParError(3)/v3;
                    hv3 -> SetBinContent(ipt+1, v3);
                    hv3 -> SetBinError(ipt+1, v3_err);
		}
	}
}


void TwoPCorr::Plot() {
	gStyle -> SetOptStat(0);
	TCanvas *c1 = new TCanvas();
	TH1F* hempty = new TH1F("","",12,0,6);
	hempty -> SetMaximum(0.20);
	hempty -> SetTitle("1k events TwoPCorr method STEG");
	hempty -> GetXaxis() -> SetTitle("p_{T}");
	hempty -> GetYaxis() -> SetTitle(Form("v_{2}"));
	hempty -> Draw();
        TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        TF1 *V3vsPt = new TF1("V3vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        TF1 *V4vsPt = new TF1("V4vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        V2vsPt -> SetParameters(3.31699,2.35142,3.49188,3.54429,.00005,1.50600);
        V3vsPt -> SetParameters(3.2,2.3,3.4,2.1,.00005,1.4);
        V4vsPt -> SetParameters(4.8,2.1,3.4,2.1,.00005,1.4);
	hv2 -> SetMarkerStyle(20);
	hv2 -> SetMarkerSize(0.8);
	hv2 -> SetMarkerColor(1);
	hv2 -> Draw("Psame");
	V2vsPt -> SetLineStyle(1);
	V2vsPt -> Draw("same");
	TLegend *leg = new TLegend(0.1,0.7,0.45,0.88);
	leg -> SetTextSize(0.05);
	leg -> SetBorderSize(0);
	leg -> SetFillStyle(0);
	leg -> AddEntry(V2vsPt, Form("input v_{2}"),"l");
	leg -> AddEntry(hv2, Form("v_{2} from TwoPCorr method"),"p");
	leg -> Draw("same");
	c1 -> Print(Form("V2vsPT.png"));

	c1 = new TCanvas();
	hempty -> GetYaxis() -> SetTitle(Form("v_{3}"));
	hempty -> Draw();
	hv3 -> SetMarkerStyle(20);
	hv3 -> SetMarkerSize(0.8);
	hv3 -> SetMarkerColor(1);
	hv3 -> Draw("Psame");
	V3vsPt -> SetLineStyle(1);
	V3vsPt -> Draw("same");
	leg = new TLegend(0.1,0.75,0.45,0.88);
	leg -> SetTextSize(0.05);
	leg -> SetBorderSize(0);
	leg -> SetFillStyle(0);
	leg -> AddEntry(V3vsPt, Form("input v_{3}"),"l");
	leg -> AddEntry(hv3, Form("v_{3} from TwoPCorr method"),"p");
	leg -> Draw("same");
	c1 -> Print(Form("V3vsPT.png"));
	
        int ipt = 3; //only draw 0.8 < pT < 1.5
	c1 = new TCanvas();
        c1->SetTheta(60);
	c1->SetPhi(30);
	TLatex latex;
	latex.SetNDC();
	hCorr[ipt]->SetTitle("1k Two Particle Corr method STEG");
	hCorr[ipt]->GetXaxis()->SetTitle("#Delta#eta");
	hCorr[ipt]->GetYaxis()->SetTitle("#Delta#phi");
	hCorr[ipt]->GetZaxis()->SetTitle("C(#Delta#eta,#Delta#phi)");
	hCorr[ipt] -> Draw("surf1");
	latex.SetNDC();
	latex.DrawLatex(0.12,0.7,Form("%.1f<p_{T}^{trig}<%.1f",ptbin[ipt],ptbin[ipt+1]));
	latex.DrawLatex(0.12,0.77,Form("%.1f<p_{T}^{asso}<%.1f",ptbin[ipt],ptbin[ipt+1]));
	c1 -> Print(Form("hCorr_pt%d.png",ipt));

	c1 = new TCanvas();
	hempty = new TH1F("","",10,-0.5*pi,1.5*pi);
	hempty -> SetMaximum(1.02);
	hempty -> SetMinimum(0.96);
	hempty -> SetTitle("1k events TwoPCorr method STEG");
	hempty -> GetXaxis() -> SetTitle("#Delta#phi");
	hempty -> GetYaxis() -> SetTitle(Form("C(#Delta#phi)"));
	hempty -> Draw();
	TF1* fcos1 = new TF1("fcos1","[0]*(1+2*[1]*cos(1*x))",-0.5*pi,1.5*pi);
	fcos1->SetParameters(fcos[ipt]->GetParameter(0),fcos[ipt]->GetParameter(1));
	TF1* fcos2 = new TF1("fcos2","[0]*(1+2*[1]*cos(2*x))",-0.5*pi,1.5*pi);
	fcos2->SetParameters(fcos[ipt]->GetParameter(0),fcos[ipt]->GetParameter(2));
	TF1* fcos3 = new TF1("fcos3","[0]*(1+2*[1]*cos(3*x))",-0.5*pi,1.5*pi);
	fcos3->SetParameters(fcos[ipt]->GetParameter(0),fcos[ipt]->GetParameter(3));
	TF1* fcos4 = new TF1("fcos4","[0]*(1+2*[1]*cos(4*x))",-0.5*pi,1.5*pi);
	fcos4->SetParameters(fcos[ipt]->GetParameter(0),fcos[ipt]->GetParameter(4));
	hIntCorr[ipt] -> SetMarkerStyle(20);
	hIntCorr[ipt] -> SetMarkerColor(1);
	hIntCorr[ipt] -> SetMarkerSize(0.8);
	hIntCorr[ipt] -> Draw("Psame");
	fcos1->SetLineColor(2);fcos1->SetLineStyle(2);fcos1->SetLineWidth(2);
	fcos2->SetLineColor(3);fcos2->SetLineStyle(2);fcos2->SetLineWidth(2);
	fcos3->SetLineColor(4);fcos3->SetLineStyle(2);fcos3->SetLineWidth(2);
	fcos4->SetLineColor(5);fcos4->SetLineStyle(2);fcos4->SetLineWidth(2);
	latex.SetNDC();
	latex.DrawLatex(0.12,0.7,Form("%.1f<p_{T}^{trig}<%.1f",ptbin[ipt],ptbin[ipt+1]));
	latex.DrawLatex(0.12,0.77,Form("%.1f<p_{T}^{asso}<%.1f",ptbin[ipt],ptbin[ipt+1]));
	fcos[ipt]->SetLineColor(1);
	fcos[ipt]->Draw("same");
	fcos1->Draw("same");
	fcos2->Draw("same");
	fcos3->Draw("same");
	fcos4->Draw("same");
	leg = new TLegend(0.6,0.6,0.85,0.85);
	leg->SetTextSize(0.04);
	leg->SetBorderSize(0);
	leg->SetFillColor(0);
	leg->AddEntry(fcos[ipt],"Fourier Fit","l");
	leg->AddEntry(fcos1,"A*(1+c_{1}cos(x))","l");
	leg->AddEntry(fcos2,"A*(1+c_{2}cos(2x))","l");
	leg->AddEntry(fcos3,"A*(1+c_{3}cos(3x))","l");
	leg->AddEntry(fcos4,"A*(1+c_{4}cos(4x))","l");
	leg->Draw("same");

	c1 -> Print(Form("hIntCorr_pt%d.png",ipt));
}
