void plot() {
        TFile *fin = TFile::Open("test_vnres_run14auau.root");
        TProfile *hResoVsCent = (TProfile*) fin -> Get("hResoVsCent");
        TH1F* hReso = (TH1F*)hResoVsCent -> ProjectionX();
        gStyle -> SetOptStat(0);
        int nharm = 2;
	TCanvas *c1 = new TCanvas();
	TH1F* hempty = new TH1F("","",100,0,100);
	hempty -> SetMaximum(0.80);
	hempty -> SetTitle("Run14 Au+Au event plane resolution");
	hempty -> GetXaxis() -> SetTitle("Centrality(%)");
	hempty -> GetYaxis() -> SetTitle(Form("Resolution",nharm));
        hempty -> Draw();
        for (int ibin = 1; ibin <= hReso-> GetNbinsX(); ibin ++) {
            hReso -> SetBinContent(ibin, sqrt(hReso->GetBinContent(ibin)));
            hReso -> SetBinError(ibin, hReso->GetBinError(ibin)/2./hReso->GetBinContent(ibin));
        }
        hReso -> Draw("same");
        c1 -> Print("hReso.png");
	float reso20_30 = hReso -> GetBinContent(3);
        hpTvnCorr = (TH1F*) hpTvnRaw -> ProjectionX();
        hpTvnCorr -> Scale(1./reso20_30); //Scale by resolution
	TCanvas *c1 = new TCanvas();
	hempty = new TH1F("","",12,0,4);
	hempty -> SetMaximum(0.30);
	hempty -> SetTitle("Run14 Au+Au event plane method");
	hempty -> GetXaxis() -> SetTitle("p_{T}");
	hempty -> GetYaxis() -> SetTitle(Form("v_{%d}",nharm));
	hempty -> Draw();
	hpTvnCorr -> SetMarkerStyle(20);
	hpTvnCorr -> SetMarkerSize(0.8);
	hpTvnCorr -> SetMarkerColor(1);
	hpTvnCorr -> Draw("Psame");
        TLatex latex;
        latex.SetNDC();
	TLegend *leg = new TLegend(0.1,0.75,0.45,0.88);
	leg -> SetTextSize(0.05);
	leg -> SetBorderSize(0);
	leg -> SetFillStyle(0);
	leg -> AddEntry(hpTvnCorr, Form("v_{%d} from EP method",nharm),"p");
	leg -> Draw("same");
        latex.DrawLatex(0.6,0.85,"20 < centrality < 30%");
	c1 -> Print(Form("V%dvsPT.png",nharm));
}
