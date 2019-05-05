#include "QCumulant.h"

QCumulant::QCumulant():
	N(0),
	ptg(),
	etag(),
	phig(),
	M_(0),
	Mpt_()
	{
	TH1::SetDefaultSumw2();
	float ptarr[] = {0, 0.2, 0.5, 0.8, 1.5, 3.0, 6.0};
        ptbin = std::vector<float>(ptarr, ptarr + sizeof(ptarr)/sizeof(ptarr[0]));
        npt = ptbin.size()-1;
	hIntcumu2 = new TH1F("hcumu2","",100, -1.1, 1.1);
	hIntcumu4 = new TH1F("hcumu4","",100, -1.1, 1.1);
	for(int ipt = 0; ipt < npt; ipt ++) {
		hcumu2[ipt] = new TH1F(Form("hcumu2_%d",ipt),"",100, -1.1, 1.1);
		hcumu4[ipt] = new TH1F(Form("hcumu4_%d",ipt),"",100, -1.1, 1.1);
	}
	nharm = 2;
	hv2 = new TH1F(Form("hv%d2",nharm),"",npt,ptbin.data());
	hv4 = new TH1F(Form("hv%d4",nharm),"",npt,ptbin.data());
}

QCumulant::~QCumulant() {
	delete tree;
}

void QCumulant::Init() {
	tree = new TChain("tree");
	for(unsigned int ifile = 0; ifile < filelist.size(); ifile++){
		tree -> Add(filelist[ifile].c_str());
	}
	tree -> SetBranchAddress("n", &N);
	tree -> SetBranchAddress("ptg", ptg);
	tree -> SetBranchAddress("etag", etag);
	tree -> SetBranchAddress("phig", phig);
}

void QCumulant::ProcessEvents() {
	for (int ievent = 0; ievent < tree -> GetEntries(); ievent ++) {
		if(ievent % 1000 == 0) std::cout << "processed " << ievent << " events"<< std::endl;
		tree -> GetEntry(ievent);
		float Qx = 0;
		float Qy = 0;		
		float Q2x = 0;
		float Q2y = 0;
		float Qxpt[10]={};
		float Qypt[10]={};
		float Q2xpt[10]={};
		float Q2ypt[10]={};
		float M = 0;
                float Mpt[10] = {};
		for (int ipart = 0; ipart < N; ipart ++) {
			float phi = phig[ipart];
			float eta = etag[ipart];
			float pt = ptg[ipart];
			float weight = 1.;
			Qx += weight * cos(nharm * phi);
			Qy += weight * sin(nharm * phi);			
			Q2x += weight * cos(2.*nharm * phi);
			Q2y += weight * sin(2.*nharm * phi);
			M += weight;
			M_ += weight;
			int ipt = -1;
			for (ipt = 0; ipt < npt; ipt ++) {
				if (pt >= ptbin[ipt] && pt < ptbin[ipt+1])
					break;
			}
			Qxpt[ipt] += weight * cos(nharm * phi);
			Qypt[ipt] += weight * sin(nharm * phi);			
			Q2xpt[ipt] += weight * cos(2.*nharm * phi);
			Q2ypt[ipt] += weight * sin(2.*nharm * phi);
			Mpt[ipt] += weight;
			Mpt_[ipt] += weight;
		}
		TComplex Qn(Qx, Qy);
		TComplex Q2n(Q2x, Q2y);
		hIntcumu2 -> Fill((Qn.Rho()*Qn.Rho()-M)/(M*(M-1)));
		hIntcumu4 -> Fill((Qn.Rho()*Qn.Rho()*Qn.Rho()*Qn.Rho()
			+Q2n.Rho()*Q2n.Rho()-2.*(Q2n*TComplex::Conjugate(Qn)*TComplex::Conjugate(Qn)).Re()
			-2.*(2.*(M-2)*Qn.Rho()*Qn.Rho()-M*(M-3)))/(M*(M-1)*(M-2)*(M-3)));
		
		for (int ipt = 0; ipt < npt; ipt ++) {
			TComplex pn(Qxpt[ipt],Qypt[ipt]);
			TComplex qn(Qxpt[ipt],Qypt[ipt]);
			TComplex p2n(Q2xpt[ipt],Q2ypt[ipt]);
			TComplex q2n(Q2xpt[ipt],Q2ypt[ipt]);
			float Mq = Mpt[ipt];
			hcumu2[ipt] -> Fill(((pn*TComplex::Conjugate(Qn)).Re()-Mq)/(Mpt[ipt]*M-Mq));
			hcumu4[ipt] -> Fill((pn*Qn*TComplex::Conjugate(Qn)*TComplex::Conjugate(Qn)
				-q2n*TComplex::Conjugate(Qn)*TComplex::Conjugate(Qn)-pn*Qn*TComplex::Conjugate(Q2n)
				-2.*M*pn*TComplex::Conjugate(Qn)-2.*Mq*Qn.Rho()*Qn.Rho()+7.*qn*TComplex::Conjugate(Qn)
				-Qn*TComplex::Conjugate(qn)+q2n*TComplex::Conjugate(Q2n)+2.*pn*TComplex::Conjugate(Qn)
				+2.*Mq*M - 6.*Mq).Re()/((Mpt[ipt]*M-3*Mq)*(M-1)*(M-2)));
		}
	} // event loop
}

void QCumulant::End() {
	float v2Int = sqrt(hIntcumu2->GetMean());
	float v2Int_err = 1./2*hIntcumu2->GetStdDev()/v2Int/sqrt(M_);	
	float v4Int = sqrt(sqrt(-(hIntcumu4->GetMean()-2*pow(hIntcumu2->GetMean(),2))));
	float v4Int_err = 1./4*hIntcumu4->GetStdDev()*pow(v4Int,-3)/sqrt(M_);
	std::cout << "integrated v"<<nharm<<"{2} = " << v2Int << std::endl;
	std::cout << "integrated v"<<nharm<<"{4} = " << v4Int << std::endl;
	for (int ipt = 0; ipt < npt; ipt ++) {
                if(Mpt_[ipt] <= 0) continue;
		hv2 -> SetBinContent(ipt+1, hcumu2[ipt]->GetMean()/v2Int);
		hv2 -> SetBinError(ipt+1, hcumu2[ipt]->GetStdDev()/sqrt(Mpt_[ipt]/v2Int));
		hv4 -> SetBinContent(ipt+1, -(hcumu4[ipt]->GetMean()-2*hcumu2[ipt]->GetMean()*hIntcumu2->GetMean())/pow(v4Int,3));
		hv4 -> SetBinError(ipt+1, hcumu4[ipt]->GetStdDev()/sqrt(Mpt_[ipt])/pow(v4Int,3));
	}
}


void QCumulant::Plot() {
	gStyle -> SetOptStat(0);
	TCanvas *c1 = new TCanvas();
	TH1F* hempty = new TH1F("","",12,0,6);
	hempty -> SetMaximum(0.20);
	hempty -> SetTitle("1M events Multi-PCorr method STEG with non-flow");
	hempty -> GetXaxis() -> SetTitle("p_{T}");
	hempty -> GetYaxis() -> SetTitle(Form("v_{%d}",nharm));
	hempty -> Draw();
	hv2 -> SetMarkerStyle(20);
	hv4 -> SetMarkerStyle(24);
	hv2 -> Draw("Psame");
	hv4 -> Draw("Psame");

        TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        TF1 *V3vsPt = new TF1("V3vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        TF1 *V4vsPt = new TF1("V4vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*([4]+(1/x)^[5])",0.2,6.0);
        V2vsPt -> SetParameters(3.31699,2.35142,3.49188,3.54429,.00005,1.50600);
        V3vsPt -> SetParameters(3.2,2.3,3.4,2.1,.00005,1.4);
        V4vsPt -> SetParameters(4.8,2.1,3.4,2.1,.00005,1.4);

	TLegend *leg = new TLegend(0.1,0.70,0.45,0.88);
	leg -> SetTextSize(0.05);
	leg -> SetBorderSize(0);
	leg -> SetFillStyle(0);
	if (nharm == 2) {
		leg -> AddEntry(V2vsPt, Form("input v_{%d}",nharm),"l");
		V2vsPt -> Draw("same");
	}	
	if (nharm == 3) {
		leg -> AddEntry(V3vsPt, Form("input v_{%d}",nharm),"l");
		V3vsPt -> Draw("same");
	}
	leg -> AddEntry(hv2, Form("v_{%d}{2} from Multi-PCorr method",nharm),"p");
	leg -> AddEntry(hv4, Form("v_{%d}{4} from Multi-PCorr method",nharm),"p");
	leg -> Draw("same");
	c1 -> Print(Form("V%dvsPT.png",nharm));
}
