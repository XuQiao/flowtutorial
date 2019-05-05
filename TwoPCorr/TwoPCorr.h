#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLegend.h"

const int maxN = 1000;
class TwoPCorr {
	TChain *tree;
	int N;
	float phig[maxN];
	float ptg[maxN];
	float etag[maxN];

        std::vector<std::vector<float> >phivec_mix;
        std::vector<std::vector<float> >ptvec_mix;
        std::vector<std::vector<float> >etavec_mix;
        
        std::vector<float> ptbin;
        int npt;

	TH2F* hSignal[10];
	TH2F* hBackgd[10];
	TH2F* hCorr[10];
	TH1F* hIntCorr[10];

        std::vector<string> filelist;

	TF1* fcos[10];
	unsigned int nDepth;
	float pi;

	TH1F* hv2;
	TH1F* hv3;

public:
	TwoPCorr();
	~TwoPCorr();
	void push_file(string file) {filelist.push_back(file);}
	void Init();
	void ProcessEvents();
	void End();
	void Plot();
};
