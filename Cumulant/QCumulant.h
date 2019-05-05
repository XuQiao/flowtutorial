#include <vector>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"
#include "TStyle.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TComplex.h"
#include "TLegend.h"

const int maxN = 1000;
class QCumulant {
	TChain *tree;
	int N;
	float phig[maxN];
	float ptg[maxN];
	float etag[maxN];

	int nharm;
        std::vector<float> ptbin;
        int npt;
	float M_;
	float Mpt_[10];

        std::vector<std::string> filelist;

	TH1F* hIntcumu2;
	TH1F* hIntcumu4;
	TH1F* hcumu2[10];
	TH1F* hcumu4[10];

	TH1F* hv2;
	TH1F* hv4;

public:
	QCumulant();
	~QCumulant();
	void set_nharm(int _nharm) { nharm = _nharm;}
	void push_file(std::string file) {filelist.push_back(file);}
	void Init();
	void ProcessEvents();
	void End();
	void Plot();
};
