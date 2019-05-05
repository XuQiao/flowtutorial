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
#include "TLegend.h"

const int maxN = 1000;
class EventPlaneAna {
	TChain *tree;
	int N;
	float phig[maxN];
	float ptg[maxN];
	float etag[maxN];

	int nharm;

        std::vector<std::string> filelist;

	TH1F* hReso;
	TProfile* hpTvnRaw;
	TH1F* hpTvnCorr;

        std::vector<double> ptbin;
	int npt;

public:
	EventPlaneAna();
	~EventPlaneAna();
	void set_nharm(int _nharm) { nharm = _nharm;}
	void push_file(string file) {filelist.push_back(file);}
	void Init();
	void ProcessEvents();
	void End();
	void Plot();
};
