#ifndef __EPANALYZER_H__
#define __EPANALYZER_H__

#include <string>
#include "SubsysReco.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 

class EPAnalyzer: public SubsysReco
{
 public:

  EPAnalyzer(const char* output="rpanase.root");
  virtual ~EPAnalyzer();

  int Init(PHCompositeNode *topNode);
  int InitRun(PHCompositeNode *topNode);
  int process_event(PHCompositeNode *topNode);
  int End(PHCompositeNode *topNode);
  
 private:
  const std::string OutputFileName;
  int ievent;
  int RunNumber;

  TFile *d_outfile;
 
  TH1F* hCentrality;
  TH1F* hBBCqS;
  TH1F* hBBCqN;
  TH1F* hBBCq;

  TH1F* hQx_N;
  TH1F* hQy_N;
  TH1F* hQx_S;
  TH1F* hQy_S;
  TH1F* hpsi_FVTXS;
  TH1F* hpsi_FVTXN;
  TH2F* hpsi_FVTXSFVTXN;
  TProfile* hReso;
  TProfile* hpTvnRaw;
};

#endif /* __EPANALYZER_H__ */
