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
  
  TH1F* hpsi;
  TProfile* hReso;
  TProfile* hpTvnRaw;
};

#endif /* __EPANALYZER_H__ */
