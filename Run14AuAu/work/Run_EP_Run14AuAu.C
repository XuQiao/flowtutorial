#include <string>
#include <iostream>
#include <fstream>
using namespace std;

void Run_EP_Run14AuAu(const char *outFile = "output.root")
{
  gSystem->Load("/phenix/plhf/xuq/phenix/flow/Run14AuAu/install/lib/libRPCalibRun.so");
  gSystem->Load("/phenix/plhf/xuq/phenix/flow/flowtutorial/Run14AuAu/install/lib/libEPAnalyzer.so");
  // For event plane angle calibration, private version by Hachiya for Run 14 Au+Au
  //gSystem->Load("librpc_subsysreco");
  //gSystem->Load("librpc_muotrackreco");
  gSystem->Load("librecal");
  //gSystem->Load("libpicodst_object.so");
  
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  //Reaction Plane Node
  recoConsts *rc = recoConsts::instance();  
  rc->set_IntFlag("RPCALIB_READFROMDB", 0);
  rc->set_CharFlag("RPCALIB_CALIBFILENAME", "");

  rc->set_IntFlag("RP_SKIP_RECENTERING", 0); // do recentering
  rc->set_IntFlag("RP_SKIP_FLATTENING", 0); // do flattening

  RPReadCalibTree *readT = new RPReadCalibTree();
  readT->setTreeFileRecent("RP_recent_run14pro106_newcent_merge.root");
  readT->setTreeFileFlat("RP_flat_run14pro106_newcent_merge.root");
  readT->setTOADname("hachiya/15.08");
  se->registerSubsystem(readT);
 
  SubsysReco *rpana = new EPAnalyzer(outFile);
  se->registerSubsystem(rpana);
}

void 
InputData(vector<string> &indata)
{
  indata.push_back("CNT");
  indata.push_back("DST_EVE");
  return;
}
