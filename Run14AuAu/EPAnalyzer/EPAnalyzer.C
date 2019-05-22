#include "EPAnalyzer.h"
#include <stdlib.h>
#include <time.h>
#include "gsl/gsl_math.h"
#include "getClass.h"
#include "phool.h"
#include "Fun4AllHistoManager.h"
#include "PHGlobal.h"
#include "PHCentralTrack.h"
#include "PHSnglCentralTrack.h"
#include "PHCentralTrackv23.h"
#include "PHSnglCentralTrackv23.h"
#include "TrigLvl1.h"
#include "PreviousEvent.h"

//#include "PhCglList.h"
#include "RpSumXYObject.h"
#include "RpSnglSumXY.h"
#include "ReactionPlaneObject.h"
#include "ReactionPlaneSngl.h"
#include "RpConst.h"
#include "RunHeader.h"
#include "RunNumberRanges.h"

#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h" 
#include "TString.h"
#include <fstream>
#include <string>

#include <PHCompositeNode.h>
#include <recoConsts.h>
#include <TMutNode.h>
#include <Fun4AllReturnCodes.h>
#include <PHMapManager.h>
#include <TRxnpRawScintMap.h>
#include <TRxnpScintMap.h>
#include <TRxnpRawXangMap.h>

#include <MUTOO.h>
#include <FVTXOO.h>
#include <TFvtxCoordMap.h>
#include <TFvtxCoord.h>
#include <PHPoint.h>

#include "SvxClusterList.h"
#include "SvxCluster.h"
#include "SvxSegmentList.h"
#include "SvxSegment.h"
#include "VtxOut.h"

#include "TOAD.h"

#include<TFvtxCompactTrkMap.h>

using namespace std;
using namespace RP;

//_____________________________________________________________________________________________________________________________
EPAnalyzer::EPAnalyzer(const char* output) :
  SubsysReco("EPAnalyzer"), OutputFileName(output), ievent(0), RunNumber(0)
{
  d_outfile=NULL;
  hCentrality = NULL;
  hQx_S = NULL;
  hQy_S = NULL;
  hQx_N = NULL;
  hQy_N = NULL;
  hBBCqS = NULL;
  hBBCqN = NULL;
  hBBCq = NULL;
  hpsi_FVTXS = NULL;
  hpsi_FVTXN = NULL;
  hpsi_FVTXSFVTXN = NULL;
  hReso = NULL;
  hpTvnRaw = NULL;
}

//_____________________________________________________________________________________________________________________________
EPAnalyzer::~EPAnalyzer()
{
  cout << " EPAnalyzer::~EPAnalyzer " << endl;
}

//_____________________________________________________________________________________________________________________________

int EPAnalyzer::Init(PHCompositeNode *topNode)
{
  cout << " EPAnalyzer::Init " << endl;
  
  char name[200];
  float pi = acos(-1.0);

  d_outfile = new TFile(OutputFileName.c_str(),"recreate");
  sprintf(name,"hCentrality");
  hCentrality = new TH1F(name,name,100,0,100);
  sprintf(name,"hQx_N");
  hQx_N = new TH1F(name,name,100,-1,1);
  sprintf(name,"hQy_N");
  hQy_N = new TH1F(name,name,100,-1,1);
  sprintf(name,"hQx_S");
  hQx_S = new TH1F(name,name,100,-1,1);
  sprintf(name,"hQy_S");
  hQy_S = new TH1F(name,name,100,-1,1);
  sprintf(name,"hBBCqS");
  hBBCqS = new TH1F(name,name,100,0,300);
  sprintf(name,"hBBCqN");
  hBBCqN = new TH1F(name,name,100,0,300);
  sprintf(name,"hBBCq");
  hBBCq = new TH1F(name,name,100,0,600);
  sprintf(name,"hpsi_FVTXS");
  hpsi_FVTXS = new TH1F(name,name,100,-pi,pi);
  sprintf(name,"hpsi_FVTXN");
  hpsi_FVTXN = new TH1F(name,name,100,-pi,pi);
  sprintf(name,"hpsi_FVTXSFVTXN");
  hpsi_FVTXSFVTXN = new TH2F(name,name,20,-pi,pi,20,-pi,pi);
  sprintf(name,"hResoVsCent");
  double centbin[] = {0,10,20,30,40,50,60,100};
  hReso = new TProfile(name,name,7,centbin,-1.1,1.1);
  sprintf(name,"hpTvnRaw");
  double ptbin[] = {0,0.2,0.5,1.0,1.6,2.5,3.5,6.0};
  hpTvnRaw = new TProfile(name,name,7,ptbin,-1.1,1.1);

  cout<<"finish of initialize"<<endl;
  return 0;
}

//_____________________________________________________________________________________________________________________________
int EPAnalyzer::InitRun(PHCompositeNode *topNode)
{
  RunHeader* d_runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!d_runheader){
    cout << PHWHERE << " can't find RunHeader " << endl;
    return -1;
  }
  RunNumber = d_runheader->get_RunNumber();
  cout << "---------------------------------------------------------------" << endl;
  cout << " EPAnalyzer::InitRun() :  run number = " << RunNumber << endl;
  cout << "---------------------------------------------------------------" << endl;
  cout << " Loading Recentering Parameters " << endl;
  cout << "---------------------------------------------------------------" << endl;
  cout << "---------------------------------------------------------------" << endl;
  cout << " Loading Flattening Parameters " << endl;
  cout << "---------------------------------------------------------------" << endl;

  cout << "EPAnalyzer::InitRun, run number= " << RunNumber << endl;

  return 0;
}

//_____________________________________________________________________________________________________________________________
int EPAnalyzer::process_event(PHCompositeNode *topNode)
{
  ievent++;

  if(ievent%10000==0) {
    cout<<"EP************* ievent= "<<ievent<<"    *************"<<endl;
  }
  PHMapManager::read(topNode);

  PHGlobal* d_global = findNode::getClass<PHGlobal>(topNode, "PHGlobal");
  if(!d_global){
    cout << PHWHERE << "Could not find PHGlobal !" << endl;
    return -1;
  }

  PHCentralTrack* d_cnt = findNode::getClass<PHCentralTrack>(topNode, "PHCentralTrack");
  if(!d_cnt)
  {
    cout << PHWHERE << "Could not find PHCentralTrack !" << endl;
    d_cnt = findNode::getClass<PHCentralTrack>(topNode, "PhCglList");
    if(!d_cnt)
    {
      cout << PHWHERE << "Could not find PhCglList !" << endl;
      return 1;
    }
  }

  TrigLvl1 *triglvl1 = findNode::getClass<TrigLvl1>(topNode, "TrigLvl1");
  if(!triglvl1)
    {
      cout << PHWHERE << "Could not find TrigLvl1 !" << endl;
      return 1;
    }

  VtxOut *vtxout = findNode::getClass<VtxOut>(topNode,"VtxOut");
  if(!vtxout) {cout<<"No Vtxout information!"<<endl;return 1;}

  RpSumXYObject* d_rp = findNode::getClass<RpSumXYObject>(topNode, "RpSumXYObject");
  
  if(!d_rp){
    cout<< PHWHERE << "Could not find the RP Object!"<< endl;
    return -1;
  }
  
  ReactionPlaneObject* rp = findNode::getClass<ReactionPlaneObject>(topNode,"ReactionPlaneObject");
  if( !rp ){
    cout<<"can't find ReactionPlaneObject "<<endl;
    return DISCARDEVENT;//exit(1);
  }

  RunHeader* d_runheader = findNode::getClass<RunHeader>(topNode, "RunHeader");
  if(!d_runheader){
    cout << PHWHERE << " can't find RunHeader " << endl;
    return -1;
  }

  //global

//  int run = d_runheader->get_RunNumber();
  int trig = triglvl1->get_lvl1_trigscaled();
  if(!(trig & 0x00000010)>0) return -1; //MinBias trigger
  float bbcv = d_global->getBbcZVertex();
  int cent = (int)d_global->getCentrality();
  float bbcqN  = d_global->getBbcChargeN();   
  float bbcqS  = d_global->getBbcChargeS();   
  if(fabs(bbcv)>10) return -1;
  
  hCentrality -> Fill(cent);
  hBBCqS -> Fill(bbcqS);
  hBBCqN -> Fill(bbcqN);
  hBBCq -> Fill(bbcqN + bbcqS);
  int nharm = 2;
  
  // Get raw Q-vectors without any recalibration
    RpSnglSumXY *s_rp;
    s_rp = d_rp->getRpSumXY(RP::calcIdCode(91,40,1));
    float Qx_S = (s_rp && s_rp->Weight() > 0 && s_rp->Weight() < 999) ? s_rp->QVector(0) / s_rp->Weight() : -9999;
    s_rp = d_rp->getRpSumXY(RP::calcIdCode(91,40,1));
    float Qy_S = (s_rp && s_rp->Weight() > 0 && s_rp->Weight() < 999) ? s_rp->QVector(1) / s_rp->Weight() : -9999;
    
    s_rp = d_rp->getRpSumXY(RP::calcIdCode(91,41,1));
    float Qx_N = (s_rp && s_rp->Weight() > 0 && s_rp->Weight() < 999) ? s_rp->QVector(0) / s_rp->Weight(): -9999;
    s_rp = d_rp->getRpSumXY(RP::calcIdCode(91,41,1));
    float Qy_N = (s_rp && s_rp->Weight() > 0 && s_rp->Weight() < 999) ? s_rp->QVector(1) / s_rp->Weight(): -9999;
    
    hQx_S -> Fill(Qx_S);
    hQy_S -> Fill(Qy_S);
    hQx_N -> Fill(Qx_N);
    hQy_N -> Fill(Qy_N);
 
    //FVTX, Particle for Event Plane angle determination
    ReactionPlaneSngl* rpsngl;
    rpsngl = rp->getReactionPlane(RP::calcIdCode(91,40,1));
    float Psi_S = (rpsngl) ? rpsngl->GetPsi() : -9999;
    rpsngl = rp->getReactionPlane(RP::calcIdCode(91,41,1));
    float Psi_N = (rpsngl) ? rpsngl->GetPsi() : -9999;

    if(fabs(Psi_S) < 4 && fabs(Psi_N) < 4) {
        hpsi_FVTXS -> Fill(Psi_S);
        hpsi_FVTXN -> Fill(Psi_N);
        hpsi_FVTXSFVTXN -> Fill(Psi_S, Psi_N);
        hReso -> Fill(cent, cos(nharm*(Psi_S - Psi_N)));//South - North
    }

// Tracks, Particle of Interest
  for(unsigned int itrk=0;itrk<d_cnt->get_npart();itrk++){
    PHSnglCentralTrack *d_scnt = d_cnt->get_track(itrk);
    if (d_scnt->get_mom()<10.0 && d_scnt->get_mom()>0.02 &&
        fabs(d_scnt->get_pc3dphi())<5.0 && fabs(d_scnt->get_pc3dz())<10.0
        && (d_scnt->get_quality() ==31 || d_scnt->get_quality() ==63)
        && fabs(d_scnt->get_zed())<75 && d_scnt->get_n0()<0
        ) {
      float phi       = d_scnt->get_phi0();
      float pt        = d_scnt->get_mom()*sin(d_scnt->get_the0());
      float eta       = -log(tan(0.5*d_scnt->get_the0()));

      double sdphi = d_scnt->get_pc3sdphi();
      double sdz = d_scnt->get_pc3sdz();
      if(fabs(sdphi)<3.0 && fabs(sdz)<3.0){
          if(cent > 20 && cent < 30) {
              if (eta > 0 && fabs(Psi_S) < 4)
            hpTvnRaw -> Fill(pt, cos(nharm * (phi - Psi_S)));
              if (eta < 0 && fabs(Psi_N) < 4)
            hpTvnRaw -> Fill(pt, cos(nharm * (phi - Psi_N)));
          }
      }
 }
 }
  
  return 0;
}   

//_____________________________________________________________________________________________________________________________
int EPAnalyzer::End(PHCompositeNode* topNode)
{
  cout << "End of EPAnalyzer for Run " << RunNumber << endl;
  cout << "Total # of events = " << ievent << endl;
  cout << "OutputFileName = " << OutputFileName << endl;

  //HistoManager->dumpHistos(OutputFileName);
  if(d_outfile) {
    d_outfile->cd();
    hQx_S -> Write();
    hQy_S -> Write();
    hQy_N -> Write();
    hQx_N -> Write();
    hpsi_FVTXS ->Write();
    hpsi_FVTXN ->Write();
    hpsi_FVTXSFVTXN ->Write();
    hBBCqN -> Write();
    hBBCqS -> Write();
    hBBCq -> Write();
    hReso ->Write();
    hpTvnRaw ->Write();
    d_outfile->Close();
  }else{

    cout<<PHWHERE<<"ERROR: No output file set!"<<endl;

  }

  return 0;
}
