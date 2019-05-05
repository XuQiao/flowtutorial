#include "TMath.h"
#include "Math/DistSampler.h"
#include "Math/DistSamplerOptions.h"
#include "Math/Factory.h"
#include "TTree.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TRandom3.h"
#include "SetSeedOwn.C"

using namespace ROOT::Math;
void proSTEGvnnf(int ifile, int N, int M, int RMSM){
    TRandom3 *r = new TRandom3(0);
    /*
    float v1 = 0.;
    float v2 = 0.05;
    float v3 = 0.;
    float v4 = 0.;
    float v5 = 0.;
    float v6 = 0.;
    */
    float pi = acos(-1);
    const int maxnh = 1000;
    int b_n;
    float b_psig;
    float b_ptg[maxnh], b_etag[maxnh], b_phig[maxnh];
    TFile *fout = new TFile(Form("output/vndatanf_N%d_M%dRMS%d_ifile%d.root",N,M,RMSM,ifile),"recreate");
    //TTree* flattree = new TTree("flattree","flattree");
    TTree* tree = new TTree("tree","tree");
    float phi;
    float Psi;
    int event = 0;
    //flattree->Branch("phi", &phi, "phi/F");
    //flattree->Branch("Psi", &Psi, "Psi/F");
    //flattree->Branch("event",&event,"event/I");

    tree->Branch("n", &b_n, "n/I");
    tree->Branch("psig", &b_psig, "psig/F");
    tree->Branch("ptg", &b_ptg, "ptg[n]/F");  
    tree->Branch("etag", &b_etag, "etag[n]/F");
    tree->Branch("phig", &b_phig, "phig[n]/F");

    TF1 *EtaDistr = new TF1("EtaDistr","exp(-(x-2.1)^2/6.3)+exp(-(x+2.1)^2/6.3)",-4,4);
    //TF1 *PhiDistr = new TF1("PhiDistr","1+2*[0]*cos(x)+2*[1]*cos(2*x)+2*[2]*cos(3*x)+2*[3]*cos(4*x)+2*[4]*cos(5*x)+2*[5]*cos(6*x)",0,2.*TMath::Pi());
    TF1 *PhiDistr = new TF1("PhiDistr","1+2*[1]*cos(x-[0])+2*[2]*cos(2*(x-[0]))+2*[3]*cos(3*(x-[0]))+2*[4]*cos(4*(x-[0]))+2*[5]*cos(5*(x-[0]))+2*[6]*cos(6*(x-[0]))",-pi,pi);
    //TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/.40))+0.0015*exp (-(x/1.5))", 0.2,10);  //V~0.12
    //TF1 *PtDistr  = new TF1("PtDistr","exp (-(x/0.90))+0.15*exp (-(x/15))", 0.1,10);    //V~=0.06
    TF1 *PtDistr  = new TF1("PtDistr","0.03*(exp (-(x/0.594540))+0.00499506*exp (-(x/1.89391)))", 0.2,6.0);       //Real Data
    //  TF1 *PtDistr = new TF1("PtDistr","[0]*x*TMath::Power(1+(sqrt(x*x+[1]*[1])-[1]/[2]),-[3])",0.2,10);
        //PtDistr->SetParameters(118.836,-0.335972,0.759243,118.836);   //Real data fit with Tsallis
    //TF1 *V1vsEta = new TF1("V1vsEta","-exp(-(x+1)^2)/20-x/30+exp(-(x-1)^2)/20",-4,4); 
    TF1* V1vsEta = new TF1("V1vsEta","exp(-(x+1)^2)/20+exp(-(x-1)^2)/20", -4, 4);
    //TF1 *V2vsPt   = new TF1("V2vsPt","((x/3)^1.8/(1+(x/3)^1.8))*(.00005+(1/x)^0.8)",0.2,10);
    //TF1 *V2vsPt = new TF1("V2vsPt","((x/[0])^[1]/(1+(x/[2])^[3]))*(.00005+(1/x)^[4])",0.1,10);
    //     V2vsPt->SetParameters(4.81159,1.80783,3.69272,3.11889,0.931485);        //Real data V~0.05
    //    V2vsPt->SetParameters(5,1.8,3,1.8,0.8); //V~0.06
    TF1 *V2vsPt = new TF1("V2vsPt","((x/3.31699)^2.35142/(1+(x/3.49188)^3.54429))*(.00005+(1/x)^1.50600)",0.2,6.0);
    TF1 *V3vsPt = new TF1("V3vsPt","((x/3.2)^2.3/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.2,6.0);
    TF1 *V4vsPt = new TF1("V4vsPt","((x/4.8)^2.1/(1+(x/3.4)^2.1))*(.00005+(1/x)^1.4)",0.2,6.0);
    TF1 *V5vsPt = new TF1("V5vsPt","((x/6.0)^3.2/(1+(x/11.4)^2.1))*(.00005+(1/x)^1.4)",0.2,6.0);
    TF1 *V6vsPt = new TF1("V6vsPt","((x/5.6)^2.4/(1+(x/4.7)^2.1))*(.00005+(1/x)^1.4)",0.2,6.0);
    /*
    DistSampler * sampler = Factory::CreateDistSampler();
    if (sampler == 0) {
        Info("multidimSampling","Default sampler %s is not available try with Foam ",
           ROOT::Math::DistSamplerOptions::DefaultSampler().c_str() );
        ROOT::Math::DistSamplerOptions::SetDefaultSampler("Foam");
    }
    sampler = Factory::CreateDistSampler();
    if (sampler == 0) {
        Error("multidimSampling","Foam sampler is not available - exit ");
        return;
    }
    */
    for (int i = 0; i < N; i ++){
        if(i%10==0)
            cout << i << " events" << endl;
        UInt_t iniseed = SetSeedOwn();
        gRandom->SetSeed(iniseed);
        Psi = r->Uniform(-pi,pi);
        b_n = (int)r->Gaus(M, RMSM);
        int b_nsame = (int)b_n * 0.1;
        int b_nback = (int)b_n * 0.1;
        //PhiDistr->SetParameters(Psi,v1,v2,v3,v4,v5,v6);
        /*
        sampler->SetFunction(*fcos, 1);
        sampler->SetRange(-pi,pi);
        bool ret = sampler->Init();
        if (!ret) {
            Error("Sampler::Init","Error initializing unuran sampler");
            return;
        }
        */
        for(int j = 0; j < b_n; j++){
            float pt = PtDistr->GetRandom();
            float eta =  EtaDistr->GetRandom();

            float v1=V1vsEta->Eval(eta);
            float v2=V2vsPt->Eval(pt);
            float v3=V3vsPt->Eval(pt);
            float v4=V4vsPt->Eval(pt);
            float v5=V5vsPt->Eval(pt);
            float v6=V6vsPt->Eval(pt);
            
            b_ptg[j]  = pt;
            b_etag[j] = eta;

            PhiDistr->SetParameters(Psi,v1,v2,v3,v4,v5,v6);
            phi = PhiDistr->GetRandom();
            b_phig[j] = phi;
            b_psig = Psi;
            //phi = sampler->Sample1D();
            //flattree->Fill();
        }
        for(int j=0;j<b_nsame;j++){
            int k = r->Integer(b_n);
            b_ptg[b_n] = b_ptg[k];
            b_etag[b_n] = b_etag[k]+2*(k%2-0.5)*1e-3;
            b_phig[b_n] = b_phig[k]+2*(k%2-0.5)*1e-3;
            b_n++;
            //flattree->Fill();
        }
        for(int j=0;j<b_nback;j++){
          int k = r->Integer(b_n-b_nsame);
          b_ptg[b_n] = b_ptg[k];
          b_etag[b_n] = -(b_etag[k]+2*(k%2-0.5)*1e-3);
          b_phig[b_n] = TMath::Pi()+(b_phig[k]+2*(k%2-0.5)*1e-3);
          if (b_phig[b_n]>pi) b_phig[b_n]=b_phig[b_n]-2.*pi; // -pi ~ pi
          b_n++;
        } // End of loop over particles
        tree->Fill();
        event++;
    }
    
    //flattree->Write();
    tree->Write();
}

