#include <iostream>
#include <fstream>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include <cstdlib>
#include <math.h>       
#include <cmath>        
#include <stdlib.h>  
#include <stdio.h>
#include <string.h>
#include "TVector3.h"
#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TProfile2D.h"
#include "TProfile.h"
#include "TGraph.h"
#include "TBrowser.h"
#include "TRandom.h"
#include <sstream>
#include <algorithm>

using namespace std;


struct Proton{
  Double_t x0,y0,z0,px0,py0,pz0;
  Double_t x1,y1,z1,px1,py1,pz1;
  Double_t Einit,Estop;
  vector<double> *tracks_X = new vector<double>;
  vector<double> *tracks_E = new vector<double>;
  Int_t   Id;
};

vector<double> Energy;
vector<double> dEdXBins;


double findWET(double, double);
void ComputeSpline(Proton *, TH3D*, TH2D*, TH2D*);

int main(int argc, char** argv){
  TApplication theApp("theApp", &argc, argv);
  gROOT->ProcessLine("#include<vector>");
  Proton Point;
  char *FileName = argv[1];
  

  TFile* phaseFile = new TFile(FileName,"update");


  TTree* t = (TTree*)phaseFile->Get("phase");
  t->SetBranchAddress("x0",&Point.x0);
  t->SetBranchAddress("y0",&Point.y0);
  t->SetBranchAddress("z0",&Point.z0);
 
  t->SetBranchAddress("px0",&Point.px0);
  t->SetBranchAddress("py0",&Point.py0);
  t->SetBranchAddress("pz0",&Point.pz0);

  t->SetBranchAddress("x1",&Point.x1);
  t->SetBranchAddress("y1",&Point.y1);
  t->SetBranchAddress("z1",&Point.z1);

  t->SetBranchAddress("px1",&Point.px1);
  t->SetBranchAddress("py1",&Point.py1);
  t->SetBranchAddress("pz1",&Point.pz1);

  t->SetBranchAddress("Einit",&Point.Einit);
  t->SetBranchAddress("Estop",&Point.Estop);
  t->SetBranchAddress("tracks_E",&Point.tracks_E);
  t->SetBranchAddress("tracks_X",&Point.tracks_X);
  t->SetBranchAddress("Id",&Point.Id);
 
  t->GetEntry(0);
 

 /* double xmin  = t->GetMinimum("x0");
  double xmax  = t->GetMaximum("x1");
  double ymin  = t->GetMinimum("y1");
  double ymax  = t->GetMaximum("y1");
  double zmin  = t->GetMinimum("z1");
  double zmax  = t->GetMaximum("z1");*/


  //Initialising the theta map
  TH3D* T3DMap = new TH3D("Edep","Edep",512,0,20,200,0,15,200,0,15);
  TH3D* Theta3DMap = new TH3D("Theta3DMap","Theta binned [-.5,.5]", 200, 0.,15., 200, 0.,15., 50, -.5, .5);
  TH2D* Integral       = new TH2D("Integral","Weighted Mean Integral",200,0,15,200,0,15);
  TH2D* weights        = new TH2D("weights","Weights",200,0,15,200,0,15);
  

  std::string line;
  std::ifstream SPWater("Water_Geant4.dat");  
  double data[2];
  while(getline(SPWater, line)) {
    stringstream ss(line);
    for(int i=0;i<2;i++) ss >> data[i];
    Energy.push_back(data[0]);
    dEdXBins.push_back(data[1]);
  }
  
 

  int NEntries = t->GetEntries();
  //Loop over all protons
  for(int i=0;i<NEntries;i++){

    t->GetEntry(i); 
    if(i%100000 == 0) cout<<i<<endl;
    //if(i==100000) break;   

    //Converting to length in cm
    transform(Point.tracks_X->begin(), Point.tracks_X->end(), Point.tracks_X->begin(), bind1st(multiplies<double>(),0.1));
    transform(Point.tracks_X->begin(), Point.tracks_X->end(), Point.tracks_X->begin(), bind2nd(plus<double>(), 10.0));
    
    Point.x0 = 0.1*Point.x0+10.; Point.y0 = 0.1*Point.y0+7.5; Point.z0 = 0.1*Point.z0+7.5;
    Point.x1 = 0.1*Point.x1+10.; Point.y1 = 0.1*Point.y1+7.5; Point.z1 = 0.1*Point.z1+7.5;

    
    //Cuts
    if(abs(Point.y1-Point.y0)>5 || abs(Point.z1-Point.z0)>5) continue;
    if(Point.Estop<80) continue; 

    
    //Calculating the variance map
    double theta = atan2(Point.py1,Point.px1)-atan2(Point.py0,Point.px0);
    Theta3DMap->Fill(Point.y1, Point.z1, theta);

    ComputeSpline(&Point,T3DMap, Integral, weights);
    

  }

  Theta3DMap->FitSlicesZ(0,1,0,1,0,0,"QN");

  TH2D* h3_2 = (TH2D*)gDirectory->Get("Theta3DMap_2");
  h3_2->Multiply(h3_2);//Squaring to get variance from stdev
  h3_2->Write("ThetafromMap",TObject::kOverwrite); 
  
  Integral->Divide(weights);
  Integral->Write("Integral", TObject::kOverwrite);

  h3_2->Divide(Integral);
  h3_2->Write("HighlandValueTMap",TObject::kOverwrite);

  phaseFile->Close();
  return 0;
 
}

////////////////////////////////////////////
// Extract WET
////////////////////////////////////////////
double findWET(double Einit,double Estop){
  int it_Einit = lower_bound(Energy.begin(), Energy.end(), Einit)-Energy.begin();
  int it_Estop = lower_bound(Energy.begin(), Energy.end(), Estop)-Energy.begin();
  double WET = 0 ;
  for(int i=it_Estop;i<it_Einit;i++){
    WET += 1./dEdXBins[i];
  }
  return WET/10.;
}

////////////////////////////////////////////
// Compute Spline
////////////////////////////////////////////
 
void ComputeSpline(Proton *Point, TH3D* T3DMap, TH2D* Integral, TH2D* weights){
  TVector3 P0(Point->px0,Point->py0,Point->pz0);
  TVector3 P1(Point->px1,Point->py1,Point->pz1);
  TVector3 X0(Point->x0,Point->y0,Point->z0);
  TVector3 X1(Point->x1,Point->y1,Point->z1);

  double WEPL = 0.00244*pow(Point->Einit,1.75);
  double WET  = findWET(Point->Einit,Point->Estop);

 

  //Calculate WET and WEPL and compute the spline 
  double Einit = Point->Einit;
  double Estop = Point->Estop;
  double Lambda1 = 1.01+0.43*pow(WET/WEPL,2);
  double Lambda2 = 0.99-0.46*pow(WET/WEPL,2);
  double Length = TVector3(X1-X0).Mag();
  P0.SetMag(Length*Lambda1);
  P1.SetMag(Length*Lambda2);
  TVector3 A,B,C,D;
  TVector3 X;

  A     =  2*X0 - 2*X1 + P0 + P1;
  B     = -3*X0 + 3*X1 - 2*P0 - P1;
  C     =  P0;
  D     =  X0;

  TVector3 X_old = X0;

  double TotL      = 0;
  std::map<pair<int,int>,double> Lengthmap;
  pair<map<pair<int,int>,double>::iterator,bool> ret;
  int binx,biny,binz;

  //Calculate the Integral of the 1/(pv)^2 function 
  double tau1, tau2, pv1, pv2;
  double Length1, Length2, IntStep;
  double integral  = 0;
  double TotInt    = 0;
  double MassP     = 938.27;
   
  tau1 	  = Point->tracks_E->at(1)/MassP;
  pv1 	  = pow((tau2+1)/(Point->tracks_E->at(1)*(tau2+2)),2);

  for(int i=0;i<Point->tracks_X->size()-1;i++){

    double t    = (Point->tracks_X->at(i)-Point->x0)/(Point->x1-Point->x0);
    X           = A*(t*t*t) + B*(t*t) + C*t + D;
    int binId   = T3DMap->FindBin(X.x(), X.y(),X.z());

    T3DMap->GetBinXYZ(binId,binx,biny,binz);
    double L    = TVector3(X-X_old).Mag(); 

    X_old       = X;     

    tau2        = Point->tracks_E->at(i+1)/MassP;
    pv2	        = pow((tau2+1)/(Point->tracks_E->at(i+1)*(tau2+2)),2);
    IntStep     = Point->tracks_X->at(i+1) - Point->tracks_X->at(i);
    Length1     = Point->tracks_X->at(i) - Point->x1;
    Length2     = Point->tracks_X->at(i+1) - Point->x1;
     
    integral    = (pv2 + pv1)*IntStep/2;
    pv1 	= pv2;
    TotInt     += integral;
      
    
    //Mapping the length crossed in each column
    pair<int,int> bin2dID = pair<int,int>(biny,binz);
    ret 	= Lengthmap.insert(pair<pair<int,int>,double>(bin2dID,L));
    if ( !ret.second ) Lengthmap[bin2dID]  += L;
    TotL+=L;      

  }
    
  //Weighting the integral
  std::map<pair<int,int>,double>::iterator it;
  for(it = Lengthmap.begin(); it != Lengthmap.end(); it++){
    biny = it->first.first;
    binz = it->first.second;
    double L = it->second;
    double y = T3DMap->GetYaxis()->GetBinCenter(biny);
    double z = T3DMap->GetZaxis()->GetBinCenter(binz);
    Integral->Fill(y,z,TotInt*L/TotL);
    weights->Fill(y,z,L/TotL);
  }
}

