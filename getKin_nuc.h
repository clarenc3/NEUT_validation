#include <iostream>
#include <iomanip>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TRandom3.h"

#include "TSystem.h"

#include "neutpart.h"
#include "neutfsipart.h"
#include "neutvect.h"
#include "neutfsivert.h"
#include "neutvtx.h"
#include "neutrootTreeSingleton.h"

void getKin(std::string fileName, int fitType);

void Usage(void);

bool isCC1ppip(NeutVect *nvect, double EnuMax);
bool isCC1npip(NeutVect *nvect, double EnuMax);
bool isCC1pi0(NeutVect *nvect, double EnuMax);
bool isCC1ppim(NeutVect *nvect, double EnuMax);
bool isCC1npim(NeutVect *nvect, double EnuMax);

TH1D* FluxUnfoldedScaling(TH1D* mcHist, TH1D* fluxHist);

double Q2CCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding=0.);
double EnuCCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding=0.);
double Wrec(TLorentzVector pnu, TLorentzVector pmu, double EHad=0, double binding=0);
double MpPi(TLorentzVector pp, TLorentzVector ppi);
