#include <iostream>
#include <cmath>

#include <omp.h>

#include "TStopwatch.h"
#include "TThread.h"

#include "TH1D.h"
#include "TH2D.h"

#include "TFile.h"
#include "TChain.h"
#include "TTree.h"

#include "TLorentzVector.h"

#include "TSystem.h"

#include "neutpart.h"
#include "neutfsipart.h"
#include "neutvect.h"
#include "neutfsivert.h"
#include "neutvtx.h"
#include "neutrootTreeSingleton.h"

void getKin(std::string fileName, int fitType, double maxMom);

void Usage(void);

// Signal definitions
bool isT2KCC0pi(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted);
bool isT2KCCpip(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted);
bool isT2KCCpi0(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted);

// Reconstruction variables
double Q2QErec(TLorentzVector pnu, TLorentzVector pmu, double binding=30.);
double EnuQErec(TLorentzVector pnu, TLorentzVector pmu, double binding=30.);
double Q2CCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding=0.);
double EnuCCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding=0.);
double Wrec(TLorentzVector pnu, TLorentzVector pmu, double EHad=0, double binding=0);
