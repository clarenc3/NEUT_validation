#include <iostream>
#include <cmath>
#include <sstream>

#include "TStopwatch.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "THStack.h"

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

// Additionally need a TFile to take the weights from
void getKin_weight(std::string fileName, int fitType, double maxMom, std::string WeightFile, std::string HistoName);

void Usage(void);

// And enum to keep track of signal definiton
enum SignalDefinition {
  kCC0pi = 0,
  kCC1pip = 1,
  kCC1pi0 = 2
};

// Signal definitions
bool isT2K_CC0pi(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted);
bool isT2K_CC1pip(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted);
bool isT2K_CC1pi0(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted);

// Reconstruction variables
double Q2CCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip);
double Q2true(TLorentzVector pnu, TLorentzVector pmu);

double EnuCCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip);
double Wrec(TLorentzVector pnu, TLorentzVector pmu);
double Wtrue(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector pprim);
