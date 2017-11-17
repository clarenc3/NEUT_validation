#include "getKin_mode_validation.h"

// The main
int main(int argc, char* argv[]) {

  if (argc != 5) {
    Usage(argc, argv);
    exit(-1);
  }

  getKin(std::string(argv[1]), std::atoi(argv[2]), std::atof(argv[3]), std::string(argv[4]));

  return 0;
};

// Print usage
void Usage(int argc, char *argv[]) {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << "You gave " << argc << " arguments" << std::endl;
  for (int i = 0; i < argc; ++i) {
    std::cout << "   " << argv[i] << std::endl;
  }
  std::cout << "./getKin_mode_validation.exe ROOT_FILE SIGNAL MAX_MOM INPUT_REWEIGHT_FILE" << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "SIGNAL is 0 = CC0pi, 1 = CC1pi, 2 = CC1pi0" << std::endl;
  std::cout << "MAX_MOM is the rough momentum scale of the neutrinos" << std::endl;
  std::cout << "INPUT_REWEIGHT_FILE is the reweighting file" << std::endl;

}

// The main loop
void getKin(std::string fileName, int sigType, double maxMom, std::string inputFile) {

  if (sigType == kCC0pi)        std::cout << "Signal is CC0pi" << std::endl;
  else if (sigType == kCC1pip)  std::cout << "Signal is CC1pi+" << std::endl;
  else if (sigType == kCC1pi0)  std::cout << "Signal is CC1pi0" << std::endl;
  else {
    std::cout << "Unrecognised signal! Exiting..." << std::endl;
    char *blarb[] = {"hi", "bye"};
    Usage(0, blarb);
    exit(-1);
  }

  // Open the input file with the reweighting histogram
  TFile *inputfile = new TFile(inputFile.c_str(), "open");

  // Open the input NEUT file
  TFile *f = TFile::Open((fileName).c_str(),"open");
  //f->ls();
  TTree *tn = (TTree*)(f->Get("neuttree"));
  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch", &nvect);

  TH1D *fluxHist = NULL;
  TH1D *eventHist = NULL;

  if (f->Get("flux_numu")) {
    fluxHist = (TH1D*)(f->Get("flux_numu")->Clone());
    eventHist = (TH1D*)(f->Get("evtrt_numu")->Clone());
  } else if (f->Get("flux_numub")) {
    fluxHist = (TH1D*)(f->Get("flux_numub")->Clone());
    eventHist = (TH1D*)(f->Get("evtrt_numub")->Clone());
  } else if (f->Get("flux_nue")) {
    fluxHist = (TH1D*)(f->Get("flux_nue")->Clone());
    eventHist = (TH1D*)(f->Get("evtrt_nue")->Clone());
  } else if (f->Get("flux_nueb")) {
    fluxHist = (TH1D*)(f->Get("flux_nueb")->Clone());
    eventHist = (TH1D*)(f->Get("evtrt_nueb")->Clone());
  }
  long int nevents = tn->GetEntries();
  double ScaleFactor = eventHist->Integral("width")*1E-38/(nevents*fluxHist->Integral("width"));

  // Restricted phase space
  bool isRestricted = false;

  // Array of the modes
  std::vector<int> ModeArray;
  ModeArray.reserve(7);
  // CC
  ModeArray.push_back(11);
  ModeArray.push_back(12);
  ModeArray.push_back(13);
  // NC
  ModeArray.push_back(31);
  ModeArray.push_back(32);
  ModeArray.push_back(33);
  ModeArray.push_back(34);

  //const double MaxMom = 10.00;
  const unsigned int nModes = 50;

  // Pmu CosThetaMu for each mode
  std::vector<TH2D*> PmuCosmu;
  std::vector<TH3D*> PmuCosmuEnu;

  // The weighted histograms
  std::vector<TH2D*> PmuCosmu_weight;
  std::vector<TH3D*> PmuCosmuEnu_weight;

  PmuCosmu.reserve(nModes);
  PmuCosmuEnu.reserve(nModes);

  PmuCosmu_weight.reserve(nModes);
  PmuCosmuEnu_weight.reserve(nModes);

  std::vector<double> PmuCosmu_max;
  std::vector<double> PmuCosmuEnu_max;
  PmuCosmu_max.reserve(nModes);
  PmuCosmuEnu_max.reserve(nModes);

  //inputfile->ls();

  // Initialise the histograms
  for (size_t i = 0; i < nModes; ++i) {
    std::stringstream ss;
    ss << "PmuCosmu_" << i;
    std::stringstream ss2;
    ss2 << "PmuCosmuEnu_" << i;

    // Get the weight histogram for each mode
    std::string title = ss.str()+"_NEUT533_vs_"+ss.str()+"_Minoo";
    std::string title2 = ss2.str()+"_NEUT533_vs_"+ss2.str()+"_Minoo";
    if (inputfile->Get(title.c_str())) {
      std::cout << "Found " << title << std::endl;
      TH2D *temp = NULL;
      inputfile->GetObject(title.c_str(), temp);
      PmuCosmu_weight[i] = (TH2D*)(temp->Clone());
      delete temp;
    } else {
      PmuCosmu_weight[i] = NULL;
    }

    if (inputfile->Get(title2.c_str())) {
      std::cout << "Found " << title2 << std::endl;
      TH3D *temp = NULL;
      inputfile->GetObject(title2.c_str(), temp);
      PmuCosmuEnu_weight[i] = (TH3D*)(temp->Clone());
      delete temp;
    } else {
      PmuCosmuEnu_weight[i] = NULL;
    }

    if (PmuCosmu_weight[i] != NULL) {
      // Make PmuCosmu and PmuCosmuEnu just clones
      PmuCosmu[i] = (TH2D*)(PmuCosmu_weight[i]->Clone());
      PmuCosmu[i]->Reset();
      PmuCosmu[i]->SetNameTitle(ss.str().c_str(), ss.str().c_str());

      PmuCosmuEnu[i] = (TH3D*)(PmuCosmuEnu_weight[i]->Clone());
      PmuCosmuEnu[i]->Reset();
      PmuCosmuEnu[i]->SetNameTitle(ss2.str().c_str(), ss2.str().c_str());
    } else {
      PmuCosmu[i] = NULL;
      PmuCosmuEnu[i] = NULL;
    }

    if (PmuCosmu[i] != NULL) {
      PmuCosmu[i]->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
      PmuCosmu[i]->GetYaxis()->SetTitle("cos#theta_{#mu}");
      PmuCosmu[i]->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#mu}dcos#theta_{#mu} (cm^{2}/GeV/1/nucleon)");
      PmuCosmu_max[i] = PmuCosmu_weight[i]->GetBinContent(PmuCosmu_weight[i]->GetMaximumBin());
    }

    if (PmuCosmuEnu[i] != NULL) {
      PmuCosmuEnu[i]->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
      PmuCosmuEnu[i]->GetYaxis()->SetTitle("cos#theta_{#mu}");
      PmuCosmuEnu[i]->GetZaxis()->SetTitle("E_{#nu} (GeV)");
      PmuCosmuEnu_max[i] = PmuCosmuEnu_weight[i]->GetBinContent(PmuCosmuEnu_weight[i]->GetMaximumBin());
    }
  }

  int eventCnt = 0;
  std::cout << "Total number of events: " << nevents << std::endl;

  TStopwatch clock;
  clock.Start();


  int countwidth = int(nevents/20);
  // Loop over the events
  for (int j = 0; j < nevents; j++) {

    tn->GetEntry(j);

    if (j % countwidth == 0) {
      std::cout << "On event #" << j << "/" << nevents << " (" << int(double(j)*100./double(nevents)) << "%)" << std::endl;
    }

    // Very simple signal definition
    if (sigType == kCC0pi) {
      if (!isT2K_CC0pi(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
      // Only plot the pion modes
    } else if (sigType == kCC1pip) {
      // Check if the event is the mode we want
      bool found_mode = false;
      for (size_t i = 0; i < ModeArray.size(); ++i) {
        if (abs(nvect->Mode) != ModeArray[i]) continue;
        found_mode = true;
      }
      if (!found_mode) continue;
      //if (!isT2K_CC1pip(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
    } else if (sigType == kCC1pi0) {
      if (!isT2K_CC1pi0(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
    }

    TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
    TLorentzVector Pmu(0.0, 0.0, 0.0, 0.0);
    TLorentzVector PnuOut(0.0, 0.0, 0.0, 0.0);

    // Loop over the particle stack
    for (int k = 2; k < nvect->Npart(); ++k) {
      // Get the PID of the particle
      int PID = (nvect->PartInfo(k))->fPID;
      // Check if it's the outgoing lepton
      if (fabs(PID) == fabs(nvect->PartInfo(0)->fPID)-1) {
        if (!(nvect->PartInfo(k))->fIsAlive && (nvect->PartInfo(k))->fStatus != 0) continue;
        Pmu = (nvect->PartInfo(k))->fP;  
      } else if (PID == nvect->PartInfo(0)->fPID) {
        PnuOut = nvect->PartInfo(k)->fP;
      }
    } // Finish the for loop over particles

    // Get the muon variables
    double plep = -999.999;
    double costhlep = -999.999;
    double enu = Pnu.E()/1000.;

    int mode = abs(nvect->Mode);
    if (mode == 11 || mode == 12 || mode == 13) {
      plep = Pmu.Vect().Mag()/1000.;
      if (plep == 0.0) continue;
      costhlep = cos(Pnu.Vect().Angle(Pmu.Vect()));
    } else if (mode == 31 || mode == 32 || mode == 33 || mode == 34) {
      plep = PnuOut.Vect().Mag()/1000.;
      if (plep == 0.0) continue;
      costhlep = cos(Pnu.Vect().Angle(PnuOut.Vect()));
    }

    // Find which bin this plep costhlep enu combination lives in
    int PmuCosmu_bin = PmuCosmu_weight[abs(nvect->Mode)]->FindFixBin(plep, costhlep);
    int PmuCosmuEnu_bin = PmuCosmuEnu_weight[abs(nvect->Mode)]->FindFixBin(plep, costhlep, enu);

    // Make sure it's the same bin in the two plots
    int PmuCosmu_bin_ref = PmuCosmu[abs(nvect->Mode)]->FindFixBin(plep, costhlep);
    int PmuCosmuEnu_bin_ref = PmuCosmuEnu[abs(nvect->Mode)]->FindFixBin(plep, costhlep, enu);

    // Make sure we're picking up the correct bin
    if (PmuCosmu_bin != PmuCosmu_bin_ref) {
      std::cerr << "PmuCosmu_bin: " << PmuCosmu_bin << std::endl;
      std::cerr << "PmuCosmu_bin_ref: " << PmuCosmu_bin_ref << std::endl;
      std::cerr << "PmuCosmu global bin numbering don't match!" << std::endl;
    }

    if (PmuCosmuEnu_bin != PmuCosmuEnu_bin_ref) {
      std::cerr << "PmuCosmuEnu_bin: " << PmuCosmuEnu_bin << std::endl;
      std::cerr << "PmuCosmuEnu_bin_ref: " << PmuCosmuEnu_bin_ref << std::endl;
      std::cerr << "PmuCosmuEnu global bin numbering don't match!" << std::endl;
    }

    // Get the bin content of said global bin
    double weight_2d = PmuCosmu_weight[abs(nvect->Mode)]->GetBinContent(PmuCosmu_bin);
    double weight_3d = PmuCosmuEnu_weight[abs(nvect->Mode)]->GetBinContent(PmuCosmuEnu_bin);

    //std::cout << 1./weight_2d << " " << 1./weight_3d << std::endl;

    if (weight_2d > PmuCosmu_max[abs(nvect->Mode)] || weight_3d > PmuCosmuEnu_max[abs(nvect->Mode)]) {
      std::cout << std::endl;
      std::cout << nvect->Mode << std::endl;
      std::cout << plep << " " << costhlep << std::endl;
      std::cout << "2d " << weight_2d << " > " << PmuCosmu_max[abs(nvect->Mode)] << std::endl;
      std::cout << "3d " << weight_3d << " > " << PmuCosmuEnu_max[abs(nvect->Mode)] << std::endl;
    }

    // Find the weight to apply to the nominal histogram
    if (weight_2d != 0) {
      PmuCosmu[abs(nvect->Mode)]->Fill(plep, costhlep, 1./weight_2d);
    }
    if (weight_3d != 0) {
      PmuCosmuEnu[abs(nvect->Mode)]->Fill(plep, costhlep, enu, 1./weight_3d);
    }
  }

  std::cout << eventCnt << "/" << nevents << " events" << std::endl;

  // Make filename for selection
  if      (sigType == kCC0pi)  fileName += "_CC0pi";
  else if (sigType == kCC1pip) fileName += "_CC1pip";
  else if (sigType == kCC1pi0) fileName += "_CC1pi0";

  fileName += "_1pi_modes_outgoing_validation.root";

  TFile *output = new TFile(fileName.c_str(), "recreate");
  output->cd();
  for (size_t i = 0; i < nModes; ++i) {
    if (PmuCosmu[i] == NULL) continue;

    if (PmuCosmu[i]->Integral() > 0) {
      PmuCosmu[i]->Scale(ScaleFactor, "width");
      std::cout << "PmuCosmu integral: " << PmuCosmu[i]->Integral() << std::endl;
      PmuCosmu[i]->Write();
    }

    if (PmuCosmu_weight[i]->Integral() > 0) {
      PmuCosmu_weight[i]->Write();
    }

    if (PmuCosmuEnu[i]->Integral() > 0) {
      PmuCosmuEnu[i]->Scale(ScaleFactor, "width");
      std::cout << "PmuCosmuEnu integral: " << PmuCosmuEnu[i]->Integral() << std::endl;
      PmuCosmuEnu[i]->Write();
    }

    if (PmuCosmuEnu_weight[i]->Integral() > 0) {
      PmuCosmuEnu_weight[i]->Write();
    }

  }
  output->Close();


  std::cout << "Wrote to " << fileName << std::endl;

  return;
};


// *********************************************
bool isT2K_CC0pi(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted) {
  // *********************************************

  // Make sure it's neutrino
  if ((nvect->PartInfo(0))->fPID != 14) return false; 

  // Make sure it's CC1mu
  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false;

  int lepCnt = 0;

  TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
  TLorentzVector Pmu;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue;

    int PID = (nvect->PartInfo(j))->fPID;
    if ((abs(PID) >= 111 && abs(PID) <= 210) || (abs(PID) >= 212 && abs(PID) <= 557) || abs(PID) == 211) { 
      return false; 
    } else if (abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15 || abs(PID) == 17) {
      lepCnt++;
      Pmu = (nvect->PartInfo(j))->fP;
    }
  }

  if (lepCnt != 1) return false;

  // relatively generic CC1pi+ definition done
  // now measurement specific (p_mu > 200 MeV, p_pi > 200 MeV, cos th_mu > 0.3,
  // cos th_pi > 0.3 in TRUE AND RECONSTRUCTED!

  if (restricted) {
    double p_mu = Pmu.Vect().Mag();
    double cos_th_mu = cos(Pmu.Vect().Angle(Pnu.Vect()));

    if (p_mu <= 200 || cos_th_mu <= 0.2) return false;
  }

  return true;
};

// *********************************************
bool isT2K_CC1pip(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted) {
  // *********************************************

  if ((nvect->PartInfo(0))->fPID != 14) return false; 

  int pipCnt = 0;
  int lepCnt = 0;

  TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
  TLorentzVector Ppip;
  TLorentzVector Pmu;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue;

    int PID = (nvect->PartInfo(j))->fPID;

    if ((abs(PID) >= 111 && abs(PID) <= 210) || 
        (abs(PID) >= 212 && abs(PID) <= 557) || 
        PID == -211) {
      return false; 

    } else if ( abs(PID) == 11 ||
        abs(PID) == 13 ||
        abs(PID) == 15 ||
        abs(PID) == 17) {
      lepCnt++;
      Pmu = (nvect->PartInfo(j))->fP;

    } else if (PID == 211) {
      pipCnt++;
      Ppip = (nvect->PartInfo(j))->fP;
    }
  }

  if (pipCnt != 1) return false; 
  if (lepCnt != 1) return false;

  // relatively generic CC1pi+ definition done
  // now measurement specific (p_mu > 200 MeV, p_pi > 200 MeV, cos th_mu > 0.2, cos th_pi > 0.2

  if (restricted) {
    double p_mu = Pmu.Vect().Mag();
    double p_pi = Ppip.Vect().Mag();
    double cos_th_mu = cos(Pnu.Vect().Angle(Pmu.Vect()));
    double cos_th_pi = cos(Pnu.Vect().Angle(Ppip.Vect()));

    if (p_mu <= 200 || p_pi <= 200 || cos_th_mu <= 0.2 || cos_th_pi <= 0.2) return false;
  }


  return true;
};

// *********************************************
bool isT2K_CC1pi0(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted) {
  // *********************************************

  if ((nvect->PartInfo(0))->fPID != 14) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false;

  int pi0Cnt = 0;
  int lepCnt = 0;

  TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
  TLorentzVector Ppi0;
  TLorentzVector Pmu;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue;
    int PID = (nvect->PartInfo(j))->fPID;
    if ((abs(PID) > 111 && abs(PID) <= 210) || (abs(PID) >= 212 && abs(PID) <= 557) || abs(PID) == 211) {
      return false; 
    } else if (abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15 || abs(PID) == 17) {
      lepCnt++;
      Pmu = (nvect->PartInfo(j))->fP;
    } else if (PID == 111) {
      pi0Cnt++;
      Ppi0 = (nvect->PartInfo(j))->fP;
    }
  }

  if (pi0Cnt != 1) return false; 
  if (lepCnt != 1) return false;

  // relatively generic CC1pi+ definition done
  // now measurement specific (p_mu > 200 MeV, p_pi > 200 MeV, cos th_mu > 0.2, cos th_pi > 0.2

  if (restricted) {
    double p_mu = Pmu.Vect().Mag();
    double p_pi = Ppi0.Vect().Mag();
    double cos_th_mu = cos(Pnu.Vect().Angle(Pmu.Vect()));
    double cos_th_pi = cos(Pnu.Vect().Angle(Ppi0.Vect()));

    if (p_mu <= 200 || p_pi <= 200 || cos_th_mu <= 0.2 || cos_th_pi <= 0.2) return false;
  }


  return true;
};


