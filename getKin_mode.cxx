#include "getKin.h"

// The main
int main(int argc, char* argv[]) {

  if (argc != 4) {
    Usage();
    exit(-1);
  }

  getKin(std::string(argv[1]), std::atoi(argv[2]), std::atof(argv[3]));

  return 0;
};

// Print usage
void Usage() {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << "./getKin.exe ROOT_FILE SIGNAL MAX_MOM" << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "SIGNAL is 0 = CC0pi, 1 = CC1pi, 2 = CC1pi0" << std::endl;
  std::cout << "MAX_MOM is the rough momentum scale of the neutrinos" << std::endl;

  return;
}

// The main loop
void getKin(std::string fileName, int sigType, double maxMom) {

  if (sigType == kCC0pi)        std::cout << "Signal is CC0pi" << std::endl;
  else if (sigType == kCC1pip)  std::cout << "Signal is CC1pi+" << std::endl;
  else if (sigType == kCC1pi0)  std::cout << "Signal is CC1pi0" << std::endl;
  else {
    std::cout << "Unrecognised signal! Exiting..." << std::endl;
    Usage();
    exit(-1);
  }

  TFile *f = TFile::Open((fileName).c_str(),"open");
  f->ls();
  TTree *tn = (TTree*)(f->Get("neuttree"));
  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch",&nvect);

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

  const double MaxMom = 10.00;
  const unsigned int nModes = 50;

  // Pmu CosThetaMu for each mode
  TH2D *PmuCosmu[nModes];
  TH3D *PmuCosmuEnu[nModes];

  for (size_t i = 0; i < nModes; ++i) {
    std::stringstream ss;
    ss << "PmuCosmu_" << i;
    std::stringstream ss2;
    ss2 << "PmuCosmuEnu_" << i;
    switch (i) {
      // CC1pi+1p
      case 11:
      case 12:
      case 13:
        {
          // Make custom arrays, so boring
          std::vector<double> Pmu;
          double mom = 0.0;
          while (mom < MaxMom) {
            Pmu.push_back(mom);
            if (mom <= 0.2) {
              mom += 0.07;
            } else if (mom > 0.2 && mom <= 0.75) {
              mom += 0.05;
            } else if (mom > 0.75 && mom <= 1.50) {
              mom += 0.10;
            } else if (mom > 1.50 && mom <= 4.00) {
              mom += 0.25;
            } else if (mom > 4.00 && mom <= 10.00) {
              mom += 0.8;
            }
          }
          Pmu.push_back(MaxMom);

          // CosMu binning
          std::vector<double> Cosmu;
          double cosmu = -1.0;
          while (cosmu < 1.0) {
            Cosmu.push_back(cosmu);
            if (cosmu <= 0.2) {
              cosmu += 0.1;
            } else if (cosmu > 0.2 && cosmu <= 0.5) {
              cosmu += 0.05;
            } else if (cosmu > 0.50 && cosmu <= 0.70) {
              cosmu += 0.03;
            } else if (cosmu > 0.70 && cosmu <= 0.95) {
              cosmu += 0.02;
            } else if (cosmu > 0.95 && cosmu <= 1.00) {
              cosmu += 0.01;
            }
          }
          Cosmu.push_back(1.0);

          // Enu binning
          std::vector<double> Enu;
          double enu = 0.0;
          while (enu < 5.0) {
            Enu.push_back(enu);
            if (enu < 0.3) {
              enu += 0.30;
            } else if (enu < 0.9) {
              enu += 0.15;
            } else if (enu < 1.4) {
              enu += 0.3;
            } else if (enu >= 1.4) {
              enu += 0.6;
            }
          }
          Enu.push_back(5.0);

          PmuCosmu[i] = new TH2D(ss.str().c_str(), ss.str().c_str(), Pmu.size()-1, &Pmu[0], Cosmu.size()-1, &Cosmu[0]);
          PmuCosmuEnu[i] = new TH3D(ss2.str().c_str(), ss2.str().c_str(), Pmu.size()-1, &Pmu[0], Cosmu.size()-1, &Cosmu[0], Enu.size()-1, &Enu[0]);
          break;
        }
        // NC1pi01n
      case 31:
      case 32:
      case 33:
      case 34:
        {
          // Make custom arrays, so boring
          std::vector<double> Pmu;
          double mom = 0.0;
          while (mom < MaxMom) {
            Pmu.push_back(mom);
            if (mom <= 0.2) {
              mom += 0.08;
            } else if (mom > 0.2 && mom <= 0.75) {
              mom += 0.05;
            } else if (mom > 0.75 && mom <= 1.50) {
              mom += 0.10;
            } else if (mom > 1.50 && mom <= 4.00) {
              mom += 0.25;
            } else if (mom > 4.00 && mom <= 10.00) {
              mom += 0.8;
            }
          }
          Pmu.push_back(MaxMom);

          std::vector<double> Cosmu;
          double cosmu = -1.0;
          while (cosmu < 1.0) {
            Cosmu.push_back(cosmu);
            if (cosmu <= 0.2) {
              cosmu += 0.1;
            } else if (cosmu > 0.2 && cosmu <= 0.5) {
              cosmu += 0.04;
            } else if (cosmu > 0.50 && cosmu <= 1.00) {
              cosmu += 0.03;
            }
          }
          Cosmu.push_back(1.0);

          // Enu binning
          std::vector<double> Enu;
          double enu = 0.0;
          while (enu < 5.0) {
            Enu.push_back(enu);
            if (enu <= 0.3) {
              enu += 0.3;
            } else if (enu < 0.9) {
              enu += 0.15;
            } else if (enu < 1.4) {
              enu += 0.3;
            } else {
              enu += 0.6;
            }
          }
          Enu.push_back(5.0);

          PmuCosmu[i] = new TH2D(ss.str().c_str(), ss.str().c_str(), Pmu.size()-1, &Pmu[0], Cosmu.size()-1, &Cosmu[0]);
          PmuCosmuEnu[i] = new TH3D(ss2.str().c_str(), ss2.str().c_str(), Pmu.size()-1, &Pmu[0], Cosmu.size()-1, &Cosmu[0], Enu.size()-1, &Enu[0]);
          break;
        }
      default:
        PmuCosmu[i] = NULL;
        PmuCosmuEnu[i] = NULL;
        break;
    }

    if (PmuCosmu[i] != NULL) {
      PmuCosmu[i]->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
      PmuCosmu[i]->GetYaxis()->SetTitle("cos#theta_{#mu}");
      PmuCosmu[i]->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#mu}dcos#theta_{#mu} (cm^{2}/GeV/1/nucleon)");

      PmuCosmuEnu[i]->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
      PmuCosmuEnu[i]->GetYaxis()->SetTitle("cos#theta_{#mu}");
      PmuCosmuEnu[i]->GetZaxis()->SetTitle("E_{#nu} (GeV)");
    }
  }

  int eventCnt = 0;

  std::cout << "Total number of events: " << nevents << std::endl;

  TStopwatch clock;
  clock.Start();

  int countwidth = int(nevents/20);
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

    PmuCosmu[abs(nvect->Mode)]->Fill(plep, costhlep);
    PmuCosmuEnu[abs(nvect->Mode)]->Fill(plep, costhlep, enu);
  }

  std::cout << eventCnt << "/" << nevents << " events" << std::endl;

  // Make filename for selection
  if      (sigType == kCC0pi)  fileName += "_CC0pi";
  else if (sigType == kCC1pip) fileName += "_CC1pip";
  else if (sigType == kCC1pi0) fileName += "_CC1pi0";

  fileName += "_1pi_modes_outgoing.root";

  TFile *output = new TFile(fileName.c_str(), "recreate");
  output->cd();
  for (size_t i = 0; i < nModes; ++i) {

    if (PmuCosmu[i] == NULL) continue;

    if (PmuCosmu[i]->Integral() > 0) {
      PmuCosmu[i]->Sumw2();
      PmuCosmu[i]->Scale(ScaleFactor, "width");
      PmuCosmu[i]->Write();
    }

    if (PmuCosmuEnu[i]->Integral() > 0) {
      PmuCosmuEnu[i]->Sumw2();
      PmuCosmuEnu[i]->Scale(ScaleFactor, "width");
      PmuCosmuEnu[i]->Write();
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



// *********************************************
// Q2 reconstructed for CC1pi+ (difference is 4 vector from pion)
double Q2CCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip) {
  // *********************************************

  double E_mu = pmu.E()/1000.;// energy of lepton in GeV
  double p_mu = pmu.Vect().Mag()/1000.; // momentum of lepton
  double m_mu = sqrt(E_mu*E_mu - p_mu*p_mu); // lepton mass
  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());

  // Just use the true Enu here
  double rEnu = pnu.E()/1000.;
  double q2 = -m_mu*m_mu + 2.*rEnu*(E_mu - p_mu*cos(th_nu_mu)); 

  return q2;
};

// *********************************************
// Q2 true
double Q2true(TLorentzVector pnu, TLorentzVector pmu) {
  // *********************************************
  return -1.0*(pnu-pmu).Mag2()/1.E6;
}

// *********************************************
// Enu reconstructed for CC1pi+
double EnuCCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding){
  // *********************************************

  double E_mu = pmu.E()/1000;
  double p_mu = pmu.Vect().Mag()/1000;
  double m_mu = sqrt(E_mu*E_mu - p_mu*p_mu); 
  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());

  double E_pip = ppip.E()/1000;
  double p_pip = ppip.Vect().Mag()/1000;
  double m_pip = sqrt(E_pip*E_pip - p_pip*p_pip); 
  double th_nu_pip = pnu.Vect().Angle(ppip.Vect());

  const double V  = binding/1000.;    // binding potential 
  const double m_n = 0.93956536;       // neutron mass   
  const double m_n_eff = m_n - V;

  double th_pip_mu = ppip.Vect().Angle(pmu.Vect());

  double rEnu = (m_mu*m_mu + m_pip*m_pip - 2*m_n_eff*(E_pip + E_mu) + 2*E_pip*E_mu - 2*p_pip*p_mu*cos(th_pip_mu))/(2*(E_pip + E_mu - p_pip*cos(th_nu_pip) - p_mu*cos(th_nu_mu) - m_n_eff));

  return rEnu;
};

// *********************************************
// Reconstructed W
double Wrec(TLorentzVector pnu, TLorentzVector pmu) {
  // *********************************************
  double E_mu = pmu.E();
  double p_mu = pmu.Vect().Mag();
  double m_mu = sqrt(E_mu*E_mu - p_mu*p_mu);
  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());

  double E_nu = pnu.E();

  // proton mass (should technically be neutron mass for interactions happening
  // on the neutron! this is easy code but boring and has negligable effect)
  const double m_p = 938.27203;       

  // Enu + Ep = Emu + EHad! 
  // proton is at rest -> Ep = m_p
  // this is dodgy and will have to be checked!
  // double E_nu = E_mu + EHad;

  double q2 = 2*E_nu*(E_mu - p_mu * cos(th_nu_mu)) - m_mu*m_mu; 

  double w_rec = sqrt(m_p*m_p - q2 + 2*m_p*(E_nu - E_mu));

  return w_rec;
};

// *********************************************
// True W, which requires the primary nucleon
double Wtrue(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector pprim) {
  // *********************************************
  return sqrt((pnu - pmu + pprim).Mag2());
}

