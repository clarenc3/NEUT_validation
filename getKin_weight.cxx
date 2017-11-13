#include "getKin_weight.h"

// The main
int main(int argc, char* argv[]) {

  if (argc != 6) {
    Usage();
    exit(-1);
  }

  getKin_weight(std::string(argv[1]), std::atoi(argv[2]), std::atof(argv[3]), std::string(argv[4]), std::string(argv[5]));

  return 0;
};

// Print usage
void Usage() {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << "./getKin.exe ROOT_FILE SIGNAL MAX_MOM CORRECTED_FILE CORRECTED_HIST" << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "SIGNAL is 0 = CC0pi, 1 = CC1pi, 2 = CC1pi0" << std::endl;
  std::cout << "MAX_MOM is the rough momentum scale of the neutrinos" << std::endl;

  return;
}

// The main loop
void getKin_weight(std::string fileName, int sigType, double maxMom, std::string WeightFileName, std::string HistoName) {

  if (sigType == kCC0pi)        std::cout << "Signal is CC0pi" << std::endl;
  else if (sigType == kCC1pip)  std::cout << "Signal is CC1pi+" << std::endl;
  else if (sigType == kCC1pi0)  std::cout << "Signal is CC1pi0" << std::endl;
  else {
    std::cout << "Unrecognised signal! Exiting..." << std::endl;
    Usage();
    exit(-1);
  }

  // First get the weights histogram
  TFile *WeightFile = new TFile(WeightFileName.c_str(), "OPEN");
  TH3D* WeightSpectrum = (TH3D*)(WeightFile->Get(HistoName.c_str()));
  WeightSpectrum->SetDirectory(0);
  WeightSpectrum->SetNameTitle("WeightSpectrum", "WeightSpectrum");
  WeightFile->Close();

  TFile *f = TFile::Open((fileName).c_str(),"open");
  TTree *tn = (TTree*)(f->Get("neuttree"));
  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch",&nvect);

  TH1D *fluxHist = (TH1D*)f->Get("flux_numu");
  f->ls();
  TH1D *eventHist = (TH1D*)f->Get("evtrt_numu");
  long int nevents = tn->GetEntries();

  double ScaleFactor = eventHist->Integral("width")*1E-38/(nevents*fluxHist->Integral("width"));

  int nBins = 50;

  // restricted T2K phase space?
  bool isRestricted = false;
  // Minimum costheta for plots
  double minAng = -1;
  // if we're in restricted costheta plots should be > 0.2
  if (isRestricted) minAng = 0.2;

  // Muon momentum 1D plot
  TH1D *hPmu = new TH1D("hPmu","hPmu", nBins, 0, 1.5);
  hPmu->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPmu->GetYaxis()->SetTitle("d#sigma/dp_{#mu} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPprot = new TH1D("hPprot","hPprot", nBins, 0, 2);
  hPprot->GetXaxis()->SetTitle("p_{p} (GeV/c)");
  hPprot->GetYaxis()->SetTitle("d#sigma/dp_{p} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPneut = new TH1D("hPneut","hPneut", nBins, 0, 2);
  hPneut->GetXaxis()->SetTitle("p_{n} (GeV/c)");
  hPneut->GetYaxis()->SetTitle("d#sigma/dp_{n} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hThmu = new TH1D("hThmu","hThmu", nBins, minAng, 1);
  hThmu->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hThmu->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu} (cm^{2}/nucleon)");

  TH1D *hThprot = new TH1D("hThprot","hThprot", nBins, 0, 1);
  hThprot->GetXaxis()->SetTitle("cos#theta_{p}");
  hThprot->GetYaxis()->SetTitle("d#sigma/dcos#theta_{p} (cm^{2}/nucleon)");

  TH1D *hThneut = new TH1D("hThneut","hThneut", nBins, -1, 1);
  hThneut->GetXaxis()->SetTitle("cos#theta_{n}");
  hThneut->GetYaxis()->SetTitle("d#sigma/dcos#theta_{n} (cm^{2}/nucleon)");

  TH1D *hQ2 = new TH1D("hQ2","hQ2", nBins, 0, 1.5);
  hQ2->GetXaxis()->SetTitle("Q^{2}_{true} (GeV^{2})");
  hQ2->GetYaxis()->SetTitle("d#sigma/dQ_{true}^{2} (cm^{2}/nucleon/GeV^{2})");

  // Make a mode array for Q2
  TH1D *hQ2_mode[52];
  for (int i = 0; i < 52; ++i) {
    std::stringstream ss;
    ss << i;

    hQ2_mode[i] = new TH1D((std::string("hQ2")+ss.str()).c_str(), (std::string("hQ2")+ss.str()).c_str(), nBins, 0, 1.5);
    hQ2_mode[i]->GetXaxis()->SetTitle("Q^{2}_{true} (GeV^{2})");
    hQ2_mode[i]->GetYaxis()->SetTitle("d#sigma/dQ_{true}^{2} (cm^{2}/nucleon/GeV^{2})");
    hQ2_mode[i]->SetTitle((ss.str()).c_str());
    hQ2_mode[i]->SetFillColor(i);
    hQ2_mode[i]->SetLineColor(i);
    hQ2_mode[i]->SetFillStyle(1001);
  }

  TH1D *hW = new TH1D("hW","hW", nBins, 1.1, 1.6);
  hW->GetXaxis()->SetTitle("W_{true} (GeV/c^{2})");
  hW->GetYaxis()->SetTitle("d#sigma/dW_{true} (cm^{2}/nucleon/(GeV/c^{2}))");

  TH1D *hW_rec = new TH1D("hW_rec","hW_rec", nBins, 1.1, 1.6);
  hW_rec->GetXaxis()->SetTitle("W_{rec} (GeV/c^{2})");
  hW_rec->GetYaxis()->SetTitle("d#sigma/dW_{rec} (cm^{2}/nucleon/(GeV/c^{2}))");

  TH2D *hWQ2 = new TH2D("hWQ2", "hWQ2", nBins, 1.1, 1.6, nBins, 0, 1.5);
  hWQ2->GetXaxis()->SetTitle("W (GeV/c^{2})");
  hWQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2}");
  hWQ2->GetZaxis()->SetTitle("d^{2}#sigma/dWdQ^{2} (cm^{2}/nucleon/(GeV/c^{2})/GeV^{2})");

  TH1D *hEnu = new TH1D("hEnu","hEnu",nBins, 0.4, 2.0);
  hEnu->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hEnu->GetYaxis()->SetTitle("#sigma (E_{#nu}) (flux averaged!)");

  TH2D *hPthetaMu = new TH2D("hPthetaMu", "hPthetaMu", nBins, 0, maxMom, nBins, minAng, 1);
  hPthetaMu->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPthetaMu->GetYaxis()->SetTitle("cos#theta_{#mu}");
  hPthetaMu->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#mu}dcos#theta_{#mu} (cm^{2}/nucleon/(GeV/c))");

  TH2D *hPmuQ2 = new TH2D("hPmuQ2", "hPmuQ2", nBins, 0, maxMom, nBins, 0, 1.5);
  hPmuQ2->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPmuQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  hPmuQ2->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#mu}dQ^{2} (cm^{2}/nucleon/(GeV/c)/GeV^{2})");

  TH2D *hThetaQ2 = new TH2D("hThetaQ2", "hThetaQ2", nBins, minAng, 1, nBins, 0, 1.5);
  hThetaQ2->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hThetaQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  hThetaQ2->GetZaxis()->SetTitle("d^{2}#sigma/dQ^{2}dcos#theta_{#mu} (cm^{2}/nucleon/(GeV/c)/GeV^{2})");

  TH2D *hPthetaProt = new TH2D("hPthetaProt", "hPthetaProt", nBins, 0, 2.0, nBins, -1, 1);
  hPthetaProt->GetXaxis()->SetTitle("p_{p} (GeV/c)");
  hPthetaProt->GetYaxis()->SetTitle("cos#theta_{p}");
  hPthetaProt->GetZaxis()->SetTitle("d^{2}#sigma/dp_{p}dcos#theta_{p} (cm^{2}/nucleon/(GeV/c))");

  TH2D *hPthetaNeut = new TH2D("hPthetaNeut", "hPthetaNeut", nBins, 0, 2.0, nBins, -1, 1);
  hPthetaNeut->GetXaxis()->SetTitle("p_{n} (GeV/c)");
  hPthetaNeut->GetYaxis()->SetTitle("cos#theta_{n}");
  hPthetaNeut->GetZaxis()->SetTitle("d^{2}#sigma/dp_{n}dcos#theta_{n} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPpi = new TH1D("hPpi","hPpi", nBins, 0, 1);
  hPpi->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
  hPpi->GetYaxis()->SetTitle("d#sigma/dp_{#pi} (cm^{2}/nucleon/(GeV/c))");
  // Make a mode array for Q2
  TH1D *hPpi_mode[52];
  for (int i = 0; i < 52; ++i) {
    std::stringstream ss;
    ss << i;

    hPpi_mode[i] = new TH1D((std::string("hPpi")+ss.str()).c_str(), (std::string("hPpi")+ss.str()).c_str(), nBins, 0, 1);
    hPpi_mode[i]->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
    hPpi_mode[i]->GetYaxis()->SetTitle("d#sigma/dp_{#pi} (cm^{2}/nucleon/(GeV/c))");
    hPpi_mode[i]->SetTitle((ss.str()).c_str());
    hPpi_mode[i]->SetFillColor(i);
    hPpi_mode[i]->SetLineColor(i);
    hPpi_mode[i]->SetFillStyle(1001);
  }

  TH1D *hThpi = new TH1D("hThpi","hThpi", nBins, minAng, 1);
  hThpi->GetXaxis()->SetTitle("cos#theta_{#pi}");
  hThpi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi} (cm^{2}/nucleon)");

  TH1D *hThprotpi = new TH1D("hThprotpi","hThprotpi", nBins, -1, 1);
  hThprotpi->GetXaxis()->SetTitle("cos#theta_{p,#pi}");
  hThprotpi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{p#pi} (cm^{2}/nucleon)");

  TH1D *hThmupi = new TH1D("hThmupi","hThmupi", nBins, -1, 1);
  hThmupi->GetXaxis()->SetTitle("cos#theta_{#mu,#pi}");
  hThmupi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu,#pi} (cm^{2}/nucleon)");

  TH2D *hPthetaPi = new TH2D("hPthetaPi", "hPthetaPi", nBins, 0, 1.0, nBins, minAng, 1);
  hPthetaPi->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
  hPthetaPi->GetYaxis()->SetTitle("cos#theta_{#pi}");
  hPthetaPi->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#pi}dcos#theta_{#pi} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hNp = new TH1D("hNp","hNp", 6, 0, 6);
  hNp->GetXaxis()->SetTitle("N_{p}");
  hNp->GetYaxis()->SetTitle("d#sigma/dN_{p} (cm^{2}/nucleon)");

  TH1D *hNn = new TH1D("hNn","hNn", 6, 0, 6);
  hNn->GetXaxis()->SetTitle("N_{n}");
  hNn->GetYaxis()->SetTitle("d#sigma/dN_{n} (cm^{2}/nucleon)");

  // THETA PION IN RESONANCE REST FRAME
  TH1D *hCosThPiRest = new TH1D("hCosThPiRest","hCosThPiRest", 50, -1, 1);
  hCosThPiRest->GetXaxis()->SetTitle("cos#theta_{#pi}^{RES}");
  hCosThPiRest->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi}^{RES} (cm^{2}/nucleon/1)");

  // THETA PION IN RESONANCE REST FRAME
  TH1D *hPhiRest = new TH1D("hPhiRest","hPhiRest", 50, 0, M_PI);
  hPhiRest->GetXaxis()->SetTitle("#phi_{#pi}^{RES}");
  hPhiRest->GetYaxis()->SetTitle("d#sigma/d#phi_{#pi}^{RES} (cm^{2}/nucleon/1)");

  int eventCnt = 0;

  std::cout << "number of events: " << nevents << std::endl;

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
    } else if (sigType == kCC1pip) {
      // Only plot the pion modes
      if (nvect->Mode != 11 && nvect->Mode != 12 && nvect->Mode != 13) continue;
      eventCnt++;
    } else if (sigType == kCC1pi0) {
      if (!isT2K_CC1pi0(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
    }

    TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
    TLorentzVector Pp_init = (nvect->PartInfo(1))->fP;
    TLorentzVector Pp;
    TLorentzVector Pn;
    TLorentzVector Ppip;
    TLorentzVector Pmu;
    TLorentzVector PnucPrimary;
    TLorentzVector PInitialState;
    TLorentzVector PpipPrimary;

    int np = 0;
    int nn = 0;

    // Loop over the particle stack
    for (int k = 2; k < nvect->Npart(); ++k) {
      
      // Get the PID of the particle
      int PID = (nvect->PartInfo(k))->fPID;

      // Pick out the pre-FSI particles for CC1pi+ or CC1pi0
      // Need these for Adler angles and W_true
      if (sigType == kCC1pip || sigType == kCC1pi0) {
        if (k < 5) {
          if ((PID == 2212) || (PID == 2112)) {
            PnucPrimary = nvect->PartInfo(k)->fP;
          } else if (PID == 211) {
            PpipPrimary = nvect->PartInfo(k)->fP;
          }
        }
      }

      if (!(nvect->PartInfo(k))->fIsAlive && (nvect->PartInfo(k))->fStatus != 0) continue;

      // Pick out the pion
      if ((PID == 211 && sigType == kCC1pip) || (PID == 111 && sigType == kCC1pi0)) {
        Ppip = nvect->PartInfo(k)->fP;

        // Pick out the highest momentum proton
      } else if (PID == 2212) {
        np++;
        if (nvect->PartInfo(k)->fP.Vect().Mag() > Pp.Vect().Mag()) {
          Pp = nvect->PartInfo(k)->fP;
        }

        // Pick out the highest momentum neutron
      } else if (PID == 2112) {
        nn++;
        if (nvect->PartInfo(k)->fP.Vect().Mag() > Pn.Vect().Mag()) {
          Pn = nvect->PartInfo(k)->fP;
        }

        // Pick out the muon
      } else if (PID == 13) {
        Pmu = (nvect->PartInfo(k))->fP;  
      }

    } // Finish the for loop over particles

    // Also get the initial state
    PInitialState = nvect->PartInfo(1)->fP;

    // Get the muon variables
    double pmu = Pmu.Vect().Mag()/1000.;
    double thmu = cos(Pnu.Vect().Angle(Pmu.Vect()));

    // Get the proton variables
    double pprot = Pp.Vect().Mag()/1000.;
    double thprot = cos(Pnu.Vect().Angle(Pp.Vect()));

    // Get the neutron variabls
    double pneut = Pn.Vect().Mag()/1000.;
    double thneut = cos(Pnu.Vect().Angle(Pn.Vect()));

    // Derived kinematic variables
    double Enu = Pnu.E()/1000.;
    double Q2 = Q2true(Pnu, Pmu);
    double W;
    double W_rec;

    // Pion kinematics
    double ppi;
    double thpi;
    double thprotpi;
    double thmupi;

    // Adler angles
    double costhAdler;
    double phiAdler;

    double weight = 1.;
    // If we have a pion
    if (sigType == kCC1pi0 || sigType == kCC1pip) {

      ppi = Ppip.Vect().Mag()/1000.;
      thpi = cos(Pnu.Vect().Angle(Ppip.Vect()));
      thprotpi = cos(Ppip.Vect().Angle(Pp.Vect()));
      thmupi = cos(Pmu.Vect().Angle(Ppip.Vect()));

      TLorentzVector PResRest = PnucPrimary + PpipPrimary;
      Ppip.Boost(-PResRest.BoostVector());
      PnucPrimary.Boost(-PResRest.BoostVector());
      costhAdler = cos(Ppip.Angle(PResRest.Vect()));

      W_rec = Wrec(Pnu, Pmu)/1000.;
      W = Wtrue(Pnu, Pmu, PInitialState)/1000.;

      weight = WeightSpectrum->GetBinContent(WeightSpectrum->FindBin(W, Q2, Enu));
    }

    // Now that we know it is signal and have the relevant kinematics, find the weight for this event
    //std::cout << weight << std::endl;


    if (pmu > 0) {
      hPmu->Fill(pmu, weight);
      hThmu->Fill(thmu, weight);
      hPthetaMu->Fill(pmu, thmu, weight);

      hPmuQ2->Fill(pmu, Q2, weight);
      hThetaQ2->Fill(thmu, Q2, weight);

      hEnu->Fill(Enu, weight);
      hQ2->Fill(Q2, weight);
      hQ2_mode[nvect->Mode]->Fill(Q2, weight);
    }

    if (pprot > 0) {
      hPprot->Fill(pprot, weight);
      hThprot->Fill(thprot, weight);
      hPthetaProt->Fill(pprot, thprot, weight);
    }

    if (pneut > 0) {
      hPneut->Fill(pneut, weight);
      hThneut->Fill(thneut, weight);
      hPthetaNeut->Fill(pneut, thneut, weight);
    }

    if (sigType == kCC1pip || sigType == kCC1pi0) {

      if (ppi > 0) {

        hPpi->Fill(ppi,weight);
        hPpi_mode[nvect->Mode]->Fill(ppi,weight);
        hThpi->Fill(thpi,weight);

        hPthetaPi->Fill(ppi, thpi, weight);
        if (pprot > 0) {
          hThprotpi->Fill(thprotpi,weight);
          hThmupi->Fill(thmupi,weight);

          hCosThPiRest->Fill(costhAdler,weight);
        }
      }

      if (pmu > 0) {
        hW->Fill(W,weight);
        hWQ2->Fill(W, Q2,weight);
        hW_rec->Fill(W_rec,weight);
      }
    }

    // Always fill these
    hNp->Fill(np,weight);
    hNn->Fill(nn,weight);

  }

  std::cout << eventCnt << "/" << nevents << " events" << std::endl;

  hPmu->Sumw2();
  hPmu->Scale(ScaleFactor,"width");
  hThmu->Sumw2();
  hThmu->Scale(ScaleFactor, "width");

  hPprot->Sumw2();
  hPprot->Scale(ScaleFactor, "width");
  hThprot->Sumw2();
  hThprot->Scale(ScaleFactor, "width");

  hPneut->Sumw2();
  hPneut->Scale(ScaleFactor, "width");
  hThneut->Sumw2();
  hThneut->Scale(ScaleFactor,"width");

  hPthetaMu->Sumw2();
  hPthetaMu->Scale(ScaleFactor, "width");
  hPthetaProt->Sumw2();
  hPthetaProt->Scale(ScaleFactor, "width");
  hPthetaNeut->Sumw2();
  hPthetaNeut->Scale(ScaleFactor, "width");

  hPmuQ2->Sumw2();
  hPmuQ2->Scale(ScaleFactor, "width");
  hThetaQ2->Sumw2();
  hThetaQ2->Scale(ScaleFactor, "width");

  hQ2->Sumw2();
  hQ2->Scale(ScaleFactor, "width");
  for (int i = 0; i < 52; ++i) {
    for (int j = 0; j < hQ2_mode[i]->GetNbinsX()+1; ++j) {
      hQ2_mode[i]->SetBinError(j+1, 0);
    }

    hQ2_mode[i]->Scale(ScaleFactor, "width");
  }

  hEnu->Sumw2();
  hEnu->Scale(ScaleFactor, "width");

  if (sigType == kCC1pi0 || sigType == kCC1pip) {
    hPpi->Sumw2();
    hPpi->Scale(ScaleFactor, "width");

    for (int i = 0; i < 52; ++i) {
      for (int j = 0; j < hPpi_mode[i]->GetNbinsX()+1; ++j) {
        hPpi_mode[i]->SetBinError(j+1, 0);
      }

      hPpi_mode[i]->Scale(ScaleFactor, "width");
    }

    hThpi->Sumw2();
    hThpi->Scale(ScaleFactor, "width");

    hThprotpi->Sumw2();
    hThprotpi->Scale(ScaleFactor, "width");

    hThmupi->Sumw2();
    hThmupi->Scale(ScaleFactor, "width");

    hCosThPiRest->Sumw2();
    hCosThPiRest->Scale(ScaleFactor, "width");

    hPthetaPi->Sumw2();
    hPthetaPi->Scale(ScaleFactor, "width");

    hWQ2->Sumw2();
    hWQ2->Scale(ScaleFactor, "width");

    hW->Sumw2();
    hW->Scale(ScaleFactor, "width");

    hW_rec->Sumw2();
    hW_rec->Scale(ScaleFactor, "width");
  }

  hNp->Sumw2();
  hNp->Scale(ScaleFactor, "width");

  hNn->Sumw2();
  hNn->Scale(ScaleFactor, "width");

  // Make filename for selection
  if      (sigType == kCC0pi)  fileName += "_CC0pi";
  else if (sigType == kCC1pip) fileName += "_CC1pip";
  else if (sigType == kCC1pi0) fileName += "_CC1pi0";

  // Make filename for not restricted
  if (!isRestricted) fileName += "_noKinCuts";

  fileName+="_outgoing_cc1pi_corr"+HistoName+".root";

  TFile *output = new TFile(fileName.c_str(), "recreate");

  output->cd();

  hPmu->Write();
  hThmu->Write();

  hPprot->Write();
  hThprot->Write();

  hPneut->Write();
  hThneut->Write();

  hPthetaMu->Write();
  hPthetaProt->Write();
  hPthetaNeut->Write();

  hPmuQ2->Write();
  hThetaQ2->Write();

  hEnu->Write();
  hQ2->Write();

  std::vector<int> nModes;
  // Count how many contributing modes we have
  for (int i = 0; i < 52; ++i) {
    if (hQ2_mode[i]->Integral() > 0) {
      nModes.push_back(i);
    }
  }
  // Now make the THStack in Q2

  THStack *q2_stack = new THStack("hQ2_stack", "hQ2_stack");
  unsigned int nModes_size = nModes.size();
  for (unsigned int i = 0; i < nModes_size; ++i) {
    hQ2_mode[nModes.at(i)]->SetFillColor(i+1);
    q2_stack->Add(hQ2_mode[nModes.at(i)]);
  }

  q2_stack->Write();


  if (sigType == kCC1pip || sigType == kCC1pi0) {
    hPpi->Write();
    THStack *ppi_stack = new THStack("hPpi_stack","hPpi_stack");
    for (unsigned int i = 0; i < nModes_size; ++i) {
      hPpi_mode[nModes.at(i)]->SetFillColor(i+1);
      ppi_stack->Add(hPpi_mode[nModes.at(i)]);
    }
    ppi_stack->Write();
    hThpi->Write();
    hThprotpi->Write();
    hThmupi->Write();
    hCosThPiRest->Write();

    hPthetaPi->Write();

    hW->Write();
    hW_rec->Write();
    hWQ2->Write();
  }

  hNp->Write();
  hNn->Write();

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

