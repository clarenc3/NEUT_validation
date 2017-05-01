#include "getKin.h"

int main(int argc, char* argv[]) {

  if (argc != 4) {
    Usage();
    exit(-1);
  }

  getKin(std::string(argv[1]), std::atoi(argv[2]), std::atof(argv[3]));

  return 0;
};

void Usage() {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << "./getKin.exe ROOT_FILE SIGNAL MAX_MOM" << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "SIGNAL is 0 = CC0pi, 1 = CC1pi, 2 = CC1pi0" << std::endl;
  std::cout << "MAX_MOM is the rough momentum scale of the neutrinos" << std::endl;

  return;
}

void getKin(std::string fileName, int sigType, double maxMom) {

  if (sigType == 0)       std::cout << "Signal is CC0pi" << std::endl;
  else if (sigType == 1)  std::cout << "Signal is CC1pi+" << std::endl;
  else if (sigType == 2)  std::cout << "Signal is CC1pi0" << std::endl;
  else {
    std::cout << "Unrecognised signal! Exiting..." << std::endl;
    exit(-1);
  }

  TFile *f = TFile::Open((fileName).c_str(),"open");
  TTree *tn = (TTree*)(f->Get("neuttree"));
  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch",&nvect);

  TH1D *fluxHist = (TH1D*)f->Get("flux_numu");
  TH1D *eventHist = (TH1D*)f->Get("evtrt_numu");
  long int nevents = tn->GetEntries();

  double scaleFactor = eventHist->Integral("width")*1E-38/(nevents*fluxHist->Integral("width"));

  int nBins = 50;

  // restricted T2K phase space?
  bool isRestricted = false;
  // Minimum costheta for plots
  double minAng = -1;
  // if we're in restricted costheta plots should be > 0.2
  if (isRestricted) minAng = 0.2;

  TH1D *hPmu = new TH1D("hPmu","hPmu", nBins, 0, maxMom);
  hPmu->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPmu->GetYaxis()->SetTitle("d#sigma/dp_{#mu} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPprot = new TH1D("hPprot","hPprot", nBins, 0, maxMom);
  hPprot->GetXaxis()->SetTitle("p_{p} (GeV/c)");
  hPprot->GetYaxis()->SetTitle("d#sigma/dp_{p} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPneut = new TH1D("hPneut","hPneut", nBins, 0, maxMom);
  hPneut->GetXaxis()->SetTitle("p_{n} (GeV/c)");
  hPneut->GetYaxis()->SetTitle("d#sigma/dp_{n} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hThmu = new TH1D("hThmu","hThmu", nBins, minAng, 1);
  hThmu->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hThmu->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu} (cm^{2}/nucleon)");

  TH1D *hThprot = new TH1D("hThprot","hThprot", nBins, -1, 1);
  hThprot->GetXaxis()->SetTitle("cos#theta_{p}");
  hThprot->GetYaxis()->SetTitle("d#sigma/dcos#theta_{p} (cm^{2}/nucleon)");

  TH1D *hThneut = new TH1D("hThneut","hThneut", nBins, -1, 1);
  hThneut->GetXaxis()->SetTitle("cos#theta_{n}");
  hThneut->GetYaxis()->SetTitle("d#sigma/dcos#theta_{n} (cm^{2}/nucleon)");

  TH1D *hQ2 = new TH1D("hQ2","hQ2", nBins, 0, 3);
  hQ2->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
  hQ2->GetYaxis()->SetTitle("d#sigma/dQ^{2} (cm^{2}/nucleon/GeV^{2})");

  TH1D *hW = new TH1D("hW","hW", nBins, 1, 2);
  hW->GetXaxis()->SetTitle("W (GeV/c^{2})");
  hW->GetYaxis()->SetTitle("d#sigma/dW (cm^{2}/nucleon/(GeV/c^{2}))");

  TH2D *hWQ2 = new TH2D("hWQ2", "hWQ2", nBins, 1, 2, nBins, 0, 3);
  hWQ2->GetXaxis()->SetTitle("W (GeV/c^{2})");
  hWQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2}");
  hWQ2->GetZaxis()->SetTitle("d^{2}#sigma/dWdQ^{2} (cm^{2}/nucleon/(GeV/c^{2})/GeV^{2})");

  TH1D *hEnu = new TH1D("hEnu","hEnu",nBins, 0, maxMom);
  hEnu->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hEnu->GetYaxis()->SetTitle("#sigma (E_{#nu}) (flux averaged!)");

  TH2D *hPthetaMu = new TH2D("hPthetaMu", "hPthetaMu", nBins, 0, maxMom, nBins, minAng, 1);
  hPthetaMu->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPthetaMu->GetYaxis()->SetTitle("cos#theta_{#mu}");
  hPthetaMu->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#mu}dcos#theta_{#mu} (cm^{2}/nucleon/(GeV/c))");

  TH2D *hPmuQ2 = new TH2D("hPmuQ2", "hPmuQ2", nBins, 0, maxMom, nBins, 0, 3);
  hPmuQ2->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPmuQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  hPmuQ2->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#mu}dQ^{2} (cm^{2}/nucleon/(GeV/c)/GeV^{2})");

  TH2D *hThetaQ2 = new TH2D("hThetaQ2", "hThetaQ2", nBins, minAng, 1, nBins, 0, 3);
  hThetaQ2->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hThetaQ2->GetYaxis()->SetTitle("Q^{2} (GeV^{2})");
  hThetaQ2->GetZaxis()->SetTitle("d^{2}#sigma/dQ^{2}dcos#theta_{#mu} (cm^{2}/nucleon/(GeV/c)/GeV^{2})");

  TH2D *hPthetaProt = new TH2D("hPthetaProt", "hPthetaProt", nBins, 0, maxMom, nBins, -1, 1);
  hPthetaProt->GetXaxis()->SetTitle("p_{p} (GeV/c)");
  hPthetaProt->GetYaxis()->SetTitle("cos#theta_{p}");
  hPthetaProt->GetZaxis()->SetTitle("d^{2}#sigma/dp_{p}dcos#theta_{p} (cm^{2}/nucleon/(GeV/c))");

  TH2D *hPthetaNeut = new TH2D("hPthetaNeut", "hPthetaNeut", nBins, 0, maxMom, nBins, -1, 1);
  hPthetaNeut->GetXaxis()->SetTitle("p_{n} (GeV/c)");
  hPthetaNeut->GetYaxis()->SetTitle("cos#theta_{n}");
  hPthetaNeut->GetZaxis()->SetTitle("d^{2}#sigma/dp_{n}dcos#theta_{n} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPpi = new TH1D("hPpi","hPpi", nBins, 0, maxMom);
  hPpi->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
  hPpi->GetYaxis()->SetTitle("d#sigma/dp_{#pi} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hThpi = new TH1D("hThpi","hThpi", nBins, minAng, 1);
  hThpi->GetXaxis()->SetTitle("cos#theta_{#pi}");
  hThpi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi} (cm^{2}/nucleon)");

  TH1D *hThprotpi = new TH1D("hThprotpi","hThprotpi", nBins, -1, 1);
  hThprotpi->GetXaxis()->SetTitle("cos#theta_{p,#pi}");
  hThprotpi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{p#pi} (cm^{2}/nucleon)");

  TH1D *hThmupi = new TH1D("hThmupi","hThmupi", nBins, -1, 1);
  hThmupi->GetXaxis()->SetTitle("cos#theta_{#mu,#pi}");
  hThmupi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu,#pi} (cm^{2}/nucleon)");

  TH2D *hPthetaPi = new TH2D("hPthetaPi", "hPthetaPi", nBins, 0, maxMom, nBins, minAng, 1);
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

  for (int j = 0; j < nevents; j++) {

    tn->GetEntry(j);

    if (j%50000 == 0) std::cout << "On event #" << j << std::endl;

    // Very simple signal definition
    if (sigType == 0) {
      if (!isT2KCC0pi(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
    } else if (sigType == 1) {
      if (!isT2KCCpip(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
    } else if (sigType == 2) {
      if (!isT2KCCpi0(nvect, 0, maxMom*3, isRestricted)) continue;
      eventCnt++;
    }

    /*
    if (nvect->Mode != 11) {
      continue;
    } else {
      eventCnt++;
    }
    */

    TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
    TLorentzVector Pp_init = (nvect->PartInfo(1))->fP;
    TLorentzVector Pp;
    TLorentzVector Pn;
    TLorentzVector Ppip;
    TLorentzVector Pmu;
    TLorentzVector PnucPrimary;
    TLorentzVector PpipPrimary;

    int np = 0;
    int nn = 0;

    // Loop over the particle stack
    for (int k = 2; k < nvect->Npart(); ++k) {
      int PID = (nvect->PartInfo(k))->fPID;
      if (k < 5) {
        if ((PID == 2212) || (PID == 2112)) {
          PnucPrimary = nvect->PartInfo(k)->fP;
        } else if (PID == 211) {
          PpipPrimary = nvect->PartInfo(k)->fP;
        }
      }
      
      if (!(nvect->PartInfo(k))->fIsAlive && (nvect->PartInfo(k))->fStatus != 0) continue;

      if ((PID == 211 && sigType == 1) || (PID == 111 && sigType == 2)) {
        Ppip = nvect->PartInfo(k)->fP;
      // select highest momentum proton
      } else if (PID == 2212) {
        np++;
        if (nvect->PartInfo(k)->fP.Vect().Mag() > Pp.Vect().Mag())
          Pp = nvect->PartInfo(k)->fP;
      } else if (PID == 13) {
        Pmu = (nvect->PartInfo(k))->fP;  
      } else if (PID == 2112) {
        nn++;
        if (nvect->PartInfo(k)->fP.Vect().Mag() > Pn.Vect().Mag())
          Pn = nvect->PartInfo(k)->fP;
      }
    }

    double pmu = Pmu.Vect().Mag()/1000.;
    double thmu = cos(Pnu.Vect().Angle(Pmu.Vect()));

    double pprot = Pp.Vect().Mag()/1000.;
    double thprot = cos(Pnu.Vect().Angle(Pp.Vect()));
    
    double pneut = Pn.Vect().Mag()/1000.;
    double thneut = cos(Pnu.Vect().Angle(Pn.Vect()));

    // Derived kinematic variables
    double Enu;
    double Q2;
    double W;

    // Pion kinematics
    double ppi;
    double thpi;
    double thprotpi;
    double thmupi;

    // Adler angles
    double costhAdler;
    double phiAdler;

    if (sigType == 1 || sigType == 2) {
      ppi = Ppip.Vect().Mag()/1000.;
      thpi = cos(Pnu.Vect().Angle(Ppip.Vect()));
      thprotpi = cos(Ppip.Vect().Angle(Pp.Vect()));
      thmupi = cos(Pmu.Vect().Angle(Ppip.Vect()));

      TLorentzVector PResRest = PnucPrimary + PpipPrimary;
      Ppip.Boost(-PResRest.BoostVector());
      PnucPrimary.Boost(-PResRest.BoostVector());
      costhAdler = cos(Ppip.Angle(PResRest.Vect()));

      Enu = EnuCCpiprec(Pnu, Pmu, Ppip);
      Q2 = Q2CCpiprec(Pnu, Pmu, Ppip);

      W = Wrec(Pnu, Pmu)/1000.;
    }
    else {
      Enu = EnuQErec(Pnu, Pmu);
      Q2 = Q2QErec(Pnu, Pmu);
    }

    hPmu->Fill(pmu);
    hThmu->Fill(thmu);

    hPprot->Fill(pprot);
    hThprot->Fill(thprot);

    hPneut->Fill(pneut);
    hThneut->Fill(thneut);

    hPthetaMu->Fill(pmu, thmu);
    hPthetaProt->Fill(pprot, thprot);
    hPthetaNeut->Fill(pneut, thneut);

    hPmuQ2->Fill(pmu, Q2);
    hThetaQ2->Fill(thmu, Q2);

    hEnu->Fill(Enu);
    hQ2->Fill(Q2);

    if (sigType == 1 || sigType == 2) {
      hPpi->Fill(ppi);
      hThpi->Fill(thpi);
      hThprotpi->Fill(thprotpi);
      hThmupi->Fill(thmupi);

      hCosThPiRest->Fill(costhAdler);

      hPthetaPi->Fill(ppi, thpi);

      hW->Fill(W);

      hWQ2->Fill(W, Q2);
    }

    hNp->Fill(np);
    hNn->Fill(nn);

  }

  // Stop clock
  clock.Stop();
  std::cout << clock.RealTime() << " seconds " << std::endl;

  std::cout << eventCnt << "/" << nevents << " events" << std::endl;

  for (int i = 0; i < hPmu->GetNbinsX()+1; i++) {
    double width = hPmu->GetBinLowEdge(i+2) - hPmu->GetBinLowEdge(i+1);
    if (hPmu->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hPmu->GetBinContent(i+1));
    hPmu->SetBinContent(i+1, hPmu->GetBinContent(i+1)*scaleFactor/width);
    hPmu->SetBinError(i+1, fracErr*hPmu->GetBinContent(i+1));
  }

  for (int i = 0; i < hThmu->GetNbinsX()+1; i++) {
    double width = hThmu->GetBinLowEdge(i+2) - hThmu->GetBinLowEdge(i+1);
    if (hThmu->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hThmu->GetBinContent(i+1));
    hThmu->SetBinContent(i+1, hThmu->GetBinContent(i+1)*scaleFactor/width);
    hThmu->SetBinError(i+1, fracErr*hThmu->GetBinContent(i+1));
  }

  for (int i = 0; i < hPprot->GetNbinsX()+1; i++) {
    double width = hPprot->GetBinLowEdge(i+2) - hPprot->GetBinLowEdge(i+1);
    if (hPprot->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hPprot->GetBinContent(i+1));
    hPprot->SetBinContent(i+1, hPprot->GetBinContent(i+1)*scaleFactor/width);
    hPprot->SetBinError(i+1, fracErr*hPprot->GetBinContent(i+1));
  }

  for (int i = 0; i < hThprot->GetNbinsX()+1; i++) {
    double width = hThprot->GetBinLowEdge(i+2) - hThprot->GetBinLowEdge(i+1);
    if (hThprot->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hThprot->GetBinContent(i+1));
    hThprot->SetBinContent(i+1, hThprot->GetBinContent(i+1)*scaleFactor/width);
    hThprot->SetBinError(i+1, fracErr*hThprot->GetBinContent(i+1));
  }

  for (int i = 0; i < hPneut->GetNbinsX()+1; i++) {
    double width = hPneut->GetBinLowEdge(i+2) - hPneut->GetBinLowEdge(i+1);
    if (hPneut->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hPneut->GetBinContent(i+1));
    hPneut->SetBinContent(i+1, hPneut->GetBinContent(i+1)*scaleFactor/width);
    hPneut->SetBinError(i+1, fracErr*hPneut->GetBinContent(i+1));
  }

  for (int i = 0; i < hThneut->GetNbinsX()+1; i++) {
    double width = hThneut->GetBinLowEdge(i+2) - hThneut->GetBinLowEdge(i+1);
    if (hThneut->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hThneut->GetBinContent(i+1));
    hThneut->SetBinContent(i+1, hThneut->GetBinContent(i+1)*scaleFactor/width);
    hThneut->SetBinError(i+1, fracErr*hThneut->GetBinContent(i+1));
  }

  //hPmu->Scale(scaleFactor, "width");
  //hThmu->Scale(scaleFactor, "width");

  //hPprot->Scale(scaleFactor, "width");
  //hThprot->Scale(scaleFactor, "width");

  hPthetaMu->Scale(scaleFactor, "width");
  hPthetaProt->Scale(scaleFactor, "width");
  hPthetaNeut->Scale(scaleFactor, "width");

  hPmuQ2->Scale(scaleFactor, "width");
  hThetaQ2->Scale(scaleFactor, "width");
  //hQ2->Scale(scaleFactor, "width");

  for (int i = 0; i < hQ2->GetNbinsX()+1; i++) {
    double width = hQ2->GetBinLowEdge(i+2) - hQ2->GetBinLowEdge(i+1);
    if (hQ2->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hQ2->GetBinContent(i+1));
    hQ2->SetBinContent(i+1, hQ2->GetBinContent(i+1)*scaleFactor/width);
    hQ2->SetBinError(i+1, fracErr*hQ2->GetBinContent(i+1));
  }

  for (int i = 0; i < hEnu->GetNbinsX()+1; i++) {
    double width = hEnu->GetBinLowEdge(i+2) - hEnu->GetBinLowEdge(i+1);
    if (hEnu->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hEnu->GetBinContent(i+1));
    hEnu->SetBinContent(i+1, hEnu->GetBinContent(i+1)*scaleFactor/*/width*/);
    hEnu->SetBinError(i+1, fracErr*hEnu->GetBinContent(i+1));
  }

  
  if (sigType == 1 || sigType == 2) {
    /*
    hPpi->Scale(scaleFactor, "width");
    hThpi->Scale(scaleFactor, "width");
    hThprotpi->Scale(scaleFactor, "width");
    hThmupi->Scale(scaleFactor, "width");
    */
    for (int i = 0; i < hPpi->GetNbinsX()+1; i++) {
      double width = hPpi->GetBinLowEdge(i+2) - hPpi->GetBinLowEdge(i+1);
      if (hPpi->GetBinContent(i+1) == 0 || width == 0) continue;
      double fracErr = sqrt(1/hPpi->GetBinContent(i+1));
      hPpi->SetBinContent(i+1, hPpi->GetBinContent(i+1)*scaleFactor/width);
      hPpi->SetBinError(i+1, fracErr*hPpi->GetBinContent(i+1));
    }

    for (int i = 0; i < hThpi->GetNbinsX()+1; i++) {
      double width = hThpi->GetBinLowEdge(i+2) - hThpi->GetBinLowEdge(i+1);
      if (hThpi->GetBinContent(i+1) == 0 || width == 0) continue;
      double fracErr = sqrt(1/hThpi->GetBinContent(i+1));
      hThpi->SetBinContent(i+1, hThpi->GetBinContent(i+1)*scaleFactor/width);
      hThpi->SetBinError(i+1, fracErr*hThpi->GetBinContent(i+1));
    }

    for (int i = 0; i < hThprotpi->GetNbinsX()+1; i++) {
      double width = hThprotpi->GetBinLowEdge(i+2) - hThprotpi->GetBinLowEdge(i+1);
      if (hThprotpi->GetBinContent(i+1) == 0 || width == 0) continue;
      double fracErr = sqrt(1/hThprotpi->GetBinContent(i+1));
      hThprotpi->SetBinContent(i+1, hThprotpi->GetBinContent(i+1)*scaleFactor/width);
      hThprotpi->SetBinError(i+1, fracErr*hThprotpi->GetBinContent(i+1));
    }

    for (int i = 0; i < hThmupi->GetNbinsX()+1; i++) {
      double width = hThmupi->GetBinLowEdge(i+2) - hThmupi->GetBinLowEdge(i+1);
      if (hThmupi->GetBinContent(i+1) == 0 || width == 0) continue;
      double fracErr = sqrt(1/hThmupi->GetBinContent(i+1));
      hThmupi->SetBinContent(i+1, hThmupi->GetBinContent(i+1)*scaleFactor/width);
      hThmupi->SetBinError(i+1, fracErr*hThmupi->GetBinContent(i+1));
    }

    for (int i = 0; i < hCosThPiRest->GetNbinsX()+1; i++) {
      double width = hCosThPiRest->GetBinLowEdge(i+2) - hCosThPiRest->GetBinLowEdge(i+1);
      if (hCosThPiRest->GetBinContent(i+1) == 0 || width == 0) continue;
      double fracErr = sqrt(1/hCosThPiRest->GetBinContent(i+1));
      hCosThPiRest->SetBinContent(i+1, hCosThPiRest->GetBinContent(i+1)*scaleFactor/width);
      hCosThPiRest->SetBinError(i+1, fracErr*hCosThPiRest->GetBinContent(i+1));
    }

    hPthetaPi->Scale(scaleFactor, "width");
    hWQ2->Scale(scaleFactor, "width");

    for (int i = 0; i < hW->GetNbinsX()+1; i++) {
      double width = hW->GetBinLowEdge(i+2) - hW->GetBinLowEdge(i+1);
      if (hW->GetBinContent(i+1) == 0 || width == 0) continue;
      double fracErr = sqrt(1/hW->GetBinContent(i+1));
      hW->SetBinContent(i+1, hW->GetBinContent(i+1)*scaleFactor/width);
      hW->SetBinError(i+1, fracErr*hW->GetBinContent(i+1));
    }
  }

  for (int i = 0; i < hNp->GetNbinsX()+1; i++) {
    double width = hNp->GetBinLowEdge(i+2) - hNp->GetBinLowEdge(i+1);
    if (hNp->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hNp->GetBinContent(i+1));
    hNp->SetBinContent(i+1, hNp->GetBinContent(i+1)*scaleFactor/width);
    hNp->SetBinError(i+1, fracErr*hNp->GetBinContent(i+1));
  }

  for (int i = 0; i < hNn->GetNbinsX()+1; i++) {
    double width = hNn->GetBinLowEdge(i+2) - hNn->GetBinLowEdge(i+1);
    if (hNn->GetBinContent(i+1) == 0 || width == 0) continue;
    double fracErr = sqrt(1/hNn->GetBinContent(i+1));
    hNn->SetBinContent(i+1, hNn->GetBinContent(i+1)*scaleFactor/width);
    hNn->SetBinError(i+1, fracErr*hNn->GetBinContent(i+1));
  }
  
  // Make filename for selection
  if (sigType == 0) fileName += "_CC0pi";
  else if (sigType == 1) fileName += "_CC1pip";
  else if (sigType == 2) fileName += "_CC1pi0";

  // Make filename for not restricted
  if (!isRestricted) fileName += "_noKinCuts";

  fileName+="_outgoing.root";

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

  if (sigType == 1 || sigType == 2) {
    hPpi->Write();
    hThpi->Write();
    hThprotpi->Write();
    hThmupi->Write();
    hCosThPiRest->Write();

    hPthetaPi->Write();

    hW->Write();
    hWQ2->Write();
  }

  hNp->Write();
  hNn->Write();

  std::cout << "Wrote to " << fileName << std::endl;

  return;
};

// T2K not unfolded phase space restrictions
bool isT2KCCpip(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted) {

  if ((nvect->PartInfo(0))->fPID != 14) return false; 

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false;

  int pipCnt = 0; // counts number of pions
  int lepCnt = 0;

  TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
  TLorentzVector Ppip;
  TLorentzVector Pmu;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move on if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if ((abs(PID) >= 111 && abs(PID) <= 210) || (abs(PID) >= 212 && abs(PID) <= 557) || PID == -211) return false; 
    //else if (abs(PID) == 1114 || abs(PID) == 2114 || (abs(PID) >= 2214 && abs(PID) <= 5554)) return false; PHOTONS, NUCLEON, MULTINUCLEON OK
    else if (abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15 || abs(PID) == 17) {
      lepCnt++;
      Pmu = (nvect->PartInfo(j))->fP;
    }
    else if (PID == 211) {
      pipCnt++;
      Ppip = (nvect->PartInfo(j))->fP;
    }
  }

  if (pipCnt != 1) return false; 
  if (lepCnt != 1) return false;

  // relatively generic CC1pi+ definition done
  // now measurement specific (p_mu > 200 MeV, p_pi > 200 MeV, cos th_mu > 0.3,
  // cos th_pi > 0.3 in TRUE AND RECONSTRUCTED!

  if (restricted) {
    double p_mu = Pmu.Vect().Mag();
    double p_pi = Ppip.Vect().Mag();
    double cos_th_mu = cos(Pnu.Vect().Angle(Pmu.Vect()));
    double cos_th_pi = cos(Pnu.Vect().Angle(Ppip.Vect()));

    if (p_mu <= 200 || p_pi <= 200 || cos_th_mu <= 0.2 || cos_th_pi <= 0.2) return false;
  }


  return true;
};

bool isT2KCCpi0(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted) {

  if ((nvect->PartInfo(0))->fPID != 14) return false; 

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false;

  int pi0Cnt = 0; // counts number of pions
  int lepCnt = 0;

  TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
  TLorentzVector Ppi0;
  TLorentzVector Pmu;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move on if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if ((abs(PID) > 111 && abs(PID) <= 210) || (abs(PID) >= 212 && abs(PID) <= 557) || abs(PID) == 211) return false; 
    //else if (abs(PID) == 1114 || abs(PID) == 2114 || (abs(PID) >= 2214 && abs(PID) <= 5554)) return false; PHOTONS, NUCLEON, MULTINUCLEON OK
    else if (abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15 || abs(PID) == 17) {
      lepCnt++;
      Pmu = (nvect->PartInfo(j))->fP;
    }
    else if (PID == 111) {
      pi0Cnt++;
      Ppi0 = (nvect->PartInfo(j))->fP;
    }
  }

  if (pi0Cnt != 1) return false; 
  if (lepCnt != 1) return false;

  // relatively generic CC1pi+ definition done
  // now measurement specific (p_mu > 200 MeV, p_pi > 200 MeV, cos th_mu > 0.3,
  // cos th_pi > 0.3 in TRUE AND RECONSTRUCTED!

  if (restricted) {
    double p_mu = Pmu.Vect().Mag();
    double p_pi = Ppi0.Vect().Mag();
    double cos_th_mu = cos(Pnu.Vect().Angle(Pmu.Vect()));
    double cos_th_pi = cos(Pnu.Vect().Angle(Ppi0.Vect()));

    if (p_mu <= 200 || p_pi <= 200 || cos_th_mu <= 0.2 || cos_th_pi <= 0.2) return false;
  }


  return true;
};

// T2K not unfolded phase space restrictions
bool isT2KCC0pi(NeutVect *nvect, double EnuMin, double EnuMax, bool restricted) {

  if ((nvect->PartInfo(0))->fPID != 14) return false; 

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false;

  int lepCnt = 0;

  TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
  TLorentzVector Pmu;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move on if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if ((abs(PID) >= 111 && abs(PID) <= 210) || (abs(PID) >= 212 && abs(PID) <= 557) || abs(PID) == 211) return false; 
    //else if (abs(PID) == 1114 || abs(PID) == 2114 || (abs(PID) >= 2214 && abs(PID) <= 5554)) return false; PHOTONS, NUCLEON, MULTINUCLEON OK
    else if (abs(PID) == 11 || abs(PID) == 13 || abs(PID) == 15 || abs(PID) == 17) {
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

// Q2 reconstructed for CCQE
double Q2QErec(TLorentzVector pnu, TLorentzVector pmu, double binding){
  double el = pmu.E()/1000.;// energy of lepton in GeV
  double pl = pmu.Vect().Mag()/1000.; // momentum of lepton
  double ml = sqrt(el*el - pl*pl); // lepton mass

  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());

  double rEnu = EnuQErec(pnu, pmu, binding); //reconstructed neutrino energy
  double q2 = -ml*ml+2.*rEnu*(el-pl*cos(th_nu_mu));
  
  return q2;
};

// Enu reconstructed for CCQE
double EnuQErec(TLorentzVector pnu, TLorentzVector pmu, double binding){
  const double V  = binding/1000.;    // binding potential 
  const double mn = 0.93956536;       // neutron mass   
  const double mp = 0.93827203;       // proton mass    
  const double mn_eff = mn - V;

  double el = pmu.E()/1000.;
  double pl = pmu.Vect().Mag()/1000.; // momentum of lepton
  double ml = sqrt(el*el - pl*pl); // lepton mass

  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());
  
  double rEnu = (2*mn_eff*el - ml*ml + mp*mp - mn_eff*mn_eff)/
    (2*(mn_eff - el + pl*cos(th_nu_mu)));
  
  return rEnu;
};

// Q2 reconstructed for CC1pi+ (difference is 4 vector from pion)
double Q2CCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding) {
  
  double E_mu = pmu.E()/1000.;// energy of lepton in GeV
  double p_mu = pmu.Vect().Mag()/1000.; // momentum of lepton
  double m_mu = sqrt(E_mu*E_mu - p_mu*p_mu); // lepton mass
  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());

  // Just use the true Enu here
  //double rEnu = EnuCCpiprec(pnu, pmu, ppip, binding); //reconstructed neutrino energy
  double rEnu = pnu.E()/1000.;
  double q2 = -m_mu*m_mu + 2.*rEnu*(E_mu - p_mu*cos(th_nu_mu)); 
  
  return q2;
};

// Enu reconstructed for CC1pi+
double EnuCCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding){

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

double Wrec(TLorentzVector pnu, TLorentzVector pmu, double EHad, double binding) {
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
