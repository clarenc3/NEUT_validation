#include "getKin_nuc.h"

int main(int argc, char* argv[]) {

  if (argc != 3) {
    Usage();
    return -1;
  }

  getKin(std::string(argv[1]), std::atoi(argv[2]));

  return 0;
};

void Usage() {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << "./getKin_nuc.exe ROOT_FILE SIGNAL" << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "SIGNAL is 0 = CC1pi+1p, 1 = CC1pi+1n, 2 = CC1pi0, 3 = CC1pi-1p, 4 = CC1pi-1n" << std::endl;

  return;
}

void getKin(std::string fileName, int sigType) {

  if (sigType == 0) std::cout << "Signal is CC1pi+1p" << std::endl;
  else if (sigType == 1) std::cout << "Signal is CC1pi+1n" << std::endl;
  else if (sigType == 2) std::cout << "Signal is CC1pi0" << std::endl;
  else if (sigType == 3) std::cout << "Signal is CC1pi-1p" << std::endl;
  else if (sigType == 4) std::cout << "Signal is CC1pi-1n" << std::endl;
  else {
    std::cout << "Unrecognised signal! Exiting..." << std::endl;
    exit(-1);
  }

  double EnuMax = 1000.;
  int nBins = 30;

  double minAng = -1;

  TFile *f = TFile::Open((fileName).c_str(),"open");
  TTree *tn = (TTree*)(f->Get("neuttree"));

  TH1D *fluxHist = NULL;
  TH1D *eventHist = NULL;

  if (sigType < 3) {
    fluxHist = (TH1D*)f->Get("flux_numu");
    eventHist = (TH1D*)f->Get("evtrt_numu");
  } else {
    fluxHist = (TH1D*)f->Get("flux_numub");
    eventHist = (TH1D*)f->Get("evtrt_numub");
  }

  if (!fluxHist || !eventHist) {
    std::cerr << "***********************************\n";
    std::cerr << "COULDN'T FIND FLUXHIST OR EVENTHIST\n";
    std::cerr << fileName << std::endl;
    std::cerr << "***********************************" << std::endl;
    exit(-1);
  }

  double maxMom = eventHist->GetBinLowEdge(eventHist->GetNbinsX()+1);
  maxMom = 2;

  double nevents = tn->GetEntries();

  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch",&nvect);

  //int minBin = fluxHist->GetXaxis()->FindBin(Double_t(0));
  //int maxBin = fluxHist->GetXaxis()->FindBin(Double_t(8));
  // Get integral over custom range
  //double integral = fluxHist->Integral(minBin, maxBin+1, "width");
  double scaleFactor = (eventHist->Integral("width")*1.E-38)/(nevents*fluxHist->Integral("width"))*(16./8.);
  //double scaleFactor = eventHist->Integral("width")*1.E-38/(nevents*integral)*(16./8.);


  TH1D *hPmu = new TH1D("hPmu","hPmu", nBins, 0, maxMom);
  hPmu->GetXaxis()->SetTitle("p_{#mu} (GeV/c)");
  hPmu->GetYaxis()->SetTitle("d#sigma/dp_{#mu} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPnuc = new TH1D("hPnuc","hPnuc", nBins, 0, maxMom);
  hPnuc->GetXaxis()->SetTitle("p_{p} (GeV/c)");
  hPnuc->GetYaxis()->SetTitle("d#sigma/dp_{p} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hThmu = new TH1D("hThmu","hThmu", nBins, minAng, 1);
  hThmu->GetXaxis()->SetTitle("cos#theta_{#mu}");
  hThmu->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu} (cm^{2}/nucleon)");

  TH1D *hThnuc = new TH1D("hThnuc","hThnuc", nBins, -1, 1);
  hThnuc->GetXaxis()->SetTitle("cos#theta_{p}");
  hThnuc->GetYaxis()->SetTitle("d#sigma/dcos#theta_{p} (cm^{2}/nucleon)");

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

  TH1D *hEnu = new TH1D("hEnu","hEnu",nBins, 0, fluxHist->GetBinLowEdge(fluxHist->GetNbinsX()));
  hEnu->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  hEnu->GetYaxis()->SetTitle("#sigma (E_{#nu})");

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

  TH2D *hPthetaNuc = new TH2D("hPthetaNuc", "hPthetaNuc", nBins, 0, maxMom, nBins, -1, 1);
  hPthetaNuc->GetXaxis()->SetTitle("p_{p} (GeV/c)");
  hPthetaNuc->GetYaxis()->SetTitle("cos#theta_{p}");
  hPthetaNuc->GetZaxis()->SetTitle("d^{2}#sigma/dp_{p}dcos#theta_{p} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hPpi = new TH1D("hPpi","hPpi", nBins, 0, maxMom);
  hPpi->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
  hPpi->GetYaxis()->SetTitle("d#sigma/dp_{#pi} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hThpi = new TH1D("hThpi","hThpi", nBins, minAng, 1);
  hThpi->GetXaxis()->SetTitle("cos#theta_{#pi}");
  hThpi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi} (cm^{2}/nucleon)");

  TH1D *hThnucpi = new TH1D("hThnucpi","hThnucpi", nBins, -1, 1);
  hThnucpi->GetXaxis()->SetTitle("cos#theta_{p,#pi}");
  hThnucpi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{p#pi} (cm^{2}/nucleon)");

  TH1D *hThmupi = new TH1D("hThmupi","hThmupi", nBins, -1, 1);
  hThmupi->GetXaxis()->SetTitle("cos#theta_{#mu,#pi}");
  hThmupi->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#mu,#pi} (cm^{2}/nucleon)");

  TH2D *hPthetaPi = new TH2D("hPthetaPi", "hPthetaPi", nBins, 0, maxMom, nBins, minAng, 1);
  hPthetaPi->GetXaxis()->SetTitle("p_{#pi} (GeV/c)");
  hPthetaPi->GetYaxis()->SetTitle("cos#theta_{#pi}");
  hPthetaPi->GetZaxis()->SetTitle("d^{2}#sigma/dp_{#pi}dcos#theta_{#pi} (cm^{2}/nucleon/(GeV/c))");

  TH1D *hNNuc = new TH1D("hNNuc","hNNuc", 6, 0, 6);
  hNNuc->GetXaxis()->SetTitle("N_{p}");
  hNNuc->GetYaxis()->SetTitle("d#sigma/dN_{p} (cm^{2}/nucleon)");

  // THETA PION IN RESONANCE REST FRAME
  TH1D *hCosThPiRest = new TH1D("hCosThPiRest","hCosThPiRest", 25, -1, 1);
  hCosThPiRest->GetXaxis()->SetTitle("cos#theta_{#pi}^{RES}");
  hCosThPiRest->GetYaxis()->SetTitle("d#sigma/dcos#theta_{#pi}^{RES} (cm^{2}/nucleon/1)");

  // THETA PION IN RESONANCE REST FRAME
  TH1D *hPhiRest = new TH1D("hPhiRest","hPhiRest", 25, 0, 2*M_PI);
  hPhiRest->GetXaxis()->SetTitle("#phi_{#pi}^{RES}");
  hPhiRest->GetYaxis()->SetTitle("d#sigma/d#phi_{#pi}^{RES} (cm^{2}/nucleon/1)");

  int eventCnt = 0;
  int printwidth = double(nevents)/double(25);

  TStopwatch clock;
  clock.Start();

  // Needed for Adler angle
  TRandom3 *rand = new TRandom3();

  for (int j = 0; j < nevents; j++) {

    tn->GetEntry(j);

    if (j%printwidth == 0) {
      std::cout << "On event #" << std::setw(8) << j << "/" << std::left << nevents << " (" << double(j)/double(nevents) * 100. << "%)" << std::endl;
    }

    // Very simple signal definition
    if (sigType == 0) {
      if (!isCC1ppip(nvect, EnuMax)) continue;
    } else if (sigType == 1) {
      if (!isCC1npip(nvect, EnuMax)) continue;
    } else if (sigType == 2) {
      if (!isCC1pi0(nvect, EnuMax)) continue;
    } else if (sigType == 3) {
      if (!isCC1ppim(nvect, EnuMax)) continue;
    } else if (sigType == 4) {
      if (!isCC1npim(nvect, EnuMax)) continue;
    }

    eventCnt++;

    TLorentzVector Pnu = (nvect->PartInfo(0))->fP;
    TLorentzVector Pp_init = (nvect->PartInfo(1))->fP;
    TLorentzVector Pnuc;
    TLorentzVector Ppi;
    TLorentzVector Pmu;

    int nNuc = 0;

    // Loop over the particle stack
    for (int k = 2; k < nvect->Npart(); ++k) {

      int PID = (nvect->PartInfo(k))->fPID;

      if (!(nvect->PartInfo(k))->fIsAlive && (nvect->PartInfo(k))->fStatus != 0) continue;

      // select the charged pion
      if ((PID == 211 && (sigType == 0 || sigType == 1)) || (PID == 111 && sigType == 2) || (PID == -211 && (sigType == 3 || sigType == 4))) {
        Ppi = nvect->PartInfo(k)->fP;
        // select highest momentum nucleon
      } else if (PID == 2212 || PID == 2112) {
        nNuc++;
        if (nvect->PartInfo(k)->fP.Vect().Mag() > Pnuc.Vect().Mag()) {
          Pnuc = nvect->PartInfo(k)->fP;
        }
      } else if ((PID == -13 && sigType > 2) || (PID == 13 && sigType < 3)) {
        Pmu = (nvect->PartInfo(k))->fP;  
      }
    }

    double W = MpPi(Pnuc, Ppi)/1.E3;
    // Minoo wanted a W cut
    if (W > 2.0) continue;

    double pmu = Pmu.Vect().Mag()/1000.;
    double thmu = cos(Pnu.Vect().Angle(Pmu.Vect()));

    double pnuc = Pnuc.Vect().Mag()/1000.;
    double thnuc = cos(Pnu.Vect().Angle(Pnuc.Vect()));

    //Enu = EnuCCpiprec(Pnu, Pmu, Ppi);
    //Q2 = Q2CCpiprec(Pnu, Pmu, Ppi);
    // Only want MC true quantities, not reconstructed
    double Enu = Pnu.E()/1.E3;
    double Q2 = -1*((Pmu-Pnu).Mag2())/1.E6;
    TLorentzVector q = Pmu-Pnu;

    double ppi;
    double thpi;
    double thnucpi;
    double thmupi;

    double costhAdler;
    double phiAdler;

    // Do the simple angles etc
    ppi = Ppi.Vect().Mag()/1000.;
    thpi = cos(Pnu.Vect().Angle(Ppi.Vect()));
    thnucpi = cos(Ppi.Vect().Angle(Pnuc.Vect()));
    thmupi = cos(Pmu.Vect().Angle(Ppi.Vect()));

    // Do the Adler angles
    // The resonance rest frame
    TLorentzVector PResRest = Pnuc + Ppi;
    // The particles in the rest frame
    TLorentzVector PnuCopy = Pnu;
    TLorentzVector PmuCopy = Pmu;
    TLorentzVector PpiCopy = Ppi;

    // Boost into the resonance rest frame
    PnuCopy.Boost(PResRest.BoostVector());
    PmuCopy.Boost(-PResRest.BoostVector());
    PpiCopy.Boost(-PResRest.BoostVector());
    costhAdler = cos(PpiCopy.Angle(PResRest.Vect()));

    // Get the vectors from the 4-vector
    TVector3 PmuVect = PmuCopy.Vect();
    TVector3 PnuVect = PnuCopy.Vect();
    TVector3 PresVect = PResRest.Vect();
    TVector3 PpiVect = PpiCopy.Vect();

    // Define the z-direction
    TVector3 zVect = (PnuVect-PmuVect).Unit();
    // Define y direction as being z (resonance direction) x pmu*
    TVector3 yVect = (zVect.Cross(PmuVect)).Unit();
    // define x direction as being y X z
    TVector3 xVect = (yVect.Cross(zVect)).Unit();

    // Project pion onto z axis
    TVector3 PpiVectZ = zVect * PpiVect.Dot(zVect);
    // Then subtract this vector off the pion vector
    TVector3 PpiVectPlane = PpiVect - PpiVectZ;

    if (PpiVectPlane.Y() > 0) {
      phiAdler = PpiVectPlane.Angle(xVect);
    } else if (PpiVectPlane.Y() < 0) {
      phiAdler = (2*M_PI-PpiVectPlane.Angle(xVect));
    } else if (PpiVectPlane.Y() == 0) {
      double randNo = rand->Rndm();
      if (randNo > 0.5) {
        phiAdler = PpiVectPlane.Angle(xVect);
      } else {
        phiAdler = (2*M_PI-PpiVectPlane.Angle(xVect));
      }
    }

    hPmu->Fill(pmu);
    hThmu->Fill(thmu);

    hPnuc->Fill(pnuc);
    hThnuc->Fill(thnuc);

    hPthetaMu->Fill(pmu, thmu);
    hPthetaNuc->Fill(pnuc, thnuc);

    hPmuQ2->Fill(pmu, Q2);
    hThetaQ2->Fill(thmu, Q2);

    hEnu->Fill(Enu);
    hQ2->Fill(Q2);

    hPpi->Fill(ppi);
    hThpi->Fill(thpi);
    hThnucpi->Fill(thnucpi);
    hThmupi->Fill(thmupi);

    hPthetaPi->Fill(ppi, thpi);

    hPhiRest->Fill(phiAdler);
    hCosThPiRest->Fill(costhAdler);

    hW->Fill(W);

    hWQ2->Fill(W, Q2);

    hNNuc->Fill(nNuc);

  }

  clock.Stop();

  std::cout << clock.RealTime() << " seconds " << std::endl;
  std::cout << eventCnt << "/" << nevents << " events" << std::endl;

  hPmu->Sumw2();
  hPmu->Scale(scaleFactor, "width");
  hThmu->Sumw2();
  hThmu->Scale(scaleFactor, "width");

  hPnuc->Sumw2();
  hPnuc->Scale(scaleFactor, "width");
  hThnuc->Sumw2();
  hThnuc->Scale(scaleFactor, "width");

  hPthetaMu->Sumw2();
  hPthetaMu->Scale(scaleFactor, "width");
  hPthetaNuc->Sumw2();
  hPthetaNuc->Scale(scaleFactor, "width");

  hPmuQ2->Sumw2();
  hPmuQ2->Scale(scaleFactor, "width");
  hThetaQ2->Sumw2();
  hThetaQ2->Scale(scaleFactor, "width");

  TH1D* fineFlux = FluxUnfoldedScaling(hEnu, fluxHist);
  hEnu->Scale(scaleFactor*fluxHist->Integral("width"));
  hQ2->Sumw2();
  hQ2->Scale(scaleFactor, "width");

  hPpi->Sumw2();
  hPpi->Scale(scaleFactor, "width");
  hThpi->Sumw2();
  hThpi->Scale(scaleFactor, "width");
  hThnucpi->Sumw2();
  hThnucpi->Scale(scaleFactor, "width");
  hThmupi->Sumw2();
  hThmupi->Scale(scaleFactor, "width");

  hCosThPiRest->Sumw2();
  hCosThPiRest->Scale(scaleFactor, "width");
  hPhiRest->Sumw2();
  hPhiRest->Scale(scaleFactor, "width");

  hW->Sumw2();
  hW->Scale(scaleFactor, "width");
  hWQ2->Sumw2();
  hWQ2->Scale(scaleFactor, "width");

  hPthetaPi->Sumw2();
  hPthetaPi->Scale(scaleFactor, "width");

  hNNuc->Sumw2();
  hNNuc->Scale(scaleFactor, "width");


  // Make filename for selection
  // fileName is now the output filename; strip the .root ending
  while (fileName.find(".root") != std::string::npos) {
    fileName = fileName.substr(0, fileName.find(".root"));
  }
  if (sigType == 0) fileName += "_CC1ppip";
  else if (sigType == 1) fileName += "_CC1npip";
  else if (sigType == 2) fileName += "_CC1pi0";
  else if (sigType == 3) fileName += "_CC1ppim";
  else if (sigType == 4) fileName += "_CC1npim";

  fileName+="_outgoing.root";

  TFile *output = new TFile(fileName.c_str(), "recreate");

  output->cd();

  fluxHist->Write();
  fineFlux->Write();
  eventHist->Write();

  hPmu->Write();
  hThmu->Write();

  hPnuc->Write();
  hThnuc->Write();

  hPthetaMu->Write();
  hPthetaNuc->Write();

  hPmuQ2->Write();
  hThetaQ2->Write();

  hEnu->Write();
  hQ2->Write();

  hPpi->Write();
  hThpi->Write();
  hThnucpi->Write();
  hThmupi->Write();

  hCosThPiRest->Write();
  hPhiRest->Write();

  hPthetaPi->Write();

  hW->Write();
  hWQ2->Write();

  hNNuc->Write();

  std::cout << "Wrote to " << fileName << ", bye!" << std::endl;

  return;
};

// This interpolates the flux by a TGraph instead of requiring the flux and MC flux to have the same binning
//******************************************************************** 
TH1D* FluxUnfoldedScaling(TH1D* mcHist, TH1D* fluxHist) {
  //******************************************************************** 

  if (true) {

    // Make a temporary TGraph which holds the points from the flux (essentially copying the TH1D to a TGraph)
    TGraph* fluxGraph = new TGraph(fluxHist->GetNbinsX());
    for (int i = 0; i < fluxHist->GetNbinsX()+1; ++i){
      fluxGraph->SetPoint(i, fluxHist->GetXaxis()->GetBinCenter(i+1), fluxHist->GetBinContent(i+1));
    }

    // Increase the resolution for the interpolation used for the flux
    int resolution = 50*fluxHist->GetXaxis()->GetNbins();
    double range = fluxHist->GetXaxis()->GetBinLowEdge(fluxHist->GetNbinsX()) - fluxHist->GetXaxis()->GetBinLowEdge(1);
    double criterion = range/double(resolution);
    std::cout << "Range in fluxHist " << fluxHist->GetXaxis()->GetBinLowEdge(fluxHist->GetNbinsX()) << " - " << fluxHist->GetXaxis()->GetBinLowEdge(1) << " = " << range << " with " << resolution << " points: " << criterion << std::endl;

    // The new interpolated flux histogram with fine binning
    TH1D* fineFlux = new TH1D("fineFlux", "fineFlux", resolution, fluxHist->GetXaxis()->GetBinLowEdge(1), fluxHist->GetXaxis()->GetBinLowEdge(fluxHist->GetNbinsX()+1));

    // Set the new TH1D with the TGraph interpolated bin content
    // Set the first bin to the fluxHist first bin content; seems like there's strange interpolation at the edge
    for (int i = 0; i < fineFlux->GetNbinsX()+1; i++) {

      // The first few bins (those between the first bin's low edge and that bin's center) will be skewed and muck up interpolation; 
      if (fluxHist->GetBinCenter(1) > fineFlux->GetXaxis()->GetBinCenter(i+1)) {
        fineFlux->SetBinContent(i+1, fluxHist->GetBinContent(1));
      } else {
        fineFlux->SetBinContent(i+1, fluxGraph->Eval(fineFlux->GetXaxis()->GetBinCenter(i+1), 0, "s"));
      }
    }

    // Loop over the mcHist and match up the bins with the new fluxHist
    for (int i = 1; i < mcHist->GetNbinsX()+1; i++) {
      // Get the low edge of the ith bin
      Double_t binLowEdge = mcHist->GetBinLowEdge(i);
      // Get the high edge of the ith bin
      Double_t binHighEdge = mcHist->GetBinLowEdge(i+1);

      // Find the correpsonding bin in the interpolated flux histogram
      // Start by finding the low edge of the bin
      Int_t fluxLow = 0;
      Int_t fluxHigh = 0;

      // Find the binLowEdge in the new interpolated flux histogram which matches the mc binning's edge
      double diff = 999; // Set the initial difference to be very large
      for (; fluxLow < fineFlux->GetNbinsX()+1; fluxLow++) {
        // the difference between the mc bin edge and the flux bin edge
        double temp = fabs(binLowEdge - fineFlux->GetBinLowEdge(fluxLow));
        // if difference is larger than previous
        if (temp < diff) {
          diff = temp;
        } else {
          break;
        }
      }

      // 0.5 GeV doesn't sound right!
      if (diff > criterion) {
        // This is a known issue for all BEBC 1pi 1DEnu classes in the first bin, where the flux histogram starts at 8.9GeV but MC binning starts at 5GeV
        std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << std::endl;
        std::cerr << "Couldn't find good low-edge bin match for flux histogram " << fluxHist->GetName() << " and MC " << mcHist->GetName() << std::endl;
        std::cout << "fluxLow = " << fluxLow << std::endl;
        std::cout << "binLowEdge - fineFlux->GetBinLowEdge(fluxLow) = " << binLowEdge << " - " << fineFlux->GetBinLowEdge(fluxLow) << " = " << binLowEdge - fineFlux->GetBinLowEdge(fluxLow) << std::endl;
      }

      diff = 999;
      for (; fluxHigh < fineFlux->GetNbinsX()+1; fluxHigh++) {
        double temp = fabs(binHighEdge - fineFlux->GetBinLowEdge(fluxHigh));
        if (temp < diff) {
          diff = temp;
        } else {
          break;
        }
      }

      if (diff > criterion) {
        // This is a known issue for anti-nu BEBC 1pi 1DEnu classes in the last bin, where the flux histogram ends at 180 GeV but MC binning continues to 200 GeV
        std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << std::endl;
        std::cerr << "Couldn't find good high-edge bin match for flux histogram " << fluxHist->GetName() << " and MC " << mcHist->GetName() << std::endl;
        std::cout << "fluxHigh = " << fluxHigh << std::endl;
        std::cout << "binHighEdge - fineFlux->GetBinLowEdge(fluxHigh) = " << binHighEdge << " - " << fineFlux->GetBinLowEdge(fluxHigh) << " = " << binHighEdge - fineFlux->GetBinLowEdge(fluxHigh) << std::endl;
      }

      // fluxHigh - 1 because Integral takes the binLowEdge into account, and fluxHigh is our high edge
      double fluxInt = fineFlux->Integral(fluxLow, fluxHigh - 1, "width");

      // Scale the bin content in bin i by the flux integral in that bin
      if (fluxInt == 0) continue;
      mcHist->SetBinContent(i, mcHist->GetBinContent(i)/fluxInt);
      mcHist->SetBinError(i, mcHist->GetBinError(i)/fluxInt);
    }

    delete fluxGraph;

    return fineFlux;

  } else {

    std::cout << "hello" << std::endl;
    double range = fluxHist->GetXaxis()->GetBinLowEdge(fluxHist->GetNbinsX()) - fluxHist->GetXaxis()->GetBinLowEdge(1);
    double criterion = range/double(fluxHist->GetXaxis()->GetNbins());

    for (int i = 1; i <= mcHist->GetNbinsX()+2; i++) {
      // Get the low edge of the ith bin
      Double_t binLowEdge = mcHist->GetBinLowEdge(i);
      // Get the high edge of the ith bin
      Double_t binHighEdge = mcHist->GetBinLowEdge(i+1);

      // Find the correpsonding bin in the interpolated flux histogram
      // Start by finding the low edge of the bin
      Int_t fluxLow = 0;
      Int_t fluxHigh = 0;

      // Find the binLowEdge in the new interpolated flux histogram which matches the mc binning's edge
      double diff = 999; // Set the initial difference to be very large
      for (; fluxLow < fluxHist->GetNbinsX()+1; fluxLow++) {
        // the difference between the mc bin edge and the flux bin edge
        double temp = fabs(binLowEdge - fluxHist->GetBinLowEdge(fluxLow));
        // if difference is larger than previous
        if (temp < diff) {
          diff = temp;
        } else {
          std::cout << binLowEdge << " " << fluxHist->GetBinLowEdge(fluxLow) << std::endl;
          break;
        }
      }

      // 0.5 GeV doesn't sound right!
      if (diff > criterion) {
        // This is a known issue for all BEBC 1pi 1DEnu classes in the first bin, where the flux histogram starts at 8.9GeV but MC binning starts at 5GeV
        std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << std::endl;
        std::cerr << "Couldn't find good low-edge bin match for flux histogram " << fluxHist->GetName() << " and MC " << mcHist->GetName() << std::endl;
        std::cout << "fluxLow = " << fluxLow << std::endl;
        std::cout << "binLowEdge - fluxHist->GetBinLowEdge(fluxLow) = " << binLowEdge << " - " << fluxHist->GetBinLowEdge(fluxLow) << " = " << binLowEdge - fluxHist->GetBinLowEdge(fluxLow) << std::endl;
      }

      diff = 999;
      for (; fluxHigh < fluxHist->GetNbinsX()+1; fluxHigh++) {
        double temp = fabs(binHighEdge - fluxHist->GetBinLowEdge(fluxHigh));
        if (temp < diff) {
          diff = temp;
        } else {
          std::cout << binHighEdge << " " << fluxHist->GetBinLowEdge(fluxHigh) << std::endl;
          break;
        }
      }

      if (diff > criterion) {
        // This is a known issue for anti-nu BEBC 1pi 1DEnu classes in the last bin, where the flux histogram ends at 180 GeV but MC binning continues to 200 GeV
        std::cerr << "Warning " << __FILE__ << ":" << __LINE__ << std::endl;
        std::cerr << "Couldn't find good high-edge bin match for flux histogram " << fluxHist->GetName() << " and MC " << mcHist->GetName() << std::endl;
        std::cout << "fluxHigh = " << fluxHigh << std::endl;
        std::cout << "binHighEdge - fluxHist->GetBinLowEdge(fluxHigh) = " << binHighEdge << " - " << fluxHist->GetBinLowEdge(fluxHigh) << " = " << binHighEdge - fluxHist->GetBinLowEdge(fluxHigh) << std::endl;
      }

      // fluxHigh - 1 because Integral takes the binLowEdge into account, and fluxHigh is our high edge
      double fluxInt = fluxHist->Integral(fluxLow, fluxHigh-1, "width");

      // Scale the bin content in bin i by the flux integral in that bin
      if (fluxInt == 0) continue;
      mcHist->SetBinContent(i, mcHist->GetBinContent(i)/fluxInt);
      mcHist->SetBinError(i, mcHist->GetBinError(i)/fluxInt);
      std::cout << i << " " << mcHist->GetBinContent(i) << std::endl;
    }
    mcHist->SetBinContent(mcHist->GetNbinsX(), 0);
    mcHist->SetBinError(mcHist->GetNbinsX(), 0);

    return fluxHist;

  }

  return NULL;


};

bool isCC1ppip(NeutVect *nvect, double EnuMax) {

  double EnuMin = 0.;

  if ((nvect->PartInfo(0))->fPID != 14) return false;

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false; 

  int pipCnt = 0;
  int lepCnt = 0;
  int protonCnt = 0;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move to next particle if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if (PID == 13)
      lepCnt++;
    else if (PID == 211)
      pipCnt++;
    else if (PID == 2212 /*&& nvect->Ibound == 0*/)
      protonCnt++;
    else
      return false; // require only three prong events! (allow photons?)
  }

  // don't think there's away of implementing spectator proton cuts in NEUT?
  // 100 MeV or larger protons

  if (pipCnt != 1) return false;
  if (lepCnt != 1) return false;
  if (protonCnt != 1) return false;

  return true;
}

bool isCC1npip(NeutVect *nvect, double EnuMax) {
  double EnuMin = 0.;

  if ((nvect->PartInfo(0))->fPID != 14) return false;

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false; 

  int pipCnt = 0;
  int lepCnt = 0;
  int neutronCnt = 0;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move to next particle if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if (PID == 13)
      lepCnt++;
    else if (PID == 211)
      pipCnt++;
    else if (PID == 2112 /*&& nvect->Ibound == 0*/)
      neutronCnt++;
    else
      return false; // require only three prong events! (allow photons?)
  }

  // don't think there's away of implementing spectator proton cuts in NEUT?
  // 100 MeV or larger protons

  if (pipCnt != 1) return false;
  if (lepCnt != 1) return false;
  if (neutronCnt != 1) return false;

  return true;
}

bool isCC1pi0(NeutVect *nvect, double EnuMax) {
  double EnuMin = 0.;

  if ((nvect->PartInfo(0))->fPID != 14) return false;

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != 13) && ((nvect->PartInfo(3))->fPID != 13)) return false; 

  int pi0Cnt = 0;
  int lepCnt = 0;
  int protonCnt = 0;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move to next particle if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if (PID == 13)
      lepCnt++;
    else if (PID == 111)
      pi0Cnt++;
    else if (PID == 2212 /*&& nvect->Ibound == 0*/)
      protonCnt++;
    else
      return false; // require only three prong events! (allow photons?)
  }

  // don't think there's away of implementing spectator proton cuts in NEUT?
  // 100 MeV or larger protons

  if (pi0Cnt != 1) return false;
  if (lepCnt != 1) return false;
  if (protonCnt != 1) return false;

  return true;
}

bool isCC1ppim(NeutVect *nvect, double EnuMax) {
  double EnuMin = 0.;

  if ((nvect->PartInfo(0))->fPID != -14) return false;

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != -13) && ((nvect->PartInfo(3))->fPID != -13)) return false; 

  int pimCnt = 0;
  int lepCnt = 0;
  int protonCnt = 0;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move to next particle if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if (PID == -13)
      lepCnt++;
    else if (PID == -211)
      pimCnt++;
    else if (PID == 2212 /*&& nvect->Ibound == 0*/)
      protonCnt++;
    else
      return false; // require only three prong events! (allow photons?)
  }

  // don't think there's away of implementing spectator proton cuts in NEUT?
  // 100 MeV or larger protons

  if (pimCnt != 1) return false;
  if (lepCnt != 1) return false;
  if (protonCnt != 1) return false;

  return true;
}

bool isCC1npim(NeutVect *nvect, double EnuMax) {
  double EnuMin = 0.;

  if ((nvect->PartInfo(0))->fPID != -14) return false;

  if (((nvect->PartInfo(0))->fP.E() < EnuMin*1000.) || ((nvect->PartInfo(0))->fP.E() > EnuMax*1000.)) return false; 

  if (((nvect->PartInfo(2))->fPID != -13) && ((nvect->PartInfo(3))->fPID != -13)) return false; 

  int pimCnt = 0;
  int lepCnt = 0;
  int neutronCnt = 0;

  for (int j = 2; j < nvect->Npart(); j++) {
    if (!((nvect->PartInfo(j))->fIsAlive) && (nvect->PartInfo(j))->fStatus != 0) continue; //move to next particle if NOT ALIVE and NOT NORMAL
    int PID = (nvect->PartInfo(j))->fPID;
    if (PID == -13)
      lepCnt++;
    else if (PID == -211)
      pimCnt++;
    else if (PID == 2112 /*&& nvect->Ibound == 0*/)
      neutronCnt++;
    else
      return false; // require only three prong events! (allow photons?)
  }

  // don't think there's away of implementing spectator proton cuts in NEUT?
  // 100 MeV or larger protons

  if (pimCnt != 1) return false;
  if (lepCnt != 1) return false;
  if (neutronCnt != 1) return false;

  return true;
}

// Q2 reconstructed for CC1pi+ (difference is 4 vector from pion)
double Q2CCpiprec(TLorentzVector pnu, TLorentzVector pmu, TLorentzVector ppip, double binding) {

  double E_mu = pmu.E()/1000.;// energy of lepton in GeV
  double p_mu = pmu.Vect().Mag()/1000.; // momentum of lepton
  double m_mu = sqrt(E_mu*E_mu - p_mu*p_mu); // lepton mass
  double th_nu_mu = pnu.Vect().Angle(pmu.Vect());

  double rEnu = EnuCCpiprec(pnu, pmu, ppip, binding); //reconstructed neutrino energy
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

double MpPi(TLorentzVector pp, TLorentzVector ppi) {
  double E_p = pp.E();
  double p_p = pp.Vect().Mag();
  double m_p = sqrt(E_p*E_p - p_p*p_p);

  double E_pi = ppi.E();
  double p_pi = ppi.Vect().Mag();
  double m_pi = sqrt(E_pi*E_pi - p_pi*p_pi);

  double th_p_pi = pp.Vect().Angle(ppi.Vect());

  // fairly easy thing to derive since bubble chambers measure the proton!
  double invMass = sqrt(m_p*m_p + m_pi*m_pi + 2*E_p*E_pi - 2*p_pi*p_p*cos(th_p_pi));

  return invMass;
};
