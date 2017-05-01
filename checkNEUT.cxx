#include <iostream>
#include <iomanip>

// ROOT includes
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TLorentzVector.h"

// NEUT includes
#include "neutpart.h"
#include "neutfsipart.h"
#include "neutvect.h"
#include "neutfsivert.h"
#include "neutvtx.h"
#include "neutrootTreeSingleton.h"


void Usage() {

  std::cout << "Wrong number of arguments!" << std::endl;
  std::cout << std::endl;
  std::cout << "./getKin.exe ROOT_FILE" << std::endl;
  std::cout << "./getKin.exe ROOT_FILE max_events" << std::endl;
  std::cout << std::endl;
  std::cout << "ROOT_FILE is the NEUT output vector file" << std::endl;
  std::cout << "max_events is the maximum number of events (if blank default to all)" << std::endl;

  return;
}

void checkNEUT(std::string fileName, int maxevents = -1) {

  // Open the TFile
  TFile *f = TFile::Open((fileName).c_str(),"open");
  TTree *tn = (TTree*)(f->Get("neuttree"));

  // Make an empty NEUT vector
  NeutVect *nvect = new NeutVect();
  tn->SetBranchAddress("vectorbranch", &nvect);

  // Make an empty NEUT vtx
  NeutVtx *nvtx = new NeutVtx();
  tn->SetBranchAddress("vertexbranch", &nvtx);

  // Some nice constants
  long unsigned int nevents = 0;
  if (maxevents == -1) {
    nevents = tn->GetEntries();
  } else {
    nevents = maxevents;
  }

  int targetA = -1;
  int targetZ = -1;
  double vnuclin = -1;
  double vnuclfi = -1;
  double pfsurf = -1;
  double pfmax = -1;
  int fluxid = -1;

  // Decides how many events we print
  int countwidth = nevents/1;

  // Holds the modes
  const int nmodes = 60;
  std::vector<int> modeVec(nmodes, 0);

  for (long unsigned int j = 0; j < nevents; j++) {

    tn->GetEntry(j);
    modeVec[nvect->Mode] += 1;

    // Print every 10000 event
    if (j % countwidth == 0) {
      std::cout << "==============================================" << std::endl;
      std::cout << "Event " << j << "/" << nevents << "(" << double(j)/ double(nevents)*100.0 << "%):" << std::endl;

      targetA = nvect->TargetA;
      targetZ = nvect->TargetA;
      vnuclin = nvect->VNuclIni;
      vnuclfi = nvect->VNuclFin;
      pfsurf = nvect->PFSurf;
      pfmax = nvect->PFMax;
      fluxid = nvect->FluxID;

      std::cout << std::setw(18) << std::left << "Intr. mode" << nvect->Mode << "\n";
      std::cout << "----" << std::endl;
      std::cout << std::setw(18) << std::left << "Nparts" << nvect->Npart() << "\n";

      // Print the particles
      for (int i = 0; i < nvect->Npart(); ++i) {
        std::cout << " Particle " << i+1 << "/" << nvect->Npart() << "\n";

        // Print PDG, mass and momentum
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "pdg" << (nvect->PartInfo(i))->fPID   << "\n";
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "mass" << (nvect->PartInfo(i))->fMass   << "\n";
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "4mom." << "(" << (nvect->PartInfo(i))->fP.Px() << ","
          << (nvect->PartInfo(i))->fP.Py() << "," 
          << (nvect->PartInfo(i))->fP.Pz() << ","
          << (nvect->PartInfo(i))->fP.E()  << ")"
          << "\n";

        // Print some vertex information
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "vertex" << nvect->VertexID(i) << "\n";
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "parent" << nvect->ParentIdx(i) << "\n";
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "parent pdg" << nvect->PartInfo(nvect->ParentIdx(i))->fPID << "\n";

        // Print alive and status
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "alive" << (nvect->PartInfo(i))->fIsAlive << "\n";
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "status" << (nvect->PartInfo(i))->fStatus  << "\n";

        // Print some particle position information
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "init." << "(" << (nvect->PartInfo(i))->fPosIni.X() << ","
          << (nvect->PartInfo(i))->fPosIni.Y() << "," 
          << (nvect->PartInfo(i))->fPosIni.Z() << ","
          << (nvect->PartInfo(i))->fPosIni.T()  << ")"
          << "\n";
        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "fin." << "(" << (nvect->PartInfo(i))->fPosFin.Px() << ","
          << (nvect->PartInfo(i))->fPosFin.Y() << "," 
          << (nvect->PartInfo(i))->fPosFin.Z() << ","
          << (nvect->PartInfo(i))->fPosFin.T()  << ")"
          << "\n";
      } // Finish looping over the normal particles

      std::cout << "----" << std::endl;
      std::cout << std::setw(18) << std::left << "Nvtx" << nvtx->Nvtx() << "\n";
      // Print some vertex information
      for (int i = 0 ; i < nvtx->Nvtx(); i++) {
        std::cout << " Vertex " << i+1 << "/" << nvtx->Nvtx() << "\n";

        std::cout << std::setw(3) << std::left << " " << std::setw(15) << "vtx pos" << "(" << (nvtx->Pos(i))->X() << ","
          << (nvtx->Pos(i))->Y() << "," 
          << (nvtx->Pos(i))->Z() << ","
          << (nvtx->Pos(i))->T()  << ")"
          << "\n";
      }

      // Print some FSI information

    } // Finish if 

  } // Finish event foor loop
  std::cout << "==============================================" << std::endl;

  std::cout << "Mode summary: " << std::endl;
  std::cout << std::right << std::setw(12) << "Mode" << std::setw(10) << "Counts" << std::setw(12) << "Percent" << std::endl;
  for (size_t i = 0; i < modeVec.size(); ++i) {
    std::cout << std::right << std::setw(10) << i << std::setw(10) << modeVec[i] << std::setw(10) << "(" << std::setprecision(2) << double(modeVec[i])/double(nevents)*100. << "%)" << std::endl;
  }
  std::cout << std::endl;
  std::cout << "==============================================" << std::endl;

  std::cout << "Nuclear summary: " << std::endl;
  std::cout << std::setw(10) << std::left << "TargetA " << targetA << std::endl;
  std::cout << std::setw(10) << std::left << "TargetZ " << targetZ << std::endl;
  std::cout << std::setw(10) << std::left << "VNuclIn " << vnuclin << std::endl;
  std::cout << std::setw(10) << std::left << "VNuclFi " << vnuclfi << std::endl;
  std::cout << std::setw(10) << std::left << "pF surf " << pfsurf << std::endl;
  std::cout << std::setw(10) << std::left << "pF max " << pfmax << std::endl;
  std::cout << std::setw(10) << std::left << "flux ID " << fluxid << std::endl;
  std::cout << "==============================================" << std::endl;

};

int main(int argc, char* argv[]) {

  if (argc != 2 && argc != 3) {
    Usage();
    exit(-1);
  }

  if (argc == 2) {
    checkNEUT(std::string(argv[1]));
  } else if (argc == 3) {
    checkNEUT(std::string(argv[1]), std::atoi(argv[2]));
  }

  return 0;
};
