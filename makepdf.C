#include "TH1D.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TKey.h"
#include "TROOT.h"
#include "TGaxis.h"
#include <iostream>
#include <algorithm>

#include <TF1.h>
#include <TLegend.h>

void loopdir(std::string file1) {


  TFile *f = TFile::Open(file1.c_str());

  TIter next(f->GetListOfKeys());
  TKey *key;
  std::string name;

  TCanvas *canv = new TCanvas("canv", "canv", 1280, 1024);
  canv->Print("asdf.pdf[");

  // Loop through all entries
  while ((key = (TKey*)next())) {

    std::string ClassName = key->GetClassName();
    TClass *cl = gROOT->GetClass(ClassName.c_str());

    // Get name of object
    name = std::string(key->GetName());
    std::cout << name << std::endl;
    std::cout << ClassName << std::endl;

    if (ClassName == "TH2D") {
      TH2D* plot2 = (TH2D*)(f->Get(name.c_str())->Clone());
      canv->cd();
      plot2->Draw("colz");
      canv->Print("asdf.pdf");
      delete plot2;
    }

    std::cout << "Done" << std::endl;

  } // end while

  canv->Print("asdf.pdf]");

  Output->Close();
}
