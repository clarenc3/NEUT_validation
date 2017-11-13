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

void makeratio(std::string file1, std::string file2) {

  std::cout << "Minoo, NEUT 5.3.3" << std::endl;

  TFile *f = TFile::Open(file1.c_str());
  TFile *f2 = TFile::Open(file2.c_str());

  // Strip out leading underscores from file1 and file2
  file1.erase(0, 3);
  file2.erase(0, 3);

  while (file1.find("/") != std::string::npos) {
    file1.replace(file1.find("/"), 1, std::string("_"));
  }

  while (file2.find("/") != std::string::npos) {
    file2.replace(file2.find("/"), 1, std::string("_"));
  }

  while (file1.find(".root") != std::string::npos) {
    file1.erase(file1.find(".root"), 5);
  }

  while (file2.find(".root") != std::string::npos) {
    file2.erase(file2.find(".root"), 5);
  }

  std::cout << file1 << " " << file2 << std::endl;

  TIter next(f->GetListOfKeys());
  TKey *key;
  std::string name;

  TFile *Output = new TFile((file1+"_"+file2+".root").c_str(), "recreate");

  // Loop through all entries
  while ((key = (TKey*)next())) {

    std::string ClassName = key->GetClassName();
    TClass *cl = gROOT->GetClass(ClassName.c_str());

    // Get name of object
    name = std::string(key->GetName());
    std::cout << name << std::endl;

    if (ClassName == "TH2D") {

      TH2D* minoo       = (TH2D*)(f->Get(name.c_str())->Clone());
      TH2D* neut_533    = (TH2D*)(f2->Get(name.c_str())->Clone());

      std::string MinooTitle = minoo->GetTitle();
      MinooTitle+="_Minoo";
      std::string NEUTTitle = minoo->GetTitle();
      NEUTTitle+="_NEUT533";
      minoo->SetNameTitle(MinooTitle.c_str(), MinooTitle.c_str());
      neut_533->SetNameTitle(NEUTTitle.c_str(), NEUTTitle.c_str());

      minoo->GetZaxis()->SetTitleOffset(1.5);
      neut_533->GetZaxis()->SetTitleOffset(1.5);

      Output->cd();
      minoo->Write();
      neut_533->Write();

      TH2D *Ratio = (TH2D*)(neut_533->Clone());
      Ratio->Divide(minoo);

      std::string RatioTitle = (NEUTTitle+"_vs_"+MinooTitle);
      Ratio->SetNameTitle(RatioTitle.c_str(), RatioTitle.c_str());
      Ratio->Write();

      // Now make some 1D projections
      TH1D *RatioPmu_Minoo    = minoo->ProjectionX((std::string(minoo->GetName())+"_px").c_str(), 1, minoo->GetXaxis()->GetNbins());
      TH1D *RatioCosmu_Minoo  = minoo->ProjectionY((std::string(minoo->GetName())+"_py").c_str(), 1, minoo->GetYaxis()->GetNbins());
      RatioPmu_Minoo->GetYaxis()->SetTitle("Minoo_px");
      RatioCosmu_Minoo->GetYaxis()->SetTitle("Minoo_py");

      TH1D *RatioPmu_533    = neut_533->ProjectionX((std::string(neut_533->GetName())+"_vs_"+std::string(minoo->GetName())+"_px").c_str(), 1, neut_533->GetXaxis()->GetNbins());
      TH1D *RatioCosmu_533  = neut_533->ProjectionY((std::string(neut_533->GetName())+"_vs_"+std::string(minoo->GetName())+"_py").c_str(), 1, neut_533->GetYaxis()->GetNbins());
      RatioPmu_533->GetYaxis()->SetTitle("533_px");
      RatioCosmu_533->GetYaxis()->SetTitle("533_py");

      RatioPmu_533->Divide(RatioPmu_Minoo);
      RatioPmu_533->SetTitle("NEUT/Minoo");
      RatioPmu_533->GetYaxis()->SetTitle("NEUT/Minoo");

      RatioCosmu_533->Divide(RatioCosmu_Minoo);
      RatioCosmu_533->SetTitle("NEUT/Minoo");
      RatioCosmu_533->GetYaxis()->SetTitle("NEUT/Minoo");

      RatioPmu_533->Write();
      RatioCosmu_533->Write();

      delete minoo;
      delete neut_533;
      delete RatioPmu_Minoo;
      delete RatioCosmu_Minoo;
      delete RatioPmu_533;
      delete RatioCosmu_533;
    }

    if (ClassName == "TH3D") {

      TH3D* minoo_3     = (TH3D*)(f->Get(name.c_str())->Clone());
      minoo_3->SetDirectory(0);
      TH3D* neut_533_3  = (TH3D*)(f2->Get(name.c_str())->Clone());
      neut_533_3->SetDirectory(0);

      std::string MinooTitle = minoo_3->GetTitle();
      MinooTitle+="_Minoo";
      std::string NEUTTitle = minoo_3->GetTitle();
      NEUTTitle+="_NEUT533";
      minoo_3->SetNameTitle(MinooTitle.c_str(), MinooTitle.c_str());
      neut_533_3->SetNameTitle(NEUTTitle.c_str(), NEUTTitle.c_str());

      minoo_3->GetZaxis()->SetTitleOffset(1.5);
      neut_533_3->GetZaxis()->SetTitleOffset(1.5);

      Output->cd();
      minoo_3->Write();
      neut_533_3->Write();

      TH3D *Ratio_3 = (TH3D*)(neut_533_3->Clone());
      Ratio_3->Divide(minoo_3);

      std::string RatioTitle = (NEUTTitle+"_vs_"+MinooTitle);
      Ratio_3->SetNameTitle(RatioTitle.c_str(), RatioTitle.c_str());
      Ratio_3->Write();

      delete minoo_3;
      delete neut_533_3;
    }

  } // end while

  std::cout << "Wrote to " << Output->GetName() << std::endl;
  Output->Close();

}
