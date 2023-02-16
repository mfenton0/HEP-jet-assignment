/*
 * Copyright (c) 2009--2018, the KLFitter developer team
 *
 * This file is part of KLFitter.
 *
 * KLFitter is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at
 * your option) any later version.
 *
 * KLFitter is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with KLFitter. If not, see <http://www.gnu.org/licenses/>.
 */

//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Feb  6 18:16:58 2018 by ROOT version 6.10/08
// from TTree TreeReaderTopLJets/Example input tree
// found on file: top-ljets-input.root
//////////////////////////////////////////////////////////

#ifndef TREEREADERTOPLJETS_H_
#define TREEREADERTOPLJETS_H_

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.
#include "vector"


class TreeReaderTopLJets {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
   Double_t        Lepton_pt;
   Double_t        Lepton_eta;
   Double_t        Lepton_cl_eta;
   Double_t        Lepton_phi;
   Double_t        Lepton_m;
   Float_t         Met_met;
   Float_t         sumet;
   Float_t         Met_phi;
   Char_t          Lepton_is_e;
   Char_t          Lepton_is_mu;
   Float_t         Jet_pt[18];
   Float_t         Jet_eta[18];
   Float_t         Jet_phi[18];
   Float_t         Jet_m[18];
   Float_t         Jet_btag_weight[18];
   Float_t         Jet_has_btag[18];
   Float_t         Jet_barcode[18];
   Int_t           Njet;

   // List of branches
   TBranch        *b_Lepton_pt;   //!
   TBranch        *b_Lepton_eta;   //!
   TBranch        *b_Lepton_cl_eta;   //!
   TBranch        *b_Lepton_phi;   //!
   TBranch        *b_Lepton_e;   //!
   TBranch        *b_Met_met;   //!
   TBranch        *b_Met_phi;   //!
   TBranch        *b_sumet;   //!
   TBranch        *b_Lepton_is_e;   //!
   TBranch        *b_Lepton_is_mu;   //!
   TBranch        *b_Jet_pt;   //!
   TBranch        *b_Jet_eta;   //!
   TBranch        *b_Jet_phi;   //!
   TBranch        *b_Jet_m;   //!
   TBranch        *b_Jet_btag_weight;   //!
   TBranch        *b_Jet_has_btag;   //!
   TBranch        *b_Jet_barcode;   //!
   TBranch        *b_Njet;   //!

  TreeReaderTopLJets(TTree *tree=0);
  virtual ~TreeReaderTopLJets();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     Init(TTree *tree);
};

TreeReaderTopLJets::TreeReaderTopLJets(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SimuInput_4Jet.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("SimuInput_4Jet.root");
    }
    f->GetObject("nominal",tree);

  }
  Init(tree);
}

TreeReaderTopLJets::~TreeReaderTopLJets() {
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

Int_t TreeReaderTopLJets::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void TreeReaderTopLJets::Init(TTree *tree) {
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).

  // Set object pointer
  
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
   fChain->SetBranchAddress("Lepton_pt", &Lepton_pt, &b_Lepton_pt);
   fChain->SetBranchAddress("Lepton_eta", &Lepton_eta, &b_Lepton_eta);
   fChain->SetBranchAddress("Lepton_cl_eta", &Lepton_cl_eta, &b_Lepton_cl_eta);
   fChain->SetBranchAddress("Lepton_phi", &Lepton_phi, &b_Lepton_phi);
   fChain->SetBranchAddress("Lepton_m", &Lepton_m, &b_Lepton_e);
   fChain->SetBranchAddress("Met_met", &Met_met, &b_Met_met);
   fChain->SetBranchAddress("Met_phi", &Met_phi, &b_Met_phi);
   fChain->SetBranchAddress("sumet", &sumet, &b_sumet);
   fChain->SetBranchAddress("Lepton_is_e", &Lepton_is_e, &b_Lepton_is_e);
   fChain->SetBranchAddress("Lepton_is_mu", &Lepton_is_mu, &b_Lepton_is_mu);
   fChain->SetBranchAddress("Jet_pt", Jet_pt, &b_Jet_pt);
   fChain->SetBranchAddress("Jet_eta", Jet_eta, &b_Jet_eta);
   fChain->SetBranchAddress("Jet_phi", Jet_phi, &b_Jet_phi);
   fChain->SetBranchAddress("Jet_m", Jet_m, &b_Jet_m);
   fChain->SetBranchAddress("Jet_btag_weight", Jet_btag_weight, &b_Jet_btag_weight);
   fChain->SetBranchAddress("Jet_has_btag", Jet_has_btag, &b_Jet_has_btag);
   fChain->SetBranchAddress("Jet_barcode", Jet_barcode, &b_Jet_barcode);
   fChain->SetBranchAddress("Njet", &Njet, &b_Njet);
}

#endif  // TREEREADERTOPLJETS_H_
