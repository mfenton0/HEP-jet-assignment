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

// Header file for the classes stored in the TTree if any.
#include "vector"

class TreeReaderTopLJets {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Fixed size dimensions of array or collections stored in the TTree if any.

  // Declaration of leaf types
  Float_t         lepton_pt;
  Float_t         lepton_eta;
  Float_t         lepton_cl_eta;
  Float_t         lepton_phi;
  Float_t         lepton_m;
  Long64_t        lepton_pid;
  Float_t         met_met;
  Float_t         sumet;
  Float_t         met_phi;
  Char_t          lepton_is_e;
  Char_t          lepton_is_mu;
  Long64_t        njet;
  /*std::vector<float>   *jet_pt;
  std::vector<float>   *jet_eta;
  std::vector<float>   *jet_phi;
  std::vector<float>   *jet_m;
  std::vector<float>   *jet_btag_weight;
  std::vector<char>    *jet_has_btag; */
  Float_t         jet_pt[12];
  Float_t         jet_eta[12];
  Float_t         jet_phi[12];
  Float_t         jet_m[12];
  Float_t         jet_has_btag[12];
  Float_t         jet_btag_weight[12];
  Long64_t        jet_barcode[12];

  // List of branches
  TBranch        *b_lepton_pt;   //!
  TBranch        *b_lepton_eta;   //!
  TBranch        *b_lepton_cl_eta;   //!
  TBranch        *b_lepton_phi;   //!
  TBranch        *b_lepton_m;   //!
  TBranch        *b_lepton_pid;   //!
  TBranch        *b_met_met;   //!
  TBranch        *b_sumet;   //!
  TBranch        *b_met_phi;   //!
  TBranch        *b_lepton_is_e;   //!
  TBranch        *b_lepton_is_mu;   //!
  TBranch        *b_njet;   //!
  TBranch        *b_jet_pt;   //!
  TBranch        *b_jet_eta;   //!
  TBranch        *b_jet_phi;   //!
  TBranch        *b_jet_m;   //!
  TBranch        *b_jet_btag_weight;   //!
  TBranch        *b_jet_has_btag;   //!

  TreeReaderTopLJets(TTree *tree=0);
  virtual ~TreeReaderTopLJets();
  virtual Int_t    GetEntry(Long64_t entry);
  virtual void     Init(TTree *tree);
};

TreeReaderTopLJets::TreeReaderTopLJets(TTree *tree) : fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("top-ljets-input.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("top-ljets-input.root");
    }
    f->GetObject("Event",tree);

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
  //jet_pt = 0;
  //jet_eta = 0;
  //jet_phi = 0;
  //jet_m = 0;
  //jet_btag_weight = 0;
  //jet_has_btag = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  fChain->SetBranchAddress("lepton_pt", &lepton_pt, &b_lepton_pt);
  fChain->SetBranchAddress("lepton_eta", &lepton_eta, &b_lepton_eta);
  fChain->SetBranchAddress("lepton_cl_eta", &lepton_cl_eta, &b_lepton_cl_eta);
  fChain->SetBranchAddress("lepton_phi", &lepton_phi, &b_lepton_phi);
  fChain->SetBranchAddress("lepton_m", &lepton_m, &b_lepton_m);
  fChain->SetBranchAddress("lepton_pid", &lepton_pid, &b_lepton_pid);
  fChain->SetBranchAddress("met_met", &met_met, &b_met_met);
  fChain->SetBranchAddress("sumet", &sumet, &b_sumet);
  fChain->SetBranchAddress("met_phi", &met_phi, &b_met_phi);
  fChain->SetBranchAddress("lepton_is_e", &lepton_is_e, &b_lepton_is_e);
  fChain->SetBranchAddress("lepton_is_mu", &lepton_is_mu, &b_lepton_is_mu);
  fChain->SetBranchAddress("njet", &njet, &b_njet);
  fChain->SetBranchAddress("jet_pt", jet_pt, &b_jet_pt);
  fChain->SetBranchAddress("jet_eta", jet_eta, &b_jet_eta);
  fChain->SetBranchAddress("jet_phi", jet_phi, &b_jet_phi);
  fChain->SetBranchAddress("jet_m", jet_m, &b_jet_m);
  fChain->SetBranchAddress("jet_btag_weight", jet_btag_weight, &b_jet_btag_weight);
  fChain->SetBranchAddress("jet_has_btag", jet_has_btag, &b_jet_has_btag);
}

#endif  // TREEREADERTOPLJETS_H_
