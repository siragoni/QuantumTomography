/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// c++ headers
#include <iostream>
#include <fstream>
// #include <vector>
// #include <algorithm>


// root headers
#include <TMath.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH1.h"
#include "THnSparse.h"
#include <TFile.h>
#include "TTree.h"
#include <TF2.h>
#include <TF1.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph2D.h>
#include <TStopwatch.h>
#include <TMatrixDSym.h>
#include <TFitResult.h>
#include <TLatex.h>
#include "TClonesArray.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TList.h"
#include "TChain.h"


// aliroot headers
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliMuonTrackCuts.h"
#include "AliAODVertex.h"


// my headers
#include "AliAnalysisTaskPbPbTree.h"


// headers for MC
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliAODMCParticle.h"



class AliAnalysisTaskPbPbTree;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskPbPbTree) // classimp: necessary for root

//_____________________________________________________________________________
AliAnalysisTaskPbPbTree::AliAnalysisTaskPbPbTree()
    : AliAnalysisTaskSE(),
      fAOD(0),
      fTree(0),
      fOutputList(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fL0inputs(0),
      fL1inputs(0),
      fRapiditySingleMuon_0(-999.),
      fPtSingleMuon_0(-999.),
      fPhiSingleMuon_0(-999.),
      fRapiditySingleMuon_1(-999.),
      fPtSingleMuon_1(-999.),
      fPhiSingleMuon_1(-999.),
      fPxSingleMuon_0(-999.),
      fPySingleMuon_0(-999.),
      fPzSingleMuon_0(-999.),
      fEnergySingleMuon_0(-999.),
      fPxSingleMuon_1(-999.),
      fPySingleMuon_1(-999.),
      fPzSingleMuon_1(-999.),
      fEnergySingleMuon_1(-999.),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      isZNAfired(0),
      isZNCfired(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0)
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
}

//_____________________________________________________________________________
AliAnalysisTaskPbPbTree::AliAnalysisTaskPbPbTree( const char* name )
    : AliAnalysisTaskSE(name),
      fAOD(0),
      fTree(0),
      fOutputList(0),
      fMuonTrackCuts(0x0),
      fRunNum(0),
      fTracklets(0),
      fL0inputs(0),
      fL1inputs(0),
      fRapiditySingleMuon_0(-999.),
      fPtSingleMuon_0(-999.),
      fPhiSingleMuon_0(-999.),
      fRapiditySingleMuon_1(-999.),
      fPtSingleMuon_1(-999.),
      fPhiSingleMuon_1(-999.),
      fPxSingleMuon_0(-999.),
      fPySingleMuon_0(-999.),
      fPzSingleMuon_0(-999.),
      fEnergySingleMuon_0(-999.),
      fPxSingleMuon_1(-999.),
      fPySingleMuon_1(-999.),
      fPzSingleMuon_1(-999.),
      fEnergySingleMuon_1(-999.),
      fZem1Energy(0),
      fZem2Energy(0),
      fZNCEnergy(0),
      fZNAEnergy(0),
      fZPCEnergy(0),
      fZPAEnergy(0),
      fZNATime(0),
      fZNCTime(0),
      isZNAfired(0),
      isZNCfired(0),
      fV0ADecision(-10),
      fV0CDecision(-10),
      fADADecision(-10),
      fADCDecision(-10),
      fIR1Map(0),
      fIR2Map(0),
      fZNATDC{0, 0, 0, 0},
      fZNCTDC{0, 0, 0, 0},
      fZPATDC{0, 0, 0, 0},
      fZPCTDC{0, 0, 0, 0},
      fV0Hits{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fV0TotalNCells(0),
      fBBFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlag{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
               0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlags(0),
      fBBCFlags(0),
      fBGAFlags(0),
      fBGCFlags(0),
      fBBFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBGFlagAD{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      fBBAFlagsAD(0),
      fBBCFlagsAD(0),
      fBGAFlagsAD(0),
      fBGCFlagsAD(0)
{
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it,
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
                                        // make changes to your AddTask macro!
    DefineOutput(2, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskPbPbTree::~AliAnalysisTaskPbPbTree()
{
    // destructor
    if(fOutputList)    {delete fOutputList;}     	// at the end of your task, it is deleted
    if(fTree)          {delete fTree;}     	// at the end of your task, it is deleted
    if(fMuonTrackCuts) {delete fMuonTrackCuts;}   // from memory by calling this function
}
//_____________________________________________________________________________
void AliAnalysisTaskPbPbTree::UserCreateOutputObjects()
{
  // create output objects
  //
  // this function is called ONCE at the start of your analysis (RUNTIME)
  // here you ceate the histograms that you want to use
  //
  // the histograms are in this case added to a tlist, this list is in the end saved
  // to an output file
  //

  //muon track cuts
  fMuonTrackCuts = new AliMuonTrackCuts("StdMuonCuts", "StdMuonCuts");
  fMuonTrackCuts->SetFilterMask(    // AliMuonTrackCuts::kMuEta     |
                                    // AliMuonTrackCuts::kMuThetaAbs|
                                    AliMuonTrackCuts::kMuPdca    |
                                    AliMuonTrackCuts::kMuMatchLpt   );
  fMuonTrackCuts->SetAllowDefaultParams(kTRUE);
  fMuonTrackCuts->Print("mask");

  fTree = new TTree("Data","Data");
	//define branches
	// fTree->Branch( "fInvariantMassDimuon", &fInvariantMassDimuon, "fInvariantMassDimuon/D");
  // fTree->Branch( "fRapidityDimuon",      &fRapidityDimuon,      "fRapidityDimuon/D");
  // fTree->Branch( "fPtDimuon",            &fPtDimuon,            "fPtDimuon/D");
  // fTree->Branch( "fTracklets",           &fTracklets,           "fTracklets/D");
  // fTree->Branch( "fV0TotalNCells",       &fV0TotalNCells,       "fV0TotalNCells/I");
  // fTree->Branch( "fZNCEnergy",           &fZNCEnergy,           "fZNCEnergy/D");
  // fTree->Branch( "fZNAEnergy",           &fZNAEnergy,           "fZNAEnergy/D");
  // fTree->Branch( "fRunNum",              &fRunNum,              "fRunNum/I");
  // fTree->Branch( "fV0ADecision",         &fV0ADecision,         "fV0ADecision/I");
  // fTree->Branch( "fV0CDecision",         &fV0CDecision,         "fV0CDecision/I");
  // fTree->Branch( "fADADecision",         &fADADecision,         "fADADecision/I");
  // fTree->Branch( "fADCDecision",         &fADCDecision,         "fADCDecision/I");
  // fTree->Branch( "isZNAfired",           &isZNAfired,           "isZNAfired/I");
  // fTree->Branch( "isZNCfired",           &isZNCfired,           "isZNCfired/I");
  fTree->Branch( "fRapiditySingleMuon_0",   &fRapiditySingleMuon_0,"fRapiditySingleMuon_0/D");
  fTree->Branch( "fPtSingleMuon_0",         &fPtSingleMuon_0,      "fPtSingleMuon_0/D");
  fTree->Branch( "fPhiSingleMuon_0",        &fPhiSingleMuon_0,     "fPhiSingleMuon_0/D");
  fTree->Branch( "fRapiditySingleMuon_1",   &fRapiditySingleMuon_1,"fRapiditySingleMuon_1/D");
  fTree->Branch( "fPtSingleMuon_1",         &fPtSingleMuon_1,      "fPtSingleMuon_1/D");
  fTree->Branch( "fPhiSingleMuon_1",        &fPhiSingleMuon_1,     "fPhiSingleMuon_1/D");
  fTree->Branch( "fPxSingleMuon_0",         &fPxSingleMuon_0,      "fPxSingleMuon_0/D");
  fTree->Branch( "fPySingleMuon_0",         &fPySingleMuon_0,      "fPySingleMuon_0/D");
  fTree->Branch( "fPzSingleMuon_0",         &fPzSingleMuon_0,      "fPzSingleMuon_0/D");
  fTree->Branch( "fEnergySingleMuon_0",     &fEnergySingleMuon_0,  "fEnergySingleMuon_0/D");
  fTree->Branch( "fPxSingleMuon_1",         &fPxSingleMuon_1,      "fPxSingleMuon_1/D");
  fTree->Branch( "fPySingleMuon_1",         &fPySingleMuon_1,      "fPySingleMuon_1/D");
  fTree->Branch( "fPzSingleMuon_1",         &fPzSingleMuon_1,      "fPzSingleMuon_1/D");
  fTree->Branch( "fEnergySingleMuon_1",     &fEnergySingleMuon_1,  "fEnergySingleMuon_1/D");


  fOutputList = new TList();          // this is a list which will contain all
                                      // of your histograms at the end of the
                                      // analysis, the contents of this list
                                      // are written to the output file

  fOutputList->SetOwner(kTRUE);       // memory management: the list is owner
                                      // of all objects it contains and will
                                      // delete them if requested




  //_______________________________
  // - End of the function
  PostData(1, fOutputList);
  PostData(2, fTree);
}
//_____________________________________________________________________________
void AliAnalysisTaskPbPbTree::NotifyRun()
{
  /// Set run number for cuts
  fMuonTrackCuts->SetRun(fInputHandler);
}
//_____________________________________________________________________________
void AliAnalysisTaskPbPbTree::UserExec(Option_t *)
{
  /* - This iSelectionCounter is used as a token. So at every passing step it is
     - increased by one. All the events are supposed to pass the first step
     - obviously, but only a few get to the end. This is effect is clearly
     - noticeable in fCounterH event with the small trial local sample.
     - Almost 160k possible events at the 0-th step, while only 2k at the 4th step.
   */
  Int_t iSelectionCounter = 0; // no selection applied yet
  // fCounterH->Fill(iSelectionCounter); // entering UserExec
  iSelectionCounter++;


  // get AOD event
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) {
      PostData(1, fOutputList);
      PostData(2, fTree);
      return;
  }
  TString trigger = fAOD->GetFiredTriggerClasses();
  if (    !(
            trigger.Contains("CMUP10")   ||
            trigger.Contains("CMUP11")   ||
            trigger.Contains("CMUP13")   ||
            trigger.Contains("CMUP6")    ||
            trigger.Contains("CMUP26")

            )
          )  {
                    PostData(1, fOutputList);
                    PostData(2, fTree);
                    return;
  }
  // fCounterH->Fill(iSelectionCounter); // right trigger found
  iSelectionCounter++;







  fRapiditySingleMuon_0 = -999.;
  fRapiditySingleMuon_1 = -999.;
  fPtSingleMuon_0       = -999.;
  fPtSingleMuon_1       = -999.;
  fPhiSingleMuon_0      = -999.;
  fPhiSingleMuon_1      = -999.;
  fPxSingleMuon_0       = -999.;
  fPxSingleMuon_1       = -999.;
  fPySingleMuon_0       = -999.;
  fPySingleMuon_1       = -999.;
  fPzSingleMuon_0       = -999.;
  fPzSingleMuon_1       = -999.;
  fEnergySingleMuon_0   = -999.;
  fEnergySingleMuon_1   = -999.;



  /* - We are now checking if there were any tracks. If there were at least one,
     - then the histogram gets filled again. If not we are returning. There
     - would be no point in going further.
   */
  Int_t nTracks(fAOD->GetNumberOfTracks());
  if(nTracks<1) {
        PostData(1, fOutputList);
        PostData(2, fTree);
        return;
  }
  // fCounterH->Fill(iSelectionCounter); // At least one track
  iSelectionCounter++;


  //_______________________________
  // EVENT DATA EXTRACTION
  /* - Eugeny Krishen's event data extraction. I am trying to implement it.
     - The only thing I am a bit worried about is whether it should go before or
     - after the "nTracks<1" check... I will try and switch it if it sounds
     - better. These data are used for the event selection and maybe later on
     - for track selection, but I did not get to that part yet. If after all of
     - this I remember to do so, I will come back to this point and correct this
     - statement. If you find this part, please, keep in mind to check the
     - following.
   */

  /* - Event information:
     - Run Number, maybe to select the GOOD Runs and discard the others;
     - Number of Tracklets, these are in this case the SPD tracklets, so the
     - almost unit vector roughly 2 cm between two pixels of the SPD in different
     - layers.
   */
  fRunNum    = fAOD->GetRunNumber();
  fTracklets = fAOD->GetTracklets()->GetNumberOfTracklets();

  /* - Trigger Inputs:
     - L0: ..... ;
     - L1: ..... .
   */
  // fL0inputs = fAOD->GetHeader()->GetL0TriggerInputs();
  // fL1inputs = fAOD->GetHeader()->GetL1TriggerInputs();

  /* - Past-future protection maps:
     - IR1: .... ;
     - IR2: .... .
   */
  // fIR1Map = fAOD->GetHeader()->GetIRInt1InteractionMap();
  // fIR2Map = fAOD->GetHeader()->GetIRInt2InteractionMap();

  /* - ZDC: we try to find the ZDC object data in the nano-AOD. If we cannot,
     - we return, because there would be no way to actually select the events
     - otherwise! We are here, so we could even check if there is a discrepancy
     - between good events with and without ZDC's information. Or at least, this
     - is my impression of it (filling fCounterH). ZDC information:
     - fZem1Energy:
     - fZem2Energy:
     - fZNAEnergy:
     - fZNCEnergy:
     - fZPAEnergy:
     - fZPCEnergy:
     - fZNATime:
     - fZNCTime:
     - fZNATDC[i]:
     - fZNCTDC[i]:
     - fZPATDC[i]:
     - fZPCTDC[i]:
   */
  AliAODZDC *dataZDC = dynamic_cast<AliAODZDC*>(fAOD->GetZDCData());
  if(!dataZDC) {
        PostData(1, fOutputList);
        PostData(2, fTree);
        return;
  }
  // fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fZem1Energy = dataZDC->GetZEM1Energy();
  fZem2Energy = dataZDC->GetZEM2Energy();
  fZNAEnergy  = dataZDC->GetZNATowerEnergy()[0];
  fZNCEnergy  = dataZDC->GetZNCTowerEnergy()[0];
  fZPAEnergy  = dataZDC->GetZPATowerEnergy()[0];
  fZPCEnergy  = dataZDC->GetZPCTowerEnergy()[0];

  fZNATime    = dataZDC->GetZNATime();
  fZNCTime    = dataZDC->GetZNCTime();

  for (Int_t i=0;i<4;i++) fZNATDC[i] = dataZDC->GetZNATDCm(i);
  for (Int_t i=0;i<4;i++) fZNCTDC[i] = dataZDC->GetZNCTDCm(i);
  for (Int_t i=0;i<4;i++) fZPATDC[i] = dataZDC->GetZPATDCm(i);
  for (Int_t i=0;i<4;i++) fZPCTDC[i] = dataZDC->GetZPCTDCm(i);


  // V0
  AliVVZERO *dataVZERO = dynamic_cast<AliVVZERO*>(fAOD->GetVZEROData());
  if(!dataVZERO) {
        PostData(1, fOutputList);
        PostData(2, fTree);
        return;
  }
  // fCounterH->Fill(iSelectionCounter);
  iSelectionCounter++;

  fV0ADecision = dataVZERO->GetV0ADecision();
  fV0CDecision = dataVZERO->GetV0CDecision();

  // Reset event info
  fBBCFlags = 0;
  fBGCFlags = 0;
  fBBAFlags = 0;
  fBGAFlags = 0;
  for (Int_t i=0; i<64; i++) {
    // get array of fired cells
    fBBFlag[i] = dataVZERO->GetBBFlag(i);
    fBGFlag[i] = dataVZERO->GetBGFlag(i);
  }

  for (Int_t i=0; i<32; i++){ // loop over cells
    fBBCFlags+=fBBFlag[i];
    fBGCFlags+=fBGFlag[i];
    fBBAFlags+=fBBFlag[i+32];
    fBGAFlags+=fBGFlag[i+32];
  }


  //_____________________________________
  // RUN SELECTION

  //_____________________________________
  // RUN SELECTION
  /* - NOTE: total run selection.
   * -
   */
  // fCounterH->Fill(15);
  Int_t listOfGoodRunNumbersLHC18q[] = { 295585, 295586, 295587, 295588, 295589, 295612,
                                         295615, 295665, 295666, 295667, 295668, 295671,
                                         295673, 295675, 295676, 295677, 295714, 295716,
                                         295717, 295718, 295719, 295723, 295725, 295753,
                                         295754, 295755, 295758, 295759, 295762, 295763,
                                         295786, 295788, 295791, 295816, 295818, 295819,
                                         295822, 295825, 295826, 295829, 295831, 295854,
                                         295855, 295856, 295859, 295860, 295861, 295863,
                                         295881, 295908, 295909, 295910, 295913, 295936,
                                         295937, 295941, 295942, 295943, 295945, 295947,
                                         296061, 296062, 296063, 296065, 296066, 296068,
                                         296123, 296128, 296132, 296133, 296134, 296135,
                                         296142, 296143, 296191, 296192, 296194, 296195,
                                         296196, 296197, 296198, 296241, 296242, 296243,
                                         296244, 296246, 296247, 296269, 296270, 296273,
                                         296279, 296280, 296303, 296304, 296307, 296309,
                                         296312, /*296376,*/ 296377, 296378, 296379, 296380,
                                         296381, 296383, 296414, 296419, 296420, 296423,
                                         296424, 296433, 296472, 296509, 296510, 296511,
                                         296514, 296516, 296547, 296548, 296549, 296550,
                                         296551, 296552, 296553, 296615, 296616, 296618,
                                         296619, 296622, 296623 };
  Int_t listOfGoodRunNumbersLHC18r[] = { 296690, 296691, 296694, 296749, 296750, 296781,
                                         296784, 296785, 296786, 296787, 296791, 296793,
                                         296794, 296799, 296836, 296838, 296839, 296848,
                                         296849, 296850, 296851, 296852, 296890, 296894,
                                         296899, 296900, 296903, 296930, 296931, 296932,
                                         296934, 296935, 296938, 296941, 296966, 296967,
                                         296968, 296969, 296971, 296975, 296976, /*296977,*/
                                         296979, 297029, 297031, 297035, 297085, 297117,
                                         297118, 297119, 297123, 297124, 297128, 297129,
                                         297132, 297133, 297193, 297194, 297196, 297218,
                                         297219, 297221, 297222, 297278, 297310, 297312,
                                         297315, 297317, 297363, 297366, 297367, 297372,
                                         297379, 297380, 297405, 297408, 297413, 297414,
                                         297415, 297441, 297442, 297446, 297450, 297451,
                                         297452, 297479, 297481, 297483, 297512, 297537,
                                         297540, 297541, 297542, 297544, 297558, 297588,
                                         297590, 297595/*, 297623, 297624*/ };
  /* - This good run number list has been taken from the analysis
     - note of Kay's talk for DIS 2017, see:
     - https://alice-notes.web.cern.ch/system/files/notes/analysis/596/2017-Feb-08-analysis_note-2017-Feb-08-analysis-note.pdf
     -
   */
  Int_t listOfGoodRunNumbersLHC15o[] = { /*244918,*/ 244980, 244982, 244983, 245064, 245066, 245068, 245145, 245146, 245151,
                                         245152, 245231, 245232, 245233, 245253, 245259, 245343, 245345, 245346, 245347,
                                         245353, 245401, 245407, 245409, 245410, 245446, 245450, 245496, 245501, 245504,
                                         245505, 245507, 245535, 245540, 245542, 245543, 245554, 245683, 245692, 245700,
                                         245705, 245729, 245731, 245738, 245752, 245759, 245766, 245775, 245785, 245793,
                                         245829, 245831, 245833, 245949, 245952, 245954, 245963, 245996, 246001, 246003,
                                         246012, 246036, 246037, 246042, 246048, 246049, 246053, 246087, 246089, 246113,
                                         246115, 246148, 246151, 246152, 246153, 246178, 246181, 246182, 246217, 246220,
                                         246222, 246225, 246272, 246275, 246276, 246390, 246391, 246392, 246424, 246428,
                                         246431, 246433, 246434, 246487, 246488, 246493, 246495, 246675, 246676, 246750,
                                         246751, 246755, 246757, 246758, 246759, 246760, 246763, 246765, 246804, 246805,
                                         246806, 246807, 246808, 246809, 246844, 246845, 246846, 246847, 246851, 246855,
                                         246859, 246864, 246865, 246867, 246871, 246930, 246937, 246942, 246945, 246948,
                                         246949, 246980, 246982, 246984, 246989, 246991, 246994
                                       };
  Bool_t checkIfGoodRun = kFALSE;
  for( Int_t iRunLHC18q = 0; iRunLHC18q < 128; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 129; iRunLHC18q++){
  // for( Int_t iRunLHC18q = 0; iRunLHC18q < 125; iRunLHC18q++){
    if( fRunNum == listOfGoodRunNumbersLHC18q[iRunLHC18q] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC18r = 0; iRunLHC18r <  97; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  98; iRunLHC18r++){
  // for( Int_t iRunLHC18r = 0; iRunLHC18r <  82; iRunLHC18r++){
    if( fRunNum == listOfGoodRunNumbersLHC18r[iRunLHC18r] ) checkIfGoodRun = kTRUE;
  }
  for( Int_t iRunLHC15o = 0; iRunLHC15o < 136/*137*/; iRunLHC15o++){
  // for( Int_t iRunLHC15o = 0; iRunLHC15o < 134; iRunLHC15o++){
    if( fRunNum == listOfGoodRunNumbersLHC15o[iRunLHC15o] ) checkIfGoodRun = kTRUE;
  }
  if(checkIfGoodRun != 1) {
       PostData(1, fOutputList);
       PostData(2, fTree);
       // cout << "OPS!" << endl;
       return;
  }
  // fCounterH->Fill(17);

  // END RUN SELECTION
  //_____________________________________



  /* - We have to get the number of fired V0C cells. So firstly, we get the
     - boolean information about the hit cells for all V0. This is done through
     - the GetBBFlag(i) method, where 0<i<32 stands for the V0C cells and
     - 32<i<64 for the V0A cells. Then I thought the easiest way to check
     - whether the number of fired V0C cells is above 2 is just to add up the
     - boolean numbers for 0<i<32. Let's see.
   */
  fV0TotalNCells = 0;
  for(Int_t iV0Hits = 0; iV0Hits < 64; iV0Hits++) {
        fV0Hits[iV0Hits] = dataVZERO->GetBBFlag(iV0Hits);
        if(fV0Hits[iV0Hits] == kTRUE) {
              if(iV0Hits < 32) fV0TotalNCells += 1;
        }
  }

  // AD
  AliVAD *dataAD = dynamic_cast<AliVAD*>(fAOD->GetADData());
  if(dataAD) {
        // fCounterH->Fill(iSelectionCounter);
        iSelectionCounter++;

        fADADecision = dataAD->GetADADecision();
        fADCDecision = dataAD->GetADCDecision();

        // Reset event info
        fBBCFlagsAD = 0;
        fBGCFlagsAD = 0;
        fBBAFlagsAD = 0;
        fBGAFlagsAD = 0;
        for(Int_t i=0; i<16; i++) {
          // get array of fired pads
          fBBFlagAD[i] = dataAD->GetBBFlag(i);
          fBGFlagAD[i] = dataAD->GetBGFlag(i);
        }

        for(Int_t i=0; i<4; i++) { // loop over pairs of pads
          if ( fBBFlagAD[i]   && fBBFlagAD[i+4]  ) fBBCFlagsAD++;
          if ( fBGFlagAD[i]   && fBGFlagAD[i+4]  ) fBGCFlagsAD++;
          if ( fBBFlagAD[i+8] && fBBFlagAD[i+12] ) fBBAFlagsAD++;
          if ( fBGFlagAD[i+8] && fBGFlagAD[i+12] ) fBGAFlagsAD++;
        }
  }
  // END EVENT DATA EXTRACTION
  //_______________________________
  // EVENT SELECTION
  /* - This is Eugeny Krishen's event selection from the talk in 14/1/2019 for
     - the PWG-UD (UPC oriented) meeting. The event selection requires:
     - CMUP11-B triggers;
     - Maximum 2 V0C cells fired;
     - Empty V0A decision;
     - Empty ADA decision;
     - Empty ADC decision;
     - 0 tracklets in SPD;
     - Exactly 2 unlike-sign muons;
   */
  /* - Empty V0A decision
     - Empty ADA decision
     - Empty ADC decision
   */
  if(fV0ADecision != 0) {
       PostData(1, fOutputList);
       PostData(2, fTree);
       return;
  }
  if(fADADecision != 0) {
       PostData(1, fOutputList);
       PostData(2, fTree);
       return;
  }
  if(fADCDecision != 0) {
       PostData(1, fOutputList);
       PostData(2, fTree);
       return;
  }
  /* - Empty V0C decision
   * - or at least in beam timing.
   */
  if( !(fV0CDecision == 0 || fV0CDecision == 1) ) {
       PostData(1, fOutputList);
       PostData(2, fTree);
       return;
  }
  /* - 0 tracklets in SPD
   */
  // if(fTracklets != 0) {
  //      PostData(1, fOutputList);
  //      PostData(2, fTree);
  //      return;
  // }
  /* - Maximum 2 V0C cells fired.
   */
  if( fV0TotalNCells > 2 ) {
       PostData(1, fOutputList);
       PostData(2, fTree);
       return;
  }

  /* - We are finally at the starting point. We loop over the tracks and select
     - the good muons. Later on everything should happen in this loop. Let us
     - see what the future has in hold.
     -
     - Saturday: I moved the creation of the AliAODTrack* track outside of the
     - loop as it would have been otherwise created for each single iteration.
     - This could have caused massive memory issues especially to grid. I have
     - added a second AliAODTrack* track[2] to hold the second supposed muon.
     - Now this is ready to send the information to two TLorentzVectors to
     - obtain the invariant mass of the J/Psi through the Mag() method of the
     - class. Hope for the best.
   */
  Int_t nGoodMuons = 0;
  AliAODTrack* track[2];
  track[0]         = 0x0;
  track[1]         = 0x0;
  for(Int_t iTrack(0); iTrack < nTracks; iTrack++) {
    /* - This should be another form of event selection.
       - I am basically requesting the presence of TWO good muons only.
       - Later I will be checking whether of they are likesign or unlikesign.
     */
    if(nGoodMuons > 2) {
         PostData(1, fOutputList);
         PostData(2, fTree);
         return;
    }
    track[nGoodMuons] = static_cast<AliAODTrack*>(fAOD->GetTrack(iTrack));
    if(!track[nGoodMuons]) return;

    // is it a good muon track?
    if(!track[nGoodMuons]->IsMuonTrack()) {
        // track[nGoodMuons] = 0x0;
        continue;
    }
    if(!fMuonTrackCuts->IsSelected(track[nGoodMuons])) {
        // track[nGoodMuons] = 0x0;
        continue;
    }

    // MUON SELECTION
    /* - This is Eugeny Krishen's MUON selection from the talk in 14/1/2019 for
       - the PWG-UD (UPC oriented) meeting. The event selection requires:
       - Muon trigger matching >=2 (1 GeV/c threshold);
       - (-4) < eta < (-2.5);
       - (17.5 cm) < R_{abs} < (89.5 cm);
       - p \times DCA cut;
    */

    // increase counter
    nGoodMuons++;
  }
  /* - We need EXACTLY 2 good muons !!!!!
     -
   */
  if( nGoodMuons != 2 ) {
        PostData(1, fOutputList);
        PostData(2, fTree);
        return;
  }
  /* - Implementing the track cut on the unlike muons
   * -
   */
  if( (track[0]->Charge()) == (track[1]->Charge()) ) {
        PostData(1, fOutputList);
        PostData(2, fTree);
        return;
  }



  fRapiditySingleMuon_0 = track[0]->Eta();
  fRapiditySingleMuon_1 = track[1]->Eta();
  fPtSingleMuon_0       = track[0]->Pt();
  fPtSingleMuon_1       = track[1]->Pt();
  fPhiSingleMuon_0      = track[0]->Phi();
  fPhiSingleMuon_1      = track[1]->Phi();
  fPxSingleMuon_0       = track[0]->Px();
  fPxSingleMuon_1       = track[1]->Px();
  fPySingleMuon_0       = track[0]->Py();
  fPySingleMuon_1       = track[1]->Py();
  fPzSingleMuon_0       = track[0]->Pz();
  fPzSingleMuon_1       = track[1]->Pz();
  fEnergySingleMuon_0   = track[0]->E( TDatabasePDG::Instance()->GetParticle(13)->Mass() );
  fEnergySingleMuon_1   = track[1]->E( TDatabasePDG::Instance()->GetParticle(13)->Mass() );




  /* - Finally the core!!!
   * - What will be happening is that we will instantiate TLorentzVectors to
   * - obtain the invariant mass of the dimuon system. If everything goes fine
   * - after this we should be able to obtain the peak of the J/Psi. But
   * - things never go as expected, so who knows!
   */
  TLorentzVector muons[2];
  TLorentzVector possibleJPsi;
  Double_t       chargeOfMuons[2];
  for(int indexMuon = 0; indexMuon < 2; indexMuon++) {
        muons[indexMuon].SetPtEtaPhiM(   track[indexMuon]->Pt(),
                                         track[indexMuon]->Eta(),
                                         track[indexMuon]->Phi(),
                                         TDatabasePDG::Instance()->GetParticle(13)->Mass()
                                       );
        possibleJPsi += muons[indexMuon];
        chargeOfMuons[indexMuon] = track[indexMuon]->Charge();
  }




  fInvariantMassDimuon = possibleJPsi.Mag();
  fRapidityDimuon      = possibleJPsi.Rapidity();
  fPtDimuon            = possibleJPsi.Pt();


  isZNAfired = kFALSE;
  isZNCfired = kFALSE;
  Bool_t isZNAfiredStrict = kFALSE;
  Bool_t isZNCfiredStrict = kFALSE;
  Int_t  counterZNA = 0;
  Int_t  counterZNC = 0;
  /* - Note that in C++ the && and || operators "short-circuit". That means that
     - they only evaluate a parameter if required. If the first parameter to &&
     - is false, or the first to || is true, the rest will not be evaluated.
     - That means that writing:
     - if ( (isZNAfired == 0) && (...) )
     - should mean effectively
     - if ( isZNAfired != 0 ) continue;
     - hence it should be *at least* one hit!!!
     -
   */
  for(Int_t iZDC = 0; iZDC < 4 ; iZDC++) {
    if ( (isZNAfired == 0) && (fZNATDC[iZDC] > -2.) && (fZNATDC[iZDC] < 2.) ) {
      isZNAfired = kTRUE;
      /* - After mail with Chiara Oppedisano, it seems like the best way
         - to proceed is to firstly call the IsZNAfired() and then filling...
         -
         - If this doesn't appear in later pulls it is because this
         - doesn't seem to suit my case...
         -
       */
      // if( dataZDC->IsZNAfired() ) {
      //   if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
      //     fZNATimeAgainstEntriesH->Fill(fZNATDC[iZDC]);
      //   }
      // }
      // fCounterZNAH->Fill(counterZNA);
    }
    if ( (isZNCfired == 0) && (fZNCTDC[iZDC] > -2.) && (fZNCTDC[iZDC] < 2.) ) {
      isZNCfired = kTRUE;
      // if( dataZDC->IsZNCfired() ) {
      //   if ( (possibleJPsi.Mag() > 2.85) && (possibleJPsi.Mag() < 3.35) ){
      //     fZNCTimeAgainstEntriesH->Fill(fZNCTDC[iZDC]);
      //   }
      // }
      // fCounterZNCH->Fill(counterZNC);
    }
    counterZNA++;
    counterZNC++;
  }







  fTree->Fill();







  // post the data
  PostData(1, fOutputList);
  PostData(2, fTree);
}
//_____________________________________________________________________________
void AliAnalysisTaskPbPbTree::Terminate(Option_t *)
{
    cout << endl;
    // terminate
    // called at the END of the analysis (when all events are processed)
}
