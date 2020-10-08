/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKPPBPYTHON_H
#define ALIANALYSISTASKPPBPYTHON_H

/**
 * @file   AliAnalysisTaskPbPbTree.h
 * @author Simone Ragoni <Simone.Ragoni@cern.ch>
 * @date   February 2019
 */

#include "AliAnalysisTaskSE.h"
#include "TBits.h"

class AliMuonTrackCuts; 					// Include class for standard muon tack cuts
class TTree;

/**
 * \file AliAnalysisTaskPbPbTree.h
 * \brief Contains the declaration of the AliAnalysisTaskPbPbTree class
 */

/**
 * \class AliAnalysisTaskPbPbTree
 * \brief Contains the tools to every part of my analysis
 */
class AliAnalysisTaskPbPbTree : public AliAnalysisTaskSE
{
    public:
                                /**
                                 * Create a new AliAnalysisTaskPbPbTree with
                                 * default constructor. Based on my understanding
                                 * this is mostly related to input/output
                                 * processes.
                                 */
                                AliAnalysisTaskPbPbTree();

                                /**
                                 * Create a new AliAnalysisTaskPbPbTree with
                                 * a specific name. This serves to instantiate
                                 * the AliAnalysisTaskSE object built with the
                                 * constructor, that will take this specific
                                 * name.
                                 *
                                 * \param name , the name taken by the AliAnalysisTaskSE object.
                                 * \param isMC , 0 if Data, 1 if MC(look at the AddTask please).
                                 */
                                AliAnalysisTaskPbPbTree(const char *name);

                                /**
                                 * Virtual destructor for the class. It will
                                 * take care of all the particular objects to be
                                 * destroyed for the class.
                                 */
        virtual                 ~AliAnalysisTaskPbPbTree();

                                /**
                                 * The function related to the instantiation of
                                 * all the histograms and the output list.
                                 */
        virtual void            UserCreateOutputObjects();

                                /**
                                 * Everything happens here. Here, the cuts are
                                 * applied, the histograms are filled and the
                                 * J/Psi peak is manifested inside the histograms.
                                 *
                                 * \param option , actually it is not used for now...
                                 */
        virtual void            UserExec(Option_t* option);

                                /**
                                 * Called at the END of the analysis (when all
                                 * events are processed). But it is not actually
                                 * doing anything! I guess it is mostly needed
                                 * for I/O purposes and GRID interfacing...
                                 */
        virtual void            Terminate(Option_t* option);

                                /**
                                 * Implement the NotifyRun to search for the new
                                 * parameters at each new runs. Sets the run
                                 * number for the successive cuts.
                                 */
        virtual void            NotifyRun();

                                /**
                                 * Use the class as a data member. It contains
                                 * the cuts for the muon track.
                                 */
        AliMuonTrackCuts*       fMuonTrackCuts;


    private:

                                /// The input events for the analysis.
        AliAODEvent*            fAOD;               //!

                                /**
                                 * The output tree containing all the
                                 * information required for a machine
                                 * learning approach.
                                 */
        TTree*                  fTree;            //!


        Double_t                fInvariantMassDimuon; //!
        Double_t                fRapidityDimuon;      //!
        Double_t                fPtDimuon;            //!
        Bool_t                  isZNAfired;           //!
        Bool_t                  isZNCfired;           //!


        Double_t                fRapiditySingleMuon_0;      //!
        Double_t                fPtSingleMuon_0;            //!
        Double_t                fPhiSingleMuon_0;           //!
        Double_t                fRapiditySingleMuon_1;      //!
        Double_t                fPtSingleMuon_1;            //!
        Double_t                fPhiSingleMuon_1;           //!
        Double_t                fPxSingleMuon_0;            //!
        Double_t                fPySingleMuon_0;            //!
        Double_t                fPzSingleMuon_0;            //!
        Double_t                fEnergySingleMuon_0;        //!
        Double_t                fPxSingleMuon_1;            //!
        Double_t                fPySingleMuon_1;            //!
        Double_t                fPzSingleMuon_1;            //!
        Double_t                fEnergySingleMuon_1;        //!


                                /**
                                 * The output list containing all the histograms
                                 * required for the analysis. In a second time
                                 * I will probably make it so to include every
                                 * possible cut variation to better compute the
                                 * systematics.
                                 */
        TList*                  fOutputList;        //!


        //_______________________________
        // CUTS
        /*
         * The following is all the possible checks for the event selections
         * and the track selection as well. Enjoy.
         */
        Int_t                   fRunNum;        //!
        Int_t                   fTracklets;     //!
        Double_t                fLumiPerRun;    //!
        Double_t*               fEtaAndPhi;     //!

        UInt_t                  fL0inputs;      //!
      	UInt_t                  fL1inputs;      //!

      	Double_t                fZem1Energy;    //!
      	Double_t                fZem2Energy;    //!

      	Double_t                fZNCEnergy;     //!
      	Double_t                fZNAEnergy;     //!
      	Double_t                fZPCEnergy;     //!
      	Double_t                fZPAEnergy;     //!
      	Double_t                fZNATDC[4];     //!
      	Double_t                fZNCTDC[4];     //!
      	Double_t                fZPATDC[4];     //!
      	Double_t                fZPCTDC[4];     //!
      	Double_t                fZNATime;       //!
      	Double_t                fZNCTime;       //!
      	Int_t                   fV0ADecision;   //!
      	Int_t                   fV0CDecision;   //!
      	Int_t                   fADADecision;   //!
      	Int_t                   fADCDecision;   //!
        TBits                   fIR1Map;        //!
        TBits                   fIR2Map;        //!


        Bool_t                  fV0Hits[64];    //!
        Int_t                   fV0TotalNCells; //!
        //_______________________________


        //_______________________________
        // TRIGGER INPUTS for MC

        // V0 inputs
        Bool_t                  fBBFlag[64];    //!
        Bool_t                  fBGFlag[64];    //!
        UInt_t                  fBBAFlags;      //!
        UInt_t                  fBBCFlags;      //!
        UInt_t                  fBGAFlags;      //!
        UInt_t                  fBGCFlags;      //!

        // AD inputs
        Bool_t                  fBBFlagAD[16];  //!
        Bool_t                  fBGFlagAD[16];  //!
        UInt_t                  fBBAFlagsAD;    //!
        UInt_t                  fBBCFlagsAD;    //!
        UInt_t                  fBGAFlagsAD;    //!
        UInt_t                  fBGCFlagsAD;    //!
        Double_t                fCounterGeneratedLevel[60000];    //!

        // FINISHED TRIGGER INPUTS for MC
        //_______________________________


        /**
         * Not implemented yet...
         */
        AliAnalysisTaskPbPbTree(const AliAnalysisTaskPbPbTree&);

        /**
         * Not implemented yet...
         */
        AliAnalysisTaskPbPbTree& operator=(const AliAnalysisTaskPbPbTree&);

        /**
         * This is important for ROOT only. I do not remember the reason anymore.
         * If I happen to encounter it again in the future, I will make sure to
         * record it!
         */
        ClassDef(AliAnalysisTaskPbPbTree, 1);
};

#endif
