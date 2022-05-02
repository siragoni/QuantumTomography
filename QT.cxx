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



using namespace std;            // std namespace: so you can do things like 'cout'

//_____________________________________________________________________________
/* -
   - Quantum Tomography for CS.
   -
 */
Double_t  CosThetaQuantumTomCS( TLorentzVector muonPositive,
                                TLorentzVector muonNegative,
                                TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  // TLorentzVector pProjCM(0.,0., -1, 1); // projectile
  // TLorentzVector pTargCM(0.,0.,  1, 1); // target

  TLorentzVector VectorMeson = possibleJPsi;

  // // TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) ).Unit();
  // TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) );
  // // TLorentzVector QuantumZ  =  VectorMeson;
  // // QuantumZ *= ( 1./QuantumZ.Mag2() );
  //
  // Double_t       xQuantumY =  pProjCM.Pz()*pTargCM.E()*VectorMeson.Py()-pProjCM.E()*pTargCM.Pz()*VectorMeson.Py();
  // Double_t       yQuantumY = -pProjCM.Pz()*pTargCM.E()*VectorMeson.Px()+pProjCM.E()*pTargCM.Pz()*VectorMeson.Px();
  // TLorentzVector QuantumYnotnormalised(xQuantumY, yQuantumY, 0., 0.);
  // // TLorentzVector QuantumY = QuantumYnotnormalised.Unit();
  // TLorentzVector QuantumY = QuantumYnotnormalised;
  // // QuantumY *= ( 1./QuantumY.Mag2() );
  // // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) ).Unit();
  // // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) );
  // // QuantumX *= ( 1./QuantumX.Mag2() );
  // TLorentzVector QuantumX = VectorMeson - pProjCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pProjCM)) )) - pTargCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pTargCM)) ));
  // // QuantumX *= ( 1./QuantumX.Mag2() );
  //
  // cout << "XZ dot product = " << QuantumX.Dot(QuantumZ) << endl;
  // cout << "YZ dot product = " << QuantumY.Dot(QuantumZ) << endl;
  // cout << "XY dot product = " << QuantumX.Dot(QuantumY) << endl;
  // cout << "XQ dot product = " << QuantumX.Dot(VectorMeson) << endl;
  //
  // TLorentzVector DifferenceMuons   = muonPositive - muonNegative;
  // Double_t       xFinalVector      = QuantumX.Dot(DifferenceMuons);
  // Double_t       yFinalVector      = QuantumY.Dot(DifferenceMuons);
  // Double_t       zFinalVector      = QuantumZ.Dot(DifferenceMuons);
  // TVector3       FinalVector( xFinalVector, yFinalVector, zFinalVector );
  //
  // Double_t       CosThetaQuantumCS = ( FinalVector.Unit() ).Z();
  // Double_t       SinThetaQuantumCS = TMath::Sqrt( 1 - CosThetaQuantumCS*CosThetaQuantumCS );
  // Double_t       PhiQuantumCS      = TMath::ASin( ( (FinalVector.Unit()).Y() )/ SinThetaQuantumCS );




  // Daniel



  //calculate z;
  TLorentzVector z;
  Float_t part1 = 0.;
  Float_t part2 = 0.;

  //Dot product: v1*v2 = t1*t2-x1*x2-y1*y2-z1*z2

  part1 = VectorMeson.Dot(pTargCM);
  part2 = VectorMeson.Dot(pProjCM);

  Float_t part3x = pProjCM.X()*part1;
  Float_t part3y = pProjCM.Y()*part1;
  Float_t part3z = pProjCM.Z()*part1;
  Float_t part3e = pProjCM.T()*part1;

  Float_t part4x = pTargCM.X()*part2;
  Float_t part4y = pTargCM.Y()*part2;
  Float_t part4z = pTargCM.Z()*part2;
  Float_t part4e = pTargCM.T()*part2;

  TLorentzVector part3(TVector3(part3x,part3y,part3z), part3e);
  TLorentzVector part4(TVector3(part4x,part4y,part4z),part4e);

  // Q=Q; pb=Pbar; pa=P; from paper
  // Un-normalized Z
  // z = part3 - part4;
  z = part4 - part3;

   //Normalized z
   Float_t normz = TMath::Sqrt(-z*z);
   Float_t znx = z.X()/normz;
   Float_t zny = z.Y()/normz;
   Float_t znz = z.Z()/normz;
   Float_t zne = z.E()/normz;

   //Normalized z
   TLorentzVector zhat(TVector3(znx,zny,znz),zne);

// calculate x
TLorentzVector x;

Float_t constant1 = (VectorMeson.Dot(VectorMeson))/(2*(VectorMeson.Dot(pProjCM)));
Float_t constant2 = (VectorMeson.Dot(VectorMeson))/(2*(VectorMeson.Dot(pTargCM)));

Float_t comp1x = pProjCM.X()*constant1;
Float_t comp1y = pProjCM.Y()*constant1;
Float_t comp1z = pProjCM.Z()*constant1;
Float_t comp1e = pProjCM.T()*constant1;

TLorentzVector comp1(TVector3(comp1x,comp1y,comp1z),comp1e);

Float_t comp2x = pTargCM.X()*constant2;
Float_t comp2y = pTargCM.Y()*constant2;
Float_t comp2z = pTargCM.Z()*constant2;
Float_t comp2e = pTargCM.T()*constant2;

TLorentzVector comp2(TVector3(comp2x,comp2y, comp2z),comp2e);

//Un-normalized x
x = VectorMeson - comp1 - comp2;


  //normalize x
  Float_t normx = TMath::Sqrt(-x*x);
  Float_t xnx = x.X()/normx;
  Float_t xny = x.Y()/normx;
  Float_t xnz = x.Z()/normx;
  Float_t xne = x.E()/normx;

   //Normalized x
  TLorentzVector xhat(TVector3(xnx,xny,xnz),xne);






// calculate y
//TLorentzVector y;
Float_t yone = pProjCM.Y()*pTargCM.Z()*VectorMeson.E() - pProjCM.Z()*pTargCM.Y()*VectorMeson.E() + pProjCM.Z()*pTargCM.E()*VectorMeson.Y() + pProjCM.E()*pTargCM.Y()*VectorMeson.Z() - pProjCM.Y()*pTargCM.E()*VectorMeson.Z() - pProjCM.E()*pTargCM.Z()*VectorMeson.Y();
Float_t ytwo = -pProjCM.Z()*pTargCM.E()*VectorMeson.X() + pProjCM.Z()*pTargCM.X()*VectorMeson.E() - pProjCM.X()*pTargCM.Z()*VectorMeson.E() + pProjCM.X()*pTargCM.E()*VectorMeson.Z() - pProjCM.E()*pTargCM.X()*VectorMeson.Z() + pProjCM.E()*pTargCM.Z()*VectorMeson.X();
Float_t ythree = pProjCM.X()*pTargCM.Y()*VectorMeson.E() - pProjCM.Y()*pTargCM.X()*VectorMeson.E() + pProjCM.Y()*pTargCM.E()*VectorMeson.X() - pProjCM.X()*pTargCM.E()*VectorMeson.Y() + pProjCM.E()*pTargCM.X()*VectorMeson.Y() - pProjCM.E()*pTargCM.Y()*VectorMeson.X();
Float_t yfour = -pProjCM.X()*pTargCM.Y()*VectorMeson.Z() + pProjCM.X()*pTargCM.Z()*VectorMeson.Y() - pProjCM.Z()*pTargCM.X()*VectorMeson.Y() + pProjCM.Z()*pTargCM.Y()*VectorMeson.X() - pProjCM.Y()*pTargCM.Z()*VectorMeson.X() + pProjCM.Y()*pTargCM.X()*VectorMeson.Z();

//Un-normalized y
TLorentzVector y(TVector3(yone,ytwo,ythree),yfour);

 //normalize y
 Float_t normy = TMath::Sqrt(-y*y);
 Float_t ynx = y.X()/normy;
 Float_t yny = y.Y()/normy;
 Float_t ynz = y.Z()/normy;
 Float_t yne = y.E()/normy;

//normalized y
 TLorentzVector yhat(TVector3(ynx,yny,ynz),yne);



 cout << "XZ dot product = " << x.Dot(z) << endl;
 cout << "YZ dot product = " << y.Dot(z) << endl;
 cout << "XY dot product = " << x.Dot(y) << endl;
 cout << "XQ dot product = " << x.Dot(VectorMeson) << endl;
 cout << "ZQ dot product = " << z.Dot(VectorMeson) << endl;



   //Lepton momentum difference
   TLorentzVector diff;
   diff = (muonPositive - muonNegative);
   Float_t diff2x = diff.X()/2.;
   Float_t diff2y = diff.Y()/2.;
   Float_t diff2z = diff.Z()/2.;
   Float_t diff2e = diff.E()/2.;
   TLorentzVector diff2(TVector3(diff2x,diff2y,diff2z),diff2e);

   //Normalize diff2
   Float_t norm2 = TMath::Sqrt(-diff2*diff2);
   Float_t diff3x = diff2.X()/norm2;
   Float_t diff3y = diff2.Y()/norm2;
   Float_t diff3z = diff2.Z()/norm2;
   Float_t diff3e = diff2.E()/norm2;

   TLorentzVector diff3(TVector3(diff3x,diff3y,diff3z),diff3e);

   //computing the angles
   Float_t cosThetaCS =  zhat*diff3;
   Float_t SinThetaCosPhiCS = xhat*diff3;
   Float_t SinThetaSinPhiCS = yhat*diff3;
  //**************************************








  return cosThetaCS;

}
//_____________________________________________________________________________
Double_t  PhiQuantumTomogrCS( TLorentzVector muonPositive,
                              TLorentzVector muonNegative,
                              TLorentzVector possibleJPsi )
{
  /* - This function computes the Collins-Soper cos(theta) for the
     - helicity of the J/Psi.
     - The idea should be to get back to a reference frame where it
     - is easier to compute and to define the proper z-axis.
     -
   */

  /* - Half of the energy per pair of the colliding nucleons.
     -
   */
  Double_t HalfSqrtSnn   = 2510.;
  Double_t MassOfLead208 = 193.6823;
  Double_t MomentumBeam  = TMath::Sqrt( HalfSqrtSnn*HalfSqrtSnn*208*208 - MassOfLead208*MassOfLead208 );
  /* - Fill the Lorentz vector for projectile and target.
     - For the moment we do not consider the crossing angle.
     - Projectile runs towards the MUON arm.
     -
   */
  // TLorentzVector pProjCM(0.,0., -MomentumBeam, HalfSqrtSnn*208); // projectile
  // TLorentzVector pTargCM(0.,0.,  MomentumBeam, HalfSqrtSnn*208); // target
  TLorentzVector pProjCM(0.,0., -1, 1); // projectile
  TLorentzVector pTargCM(0.,0.,  1, 1); // target

  TLorentzVector VectorMeson = possibleJPsi;

  // // TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) ).Unit();
  // TLorentzVector QuantumZ  = ( pProjCM*( VectorMeson.Dot(pTargCM) ) - pTargCM*( VectorMeson.Dot(pProjCM) ) );
  // // TLorentzVector QuantumZ  =  VectorMeson;
  // // QuantumZ *= ( 1./QuantumZ.Mag2() );
  //
  // Double_t       xQuantumY =  pProjCM.Pz()*pTargCM.E()*VectorMeson.Py()-pProjCM.E()*pTargCM.Pz()*VectorMeson.Py();
  // Double_t       yQuantumY = -pProjCM.Pz()*pTargCM.E()*VectorMeson.Px()+pProjCM.E()*pTargCM.Pz()*VectorMeson.Px();
  // TLorentzVector QuantumYnotnormalised(xQuantumY, yQuantumY, 0., 0.);
  // // TLorentzVector QuantumY = QuantumYnotnormalised.Unit();
  // TLorentzVector QuantumY = QuantumYnotnormalised;
  // // QuantumY *= ( 1./QuantumY.Mag2() );
  // // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) ).Unit();
  // // TLorentzVector QuantumX = ( QuantumY.Dot(QuantumZ) );
  // // QuantumX *= ( 1./QuantumX.Mag2() );
  // TLorentzVector QuantumX = VectorMeson - pProjCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pProjCM)) )) - pTargCM*(VectorMeson.Mag2()/( 2*(VectorMeson.Dot(pTargCM)) ));
  // // QuantumX *= ( 1./QuantumX.Mag2() );
  //
  // cout << "XZ dot product = " << QuantumX.Dot(QuantumZ) << endl;
  // cout << "YZ dot product = " << QuantumY.Dot(QuantumZ) << endl;
  // cout << "XY dot product = " << QuantumX.Dot(QuantumY) << endl;
  // cout << "XQ dot product = " << QuantumX.Dot(VectorMeson) << endl;
  //
  // TLorentzVector DifferenceMuons   = muonPositive - muonNegative;
  // Double_t       xFinalVector      = QuantumX.Dot(DifferenceMuons);
  // Double_t       yFinalVector      = QuantumY.Dot(DifferenceMuons);
  // Double_t       zFinalVector      = QuantumZ.Dot(DifferenceMuons);
  // TVector3       FinalVector( xFinalVector, yFinalVector, zFinalVector );
  //
  // Double_t       CosThetaQuantumCS = ( FinalVector.Unit() ).Z();
  // Double_t       SinThetaQuantumCS = TMath::Sqrt( 1 - CosThetaQuantumCS*CosThetaQuantumCS );
  // Double_t       PhiQuantumCS      = TMath::ASin( ( (FinalVector.Unit()).Y() )/ SinThetaQuantumCS );




  // Daniel



  //calculate z;
  TLorentzVector z;
  Float_t part1 = 0.;
  Float_t part2 = 0.;

  //Dot product: v1*v2 = t1*t2-x1*x2-y1*y2-z1*z2

  part1 = VectorMeson.Dot(pTargCM);
  part2 = VectorMeson.Dot(pProjCM);

  Float_t part3x = pProjCM.X()*part1;
  Float_t part3y = pProjCM.Y()*part1;
  Float_t part3z = pProjCM.Z()*part1;
  Float_t part3e = pProjCM.T()*part1;

  Float_t part4x = pTargCM.X()*part2;
  Float_t part4y = pTargCM.Y()*part2;
  Float_t part4z = pTargCM.Z()*part2;
  Float_t part4e = pTargCM.T()*part2;

  TLorentzVector part3(TVector3(part3x,part3y,part3z), part3e);
  TLorentzVector part4(TVector3(part4x,part4y,part4z),part4e);

  // Q=Q; pb=Pbar; pa=P; from paper
  // Un-normalized Z
  // z = part3 - part4;
  z = part4 - part3;

   //Normalized z
   Float_t normz = TMath::Sqrt(-z*z);
   Float_t znx = z.X()/normz;
   Float_t zny = z.Y()/normz;
   Float_t znz = z.Z()/normz;
   Float_t zne = z.E()/normz;

   //Normalized z
   TLorentzVector zhat(TVector3(znx,zny,znz),zne);

// calculate x
TLorentzVector x;

Float_t constant1 = (VectorMeson.Dot(VectorMeson))/(2*(VectorMeson.Dot(pProjCM)));
Float_t constant2 = (VectorMeson.Dot(VectorMeson))/(2*(VectorMeson.Dot(pTargCM)));

Float_t comp1x = pProjCM.X()*constant1;
Float_t comp1y = pProjCM.Y()*constant1;
Float_t comp1z = pProjCM.Z()*constant1;
Float_t comp1e = pProjCM.T()*constant1;

TLorentzVector comp1(TVector3(comp1x,comp1y,comp1z),comp1e);

Float_t comp2x = pTargCM.X()*constant2;
Float_t comp2y = pTargCM.Y()*constant2;
Float_t comp2z = pTargCM.Z()*constant2;
Float_t comp2e = pTargCM.T()*constant2;

TLorentzVector comp2(TVector3(comp2x,comp2y, comp2z),comp2e);

//Un-normalized x
x = VectorMeson - comp1 - comp2;


  //normalize x
  Float_t normx = TMath::Sqrt(-x*x);
  Float_t xnx = x.X()/normx;
  Float_t xny = x.Y()/normx;
  Float_t xnz = x.Z()/normx;
  Float_t xne = x.E()/normx;

   //Normalized x
  TLorentzVector xhat(TVector3(xnx,xny,xnz),xne);






// calculate y
//TLorentzVector y;
Float_t yone = pProjCM.Y()*pTargCM.Z()*VectorMeson.E() - pProjCM.Z()*pTargCM.Y()*VectorMeson.E() + pProjCM.Z()*pTargCM.E()*VectorMeson.Y() + pProjCM.E()*pTargCM.Y()*VectorMeson.Z() - pProjCM.Y()*pTargCM.E()*VectorMeson.Z() - pProjCM.E()*pTargCM.Z()*VectorMeson.Y();
Float_t ytwo = -pProjCM.Z()*pTargCM.E()*VectorMeson.X() + pProjCM.Z()*pTargCM.X()*VectorMeson.E() - pProjCM.X()*pTargCM.Z()*VectorMeson.E() + pProjCM.X()*pTargCM.E()*VectorMeson.Z() - pProjCM.E()*pTargCM.X()*VectorMeson.Z() + pProjCM.E()*pTargCM.Z()*VectorMeson.X();
Float_t ythree = pProjCM.X()*pTargCM.Y()*VectorMeson.E() - pProjCM.Y()*pTargCM.X()*VectorMeson.E() + pProjCM.Y()*pTargCM.E()*VectorMeson.X() - pProjCM.X()*pTargCM.E()*VectorMeson.Y() + pProjCM.E()*pTargCM.X()*VectorMeson.Y() - pProjCM.E()*pTargCM.Y()*VectorMeson.X();
Float_t yfour = -pProjCM.X()*pTargCM.Y()*VectorMeson.Z() + pProjCM.X()*pTargCM.Z()*VectorMeson.Y() - pProjCM.Z()*pTargCM.X()*VectorMeson.Y() + pProjCM.Z()*pTargCM.Y()*VectorMeson.X() - pProjCM.Y()*pTargCM.Z()*VectorMeson.X() + pProjCM.Y()*pTargCM.X()*VectorMeson.Z();

//Un-normalized y
TLorentzVector y(TVector3(yone,ytwo,ythree),yfour);

 //normalize y
 Float_t normy = TMath::Sqrt(-y*y);
 Float_t ynx = y.X()/normy;
 Float_t yny = y.Y()/normy;
 Float_t ynz = y.Z()/normy;
 Float_t yne = y.E()/normy;

//normalized y
 TLorentzVector yhat(TVector3(ynx,yny,ynz),yne);



 cout << "XZ dot product = " << x.Dot(z) << endl;
 cout << "YZ dot product = " << y.Dot(z) << endl;
 cout << "XY dot product = " << x.Dot(y) << endl;
 cout << "XQ dot product = " << x.Dot(VectorMeson) << endl;
 cout << "ZQ dot product = " << z.Dot(VectorMeson) << endl;



   //Lepton momentum difference
   TLorentzVector diff;
   diff = (muonPositive - muonNegative);
   Float_t diff2x = diff.X()/2.;
   Float_t diff2y = diff.Y()/2.;
   Float_t diff2z = diff.Z()/2.;
   Float_t diff2e = diff.E()/2.;
   TLorentzVector diff2(TVector3(diff2x,diff2y,diff2z),diff2e);

   //Normalize diff2
   Float_t norm2 = TMath::Sqrt(-diff2*diff2);
   Float_t diff3x = diff2.X()/norm2;
   Float_t diff3y = diff2.Y()/norm2;
   Float_t diff3z = diff2.Z()/norm2;
   Float_t diff3e = diff2.E()/norm2;

   TLorentzVector diff3(TVector3(diff3x,diff3y,diff3z),diff3e);

   //computing the angles
   Float_t cosThetaCS       = zhat*diff3;
   Float_t SinThetaCosPhiCS = xhat*diff3;
   Float_t SinThetaSinPhiCS = yhat*diff3;
  //**************************************



  Double_t SinTheta = TMath::Sqrt( 1 - cosThetaCS*cosThetaCS );
  Double_t Phi      = TMath::ASin( SinThetaSinPhiCS / SinTheta );




  // Double_t       CosThetaQuantumCS = ( FinalVector.Unit() ).Z();
  // Double_t       SinThetaQuantumCS = TMath::Sqrt( 1 - CosThetaQuantumCS*CosThetaQuantumCS );
  // Double_t       PhiQuantumCS      = TMath::ASin( ( (FinalVector.Unit()).Y() )/ SinThetaQuantumCS );


  return Phi;

}
//_____________________________________________________________________________
void ComputeIt()
{
  TLorentzVector muonPositive( 1.201056,  0.586901, 21.715403, 21.756766 );
  TLorentzVector muonNegative(-1.260339, -0.528812, -9.929302, 10.023487 );
  TLorentzVector VectorMeson = muonPositive + muonNegative;

  Double_t CosThetaQuantumCSvalue = CosThetaQuantumTomCS(muonPositive, muonNegative, VectorMeson);
  std::cout << "CosThetaQuantumCSvalue = " << CosThetaQuantumCSvalue <<'\n';





  //  ROOT USING LIGHT CONE VECTORS:
  // XZ dot product = 6.42065e-06
  // YZ dot product = 0
  // XY dot product = 3.01748e-11
  // XQ dot product = 1.96022e-06
  // ZQ dot product = -2.61186e-05
  // CosThetaQuantumCSvalue = -0.995799
  //
  //
  // PYTHON LIGHT CONES
  // X  =  (-0.0592830000000002, 0.05808900000000006, 9.320944201007819e-05, 0.0002513316022891843)
  // Y  =  (-31666347903.0162, -32317239111.2691, 0, 0)
  // Z  =  (0.0, 0.0, -63.560506, -23.572201999999997)
  // Xh =  (-0.714265893086271, 0.699880091484714, 0.0011230255780139546, 0.0030281462032941495)
  // Yh =  (-0.699877324006233, -0.714263068723194, 0, 0)
  // Zh =  (0.0, 0.0, -1.0767881005148563, -0.39934022263026814)
  // XY =  1.11022302462516e-16
  // XZ =  -1.2420620087993939e-15
  // YZ =  0
  // ZQ =  0.0
  // CosThetaQuantumCSvalue = 0.995798810540499
  //
  // PYTHON TRUE BEAM COORDINATES
  // X  =  (-0.0592830000000002, 0.05808900000000006, 9.509020021170045e-05, 0.00025202910479649177)
  // Y  =  (-31666347903.0162, -32317239111.2691, 0, 0)
  // Z  =  (0.0, 0.0, -17324528704984.98, -6425016348811.073)
  // Xh =  (-0.7142658929281116, 0.69988009132974, 0.0011456857237812342, 0.003036549998841796)
  // Yh =  (-0.699877324006233, -0.714263068723194, 0, 0)
  // Zh =  (0.0, 0.0, -1.0767881005148563, -0.3993402226302681)
  // XY =  5.55111512312578e-17
  // XZ =  2.104420173196014e-05
  // YZ =  0
  // ZQ =  1.7763568394002505e-15
  // CosThetaQuantumCSvalue = 0.995798810540499   // basically same but with inverted sign
  //
  //
  // ROOT WITH BEAM COORDINATES
  // XZ dot product = -3.74212e+07
  // YZ dot product = 0
  // XY dot product = 3.08428
  // XQ dot product = -5.0828e-05
  // ZQ dot product = -7.79453e+06
  // CosThetaQuantumCSvalue = -0.995799


}
