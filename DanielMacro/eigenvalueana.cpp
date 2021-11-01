#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLeaf.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TF1.h"
#include "TLatex.h"
using namespace std;
#include <math.h>
#include <vector>

TArrayD data(9);

//_______________________________
void input_data(Int_t i){
  for (size_t i = 0; i < 9; i++) {
    data[i] = 0.;
  }
  if(         i == 0 ){
    data[0] = 0.000787;
    data[4] = 0.316639;
    data[6] = 0.015815;
    data[8] = 0.317957;
  } else if ( i == 1 ){
    data[0] = 0.008094;
    data[4] = 0.269297;
    data[6] = 0.054889;
    data[8] = 0.372213;
  } else if ( i == 2 ){
    data[0] = 0.000815;
    data[4] = 0.349335;
    data[6] = 0.015569;
    data[8] = 0.297438;
  } else if ( i == 3 ){
    data[0] = 0.000001;
    data[4] = 0.362348;
    data[6] = 0.000675;
    data[8] = 0.349528;
  } else if ( i == 4 ){
    data[0] = 0.009978;
    data[4] = 0.374140;
    data[6] =-0.061210;
    data[8] = 0.375516;
  } else if ( i == 5 ){
    data[0] = 0.055535;
    data[4] = 0.248970;
    data[6] =-0.166230;
    data[8] = 0.497566;
  } else if ( i == 6 ){
    data[0] = 0.002244;
    data[4] = 0.313791;
    data[6] =-0.027315;
    data[8] = 0.332445;
  } else if ( i == 7 ){
    data[0] = 0.007446;
    data[4] = 0.287387;
    data[6] = 0.050184;
    data[8] = 0.338240;
  } else if ( i == 8 ){
    data[0] = 0.000048;
    data[4] = 0.320893;
    data[6] =-0.003999;
    data[8] = 0.332889;
  } else if ( i == 9 ){
    data[0] = 0.002459;
    data[4] = 0.333986;
    data[6] = 0.027424;
    data[8] = 0.305909;
  } else if ( i == 10){
    data[0] = 0.000124;
    data[4] = 0.370709;
    data[6] = 0.005669;
    data[8] = 0.259843;
  } else if ( i == 11){
    data[0] = 0.005556;
    data[4] = 0.293723;
    data[6] = 0.043639;
    data[8] = 0.342732;
  }
}
//_____________________________________________________________________________
/* - Eigenvalue ana.
 * -
 */
Double_t EigenvalueEntropy(Int_t i, Int_t returnValue){
  TMatrixD h(3,3);
  input_data(i);
  // TArrayD data(9);
  // data[0] = 0.000787;
  // data[1] = 0.;
  // data[2] = 0.;
  // data[3] = 0.;
  // data[4] = 0.316639;
  // data[5] = 0.;
  // data[6] = 0.015815;
  // data[7] = 0.;
  // data[8] = 0.317957;
  h.SetMatrixArray(data.GetArray());
  h.Print();
  TMatrixDEigen eigen(h);
  auto eigenvalues = eigen.GetEigenValues();
  eigenvalues.Print();
  auto eigenvectors = eigen.GetEigenVectors();
  eigenvectors.Print();



  auto h2 = h*h;
  h2.Print();

  cout << h2(0, 0);
  cout << h2(1, 1);

  Double_t Trace = h2(0, 0)+h2(1, 1)+h2(2, 2);

  Double_t Degree = TMath::Sqrt(1.5*Trace - 0.5);
  cout << "Degree = " << Degree;




  Double_t Entropy = -999.;
  TMatrixD logeigenvalues(3,3);
  logeigenvalues(0,0) = TMath::Log(eigenvalues(0,0));
  logeigenvalues(1,1) = TMath::Log(eigenvalues(1,1));
  logeigenvalues(2,2) = TMath::Log(eigenvalues(2,2));


  TMatrixD eigenTimesLogeigen = eigenvalues*logeigenvalues;
  Entropy = (-1.)*(eigenTimesLogeigen(0,0) + eigenTimesLogeigen(1,1) + eigenTimesLogeigen(2,2));
  cout << "Entropy = " << Entropy;


  if (       returnValue == 0){
    return Entropy;
  } else if (returnValue == 1){
    return eigenvalues(0,0);
  } else if (returnValue == 2){
    return eigenvalues(1,1);
  } else if (returnValue == 3){
    return eigenvalues(2,2);
  } else {
    return Entropy;
  }

}
//_____________________________________________________________________________
void Plot() {
  Double_t EntropyValues[12];
  Double_t EigenZero[12];
  Double_t EigenOne[12];
  Double_t EigenTwo[12];
  for (size_t i = 0; i < 12; i++) {
    EntropyValues[i] = EigenvalueEntropy(i,0);
    EigenZero[i]     = EigenvalueEntropy(i,1);
    EigenOne[i]      = EigenvalueEntropy(i,2);
    EigenTwo[i]      = EigenvalueEntropy(i,3);
  }






  const Int_t n = 6;
  Double_t x[6] = {2., 3., 4., 5., 7., 10.};
  Double_t EntropyHE[6] = {EntropyValues[0], EntropyValues[1], EntropyValues[2],
                           EntropyValues[3], EntropyValues[4], EntropyValues[5]};
  Double_t EntropyCS[6] = {EntropyValues[6], EntropyValues[7], EntropyValues[8],
                           EntropyValues[9], EntropyValues[10],EntropyValues[11]};
  Double_t Eigen0HE[6] = { EigenZero[0], EigenZero[1], EigenZero[2],
                           EigenZero[3], EigenZero[4], EigenZero[5]};
  Double_t Eigen0CS[6] = { EigenZero[6], EigenZero[7], EigenZero[8],
                           EigenZero[9], EigenZero[10],EigenZero[11]};
  Double_t Eigen1HE[6] = { EigenOne[0], EigenOne[1], EigenOne[2],
                           EigenOne[3], EigenOne[4], EigenOne[5]};
  Double_t Eigen1CS[6] = { EigenOne[6], EigenOne[7], EigenOne[8],
                           EigenOne[9], EigenOne[10],EigenOne[11]};
  Double_t Eigen2HE[6] = { EigenTwo[0], EigenTwo[1], EigenTwo[2],
                           EigenTwo[3], EigenTwo[4], EigenTwo[5]};
  Double_t Eigen2CS[6] = { EigenTwo[6], EigenTwo[7], EigenTwo[8],
                           EigenTwo[9], EigenTwo[10],EigenTwo[11]};


  TCanvas *c1 = new TCanvas("c1","Entanglement entropy",200,10,700,500);
  c1->SetGrid();
  TGraph *grEntropyHE = new TGraph(n,x,EntropyHE);
  grEntropyHE->GetHistogram()->SetMaximum(1.);   // along
  grEntropyHE->GetHistogram()->SetMinimum(0.5);  //   Y
  grEntropyHE->SetLineColor(2);
  grEntropyHE->SetLineWidth(4);
  grEntropyHE->SetMarkerColor(4);
  grEntropyHE->SetMarkerStyle(21);
  grEntropyHE->SetTitle("Entanglement entropy");
  grEntropyHE->GetXaxis()->SetTitle("p_{T}");
  grEntropyHE->GetYaxis()->SetTitle("Entanglement Entropy [a.u.]");
  grEntropyHE->Draw("ACP");
  auto leg = new TLatex(4.,0.8,"HELICITY");
  leg->SetTextSize(0.08);
  leg->SetLineWidth(2);
  leg->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c1->Update();
  c1->GetFrame()->SetBorderSize(12);
  c1->Modified();
  c1->SaveAs("plots/entropyHE.pdf");








  TCanvas *c12 = new TCanvas("c12","Entanglement entropy",200,10,700,500);
  c12->SetGrid();
  TGraph *grEntropyCS = new TGraph(n,x,EntropyCS);
  grEntropyCS->GetHistogram()->SetMaximum(1.);   // along
  grEntropyCS->GetHistogram()->SetMinimum(0.5);  //   Y
  grEntropyCS->SetLineColor(2);
  grEntropyCS->SetLineWidth(4);
  grEntropyCS->SetMarkerColor(4);
  grEntropyCS->SetMarkerStyle(21);
  grEntropyCS->SetTitle("Entanglement entropy");
  grEntropyCS->GetXaxis()->SetTitle("p_{T}");
  grEntropyCS->GetYaxis()->SetTitle("Entanglement Entropy [a.u.]");
  grEntropyCS->Draw("ACP");
  auto leg2 = new TLatex(4.,0.8,"COLLINS-SOPER");
  leg2->SetTextSize(0.08);
  leg2->SetLineWidth(2);
  leg2->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c12->Update();
  c12->GetFrame()->SetBorderSize(12);
  c12->Modified();
  c12->SaveAs("plots/entropyCS.pdf");










  TCanvas *c13 = new TCanvas("c13","Eigenvalue",200,10,700,500);
  c13->SetGrid();
  TGraph *grEigen0HE = new TGraph(n,x,Eigen0HE);
  grEigen0HE->GetHistogram()->SetMaximum(1.);   // along
  grEigen0HE->GetHistogram()->SetMinimum(0.);  //   Y
  grEigen0HE->SetLineColor(2);
  grEigen0HE->SetLineWidth(4);
  grEigen0HE->SetMarkerColor(4);
  grEigen0HE->SetMarkerStyle(21);
  grEigen0HE->SetTitle("First Eigenvalue");
  grEigen0HE->GetXaxis()->SetTitle("p_{T}");
  grEigen0HE->GetYaxis()->SetTitle("First Eigenvalue [a.u.]");
  grEigen0HE->Draw("ACP");
  auto leg3 = new TLatex(4.,0.8,"HELICITY");
  leg3->SetTextSize(0.08);
  leg3->SetLineWidth(2);
  leg3->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c13->Update();
  c13->GetFrame()->SetBorderSize(12);
  c13->Modified();
  c13->SaveAs("plots/eigen0HE.pdf");








  TCanvas *c14 = new TCanvas("c14","Eigenvalue",200,10,700,500);
  c14->SetGrid();
  TGraph *grEigen0CS = new TGraph(n,x,Eigen0CS);
  grEigen0CS->GetHistogram()->SetMaximum(1.);   // along
  grEigen0CS->GetHistogram()->SetMinimum(0.);  //   Y
  grEigen0CS->SetLineColor(2);
  grEigen0CS->SetLineWidth(4);
  grEigen0CS->SetMarkerColor(4);
  grEigen0CS->SetMarkerStyle(21);
  grEigen0CS->SetTitle("First Eigenvalue");
  grEigen0CS->GetXaxis()->SetTitle("p_{T}");
  grEigen0CS->GetYaxis()->SetTitle("First Eigenvalue [a.u.]");
  grEigen0CS->Draw("ACP");
  auto leg4 = new TLatex(4.,0.8,"COLLINS-SOPER");
  leg4->SetTextSize(0.08);
  leg4->SetLineWidth(2);
  leg4->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c14->Update();
  c14->GetFrame()->SetBorderSize(12);
  c14->Modified();
  c14->SaveAs("plots/eigen0CS.pdf");
























  TCanvas *c15 = new TCanvas("c15","Eigenvalue",200,10,700,500);
  c15->SetGrid();
  TGraph *grEigen1HE = new TGraph(n,x,Eigen1HE);
  grEigen1HE->GetHistogram()->SetMaximum(1.);   // along
  grEigen1HE->GetHistogram()->SetMinimum(0.);  //   Y
  grEigen1HE->SetLineColor(2);
  grEigen1HE->SetLineWidth(4);
  grEigen1HE->SetMarkerColor(4);
  grEigen1HE->SetMarkerStyle(21);
  grEigen1HE->SetTitle("Second Eigenvalue");
  grEigen1HE->GetXaxis()->SetTitle("p_{T}");
  grEigen1HE->GetYaxis()->SetTitle("Second Eigenvalue [a.u.]");
  grEigen1HE->Draw("ACP");
  auto leg5 = new TLatex(4.,0.8,"HELICITY");
  leg5->SetTextSize(0.08);
  leg5->SetLineWidth(2);
  leg5->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c15->Update();
  c15->GetFrame()->SetBorderSize(12);
  c15->Modified();
  c15->SaveAs("plots/eigen1HE.pdf");









  TCanvas *c16 = new TCanvas("c16","Eigenvalue",200,10,700,500);
  c16->SetGrid();
  TGraph *grEigen1CS = new TGraph(n,x,Eigen1CS);
  grEigen1CS->GetHistogram()->SetMaximum(1.);   // along
  grEigen1CS->GetHistogram()->SetMinimum(0.);  //   Y
  grEigen1CS->SetLineColor(2);
  grEigen1CS->SetLineWidth(4);
  grEigen1CS->SetMarkerColor(4);
  grEigen1CS->SetMarkerStyle(21);
  grEigen1CS->SetTitle("Second Eigenvalues");
  grEigen1CS->GetXaxis()->SetTitle("p_{T}");
  grEigen1CS->GetYaxis()->SetTitle("Second Eigenvalue [a.u.]");
  grEigen1CS->Draw("ACP");
  auto leg6 = new TLatex(4.,0.8,"COLLINS-SOPER");
  leg6->SetTextSize(0.08);
  leg6->SetLineWidth(2);
  leg6->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c16->Update();
  c16->GetFrame()->SetBorderSize(12);
  c16->Modified();
  c16->SaveAs("plots/eigen1CS.pdf");














  TCanvas *c17 = new TCanvas("c17","Eigenvalue",200,10,700,500);
  c17->SetGrid();
  TGraph *grEigen2HE = new TGraph(n,x,Eigen2HE);
  grEigen2HE->GetHistogram()->SetMaximum(1.);   // along
  grEigen2HE->GetHistogram()->SetMinimum(0.);  //   Y
  grEigen2HE->SetLineColor(2);
  grEigen2HE->SetLineWidth(4);
  grEigen2HE->SetMarkerColor(4);
  grEigen2HE->SetMarkerStyle(21);
  grEigen2HE->SetTitle("Third Eigenvalue");
  grEigen2HE->GetXaxis()->SetTitle("p_{T}");
  grEigen2HE->GetYaxis()->SetTitle("Third Eigenvalue [a.u.]");
  grEigen2HE->Draw("ACP");
  auto leg7 = new TLatex(4.,0.8,"HELICITY");
  leg7->SetTextSize(0.08);
  leg7->SetLineWidth(2);
  leg7->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c17->Update();
  c17->GetFrame()->SetBorderSize(12);
  c17->Modified();
  c17->SaveAs("plots/eigen2HE.pdf");









  TCanvas *c18 = new TCanvas("c18","Eigenvalue",200,10,700,500);
  c18->SetGrid();
  TGraph *grEigen2CS = new TGraph(n,x,Eigen2CS);
  grEigen2CS->GetHistogram()->SetMaximum(1.);   // along
  grEigen2CS->GetHistogram()->SetMinimum(0.);  //   Y
  grEigen2CS->SetLineColor(2);
  grEigen2CS->SetLineWidth(4);
  grEigen2CS->SetMarkerColor(4);
  grEigen2CS->SetMarkerStyle(21);
  grEigen2CS->SetTitle("Third Eigenvalues");
  grEigen2CS->GetXaxis()->SetTitle("p_{T}");
  grEigen2CS->GetYaxis()->SetTitle("Third Eigenvalue [a.u.]");
  grEigen2CS->Draw("ACP");
  auto leg8 = new TLatex(4.,0.8,"COLLINS-SOPER");
  leg8->SetTextSize(0.08);
  leg8->SetLineWidth(2);
  leg8->Draw("same");
  // TCanvas::Update() draws the frame, after which one can change it
  c18->Update();
  c18->GetFrame()->SetBorderSize(12);
  c18->Modified();
  c18->SaveAs("plots/eigen2CS.pdf");

}
