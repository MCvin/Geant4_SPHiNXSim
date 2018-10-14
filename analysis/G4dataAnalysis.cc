#include <iostream>
#include <fstream>
#include <TFile.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TPaveStats.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLegend.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TVector3.h>

#include "G4data.h"

// user settings
bool vetoON = true;
float hitThreshold = 2;   // in keV, now ADC = keV
float trigThreshold = 30; // in keV, now ADC = keV
float hitThreshold_BGO = 50;   // in keV, now ADC = keV
float trigThreshold_BGO = 50; // in keV, now ADC = keV
float upperD = 300;       // in keV, now ADC = keV
bool saveHistos = false;

// global variables
float nPhE_keV[88];   // number of photo-electrons per keV per unit
float nADC_PhE[88];   // number of ADC channels per photo-electrons per unit
float detpos[88][2];  // position of the detector units (Slow, Fast and HighZ)
float Angle[88][88];  // angle between 2 detector units (0deg is x, 90deg is -y)
float weight[28][28]; // weight to deal with the plastic anisotropy

TH1F *Sn, *S1, *S2;
TH1F *Dn, *D1, *D2T, *D2P, *D2C;
TH1F *M2;

void createHistos()
{
  // **************** Nhits ****************
  Sn = new TH1F("Snhits", "G4 spectrum, all events", 500, 0, 500);
  Sn->GetXaxis()->SetTitle("Energy (keV)");
  Sn->GetYaxis()->SetTitle("Counts/s");
  Sn->GetXaxis()->CenterTitle();
  Sn->GetYaxis()->CenterTitle();
  Dn = new TH1F("Dnhits", "G4 distribution, all events", 88, -0.5, 87.5);
  Dn->GetXaxis()->SetTitle("Detector id");
  Dn->GetYaxis()->SetTitle("Counts");
  Dn->GetXaxis()->CenterTitle();
  Dn->GetYaxis()->CenterTitle();

  // **************** 1hit *****************
  S1 = (TH1F*)Sn->Clone("S1hit");
  S1->SetTitle("G4 spectrum, 1 hit events");
  D1 = (TH1F*)Dn->Clone("D1hit");
  D1->SetTitle("G4 distribution, 1hit events");

  // **************** 2hits ****************
  S2 = (TH1F*)Sn->Clone("S2hits");
  S2->SetTitle("G4 spectrum, 2 hit events");
  D2T = (TH1F*)Dn->Clone("D2hitsTotal");
  D2T->SetTitle("G4 distribution, 2hit events");
  D2P = (TH1F*)D2T->Clone("D2hitsPhoto");
  D2C = (TH1F*)D2T->Clone("D2hitsCompton");
  M2 = new TH1F("M2hits", "G4 modulation curve, 2 hit events", 12, 0, 360);
  M2->GetXaxis()->SetTitle("Azimutal Angle (deg)");
  M2->GetYaxis()->SetTitle("Counts");
  M2->GetXaxis()->CenterTitle();
  M2->GetYaxis()->CenterTitle();
}

void displayHistos()
{
  // gStyle->SetOptFit();
  // // **************** Nhits ****************
  TCanvas* cSn = new TCanvas("G4spectrumNhits","G4spectrumNhits");
  cSn->SetGrid();
  cSn->SetLogy();
  //TFile fileBS("../IRF/BackgroundP30-300keVBGO50-300keV/SnhitsB1000s.root");
  //TH1F *SnB = (TH1F*)fileBS.Get("Snhits");
  //SnB->Scale(50.0/1000.0);
  //Sn->Add(SnB,-1);
  Sn->Draw("");
  // TFile fileSn("SnhitsB1000s.root", "new");
  // Sn->Write();
  TF1* Sf = new TF1("Sf", "[0]*TMath::Power(x,[1])", 60,180);
  Sf->SetLineColor(kBlue);
  Sn->Fit("Sf", "R");
  // std::cout << "alpha " << Sf->GetParameter(1) << " Nhits " << Sn->GetEntries() << std::endl;
  // TCanvas* cDn = new TCanvas("G4distributionNhits","G4distributionNhits");
  // cDn->SetGrid();
  // cDn->SetLogy();
  // Dn->Draw("histe");

  // // **************** 1hit *****************
  TCanvas* cS1 = new TCanvas("G4spectrum1hit","G4spectrum1hit");
  cS1->SetGrid();
  cS1->SetLogy();
  S1->Draw("");
  //S1->Fit("Sf", "R");
  // TCanvas* cD1 = new TCanvas("G4distribution1hit","G4distribution1hit");
  // cD1->SetGrid();
  // cD1->SetLogy();
  // D1->Draw("histe");
  // TFile fileD1("D1hitB1000s.root", "new");
  // D1->Write();

  // // **************** 2hits ****************
  // TCanvas* cS2 = new TCanvas("G4spectrum2hits","G4spectrum2hits");
  // cS2->SetGrid();
  // cS2->SetLogy();
  // S2->Draw("");
  // //S2->Fit("Sf", "R");
  // TCanvas* cD2 = new TCanvas("G4distribution2hits","G4distribution2hits");
  // cD2->SetGrid();
  // cD2->SetLogy();
  // D2T->Draw("histe");
  // D2P->SetLineColor(kGreen);
  // D2P->Draw("sames");
  // D2C->SetLineColor(kRed);
  // D2C->Draw("sames");
  // TLegend* lD2 = new TLegend(0.45, 0.6, 0.85, 0.75);
  // lD2->AddEntry(D2T, "Total hits", "l");
  // lD2->AddEntry(D2P, "Photo hits (Emax)", "l");
  // lD2->AddEntry(D2C, "Compton hits (Emin)", "l");
  // lD2->Draw();
  // modulation curve  
  //To write the histogram of the unpolarised distribution
  // TFile fileM2("M2hitsB1000s.root", "new");
  // M2->Write();

  //Sn->Draw();

  //TCanvas* lots = new TCanvas("Sn_plots","Sn_plots");
  //TFile fileM1("Snhits_alpha_1m.root");
  //TH1F *M1 = (TH1F*)fileM1.Get("Snhits");
  //M1->Draw();
  
  
//   std::cout << "number of background 2-hits " << M2B->GetEntries() << std::endl;
  //std::cout << "number of signal 2-hits " << M2->GetEntries() << std::endl;
  //To read this histogram in another Root session, do:
  //TFile file("M2hitsNP.root");
  //TFile fileB("../IRF/BackgroundP30-300keVBGO50-300keV/M2hitsB1000s.root");
  //TH1F *M2B = (TH1F*)fileB.Get("M2hits");
  //TH1F *M2NP = (TH1F*)file.Get("M2hits");
  //M2NP->Scale(12/M2NP->Integral());
  //M2B->Scale(50.0/1000.0);
  //M2->Divide(M2NP);
  //M2->Add(M2B,-1);
  TCanvas* cM2 = new TCanvas("G4modulation2hits","G4modulation2hits");
  M2->Draw("E");
  // //One 180-degree fit only ("normal" modulation curve)
  TF1* Mf = new TF1("Mf", "[0]*(1+[1]*cos(2*3.1415/180.0*(x-[2])))", -10,370);
  Mf->SetLineColor(kBlue);
  //Combined 180-degree and 360-degree fit (useful for fitting asymmetric signal)
  //TF1* Mf = new TF1("Mf", "[0]*(1+[1]*cos(2*3.1415/180.0*(x-[2]))+[3]*cos(3.1415/180.0*(x-[4])))", -10,370);
  //M2->Fit("Mf", "RN"); // for batch mode
  M2->Fit("Mf", "R", "E");
  std::cout << "MF " << Mf->GetParameter(1) << " N2hits " << M2->GetEntries() << " PA " << Mf->GetParameter(2) << " theta " << Mf->GetParameter(3) << " phi " << Mf->GetParameter(4) << " alpha " << Sf->GetParameter(1) << std::endl;

  if (saveHistos) {
  }
}

void ReadCalibrationParameters()
{
  std::ifstream lightyield("light_yields.dat");
  std::ifstream pmtgain("PMT_gains.dat");
  std::ifstream detpositions("det_positions.dat");
  std::ifstream angmatrix("angle_matrix.dat");
  std::ifstream weightmatrix("weight_matrix.dat");
  for (int x = 0; x < 88; x++) {
    lightyield >> nPhE_keV[x];
    pmtgain >> nADC_PhE[x];
  }
  for (int i = 0; i < 88; i++) {
    detpositions >> detpos[i][0] >> detpos[i][1];
    for (int j = 0; j < 88; j++) {
      angmatrix >> Angle[i][j];
    }
  }
  for (int i = 0; i < 28; i++) {
    for (int j = 0; j < 28; j++) {
      weightmatrix >> weight[i][j];
    }
  }
  lightyield.close();
  pmtgain.close();
  detpositions.close();
  angmatrix.close();
  weightmatrix.close();
}

void EnergyConversion_keV2ADC(G4data* G4)
{
  TRandom3 random(0);   // "TRandom3 is the best random engine in ROOT. 0 means seed from system clock."

  for (int x = 0; x < 60; x++) {
    float E_fast = G4->Fast->EdepTot[x];
    float E_highZ = G4->HighZ->EdepTot[x];
    // correction for scintillator non-linearity response
      // slow scintillator
      // fast scintillator
      E_fast = E_fast * (1.001 - 0.486 * exp(-0.0902 * E_fast));
      // highZ
      E_highZ = E_highZ * (1.001 - 0.486 * exp(-0.0902 * E_highZ));

    // correction for scintillator position dependance (need to output everything from G4)
      // slow scintillator
      // fast scintillator
      // highZ

    // correction for PMT photon-electron production (low energy cutoff)
      // slow scintillator
      // fast scintillator
      int nPhE_fast = random.Poisson(E_fast * nPhE_keV[x]);
      // highZ
      int nPhE_highZ = random.Poisson(E_highZ * nPhE_keV[x+28]);

    // correction for PMT gain fluctuation
      // slow scintillator
      // fast scintillator
      int ADC_fast = random.Gaus(nPhE_fast, 0.35*sqrt(nPhE_fast)) * nADC_PhE[x];
      // highZ
      int ADC_highZ = random.Gaus(nPhE_highZ, 0.35*sqrt(nPhE_highZ)) * nADC_PhE[x+28];

    G4->Fast->EdepTot[x] = ADC_fast;
    G4->HighZ->EdepTot[x] = ADC_highZ;
  }
}

void EnergyReconstruction_ADC2keV(G4data* G4)
{
  for (int x = 0; x < 60; x++) {
    int ADC_fast = G4->Fast->EdepTot[x];
    int ADC_highZ = G4->HighZ->EdepTot[x];
    // ADC to keV
      // slow scintillator
      // fast scintillator
      float E_fast = ADC_fast / nADC_PhE[x] / nPhE_keV[x];
      // highZ
      float E_highZ = ADC_highZ / nADC_PhE[x] / nPhE_keV[x];

    // correction for scintillator non-linearity response
      // slow scintillator
      // fast scintillator
      E_fast = E_fast / (1.001 - 0.486 * exp(-0.0902 * E_fast)); // not exactly correct
      // highZ
      E_highZ = E_highZ / (1.001 - 0.486 * exp(-0.0902 * E_highZ)); // not exactly correct

    G4->Fast->EdepTot[x] = E_fast;
    G4->HighZ->EdepTot[x] = E_highZ;
  }
}

bool EventsSelection(G4data* G4)
{
  bool valid = false;
  for (int x = 0; x < 60; x++) {
    // upper discriminator
    if (G4->Fast->EdepTot[x] > upperD) return false;  // if UD the function stops and return false
    if (G4->HighZ->EdepTot[x] > upperD) return false; // if UD the function stops and return false
    // trigger threshold
    if (G4->Fast->EdepTot[x] > trigThreshold) valid = true;
    if (G4->HighZ->EdepTot[x] > trigThreshold_BGO) valid = true;
    // hit threshold
    if (G4->Fast->EdepTot[x] < hitThreshold) G4->Fast->EdepTot[x] = 0;
    if (G4->HighZ->EdepTot[x] < hitThreshold_BGO) G4->HighZ->EdepTot[x] = 0;
  }
  // update the number of hitten units after the hit threshold
  G4->Fast->nUnits = 0;
  G4->HighZ->nUnits = 0;
  for (int x = 0; x < 60; x++) {
    if (G4->Fast->EdepTot[x] > 0) G4->Fast->nUnits += 1;
    if (G4->HighZ->EdepTot[x] > 0) G4->HighZ->nUnits += 1;
  }

  return valid;
}

void fillnhits(G4data* G4)
{
  for (int x = 0; x < 60; x++) {
    // Spectra
    if (G4->Fast->EdepTot[x] > 0) Sn->Fill(G4->Fast->EdepTot[x]);
    if (G4->HighZ->EdepTot[x] > 0) Sn->Fill(G4->HighZ->EdepTot[x]);

    // Distribution
    if (G4->Fast->EdepTot[x] > 0) Dn->Fill(x);
    if (G4->HighZ->EdepTot[x] > 0) Dn->Fill(x+28);
  }
}

void fill1hit(G4data* G4)
{
  if ((G4->Fast->nUnits + G4->HighZ->nUnits) == 1) {
    for (int x = 0; x < 60; x++) {
      // Spectra
      if (G4->Fast->EdepTot[x] > 0) S1->Fill(G4->Fast->EdepTot[x]);
      if (G4->HighZ->EdepTot[x] > 0) S1->Fill(G4->HighZ->EdepTot[x]);

      // Distribution
      if (G4->Fast->EdepTot[x] > 0) D1->Fill(x);
      if (G4->HighZ->EdepTot[x] > 0) D1->Fill(x+28);
    }
  }
}

void fill2hits(G4data* G4)
{
  int h1 = -1;
  int h2 = -1;
  int htemp = -1;
  float E1 = -1;
  float E2 = -1;
  float Etemp = -1;
  if ((G4->Fast->nUnits + G4->HighZ->nUnits) == 2) {
    // Spectra
    for (int x = 0; x < 60; x++) {
      if (G4->Fast->EdepTot[x] > 0) S2->Fill(G4->Fast->EdepTot[x]);
      if (G4->HighZ->EdepTot[x] > 0) S2->Fill(G4->HighZ->EdepTot[x]);
    }
    // Distribution
      // Total energy
      for (int x = 0; x < 60; x++) {
	if (G4->Fast->EdepTot[x] > 0) {
	  D2T->Fill(x);
	  if (h1 == -1) {h1 = x; E1 = G4->Fast->EdepTot[x];}
	  else if (h2 == -1) {h2 = x; E2 = G4->Fast->EdepTot[x];}
	}
	if (G4->HighZ->EdepTot[x] > 0) {
	  D2T->Fill(x+28);
	  if (h1 == -1) {h1 = x+28; E1 = G4->HighZ->EdepTot[x];}
	  else if (h2 == -1) {h2 = x+28; E2 = G4->HighZ->EdepTot[x];}
	}
      }
      // to always have E1 < E2
      if (E2 < E1) {htemp = h1; h1 = h2; h2 = htemp; Etemp = E1; E1 = E2; E2 = Etemp;}

      // Photopeak
      int pP = TMath::LocMax(60, G4->Fast->EdepTot);
      int pHZ = TMath::LocMax(60, G4->HighZ->EdepTot);
      if (G4->Fast->EdepTot[pP] > G4->HighZ->EdepTot[pHZ]) {
	D2P->Fill(pP);
	pHZ = -1;
      }
      else {
	D2P->Fill(pHZ+28);
	pP = -1;
      }

      // Compton
      for (int x = 0; x < 60; x++) {
	if ((G4->Fast->EdepTot[x] > 0) && (x != pP)) D2C->Fill(x);
	if ((G4->HighZ->EdepTot[x] > 0) && (x != pHZ)) D2C->Fill(x+28);
      }

    // Angular distribution
      if ((h1 < 28) && (h2 < 28)) M2->Fill(Angle[h1][h2],weight[h1][h2]); // Fast to Fast, corrected
      if ((h1 < 28) && (h2 > 27)) M2->Fill(Angle[h1][h2],weight[h1][h1]); // Fast to HighZ, corrected
      if (h1 > 27) M2->Fill(Angle[h1][h2]); // HighZ to any
  }
  // if ((G4->Fast->nUnits + G4->HighZ->nUnits) == 3) {
  //   M2->Fill(0);
  // }
}

void printAngleMatrix()
{
  // compute angles between detector units
  double theta;
  double deltaX, deltaY;
  for (int i = 0; i < 88; i++) {
    for (int j = i+1; j < 88; j++) {
      deltaX = detpos[j][0]-detpos[i][0];
      deltaY = detpos[j][1]-detpos[i][1];
      theta = std::acos(deltaX/std::sqrt(deltaX*deltaX+deltaY*deltaY))*180/std::acos(-1.0);
      if(deltaY > 0){ //Vector pointing towards y
	Angle[i][j] = 360 - theta;
      }
      else{ //Vector pointing towards -y
	Angle[i][j] = theta;
      }
      if (Angle[i][j] < 180) Angle[j][i] = Angle[i][j] + 180;
      else Angle[j][i] = Angle[i][j] - 180;
    }
  }

  // print angles to a file
  std::ofstream Amatrix("angle_matrix.dat");
  Amatrix.precision(4);
  for (int i = 0; i < 88; i++) {
    for (int j = 0; j < 88; j++) {
      Amatrix << Angle[i][j] << " ";
    }
    Amatrix << std::endl;
  }
  Amatrix.close();
}

void getLocalisation()
{
  char *filename;
  TH1F *D1NP, *D1copy;
  TFile fileB("IRF/BackgroundEcut30-300keV/D1hitB1000s.root");
  TH1F *D1B = (TH1F*)fileB.Get("D1hit");
  D1B->Scale(50.0/1000.0);
  TF1 *Df = new TF1("Df", "[0]", 0,88);

  for (int i = 0; i < 91; i++) {
    D1copy = (TH1F*)D1->Clone("D1copy");
    filename = Form("IRF/Alpha-1.05_E5-1000keV_ECut30-300keV/NP0-90deg/D1hitNP%d.root",i);
    TFile fileNP(filename);
    D1NP = (TH1F*)fileNP.Get("D1hit");
    D1NP->Scale(88.0/D1NP->Integral());
    D1copy->Add(D1B,-1);
    D1copy->Divide(D1NP);
    D1copy->Fit("Df", "R", "E");
    std::cout << "OffAxis " << i  << " Chi2 " << Df->GetChisquare()/88 << std::endl;
    fileNP.Close("R");
  }
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", 0, 0); // load every root library in memory

  // open the Tree and create objects to contain the branches
  TFile *f;
  if (argc == 2) f = new TFile(argv[1],"r");
  else f = new TFile("G4data.root","r");
  TTree *tree = (TTree*) f->Get("G4data");
  G4data* G4 = new G4data; //Declare the object
  G4->SetAddress(tree);    //Set every address explicitly (IMPORTANT!)

  // read the characteristics of the PDCs from ascii files
  ReadCalibrationParameters();

  // compute and print the angle matrix to a file
  //printAngleMatrix();

  // create histograms
  createHistos();

  double sumEp = 0;
  double sumEb = 0;

  // loop over all entries in the tree
  for (long i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);

    for (int x = 0; x < 60; x++) {
      sumEp += G4->Fast->EdepTot[x];
      sumEb += G4->HighZ->EdepTot[x];
    }

    // energy conversion: keV to ADC
    //EnergyConversion_keV2ADC(G4);

    // energy reconstruction: ADC to keV
    //EnergyReconstruction_ADC2keV(G4);

    // events selection
    bool validEvent = EventsSelection(G4);

    if (validEvent) {
      // Nhits
      fillnhits(G4);
      // 1hit
      fill1hit(G4);
      // 2hits
      fill2hits(G4);
    }
  }

  std::cout << "Total energy in plastic: " << sumEp << std::endl;
  std::cout << "Total energy in BGO: " << sumEb << std::endl;

  // get GRB localisation
  //getLocalisation();

  // display histograms
  displayHistos();

  theApp.Run();
  return 0;
}
