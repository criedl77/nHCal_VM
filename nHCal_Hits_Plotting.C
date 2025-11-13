#include "MyConstants.h"

void plot_nHCalRecHitsPosZ(TString strang, TH1F *nHCalRecHitsPosZ_all);
void plot_nHCalRecHitsPosXY(TString strang, TH2F *nHCalRecHitsPosXY_all);

void nHCal_Hits_Plotting(TString strang){
  
  gSystem->Exec("date");
  TString flavor = "nHCal_VM";
  
  if (!strang.IsNull()) {
    // Directory does not exist (): try to make it - this doesn't work as it should...
    gSystem->mkdir(strang.Data(), kTRUE);
  }
  cout << "Output pdf directory = type of analyzed data:\n " << strang << " .\n";

  // define and open input file:
  TString infile= TString("out.") + strang + TString("-") + flavor + TString(".root");
  TFile *ifile = TFile::Open(infile,"READ");
  
  cout << "Reading infile:\n " << infile << " .\n";

  ///////////////////////////////////////////////////////////
  //Get the histograms to plot:
        
    // Clusters
    //// all
    TH1F *nHCalClustersnHits_all = (TH1F*)ifile->Get("HCalClustersnHits_all");
    TH1F *nHCalClustersEnergy_all = (TH1F*)ifile->Get("nHCalClustersEnergy_all");
    TH2F *nHCalClustersPosXY_all = (TH2F*)ifile->Get("nHCalClustersPosXY_all");
    TH1F *nHCalClustersPosZ_all = (TH1F*)ifile->Get("nHCalClustersPosZ_all");
    TH2F *nHCalClustersPosZX_all = (TH2F*)ifile->Get("nHCalClustersPosZX_all");
    TH2F *nHCalClustersPosZY_all = (TH2F*)ifile->Get("nHCalClustersPosZY_all");
    TH3F *nHCalClustersPosXYZ_all = (TH3F*)ifile->Get("nHCalClustersPosXYZ_all");
    //// XXX muons: all-->muons
    
    //// XXX charged kaons: all-->kaons
    
    //// XXX electrons (+-): all-->electrons
    

    // Hits
    //// all
    TH1F *nHCalRecHitsPosZ_all = (TH1F*)ifile->Get("nHCalRecHitsPosZ_all");
    TH2F *nHCalRecHitsPosXY_all = (TH2F*)ifile->Get("nHCalRecHitsPosXY_all");
    TH3F *nHCalRecHitsPosXYZ_all = (TH3F*)ifile->Get("nHCalRecHitsPosXYZ_all");
    TH1F *nHCalRecHitsE_all = (TH1F*)ifile->Get("nHCalRecHitsE_all");
      
  // Create a vector to store the histogram pointers:
    std::vector<TH1F*> nHCalRecHitsE_L_all;

    for (int i = 1; i <= 10; ++i) {
          std::ostringstream name;
          name << "nHCalRecHitsE_L" << i << "_all";
          TH1F* hist = (TH1F*)ifile->Get(name.str().c_str());
          nHCalRecHitsE_L_all.push_back(hist);
      }
      
    TH2F *nHCalRecHitsE_Vs_PosZ_all = (TH2F*)ifile->Get("nHCalRecHitsE_Vs_PosZ_all");

  
  
  ///////////////////////////////////////////////////////////
  // Plot:
  
    // Hit distributions:
    plot_nHCalRecHitsPosZ(strang, nHCalRecHitsPosZ_all);
    plot_nHCalRecHitsPosXY(strang, nHCalRecHitsPosXY_all);
  
  
  ///////////////////////////////////////////////////////////
  
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of main macro
//////////////////

/////////////////
// sub-macros:

void plot_nHCalRecHitsPosZ(TString strang, TH1F *nHCalRecHitsPosZ_all){
  TString name = TString("nHCalRecHitsPosZ");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  //gStyle->SetOptStat(0); //no stats box

  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
    nHCalRecHitsPosZ_all->SetLineStyle(1);
    nHCalRecHitsPosZ_all->SetTitle(strang);
    nHCalRecHitsPosZ_all->SetLineColor(kRed);
    nHCalRecHitsPosZ_all->Draw();
    canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetHeader("Reco hits in nHCal", "C");
  //leg->AddEntry(nHCalRecHitsPosZ_all,"all","l");
  leg->Draw();

  // add vertical lines for nHCal range
  Int_t binmax = nHCalRecHitsPosZ_all->GetMaximumBin();
  Double_t y_max = 0.6*nHCalRecHitsPosZ_all->GetBinContent(binmax);
  TLine *z_min_nhcal_line= new TLine(10*z_min_nhcal,0.,10*z_min_nhcal,y_max);
  z_min_nhcal_line->SetLineColor(kBlack);
  z_min_nhcal_line->SetLineWidth(2);
  z_min_nhcal_line->SetLineStyle(kDashed);
  z_min_nhcal_line->Draw("same");
  TLine *z_max_nhcal_line= new TLine(10*z_max_nhcal,0.,10*z_max_nhcal,y_max);
  z_max_nhcal_line->SetLineColor(kBlack);
  z_max_nhcal_line->SetLineWidth(2);
  z_max_nhcal_line->SetLineStyle(kDashed);
  z_max_nhcal_line->Draw("same");
  canvas->Print(filename, "pdf");
 
} // end of plot_nHCalRecHitsPosZ()

void plot_nHCalRecHitsPosXY(TString strang, TH2F *nHCalRecHitsPosXY_all){
  TString name = TString("nHCalRecHitsPosXY");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box

  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
    nHCalRecHitsPosXY_all->SetLineStyle(1);
    nHCalRecHitsPosXY_all->SetTitle(strang);
    nHCalRecHitsPosZ_all->Draw("colz");
    canvas->Draw();

    // Create a circle (TEllipse) centered at (0,0) with nHCal radius 
    TEllipse *circle = new TEllipse(0, 0, hx_max_nhcal, hx_max_nhcal);
    circle->SetFillStyle(0);     // no fill, just outline
    circle->SetLineColor(kRed);  // red border
    circle->SetLineWidth(2);     // thicker line

    // Draw it on top of the histogram
    circle->Draw("same");
    
    
  canvas->Print(filename, "pdf");
 
} // end of plot_nHCalRecHitsPosZ()
