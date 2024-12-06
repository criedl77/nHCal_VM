#include "MyConstants.h"

void plot_trueEta_species(TString strang, TH1F *electronEta, TH1F *muonEta, TH1F *protonEta, TH1F *pionEta, TH1F *kaonEta, TH1F *rho0Eta, TH1F *jpsiEta, TH1F *phiEta);
void plot_recEta_species(TString strang, TH1F *electronRecEta, TH1F *muonRecEta, TH1F *protonRecEta, TH1F *pionRecEta, TH1F *kaonRecEta);
void plot_genVSrecEta_species(TString strang, TH1F *electronEta, TH1F *electronRecEta, TH1F *muonEta, TH1F *muonRecEta, TH1F *pionEta, TH1F *pionRecEta, TH1F *protonEta, TH1F *protonRecEta, TH1F *kaonEta, TH1F *kaonRecEta);

void plot_Eta_decay_rho0_pipi( TString strang, TH1F *rho0Eta, TH1F *pipmfromrho0RecEta);
void plot_Eta_decay_phi_KK(TString strang, TH1F *phiEta, TH1F *kpmfromphiRecEta);
void plot_Eta_decay_jpsi_ee(TString strang, TH1F *jpsiEta, TH1F *epmfromjpsiRecEta);

void plot_kpmfromphi_momentum(TString strang, TH1F *kpmfromphiRecMom, TH1F *kpmfromphiRecMom_nHCal);
void plot_kpmfromphi_decaylength(TString strang, TH1F *kpmfromphiRecDecayLength, TH1F *kpmfromphiRecDecayLength_nHCal);
void plot_kpmfromphi_zdecay(TString strang, TH1F *kpmfromphiRecZdecay, TH1F *kpmfromphiRecZdecay_nHCal);

void nHCal_VM_Plotting(){
  
  gSystem->Exec("date");
  TString flavor = "nHCal_VM"; 
  
  // Define name of input file (= output file of nHCal_VM_Analysis.C):
  //TString strang = "pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_5.hepmc3.tree"; // has no tree "events"
  //TString strang = "sartre_bnonsat_Au_phi_ab_eAu_1.3998.eicrecon.tree.edm4eic";
  //TString strang = "EpIC1.0.0-1.0_DVMP_10x100_1.0065.eicrecon.tree.edm4eic";
  //TString strang = "rho_10x100_uChannel_Q2of0to10_hiDiv.0047.eicrecon.tree.edm4eic";
  //TString strang = "pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run39.ab.0606.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_5x41_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.0503.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_2.0256.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_3.0294.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0998.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_3.0294.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0015.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_1.0019.eicrecon.tree.edm4eic";
  //TString strang = "pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run39.ab.0606.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_5x41_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.0503.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_2.0256.eicrecon.tree.edm4eic";      
  //TString strang = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0998.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_3.0294.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0015.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_1.0019.eicrecon.tree.edm4eic";
  //TString strang = "pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1";
  //TString strang = "pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run39_10runs";
  //TString strang = "pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run39";
  //TString strang = "rho_10x100_uChannel_Q2of0to10_hiDiv";
  //TString strang = "sartre_bnonsat_Au_phi_ab_eAu_q2_15_1_1000runs";
  //TString strang = "pythia8CCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1_1000runs";
  //TString strang = "podio_output_100events";
  //TString strang = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1_100files";   

  //TString strang = "Sartre_Au_phi_10runs";

  //TString strang = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1_10files";
  //TString strang = "pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1_1000runs";
  //TString strang = "pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run39_10runs";

  // Link instead of having to change here all the time: (doesn't work conveniently like this - then I overwrite all my plots in the current directory...)
  TString strang = "current";
  
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
  
  TH1F *electronEta = (TH1F*)ifile->Get("electronEta");
  TH1F *muonEta = (TH1F*)ifile->Get("muonEta");
  TH1F *protonEta = (TH1F*)ifile->Get("protonEta");
  TH1F *pionEta = (TH1F*)ifile->Get("pionEta");
  TH1F *kaonEta = (TH1F*)ifile->Get("kaonEta");
  TH1F *rho0Eta = (TH1F*)ifile->Get("rho0Eta");
  TH1F *phiEta = (TH1F*)ifile->Get("phiEta");
  TH1F *jpsiEta = (TH1F*)ifile->Get("jpsiEta");

  TH1F *electronRecEta = (TH1F*)ifile->Get("electronRecEta");
  TH1F *muonRecEta = (TH1F*)ifile->Get("muonRecEta");
  TH1F *protonRecEta = (TH1F*)ifile->Get("protonRecEta");
  TH1F *pionRecEta = (TH1F*)ifile->Get("pionRecEta");
  TH1F *kaonRecEta = (TH1F*)ifile->Get("kaonRecEta");

  TH1F *pipmfromrho0Eta = (TH1F*)ifile->Get("pipmfromrho0Eta");
  TH1F *pipmfromrho0RecEta = (TH1F*)ifile->Get("pipmfromrho0RecEta");

  TH1F *kpmfromphiEta = (TH1F*)ifile->Get("kpmfromphiEta");
  TH1F *kpmfromphiRecEta = (TH1F*)ifile->Get("kpmfromphiRecEta");

  TH1F *kpmfromphiRecMom = (TH1F*)ifile->Get("kpmfromphiRecMom"); 
  TH1F *kpmfromphiRecMom_nHCal = (TH1F*)ifile->Get("kpmfromphiRecMom_nHCal"); 
  TH1F *kpmfromphiRecDecayLength = (TH1F*)ifile->Get("kpmfromphiRecDecayLength"); 
  TH1F *kpmfromphiRecDecayLength_nHCal = (TH1F*)ifile->Get("kpmfromphiRecDecayLength_nHCal");  
  TH1F *kpmfromphiRecZdecay = (TH1F*)ifile->Get("kpmfromphiRecZdecay");  
  TH1F *kpmfromphiRecZdecay_nHCal = (TH1F*)ifile->Get("kpmfromphiRecZdecay_nHCal"); 
  
  TH1F *epmfromjpsiEta = (TH1F*)ifile->Get("epmfromjpsiEta");
  TH1F *epmfromjpsiRecEta = (TH1F*)ifile->Get("epmfromjpsiRecEta");
  
  ///////////////////////////////////////////////////////////
  // Plot:
  
  //plot_trueEta_species(strang, electronEta, muonEta, protonEta, pionEta, kaonEta, rho0Eta, jpsiEta, phiEta);
  //plot_recEta_species(strang, electronRecEta, muonRecEta, protonRecEta, pionRecEta, kaonRecEta);
  //plot_genVSrecEta_species(strang, electronEta, electronRecEta, muonEta, muonRecEta, pionEta, pionRecEta, protonEta, protonRecEta, kaonEta, kaonRecEta);
  
  //plot_Eta_decay_rho0_pipi( strang, rho0Eta, pipmfromrho0RecEta);
  //plot_Eta_decay_phi_KK(strang, phiEta, kpmfromphiRecEta);
  //plot_Eta_decay_jpsi_ee(strang, jpsiEta, epmfromjpsiRecEta);
  
  plot_kpmfromphi_momentum(strang, kpmfromphiRecMom, kpmfromphiRecMom_nHCal);
  plot_kpmfromphi_decaylength(strang, kpmfromphiRecDecayLength, kpmfromphiRecDecayLength_nHCal);
  plot_kpmfromphi_zdecay(strang, kpmfromphiRecZdecay, kpmfromphiRecZdecay_nHCal);
  
  ///////////////////////////////////////////////////////////
  
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of main macro  
//////////////////

/////////////////
// sub-macros:
void plot_trueEta_species(TString strang, TH1F *electronEta, TH1F *muonEta, TH1F *protonEta, TH1F *pionEta, TH1F *kaonEta, TH1F *rho0Eta, TH1F *jpsiEta, TH1F *phiEta){

  TString name = TString("trueEta_species");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  electronEta->SetTitle(strang);
  electronEta->SetLineColor(kBlue);
  electronEta->Draw();
  muonEta->SetLineColor(kMagenta);
  muonEta->Draw("same");
  protonEta->SetLineColor(kRed);
  protonEta->Draw("same");
  pionEta->SetLineColor(kOrange+1);
  pionEta->Draw("same");
  kaonEta->SetLineColor(kCyan);
  kaonEta->Draw("same");
  rho0Eta->SetLineColor(kBlack);
  rho0Eta->Draw("same");
  jpsiEta->SetLineColor(kGreen+3);
  jpsiEta->Draw("same");
  phiEta->SetLineColor(kCyan+3);
  phiEta->Draw("same");
  canvas->Draw();

  auto leg = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header  
  leg->SetHeader("Generated particles", "C"); // option "C" allows to center the header
  leg->SetFillStyle(0);
  leg->AddEntry(electronEta,"electrons (#pm)","l");
  leg->AddEntry(muonEta,"muons (#pm)","l");
  leg->AddEntry(pionEta,"pions (#pm)","l");
  leg->AddEntry(kaonEta,"kaons (#pm)","l");
  leg->AddEntry(protonEta,"protons (#pm)","l");
  leg->AddEntry(rho0Eta,"#rho^{0} (770)","l");
  leg->AddEntry(jpsiEta,"J/#psi","l");
  leg->AddEntry(phiEta,"#phi(1020)","l");
  leg->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax = electronEta->GetMaximumBin();
  Double_t y_max = electronEta->GetBinContent(binmax);
  TLine *eta_min_nhcal_line= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max);  // (x1,y1,x2,y2)
  eta_min_nhcal_line->SetLineColor(kBlack);
  eta_min_nhcal_line->SetLineWidth(2);
  eta_min_nhcal_line->SetLineStyle(kDashed);
  eta_min_nhcal_line->Draw("same");
  TLine *eta_max_nhcal_line= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max);  // (x1,y1,x2,y2)
  eta_max_nhcal_line->SetLineColor(kBlack);
  eta_max_nhcal_line->SetLineWidth(2);
  eta_max_nhcal_line->SetLineStyle(kDashed);
  eta_max_nhcal_line->Draw("same");
  canvas->Print(filename, "pdf");          
  
}// end of plot_trueEta_species

void plot_recEta_species(TString strang, TH1F *electronRecEta, TH1F *muonRecEta, TH1F *protonRecEta, TH1F *pionRecEta, TH1F *kaonRecEta){

  TString name = TString("recEta_species");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box    

  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  electronRecEta->SetTitle(strang);
  electronRecEta->SetLineColor(kBlue);
  electronRecEta->Draw();
  muonRecEta->SetLineColor(kMagenta);
  muonRecEta->Draw("same");
  protonRecEta->SetLineColor(kRed);
  protonRecEta->Draw("same");
  pionRecEta->SetLineColor(kOrange+1);
  pionRecEta->Draw("same");
  kaonRecEta->SetLineColor(kCyan);
  kaonRecEta->Draw("same");
  canvas->Draw();

  auto leg = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header 
  leg->SetHeader("Reconstructed particles", "C"); // option "C" allows to center the header
  leg->SetFillStyle(0);
  leg->AddEntry(electronRecEta,"electrons (#pm)","l");
  leg->AddEntry(muonRecEta,"muons (#pm)","l");
  leg->AddEntry(pionRecEta,"pions (#pm)","l");
  leg->AddEntry(kaonRecEta,"kaons (#pm)","l");
  leg->AddEntry(protonRecEta,"protons (#pm)","l");
  leg->Draw();

  // add vertical lines for nHCal acceptance:                                                                                       
  Int_t binmax = electronRecEta->GetMaximumBin();
  Double_t y_max = electronRecEta->GetBinContent(binmax);
  TLine *eta_min_nhcal_line= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max);  // (x1,y1,x2,y2)
  eta_min_nhcal_line->SetLineColor(kBlack);
  eta_min_nhcal_line->SetLineWidth(2);
  eta_min_nhcal_line->SetLineStyle(kDashed);
  eta_min_nhcal_line->Draw("same");
  TLine *eta_max_nhcal_line= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max);  // (x1,y1,x2,y2)
  eta_max_nhcal_line->SetLineColor(kBlack);
  eta_max_nhcal_line->SetLineWidth(2);
  eta_max_nhcal_line->SetLineStyle(kDashed);
  eta_max_nhcal_line->Draw("same");
  canvas->Print(filename, "pdf");
  
} // end of plot_recEta_species()

void plot_genVSrecEta_species(TString strang, TH1F *electronEta, TH1F *electronRecEta, TH1F *muonEta, TH1F *muonRecEta, TH1F *pionEta, TH1F *pionRecEta, TH1F *protonEta, TH1F *protonRecEta, TH1F *kaonEta, TH1F *kaonRecEta){

  TString name = TString("gen-recEta_species");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");
  TCanvas *canvas = new TCanvas(name, strang, 1200, 600);

  //devide the canvas into several pads (x,y):
  canvas->Divide(3,2);
  gStyle->SetOptStat(1111); //now we want some stat boxes  
  
  // go to the first one and do your thing
  canvas->cd(1);
  electronEta->SetTitle(strang);
  electronEta->SetLineColor(kBlue);
  electronEta->Draw();
  electronRecEta->SetLineStyle(2);
  electronRecEta->Draw("same");

  auto leg1 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
  leg1->SetBorderSize(0);
  leg1->SetTextSize(0.05);
  leg1->SetFillStyle(0);
  leg1->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg1->AddEntry(electronEta,"electrons (#pm) gen","l");
  leg1->AddEntry(electronRecEta,"electrons (#pm) rec","l");
  leg1->Draw();

  // add vertical lines for nHCal acceptance 
  Int_t binmax1 = electronEta->GetMaximumBin();
  Double_t y_max1 = electronEta->GetBinContent(binmax1);
  TLine *eta_min_nhcal_line1= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max1);  // (x1,y1,x2,y2)
  eta_min_nhcal_line1->SetLineColor(kBlack);
  eta_min_nhcal_line1->SetLineWidth(2);
  eta_min_nhcal_line1->SetLineStyle(kDashed);
  eta_min_nhcal_line1->Draw("same");
  TLine *eta_max_nhcal_line1= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max1);  // (x1,y1,x2,y2)
  eta_max_nhcal_line1->SetLineColor(kBlack);
  eta_max_nhcal_line1->SetLineWidth(2);
  eta_max_nhcal_line1->SetLineStyle(kDashed);
  eta_max_nhcal_line1->Draw("same");
  
  // pad 2
  canvas->cd(2);
  muonEta->SetTitle(strang);
  muonEta->SetLineColor(kMagenta);
  muonEta->Draw();
  muonRecEta->SetLineStyle(2);
  muonRecEta->Draw("same");

  auto leg2 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header    
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.05);
  leg2->SetFillStyle(0);
  leg2->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg2->AddEntry(muonEta,"muons (#pm) gen","l");
  leg2->AddEntry(muonRecEta,"muons (#pm) rec","l");
  leg2->Draw();

// add vertical lines for nHCal acceptance
  Int_t binmax_2 = muonEta->GetMaximumBin();
  Double_t y_max2 = muonEta->GetBinContent(binmax_2);
  TLine *eta_min_nhcal_line2= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max2);
  eta_min_nhcal_line2->SetLineColor(kBlack);
  eta_min_nhcal_line2->SetLineWidth(2);
  eta_min_nhcal_line2->SetLineStyle(kDashed);
  eta_min_nhcal_line2->Draw("same");
  TLine *eta_max_nhcal_line2= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max2);
  eta_max_nhcal_line2->SetLineColor(kBlack);
  eta_max_nhcal_line2->SetLineWidth(2);
  eta_max_nhcal_line2->SetLineStyle(kDashed);
  eta_max_nhcal_line2->Draw("same");

  // pad 3
  canvas->cd(3);
  pionEta->SetTitle(strang);
  pionEta->SetLineColor(kOrange+1);
  pionEta->Draw();
  pionRecEta->SetLineStyle(2);
  pionRecEta->Draw("same");

  auto leg3 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
  leg3->SetBorderSize(0);
  leg3->SetTextSize(0.05);
  leg3->SetFillStyle(0);
  leg3->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header  
  leg3->AddEntry(pionEta,"pions (#pm) gen","l");
  leg3->AddEntry(pionRecEta,"pions (#pm) rec","l");
  leg3->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_3 = pionEta->GetMaximumBin();
  Double_t y_max3 = pionEta->GetBinContent(binmax_3);
  TLine *eta_min_nhcal_line3= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max3);  
  eta_min_nhcal_line3->SetLineColor(kBlack);
  eta_min_nhcal_line3->SetLineWidth(2);
  eta_min_nhcal_line3->SetLineStyle(kDashed);
  eta_min_nhcal_line3->Draw("same");
  TLine *eta_max_nhcal_line3= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max3);  
  eta_max_nhcal_line3->SetLineColor(kBlack);
  eta_max_nhcal_line3->SetLineWidth(2);
  eta_max_nhcal_line3->SetLineStyle(kDashed);
  eta_max_nhcal_line3->Draw("same");
  
  // pad 4
  canvas->cd(4);
  protonEta->SetTitle(strang);
  protonEta->SetLineColor(kRed);
  protonEta->Draw();
  protonRecEta->SetLineStyle(2);
  protonRecEta->Draw("same");

  auto leg4 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header					      
  leg4->SetBorderSize(0);
  leg4->SetTextSize(0.05);
  leg4->SetFillStyle(0);
  leg4->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg4->AddEntry(protonEta,"protons (#pm) gen","l");
  leg4->AddEntry(protonRecEta,"protons (#pm) rec","l");
  leg4->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_4 = protonEta->GetMaximumBin();
  Double_t y_max4 = protonEta->GetBinContent(binmax_4);
  TLine *eta_min_nhcal_line4= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max4);
  eta_min_nhcal_line4->SetLineColor(kBlack);
  eta_min_nhcal_line4->SetLineWidth(2);
  eta_min_nhcal_line4->SetLineStyle(kDashed);
  eta_min_nhcal_line4->Draw("same");
  TLine *eta_max_nhcal_line4= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max4);
  eta_max_nhcal_line4->SetLineColor(kBlack);
  eta_max_nhcal_line4->SetLineWidth(2);
  eta_max_nhcal_line4->SetLineStyle(kDashed);
  eta_max_nhcal_line4->Draw("same");

  // pad 5
  canvas->cd(5);
  kaonEta->SetTitle(strang);
  kaonEta->SetLineColor(kCyan);
  kaonEta->Draw();
  kaonRecEta->SetLineStyle(2);
  kaonRecEta->Draw("same");

  auto leg5 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header						  
  leg5->SetBorderSize(0);
  leg5->SetTextSize(0.05);
  leg5->SetFillStyle(0);
  leg5->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg5->AddEntry(kaonEta,"kaons (#pm) gen","l");
  leg5->AddEntry(kaonRecEta,"kaons (#pm) rec","l");
  leg5->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_5 = kaonEta->GetMaximumBin();
  Double_t y_max5 = kaonEta->GetBinContent(binmax_5);
  TLine *eta_min_nhcal_line5= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max5);
  eta_min_nhcal_line5->SetLineColor(kBlack);
  eta_min_nhcal_line5->SetLineWidth(2);
  eta_min_nhcal_line5->SetLineStyle(kDashed);
  eta_min_nhcal_line5->Draw("same");
  TLine *eta_max_nhcal_line5= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max5);
  eta_max_nhcal_line5->SetLineColor(kBlack);
  eta_max_nhcal_line5->SetLineWidth(2);
  eta_max_nhcal_line5->SetLineStyle(kDashed);
  eta_max_nhcal_line5->Draw("same");
  
  canvas->Draw();
  canvas->Print(filename, "pdf");
     
} // end of plot_genVSrecEta_species()

void plot_Eta_decay_rho0_pipi( TString strang, TH1F *rho0Eta, TH1F *pipmfromrho0RecEta){

  TString name = TString("Eta_decay_rho0");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  rho0Eta->SetLineColor(kBlack);
  rho0Eta->SetLineStyle(1);
  rho0Eta->Draw();
  pipmfromrho0RecEta->SetLineStyle(2);
  pipmfromrho0RecEta->SetTitle(strang);
  pipmfromrho0RecEta->SetLineColor(kRed);
  pipmfromrho0RecEta->Draw("same"); 
  canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetHeader("generated #rho^{0}(770) and decay pions", "C"); 
  leg->AddEntry(pipmfromrho0RecEta,"reco pions (#pm)","l");
  leg->AddEntry(rho0Eta,"gen #rho^{0}","l");
  leg->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax = pipmfromrho0RecEta->GetMaximumBin();
  Double_t y_max = 0.6*pipmfromrho0RecEta->GetBinContent(binmax);
  TLine *eta_min_nhcal_line= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max); 
  eta_min_nhcal_line->SetLineColor(kBlack);
  eta_min_nhcal_line->SetLineWidth(2);
  eta_min_nhcal_line->SetLineStyle(kDashed);
  eta_min_nhcal_line->Draw("same");
  TLine *eta_max_nhcal_line= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max);
  eta_max_nhcal_line->SetLineColor(kBlack);
  eta_max_nhcal_line->SetLineWidth(2);
  eta_max_nhcal_line->SetLineStyle(kDashed);
  eta_max_nhcal_line->Draw("same"); 
  canvas->Print(filename, "pdf");
  
}// end of plot_Eta_decay_rho0_pipi()


void plot_Eta_decay_phi_KK(TString strang, TH1F *phiEta, TH1F *kpmfromphiRecEta){

  TString name = TString("Eta_decay_phi");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  // meager channels:
  //phiEta->SetLineColor(kBlack);
  //phiEta->SetLineStyle(1);
  //phiEta->Draw();
  //kpmfromphiRecEta->SetLineStyle(2);
  //kpmfromphiRecEta->SetTitle(strang);
  //kpmfromphiRecEta->SetLineColor(kRed);
  //kpmfromphiRecEta->Draw("same");
  // fat channels:
  kpmfromphiRecEta->SetLineStyle(2);
  kpmfromphiRecEta->SetTitle(strang);
  kpmfromphiRecEta->SetLineColor(kRed);
  kpmfromphiRecEta->Draw();
  phiEta->SetLineStyle(1);
  phiEta->SetLineColor(kBlack);
  phiEta->SetLineStyle(1);
  phiEta->Draw("same");
  canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetHeader("generated #phi(1020) and decay kaons", "C"); 
  leg->AddEntry(kpmfromphiRecEta,"reco kaons (#pm)","l");
  leg->AddEntry(phiEta,"gen #phi","l");
  leg->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax = kpmfromphiRecEta->GetMaximumBin();
  Double_t y_max = 0.6*kpmfromphiRecEta->GetBinContent(binmax);
  TLine *eta_min_nhcal_line= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max);
  eta_min_nhcal_line->SetLineColor(kBlack);
  eta_min_nhcal_line->SetLineWidth(2);
  eta_min_nhcal_line->SetLineStyle(kDashed);
  eta_min_nhcal_line->Draw("same");
  TLine *eta_max_nhcal_line= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max);
  eta_max_nhcal_line->SetLineColor(kBlack);
  eta_max_nhcal_line->SetLineWidth(2);
  eta_max_nhcal_line->SetLineStyle(kDashed);
  eta_max_nhcal_line->Draw("same");   
  canvas->Print(filename, "pdf");
  
} // end of plot_Eta_decay_phi_KK()

void plot_Eta_decay_jpsi_ee(TString strang, TH1F *jpsiEta, TH1F *epmfromjpsiRecEta){
                                                                                                         
  TString name = TString("Eta_decay_jpsi");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  // fat channels:
  epmfromjpsiRecEta->SetLineStyle(2);
  epmfromjpsiRecEta->SetTitle(strang);
  epmfromjpsiRecEta->SetLineColor(kRed);
  epmfromjpsiRecEta->Draw();
  jpsiEta->SetLineStyle(1);
  jpsiEta->SetLineColor(kBlack);
  jpsiEta->SetLineStyle(1);
  jpsiEta->Draw("same");
  canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->SetHeader("generated J/#psi and decay electrons", "C"); 
  leg->AddEntry(epmfromjpsiRecEta,"reco electrons (#pm)","l");
  leg->AddEntry(jpsiEta,"gen J/#psi","l");
  leg->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax = epmfromjpsiRecEta->GetMaximumBin();
  Double_t y_max = 0.6*epmfromjpsiRecEta->GetBinContent(binmax);
  TLine *eta_min_nhcal_line= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max);
  eta_min_nhcal_line->SetLineColor(kBlack);
  eta_min_nhcal_line->SetLineWidth(2);
  eta_min_nhcal_line->SetLineStyle(kDashed);
  eta_min_nhcal_line->Draw("same");
  TLine *eta_max_nhcal_line= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max);
  eta_max_nhcal_line->SetLineColor(kBlack);
  eta_max_nhcal_line->SetLineWidth(2);
  eta_max_nhcal_line->SetLineStyle(kDashed);
  eta_max_nhcal_line->Draw("same");       
  canvas->Print(filename, "pdf");
 
} // end of plot_Eta_decay_jpsi_ee()

void plot_kpmfromphi_momentum(TString strang, TH1F *kpmfromphiRecMom, TH1F *kpmfromphiRecMom_nHCal){
  
  TString name = TString("kpmfromphi_momentum");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  kpmfromphiRecMom->SetTitle(strang);
  kpmfromphiRecMom->SetLineColor(kBlack);
  kpmfromphiRecMom->Draw();
  kpmfromphiRecMom_nHCal->SetLineColor(kRed);
  kpmfromphiRecMom_nHCal->Draw("same");  
  canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header  
  leg->SetHeader("Kaons from #phi(1020) decay - momentum", "C"); // option "C" allows to center the header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(kpmfromphiRecMom,"all","l");
  leg->AddEntry(kpmfromphiRecMom_nHCal,"in nHCal acceptance","l");
  leg->Draw();
  canvas->Print(filename, "pdf");
  
} //end of plot_kpmfromphi_momentum()

void plot_kpmfromphi_decaylength(TString strang, TH1F *kpmfromphiRecDecayLength, TH1F *kpmfromphiRecDecayLength_nHCal){
  
  TString name = TString("kpmfromphi_decaylength");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  kpmfromphiRecDecayLength->SetTitle(strang);
  kpmfromphiRecDecayLength->SetLineColor(kBlack);
  kpmfromphiRecDecayLength->Draw();
  kpmfromphiRecDecayLength_nHCal->SetLineColor(kRed);
  kpmfromphiRecDecayLength_nHCal->Draw("same");  
  canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header  
  leg->SetHeader("Kaons from #phi(1020) decay - decay length", "C"); // option "C" allows to center the header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(kpmfromphiRecDecayLength,"all","l");
  leg->AddEntry(kpmfromphiRecDecayLength_nHCal,"in nHCal acceptance","l");
  leg->Draw();
  canvas->Print(filename, "pdf");
  
} // end of plot_kpmfromphi_decaylength()

void plot_kpmfromphi_zdecay(TString strang, TH1F *kpmfromphiRecZdecay, TH1F *kpmfromphiRecZdecay_nHCal){
  
  TString name = TString("kpmfromphi_zdecay");
  TString filename = strang + TString("/") + TString(name) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas = new TCanvas(name, strang, 800, 600);
  kpmfromphiRecZdecay->SetTitle(strang);
  kpmfromphiRecZdecay->SetLineColor(kBlack);
  kpmfromphiRecZdecay->Draw();
  kpmfromphiRecZdecay_nHCal->SetLineColor(kRed);
  kpmfromphiRecZdecay_nHCal->Draw("same");  
  canvas->Draw();

  auto leg = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header  
  leg->SetHeader("Kaons from #phi(1020) decay - z-location of kaon decay", "C"); // option "C" allows to center the header
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.05);
  leg->AddEntry(kpmfromphiRecZdecay,"all","l");
  leg->AddEntry(kpmfromphiRecZdecay_nHCal,"in nHCal acceptance","l");
  leg->Draw();

  // add vertical lines for nHCal z-min and z-max:
  Int_t binmax = kpmfromphiRecZdecay->GetMaximumBin();
  Double_t y_max = 0.45*kpmfromphiRecZdecay->GetBinContent(binmax);
  TLine *z_min_nhcal_line= new TLine(z_nhcal_min,0.,z_nhcal_min,y_max);  // (x1,y1,x2,y2)
  z_min_nhcal_line->SetLineColor(kBlack);
  z_min_nhcal_line->SetLineWidth(2);
  z_min_nhcal_line->SetLineStyle(kDashed);
  z_min_nhcal_line->Draw("same");
  TLine *z_max_nhcal_line= new TLine(z_nhcal_max,0.,z_nhcal_max,y_max);  // (x1,y1,x2,y2)
  z_max_nhcal_line->SetLineColor(kBlack);
  z_max_nhcal_line->SetLineWidth(2);
  z_max_nhcal_line->SetLineStyle(kDashed);
  z_max_nhcal_line->Draw("same");
  
  canvas->Print(filename, "pdf");          

} // end of plot_kpmfromphi_zdecay()




  //TPaveText *t = new TPaveText(.05,.3,.95,.6, "NDC");                                                                                   
  //t->AddText("This line is blue"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);                                          
  //t->AddText("This line is red");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kRed);                                           
  //t->Draw();     
