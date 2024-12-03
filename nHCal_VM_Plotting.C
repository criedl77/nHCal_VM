void nHCal_VM_Plotting()
{
  gSystem->Exec("date");

  //
  double eta_min_nhcal = -4.05; 
  double eta_max_nhcal = -1.2; 
  //
  double eta_min_bhcal = -1.2;
  double eta_max_bhcal = 1.18;
 //
  double eta_min_lfhcal = 1.18;			    
  double eta_max_lfhcal = 4.2;
  //
  double z_nhcal_min = -3.95; // (2024-12-03) assumed start of nHCal in z-direction [m], from $DETECTOR_PATH/compact/definitions.xml
  double z_nhcal_thickness = 0.45; // nHCal thickness in z [m]
  double z_nhcal_max = z_nhcal_min - z_nhcal_thickness;
  //

  // Define name of input file (= output file of nHCal_VM_Analysis.C):
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_5.0001.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_vtxfix_5.hepmc3.tree"); // has no tree "events"
  //TString strang_ram("sartre_bnonsat_Au_phi_ab_eAu_1.3998.eicrecon.tree.edm4eic");
  //TString strang_ram("EpIC1.0.0-1.0_DVMP_10x100_1.0065.eicrecon.tree.edm4eic");
  //TString strang_ram("rho_10x100_uChannel_Q2of0to10_hiDiv.0047.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run39.ab.0606.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_5x41_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.0503.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_2.0256.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_3.0294.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0998.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_3.0294.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0015.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_1.0019.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia_ep_noradcor_10x100_q2_0.000000001_1.0_run39.ab.0606.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_5x41_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_5.0503.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_10x100_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_2.0256.eicrecon.tree.edm4eic");                         
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1.0998.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=10_beamEffects_xAngle=-0.025_hiDiv_3.0294.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1.0015.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=1000_beamEffects_xAngle=-0.025_hiDiv_1.0019.eicrecon.tree.edm4eic");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1");
  //TString strang_ram("pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run39_10runs");
  //TString strang_ram("pythia_ep_noradcor_18x275_q2_0.000000001_1.0_run39");
  //TString strang_ram("rho_10x100_uChannel_Q2of0to10_hiDiv");
  //TString strang_ram("sartre_bnonsat_Au_phi_ab_eAu_q2_15_1_1000runs");
  TString strang_ram("Sartre_Au_phi_10runs");
  //TString strang_ram("pythia8CCDIS_18x275_minQ2=100_beamEffects_xAngle=-0.025_hiDiv_1_1000runs");
  //TString strang_ram("podio_output_100events");
  //TString strang_ram("pythia8NCDIS_18x275_minQ2=1_beamEffects_xAngle=-0.025_hiDiv_1_100files");   
  
  const char *strang=strang_ram.Data();

  cout << "Analyzed data will be of the type:\n " << strang << " .\n";

  // define flavor of this plotting macro:
  
  TString flavor_ram("nHCal_VM");
  const char *flavor=flavor_ram.Data();

  //define and create, if not existing, pdf output directory:  
  TString pdfdir = strang;
  if (!pdfdir.IsNull()) {
    // Directory does not exist: try to make it
    gSystem->mkdir(pdfdir.Data(), kTRUE);
    cout << "Created output pdf directory:\n " << pdfdir << " .\n";
  }
  
  // define and open input file:
  TString infile_ram= TString("out.") + strang + TString("-") + flavor + TString(".root");
  const char *infile=infile_ram.Data();
  TFile *ifile = TFile::Open(infile,"READ");
  
  cout << "Reading infile:\n " << infile << " .\n";


  ///////////////////////////////////////////////////////////
  //Histograms to plot: 
  
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
  
  // FILE 7 - kpmfromphidecay momentum
  
  TString name7 = TString("kpmfromphi_momentum");
  TString filename7 = pdfdir + TString("/") + TString(name7) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas7 = new TCanvas(name7, strang, 800, 600);
  kpmfromphiRecMom->SetTitle(strang);
  kpmfromphiRecMom->SetLineColor(kBlack);
  kpmfromphiRecMom->Draw();
  kpmfromphiRecMom_nHCal->SetLineColor(kRed);
  kpmfromphiRecMom_nHCal->Draw("same");  
  canvas7->Draw();

  auto leg7 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header  
  leg7->SetHeader("Kaon momenta from #phi decay", "C"); // option "C" allows to center the header
  leg7->SetFillStyle(0);
  leg7->AddEntry(kpmfromphiRecMom,"all","l");
  leg7->AddEntry(kpmfromphiRecMom_nHCal,"in nHCal acceptance","l");
  leg7->Draw();
  canvas7->Print(filename7, "pdf");          
  // end file 7
  

  
  // FILE 1 - generated eta //
  
  // Define the name of the plot:
  TString name1 = TString("trueEta_species");
  // Define the name of the pdf file:
  TString filename1 = pdfdir + TString("/") + TString(name1) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box
  
  TCanvas *canvas = new TCanvas(name1, strang, 800, 600);
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
  Int_t binmax_1 = electronEta->GetMaximumBin();
  Double_t y_max1 = electronEta->GetBinContent(binmax_1);
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

  //
  canvas->Print(filename1, "pdf");          
  // end of generated eta

  // FILE 2 - reconstructed eta //
  
  // Define the name of the plot:
  TString name2 = TString("recEta_species");
  // Define the name of the pdf file:
  TString filename2 = pdfdir + TString("/") + TString(name2) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box    

  TCanvas *canvas2 = new TCanvas(name2, strang, 800, 600);
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

  canvas2->Draw();

  auto leg2 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header                                                                         
  leg2->SetHeader("Reconstructed particles", "C"); // option "C" allows to center the header
  leg2->SetFillStyle(0);
  leg2->AddEntry(electronRecEta,"electrons (#pm)","l");
  leg2->AddEntry(muonRecEta,"muons (#pm)","l");
  leg2->AddEntry(pionRecEta,"pions (#pm)","l");
  leg2->AddEntry(kaonRecEta,"kaons (#pm)","l");
  leg2->AddEntry(protonRecEta,"protons (#pm)","l");
  leg2->Draw();

  // add vertical lines for nHCal acceptance                                                                                              
  Int_t binmax_2 = electronRecEta->GetMaximumBin();
  Double_t y_max2 = electronRecEta->GetBinContent(binmax_2);
  TLine *eta_min_nhcal_line2= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max2);  // (x1,y1,x2,y2)
  eta_min_nhcal_line2->SetLineColor(kBlack);
  eta_min_nhcal_line2->SetLineWidth(2);
  eta_min_nhcal_line2->SetLineStyle(kDashed);
  eta_min_nhcal_line2->Draw("same");
  TLine *eta_max_nhcal_line2= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max2);  // (x1,y1,x2,y2)
  eta_max_nhcal_line2->SetLineColor(kBlack);
  eta_max_nhcal_line2->SetLineWidth(2);
  eta_max_nhcal_line2->SetLineStyle(kDashed);
  eta_max_nhcal_line2->Draw("same");
  
  //
  canvas2->Print(filename2, "pdf");
  // end of reconstructed eta  

  // FILE 3 - gen & rec eta: //                                                                                                         

  // Define the name of the plot:
  TString name3 = TString("gen-recEta_species");
  // Define the name of the pdf file:
  TString filename3 = pdfdir + TString("/") + TString(name3) + TString(".pdf");

  TCanvas *canvas3 = new TCanvas(name3, strang, 1200, 600);

  //devide the canvas into several pads (x,y):
  canvas3->Divide(3,2);

  gStyle->SetOptStat(1111); //now we want some stat boxes  
  
  // go to the first one and do your thing
  canvas3->cd(1);
  electronEta->SetTitle(strang);
  electronEta->SetLineColor(kBlue);
  electronEta->Draw();
  electronRecEta->SetLineStyle(2);
  electronRecEta->Draw("same");

  auto leg31 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
  leg31->SetBorderSize(0);
  leg31->SetTextSize(0.05);
  leg31->SetFillStyle(0);
  leg31->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg31->AddEntry(electronEta,"electrons (#pm) gen","l");
  leg31->AddEntry(electronRecEta,"electrons (#pm) rec","l");
  leg31->Draw();

  // add vertical lines for nHCal acceptance (use existing lines for generated eta)
  eta_min_nhcal_line1->Draw("same");
  eta_max_nhcal_line1->Draw("same");
  
  // pad 2
  canvas3->cd(2);
  muonEta->SetTitle(strang);
  muonEta->SetLineColor(kMagenta);
  muonEta->Draw();
  muonRecEta->SetLineStyle(2);
  muonRecEta->Draw("same");

  auto leg32 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header                                                                       
  leg32->SetBorderSize(0);
  leg32->SetTextSize(0.05);
  leg32->SetFillStyle(0);
  leg32->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg32->AddEntry(muonEta,"muons (#pm) gen","l");
  leg32->AddEntry(muonRecEta,"muons (#pm) rec","l");
  leg32->Draw();

// add vertical lines for nHCal acceptance
  Int_t binmax_32 = muonEta->GetMaximumBin();
  Double_t y_max32 = muonEta->GetBinContent(binmax_32);
  TLine *eta_min_nhcal_line32= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max32);
  eta_min_nhcal_line32->SetLineColor(kBlack);
  eta_min_nhcal_line32->SetLineWidth(2);
  eta_min_nhcal_line32->SetLineStyle(kDashed);
  eta_min_nhcal_line32->Draw("same");
  TLine *eta_max_nhcal_line32= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max32);
  eta_max_nhcal_line32->SetLineColor(kBlack);
  eta_max_nhcal_line32->SetLineWidth(2);
  eta_max_nhcal_line32->SetLineStyle(kDashed);
  eta_max_nhcal_line32->Draw("same");

  // pad 3
  canvas3->cd(3);
  pionEta->SetTitle(strang);
  pionEta->SetLineColor(kOrange+1);
  pionEta->Draw();
  pionRecEta->SetLineStyle(2);
  pionRecEta->Draw("same");

  auto leg33 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
                                                                                                                                          
  leg33->SetBorderSize(0);
  leg33->SetTextSize(0.05);
  leg33->SetFillStyle(0);
  leg33->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header                                          
  leg33->AddEntry(pionEta,"pions (#pm) gen","l");
  leg33->AddEntry(pionRecEta,"pions (#pm) rec","l");
  leg33->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_33 = pionEta->GetMaximumBin();
  Double_t y_max33 = pionEta->GetBinContent(binmax_33);
  TLine *eta_min_nhcal_line33= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max33);  
  eta_min_nhcal_line33->SetLineColor(kBlack);
  eta_min_nhcal_line33->SetLineWidth(2);
  eta_min_nhcal_line33->SetLineStyle(kDashed);
  eta_min_nhcal_line33->Draw("same");
  TLine *eta_max_nhcal_line33= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max33);  
  eta_max_nhcal_line33->SetLineColor(kBlack);
  eta_max_nhcal_line33->SetLineWidth(2);
  eta_max_nhcal_line33->SetLineStyle(kDashed);
  eta_max_nhcal_line33->Draw("same");
  
  // pad 4
  canvas3->cd(4);
  protonEta->SetTitle(strang);
  protonEta->SetLineColor(kRed);
  protonEta->Draw();
  protonRecEta->SetLineStyle(2);
  protonRecEta->Draw("same");

  auto leg34 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
					      
  leg34->SetBorderSize(0);
  leg34->SetTextSize(0.05);
  leg34->SetFillStyle(0);
  leg34->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg34->AddEntry(protonEta,"protons (#pm) gen","l");
  leg34->AddEntry(protonRecEta,"protons (#pm) rec","l");
  leg34->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_34 = protonEta->GetMaximumBin();
  Double_t y_max34 = protonEta->GetBinContent(binmax_34);
  TLine *eta_min_nhcal_line34= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max34);
  eta_min_nhcal_line34->SetLineColor(kBlack);
  eta_min_nhcal_line34->SetLineWidth(2);
  eta_min_nhcal_line34->SetLineStyle(kDashed);
  eta_min_nhcal_line34->Draw("same");
  TLine *eta_max_nhcal_line34= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max34);
  eta_max_nhcal_line34->SetLineColor(kBlack);
  eta_max_nhcal_line34->SetLineWidth(2);
  eta_max_nhcal_line34->SetLineStyle(kDashed);
  eta_max_nhcal_line34->Draw("same");

  // pad 5
  canvas3->cd(5);
  kaonEta->SetTitle(strang);
  kaonEta->SetLineColor(kCyan);
  kaonEta->Draw();
  kaonRecEta->SetLineStyle(2);
  kaonRecEta->Draw("same");

  auto leg35 = new TLegend(0.48,0.6,0.68,0.88); //x1,y1,x2,y2,header
  						  
  leg35->SetBorderSize(0);
  leg35->SetTextSize(0.05);
  leg35->SetFillStyle(0);
  leg35->SetHeader("Generated vs. reconstructed", "C"); // option "C" allows to center the header
  leg35->AddEntry(kaonEta,"kaons (#pm) gen","l");
  leg35->AddEntry(kaonRecEta,"kaons (#pm) rec","l");
  leg35->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_35 = kaonEta->GetMaximumBin();
  Double_t y_max35 = kaonEta->GetBinContent(binmax_35);
  TLine *eta_min_nhcal_line35= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max35);
  eta_min_nhcal_line35->SetLineColor(kBlack);
  eta_min_nhcal_line35->SetLineWidth(2);
  eta_min_nhcal_line35->SetLineStyle(kDashed);
  eta_min_nhcal_line35->Draw("same");
  TLine *eta_max_nhcal_line35= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max35);
  eta_max_nhcal_line35->SetLineColor(kBlack);
  eta_max_nhcal_line35->SetLineWidth(2);
  eta_max_nhcal_line35->SetLineStyle(kDashed);
  eta_max_nhcal_line35->Draw("same");
  
  // draw canvas and print to pdf
  canvas3->Draw();
  canvas3->Print(filename3, "pdf");
  // end of gen & rec eta       

  // FILE 4 - eta decay rho0 to pi+ pi- //                                                                                                           
  // Define the name of the plot:
  TString name4 = TString("Eta_decay_rho0");
  // Define the name of the pdf file:
  TString filename4 = pdfdir + TString("/") + TString(name4) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas4 = new TCanvas(name4, strang, 800, 600);
  rho0Eta->SetLineColor(kBlack);
  rho0Eta->SetLineStyle(1);
  rho0Eta->Draw();
  pipmfromrho0RecEta->SetLineStyle(2);
  pipmfromrho0RecEta->SetTitle(strang);
  pipmfromrho0RecEta->SetLineColor(kRed);
  pipmfromrho0RecEta->Draw("same");
  
  canvas4->Draw();

  auto leg4 = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg4->SetBorderSize(0);
  leg4->SetFillStyle(0);
  leg4->SetTextSize(0.05);
  leg4->SetHeader("generated #rho^{0}(770) and decay pions", "C"); 
  leg4->AddEntry(pipmfromrho0RecEta,"reco pions (#pm)","l");
  leg4->AddEntry(rho0Eta,"gen #rho^{0}","l");
  leg4->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_4 = pipmfromrho0Eta->GetMaximumBin();
  Double_t y_max4 = 0.6*pipmfromrho0Eta->GetBinContent(binmax_4);
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
  //                                                                                                                                        
  canvas4->Print(filename4, "pdf");
  // end of decay rho0 to pipi eta   

  // FILE 5 - eta decay phi to K+ K- //                                                                                                           
  // Define the name of the plot:
  TString name5 = TString("Eta_decay_phi");
  // Define the name of the pdf file:
  TString filename5 = pdfdir + TString("/") + TString(name5) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas5 = new TCanvas(name5, strang, 800, 600);
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
  
  canvas5->Draw();

  auto leg5 = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg5->SetBorderSize(0);
  leg5->SetFillStyle(0);
  leg5->SetTextSize(0.05);
  leg5->SetHeader("generated #phi(1020) and decay kaons", "C"); 
  leg5->AddEntry(kpmfromphiRecEta,"reco kaons (#pm)","l");
  leg5->AddEntry(phiEta,"gen #phi","l");
  leg5->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_5 = kpmfromphiEta->GetMaximumBin();
  Double_t y_max5 = 0.6*kpmfromphiEta->GetBinContent(binmax_5);
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
  //                                                                                                                                        
  canvas5->Print(filename5, "pdf");
  // end of decay phi to KK eta

  // FILE 6 - eta decay jpsi to e+ e- //                                                                                                           
  // Define the name of the plot:
  TString name6 = TString("Eta_decay_jpsi");
  // Define the name of the pdf file:
  TString filename6 = pdfdir + TString("/") + TString(name6) + TString(".pdf");

  gStyle->SetOptStat(0); //no stats box                 

  TCanvas *canvas6 = new TCanvas(name6, strang, 800, 600);
  // fat channels:
  epmfromjpsiRecEta->SetLineStyle(2);
  epmfromjpsiRecEta->SetTitle(strang);
  epmfromjpsiRecEta->SetLineColor(kRed);
  epmfromjpsiRecEta->Draw();
  jpsiEta->SetLineStyle(1);
  jpsiEta->SetLineColor(kBlack);
  jpsiEta->SetLineStyle(1);
  jpsiEta->Draw("same");
  
  canvas6->Draw();

  auto leg6 = new TLegend(0.25,0.6,0.75,0.88); //x1,y1,x2,y2,header
  leg6->SetBorderSize(0);
  leg6->SetFillStyle(0);
  leg6->SetTextSize(0.05);
  leg6->SetHeader("generated J/#psi and decay electrons", "C"); 
  leg6->AddEntry(epmfromjpsiRecEta,"reco electrons (#pm)","l");
  leg6->AddEntry(jpsiEta,"gen J/#psi","l");
  leg6->Draw();

  // add vertical lines for nHCal acceptance
  Int_t binmax_6 = epmfromjpsiEta->GetMaximumBin();
  Double_t y_max6 = 0.6*epmfromjpsiEta->GetBinContent(binmax_6);
  TLine *eta_min_nhcal_line6= new TLine(eta_min_nhcal,0.,eta_min_nhcal,y_max6);
  eta_min_nhcal_line6->SetLineColor(kBlack);
  eta_min_nhcal_line6->SetLineWidth(2);
  eta_min_nhcal_line6->SetLineStyle(kDashed);
  eta_min_nhcal_line6->Draw("same");
  TLine *eta_max_nhcal_line6= new TLine(eta_max_nhcal,0.,eta_max_nhcal,y_max6);
  eta_max_nhcal_line6->SetLineColor(kBlack);
  eta_max_nhcal_line6->SetLineWidth(2);
  eta_max_nhcal_line6->SetLineStyle(kDashed);
  eta_max_nhcal_line6->Draw("same");
  //                                                                                                                                        
  canvas6->Print(filename6, "pdf");
  // end of decay jspi to ee eta   
  
  //
  cout << "Thank you for running Caro's macro.\n";
  gSystem->Exec("date");
  
} // end of macro  
//////////////////

  //TPaveText *t = new TPaveText(.05,.3,.95,.6, "NDC");                                                                                   
  //t->AddText("This line is blue"); ((TText*)t->GetListOfLines()->Last())->SetTextColor(kBlue);                                          
  //t->AddText("This line is red");  ((TText*)t->GetListOfLines()->Last())->SetTextColor(kRed);                                           
  //t->Draw();     
