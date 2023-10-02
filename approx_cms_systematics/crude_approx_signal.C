#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TRandom3.h>
#include <cmath>

void GenerateSystematicHistogram()
{
  // Open the input file and retrieve the histogram
  TFile* inputFile = new TFile("cms_tau_mu3prong_orig.root", "READ");
  TH1D* inputHist = (TH1D*)inputFile->Get("cms_SR");

  // Clone the input histogram to the output file
  TH1D* outputHist = (TH1D*)inputHist->Clone();
  outputHist->SetName("cms_SR");

  // Create the systematic histograms
  TH1D* systHist_muon_eff = new TH1D("cms_SR_syst_muon_eff", "CMS Signal Region with Systematic (Muon Efficiency)", 17, 2.5, 19.5);
  TH1D* systHist_pion_eff = new TH1D("cms_SR_syst_pion_eff", "CMS Signal Region with Systematic (Pion Efficiency)", 17, 2.5, 19.5);
  TH1D* systHist_simSampleSizeBinByBinVariation = new TH1D("cms_SR_syst_simSampleSizeBinByBinVariation", "CMS Signal Region with Systematic (Simulated Sample Size Bin-By-Bin Variation)", 17, 2.5, 19.5);
  TH1D* systHist_tLeptonBRVariation = new TH1D("cms_SR_syst_tLeptonBRVariation", "CMS Signal Region with Systematic (t-Lepton Branching Ratio Variation)", 17, 2.5, 19.5);
  TH1D* systHist_simSampleSizeEfficiencyVariation = new TH1D("cms_SR_syst_simSampleSizeEfficiencyVariation", "CMS Signal Region with Systematic (Simulated Sample Size Efficiency Variation)", 17, 2.5, 19.5);
  TH1D* systHist_lumi = new TH1D("cms_SR_syst_lumi", "CMS Signal Region with Systematic (Luminosity)", 17, 2.5, 19.5);

  // Set the random seed for reproducibility
  // gRandom->SetSeed(0);

  // Loop over the bins of the input histogram
  for (int i = 1; i <= inputHist->GetNbinsX(); ++i)
  {
    // Get the content of the input histogram bin
    double content = inputHist->GetBinContent(i);

    if (content == 0.0){
      continue;
    }

    double contentError = inputHist->GetBinError(i);

    // Calculate the standard deviations for each variation
    double stdDev_muon_eff = content * 0.067;
    double stdDev_pion_eff = content * 0.036;
    double stdDev_simSampleSizeBinByBinVariation = content * 0.03;
    double stdDev_tLeptonBRVariation = content * 0.006;
    double stdDev_simSampleSizeEfficiencyVariation = content * 0.011;

    // Generate 1000 systematic variations for each bin
    double sumSquares_muon_eff = 0.0;
    double sumSquares_pion_eff = 0.0;
    double sumSquares_simSampleSizeBinByBinVariation = 0.0;
    double sumSquares_tLeptonBRVariation = 0.0;
    double sumSquares_simSampleSizeEfficiencyVariation = 0.0;

    for (int j = 0; j < 10000; ++j)
    {
      // Generate random numbers from Gaussian distributions using the means and standard deviations
      double variation_muon_eff = gRandom->Gaus(content, stdDev_muon_eff);
      double variation_pion_eff = gRandom->Gaus(content, stdDev_pion_eff);
      double variation_simSampleSizeBinByBinVariation = gRandom->Gaus(content, stdDev_simSampleSizeBinByBinVariation);
      double variation_tLeptonBRVariation = gRandom->Gaus(content, stdDev_tLeptonBRVariation);
      double variation_simSampleSizeEfficiencyVariation = gRandom->Gaus(content, stdDev_simSampleSizeEfficiencyVariation);

      // Add the squares of the variations to the sum of squares for each variation
      sumSquares_muon_eff += (variation_muon_eff - content) * (variation_muon_eff - content);
      sumSquares_pion_eff += (variation_pion_eff - content) * (variation_pion_eff - content);
      sumSquares_simSampleSizeBinByBinVariation += (variation_simSampleSizeBinByBinVariation - content) * (variation_simSampleSizeBinByBinVariation - content);
      sumSquares_tLeptonBRVariation += (variation_tLeptonBRVariation - content) * (variation_tLeptonBRVariation - content);
      sumSquares_simSampleSizeEfficiencyVariation += (variation_simSampleSizeEfficiencyVariation - content) * (variation_simSampleSizeEfficiencyVariation - content);
    }

    // Calculate the RMS for each variation as the square root of the average of the sum of squares
    double rms_muon_eff = std::sqrt(sumSquares_muon_eff / 10000.0);
    double rms_pion_eff = std::sqrt(sumSquares_pion_eff / 10000.0);
    double rms_simSampleSizeBinByBinVariation = std::sqrt(sumSquares_simSampleSizeBinByBinVariation / 10000.0);
    double rms_tLeptonBRVariation = std::sqrt(sumSquares_tLeptonBRVariation / 10000.0);
    double rms_simSampleSizeEfficiencyVariation = std::sqrt(sumSquares_simSampleSizeEfficiencyVariation / 10000.0);

    // Set the bin contents for each variation
    systHist_muon_eff->SetBinContent(i, rms_muon_eff + content);
    systHist_pion_eff->SetBinContent(i, rms_pion_eff + content);
    systHist_simSampleSizeBinByBinVariation->SetBinContent(i, rms_simSampleSizeBinByBinVariation + content);
    systHist_tLeptonBRVariation->SetBinContent(i, rms_tLeptonBRVariation + content);
    systHist_simSampleSizeEfficiencyVariation->SetBinContent(i, rms_simSampleSizeEfficiencyVariation + content);

    // Handle systHist_lumi differently by multiplying inputHist by 5% and adding it
    double lumi_variation = content * 0.05; // 5% variation
    systHist_lumi->SetBinContent(i, content + lumi_variation);

    // Set the bin errors for each variation
    if (content == 0.0){
      systHist_muon_eff->SetBinError(i, 0.0);
      systHist_pion_eff->SetBinError(i, 0.0);
      systHist_simSampleSizeBinByBinVariation->SetBinError(i, 0.0);
      systHist_tLeptonBRVariation->SetBinError(i, 0.0);
      systHist_simSampleSizeEfficiencyVariation->SetBinError(i, 0.0);
      systHist_lumi->SetBinError(i, 0.0);
    } else {
      systHist_muon_eff->SetBinError(i, (contentError / content) * systHist_muon_eff->GetBinContent(i));
      systHist_pion_eff->SetBinError(i, (contentError / content) * systHist_pion_eff->GetBinContent(i));
      systHist_simSampleSizeBinByBinVariation->SetBinError(i, (contentError / content) * systHist_simSampleSizeBinByBinVariation->GetBinContent(i));
      systHist_tLeptonBRVariation->SetBinError(i, (contentError / content) * systHist_tLeptonBRVariation->GetBinContent(i));
      systHist_simSampleSizeEfficiencyVariation->SetBinError(i, (contentError / content) * systHist_simSampleSizeEfficiencyVariation->GetBinContent(i));
      systHist_lumi->SetBinError(i, (contentError / content) * systHist_lumi->GetBinContent(i));
    }
  }

  // Write histograms to the output file
  TFile* outputFile = new TFile("cms_tau_mu3prong_syst.root", "RECREATE");
  outputHist->Write();
  systHist_muon_eff->Write();
  systHist_pion_eff->Write();
  systHist_simSampleSizeBinByBinVariation->Write();
  systHist_tLeptonBRVariation->Write();
  systHist_simSampleSizeEfficiencyVariation->Write();
  systHist_lumi->Write();

  // Clean up
  delete inputFile;
  delete outputFile;
}

int main()
{
  GenerateSystematicHistogram();
  return 0;
}


// g++ -o crude_approx_signal crude_approx_signal.C $(root-config --cflags --libs)
