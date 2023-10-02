#include <iostream>
#include <TFile.h>
#include <TH1.h>
#include <TRandom3.h>
#include <cmath>

void GenerateSystematicHistogram()
{
  // Open the input file and retrieve the histogram
  TFile* inputFile = new TFile("cms_bckg.root", "READ");
  TH1D* inputHist = (TH1D*)inputFile->Get("cms_SR");

  // Clone the input histogram to the output file
  TH1D* outputHist = (TH1D*)inputHist->Clone();
  outputHist->SetName("cms_SR");

  // Create the systematic histograms
  TH1D* systHist_hf_scale = new TH1D("cms_SR_HF_SCALE_SYS", "CMS Signal Region with Systematic (HF scale)", 17, 2.5, 19.5);
  TH1D* systHist_nch_eff = new TH1D("cms_SR_Nch_SYS", "CMS Signal Region with Systematic (Effect of nch)", 17, 2.5, 19.5);
  TH1D* systHist_lumi = new TH1D("cms_SR_LUMI_SYS", "CMS Signal Region with Systematic (Luminosity)", 17, 2.5, 19.5);

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
    double stdDev_hf_scale = content * 0.067;
    double stdDev_nch_eff = content * 0.036;
    double stdDev_lumi = content * 0.05;

    // Generate 1000 systematic variations for each bin
    double sumSquares_hf_scale = 0.0;
    double sumSquares_nch_eff = 0.0;
    double sumSquares_lumi = 0.0;

    for (int j = 0; j < 10000; ++j)
    {
      // Generate random numbers from Gaussian distributions using the means and standard deviations
      double variation_hf_scale = gRandom->Gaus(content, stdDev_hf_scale);
      double variation_nch_eff = gRandom->Gaus(content, stdDev_nch_eff);
      double variation_lumi = gRandom->Gaus(content, stdDev_lumi);

      // Add the squares of the variations to the sum of squares for each variation
      sumSquares_hf_scale += (variation_hf_scale - content) * (variation_hf_scale - content);
      sumSquares_nch_eff += (variation_nch_eff - content) * (variation_nch_eff - content);
  
      sumSquares_lumi += (variation_lumi - content) * (variation_lumi - content);
    }

    // Calculate the RMS for each variation as the square root of the average of the sum of squares
    double rms_hf_scale = std::sqrt(sumSquares_hf_scale / 10000.0);
    double rms_nch_eff = std::sqrt(sumSquares_nch_eff / 10000.0);

    double rms_lumi = std::sqrt(sumSquares_lumi / 10000.0);

    // Set the bin contents for each variation
    systHist_hf_scale->SetBinContent(i, rms_hf_scale + content);
    systHist_nch_eff->SetBinContent(i, rms_nch_eff + content);
    systHist_lumi->SetBinContent(i, rms_lumi + content);

    // Set the bin errors for each variation
    if (content == 0.0){
      systHist_hf_scale->SetBinError(i, 0.0);
      systHist_nch_eff->SetBinError(i, 0.0);
      systHist_lumi->SetBinError(i, 0.0);
    } else {
      systHist_hf_scale->SetBinError(i, (contentError / content) * systHist_hf_scale->GetBinContent(i));
      systHist_nch_eff->SetBinError(i, (contentError / content) * systHist_nch_eff->GetBinContent(i));

      systHist_lumi->SetBinError(i, (contentError / content) * systHist_lumi->GetBinContent(i));
    }

  }

  // // Write histograms to the output file

  // outputHist->Write();
  // systHist_hf_scale->Write();
  // systHist_nch_eff->Write();
  // systHist_lumi->Write();
  
  // outputFile->mkdir("nominal");
  // outputFile->cd("nominal");
  // outputHist->Write();
  // systHist_hf_scale->Write();
  // systHist_nch_eff->Write();
  // systHist_lumi->Write();
  // outputFile->cd("..");

    // Create a new histogram for the systematic variations
  TFile* outputFile = new TFile("cms_bckg_syst.root", "RECREATE");

    // Write histograms to the output file
  outputHist->Write();

  gDirectory->mkdir ("cms_SR_LUMI_SYS");
  gDirectory->cd ("cms_SR_LUMI_SYS");
  systHist_lumi->Write();
  gDirectory->cd ("..");

  gDirectory->mkdir ("cms_SR_HF_SCALE_SYS");
  gDirectory->cd ("cms_SR_HF_SCALE_SYS");
  systHist_hf_scale->Write();
  gDirectory->cd (".."); 

  gDirectory->mkdir ("cms_SR_Nch_SYS");
  gDirectory->cd ("cms_SR_Nch_SYS");
  systHist_nch_eff->Write();
  gDirectory->cd ("..");    

  // Clean up
  delete inputFile;
  delete outputFile;
}

int main()
{
  GenerateSystematicHistogram();
  return 0;
}


// g++ -o crude_approx_bckg crude_approx_bckg.C $(root-config --cflags --libs)
