void PileupHistograms_0905_63mb_pm5_hMC25ns()
{
//=========Macro generated from canvas: c1/c1
//=========  (Mon Mar 13 16:24:38 2017) by ROOT version6.06/01
   TCanvas *c1 = new TCanvas("c1", "c1",65,52,700,500);
   c1->Range(-12.5,-0.008903069,112.5,0.08012762);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *hMC25ns__2 = new TH1F("hMC25ns__2","",100,0,100);
   hMC25ns__2->SetBinContent(1,0.0008293129);
   hMC25ns__2->SetBinContent(2,0.001242761);
   hMC25ns__2->SetBinContent(3,0.003393292);
   hMC25ns__2->SetBinContent(4,0.004082247);
   hMC25ns__2->SetBinContent(5,0.003830366);
   hMC25ns__2->SetBinContent(6,0.006591593);
   hMC25ns__2->SetBinContent(7,0.008160227);
   hMC25ns__2->SetBinContent(8,0.009436408);
   hMC25ns__2->SetBinContent(9,0.01377774);
   hMC25ns__2->SetBinContent(10,0.01705939);
   hMC25ns__2->SetBinContent(11,0.0213193);
   hMC25ns__2->SetBinContent(12,0.02473432);
   hMC25ns__2->SetBinContent(13,0.02808488);
   hMC25ns__2->SetBinContent(14,0.03233085);
   hMC25ns__2->SetBinContent(15,0.03703943);
   hMC25ns__2->SetBinContent(16,0.04569177);
   hMC25ns__2->SetBinContent(17,0.05587629);
   hMC25ns__2->SetBinContent(18,0.05769562);
   hMC25ns__2->SetBinContent(19,0.06253253);
   hMC25ns__2->SetBinContent(20,0.05916037);
   hMC25ns__2->SetBinContent(21,0.06566508);
   hMC25ns__2->SetBinContent(22,0.0678329);
   hMC25ns__2->SetBinContent(23,0.06251422);
   hMC25ns__2->SetBinContent(24,0.05480684);
   hMC25ns__2->SetBinContent(25,0.05038933);
   hMC25ns__2->SetBinContent(26,0.04020982);
   hMC25ns__2->SetBinContent(27,0.0374447);
   hMC25ns__2->SetBinContent(28,0.02996616);
   hMC25ns__2->SetBinContent(29,0.02720248);
   hMC25ns__2->SetBinContent(30,0.02193284);
   hMC25ns__2->SetBinContent(31,0.01795866);
   hMC25ns__2->SetBinContent(32,0.01429267);
   hMC25ns__2->SetBinContent(33,0.008399417);
   hMC25ns__2->SetBinContent(34,0.005223664);
   hMC25ns__2->SetBinContent(35,0.00224458);
   hMC25ns__2->SetBinContent(36,0.000779275);
   hMC25ns__2->SetBinContent(37,0.0001970666);
   hMC25ns__2->SetBinContent(38,7.160318e-05);
   hMC25ns__2->SetEntries(100);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *AText = ptstats->AddText("hMC25ns");
   AText->SetTextSize(0.0368);
   AText = ptstats->AddText("Entries = 100    ");
   AText = ptstats->AddText("Mean  =  19.82");
   AText = ptstats->AddText("Std Dev   =  6.268");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hMC25ns__2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hMC25ns__2);

   Int_t ci;      // for color index setting
   TColor *color; // for color definition with alpha
   ci = TColor::GetColor("#000099");
   hMC25ns__2->SetLineColor(ci);
   hMC25ns__2->GetXaxis()->SetLabelFont(42);
   hMC25ns__2->GetXaxis()->SetLabelSize(0.035);
   hMC25ns__2->GetXaxis()->SetTitleSize(0.035);
   hMC25ns__2->GetXaxis()->SetTitleFont(42);
   hMC25ns__2->GetYaxis()->SetLabelFont(42);
   hMC25ns__2->GetYaxis()->SetLabelSize(0.035);
   hMC25ns__2->GetYaxis()->SetTitleSize(0.035);
   hMC25ns__2->GetYaxis()->SetTitleFont(42);
   hMC25ns__2->GetZaxis()->SetLabelFont(42);
   hMC25ns__2->GetZaxis()->SetLabelSize(0.035);
   hMC25ns__2->GetZaxis()->SetTitleSize(0.035);
   hMC25ns__2->GetZaxis()->SetTitleFont(42);
   hMC25ns__2->Draw("");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
