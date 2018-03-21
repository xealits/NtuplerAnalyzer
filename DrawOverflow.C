TH1F *DrawOverflow(TH1F* h)
{
//function to paint the histogram h with an extra bin for overflows
UInt_t nx = h->GetNbinsX()+1; Double_t *xbins= new Double_t[nx+1];
for (UInt_t i=0;i<nx;i++)
        xbins[i]=h->GetBinLowEdge(i+1);
xbins[nx]=xbins[nx-1]+h->GetBinWidth(nx);
//book a temporary histogram having extra bins for overflows
TH1F *htmp = new TH1F(h->GetName(), h->GetTitle(), nx, xbins);
htmp->Sumw2();
//fill the new histogram including the overflows
for (UInt_t i=1; i<=nx; i++)
        {
        htmp->SetBinContent(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinContent(i)); htmp->SetBinError(htmp->FindBin(htmp->GetBinCenter(i)),h->GetBinError(i));
        }
htmp->SetBinContent(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinContent(0)); htmp->SetBinError(htmp->FindBin(h->GetBinLowEdge(1)-1), h->GetBinError(0));
// Restore the number of entries
htmp->SetEntries(h->GetEffectiveEntries());
return htmp;
}


