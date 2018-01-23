
void grshade() {
   gROOT->SetBatch();
   gStyle->SetOptStat(0);

   TCanvas *c1 = new TCanvas("c1","Expected Bias Band Graph",200,10,800,800);

   c1->SetGrid();
   c1->DrawFrame(0.85,0.75,1.15,1.25);
   
   const Int_t n = 3;
   Double_t x[n], y[n],
       ymin[n], ymax[n],
       ymin2[n], ymax2[n];
   Int_t i;
   /*
   for (i=0;i<n;i++) {
     x[i] = 0.1+i*0.1;
     ymax[i] = 10*sin(x[i]+0.2);
     ymin[i] = 8*sin(x[i]+0.1);
     y[i] = 9*sin(x[i]+0.15);
   }
   */
   x    [0] = 0.9;
   y    [0] = 0.8999; // median
   ymin [0] = (0.8395 + 0.8405) * 0.5; // sigma dev
   ymax [0] = 0.9645;
   ymin2[0] = (0.7835 + 0.7845) * 0.5; // 2 sigma dev
   ymax2[0] = (2*1.0345 + 4*1.0335) / (2+4);

   x    [1] = 1.0;
   y    [1] = 1.005;
   ymin [1] = 0.934;
   ymax [1] = 1.07;
   ymin2[1] = 0.8735;
   ymax2[1] = (1*1.1465 + 4*1.1455) / (1+4);

   x    [2] = 1.1;
   y    [2] = 1.0999;
   ymin [2] = (1.0285 + 1.0295) * 0.5;
   ymax [2] = 1.1764;
   ymin2[2] = 0.9625;
   ymax2[2] = (2*1.2585 + 4*1.2575) / (2+4);

   TGraph *grmin = new TGraph(n,x,ymin);
   TGraph *grmax = new TGraph(n,x,ymax);
   TGraph *gr    = new TGraph(n,x,y);

   TGraph *grshade  = new TGraph(2*n);
   TGraph *grshade2 = new TGraph(2*n);
   for (i=0;i<n;i++) {
      grshade->SetPoint(i,x[i],ymax[i]);
      grshade->SetPoint(n+i,x[n-i-1],ymin[n-i-1]);
      grshade2->SetPoint(i,x[i],ymax2[i]);
      grshade2->SetPoint(n+i,x[n-i-1],ymin2[n-i-1]);
   }

   //grshade->SetFillStyle(3013);
   //grshade2->SetFillColor(20);
   //grshade2->SetFillColor(18);
   grshade2->SetFillColor(5); // brazil plot is about 5 yellow 2sigma 3 green 1sigma
   grshade2->Draw("f");
   //grshade->SetFillColor(25);
   grshade->SetFillColor(3);
   grshade->Draw("same f");

   //grmin->Draw("l");
   //grmax->Draw("l");

   gr->SetLineWidth(0);
   //gr->SetLineStyle(10);
   gr->SetMarkerColor(4);
   gr->SetMarkerStyle(21);
   gr->Draw("CP");

   TLine* l = new TLine(0.86, 0.86, 1.14, 1.14);
   l->SetLineWidth(1);
   //l->SetLineStyle(4);
   //l->SetLineColor(kRed);
   l->Draw("same");

   c1->SaveAs("expected-limits-band.png");
}

