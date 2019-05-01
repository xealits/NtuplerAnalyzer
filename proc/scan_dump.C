int scan_dump(TString filename, TString out_file_name){

TFile tfile(filename);
TTree* ttree = (TTree*) tfile.Get("ntupler/reduced_ttree");

((TTreePlayer*)(ttree->GetPlayer()))->SetScanRedirect(true);
((TTreePlayer*)(ttree->GetPlayer()))->SetScanFileName(out_file_name);

ttree->Scan("indexevents:lep_p4[0].pt()", "HLT_mu && leps_ID == -11*13 && ((abs(lep_id[0]) == 13 && lep_matched_HLT[0]) || (abs(lep_id[1]) == 13 && lep_matched_HLT[1])) && lep_p4[0].pt()>30 && lep_p4[1].pt()>30 && abs(lep_p4[0].eta()) < 2.4 && abs(lep_p4[1].eta()) < 2.4 && lep_dxy[0] <0.01 && lep_dz[0]<0.02");

return 0;
}

