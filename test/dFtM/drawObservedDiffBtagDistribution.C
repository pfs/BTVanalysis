#include "TFile.h"
#include "TCanvas.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "TLegend.h"
#include "THStack.h"
#include "TSystem.h"

#include <iostream>
#include <vector>
#include <map>

using namespace std;

Int_t nTagsPerCat=5;
TString outDir="diffbtag";
TString catLabel="";
TObjArray toSave;

void reweightPredictionToExpectation(TH1F *targetBtag, TH1F *modelBtag);
TH1F *showBtagMultiplicity(THStack *expectations, TH1F *data, TString title, TString name);
TH1F *dyTemplateClosure(TFile *inF,TString dist);
TH1F *dyDataDrivenTemplate(TFile *inF,TString dist);
void addCompatibilityPave(TH1 *data,TH1 *model);
void drawDistributionsFrom(TFile *inF,TString corrDist,TString rawDist="");
void drawObservedBtagDistribution(TString url,TString tagger);
void drawObservedDiffBtagDistribution(TString url,TString tagger);
void saveDistributions(TString tagger);


//
void saveDistributions(TString tagger)
{
  if(toSave.GetEntriesFast()==0) return;
  TFile *fOut=TFile::Open(outDir+"/"+tagger+"_btags.root","RECREATE");
  for(Int_t i=0; i<toSave.GetEntriesFast(); i++)
    {
      TH1 *h=(TH1 *)toSave.At(i)->Clone();
      h->SetDirectory(fOut);
      h->Write();
    }
  fOut->Close();
  cout << "Saved " << toSave.GetEntriesFast() << " distributions for analysis @ " << outDir+"/btagMultForAnalysis.root" << endl;
}

//
void addCompatibilityPave(TH1 *data,TH1 *model)
{
  TPaveText *pave = new TPaveText(0.94,0.4,1.0,0.78,"brNDC");
  pave->SetBorderSize(0);
  pave->SetFillStyle(0);
  pave->SetTextAlign(21);
  pave->SetTextFont(42);
  pave->SetTextColor(kGray+2);
  pave->SetTextSize(0.04);
  char buf[1000];
  sprintf(buf,"#chi^{2}/ndof : %3.2f , K-S prob : %3.2f",data->Chi2Test(model,"WWCHI2/NDF"),data->KolmogorovTest(model,"") );
  pave->AddText(buf)->SetTextAngle(-90);
  pave->Draw();
}

//
void reweightPredictionToExpectation(TH1F *targetBtag, TH1F *modelBtag)
{
  if(targetBtag==0 || modelBtag==0) return;
  for(Int_t ibin=1; ibin<=targetBtag->GetXaxis()->GetNbins()/nTagsPerCat; ibin++)
    {
      Float_t targetTotal(0),modelTotal(0);
      for(Int_t jbin=1; jbin<=nTagsPerCat; jbin++){
	Int_t idx=jbin+(ibin-1)*nTagsPerCat;
	targetTotal += targetBtag->GetBinContent(idx);
	modelTotal  += modelBtag->GetBinContent(idx);
      }

      //if not counts in model, assume the target
      if(modelTotal==0) {
	for(Int_t jbin=1; jbin<=nTagsPerCat; jbin++){
	  Int_t idx=jbin+(ibin-1)*nTagsPerCat;
	  modelBtag->SetBinContent(idx,targetBtag->GetBinContent(idx));
	  modelBtag->SetBinError(idx,targetBtag->GetBinError(idx));
	}
      }
      else{
	Float_t sf=targetTotal/modelTotal;
	cout << ibin << " " << sf << endl;
	for(Int_t jbin=1; jbin<=nTagsPerCat; jbin++){
	  Int_t idx=jbin+(ibin-1)*nTagsPerCat;
	  Float_t cts=modelBtag->GetBinContent(idx);
	  Float_t err=modelBtag->GetBinError(idx);
	  modelBtag->SetBinContent(idx,cts*sf);
	  modelBtag->SetBinError(idx,err*sf);
	}
      }
    }
}

//
TH1F *dyTemplateClosure(TFile *inF,TString dist)
{
  //get plots from file
  TH1F *target=0, *templ=0;
  TString ch[]={"ee","mumu"};
  TString proc("DY");
  //TString proc("Z#rightarrow ll");
  for(size_t i=0; i<2; i++)
    {
      TH1F *itarget=(TH1F *)inF->Get(proc+"/"+ch[i]+"_"+dist);              
      TH1F *itempl=(TH1F *)inF->Get(proc+"/"+ch[i]+"zlowmet_"+dist);
      if(itarget==0 && itempl==0) continue;
      if(target==0){
		target=(TH1F *)itarget->Clone("dy_"+dist);
		target->SetDirectory(0);
		target->SetTitle("Signal reg.");
		templ=(TH1F *)itempl->Clone("dytempl_"+dist);
		templ->SetTitle("Control reg.");
		templ->SetDirectory(0); 
		templ->SetMarkerStyle(20);
      }else{
		target->Add(itarget);
		templ->Add(itempl);
      }
    }
  if(target==0 || templ==0) return 0;

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //perform the reweighting
  reweightPredictionToExpectation(target,templ);
  
  //show result
  THStack *templStack=new THStack;
  templStack->Add(templ,"hist");
  //return showBtagMultiplicity(templStack,target,"CMS simulation, #sqrt{s}=8 TeV, " + proc,"dyclosure"+dist);
  return showBtagMultiplicity(templStack,target,"#bf{CMS} #it{simulation, " + proc +"}   (7 TeV)","dyclosure"+dist);
}

//
TH1F *dyDataDrivenTemplate(TFile *inF,TString dist)
{
  TH1F *data=0,*nonDy=0;
  TList *dirs=inF->GetListOfKeys();
  for(int iproc=0; iproc<dirs->GetEntries(); iproc++)
    {
      TString iDir(dirs->At(iproc)->GetName());
  cout << __LINE__ << endl;

      TH1F *h=(TH1F *)inF->Get(iDir+"/allzlowmet_"+dist);
  cout << __LINE__ << endl;
  	  if(h==0) continue;
      //if(iDir.Contains("Z#rightarrow ll")) continue;
      if(iDir.Contains("DY")) continue;
      if(iDir=="Data"){
		data=(TH1F *)h->Clone("dydatadriven_"+dist);
		data->SetDirectory(0);
      }else if(nonDy==0){
		nonDy=(TH1F *)h->Clone("dydatadrivenres_"+dist);
		nonDy->SetDirectory(0);
      }else{
		nonDy->Add(h);
      }
    }
  cout << __LINE__ << endl;

  //add closure test error in quadrature
  TH1F *dyClosureH=dyTemplateClosure(inF,dist);
  cout << __LINE__ << endl;
  if(data) 
  {
  	data->Add(nonDy,-1);
  	if(dyClosureH)
  	{
  		for(int ibin=1; ibin<=data->GetXaxis()->GetNbins(); ibin++){
    		float unc=fabs(1.0-dyClosureH->GetBinContent(ibin,1));
    		unc=sqrt( pow(data->GetBinError(ibin),2) + pow(data->GetBinContent(ibin)*unc,2) );
    		data->SetBinError(ibin,unc);
  		}
  	}
  	else cout << "Unable to close DY template for " << dist << endl;
  }
  else cout << "Unable to produce DY template for " << dist << endl;

  cout << __LINE__ << endl;
  
  //all done here
  return data;
}

//
TH1F* showBtagMultiplicity(THStack *expectations, TH1F *data, TString title, TString name)
{
  if(outDir!="" && outDir!="./") gSystem->Exec("mkdir -p " + outDir);

  TCanvas *c=new TCanvas(name,name,970,600);
 
  TH1F *totalExpectations=(TH1F *)expectations->GetStack()->At( expectations->GetStack()->GetEntriesFast()-1);
  
  //main plot
  TPad* t1 = new TPad("t1"+name,"t1"+name, 0.0, 0.3, 1.0, 1.0);
  t1->Draw();
  t1->cd();
  t1->SetLeftMargin(0.15);
  t1->SetRightMargin(0.05);
  t1->SetTopMargin(0.1);
  t1->SetBottomMargin(0.05);
  
  TH1F *frame=(TH1F *)totalExpectations->Clone("frame"+name);
  frame->SetDirectory(0);
  frame->Reset("ICE");
  frame->Draw();
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleOffset(1.);
  Float_t maxY(max(totalExpectations->GetMaximum(),data->GetMaximum())*1.1);
  frame->GetYaxis()->SetRangeUser(1e-1,maxY);
  frame->GetXaxis()->SetTitle(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetXaxis()->SetTitleSize(0);
  
  expectations->Draw("histsame");

  TGraphErrors *gr=new TGraphErrors( totalExpectations );
  gr->SetFillStyle(1001);
  gr->SetFillColor(kGray);
  gr->Draw("e2");
  
  data->Draw("e1same");

  TPaveText *pt=new TPaveText(0.13,0.9,0.4,0.99,"brNDC");
  pt->SetBorderSize(0);
  pt->SetFillStyle(0);
  pt->SetTextAlign(12);
  pt->SetTextFont(42);
  pt->SetTextSize(0.04);
  pt->AddText(title);
  pt->Draw();

  if(catLabel!="")
    {
      pt=new TPaveText(0.85,0.9,0.94,0.95,"brNDC");
      pt->SetBorderSize(0);
      pt->SetFillStyle(0);
      pt->SetTextFont(42);
      pt->SetTextAlign(12);
      pt->SetTextSize(0.03);
      pt->SetTextColor(kBlue);
      pt->AddText(catLabel);
      pt->Draw();
    }

  TLegend *leg=new TLegend(0.5,0.9,0.9,0.99,"","brNDC");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  Int_t nentries(expectations->GetStack()->GetEntriesFast());
  for(Int_t i=0; i<nentries; i++)
    {
      TH1F *ih=(TH1F *)expectations->GetStack()->At(i);
      leg->AddEntry(ih,ih->GetTitle(),"f");
    }
  leg->AddEntry(data,data->GetTitle(),"p");
  if(nentries>3) leg->SetNColumns(4);
  else           leg->SetNColumns(nentries+1);
  leg->Draw();

  //auxiliary lines and labels
  TLine *l=new TLine;
  //l->SetLineColor(kGray); 
  l->SetLineStyle(7);
  if(nTagsPerCat==5){
    l->DrawLine(15,0,15,maxY); 
    l->DrawLine(30,0,30,maxY);
    for(int i=0; i<9; i++)
      {
	TString label("="); label += (i%3+2); label+=" jets";
	TPaveText *labelt=new TPaveText(2+5*i,maxY*0.95,3+5*i,maxY,"");
	labelt->SetBorderSize(0);
	labelt->SetFillStyle(0);
	labelt->AddText(label);
	labelt->SetTextFont(52);
	labelt->SetTextSize(0.035);
	labelt->Draw();
      }
  }
  else{
    l->DrawLine(3,0,3,maxY); 
    l->DrawLine(6,0,6,maxY);
  }

  for(int i=0; i<3; i++)
    {
      TString label("ee events");
      if(i==1) label="#mu#mu events";
      else if(i==2)    label="e#mu events";
      TPaveText *labelt=nTagsPerCat==5 ? 
	new TPaveText(6+15*i,maxY*0.9,7+15*i,maxY*0.95,"") :
	new TPaveText(0.5+3*i,maxY*0.95,1.5+5*i,maxY,"") ;
      labelt->SetBorderSize(0);
      labelt->SetFillStyle(0);
      labelt->SetTextFont(42);
      labelt->SetTextAlign(12);
      labelt->SetTextSize(0.035);
      labelt->AddText(label);
      labelt->Draw();
    }
  

  addCompatibilityPave(data,totalExpectations);

  //ratio
  c->cd();
  TPad* t2 = new TPad("t2","t2", 0.0, 0.0, 1.0, 0.3);
  t2->Draw();
  t2->cd();
  t2->SetPad(0,0.0,1.0,0.3);
  t2->SetTopMargin(0.05);
  t2->SetLeftMargin(0.15);
  t2->SetRightMargin(0.05);
  t2->SetBottomMargin(0.4);
  
  TH1F *denRelUncH=(TH1F *) totalExpectations->Clone("relunc"+name);
  denRelUncH->SetDirectory(0);
  for(int xbin=1; xbin<=totalExpectations->GetXaxis()->GetNbins(); xbin++)
    {
      if(totalExpectations->GetBinContent(xbin)==0) continue;
      Double_t err=totalExpectations->GetBinError(xbin)/totalExpectations->GetBinContent(xbin);
      denRelUncH->SetBinContent(xbin,1);
      denRelUncH->SetBinError(xbin,err);
    }
  TGraphErrors *denRelUnc=new TGraphErrors(denRelUncH);
  denRelUnc->SetLineColor(1);
  denRelUnc->SetFillStyle(3001);
  denRelUnc->SetFillColor(kGray);
  denRelUnc->SetMarkerColor(1);
  denRelUnc->SetMarkerStyle(1);
  denRelUncH->Reset("ICE");       
  denRelUncH->Draw();
  denRelUnc->Draw("2");
  float yscale = (1.0-0.2)/(0.28-0);       
  denRelUncH->GetYaxis()->SetTitle("Obs./Exp.");
  denRelUncH->SetMinimum(0.6);
  denRelUncH->SetMaximum(1.5);
  denRelUncH->GetXaxis()->SetTitleOffset(1.3);
  denRelUncH->GetXaxis()->SetLabelSize(0.04*yscale);
  denRelUncH->GetXaxis()->SetTitleSize(0.04*yscale);
  denRelUncH->GetXaxis()->SetTickLength(0.03*yscale);
  denRelUncH->GetYaxis()->SetTitleOffset(0.4);
  denRelUncH->GetYaxis()->SetNdivisions(5);
  denRelUncH->GetYaxis()->SetLabelSize(0.033*yscale);
  denRelUncH->GetYaxis()->SetTitleSize(0.04*yscale);
  
  TH1F *dataToExpectH = (TH1F*) data->Clone("data2exp"+name);
  dataToExpectH->SetDirectory(0);
  dataToExpectH->Divide(totalExpectations);
  dataToExpectH->Draw("same");

  l->DrawLine(15,0.6,15,1.5);
  l->DrawLine(30,0.6,30,1.5);

  c->cd();
  c->Modified();
  c->Update();
  c->SaveAs(outDir+"/"+name+".pdf");
  c->SaveAs(outDir+"/"+name+".png");
  c->SaveAs(outDir+"/"+name+".C");

  //return the closure test
  return dataToExpectH;
}

//
void drawDistributionsFrom(TFile *inF,TString corrDist,TString rawDist)
{
  std::map<TString, TH1F *> rawHistos,corrHistos;
  TList *dirs=inF->GetListOfKeys();
  for(int iproc=0; iproc<dirs->GetEntries(); iproc++)
    {
      TString iDir(dirs->At(iproc)->GetName());
      TH1F *h=(TH1F *)inF->Get(iDir+"/"+rawDist);
      if(h!=0)
      {
	  	h->SetDirectory(0);
	  	rawHistos[iDir]=h;
	  }

      h=(TH1F *)inF->Get(iDir+"/"+corrDist);
      if(h!=0)
      {
	  	h->SetDirectory(0);
	  	corrHistos[iDir]=h;
      }
    }
  if(rawHistos.size()==0) cout << "No raw histos found" << endl;
  if(corrHistos.size()==0)cout << "No corrected histos found" << endl;
  if(rawHistos.size()+corrHistos.size()==0) return;

  TH1F *data=0;
  if(rawHistos.find("Data")!=rawHistos.end())        data=rawHistos["Data"];
  else if(corrHistos.find("Data")!=corrHistos.end()) data=corrHistos["Data"];
  if(data==0) { cout << "No data histo found?" << endl; return; }
  data->SetTitle("Data");
  toSave.Add(data);

  THStack *rawStack        = rawHistos.size() ? new THStack : 0;
  THStack *corrStack       = corrHistos.size() ? new THStack : 0;
  THStack *corrStackWithDD = corrHistos.size() ? new THStack : 0;
  TString procs[]   = {"VV","W","other t#bar{t}","Single top","DY","t#bar{t}"};
  //TString procs[]   = {"VV","W,multijets","other t#bar{t}","Single top","Z#rightarrow ll","t#bar{t}"};
  Int_t   fcolors[] = {592,  809,         822,             824,         831,               614      };
  Int_t   nprocs(sizeof(procs)/sizeof(TString));
  TH1F *otherBkgH=0;
  for(Int_t iproc=0; iproc<nprocs; iproc++)
    {
  cout << __LINE__ << endl;

      if(rawHistos.find(procs[iproc]) != rawHistos.end())
      {
	  	TH1F *h=rawHistos[ procs[iproc] ];
	  	h->SetFillStyle(1001);
	  	h->SetFillColor(fcolors[iproc]);
	  	h->SetLineColor(1);
	  	h->SetMarkerStyle(1);
	  	h->SetMarkerColor(1);
	  	h->SetTitle(procs[iproc]);
	  	rawStack->Add(h,"hist");
	  }
  cout << __LINE__ << endl;

      if(corrHistos.find(procs[iproc]) != corrHistos.end())
      {
	  	TH1F *h=corrHistos[ procs[iproc] ];
	  	h->SetFillStyle(1001);
	  	h->SetFillColor(fcolors[iproc]);
	  	h->SetLineColor(1);
	  	h->SetMarkerStyle(1);
	  	h->SetMarkerColor(1);
	  	h->SetTitle(procs[iproc]);
	  	corrStack->Add(h,"hist");
  cout << __LINE__ << endl;

	  	if(procs[iproc]=="DY")
	    {
	       cout << procs[iproc] << " " << corrDist << " " << h << endl;
	      TH1F *ddh=dyDataDrivenTemplate(inF,corrDist);
  cout << __LINE__ << endl;

	      if(ddh)
	      {
	      	ddh->SetFillStyle(1001);
	      	ddh->SetFillColor(fcolors[iproc]);
	      	ddh->SetLineColor(1);
	      	ddh->SetMarkerStyle(1);
	      	ddh->SetMarkerColor(1);
	      	ddh->SetTitle(procs[iproc] + " (data)");
	      	reweightPredictionToExpectation(h,ddh);
	      	corrStackWithDD->Add(ddh,"hist");
	      	toSave.Add(ddh);
	      }
	      else
	      {
	      	corrStackWithDD->Add(h,"hist");
	      	toSave.Add(h);
	      }
	    }
	  	else
	    {
	      //save other backgrounds estimated from MC
	      if(procs[iproc]!="t#bar{t}")
	      {
			if(otherBkgH==0)
			{
		  		otherBkgH=(TH1F *)h->Clone("othersmc_"+corrDist);
		  		otherBkgH->SetTitle("Other processes");
		  		otherBkgH->SetDirectory(0);
		  		toSave.Add(otherBkgH);
		  	}
		  	else
		  	{
		  		otherBkgH->Add(h);
		  	}
		  }
		  else
		  {
			h->SetName("ttbarmc_"+corrDist);
			toSave.Add(h);
		  }
  cout << __LINE__ << endl;

	     corrStackWithDD->Add(h,"hist");
		}
    }
  }
  cout << __LINE__ << endl;

  //TString header("CMS preliminary, #sqrt{s}=8 TeV, #scale[0.5]{#int}L=19.7 fb^{-1}");
  TString header("#bf{CMS} #it{preliminary}       #scale[0.5]{#int}L=5.0 fb^{-1} (7 TeV)");
  if(rawStack && data) showBtagMultiplicity(rawStack,        data, header, "rawpred"+rawDist);
    cout << __LINE__ << endl;

  if(corrStack && data) showBtagMultiplicity(corrStack,       data, header, "corpred"+corrDist);
  cout << __LINE__ << endl;

  if(corrStackWithDD && data) showBtagMultiplicity(corrStackWithDD, data, header, "ddcorpred"+corrDist);
    cout << __LINE__ << endl;

}


//
void drawObservedBtagDistribution(TString url,TString tagger)
{
  toSave.Clear();
  nTagsPerCat=5;
  TFile *inF=TFile::Open(url);
  drawDistributionsFrom(inF,tagger+"btagsextendedcorr",tagger+"btagsextended");
  inF->Close();
  saveDistributions(tagger);
}


//
void drawObservedDiffBtagDistribution(TString url,TString tagger)
{
  outDir="diffbtag/"+tagger;

  toSave.Clear();
  nTagsPerCat=3;
  TFile *inF=TFile::Open(url);
  TH1F *kincats1H=(TH1F *)inF->Get("t#bar{t}/btvkincats1");
  TH1F *kincats2H=(TH1F *)inF->Get("t#bar{t}/btvkincats2");
  if(kincats1H==0){
    cout << "File does not seem to have the KIN categories defined" << endl;
    inF->Close();
    return;
  }
  
  Float_t rescale=kincats1H->GetBinContent(1)/30.; //a "feature" to be corrected in the filling of these histos

  Int_t ncats=kincats1H->GetXaxis()->GetNbins();
  for(Int_t icat=1; icat<=ncats; icat++)
    {
      char buf[200];
      sprintf(buf,"[p_{T}^{1}>%3.0f,p_{T}^{2}>%3.0f]",kincats1H->GetBinContent(icat)/rescale,kincats2H->GetBinContent(icat)/rescale);
      catLabel=buf;

      TString tagCountHisto("btvkin"); tagCountHisto+=icat; tagCountHisto+=tagger;
      drawDistributionsFrom(inF,tagCountHisto);
    }
  inF->Close();
  saveDistributions(tagger);
}



