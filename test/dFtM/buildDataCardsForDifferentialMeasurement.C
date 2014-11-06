#include "TFile.h"
#include "TString.h"
#include "TH1F.h"
#include "TSystem.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

using namespace std;

TString outDir="diffbtag/";

void buildDataCardsForDifferentialMeasurement(TString dataUrl, TString mcUrl, TString systUrl,TString tagger);

//
void buildDataCardsForDifferentialMeasurement(TString dataUrl, TString mcUrl, TString systUrl, TString tagger)
{

  outDir="diffbtag/"+tagger;
  if(outDir!="" && outDir!="./") gSystem->Exec("mkdir -p " + outDir);

  TFile *mcF=TFile::Open(mcUrl);
  TH1F *kincats1H=(TH1F *)mcF->Get("t#bar{t}/btvkincats1");
  TH1F *kincats2H=(TH1F *)mcF->Get("t#bar{t}/btvkincats2");
  if(kincats1H==0){
    cout << "File does not seem to have the KIN categories defined" << endl;
    mcF->Close();
    return;
  }
  Int_t nCats(kincats1H->GetXaxis()->GetNbins());
  //Int_t nCats=4;

  Int_t curCat(0);
  Float_t lastKinVal(-1);
  std::map<Float_t,Int_t> jetCats;
  for(Int_t xbin=1; xbin<=nCats; xbin++){
    Float_t kinVal=kincats1H->GetBinContent(xbin);
    if(lastKinVal<0) { lastKinVal=kinVal; curCat++; }
    else if(kinVal != lastKinVal) { lastKinVal=kinVal; curCat++; }
    jetCats[kinVal]=curCat;
  }
  std::vector< std::pair<Int_t,Int_t> > evCatToJetCat;
  for(Int_t xbin=1; xbin<=nCats; xbin++){
    Float_t j1KinVal=kincats1H->GetBinContent(xbin);
    Float_t j2KinVal=kincats2H->GetBinContent(xbin);
    evCatToJetCat.push_back( std::pair<Int_t,Int_t>( jetCats[j1KinVal], jetCats[j2KinVal] ) );
  }
  
  cout << "Jet kinematics will be mapped in categories as: " << endl;
  for(std::vector<std::pair<Int_t,Int_t> >::iterator jIt=evCatToJetCat.begin(); jIt!=evCatToJetCat.end(); jIt++)
      cout << "\t(" << jIt->first << "," << jIt->second << ")" << endl;

  cout << "Building datacards for " << tagger  << " with " << nCats << " categories" << endl;

  TString channels[]={"ee","mumu","emu"};
  const size_t nchannels(sizeof(channels)/sizeof(TString));

  //
  TFile *systF=TFile::Open(systUrl);
  std::map<TString,std::pair<TString,TString> > signalSysts;
  if(systF){
    signalSysts["q2"]   = std::pair<TString,TString>("t#bar{t}systq2up","t#bar{t}systq2down");
    //signalSysts["meps"] = std::pair<TString,TString>("t#bar{t}systmepsup","t#bar{t}systmepsdown");
    //signalSysts["sig"]  = std::pair<TString,TString>("t#bar{t}systpowhegpy","");
    //signalSysts["ue"]   = std::pair<TString,TString>("t#bar{t}systtunep11","t#bar{t}systtunep11tev");
    //signalSysts["cr"]   = std::pair<TString,TString>("t#bar{t}systtunep11","t#bar{t}systtunep11nocr");
    //signalSysts["had"]  = std::pair<TString,TString>("t#bar{t}systpowhegpy","t#bar{t}systpowheghw");
  }

  //
  // TAGGING EFFICIENCIES
  //
  ofstream tagEffJson;
  tagEffJson.open ((outDir+"/"+tagger+"_eff.json").Data());
  tagEffJson << std::setprecision(3) << "{" << endl
	     << "  \"tagger\":\"" << tagger << "\","<< endl;
  std::map<TString,std::stringstream *> absepsCards;
  absepsCards["b"]=new std::stringstream; *(absepsCards["b"]) << "  \"absepsb\":{\n";
  absepsCards["c"]=new std::stringstream; *(absepsCards["c"]) << "  \"absepsc\":{\n";
  absepsCards["l"]=new std::stringstream; *(absepsCards["l"]) << "  \"absepsl\":{\n";
  Int_t kinCat(0);
  for(Int_t icat=1; icat<nCats; icat++)
    {
      //histogram name
      TString hname("btvflavcountskin"); hname += (icat-1);

      //nominal prediction
      TH1F *centralEff=(TH1F *) mcF->Get("t#bar{t}/"+hname+tagger);
      if(centralEff==0) continue;
      centralEff->Divide( (TH1F *) mcF->Get("t#bar{t}/"+hname) );

      //add signal sytematics
      std::map<TString, TH1F *> signalSystsH;
      for(std::map<TString,std::pair<TString,TString> >::iterator systIt=signalSysts.begin();
	  systIt != signalSysts.end();
	  systIt++)
	{
	  TH1F * effVarUp=(TH1F *)systF->Get(systIt->second.first+"/"+hname+tagger)->Clone("effvarup");
	  if(effVarUp==0){
	    cout << systIt->first << " does not seem to have valid " << hname << " for up unc." << endl;
	    continue;
	  }
	  effVarUp->Divide( (TH1F *)systF->Get(systIt->second.first+"/"+hname) );

	  TH1F * effVarDown=0;
	  if(systIt->second.second!=""){
	    effVarDown=(TH1F *)systF->Get(systIt->second.second+"/"+hname+tagger)->Clone("effvardown");
	    effVarDown->Divide( (TH1F *)systF->Get(systIt->second.second+"/"+hname) );
	  }else{
	    effVarDown=(TH1F *)centralEff->Clone("effvardown");
	  }
	  if(effVarDown==0){
	    cout << systIt->first << " does not seem to have valid " << hname << " for down unc." << endl;
	    continue;
	  }

	  TH1F *effDiff=(TH1F *)effVarUp->Clone("effdiff");
	  effDiff->Add(effVarDown,-1);
	  TH1F *effSum=(TH1F *)effVarUp->Clone("effsum");
	  effSum->Add(effVarDown);

	  TH1F *relDiff=(TH1F *)effDiff->Clone(systIt->first);
	  relDiff->Divide(effSum);
	  relDiff->SetDirectory(0);
	  signalSystsH[systIt->first]=relDiff;

	  //free up unused mem
	  delete effVarUp;
	  delete effVarDown;
	  delete effDiff;
	  delete effSum;
	}

      //now report
      kinCat++;
      for(Int_t iflav=1; iflav<=3; iflav++)
	{
	  
	  Float_t val=centralEff->GetBinContent(iflav);
	  Float_t mcstat=centralEff->GetBinError(iflav)/val;
	  TString flav(centralEff->GetXaxis()->GetBinLabel(iflav));
	  *(absepsCards[flav]) << "    \"k_" << kinCat << "\":{"
			       << "\"val\":" << val
			       << ",\"mcstat\":" << mcstat << flush;  
	  for(std::map<TString, TH1F *>::iterator systIt=signalSystsH.begin();
	      systIt!=signalSystsH.end();
	      systIt++)
	    {
	      //Float_t relDiff=TMath::Abs(systIt->second->GetBinContent(iflav));
	      Float_t relDiff=systIt->second->GetBinContent(iflav);
	      if(fabs(relDiff)<1e-3 || isnan(relDiff)) continue;
	      *(absepsCards[flav]) << ",\""+systIt->first << "\":" << relDiff << flush;
	    }
	  
	  *(absepsCards[flav]) << "}," << endl;
	}
    }
  tagEffJson << absepsCards["b"]->str() << "  }," << endl 
	     << absepsCards["c"]->str() << "  }," << endl 
	     << absepsCards["l"]->str() << "  }" << endl 
	     << "}" << endl;
  tagEffJson.close();

  cout << "Expected tagging efficiencies have been stored in " << outDir << "/" + tagger << "_eff.json" << endl;


  //
  // FLAVOUR COMPONENTS
  //
  cout << "Building datacards for " << tagger  << " with " << nCats << " categories" << endl;
  ofstream flavCompJson;
  flavCompJson.open ((outDir+"/flavbreakup.json").Data());
  flavCompJson << std::setprecision(3) << "{" << endl;
  
  //map the categories to the jet kinematics
  Int_t icat(1);
  flavCompJson << "  \"catmap\":{" << endl
	       << "    " << flush;
  for(std::vector<std::pair<Int_t,Int_t> >::iterator jIt=evCatToJetCat.begin(); jIt!=evCatToJetCat.end(); jIt++,icat++)
    flavCompJson << "\"k" << icat << "\":\"(" << jIt->first << "," << jIt->second << ")\"," << flush;
  flavCompJson << endl << "  }," << endl;

  //now save the f_ij coefficients
  std::map<TString,std::stringstream *> fijCards;
  for(size_t ich=0; ich<nchannels; ich++)
    {
      for(icat=1; icat<nCats; icat++)
	{
	  //histogram name
	  TString hname(channels[ich]+"_btvfijkin"); hname += (icat-1);
	  
	  //nominal prediction
	  TH1F *centralExp=(TH1F *) mcF->Get("t#bar{t}/"+hname);
	  if(centralExp==0) continue;
	  Float_t totalExp=centralExp->Integral();
	  if(totalExp==0) continue;
	  centralExp->Scale(1./totalExp);
	  
	  //add signal sytematics
	  std::map<TString, TH1F *> signalSystsH;
	  for(std::map<TString,std::pair<TString,TString> >::iterator systIt=signalSysts.begin();
	      systIt != signalSysts.end();
	      systIt++)
	    {
	      TH1F * varUp=(TH1F *)systF->Get(systIt->second.first+"/"+hname);
	      if(varUp==0){
		cout << systIt->first << " does not seem to have valid " << hname << " for up unc." << endl;
		continue;
	      }
	      varUp=(TH1F *)varUp->Clone("varup");
	      Float_t totalExpUp(varUp->Integral());
	      if(totalExpUp>0) varUp->Scale(1./totalExpUp);
 
	      TH1F * varDown=0;
	      if(systIt->second.second!=""){
		varDown=(TH1F *)systF->Get(systIt->second.second+"/"+hname);
		if(varDown!=0){
		  varDown=(TH1F *)varDown->Clone("vardown");
		  Float_t totalExpDown(varDown->Integral());
		  if(totalExpDown>0) varDown->Scale(1./totalExpDown);
		}
	      }else{
		varDown=(TH1F *)centralExp->Clone("vardown");
	      }
	      if(varDown==0){
		cout << systIt->first << " does not seem to have valid " << hname << " for down unc." << endl;
		continue;
	      }

	      TH1F *diff=(TH1F *)varUp->Clone("diff");
	      diff->Add(varDown,-1);
	      TH1F *sum=(TH1F *)varUp->Clone("sum");
	      sum->Add(varDown);
	      
	      TH1F *relDiff=(TH1F *)diff->Clone(systIt->first);
	      relDiff->Divide(sum);
	      relDiff->SetDirectory(0);
	      signalSystsH[systIt->first]=relDiff;
	      
	      //free up unused mem
	      delete varUp;
	      delete varDown;
	      delete diff;
	      delete sum;
	    }

	  //now report
	  for(Int_t iflav=1; iflav<=centralExp->GetXaxis()->GetNbins(); iflav++)
	    {
	      Float_t val=centralExp->GetBinContent(iflav);
	      Float_t mcstat=centralExp->GetBinError(iflav)/val;
	      
	      TString tag(centralExp->GetXaxis()->GetBinLabel(iflav));
	      if(fijCards.find(tag)==fijCards.end()) {
		fijCards[tag]=new std::stringstream;
		*(fijCards[tag]) << "  \"f" << tag << "\":{" << endl;
	      }

	      //central value
	      *(fijCards[tag]) << "    \"" << channels[ich] << "k" << icat << "\":{"
			       << "\"val\":" << val << flush;

	      //if not null add uncertainties as well
	      if(val!=0){
		*(fijCards[tag]) << ",\"mcstat\":" << mcstat << flush;  
		
		for(std::map<TString, TH1F *>::iterator systIt=signalSystsH.begin();
		    systIt!=signalSystsH.end();
		    systIt++)
		  {
		    //Float_t relDiff=TMath::Abs(systIt->second->GetBinContent(iflav));
		    Float_t relDiff=systIt->second->GetBinContent(iflav);
		    if(fabs(relDiff)<1e-3 || isnan(relDiff)) continue;
		    *(fijCards[tag]) << ",\""+systIt->first << "\":" << relDiff << flush;
		  }
	      }
	      
	      *(fijCards[tag]) << "}," << endl;
	    }
	}
    }
  for(std::map<TString,std::stringstream *>::iterator it = fijCards.begin();
      it!=fijCards.end();
      it++)
    flavCompJson << it->second->str() << "  }," << endl;
  flavCompJson << "}" << endl;
  flavCompJson.close();

  cout << "Expected flavour breakup per channel has been stored in " << outDir << "/flavbreakup.json" << endl;
}

