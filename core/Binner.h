////////////////////////////////////////////////////////////////
///
///Class:               Binner
///Description:
///           

#pragma once

#include "Bins.h"
#include "Setup.h"
#include <map>
#include <utility>

//#pragma link C++ class std::map<TString, strings_t>+;
//#pragma link C++ class std::map<TString, TString>+;

namespace HS{
  namespace FIT{

    using strings_t = std::vector<TString>;

    class Binner  {
      
    public:
      Binner(Setup& setup);
      Binner() =default;
      Binner(const Binner&)=default;
      Binner(Binner&&)=default;
      virtual ~Binner()=default;
      Binner& operator=(const Binner& other)=default;
      Binner& operator=(Binner&& other) = default;

      void ReloadData(const TString& fname,const TString& name);
      void SplitData(TTree* tree,const TString& name="Data");
      void SplitData(const TString& tname,TString fname,const TString& name="Data");
      // void LoadAuxVar(TString vname);
      void LoadBinVar(TString opt,Int_t nbins,Double_t min,Double_t max);
      void LoadBinVar(const TString& opt,Int_t nbins,Double_t* xbins);

      void AddCut(TString cut){fSelection=std::move(cut);}
      TString GetCut(){return fSelection;}
      void RemoveAllCuts(){fSelection=TString();}
      
      void KeepBranch(TString name){fKeepBranches.push_back(name);}

      void LoadSetup(Setup &setup);

      const TString TreeName(const TString& name="Data") const {
	if(fNameToTree.count(name))
	  return fNameToTree.at(name);
	std::cout<<"Binner::TreeName files  named "<<name<<" not found"<<std::endl;
	return TString();

      }
      const strings_t FileNames(const TString& name="Data") const {
	if(fNameToFiles.count(name))
	  return fNameToFiles.at(name);
	std::cout<<"Binner::FileNames files  named "<<name<<" not found"<<fNameToFiles.size()<<std::endl;
	return strings_t();
      }

      const TString BinName(UInt_t i);
      
      UInt_t GetSize() const{return fBinNames.size();}
      void InitBins();
      
      Bool_t IsSetup() const{return fIsSetup;}

      Bins& GetBins(){return fBins;}

    protected:
      
    private:
      Bins fBins;

      strings_t fBinNames;
      strings_t fVarNames;
      std::map<TString, strings_t> fNameToFiles;
      std::map<TString, TString> fNameToTree;
      std::vector<TString> fKeepBranches;

      TString fOutDir;
      TString fSelection;
      
      Bool_t fIsSetup=kFALSE;
      
    };
    
  }//namespace FIT
}//namespace HS

