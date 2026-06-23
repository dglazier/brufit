/**
 * @file BruEventsPDF.cpp
 */

#include "BruEventsPDF.h"
#include <RooRealVar.h>
#include <RooCategory.h> 
#include <RooRandom.h>
#include <RooDataSet.h> 
#include <RooMsgService.h> 
#include <RooGaussian.h>
#include <RooPlot.h>
#include <RooHist.h>
#include <RooPullVar.h> 
#include <RooFitResult.h>

#include <TMath.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TEntryList.h>
#include <TError.h>

#include <algorithm> 
#include <random>
#include <iostream>

namespace bru {

    Bool_t BruEventsPDF::BruEventsPDF_IsPlotting = kFALSE;

    void BruEventsPDF::SetIsPlotting(Bool_t is) {
        BruEventsPDF_IsPlotting = is;	
    }

    BruEventsPDF::BruEventsPDF(const BruEventsPDF& other, const char* name) : RooAbsPdf(other, name) {
        _IsClone = kTRUE;
        _Parent = const_cast<BruEventsPDF*>(&other);

        _DataCache = other._DataCache; 
        
        _AssertPosDataReal = other._AssertPosDataReal;
        _AssertPosDataCats = other._AssertPosDataCats;
        _Napd = other._Napd;

        _NInt = other._NInt;
        _Geni = other._Geni;
        _TruthPrefix = other._TruthPrefix;
        _ForceConstInt = other._ForceConstInt;
        _ForceNumInt = other._ForceNumInt;
        _ConstInt = other._ConstInt;
        _CheckInt = other._CheckInt;
        _UseWeightsGen = other._UseWeightsGen;
        _Cut = other._Cut;
        _InWeightCut = other._InWeightCut;
        _IsValid = other._IsValid;
        _WgtsConf = other._WgtsConf;
        _HistIntegrals = other._HistIntegrals;
        _MaxValue = other._MaxValue;
        _IntRangeLow = other._IntRangeLow;
        _IntRangeHigh = other._IntRangeHigh;

        _ProtoRealVars = other._ProtoRealVars;
        _ProtoCatVars = other._ProtoCatVars;

        _GeneratedIndices = other._GeneratedIndices;
        _GeneratedWeights = other._GeneratedWeights;

        _Last = other._Last;
        _LastLength = other._LastLength;
    }

    BruEventsPDF::~BruEventsPDF() {
        if (!_IsClone) {
            if (!_GeneratedIndices.empty()) {
                TFile entryFile(TString("entryFile_") + GetName() + ".root", "recreate");
                TEntryList elist("GenEvents", "Generated Event Indices");
                for (auto idx : _GeneratedIndices) elist.Enter(idx);
                elist.Write();
                entryFile.Close();
            }
            if (!_GeneratedWeights.empty()) {
                HS::FIT::Weights wgts("genWeights");
                wgts.SetSpecies(GetName());
                wgts.SetFile(TString(GetName()) + "Weights.root");
                for (size_t i = 0; i < _GeneratedWeights.size(); ++i) {
                    wgts.FillWeight(i, _GeneratedWeights[i]);
                }
                wgts.Save();
            }
        }
        for (auto &i : _VarSet) delete i;
        _VarSet.clear();
    }

    void BruEventsPDF::InitSets() {
        _Npars = _ParSet.size();
        _Nvars = _ProxSet.size();
        _Ncats = _CatSet.size();
        _LastLength = _Npars + 1;
        _Last.assign(_Npars + 1, 100.0);
    }

    RooArgSet BruEventsPDF::VarSet(Int_t iset) const {
        RooArgSet aset(Form("VarSet_%d", iset));
        if (iset == 0) {
            for (auto i : _ProxSet) aset.add(i->arg());
            for (auto i : _CatSet) aset.add(i->arg());
        } else {
            for (auto j : _ProxSet) {
                if (_ProxSet[iset - 1]->GetName() != j->GetName()) aset.add(j->arg());
            }
            for (auto i : _CatSet) aset.add(i->arg());
        }
        return aset;
    }

    void BruEventsPDF::SetCache(std::shared_ptr<const MCEventCache> cache) {
        _DataCache = cache;
        if (_DataCache) _ConstInt = _DataCache->_NTreeEntries;
    }

    Bool_t BruEventsPDF::SetEvTree(TTree* tree, TString cut, TTree* MCGenTree) {
      TDirectory* saveDir = gDirectory; // 1. LOCK DIRECTORY
        if (!tree || !tree->GetEntries()) return kFALSE;
        if (cut == TString()) _Cut = _InWeightCut;
        else if (_InWeightCut == TString()) _Cut = cut;
        else _Cut = cut + "&&" + _InWeightCut;

        std::vector<TString> varNames, catNames;
        
        _ProtoRealVars.clear();
        for (auto p : _ProxSet) {
            varNames.push_back(p->GetName());
            if (!tree->GetBranch(p->GetName())) _ProtoRealVars.push_back(p->GetName());
        }
        
        _ProtoCatVars.clear();
        for (auto c : _CatSet) {
            catNames.push_back(c->GetName());
            if (!tree->GetBranch(c->GetName())) _ProtoCatVars.push_back(c->GetName());
        }

        auto newCache = std::make_shared<MCEventCache>();
        _BranchStatus = newCache->LoadTree(tree, _Cut, varNames, catNames, _TruthPrefix, _WgtsConf, MCGenTree);

        _ConstInt = newCache->_NTreeEntries;
        _DataCache = newCache; 
        _IsValid = _BranchStatus;
	if (saveDir) saveDir->cd(); // 2. RESTORE DIRECTORY
	
        return _BranchStatus;
    }

    Int_t BruEventsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const {
        if (!_DataCache) return 0; 
        if (matchArgs(directVars, generateVars, VarSet(0))) return 1;
        return 0;
    }

    void BruEventsPDF::initIntegrator() const {}

    void BruEventsPDF::initGenerator(Int_t code) {
        if (!_DataCache) return;
        if (_Parent->GetMaxValue() == 0 || _Parent->CheckChange()) {	
            Double_t value = 0;		
            if (code == 1) {	
                _MaxValue = 0;
                for (Int_t i = 0; i < _DataCache->_NTreeEntries; i++) {
                    _TreeEntry = i;
                    value = evaluateMC(&_DataCache->_vecRealGen, &_DataCache->_vecCatGen);
                    if (value < 0) { std::cout << " BruEventsPDF::initGenerator -ve intensity !!! " << value << std::endl; exit(0); }
                    if (value > _MaxValue) _MaxValue = value * 1.01;
                }
                _Parent->SetMaxValue(_MaxValue);
            }
        }		
        _GeneratedIndices.clear();
        _GeneratedWeights.clear();
    }

    void BruEventsPDF::generateEvent(Int_t code) {
        Double_t value = 0;
        if (!_DataCache) return;
        
        if (!_UseWeightsGen) {
            while (_Geni < _DataCache->_NTreeEntries) {
                _TreeEntry = IncrementGeni();
                if (!CheckRange("")) continue;
                value = evaluateMC(&_DataCache->_vecRealGen, &_DataCache->_vecCatGen); 
                if (value > _MaxValue * RooRandom::uniform()) {
                    for (Int_t i = 0; i < _Nvars; i++) *(_ProxSet[i]) = _DataCache->_vecReal[_TreeEntry * _Nvars + i]; 
                    for (Int_t i = 0; i < _Ncats; i++) *(_CatSet[i]) = _DataCache->_vecCat[_TreeEntry * _Ncats + i];
                    _GeneratedIndices.push_back(_DataCache->_TreeEntryNumber[_TreeEntry]);
                    return;
                }
            }
        } else {
            while (_Geni < _DataCache->_NTreeEntries) {
                _Parent->SetGeni(_Geni);
                _TreeEntry = _Geni++;
                if (!CheckRange("")) continue;
                value = evaluateMC(&_DataCache->_vecRealGen, &_DataCache->_vecCatGen); 
                for (Int_t i = 0; i < _Nvars; i++) *(_ProxSet[i]) = _DataCache->_vecReal[_TreeEntry * _Nvars + i]; 
                for (Int_t i = 0; i < _Ncats; i++) *(_CatSet[i]) = _DataCache->_vecCat[_TreeEntry * _Ncats + i];
                _GeneratedWeights.push_back(value);
                _GeneratedIndices.push_back(_Geni - 1);
                return;
            }
        }
        Fatal("BruEventsPDF::generateEvent", "Ran out of events at %lld", _Geni);
    }

    // ORIGINAL LOGIC EXACTLY PRESERVED
    Int_t BruEventsPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
        if (_ForceNumInt) return 0;
        if (!_DataCache && !_ForceConstInt) return 0;

        if (_ProxSet.size() == 1 && _CatSet.size() == 0) {
            if (matchArgs(allVars, analVars, VarSet(0))) {
                if (BruEventsPDF_IsPlotting && _HistIntegrals.size() == 0)
                    HistIntegrals(rangeName);
                return 1;
            }
        } else {
            for (UInt_t i = 0; i < 1 + _ProxSet.size(); i++) {
                if (!_DataCache && _ForceConstInt && i == 0) { return 1; }
                else if (!_DataCache) return 0;
                if (matchArgs(allVars, analVars, VarSet(i))) { return i + 1; }
            }
        }
        return 0; 
    }

    Double_t BruEventsPDF::analyticalIntegralForSampling(const char* rangeName) const {
      Double_t integral=0;
      Long64_t accepted=0;
      Long64_t all=0;
      Long64_t ilow=0;
      Long64_t ihigh=0;

      SetLowHighVals(ilow,ihigh); 
      
      if(CheckChange()){
        std::vector<double> values(ihigh-ilow);
        for(Long64_t ie=ilow;ie<ihigh;ie++){
          _TreeEntry=ie;
          if(!CheckRange(rangeName)){
            values[all]=0;
            ++all;
            continue;
          }
          values[all]=evaluateMC(&_DataCache->_vecReal,&_DataCache->_vecCat)*GetIntegralWeight(ie);
          integral+=values[all];
          ++accepted; 
          ++all; 
        }
        
        integral/=accepted;
        
        double sum_of_diffs = 0.;
        std::for_each(values.begin(), values.end(), [&sum_of_diffs,&integral] (double n) {
            double term = (n-integral);
            sum_of_diffs += term*term;
        });
        
        if (accepted > 1) {
            _SigmaIntegral = TMath::Sqrt(sum_of_diffs / (accepted - 1)) / TMath::Sqrt(accepted);
        } else {
            _SigmaIntegral = 0.0;
        }
        _Last[0]= integral;
      }      

      return sampleIntegral(_Last[0],_SigmaIntegral);
    }

    Double_t BruEventsPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
        if(code==1&&_ForceConstInt&&!_DataCache) {_Last[0]=1;return _Last[0];}
        Long64_t NEv=0;
        Double_t integral=0.;
       
        if(code==1)
            if(!CheckChange()) return _Last[0];
        
        if(code==1){
            auto check= AssertPositivePDF();
            if(check==kFALSE) return _Last[0]=0;

            Long64_t accepted=0;
            Long64_t ilow=0;
            Long64_t ihigh=0;

            SetLowHighVals(ilow,ihigh); 
            for(Long64_t ie=ilow;ie<ihigh;ie++){
                _TreeEntry=ie;
                if(!CheckRange(rangeName)) continue;
                accepted++;
                integral+=evaluateMC(&_DataCache->_vecReal,&_DataCache->_vecCat)*GetIntegralWeight(ie);
            }
        
            if (accepted > 0) integral/=accepted;
            else integral = 0;
        }
        else {
            if(_HistIntegrals.size()==0)
                HistIntegrals(rangeName);
        
            Int_t vindex=code-2;
            Double_t vval=*(_ProxSet[vindex]);
            integral=_HistIntegrals[vindex].Interpolate(vval);
            if(integral<0) integral=0;
            return integral;
        }
        _Last[0]=integral;
        return _Last[0];
    }

    Double_t BruEventsPDF::unnormalisedIntegral(Int_t code, const char* rangeName) const {
      Double_t integral=0;
      Double_t nev=0;
      Double_t nMC=0;
      if(!_DataCache) return 0;
      
      if(code==1){
        for(Long64_t ie=0;ie<_DataCache->_NTreeEntries;ie++){
          _TreeEntry=ie;
          if(!CheckRange(TString(rangeName).Data())) continue;
          integral+=evaluateMC(&_DataCache->_vecReal,&_DataCache->_vecCat)*GetIntegralWeight(ie);
          nev++;
        }
      }
      else if(code==2 && _DataCache->_HasMCGenTree){
        for(Long64_t ie=0;ie<_DataCache->_NMCGenTreeEntries;ie++){
          _TreeEntry=ie;
          integral+=evaluateMC(&_DataCache->_vecRealMCGen,&_DataCache->_vecCatMCGen);
          nMC++;
        }
      }
      else{
        return 0;
      }
      return integral;
    }
    
    void BruEventsPDF::HistIntegrals(const char* rangeName) const {
      if(!_DataCache) return;
      Long64_t ilow=0;
      Long64_t ihigh=0;
      SetLowHighVals(ilow,ihigh);
      for(Int_t i=0;i<_Nvars;i++){
        auto  arg=dynamic_cast<const RooRealVar*>(&_ProxSet[i]->arg());
        if(arg)
          _HistIntegrals.emplace_back(arg->GetName(),arg->GetName(),arg->getBins(),arg->getMin(),arg->getMax());
      }
      Long64_t accepted=0;
      for(Int_t ie=ilow;ie<ihigh;ie++){
        _TreeEntry=ie;
        if(!CheckRange(TString(rangeName).Data())){continue;}
        accepted++;
        Double_t value=evaluateMC(&_DataCache->_vecReal,&_DataCache->_vecCat)*GetIntegralWeight(ie);
        for(Int_t vindex=0;vindex<_Nvars;vindex++){
          _HistIntegrals[vindex].Fill(_DataCache->_vecReal[_TreeEntry*_Nvars+vindex],value/_HistIntegrals[vindex].GetBinWidth(1));
        }
      }
      for(Int_t vindex=0;vindex<_Nvars;vindex++) {
        if(accepted > 0) _HistIntegrals[vindex].Scale(1./accepted);
      }
  
      _Parent->SetHistIntegrals(_HistIntegrals);
    }

    void BruEventsPDF::SetLowHighVals(Long64_t& ilow, Long64_t& ihigh) const {
        ilow = 0; ihigh = 0;
        if (_Parent) {
            ilow = _Parent->GetIntRangeLow();
            ihigh = _Parent->GetIntRangeHigh();
        } else {
            ilow = GetIntRangeLow();
            ihigh = GetIntRangeHigh();
        }
        if (ihigh == 0 && _NInt > -1) ihigh = _NInt;
        else if (ihigh == 0 && _DataCache) ihigh = _DataCache->_NTreeEntries; 
        if (_DataCache && ihigh > (Long64_t)_DataCache->_NTreeEntries) ihigh = _DataCache->_NTreeEntries;
    }

    Bool_t BruEventsPDF::CheckRange(const char* rangeName) const {
      if(!_DataCache) return kFALSE;
      for(UInt_t i=0;i<_ProxSet.size();i++){
        auto var=(dynamic_cast<const RooRealVar*>(&(_ProxSet[i]->arg())));
        if(!var->inRange(_DataCache->_vecReal[_TreeEntry*_Nvars+i],TString(rangeName).Data())){return kFALSE;}
      }
      return kTRUE;
    }

    Bool_t BruEventsPDF::CheckChange() const {
      Bool_t hasChanged=false;
      for(Int_t i=1;i<_Npars+1;i++)
        if(_Last[i]!=(*(_ParSet[i-1]))){
          hasChanged=true;
        }
      if(hasChanged){
        for(Int_t i=1;i<_Npars+1;i++){
          _Last[i]=*(_ParSet[i-1]);
        }
      }
      return hasChanged;
    }

    void BruEventsPDF::CheckIntegralParDep(Int_t Ntests) {
        _CheckInt = Ntests;
        if (!_DataCache) return; 
    
        Long64_t saveNint = _NInt;
        _NInt = _DataCache->_NTreeEntries;
        Ntests = (Ntests * _ParSet.size());
    
        RooRealVar integral("integral", "integral", 0, 0, 2);
        if (_NInt > 0) integral.setError(sqrt(_NInt) / _NInt); 
        RooDataSet ds("intds", "intds", RooArgSet(integral));
        std::vector<Double_t> SavedPars;
        
        for (auto &ip : _ParSet) {
            auto par = (dynamic_cast<const RooRealVar*>(&(ip->arg())));
            SavedPars.push_back(par->getValV());
        }
        
        for (Int_t ir = 0; ir < Ntests; ir++) { 
            for (auto &ip : _ParSet) {
                auto par = (RooRealVar*)(&(ip->arg()));
                par->setVal((par->getMax("") - par->getMin("")) * RooRandom::uniform() + par->getMin(""));
            }
            integral.setVal(analyticalIntegral(1, ""));
            ds.add(RooArgSet(integral));
        }
    
        Double_t low = 0;
        Double_t high = 0;
        ds.getRange(integral, low, high);
        integral.setRange(low, high);
        RooPlot *frame = integral.frame();
        ds.plotOn(frame);
    
        frame->Draw();
    
        RooRealVar mean("mean", "mean", ds.mean(integral));
        RooRealVar pvar("IntPull", "Integral Pull Dist.", -5, 5);
        RooPullVar pull("IntPull", "Integral Pull Dist.", integral, mean);
        ds.addColumn(pull, kFALSE);
    
        ds.getRange(pvar, low, high);
        pvar.setRange(low, high);
        RooPlot *framePull = pvar.frame();
        ds.plotOn(framePull);
        RooRealVar mp("mp", "mp", 0, -5, 5);
        RooRealVar sp("sp", "sp", 1, 0, 100);
        RooGaussian gp("gp", "gp", pvar, mp, sp);
        gp.fitTo(ds);
        gp.paramOn(framePull);
        gp.plotOn(framePull);
    
        new TCanvas();
        framePull->Draw();
    
        _ConstInt = mean.getVal();
        _NInt = saveNint;
        if (sp.getVal() < 2) SetConstInt();
        for (UInt_t ip = 0; ip < _ParSet.size(); ip++) {
            auto par = (RooRealVar*)(&(_ParSet[ip]->arg()));
            par->setVal(SavedPars[ip]);
        }
        _CheckInt = kFALSE; 
    }

    Bool_t BruEventsPDF::AddProtoData(const RooDataSet* data) {
        if (!_DataCache || !_DataCache->_NTreeEntries) return kFALSE;
	TDirectory* saveDir = gDirectory; // 1. LOCK DIRECTORY
       
        auto mutableCache = std::make_shared<MCEventCache>(*_DataCache);
        const RooArgSet *dataVars = data->get();
        Long64_t Nentries = data->numEntries();
        std::vector<Long64_t> vrandom(Nentries);
        for (Long64_t ir = 0; ir < Nentries; ir++) vrandom[ir] = ir;
        std::shuffle(vrandom.begin(), vrandom.end(), std::mt19937(std::random_device()()));
    
        std::vector<Short_t> protoDataForVar, protoDataForCat;
        for (auto* arg : *dataVars) {
            if (TString("UID") == arg->GetName()) continue; 
            
            for (Int_t ip = 0; ip < _Nvars; ip++) {
                if (TString(arg->GetName()) == TString(_ProxSet[ip]->GetName())) {
                    if (std::find(_ProtoRealVars.begin(), _ProtoRealVars.end(), TString(arg->GetName())) != _ProtoRealVars.end()) {
                        protoDataForVar.push_back(ip);
                    }
                }
            }
            for (Int_t ip = 0; ip < _Ncats; ip++) {
                if (TString(arg->GetName()) == TString(_CatSet[ip]->GetName())) {
                    if (std::find(_ProtoCatVars.begin(), _ProtoCatVars.end(), TString(arg->GetName())) != _ProtoCatVars.end()) {
                        protoDataForCat.push_back(ip);
                    }
                }
            }
        }
    
        Long64_t idata = 0;
        if (!(protoDataForVar.size() + protoDataForCat.size())) {
            return kTRUE; 
        }
    
        for (Long64_t id = 0; id < mutableCache->_NTreeEntries; id++) {
            dataVars = data->get(vrandom[idata]);
            for (short ip : protoDataForVar) {
                Double_t val = dataVars->getRealValue(_ProxSet[ip]->GetName());
                mutableCache->_vecReal[id * _Nvars + ip] = val;
                mutableCache->_vecRealGen[id * _Nvars + ip] = val;
            }  
            for (short ip : protoDataForCat) {
                Int_t val = dataVars->getCatIndex(_CatSet[ip]->GetName());
                mutableCache->_vecCat[id * _Ncats + ip] = val;
                mutableCache->_vecCatGen[id * _Ncats + ip] = val;     
            }
            
            if(idata==(Long64_t)vrandom.size()-1){
                std::shuffle(vrandom.begin(),vrandom.end(), std::mt19937(std::random_device()()));
                idata=0;
            } else {
                idata++;
            }
        }
        
        _DataCache = mutableCache;
	if (saveDir) saveDir->cd(); // 2. RESTORE DIRECTORY
       
	return kTRUE;  
    }

    void BruEventsPDF::SetNextRange(Int_t ir) {
        Long64_t Nentries = _DataCache ? _DataCache->_NTreeEntries : 0;
        Int_t range = ((Double_t)Nentries) / _NRanges;
        _IntRangeLow = ir * range;
        _IntRangeHigh = (ir + 1) * range;
    }

    void BruEventsPDF::MakeAssertPostiveData() {
        auto saveTreeEntry = _TreeEntry;
        _TreeEntry = 0;

        auto NVars = _ProxSet.size();
        _AssertPosDataReal.resize(_Napd * NVars);

        auto NCats = _CatSet.size();
        _AssertPosDataCats.resize(_Napd * NCats);

        for (Long64_t iapd = 0; iapd < _Napd; ++iapd) {
            UInt_t ivar = 0;
            for (auto v : _ProxSet) {
                auto vr = dynamic_cast<const RooRealVar*>(&v->arg());
                if (vr != nullptr) {
                    _AssertPosDataReal[_TreeEntry * NVars + ivar] = gRandom->Uniform(vr->getMin(""), vr->getMax(""));
                    ++ivar;
                }
            }
            UInt_t icat = 0;
            for (auto v : _CatSet) {
                auto vc = dynamic_cast<const RooCategory*>(&v->arg());
                if (vc != nullptr) {
                    auto catstate = gRandom->Integer(vc->size());
                    auto val = vc->getOrdinal(catstate).second;
                    _AssertPosDataCats[_TreeEntry * NCats + icat] = val;
                    ++icat;
                }
            }
            _TreeEntry++;
        }
        _TreeEntry = saveTreeEntry;
    }
   
    Bool_t BruEventsPDF::AssertPositivePDF() const {
        InitAssertPositiveCheck();
        auto saveTreeEntry = _TreeEntry;
        _TreeEntry = 0;
        for (Long64_t iapd = 0; iapd < _Napd; ++iapd) {
            auto val = evaluateMC(&_AssertPosDataReal, &_AssertPosDataCats);
            ++_TreeEntry;
            if (val < -1E-4) { 
                logEvalError("BruEventsPDF::AssertPositivePDF() PDF cannot be -ve...");
                _TreeEntry = saveTreeEntry;
                FinishAssertPositiveCheck(); 
                return kFALSE;
            }
        }
        _TreeEntry = saveTreeEntry;
        FinishAssertPositiveCheck();
        return kTRUE;
    }

} // namespace bru
// /**
//  * @file BruEventsPDF.cpp
//  */

// #include "BruEventsPDF.h"
// #include <RooRealVar.h>
// #include <RooCategory.h> 
// #include <RooRandom.h>
// #include <RooDataSet.h> 
// #include <RooMsgService.h> 
// #include <RooGaussian.h>
// #include <RooPlot.h>
// #include <RooHist.h>
// #include <RooPullVar.h> 
// #include <RooFitResult.h>

// #include <TMath.h>
// #include <TCanvas.h>
// #include <TFile.h>
// #include <TEntryList.h>
// #include <TError.h>

// #include <algorithm> 
// #include <random>
// #include <iostream>

// namespace bru {

//     Bool_t BruEventsPDF::BruEventsPDF_IsPlotting = kFALSE;

//     void BruEventsPDF::SetIsPlotting(Bool_t is) {
//         BruEventsPDF_IsPlotting = is;	
//     }

//     BruEventsPDF::BruEventsPDF(const BruEventsPDF& other, const char* name) : RooAbsPdf(other, name) {
//         _IsClone = kTRUE;
//         _Parent = const_cast<BruEventsPDF*>(&other);

//         _DataCache = other._DataCache; 

//         _AssertPosDataReal = other._AssertPosDataReal;
//         _AssertPosDataCats = other._AssertPosDataCats;
//         _Napd = other._Napd;

//         _NInt = other._NInt;
//         _Geni = other._Geni;
//         _ForceConstInt = other._ForceConstInt;
//         _ForceNumInt = other._ForceNumInt;
//         _ConstInt = other._ConstInt;
//         _CheckInt = other._CheckInt;
//         _UseWeightsGen = other._UseWeightsGen;
//         _HistIntegrals = other._HistIntegrals;
//         _MaxValue = other._MaxValue;
//         _IntRangeLow = other._IntRangeLow;
//         _IntRangeHigh = other._IntRangeHigh;

//         _TruthPrefix = other._TruthPrefix;
//         _WgtsConf = other._WgtsConf;
//         _Cut = other._Cut;
//         _InWeightCut = other._InWeightCut;
//         _IsValid = other._IsValid;
//         _BranchStatus = other._BranchStatus;
        
//         _ProtoRealVars = other._ProtoRealVars;
//         _ProtoCatVars = other._ProtoCatVars;

//         _GeneratedIndices = other._GeneratedIndices;
//         _GeneratedWeights = other._GeneratedWeights;

//         _Last = other._Last;
//         _LastLength = other._LastLength;
//     }

//     BruEventsPDF::~BruEventsPDF() {
//         if (!_IsClone) {
//             if (!_GeneratedIndices.empty()) {
//                 TFile entryFile(TString("entryFile_") + GetName() + ".root", "recreate");
//                 TEntryList elist("GenEvents", "Generated Event Indices");
//                 for (auto idx : _GeneratedIndices) elist.Enter(idx);
//                 elist.Write();
//                 entryFile.Close();
//                 std::cout << "BruEventsPDF saved " << _GeneratedIndices.size() << " generated events to entryFile_" << GetName() << ".root\n";
//             }
//             if (!_GeneratedWeights.empty()) {
//                 HS::FIT::Weights wgts("genWeights");
//                 wgts.SetSpecies(GetName());
//                 wgts.SetFile(TString(GetName()) + "Weights.root");
//                 for (size_t i = 0; i < _GeneratedWeights.size(); ++i) {
//                     wgts.FillWeight(i, _GeneratedWeights[i]);
//                 }
//                 wgts.Save();
//                 std::cout << "BruEventsPDF saved " << _GeneratedWeights.size() << " generated weights to " << GetName() << "Weights.root\n";
//             }
//         }
//         for (auto &i : _VarSet) delete i;
//         _VarSet.clear();
//     }

//     void BruEventsPDF::InitSets() {
//         _Npars = _ParSet.size();
//         _Nvars = _ProxSet.size();
//         _Ncats = _CatSet.size();
//         _LastLength = _Npars + 1;
//         _Last.assign(_Npars + 1, 100.0);
//     }

//     RooArgSet BruEventsPDF::VarSet(Int_t iset) const {
//         RooArgSet aset(Form("VarSet_%d", iset));
//         if (iset == 0) {
//             for (auto i : _ProxSet) aset.add(i->arg());
//             for (auto i : _CatSet) aset.add(i->arg());
//         } else {
//             for (auto j : _ProxSet) {
//                 if (_ProxSet[iset - 1]->GetName() != j->GetName()) aset.add(j->arg());
//             }
//             for (auto i : _CatSet) aset.add(i->arg());
//         }
//         return aset;
//     }

//     void BruEventsPDF::SetCache(std::shared_ptr<const MCEventCache> cache) {
//         _DataCache = cache;
//         if (_DataCache) _ConstInt = _DataCache->_NTreeEntries;
//     }

//     Bool_t BruEventsPDF::SetEvTree(TTree* tree, TString cut, TTree* MCGenTree) {
//         if (!tree || !tree->GetEntries()) return kFALSE;
//         if (cut == TString()) _Cut = _InWeightCut;
//         else if (_InWeightCut == TString()) _Cut = cut;
//         else _Cut = cut + "&&" + _InWeightCut;

//         std::vector<TString> varNames, catNames;
        
//         _ProtoRealVars.clear();
//         for (auto p : _ProxSet) {
//             varNames.push_back(p->GetName());
//             if (!tree->GetBranch(p->GetName())) _ProtoRealVars.push_back(p->GetName());
//         }
        
//         _ProtoCatVars.clear();
//         for (auto c : _CatSet) {
//             catNames.push_back(c->GetName());
//             if (!tree->GetBranch(c->GetName())) _ProtoCatVars.push_back(c->GetName());
//         }

//         auto newCache = std::make_shared<MCEventCache>();
//         _BranchStatus = newCache->LoadTree(tree, _Cut, varNames, catNames, _TruthPrefix, _WgtsConf, MCGenTree);

//         _ConstInt = newCache->_NTreeEntries;
//         _DataCache = newCache; 
//         _IsValid = _BranchStatus;
//         return _BranchStatus;
//     }

//     Int_t BruEventsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const {
//         if (!_DataCache) return 0; 
//         if (matchArgs(directVars, generateVars, VarSet(0))) return 1;
//         return 0;
//     }

//     void BruEventsPDF::initIntegrator() {}

//     void BruEventsPDF::initGenerator(Int_t code) {
//         if (!_DataCache) return;
//         if (_Parent->GetMaxValue() == 0 || _Parent->CheckChange()) {	
//             Double_t value = 0;		
//             if (code == 1) {	
//                 _MaxValue = 0;
//                 for (Int_t i = 0; i < _DataCache->_NTreeEntries; i++) {
//                     _TreeEntry = i;
//                     value = evaluateMC(&_DataCache->_vecRealGen, &_DataCache->_vecCatGen);
//                     if (value < 0) { std::cout << " BruEventsPDF::initGenerator -ve intensity !!! " << value << std::endl; exit(0); }
//                     if (value > _MaxValue) _MaxValue = value * 1.01;
//                 }
//                 _Parent->SetMaxValue(_MaxValue);
//             }
//         }		
//         _GeneratedIndices.clear();
//         _GeneratedWeights.clear();
//     }

//     void BruEventsPDF::generateEvent(Int_t code) {
//         Double_t value = 0;
//         if (!_DataCache) return;
        
//         while (_Geni < _DataCache->_NTreeEntries) {
//             _TreeEntry = _UseWeightsGen ? _Geni++ : IncrementGeni();
//             if (_UseWeightsGen) _Parent->SetGeni(_Geni);
            
//             if (!CheckRange("")) continue;
//             value = evaluateMC(&_DataCache->_vecRealGen, &_DataCache->_vecCatGen); 
            
//             if (_UseWeightsGen || (value > _MaxValue * RooRandom::uniform())) {
//                 for (Int_t i = 0; i < _Nvars; i++) *(_ProxSet[i]) = _DataCache->_vecReal[_TreeEntry * _Nvars + i]; 
//                 for (Int_t i = 0; i < _Ncats; i++) *(_CatSet[i]) = _DataCache->_vecCat[_TreeEntry * _Ncats + i];
                
//                 if (_UseWeightsGen) {
//                     _GeneratedWeights.push_back(value);
//                     _GeneratedIndices.push_back(_Geni - 1);
//                 } else {
//                     _GeneratedIndices.push_back(_DataCache->_TreeEntryNumber[_TreeEntry]);
//                 }
//                 return;
//             }
//         }
//         Fatal("BruEventsPDF::generateEvent", "Ran out of events at %lld", _Geni);
//     }

//     Int_t BruEventsPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
//         if (_ForceNumInt) return 0; 

//         Int_t unintegrated_count = 0;
//         Int_t plot_var_index = -1;

//         for (UInt_t i = 0; i < _ProxSet.size(); i++) {
//             if (!allVars.find(_ProxSet[i]->GetName())) {
//                 unintegrated_count++;
//                 plot_var_index = i;
//             }
//         }

//         if (unintegrated_count == 0) {
//             for (auto p : _ProxSet) {
//                 RooAbsArg* arg = allVars.find(p->GetName());
//                 if (arg) analVars.add(*arg);
//             }
//             for (auto c : _CatSet) {
//                 RooAbsArg* arg = allVars.find(c->GetName());
//                 if (arg) analVars.add(*arg);
//             }
//             if (!_DataCache && !_ForceConstInt) {
//                 Fatal("BruEventsPDF::getAnalyticalIntegral", "RooFit attempted SILENT numeric integration! No MC cache is loaded. Call SetNumInt(kTRUE).");
//             }
//             if (BruEventsPDF_IsPlotting && _ProxSet.size() == 1 && _CatSet.size() == 0 && _HistIntegrals.empty()) {
//                 HistIntegrals(rangeName);
//             }
//             return 1;
//         }

//         if (unintegrated_count == 1 && _DataCache) {
//             for (UInt_t i = 0; i < _ProxSet.size(); i++) {
//                 if (i != (UInt_t)plot_var_index) {
//                     RooAbsArg* arg = allVars.find(_ProxSet[i]->GetName());
//                     if (arg) analVars.add(*arg);
//                 }
//             }
//             for (auto c : _CatSet) {
//                 RooAbsArg* arg = allVars.find(c->GetName());
//                 if (arg) analVars.add(*arg);
//             }
//             return plot_var_index + 2; 
//         }

//         if (allVars.getSize() > 0) {
//             std::cout << "\n\n[BruEventsPDF ERROR] Failed to match analytical integral for PDF: " << GetName() << std::endl;
//             std::cout << "RooFit asked to integrate over: "; allVars.Print("");
//             Fatal("BruEventsPDF::getAnalyticalIntegral", "RooFit requested an integration over a subset of variables this PDF cannot handle analytically. Call SetNumInt(kTRUE).");
//         }
//         return 0; 
//     }

//     Double_t BruEventsPDF::analyticalIntegralForSampling(const char* rangeName) const {
//         Double_t integral = 0;
//         Long64_t accepted = 0;
//         Long64_t all = 0;
//         Long64_t ilow = 0;
//         Long64_t ihigh = 0;

//         SetLowHighVals(ilow, ihigh); 
        
//         if (CheckChange()) {
//             std::vector<double> values(ihigh - ilow);
//             for (Long64_t ie = ilow; ie < ihigh; ie++) {
//                 _TreeEntry = ie;
//                 if (!CheckRange(rangeName)) {
//                     values[all] = 0;
//                     ++all;
//                     continue;
//                 }
//                 values[all] = evaluateMC(&_DataCache->_vecReal, &_DataCache->_vecCat) * _DataCache->GetWeight(ie);
//                 integral += values[all];
//                 ++accepted; 
//                 ++all; 
//             }
//             if (accepted > 0) integral /= accepted;
            
//             double sum_of_diffs = 0.;
//             std::for_each(values.begin(), values.end(), [&sum_of_diffs, &integral] (double n) {
//                 double term = (n - integral);
//                 sum_of_diffs += term * term;
//             });
            
//             if (accepted > 1) {
//                 _SigmaIntegral = TMath::Sqrt(sum_of_diffs / (accepted - 1)) / TMath::Sqrt(accepted);
//             } else {
//                 _SigmaIntegral = 0.0;
//             }
//             _Last[0] = integral;
//         }      
//         return sampleIntegral(_Last[0], _SigmaIntegral);
//     }

//     Double_t BruEventsPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
//         if (code == 1) {
//             // BUG FIX: Only force constant integral during the full normalization step
//             if (_ForceConstInt) {
//                 _Last[0] = 1.0; 
//                 return _Last[0];
//             }

//             if (!CheckChange()) return _Last[0];
//             if (!_DataCache) Fatal("BruEventsPDF::analyticalIntegral", "Attempted MC Integration but _DataCache is null!");
            
//             auto check = AssertPositivePDF();
//             if (check == kFALSE) return _Last[0] = 0;

//             Long64_t accepted = 0;
//             Long64_t ilow = 0;
//             Long64_t ihigh = 0;
//             SetLowHighVals(ilow, ihigh); 
//             Double_t integral = 0.;

//             bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);

//             // ORIGINAL CLEAN LOOP RESTORED
//             for (Long64_t ie = ilow; ie < ihigh; ie++) {
//                 _TreeEntry = ie;
//                 if (hasRange && !CheckRange(rangeName)) continue; 
                
//                 accepted++;
//                 integral += evaluateMC(&_DataCache->_vecReal, &_DataCache->_vecCat) * _DataCache->GetWeight(ie);
//             }

//             if (accepted > 0) integral /= accepted;
//             else integral = 0;
            
//             _Last[0] = integral;
//             return _Last[0];
//         } 
//         else {
//             if (_HistIntegrals.empty()) HistIntegrals(rangeName);
            
//             Int_t vindex = code - 2;
//             if (vindex < 0 || vindex >= (Int_t)_ProxSet.size()) Fatal("BruEventsPDF::analyticalIntegral", "Invalid integration code");
            
//             Double_t vval = *(_ProxSet[vindex]);
//             Double_t integral = _HistIntegrals[vindex].Interpolate(vval);
//             if (integral < 0) integral = 0;
//             return integral;
//         }
//     }

//     Double_t BruEventsPDF::unnormalisedIntegral(Int_t code, const char* rangeName) const {
//         Double_t integral = 0;
//         Double_t nev = 0;
//         Double_t nMC = 0;
//         if (!_DataCache) return 0;

//         if (code == 1) {
//             bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);
            
//             // ORIGINAL CLEAN LOOP RESTORED
//             for (Long64_t ie = 0; ie < _DataCache->_NTreeEntries; ie++) {
//                 _TreeEntry = ie;
//                 if (hasRange && !CheckRange(rangeName)) continue; 
                
//                 integral += evaluateMC(&_DataCache->_vecReal, &_DataCache->_vecCat) * _DataCache->GetWeight(ie);
//                 nev++;
//             }
//         } else if (code == 2 && _DataCache->_HasMCGenTree) {
//             for (Long64_t ie = 0; ie < _DataCache->_NMCGenTreeEntries; ie++) {
//                 _TreeEntry = ie;
//                 integral += evaluateMC(&_DataCache->_vecRealMCGen, &_DataCache->_vecCatMCGen);
//                 nMC++;
//             }
//         }
//         return integral;
//     }
    
//     void BruEventsPDF::HistIntegrals(const char* rangeName) const {
//         if (!_DataCache) return;
//         Long64_t ilow = 0;
//         Long64_t ihigh = 0;
//         SetLowHighVals(ilow, ihigh);
      
//         for (Int_t i = 0; i < _Nvars; i++) {
//             auto arg = dynamic_cast<const RooRealVar*>(&_ProxSet[i]->arg());
//             if (arg) _HistIntegrals.emplace_back(arg->GetName(), arg->GetName(), arg->getBins(), arg->getMin(), arg->getMax());
//         }
      
//         bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);
//         Long64_t accepted = 0;
        
//         // ORIGINAL CLEAN LOOP RESTORED
//         for (Int_t ie = ilow; ie < ihigh; ie++) {
//             _TreeEntry = ie;
//             if (hasRange && !CheckRange(rangeName)) continue; 
            
//             accepted++;
//             Double_t w = _DataCache->GetWeight(ie);
//             Double_t value = evaluateMC(&_DataCache->_vecReal, &_DataCache->_vecCat) * w;
            
//             for (Int_t vindex = 0; vindex < _Nvars; vindex++) {
//                 _HistIntegrals[vindex].Fill(_DataCache->_vecReal[_TreeEntry * _Nvars + vindex], value / _HistIntegrals[vindex].GetBinWidth(1));
//             }
//         }
//         for (Int_t vindex = 0; vindex < _Nvars; vindex++) {
//             if (accepted > 0) _HistIntegrals[vindex].Scale(1. / accepted);
//         }
//         _Parent->SetHistIntegrals(_HistIntegrals);
//     }

//     void BruEventsPDF::SetLowHighVals(Long64_t& ilow, Long64_t& ihigh) const {
//         ilow = 0; ihigh = 0;
//         if (_Parent) {
//             ilow = _Parent->GetIntRangeLow();
//             ihigh = _Parent->GetIntRangeHigh();
//         } else {
//             ilow = GetIntRangeLow();
//             ihigh = GetIntRangeHigh();
//         }
//         if (ihigh == 0 && _NInt > -1) ihigh = _NInt;
//         else if (ihigh == 0 && _DataCache) ihigh = _DataCache->_NTreeEntries; 
//         if (_DataCache && ihigh > (Long64_t)_DataCache->_NTreeEntries) ihigh = _DataCache->_NTreeEntries;
//     }

   
//   Bool_t BruEventsPDF::CheckRange(const char* rangeName) const {
//     if (!_DataCache) return kFALSE;
    
//     const auto& vecReal = _DataCache->_vecReal; 
    
//     for (UInt_t i = 0; i < _ProxSet.size(); i++) {
//       auto var = dynamic_cast<const RooRealVar*>(&(_ProxSet[i]->arg()));
//              if (!var->inRange(vecReal[_TreeEntry * _Nvars + i], rangeName)) return kFALSE;
//     }
//     return kTRUE;
//   }
  
//     Bool_t BruEventsPDF::CheckChange() const {
//         Bool_t hasChanged = false;
//         for (Int_t i = 1; i < _Npars + 1; i++) {
//             if (_Last[i] != *(_ParSet[i - 1])) hasChanged = true;
//         }
//         if (hasChanged) {
//             for (Int_t i = 1; i < _Npars + 1; i++) _Last[i] = *(_ParSet[i - 1]);
//         }
//         return hasChanged;
//     }

//     void BruEventsPDF::CheckIntegralParDep(Int_t Ntests) {
//         _CheckInt = Ntests;
//         if (!_DataCache) return; 
    
//         Long64_t saveNint = _NInt;
//         _NInt = _DataCache->_NTreeEntries;
//         Ntests = (Ntests * _ParSet.size());
    
//         RooRealVar integral("integral", "integral", 0, 0, 2);
//         if (_NInt > 0) integral.setError(sqrt(_NInt) / _NInt); 
//         RooDataSet ds("intds", "intds", RooArgSet(integral));
//         std::vector<Double_t> SavedPars;
        
//         for (auto &ip : _ParSet) {
//             auto par = (dynamic_cast<const RooRealVar*>(&(ip->arg())));
//             SavedPars.push_back(par->getValV());
//         }
        
//         for (Int_t ir = 0; ir < Ntests; ir++) { 
//             for (auto &ip : _ParSet) {
//                 auto par = (RooRealVar*)(&(ip->arg()));
//                 par->setVal((par->getMax("") - par->getMin("")) * RooRandom::uniform() + par->getMin(""));
//             }
//             integral.setVal(analyticalIntegral(1, ""));
//             ds.add(RooArgSet(integral));
//         }
    
//         Double_t low = 0;
//         Double_t high = 0;
//         ds.getRange(integral, low, high);
//         integral.setRange(low, high);
//         RooPlot *frame = integral.frame();
//         ds.plotOn(frame);
    
//         frame->Draw();
    
//         RooRealVar mean("mean", "mean", ds.mean(integral));
//         RooRealVar pvar("IntPull", "Integral Pull Dist.", -5, 5);
//         RooPullVar pull("IntPull", "Integral Pull Dist.", integral, mean);
//         ds.addColumn(pull, kFALSE);
    
//         ds.getRange(pvar, low, high);
//         pvar.setRange(low, high);
//         RooPlot *framePull = pvar.frame();
//         ds.plotOn(framePull);
//         RooRealVar mp("mp", "mp", 0, -5, 5);
//         RooRealVar sp("sp", "sp", 1, 0, 100);
//         RooGaussian gp("gp", "gp", pvar, mp, sp);
//         gp.fitTo(ds);
//         gp.paramOn(framePull);
//         gp.plotOn(framePull);
    
//         new TCanvas();
//         framePull->Draw();
    
//         _ConstInt = mean.getVal();
//         _NInt = saveNint;
//         if (sp.getVal() < 2) SetConstInt();
//         for (UInt_t ip = 0; ip < _ParSet.size(); ip++) {
//             auto par = (RooRealVar*)(&(_ParSet[ip]->arg()));
//             par->setVal(SavedPars[ip]);
//         }
//         _CheckInt = kFALSE; 
//     }

//     Bool_t BruEventsPDF::AddProtoData(const RooDataSet* data) {
//         if (!_DataCache || !_DataCache->_NTreeEntries) return kFALSE;
    
//         auto mutableCache = std::make_shared<MCEventCache>(*_DataCache);
//         const RooArgSet *dataVars = data->get();
//         Long64_t Nentries = data->numEntries();
//         std::vector<Long64_t> vrandom(Nentries);
//         for (Long64_t ir = 0; ir < Nentries; ir++) vrandom[ir] = ir;
//         std::shuffle(vrandom.begin(), vrandom.end(), std::mt19937(std::random_device()()));
    
//         std::vector<Short_t> protoDataForVar, protoDataForCat;
//         for (auto* arg : *dataVars) {
//             if (TString("UID") == arg->GetName()) continue; 
            
//             for (Int_t ip = 0; ip < _Nvars; ip++) {
//                 if (TString(arg->GetName()) == TString(_ProxSet[ip]->GetName())) {
//                     if (std::find(_ProtoRealVars.begin(), _ProtoRealVars.end(), TString(arg->GetName())) != _ProtoRealVars.end()) {
//                         protoDataForVar.push_back(ip);
//                     }
//                 }
//             }
//             for (Int_t ip = 0; ip < _Ncats; ip++) {
//                 if (TString(arg->GetName()) == TString(_CatSet[ip]->GetName())) {
//                     if (std::find(_ProtoCatVars.begin(), _ProtoCatVars.end(), TString(arg->GetName())) != _ProtoCatVars.end()) {
//                         protoDataForCat.push_back(ip);
//                     }
//                 }
//             }
//         }
    
//         Long64_t idata = 0;
//         if (!(protoDataForVar.size() + protoDataForCat.size())) {
//             return kTRUE; 
//         }
    
//         for (Long64_t id = 0; id < mutableCache->_NTreeEntries; id++) {
//             dataVars = data->get(vrandom[idata]);
//             for (short ip : protoDataForVar) {
//                 Double_t val = dataVars->getRealValue(_ProxSet[ip]->GetName());
//                 mutableCache->_vecReal[id * _Nvars + ip] = val;
//                 mutableCache->_vecRealGen[id * _Nvars + ip] = val;
//             }  
//             for (short ip : protoDataForCat) {
//                 Int_t val = dataVars->getCatIndex(_CatSet[ip]->GetName());
//                 mutableCache->_vecCat[id * _Ncats + ip] = val;
//                 mutableCache->_vecCatGen[id * _Ncats + ip] = val;     
//             }
            
//             // BUG FIX: Accurate increment before shuffle check
//             idata++;
//             if (idata >= (Long64_t)vrandom.size()) {
//                 std::shuffle(vrandom.begin(), vrandom.end(), std::mt19937(std::random_device()()));
//                 idata = 0;
//             }
//         }
        
//         _DataCache = mutableCache;
//         return kTRUE;  
//     }

//     void BruEventsPDF::SetNextRange(Int_t ir) {
//         Long64_t Nentries = _DataCache ? _DataCache->_NTreeEntries : 0;
//         Int_t range = ((Double_t)Nentries) / _NRanges;
//         _IntRangeLow = ir * range;
//         _IntRangeHigh = (ir + 1) * range;
//     }

//     void BruEventsPDF::MakeAssertPostiveData() {
//         auto saveTreeEntry = _TreeEntry;
//         _TreeEntry = 0;

//         auto NVars = _ProxSet.size();
//         _AssertPosDataReal.resize(_Napd * NVars);

//         auto NCats = _CatSet.size();
//         _AssertPosDataCats.resize(_Napd * NCats);

//         for (Long64_t iapd = 0; iapd < _Napd; ++iapd) {
//             UInt_t ivar = 0;
//             for (auto v : _ProxSet) {
//                 auto vr = dynamic_cast<const RooRealVar*>(&v->arg());
//                 if (vr != nullptr) {
//                     _AssertPosDataReal[_TreeEntry * NVars + ivar] = gRandom->Uniform(vr->getMin(""), vr->getMax(""));
//                     ++ivar;
//                 }
//             }
//             UInt_t icat = 0;
//             for (auto v : _CatSet) {
//                 auto vc = dynamic_cast<const RooCategory*>(&v->arg());
//                 if (vc != nullptr) {
//                     auto catstate = gRandom->Integer(vc->size());
//                     auto val = vc->getOrdinal(catstate).second;
//                     _AssertPosDataCats[_TreeEntry * NCats + icat] = val;
//                     ++icat;
//                 }
//             }
//             _TreeEntry++;
//         }
//         _TreeEntry = saveTreeEntry;
//     }
   
//     Bool_t BruEventsPDF::AssertPositivePDF() const {
//         InitAssertPositiveCheck();
//         auto saveTreeEntry = _TreeEntry;
//         _TreeEntry = 0;
//         for (Long64_t iapd = 0; iapd < _Napd; ++iapd) {
//             auto val = evaluateMC(&_AssertPosDataReal, &_AssertPosDataCats);
//             ++_TreeEntry;
//             if (val < -1E-4) { 
//                 logEvalError("BruEventsPDF::AssertPositivePDF() PDF cannot be -ve...");
//                 _TreeEntry = saveTreeEntry;
//                 FinishAssertPositiveCheck(); 
//                 return kFALSE;
//             }
//         }
//         _TreeEntry = saveTreeEntry;
//         FinishAssertPositiveCheck();
//         return kTRUE;
//     }

// } // namespace bru
