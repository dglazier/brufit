/**
 * @file BruComponentsPDF.cpp
 */

#include "BruComponentsPDF.h" 
#include <RooAbsReal.h> 
#include <RooAbsArg.h>
#include <RooAbsCategory.h> 
#include <cmath> 
#include "TMath.h" 

namespace bru {
    
    BruComponentsPDF::BruComponentsPDF(const char *name, const char *title, Double_t base, const RooArgList& obsList, const std::vector<RooArgList> compList)
      : BruEventsPDF(name, title), _BaseLine(base), _ActualObs{"ActualObs", "Actual observables", this}, _ActualCats{"ActualCats", "Actual categories", this}, _ActualComps{"ActualComps", "Actual components", this}
    {
        _NObs = 0;
        _NCats = 0;
        for (Int_t i = 0; i < obsList.getSize(); i++) {
            if (dynamic_cast<RooRealVar*>(&obsList[i])) {
                std::unique_ptr<RooRealProxy> tempR{new RooRealProxy(obsList[i].GetName(), obsList[i].GetName(), this, dynamic_cast<RooAbsReal&>(obsList[i]))};
                _NObs++;
                _ActualObs.add(dynamic_cast<RooAbsReal&>(obsList[i]));
                _Observables.push_back(std::move(tempR));
                continue;
            }
            if (dynamic_cast<RooCategory*>(&obsList[i])) {
                std::unique_ptr<RooCategoryProxy> tempC{new RooCategoryProxy(obsList[i].GetName(), obsList[i].GetName(), this, dynamic_cast<RooAbsCategory&>(obsList[i]))};
                _NCats++;
                _ActualCats.add(dynamic_cast<RooAbsCategory&>(obsList[i]));
                _Categories.push_back(std::move(tempC));
                continue;
            }
        }

        _NComps = compList.size();
        for (auto& comp: compList) {
            vecUPtrReal vterms;
            for (Int_t i = 0; i < comp.getSize(); i++) {
                _ActualComps.add(dynamic_cast<RooAbsReal&>(comp[i]));
                std::unique_ptr<RooRealProxy> temp{new RooRealProxy(comp[i].GetName(), comp[i].GetName(), this, dynamic_cast<RooAbsReal&>(comp[i]))};
                vterms.push_back(std::move(temp));
            }
            _Components.push_back(std::move(vterms));
        }
      
        MakeSets();
     }

    BruComponentsPDF::BruComponentsPDF(const BruComponentsPDF& other, const char* name) :
      BruEventsPDF(other, name),
      _BaseLine(other._BaseLine),
      _ActualComps("AllComponents", this, other._ActualComps),
      _ActualCats("AllCategories", this, other._ActualCats),
      _ActualObs("AllObservables", this, other._ActualObs)
    {
        _WeightedBaseLine = other._WeightedBaseLine;

        for (const auto& fObservable : other._Observables) {
            unique_ptr<RooRealProxy> temp{new RooRealProxy(fObservable->GetName(), this, *fObservable)};
            _Observables.push_back(std::move(temp));
        }
        for (const auto& fCategorie : other._Categories) {
            unique_ptr<RooCategoryProxy> temp{new RooCategoryProxy(fCategorie->GetName(), this, *fCategorie)};
            _Categories.push_back(std::move(temp));
        }
      
        for (auto& comp: other._Components) {
            vecUPtrReal vterms;
            for (const auto& i : comp) {
                unique_ptr<RooRealProxy> temp{new RooRealProxy(i->GetName(), this, *i)};
                vterms.push_back(std::move(temp));
            }
            _Components.push_back(std::move(vterms));
        }

        _NObs = other._NObs;
        _NCats = other._NCats;
        _NComps = other._NComps;
      
        MakeSets(); 
        
        _Last = other._Last;
        _LastLength = other._LastLength;
        _CacheCompDepIntegral = other._CacheCompDepIntegral;
        _CacheCompDepSigmaIntegral = other._CacheCompDepSigmaIntegral;
        _MCAPDepTerm = other._MCAPDepTerm;
    } 

    void BruComponentsPDF::MakeSets() {
        for (auto& obs: _Observables) {
            _ProxSet.push_back(obs.get());
            _IntegrateObs.push_back(dynamic_cast<RooRealVar*>(_IntegrateSet.addClone(obs->arg())));
        }
        for (auto& cat: _Categories) {
            _CatSet.push_back(cat.get());
            _IntegrateCats.push_back(dynamic_cast<RooCategory*>(_IntegrateSet.addClone(cat->arg())));
        }
 
        for (auto& comp: _Components) {
            for (auto& term: comp) {
                auto argTerm = _ActualComps.find(term->GetName());
                auto vars = argTerm->getVariables();
          
                for (auto* arg: *vars) {
                    if (!_ActualObs.contains(*arg) && !_ActualCats.contains(*arg) && !_Parameters.contains(*arg)) {
                        _Parameters.add(*arg);
                        _myVarProxies.push_back(std::unique_ptr<RooRealProxy>{new RooRealProxy(arg->GetName(), arg->GetName(), this, *dynamic_cast<RooAbsReal*>(arg))});
                        _ParSet.push_back(_myVarProxies.back().get());
                    }
                }
          
                if (vars->getSize() == 0) {
                    if (!_ActualObs.contains(*argTerm) && !_ActualCats.contains(*argTerm) && !_Parameters.contains(*argTerm)) {
                        _Parameters.add(*argTerm);
                        _myVarProxies.push_back(std::unique_ptr<RooRealProxy>{new RooRealProxy(argTerm->GetName(), argTerm->GetName(), this, *dynamic_cast<RooAbsReal*>(argTerm))});
                        _ParSet.push_back(_myVarProxies.back().get());
                    }
                }
            }
        }
        InitSets();
    }

    Int_t BruComponentsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const {	
        if (!_DataCache) return 0; 
        if (matchArgs(directVars, generateVars, VarSet(0))) return 1;
        return 0;
    }

    Double_t BruComponentsPDF::evaluateData() const {
        Double_t val = _BaseLine;
        for (auto& comp: _Components) {
            Double_t product = 1;
            for (auto& term: comp) {
                product *= *term.get(); 
            }
            val += product; 
        }
        return val;
    }
    
    void BruComponentsPDF::RedirectServersToPdf() const {
        for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_IntegrateSet);
            }	
        }	
    }

    void BruComponentsPDF::RedirectServersToData() const {
        for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_ActualObs);
            }	
        }	
    }

    void BruComponentsPDF::HistIntegrals(const char* rangeName) const {
      initIntegrator();
      
        for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_IntegrateSet);
            }
        }
     
        BruEventsPDF::HistIntegrals(rangeName);
     
        for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_ActualObs);
            }
        }
    }

    Double_t BruComponentsPDF::cacheMCAP(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const {
        _MCAPDepTerm.resize(_Napd);
	for(auto& row : _MCAPDepTerm) row.resize(_NComps);
	
        for (Long64_t iev = 0; iev < _Napd; ++iev) {
            for (Int_t ii = 0; ii < _Nvars; ii++) {
                _IntegrateObs[ii]->setVal(vars->at(iev * _Nvars + ii));
            }
            for (Int_t ii = 0; ii < _Ncats; ii++) {
                _IntegrateCats[ii]->setIndex(cats->at(iev * _Ncats + ii));
            }
      
            Double_t val = _BaseLine;
            Int_t icomp = 0;
            for (auto& comp : _Components) {
                Double_t depprod = 1; 
                for (auto& term : _DependentTermProxy[icomp]) {
                    depprod *= *term;
                }
  		_MCAPDepTerm[iev][icomp] = depprod;
                Double_t product = 1;
                for (auto& term : comp) {
                    product *= *term.get(); 
                }
                val += product; 
                ++icomp;
            }
        }
        return 0.;
    }

    Double_t BruComponentsPDF::evaluateMCAP() const {
        Double_t val = _BaseLine;
        Int_t icomp = 0;
     
        for (auto& comp : _Components) {
            Double_t product = _MCAPDepTerm[_TreeEntry][icomp];
            for (auto& term : _IndependentTermProxy[icomp]) {
                product *= *term;
            }
            val += product; 
            ++icomp;
        }
        return val;
    }
    
    Double_t BruComponentsPDF::evaluateMC(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const {
        if (_assertPostive) return evaluateMCAP();
      
        for (Int_t ii = 0; ii < _Nvars; ii++) {
            _IntegrateObs[ii]->setVal(vars->at(_TreeEntry * _Nvars + ii));
        }
        for (Int_t ii = 0; ii < _Ncats; ii++) {
            _IntegrateCats[ii]->setIndex(cats->at(_TreeEntry * _Ncats + ii));
        }
        return evaluateData();
    }

    bool BruComponentsPDF::isDirectGenSafe(const RooAbsArg& arg) const {
        if (_ActualObs.find(arg.GetName())) return kTRUE;
        if (_ActualCats.find(arg.GetName())) return kTRUE;
        return kFALSE;
    }

    void BruComponentsPDF::initGenerator(Int_t code) {
      initIntegrator();
        RedirectServersToPdf();
        BruEventsPDF::initGenerator(code);
    }

    void BruComponentsPDF::initIntegrator() const {
      if (!_once) return; // GUARD: Ensure this only ever runs once
      _once = kFALSE;     // Flip the flag
       
      BruEventsPDF::initIntegrator();
     
        _DependentTermProxy.resize(_NComps);
        _DependentTermParams.resize(_NComps);
        _PrevParVals.resize(_NComps);
        _CacheCompDepIntegral.resize(_NComps);
        _CacheCompDepSigmaIntegral.resize(_NComps);
        _IndependentTermProxy.resize(_NComps);
    
        UInt_t icomp = 0;
        for (auto& comp : _Components) {
            _CacheCompDepIntegral[icomp] = 1;
            _CacheCompDepSigmaIntegral[icomp] = 0;
            UInt_t iterm = 0;
            Double_t product = 1;
            for (auto& term : comp) {
                auto arg = _ActualComps.find(term->GetName());
                auto deps = arg->getObservables(VarSet(0));
        
                if (deps->getSize()) {
                    _DependentTermProxy[icomp].push_back(term.get());
                    auto parDeps = arg->getObservables(_Parameters);
            
                    if (parDeps->getSize()) {
                        for (auto* p_arg : *parDeps) {
                            auto* rarg = dynamic_cast<RooRealVar*>(p_arg);      
                            if (!vecContains(rarg, _DependentTermParams[icomp])) {
                                _DependentTermParams[icomp].push_back(rarg);
                                Double_t initf = -1E6;
                                _PrevParVals[icomp].push_back(initf);
                            }
                        }
                    }
                } else {
                    _IndependentTermProxy[icomp].push_back(term.get());
                }
                iterm++;
            }
            icomp++;
        }
    }

    void BruComponentsPDF::DoFirstIntegrations(const char* rangeName) const {
        for (UInt_t icomp = 0; icomp < _NComps; icomp++)
            _RecalcComponent.push_back(icomp);
        
        RecalcComponentIntegrals(0, rangeName);

        if (_UseSamplingIntegral == kTRUE)
            RecalcComponentIntegralsSampling(0, rangeName);
        
        _FirstCalculation = kFALSE;
    }
    
    Double_t BruComponentsPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
        if (code != 1) return BruEventsPDF::analyticalIntegral(code, rangeName);
        if (code == 1 && _ForceConstInt && !_DataCache) { _Last[0] = 1; return _Last[0]; }
  
        if (!CheckChange()) return _Last[0];
	initIntegrator();
	
        auto check = AssertPositivePDF();
        if (check == kFALSE) return _Last[0] = 0;
       
        if (_FirstCalculation == kTRUE) DoFirstIntegrations();
       
        if (_WeightedBaseLine == 0 && _BaseLine != 0 && _DataCache->_UseEvWeights)
            CalcWeightedBaseLine(rangeName);
        else
            _WeightedBaseLine = _BaseLine;

        Bool_t needRecalc = kFALSE;
        _RecalcComponent.clear();

        for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
            if (_DependentTermProxy[icomp].size()) {
                UInt_t ipar = 0;
                for (auto par : _DependentTermParams[icomp]) {
                    Double_t previous = _PrevParVals[icomp][ipar];
                    Double_t pval = par->getVal();
             
                    if (pval != previous) { 
                        needRecalc = kTRUE;
                        if (!vecContains(icomp, _RecalcComponent)) _RecalcComponent.push_back(icomp);
                    }
                    ipar++;
                }
            }
        }
     
        if (needRecalc) RecalcComponentIntegrals(code, rangeName);

        Double_t integral = _WeightedBaseLine;
        for (UInt_t icomp = 0; icomp < _NComps; icomp++)
            integral += componentIntegral(icomp);

        _Last[0] = integral;
        return integral;
    }
     
    void BruComponentsPDF::CalcWeightedBaseLine(const char* rangeName) const {
        Long64_t ilow, ihigh = 0;
        SetLowHighVals(ilow, ihigh);
      
        for (const auto& icomp : _RecalcComponent) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_IntegrateSet);
            }
        }
     
        _WeightedBaseLine = 0;
        Long64_t accepted = 0;   
      
        for (Long64_t ie = ilow; ie < ihigh; ie++) {
            _TreeEntry = ie;
            if (!CheckRange(TString(rangeName).Data())) continue;
            accepted++;
            _WeightedBaseLine += GetIntegralWeight(ie);
        }
        if (accepted > 0) _WeightedBaseLine /= accepted;
    }

    void BruComponentsPDF::RecalcComponentIntegrals(Int_t code, const char* rangeName) const {
        Long64_t ilow, ihigh = 0;
        SetLowHighVals(ilow, ihigh);
      
        for (const auto& icomp : _RecalcComponent) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_IntegrateSet);
            }
        }

        Long64_t accepted = 0;
      
        for (Long64_t ie = ilow; ie < ihigh; ie++) {
            _TreeEntry = ie;
            if (!CheckRange(TString(rangeName).Data())) continue;
            accepted++;
        
            for (Int_t ii = 0; ii < _Nvars; ii++)
                _IntegrateObs[ii]->setVal(_DataCache->_vecReal[_TreeEntry * _Nvars + ii]);
            for (Int_t ii = 0; ii < _Ncats; ii++)
                _IntegrateCats[ii]->setIndex(_DataCache->_vecCat[_TreeEntry * _Ncats + ii]);
          
            for (const auto& icomp : _RecalcComponent) {
                Double_t product = 1.;
                for (const auto& term : _DependentTermProxy[icomp]) {
                    product *= *term;
                }
                product *= GetIntegralWeight(ie);
                _CacheCompDepIntegral[icomp] += product;
            }
        }
       
        for (const auto& icomp : _RecalcComponent) {
            _CacheCompDepIntegral[icomp] = _CacheCompDepIntegral[icomp] / accepted;
        }
        for (const auto& icomp : _RecalcComponent) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_ActualObs);
            }
        }
    }
    
    void BruComponentsPDF::RecalcComponentIntegralsSampling(Int_t code, const char* rangeName) const {
        if (_RecalcComponent.empty() == kTRUE) return;
      
        Long64_t ilow, ihigh = 0;
        SetLowHighVals(ilow, ihigh);
      
        for (const auto& icomp : _RecalcComponent) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_IntegrateSet);
            }
        }

        Long64_t accepted = 0;
        Long64_t all = 0;
        std::vector<Double_t> sumSquares(_RecalcComponent.size());
      
        for (Long64_t ie = ilow; ie < ihigh; ie++) {
            _TreeEntry = ie;
            if (!CheckRange(rangeName)) { ++all; continue; }
        
            for (Int_t ii = 0; ii < _Nvars; ii++)
                _IntegrateObs[ii]->setVal(_DataCache->_vecReal[_TreeEntry * _Nvars + ii]);
            for (Int_t ii = 0; ii < _Ncats; ii++)
                _IntegrateCats[ii]->setIndex(_DataCache->_vecCat[_TreeEntry * _Ncats + ii]);
          
            for (const auto& icomp : _RecalcComponent) {
                Double_t product = 1;
                for (const auto& term : _DependentTermProxy[icomp]) {
                    product *= *term;
                }
                product *= GetIntegralWeight(ie);
          
                _CacheCompDepIntegral[icomp] += product;
                sumSquares[icomp] += product * product;
            }
            ++accepted;
            ++all;
        }
      
        for (const auto& icomp : _RecalcComponent) {
            _CacheCompDepIntegral[icomp] = _CacheCompDepIntegral[icomp] / (accepted - 1);
            _CacheCompDepSigmaIntegral[icomp] = sumSquares[icomp] / accepted;
        }
        
        _NUsedForIntegral = accepted;
      
        for (const auto& icomp : _RecalcComponent) {
            for (const auto& term : _DependentTermProxy[icomp]) {
                auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
                unconstTerm->recursiveRedirectServers(_ActualObs);
            }
        }
    }

    Double_t BruComponentsPDF::componentIntegral(Int_t icomp) const {
        Double_t product = 1;
        product *= _CacheCompDepIntegral[icomp];
        for (auto& term : _IndependentTermProxy[icomp]) {
            product *= *term;
        }
        return product; 
    }

    Double_t BruComponentsPDF::componentVariance(Int_t icomp) const {
        Double_t product = 1;
        product *= _CacheCompDepSigmaIntegral[icomp];
        for (auto& term : _IndependentTermProxy[icomp]) {
            product *= (*term) * (*term);
        }
        return product; 
    }
    
    Bool_t BruComponentsPDF::CheckChange() const {
        Bool_t hasChanged = false;
        for (Int_t i = 1; i < _Npars + 1; i++) {
            if (_Last[i] != *(_ParSet[i - 1])) {
                hasChanged = true;
            }
        }
        if (hasChanged) {
            for (Int_t i = 1; i < _Npars + 1; i++) {
                _Last[i] = *(_ParSet[i - 1]);
            }
        }
        return hasChanged;
    }   

} // namespace bru

// /**
//  * @file BruComponentsPDF.cpp
//  */

// #include "BruComponentsPDF.h" 
// #include <RooAbsReal.h> 
// #include <RooAbsArg.h>
// #include <RooAbsCategory.h> 
// #include <cmath> 
// #include "TMath.h" 

// namespace bru {
    
//     BruComponentsPDF::BruComponentsPDF(const char *name, const char *title, Double_t base, const RooArgList& obsList, const std::vector<RooArgList> compList)
//       : BruEventsPDF(name, title), _BaseLine(base), _ActualObs{"ActualObs", "Actual observables", this}, _ActualCats{"ActualCats", "Actual categories", this}, _ActualComps{"ActualComps", "Actual components", this}
//     {
//         _NObs = 0;
//         _NCats = 0;
//         for (Int_t i = 0; i < obsList.getSize(); i++) {
//             if (dynamic_cast<RooRealVar*>(&obsList[i])) {
//                 std::unique_ptr<RooRealProxy> tempR{new RooRealProxy(obsList[i].GetName(), obsList[i].GetName(), this, dynamic_cast<RooAbsReal&>(obsList[i]))};
//                 _NObs++;
//                 _ActualObs.add(dynamic_cast<RooAbsReal&>(obsList[i]));
//                 _Observables.push_back(std::move(tempR));
//                 continue;
//             }
//             if (dynamic_cast<RooCategory*>(&obsList[i])) {
//                 std::unique_ptr<RooCategoryProxy> tempC{new RooCategoryProxy(obsList[i].GetName(), obsList[i].GetName(), this, dynamic_cast<RooAbsCategory&>(obsList[i]))};
//                 _NCats++;
//                 _ActualCats.add(dynamic_cast<RooAbsCategory&>(obsList[i]));
//                 _Categories.push_back(std::move(tempC));
//                 continue;
//             }
//         }

//         _NComps = compList.size();
//         for (auto& comp: compList) {
//             vecUPtrReal vterms;
//             for (Int_t i = 0; i < comp.getSize(); i++) {
//                 _ActualComps.add(dynamic_cast<RooAbsReal&>(comp[i]));
//                 std::unique_ptr<RooRealProxy> temp{new RooRealProxy(comp[i].GetName(), comp[i].GetName(), this, dynamic_cast<RooAbsReal&>(comp[i]))};
//                 vterms.push_back(std::move(temp));
//             }
//             _Components.push_back(std::move(vterms));
//         }
      
//         MakeSets();
//         _IntegratorInit = kFALSE;
//     }

//     BruComponentsPDF::BruComponentsPDF(const BruComponentsPDF& other, const char* name) :
//       BruEventsPDF(other, name),
//       _BaseLine(other._BaseLine),
//       _ActualComps("AllComponents", this, other._ActualComps),
//       _ActualCats("AllCategories", this, other._ActualCats),
//       _ActualObs("AllObservables", this, other._ActualObs)
//     {
//         _WeightedBaseLine = other._WeightedBaseLine;

//         // BUG FIX: Creating proxies from the _ActualObs correctly reconnects clones for plotting
//         for (Int_t i = 0; i < _ActualObs.getSize(); i++) {
//             auto& clonedArg = static_cast<RooAbsReal&>(_ActualObs[i]);
//             unique_ptr<RooRealProxy> temp{new RooRealProxy(clonedArg.GetName(), clonedArg.GetName(), this, clonedArg)};
//             _Observables.push_back(std::move(temp));
//         }
//         for (Int_t i = 0; i < _ActualCats.getSize(); i++) {
//             auto& clonedCat = static_cast<RooAbsCategory&>(_ActualCats[i]);
//             unique_ptr<RooCategoryProxy> temp{new RooCategoryProxy(clonedCat.GetName(), clonedCat.GetName(), this, clonedCat)};
//             _Categories.push_back(std::move(temp));
//         }
      
//         Int_t flat_idx = 0;
//         for (auto& comp: other._Components) {
//             vecUPtrReal vterms;
//             for (const auto& i : comp) {
//                 auto& clonedComp = static_cast<RooAbsReal&>(_ActualComps[flat_idx]);
//                 unique_ptr<RooRealProxy> temp{new RooRealProxy(clonedComp.GetName(), clonedComp.GetName(), this, clonedComp)};
//                 vterms.push_back(std::move(temp));
//                 flat_idx++;
//             }
//             _Components.push_back(std::move(vterms));
//         }

//         _NObs = other._NObs;
//         _NCats = other._NCats;
//         _NComps = other._NComps;
      
//         MakeSets(); 
        
//         _Last = other._Last;
//         _LastLength = other._LastLength;
//         _CacheCompDepIntegral = other._CacheCompDepIntegral;
//         _CacheCompDepSigmaIntegral = other._CacheCompDepSigmaIntegral;
//         _MCAPDepTerm = other._MCAPDepTerm;
//         _IntegratorInit = kFALSE;
//     } 

//     void BruComponentsPDF::MakeSets() {
//         for (auto& obs: _Observables) {
//             _ProxSet.push_back(obs.get());
//             _IntegrateObs.push_back(dynamic_cast<RooRealVar*>(_IntegrateSet.addClone(obs->arg())));
//         }
//         for (auto& cat: _Categories) {
//             _CatSet.push_back(cat.get());
//             _IntegrateCats.push_back(dynamic_cast<RooCategory*>(_IntegrateSet.addClone(cat->arg())));
//         }
 
//         for (auto& comp: _Components) {
//             for (auto& term: comp) {
//                 auto argTerm = _ActualComps.find(term->GetName());
//                 auto vars = argTerm->getVariables();
          
//                 for (auto* arg: *vars) {
//                     if (!_ActualObs.contains(*arg) && !_ActualCats.contains(*arg) && !_Parameters.contains(*arg)) {
//                         _Parameters.add(*arg);
//                         _myVarProxies.push_back(std::unique_ptr<RooRealProxy>{new RooRealProxy(arg->GetName(), arg->GetName(), this, *dynamic_cast<RooAbsReal*>(arg))});
//                         _ParSet.push_back(_myVarProxies.back().get());
//                     }
//                 }
          
//                 if (vars->getSize() == 0) {
//                     if (!_ActualObs.contains(*argTerm) && !_ActualCats.contains(*argTerm) && !_Parameters.contains(*argTerm)) {
//                         _Parameters.add(*argTerm);
//                         _myVarProxies.push_back(std::unique_ptr<RooRealProxy>{new RooRealProxy(argTerm->GetName(), argTerm->GetName(), this, *dynamic_cast<RooAbsReal*>(argTerm))});
//                         _ParSet.push_back(_myVarProxies.back().get());
//                     }
//                 }
//             }
//         }
//         InitSets();
//     }

//     void BruComponentsPDF::SetupIntegrator() const {
//         if (_IntegratorInit) return;
        
//         auto self = const_cast<BruComponentsPDF*>(this);
     
//         self->_DependentTermProxy.resize(_NComps);
//         self->_DependentTermParams.resize(_NComps);
//         self->_PrevParVals.resize(_NComps);
//         self->_CacheCompDepIntegral.resize(_NComps);
//         self->_CacheCompDepSigmaIntegral.resize(_NComps);
//         self->_IndependentTermProxy.resize(_NComps);
    
//         UInt_t icomp = 0;
//         for (auto& comp : _Components) {
//             // BUG FIX: Prevent infinite sum accumulation
//             self->_CacheCompDepIntegral[icomp] = 0.0;
//             self->_CacheCompDepSigmaIntegral[icomp] = 0.0;
            
//             for (auto& term : comp) {
//                 auto arg = _ActualComps.find(term->GetName());
//                 auto deps = arg->getObservables(VarSet(0)); 
        
//                 if (deps->getSize() > 0) {
//                     self->_DependentTermProxy[icomp].push_back(term.get());
//                     auto parDeps = arg->getObservables(_Parameters);
            
//                     if (parDeps->getSize() > 0) {
//                         for (auto* p_arg : *parDeps) {
//                             auto* rarg = dynamic_cast<RooRealVar*>(p_arg);      
//                             if (!vecContains(rarg, self->_DependentTermParams[icomp])) {
//                                 self->_DependentTermParams[icomp].push_back(rarg);
//                                 self->_PrevParVals[icomp].push_back(-1E6);
//                             }
//                         }
//                     }
//                     delete parDeps; 
//                 } else {
//                     self->_IndependentTermProxy[icomp].push_back(term.get());
//                 }
//                 delete deps; 
//             }
//             icomp++;
//         }
//         self->_IntegratorInit = kTRUE;
//     }

//     Int_t BruComponentsPDF::getGenerator(const RooArgSet& directVars, RooArgSet &generateVars, Bool_t staticInitOK) const {	
//         if (!_DataCache) return 0; 
//         if (matchArgs(directVars, generateVars, VarSet(0))) return 1;
//         return 0;
//     }

 
  
//   Double_t BruComponentsPDF::evaluateData() const {
//         Double_t val = _BaseLine;
//         for (auto& comp: _Components) {
//             Double_t product = 1;
//             for (auto& term: comp) {
//                 product *= *term.get(); 
//             }
//             val += product; 
//         }
//         return val;
//     }
    
//     void BruComponentsPDF::RedirectServersToPdf() const {
//         for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
//             for (const auto& term : _DependentTermProxy[icomp]) {
//                 auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
//                 unconstTerm->recursiveRedirectServers(_IntegrateSet);
//             }	
//         }	
//     }

//     void BruComponentsPDF::RedirectServersToData() const {
//         for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
//             for (const auto& term : _DependentTermProxy[icomp]) {
//                 auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
//                 unconstTerm->recursiveRedirectServers(_ActualObs);
//             }	
//         }	
//     }

//     void BruComponentsPDF::HistIntegrals(const char* rangeName) const {
//         SetupIntegrator(); 
//         for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
//             for (const auto& term : _DependentTermProxy[icomp]) {
//                 auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
//                 unconstTerm->recursiveRedirectServers(_IntegrateSet);
//             }
//         }
     
//         BruEventsPDF::HistIntegrals(rangeName);
     
//         for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
//             for (const auto& term : _DependentTermProxy[icomp]) {
//                 auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
//                 unconstTerm->recursiveRedirectServers(_ActualObs);
//             }
//         }
//     }

//     Double_t BruComponentsPDF::cacheMCAP(const std::vector<Double_t> *vars, const std::vector<Int_t> *cats) const {
//         _MCAPDepTerm.resize(_Napd);
//         const Double_t* rawVars = vars->data();
//         const Int_t* rawCats = cats->data();
        
//         for (Long64_t iev = 0; iev < _Napd; ++iev) {
//             for (Int_t ii = 0; ii < _Nvars; ii++) {
//                 _IntegrateObs[ii]->setVal(rawVars[iev * _Nvars + ii]);
//             }
//             for (Int_t ii = 0; ii < _Ncats; ii++) {
//                 _IntegrateCats[ii]->setIndex(rawCats[iev * _Ncats + ii]);
//             }
      
//             Double_t val = _BaseLine;
//             Int_t icomp = 0;
//             for (auto& comp : _Components) {
//                 Double_t depprod = 1; 
//                 for (auto& term : _DependentTermProxy[icomp]) {
//                     depprod *= *term;
//                 }
//                 _MCAPDepTerm[iev].push_back(depprod);

//                 Double_t product = 1;
//                 for (auto& term : comp) {
//                     product *= *term.get(); 
//                 }
//                 val += product; 
//                 ++icomp;
//             }
//         }
//         return 0.;
//     }

//     Double_t BruComponentsPDF::evaluateMCAP() const {
//         Double_t val = _BaseLine;
//         Int_t icomp = 0;
     
//         for (auto& comp : _Components) {
//             Double_t product = _MCAPDepTerm[_TreeEntry][icomp];
//             for (auto& term : _IndependentTermProxy[icomp]) {
//                 product *= *term;
//             }
//             val += product; 
//             ++icomp;
//         }
//         return val;
//     }
    
//   Double_t BruComponentsPDF::evaluateMC(const std::vector<Double_t> *vars, const std::vector<Int_t> *cats) const {
//     if (_assertPostive) return evaluateMCAP();
    
//     const Double_t* rawVars = vars->data();
//     const Int_t* rawCats = cats->data();
    
//       for (Int_t ii = 0; ii < _Nvars; ii++) {
//       static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->setVal(rawVars[_TreeEntry * _Nvars + ii]);
//     }
//     for (Int_t ii = 0; ii < _Ncats; ii++) {
//       static_cast<RooCategory*>(&_CatSet[ii]->arg())->setIndex(rawCats[_TreeEntry * _Ncats + ii]);
//     }
//     return evaluateData();
//   }

 

//     bool BruComponentsPDF::isDirectGenSafe(const RooAbsArg& arg) const {
//         if (_ActualObs.find(arg.GetName())) return kTRUE;
//         if (_ActualCats.find(arg.GetName())) return kTRUE;
//         return kFALSE;
//     }

//     void BruComponentsPDF::initGenerator(Int_t code) {
//         SetupIntegrator();
//         RedirectServersToPdf();
//         BruEventsPDF::initGenerator(code);
//     }

//     void BruComponentsPDF::DoFirstIntegrations(const char* rangeName) const {
//         SetupIntegrator(); 
//         for (UInt_t icomp = 0; icomp < _NComps; icomp++)
//             _RecalcComponent.push_back(icomp);
        
//         RecalcComponentIntegrals(0, rangeName);

//         if (_UseSamplingIntegral == kTRUE)
//             RecalcComponentIntegralsSampling(0, rangeName);
        
//         _FirstCalculation = kFALSE;
//     }
    
//     Double_t BruComponentsPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
//         // BUG FIX: Safely delegate plotting (code > 1) before the constant override
//         if (code != 1) return BruEventsPDF::analyticalIntegral(code, rangeName);
//         if (_ForceConstInt) { _Last[0] = 1.0; return _Last[0]; }
//         if (!CheckChange()) return _Last[0];
//         if (!_DataCache) { Fatal("BruComponentsPDF::analyticalIntegral", "MC Cache is null!"); }

//         SetupIntegrator();

//         auto check = AssertPositivePDF();
//         if (check == kFALSE) return _Last[0] = 0;
       
//         if (_FirstCalculation == kTRUE) DoFirstIntegrations();
       
//         if (_WeightedBaseLine == 0 && _BaseLine != 0 && _DataCache->_UseEvWeights)
//             CalcWeightedBaseLine(rangeName);
//         else
//             _WeightedBaseLine = _BaseLine;

//         Bool_t needRecalc = kFALSE;
//         _RecalcComponent.clear();

//         for (UInt_t icomp = 0; icomp < _NComps; icomp++) {
//             if (_DependentTermProxy[icomp].size()) {
//                 UInt_t ipar = 0;
//                 for (auto par : _DependentTermParams[icomp]) {
//                     Double_t previous = _PrevParVals[icomp][ipar];
//                     Double_t pval = par->getVal();
             
//                     if (pval != previous) { 
//                         needRecalc = kTRUE;
//                         if (!vecContains(icomp, _RecalcComponent)) _RecalcComponent.push_back(icomp);
//                     }
//                     ipar++;
//                 }
//             }
//         }
     
//         if (needRecalc) RecalcComponentIntegrals(code, rangeName);

//         Double_t integral = _WeightedBaseLine;
//         for (UInt_t icomp = 0; icomp < _NComps; icomp++)
//             integral += componentIntegral(icomp);

//         _Last[0] = integral;
//         return integral;
//     }
     
//  void BruComponentsPDF::CalcWeightedBaseLine(const char* rangeName) const {
//         if (!_DataCache) return;
//         Long64_t ilow, ihigh = 0;
//         SetLowHighVals(ilow, ihigh);
      
//          _WeightedBaseLine = 0;
//         Long64_t accepted = 0;   
//         bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);

//         for (Long64_t ie = ilow; ie < ihigh; ie++) {
//             _TreeEntry = ie;
//             if (hasRange && !CheckRange(rangeName)) continue; 
            
//             accepted++;
//             _WeightedBaseLine += _DataCache->GetWeight(ie);
//         }

//         if (accepted > 0) _WeightedBaseLine /= accepted;
//     }

//  void BruComponentsPDF::RecalcComponentIntegrals(Int_t code, const char* rangeName) const {
//         if (!_DataCache) return;
//         Long64_t ilow, ihigh = 0;
//         SetLowHighVals(ilow, ihigh);
      
//         for (const auto& icomp : _RecalcComponent) {
//             _CacheCompDepIntegral[icomp] = 0.0;
//             // FUNDAMENTAL FIX: Deleted recursiveRedirectServers!
//         }

//         Long64_t accepted = 0;
//         bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);
//         const auto& vecReal = _DataCache->_vecReal;
//         const auto& vecCat = _DataCache->_vecCat;

//         // 1. SAVE OBSERVABLE STATE
//         std::vector<Double_t> savedObs(_Nvars);
//         std::vector<Int_t> savedCats(_Ncats);
//         for(Int_t ii=0; ii<_Nvars; ii++) savedObs[ii] = static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->getVal();
//         for(Int_t ii=0; ii<_Ncats; ii++) savedCats[ii] = static_cast<RooCategory*>(&_CatSet[ii]->arg())->getIndex();

//         // 2. RUN FAST LOOP
//         for (Long64_t ie = ilow; ie < ihigh; ie++) {
//             _TreeEntry = ie;
//             if (hasRange && !CheckRange(rangeName)) continue; 
            
//             accepted++;
        
//             for (Int_t ii = 0; ii < _Nvars; ii++)
//                 static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->setVal(vecReal[ie * _Nvars + ii]);
//             for (Int_t ii = 0; ii < _Ncats; ii++)
//                 static_cast<RooCategory*>(&_CatSet[ii]->arg())->setIndex(vecCat[ie * _Ncats + ii]);
          
//             Double_t weight = _DataCache->GetWeight(ie);
            
//             for (const auto& icomp : _RecalcComponent) {
//                 Double_t product = 1.0;
//                 for (const auto& term : _DependentTermProxy[icomp]) { product *= *term; }
//                 _CacheCompDepIntegral[icomp] += (product * weight);
//             }
//         }
       
//         for (const auto& icomp : _RecalcComponent) {
//             if (accepted > 0) _CacheCompDepIntegral[icomp] = _CacheCompDepIntegral[icomp] / accepted;
        
//             UInt_t ipar = 0;
//             for (auto par : _DependentTermParams[icomp]) {
//                 _PrevParVals[icomp][ipar] = par->getVal();
//                 ipar++;
//             }
//         }
      
//         // 3. RESTORE OBSERVABLE STATE
//         for(Int_t ii=0; ii<_Nvars; ii++) static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->setVal(savedObs[ii]);
//         for(Int_t ii=0; ii<_Ncats; ii++) static_cast<RooCategory*>(&_CatSet[ii]->arg())->setIndex(savedCats[ii]);
//     }
    
//     void BruComponentsPDF::RecalcComponentIntegralsSampling(Int_t code, const char* rangeName) const {
//         if (_RecalcComponent.empty() == kTRUE || !_DataCache) return;
      
//         Long64_t ilow, ihigh = 0;
//         SetLowHighVals(ilow, ihigh);
      
//         for (const auto& icomp : _RecalcComponent) {
//             _CacheCompDepIntegral[icomp] = 0.0;
//             _CacheCompDepSigmaIntegral[icomp] = 0.0;
//             // FUNDAMENTAL FIX: Deleted recursiveRedirectServers!
//         }

//         Long64_t accepted = 0, all = 0;
//         std::vector<Double_t> sumSquares(_NComps, 0.0);
//         bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);
//         const auto& vecReal = _DataCache->_vecReal;
//         const auto& vecCat = _DataCache->_vecCat;

//         // 1. SAVE OBSERVABLE STATE
//         std::vector<Double_t> savedObs(_Nvars);
//         std::vector<Int_t> savedCats(_Ncats);
//         for(Int_t ii=0; ii<_Nvars; ii++) savedObs[ii] = static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->getVal();
//         for(Int_t ii=0; ii<_Ncats; ii++) savedCats[ii] = static_cast<RooCategory*>(&_CatSet[ii]->arg())->getIndex();

//         // 2. RUN FAST LOOP
//         for (Long64_t ie = ilow; ie < ihigh; ie++) {
//             _TreeEntry = ie;
//             if (hasRange && !CheckRange(rangeName)) { ++all; continue; }
        
//             for (Int_t ii = 0; ii < _Nvars; ii++)
//                 static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->setVal(vecReal[ie * _Nvars + ii]);
//             for (Int_t ii = 0; ii < _Ncats; ii++)
//                 static_cast<RooCategory*>(&_CatSet[ii]->arg())->setIndex(vecCat[ie * _Ncats + ii]);
          
//             Double_t weight = _DataCache->GetWeight(ie);

//             for (const auto& icomp : _RecalcComponent) {
//                 Double_t product = 1.0;
//                 for (const auto& term : _DependentTermProxy[icomp]) { product *= *term; }
//                 product *= weight;
          
//                 _CacheCompDepIntegral[icomp] += product;
//                 sumSquares[icomp] += product * product;
//             }
//             ++accepted;
//             ++all;
//         }
      
//         for (const auto& icomp : _RecalcComponent) {
//             if (accepted > 1) {
//                 _CacheCompDepIntegral[icomp] = _CacheCompDepIntegral[icomp] / accepted;
//                 _CacheCompDepSigmaIntegral[icomp] = sumSquares[icomp] / accepted;
//             }
//         }
        
//         _NUsedForIntegral = accepted;
      
//         // 3. RESTORE OBSERVABLE STATE
//         for(Int_t ii=0; ii<_Nvars; ii++) static_cast<RooRealVar*>(&_ProxSet[ii]->arg())->setVal(savedObs[ii]);
//         for(Int_t ii=0; ii<_Ncats; ii++) static_cast<RooCategory*>(&_CatSet[ii]->arg())->setIndex(savedCats[ii]);
//     }
  
 
//   void BruComponentsPDF::RecalcComponentIntegralsSampling(Int_t code, const char* rangeName) const {
//         if (_RecalcComponent.empty() == kTRUE || !_DataCache) return;
      
//         Long64_t ilow, ihigh = 0;
//         SetLowHighVals(ilow, ihigh);
      
//         for (const auto& icomp : _RecalcComponent) {
//             _CacheCompDepIntegral[icomp] = 0.0;
//             _CacheCompDepSigmaIntegral[icomp] = 0.0;
//             for (const auto& term : _DependentTermProxy[icomp]) {
//                 auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
//                 unconstTerm->recursiveRedirectServers(_IntegrateSet);
//             }
//         }

//         Long64_t accepted = 0;
//         Long64_t all = 0;
//         std::vector<Double_t> sumSquares(_NComps, 0.0);
//         bool hasRange = (rangeName != nullptr && strlen(rangeName) > 0);

//         for (Long64_t ie = ilow; ie < ihigh; ie++) {
//             _TreeEntry = ie;
//             if (hasRange && !CheckRange(rangeName)) { ++all; continue; }
        
//             for (Int_t ii = 0; ii < _Nvars; ii++)
//                 _IntegrateObs[ii]->setVal(_DataCache->_vecReal[ie * _Nvars + ii]);
//             for (Int_t ii = 0; ii < _Ncats; ii++)
//                 _IntegrateCats[ii]->setIndex(_DataCache->_vecCat[ie * _Ncats + ii]);
          
//             Double_t weight = _DataCache->GetWeight(ie);

//             for (const auto& icomp : _RecalcComponent) {
//                 Double_t product = 1.0;
//                 for (const auto& term : _DependentTermProxy[icomp]) {
//                     product *= *term;
//                 }
//                 product *= weight;
          
//                 _CacheCompDepIntegral[icomp] += product;
//                 sumSquares[icomp] += product * product;
//             }
//             ++accepted;
//             ++all;
//         }
      
//         for (const auto& icomp : _RecalcComponent) {
//             if (accepted > 1) {
//                 _CacheCompDepIntegral[icomp] = _CacheCompDepIntegral[icomp] / accepted;
//                 _CacheCompDepSigmaIntegral[icomp] = sumSquares[icomp] / accepted;
//             }
//         }
        
//         _NUsedForIntegral = accepted;
      
//         for (const auto& icomp : _RecalcComponent) {
//             for (const auto& term : _DependentTermProxy[icomp]) {
//                 auto unconstTerm = const_cast<RooAbsReal*>(&term->arg());
//                 unconstTerm->recursiveRedirectServers(_ActualObs);
//             }
//         }
//     }
//     Double_t BruComponentsPDF::componentIntegral(Int_t icomp) const {
//         Double_t product = 1;
//         product *= _CacheCompDepIntegral[icomp];
//         for (auto& term : _IndependentTermProxy[icomp]) {
//             product *= *term;
//         }
//         return product; 
//     }

//     Double_t BruComponentsPDF::componentVariance(Int_t icomp) const {
//         Double_t product = 1;
//         product *= _CacheCompDepSigmaIntegral[icomp];
//         for (auto& term : _IndependentTermProxy[icomp]) {
//             product *= (*term) * (*term);
//         }
//         return product; 
//     }
    
//     Bool_t BruComponentsPDF::CheckChange() const {
//         Bool_t hasChanged = false;
//         for (Int_t i = 1; i < _Npars + 1; i++) {
//             if (_Last[i] != *(_ParSet[i - 1])) {
//                 hasChanged = true;
//             }
//         }
//         if (hasChanged) {
//             for (Int_t i = 1; i < _Npars + 1; i++) {
//                 _Last[i] = *(_ParSet[i - 1]);
//             }
//         }
//         return hasChanged;
//     }   

// } // namespace bru
