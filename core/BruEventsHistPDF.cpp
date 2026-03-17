/**
 * @file BruEventsHistPDF.cpp
 */

#include <Riostream.h> 

#include "BruEventsHistPDF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <cmath> 
#include <TMath.h> 
#include <TF1.h>
#include <TH1.h>
#include <TRandom3.h>

namespace bru {

    BruEventsHistPDF::BruEventsHistPDF(const char *name, const char *title, RooAbsReal& in_x, RooAbsReal& in_alpha, RooAbsReal& in_offset, RooAbsReal& in_scale, Int_t applySmooth, Int_t interp, Int_t xbins, Int_t nsamp, Int_t abins) :
      BruEventsPDF(name, title),
      x("x", "x", this, in_x),
      offset("offset", "offset", this, in_offset),
      scale("scale", "scale", this, in_scale),
      alpha("alpha", "alpha", this, in_alpha),
      _applySmooth(applySmooth),
      _Interpolate(interp),
      _NAlphaBins(abins),
      _NXBins0(xbins),
      _NIntSamples(nsamp)
    {
        MakeSets();
        x.SetName(in_x.GetName());
        offset.SetName(in_offset.GetName());
        scale.SetName(in_scale.GetName());
        alpha.SetName(in_alpha.GetName());
      
        auto *rx = dynamic_cast<RooRealVar*>(&in_x);
        auto *ra = dynamic_cast<RooRealVar*>(&in_alpha);
        auto *rs = dynamic_cast<RooRealVar*>(&in_scale);
        auto *ro = dynamic_cast<RooRealVar*>(&in_offset);
      
        // std::cout << "BruEventsHistPDF::BruEventsHistPDF using fudge parameters:" << std::endl;
        std::vector<RooRealVar*> vars = {ra, rs, ro};
        for (auto& var : vars) {
            if (var->getMax() == var->getMin()) {  
                var->setConstant();
            } else if (var->isConstant()) {  
                var->setMin(var->getVal());
                var->setMax(var->getVal());
            }
	    // var->Print();
        }
      
        Int_t NBins0 = _NXBins0;
        Double_t rsmin = 1;
        if (rs->getMin()) rsmin = rs->getMin();
        else std::cout << "BruEventsHistPDF::BruEventsHistPDF Warning no scale minimum set take = 1" << std::endl;
      
        Double_t mid = (rx->getMax() + rx->getMin()) / 2;
        Double_t diff = (rx->getMax() - rx->getMin()) / 2;
        Double_t rMin = mid - diff / rsmin + ro->getMin(); 
        Double_t rMax = mid + diff / rsmin + ro->getMax();
        Int_t NbinX = NBins0 / rsmin;
        if (NbinX < 10) NbinX = 10; 
      
	//  std::cout << "BruEventsHistPDF::BruEventsHistPDF using binning (" << NbinX << ", " << rMin << ", " << rMax << ") for x-axis variable '" << in_x.GetName() << "' of PDF '" << name << "'" << std::endl;

        Construct2DHist(TAxis(NbinX, rMin, rMax), TAxis(_NAlphaBins, ra->getMin(), ra->getMax()), ra->isConstant());

        _x_off = new RooRealVar(TString("off") + in_x.GetName(), "Vx_off", mid, rMin, rMax);
        _alphaVar = new RooRealVar("Valpha", "Valpha", 0, alpha.min(), alpha.max());

        if (not ra->isConstant()) _AlphaConstr = new RooGaussian(TString("AlphaConstr") + GetName(), "AlphaConstr", in_alpha, RooFit::RooConst(ra->getVal()), RooFit::RooConst(ra->getMax() / 5));
        if (not ro->isConstant()) _OffConstr = new RooGaussian(TString("OffConstr") + GetName(), "OffConstr", in_offset, RooFit::RooConst(ro->getVal()), RooFit::RooConst((ro->getMax() - ro->getMin()) / 5 / 2));
        if (not rs->isConstant()) _ScaleConstr = new RooGaussian(TString("ScConstr") + GetName(), "ScConstr", in_scale, RooFit::RooConst(rs->getVal()), RooFit::RooConst((rs->getMax() - rs->getMin()) / 5 / 2));
    } 

    void BruEventsHistPDF::Construct2DHist(TAxis xaxis, TAxis AlphAxis, Bool_t isAlphaConst) {
        if (_RHist != nullptr) { delete _RHist; _RHist = nullptr; }
      
        Int_t NAlphBins = AlphAxis.GetNbins();
        Int_t NbinX = xaxis.GetNbins();
        auto xedges = GetBinVector(xaxis);    
        
        if (isAlphaConst) {
            _RHist = new TH2D(TString("hmc_model_") + x.GetName() + GetName(), TString("MC model for ") + x.GetName(), NbinX, xedges.data(), 1, AlphAxis.GetXmin() - 1, AlphAxis.GetXmax() + 1);
        } else {
            _RHist = new TH2D(TString("hmc_model_") + x.GetName() + GetName(), TString("MC model for ") + x.GetName(), NbinX, xedges.data(), NAlphBins + 1, AlphAxis.GetXmin() - AlphAxis.GetBinCenter(1), AlphAxis.GetXmax() + AlphAxis.GetBinCenter(1));
        }
	_RHist->SetDirectory(nullptr);
    }
    
    BruEventsHistPDF::BruEventsHistPDF(const BruEventsHistPDF& other, const char* name) :  
        BruEventsPDF(other, name),
        x("x", this, other.x),
        offset("offset", this, other.offset),
        scale("scale", this, other.scale),
        alpha("alpha", this, other.alpha)
    {
        x.SetName(other.x.GetName());
        offset.SetName(other.offset.GetName());
        scale.SetName(other.scale.GetName());
        alpha.SetName(other.alpha.GetName());
      
        if (other._x_off) _x_off = dynamic_cast<RooRealVar*>(other._x_off->Clone());
        if (other._alphaVar) _alphaVar = dynamic_cast<RooRealVar*>(other._alphaVar->Clone());
        if (other._Hist) _Hist = dynamic_cast<RooDataHist*>(other._Hist->Clone(other._Hist->GetName()));
        if (other._RHist) _RHist = dynamic_cast<TH2D*>(other._RHist->Clone(other._RHist->GetName()));
        if (other._AlphaConstr) _AlphaConstr = dynamic_cast<RooGaussian*>(other._AlphaConstr->Clone());
        if (other._OffConstr) _OffConstr = dynamic_cast<RooGaussian*>(other._OffConstr->Clone());
        if (other._ScaleConstr) _ScaleConstr = dynamic_cast<RooGaussian*>(other._ScaleConstr->Clone());
      
        _VarMax = other._VarMax;
        _applySmooth = other._applySmooth;
        _Interpolate = other._Interpolate;
        _NAlphaBins = other._NAlphaBins;
        _NXBins0 = other._NXBins0;
        _NIntSamples = other._NIntSamples;
        _UseHistGenerator = other._UseHistGenerator;

        if (_RHist) _RHist->SetDirectory(nullptr);
      
        MakeSets(); 

        _Last = other._Last;
        _LastLength = other._LastLength;
    }

    BruEventsHistPDF::~BruEventsHistPDF() {
        if (_Hist) delete _Hist;
        if (_RHist) delete _RHist; 
        if (_x_off) delete _x_off;
        if (_alphaVar) delete _alphaVar;
        if (_AlphaConstr) delete _AlphaConstr;
        if (_OffConstr) delete _OffConstr;
        if (_ScaleConstr) delete _ScaleConstr;
    }

    void BruEventsHistPDF::MakeSets() {
        _ProxSet.clear();
        _VarSet.clear();
        _ParSet.clear();
        _ProxSet.push_back(&x);
        _ParSet.push_back(&offset);
        _ParSet.push_back(&scale);
        _ParSet.push_back(&alpha);
        InitSets();
    }

    Double_t BruEventsHistPDF::evaluate() const {
        if (!_Hist) return 1;
        Double_t arg = (x - _VarMax) * scale + _VarMax;
        arg = arg - offset;
        _x_off->setVal(arg);
        _alphaVar->setVal(Double_t(alpha));
        auto result = _Hist->weight(RooArgSet(*_x_off, *_alphaVar), _Interpolate, kFALSE);
        return result > 0 ? result : 0;
    } 

    Double_t BruEventsHistPDF::evaluateMC(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const {
        Double_t mcx = (*vars)[0];
        return evaluateMC(mcx);  
    }

    Double_t BruEventsHistPDF::evaluateMC(Double_t mcx) const {
        Double_t arg = (mcx - _VarMax) * scale + _VarMax;
        arg = arg - offset;
        _x_off->setVal(arg);
        _alphaVar->setVal(Double_t(alpha));
        auto result = _Hist->weight(RooArgSet(*_x_off, *_alphaVar), _Interpolate, kFALSE);
        return result > 0 ? result : 0;
    }

    Bool_t BruEventsHistPDF::SetEvTree(TTree* tree, TString cut, TTree* MCGenTree) {
        Bool_t OK = BruEventsPDF::SetEvTree(tree, cut, MCGenTree);
        if (!_Hist) CreateHistPdf();
        return (OK && _Hist) ? kTRUE : kFALSE;
    }

    void BruEventsHistPDF::CreateHistPdf() {
        _ConstInt = _DataCache ? _DataCache->_NTreeEntries : 0;
        TH1D his1("his1D", "his1D", _RHist->GetXaxis()->GetNbins(), _RHist->GetXaxis()->GetXmin(), _RHist->GetXaxis()->GetXmax());
        FillBase1DHist(his1);
  
        CheckForNegativeBins(his1);
        if (_applySmooth) his1.Smooth(_applySmooth);

        for (Int_t jtemp = 1; jtemp <= _RHist->GetNbinsX(); jtemp++) {
            _RHist->Fill(_RHist->GetXaxis()->GetBinCenter(jtemp), _RHist->GetYaxis()->GetBinCenter(1), his1.GetBinContent(jtemp));
        }

        TF1 gausnX("gausnX", "gausn(0)", _RHist->GetXaxis()->GetXmin(), _RHist->GetXaxis()->GetXmax());
        const auto NbinsX = _RHist->GetNbinsX();
        for (Int_t ia = 2; ia <= _RHist->GetNbinsY(); ia++) {
            Double_t vAlphb = _RHist->GetYaxis()->GetBinCenter(ia);
            for (Int_t itemp = 1; itemp <= NbinsX; itemp++) {
                Double_t vari = _RHist->GetXaxis()->GetBinCenter(itemp);
                Double_t NX = his1.GetBinContent(itemp);
                if (!NX) continue;
                gausnX.SetParameters(NX, vari, vAlphb);
                for (Int_t jtemp = 1; jtemp <= NbinsX; jtemp++) {
                    Double_t varj = _RHist->GetXaxis()->GetBinCenter(jtemp);
                    _RHist->Fill(varj, vAlphb, gausnX.Eval(varj));
                }
            }
        }

        for (Int_t ia = 1; ia <= _RHist->GetNbinsY(); ia++) {
            Double_t max_cont = 0;
            for (Int_t jtemp = 1; jtemp <= _RHist->GetNbinsX(); jtemp++) {
                Double_t cont = _RHist->GetBinContent(jtemp, ia);
                if (cont > max_cont) max_cont = cont;
            }
            for (Int_t jtemp = 1; jtemp <= _RHist->GetNbinsX(); jtemp++) {
                Double_t cont = _RHist->GetBinContent(jtemp, ia);
                if (cont == 0) _RHist->SetBinContent(jtemp, ia, 1E-10 * max_cont);
            }
        }

        Int_t bx, by, bz;
        _RHist->GetBinXYZ(_RHist->GetMaximumBin(), bx, by, bz);
        _VarMax = _RHist->GetXaxis()->GetBinCenter(bx);
      
        _Hist = new RooDataHist(_RHist->GetName(), _RHist->GetName(), RooArgSet(*_x_off, *_alphaVar), RooFit::Import(*_RHist));
    }

    Int_t BruEventsHistPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
        return BruEventsPDF::getAnalyticalIntegral(allVars, analVars, rangeName);
    }

    Double_t BruEventsHistPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
        if (code == 1) {
            if (!CheckChange()) return _Last[0];
            Double_t integral = 0;
            Double_t min = _RHist->GetXaxis()->GetXmin();
            Double_t max = _RHist->GetXaxis()->GetXmax();
            Double_t delta = (max - min) / _NIntSamples;
            auto var = (RooRealVar*)(&(x.arg()));
        
            Double_t rangeMin = var->getMin(rangeName);
            Double_t rangeMax = var->getMax(rangeName);
        
            for (Int_t ie = 1; ie <= _NIntSamples; ie++) {
                Double_t val = min + delta * ie;
                if (val < rangeMin || val > rangeMax) continue;
                integral += evaluateMC(val);
            }
            _Last[0] = integral * delta;
            return _Last[0];
        }
        return 1; 
    }

    void BruEventsHistPDF::FillBase1DHist(TH1D& his1) {
        Long64_t NFT = _DataCache ? _DataCache->_NTreeEntries : 0;
        his1.SetDirectory(nullptr);
        for (Int_t itr = 0; itr < NFT; itr++) {
            _TreeEntry = itr;
            if (_DataCache && _DataCache->_vecReal.size() > 0) {
                Double_t tvar = _DataCache->_vecReal[_TreeEntry * _Nvars + 0];
                his1.Fill(tvar, GetIntegralWeight(itr));
            }
        }
        _TreeEntry = 0;
    }

    void BruEventsHistPDF::CheckForNegativeBins(TH1D& his1) {
        for (Int_t ix = 0; ix < his1.GetEntries(); ix++) {
            if (his1.GetBinContent(ix) < 0) {
                std::cout << " BruEventsHistPDF::CreateHistPdf() weights have resulted in some -ve bins. Will try averaging..." << std::endl;
                Double_t b0, b1, b2;
                b0 = (ix == 1) ? his1.GetBinContent(ix) : his1.GetBinContent(ix - 1);
                b1 = his1.GetBinContent(ix);
                b2 = (ix == his1.GetEntries()) ? his1.GetBinContent(ix) : his1.GetBinContent(ix + 1);
      
                Double_t bmean = (b0 + b1 + b2) / 3;
                if (bmean < 0) bmean = 0;
                his1.SetBinContent(ix - 1, bmean);
                his1.SetBinContent(ix, bmean);
                his1.SetBinContent(ix + 1, bmean);
            }
        }
    }
    
    void BruEventsHistPDF::ResetTree() {
        BruEventsPDF::ResetTree();
        if (_Hist) {
            delete _Hist;
            _Hist = nullptr;
            _RHist->Reset();
        }
    }

    void BruEventsHistPDF::generateEvent(Int_t code) {
        if (_UseHistGenerator) {
            Bool_t inRange = kFALSE;
            Double_t genx = 0;
            auto var = static_cast<const RooRealVar*>(&x.arg());
        
            Double_t rangeMin = var->getMin("");
            Double_t rangeMax = var->getMax("");
        
            while (inRange == kFALSE) {
                genx = _GenHist.GetRandom();
                Double_t arg = (genx - _VarMax) / scale + _VarMax;
                arg = arg + offset / scale;
          
                if (arg < rangeMin || arg > rangeMax) { continue; }
                inRange = kTRUE;
                x = arg; 
            }
        } else {
            BruEventsPDF::generateEvent(code);
        }
    }
    
    void BruEventsHistPDF::initGenerator(Int_t code) {
        if (_UseHistGenerator == kFALSE) {
            BruEventsPDF::initGenerator(code);
            return;
        }
        if (_GenHist.Integral() > 0) return;
      
        Double_t arg = (x - _VarMax) * scale + _VarMax;
        arg = arg - offset;
        _x_off->setVal(arg);
        _alphaVar->setVal(Double_t(alpha));
        auto gbin = _RHist->FindFixBin(arg, Double_t(alpha));
        Int_t xbin = 0, abin = 0, zbin = 0;
        _RHist->GetBinXYZ(gbin, xbin, abin, zbin);
        _GenHist = *_RHist->ProjectionX(TString("projX_") + _RHist->GetName(), abin, abin);
    }

} // namespace bru
