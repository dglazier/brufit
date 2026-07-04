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
#include <algorithm>
#include <stdexcept>

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
      
        std::vector<RooRealVar*> vars = {ra, rs, ro};
        for (auto& var : vars) {
            if (var->getMax() == var->getMin()) {  
                var->setConstant();
            } else if (var->isConstant()) {  
                var->setMin(var->getVal());
                var->setMax(var->getVal());
            }
        }
      
        Double_t mid = (rx->getMax() + rx->getMin()) / 2;
        Double_t diff = (rx->getMax() - rx->getMin()) / 2;
        Double_t rsmin = (rs->getMin() > 0) ? rs->getMin() : 1.0;
        Double_t rMin = mid - diff / rsmin + ro->getMin(); 
        Double_t rMax = mid + diff / rsmin + ro->getMax();
        Int_t NbinX = std::max(10, (int)(_NXBins0 / rsmin));
      
        Construct2DHist(TAxis(NbinX, rMin, rMax), TAxis(_NAlphaBins, ra->getMin(), ra->getMax()), ra->isConstant());

        _x_off = new RooRealVar(TString("off") + in_x.GetName(), "Vx_off", mid, rMin, rMax);
        _alphaVar = new RooRealVar("Valpha", "Valpha", 0, alpha.min(), alpha.max());

        if (not ra->isConstant()) _AlphaConstr = new RooGaussian(TString("AlphaConstr") + GetName(), "AlphaConstr", in_alpha, RooFit::RooConst(ra->getVal()), RooFit::RooConst(ra->getMax() / 5));
        if (not ro->isConstant()) _OffConstr = new RooGaussian(TString("OffConstr") + GetName(), "OffConstr", in_offset, RooFit::RooConst(ro->getVal()), RooFit::RooConst((ro->getMax() - ro->getMin()) / 5 / 2));
        if (not rs->isConstant()) _ScaleConstr = new RooGaussian(TString("ScConstr") + GetName(), "ScConstr", in_scale, RooFit::RooConst(rs->getVal()), RooFit::RooConst((rs->getMax() - rs->getMin()) / 5 / 2));
    } 

    void BruEventsHistPDF::Construct2DHist(TAxis xaxis, TAxis AlphAxis, Bool_t isAlphaConst) {
        if (_RHist != nullptr) { delete _RHist; _RHist = nullptr; }
        Int_t NbinX = xaxis.GetNbins();
        auto xedges = GetBinVector(xaxis);    
        if (isAlphaConst) {
            _RHist = new TH2D(TString("hmc_model_") + x.GetName() + GetName(), TString("MC model"), NbinX, xedges.data(), 1, AlphAxis.GetXmin() - 1, AlphAxis.GetXmax() + 1);
        } else {
            _RHist = new TH2D(TString("hmc_model_") + x.GetName() + GetName(), TString("MC model"), NbinX, xedges.data(), AlphAxis.GetNbins() + 1, AlphAxis.GetXmin() - AlphAxis.GetBinCenter(1), AlphAxis.GetXmax() + AlphAxis.GetBinCenter(1));
        }
        _RHist->SetDirectory(nullptr);
    }
    
    BruEventsHistPDF::BruEventsHistPDF(const BruEventsHistPDF& other, const char* name) :  
        BruEventsPDF(other, name),
        x("x", this, other.x),
        offset("offset", this, other.offset),
        scale("scale", this, other.scale),
        alpha("alpha", this, other.alpha),
        _cachedBins(other._cachedBins),
        _xMin(other._xMin),
        _xWidth(other._xWidth),
        _nx(other._nx),
        _yMin(other._yMin),
        _yWidth(other._yWidth),
        _ny(other._ny)
    {
        x.SetName(other.x.GetName());
        offset.SetName(other.offset.GetName());
        scale.SetName(other.scale.GetName());
        alpha.SetName(other.alpha.GetName());

        if (other._x_off) _x_off = (RooRealVar*)other._x_off->Clone();
        if (other._alphaVar) _alphaVar = (RooRealVar*)other._alphaVar->Clone();
        if (other._Hist) _Hist = (RooDataHist*)other._Hist->Clone(other._Hist->GetName());
        if (other._RHist) _RHist = (TH2D*)other._RHist->Clone(other._RHist->GetName());
        if (other._AlphaConstr) _AlphaConstr = (RooGaussian*)other._AlphaConstr->Clone();
        if (other._OffConstr) _OffConstr = (RooGaussian*)other._OffConstr->Clone();
        if (other._ScaleConstr) _ScaleConstr = (RooGaussian*)other._ScaleConstr->Clone();
      
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
        _ProxSet.clear(); _VarSet.clear(); _ParSet.clear();
        _ProxSet.push_back(&x);
        _ParSet.push_back(&offset); _ParSet.push_back(&scale); _ParSet.push_back(&alpha);
        InitSets();
    }

    void BruEventsHistPDF::initializeCache() {
        if (!_Hist) return;
        auto* vars = (RooArgSet*)_Hist->get();
        RooRealVar* xVar = dynamic_cast<RooRealVar*>((*vars)[0]);
        RooRealVar* yVar = dynamic_cast<RooRealVar*>((*vars)[1]);
        
        if (!xVar || !yVar) {
            std::cout << "BruEventsHistPDF ERROR: Axis mismatch in " << GetName() << std::endl;
            return;
        }

        _nx = xVar->numBins(); _xMin = xVar->getMin(); _xWidth = (xVar->getMax() - _xMin) / _nx;
        _ny = yVar->numBins(); _yMin = yVar->getMin(); _yWidth = (yVar->getMax() - _yMin) / _ny;
        _cachedBins.resize(_nx * _ny);

        for (int iy = 0; iy < _ny; ++iy) {
            yVar->setBin(iy);
            for (int ix = 0; ix < _nx; ++ix) {
                xVar->setBin(ix);
                _cachedBins[iy * _nx + ix] = _Hist->weight(*vars, 0, kFALSE);
            }
        }
    }

    Double_t BruEventsHistPDF::evaluate() const {
        if (!_Hist) return 1;
        Double_t arg = (x - _VarMax) * scale + _VarMax - offset;
        _x_off->setVal(arg); _alphaVar->setVal((double)alpha);
        auto result = _Hist->weight(RooArgSet(*_x_off, *_alphaVar), _Interpolate, kFALSE);
        return result > 0 ? result : 0;
    } 

    void BruEventsHistPDF::doEval(RooFit::EvalContext & ctx) const {
        if (!_Hist || _cachedBins.empty()) {
            std::fill(ctx.output().begin(), ctx.output().end(), 1.0); return;
        }

        auto xData = ctx.at(x); auto output = ctx.output();
        const double cSc = scale, cOff = offset, vMax = _VarMax, cAl = alpha;
        const double* bins = _cachedBins.data(); const int nx = _nx;

        if (_Interpolate) {
            // --- Y-AXIS EXTRAPOLATION ---
            double cY = (cAl - _yMin) / _yWidth - 0.5;
            cY = std::max(-0.5, std::min(cY, _ny - 0.5)); // Allow extrapolation to physical edges
            int iy = std::floor(cY);
            iy = std::max(0, std::min(iy, _ny - 2));      // Safe array index
            double yFr = cY - iy;
            
            const double wY0 = 1.0 - yFr, wY1 = yFr;
            const int r0 = iy * nx;
            const int r1 = std::min(iy + 1, _ny - 1) * nx;

            for (size_t i = 0; i < xData.size(); ++i) {
                double arg = (xData[i] - vMax) * cSc + vMax - cOff;
                
                // --- X-AXIS EXTRAPOLATION ---
                double cX = (arg - _xMin) / _xWidth - 0.5;
                cX = std::max(-0.5, std::min(cX, nx - 0.5)); 
                int ix = std::floor(cX);
                ix = std::max(0, std::min(ix, nx - 2));      
                double xFr = cX - ix;

                int ix0 = ix;
                int ix1 = std::min(ix + 1, nx - 1);

                double b00 = bins[r0 + ix0], b10 = bins[r0 + ix1];
                double b01 = bins[r1 + ix0], b11 = bins[r1 + ix1];

                double vL = b00 * wY0 + b01 * wY1;
                double vR = b10 * wY0 + b11 * wY1;
                double res = vL + xFr * (vR - vL);
                
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        } else {
            int iy_nearest = std::max(0, std::min((int)((cAl - _yMin) / _yWidth), _ny - 1));
            const int r_nearest = iy_nearest * nx;

            for (size_t i = 0; i < xData.size(); ++i) {
                double arg = (xData[i] - vMax) * cSc + vMax - cOff;
                int ix_nearest = std::max(0, std::min((int)((arg - _xMin) / _xWidth), nx - 1));
                
                double res = bins[r_nearest + ix_nearest];
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        }
    }

    Double_t BruEventsHistPDF::evaluateMC(const std::vector<Float_t> *vars, const std::vector<Int_t> *cats) const {
        return evaluateMC((*vars)[0]);  
    }

    Double_t BruEventsHistPDF::evaluateMC(Double_t mcx) const {
        Double_t arg = (mcx - _VarMax) * scale + _VarMax - offset;
        _x_off->setVal(arg); _alphaVar->setVal((double)alpha);
        auto result = _Hist->weight(RooArgSet(*_x_off, *_alphaVar), _Interpolate, kFALSE);
        return result > 0 ? result : 0;
    }

    void BruEventsHistPDF::evaluateMCBatch(const std::vector<Double_t>& mcx_array, std::vector<Double_t>& output) const {
        size_t nEvents = mcx_array.size();
        if (output.size() != nEvents) output.resize(nEvents);
        
        if (!_Hist || _cachedBins.empty()) {
            std::fill(output.begin(), output.end(), 1.0); return;
        }

        const double cSc = scale, cOff = offset, vMax = _VarMax, cAl = alpha;
        const double* bins = _cachedBins.data(); 
        const int nx = _nx;

        if (_Interpolate) {
            double cY = (cAl - _yMin) / _yWidth - 0.5;
            cY = std::max(-0.5, std::min(cY, _ny - 0.5));
            int iy = std::floor(cY);
            iy = std::max(0, std::min(iy, _ny - 2));
            double yFr = cY - iy;
            
            const double wY0 = 1.0 - yFr, wY1 = yFr;
            const int r0 = iy * nx;
            const int r1 = std::min(iy + 1, _ny - 1) * nx;

            for (size_t i = 0; i < nEvents; ++i) {
                double arg = (mcx_array[i] - vMax) * cSc + vMax - cOff;
                
                double cX = (arg - _xMin) / _xWidth - 0.5;
                cX = std::max(-0.5, std::min(cX, nx - 0.5));
                int ix = std::floor(cX);
                ix = std::max(0, std::min(ix, nx - 2));
                double xFr = cX - ix;

                int ix0 = ix;
                int ix1 = std::min(ix + 1, nx - 1);

                double b00 = bins[r0 + ix0], b10 = bins[r0 + ix1];
                double b01 = bins[r1 + ix0], b11 = bins[r1 + ix1];

                double vL = b00 * wY0 + b01 * wY1;
                double vR = b10 * wY0 + b11 * wY1;
                double res = vL + xFr * (vR - vL);
                
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        } else {
            int iy_nearest = std::max(0, std::min((int)((cAl - _yMin) / _yWidth), _ny - 1));
            const int r_nearest = iy_nearest * nx;

            for (size_t i = 0; i < nEvents; ++i) {
                double arg = (mcx_array[i] - vMax) * cSc + vMax - cOff;
                int ix_nearest = std::max(0, std::min((int)((arg - _xMin) / _xWidth), nx - 1));
                
                double res = bins[r_nearest + ix_nearest];
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        }
    }

    Double_t BruEventsHistPDF::analyticalIntegral(Int_t code, const char* rangeName) const {
        if (code == 1) {
            if (!CheckChange()) return _Last[0];
            Double_t min = _RHist->GetXaxis()->GetXmin(), max = _RHist->GetXaxis()->GetXmax();
            Double_t delta = (max - min) / _NIntSamples;
            auto var = (RooRealVar*)(&(x.arg()));
            Double_t rMin = var->getMin(rangeName), rMax = var->getMax(rangeName);
        
            std::vector<Double_t> pts; pts.reserve(_NIntSamples);
            for (Int_t i = 1; i <= _NIntSamples; ++i) {
                double v = min + delta * i;
                if (v >= rMin && v <= rMax) pts.push_back(v);
            }
            std::vector<Double_t> res; evaluateMCBatch(pts, res);
            Double_t integral = 0;
            for (auto r : res) integral += std::max(r, 1e-18);
            
            _Last[0] = integral * delta; 
            return (_Last[0] > 1e-18) ? _Last[0] : 1e-18;
        }
        return 1; 
    }

    void BruEventsHistPDF::CreateHistPdf() {
        _ConstInt = _DataCache ? _DataCache->_NTreeEntries : 0;
        TH1D his1("his1D", "his1D", _RHist->GetXaxis()->GetNbins(), _RHist->GetXaxis()->GetXmin(), _RHist->GetXaxis()->GetXmax());
        FillBase1DHist(his1); CheckForNegativeBins(his1);
        if (_applySmooth) his1.Smooth(_applySmooth);

        for (Int_t j = 1; j <= _RHist->GetNbinsX(); ++j)
            _RHist->Fill(_RHist->GetXaxis()->GetBinCenter(j), _RHist->GetYaxis()->GetBinCenter(1), his1.GetBinContent(j));

        TF1 gX("gX", "gausn(0)", _RHist->GetXaxis()->GetXmin(), _RHist->GetXaxis()->GetXmax());
        for (Int_t ia = 2; ia <= _RHist->GetNbinsY(); ++ia) {
            Double_t vA = _RHist->GetYaxis()->GetBinCenter(ia);
            for (Int_t i = 1; i <= _RHist->GetNbinsX(); ++i) {
                Double_t NX = his1.GetBinContent(i); if (!NX) continue;
                gX.SetParameters(NX, _RHist->GetXaxis()->GetBinCenter(i), vA);
                for (Int_t j = 1; j <= _RHist->GetNbinsX(); ++j)
                    _RHist->Fill(_RHist->GetXaxis()->GetBinCenter(j), vA, gX.Eval(_RHist->GetXaxis()->GetBinCenter(j)));
            }
        }
        for (Int_t ia = 1; ia <= _RHist->GetNbinsY(); ++ia) {
            Double_t mx = 0;
            for (Int_t j = 1; j <= _RHist->GetNbinsX(); ++j) mx = std::max(mx, _RHist->GetBinContent(j, ia));
            for (Int_t j = 1; j <= _RHist->GetNbinsX(); ++j) if (_RHist->GetBinContent(j,ia) == 0) _RHist->SetBinContent(j,ia,1E-10*mx);
        }

        Int_t bx, by, bz; _RHist->GetBinXYZ(_RHist->GetMaximumBin(), bx, by, bz);
        _VarMax = _RHist->GetXaxis()->GetBinCenter(bx);
        _Hist = new RooDataHist(_RHist->GetName(), _RHist->GetName(), RooArgSet(*_x_off, *_alphaVar), RooFit::Import(*_RHist));
        
        initializeCache(); 
    }

    Bool_t BruEventsHistPDF::SetEvTree(TTree* tree, TString cut, TTree* MCGenTree) {
        Bool_t OK = BruEventsPDF::SetEvTree(tree, cut, MCGenTree);
        if (!_Hist) CreateHistPdf(); return OK && _Hist;
    }

    void BruEventsHistPDF::ResetTree() {
        BruEventsPDF::ResetTree();
        if (_Hist) { delete _Hist; _Hist = nullptr; _RHist->Reset(); }
        _cachedBins.clear();
    }

    void BruEventsHistPDF::FillBase1DHist(TH1D& his1) {
        Long64_t NFT = _DataCache ? _DataCache->_NTreeEntries : 0;
        his1.SetDirectory(nullptr);
        for (Int_t itr = 0; itr < NFT; ++itr) {
            _TreeEntry = itr;
            if (_DataCache && !_DataCache->_vecReal.empty())
                his1.Fill(_DataCache->_vecReal[_TreeEntry * _Nvars], GetIntegralWeight(itr));
        }
        _TreeEntry = 0;
    }

    void BruEventsHistPDF::CheckForNegativeBins(TH1D& his1) {
        for (Int_t ix = 1; ix <= his1.GetNbinsX(); ++ix) {
            if (his1.GetBinContent(ix) < 0) {
                double b0 = his1.GetBinContent(std::max(1, ix-1));
                double b1 = his1.GetBinContent(ix);
                double b2 = his1.GetBinContent(std::min((int)his1.GetNbinsX(), ix+1));
                double avg = std::max(0.0, (b0+b1+b2)/3.0);
                his1.SetBinContent(ix, avg);
            }
        }
    }

    void BruEventsHistPDF::generateEvent(Int_t code) {
        if (_UseHistGenerator) {
            auto var = (RooRealVar*)&x.arg();
            double rMin = var->getMin(), rMax = var->getMax();
            while (true) {
                double gx = _GenHist.GetRandom();
                double arg = (gx - _VarMax) / scale + _VarMax + offset / scale;
                if (arg >= rMin && arg <= rMax) { x = arg; break; }
            }
        } else BruEventsPDF::generateEvent(code);
    }
    
    void BruEventsHistPDF::initGenerator(Int_t code) {
        if (!_UseHistGenerator) { BruEventsPDF::initGenerator(code); return; }
        if (_GenHist.Integral() > 0) return;
        _x_off->setVal((double)x - _VarMax * scale + _VarMax - offset);
        _alphaVar->setVal((double)alpha);
        int bin = _RHist->FindFixBin(_x_off->getVal(), _alphaVar->getVal());
        int bx, by, bz; _RHist->GetBinXYZ(bin, bx, by, bz);
        _GenHist = *_RHist->ProjectionX(TString("proj_")+_RHist->GetName(), by, by);
    }

    Int_t BruEventsHistPDF::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName) const {
        return BruEventsPDF::getAnalyticalIntegral(allVars, analVars, rangeName);
    }
  Double_t BruEventsHistPDF::GetRelativeVariance() const {
        if (!_TrackMCVariance || !_DataCache || _DataCache->_NTreeEntries <= 1) return 1E-12;
        
        Double_t sumW = 0;
        Double_t sumW2 = 0;
        
        for (Long64_t ie = 0; ie < _DataCache->_NTreeEntries; ie++) {
            Double_t w = GetIntegralWeight(ie);
            sumW += w;
            sumW2 += w * w;
        }
        
        if (sumW <= 0 || sumW2 <= 0) return 1E-12;
        
        Double_t nEff = (sumW * sumW) / sumW2;
        
        // Relative Standard Error of the generated template
        Double_t relError = std::sqrt(1.0 / nEff);
        return relError > 0 ? relError : 1E-12;
    }

} // namespace bru
