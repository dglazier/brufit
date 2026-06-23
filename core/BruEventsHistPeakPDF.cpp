/**
 * @file BruEventsHistPeakPDF.cpp
 */

#include <Riostream.h> 
#include "BruEventsHistPeakPDF.h" 
#include <RooAbsReal.h> 
#include <RooAbsCategory.h> 
#include <cmath> 
#include <TMath.h> 
#include <TF1.h>
#include <TH1.h>
#include <TRandom3.h>
#include <vector>
#include <algorithm>
#include <iostream>

namespace bru {

    BruEventsHistPeakPDF::BruEventsHistPeakPDF(const char *name, const char *title, RooAbsReal& in_x, RooAbsReal& in_alpha, RooAbsReal& in_offset, RooAbsReal& in_scale, Int_t applySmooth, Int_t interp, Int_t xbins, Int_t nsamp, Int_t abins) :
        BruEventsHistPDF(name, title, in_x, in_alpha, in_offset, in_scale, applySmooth, interp, xbins, nsamp, abins)
    {
    }

    BruEventsHistPeakPDF::BruEventsHistPeakPDF(const BruEventsHistPeakPDF& other, const char* name) :  
        BruEventsHistPDF(other, name),
        _lut(other._lut),
        _varCenters(other._varCenters),
        _nLut(other._nLut),
        _lutWidth(other._lutWidth)
    {
    }

    void BruEventsHistPeakPDF::initializeCache() {
        BruEventsHistPDF::initializeCache();
        
        if (!_Hist) return;
        auto* vars = (RooArgSet*)_Hist->get();
        RooRealVar* xVar = dynamic_cast<RooRealVar*>((*vars)[0]);
        
        _varCenters.resize(_nx);
        double minBinWidth = 1e9; 
        
        for (int ix = 0; ix < _nx; ++ix) {
            xVar->setBin(ix);
            _varCenters[ix] = xVar->getVal();
            
            double w = xVar->getBinWidth(ix);
            if (w < minBinWidth) minBinWidth = w; 
        }

        double totalRange = _nx * _xWidth; 
        _nLut = std::ceil(totalRange / (minBinWidth / 2.0));
        if (_nLut > 500000) _nLut = 500000; 
        if (_nLut < _nx * 2) _nLut = _nx * 2;
        
        _lut.resize(_nLut);
        _lutWidth = totalRange / _nLut; 
        
        for (int i = 0; i < _nLut; ++i) {
            double lut_x = _xMin + (i + 0.5) * _lutWidth;
            
            auto it = std::upper_bound(_varCenters.begin(), _varCenters.end(), lut_x);
            int left_center_idx = std::distance(_varCenters.begin(), it) - 1;
            
            _lut[i] = std::max(0, std::min(left_center_idx, _nx - 2)); 
        }
    }

    void BruEventsHistPeakPDF::doEval(RooFit::EvalContext & ctx) const {
        if (!_Hist || _cachedBins.empty()) {
            std::fill(ctx.output().begin(), ctx.output().end(), 1.0); return;
        }

        auto xData = ctx.at(x); 
        auto output = ctx.output();
        
        const double cSc = scale, cOff = offset, vMax = _VarMax, cAl = alpha;
        const double* bins = _cachedBins.data(); 
        const int nx = _nx;
        const int* lut = _lut.data();
        const double* centers = _varCenters.data();

        if (_Interpolate) {
            // Y-AXIS EXTRAPOLATION
            double cY = (cAl - _yMin) / _yWidth - 0.5;
            cY = std::max(-0.5, std::min(cY, _ny - 0.5));
            int iy = std::floor(cY);
            iy = std::max(0, std::min(iy, _ny - 2));
            double yFr = cY - iy;
            const double wY0 = 1.0 - yFr, wY1 = yFr;
            const int r0 = iy * nx, r1 = std::min(iy + 1, _ny - 1) * nx;

            for (size_t i = 0; i < xData.size(); ++i) {
                double arg = (xData[i] - vMax) * cSc + vMax - cOff;

                // 1. O(1) LUT Guess
                double lut_continuous = (arg - _xMin) / _lutWidth;
                int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
                int ix = lut[lut_idx];

                // 2. EXACT ALIGNMENT FIX (Prevents "LUT Cliffs" for Minuit)
                while (ix < nx - 2 && arg >= centers[ix + 1]) ix++;
                while (ix > 0 && arg < centers[ix]) ix--;

                double c0 = centers[ix];
                double c1 = centers[ix + 1];
                
                // 3. LINEAR EXTRAPOLATION FIX (Prevents flat plateaus at edges)
                double xFrac = (arg - c0) / (c1 - c0);
                xFrac = std::max(-0.5, std::min(xFrac, 1.5));

                double b00 = bins[r0 + ix], b10 = bins[r0 + ix + 1];
                double b01 = bins[r1 + ix], b11 = bins[r1 + ix + 1];

                double vL = b00 * wY0 + b01 * wY1;
                double vR = b10 * wY0 + b11 * wY1;
                double res = vL + xFrac * (vR - vL);
                
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        } else {
            int iy = std::max(0, std::min((int)((cAl - _yMin) / _yWidth), _ny - 1));
            const int r_nearest = iy * nx;

            for (size_t i = 0; i < xData.size(); ++i) {
                double arg = (xData[i] - vMax) * cSc + vMax - cOff;
                
                double lut_continuous = (arg - _xMin) / _lutWidth;
                int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
                int ix = lut[lut_idx];

                // EXACT ALIGNMENT FIX
                while (ix < nx - 2 && arg >= centers[ix + 1]) ix++;
                while (ix > 0 && arg < centers[ix]) ix--;

                double c0 = centers[ix];
                double c1 = centers[ix + 1];
                int nearest_ix = (arg - c0) < (c1 - arg) ? ix : (ix + 1);

                double res = bins[r_nearest + nearest_ix];
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        }
    }

    void BruEventsHistPeakPDF::evaluateMCBatch(const std::vector<Double_t>& mcx_array, std::vector<Double_t>& output) const {
        size_t nEvents = mcx_array.size();
        if (output.size() != nEvents) output.resize(nEvents);
        
        if (!_Hist || _cachedBins.empty()) {
            std::fill(output.begin(), output.end(), 1.0); return;
        }

        const double cSc = scale, cOff = offset, vMax = _VarMax, cAl = alpha;
        const double* bins = _cachedBins.data(); 
        const int nx = _nx;
        const int* lut = _lut.data();
        const double* centers = _varCenters.data();

        if (_Interpolate) {
            // Y-AXIS EXTRAPOLATION
            double cY = (cAl - _yMin) / _yWidth - 0.5;
            cY = std::max(-0.5, std::min(cY, _ny - 0.5));
            int iy = std::floor(cY);
            iy = std::max(0, std::min(iy, _ny - 2));
            double yFr = cY - iy;
            const double wY0 = 1.0 - yFr, wY1 = yFr;
            const int r0 = iy * nx, r1 = std::min(iy + 1, _ny - 1) * nx;

            for (size_t i = 0; i < nEvents; ++i) {
                double arg = (mcx_array[i] - vMax) * cSc + vMax - cOff;

                // 1. O(1) LUT Guess
                double lut_continuous = (arg - _xMin) / _lutWidth;
                int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
                int ix = lut[lut_idx];

                // 2. EXACT ALIGNMENT FIX
                while (ix < nx - 2 && arg >= centers[ix + 1]) ix++;
                while (ix > 0 && arg < centers[ix]) ix--;

                double c0 = centers[ix];
                double c1 = centers[ix + 1];
                
                // 3. LINEAR EXTRAPOLATION FIX
                double xFrac = (arg - c0) / (c1 - c0);
                xFrac = std::max(-0.5, std::min(xFrac, 1.5));

                double b00 = bins[r0 + ix], b10 = bins[r0 + ix + 1];
                double b01 = bins[r1 + ix], b11 = bins[r1 + ix + 1];

                double vL = b00 * wY0 + b01 * wY1;
                double vR = b10 * wY0 + b11 * wY1;
                double res = vL + xFrac * (vR - vL);
                
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        } else {
            int iy = std::max(0, std::min((int)((cAl - _yMin) / _yWidth), _ny - 1));
            const int r_nearest = iy * nx;

            for (size_t i = 0; i < nEvents; ++i) {
                double arg = (mcx_array[i] - vMax) * cSc + vMax - cOff;
                
                double lut_continuous = (arg - _xMin) / _lutWidth;
                int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
                int ix = lut[lut_idx];

                // EXACT ALIGNMENT FIX
                while (ix < nx - 2 && arg >= centers[ix + 1]) ix++;
                while (ix > 0 && arg < centers[ix]) ix--;

                double c0 = centers[ix];
                double c1 = centers[ix + 1];
                int nearest_ix = (arg - c0) < (c1 - arg) ? ix : (ix + 1);

                double res = bins[r_nearest + nearest_ix];
                output[i] = res > 1e-15 ? res : 1e-15;
            }
        }
    }

    void BruEventsHistPeakPDF::FillBase1DHist(TH1D& his1) {
        std::cout << "BruEventsHistPeakPDF::FillBase1DHist " << his1.GetName() << std::endl;
        
        BruEventsHistPDF::FillBase1DHist(his1);
        
        int Nbins0 = his1.GetNbinsX();
        double xmin0 = his1.GetBinLowEdge(0);
        double xmax0 = his1.GetBinLowEdge(his1.GetNbinsX() + 1);
        double xmin = xmin0;
        double xmax = xmax0;
        double range = xmax - xmin;
        double binwidth0 = range / Nbins0;
        double binwidth = binwidth0;
 
        auto maxBin = his1.GetMaximumBin();
        int deltaBin = his1.GetRMS() / binwidth; 
        
        if (maxBin - deltaBin < 0) deltaBin = maxBin;
        if (maxBin + deltaBin > Nbins0) deltaBin = Nbins0 - maxBin;
        
        xmin = his1.GetBinLowEdge(maxBin - deltaBin);
        xmax = his1.GetBinCenter(maxBin + deltaBin) + binwidth;
        
        range = xmax - xmin;
        binwidth = range / Nbins0;
        
        TH1D peakhist("peakhist", "peakhist", Nbins0, xmin, xmax);
        peakhist.SetDirectory(nullptr);
       
        BruEventsHistPDF::FillBase1DHist(peakhist);

        auto lastEdge = xmin;
        std::vector<double> bins;
        std::vector<double> fracs;
        bins.push_back(lastEdge);
        
        auto hmax = peakhist.GetMaximum();
        double sumFracs = 0.;
        double minFrac = 10;
        
        for (auto i = 1; i <= peakhist.GetNbinsX(); i++) {
            Double_t mainpeak = peakhist.GetBinContent(i);
            Double_t leftpeak = (i == 1) ? 0 : peakhist.GetBinContent(i - 1);
            Double_t rightpeak = (i == peakhist.GetNbinsX()) ? 0 : peakhist.GetBinContent(i + 1);
            
            auto binFrac = (4 * 1.1 * hmax - 2 * mainpeak - leftpeak - rightpeak);
            if (binFrac < 0) binFrac = 2;
            
            fracs.push_back(binFrac);
            if (binFrac < minFrac) minFrac = binFrac;
        }
        
        for (auto& frac : fracs) {
            if (frac == 2) frac = minFrac; 
            sumFracs += frac;
        }
        for (auto& frac : fracs) {
            frac /= sumFracs;
            lastEdge += frac * range;
            bins.push_back(lastEdge);
        }
     
        double edge = bins[0];
        Int_t gradual = 3;
        while (edge > xmin0) {
            double gradfactor = gradual > 0 ? 1. / gradual-- : 1;
            edge -= binwidth0 * gradfactor;
            bins.push_back(edge);
            std::rotate(bins.rbegin(), bins.rbegin() + 1, bins.rend());
        }
        
        edge = bins.back();
        gradual = 3;
        while (edge < xmax0) {
            double gradfactor = gradual > 0 ? 1. / gradual-- : 1;
            edge += binwidth0 * gradfactor;
            bins.push_back(edge);
        }
        
        his1 = TH1D("adaptHist", "adaptHist", bins.size() - 1, bins.data());
        his1.SetDirectory(nullptr);
        BruEventsHistPDF::FillBase1DHist(his1);
        for (auto i = 1; i <= his1.GetNbinsX(); i++) {
            his1.SetBinContent(i, his1.GetBinContent(i) / his1.GetBinWidth(i));
        }

        auto *ra = dynamic_cast<const RooRealVar*>(&alpha.arg());
        Construct2DHist(*his1.GetXaxis(), TAxis(_NAlphaBins, ra->getMin(), ra->getMax()), ra->isConstant());
    }

} // namespace bru


// /**
//  * @file BruEventsHistPeakPDF.cpp
//  */

// #include <Riostream.h> 
// #include "BruEventsHistPeakPDF.h" 
// #include <RooAbsReal.h> 
// #include <RooAbsCategory.h> 
// #include <cmath> 
// #include <TMath.h> 
// #include <TF1.h>
// #include <TH1.h>
// #include <TRandom3.h>
// #include <vector>
// #include <algorithm>
// #include <iostream>

// namespace bru {

//     BruEventsHistPeakPDF::BruEventsHistPeakPDF(const char *name, const char *title, RooAbsReal& in_x, RooAbsReal& in_alpha, RooAbsReal& in_offset, RooAbsReal& in_scale, Int_t applySmooth, Int_t interp, Int_t xbins, Int_t nsamp, Int_t abins) :
//         BruEventsHistPDF(name, title, in_x, in_alpha, in_offset, in_scale, applySmooth, interp, xbins, nsamp, abins)
//     {
//     }

//     BruEventsHistPeakPDF::BruEventsHistPeakPDF(const BruEventsHistPeakPDF& other, const char* name) :  
//         BruEventsHistPDF(other, name),
//         _lut(other._lut),
//         _varCenters(other._varCenters),
//         _nLut(other._nLut),
//         _lutWidth(other._lutWidth)
//     {
//     }

//     // =========================================================================
//     // LOOKUP TABLE INITIALIZATION (Dynamic Sizing + Center Search)
//     // =========================================================================
//     void BruEventsHistPeakPDF::initializeCache() {
//         BruEventsHistPDF::initializeCache();
        
//         if (!_Hist) return;
//         auto* vars = (RooArgSet*)_Hist->get();
//         RooRealVar* xVar = dynamic_cast<RooRealVar*>((*vars)[0]);
        
//         _varCenters.resize(_nx);
//         double minBinWidth = 1e9; 
        
//         for (int ix = 0; ix < _nx; ++ix) {
//             xVar->setBin(ix);
//             _varCenters[ix] = xVar->getVal();
            
//             double w = xVar->getBinWidth(ix);
//             if (w < minBinWidth) minBinWidth = w; 
//         }

//         double totalRange = _nx * _xWidth; 
//         _nLut = std::ceil(totalRange / (minBinWidth / 2.0));
//         if (_nLut > 500000) _nLut = 500000; 
//         if (_nLut < _nx * 2) _nLut = _nx * 2;
        
//         _lut.resize(_nLut);
//         _lutWidth = totalRange / _nLut; 
        
//         for (int i = 0; i < _nLut; ++i) {
//             double lut_x = _xMin + (i + 0.5) * _lutWidth;
            
//             auto it = std::upper_bound(_varCenters.begin(), _varCenters.end(), lut_x);
//             int left_center_idx = std::distance(_varCenters.begin(), it) - 1;
            
//             _lut[i] = std::max(0, std::min(left_center_idx, _nx - 2)); 
//         }
//     }

//     // =========================================================================
//     // LUT-ACCELERATED AVX2 EVALUATION (With Interpolation Switch)
//     // =========================================================================
//     void BruEventsHistPeakPDF::doEval(RooFit::EvalContext & ctx) const {
//         if (!_Hist || _cachedBins.empty()) {
//             std::fill(ctx.output().begin(), ctx.output().end(), 1.0); return;
//         }

//         auto xData = ctx.at(x); 
//         auto output = ctx.output();
        
//         const double cSc = scale, cOff = offset, vMax = _VarMax, cAl = alpha;
//         const double* bins = _cachedBins.data(); 
//         const int nx = _nx;
//         const int* lut = _lut.data();
//         const double* centers = _varCenters.data();

//         if (_Interpolate) {
//             // --- PATH A: SMOOTH BILINEAR INTERPOLATION ---
//             double cY = (cAl - _yMin) / _yWidth - 0.5;
//             double clY = std::max(0.0, std::min(cY, (double)(_ny - 1)));
//             int iy = std::min((int)clY, _ny - 2);
//             double yFr = clY - iy;
//             const double wY0 = 1.0 - yFr, wY1 = yFr;
//             const int r0 = iy * nx, r1 = (iy + 1) * nx;

//             for (size_t i = 0; i < xData.size(); ++i) {
//                 double arg = (xData[i] - vMax) * cSc + vMax - cOff;

//                 double lut_continuous = (arg - _xMin) / _lutWidth;
//                 int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
//                 int ix = lut[lut_idx];

//                 double c0 = centers[ix];
//                 double c1 = centers[ix + 1];
//                 double xFrac = (arg - c0) / (c1 - c0);
//                 xFrac = std::max(0.0, std::min(xFrac, 1.0));

//                 double b00 = bins[r0 + ix], b10 = bins[r0 + ix + 1];
//                 double b01 = bins[r1 + ix], b11 = bins[r1 + ix + 1];

//                 double vL = b00 * wY0 + b01 * wY1;
//                 double vR = b10 * wY0 + b11 * wY1;
//                 double res = vL + xFrac * (vR - vL);
                
//                 output[i] = res > 1e-15 ? res : 1e-15;
//             }
//         } else {
//             // --- PATH B: STEP FUNCTION (Nearest Bin, No Interpolation) ---
//             int iy = std::max(0, std::min((int)((cAl - _yMin) / _yWidth), _ny - 1));
//             const int r_nearest = iy * nx;

//             for (size_t i = 0; i < xData.size(); ++i) {
//                 double arg = (xData[i] - vMax) * cSc + vMax - cOff;
                
//                 double lut_continuous = (arg - _xMin) / _lutWidth;
//                 int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
//                 int ix = lut[lut_idx];

//                 double c0 = centers[ix];
//                 double c1 = centers[ix + 1];
                
//                 // Branchless math to find closest center
//                 int nearest_ix = (arg - c0) < (c1 - arg) ? ix : (ix + 1);

//                 double res = bins[r_nearest + nearest_ix];
//                 output[i] = res > 1e-15 ? res : 1e-15;
//             }
//         }
//     }

//     // =========================================================================
//     // MC BATCH INTEGRATION (With Interpolation Switch)
//     // =========================================================================
//     void BruEventsHistPeakPDF::evaluateMCBatch(const std::vector<Double_t>& mcx_array, std::vector<Double_t>& output) const {
//         size_t nEvents = mcx_array.size();
//         if (output.size() != nEvents) output.resize(nEvents);
        
//         if (!_Hist || _cachedBins.empty()) {
//             std::fill(output.begin(), output.end(), 1.0); return;
//         }

//         const double cSc = scale, cOff = offset, vMax = _VarMax, cAl = alpha;
//         const double* bins = _cachedBins.data(); 
//         const int nx = _nx;
//         const int* lut = _lut.data();
//         const double* centers = _varCenters.data();

//         if (_Interpolate) {
//             // --- PATH A: SMOOTH BILINEAR INTERPOLATION ---
//             double cY = (cAl - _yMin) / _yWidth - 0.5;
//             double clY = std::max(0.0, std::min(cY, (double)(_ny - 1)));
//             int iy = std::min((int)clY, _ny - 2);
//             double yFr = clY - iy;
//             const double wY0 = 1.0 - yFr, wY1 = yFr;
//             const int r0 = iy * nx, r1 = (iy + 1) * nx;

//             for (size_t i = 0; i < nEvents; ++i) {
//                 double arg = (mcx_array[i] - vMax) * cSc + vMax - cOff;

//                 double lut_continuous = (arg - _xMin) / _lutWidth;
//                 int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
//                 int ix = lut[lut_idx];

//                 double c0 = centers[ix];
//                 double c1 = centers[ix + 1];
//                 double xFrac = (arg - c0) / (c1 - c0);
//                 xFrac = std::max(0.0, std::min(xFrac, 1.0));

//                 double b00 = bins[r0 + ix], b10 = bins[r0 + ix + 1];
//                 double b01 = bins[r1 + ix], b11 = bins[r1 + ix + 1];

//                 double vL = b00 * wY0 + b01 * wY1;
//                 double vR = b10 * wY0 + b11 * wY1;
//                 double res = vL + xFrac * (vR - vL);
                
//                 output[i] = res > 1e-15 ? res : 1e-15;
//             }
//         } else {
//             // --- PATH B: STEP FUNCTION (Nearest Bin, No Interpolation) ---
//             int iy = std::max(0, std::min((int)((cAl - _yMin) / _yWidth), _ny - 1));
//             const int r_nearest = iy * nx;

//             for (size_t i = 0; i < nEvents; ++i) {
//                 double arg = (mcx_array[i] - vMax) * cSc + vMax - cOff;
                
//                 double lut_continuous = (arg - _xMin) / _lutWidth;
//                 int lut_idx = std::max(0, std::min((int)lut_continuous, _nLut - 1));
//                 int ix = lut[lut_idx];

//                 double c0 = centers[ix];
//                 double c1 = centers[ix + 1];
                
//                 // Branchless math to find closest center
//                 int nearest_ix = (arg - c0) < (c1 - arg) ? ix : (ix + 1);

//                 double res = bins[r_nearest + nearest_ix];
//                 output[i] = res > 1e-15 ? res : 1e-15;
//             }
//         }
//     }

//     // =========================================================================
//     // PEAK BINNING LOGIC (Unchanged)
//     // =========================================================================
//     void BruEventsHistPeakPDF::FillBase1DHist(TH1D& his1) {
//         std::cout << "BruEventsHistPeakPDF::FillBase1DHist " << his1.GetName() << std::endl;
        
//         BruEventsHistPDF::FillBase1DHist(his1);
        
//         int Nbins0 = his1.GetNbinsX();
//         double xmin0 = his1.GetBinLowEdge(0);
//         double xmax0 = his1.GetBinLowEdge(his1.GetNbinsX() + 1);
//         double xmin = xmin0;
//         double xmax = xmax0;
//         double range = xmax - xmin;
//         double binwidth0 = range / Nbins0;
//         double binwidth = binwidth0;
 
//         auto maxBin = his1.GetMaximumBin();
//         int deltaBin = his1.GetRMS() / binwidth; 
        
//         if (maxBin - deltaBin < 0) deltaBin = maxBin;
//         if (maxBin + deltaBin > Nbins0) deltaBin = Nbins0 - maxBin;
        
//         xmin = his1.GetBinLowEdge(maxBin - deltaBin);
//         xmax = his1.GetBinCenter(maxBin + deltaBin) + binwidth;
        
//         range = xmax - xmin;
//         binwidth = range / Nbins0;
        
//         TH1D peakhist("peakhist", "peakhist", Nbins0, xmin, xmax);
//         peakhist.SetDirectory(nullptr);
       
//         BruEventsHistPDF::FillBase1DHist(peakhist);

//         auto lastEdge = xmin;
//         std::vector<double> bins;
//         std::vector<double> fracs;
//         bins.push_back(lastEdge);
        
//         auto hmax = peakhist.GetMaximum();
//         double sumFracs = 0.;
//         double minFrac = 10;
        
//         for (auto i = 1; i <= peakhist.GetNbinsX(); i++) {
//             Double_t mainpeak = peakhist.GetBinContent(i);
//             Double_t leftpeak = (i == 1) ? 0 : peakhist.GetBinContent(i - 1);
//             Double_t rightpeak = (i == peakhist.GetNbinsX()) ? 0 : peakhist.GetBinContent(i + 1);
            
//             auto binFrac = (4 * 1.1 * hmax - 2 * mainpeak - leftpeak - rightpeak);
//             if (binFrac < 0) binFrac = 2;
            
//             fracs.push_back(binFrac);
//             if (binFrac < minFrac) minFrac = binFrac;
//         }
        
//         for (auto& frac : fracs) {
//             if (frac == 2) frac = minFrac; 
//             sumFracs += frac;
//         }
//         for (auto& frac : fracs) {
//             frac /= sumFracs;
//             lastEdge += frac * range;
//             bins.push_back(lastEdge);
//         }
     
//         double edge = bins[0];
//         Int_t gradual = 3;
//         while (edge > xmin0) {
//             double gradfactor = gradual > 0 ? 1. / gradual-- : 1;
//             edge -= binwidth0 * gradfactor;
//             bins.push_back(edge);
//             std::rotate(bins.rbegin(), bins.rbegin() + 1, bins.rend());
//         }
        
//         edge = bins.back();
//         gradual = 3;
//         while (edge < xmax0) {
//             double gradfactor = gradual > 0 ? 1. / gradual-- : 1;
//             edge += binwidth0 * gradfactor;
//             bins.push_back(edge);
//         }
        
//         his1 = TH1D("adaptHist", "adaptHist", bins.size() - 1, bins.data());
//         his1.SetDirectory(nullptr);
//         BruEventsHistPDF::FillBase1DHist(his1);
//         for (auto i = 1; i <= his1.GetNbinsX(); i++) {
//             his1.SetBinContent(i, his1.GetBinContent(i) / his1.GetBinWidth(i));
//         }

//         auto *ra = dynamic_cast<const RooRealVar*>(&alpha.arg());
//         Construct2DHist(*his1.GetXaxis(), TAxis(_NAlphaBins, ra->getMin(), ra->getMax()), ra->isConstant());
//     }

// } // namespace bru
