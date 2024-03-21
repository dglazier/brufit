
#include <Riostream.h> 

#include "BruEventsHistPeakPDF.h" 
#include "RooAbsReal.h" 
#include "RooAbsCategory.h" 
#include <cmath> 
#include <TMath.h> 
#include <TF1.h>
#include <TH1.h>
#include <TRandom3.h>

namespace HS{
  namespace FIT{

    BruEventsHistPeakPDF::BruEventsHistPeakPDF(const char *name, const char *title, RooAbsReal& _x,RooAbsReal& _alpha,RooAbsReal& _offset, RooAbsReal& _scale, Int_t applySmooth, Int_t interp, Int_t xbins, Int_t nsamp, Int_t abins) :
      RooHSEventsHistPDF(name,title,_x,_alpha,_offset,_scale, applySmooth, interp,xbins,nsamp,abins)
    {
  
     
    
    }

    BruEventsHistPeakPDF::BruEventsHistPeakPDF(const BruEventsHistPeakPDF& other, const char* name) :  
      RooHSEventsHistPDF(other,name)
    {
    }


    void  BruEventsHistPeakPDF::FillBase1DHist(TH1D& his1){
      cout<<"BruEventsHistPeakPDF::FillBase1DHist "<<his1.GetName()<<endl;
      //get standard 1D histogram
      RooHSEventsHistPDF::FillBase1DHist(his1);
      //now adapt the binning around the peak
     int Nbins0 = his1.GetNbinsX();
     double xmin0 = his1.GetBinLowEdge(0);
     double xmax0 = his1.GetBinLowEdge(his1.GetNbinsX()+1);
     double xmin=xmin0;
     double xmax=xmax0;
     double range = xmax-xmin;
     double binwidth0 = range/Nbins0;
     double binwidth = binwidth0;
 
     //select symmetric region around peak to use adpative binning
     auto maxBin = his1.GetMaximumBin();
     //need to center peak to get binning scheme correct
     int deltaBin = his1.GetRMS()/binwidth; //range taken as +- 1 sigma
     //check for boundaries
     if(maxBin - deltaBin < 0){
       deltaBin = maxBin;
     }
     if(maxBin + deltaBin > Nbins0){
       deltaBin = Nbins0-maxBin;
     }
     xmin = his1.GetBinLowEdge(maxBin - deltaBin);
     xmax = his1.GetBinCenter(maxBin+deltaBin)+binwidth;
     cout<<maxBin<<" "<<his1.GetBinCenter(maxBin)<<" "<<xmin<<" "<<xmax<<endl;;
     //Calculate new range and bin width for around the peak
     range = xmax-xmin;
     binwidth = range/Nbins0;
     //remake 1D histogram with new binning
     TH1D peakhist("peakhist","peakhist",Nbins0,xmin,xmax);
     RooHSEventsHistPDF::FillBase1DHist(peakhist);

     //create variable bins from xmin to xmax
     //these are depedent on the bin content
     //smaller bins for larger yield
     auto lastEdge = xmin;
     vector<double> bins;
     vector<double> fracs;
     bins.push_back(lastEdge);
     //get peak yield 
     auto hmax = peakhist.GetMaximum();
     double sumFracs =0.;
     double minFrac = 10;
     //loopp over old bins and create new bins array
     for(auto i = 1 ; i <= peakhist.GetNbinsX(); i++){
       //This linescales the bin widths dependent on their content
       //We also include the nearest neighbours weighted by 0.5
       //This stops the smallest bin being too small
       //similar the factor 1.01 stops bincontent = hmax
       auto binFrac = (4*1.01*hmax-2*peakhist.GetBinContent(i)-peakhist.GetBinContent(i-1)-peakhist.GetBinContent(i-1));
       //protect from nonsense
       if(binFrac<0)binFrac=2;
       //      cout<<binFrac<<" "<<hist->GetBinContent(i)<<" "<<hist->GetBinCenter(i)<<endl;
       fracs.push_back(binFrac);
       if(binFrac<minFrac) minFrac = binFrac;
     }
     //now we need to rescale all of the bin widths to cover the full range
     for(auto& frac:fracs){//correct in case of zero bin contents
       if(frac==2) frac = minFrac;
       sumFracs+=frac;
     }
     for(auto& frac:fracs){
       frac/=sumFracs;
       lastEdge+=frac*range;
       bins.push_back(lastEdge);
     }
     //now need to fill out rest with inital binning scheme.
     double edge = bins[0];
     while(edge>xmin0){
       edge-=binwidth0;
       //this inserts element at start
       bins.push_back(edge);
       std::rotate(bins.rbegin(), bins.rbegin() + 1, bins.rend());
     }
     edge = bins.back();
     cout<<"edge "<<edge<<" "<<xmax0<<endl;
     while(edge<xmax0){
       edge+=binwidth0*2;
       //this inserts element at back
       bins.push_back(edge);
     }
     //replace original histogram
     his1 = TH1D("adaptHist","adaptHist",bins.size()-1,bins.data());
     RooHSEventsHistPDF::FillBase1DHist(his1);
     for(auto i = 1 ; i <= his1.GetNbinsX(); i++){
       his1.SetBinContent(i,his1.GetBinContent(i)/his1.GetBinWidth(i));
     }

     //reconstruct 2D hist scheme
     auto *ra=dynamic_cast<const RooRealVar*>(&alpha.arg());
     cout<<"BruEventsHistPeakPDF::FillBase1DHist "<<ra<<" "<<his1.GetNbinsX()<<endl;
     Construct2DHist(*his1.GetXaxis(),TAxis(fNAlphaBins,ra->getMin(),ra->getMax()),ra->isConstant());
     //all done
    }
  }
}
