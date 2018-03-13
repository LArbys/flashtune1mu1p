#include <iostream>
#include "Minuit2/Minuit2Minimizer.h"

#include "HandScanTable.h"
#include "FCNFlashMatch.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TText.h"


int main( int nargs, char ** argv ) {

  bool printHists = true;
  std::string pngfolder="./png";
  
  // Minimizer
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.01);

  // Function
  ftune::FCNFlashMatch llflash( "goodrecohandscan_1mu1p.tab", "flashdata.root" );

  double step = 0.01;
  double variable[33] = { 1.0 };
  for (int i=0; i<33; i++)
    variable[i] = 1.0;

  min.SetFunction( llflash );

  // Set the free variables to be minimized!
  for (int ipar=0; ipar<33; ipar++) {
    if ( ipar==0 )
      min.SetLimitedVariable(ipar,"global", variable[ipar], step, 0.5, 1.5 );
    else {
      char parname[10];
      sprintf( parname, "opch%02d", ipar-1 );
      min.SetLimitedVariable(ipar,parname,1.0,step, 0.5, 1.5);
    }
  }
  
  min.Minimize();

  min.PrintResults();
  const double* pars = min.X();
  for (int i=0; i<33; i++) {
    std::cout << "[" << i << "] " << min.VariableName( i ) << " : " << pars[i] << std::endl;
  }

  llflash.EvalChi2( pars );
  
  // print out hists ...
  if ( printHists ) {

    TCanvas c("c","c",800,600);

    TH1D hchi2("hchi2", "", 500, 0, 100);
    TH1D hratio("hratio", "", 100, -5, 5 );
    
    for (int i=0; i<llflash.GetEntries(); i++) {
      std::vector<float> datape = llflash.GetEntryDataPE(i);
      std::vector<float> hypope = llflash.GetEntryHypoPE(i);

      TH1D hdata("datape","",32,0,32);
      TH1D hhypo("hypope","",32,0,32);

      for (int ich=0; ich<32; ich++) {
	hdata.SetBinContent( ich+1, datape[ich] );
	hhypo.SetBinContent( ich+1, pars[0]*pars[ich+1]*hypope[ich] );
      }
      float totdata = hdata.Integral();
      float totmc   = hhypo.Integral();

      float max = hdata.GetMaximum();
      if ( max < hhypo.GetMaximum() )
	max = hhypo.GetMaximum();

      hdata.SetMaximum( 1.1*max );

      c.cd();
      hdata.Draw("E1");
      hhypo.SetLineColor(kRed);
      hhypo.Draw("samehist");

      char zchi2[50];
      sprintf( zchi2, "chi^{2}: %.2f", llflash.GetEntryLastChi2( i ) );
      TText tchi( 20, 12, zchi2 );
      tchi.Draw();
      hchi2.Fill( llflash.GetEntryLastChi2(i)/32.0 );

      char zratio[50];
      sprintf( zratio, "mc/data: %.2f", totmc/totdata );
      TText tratio( 20, 8, zratio );
      tratio.Draw();
      hratio.Fill( (totmc-totdata)/totdata );

      c.Update();

      int run, subrun, event, vtxid;
      llflash.GetEntryRSEV( i, run, subrun, event, vtxid );
      char zcanvas[200];
      sprintf( zcanvas, "png/pe_entry%d_run%d_subrun%d_event%d_vtx%d.png", i, run, subrun, event, vtxid );
      c.SaveAs( zcanvas );

      c.cd();
      hchi2.Draw();
      c.SaveAs( "png/chi2.png");

      c.cd();
      hratio.Draw();
      c.SaveAs( "png/mc_data_ratio.png" );
      
    }//end of entry loop
    
  }//if print Hists

  
  return 0;
}
