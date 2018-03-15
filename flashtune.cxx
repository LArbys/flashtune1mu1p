#include <iostream>
#include "Minuit2/Minuit2Minimizer.h"

#include "HandScanTable.h"
#include "FCNFlashMatch.h"

#include "TCanvas.h"
#include "TH1D.h"
#include "TText.h"
#include "TStyle.h"

int main( int nargs, char ** argv ) {

  gStyle->SetOptStat(0);
  
  bool printHists = true;
  std::string pngfolder="./png";
  bool fShapeOnly = true;
  
  // Minimizer
  ROOT::Minuit2::Minuit2Minimizer min ( ROOT::Minuit2::kMigrad );
  min.SetMaxFunctionCalls(1000000);
  min.SetMaxIterations(100000);
  min.SetTolerance(0.01);

  // Function
  //ftune::FCNFlashMatch llflash( "goodrecohandscan_1mu1p.tab", "flashdata.root" );
  //ftune::FCNFlashMatch llflash( "", "flashdata_extbnb.root" );
  ftune::FCNFlashMatch llflash( "", "flashdata_shower.root", ftune::FCNFlashMatch::kElectron );
  llflash.fitShapeOnly(fShapeOnly);

  double step = 0.01;
  double variable[33] = { 1.0 };
  variable[32] = 0.5;

  min.SetFunction( llflash );

  // Set the free variables to be minimized!
  for (int ipar=0; ipar<33; ipar++) {
    if (ipar<32) {
      // pmt parameters
      char parname[10];
      sprintf( parname, "opch%02d", ipar );
      min.SetLimitedVariable(ipar,parname,1.0,step, 0.5, 1.5);
    }
    else if ( ipar==32 ) {
      min.SetLimitedVariable(ipar,"global",variable[ipar],step, 0.0, 10.0);
    }
    
  }
  
  //min.Minimize();

  min.PrintResults();
  const double* pars = min.X();
  for (int i=0; i<33; i++) {
    std::cout << "[" << i << "] " << min.VariableName( i ) << " : " << pars[i] << std::endl;
  }

  llflash.EvalChi2( pars );

  int chi2pass = 0;
  int allchi2 = 0;
  for (int i=0; i<llflash.GetEntries(); i++) {
    float chi2 = llflash.GetEntryLastChi2(i);
    if ( chi2<3.0 )
      chi2pass++;
    allchi2++;
  }

  std::cout << "Vertices passing: " << chi2pass << " (" << float(chi2pass)/allchi2*100.0 << "%)" << std::endl;
  
  // print out hists ...
  if ( printHists ) {

    TCanvas c("c","c",800,600);

    int maxbin = 5000;
    if ( fShapeOnly )
      maxbin = 20;
    TH1D hchi2("hchi2", "", 100, 0, maxbin);
    TH1D hratio("hratio", "", 100, -5, 5 );
    
    for (int i=0; i<llflash.GetEntries(); i++) {
      std::vector<float> datape = llflash.GetEntryDataPE(i);
      std::vector<float> hypope = llflash.GetEntryHypoPE(i,fShapeOnly,pars);

      TH1D hdata("datape","",32,0,32);
      TH1D hhypo("hypope","",32,0,32);

      for (int ich=0; ich<32; ich++) {
	hdata.SetBinContent( ich+1, datape[ich] );
	hhypo.SetBinContent( ich+1, hypope[ich] );
      }
      float totdata = hdata.Integral();
      float totmc   = hhypo.Integral();

      float max = hdata.GetMaximum();
      if ( max < hhypo.GetMaximum() )
	max = hhypo.GetMaximum();

      hdata.SetMaximum( 1.1*max );

      c.cd();
      //hdata.Draw("E1");
      hdata.Draw("hist");
      hhypo.SetLineColor(kRed);
      hhypo.Draw("samehist");

      char zchi2[50];
      sprintf( zchi2, "chi2: %.2f", llflash.GetEntryLastChi2( i ) );
      TText tchi( 0.6, 0.8, zchi2 );
      tchi.SetNDC();
      tchi.Draw();
      hchi2.Fill( llflash.GetEntryLastChi2(i) );

      char zratio[50];
      sprintf( zratio, "mc/data: %.2f", totmc/totdata );
      TText tratio( 0.6, 0.7, zratio );
      tratio.SetNDC();      
      tratio.Draw();
      hratio.Fill( (totmc-totdata)/totdata );

      char zvertex[50];
      std::vector<float> vertex = llflash.GetEntryVertex(i);
      sprintf( zvertex, "pos: (%.1f,%.1f,%.1f)", vertex[0], vertex[1], vertex[2] );
      TText tvertex( 0.6, 0.6, zvertex );
      tvertex.SetNDC();
      tvertex.Draw();

      char zdir[50];
      std::vector<float> ldir = llflash.GetEntryDir(i);
      sprintf( zdir, "dir: (%.1f,%.1f,%.1f)", ldir[0], ldir[1], ldir[2] );
      TText tdir( 0.6, 0.5, zdir );
      tdir.SetNDC();      
      tdir.Draw();

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
