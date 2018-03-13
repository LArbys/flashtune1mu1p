#include "FCNFlashMatch.h"

#include <iostream>

#include "TFile.h"
#include "TTree.h"

#include "HandScanTable.h"

namespace ftune {

  FCNFlashMatch::FCNFlashMatch( std::string handscantable, std::string flashdata ) {

    // load hand-scan table
    m_goodreco_info = new HandScanTable( handscantable );

    fNEntries = m_goodreco_info->numEntries();
    
    // data tables
    m_data_array.resize(32*fNEntries,0);

    // hypo thesis tables
    m_hypo_array.resize(32*fNEntries,0);
    
    // load flash data
    int run;
    int subrun;
    int event;
    int vtxid;

    std::vector<float>* pe_data = 0;
    std::vector<float>* pe_hypo = 0;

    TFile froot( flashdata.c_str() );
    TTree* tree = (TTree*)froot.Get("ffmatch");
    tree->SetBranchAddress( "run", &run );
    tree->SetBranchAddress( "subrun", &subrun );
    tree->SetBranchAddress( "event", &event );
    tree->SetBranchAddress( "vertexid", &vtxid );
    tree->SetBranchAddress( "data_pe", pe_data );
    tree->SetBranchAddress( "hypo_pe_1mu1p", pe_hypo );

    unsigned long ientry = 0;
    unsigned long bytes = tree->GetEntry(ientry);

    int idata_entry = 0;
    while ( bytes>0 ) {

      // we using this vertex?
      int entry_goodvtx = m_goodreco_info->GetVertexID( run, subrun, event );
      if ( entry_goodvtx==vtxid ) {
	// a good vertex. copy into data arrays
	for (int ich=0; ich<32; ich++) {
	  m_data_array[ 32*idata_entry + ich ] = (*pe_data)[ich];
	  m_hypo_array[ 32*idata_entry + ich ] = (*pe_hypo)[ich];
	}
	idata_entry++;
      }

      if ( idata_entry==fNEntries ) {
	std::cout << "we've filled all the allotted data positions (" << fNEntries << ")" << std::endl;
	break;
      }
      
      ientry++;
      tree->GetEntry(ientry);
    }

    // don't need this anymore
    delete m_goodreco_info;
  }

  FCNFlashMatch::~FCNFlashMatch() {
    m_data_array.clear();
    m_hypo_array.clear();
  }
  
  double FCNFlashMatch::operator() (const std::vector<double>& x) const {

    float globalfactor = x[0];

    float chi2 = 0;
    for (int ientry=0; ientry<fNEntries; ientry++) {

      for (int ich=0; ich<32; ich++) {
	float obs  = m_data_array[ ientry*32 + ich ];
	float pred = m_hypo_array[ ientry*32 + ich ]*globalfactor*x[ich+1]; // we allow for global correction and per-pmt correction
	if ( pred<1.0e-3)
	  pred = 1.0e-3;
	float diff = obs-pred;
	chi2 += diff*diff/pred;
      }
    }

    float ll = -0.5*chi2;

    return ll;
  };


}
