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

    std::cout << "Number of Good Reco Entries: " << fNEntries << std::endl;
    
    // data tables
    m_data_array.resize(32*fNEntries,0);

    // hypo thesis tables
    m_hypo_array.resize(32*fNEntries,0);
    
    // load flash data
    int run;
    int subrun;
    int event;
    int vtxid;

    // std::vector<float>* pe_data = 0;
    // std::vector<float>* pe_hypo = 0;
    // TBranch* bpe_data = 0;
    // TBranch* bpe_hypo = 0;
    float pe_data[32];
    float pe_hypo[32];

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

      // std::cout << "entry " << ientry << ": "
      // 		<< "(rse)=(" << run << "," << subrun << "," << event << ") "
      // 		<< "vtxid=" << vtxid << " "
      // 		<< "goodvtx=" << entry_goodvtx
      // 		<< std::endl;
      
      if ( entry_goodvtx==vtxid ) {
	// a good vertex. copy into data arrays
	for (int ich=0; ich<32; ich++) {
	  m_data_array[ 32*idata_entry + ich ] = pe_data[ich];
	  m_hypo_array[ 32*idata_entry + ich ] = pe_hypo[ich];
	}
	std::vector<int> rsev(4);
	rsev[0] = run;
	rsev[1] = subrun;
	rsev[2] = event;
	rsev[3] = vtxid;
	m_rse.push_back( rsev );
	idata_entry++;
      }
      
      if ( idata_entry==fNEntries ) {
	std::cout << "we've filled all the allotted data positions (" << fNEntries << ")" << std::endl;
	break;
      }
      
      ientry++;
      bytes = tree->GetEntry(ientry);
    }

    std::cout << "Filled " << idata_entry << " entries"  << std::endl;
    fNEntries = idata_entry;
    m_last_chi2.resize(fNEntries,0.0);

    // don't need this anymore
    delete m_goodreco_info;
  }

  FCNFlashMatch::~FCNFlashMatch() {
    m_data_array.clear();
    m_hypo_array.clear();
  }
  
  //double FCNFlashMatch::operator() (const std::vector<double>& x) const {
  double FCNFlashMatch::DoEval( const double* x ) const {
    
    float globalfactor = x[0];

    // std::cout << "pars [ ";
    // for (int i=0; i<33; i++)
    //   std::cout << x[i] << " ";
    // std::cout << "]" << std::endl;
    
    float chi2 = 0;
    for (int ientry=0; ientry<fNEntries; ientry++) {
      float chi2_entry = 0.;
      for (int ich=0; ich<32; ich++) {
	float obs  = m_data_array[ ientry*32 + ich ];
	float pred = m_hypo_array[ ientry*32 + ich ]*globalfactor*x[ich+1]; // we allow for global correction and per-pmt correction
	if ( pred<1.0e-3)
	  pred = 1.0e-3;
	float diff = obs-pred;
	float err = sqrt( obs*obs + pred*pred );
	chi2_entry += diff*diff/err;
      }
      //std::cout << "  " << chi2_entry << std::endl;
      chi2 += chi2_entry;
    }

    chi2 /= (fNEntries*32-33); // number of bins we're fitting - dof
    
    std::cout << "total chi2=" << chi2 << std::endl;
    //std::cin.get();
    return chi2;
  };

  void FCNFlashMatch::EvalChi2( const double* x ) {
    
    float globalfactor = x[0];

    // std::cout << "pars [ ";
    // for (int i=0; i<33; i++)
    //   std::cout << x[i] << " ";
    // std::cout << "]" << std::endl;
    
    for (int ientry=0; ientry<fNEntries; ientry++) {
      float chi2_entry = 0.;
      for (int ich=0; ich<32; ich++) {
	float obs  = m_data_array[ ientry*32 + ich ];
	float pred = m_hypo_array[ ientry*32 + ich ]*globalfactor*x[ich+1]; // we allow for global correction and per-pmt correction
	if ( pred<1.0e-3)
	  pred = 1.0e-3;
	float diff = obs-pred;
	float err = sqrt( obs*obs + pred*pred );
	chi2_entry += diff*diff/err;
      }
      //std::cout << "  " << chi2_entry << std::endl;
      m_last_chi2[ientry] = chi2_entry;
    }

  };


  // copy constructor
  ROOT::Math::IBaseFunctionMultiDim * FCNFlashMatch::Clone() const {

    FCNFlashMatch* clone = new FCNFlashMatch;
    clone->fNEntries    = this->fNEntries;
    clone->m_data_array = this->m_data_array;
    clone->m_hypo_array = this->m_hypo_array;
    
  }

  std::vector<float> FCNFlashMatch::GetEntryDataPE( int ientry ) {
    std::vector<float> pe(32);
    for ( int i=0; i<32; i++) {
      pe[i] = m_data_array[ 32*ientry + i ];
    }
    return pe;
  }
  
  std::vector<float> FCNFlashMatch::GetEntryHypoPE( int ientry ) {
    std::vector<float> pe(32);
    for ( int i=0; i<32; i++) {
      pe[i] = m_hypo_array[ 32*ientry + i ];
    }
    return pe;
  }

  void FCNFlashMatch::GetEntryRSEV( int ientry, int& run, int& subrun, int& event, int& vtxid ) {
    const std::vector<int>& rsev = m_rse[ientry];
    run = rsev[0];
    subrun = rsev[1];
    event = rsev[2];
    vtxid = rsev[3];
  }
  
}
