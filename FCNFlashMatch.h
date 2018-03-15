#ifndef __FCNFLASHMATCH_H__
#define __FCNFLASHMATCH_H__

//#include "Minuit2/FCNBase.h"
#include "Math/Minimizer.h"
#include <vector>
#include <string>

namespace ftune {

  class HandScanTable;
  
  //class FCNFlashMatch : public ROOT::Minuit2::FCNBase {
  class FCNFlashMatch : public ROOT::Math::IMultiGenFunction {

  public:

    typedef enum { kMuon, kElectron } HypoType; 
    
    FCNFlashMatch() {};
    FCNFlashMatch( std::string handscantable, std::string flashdata, FCNFlashMatch::HypoType hypotype=kMuon );
    virtual ~FCNFlashMatch();

    //virtual double operator() (const std::vector<double>& x) const;
    //virtual double Up() const { return 1.0; };
    
    virtual double DoEval( const double* x) const ;
    virtual unsigned int NDim() const {
      if ( fShapeOnly ) return 33;
      else return 33;
    };
    virtual ROOT::Math::IBaseFunctionMultiDim * Clone() const;
    void EvalChi2( const double* x );
    void fitShapeOnly( bool shapeonly ) { fShapeOnly = shapeonly; };
    
  protected:

    HandScanTable* m_goodreco_info;
    int fNEntries;
    bool fShapeOnly;
    std::vector<float> m_data_array;
    std::vector<float> m_hypo_array;
    std::vector< std::vector<int> > m_rse;
    std::vector< std::vector<float> > m_pos;
    std::vector< std::vector<float> > m_dir;    
    std::vector< float > m_last_chi2;

  public:
    
    int GetEntries() const { return fNEntries; };
    std::vector<float> GetEntryDataPE( int ientry );
    std::vector<float> GetEntryHypoPE( int ientry, bool shapeonly=false, const double *pars=NULL );
    void GetEntryRSEV( int ientry, int& run, int& subrun, int& event, int& vtxid );
    float GetEntryLastChi2( int ientry ) { return m_last_chi2[ientry]; };
    std::vector<float> GetEntryVertex( int ientry ) { return m_pos[ientry]; };
    std::vector<float> GetEntryDir( int ientry ) { return m_dir[ientry]; };    

  };

}


#endif
