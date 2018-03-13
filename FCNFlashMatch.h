#ifndef __FCNFLASHMATCH_H__
#define __FCNFLASHMATCH_H__

#include "Minuit2/FCNBase.h"
#include <vector>
#include <string>

namespace ftune {

  class HandScanTable;
  
  class FCNFlashMatch : public ROOT::Minuit2::FCNBase {

  public:

    FCNFlashMatch() {};
    FCNFlashMatch( std::string handscantable, std::string flashdata );
    virtual ~FCNFlashMatch();

    virtual double operator() (const std::vector<double>& x) const;

  protected:

    HandScanTable* m_goodreco_info;
    int fNEntries;
    std::vector<float> m_data_array;
    std::vector<float> m_hypo_array;
    
  };

}


#endif
