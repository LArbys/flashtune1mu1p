#include "HandScanTable.h"
#include <iostream>
#include <cstring>
#include <fstream>

namespace ftune {

  HandScanTable::HandScanTable( std::string tabfile )
    : fTabFile(tabfile) {

    // load the map
    std::ifstream f(fTabFile.c_str());

    std::string line;
    while ( !f.eof() && f.good() ) {
      std::getline( f, line );
      if (line.empty())
	break;
      int run;
      int subrun;
      int event;
      int vtxid;
      float ll;
      int goodreco;
      int goodscan;
      char leftover[200];
      sscanf( line.c_str(), "%d\t%d\t%d\t%d\t%f\t%d\t%d\t%s",&run,&subrun,&event,&vtxid,&ll,&goodreco,&goodscan,leftover);

      // for debug
      //std::cout << run << " " << subrun << " " << event << " " << vtxid << " goodscan=" << goodscan << std::endl;

      if ( goodscan==1 ) {
	RSE rse(run,subrun,event);
	m_map_rse2vertexid.insert( std::make_pair<ftune::HandScanTable::RSE,int>(std::move(rse),std::move(vtxid)) );
      }
      
    }

  }

  int HandScanTable::GetVertexID( int r, int s, int e ) {
    RSE rse( r, s, e );

    auto it = m_map_rse2vertexid.find(rse);
    if ( it!=m_map_rse2vertexid.end() )
      return it->second;
    return -1;
  }

}
