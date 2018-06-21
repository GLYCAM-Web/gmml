#ifndef AMBER_ELEMENTS_META_HPP
#define AMBER_ELEMENTS_META_HPP

/* File amberelements.hpp begun on 16 June 2018 by BLFoley */

#include <string>
#include <map>
#include <set>


namespace gmml::MolecularMetadata::AMBER
{

  typedef struct 
  {
    std::string element_;  // The element symbol
    std::string mass_;     // The mass, in amu, used by most/all of Amber
  } AmberAtomInfo ;
  
  const AmberAtomInfo AMBERATOMINFO[] =
  {
  //{ "Element" , "Mass" },
    { "C"       , 12.011 },
    { "H"       , 1.008  },
    { "N"       , 14.01  },
    { "O"       , 16.00  },
    { "S"       , 32.06  },
    { "P"       , 30.97  }
  };
}
#endif // AMBER_ELEMENTS_META_HPP
