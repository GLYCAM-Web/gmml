#include "../../../includes/MolecularMetadata/AMBER/amberelements.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::AMBER::AmberElementContainer;

AmberElementContainer::AmberElementContainer()
{            // const AmberAtomInfo AMBERATOMINFO[] =
    amberElementVector_ =
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
