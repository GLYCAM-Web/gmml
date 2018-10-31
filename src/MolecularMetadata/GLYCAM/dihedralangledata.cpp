#include "../../../../includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer;

//struct DihedralAngleData
//{
//    std::string linking_atom1_ ;
//    std::string linking_atom2_ ;
//    std::string dihedral_angle_name_ ;
//    double default_angle_value_ ;
//    double lower_deviation_ ;
//    double upper_deviation_ ;
//    std::string name_ ;
//    int index_ ; // if two entries match the criteria, and have the same index, the later entry should overwrite the earlier.
//    std::string residue1_condition_ ;
//    std::string residue2_condition_ ;
//    std::string atom4_ ; // I always want rotation to be in the direction of protein->glycan and within glycan from reducing terminal to non-reducing.
//    std::string atom3_ ;
//    std::string atom2_ ;
//    std::string atom1_ ;
//} ;

// A note on index_. When two entries match the criteria, sometimes you want the second entry to override the first, in that case give them the same index_ number.
// Further, if two entries are mutually exclusive, but apply to the same dihedral angle (e.g. phi), give them the same index_ number.
// If multiple entries both apply to phi, but have different index_ numbers, then two rotamers exist

DihedralAngleDataContainer::DihedralAngleDataContainer()
{   // const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
    dihedralAngleDataVector_ =
    {
        { "C1", ".."     , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  , 1  , "none"   , "none"                  ,  "C2" , "C1" , "O." , "C."  }, // phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
        { "C2", ".."     , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  , 1  , "none"   , "none"                  ,  "C1" , "C2" , "O." , "C."  }, // phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
        { "C2", ".."     , "phi"  , -60.0  ,  20.0  ,  20.0  , "-g" , 2  , "ketose" , "none"                  ,  "C1" , "C2" , "O." , "C."  },

        { "C.", "O[1-5]" , "psi"  ,   0.0  ,  20.0  ,  20.0  , "t"  , 3  , "none"   , "none"                  ,  "C." , "O." , "C." , "H."  }, // psi should be C(ano)-Ox-Cx-Hx, if Cx is ring, otherwise, C(ano)-Ox-Cx-C(x-1)
        { "C.", "O[6-9]" , "psi"  , 180.0  ,  20.0  ,  20.0  , "t"  , 3  , "none"   , "none"                  ,  "C." , "O." , "C." , "C."  },

        { "C.", "O6"     , "omg"  , -60.0  ,  20.0  ,  20.0  , "gg" , 4  , "none"   , "none"                  ,  "O6" , "C6" , "C5" , "O5"  }, // omg is O6-C5-C5-O5
        { "C.", "O6"     , "omg"  ,  60.0  ,  20.0  ,  20.0  , "gt" , 5  , "none"   , "none"                  ,  "O6" , "C6" , "C5" , "O5"  },
        { "C.", "O6"     , "omg"  , 180.0  ,  20.0  ,  20.0  , "tg" , 6  , "none"   , "gauche-effect=galacto" ,  "O6" , "C6" , "C5" , "O5"  }, // need to add this tag
        { "C.", ".."     , "chi1" , -60.0  ,  20.0  ,  20.0  , ""   , 7  , "none"   , "protein"               ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", ".."     , "chi1" ,  60.0  ,  20.0  ,  20.0  , ""   , 8  , "none"   , "protein"               ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", ".."     , "chi1" , 180.0  ,  20.0  ,  20.0  , ""   , 9  , "none"   , "protein"               ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "ND2"    , "chi2" , -60.0  ,  60.0  ,  60.0  , ""   , 10 , "none"   , "protein"               ,  "ND2", "CG" , "CB" , "CA"  }, // Asn
        { "C.", "OG"     , "chi2" , -60.0  ,  60.0  ,  60.0  , ""   , 10 , "none"   , "protein"               ,  "OG" , "CG" , "CB" , "CA"  }, // Ser/Thr
        { "C.", "OH"     , "chi2" , -60.0  ,  60.0  ,  60.0  , ""   , 10 , "none"   , "protein"               ,  "OH" , "CG" , "CB" , "CA"  }, // Tyr
        //{ "C2", "O8"   , "omg7" ,
    };
}
