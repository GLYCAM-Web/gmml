#include "../../../../includes/MolecularMetadata/GLYCAM/dihedralangledata.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::DihedralAngleDataContainer;

//struct DihedralAngleData
//{
//std::string linking_atom1_ ;
//std::string linking_atom2_ ;
//std::string dihedral_angle_name_ ;
//double default_angle_value_ ;
//double lower_deviation_ ;
//double upper_deviation_ ;
//std::string rotamer_name_ ;
//double index_ ; // if two entries match the criteria, and have the same index, the later entry should overwrite the earlier.
//std::string residue1_condition_ ;
//std::string residue2_condition_ ;
//std::string atom1_ ;
//std::string atom2_ ;
//std::string atom3_ ;
//std::string atom4_ ;
//} ;

// A note on index_. When two entries match the criteria, sometimes you want the second entry to override the first, in that case give them the same index_ number.
// Further, if two entries are mutually exclusive, but apply to the same dihedral angle (e.g. phi), give them the same index_ number.
// If multiple entries both apply to phi, but have different index_ numbers, then two rotamers exist

DihedralAngleDataContainer::DihedralAngleDataContainer()
{   // const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
    dihedralAngleDataVector_ =
    {    // Regex:       , name                                                                                 // Atom names this applies to
        { "C1", "O."     , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  1.1 , "none"   , "none"                  ,  "C2" , "C1" , "O." , "C."  }, // phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
        { "C2", "O."     , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  1.1 , "none"   , "none"                  ,  "C1" , "C2" , "O." , "C."  }, // phi should be C2-C1(ano)-Ox-Cx, or C1-C2(ano)-Ox-Cx
        { "C2", "O."     , "phi"  , -60.0  ,  20.0  ,  20.0  , "-g" ,  1.2 , "ulosonate"  , "none"                  ,  "C1" , "C2" , "O." , "C."  },

        { "C.", "O[1-5]" , "psi"  ,   0.0  ,  20.0  ,  20.0  , "t"  ,  2.1 , "none"   , "none"                  ,  "C." , "O." , "C." , "H."  }, // psi should be C(ano)-Ox-Cx-Hx, if Cx is ring, otherwise, C(ano)-Ox-Cx-C(x-1)
        { "C.", "O[6-9]" , "psi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  2.1 , "none"   , "none"                  ,  "C." , "O." , "C." , "C."  },

        { "C.", "O6"     , "omg"  , -60.0  ,  20.0  ,  20.0  , "gg" ,  3.1 , "none"   , "none"                  ,  "O6" , "C6" , "C5" , "O5"  }, // omg is O6-C5-C5-O5
        { "C.", "O6"     , "omg"  ,  60.0  ,  20.0  ,  20.0  , "gt" ,  3.2 , "none"   , "none"                  ,  "O6" , "C6" , "C5" , "O5"  },
        { "C.", "O6"     , "omg"  , 180.0  ,  20.0  ,  20.0  , "tg" ,  3.3 , "none"   , "gauche-effect=galacto" ,  "O6" , "C6" , "C5" , "O5"  }, // need to add this tag

         // 2-8 linkages
      //{ "C2", "O8"   , "omg7" ,    };

        // Protein linkages
        // ASN // Values not checked
        { "C.", "ND2"    , "chi1" , -60.0  ,  20.0  ,  20.0  , ""   ,  3.1 , "none"   , "amino-acid"            ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "ND2"    , "chi1" ,  60.0  ,  20.0  ,  20.0  , ""   ,  3.2 , "none"   , "amino-acid"            ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "ND2"    , "chi1" , 180.0  ,  20.0  ,  20.0  , ""   ,  3.3 , "none"   , "amino-acid"            ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "ND2"    , "chi2" , -60.0  ,  60.0  ,  60.0  , ""   ,  4.1 , "none"   , "amino-acid"            ,  "ND2", "CG" , "CB" , "CA"  },
        { "C.", "ND2"    , "psi"  , -60.0  ,  60.0  ,  60.0  , ""   ,  2.1 , "none"   , "amino-acid"            ,  "C." , "ND2", "CG" , "CB"  },
        { "C1", "ND2"    , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  1.1 , "none"   , "amino-acid"            ,  "C." , "C." , "ND2", "CG"  },

        // THR // Values not checked
        { "C.", "OG1"    , "chi1" , -60.0  ,  20.0  ,  20.0  , ""   ,  3.1 , "none"   , "amino-acid"            ,  "OG1", "CB" , "CA" , "N"   },
        { "C.", "OG1"    , "chi1" ,  60.0  ,  20.0  ,  20.0  , ""   ,  3.2 , "none"   , "amino-acid"            ,  "OG1", "CB" , "CA" , "N"   },
        { "C.", "OG1"    , "chi1" , 180.0  ,  20.0  ,  20.0  , ""   ,  3.3 , "none"   , "amino-acid"            ,  "OG1", "CB" , "CA" , "N"   },
        { "C.", "OG1"    , "psi"  , -60.0  ,  60.0  ,  60.0  , ""   ,  2.1 , "none"   , "amino-acid"            ,  "C." , "OG1", "CB" , "CA"  },
        { "C.", "OG1"    , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  1.1 , "none"   , "amino-acid"            ,  "C." , "C." , "OG1", "CB"  },

         // SER // Values not checked
        { "C.", "OG"     , "chi1" , -60.0  ,  20.0  ,  20.0  , ""   ,  3.1 , "none"   , "amino-acid"            ,  "OG" , "CB" , "CA" , "N"   },
        { "C.", "OG"     , "chi1" ,  60.0  ,  20.0  ,  20.0  , ""   ,  3.2 , "none"   , "amino-acid"            ,  "OG" , "CB" , "CA" , "N"   },
        { "C.", "OG"     , "chi1" , 180.0  ,  20.0  ,  20.0  , ""   ,  3.3 , "none"   , "amino-acid"            ,  "OG" , "CB" , "CA" , "N"   },
        { "C.", "OG"     , "psi"  , -60.0  ,  60.0  ,  60.0  , ""   ,  2.1 , "none"   , "amino-acid"            ,  "C." , "OG" , "CB" , "CA"  },
        { "C.", "OG"     , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  1.1 , "none"   , "amino-acid"            ,  "C." , "C." , "OG1", "CB"  },

         // TYR // Values not checked
        { "C.", "OH"     , "chi1" , -60.0  ,  20.0  ,  20.0  , ""   ,  6.1 , "none"   , "amino-acid"            ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "OH"     , "chi1" ,  60.0  ,  20.0  ,  20.0  , ""   ,  6.2 , "none"   , "amino-acid"            ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "OH"     , "chi1" , 180.0  ,  20.0  ,  20.0  , ""   ,  6.3 , "none"   , "amino-acid"            ,  "CG" , "CB" , "CA" , "N"   },
        { "C.", "OH"     , "chi2" , -60.0  ,  60.0  ,  60.0  , ""   ,  7.1 , "none"   , "amino-acid"            ,  "CD1", "CG" , "CB" , "CA"  },
        { "C.", "OH"     , "psi"  , -60.0  ,  60.0  ,  60.0  , ""   ,  2.1 , "none"   , "amino-acid"            ,  "C." , "OH" , "CZ" , "CE1" },
        { "C.", "OH"     , "phi"  , 180.0  ,  20.0  ,  20.0  , "t"  ,  1.1 , "none"   , "amino-acid"            ,  "C." , "C." , "OH ", "CZ"  },

    };
}

//    Statistical analysis of the protein environment of N-glycosylation sites: implications for occupancy, structure, and folding
//    Andrei-J. Petrescu  Adina-L. Milac  Stefana M. Petrescu  Raymond A. Dwek Mark R. Wormald
//    Glycobiology, Volume 14, Issue 2, 1 February 2004, Pages 103–114,
