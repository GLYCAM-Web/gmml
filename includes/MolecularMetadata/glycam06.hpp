#ifndef GLYCAM06_META_HPP
#define GLYCAM06_META_HPP

#include <string>
#include <vector>
#include <set>

namespace gmml::MolecularMetadata::GLYCAM
{
  // Move this to an amber space once the proper design is apparent
  typedef struct 
  {
    std::string type_ ;                 // The atom type
    std::string element_ ;              // The element symbol 
    std::string amber_hybridization_ ;  // Hybridization expected by AMBER
    std::string hybridization_ ;        // Hybridization a chemist would specify
    std::string vdw_radius_ ;           // Default non-bonded atom radius
    std::string note_;                  // Note about this atom type
  } AmberAtomTypeInfo;

  const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
  {
    { "Type" , "Element" , "AMBER_Hybridization" ,"Hybridization" , "VDW_Radius" , "Note" },
  // Carbons
    { "C"    , "C"       , "sp2"                 , "sp2"          , "1.9080"           , "carbonyl group" },
    { "Cg"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "aliphatic glycan" },
    { "Cy"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "sialic acid only!" },
    { "Ck"   , "C"       , "sp2"                 , "sp2"          , "1.9080"           , "alkenes" },
    { "CT"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "aliphatic general" },
    { "Cj"   , "C"       , "sp2"                 , "sp2"          , "1.9080"           , "alkenes for adjacent double bonds" },
    { "Cp"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "bonded to an oxygen atom bonded to a phosphorus atom" },
    { "2C"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "aliphatic with two (duo) heavy atoms copied from ff12SB" },
    { "3C"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "aliphatic with three (tres) heavy atoms copied from ff12SB" },
    { "CX"   , "C"       , "sp3"                 , "sp3"          , "1.9080"           , "protein C-alpha (new to ff10)  copied from parm10" },
  // Hydrogens
    { "H"    , "H"       , "sp3"                 , "s"            , "0.6000"           , "H Bonded to nitrogen atoms" },
    { "H1"   , "H"       , "sp3"                 , "s"            , "1.3870"           , "H aliph. bond. to C with 1 electrwd. groups" },
    { "H2"   , "H"       , "sp3"                 , "s"            , "1.2870"           , "H aliph. bond. to C with 2 electrwd. groups" },
    { "Ha"   , "H"       , "sp3"                 , "s"            , "1.4590"           , "H aliph. bond. to C in alkenes eg Ck" },
    { "Hp"   , "H"       , "sp3"                 , "s"            , "1.1000"           , "H bonded to C next to positively charged group" },
    { "Hc"   , "H"       , "sp3"                 , "s"            , "1.4870"           , "H aliph. bond. to C without electrwd. groups" },
    { "Ho"   , "H"       , "sp3"                 , "s"            , "0.2000"           , "H hydroxyl group" },
    { "HW"   , "H"       , "sp3"                 , "s"            , "0.0000"           , "H TIP3P water" },
  // Nitrogens
    { "Ng"   , "N"       , "sp2"                 , "sp2"          , "1.8240"           , "sp2 N amide group" },
    { "NT"   , "N"       , "sp3"                 , "sp3"          , "1.8240"           , "sp3 N amine group" },
    { "N3"   , "N"       , "sp3"                 , "sp3"          , "1.8240"           , "sp3 N for charged amino groups (Lys, phospholipids, etc)" },
  // Oxygens
    { "Oh"   , "O"       , "sp3"                 , "sp3"          , "1.7210"           , "hydroxyl group" },
    { "Os"   , "O"       , "sp3"                 , "sp3"          , "1.6837"           , "ether" },
    { "O"    , "O"       , "sp2"                 , "sp2"          , "1.6612"           , "carbonyl group" },
    { "O2"   , "O"       , "sp2"                 , "sp2"          , "1.6612"           , "carboxyl group" },
    { "OW"   , "O"       , "sp3"                 , "sp3"          , "1.7683"           , "TIP3P water" },
    { "Oy"   , "O"       , "sp3"                 , "sp3"          , "1.6837"           , "ether - for sialic acid only!" },
  // Other
    { "S"    , "S"       , "sp3"                 , "sp3"          , "2.0000"           , "sulphates" },
    { "Sm"   , "S"       , "sp3"                 , "sp3"          , "1.7210"           , "sulfane carbohydrate linkage (-CH2-S-CH2-)" },
    { "P"    , "P"       , "sp3"                 , "sp3"          , "2.1000"           , "phosphates" },
  };
  
  // Fix the following to be in C++ syntax....
  // alias GLYCAM06ATOMTYPES=Glycam06j1AtomTypes;

  typedef struct {
    std::string type1_;  // One of the atom types
    std::string type2_;  // The other atom type
    std::string length_; // Default length in Angstroms
    std::string note_;   // Note about this bond.
  } BondLengthByTypePair;

  const BondLengthByTypePair Glycam06j1BondLengths[] =
  {
    { "Type_1" , "Type_2" , "Length" , "Note" },
    { "Cg"     , "N"      , "337.0"  , "Copy of Cg-Ng from GLYCAM06" },
    { "CT"     , "Os"     , "320.0"  , "Copy of CT-OS from parm10.dat" },
    { "2C"     , "Os"     , "320.0"  , "Copy of CT-OS from parm10.dat" },
    { "3C"     , "Os"     , "320.0"  , "Copy of CT-OS from parm10.dat " },
    { "S"      , "Ng"     , "238.0"  , "N-Sulfate - Using avg value from ZULPIF and ZULPIF01 (CSD 1.638 A)" },
    { "Cg"     , "Sm"     , "237.0"  , "Changed from 222.0 based on methanethiol nmodes" },
    { "N3"     , "H"      , "434.0"  , "Parm10" },
    { "N3"     , "Cg"     , "355.0"  , "Methanaminium (eqm value from crystal avg)" },
    { "Cj"     , "Ck"     , "350.9"  , "CRC manual for 1,3-butadiene" },
    { "Os"     , "Cj"     , "454.5"  , "copy of Os-Ck" },
    { "Os"     , "Ck"     , "454.5"  , "Methoxyethene (JACS 1993, 115, 11921)" },
    { "Cj"     , "Cj"     , "629.0"  , "copy of Ck-Ck" },
    { "Ck"     , "Ck"     , "629.0"  , "JCC 1996, 17 (5&6),669" },
    { "Cg"     , "Cj"     , "324.0"  , "copy of Cg-Ck" },
    { "Cg"     , "Ck"     , "324.0"  , "JCC 1996, 17 (5&6),669" },
    { "Cj"     , "Ha"     , "360.0"  , "copy of Ck-Ha" },
    { "Ck"     , "Ha"     , "360.0"  , "Ethane for alkenes" },
    { "Cg"     , "Hp"     , "360.0"  , "Copy of Cg-Hc" },
    { "NT"     , "H"      , "434.0"  , "Parm99" },
    { "NT"     , "Cg"     , "352.0"  , "K calculated from methyl amine (eqm value from crystal average)" },
    { "Cg"     , "Cp"     , "310.0"  , "Copy of Cg-Cg" },
    { "Cp"     , "H2"     , "340.0"  , "Copy of Cg-H1" },
    { "Cp"     , "H1"     , "340.0"  , "Copy of Cg-H1" },
    { "Cp"     , "Os"     , "285.0"  , "Copy of Cg-Os" },
    { "P"      , "Os"     , "230.0"  , "Parm94" },
    { "P"      , "O2"     , "525.0"  , "Parm10" },
    { "S"      , "Os"     , "206.0"  , "K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)" },
    { "S"      , "O2"     , "620.0"  , "K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)" },
    { "C"      , "Os"     , "450.0"  , "Parm99" },
    { "OW"     , "HW"     , "553.0"  , "TIP3P water" },
    { "HW"     , "HW"     , "553.0"  , "TIP3P water" },
    { "Cg"     , "Cg"     , "310.0"  , "Butane (gauche, and trans)" },
    { "Cg"     , "Hc"     , "340.0"  , "Parm94" },
    { "Cg"     , "H1"     , "340.0"  , "Parm94" },
    { "Cg"     , "H2"     , "340.0"  , "Parm94" },
    { "Cg"     , "Oh"     , "320.0"  , "Methanol" },
    { "Cg"     , "C"      , "220.0"  , "2-Methylpropanoate" },
    { "Cg"     , "Ng"     , "337.0"  , "N-Methylethanamide" },
    { "C"      , "O"      , "570.0"  , "Parm10" },
    { "C"      , "Ng"     , "490.0"  , "Parm94" },
    { "C"      , "O2"     , "656.0"  , "Parm10" },
    { "C"      , "Hc"     , "331.0"  , "Parm91" },
    { "C"      , "H1"     , "410.0"  , "Methanol" },
    { "Oh"     , "Ho"     , "553.0"  , "Methanol" },
    { "Ng"     , "H"      , "434.0"  , "Parm94" },
    { "Cy"     , "Oh"     , "320.0"  , "Parm94 - for sialic acid only!" },
    { "Cg"     , "Os"     , "285.0"  , "Parm94   K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)" },
    { "Cg"     , "Oy"     , "320.0"  , "Parm94 - for sialic acid only!" },
    { "Cy"     , "Cg"     , "310.0"  , "Butane (gauche, and trans) - for sialic acid only!" },
    { "Cy"     , "Os"     , "320.0"  , "Parm94 - for sialic acid only!" },
    { "Cy"     , "Oy"     , "320.0"  , "Parm94 - for sialic acid only!" },
    { "Cy"     , "C"      , "220.0"  , "2-Methylpropanoate - for sialic acid only!" }
};  

}
#endif // GLYCAM06META_HPP
