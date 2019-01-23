#ifndef GLYCAM06_META_HPP
#define GLYCAM06_META_HPP

/* File glycam06.hpp begun on 16 June 2018 by BLFoley */

#include <string>
#include <vector>
#include <set>

namespace gmml
{
	namespace MolecularMetadata
	{
		namespace GLYCAM
		{
		  // Move this to an amber space once the proper design is apparent
		  typedef struct
		  {
		    std::string type_ ;                 // The atom type
		    std::string element_ ;              // The element symbol
		    std::string amber_hybridization_ ;  // Hybridization expected by AMBER
		    std::string hybridization_ ;        // Hybridization a chemist would specify
		    double vdw_radius_ ;                // Default non-bonded atom radius
		    std::string note_;                  // Note about this atom type
		  } AmberAtomTypeInfo;

		  const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
		  {
		  // The '1.0' in the 5th column is because I can't have column headers if the numbers aren't strings.
		  //{ "Type" , "Element" , "AMBER_Hybridization" ,"Hybridization" , "VDW_Radius"     , "Note" },
		  // Carbons
		    { "C"    , "C"       , "sp2"                 , "sp2"          , 1.9080           , "carbonyl group" },
		    { "Cg"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "aliphatic glycan" },
		    { "Cy"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "sialic acid only!" },
		    { "Ck"   , "C"       , "sp2"                 , "sp2"          , 1.9080           , "alkenes" },
		    { "CT"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "aliphatic general" },
		    { "Cj"   , "C"       , "sp2"                 , "sp2"          , 1.9080           , "alkenes for adjacent double bonds" },
		    { "Cp"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "bonded to an oxygen atom bonded to a phosphorus atom" },
		    { "2C"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "aliphatic with two (duo) heavy atoms copied from ff12SB" },
		    { "3C"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "aliphatic with three (tres) heavy atoms copied from ff12SB" },
		    { "CX"   , "C"       , "sp3"                 , "sp3"          , 1.9080           , "protein C-alpha (new to ff10)  copied from parm10" },
		  // Hydrogens
		    { "H"    , "H"       , "sp3"                 , "s"            , 0.6000           , "H Bonded to nitrogen atoms" },
		    { "H1"   , "H"       , "sp3"                 , "s"            , 1.3870           , "H aliph. bond. to C with 1 electrwd. groups" },
		    { "H2"   , "H"       , "sp3"                 , "s"            , 1.2870           , "H aliph. bond. to C with 2 electrwd. groups" },
		    { "Ha"   , "H"       , "sp3"                 , "s"            , 1.4590           , "H aliph. bond. to C in alkenes eg Ck" },
		    { "Hp"   , "H"       , "sp3"                 , "s"            , 1.1000           , "H bonded to C next to positively charged group" },
		    { "Hc"   , "H"       , "sp3"                 , "s"            , 1.4870           , "H aliph. bond. to C without electrwd. groups" },
		    { "Ho"   , "H"       , "sp3"                 , "s"            , 0.2000           , "H hydroxyl group" },
		    { "HW"   , "H"       , "sp3"                 , "s"            , 0.0000           , "H TIP3P water" },
		  // Nitrogens
		    { "Ng"   , "N"       , "sp2"                 , "sp2"          , 1.8240           , "sp2 N amide group" },
		    { "NT"   , "N"       , "sp3"                 , "sp3"          , 1.8240           , "sp3 N amine group" },
		    { "N3"   , "N"       , "sp3"                 , "sp3"          , 1.8240           , "sp3 N for charged amino groups (Lys, phospholipids, etc)" },
		  // Oxygens
		    { "Oh"   , "O"       , "sp3"                 , "sp3"          , 1.7210           , "hydroxyl group" },
		    { "Os"   , "O"       , "sp3"                 , "sp3"          , 1.6837           , "ether" },
		    { "O"    , "O"       , "sp2"                 , "sp2"          , 1.6612           , "carbonyl group" },
		    { "O2"   , "O"       , "sp2"                 , "sp2"          , 1.6612           , "carboxyl group" },
		    { "OW"   , "O"       , "sp3"                 , "sp3"          , 1.7683           , "TIP3P water" },
		    { "Oy"   , "O"       , "sp3"                 , "sp3"          , 1.6837           , "ether - for sialic acid only!" },
		  // Other
		    { "S"    , "S"       , "sp3"                 , "sp3"          , 2.0000           , "sulphates" },
		    { "Sm"   , "S"       , "sp3"                 , "sp3"          , 1.7210           , "sulfane carbohydrate linkage (-CH2-S-CH2-)" },
		    { "P"    , "P"       , "sp3"                 , "sp3"          , 2.1000           , "phosphates" },
		  };

		  // Fix the following to be in C++ syntax....
		  // alias GLYCAM06ATOMTYPES=Glycam06j1AtomTypes;

		  typedef struct {
		    std::string type1_;  // One of the atom types
		    std::string type2_;  // The other atom type
		    double length_;      // Default length in Angstroms
		    std::string note_;   // Note about this bond.
		  } BondLengthByTypePair;

		  const BondLengthByTypePair Glycam06j1BondLengths[] =
		  {
		  //{ "Type_1" , "Type_2" , "Length" , "Note" },
		    { "Cg"     , "N"      , 1.450    , "Copy of Cg-Ng from GLYCAM06" },
		    { "CT"     , "Os"     , 1.410    , "Copy of CT-OS from parm10.dat" },
		    { "2C"     , "Os"     , 1.410    , "Copy of CT-OS from parm10.dat" },
		    { "3C"     , "Os"     , 1.410    , "Copy of CT-OS from parm10.dat " },
		    { "S"      , "Ng"     , 1.675    , "N-Sulfate - Using avg value from ZULPIF and ZULPIF01 (CSD 1.638 A)" },
		    { "Cg"     , "Sm"     , 1.810    , "Changed from 222.0 based on methanethiol nmodes" },
		    { "N3"     , "H"      , 1.010    , "Parm10" },
		    { "N3"     , "Cg"     , 1.490    , "Methanaminium (eqm value from crystal avg)" },
		    { "Cj"     , "Ck"     , 1.467    , "CRC manual for 1,3-butadiene" },
		    { "Os"     , "Cj"     , 1.359    , "copy of Os-Ck" },
		    { "Os"     , "Ck"     , 1.359    , "Methoxyethene (JACS 1993, 115, 11921)" },
		    { "Cj"     , "Cj"     , 1.337    , "copy of Ck-Ck" },
		    { "Ck"     , "Ck"     , 1.337    , "JCC 1996, 17 (5&6),669" },
		    { "Cg"     , "Cj"     , 1.514    , "copy of Cg-Ck" },
		    { "Cg"     , "Ck"     , 1.514    , "JCC 1996, 17 (5&6),669" },
		    { "Cj"     , "Ha"     , 1.095    , "copy of Ck-Ha" },
		    { "Ck"     , "Ha"     , 1.095    , "Ethane for alkenes" },
		    { "Cg"     , "Hp"     , 1.095    , "Copy of Cg-Hc" },
		    { "NT"     , "H"      , 1.010    , "Parm99" },
		    { "NT"     , "Cg"     , 1.470    , "K calculated from methyl amine (eqm value from crystal average)" },
		    { "Cg"     , "Cp"     , 1.520    , "Copy of Cg-Cg" },
		    { "Cp"     , "H2"     , 1.090    , "Copy of Cg-H1" },
		    { "Cp"     , "H1"     , 1.090    , "Copy of Cg-H1" },
		    { "Cp"     , "Os"     , 1.460    , "Copy of Cg-Os" },
		    { "P"      , "Os"     , 1.610    , "Parm94" },
		    { "P"      , "O2"     , 1.480    , "Parm10" },
		    { "S"      , "Os"     , 1.589    , "K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)" },
		    { "S"      , "O2"     , 1.440    , "K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)" },
		    { "C"      , "Os"     , 1.323    , "Parm99" },
		    { "OW"     , "HW"     , 0.9572   , "TIP3P water" },
		    { "HW"     , "HW"     , 1.5136   , "TIP3P water" },
		    { "Cg"     , "Cg"     , 1.520    , "Butane (gauche, and trans)" },
		    { "Cg"     , "Hc"     , 1.090    , "Parm94" },
		    { "Cg"     , "H1"     , 1.090    , "Parm94" },
		    { "Cg"     , "H2"     , 1.090    , "Parm94" },
		    { "Cg"     , "Oh"     , 1.430    , "Methanol" },
		    { "Cg"     , "C"      , 1.530    , "2-Methylpropanoate" },
		    { "Cg"     , "Ng"     , 1.450    , "N-Methylethanamide" },
		    { "C"      , "O"      , 1.229    , "Parm10" },
		    { "C"      , "Ng"     , 1.335    , "Parm94" },
		    { "C"      , "O2"     , 1.250    , "Parm10" },
		    { "C"      , "Hc"     , 1.090    , "Parm91" },
		    { "C"      , "H1"     , 1.092    , "Methanol" },
		    { "Oh"     , "Ho"     , 0.960    , "Methanol" },
		    { "Ng"     , "H"      , 1.010    , "Parm94" },
		    { "Cy"     , "Oh"     , 1.410    , "Parm94 - for sialic acid only!" },
		    { "Cg"     , "Os"     , 1.460    , "Parm94   K calculated from methyl sulphate (eqm value from THEOCHEM 395/396 (1997) 107-122)" },
		    { "Cg"     , "Oy"     , 1.410    , "Parm94 - for sialic acid only!" },
		    { "Cy"     , "Cg"     , 1.520    , "Butane (gauche, and trans) - for sialic acid only!" },
		    { "Cy"     , "Os"     , 1.410    , "Parm94 - for sialic acid only!" },
		    { "Cy"     , "Oy"     , 1.410    , "Parm94 - for sialic acid only!" },
		    { "Cy"     , "C"      , 1.530    , "2-Methylpropanoate - for sialic acid only!" }
		  };
		}
	}
}
#endif // GLYCAM06META_HPP
