#include "../../../includes/MolecularMetadata/GLYCAM/amberatomtypeinfo.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::AmberAtomTypeInfoContainer;

AmberAtomTypeInfoContainer::AmberAtomTypeInfoContainer()
{            // const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
    amberAtomTypeInfoVector_ =
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
}
