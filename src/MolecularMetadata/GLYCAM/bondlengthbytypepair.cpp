#include "../../../includes/MolecularMetadata/GLYCAM/bondlengthbytypepair.hpp"

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::BondLengthByTypePairContainer;

BondLengthByTypePairContainer::BondLengthByTypePairContainer()
{
    //  const BondLengthByTypePair Glycam06j1BondLengths[] =
    bondLengthByTypePairVector_ =
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
