#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp"
#include <iostream> // for cout, can remove after debug
#include <locale> // for isLower()
#include <sstream> // for stringstream

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNamesToCodesLookupContainer;

//////////////////////////////////////////////////////////
//                      QUERY FUNCTIONS                 //
//////////////////////////////////////////////////////////

//std::string Glycam06ResidueNamesToCodesLookup::GetCodeForResidue(std::string query)
//{
//    std::regex regex1(query, std::regex_constants::ECMAScript);
//    std::string match;
//    //            std::string str("1231");
//    //            std::regex r("^(\\d)\\d"); // entire match will be 2 numbers
//    //            std::smatch m;
//    //            std::regex_search(str, m, r);
//            //for(auto v: m) std::cout << v << std::endl;
//    //            std::smatch results1;
//    //            std::string string1 = linking_atom2->GetName();
//    //            bool well = std::regex_search(string1, results1, regex1);
//    //            std::cout << well << "\n";
//    //            for(auto v: results1) std::cout << v << std::endl;
//    // Iterate over the multimap using range based for loop
//    for (std::pair<std::string, std::string> elem : glycam06ResidueNamesToCodesLookupMap_)
//    {
//        if (elem.first.compare(query)==0)
//        {
//            match.push_back(elem.second);
//        }
//        //std::cout << elem.first << " :: " << elem.second << std::endl;
//    }
//    return match;
//}

// This won't be correct all the time, but should be ok for displaying labels to humans.
std::string Glycam06ResidueNamesToCodesLookupContainer::GetResidueForCode(std::string residueNameInGLYCAMFormat)
{   // residueNameInGLYCAMFormat. E.g. 0GB, 3YB, 3BC.
    std::cout << "Input name is " << residueNameInGLYCAMFormat << std::endl;
    if (residueNameInGLYCAMFormat.length() < 3)
        return "PrObLeM";

    for (const auto& entry : ResidueNamesCodesTypesLookupTable_)
    { // residueName_, glycamCode_ , residueType_
        if ( (residueNameInGLYCAMFormat == entry.glycamCode_) && (entry.residueType_ != "Saccharide") )
            return entry.residueName_;
    }

    char secondLetter = residueNameInGLYCAMFormat.at(1);
    char thirdLetter = residueNameInGLYCAMFormat.at(2);
    //char secondLetter = residueNameInGLYCAMFormat[residueNameInGLYCAMFormat.size() - 2]; // sush.
    std::cout << "Letters are : " << secondLetter << " and " << thirdLetter << std::endl;

    char configurationDvsL = 'D';
    if (islower(secondLetter)) // Lower case means L type sugar. E.g. 0fA is LFucpA
        configurationDvsL = 'L';
    secondLetter = toupper(secondLetter); // Make it upper for matching to resname below.

    bool addRingformAndAnomer = true;
    char ringform = 'p';                               // Default, may change below if furanose
    char anomericConfiguration = tolower(thirdLetter);  // Default, may change below if furanose
    std::string query(1, secondLetter);                // Default, may change below
    if (thirdLetter == 'U' || thirdLetter == 'D' || thirdLetter == 'A' || thirdLetter == 'B' || thirdLetter == 'X')
    {
        if (thirdLetter == 'U' || thirdLetter == 'D') // E.g. 0GU is DGalfb
            ringform = 'f'; // Furanose
        if (thirdLetter == 'U') // Furanoses: U for up means b. D for down means a.
            anomericConfiguration = 'b';
        if (thirdLetter == 'D')
            anomericConfiguration = 'a';
    }
    else // lastLetter isn't one of the above.
    {
        std::stringstream ss;
        ss << secondLetter << thirdLetter;
        query = ss.str(); // Weirdos like BC for Bac
        addRingformAndAnomer = false;
    }
    std::cout << "Query is " << query << std::endl;
    std::string match = configurationDvsL + "Hex" + anomericConfiguration; // Default is Hex?
    for (const auto& entry : ResidueNamesCodesTypesLookupTable_)
    { // residueName_, glycamCode_ , residueType_
        if (entry.glycamCode_ == query)
        {
             if (addRingformAndAnomer)
                 match = configurationDvsL + entry.residueName_ + ringform + anomericConfiguration;
             else
                 match = configurationDvsL + entry.residueName_;
             std::cout << "Match is " << match << std::endl;
             return match;
        }
    }
    return match;
}

//std::string Glycam06ResidueNamesToCodesLookupContainer::GetResidueForCode(std::string residueNameInGLYCAMFormat)
//{ // lettersFromGLYCAMFormatResidueName. E.g. 0GB, 3YB, 3BC.
//    if (residueNameInGLYCAMFormat.length() != 3)
//    { // Problem.
//        std::cout << "Problem" << std::endl;
//        return "PrObLeM";
//    }

//    char firstLetter = *residueNameInGLYCAMFormat.begin();
//    char lastLetter = *residueNameInGLYCAMFormat.rbegin(); //points to the last character
//    char middleLetter = residueNameInGLYCAMFormat[residueNameInGLYCAMFormat.size() - 2]; // sush.
//    std::cout << "Letters are : " << firstLetter << ", " << middleLetter << " and " << lastLetter << std::endl;

//    char configurationDvsL = 'D';
//    if (islower(middleLetter)) // Lower case means L type sugar. E.g. 0fA is LFucpA
//        configurationDvsL = 'L';
//    middleLetter = toupper(middleLetter); // Make it upper for matching to resname below.

//    bool addAnomer = true;                             // Default, may change below if derivative or aglycone
//    char ringform = 'p';                               // Default, may change below if furanose
//    char anomericConfiguration = tolower(lastLetter);  // Default, may change below if furanose
//    std::string query(1, middleLetter);                // Default, may change below
//    if (lastLetter == 'U' || lastLetter == 'D' || lastLetter == 'A' || lastLetter == 'B' || lastLetter == 'X')
//    {
//        if (lastLetter == 'U' || lastLetter == 'D') // E.g. 0GU is DGalfb
//            ringform = 'f'; // Furanose
//        if (lastLetter == 'U') // Furanoses: U for up means b. D for down means a.
//            anomericConfiguration = 'b';
//        if (lastLetter == 'D')
//            anomericConfiguration = 'a';
//    }
//    else // lastLetter isn't one of the above.
//    {
//        std::stringstream ss;
//        ss << middleLetter << lastLetter;
//        query = ss.str(); // Weirdos like BC for Bac
//        addAnomer = false;
//    }
//    std::cout << "Query is " << query << std::endl;
//    std::string match = configurationDvsL + "Hex" + anomericConfiguration; // Default is Hex?
//    for (std::pair<std::string, std::string> elem : ResidueNamesCodesTypesLookupTable_)
//    {
//        if (elem.second.compare(query)==0)
//        {
//             if (addAnomer)
//                 match = configurationDvsL + elem.first + ringform + anomericConfiguration;
//             else
//                 match = configurationDvsL + elem.first + ringform;
//             std::cout << "Match is " << match << std::endl;
//             return match;
//        }
//    }
//    return match;
//}


//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                    INITIALIZER                       //
//////////////////////////////////////////////////////////


//DihedralAngleDataContainer::DihedralAngleDataContainer()
//{   // const AmberAtomTypeInfo Glycam06j1AtomTypes[] =
//    dihedralAngleDataVector_ =
//    { //Regex1, Regex2   , Name   , Angle  , Upper  , Lower  , Weight, Entry Type    , Tag  , D , R , Residue 1 type , Residue 2 type,           , Atom names                                                               // Atom names this applies to
//        { "C1"

Glycam06ResidueNamesToCodesLookupContainer::Glycam06ResidueNamesToCodesLookupContainer()
{
    ResidueNamesCodesTypesLookupTable_ =
    {// residueName_, glycamCode_ , residueType_
        {"All"         , "N"  , "Saccharide"},
        {"Alt"         , "E"  , "Saccharide"},
        {"Ara"         , "A"  , "Saccharide"},
        {"Fru"         , "C"  , "Saccharide"},
        {"Fuc"         , "F"  , "Saccharide"},
        {"Gal"         , "L"  , "Saccharide"},
        {"GalA"        , "O"  , "Saccharide"},
        {"GalNAc"      , "V"  , "Saccharide"},
        {"Glc"         , "G"  , "Saccharide"},
        {"GlcA"        , "Z"  , "Saccharide"},
        {"GlcNAc"      , "Y"  , "Saccharide"},
        {"Gul"         , "K"  , "Saccharide"},
        {"Ido"         , "I"  , "Saccharide"},
        {"IdoA"        , "U"  , "Saccharide"},
        {"Lyx"         , "D"  , "Saccharide"},
        {"Man"         , "M"  , "Saccharide"},
        {"ManNAc"      , "W"  , "Saccharide"},
        {"Neu5Ac"      , "S"  , "Saccharide"},
        {"NeuNAc"      , "S"  , "Saccharide"},
        {"Psi"         , "P"  , "Saccharide"},
        {"Qui"         , "Q"  , "Saccharide"},
        {"Rha"         , "H"  , "Saccharide"},
        {"Rib"         , "R"  , "Saccharide"},
        {"Sor"         , "B"  , "Saccharide"},
        {"Tag"         , "J"  , "Saccharide"},
        {"Tal"         , "T"  , "Saccharide"},
        {"Xyl"         , "X"  , "Saccharide"},
        {"GlcNS"       , "Y"  , "Saccharide"},
        {"Tyvp"        , "TV" , "Saccharide"}, // Can be alpha or beta. Feck.
        {"45UIdopa"    , "45" , "Saccharide"}, // Unsaturated 4,5-unsaturated uronate
        {"LIdoA(1C4)pa", "UA1", "Saccharide"}, // e.g 4uA1 This is a mistake imo. Should be U1, not UA1.
        {"LIdoA(2SO)pa", "UA2", "Saccharide"},
        {"LIdoA(4C1)pa", "UA3", "Saccharide"},
        {"Neu5Gcpa"    , "GL" , "Saccharide"},
        {"NeuNGcpa"    , "GL" , "Saccharide"},
        {"KDNpa"       , "KN" , "Saccharide"},
        {"KDOpa"       , "KO" , "Saccharide"},
        {"Bacpb"       , "BC" , "Saccharide"},
        {"ROH"         , "ROH", "Aglycone"  },
        {"OME"         , "OME", "Aglycone"  },
        {"OtBu"        , "TBT", "Aglycone"  },
        {"OThr"        , "OLT", "Amino-acid"},
        {"OSer"        , "OST", "Amino-acid"},
        {"OTyr"        , "OLY", "Amino-acid"},
        {"NAsn"        , "NLN", "Amino-acid"},
        {"Sulpho"      , "SO3", "Derivative"}, // Phosfo? Phosno.
        {"Phospho"     , "PO3", "Derivative"},
        {"Methyl"      , "MEX", "Derivative"},
        {"Acetyl"      , "ACX", "Derivative"},
    };
}
