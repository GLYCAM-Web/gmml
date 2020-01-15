#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residuecodes.hpp"
#include <iostream> // for cout, can remove after debug
#include <locale> // for isLower()
#include <sstream> // for stringstream

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////

using gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNamesToCodesLookup;
//using gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNamesToCodesLookupContainerVector;

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

std::string Glycam06ResidueNamesToCodesLookup::GetResidueForCode(std::string residueNameInGLYCAMFormat)
{ // lettersFromGLYCAMFormatResidueName. E.g. 0GB, 3YB, 3BC.
    if (residueNameInGLYCAMFormat.length() < 2)
    { // Problem.
        std::cout << "Problem" << std::endl;
        return "Problem";
    }
    std::cout << "Input name is " << residueNameInGLYCAMFormat << std::endl;
    char lastLetter = *residueNameInGLYCAMFormat.rbegin(); //points to the last character
    char middleLetter = residueNameInGLYCAMFormat[residueNameInGLYCAMFormat.size() - 2];
    std::cout << "last and midd letter are : " << lastLetter << " and " << middleLetter << std::endl;

    char configurationDvsL = 'D';
    if (islower(middleLetter)) // Lower case means L type sugar. E.g. 0fA is LFucpA
        configurationDvsL = 'L';
    middleLetter = toupper(middleLetter); // Make it upper for matching below.

    bool addAnomer = true;                             // Default, may change below
    char ringform = 'p';                               // Default, may change below
    char anomericConfiguration = tolower(lastLetter);  // Default, may change below
    std::string query(1, middleLetter);                // Default, may change below
    if (lastLetter == 'U' || lastLetter == 'D' || lastLetter == 'A' || lastLetter == 'B' || lastLetter == 'X')
    {
        if (lastLetter == 'U' || lastLetter == 'D') // E.g. 0GU is DGalfb
            ringform = 'f'; // Furanose
        if (lastLetter == 'U') // Furanoses: U for up means b. D for down means a.
            anomericConfiguration = 'b';
        if (lastLetter == 'D')
            anomericConfiguration = 'a';
    }
    else // lastLetter isn't one of the above.
    {
        std::stringstream ss;
        ss << middleLetter << lastLetter;
        query = ss.str(); // Weirdos like BC for Bac
        addAnomer = false;
    }
    std::cout << "Query is " << query << std::endl;
    std::string match = configurationDvsL + "Hex" + anomericConfiguration; // Default is Hex?
    for (std::pair<std::string, std::string> elem : glycam06ResidueNamesToCodesLookupMap_)
    {
        if (elem.second.compare(query)==0)
        {
             if (addAnomer)
                 match = configurationDvsL + elem.first + ringform + anomericConfiguration;
             else
                 match = configurationDvsL + elem.first + ringform;
             std::cout << "Match is " << match << std::endl;
             return match;
        }
    }
    return match;
}


//////////////////////////////////////////////////////////
//                    PRIVATE FUNCTIONS                 //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                    INITIALIZER                       //
//////////////////////////////////////////////////////////


Glycam06ResidueNamesToCodesLookup::Glycam06ResidueNamesToCodesLookup()
{
    glycam06ResidueNamesToCodesLookupMap_ =
    {
        {"All", "N"},
        {"Alt", "E"},
        {"Ara", "A"},
        {"Bac", "BC"},
        {"Fru", "C"},
        {"Fuc", "F"},
        {"Gal", "L"},
        {"GalA", "O"},
        {"GalNAc", "V"},
        {"Glc", "G"},
        {"GlcA", "Z"},
        {"GlcNAc", "Y"},
        {"Gul", "K"},
        {"Ido", "I"},
        {"IdoA", "U"},
        {"IdoA(1C4)", "U1"},
        {"IdoA(2SO)", "U2"},
        {"IdoA(4C1)", "U3"},
        {"Lyx", "D"},
        {"Man", "M"},
        {"ManNAc", "W"},
        {"Neu5Ac", "S"},
        {"Neu5Gc", "GL"},
        {"NeuNAc", "S"},
        {"NeuNGc", "GL"},
        {"Psi", "P"},
        {"Qui", "Q"},
        {"Rha", "H"},
        {"Rib", "R"},
        {"Sor", "B"},
        {"Tag", "J"},
        {"Tal", "T"},
        {"Tyv", "TV"},
        {"Xyl", "X"},
        {"GlcNS", "Y"},
        {"U", "045"}, // ? I don't know what these are. OG.
        {"S", "245"}, // ? I don't know what these are. OG.
        {"KDN", "KN"},
        {"KDO", "KO"},
        {"ROH", "OH"},
        {"OME", "ME"},
        {"OtBu", "BT"},
        {"OThr", "LT"},
        {"OSer", "ST"},
        {"OTyr", "LY"},
        {"NAsn", "LN"},
        {"Sulfo", "O3"},
        {"Methyl", "EX"},
    };
}
