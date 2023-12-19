#include "includes/MolecularMetadata/GLYCAM/glycam06Functions.hpp"
#include "includes/CodeUtils/strings.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <algorithm> // sort

std::string GlycamMetadata::GetGlycam06ResidueLinkageCode(const std::string query)
{ // If something like 3,2,6 or 6,4,2 comes in, we need to sort it first.
    std::vector<std::string> queryLinkages = codeUtils::split(query, ',');
    std::sort(queryLinkages.begin(), queryLinkages.end());
    std::string sortedQuery = "";
    std::string deliminator = "";
    for (auto& element : queryLinkages)
    {
        sortedQuery += (deliminator + element);
        deliminator = ",";
    }
    static const std::unordered_map<std::string, std::string> Glycam06LinkageCodeLookup = {
  // linkage    , glycamLinkageCode
        { "Terminal", "0"},
        {        "1", "1"},
        {        "2", "2"},
        {        "3", "3"},
        {        "4", "4"},
        {        "5", "5"},
        {        "6", "6"},
        {        "7", "7"},
        {        "8", "8"},
        {        "9", "9"},
        {      "2,3", "Z"},
        {      "2,4", "Y"},
        {      "2,6", "X"},
        {      "3,4", "W"},
        {      "3,6", "V"},
        {      "4,6", "U"},
        {    "2,3,4", "T"},
        {    "2,3,6", "S"},
        {    "2,4,6", "R"},
        {    "3,4,6", "Q"},
        {  "2,3,4,6", "P"},
        {      "4,7", "K"},
        {      "4,8", "J"},
        {      "4,9", "I"},
        {      "7,8", "H"},
        {      "7,9", "G"},
        {      "8,9", "F"},
        {    "4,7,8", "E"},
        {    "4,7,9", "D"},
        {    "4,8,9", "C"},
        {    "7,8,9", "B"},
        {  "4,7,8,9", "A"},
        {      "2,5", "u"},
        {      "3,5", "t"},
        {      "4,5", "s"},
        {      "5,6", "r"},
        {      "5,7", "q"},
        {      "5,8", "p"},
        {      "5,9", "o"},
        {    "2,3,5", "n"},
        {    "2,5,6", "m"},
        {    "3,5,6", "l"},
        {    "4,5,7", "k"},
        {    "4,5,8", "j"},
        {    "4,5,9", "i"},
        {    "5,7,8", "h"},
        {    "5,7,9", "g"},
        {    "5,8,9", "f"},
        {  "2,3,5,6", "e"},
        {  "4,5,7,8", "d"},
        {  "4,5,7,9", "c"},
        {  "5,7,8,9", "b"},
        {"4,5,7,8,9", "a"},
 // These are ones the generation code found were missing from Xiao's table:
        {    "1,4,5", "v"},
        {    "2,4,5", "w"},
        {    "3,4,5", "x"},
        {  "1,3,4,5", "y"},
        {  "4,5,8,9", "z"},
 // These are repeats of the 2 version. Not allowing branches from the anomeric oxygen version, e.g. 1GA, means
  // we can do this for C2 anomeric sugars:
        {      "1,3", "Z"},
        {      "1,4", "Y"},
        {      "1,6", "X"},
        {    "1,3,4", "T"},
        {    "1,3,6", "S"},
        {    "1,4,6", "R"},
        {  "1,3,4,6", "P"},
        {      "1,5", "u"},
        {    "1,3,5", "n"},
        {    "1,5,6", "m"},
        {  "1,3,5,6", "e"},
 // Unused: L M N O
    };
    return codeUtils::FindStringInStringMap(sortedQuery, Glycam06LinkageCodeLookup);
}

struct NamesCodesTypes
{
    std::string residueName_;
    std::string glycamCode_;
    std::string residueType_;
};

static const std::vector<NamesCodesTypes> ResidueNameCodeLookup = {
  // residueName_, glycamCode_ , residueType_
    {        "All",   "N", "Saccharide"},
    {        "Alt",   "E", "Saccharide"},
    {        "Ara",   "A", "Saccharide"},
    {        "Fru",   "C", "Saccharide"},
    {        "Fuc",   "F", "Saccharide"},
    {        "Gal",   "L", "Saccharide"},
    {       "GalA",   "O", "Saccharide"},
    {     "GalNAc",   "V", "Saccharide"},
    {        "Glc",   "G", "Saccharide"},
    {       "GlcA",   "Z", "Saccharide"},
    {     "GlcNAc",   "Y", "Saccharide"},
    {        "Gul",   "K", "Saccharide"},
    {        "Ido",   "I", "Saccharide"},
    {       "IdoA",   "U", "Saccharide"},
    {        "Lyx",   "D", "Saccharide"},
    {        "Man",   "M", "Saccharide"},
    {     "ManNAc",   "W", "Saccharide"},
    {     "Neu5Ac",   "S", "Saccharide"},
    {     "NeuNAc",   "S", "Saccharide"},
    {        "Psi",   "P", "Saccharide"},
    {        "Qui",   "Q", "Saccharide"},
    {        "Rha",   "H", "Saccharide"},
    {        "Rib",   "R", "Saccharide"},
    {        "Sor",   "B", "Saccharide"},
    {        "Tag",   "J", "Saccharide"},
    {        "Tal",   "T", "Saccharide"},
    {        "Xyl",   "X", "Saccharide"},
    {    "GlcpNSa",  "YS", "Saccharide"},
    {    "GlcpNSb",  "Ys", "Saccharide"},
    {     "GlcpNa",  "YN", "Saccharide"},
    {     "GlcpNb",  "Yn", "Saccharide"},
    {    "GlcpNPa", "YNP", "Saccharide"},
    {    "GlcpNPb", "YnP", "Saccharide"},
    {      "Tyvpa",  "TV", "Saccharide"},
    {      "Tyvpb",  "Tv", "Saccharide"},
    {        "dUA",  "45", "Saccharide"}, // Unsaturated 4,5-unsaturated uronate.
    {"LIdopA(1C4)", "uA1", "Saccharide"}, // e.g 0uA1 with the 1 over-running.
    {"LIdopA(2SO)", "uA2", "Saccharide"},
    {"LIdopA(4C1)", "uA3", "Saccharide"},
    {   "Neup5Gca",  "GL", "Saccharide"},
    {   "NeupNGca",  "GL", "Saccharide"},
    {      "KDNpa",  "KN", "Saccharide"},
    {      "KDOpa",  "KO", "Saccharide"},
    {      "KDNpb",  "Kn", "Saccharide"},
    {      "KDOpb",  "Ko", "Saccharide"},
    {      "Bacpb",  "BC", "Saccharide"},
    {        "ROH", "ROH",   "Aglycone"},
    {         "OH", "ROH",   "Aglycone"},
    {        "OME", "OME",   "Aglycone"},
    {       "OtBu", "TBT",   "Aglycone"},
    {       "OThr", "OLT", "Amino-acid"},
    {       "OSer", "OST", "Amino-acid"},
    {       "OTyr", "OLY", "Amino-acid"},
    {       "NAsn", "NLN", "Amino-acid"},
    {     "Sulpho", "SO3", "Derivative"}, // Phosfo? Phosno.
    {    "Phospho", "PO3", "Derivative"},
    {     "Methyl", "MEX", "Derivative"},
    {     "Acetyl", "ACX", "Derivative"},
    {          "S", "SO3", "Derivative"},
    {          "P", "PO3", "Derivative"},
    {         "Me", "MEX", "Derivative"},
    {         "Ac", "ACX", "Derivative"},
    {          "A", "ACX", "Derivative"},
};

#include <iostream>

std::string GlycamMetadata::GetNameForCode(const std::string query)
{
    for (auto& entry : ResidueNameCodeLookup)
    {
        if (entry.glycamCode_ == query)
        {
            return entry.residueName_;
        }
    }
    return "";
}

std::string GlycamMetadata::GetCodeForName(const std::string query)
{
    // e.g. input = Fuc, output = F.
    for (auto& entry : ResidueNameCodeLookup)
    {
        if (entry.residueName_ == query)
        {
            return entry.glycamCode_;
        }
    }
    return "";
}

std::string GlycamMetadata::GetTypeForCode(const std::string query)
{
    // e.g. input = Fuc, output = F.
    for (auto& entry : ResidueNameCodeLookup)
    {
        if (entry.glycamCode_ == query)
        {
            return entry.residueType_;
        }
    }
    return "";
}

std::string GlycamMetadata::GetDescriptiveNameForGlycamResidueName(const std::string residueNameInGLYCAMFormat)
{
    if (residueNameInGLYCAMFormat.length() < 3)
    {
        std::string message = "Residue name is too short: " + residueNameInGLYCAMFormat;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    std::string type = GlycamMetadata::GetTypeForCode(residueNameInGLYCAMFormat);
    if (!type.empty() && type != "Saccharide")
    { // This just works for e.g. ROH, NLN
        return GlycamMetadata::GetNameForCode(residueNameInGLYCAMFormat);
    }
    // For glycam we made it harder:
    char secondLetter = residueNameInGLYCAMFormat.at(1);
    char thirdLetter  = residueNameInGLYCAMFormat.at(2);
    // char secondLetter = residueNameInGLYCAMFormat[residueNameInGLYCAMFormat.size() - 2]; // sush.
    //    std::cout << "Letters are : " << secondLetter << " and " << thirdLetter << std::endl;

    char configurationDvsL = 'D';
    if (islower(secondLetter))
    { // Lower case means L type sugar. E.g. 0fA is LFucpA
        configurationDvsL = 'L';
    }
    secondLetter = toupper(secondLetter); // Make it upper for matching to resname below.

    bool addRingformAndAnomer  = true;
    char ringform              = 'p';                  // Default, may change below if furanose
    char anomericConfiguration = tolower(thirdLetter); // Default, may change below if furanose
    std::string query(1, secondLetter);                // Default, may change below
    if (thirdLetter == 'U' || thirdLetter == 'D' || thirdLetter == 'A' || thirdLetter == 'B' || thirdLetter == 'X')
    {
        if (thirdLetter == 'U' || thirdLetter == 'D')
        {                   // E.g. 0GU is DGalfb
            ringform = 'f'; // Furanose
        }
        if (thirdLetter == 'U')
        { // Furanoses: U for up means b. D for down means a.
            anomericConfiguration = 'b';
        }
        if (thirdLetter == 'D')
        {
            anomericConfiguration = 'a';
        }
    }
    else // thirdLetter isn't one of the above.
    {
        std::stringstream ss;
        ss << secondLetter << thirdLetter;
        query                = ss.str(); // Weirdos like BC for Bac
        addRingformAndAnomer = false;
    }
    // std::cout << "Query is " << query << std::endl;
    // std::string match = configurationDvsL + "Hex" + anomericConfiguration; // Default is Hex?
    std::string name = GlycamMetadata::GetNameForCode(query);
    if (addRingformAndAnomer)
    {
        return configurationDvsL + name + ringform + anomericConfiguration;
    }
    return configurationDvsL + name;
}

struct ChargeAdjustmentInfo
{
    std::string glycamResidueCode_;
    std::string adjustmentAtom_;
    double charge_;
};

std::vector<ChargeAdjustmentInfo> ChargeAdjustmentInfo = {
  // residueCode_, adjustmentAtom_ , charge_
    {"ACX", "C",  0.008},
    {"MEX", "C", -0.039},
    {"SO3", "O", +0.031},
};

double GlycamMetadata::GetAdjustmentCharge(std::string queryResidueName)
{ // e.g. input = SO3, output = +0.008.
    for (const auto& entry : ChargeAdjustmentInfo)
    {
        if (entry.glycamResidueCode_ == queryResidueName)
        {
            return entry.charge_;
        }
    }
    return 0.0;
}

std::string GlycamMetadata::GetAdjustmentAtom(std::string queryResidueName)
{ // e.g. input = SO3, output = O.
    for (const auto& entry : ChargeAdjustmentInfo)
    {
        if (entry.glycamResidueCode_ == queryResidueName)
        {
            return entry.adjustmentAtom_;
        }
    }
    return "O";
}

static const std::unordered_map<std::string, std::string> glycam06DerivativeAglyconeConnectionAtomLookup_ = {
  // Glycam06 Residue Name , Linking atom
  // Aglycones
    {"ROH",  "O1"},
    {"OME",   "O"},
    {"TBT",  "O1"},
 // Aglycones amino acid
    {"NLN", "ND2"},
    {"OLT", "OG1"},
    {"OLS",  "OG"},
 // Derivatives
    {"SO3",  "S1"},
    {"MEX", "CH3"},
    {"ACX", "C1A"},
};

std::string GlycamMetadata::GetConnectionAtomForResidue(const std::string query)
{
    std::string result = codeUtils::FindStringInStringMap(query, glycam06DerivativeAglyconeConnectionAtomLookup_);
    if (result.empty())
    {
        std::string message = "This derivative or aglycone residue is not supported: " + query;
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return result;
}
