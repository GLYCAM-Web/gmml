#include "includes/CentralDataStructure/CondensedSequence/parsedResidue.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <sstream>

using cds::ResidueType;
using cdsCondensedSequence::ParsedResidue;

ParsedResidue::ParsedResidue(std::string residueString, ResidueType specifiedType) : fullResidueString_(residueString)
//: Node(residueString), fullResidueString_ (residueString)
{
    this->setName(residueString);
    this->ParseResidueStringIntoComponents(residueString, specifiedType);
}

ParsedResidue::ParsedResidue(std::string residueString, ParsedResidue* neighbor, ResidueType specifiedType)
    //: Node(residueString), fullResidueString_ (residueString)
    : fullResidueString_(residueString)
{
    this->setName(residueString);
    this->ParseResidueStringIntoComponents(residueString, specifiedType);
    this->AddLinkage(neighbor);
}

void ParsedResidue::AddLinkage(ParsedResidue* otherRes)
{
    if (this->GetType() == ResidueType::Sugar)
    {
        this->addParent(this->GetConfiguration() + this->GetLinkage(), otherRes);
    }
    else
    {
        this->addParent(this->GetLinkage(), otherRes);
    }
    // std::cout << this->getName() << " with linkage " << this->GetLinkage() << " has parent " << otherRes->getName()
    // << std::endl;
    return;
}

std::string ParsedResidue::GetLink() const
{
    std::string linkage = this->GetLinkage();
    switch (this->GetType())
    {
        case (ResidueType::Sugar):
            return linkage.substr(linkage.size() - 1, 1);
        case (ResidueType::Derivative):
            return linkage.substr(linkage.size() - 1, 1);
        case (ResidueType::Deoxy):
            return linkage.substr(0, 1);
        default:
            return "0";
    }
}

std::vector<ParsedResidue*> ParsedResidue::GetChildren() const
{
    std::vector<ParsedResidue*> resRet;
    for (auto& currNodeRes : this->getChildren())
    {
        resRet.push_back(static_cast<ParsedResidue*>(currNodeRes));
    }
    return resRet;
}

std::vector<ParsedResidue*> ParsedResidue::GetParents() const
{
    std::vector<ParsedResidue*> resRet;
    for (auto& currNodeRes : this->getParents())
    {
        resRet.push_back(static_cast<ParsedResidue*>(currNodeRes));
    }
    return resRet;
}

ParsedResidue* ParsedResidue::GetParent() const
{
    std::vector<cds::Residue*> parents = this->getParents();
    if (parents.empty())
    {
        std::string message = "Bad situation: Parsed residue named " + this->getName() + " has no parent\n";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    return static_cast<ParsedResidue*>(parents.front());
}

std::string ParsedResidue::GetChildLinkagesForGlycamResidueNaming() const
{
    std::string linkages;
    for (auto& child : this->GetChildren())
    { // For glycam residue name, e.g. 2YB, do not want deoxy linkages to impact the residue name.
        if (child->GetType() != ResidueType::Deoxy)
        {
            linkages += (child->GetLink() + ",");
        }
    }
    if (linkages.empty())
    {
        linkages = "Terminal";
    }
    else // Erase the last ","
    {
        linkages.erase(std::prev(linkages.end()));
    }
    return linkages;
}

std::string ParsedResidue::GetName(const bool withLabels) const
{
    if (withLabels)
    {
        return this->findLabelContaining("&Label=");
    }
    return this->GetIsomer() + this->GetResidueName() + this->GetRingType() + this->GetResidueModifier() +
           this->GetRingShape();
}

std::string ParsedResidue::GetLinkageName(const bool withLabels) const
{ // Should only ever be zero or one inEdges in my current design.
    // e.g  Galb1-4[Glcb1-3]Manb1-OH. OH is the parent to Manb, the linkage is b1-.
    // Manb is parent (has out edges) to both Glc and Gal, but they only have one in edge each from Manb. The inedge
    // from Man to Gal has the linkage name 1-4.
    for (auto& linkage : this->getInEdges())
    {
        if (withLabels)
        {
            return linkage->findLabelContaining("&Label=");
        }
        else
        {
            return linkage->getLabel();
        }
    }
    return ""; // aglycone/reducing terminal will not have linkage.
}

void ParsedResidue::ParseResidueStringIntoComponents(std::string residueString, ResidueType specifiedType)
{
    // gmml::log(__LINE__, __FILE__, gmml::INF, "PARSING RESIDUE: " + residueString);
    // Set defaults:
    this->SetIsomer("");
    this->SetResidueName("");
    this->SetRingType("");
    this->SetRingShape("");
    this->SetResidueModifier("");
    this->SetConfiguration("");
    this->SetLinkage("");
    this->SetType(specifiedType);
    if ((residueString.find('-') != std::string::npos) || (specifiedType == ResidueType::Sugar))
    { // E.g. DManpNAca1-4 . Isomer (D or L), residueName (ManNAc), ring type (f or p), configuration (a or b), linkage
      // (1-4)
        // Reading from front.
        this->SetType(ResidueType::Sugar);
        // Assumptions
        size_t residueStart  = 1; // e.g. Gal, Glc, Ido
        size_t modifierStart = 5; // E.g. NAc, A, A(1C4)
        // Checks
        std::string isomer   = residueString.substr(0, 1);
        if ((isomer == "D") || (isomer == "L"))
        {
            this->SetIsomer(isomer);
        }
        else
        {
            residueStart--;
            modifierStart--;
        }
        this->SetResidueName(residueString.substr(residueStart, 3));
        size_t ringPosition  = (residueStart + 3);
        std::string ringType = residueString.substr(ringPosition, 1);
        if ((ringType == "p") || (ringType == "f"))
        {
            this->SetRingType(ringType);
        }
        else
        {
            modifierStart--;
        }
        // Find the dash, read around it.
        size_t dashPosition = residueString.find('-');
        if (dashPosition == std::string::npos) // There is no -. e.g. Fru in DGlcpa1-2DFrufb
        {
            dashPosition = residueString.size() + 1;
        }
        else
        {
            this->SetLinkage(residueString.substr((dashPosition - 1), 3));
        }
        std::string configuration = residueString.substr(dashPosition - 2, 1);
        if ((configuration == "a" || configuration == "b"))
        {
            this->SetConfiguration(residueString.substr(dashPosition - 2, 1));
        }
        // Find any special modifiers e.g. NAc, Gc, A in IdoA
        size_t modifierLength = (dashPosition - modifierStart - 2); // They are 2 apart if no modifier
        // logss << "modifierLength is " << modifierLength << ", dashPosition was " << dashPosition << ", ringPosition
        // was " << ringPosition << std::endl;
        if (modifierLength > 100)
        {
            gmml::log(__LINE__, __FILE__, gmml::WAR,
                      "This is a non-standard residue string that may cause issues: " + residueString);
        }
        if (modifierLength > 0 && modifierLength < 100)
        {
            this->SetResidueModifier(residueString.substr(modifierStart, modifierLength));
            this->ExciseRingShapeFromModifier();
        }
        else
        {
            this->SetResidueModifier("");
        }
    }
    else if (isdigit(residueString[0]))
    { // A derivative e.g. 3S, 6Me. Linkage followed by residue name. No configuration.
        this->SetType(ResidueType::Derivative);
        this->SetLinkage(residueString.substr(0, 1));
        this->SetResidueName(residueString.substr(1));                          // From position 1 to the end.
        if ((this->GetResidueName() == "D") || (this->GetResidueName() == "H")) // Supporting both for now.
        {
            this->SetType(ResidueType::Deoxy);
        }
    }
    else if (specifiedType == ResidueType::Aglycone)
    { // A terminal
        this->SetResidueName(residueString);
    }
    else
    { // Dunno.
        std::string message = "Error: we can't parse this residue: \"" + residueString + "\"";
        gmml::log(__LINE__, __FILE__, gmml::ERR, message);
        throw std::runtime_error(message);
    }
    // gmml::log(__LINE__, __FILE__, gmml::INF, this->PrintToString());
}

std::string ParsedResidue::PrintToString() const
{
    std::stringstream ss;
    ss << this->GetIsomer() << "_" << this->GetResidueName() << "_" << this->GetRingType() << "_"
       << this->GetResidueModifier() << "_" << this->GetRingShape() << "_" << this->GetConfiguration() << "_"
       << this->GetLinkage() << ". Type: " << this->GetType() << ".\n";
    return ss.str();
}

std::string ParsedResidue::GetMonosaccharideName() const
{
    return this->GetIsomer() + this->GetResidueName() + this->GetResidueModifier() + this->GetRingShape();
}

void ParsedResidue::ExciseRingShapeFromModifier()
{ // E.g. LIdopA(4C1)a1-4 with modifier "A(4C1)", which here gets broken into ring shape "4C1" and modifier "A".
    std::string modifier            = this->GetResidueModifier();
    size_t leftParenthesisPosition  = modifier.find('(');
    size_t rightParenthesisPosition = modifier.find(')');

    if ((leftParenthesisPosition == std::string::npos) || (rightParenthesisPosition == std::string::npos))
    { // If there isn't a ring shape declared.
        return;
    }
    else
    { // Assumes it's always at end of modifiers
        this->SetRingShape(modifier.substr(leftParenthesisPosition));
        this->SetResidueModifier(modifier.substr(0, leftParenthesisPosition));
    }
    return;
}
