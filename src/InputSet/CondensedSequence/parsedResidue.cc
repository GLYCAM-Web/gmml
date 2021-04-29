
#include <sstream>
#include "includes/InputSet/CondensedSequence/parsedResidue.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp"

using CondensedSequence::ParsedResidue;

ParsedResidue::ParsedResidue(std::string residueString, ParsedResidue::Type specifiedType) 
: Node(this, residueString), fullResidueString_ (residueString)  
{
    this->ParseResidueStringIntoComponents(residueString, specifiedType);
}

ParsedResidue::ParsedResidue(std::string residueString, ParsedResidue* neighbor, ParsedResidue::Type specifiedType) 
: Node(this, residueString), fullResidueString_ (residueString) 
{
    this->ParseResidueStringIntoComponents(residueString, specifiedType);
	this->AddLinkage(neighbor);
}

void ParsedResidue::AddLinkage(ParsedResidue* otherRes) 
{
    //std::string label = this->GetName() + "->" + otherRes->GetName();
    //std::cout << "Adding Edge: " << label << std::endl;
    this->AddEdge(otherRes, this->GetConfiguration() + this->GetLinkage());
   // std::cout << "Is it getting destroyed now?\n";
}

char ParsedResidue::GetLink()
{
    switch (this->GetType())
    {
        case (Type::Sugar):
            return this->GetLinkage().back();
        case (Type::Derivative):
            return this->GetLinkage().front();
        case (Type::Deoxy):
            return this->GetLinkage().front();
        default:
            return '0';
    }
}


std::vector<ParsedResidue*> ParsedResidue::GetChildren()
{
    return this->GetIncomingNeighborObjects();
}

std::string ParsedResidue::GetChildLinkages()
{
    std::string linkages;
    for (auto &child : this->GetChildren())
    {
        std::string link(1, child->GetLink()); // convert from char via string constructor.
        linkages += (link + ",");
    }
    if (linkages.empty())
    {
        linkages = "Terminal";
    }
    else     // Erase the last ","
    {
        linkages.erase(std::prev(linkages.end()));
    }
    return linkages;
}

std::string ParsedResidue::GetName(const bool withLabels)
{
    if (withLabels)
    {
        return FindLabelContaining("&Label=");
    }
    return this->GetIsomer() + this->GetResidueName() + this->GetRingType() + this->GetResidueModifier();
}

std::string ParsedResidue::GetLinkageName(const bool withLabels)
{   // Should only ever be zero or one outEdges in my current design.
    for (auto &linkage : this->GetOutEdges())
    {
        if (withLabels)
        {
            return linkage->FindLabelContaining("&Label=");
        }
        else
        {
            return linkage->GetLabel();
        }
    }    
    return ""; // aglycone/reducing terminal will not have linkage.
}

void ParsedResidue::ParseResidueStringIntoComponents(std::string residueString, ParsedResidue::Type specifiedType)
{
	//std::cout << "PARSING RESIDUE: " << residueString << std::endl;
    // Set defaults:
    this->SetIsomer('\0');
    this->SetResidueName("");
    this->SetRingType('\0');
    this->SetResidueModifier("");
    this->SetConfiguration('\0');
    this->SetLinkage("");
    this->SetType(specifiedType);
	if ( (residueString.find('-') != std::string::npos) || (specifiedType == Type::Sugar) )
    { // E.g. DManpNAca1-4 . Isomer (D or L), residueName (ManNAc), ring type (f or p), configuration (a or b), linkage (1-4)
    	// Reading from front.
        this->SetType(Type::Sugar);
        // Assumptions
        size_t residueStart = 1; // e.g. Gal, Glc, Ido
        size_t modifierStart = 5; // E.g. NAc, A, A(1C4)
        // Checks
        char isomer = residueString[0];
        if ((isomer == 'D') || (isomer == 'L')) 
        {
            this->SetIsomer(isomer);
        }
        else
        {
            residueStart--;
            modifierStart--;
        }
        this->SetResidueName(residueString.substr(residueStart, 3));
        size_t ringPosition = (residueStart + 3);
        char ringType = residueString[ringPosition];
        if (( ringType == 'p') || (ringType == 'f'))
        {
            this->SetRingType(ringType);
        }
        else
        {
            modifierStart--;
        }
        // Find the dash, read around it.
        size_t dashPosition = residueString.find('-');
        if (dashPosition == std::string::npos) // There is no -
        {
            dashPosition = residueString.size(); 
        }
        this->SetLinkage(residueString.substr((dashPosition - 1), 3 ));
        char configuration = residueString[dashPosition - 2];
        if (( configuration == 'a' || configuration == 'b'))
        {
            this->SetConfiguration(residueString[dashPosition - 2]);
        }
        else
        {
            modifierStart--;
        }
    	// Find any special modifiers e.g. NAc, Gc, A in IdoA
    	size_t modifierLength = (dashPosition - modifierStart - 2); // They are 3 apart if no modifier
        //std::cout << "modifierLength is " << modifierLength << ", dashPosition was " << dashPosition << ", ringPosition was " << ringPosition << std::endl;
        if (modifierLength > 100)
        {
            std::string message = "Non standard glycam residue string: " + residueString;
            //throw message;
            std::cout << message << std::endl;
        }
    	if (modifierLength > 0 && modifierLength < 100)
        {
    		this->SetResidueModifier(residueString.substr(modifierStart, modifierLength));
            //std::cout << "Modifier is " << this->GetResidueModifier() << std::endl;
        }
    	else
        {
    		this->SetResidueModifier("");
        }
    }
    else if ( isdigit(residueString[0]) )
    { // A derivative e.g. 3S, 6Me. Linkage followed by residue name. No configuration.
        this->SetType(Type::Derivative);
    	std::string linkage(1,residueString[0]);
    	this->SetLinkage(linkage);
    	this->SetResidueName(residueString.substr(1)); // From position 1 to the end.
        if (this->GetResidueName() == "D")
        {
            this->SetType(Type::Deoxy);
        }
    }
    else if (specifiedType == Type::Aglycone)
    { // A terminal
    	this->SetResidueName(residueString);
    }
    else
    { // Dunno.
        std::string message = "Error: we can't parse this residue: \"" + residueString + "\""; 
        throw message;
    }
    //std::cout << this->Print();
}

std::string ParsedResidue::Print()
{
	std::stringstream ss;
    ss << this->GetIsomer() << "_" 
				<< this->GetResidueName() << "_"
				<< this->GetRingType() << "_"
				<< this->GetResidueModifier() << "_"
				<< this->GetConfiguration() << "_"
				<< this->GetLinkage() << ".\n";
    return ss.str();
}

std::string ParsedResidue::GetGlycamResidueName()
{
    std::string linkages = "";
    if (this->GetType() == Type::Sugar)
    {
        linkages = this->GetChildLinkages();
    }
    try
    {
        auto code = gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNameGenerator(linkages, this->GetIsomer(), this->GetResidueName(), 
                                                                            this->GetRingType(), this->GetResidueModifier(), this->GetConfiguration() );
        return code;
    }
    catch (const std::string exception)
    {
        std::cerr << "Error: " << exception << std::endl;
    }
    return "";
}
