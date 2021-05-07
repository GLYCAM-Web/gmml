
#include <sstream>
#include <sys/stat.h> // for checking if file exists
#include "includes/InputSet/CondensedSequence/parsedResidue.hpp"
#include "includes/MolecularMetadata/GLYCAM/glycam06ResidueNameGenerator.hpp"

using CondensedSequence::ParsedResidue;
//using std::filesystem::exists;

bool file_exists (const char *filename)
{
    struct stat buffer;
    std::cout << filename << "\n";
    return (stat (filename, &buffer) == 0);
}

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
    if ( this->GetType() == Type::Sugar ) 
    {
        this->AddEdge(otherRes, this->GetConfiguration() + this->GetLinkage());
    }
    else
    {
        this->AddEdge(otherRes, this->GetLinkage());
    }
    return;
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

std::vector<ParsedResidue*> ParsedResidue::GetParents()
{
    return this->GetOutgoingNeighborObjects();
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
        if (dashPosition == std::string::npos) // There is no -. e.g. Fru in DGlcpa1-2DFrufb
        {
            dashPosition = residueString.size() + 1;
        }
        else
        {
            this->SetLinkage(residueString.substr((dashPosition - 1), 3 ));
        }
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
    	size_t modifierLength = (dashPosition - modifierStart - 2); // They are 2 apart if no modifier
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
    std::cout << this->Print();
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
        std::string code = gmml::MolecularMetadata::GLYCAM::Glycam06ResidueNameGenerator(linkages, this->GetIsomer(), this->GetResidueName(), 
                                                                            this->GetRingType(), this->GetResidueModifier(), this->GetConfiguration() );
        return code;
    }
    catch (const std::string exception)
    {
        std::cerr << "Error: " << exception << std::endl;
    }
    return "";
}

std::string ParsedResidue::GetGraphVizLine()
{
    std::string SnfgFilePath = "/home/oliver/Programs/GLYCAM_Dev_Env/V_2/Web_Programs/gems/gmml/includes/MolecularMetadata/Sugars/CFG_Symbol_Images/";
    std::cout << "Getting GraphVizLine for " << this->GetName() << "\n"; 
    std::stringstream ss;
    ss << this->GetIndex() << " [";
    // Aglycone
    if (this->GetType() == Type::Aglycone)
    {
        ss << "shape=box label=\"" << this->GetSimpleName() << "\"]";
        return ss.str();
    }
    // Sugar
    std::string label = "";
    std::string imageFile = SnfgFilePath + this->GetImageFileName() + ".svg";
    
    std::cout << "Searching for image: " << imageFile << "\n";
    if(file_exists(imageFile.c_str()))
    {
        std::cout << "FOUND IT\n";
        (this->GetRingType() == 'f') ? label = "f" : label = "";
        ss << "label=\"" << label << "\" height=\"0.7\" image=\"" << imageFile << "\"];\n";
    }
    else
    {
        ss << "shape=circle height=\"0.7\" label=\"" << this->GetSimpleName() << "\"];\n";
    }
    // Derivatives
    std::string derivativeStr = "";
    for (auto &childLink : this->GetChildren())
    {
        if (childLink->GetType() == Type::Derivative) 
        {
            derivativeStr += childLink->GetLinkageName() + childLink->GetName() + " ";
        }
    }
    if (! derivativeStr.empty())
    {
        ss << "\n" << "b" << this->GetIndex(); 
        ss << "[ shape=\"plaintext\",fontsize=\"12\",forcelabels=\"true\"; height = \"0.3\"; labelloc = b;  label=\""; 
        ss << derivativeStr << "\"]\n";
        ss << "{ rank=\"same\"; b" << this->GetIndex() << " " << this->GetIndex() << "}\n";
        ss << "{nodesep=\"0.02\";b" << this->GetIndex() << ";" << this->GetIndex() << "}\n";
    }
    // Linkage
    for (auto &parent : this->GetParents())
    { // There is either 1 or 0, this covers both cases 
        ss << this->GetIndex() << "--" << parent->GetIndex() << "[label=\"" << this->GetLinkageName() << "\" ];\n";
    }
    return ss.str();
}

std::string ParsedResidue::GetSimpleName()
{
    return this->GetResidueName() + this->GetResidueModifier();
}

std::string ParsedResidue::GetImageFileName()
{
    std::stringstream ss;
    ss << this->GetIsomer() << this->GetResidueName() << this->GetResidueModifier();
    return ss.str();
}


