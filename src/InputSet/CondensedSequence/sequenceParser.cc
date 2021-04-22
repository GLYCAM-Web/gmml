#include <string>
#include <algorithm> // Reverse function.

#include "includes/InputSet/CondensedSequence/sequenceParser.hpp"
#include "includes/MolecularModeling/Graph/Graph.hpp" // For print

using CondensedSequence::SequenceParser;
using CondensedSequence::ParsedResidue;
using TemplateGraph::Graph;

SequenceParser::SequenceParser (std::string inputSequence)
{
	if (inputSequence.find(';') != std::string::npos)
	{
		std::cout << "Found labels in input\n";
		//this->TokenizeLabelledInput(inputSequence);
	}
	else
	{
		std::cout << "Parsing unlabelled input sequence: " << inputSequence << "\n";
		this->ParseCondensedSequence(inputSequence);
	}
}

std::string SequenceParser::Print()
{
    std::string output;
    for (auto &residue : this->GetParsedResidues())
    {
        output += residue->Print();
    }
    return output;
}

ParsedResidue* SequenceParser::FindTerminalResidue()
{
    for (auto &residue : this->GetParsedResidues())
    {
        if (residue->GetType() == ParsedResidue::Type::Aglycone)
        {
            return residue;       
        }
    }
    //throw would be better? How do I handle Ano-Ano or cycles?
    return this->GetParsedResidues().at(0);
}

std::vector<ParsedResidue*> SequenceParser::GetParsedResidues()
{
    std::vector<ParsedResidue*> rawResidues;
    for(auto &residue : parsedResidues_)
    {
        rawResidues.push_back(residue.get());
    }
    return rawResidues;
}

std::vector<ParsedResidue*> SequenceParser::GetParsedResiduesOrderedByConnectivity()
{
    std::vector<ParsedResidue*> rawResidues;
    // Go via Graph so order decided by connectivity, depth first traversal:
    TemplateGraph::Graph<ParsedResidue> sequenceGraph(this->FindTerminalResidue());
    for(auto &node : sequenceGraph.GetNodes())
    {
        rawResidues.push_back(node->GetObjectPtr());
    }
    return rawResidues;
}


// void SequenceParser::TokenizeLabelledInput(std::string inString)
// {
// 	// char delimiter = ';';
// 	// std::vector<std::string> tokens = gmml::splitStringByDelimiter(inString, delimiter);
// 	// delimiter = '&';
// 	// for (auto &element : tokens)
// 	// {
// 	// 	std::vector<std::string> labelAndLabelee = gmml::splitStringByDelimiter(element, delimiter);
// 	// 	for (auto &subtoken : labelAndLabelee)
// 	// 	{
// 	// 		std::cout << " " << subtoken << " \n";
// 	// 	}
// 	// }
// }

bool SequenceParser::ParseCondensedSequence(const std::string sequence)
{
    // Reading from the end of the string.
	size_t i = (sequence.find_last_of('-') + 1);
    if (isdigit(sequence[i]))
    { // e.g. DGlcpa1-2DFrufb 
        ++i; // ano-ano
        parsedResidues_.push_back(std::make_unique<ParsedResidue>(sequence.substr(i), ParsedResidue::Type::Sugar));
    }
    else
    { // e.g. DGlcpa1-OH
        parsedResidues_.push_back(std::make_unique<ParsedResidue>(sequence.substr(i), ParsedResidue::Type::Aglycone));
    }
    auto terminal = parsedResidues_.back().get();
    this->RecurveParse(i, sequence, terminal); 
    return true;
}

void SequenceParser::RecurveParse(size_t &i, const std::string sequence, ParsedResidue* parent)
{
    bool branchStart = true;
    size_t windowEnd = i;
    while (i > 0)
    {
        i--;
        if ( sequence[i] == '-' )
        {
            if (branchStart)
            {
                branchStart = false;
            }
            else
            {
                parent = this->SaveResidue(i + 2, windowEnd, sequence, parent);
                windowEnd = i + 2; // Get to the right side of the number e.g. 1-4
            }
        }
        else if ( sequence[i] == '[' )
        {
            parent = this->SaveResidue(i + 1, windowEnd, sequence, parent);
            return;
        }
        else if ( sequence[i] == ']' )
        {
            if ((windowEnd - i) > 7) 
            {   // if not a derivative, save
                parent = this->SaveResidue(i + 1, windowEnd, sequence, parent);
                this->RecurveParse(i, sequence, parent);
                windowEnd = i;
                branchStart = true; // reset this when you fall out a level
            }
            else // Derivative
            {   
                this->RecurveParse(i, sequence, parent);
            }
        }
        else if ( sequence[i] == ',' )
        {
            parent = this->SaveResidue(i + 1, windowEnd, sequence, parent);
            windowEnd = i;
        }
        else if ( i == 0)
        {
            parent = this->SaveResidue(i, windowEnd, sequence, parent);
        }
    }
    return; 
}

ParsedResidue* SequenceParser::SaveResidue(const size_t windowStart, const size_t windowEnd, const std::string sequence, ParsedResidue* parent)
{
    std::string residueString = sequence.substr(windowStart, (windowEnd - windowStart));
    //std::cout << "At start of save: " << residueString << std::endl;
    // Splice out anything within [ and ].
    if (residueString.find('[') != std::string::npos) 
    {
        size_t branch_start = residueString.find_first_of('[');
        size_t branch_finish = residueString.find_last_of(']') + 1;
        std::string firstPart = residueString.substr(0, branch_start);
        std::string lastPart = residueString.substr(branch_finish);
        residueString = firstPart + lastPart;
        //std::cout << firstPart << " + " << lastPart << std::endl;
    }
    if (residueString.find('-') != std::string::npos)
    {
        std::cout << "Saving " << residueString << " with parent " << parent->GetName() <<  std::endl;
        parsedResidues_.push_back(std::make_unique<ParsedResidue>(residueString, parent));
        auto newRes = parsedResidues_.back().get();
        if(this->DerivativesExist())
        {
            for(auto &derivative : this->ExtractDerivatives())
            {
                std::cout << "Saving " << derivative << " with parent " << newRes->GetName() <<  std::endl;
                parsedResidues_.push_back(std::make_unique<ParsedResidue>(derivative, newRes));
            }
        }
        return newRes;
    }
    else // A derivatve. The parent residue doesn't exist yet, so save it.
    {
        //std::cout << "Temporarily holding derivative: " << residueString << "\n";
        this->SaveDerivative(residueString);
        return parent;
    }
}

std::vector<std::string> SequenceParser::ExtractDerivatives()
{
    auto derivativesCopy = savedDerivatives_;
    savedDerivatives_.clear();
    return derivativesCopy;
}


