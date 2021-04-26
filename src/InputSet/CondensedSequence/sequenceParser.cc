#include <string>
#include <algorithm> // Reverse function.

#include "includes/InputSet/CondensedSequence/sequenceParser.hpp"

using CondensedSequence::SequenceParser;
using CondensedSequence::ParsedResidue;

SequenceParser::SequenceParser (std::string inputSequence)
{
    
	if (inputSequence.find(';') != std::string::npos)
	{
		std::cout << "Found labels in input\n";
		this->ParseLabelledInput(inputSequence);
	}
	else
	{
		std::cout << "Parsing unlabelled input sequence:\n" << inputSequence << "\n";
        try 
        {
            if (this->CheckSequenceSanity(inputSequence))
            {
                std::cout << "Sequence is sane\n";
                this->ParseCondensedSequence(inputSequence);
            }
        }
        catch (std::string str)
        {
            std::cerr << str;
        }
	}
}

std::string SequenceParser::Print()
{
    std::string output = "";
    for (auto &residue : this->GetParsedResidues())
    {
        output += residue->Print();
    }
    return output;
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

void SequenceParser::ParseLabelledInput(std::string inString)
{
	// char delimiter = ';';
	// std::vector<std::string> tokens = gmml::splitStringByDelimiter(inString, delimiter);
	// delimiter = '&';
	// for (auto &element : tokens)
	// {
	// 	// std::vector<std::string> labelAndLabelee = gmml::splitStringByDelimiter(element, delimiter);
	// 	// for (auto &subtoken : labelAndLabelee)
	// 	// {
	// 		std::cout << " " << subtoken << " \n";
	// 	}
	// }
    std::cout << "We can't handle labeled stuff yet: " << inString << "\n";
}

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

bool SequenceParser::CheckSequenceSanity(std::string sequence)
{
    try
    {
        if ( sequence.empty() )
        {
            throw "Error: sequence is empty: ";
        }
        if ( sequence.find(" ") != std::string::npos )
        {
            throw "Error: sequence contains a space: ";
        }
        if ( sequence.find("cake") != std::string::npos )
        {
            throw "Error: the cake is a lie: ";
        }
        size_t a = std::count(sequence.begin(), sequence.end(), '[');
        size_t b = std::count(sequence.begin(), sequence.end(), ']');
        if (a != b)
        {
            throw "Error: the number of [ doesn't match the number of ]. Bad branch in : ";
        }
    }
    catch (char const* exception)
    {
        std::cerr << exception;
        std::cerr << "\n>>>" << sequence << "<<<\n";
        //throw message; 
        return false;
    }
    return true;
}


