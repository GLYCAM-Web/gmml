#include "includes/CentralDataStructure/CondensedSequence/sequenceParser.hpp"
#include "includes/CodeUtils/logging.hpp"
#include <algorithm> // Reverse function.
#include <sstream>
#include <cctype> // isdigit()

using cdsCondensedSequence::ParsedResidue;
using cdsCondensedSequence::SequenceParser;

SequenceParser::SequenceParser(std::string inputSequence)
{
    if (inputSequence.find('<') != std::string::npos)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Found repeating unit in input\n");
        inputSequence = this->parseRepeatingUnits(inputSequence);
    }
    if (inputSequence.find(';') != std::string::npos)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Found labels in input\n");
        this->ParseLabelledInput(inputSequence);
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Parsing unlabelled input sequence:\n" + inputSequence + "\n");
        if (this->CheckSequenceSanity(inputSequence))
        {
            gmml::log(
                __LINE__, __FILE__, gmml::INF,
                "Sequence passed initial sanity checks for things like special characters or incorrect branching.\n");
            this->ParseCondensedSequence(inputSequence);
        }
    }
    gmml::log(__LINE__, __FILE__, gmml::INF, "SequenceParser constructor complete");
    return;
}

// Note Rob waved the wand and changed the format as of 2023-05-08. Changes below reflect that. Tails didn't change,
// Heads yes and have become optional. Examples: DGlcpa1-2[4DGlcpa1-]<4>OH becomes
// DGlcpa1-2DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH [4DGlcpa1-]<4>OH becomes DGlcpa1-4DGlcpa1-4DGlcpa1-4DGlcpa1-OH
// DGlcpa1-3[4DGlcpa1-]<9>2DManpa1-[4DGalpNAca1-]<4>OH // Multiple repeats
// DGlcpa1-2[4DGlcpa1-3DManpa1-]<9>OH // Disacc repeats
// DGlcpa1-2[4DGlcpa1-3[DAllpb1-2]DManpa1-]<9>OH // Branched repeats, unavailable on legacy
// Repeats within repeats?
// DGlcpa1-4[4[DGalpa1-3]DGlcpa1-3[DAllpb1-2]DManpa1-]<3>OH Repeats with branches on the leftmost Residue?
// DGlcpa1-[4DGlcpa1-3DManpa1-]<3>4DGalpa1-OH // Tails that aren't OH
std::string SequenceParser::parseRepeatingUnits(const std::string inputSequence)
{
    unsigned int repeatCharacterEndLocation   = inputSequence.find_first_of('>');
    unsigned int repeatCharacterStartLocation = inputSequence.find_first_of('<');
    if (repeatCharacterStartLocation == std::string::npos)
    {
        throw std::runtime_error("No '<' found in sequence with repeating syntax symbol '>' : " + inputSequence);
    }
    // Get the number of repeats:
    int numberRepeats = 0;
    try
    {
        unsigned int numberStart = repeatCharacterStartLocation + 1;
        std::string stringNumber = inputSequence.substr(numberStart, (repeatCharacterEndLocation - numberStart));
        numberRepeats            = std::stoi(stringNumber);
    }
    catch (...) // if e.g. stoi throws
    {
        throw std::runtime_error("Number of repeating units not specified correctly in repeating unit: " +
                                 inputSequence);
    }
    unsigned int i = repeatCharacterStartLocation;
    // Ensure next char is the ] of the repeating unit
    if (inputSequence[--i] != ']')
    {
        throw std::runtime_error("Missing or incorrect usage of ']' in repeating sequence: " + inputSequence);
    }
    // Ok now go find the position of the start of the repeating unit, considering branches
    unsigned int repeatEnd   = i;
    unsigned int repeatStart = this->seekRepeatStart(inputSequence, i);
    // std::cout << "Repeat starts at position " << repeatStart << " and ends here: " << repeatEnd << std::endl;
    //  std::cout << "Number of repeats: " << numberRepeats << "\n";
    std::string beforeRepeat = inputSequence.substr(0, repeatStart);
    // Check if using the old nomenclature. i.e. DGlcpa1-[4DGlcpa1-]<4> was old, DGlcpa1-4[4DGlcpa1-]<3> is new.
    if (inputSequence.find("-[") != std::string::npos)
    {
        std::string message = "Old repeat syntax detected as \"-[\" found in " + inputSequence;
        gmml::log(__LINE__, __FILE__, gmml::INF, message);
        //  I think I'll just copy the e.g. 4 to make it look like the new syntax.
        std::string numberToInsert = inputSequence.substr(repeatStart + 1, 1);
        gmml::log(__LINE__, __FILE__, gmml::INF,
                  "Adding this before repeat unit to make it be new syntax: " + numberToInsert);
        beforeRepeat += numberToInsert;
    }

    std::string firstRepeat = inputSequence.substr(
        repeatStart + 2, (repeatEnd - repeatStart - 2)); // firstRepeat does not have the e.g. 4 in 4DGlcpa1-
    std::string repeat = inputSequence.substr(repeatStart + 1, (repeatEnd - repeatStart - 1));
    std::string after  = inputSequence.substr(repeatCharacterEndLocation + 1);
    std::string newInputString;
    newInputString += beforeRepeat;
    newInputString += firstRepeat;
    numberRepeats--; // firstRepeat added already.
    for (int j = 1; j <= numberRepeats; ++j)
    {
        newInputString += repeat;
    }
    newInputString += after;
    // Check if there are more repeating units and deal with them recursively:
    if (after.find('>') != std::string::npos)
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Sequence with some repeats processed: " + newInputString);
        //        std::cout << "Sequence with some repeats processed: " << newInputString << "\n";
        newInputString = this->parseRepeatingUnits(newInputString);
    }
    else
    {
        gmml::log(__LINE__, __FILE__, gmml::INF, "Sequence with all repeats processed: " + newInputString);
        //        std::cout << "Sequence with all repeats processed: " << newInputString << "\n";
    }
    return newInputString;
}

// There may be branches which also use [], so need to check for those and find the [ that starts the repeat
unsigned int SequenceParser::seekRepeatStart(const std::string& inputSequence, unsigned int i)
{
    int branches = 0;

    while (i > 0)
    {
        --i; // skip the initial ]
        // std::cout << inputSequence[i];
        if (inputSequence[i] == ']')
        {
            ++branches;
        }
        if (inputSequence[i] == '[')
        {
            if (branches == 0)
            {
                return i;
            }
            else
            {
                --branches;
            }
        }
    }
    //    std::cout << "\n";
    throw std::runtime_error("Did not find corresponding '[' in repeat unit of repeating sequence: " + inputSequence);
}

std::string SequenceParser::Print() const
{
    std::string output = "";
    for (auto& residue : this->getResidues())
    {
        output += static_cast<ParsedResidue*>(residue)->PrintToString();
    }
    return output;
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
    throw "Error: SequenceParser can't handle labeled stuff yet: " + inString + "\n";
}

bool SequenceParser::ParseCondensedSequence(const std::string sequence)
{
    // Reading from the rightmost end of the string, get the aglycone first.
    size_t i = (sequence.find_last_of('-') + 1);
    if (isdigit(sequence[i]))   // Indicates ano-ano and not e.g. Sugar-OME or -ROH etc
    {                           // e.g. DGlcpa1-2DFrufb
        ++i;                    // ano-ano
        if (sequence[i] == ']') // e.g. DGlcpa1-2[LFucpa1-1]DFrufb also ano-ano, but with branching
        {
            ++i;
        }
        this->addResidue(std::make_unique<ParsedResidue>(sequence.substr(i), cds::ResidueType::Sugar));
    }
    else
    { // e.g. DGlcpa1-OH
        this->addResidue(std::make_unique<ParsedResidue>(sequence.substr(i), cds::ResidueType::Aglycone));
    }
    ParsedResidue* terminal = static_cast<ParsedResidue*>(this->getResidues().back());
    this->RecurveParseAlt(i, sequence, terminal);
    return true;
}

void SequenceParser::RecurveParseAlt(size_t& i, const std::string sequence, ParsedResidue* parent)
{
    // std::stringstream logss;
    // logss << "Started RecurveParseAlt with i=" << i << " and sequence[i] is " << sequence[i] << "\n";
    size_t windowEnd = i;
    while (i > 0)
    {
        i--;
        // logss << "i=" << i << " sequence[i]=" << sequence[i] << "\n";
        if ((sequence[i] == '-') && ((windowEnd - i) > 5))
        { // dash and have read enough that there is a residue to save.
          //  logss << "- with i=" << i << "\n";
            parent    = this->SaveResidue(i + 2, windowEnd, sequence, parent);
            windowEnd = i + 2; // Get to the right side of the number e.g. 1-4
        }
        else if (sequence[i] == ']')
        { // Start of branch: recurve. Maybe save if have read enough.
          //  logss << "Starting recurve down a branch as have found a ] with i=" << i << "\n";
            if ((windowEnd - i) > 5)
            { // if not a derivative start and have read enough, save.
                parent = this->SaveResidue(i + 1, windowEnd, sequence, parent);
                //     logss << "Recurving with branch\n";
                this->RecurveParseAlt(i, sequence, parent);
                //     logss << "Returned from branch\n";
                windowEnd = i;
            }
            else if ((windowEnd - i) > 1) // and not > 5
            {                             // Derivative. Note that windowEnd does not move.
                                          //    logss << "Recurving with derivative\n";
                this->RecurveParseAlt(i, sequence, parent);
                //    logss << "Returned from derivative recurve\n";
            }
            else
            { // Return from branch and find new one. e.g. [Gal][Gal]
                //    logss << "Recurving with ][ branch\n";
                this->RecurveParseAlt(i, sequence, parent);
                //    logss << "Returned from ][ branch\n";
                windowEnd = i;
            }
        }
        else if (sequence[i] == '[')
        { // End of branch
            // logss << "Saving and returning from recurve due to [ with i=" << i << "\n";
            this->SaveResidue(i + 1, windowEnd, sequence, parent);
            // gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
            return;
        }
        else if (sequence[i] == ',')
        { // , in derivative list
            // logss << ", so saving derivative with i=" << i << "\n";
            parent    = this->SaveResidue(i + 1, windowEnd, sequence, parent);
            windowEnd = i;
        }
        else if (i == 0)
        { // Fin.
            // logss << "i == 0 so need to save and then I'm finished!" << "\n";
            parent = this->SaveResidue(i, windowEnd, sequence, parent);
        }
    }
    // gmml::log(__LINE__, __FILE__, gmml::INF, logss.str());
    return;
}

ParsedResidue* SequenceParser::SaveResidue(const size_t windowStart, const size_t windowEnd, const std::string sequence,
                                           ParsedResidue* parent)
{
    std::string residueString = sequence.substr(windowStart, (windowEnd - windowStart));
    // gmml::log(__LINE__, __FILE__, gmml::INF, " Saving " + residueString + " with parent " + parent->GetName() );
    //  Splice out anything within [ and ].
    if (residueString.find('[') != std::string::npos)
    {
        size_t branch_start   = residueString.find_first_of('[');
        size_t branch_finish  = residueString.find_last_of(']') + 1;
        std::string firstPart = residueString.substr(0, branch_start);
        std::string lastPart  = residueString.substr(branch_finish);
        residueString         = firstPart + lastPart;
    }
    if (residueString.find('-') != std::string::npos)
    {
        this->addResidue(std::make_unique<ParsedResidue>(residueString, parent));
        ParsedResidue* newRes = static_cast<ParsedResidue*>(this->getResidues().back());
        if (this->DerivativesExist())
        {
            for (auto& derivative : this->ExtractDerivatives())
            {
                // gmml::log(__LINE__, __FILE__, gmml::INF, " Saving derivative " + derivative + " with parent " +
                // newRes->GetName() );
                this->addResidue(std::make_unique<ParsedResidue>(derivative, newRes));
            }
        }
        return newRes;
    }
    else // A derivatve. The parent residue doesn't exist yet, so save it.
    {
        // gmml::log(__LINE__, __FILE__, gmml::INF, "Temporarily holding derivative: " + residueString);
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
    if (sequence.empty())
    {
        throw "Error: sequence is empty:>>>" + sequence + "<<<";
    }
    if (sequence.find("cake") != std::string::npos)
    {
        throw "Error: the cake is a lie:>>>" + sequence + "<<<";
    }
    if (sequence.find(" ") != std::string::npos)
    {
        throw "Error: sequence contains a space:>>>" + sequence + "<<<";
    }
    std::vector<char> badChars = {'\'', '_', '+', '"', '`'};
    for (auto& badChar : badChars)
    {
        if (sequence.find(badChar) != std::string::npos)
        {
            std::string s(1, badChar); // convert to string
            throw "Error: sequence cannot contain this:\'" + s + "\':>>>" + sequence + "<<<";
        }
    }
    size_t a = std::count(sequence.begin(), sequence.end(), '[');
    size_t b = std::count(sequence.begin(), sequence.end(), ']');
    if (a != b)
    {
        throw "Error: the number of [ doesn't match the number of ]. Bad branch in :>>>" + sequence + "<<<";
    }
    return true;
}
