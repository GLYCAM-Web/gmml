#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequenceResidue::CondensedSequenceResidue()
{
}

CondensedSequenceResidue::CondensedSequenceResidue(string residue_string)
{
    int dash_index = residue_string.find('-');
    if(dash_index == string::npos)
    {
        this->name_ = residue_string;
        this->is_terminal_ = true;
    }
    else if(residue_string.find("Unknown") != string::npos)
    {
        this->name_ = "UNK";
    }
    else
    {
        if(residue_string.empty())
            throw CondensedSequenceProcessingException("Invalid residue in sequence");

        char isomer_letter = residue_string[0];
        if(isomer_letter == 'D' || isomer_letter == 'd')
            this->isomer_ = "D";
        else if(isomer_letter == 'L' || isomer_letter == 'l')
            this->isomer_ = "L";
        else
            throw CondensedSequenceProcessingException("Invalid isomer in residue " + residue_string);

        if(dash_index <= 2)
            throw CondensedSequenceProcessingException("Invalid residue in sequence");

        char configuration_letter = residue_string[dash_index - 2];
        if(configuration_letter == 'A' || configuration_letter == 'a')
            this->configuration_ = 'A';
        else if(configuration_letter == 'B' || configuration_letter == 'b')
            this->configuration_ = 'B';
        else if(configuration_letter == 'X' || configuration_letter == 'x')
            this->configuration_ = 'X';
        else
            throw CondensedSequenceProcessingException("Invalid configuration in residue " + residue_string);

        if(!std::isdigit(residue_string[dash_index - 1]))
            throw CondensedSequenceProcessingException("Invalid residue in sequence");

        this->anomeric_carbon_ = gmml::ConvertString<int>(gmml::ConvertT<char>(residue_string[dash_index - 1]));
        if(dash_index != residue_string.size() - 1)
        {
            if(!std::isdigit(residue_string[dash_index + 1]))
                throw CondensedSequenceProcessingException("Invalid residue in sequence");
            this->oxygen_position_ = gmml::ConvertString<int>(gmml::ConvertT<char>(residue_string[dash_index + 1]));
        }
        else
            this->oxygen_position_ = 1;

        int left_bracket = residue_string.find('[');
        int right_bracket = residue_string.find(']');
        if(left_bracket == string::npos || right_bracket == string::npos)
            this->name_ = residue_string.substr(1, dash_index - 3);
        else
        {
            this->name_ = residue_string.substr(1, left_bracket - 1);
            string derivatives = residue_string.substr(left_bracket + 1, right_bracket - left_bracket - 1);
            vector<string> derivatives_tokens = gmml::Split(derivatives, ",");
            for(vector<string>::iterator it = derivatives_tokens.begin(); it != derivatives_tokens.end(); ++it)
            {
                if(!std::isdigit(it->at(0)))
                    throw CondensedSequenceProcessingException("Invalid derivative position in sequence");
                this->derivatives_[ConvertString<int>(ConvertT<char>(it->at(0)))] = it->substr(1);
            }
        }
    }
    parent_id_ = -1;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
bool CondensedSequenceResidue::GetIsTerminal()
{
    return is_terminal_;
}
string CondensedSequenceResidue::GetIsomer()
{
    return isomer_;
}
string CondensedSequenceResidue::GetConfiguration()
{
    return configuration_;
}
string CondensedSequenceResidue::GetName()
{
    return name_;
}
int CondensedSequenceResidue::GetAnomericCarbon()
{
    return anomeric_carbon_;
}
int CondensedSequenceResidue::GetOxygenPosition()
{
    return oxygen_position_;
}
CondensedSequenceResidue::DerivativeMap CondensedSequenceResidue::GetDerivatives()
{
    return derivatives_;
}
int CondensedSequenceResidue::GetParentId()
{
    return parent_id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequenceResidue::SetIsTerminal(bool is_terminal)
{
    is_terminal_ = is_terminal;
}
void CondensedSequenceResidue::SetIsomer(string isomer)
{
    isomer_ = isomer;
}
void CondensedSequenceResidue::SetConfiguration(string configuration)
{
    configuration_ = configuration;
}
void CondensedSequenceResidue::SetName(string name)
{
    name_ = name;
}
void CondensedSequenceResidue::SetAnomericCarbon(int anomeric_carbon)
{
    anomeric_carbon_ = anomeric_carbon;
}
void CondensedSequenceResidue::SetOxygenPosition(int oxygen_position)
{
    oxygen_position_ = oxygen_position;
}
void CondensedSequenceResidue::SetDerivatives(DerivativeMap derivatives)
{
    derivatives_.clear();
    for(DerivativeMap::iterator it = derivatives.begin(); it != derivatives.end(); it++)
    {
        int index = (*it).first;
        string derivative = (*it).second;
        derivatives_[index] = derivative;
    }
}
void CondensedSequenceResidue::SetParentId(int parent_id)
{
    parent_id_ = parent_id;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequenceResidue::Print(ostream &out)
{}





