#include <stdexcept>
#include <stack>
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceamberprepresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace gmml;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequence::CondensedSequence()
{
}

CondensedSequence::CondensedSequence(string sequence)
{
    residues_ = CondensedSequenceResidueVector();
    tokens_ = CondensedSequenceTokenTypeVector();
    condensed_sequence_residue_tree_ = CondensedSequenceResidueTree();
    ParseCondensedSequence(sequence);
    BuildArrayTreeOfCondensedSequenceResidue();
    BuildArrayTreeOfCondensedSequenceAmberPrepResidue(this->condensed_sequence_residue_tree_);
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
CondensedSequence::CondensedSequenceResidueVector CondensedSequence::GetResidues()
{
    return residues_;
}
CondensedSequence::CondensedSequenceTokenTypeVector CondensedSequence::GetTokens()
{
    return tokens_;
}
CondensedSequence::CondensedSequenceResidueTree CondensedSequence::GetCondensedSequenceResidueTree()
{
    return condensed_sequence_residue_tree_;
}
CondensedSequence::CondensedSequenceAmberPrepResidueTree CondensedSequence::GetCondensedSequenceAmberPrepResidueTree()
{
    return condensed_sequence_amber_prep_residue_tree_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequence::SetResidues(CondensedSequenceResidueVector residues)
{
    residues_.clear();
    for(CondensedSequenceResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}
void CondensedSequence::AddResidue(CondensedSequenceResidue* residue)
{
    residues_.push_back(residue);
}
void CondensedSequence::SetTokens(CondensedSequenceTokenTypeVector tokens)
{
    tokens_.clear();
    for(CondensedSequenceTokenTypeVector::iterator it = tokens.begin(); it != tokens.end(); it++)
    {
        tokens_.push_back(*it);
    }
}
void CondensedSequence::AddToken(CondensedSequenceTokenType token)
{
    tokens_.push_back(token);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
int CondensedSequence::InsertNodeInCondensedSequenceResidueTree(CondensedSequenceResidue *condensed_residue, int parent_node_id)
{
    if(parent_node_id != -1 && parent_node_id >= condensed_sequence_residue_tree_.size())
        throw std::invalid_argument("ArrayTree::insert - invalid parent index(" + gmml::ConvertT(parent_node_id) + ")");
//    condensed_sequence_residue_tree_.push_back(std::make_pair(condensed_residue, parent_node_id));
    condensed_residue->SetParentId(parent_node_id);
    condensed_sequence_residue_tree_.push_back(condensed_residue);
    return condensed_sequence_residue_tree_.size() - 1;
}

int CondensedSequence::InsertNodeInCondensedSequenceAmberPrepResidueTree(CondensedSequenceAmberPrepResidue *condensed_amber_prep_residue, int parent_node_id)
{
    if(parent_node_id != -1 && parent_node_id >= condensed_sequence_amber_prep_residue_tree_.size())
        throw std::invalid_argument("ArrayTree::insert - invalid parent index(" + gmml::ConvertT(parent_node_id) + ")");
//    condensed_sequence_amber_prep_residue_tree_.push_back(std::make_pair(condensed_tree_residue, parent_node_id));
    condensed_amber_prep_residue->SetParentId(parent_node_id);
    condensed_sequence_amber_prep_residue_tree_.push_back(condensed_amber_prep_residue);
    return condensed_sequence_amber_prep_residue_tree_.size() - 1;
}

void CondensedSequence::ParseCondensedSequence(string sequence)
{
    bool reading_residue = true;
    int start_index = 0;
    for(int i = 0; i < sequence.size(); i++)
    {
        switch (sequence[i])
        {
            case '-':
            {
                bool terminal = false;
                int end_index = i + 1;
                if(!std::isdigit(sequence[i+1]))
                {
                    terminal = true;
                    end_index--;
                }
                string residue = sequence.substr(start_index, end_index - start_index + 1);
                residues_.push_back(new CondensedSequenceResidue(residue));
                if(terminal)
                    start_index = i + 1;
                else
                {
                    start_index = i + 2;
                    i++;
                }
                tokens_.push_back(CONDENSED_SEQUENCE_RESIDUE);
                reading_residue = false;
                break;
            }
            case '[':
            {
                if(!reading_residue)
                {
                    tokens_.push_back(CONDENSED_SEQUENCE_LEFT_BRACKET);
                    start_index = i + 1;
                }
                break;
            }
            case ']':
            {
                if(!reading_residue)
                {
                    tokens_.push_back(CONDENSED_SEQUENCE_RIGHT_BRACKET);
                    start_index = i + 1;
                }
                break;
            }
            default:
            {
                reading_residue = true;
                break;
            }
        }
    }
    string terminal_residue = sequence.substr(sequence.find_last_of('-') + 1);
    if(terminal_residue.compare("") != 0)
    {
        residues_.push_back(new CondensedSequenceResidue(terminal_residue));
        tokens_.push_back(CONDENSED_SEQUENCE_RESIDUE);
    }
}

void CondensedSequence::BuildArrayTreeOfCondensedSequenceResidue()
{
    CondensedSequenceTokenTypeVector::reverse_iterator current_token = tokens_.rbegin();
    CondensedSequenceResidueVector::reverse_iterator current_residue = residues_.rbegin();

    if(residues_.size() == 0)
        return;

    stack<int> residue_stack;
    residue_stack.push(this->InsertNodeInCondensedSequenceResidueTree(*current_residue));
    while(++current_token != tokens_.rend())
    {
        switch(*current_token)
        {
            case CONDENSED_SEQUENCE_LEFT_BRACKET:
            {
                if (residue_stack.empty())
                    throw CondensedSequenceProcessingException("Invalid branching in sequence");
                residue_stack.pop();
                break;
            }
            case CONDENSED_SEQUENCE_RIGHT_BRACKET:
            {
                ++current_token;
                ++current_residue;
                if(current_residue == residues_.rend())
                    throw CondensedSequenceProcessingException("Invalid sequence of residues");
                if(residue_stack.empty())
                    throw CondensedSequenceProcessingException("Invalid sequence");
                residue_stack.push(this->InsertNodeInCondensedSequenceResidueTree(*current_residue, residue_stack.top()));
                break;
            }
            case CONDENSED_SEQUENCE_RESIDUE:
            {
                current_residue++;
                if(current_residue == residues_.rend())
                    throw CondensedSequenceProcessingException("Invalid sequence of residues");
                if(residue_stack.empty())
                    throw CondensedSequenceProcessingException("Invalid sequence");
                int parent = this->InsertNodeInCondensedSequenceResidueTree(*current_residue, residue_stack.top());
                residue_stack.pop();
                residue_stack.push(parent);
                break;
            }
        }
    }
}

void CondensedSequence::BuildArrayTreeOfCondensedSequenceAmberPrepResidue(CondensedSequenceResidueTree residue_tree)
{
    condensed_sequence_amber_prep_residue_tree_ = CondensedSequenceAmberPrepResidueTree();
    vector<vector<int> > open_valences = vector<vector<int> >(residue_tree.size());
    for(int i = 0; i < residue_tree.size(); i++)
    {
        int parent = residue_tree.at(i)->GetParentId();
        CondensedSequenceResidue* residue = residue_tree.at(i);
        if(parent != -1)
        {
            int oxygen_position = residue->GetOxygenPosition();
            open_valences[parent].push_back(oxygen_position);
        }
    }

    string terminal = residue_tree.at(0)->GetName();
    this->InsertNodeInCondensedSequenceAmberPrepResidueTree(new CondensedSequenceAmberPrepResidue(this->GetAmberPrepTerminalResidueCodeOfTerminalResidue(terminal)));

    int current_derivative_count = 0;
    vector<int> derivatives = vector<int>(residue_tree.size(), 0);
    for(int i = 1; i < residue_tree.size(); i++)
    {
        derivatives[i] = current_derivative_count;
        CondensedSequenceResidue* condensed_residue = residue_tree.at(i);
        int parent = residue_tree.at(i)->GetParentId();


        string parent_name = residue_tree.at(parent)->GetName();
        string anomeric_carbon = "C" + ConvertT<int>(condensed_residue->GetAnomericCarbon());
        string oxygen_position;
        (parent_name.compare("OME") == 0) ? oxygen_position = "O" : oxygen_position = "O" + ConvertT<int>(condensed_residue->GetOxygenPosition());

        try
        {
            CondensedSequenceAmberPrepResidue* tree_residue = new CondensedSequenceAmberPrepResidue(this->GetAmberPrepResidueCodeOfCondensedResidue(
                                                                                                        condensed_residue, open_valences[i], parent_name)
                                                                                                    , anomeric_carbon, oxygen_position);

            int residue_index = this->InsertNodeInCondensedSequenceAmberPrepResidueTree(tree_residue, parent + derivatives[parent]);

            CondensedSequenceResidue::DerivativeMap condensed_residue_derivatives = condensed_residue->GetDerivatives();
            for(CondensedSequenceResidue::DerivativeMap::iterator it = condensed_residue_derivatives.begin(); it != condensed_residue_derivatives.end(); ++it)
            {
                string derivative_name = it->second;
                int derivative_index = it->first;
                this->InsertNodeInCondensedSequenceAmberPrepResidueTree(this->GetCondensedSequenceDerivativeAmberPrepResidue(derivative_name, derivative_index), residue_index);
                current_derivative_count++;
            }
        }
        catch(exception ex)
        {
            CondensedSequenceAmberPrepResidue* tree_residue = new CondensedSequenceAmberPrepResidue(condensed_residue->GetName().substr(0,3)
                                                                                                    , anomeric_carbon, oxygen_position);

            this->InsertNodeInCondensedSequenceAmberPrepResidueTree(tree_residue, parent + derivatives[parent]);

            cout << "Invalid residue in the sequence (" << condensed_residue->GetName().substr(0,3) << ")" << endl;
        }
    }
}

string CondensedSequence::GetAmberPrepTerminalResidueCodeOfTerminalResidue(string terminal_residue_name)
{
    if(terminal_residue_name.compare("OH") == 0 || terminal_residue_name.compare("ROH") == 0)
        return "ROH";
    else if(terminal_residue_name.compare("OME") == 0)
        return "OME";
    else if(terminal_residue_name.compare("OtBu") == 0 || terminal_residue_name.compare("TBT") == 0)
        return "TBT";
    else if(AmberGlycamLookup(terminal_residue_name).amber_name_.compare("") != 0 ||
            AmberGlycamLookup(terminal_residue_name).glycam_name_.compare("") != 0)
        return AmberGlycamLookup(terminal_residue_name).glycam_name_;
    throw CondensedSequenceProcessingException("Invalid aglycon " + terminal_residue_name);
}

string CondensedSequence::GetAmberPrepResidueCodeOfCondensedResidue(CondensedSequenceResidue *condensed_residue, vector<int> open_valences, string parent_name)
{
    if(condensed_residue->GetName().compare("UNK") == 0 || condensed_residue->GetName().compare("Unknown") == 0)
        return "UNK";
    vector<int> new_valences_list = open_valences;
    CondensedSequenceResidue::DerivativeMap condensed_residue_derivatives = condensed_residue->GetDerivatives();
    for(CondensedSequenceResidue::DerivativeMap::iterator it = condensed_residue_derivatives.begin(); it != condensed_residue_derivatives.end(); ++it)
        new_valences_list.push_back(it->first);
    string residue_name = condensed_residue->GetName();
    string isomer = condensed_residue->GetIsomer();
    string configuration = condensed_residue->GetConfiguration();
    string ring_type = "P";
    if(residue_name.size() > 3)
    {
        char ring_letter = residue_name[3];
        if(ring_letter == 'p' || ring_letter == 'f')
        {
            residue_name.erase(residue_name.begin() + 3);
            if(ring_letter == 'f')
                ring_type = "F";
        }
    }

    bitset<10> open_valences_check = bitset<10>();
    for(int i = 0; i < open_valences.size(); i ++)
    {
        if(open_valences[i] < 0 || open_valences[i] >= 10)
            throw CondensedSequenceProcessingException("Invalid open valence");
        open_valences_check.set(open_valences[i]);
    }

    string residue_code = this->GetFirstLetterOfAmberPrepResidueCode(open_valences_check) + this->GetSecondLetterOfAmberPrepResidueCode(residue_name, isomer);
    if(residue_code.size() < 3)
        residue_code += this->GetThirdLetterOfAmberPrepResidueCode(configuration, ring_type);
    return residue_code;
}

string CondensedSequence::GetFirstLetterOfAmberPrepResidueCode(bitset<10> open_valences_check)
{
    bitset<10> open_valences_check_temp = open_valences_check;
    if (open_valences_check_temp.count() == 4)
    {
        if (open_valences_check_temp[4] && open_valences_check_temp[7] && open_valences_check_temp[8] && open_valences_check_temp[9])
            return "A";
        if (open_valences_check_temp[2] && open_valences_check_temp[3] && open_valences_check_temp[4] && open_valences_check_temp[6])
            return "P";
    }
    else if (open_valences_check_temp.count() == 3)
    {
        if (open_valences_check_temp[4] && open_valences_check_temp[7] && open_valences_check_temp[8])
            return "E";
        if (open_valences_check_temp[4] && open_valences_check_temp[7] && open_valences_check_temp[9])
            return "D";
        if (open_valences_check_temp[4] && open_valences_check_temp[8] && open_valences_check_temp[9])
            return "C";
        if (open_valences_check_temp[7] && open_valences_check_temp[8] && open_valences_check_temp[9])
            return "B";
        if (open_valences_check_temp[3] && open_valences_check_temp[4] && open_valences_check_temp[6])
            return "Q";
        if (open_valences_check_temp[2] && open_valences_check_temp[4] && open_valences_check_temp[6])
            return "R";
        if (open_valences_check_temp[2] && open_valences_check_temp[3] && open_valences_check_temp[6])
            return "S";
        if (open_valences_check_temp[2] && open_valences_check_temp[3] && open_valences_check_temp[4])
            return "T";
    }
    else if (open_valences_check_temp.count() == 2)
    {
        if (open_valences_check_temp[4] && open_valences_check_temp[7])
            return "K";
        if (open_valences_check_temp[4] && open_valences_check_temp[8])
            return "J";
        if (open_valences_check_temp[4] && open_valences_check_temp[9])
            return "I";
        if (open_valences_check_temp[7] && open_valences_check_temp[8])
            return "H";
        if (open_valences_check_temp[7] && open_valences_check_temp[9])
            return "G";
        if (open_valences_check_temp[8] && open_valences_check_temp[9])
            return "F";
        if (open_valences_check_temp[4] && open_valences_check_temp[6])
            return "U";
        if (open_valences_check_temp[3] && open_valences_check_temp[6])
            return "V";
        if (open_valences_check_temp[3] && open_valences_check_temp[4])
            return "W";
        if (open_valences_check_temp[2] && open_valences_check_temp[6])
            return "X";
        if (open_valences_check_temp[2] && open_valences_check_temp[4])
            return "Y";
        if (open_valences_check_temp[2] && open_valences_check_temp[3])
            return "Z";
    }
    else if (open_valences_check_temp.count() == 1)
    {
        for (int i = 1; i < open_valences_check_temp.size(); i++)
            if (open_valences_check_temp[i])
                return ConvertT<int>(i);
    }
    else if (open_valences_check_temp.none())
    {
        return "0";
    }

    throw CondensedSequenceProcessingException("There is no code in the GLYCAM code set for residues with open valences at the given positions.");
}

string CondensedSequence::GetSecondLetterOfAmberPrepResidueCode(string residue_name, string isomer)
{
    gmml::ResidueCodeName residue_name_code = gmml::ResidueNameCodeLookup(residue_name);
    if(residue_name_code.name_.compare("") != 0)
    {
        string code = residue_name_code.code_;
        if(isomer.compare("L") == 0)
            std::transform(code.begin(), code.end(), code.begin(), ::tolower);
        return code;
    }
    throw CondensedSequenceProcessingException(residue_name + " is not a valid residue");
}

string CondensedSequence::GetThirdLetterOfAmberPrepResidueCode(string configuration, string ring_type)
{
    if(configuration.compare("X") == 0)
        return "X";
    if(ring_type.compare("P") == 0)
        return (configuration.compare("A") == 0) ? "A" : "B";
    else
        return (configuration.compare("A") == 0) ? "D" : "U";
}

CondensedSequenceAmberPrepResidue* CondensedSequence::GetCondensedSequenceDerivativeAmberPrepResidue(string derivative_name, int derivative_index)
{
    string oxygen_name = "O" + ConvertT<int>(derivative_index);
    if(derivative_name.compare("S") == 0)
        return new CondensedSequenceAmberPrepResidue("SO3", "S1", oxygen_name, true);
    else if(derivative_name.compare("Me") == 0)
        return new CondensedSequenceAmberPrepResidue("MEX", "CH3", oxygen_name, true);
    else if(derivative_name.compare("A") == 0)
        return new CondensedSequenceAmberPrepResidue("ACX", "C1A", oxygen_name, true);
    throw CondensedSequenceProcessingException("There is no derivative in the GLYCAM code set represented by the letter " + derivative_name);
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequence::Print(ostream &out)
{}





