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
    for(unsigned int i = 0; i < sequence.size(); i++)
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
    for(unsigned int i = 0; i < residue_tree.size(); i++)
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
    for(unsigned int i = 1; i < residue_tree.size(); i++)
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
            throw CondensedSequenceProcessingException("Invalid residue in the sequence (" + condensed_residue->GetName().substr(0,3) + ")");
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
    for(unsigned int i = 0; i < open_valences.size(); i ++)
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
    string carbon_name = "C" + ConvertT<int>(derivative_index);;
    if(derivative_name.compare("S") == 0)
        return new CondensedSequenceAmberPrepResidue("SO3", carbon_name, oxygen_name, true);
    else if(derivative_name.compare("Me") == 0)
        return new CondensedSequenceAmberPrepResidue("MEX", carbon_name, oxygen_name, true);
    else if(derivative_name.compare("A") == 0)
        return new CondensedSequenceAmberPrepResidue("ACX", carbon_name, oxygen_name, true);
    throw CondensedSequenceProcessingException("There is no derivative in the GLYCAM code set represented by the letter " + derivative_name);
}

CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo CondensedSequence::GetCondensedSequenceRotamersAndGlycosidicAnglesInfo(CondensedSequenceResidueTree residue_tree)
{
    CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles = CondensedSequenceRotamersAndGlycosidicAnglesInfo();
    int linkage_index = 0;
    for(unsigned int i = 0; i < residue_tree.size(); i++)
    {
        int parent = residue_tree.at(i)->GetParentId();
        if(parent >= 0)
        {
            linkage_index++;
            CondensedSequenceResidue* residue = residue_tree.at(i);
            string residue_absolute_name = residue->GetName().substr(0, 3) + residue->GetName().substr(4);
            char ring_letter = residue->GetName()[3];
            vector<pair<string, vector<string> > > possible_rotamers = vector<pair<string, vector<string> > >();
            vector<pair<string, vector<string> > > selected_rotamers = vector<pair<string, vector<string> > >();
            vector<pair<string, double> > enabled_glycosidic_angles = vector<pair<string, double> >();
            enabled_glycosidic_angles.push_back(make_pair<string, double>("phi", dNotSet));
            enabled_glycosidic_angles.push_back(make_pair<string, double>("psi", dNotSet));
            if(ring_letter == 'p')
            {
                if(parent > 0)
                {

                    CondensedSequenceResidue* parent_residue = residue_tree.at(parent);
                    string parent_residue_absolute_name = parent_residue->GetName().substr(0,3) + parent_residue->GetName().substr(4);
                    stringstream rotamers_name;
                    rotamers_name << residue->GetIsomer() << residue->GetName() << residue->GetConfiguration() << residue->GetAnomericCarbon() << "-" <<
                                     residue->GetOxygenPosition() << parent_residue->GetIsomer() << parent_residue->GetName() << parent_residue->GetConfiguration();
                    switch(ResidueNameIndexLookup(residue_absolute_name).index_)
                    {
                        case 10:
                        case 1:
                        case 2:
                        case 12:
                        case 21:
                        case 20:
                            if(residue->GetOxygenPosition() == 6)
                            {
                                vector<string> rot = vector<string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                possible_rotamers.push_back(make_pair("omega", rot));

                                vector<string> rot1 = vector<string>();
                                rot1.push_back("gg");
                                rot1.push_back("gt");
                                selected_rotamers.push_back(make_pair("omega", rot1));

                                enabled_glycosidic_angles.push_back(make_pair<string, double>("omega", dNotSet));
                            }
                            break;
                        case 14:
                        case 13:
                        case 33:
                        case 9:
                        case 7:
                            if(residue->GetOxygenPosition() == 6)
                            {
                                vector<string> rot = vector<string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                possible_rotamers.push_back(make_pair("omega", rot));
                                selected_rotamers.push_back(make_pair("omega", rot));

                                enabled_glycosidic_angles.push_back(make_pair<string, double>("omega", dNotSet));
                            }
                            break;
                        case 23:
                        case 24:
                            switch(ResidueNameIndexLookup(parent_residue_absolute_name).index_)
                            {
                                case 10:
                                case 1:
                                case 2:
                                case 12:
                                case 21:
                                case 20:
                                    if(residue->GetOxygenPosition() == 6)
                                    {
                                        vector<string> rot = vector<string>();
                                        rot.push_back("gg");
                                        rot.push_back("gt");
                                        rot.push_back("tg");
                                        possible_rotamers.push_back(make_pair("omega", rot));

                                        vector<string> rot1 = vector<string>();
                                        rot1.push_back("gg");
                                        rot1.push_back("gt");
                                        selected_rotamers.push_back(make_pair("omega", rot1));

                                        enabled_glycosidic_angles.push_back(make_pair<string, double>("omega", dNotSet));
                                    }
                                    break;
                                case 14:
                                case 13:
                                case 33:
                                case 9:
                                case 7:
                                    if(residue->GetOxygenPosition() == 6)
                                    {
                                        vector<string> rot = vector<string>();
                                        rot.push_back("gg");
                                        rot.push_back("gt");
                                        rot.push_back("tg");
                                        possible_rotamers.push_back(make_pair("omega", rot));
                                        selected_rotamers.push_back(make_pair("omega", rot));

                                        enabled_glycosidic_angles.push_back(make_pair<string, double>("omega", dNotSet));
                                    }
                                    break;
                            }
                            if(ResidueNameIndexLookup(parent_residue_absolute_name).index_ != 23 &&
                                    ResidueNameIndexLookup(parent_residue_absolute_name).index_ != 24)
                            {
                                vector<string> rot = vector<string>();
                                rot.push_back("t");
                                rot.push_back("g");
                                rot.push_back("-g");
                                possible_rotamers.push_back(make_pair("phi", rot));

                                vector<string> rot1 = vector<string>();
                                rot1.push_back("t");
                                rot1.push_back("-g");
                                selected_rotamers.push_back(make_pair("phi", rot1));
                            }
                            break;
                    }
                    RotamersAndGlycosidicAnglesInfo* info = new RotamersAndGlycosidicAnglesInfo(linkage_index, possible_rotamers, selected_rotamers, enabled_glycosidic_angles);
                    RotamerNameInfoPair pair_info = make_pair(rotamers_name.str(), info);

                    rotamers_glycosidic_angles.push_back(pair_info);
                }

                vector<pair<string, vector<string> > > der_possible_rotamers = vector<pair<string, vector<string> > >();
                vector<pair<string, vector<string> > > der_selected_rotamers = vector<pair<string, vector<string> > >();
                vector<pair<string, double> > der_enabled_glycosidic_angles = vector<pair<string, double> >();
                der_enabled_glycosidic_angles.push_back(make_pair<string, double>("phi", dNotSet));
                der_enabled_glycosidic_angles.push_back(make_pair<string, double>("psi", dNotSet));
                CondensedSequenceResidue::DerivativeMap derivatives = residue->GetDerivatives();
                for(CondensedSequenceResidue::DerivativeMap::iterator it = derivatives.begin(); it != derivatives.end(); it++)
                {
                    int derivative_index = (*it).first;
                    string derivative_name = (*it).second;
                    switch(ResidueNameIndexLookup(residue_absolute_name).index_)
                    {
                        case 10:
                        case 1:
                        case 2:
                        case 12:
                        case 21:
                        case 20:
                            if(derivative_index == 6 && derivative_name.compare("S") == 0)
                            {
                                vector<string> rot = vector<string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                der_possible_rotamers.push_back(make_pair("omega", rot));

                                vector<string> rot1 = vector<string>();
                                rot1.push_back("gg");
                                rot1.push_back("gt");
                                der_selected_rotamers.push_back(make_pair("omega", rot1));

                                der_enabled_glycosidic_angles.push_back(make_pair<string, double>("omega", dNotSet));
                            }
                            break;
                        case 14:
                        case 13:
                        case 33:
                        case 9:
                        case 7:
                            if(derivative_index == 6 && derivative_name.compare("S") == 0)
                            {
                                vector<string> rot = vector<string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                der_possible_rotamers.push_back(make_pair("omega", rot));
                                der_selected_rotamers.push_back(make_pair("omega", rot));

                                der_enabled_glycosidic_angles.push_back(make_pair<string, double>("omega", dNotSet));
                            }
                            break;
                    }
                    linkage_index++;
                    stringstream der_rotamers_name;
                    der_rotamers_name << residue->GetIsomer() << residue->GetName() << residue->GetConfiguration()
                                      << "[" << derivative_index << derivative_name << "]";
                    RotamersAndGlycosidicAnglesInfo* der_info = new RotamersAndGlycosidicAnglesInfo(linkage_index, der_possible_rotamers,
                                                                                                    der_selected_rotamers, der_enabled_glycosidic_angles);
                    RotamerNameInfoPair der_pair_info = make_pair(der_rotamers_name.str(), der_info);

                    rotamers_glycosidic_angles.push_back(der_pair_info);
                }
            }
        }
    }
    return rotamers_glycosidic_angles;
}

int CondensedSequence::CountAllPossibleSelectedRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info)
{
    int count = 1;
    for(unsigned int i = 0; i < rotamers_glycosidic_angles_info.size(); i++)
    {
        vector<pair<string, double> > angles = rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_;
        vector<pair<string, vector<string> > > selected_rotamers = rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_;
        int omega_rotamers = 0;
        int phi_rotamers = 0;
        for(unsigned int j = 0; j < selected_rotamers.size(); j++)
        {
            pair<string, vector<string> > rot = selected_rotamers.at(j);
            if(rot.first.compare("phi") == 0)
                phi_rotamers+= rot.second.size();
            if(rot.first.compare("omega") == 0)
                omega_rotamers += rot.second.size();
        }
        if(angles.at(0).second == dNotSet) // Phi
            count *= (phi_rotamers == 0) ? 1 : phi_rotamers;
        if(angles.size() == 3)      // Phi, Psi, Omega
        {
            if(angles.at(2).second == dNotSet)
                count *= (omega_rotamers == 0) ? 1 : omega_rotamers;
        }
    }

    return count;
}

int CondensedSequence::CountAllPossible28LinkagesRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info)
{
    int count = 1;
    for(unsigned int i = 0; i < rotamers_glycosidic_angles_info.size(); i++)
    {
        string rotamer_name = rotamers_glycosidic_angles_info.at(i).first;
        if(rotamer_name.find("DNeup5AcA2-8") != string::npos ||
                rotamer_name.find("DNeup5AcB2-8") != string::npos)
            count *= 6;
        if(rotamer_name.find("DNeup5GcA2-8") != string::npos)
            count *= 5;
    }
    return count;
}

vector<vector<int> > CondensedSequence::CreateBaseMapAllPossibleSelectedRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info)
{
    vector<vector<int> > mapper = vector<vector<int> >();
    for(unsigned int i = 0; i < rotamers_glycosidic_angles_info.size(); i++)
    {
        vector<pair<string, double> > angles = rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_;
        vector<pair<string, vector<string> > > selected_rotamers = rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_;
        int omega_rotamers = 0;
        int psi_rotamers = 1;
        int phi_rotamers = 0;
        for(unsigned int j = 0; j < selected_rotamers.size(); j++)
        {
            pair<string, vector<string> > rot = selected_rotamers.at(j);
            if(rot.first.compare("phi") == 0)
                phi_rotamers += rot.second.size();
            if(rot.first.compare("omega") == 0)
                omega_rotamers += rot.second.size();
        }
        if(angles.at(0).second == dNotSet) // Phi
            phi_rotamers = (phi_rotamers == 0) ? 1 : phi_rotamers;
        if(angles.size() == 3)      // Phi, Psi, Omega
        {
            if(angles.at(2).second == dNotSet)
                omega_rotamers = (omega_rotamers == 0) ? 1 : omega_rotamers;
        }

        vector<int> val = vector<int>();
        val.push_back(phi_rotamers);
        val.push_back(psi_rotamers);
        val.push_back(omega_rotamers);

        mapper.push_back(val);
    }

    return mapper;
}

CondensedSequence::IndexLinkageConfigurationMap CondensedSequence::CreateIndexLinkageConfigurationMap(
        CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info, IndexNameMap& names)
{
    IndexLinkageConfigurationMap mapper = IndexLinkageConfigurationMap();
    IndexConfigurationNameMap name_mapper = IndexConfigurationNameMap();
    vector<vector<int> > mapping = this->CreateBaseMapAllPossibleSelectedRotamers(rotamers_glycosidic_angles_info);
    vector<vector<vector<double> > > phi_psi_omega_vector_map = vector<vector<vector<double> > >();
    vector<vector<vector<string> > > phi_psi_omega_rt_vector_map = vector<vector<vector<string> > >();
    for(unsigned int i = 0; i < mapping.size(); i++)
    {
        vector<double> phi = vector<double>();
        vector<string> phi_rt = vector<string>();
        // phi => dNotSet in all phi cases means default value
        if(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(0).second == dNotSet) // value not set
        {
            if(rotamers_glycosidic_angles_info.at(i).second->possible_rotamers_.size() == 0) // no possible rotamer
            {
                phi.push_back(dNotSet);
                phi_rt.push_back("df");
            }
            else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.size() == 0) // no selected rotamer
            {
                phi.push_back(dNotSet);
                phi_rt.push_back("df");
            }
            else
            {
                bool phi_check = false;
                for(unsigned int j = 0; j < rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.size(); j++)
                {
                    if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).first.compare("phi") == 0)
                    {
                        if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.size() == 0)
                        {
                            phi.push_back(dNotSet);
                            phi_rt.push_back("df");
                        }
                        else
                        {
                            for(unsigned int k = 0; k < rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.size(); k++)
                            {
                                if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.at(k).compare("t") == 0)
                                {
                                    phi.push_back(180.0);
                                    phi_rt.push_back("t");
                                }
                                else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.at(k).compare("g") == 0)
                                {
                                    phi.push_back(60.0);
                                    phi_rt.push_back("g");
                                }
                                else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.at(k).compare("-g") == 0)
                                {
                                    phi.push_back(-60.0);
                                    phi_rt.push_back("-g");
                                }
                                else
                                {
                                    phi.push_back(dNotSet);
                                    phi_rt.push_back("df");
                                }
                            }
                        }
                        phi_check = true;
                    }
                }
                if(!phi_check)
                {
                    phi.push_back(dNotSet);
                    phi_rt.push_back("df");
                }
            }
        }
        else
        {
            phi.push_back(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(0).second);
            phi_rt.push_back("cu");
        }

        vector<double> psi = vector<double>();
        vector<string> psi_rt = vector<string>();
        if(rotamers_glycosidic_angles_info.at(i).first.find("[") == string::npos)
        {
            psi.push_back(dNotSet); // standard psi angle
            psi_rt.push_back("df");
        }
        else
        {
            psi.push_back(dNotSet);
            psi_rt.push_back("df");
        }

        // omega => dNotSet in all omega cases means not applicable
        // for the cases that the angle is not set by any of selecting rotamers or setting the angle, it is set to 180.0
        vector<double> omega = vector<double>();
        vector<string> omega_rt = vector<string>();
        if(mapping.at(i).at(2) == 0)
        {
            omega.push_back(dNotSet);
            omega_rt.push_back("na");
        }
        else
        {
            if(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(2).second == dNotSet) // value not set
            {
                if(rotamers_glycosidic_angles_info.at(i).second->possible_rotamers_.size() == 0) // no possible rotamer, never happens
                {
                    omega.push_back(180.0);
                    omega_rt.push_back("df");
                }
                else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.size() == 0) // no selected rotamer
                {
                    omega.push_back(180.0);
                    omega_rt.push_back("df");
                }
                else
                {
                    bool omega_check = false;
                    for(unsigned int j = 0; j < rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.size(); j++)
                    {
                        if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).first.compare("omega") == 0)
                        {
                            if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.size() == 0)
                            {
                                omega.push_back(180.0);
                                omega_rt.push_back("df");
                            }
                            else
                            {
                                for(unsigned int k = 0; k < rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.size(); k++)
                                {
                                    if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.at(k).compare("tg") == 0)
                                    {
                                        omega.push_back(180.0);
                                        omega_rt.push_back("tg");
                                    }
                                    else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.at(k).compare("gt") == 0)
                                    {
                                        omega.push_back(60.0);
                                        omega_rt.push_back("gt");
                                    }
                                    else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.at(j).second.at(k).compare("gg") == 0)
                                    {
                                        omega.push_back(-60.0);
                                        omega_rt.push_back("gg");
                                    }
                                    else
                                    {
                                        omega.push_back(180.0);
                                        omega_rt.push_back("df");
                                    }
                                }
                            }
                            omega_check = true;
                        }
                    }
                    if(!omega_check)
                    {
                        omega.push_back(180.0);
                        omega_rt.push_back("df");
                    }
                }
            }
            else
            {
                omega.push_back(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(0).second);
                omega_rt.push_back("cu");
            }
        }

        // Build all <phi, psi, omega> for each linkage
        vector<vector<double> > phi_psi_omega = vector<vector<double> >();
        vector<vector<string> > phi_psi_omega_rt = vector<vector<string> >();
        for(unsigned int j = 0; j < phi.size(); j++)
        {
            for(unsigned int k = 0; k < psi.size(); k++)
            {
                for(unsigned int l = 0; l < omega.size(); l++)
                {
                    vector<double> combination = vector<double>();
                    vector<string> combination_rt = vector<string>();
                    combination.push_back(phi.at(j));
                    combination_rt.push_back(phi_rt.at(j));
                    combination.push_back(psi.at(k));
                    combination_rt.push_back(psi_rt.at(k));
                    combination.push_back(omega.at(l));
                    combination_rt.push_back(omega_rt.at(l));
                    string rotamer_name = rotamers_glycosidic_angles_info.at(i).first;
                    if(rotamer_name.find("DNeup5AcA2-8") != string::npos ||
                            rotamer_name.find("DNeup5AcB2-8") != string::npos)
                    {
                        for(int m = 0; m < 6; m++)
                        {
                            vector<double> new_combination = combination;
                            vector<string> new_combination_rt = combination_rt;
                            for(int n = 0; n < 5; n++)
                            {
                                new_combination.push_back(EXTERNAL28LINKAGEROTAMERS[m][n]);
                                new_combination_rt.push_back("E" + ConvertT<int>(m) + "/" + ConvertT<int>(n));
                            }
                            phi_psi_omega.push_back(new_combination);
                            phi_psi_omega_rt.push_back(new_combination_rt);
                        }
                    }
                    else if(rotamer_name.find("DNeup5GcA2-8") != string::npos)
                    {
                        for(int m = 0; m < 5; m++)
                        {
                            vector<double> new_combination = combination;
                            vector<string> new_combination_rt = combination_rt;
                            for(int n = 0; n < 5; n++)
                            {
                                new_combination.push_back(INTERNAL28LINKAGEROTAMERS[m][n]);
                                new_combination_rt.push_back("I" + ConvertT<int>(m) + "/" + ConvertT<int>(n));
                            }
                            phi_psi_omega.push_back(new_combination);
                            phi_psi_omega_rt.push_back(new_combination_rt);
                        }
                    }
                    else
                    {
                        phi_psi_omega.push_back(combination);
                        phi_psi_omega_rt.push_back(combination_rt);
                    }
                }
            }
        }
        phi_psi_omega_vector_map.push_back(phi_psi_omega);
        phi_psi_omega_rt_vector_map.push_back(phi_psi_omega_rt);

    }

    // Print all rotamers for each linkage
    /*for(unsigned int i = 0; i < phi_psi_omega_vector_map.size(); i++)
    {
        vector<vector<double> > phi_psi_omega_combination = phi_psi_omega_vector_map.at(i);
        for(unsigned int j = 0; j < phi_psi_omega_combination.size(); j++)
        {
            vector<double> phi_psi_omega = phi_psi_omega_combination.at(j);
            cout << "<";
            for(unsigned int k = 0; k < phi_psi_omega.size(); k++)
            {
                cout << phi_psi_omega.at(k);
                (k < phi_psi_omega.size() - 1) ? cout << ", " : cout << "";
            }
            cout << ">";
        }
        cout << endl;
    }*/

    // Build all rotamers with combination of each phi_psi_omega value for each linkage with the others
    for(unsigned int j = 0; j < phi_psi_omega_vector_map.size(); j++)
    {
        IndexLinkageConfigurationMap new_mapper = IndexLinkageConfigurationMap();
        IndexConfigurationNameMap new_name_mapper = IndexConfigurationNameMap();
        new_mapper = mapper;
        new_name_mapper = name_mapper;
        int counter = 0;
        for(unsigned int k = 0; k < phi_psi_omega_vector_map.at(j).size(); k++)
        {
            if(j == 0)
            {
                mapper[k] = vector<vector<double> >();
                mapper[k].push_back(phi_psi_omega_vector_map.at(j).at(k));
                name_mapper[k].push_back(phi_psi_omega_rt_vector_map.at(j).at(k));
            }
            else
            {
                if(k == 0)
                    for(unsigned int i = 0; i < mapper.size(); i++)
                    {
                        mapper[i].push_back(phi_psi_omega_vector_map.at(j).at(k));
                        name_mapper[i].push_back(phi_psi_omega_rt_vector_map.at(j).at(k));
                    }
                else
                {
                    for(unsigned int i = 0; i < new_mapper.size(); i++)
                    {
                        vector<vector<double> > res = new_mapper[i];
                        vector<vector<string> > res_rt = new_name_mapper[i];
                        res.push_back(phi_psi_omega_vector_map.at(j).at(k));
                        res_rt.push_back(phi_psi_omega_rt_vector_map.at(j).at(k));
                        mapper[new_mapper.size() + counter] = res;
                        name_mapper[new_name_mapper.size() + counter] = res_rt;
                        counter++;
                    }
                }
            }
        }
    }

    for(int i = 0; i < mapper.size(); i++)
    {
        /*cout << i << ": ";
        for(unsigned int j = 0; j < mapper[i].size(); j++)
        {
            cout << "<";
            for(unsigned int k = 0; k < mapper[i].at(j).size(); k++)
            {
                cout << mapper[i].at(j).at(k);
                (k < mapper[i].at(j).size() - 1) ? cout << ", " : cout << " ";
            }
            cout << ">";            
        }*/
        stringstream ss;
        for(unsigned int j = 0; j < mapper[i].size(); j++)
        {
            for(unsigned int k = 0; k < mapper[i].at(j).size(); k++)
            {
                if(name_mapper[i].at(j).at(k).compare("df") != 0 &&
                        name_mapper[i].at(j).at(k).compare("cu") != 0 &&
                        name_mapper[i].at(j).at(k).compare("na") != 0)
                    ss << name_mapper[i].at(j).at(k) << "-";
            }
        }
        names[i] = ss.str().substr(0,ss.str().size()-1);
        /*cout << "<" << ss.str().substr(0,ss.str().size()-1) << ">";
        cout << endl;*/
    }
    return mapper;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequence::Print(ostream &out)
{

}





