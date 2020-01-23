#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/InputSet/Utilities/response.hpp"
#include "../../../includes/Glycan/note.hpp"

using CondensedSequenceSpace::CondensedSequenceResidue;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequenceResidue::CondensedSequenceResidue()
{
}

CondensedSequenceResidue::CondensedSequenceResidue(std::string residue_string, CondensedSequenceSpace::CondensedSequence* condensed_sequence)
{
    size_t dash_index = residue_string.find('-');
    if(dash_index == std::string::npos) //If there is no dashes, then this must be a terminal residue, which is either an aglycone or a sugar bearing ano-ano linkage. 
    {
	//If the first character of residue_string is not a number, for instance "OH", it is an aglycone
	if (!std::isdigit(residue_string[0])){
            this->name_ = residue_string;
            this->is_terminal_aglycone_ = true;
        }
	//Otherwise, it is a ano-ano sugar, for instance "1]DFrufb".
	else{
            //In this case, see if the 2nd character is ]. If so, ignore this character.3rd becomes isomer, 4th to (n-1)th becomes residue name, last char becomes configuration. 	    
	    if (residue_string.substr(1,1) == "]"){
	        char isomer_letter = residue_string[2];
		if (isomer_letter == 'D' || isomer_letter == 'd'){
		    this->isomer_ = "D";
		}
		else if (isomer_letter == 'L' || isomer_letter == 'l'){
		    this->isomer_ = "L";
		}
                else{
                    std::stringstream error_notice;
                    error_notice << "ERROR at residue " << residue_string << ": invalid isomer, must be either D/L";
                    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                    throw CondensedSequenceProcessingException("Invalid isomer in residue " + residue_string);
                }

		this->name_ = residue_string.substr(3, 4);

		char configuration_letter = residue_string[residue_string.size()-1];
		//std::cout << "config letter is: " << configuration_letter << std::endl;
                if(configuration_letter == 'A' || configuration_letter == 'a')
                    this->configuration_ = 'A';
                else if(configuration_letter == 'B' || configuration_letter == 'b')
                    this->configuration_ = 'B';
                else if(configuration_letter == 'X' || configuration_letter == 'x')
                    this->configuration_ = 'X';
                else{
                    std::stringstream error_notice;
                    error_notice << "ERROR at residue " << residue_string << ": invalid anomeric configuration, must be either a/b/x";
                    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                    //throw CondensedSequenceProcessingException("Invalid configuration in residue " + residue_string);
                }
	    }
	    else{ //If residue_string does not contain "]",for example 1DFrufb
	        char isomer_letter = residue_string[1];
		if (isomer_letter == 'D' || isomer_letter == 'd'){
                    this->isomer_ = "D";
                }
                else if (isomer_letter == 'L' || isomer_letter == 'l'){
                    this->isomer_ = "L";
                }
                else{
                    std::stringstream error_notice;
                    error_notice << "ERROR at residue " << residue_string << ": invalid isomer, must be either D/L";
                    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                    throw CondensedSequenceProcessingException("Invalid isomer in residue " + residue_string);
                }

		this->name_ = residue_string.substr(2, 4);
                char configuration_letter = residue_string[residue_string.size()-1];
                //std::cout << "config letter is: " << configuration_letter << std::endl;
                if(configuration_letter == 'A' || configuration_letter == 'a')
                    this->configuration_ = 'A';
                else if(configuration_letter == 'B' || configuration_letter == 'b')
                    this->configuration_ = 'B';
                else if(configuration_letter == 'X' || configuration_letter == 'x')
                    this->configuration_ = 'X';
                else{
                    std::stringstream error_notice;
                    error_notice << "ERROR at residue " << residue_string << ": invalid anomeric configuration, must be either a/b/x";
                    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                    //throw CondensedSequenceProcessingException("Invalid configuration in residue " + residue_string);
                }
	    }
	    this->is_terminal_sugar_ = true;
	}
    }
    else
    {
        if(residue_string.empty()){

	    std::string error_notice = "Empty residue token detected, but can't determine its exact location in sequence string";
	    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice));
            throw CondensedSequenceProcessingException("Invalid residue in sequence");
	}
	
	else{
            if(residue_string.find("Unknown") != std::string::npos)
                this->name_ = "UNK";
            else
            {
                char isomer_letter = residue_string[0];
                if(isomer_letter == 'D' || isomer_letter == 'd')
                    this->isomer_ = "D";
                else if(isomer_letter == 'L' || isomer_letter == 'l')
                    this->isomer_ = "L";
                else{
		    std::stringstream error_notice;
		    error_notice << "ERROR at residue " << residue_string << ": invalid isomer, must be either D/L";
		    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));

                    //throw CondensedSequenceProcessingException("Invalid isomer in residue " + residue_string);
	        }

                if(dash_index <= 2)
                    throw CondensedSequenceProcessingException("Invalid residue in sequence");

                char configuration_letter = residue_string[dash_index - 2];
                if(configuration_letter == 'A' || configuration_letter == 'a')
                    this->configuration_ = 'A';
                else if(configuration_letter == 'B' || configuration_letter == 'b')
                    this->configuration_ = 'B';
                else if(configuration_letter == 'X' || configuration_letter == 'x')
                    this->configuration_ = 'X';
                else{
		    std::stringstream error_notice;
		    error_notice << "ERROR at residue " << residue_string << ": invalid anomeric configuration, must be either a/b/x";
		    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                    //throw CondensedSequenceProcessingException("Invalid configuration in residue " + residue_string);
	        }

                if(!std::isdigit(residue_string[dash_index - 1])){

		    std::stringstream error_notice;
		    error_notice << "ERROR at residue " << residue_string << ": non-numeric anomeric position, usually 1 for aldose and 2 for ketose";
		    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                    //throw CondensedSequenceProcessingException("Invalid residue in sequence");
	        }

	        else{
                    this->anomeric_carbon_ = gmml::ConvertString<int>(gmml::ConvertT<char>(residue_string[dash_index - 1]));
	        }

                if(dash_index != residue_string.size() - 1)
                {
                    if(!std::isdigit(residue_string[dash_index + 1])){

			std::stringstream error_notice;
		        error_notice << "ERROR at residue " << residue_string << ": missing open valence position.";
		        condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                        //throw CondensedSequenceProcessingException("Invalid residue in sequence");
		    }
		    else{
                        this->oxygen_position_ = gmml::ConvertString<int>(gmml::ConvertT<char>(residue_string[dash_index + 1]));
		    }
                }

                else{

		    std::stringstream error_notice;
		    error_notice << "WARNING at residue " << residue_string << ": unresolved open valence position, since not exactly 1 character found after dash.Setting to default value 1.";
		    error_notice << "Ignore this warning if this residue is connected to the aglycone";
		    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));

                    this->oxygen_position_ = 1;
	        }

                size_t left_bracket = residue_string.find('[');
                size_t right_bracket = residue_string.find(']');
                if(left_bracket == std::string::npos || right_bracket == std::string::npos){
                    this->name_ = residue_string.substr(1, dash_index - 3);
		}
                else
                {
                    this->name_ = residue_string.substr(1, left_bracket - 1);
                    std::string derivatives = residue_string.substr(left_bracket + 1, right_bracket - left_bracket - 1);
                    std::vector<std::string> derivatives_tokens = gmml::Split(derivatives, ",");
                    for(std::vector<std::string>::iterator it = derivatives_tokens.begin(); it != derivatives_tokens.end(); ++it)
                    {
                        if(!std::isdigit(it->at(0))){

			    std::stringstream error_notice;
			    error_notice << "ERROR at residue " << residue_string << ": Illegal derivative position at " << *it << " ,derivative postion must be an integer.";
			    condensed_sequence->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, error_notice.str()));
                            //throw CondensedSequenceProcessingException("Invalid derivative position in sequence");
		        }
		        else{
                            this->derivatives_[gmml::ConvertString<int>(gmml::ConvertT<char>(it->at(0)))] = it->substr(1);
		        }
                    }
                }
	    }//else
        }
    }
    parent_id_ = gmml::iNotSet;
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
bool CondensedSequenceResidue::GetIsTerminalAglycone()
{
    return is_terminal_aglycone_;
}
bool CondensedSequenceResidue::GetIsTerminalSugar()
{
    return is_terminal_sugar_;
}
std::string CondensedSequenceResidue::GetIsomer()
{
    return isomer_;
}
std::string CondensedSequenceResidue::GetConfiguration()
{
    return configuration_;
}
std::string CondensedSequenceResidue::GetName()
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
std::vector<int> CondensedSequenceResidue::GetChildIds()
{
    return child_ids_;
}
int CondensedSequenceResidue::GetBondId()  //Added by Yao 08/03/2018
{
    return bond_id_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequenceResidue::SetIsTerminalAglycone(bool is_terminal_aglycone)
{
    is_terminal_aglycone_ = is_terminal_aglycone;
}
void CondensedSequenceResidue::SetIsTerminalSugar(bool is_terminal_sugar)
{
    is_terminal_sugar_ = is_terminal_sugar;
}
void CondensedSequenceResidue::SetIsomer(std::string isomer)
{
    isomer_ = isomer;
}
void CondensedSequenceResidue::SetConfiguration(std::string configuration)
{
    configuration_ = configuration;
}
void CondensedSequenceResidue::SetName(std::string name)
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
        std::string derivative = (*it).second;
        derivatives_[index] = derivative;
    }
}
void CondensedSequenceResidue::SetParentId(int parent_id)
{
    parent_id_ = parent_id;
}
void CondensedSequenceResidue::AddChildId(int child_id)
{
    child_ids_.push_back(child_id);
}
void CondensedSequenceResidue::RemoveChildId(int child_id){
    if (std::find(child_ids_.begin(), child_ids_.end(), child_id) != child_ids_.end()){
        child_ids_.erase(std::find(child_ids_.begin(), child_ids_.end(), child_id));
    }
}
void CondensedSequenceResidue::SetBondId(int bond_id)  //Added by Yao 08/03/2018
{
    bond_id_ = bond_id;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequenceResidue::Print(std::ostream &out)
{
    out << "";
}
