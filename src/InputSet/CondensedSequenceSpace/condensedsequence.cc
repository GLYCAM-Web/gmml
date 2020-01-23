#include <stdexcept>
#include <stack>
#include <iostream>
#include <boost/algorithm/string.hpp>
#include <map>
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/InputSet/Utilities/response.hpp"
#include "../../../includes/MolecularModeling/assembly.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceresidue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceglycam06residue.hpp"
#include "../../../includes/InputSet/CondensedSequenceSpace/condensedsequenceprocessingexception.hpp"
#include "../../../includes/common.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/MolecularMetadata/GLYCAM/glycam06residueinfo.hpp"
using CondensedSequenceSpace::CondensedSequence;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequence::CondensedSequence()
{
}

CondensedSequence::CondensedSequence(std::string sequence)
{
    input_sequence_ = sequence;
    residues_ = CondensedSequenceResidueVector();
    tokens_ = CondensedSequenceTokenTypeVector();
    condensed_sequence_residue_tree_ = CondensedSequenceResidueTree();
    bool sequence_is_sane = ParseSequenceAndCheckSanity(sequence);
    if (sequence_is_sane){  //Also if MD service is requested
        bool MD_eligible = BuildArrayTreeOfCondensedSequenceGlycam06Residue(this->condensed_sequence_residue_tree_);
	if (!MD_eligible){
    	    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "This sequence is not eligible for MD."));
	}
    }
    DetectAnomericAnomericLinkages();
    if (this->anomeric_anomeric_linkages_.size() > 1){
	//Find the anomeric-anomeric linkage, if there is one.  Name that linkage '0'.  If there is more than one, punt.
	//Throw exception?std::exit(1)?Add error notice? Somehow stop the code from doing anything else. 
    }
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
CondensedSequence::CondensedSequenceGlycam06ResidueTree CondensedSequence::GetCondensedSequenceGlycam06ResidueTree()
{
    return condensed_sequence_glycam06_residue_tree_;
}
InputOutput::Response CondensedSequence::GetResponse()  //This is for gems to obtain a copy of the response object
{
    return this->response_;
}

void CondensedSequence::WriteGraphVizDotFile(GraphVizDotConfig& configs)
{
    // Declare the LABELS vector
    // This vector will hold all of the possible Residues we have svg files for.
    // The order of the Residues in the vector will matter because it will iterate
    // one at a time and once it finds one, it will break out of the loop and stop
    // looking because the std::string::find function will find multiple versions
    // of the residue string in one name. I.E. GlcNAc and GlcN will both be found
    // in residue GlcNAc. To combat this, the searching will look for GlcNAc first
    // and if it isn't found, it will then look for GlcN.
    std::vector<std::string> svg_residue_names = {
        "",
        "All",
        "Alt",
        "Ara",
        "Fru",
        "Fuc",
		"GlcNAc", "GlcA", "GlcN"/* Has some extra code to add. */, "Glc",
		"GalA", "GalNAc", "GalN", "Gal",
        "Gul",
        "Ido",
        "IdoA",
        "KDN",
        "KDO",
        "Lyx",
		"ManNAc", "ManA", "ManN", "Man",
		"Neu5Ac", "NeuNAc",
		"Neu5Gc", "NeuNGc",
        "Psi",
        "Qui",
        "Rha",
        "Rib"
        "Sor",
        "Tag",
        "Tal",
		"Xyl",
	};

	// Declare the BOX shape vector.
    std::vector<std::string> box_shape_residue_names = {
        "",
		"OH", "OME", "OtBu",
		"R", "Fak"
    };

	// Declare the OTHER shape vector.
    std::vector<std::string> other_shape_residue_names = {
        "",
    };

    // Declare an empty std::stringstream, which is used to generate what will go into
    // the dot file.
    std::stringstream ss;

    // Start generating the information to be put into the file.
    ss << "graph G {" << std::endl;
    ss << "\tgraph [splines=false forcelabels=true  dpi=" << configs.dpi_ << "];" << std::endl;
    ss << "\tnode [shape=\"none\" fontname=Helvetica labelfontsize=12 forcelabels=\"true\"; label=\"none\" size=50 fixedsize=\"true\" scale=\"true\"];" << std::endl;
    ss << "\tedge [labelfontsize=12 fontname=Helvetica labeldistance=1.2 labelangle = 320.0];" << std::endl;
    ss << "\trankdir=RL" << std::endl;
    ss << "\tnodesep=\"0.01\" ;" << std::endl;

    // Iterate through the Residues to add their information to the dot file.
    for( unsigned int i = 0; i < this->condensed_sequence_residue_tree_.size(); i++) {
		CondensedSequenceSpace::CondensedSequenceResidue* residue = this->condensed_sequence_residue_tree_[i];
		// Need to first get some variables from the current residue.
		// This could probably be split up into the specific section in the loop,
		//  but this was easier to track for writing this function.
      	std::string residue_name = residue->GetName();
    	std::string residue_isomer = residue->GetIsomer();
		std::string residue_ring_type = "";
		std::string residue_file_name = residue_name;
		if(residue_name.length() > 3) {
			residue_ring_type = residue_name.substr(3, 1);
			residue_name = residue_name.erase(3, 1);
		}
		DerivativeMap derivatives = residue->GetDerivatives();

		// Start to construct the information for the residue.
		ss << "\t" << i << " ";

		// If the residue is in the LABELS vector then we have an svg image for it.
		if( std::find( svg_residue_names.begin(), svg_residue_names.end(), residue_name ) != svg_residue_names.end() ) {
				ss << "[label=\"\" image=\"" << configs.svg_directory_path_ << residue_isomer << residue_file_name << ".svg\" height=\"0.7\"]";
		}
		// If the residue is not in the LABELS vector, but is in the BOX shape vector.
		else if( std::find( box_shape_residue_names.begin(), box_shape_residue_names.end(), residue_name ) != box_shape_residue_names.end() ) {
			// Start to construct the shape information.
			ss << "[shape=box ";
			if( residue_name == "R" ) {
				ss << "fontsize=20 label=\"" << "R";
			} else if( residue_name == "Fak" ) {
				ss << "fontsize=22 label=\"" << "R2";
			} else {
				ss << "label=\"" << residue_name;
			}
			// Close the angle brackets.
			ss << "\"]";
		}
		// If the residue is not in the LABELS or BOX vector, but is in the OTHER shape vector.
		else if(  std::find( other_shape_residue_names.begin(), other_shape_residue_names.end(), residue_name ) != other_shape_residue_names.end() ||
                  residue_ring_type == "p" ||
                  residue_ring_type == "f" ) {
		  // If this is a pyranose then the shape is a hexagon.
		  if( residue_ring_type == "p" ) {
			  ss << "[shape=hexagon ";
	          if( residue_isomer == "D" ) {
	            ss << "label=\"";
	          } else if( residue_isomer == "L" ) {
	            ss << "style=dashed label=\"";
	          }
		  }
		  // If this is a furanose then the shape is a pentagon.
		  else if( residue_ring_type == "f" ) {
			  ss << "[shape=pentagon ";
			  if( residue_isomer == "D" ) {
				  ss << "label=\"";
			  } else if( residue_isomer == "L" ) {
				  ss << "style=dashed label=\"";
			  }
		  }
		  // Close the angle brackets.
		  ss << residue_name << "\"]";
	   }
       ss << std::endl;

	   // Need to gather some information if there are derivatives for this residue.
	   if( !derivatives.empty() ) {
		   std::stringstream derivatives_stream;
		   // Iterate through the derivatives to construct the derivatives label.
		   for( DerivativeMap::iterator it = derivatives.begin(); it != derivatives.end(); it++ ) {
			   int position = (*it).first;
			   std::string derivative = (*it).second;
			   derivatives_stream << position << derivative << " ";
		   }

		   // Turn the stream into a string.
		   std::string derivatives_label = derivatives_stream.str();

		   // Start to construct the information to be printed.
		   ss << "\tb" << i << " [shape=\"plaintext\",";

		   // If the derivatives label is too long then split it up into two lines.
		   if( derivatives_label.length() > 10 ) {
			   ss << "fontsize=\"9\",forcelabels=\"true\"; height=\"0.3\"; labelloc=b; label=\"" << derivatives_label.substr(0,8) << "\n" << derivatives_label.substr(9) << "\"]" << std::endl;
		   } else {
			   ss << "fontsize=\"12\",forcelabels=\"true\"; height=\"0.3\"; labelloc=b;  label=\"" << derivatives_label << "\"]" << std::endl;
		   }

		   // More derivative information.
		   ss << "\t{rank=\"same\"; b" << i << " " << i << "}" << std::endl;
		   ss << "\t{nodesep=\"0.02\"; b" << i << ";" << i << "}" << std::endl;
		   ss << "\tb" << i << "--" << i << " [style=invis];" << std::endl;
	   }
	}

	for( unsigned int i = 1; i < this->condensed_sequence_residue_tree_.size(); i++ ) {
		// Get the residue at the index.
		CondensedSequenceSpace::CondensedSequenceResidue* residue = this->condensed_sequence_residue_tree_[i];
		// Pull some variables from the residue.
		int residue_parent_id = residue->GetParentId();
		std::string residue_configuration = residue->GetConfiguration();
		int residue_oxygen_position = residue->GetOxygenPosition();

		// We store the CondensedSequenceResidue.configuration_ as UpperCase, so this will convert it to LowerCase.
		std::transform(residue_configuration.begin(), residue_configuration.end(), residue_configuration.begin(), ::tolower);

		// Start constructing the infromation to be printed.
		ss << "\t" << residue_parent_id << "--" << i << " [";
		if( configs.show_config_labels_ ) {
			ss << "headlabel=\"" << residue_configuration << "\"";
		}
		// Don't need to print taillabel or edge label for first Residue.
		// TODO This may need to change if the terminal residue is a sugar.
		if( i > 1 ) {
			if( configs.show_position_labels_ ) {
				ss << " taillabel=\"" << residue_oxygen_position  << "\"";
			}
			if( configs.show_edge_labels_ ) {
				ss << " label=<<B>" << (i - 1) << "</B>>";
			}
		}
		// Close the angle brackets
		ss << "]" << std::endl;
	}
	// Close the curly brace.
    ss << "}" << std::endl;

	// Declare std::ofstream Object, which opens file. This will override
	// anything previously in the file.
	std::ofstream out_file(configs.file_name_, std::fstream::trunc);
	// Push the std::stringstream into the file.
	out_file << ss.str();
	// Make sure to close the file when finished.
	out_file.close();
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequence::AddNoteToResponse(Glycan::Note* new_note)
{
    this->response_.AddNotice(new_note);
}

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
void CondensedSequence::AddToken(gmml::CondensedSequenceTokenType token)
{
    tokens_.push_back(token);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
bool CondensedSequence::ParseSequenceAndCheckSanity(std::string sequence)
{
    ParseCondensedSequence(sequence, this);  //There is error checking in the condense residue constuctor class as well, but the code construct does not allow me to let this function return.
    std::string bad_residue_notice = "Bad residue information detected.Cannot proceed.See other notices.Aborted.";
    std::string good_residue_notice = "Sequence is sane.";
    bool residues_are_sane = this->CheckResidueTokenSanity();
    if (!residues_are_sane){
	this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, bad_residue_notice));
        return false;
    }

    else {
        int condensed_residue_tree_success_status = BuildArrayTreeOfCondensedSequenceResidue();
        if (condensed_residue_tree_success_status != 0){
	    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, bad_residue_notice));
            return false;
	}

    	else{
	    bool connections_are_sane = this->CheckLinkageAndDerivativeSanity();
	    if (!connections_are_sane){
	        this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, bad_residue_notice));
	        return false;
	    }

	    else{
	        this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::COMMENT, Glycan::NoteCat::GLYCOSIDIC, good_residue_notice));
	        return true;
	    }
        }
    }
}

bool CondensedSequence::CheckResidueTokenSanity()
{
    for (CondensedSequenceResidueVector::iterator res_it = this->residues_.begin(); res_it != this->residues_.end(); res_it++){

	CondensedSequenceResidue* condensed_residue = *res_it;
	std::string isomer = condensed_residue->GetIsomer();
        std::string anomeric_configuration = condensed_residue->GetConfiguration();
        std::string residue_name = condensed_residue->GetName();
        int anomeric_position = condensed_residue->GetAnomericCarbon();

	//Read in Metadata
	//Search Metadata for residue name. If no match throw error.
	//If there exists an entry, get its likely isomers, anomeric configuration,and anomeric configuration.
	//Compare with input values. If does not match for isomer and anomeric configuration.Throw warning.For anomeric positon mismatch throw error.
    }
    return true;
}

bool CondensedSequence::CheckLinkageAndDerivativeSanity()
{
    for (CondensedSequenceResidueVector::iterator res_it = this->condensed_sequence_residue_tree_.begin(); res_it != this->condensed_sequence_residue_tree_.end(); res_it++){

        CondensedSequenceResidue* condensed_residue = *res_it;
        std::string residue_name = condensed_residue->GetName();
        //Read in Metadata
        //Search Medata with residue name(which is valid now) as key, get allowed open valence positions.
        //Make a map<int, bool> that records wheter each available position is occupied.
        CondensedSequenceResidue::DerivativeMap derivatives = condensed_residue->GetDerivatives();
        //Check each derivative position. Throw error if attached to disallowed open valence positions.
        //For legal derivatives, mark its position as occupied.
        std::vector<int> child_ids = condensed_residue->GetChildIds();
        for (std::vector<int>::iterator child_it = child_ids.begin(); child_it != child_ids.end(); child_it++){
            int child_id = *child_it;
            CondensedSequenceResidue* child_residue = this->condensed_sequence_residue_tree_[child_id];
            int child_open_valence_position = child_residue->GetOxygenPosition();
            //If child open valence does not exist as key in map, then it's accessing disallowed positions. If that position is allowed by occupied, that's also an error.
        }
    }
    return true;
}
int CondensedSequence::InsertNodeInCondensedSequenceResidueTree(CondensedSequenceResidue *condensed_residue, int parent_node_id, int bond_id)
{
    if (parent_node_id >= (int)condensed_sequence_residue_tree_.size()){
	std::stringstream note_str;
	note_str << "Invalid parent index at residue " << condensed_residue->GetName() << ", parent index is " << parent_node_id;
	this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, note_str.str()));

	parent_node_id = gmml::iNotSet;
        //throw std::invalid_argument("ArrayTree::insert - invalid parent index(" + gmml::ConvertT(parent_node_id) + ")");
    }
    //    condensed_sequence_residue_tree_.push_back(std::make_pair(condensed_residue, parent_node_id));
    condensed_residue->SetParentId(parent_node_id);
    condensed_residue->SetBondId(bond_id);

    if (parent_node_id != gmml::iNotSet){
        condensed_sequence_residue_tree_[parent_node_id]->AddChildId(condensed_sequence_residue_tree_.size());  //Not only child knows parent, parent also knows child.
    }
    condensed_sequence_residue_tree_.push_back(condensed_residue);
    return condensed_sequence_residue_tree_.size() - 1;
}

int CondensedSequence::InsertNodeInCondensedSequenceGlycam06ResidueTree(CondensedSequenceSpace::CondensedSequenceGlycam06Residue* condensed_glycam06_residue, int parent_node_id, int bond_id)
{
    if(parent_node_id != gmml::iNotSet && parent_node_id >= (int)condensed_sequence_glycam06_residue_tree_.size())
        throw std::invalid_argument("ArrayTree::insert - invalid parent index(" + gmml::ConvertT(parent_node_id) + ")");
//    condensed_sequence_glycam06_residue_tree_.push_back(std::make_pair(condensed_tree_residue, parent_node_id));
    condensed_glycam06_residue->SetParentId(parent_node_id);
    condensed_glycam06_residue->SetBondId(bond_id);
    condensed_sequence_glycam06_residue_tree_.push_back(condensed_glycam06_residue);
    return condensed_sequence_glycam06_residue_tree_.size() - 1;
}

bool CondensedSequence::ParseCondensedSequence(std::string sequence, CondensedSequence* condensed_sequence)
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
                std::string residue = sequence.substr(start_index, end_index - start_index + 1);
                residues_.push_back(new CondensedSequenceResidue(residue, this));
                if(terminal)
                    start_index = i + 1;
                else
                {
                    start_index = i + 2;
                    i++;
                }
                tokens_.push_back(gmml::CONDENSED_SEQUENCE_RESIDUE);
                reading_residue = false;
                break;
            }
            case '[':
            {
                if(!reading_residue)
                {
                    tokens_.push_back(gmml::CONDENSED_SEQUENCE_LEFT_BRACKET);
                    start_index = i + 1;
                }
                break;
            }
            case ']':
            {
                if(!reading_residue)
                {
                    tokens_.push_back(gmml::CONDENSED_SEQUENCE_RIGHT_BRACKET);
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
    std::string terminal_residue = sequence.substr(sequence.find_last_of('-') + 1);
    if(terminal_residue.compare("") != 0)
    {
        residues_.push_back(new CondensedSequenceResidue(terminal_residue, this));
        tokens_.push_back(gmml::CONDENSED_SEQUENCE_RESIDUE);
    }
    return true;
}

int CondensedSequence::BuildArrayTreeOfCondensedSequenceResidue()
{
    CondensedSequenceTokenTypeVector::reverse_iterator current_token = tokens_.rbegin();
    CondensedSequenceResidueVector::reverse_iterator current_residue = residues_.rbegin();

    int return_value = 0; //zero is success.Non zero is error.
    if(residues_.size() == 0){
        return_value = 1;
	this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "Condensed sequence does not contain any detectable residues."));
    }

    std::stack<int> residue_stack;
    residue_stack.push(this->InsertNodeInCondensedSequenceResidueTree(*current_residue, gmml::iNotSet, gmml::iNotSet));
    int bond_count = 0;  //Added by Yao 08/03/2018. Whenever a parent-child relationship is found, there is a new bond. This int count the index of the new bond to be added, and becomes its label.
    while(++current_token != tokens_.rend())
    {
        switch(*current_token)
        {
            case gmml::CONDENSED_SEQUENCE_LEFT_BRACKET:
            {
                if (residue_stack.empty()){
		    return_value = 1;
		    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "Invalid branching in sequence"));
                    //throw CondensedSequenceProcessingException("Invalid branching in sequence");
		}
                residue_stack.pop();
                break;
            }
            case gmml::CONDENSED_SEQUENCE_RIGHT_BRACKET:
            {
                ++current_token;
                ++current_residue;
                if(current_residue == residues_.rend()){
		    return_value = 1;
		    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "Invalid sequence of residues"));
                    //throw CondensedSequenceProcessingException("Invalid sequence of residues");
		}
                if(residue_stack.empty()){
		    return_value = 1;
		    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "Invalid sequence"));
                    //throw CondensedSequenceProcessingException("Invalid sequence");
		}
		if (current_residue != residues_.rend() && !residue_stack.empty()){
                    residue_stack.push(this->InsertNodeInCondensedSequenceResidueTree(*current_residue, residue_stack.top(), bond_count));
		    bond_count++;
		}
                break;
            }
            case gmml::CONDENSED_SEQUENCE_RESIDUE:
            {
                current_residue++;
                if(current_residue == residues_.rend()){
		    return_value = 1;
		    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "Invalid sequence of residues"));
                    //throw CondensedSequenceProcessingException("Invalid sequence of residues");
		}
                if(residue_stack.empty()){
		    return_value = 1;
		    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::ERROR, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, "Invalid sequence of residues"));
                    //throw CondensedSequenceProcessingException("Invalid sequence");
		}
		if(current_residue != residues_.rend() && !residue_stack.empty()){
                    int parent = this->InsertNodeInCondensedSequenceResidueTree(*current_residue, residue_stack.top(), bond_count);
                    residue_stack.pop();
                    residue_stack.push(parent);
		    bond_count++;
		}
                break;
            }
        }
    }
    return return_value;
}

int CondensedSequence::ReEvaluateParentIdentityUponAnomericAnomericLinkage()
{
    //Count the residues on either side of '0'.  If one side has more residues, then that side becomes the "left" side.  The 'left' side gets the positive integer labels/indexes.
    //If both sides have the same number of residues, the one with more branches goes on the left.
    //If both sides have the same number of residues and branches, then the lowest-numbered anomeric carbon goes on the left.
    //If all are the same still, then it goes alphabetically or by whatever else we decide.
    std::pair<int, int> anomeric_linkage = this->anomeric_anomeric_linkages_[0]; 
    int current_parent_residue_index = anomeric_linkage.first;
    int the_other_residue_index = anomeric_linkage.second;
    int num_residues_1 = 0, num_branches_1 = 0, num_residues_2 = 0, num_branches_2 = 0;
    this->RecursivelyCountNumberofDownstreamResiduesAndBranches(current_parent_residue_index, num_residues_1, num_branches_1);
    this->RecursivelyCountNumberofDownstreamResiduesAndBranches(the_other_residue_index, num_residues_2, num_branches_2);
    int reevaluated_parent_index = current_parent_residue_index;
    if (num_residues_1 < num_residues_2){
	reevaluated_parent_index = the_other_residue_index;
    }
    else if (num_residues_1 == num_residues_2){
        if (num_branches_1 < num_branches_2){
	    reevaluated_parent_index = the_other_residue_index;
	}
	else if (num_branches_1 == num_branches_2){
	    // What does "lowest-numbered anomeric carbon" mean?
	}
    }

    if (reevaluated_parent_index == the_other_residue_index){ //If we need to switch parent to the other residue
        CondensedSequenceResidue* current_parent_residue = this->condensed_sequence_residue_tree_[current_parent_residue_index];
        CondensedSequenceResidue* the_other_residue = this->condensed_sequence_residue_tree_[the_other_residue_index];
        current_parent_residue->RemoveChildId(the_other_residue_index);
        current_parent_residue->SetParentId(the_other_residue_index);
	the_other_residue->AddChildId(current_parent_residue_index);
	the_other_residue->SetParentId(gmml::iNotSet);
	current_parent_residue->SetIsTerminalSugar(false);
	the_other_residue->SetIsTerminalSugar(true);
	current_parent_residue->SetAnomericCarbon(the_other_residue->GetOxygenPosition());
	current_parent_residue->SetOxygenPosition(the_other_residue->GetAnomericCarbon());
    }
    return reevaluated_parent_index;
}

void CondensedSequence::RecursivelyCountNumberofDownstreamResiduesAndBranches(int parent_residue_index, int& num_residues, int& num_branches)
{
    CondensedSequenceResidue* parent_residue = this->condensed_sequence_residue_tree_[parent_residue_index];
    std::vector<int> child_ids = parent_residue->GetChildIds();
    //Exclude the anomeric-anomeric path. Remove the other residue (not the current parent) in the ano-ano link from child residues, if it exists.
    if (std::find(child_ids.begin(), child_ids.end(), this->anomeric_anomeric_linkages_[0].second) != child_ids.end()){
        child_ids.erase(std::find(child_ids.begin(), child_ids.end(), this->anomeric_anomeric_linkages_[0].second));
    }
    unsigned int num_childs = child_ids.size();
    num_residues += num_childs;
    if (num_childs > 1){
        num_branches += (num_childs - 1);
    }
    for (unsigned int i = 0; i < num_childs; i++){
        int child_residue_index = child_ids[i];
        this->RecursivelyCountNumberofDownstreamResiduesAndBranches(child_residue_index, num_residues, num_branches);
    }
}

void CondensedSequence::DetectAnomericAnomericLinkages()
{
    gmml::MolecularMetadata::GLYCAM::Glycam06NamesToTypesLookupContainer metadata_residueNamesToTypes;
    for (unsigned int i = 0; i < this->condensed_sequence_residue_tree_.size(); i++){
	CondensedSequenceResidue* residue = condensed_sequence_residue_tree_[i];
	std::string glycam_06_name = condensed_sequence_glycam06_residue_tree_[i]->GetName();
	std::vector<std::string> all_types = metadata_residueNamesToTypes.GetTypesForResidue(glycam_06_name);
	int anomeric_position = 0;

	if (std::find(all_types.begin(), all_types.end(), "aldose") != all_types.end()){
	    anomeric_position = 1;
	}
	else if (std::find(all_types.begin(), all_types.end(), "ketose") != all_types.end()){
	    anomeric_position = 2;
	}

	std::vector<int> child_ids = residue->GetChildIds();
	for (unsigned int j = 0; j < child_ids.size(); j++){
	    CondensedSequenceResidue* child_residue = condensed_sequence_residue_tree_[child_ids[j]];
	    int child_parent_open_valence_position = child_residue->GetOxygenPosition();
            if (child_parent_open_valence_position == anomeric_position){
		this->anomeric_anomeric_linkages_.push_back(std::make_pair(i, child_ids[j]));
	    }
	}


	
    }

}

bool CondensedSequence::BuildArrayTreeOfCondensedSequenceGlycam06Residue(CondensedSequenceResidueTree residue_tree)
{
    bool MD_eligible = false;
    condensed_sequence_glycam06_residue_tree_ = CondensedSequenceGlycam06ResidueTree();
    std::vector<std::vector<int> > open_valences = std::vector<std::vector<int> >(residue_tree.size());
    for(unsigned int i = 0; i < residue_tree.size(); i++)
    {
        int parent = residue_tree.at(i)->GetParentId();
        CondensedSequenceResidue* residue = residue_tree.at(i);
        if(parent != gmml::iNotSet)
        {
            int oxygen_position = residue->GetOxygenPosition();
            open_valences[parent].push_back(oxygen_position);
        }

        CondensedSequenceResidue::DerivativeMap condensed_residue_derivatives = residue->GetDerivatives();
        for(CondensedSequenceResidue::DerivativeMap::iterator it = condensed_residue_derivatives.begin(); it != condensed_residue_derivatives.end(); ++it)
        {
            int derivative_index = it->first;
	    std::string derivative_name = it->second;
	    if (derivative_name != "D"){  //Deoxy shouldn't be considered as open valence
	        open_valences[i].push_back(derivative_index);
	    }
	}
    }
    CondensedSequenceResidue* terminal_residue = residue_tree.at(0);
    std::string terminal = residue_tree.at(0)->GetName();

    int current_derivative_count = 0;
    std::vector<int> derivatives = std::vector<int>(residue_tree.size(), 0);
    derivatives[0] = current_derivative_count;
    if (terminal_residue->GetIsTerminalAglycone()){
        this->InsertNodeInCondensedSequenceGlycam06ResidueTree(new CondensedSequenceSpace::CondensedSequenceGlycam06Residue(this->GetGlycam06TerminalResidueCodeOfTerminalResidue(terminal)), 
			                                       gmml::iNotSet, gmml::iNotSet);
    }
    else if (terminal_residue->GetIsTerminalSugar()){
        CondensedSequenceSpace::CondensedSequenceGlycam06Residue* tree_residue = new CondensedSequenceSpace::CondensedSequenceGlycam06Residue(this->GetGlycam06ResidueCodeOfCondensedResidue(
                                                                                                        terminal_residue, open_valences[0]));
        this->InsertNodeInCondensedSequenceGlycam06ResidueTree(tree_residue, gmml::iNotSet, gmml::iNotSet); //Giving fake derivative and bond id. They will be deprecated soon anyway.  
    }

    for(unsigned int i = 1; i < residue_tree.size(); i++)
    {
        derivatives[i] = current_derivative_count;
        CondensedSequenceResidue* condensed_residue = residue_tree.at(i);
	int condensed_residue_bond_id = condensed_residue->GetBondId();
        int parent = residue_tree.at(i)->GetParentId();


        std::string parent_name = residue_tree.at(parent)->GetName();
        std::string anomeric_carbon = "C" + gmml::ConvertT<int>(condensed_residue->GetAnomericCarbon());
        std::string oxygen_position;
        (parent_name.compare("OME") == 0) ? oxygen_position = "O" : oxygen_position = "O" + gmml::ConvertT<int>(condensed_residue->GetOxygenPosition());

        try
        {
	    int glycam_06_residue_bond_id = condensed_residue_bond_id + current_derivative_count;
            CondensedSequenceSpace::CondensedSequenceGlycam06Residue* tree_residue = new CondensedSequenceSpace::CondensedSequenceGlycam06Residue(this->GetGlycam06ResidueCodeOfCondensedResidue(
                                                                                                        condensed_residue, open_valences[i])
                                                                                                    , anomeric_carbon, oxygen_position);

            int residue_index = this->InsertNodeInCondensedSequenceGlycam06ResidueTree(tree_residue, parent + derivatives[parent], glycam_06_residue_bond_id);

            CondensedSequenceResidue::DerivativeMap condensed_residue_derivatives = condensed_residue->GetDerivatives();
            for(CondensedSequenceResidue::DerivativeMap::iterator it = condensed_residue_derivatives.begin(); it != condensed_residue_derivatives.end(); ++it)
            {
		glycam_06_residue_bond_id++;
                std::string derivative_name = it->second;
                int derivative_index = it->first;
                this->InsertNodeInCondensedSequenceGlycam06ResidueTree(this->GetCondensedSequenceDerivativeGlycam06Residue(derivative_name, derivative_index), residue_index, glycam_06_residue_bond_id);
                current_derivative_count++;
            }
        }
        catch(std::exception ex)
        {
            CondensedSequenceSpace::CondensedSequenceGlycam06Residue* tree_residue = new CondensedSequenceSpace::CondensedSequenceGlycam06Residue(condensed_residue->GetName().substr(0,3)
                                                                                                    , anomeric_carbon, oxygen_position);

            this->InsertNodeInCondensedSequenceGlycam06ResidueTree(tree_residue, parent + derivatives[parent], gmml::iNotSet);

	    std::stringstream notice;
	    notice << "Residue Not eligible for MD: (" << condensed_residue->GetName().substr(0,3) << ")";
	    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, notice.str()));

	    MD_eligible = false;
            //throw CondensedSequenceProcessingException("Invalid residue in the sequence (" + condensed_residue->GetName().substr(0,3) + ")");
        }
    }
    MD_eligible = true;
    return MD_eligible;
}

std::string CondensedSequence::GetGlycam06TerminalResidueCodeOfTerminalResidue(std::string terminal_residue_name)
{
    if(terminal_residue_name.compare("OH") == 0 || terminal_residue_name.compare("ROH") == 0)
        return "ROH";
    else if(terminal_residue_name.compare("OME") == 0)
        return "OME";
    else if(terminal_residue_name.compare("OtBu") == 0 || terminal_residue_name.compare("TBT") == 0)
        return "TBT";
    else if (terminal_residue_name.compare("UNK") == 0)
	return "UNK";
    else if(gmml::AminoacidGlycamLookup(terminal_residue_name).aminoacid_name_.compare("") != 0 ||
            gmml::AminoacidGlycamLookup(terminal_residue_name).glycam_name_.compare("") != 0)
        return gmml::AminoacidGlycamLookup(terminal_residue_name).glycam_name_;
    else {
        //throw CondensedSequenceProcessingException("Invalid aglycon " + terminal_residue_name);
    }
    return "UNK"; //To prevent "Control reaches end of non-void function"
}

std::string CondensedSequence::GetGlycam06ResidueCodeOfCondensedResidue(CondensedSequenceResidue *condensed_residue, std::vector<int> open_valences)
{
    if(condensed_residue->GetName().compare("UNK") == 0 || condensed_residue->GetName().compare("Unknown") == 0)
        return "UNK";
    std::vector<int> new_valences_list = open_valences;
    CondensedSequenceResidue::DerivativeMap condensed_residue_derivatives = condensed_residue->GetDerivatives();
    for(CondensedSequenceResidue::DerivativeMap::iterator it = condensed_residue_derivatives.begin(); it != condensed_residue_derivatives.end(); ++it)
        new_valences_list.push_back(it->first);
    std::string residue_name = condensed_residue->GetName();
    std::string isomer = condensed_residue->GetIsomer();
    std::string configuration = condensed_residue->GetConfiguration();
    std::string ring_type = "P";
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

    std::bitset<10> open_valences_check = std::bitset<10>();
    for(unsigned int i = 0; i < open_valences.size(); i ++)
    {
        if(open_valences[i] < 0 || open_valences[i] >= 10){
	    std::stringstream notice;
	    notice << "Not eligible for MD: " << condensed_residue->GetName() << "unsupported open valence position: <0 or >= 10";
	    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, notice.str()));
            //throw CondensedSequenceProcessingException("Invalid open valence");
	}
        open_valences_check.set(open_valences[i]);
    }

    std::string residue_code = this->GetFirstLetterOfGlycam06ResidueCode(open_valences_check) + this->GetSecondLetterOfGlycam06ResidueCode(residue_name, isomer);
    if(residue_code.size() < 3)
        residue_code += this->GetThirdLetterOfGlycam06ResidueCode(configuration, ring_type);
    return residue_code;
}

std::string CondensedSequence::GetFirstLetterOfGlycam06ResidueCode(std::bitset<10> open_valences_check)
{
    std::bitset<10> open_valences_check_temp = open_valences_check;
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
        for (unsigned int i = 1; i < open_valences_check_temp.size(); i++)
            if (open_valences_check_temp[i])
                return gmml::ConvertT<int>(i);
    }
    else if (open_valences_check_temp.none())
    {
        return "0";
    }

    std::string notice = "Not eligible for MD: There is no code in the GLYCAM code set for residues with open valences at the given positions.";
    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, notice));

    //throw CondensedSequenceProcessingException("There is no code in the GLYCAM code set for residues with open valences at the given positions.");
    return "";
}

std::string CondensedSequence::GetSecondLetterOfGlycam06ResidueCode(std::string residue_name, std::string isomer)
{
    gmml::ResidueCodeName residue_name_code = gmml::ResidueNameCodeLookup(residue_name);
    if(residue_name_code.name_.compare("") != 0)
    {
        std::string code = residue_name_code.code_;
        if(isomer.compare("L") == 0)
            std::transform(code.begin(), code.end(), code.begin(), ::tolower);
        return code;
    }

    std::stringstream notice;
    notice << "Not eligible for MD: " << residue_name << "This residue name has no corresponding entry in glycam06 force field.";
    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, notice.str()));
    //throw CondensedSequenceProcessingException(residue_name + " is not a valid residue");
    return "";
}

std::string CondensedSequence::GetThirdLetterOfGlycam06ResidueCode(std::string configuration, std::string ring_type)
{
    if(configuration.compare("X") == 0)
        return "X";
    if(ring_type.compare("P") == 0)
        return (configuration.compare("A") == 0) ? "A" : "B";
    else
        return (configuration.compare("A") == 0) ? "D" : "U";
}

CondensedSequenceSpace::CondensedSequenceGlycam06Residue* CondensedSequence::GetCondensedSequenceDerivativeGlycam06Residue(std::string derivative_name, int derivative_index)
{
    //Oxygen names are usually OK, but why does an SO3 have a C atom? Here basicaaly oxygen_name = tail atom name, carbon_name = head atom name. Right now, the carbon_name has to be hard-coded.
    std::string oxygen_name = "O" + gmml::ConvertT<int>(derivative_index);
    //std::string carbon_name = "C" + gmml::ConvertT<int>(derivative_index); //Can't always use C* as head atom name.
    if(derivative_name.compare("S") == 0)
        return new CondensedSequenceSpace::CondensedSequenceGlycam06Residue("SO3", "S1", oxygen_name, true);
    else if(derivative_name.compare("Me") == 0)
        return new CondensedSequenceSpace::CondensedSequenceGlycam06Residue("MEX", "CH3", oxygen_name, true);
    else if(derivative_name.compare("A") == 0)
        return new CondensedSequenceSpace::CondensedSequenceGlycam06Residue("ACX", "C1A", oxygen_name, true);
    else if(derivative_name.compare("D") == 0)
	//D means deoxy.This derivative is not a template, but an action, of removing the oxygen this derivative attaches to. Here I create a false glycam06 residue for this purpose later on.
        return new CondensedSequenceSpace::CondensedSequenceGlycam06Residue("Deoxy", "Deoxy",oxygen_name, true);
    //Later PO3 might need to be added, but now I dont' now the name of its head atom yet.
    std::stringstream notice;
    notice << "Not eligible for MD: " << derivative_name << "This derivative name has no corresponding entry in glycam06 force field.";
    this->AddNoteToResponse(new Glycan::Note(Glycan::NoteType::WARNING, Glycan::NoteCat::IMPROPER_CONDENSED_SEQUENCE, notice.str()));
    //throw CondensedSequenceProcessingException("There is no derivative in the GLYCAM code set represented by the letter " + derivative_name);
    // If none of the above, return something:
//    std::cout << "WARNING: There is no derivative in the GLYCAM code set represented by the letter " << derivative_name << "\n";
    return new CondensedSequenceSpace::CondensedSequenceGlycam06Residue("UNK", "UNK", oxygen_name, true);
}

std::string CondensedSequence::BuildLabeledCondensedSequence(CondensedSequence::Reordering_Approach labeling_approach, CondensedSequence::Reordering_Approach reordering_approach, bool label)
{
    //WARNING: preserve user input option doesn't work well with the labeing function, but okay with reordering function. 
    std::string labeled_sequence = "";
    int reevaluated_parent_index = 0; //The default parent in the input sequence is 0
    int branch_depth = 0;
    std::vector<int> longest_path;  //Elements in this vector is the residue ids on the longest path. Will be used or ignored based on reordering approach. See below.
    std::map<unsigned int, std::string> residue_label_map;
    std::map<unsigned int, std::string> bond_label_map;

    if (reordering_approach != CondensedSequence::Reordering_Approach::PRESERVE_USER_INPUT && this->anomeric_anomeric_linkages_.size() == 1){
	reevaluated_parent_index = this->ReEvaluateParentIdentityUponAnomericAnomericLinkage();
    }

    if (labeling_approach == CondensedSequence::Reordering_Approach::LONGEST_CHAIN || reordering_approach == CondensedSequence::Reordering_Approach::LONGEST_CHAIN){
        this->FindLongestPath(longest_path);
    }

    if (label){
        int current_residue_label_index = 0, current_bond_label_index = 0;
        RecursivelyLabelCondensedSequence(reevaluated_parent_index, current_residue_label_index, current_bond_label_index, residue_label_map, bond_label_map, labeling_approach, longest_path); 
    }

    this->RecursivelyBuildLabeledCondensedSequence(reevaluated_parent_index, branch_depth, labeled_sequence, reordering_approach, residue_label_map, bond_label_map, longest_path, label);
    //Starts with the first condensed residue in the vector, which is the aglycone. Then go recursively down the tree.
    return labeled_sequence;
}

void CondensedSequence::FindLongestPath (std::vector<int>& longest_path)
{
    for (unsigned int i = 0; i < this->condensed_sequence_residue_tree_.size(); i++){
        if (this->condensed_sequence_residue_tree_[i]->GetChildIds().size() == 0){
            CondensedSequenceSpace::CondensedSequenceResidue* current_residue = this->condensed_sequence_residue_tree_[i];  //This finds a non reducing end residue;
            std::vector<int> current_path = std::vector<int>();
            current_path.push_back(i); //i is the index of this non reducing end residue

            while(current_residue->GetParentId() != gmml::iNotSet){ //While loop stops at a residue that has no parent, which is the aglycone.
                current_path.push_back(current_residue->GetParentId());
                current_residue = this->condensed_sequence_residue_tree_[current_residue->GetParentId()];

                if (current_residue->GetParentId() == gmml::iNotSet){  //If the aglycone is reached.
                    bool new_path_found = false;
                    if (current_path.size() > longest_path.size()){  //If a longer path is found, overwrite the current longest path.
                        new_path_found = true;
                    }
                    else if (current_path.size() == longest_path.size()){  //If the new path is exactly the same length, choose based on lower branch index at the first deverging point.
                        CondensedSequenceSpace::CondensedSequenceResidue* branching_residue = NULL;
                        for (int j=current_path.size()-1; j>=0; j--){  //reverse iterating the path to find the first branching point.
                            if (this->condensed_sequence_residue_tree_[current_path[j]]->GetChildIds().size() > 1){
                                branching_residue = this->condensed_sequence_residue_tree_[current_path[j]];
                                break;
                            }
                        }
                        if (branching_residue == NULL){
//                            std::cout << "Two paths at equal lengths are found. Can't decide which residue has lower branch index.This is a bug." << std::endl;
                        }
                        else{
                            int current_path_branching_index = 0, longest_path_branching_index = 0;
                            std::vector<int> child_ids = branching_residue->GetChildIds();
                            for (std::vector<int>::iterator it = child_ids.begin(); it != child_ids.end(); it++){
                                int child_id = *it;
                                if (std::find(current_path.begin(), current_path.end(), child_id) != current_path.end()){
                                    current_path_branching_index = child_id;
                                }
                                if (std::find(longest_path.begin(), longest_path.end(), child_id) != longest_path.end()){
                                    longest_path_branching_index = child_id;
                                }
                            }
                            if (current_path_branching_index < longest_path_branching_index){
                                new_path_found = true;
                            }
                        }
                    }
                    if (new_path_found){
                        longest_path.clear();
                        for (unsigned int k=0; k<current_path.size(); k++){
                            longest_path.push_back(current_path[k]);
                        }
                    }
                }
            }
        }
    }
}

void CondensedSequence::RecursivelyLabelCondensedSequence(int current_residue_index, int& current_residue_label_index, int& current_bond_label_index, 
		                                          std::map<unsigned int, std::string>& residue_label_map, std::map<unsigned int, std::string>& bond_label_map, 
							  CondensedSequence::Reordering_Approach labeling_approach, std::vector<int>& longest_path)
{
    CondensedSequenceSpace::CondensedSequenceResidue* current_residue = this->condensed_sequence_residue_tree_[current_residue_index];
    //Insert residue label after residue name,example: Glcp&Label_residueId=1;a1-4
    std::string residue_label; 
    residue_label.insert(0, ";"); //semicolon is the right delimiter of a label;
    residue_label.insert(0, std::to_string(current_residue_label_index)); //Add residue index to label.
    residue_label.insert(0, "&Label=residue-"); //&Label is the left delimiter of a label. 
    residue_label_map[current_residue_index] = residue_label;
    current_residue_label_index++;

    //Insert bond label. A bond label is asoociated with the child residue in the bond, not the parent. 
    //No bond label for aglycone/terminal residue
    std::string bond_label;
    if (!(current_residue->GetIsTerminalAglycone() || current_residue->GetIsTerminalSugar())){  //Only create bond label for non-terminal/aglycone residue.
        bond_label.insert(0, ";");
        bond_label.insert(0, std::to_string(current_bond_label_index));
        bond_label.insert(0, "&Label=link-");
        current_bond_label_index++;
    }
    bond_label_map[current_residue_index] = bond_label;

    std::vector<int> child_ids = current_residue->GetChildIds(); //Get the index of child residues. Sort this vector based on the reordering approach below. 

    if (labeling_approach == CondensedSequence::Reordering_Approach::LOWEST_INDEX){
        //If sequence is to be reordered based on index, rearrange child in ascending order based on parent open valene position.Use a temporary map for sorting.
        std::map<int, int, std::less<int> > openvalence_childid_map = std::map<int, int, std::less<int> >();
        for (std::vector<int>::iterator it= child_ids.begin(); it != child_ids.end(); it++){
            int child_id = *it;
            CondensedSequenceSpace::CondensedSequenceResidue* child_residue = this->condensed_sequence_residue_tree_[child_id];
            int open_valence_position = child_residue->GetOxygenPosition();
            openvalence_childid_map[open_valence_position] = child_id;
        }
        child_ids.clear();
        for (std::map<int, int>::iterator map_it = openvalence_childid_map.begin(); map_it != openvalence_childid_map.end(); map_it++){
            child_ids.push_back(map_it->second);
        }
    }

    else if (labeling_approach == CondensedSequence::Reordering_Approach::LONGEST_CHAIN){
        //If sequence is to be reordered based on longest chain, sort main chain child that's on the predetermined longest path to the front of vector, preserve the order of the rest
        for (std::vector<int>::iterator it= child_ids.begin(); it != child_ids.end(); it++){
            int child_id = *it;
            if (std::find(longest_path.begin(), longest_path.end(), child_id) != longest_path.end()){
                child_ids.erase(it);
                child_ids.insert(child_ids.begin(), child_id);
		break;
            }
        }
    }

    for (std::vector<int>::iterator it = child_ids.begin(); it != child_ids.end(); it++){
        this->RecursivelyLabelCondensedSequence(*it, current_residue_label_index, current_bond_label_index, residue_label_map, bond_label_map, labeling_approach, longest_path);
    }
}

void CondensedSequence::RecursivelyBuildLabeledCondensedSequence(int current_index, int& branch_depth, std::string& labeled_sequence, CondensedSequence::Reordering_Approach reordering_approach, std::map<unsigned int, std::string>& residue_label_map, std::map<unsigned int, std::string>& bond_label_map, std::vector<int>& longest_path, bool label)
{
    CondensedSequenceSpace::CondensedSequenceResidue* current_residue = this->condensed_sequence_residue_tree_[current_index];
    if (current_residue->GetIsTerminalAglycone()){  //If current residue is aglycone

        if (label){
	    labeled_sequence.insert(0, residue_label_map[current_index]);
        }
        //Done residue label

        labeled_sequence.insert(0, current_residue->GetName());
        labeled_sequence.insert(0, "-");
    }

    else{
        if (label){
	    labeled_sequence.insert(0, bond_label_map[current_index]);
        }
        //Done bond label
        /*if (current_residue->GetIsTerminalSugar()){ //If current residue is terminal sugar (i.e. anomeric-anomeric linkage)
	    labeled_sequence.insert(0, boost::algorithm::to_lower_copy(current_residue->GetConfiguration()));
            labeled_sequence.insert(0, current_residue->GetName());
            labeled_sequence.insert(0, current_residue->GetIsomer());
        }*/
	
        //The residues connected to the aglycone has OxygenPosition set to 1, actually there shouldn't be such a value.So ignore.
	//else if (this->condensed_sequence_residue_tree_[current_residue->GetParentId()]->GetIsTerminalSugar()){
        //}

	if (!current_residue->GetIsTerminalSugar()){
            labeled_sequence.insert(0, std::to_string(current_residue->GetOxygenPosition()));
            labeled_sequence.insert(0, "-");
            labeled_sequence.insert(0, std::to_string(current_residue->GetAnomericCarbon()));
	}

        labeled_sequence.insert(0, boost::algorithm::to_lower_copy(current_residue->GetConfiguration()));  //In residue class anomeric configuration is upper case. Convert to lower case.

        CondensedSequenceResidue::DerivativeMap derivatives = current_residue->GetDerivatives();
        if (derivatives.size() > 0){
            labeled_sequence.insert(0, "]");
            for (CondensedSequenceResidue::DerivativeMap::reverse_iterator derivative_rit = derivatives.rbegin(); derivative_rit != derivatives.rend(); derivative_rit++){
                labeled_sequence.insert(0,derivative_rit->second);
                labeled_sequence.insert(0,std::to_string(derivative_rit->first));
                labeled_sequence.insert(0,",");
            }
            labeled_sequence.erase(0,1);
            labeled_sequence.insert(0, "[");
        }

        //Insert residue label after residue name,example: Glcp&Label_residueId=1;a1-4
        if (label){
	    labeled_sequence.insert(0, residue_label_map[current_index]);
        }
        //Done residue label
	
        labeled_sequence.insert(0, current_residue->GetName());
        labeled_sequence.insert(0, current_residue->GetIsomer());

    }

    std::vector<int> child_ids = current_residue->GetChildIds();

    if (reordering_approach == CondensedSequence::Reordering_Approach::LOWEST_INDEX){
        //If sequence is to be reordered based on index, rearrange child in descending order based on parent open valene position.Use a temporary map for sorting.
        std::map<int, int, std::greater<int> > openvalence_childid_map = std::map<int, int, std::greater<int> >();
        for (std::vector<int>::iterator it= child_ids.begin(); it != child_ids.end(); it++){
            int child_id = *it;
            CondensedSequenceSpace::CondensedSequenceResidue* child_residue = this->condensed_sequence_residue_tree_[child_id];
            int open_valence_position = child_residue->GetOxygenPosition();
            openvalence_childid_map[open_valence_position] = child_id;
        }
        child_ids.clear();
        for (std::map<int, int>::iterator map_it = openvalence_childid_map.begin(); map_it != openvalence_childid_map.end(); map_it++){
            child_ids.push_back(map_it->second);
        }
    }
    else if (reordering_approach == CondensedSequence::Reordering_Approach::LONGEST_CHAIN){
        //If sequence is to be reordered based on longest chain, sort main chain child on the longest chain to the end of vector, preserve the order of the rest
        for (std::vector<int>::iterator it= child_ids.begin(); it != child_ids.end(); it++){
            int child_id = *it;
            if (std::find(longest_path.begin(), longest_path.end(), child_id) != longest_path.end()){
                child_ids.erase(it);
                child_ids.push_back(child_id);
            }
        }
    }

    if (child_ids.size() == 0 && branch_depth > 0){  //No child means we're at a non-reducing end residue. if at a branch, add left bracket [
        labeled_sequence.insert(0,"[");  //End of a branch
        branch_depth --;
    }


    for (std::vector<int>::iterator it = child_ids.begin(); it != child_ids.end(); it++){
        if (child_ids.size() > 1 && it != child_ids.end()-1){
            branch_depth++;
            labeled_sequence.insert(0,"]");
        }
        this->RecursivelyBuildLabeledCondensedSequence(*it, branch_depth, labeled_sequence, reordering_approach, residue_label_map, bond_label_map, longest_path, label);
    }

}

CondensedSequence::CondensedSequenceRotamersAndGlycosidicAnglesInfo CondensedSequence::GetCondensedSequenceRotamersAndGlycosidicAnglesInfo(CondensedSequenceResidueTree residue_tree)
{
    CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles = CondensedSequenceRotamersAndGlycosidicAnglesInfo();
    //int linkage_index = 0;
    int linkage_index = -1; //For testing front end. -- Yao
    for(unsigned int i = 0; i < residue_tree.size(); i++)
    {
        int parent = residue_tree.at(i)->GetParentId();
        if(parent >= 0)
        {
            linkage_index++;
            CondensedSequenceResidue* residue = residue_tree.at(i);
            std::string residue_absolute_name = residue->GetName().substr(0, 3) + residue->GetName().substr(4);
            std::string parent_res_name = residue_tree[parent]->GetName();
            char ring_letter = residue->GetName()[3];
            std::vector<std::pair<std::string, std::vector<std::string> > > possible_rotamers = std::vector<std::pair<std::string, std::vector<std::string> > >();
            std::vector<std::pair<std::string, std::vector<std::string> > > selected_rotamers = std::vector<std::pair<std::string, std::vector<std::string> > >();
            std::vector<std::pair<std::string, double> > enabled_glycosidic_angles = std::vector<std::pair<std::string, double> >();
            enabled_glycosidic_angles.push_back(std::make_pair("phi", gmml::dNotSet));
            enabled_glycosidic_angles.push_back(std::make_pair("psi", gmml::dNotSet));
            if(ring_letter == 'p')
            {
                if(parent >= 0) //For the first sugar residue, its parent is the aglycon,whose index is 0. So parent > 0 prevents addition of 1st sugar-aglycon bond into vector.
                {

                    CondensedSequenceResidue* parent_residue = residue_tree.at(parent);
                    std::string parent_residue_absolute_name = "";
                    //If parent is something like an aglycon, OH, for example, its residue name is too short for substr(4), which causes segfault.
                    //Below I added an if statement that if length of parent residue name is lower than 5, copy parent residue name. Otherwise, do the regular substring manipulation.
		    if (parent_residue->GetName().length() < 5){
			parent_residue_absolute_name = parent_residue->GetName();
                    }
                    else{
                        parent_residue_absolute_name = parent_residue->GetName().substr(0,3) + parent_residue->GetName().substr(4);
                    }
                    std::stringstream rotamers_name;
                    rotamers_name << residue->GetIsomer() << residue->GetName() << residue->GetConfiguration() << residue->GetAnomericCarbon() << "-" <<
                                     residue->GetOxygenPosition() << parent_residue->GetIsomer() << parent_residue->GetName() << parent_residue->GetConfiguration();
                    switch(gmml::ResidueNameIndexLookup(residue_absolute_name).index_)
                    {
                        case 10:
                        case 1:
                        case 2:
                        case 12:
                        case 21:
                        case 20:
                            if(residue->GetOxygenPosition() == 6)
                            {
                                std::vector<std::string> rot = std::vector<std::string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                possible_rotamers.push_back(std::make_pair("omega", rot));

                                std::vector<std::string> rot1 = std::vector<std::string>();
                                rot1.push_back("gg");
                                rot1.push_back("gt");
                                selected_rotamers.push_back(std::make_pair("omega", rot1));

                                enabled_glycosidic_angles.push_back(std::make_pair("omega", gmml::dNotSet));
                            }
                            break;
                        case 14:
                        case 13:
                        case 33:
                        case 9:
                        case 7:
                            if(residue->GetOxygenPosition() == 6)
                            {
                                std::vector<std::string> rot = std::vector<std::string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                possible_rotamers.push_back(std::make_pair("omega", rot));
                                selected_rotamers.push_back(std::make_pair("omega", rot));

                                enabled_glycosidic_angles.push_back(std::make_pair("omega", gmml::dNotSet));
                            }
                            break;
                        case 23:
                        case 24:
                            switch(gmml::ResidueNameIndexLookup(parent_residue_absolute_name).index_)
                            {
                                case 10:
                                case 1:
                                case 2:
                                case 12:
                                case 21:
                                case 20:
                                    if(residue->GetOxygenPosition() == 6)
                                    {
                                        std::vector<std::string> rot = std::vector<std::string>();
                                        rot.push_back("gg");
                                        rot.push_back("gt");
                                        rot.push_back("tg");
                                        possible_rotamers.push_back(std::make_pair("omega", rot));

                                        std::vector<std::string> rot1 = std::vector<std::string>();
                                        rot1.push_back("gg");
                                        rot1.push_back("gt");
                                        selected_rotamers.push_back(std::make_pair("omega", rot1));

                                        enabled_glycosidic_angles.push_back(std::make_pair("omega", gmml::dNotSet));
                                    }
                                    break;
                                case 14:
                                case 13:
                                case 33:
                                case 9:
                                case 7:
                                    if(residue->GetOxygenPosition() == 6)
                                    {
                                        std::vector<std::string> rot = std::vector<std::string>();
                                        rot.push_back("gg");
                                        rot.push_back("gt");
                                        rot.push_back("tg");
                                        possible_rotamers.push_back(std::make_pair("omega", rot));
                                        selected_rotamers.push_back(std::make_pair("omega", rot));

                                        enabled_glycosidic_angles.push_back(std::make_pair("omega", gmml::dNotSet));
                                    }
                                    break;
                            }
                            if(gmml::ResidueNameIndexLookup(parent_residue_absolute_name).index_ != 23 &&
                                    gmml::ResidueNameIndexLookup(parent_residue_absolute_name).index_ != 24)
                            {
                                std::vector<std::string> rot = std::vector<std::string>();
                                rot.push_back("t");
                                rot.push_back("g");
                                rot.push_back("-g");
                                possible_rotamers.push_back(std::make_pair("phi", rot));

                                std::vector<std::string> rot1 = std::vector<std::string>();
                                rot1.push_back("t");
                                rot1.push_back("-g");
                                selected_rotamers.push_back(std::make_pair("phi", rot1));
                            }
                            break;
                    }
                    RotamersAndGlycosidicAnglesInfo* info = new RotamersAndGlycosidicAnglesInfo(linkage_index, possible_rotamers, selected_rotamers, enabled_glycosidic_angles);
                    RotamerNameInfoPair pair_info = std::make_pair(rotamers_name.str(), info);

                    rotamers_glycosidic_angles.push_back(pair_info);
                }

                std::vector<std::pair<std::string, std::vector<std::string> > > der_possible_rotamers = std::vector<std::pair<std::string, std::vector<std::string> > >();
                std::vector<std::pair<std::string, std::vector<std::string> > > der_selected_rotamers = std::vector<std::pair<std::string, std::vector<std::string> > >();
                std::vector<std::pair<std::string, double> > der_enabled_glycosidic_angles = std::vector<std::pair<std::string, double> >();
                der_enabled_glycosidic_angles.push_back(std::make_pair("phi", gmml::dNotSet));
                der_enabled_glycosidic_angles.push_back(std::make_pair("psi", gmml::dNotSet));
                CondensedSequenceResidue::DerivativeMap derivatives = residue->GetDerivatives();
                for(CondensedSequenceResidue::DerivativeMap::iterator it = derivatives.begin(); it != derivatives.end(); it++)
                {
                    int derivative_index = (*it).first;
                    std::string derivative_name = (*it).second;
                    switch(gmml::ResidueNameIndexLookup(residue_absolute_name).index_)
                    {
                        case 10:
                        case 1:
                        case 2:
                        case 12:
                        case 21:
                        case 20:
                            if(derivative_index == 6 && derivative_name.compare("S") == 0)
                            {
                                std::vector<std::string> rot = std::vector<std::string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                der_possible_rotamers.push_back(std::make_pair("omega", rot));

                                std::vector<std::string> rot1 = std::vector<std::string>();
                                rot1.push_back("gg");
                                rot1.push_back("gt");
                                der_selected_rotamers.push_back(std::make_pair("omega", rot1));

                                der_enabled_glycosidic_angles.push_back(std::make_pair("omega", gmml::dNotSet));
                            }
                            break;
                        case 14:
                        case 13:
                        case 33:
                        case 9:
                        case 7:
                            if(derivative_index == 6 && derivative_name.compare("S") == 0)
                            {
                                std::vector<std::string> rot = std::vector<std::string>();
                                rot.push_back("gg");
                                rot.push_back("gt");
                                rot.push_back("tg");
                                der_possible_rotamers.push_back(std::make_pair("omega", rot));
                                der_selected_rotamers.push_back(std::make_pair("omega", rot));

                                der_enabled_glycosidic_angles.push_back(std::make_pair("omega", gmml::dNotSet));
                            }
                            break;
                    }
                    //linkage_index++;  //Do not increment bond index upon derivative
                    std::stringstream der_rotamers_name;
                    der_rotamers_name << residue->GetIsomer() << residue->GetName() << residue->GetConfiguration()
                                      << "[" << derivative_index << derivative_name << "]";
                    RotamersAndGlycosidicAnglesInfo* der_info = new RotamersAndGlycosidicAnglesInfo(linkage_index, der_possible_rotamers,
                                                                                                    der_selected_rotamers, der_enabled_glycosidic_angles);
                    RotamerNameInfoPair der_pair_info = std::make_pair(der_rotamers_name.str(), der_info);

                    rotamers_glycosidic_angles.push_back(der_pair_info);
                }
            }

            if(ring_letter == 'f')
            {
                if(parent > 0)
                {

                    CondensedSequenceResidue* parent_residue = residue_tree.at(parent);
                    std::string parent_residue_absolute_name = parent_residue->GetName().substr(0,3) + parent_residue->GetName().substr(4);
                    std::stringstream rotamers_name;
                    rotamers_name << residue->GetIsomer() << residue->GetName() << residue->GetConfiguration() << residue->GetAnomericCarbon() << "-" <<
                                  residue->GetOxygenPosition() << parent_residue->GetIsomer() << parent_residue->GetName() << parent_residue->GetConfiguration();

                    RotamersAndGlycosidicAnglesInfo* info = new RotamersAndGlycosidicAnglesInfo(linkage_index, possible_rotamers, selected_rotamers, enabled_glycosidic_angles);
                    RotamerNameInfoPair pair_info = std::make_pair(rotamers_name.str(), info);

                    rotamers_glycosidic_angles.push_back(pair_info);
                }




                std::vector<std::pair<std::string, std::vector<std::string> > > der_possible_rotamers = std::vector<std::pair<std::string, std::vector<std::string> > >();
                std::vector<std::pair<std::string, std::vector<std::string> > > der_selected_rotamers = std::vector<std::pair<std::string, std::vector<std::string> > >();
                std::vector<std::pair<std::string, double> > der_enabled_glycosidic_angles = std::vector<std::pair<std::string, double> >();
                der_enabled_glycosidic_angles.push_back(std::make_pair("phi", gmml::dNotSet));
                der_enabled_glycosidic_angles.push_back(std::make_pair("psi", gmml::dNotSet));
                CondensedSequenceResidue::DerivativeMap derivatives = residue->GetDerivatives();
                for(CondensedSequenceResidue::DerivativeMap::iterator it = derivatives.begin(); it != derivatives.end(); it++)
                {
                    int derivative_index = (*it).first;
                    std::string derivative_name = (*it).second;
                    linkage_index++;
                    std::stringstream der_rotamers_name;
                    der_rotamers_name << residue->GetIsomer() << residue->GetName() << residue->GetConfiguration()
                                      << "[" << derivative_index << derivative_name << "]";
                    RotamersAndGlycosidicAnglesInfo* der_info = new RotamersAndGlycosidicAnglesInfo(linkage_index, der_possible_rotamers,
                                                                                                    der_selected_rotamers, der_enabled_glycosidic_angles);
                    RotamerNameInfoPair der_pair_info = std::make_pair(der_rotamers_name.str(), der_info);

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
        std::vector<std::pair<std::string, double> > angles = rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_;
        std::vector<std::pair<std::string, std::vector<std::string> > > selected_rotamers = rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_;
        int omega_rotamers = 0;
        int phi_rotamers = 0;
        for(unsigned int j = 0; j < selected_rotamers.size(); j++)
        {
            std::pair<std::string, std::vector<std::string> > rot = selected_rotamers.at(j);
            if(rot.first.compare("phi") == 0)
                phi_rotamers+= rot.second.size();
            if(rot.first.compare("omega") == 0)
                omega_rotamers += rot.second.size();
        }
        if(angles.at(0).second == gmml::dNotSet) // Phi
            count *= (phi_rotamers == 0) ? 1 : phi_rotamers;
        if(angles.size() == 3)      // Phi, Psi, Omega
        {
            if(angles.at(2).second == gmml::dNotSet)
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
        std::string rotamer_name = rotamers_glycosidic_angles_info.at(i).first;
        if(rotamer_name.find("DNeup5AcA2-8") != std::string::npos ||
                rotamer_name.find("DNeup5AcB2-8") != std::string::npos)
            count *= 6;
        if(rotamer_name.find("DNeup5GcA2-8") != std::string::npos)
            count *= 5;
    }
    return count;
}

std::vector<std::vector<int> > CondensedSequence::CreateBaseMapAllPossibleSelectedRotamers(CondensedSequenceRotamersAndGlycosidicAnglesInfo rotamers_glycosidic_angles_info)
{
    std::vector<std::vector<int> > mapper = std::vector<std::vector<int> >();
    for(unsigned int i = 0; i < rotamers_glycosidic_angles_info.size(); i++)
    {
        std::vector<std::pair<std::string, double> > angles = rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_;
        std::vector<std::pair<std::string, std::vector<std::string> > > selected_rotamers = rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_;
        int omega_rotamers = 0;
        int psi_rotamers = 1;
        int phi_rotamers = 0;
        for(unsigned int j = 0; j < selected_rotamers.size(); j++)
        {
            std::pair<std::string, std::vector<std::string> > rot = selected_rotamers.at(j);
            if(rot.first.compare("phi") == 0)
                phi_rotamers += rot.second.size();
            if(rot.first.compare("omega") == 0)
                omega_rotamers += rot.second.size();
        }
        if(angles.at(0).second == gmml::dNotSet) // Phi
            phi_rotamers = (phi_rotamers == 0) ? 1 : phi_rotamers;
        if(angles.size() == 3)      // Phi, Psi, Omega
        {
            if(angles.at(2).second == gmml::dNotSet)
                omega_rotamers = (omega_rotamers == 0) ? 1 : omega_rotamers;
        }

        std::vector<int> val = std::vector<int>();
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
    std::vector<std::vector<int> > mapping = this->CreateBaseMapAllPossibleSelectedRotamers(rotamers_glycosidic_angles_info);
    std::vector<std::vector<std::vector<double> > > phi_psi_omega_vector_map = std::vector<std::vector<std::vector<double> > >();
    std::vector<std::vector<std::vector<std::string> > > phi_psi_omega_rt_vector_map = std::vector<std::vector<std::vector<std::string> > >();
    for(unsigned int i = 0; i < mapping.size(); i++)
    {
        std::vector<double> phi = std::vector<double>();
        std::vector<std::string> phi_rt = std::vector<std::string>();
        // phi => gmml::dNotSet in all phi cases means default value
        if(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(0).second == gmml::dNotSet) // value not set
        {
            if(rotamers_glycosidic_angles_info.at(i).second->possible_rotamers_.size() == 0) // no possible rotamer
            {
                phi.push_back(gmml::dNotSet);
                phi_rt.push_back("df");
            }
            else if(rotamers_glycosidic_angles_info.at(i).second->selected_rotamers_.size() == 0) // no selected rotamer
            {
                phi.push_back(gmml::dNotSet);
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
                            phi.push_back(gmml::dNotSet);
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
                                    phi.push_back(gmml::dNotSet);
                                    phi_rt.push_back("df");
                                }
                            }
                        }
                        phi_check = true;
                    }
                }
                if(!phi_check)
                {
                    phi.push_back(gmml::dNotSet);
                    phi_rt.push_back("df");
                }
            }
        }
        else
        {
            phi.push_back(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(0).second);
            phi_rt.push_back("cu");
        }

        std::vector<double> psi = std::vector<double>();
        std::vector<std::string> psi_rt = std::vector<std::string>();
        if(rotamers_glycosidic_angles_info.at(i).first.find("[") == std::string::npos)
        {
            psi.push_back(gmml::dNotSet); // standard psi angle
            psi_rt.push_back("df");
        }
        else
        {
            psi.push_back(gmml::dNotSet);
            psi_rt.push_back("df");
        }

        // omega => gmml::dNotSet in all omega cases means not applicable
        // for the cases that the angle is not set by any of selecting rotamers or setting the angle, it is set to 180.0
        std::vector<double> omega = std::vector<double>();
        std::vector<std::string> omega_rt = std::vector<std::string>();
        if(mapping.at(i).at(2) == 0)
        {
            omega.push_back(gmml::dNotSet);
            omega_rt.push_back("na");
        }
        else
        {
            if(rotamers_glycosidic_angles_info.at(i).second->enabled_glycosidic_angles_.at(2).second == gmml::dNotSet) // value not set
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
        std::vector<std::vector<double> > phi_psi_omega = std::vector<std::vector<double> >();
        std::vector<std::vector<std::string> > phi_psi_omega_rt = std::vector<std::vector<std::string> >();
        for(unsigned int j = 0; j < phi.size(); j++)
        {
            for(unsigned int k = 0; k < psi.size(); k++)
            {
                for(unsigned int l = 0; l < omega.size(); l++)
                {
                    std::vector<double> combination = std::vector<double>();
                    std::vector<std::string> combination_rt = std::vector<std::string>();
                    combination.push_back(phi.at(j));
                    combination_rt.push_back(phi_rt.at(j));
                    combination.push_back(psi.at(k));
                    combination_rt.push_back(psi_rt.at(k));
                    combination.push_back(omega.at(l));
                    combination_rt.push_back(omega_rt.at(l));
                    std::string rotamer_name = rotamers_glycosidic_angles_info.at(i).first;
                    if(rotamer_name.find("DNeup5AcA2-8") != std::string::npos ||
                            rotamer_name.find("DNeup5AcB2-8") != std::string::npos)
                    {
                        for(int m = 0; m < 6; m++)
                        {
                            std::vector<double> new_combination = combination;
                            std::vector<std::string> new_combination_rt = combination_rt;
                            for(int n = 0; n < 5; n++)
                            {
                                new_combination.push_back(gmml::EXTERNAL28LINKAGEROTAMERS[m][n]);
                                new_combination_rt.push_back("E" + gmml::ConvertT<int>(m) + "/" + gmml::ConvertT<int>(n));
                            }
                            phi_psi_omega.push_back(new_combination);
                            phi_psi_omega_rt.push_back(new_combination_rt);
                        }
                    }
                    else if(rotamer_name.find("DNeup5GcA2-8") != std::string::npos)
                    {
                        for(int m = 0; m < 5; m++)
                        {
                            std::vector<double> new_combination = combination;
                            std::vector<std::string> new_combination_rt = combination_rt;
                            for(int n = 0; n < 5; n++)
                            {
                                new_combination.push_back(gmml::INTERNAL28LINKAGEROTAMERS[m][n]);
                                new_combination_rt.push_back("I" + gmml::ConvertT<int>(m) + "/" + gmml::ConvertT<int>(n));
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
        std::vector<std::vector<double> > phi_psi_omega_combination = phi_psi_omega_vector_map.at(i);
        for(unsigned int j = 0; j < phi_psi_omega_combination.size(); j++)
        {
            std::vector<double> phi_psi_omega = phi_psi_omega_combination.at(j);
            std::cout << "<";
            for(unsigned int k = 0; k < phi_psi_omega.size(); k++)
            {
                std::cout << phi_psi_omega.at(k);
                (k < phi_psi_omega.size() - 1) ? std::cout << ", " : std::cout << "";
            }
            std::cout << ">";
        }
        std::cout << std::endl;
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
                mapper[k] = std::vector<std::vector<double> >();
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
                        std::vector<std::vector<double> > res = new_mapper[i];
                        std::vector<std::vector<std::string> > res_rt = new_name_mapper[i];
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

    for(unsigned int i = 0; i < mapper.size(); i++)
    {
        /*std::cout << i << ": ";
        for(unsigned int j = 0; j < mapper[i].size(); j++)
        {
            std::cout << "<";
            for(unsigned int k = 0; k < mapper[i].at(j).size(); k++)
            {
                std::cout << mapper[i].at(j).at(k);
                (k < mapper[i].at(j).size() - 1) ? std::cout << ", " : std::cout << " ";
            }
            std::cout << ">";
        }*/
        std::stringstream ss;
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
        /*std::cout << "<" << ss.str().substr(0,ss.str().size()-1) << ">";
        std::cout << std::endl;*/
    }
    return mapper;
}

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequence::Print(std::ostream &out)
{
    out << "";
}
