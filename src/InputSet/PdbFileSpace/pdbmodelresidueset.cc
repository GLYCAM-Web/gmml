#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/utils.hpp"

using PdbFileSpace::PdbModelResidueSet;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelResidueSet::PdbModelResidueSet() {}

PdbModelResidueSet::PdbModelResidueSet(std::stringstream &residue_set_block)
{
    std::string line;
    getline(residue_set_block, line);
    std::string temp = line;
    int atom_card_index = 0;
    int heterogen_atom_card_index = 0;
    while(!gmml::Trim(temp).empty())
    {
        std::stringstream atom_block;
        std::stringstream head_exceptional_atom_block;
        std::stringstream tail_exceptional_atom_block;
        std::stringstream heterogen_atom_block;

        while(line.find("HETATM") != std::string::npos)
        {
            head_exceptional_atom_block << line << std::endl;
            getline(residue_set_block, line);
            temp = line;
        }
        if(line.find("ATOM") != std::string::npos || line.find("TER") != std::string::npos)
        {
            atom_block << head_exceptional_atom_block.str();
        }
        else
        {
            heterogen_atom_block << head_exceptional_atom_block.str();
        }
        while(line.find("ATOM") != std::string::npos || line.find("ANISOU") != std::string::npos)
        {
            /// Extract ATOM section of the given residue set block
            while(line.find("ATOM") != std::string::npos)
            {
                atom_block << line << std::endl;             /// Append all lines of ATOM section to create a block of stream of each atom set
                getline(residue_set_block, line);       /// Read next line
                temp = line;
            }
            /// Discard ANISOU section
            while(line.find("ANISOU") != std::string::npos)
            {
                getline(residue_set_block,line);        /// Skip lines
                temp = line;
            }
        }
        while(line.find("HETATM") != std::string::npos || line.find("ANISOU") != std::string::npos)
        {
            tail_exceptional_atom_block << line << std::endl;
            getline(residue_set_block, line);
            temp = line;
            /// Discard ANISOU section
            while(line.find("ANISOU") != std::string::npos)
            {
                getline(residue_set_block,line);        /// Skip lines
                temp = line;
            }
        }
        if(line.find("TER") != std::string::npos)
        {
            atom_block << tail_exceptional_atom_block.str();
        }
        else
        {
            heterogen_atom_block << tail_exceptional_atom_block.str();
        }

        /// End of an ATOM section in the given residue set block
        if(line.find("TER") != std::string::npos || gmml::Trim(temp).empty())
        {
            std::stringstream card_index;
            card_index << "ATOM_" << atom_card_index;
            PdbAtomSection* atoms = new PdbAtomSection(atom_block, card_index.str());
            this->AddAtom(atoms);
            atom_card_index++;
            /// Check for the end of the block
            if(gmml::Trim(temp).empty())
            {

            }
            /// Section may have another ATOM/HETATM section
            else
            {
                getline(residue_set_block, line);
                temp = line;
                continue;
            }
        }
        /// Extract HETATM section of the given residue set block
        while(line.find("HETATM") != std::string::npos || line.find("ANISOU") != std::string::npos)
        {
            heterogen_atom_block << line << std::endl;       /// Append all lines of HETATM section to create a block of stream of each heterogen atom set
            getline(residue_set_block, line);           /// Read next line
            temp = line;
            /// Discard ANISOU section
            while(line.find("ANISOU") != std::string::npos)
            {
                getline(residue_set_block,line);        /// Skip lines
                temp = line;
            }
        }
        if(gmml::Trim(temp).empty())
        {
            std::stringstream card_index;
            card_index << "HETATOM_" << heterogen_atom_card_index;
            PdbHeterogenAtomSection* heterogen_atoms = new PdbHeterogenAtomSection(heterogen_atom_block, card_index.str());
            this->AddHeterogenAtom(heterogen_atoms);
            heterogen_atom_card_index++;
            break;
        }
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbModelResidueSet::AtomCardVector PdbModelResidueSet::GetAtomCards(){
    return atoms_;
}

PdbModelResidueSet::HeterogenAtomCardVector PdbModelResidueSet::GetHeterogenAtomCards(){
    return heterogen_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void PdbModelResidueSet::SetAtomCards(PdbModelResidueSet::AtomCardVector atoms){
    atoms_.clear();
    for(AtomCardVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_.push_back(*it);
    }
}

void PdbModelResidueSet::AddAtom(PdbAtomSection *atom)
{
    atoms_.push_back(atom);
}

void PdbModelResidueSet::SetHeterogenAtoms(PdbModelResidueSet::HeterogenAtomCardVector heterogen_atoms){
    heterogen_atoms_.clear();
    for(HeterogenAtomCardVector::iterator it = heterogen_atoms.begin(); it != heterogen_atoms.end(); it++)
    {
        heterogen_atoms_.push_back(*it);
    }
}

void PdbModelResidueSet::AddHeterogenAtom(PdbHeterogenAtomSection *heterogen_atom)
{
    heterogen_atoms_.push_back(heterogen_atom);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModelResidueSet::Print(std::ostream &out)
{
    out << "----------------- Atoms -------------" << std::endl;
    for(PdbModelResidueSet::AtomCardVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << std::endl << "-------------- Heterogen Atoms --------------" << std::endl;
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it = heterogen_atoms_.begin(); it != heterogen_atoms_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << std::endl;
}
