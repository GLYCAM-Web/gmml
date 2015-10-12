#include "../../../includes/InputSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomcard.hpp"
#include "../../../includes/utils.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbModelResidueSet::PdbModelResidueSet() {}

PdbModelResidueSet::PdbModelResidueSet(stringstream &residue_set_block)
{
    string line;
    getline(residue_set_block, line);
    string temp = line;
    int atom_card_index = 0;
    int heterogen_atom_card_index = 0;
    while(!Trim(temp).empty())
    {
        stringstream atom_block;
        stringstream head_exceptional_atom_block;
        stringstream tail_exceptional_atom_block;
        stringstream heterogen_atom_block;

        while(line.find("HETATM") != string::npos)
        {
            head_exceptional_atom_block << line << endl;
            getline(residue_set_block, line);
            temp = line;
        }
        if(line.find("ATOM") != string::npos || line.find("ANISOU") != string::npos || line.find("TER") != string::npos)
        {
            atom_block << head_exceptional_atom_block.str();
        }
        else
        {
            heterogen_atom_block << head_exceptional_atom_block.str();
        }
        while(line.find("ATOM") != string::npos || line.find("ANISOU") != string::npos)
        {
            /// Extract ATOM section of the given residue set block
            while(line.find("ATOM") != string::npos)
            {
                atom_block << line << endl;             /// Append all lines of ATOM section to create a block of stream of each atom set
                getline(residue_set_block, line);       /// Read next line
                temp = line;
            }
            /// Discard ANISOU section
            while(line.find("ANISOU") != string::npos)
            {
                getline(residue_set_block,line);        /// Skip lines
                temp = line;
            }
        }
        while(line.find("HETATM") != string::npos || line.find("ANISOU") != string::npos)
        {
            tail_exceptional_atom_block << line << endl;
            getline(residue_set_block, line);
            temp = line;
            /// Discard ANISOU section
            while(line.find("ANISOU") != string::npos)
            {
                getline(residue_set_block,line);        /// Skip lines
                temp = line;
            }
        }
        if(line.find("TER") != string::npos)
        {
            atom_block << tail_exceptional_atom_block.str();
        }
        else
        {
            heterogen_atom_block << tail_exceptional_atom_block.str();
        }
        /// End of an ATOM section in the given residue set block
        if(line.find("TER") != string::npos || Trim(temp).empty())
        {
            stringstream card_index;
            card_index << "ATOM_" << atom_card_index;
            PdbAtomCard* atoms = new PdbAtomCard(atom_block, card_index.str());
            this->AddAtom(atoms);
            atom_card_index++;
            /// Check for the end of the block
            if(Trim(temp).empty())
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
        while(line.find("HETATM") != string::npos  || line.find("ANISOU") != string::npos)
        {
            heterogen_atom_block << line << endl;       /// Append all lines of HETATM section to create a block of stream of each heterogen atom set
            getline(residue_set_block, line);           /// Read next line
            temp = line;
            /// Discard ANISOU section
            while(line.find("ANISOU") != string::npos)
            {
                getline(residue_set_block,line);        /// Skip lines
                temp = line;
            }
        }
        if(Trim(temp).empty())
        {
            stringstream card_index;
            card_index << "HETATOM_" << heterogen_atom_card_index;
            PdbHeterogenAtomCard* heterogen_atoms = new PdbHeterogenAtomCard(heterogen_atom_block, card_index.str());
            this->AddHeterogenAtom(heterogen_atoms);
            heterogen_atom_card_index++;
            break;
        }
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbModelResidueSet::AtomCardVector PdbModelResidueSet::GetAtoms(){
    return atoms_;
}

PdbModelResidueSet::HeterogenAtomCardVector PdbModelResidueSet::GetHeterogenAtoms(){
    return heterogen_atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////

void PdbModelResidueSet::SetAtoms(PdbModelResidueSet::AtomCardVector atoms){
    atoms_.clear();
    for(AtomCardVector::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        atoms_.push_back(*it);
    }
}

void PdbModelResidueSet::AddAtom(PdbAtomCard *atom)
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

void PdbModelResidueSet::AddHeterogenAtom(PdbHeterogenAtomCard *heterogen_atom)
{
    heterogen_atoms_.push_back(heterogen_atom);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbModelResidueSet::Print(ostream &out)
{
    out << "----------------- Atoms -------------" << endl;
    for(PdbModelResidueSet::AtomCardVector::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << endl << "-------------- Heterogen Atoms --------------" << endl;
    for(PdbModelResidueSet::HeterogenAtomCardVector::iterator it = heterogen_atoms_.begin(); it != heterogen_atoms_.end(); it++)
    {
        (*it)->Print(out);
    }
    out << endl;
}
