#include "../../../includes/FileSet/PdbFileSpace/pdbmodelresidueset.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/FileSet/PdbFileSpace/pdbheterogenatomcard.hpp"
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
    while(!Trim(temp).empty())
    {
        stringstream atom_block;
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
        /// End of a ATOM section in the given residue set block
        if(line.find("TER") != string::npos)
        {
            PdbAtomCard* atoms = new PdbAtomCard(atom_block);
            this->AddAtom(atoms);
            /// Section may have another ATOM/HETATM section
            getline(residue_set_block, line);
            temp = line;
        }
        stringstream heterogen_atom_block;
        /// Extract HETATM section of the given residue set block
        if(line.find("HETATM") != string::npos)
        {
            while(line.find("HETATM") != string::npos)
            {
                heterogen_atom_block << line << endl;       /// Append all lines of HETATM section to create a block of stream of each heterogen atom set
                getline(residue_set_block, line);           /// Read next line
                temp = line;
            }
            PdbHeterogenAtomCard* heterogen_atoms = new PdbHeterogenAtomCard(heterogen_atom_block);
            this->AddHeterogenAtom(heterogen_atoms);
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
