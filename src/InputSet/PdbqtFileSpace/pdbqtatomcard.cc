#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatomcard.hpp"
#include "../../../includes/InputSet/PdbqtFileSpace/pdbqtatom.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbqtFileSpace::PdbqtAtomCard;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbqtAtomCard::PdbqtAtomCard() : record_name_("ATOM"){}

PdbqtAtomCard::PdbqtAtomCard(std::ifstream &stream_block)
{
    atoms_ = PdbqtAtomMap();
    record_name_ = "ATOM";
    std::string line;
    while (getline(stream_block, line))
    {
	if (line.find("ATOM") != std::string::npos || line.find("HETATM") != std::string::npos){
            PdbqtAtom* atom = new PdbqtAtom(line);
            atoms_[atom->GetAtomSerialNumber()] = atom;
	}
	else{ //File stream has gone one line beyond atom section. Rewind by one line and exit.
	    int offset = -1*((int)line.length() +1);  //Rewind file stream postion by length of current line + 1, to go back to the last line. 
            stream_block.seekg(offset, stream_block.cur); //Go back one line    
	    break;
	}
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbqtAtomCard::GetRecordName()
{
    return record_name_;
}

PdbqtAtomCard::PdbqtAtomMap PdbqtAtomCard::GetAtoms()
{
    return atoms_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbqtAtomCard::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbqtAtomCard::SetAtoms(PdbqtAtomMap atoms)
{
    atoms_.clear();
    for(PdbqtAtomMap::iterator it = atoms.begin(); it != atoms.end(); it++)
    {
        PdbqtAtom* atom = (*it).second;
        int serial_number = (*it).first;
        atoms_[serial_number] = atom;
    }
}

void PdbqtAtomCard::AddAtom(PdbqtAtom* atom)
{
    int index = atom->GetAtomSerialNumber();
    if (atoms_.find(index) != atoms_.end()){
	std::cout << "Warning: Atom object with index: " << index << " alrady exists. " << std::endl; 
	std::cout << "Since std::map allows only one key-value pair. This add atom attempt overrides the previous one." << std::endl; 
    }

    atoms_[index] = atom;
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbqtAtomCard::Print(std::ostream &out)
{
//    out << "Record Name: " << record_name_ << std::endl <<
    out << "_________________ Atoms _______________" << std::endl;
    for(PdbqtAtomCard::PdbqtAtomMap::iterator it = atoms_.begin(); it != atoms_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != gmml::iNotSet)
            out << (it)->first << std::endl;
        else
            out << " " << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
    out << std::endl;
}
