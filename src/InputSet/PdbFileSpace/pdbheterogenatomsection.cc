#include "../../../includes/InputSet/PdbFileSpace/pdbheterogenatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbHeterogenAtomSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbHeterogenAtomSection::PdbHeterogenAtomSection() : record_name_("HETATM") {}

PdbHeterogenAtomSection::PdbHeterogenAtomSection(std::stringstream &stream_block, std::string index)
{
    heterogen_atom_cards_ = PdbHeterogenAtomCardMap();
    ordered_heterogen_atom_cards_ = PdbHeterogenAtomOrderVector();
    std::string line;
    bool is_record_name_set = false;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
//            record_name_ = line.substr(0,6);
            record_name_ = "HETATM";
            gmml::Trim(record_name_);
            is_record_name_set=true;
        }
        
        std::string this_atom_name = line.substr(12, 4);
        this_atom_name = gmml::Trim(this_atom_name);
        std::string alternate_id = line.substr(16,1);
        alternate_id = gmml::Trim(alternate_id);
          // gmml::log(__LINE__, __FILE__, gmml::INF, this_atom_name + alternate_id + line.substr(17,3));
        PdbAtomCard* atom = new PdbAtomCard(line);
        if((alternate_id.empty()) || (alternate_id == "A"))
        {
          //int ch = 97 + ConvertString<int>(Split(index,"_")[1]);
          atom->SetAtomCardIndexInResidueSet(index);
          //atom->SetAtomChainId((char)ch);
          heterogen_atom_cards_[atom->GetAtomSerialNumber()] = atom;
          ordered_heterogen_atom_cards_.push_back(atom);

          getline(stream_block, line);
          temp = line;
        }
        else //There are alternate coordinates that need to be added to the correct atom.
        {
          
          for(std::vector<PdbAtomCard*>::reverse_iterator rit = ordered_heterogen_atom_cards_.rbegin(); rit != ordered_heterogen_atom_cards_.rend(); rit++)
          {
            PdbAtomCard* thisAtomCard = *rit;
            // gmml::log(__LINE__, __FILE__, gmml::INF, thisAtomCard->GetAtomName() + " " + atom->GetAtomName());
            if(thisAtomCard->GetAtomName() == atom->GetAtomName())
            {
              // gmml::log(__LINE__, __FILE__, gmml::INF, thisAtomCard->GetAtomResidueName() + " " + atom->GetAtomResidueName());
              if(thisAtomCard->GetAtomResidueName() == atom->GetAtomResidueName())
              {
                // gmml::log(__LINE__, __FILE__, gmml::INF,thisAtomCard->GetAtomChainId()  + " " + atom->GetAtomChainId());
                // if(thisAtomCard->GetAtomChainId() == atom->GetAtomChainId())
                // {
                // gmml::log(__LINE__, __FILE__, gmml::INF,thisAtomCard->GetAtomResidueSequenceNumber()  + " " + atom->GetAtomResidueSequenceNumber());
                  // if(thisAtomCard->GetAtomResidueSequenceNumber() == atom->GetAtomResidueSequenceNumber())
                  // {
                    // gmml::log(__LINE__, __FILE__, gmml::INF, "Added card");
                    thisAtomCard->AddAlternateLocation(atom);
                    // std::stringstream test;
                    // test << thisAtomCard->GetAlternateAtomCards().size();
                    // gmml::log(__LINE__, __FILE__, gmml::INF, test.str());
                    // atom->AddAlternateLocation(thisAtomCard);
                    break;
                  // }
                // }
              }
            }
          }
          
          getline(stream_block, line);
          temp = line;
          // atom->SetAtomCardIndexInResidueSet(index);
          // //atom->SetAtomChainId((char)ch);
          // heterogen_atom_cards_[atom->GetAtomSerialNumber()] = atom;
          // ordered_heterogen_atom_cards_.push_back(atom);
        }
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbHeterogenAtomSection::GetRecordName()
{
    return record_name_;
}

PdbHeterogenAtomSection::PdbHeterogenAtomCardMap PdbHeterogenAtomSection::GetHeterogenAtomCards()
{
    return heterogen_atom_cards_;
}

PdbHeterogenAtomSection::PdbHeterogenAtomOrderVector PdbHeterogenAtomSection::GetOrderedHeterogenAtomCards()
{
    return ordered_heterogen_atom_cards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbHeterogenAtomSection::SetRecordName(const std::string record_name)
{
    record_name_ = record_name;
}

void PdbHeterogenAtomSection::SetHeterogenAtoms(PdbHeterogenAtomCardMap heterogen_atoms)
{
    heterogen_atom_cards_.clear();
    for(PdbHeterogenAtomCardMap::iterator it = heterogen_atoms.begin(); it != heterogen_atoms.end(); it++)
    {
        PdbHeterogenAtom* heterogen_atom = (*it).second;
        int serial_number = (*it).first;
        heterogen_atom_cards_[serial_number] = heterogen_atom;
    }
}

void PdbHeterogenAtomSection::SetOrderedHeterogenAtoms(PdbHeterogenAtomOrderVector ordered_heterogen_atoms)
{
    ordered_heterogen_atom_cards_.clear();
    for(PdbHeterogenAtomOrderVector::iterator it = ordered_heterogen_atoms.begin(); it != ordered_heterogen_atoms.end(); it++)
    {
        PdbAtomCard* atom = (*it);
        ordered_heterogen_atom_cards_.push_back(atom);
    }
}
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbHeterogenAtomSection::Print(std::ostream &out)
{
    out << "Record Name: " << record_name_ << std::endl <<
           "________________ Heterogen Atoms ___________________" << std::endl;
    for(PdbHeterogenAtomSection::PdbHeterogenAtomCardMap::iterator it = heterogen_atom_cards_.begin(); it != heterogen_atom_cards_.end(); it++)
    {
        out << "Atom Serial Number: ";
        if((it)->first != gmml::iNotSet)
            out << (it)->first << std::endl;
        else
            out << " " << std::endl;
        (it)->second->Print();
        out << std::endl;
    }
}
