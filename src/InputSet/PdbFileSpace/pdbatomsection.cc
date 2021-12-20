
#include "../../../includes/InputSet/PdbFileSpace/pdbatomsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbatomcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbAtomSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbAtomSection::PdbAtomSection() : recordName_("ATOM") {}

PdbAtomSection::PdbAtomSection(std::stringstream &stream_block, std::string index)
{
    atomCardMap_ = PdbAtomMap();
    orderedAtomCards_ = PdbAtomCardOrderVector();
    std::string line;
    bool is_record_name_set = false;
//    cout << stream_block.str() << std::endl;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        if(!is_record_name_set){
//            record_name_ = line.substr(0,6);
            recordName_ = "ATOM";
            gmml::Trim(recordName_);
            is_record_name_set=true;
        }
        std::string this_atom_name = line.substr(12, 4);
        this_atom_name = gmml::Trim(this_atom_name);
        std::string alternate_id = line.substr(16,1);
        alternate_id = gmml::Trim(alternate_id);
        
        PdbFileSpace::PdbAtomCard* atom_card = new PdbFileSpace::PdbAtomCard(line);
        if((alternate_id.empty()) || (alternate_id == "A"))
        {
          //int ch = 65 + ConvertString<int>(Split(index, "_")[1]);
          atom_card->SetAtomCardIndexInResidueSet(index);
          //atom->SetAtomChainId((char)ch);
          atomCardMap_[atom_card->GetAtomSerialNumber()] = atom_card;
          orderedAtomCards_.push_back(atom_card);
          
          getline(stream_block, line);
          temp = line;
        }
// OG Dec 2021. We don't want to deal with alternate coordinates in gmml for now, just take the A version and ignore the rest.
        else //These are alternate coordinates that need to be added to the correct atom
        {
//
//          for(std::vector<PdbAtomCard*>::reverse_iterator rit = ordered_atom_cards_.rbegin(); rit != ordered_atom_cards_.rend(); rit++)
//          {
//            PdbAtomCard* thisAtomCard = *rit;
//            // gmml::log(__LINE__, __FILE__, gmml::INF, thisAtomCard->GetAtomName() + " " + atom_card->GetAtomName());
//            if(thisAtomCard->GetAtomName() == atom_card->GetAtomName())
//            {
//              // gmml::log(__LINE__, __FILE__, gmml::INF, thisAtomCard->GetAtomResidueName() + " " + atom_card->GetAtomResidueName());
//              if(thisAtomCard->GetAtomResidueName() == atom_card->GetAtomResidueName())
//              {
//                // gmml::log(__LINE__, __FILE__, gmml::INF,thisAtomCard->GetAtomChainId()  + " " + atom_card->GetAtomChainId());
//                // if(thisAtomCard->GetAtomChainId() == atom_card->GetAtomChainId())
//                // {
//                // gmml::log(__LINE__, __FILE__, gmml::INF,thisAtomCard->GetAtomResidueSequenceNumber()  + " " + atom_card->GetAtomResidueSequenceNumber());
//                  // if(thisAtomCard->GetAtomResidueSequenceNumber() == atom_card->GetAtomResidueSequenceNumber())
//                  // {
//                    // gmml::log(__LINE__, __FILE__, gmml::INF, "Added card");
//                    thisAtomCard->AddAlternateLocation(atom_card);
//                    // std::stringstream test;
//                    // test << thisAtomCard->GetAlternateAtomCards().size();
//                    // gmml::log(__LINE__, __FILE__, gmml::INF, test.str());
//                    atom_card->AddAlternateLocation(thisAtomCard);
//                    break;
//                  // }
//                // }
//              }
//            }
//          }
//
          getline(stream_block, line);
          temp = line;
//          // atom_card->SetAtomCardIndexInResidueSet(index);
//          // //atom->SetAtomChainId((char)ch);
//          // atom_cards_[atom_card->GetAtomSerialNumber()] = atom_card;
//          // ordered_atom_cards_.push_back(atom_card);
        }
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
std::string PdbAtomSection::GetRecordName()
{
    return recordName_;
}
PdbAtomSection::PdbAtomMap PdbAtomSection::GetAtomCards()
{
    return atomCardMap_;
}
PdbAtomSection::PdbAtomCardOrderVector PdbAtomSection::GetOrderedAtomCards()
{
    return orderedAtomCards_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void PdbAtomSection::SetRecordName(const std::string record_name)
{
    recordName_ = record_name;
}

void PdbAtomSection::SetAtomCards(PdbAtomMap atom_cards)
{
    atomCardMap_.clear();
    for(PdbAtomMap::iterator it = atom_cards.begin(); it != atom_cards.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom_card = (*it).second;
        int serial_number = (*it).first;
        atomCardMap_[serial_number] = atom_card;
    }
}
void PdbAtomSection::SetOrderedAtomCards(PdbAtomCardOrderVector ordered_atom_cards)
{
    orderedAtomCards_.clear();
    for(PdbAtomCardOrderVector::iterator it = ordered_atom_cards.begin(); it != ordered_atom_cards.end(); it++)
    {
        PdbFileSpace::PdbAtomCard* atom_card = (*it);
        orderedAtomCards_.push_back(atom_card);
    }
}

void PdbAtomSection::InsertAtomCard(PdbAtomCard *insertionCard, PdbAtomCard* referenceCard)
{
    std::vector<PdbAtomCard*>::iterator positionOfRef = std::find(orderedAtomCards_.begin(), orderedAtomCards_.begin(), referenceCard);
    orderedAtomCards_.insert(positionOfRef, insertionCard); // Default is end, if not found. That's ok ref card defaults to nullptr so if not passed you get it on the end.
    atomCardMap_[insertionCard->GetAtomSerialNumber()] = insertionCard; // yeah well, this is a silly idea anyway.
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbAtomSection::Print(std::ostream &out)
{
    out << "Record Name: " << recordName_ << std::endl <<
           "_________________ Atoms _______________" << std::endl;
    for(PdbAtomSection::PdbAtomMap::iterator it = atomCardMap_.begin(); it != atomCardMap_.end(); it++)
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
