#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"


using PdbFileSpace::PdbSupersededEntriesSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSupersededEntriesSection::PdbSupersededEntriesSection() {}
PdbSupersededEntriesSection::PdbSupersededEntriesSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        PdbSupersededEntriesCard* superseded_entries = new PdbSupersededEntriesCard(line);
        AddSupersededEntriesCards(superseded_entries);
        this->Print();
        getline(stream_block, line);
        temp = line;
    }
    }

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbSupersededEntriesSection::SupersededEntriesCardVector PdbSupersededEntriesSection::GetSupersededEntriesCards(){
    return superseded_entries_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbSupersededEntriesSection::SetSupersededEntriesCards(SupersededEntriesCardVector superseded_entries){
    superseded_entries_.clear();
    for(SupersededEntriesCardVector::iterator it = superseded_entries.begin(); it != superseded_entries.end(); it++)
    {
        superseded_entries.push_back(*it);
    }
}

void PdbSupersededEntriesSection::AddSupersededEntriesCards(PdbSupersededEntriesCard *superseded_entries)
{
    superseded_entries_.push_back(superseded_entries);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSupersededEntriesSection::Print(std::ostream &out)
{
    for(SupersededEntriesCardVector::iterator it = superseded_entries_.begin(); it != superseded_entries_.end(); it++)
            (*it)->Print(out);
    out << std::endl;
}
