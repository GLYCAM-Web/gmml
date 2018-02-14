#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriessection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsupersededentriescard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSupersededEntriesSection::PdbSupersededEntriesSection() {}
PdbSupersededEntriesSection::PdbSupersededEntriesSection(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    SupersededEntriesCardVector superseded_entries_cards;
    while (!Trim(temp).empty())
    {
        PdbSupersededEntriesCard* superseded_entries = new PdbSupersededEntriesCard(line);
        AddSupersededEntriesCards(superseded_entries);
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
void PdbSupersededEntriesSection::Print(ostream &out)
{
    for(SupersededEntriesCardVector::iterator it = superseded_entries_.begin(); it != superseded_entries_.end(); it++)
            (*it)->Print(out);
    out << endl;
}
