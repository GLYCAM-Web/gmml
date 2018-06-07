#include "../../../includes/InputSet/PdbFileSpace/pdbrevisiondatasection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbrevisiondatacard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbRevisionDataSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbRevisionDataSection::PdbRevisionDataSection() {}
PdbRevisionDataSection::PdbRevisionDataSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    RevisionDataCardVector revision_data_cards;
    while (!gmml::Trim(temp).empty())
    {
        PdbRevisionDataCard* revision_data = new PdbRevisionDataCard(line);
        AddRevisionDataCards(revision_data);

        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbRevisionDataSection::RevisionDataCardVector PdbRevisionDataSection::GetRevisionDataCards(){
    return revision_data_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbRevisionDataSection::SetRevisionDataCards(RevisionDataCardVector revision_data){
    revision_data_.clear();
    for(RevisionDataCardVector::iterator it = revision_data.begin(); it != revision_data.end(); it++)
    {
        revision_data.push_back(*it);
    }
}

void PdbRevisionDataSection::AddRevisionDataCards(PdbRevisionDataCard *revision_data)
{
    revision_data_.push_back(revision_data);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbRevisionDataSection::Print(std::ostream &out)
{
    for(RevisionDataCardVector::iterator it = revision_data_.begin(); it != revision_data_.end(); it++)
    {
      (*it)->Print(out);
      out << std::endl;
    }

}
