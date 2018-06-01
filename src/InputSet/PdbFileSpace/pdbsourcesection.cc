#include "../../../includes/InputSet/PdbFileSpace/pdbsourcesection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsourcecard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSourceSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSourceSection::PdbSourceSection() {}
PdbSourceSection::PdbSourceSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    SourceCardVector source_cards;
    while (!gmml::Trim(temp).empty())
    {
        PdbSourceCard* source = new PdbSourceCard(line);
        AddSourceCards(source);
        getline(stream_block, line);
        temp = line;
    }
    }

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbSourceSection::SourceCardVector PdbSourceSection::GetSourceCards(){
    return source_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbSourceSection::SetSourceCards(SourceCardVector source){
    source_.clear();
    for(SourceCardVector::iterator it = source.begin(); it != source.end(); it++)
    {
        source.push_back(*it);
    }
}

void PdbSourceSection::AddSourceCards(PdbSourceCard *source)
{
    source_.push_back(source);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSourceSection::Print(std::ostream &out)
{
    for(SourceCardVector::iterator it = source_.begin(); it != source_.end(); it++)
            (*it)->Print(out);
    out << std::endl;
}
