#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using PdbFileSpace::PdbSequenceAdvancedSection;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSequenceAdvancedSection::PdbSequenceAdvancedSection() {}
PdbSequenceAdvancedSection::PdbSequenceAdvancedSection(std::stringstream &stream_block)
{
    std::string line;
    getline(stream_block, line);
    std::string temp = line;
    while (!gmml::Trim(temp).empty())
    {
        PdbSequenceAdvancedCard* sequence_advanced = new PdbSequenceAdvancedCard(line);
        AddSequenceAdvancedCards(sequence_advanced);
        getline(stream_block, line);
        temp = line;
    }
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////

PdbSequenceAdvancedSection::SequenceAdvancedCardVector PdbSequenceAdvancedSection::GetSequenceAdvancedCards()
{
    return sequence_advanced_;
}

//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void PdbSequenceAdvancedSection::SetSequenceAdvancedCards(SequenceAdvancedCardVector sequence_advanced)
{
    sequence_advanced_.clear();
    for(SequenceAdvancedCardVector::iterator it = sequence_advanced.begin(); it != sequence_advanced.end(); it++)
    {
        sequence_advanced.push_back(*it);
    }
}

void PdbSequenceAdvancedSection::AddSequenceAdvancedCards(PdbSequenceAdvancedCard *sequence_advanced)
{
    sequence_advanced_.push_back(sequence_advanced);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void PdbSequenceAdvancedSection::Print(std::ostream &out)
{
    for(SequenceAdvancedCardVector::iterator it = sequence_advanced_.begin(); it != sequence_advanced_.end(); it++)
            (*it)->Print(out);
    out << std::endl;
}
