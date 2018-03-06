#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedsection.hpp"
#include "../../../includes/InputSet/PdbFileSpace/pdbsequenceadvancedcard.hpp"
#include "../../../includes/utils.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace PdbFileSpace;
using namespace gmml;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
PdbSequenceAdvancedSection::PdbSequenceAdvancedSection() {}
PdbSequenceAdvancedSection::PdbSequenceAdvancedSection(stringstream &stream_block)
{
    string line;
    getline(stream_block, line);
    string temp = line;
    while (!Trim(temp).empty())
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
void PdbSequenceAdvancedSection::Print(ostream &out)
{
    for(SequenceAdvancedCardVector::iterator it = sequence_advanced_.begin(); it != sequence_advanced_.end(); it++)
            (*it)->Print(out);
    out << endl;
}
