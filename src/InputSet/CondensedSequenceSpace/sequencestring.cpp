#include "../../../includes/InputSet/CondensedSequenceSpace/sequencestring.hpp"

using CondensedSequenceSpace::SequenceString;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
SequenceString::SequenceString (std::string sequence)
{
    this->SetOriginalSequence (sequence);
} // end SequenceString
SequenceString::SequenceString (const SequenceString& another_sequence)
{
    this->SetOriginalSequence (another_sequence.GetOriginalSequence ());
} // end SequenceString
//////////////////////////////////////////////////////////
//                       ACCESSOR                       //
//////////////////////////////////////////////////////////
std::string SequenceString::GetOriginalSequence () const
{
    return this->original_sequence;
} // end GetOriginalSequence
//////////////////////////////////////////////////////////
//                       MUTATOR                        //
//////////////////////////////////////////////////////////
void SequenceString::SetOriginalSequence (std::string sequence)
{
    this->original_sequence = sequence;
} // end SetOriginalSequence
//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////
// Write a function to convert the original_sequence to a certain order style.
//////////////////////////////////////////////////////////
//                       OPERATORS                      //
//////////////////////////////////////////////////////////
void SequenceString::operator=(const SequenceString& another_sequence)
{
    this->SetOriginalSequence (another_sequence.GetOriginalSequence ());
} // end operator=
//////////////////////////////////////////////////////////
//                       DISPLAY FUNCTION               //
//////////////////////////////////////////////////////////
std::ostream& operator<<(std::ostream& os, const SequenceString& sequence_string)
{
    return os << sequence_string.GetOriginalSequence ();
} // operator<<