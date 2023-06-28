#include "includes/gmml.hpp"
#include <iostream>
#include <string>

int main(int argc, char* argv[])
{
    // Testing Yao's reordering thing
    std::string sequence = argv[1];
    CondensedSequenceSpace::CondensedSequence cond_seq(sequence);
    if (cond_seq.GetIsSequenceOkay())
    {
        std::cout << "Output longest is: "
                  << cond_seq.BuildLabeledCondensedSequence(
                         CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LONGEST_CHAIN,
                         CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LONGEST_CHAIN, false)
                  << std::endl;
        std::cout << "Output lowest index is: "
                  << cond_seq.BuildLabeledCondensedSequence(
                         CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LOWEST_INDEX,
                         CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LOWEST_INDEX, false)
                  << std::endl;
        std::cout << "Output lowest index Labeled is: "
                  << cond_seq.BuildLabeledCondensedSequence(
                         CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LOWEST_INDEX,
                         CondensedSequenceSpace::CondensedSequence::Reordering_Approach::LOWEST_INDEX, true)
                  << std::endl;
    }
    return 0;
}
