#include "includes/CentralDataStructure/InternalPrograms/Sequence/sequence.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"
#include "includes/CodeUtils/logging.hpp"

using cdsCondensedSequence::SequenceManipulator;
using CondensedSequence::Sequence;

Sequence::Sequence(std::string condensedSequence)
{
    SequenceManipulator manipulator(condensedSequence);
    this->setInputSequence(condensedSequence);
    this->setInterpretedSequence(manipulator.Print());
    this->setIndexOrdered(manipulator.ReorderSequence());
    manipulator.LabelSequence();
    bool withLabels = true;
    this->setIndexLabeled(manipulator.Print(withLabels));
}
