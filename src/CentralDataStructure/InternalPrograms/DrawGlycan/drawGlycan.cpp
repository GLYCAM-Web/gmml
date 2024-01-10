#include "includes/CentralDataStructure/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"

using cdsCondensedSequence::SequenceManipulator;
using CondensedSequence::DrawGlycan;

DrawGlycan::DrawGlycan(cdsCondensedSequence::GraphVizDotConfig configs, std::string condensedSequence)
{
    SequenceManipulator manipulator(condensedSequence);
    manipulator.PrintGraphViz(configs);
}
