#include "includes/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/CentralDataStructure/CondensedSequence/sequenceManipulator.hpp"

using CondensedSequence::DrawGlycan;
using cdsCondensedSequence::SequenceManipulator;

DrawGlycan::DrawGlycan(cdsCondensedSequence::GraphVizDotConfig configs, std::string condensedSequence)
{
    SequenceManipulator manipulator(condensedSequence);
    manipulator.PrintGraphViz(configs);
}

