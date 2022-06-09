#include "includes/InternalPrograms/DrawGlycan/drawGlycan.hpp"
#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"

using CondensedSequence::DrawGlycan;
using CondensedSequence::SequenceManipulator;


DrawGlycan::DrawGlycan(GraphVizDotConfig configs, std::string condensedSequence)
{
    SequenceManipulator manipulator(condensedSequence);
    manipulator.PrintGraphViz(configs);
}

