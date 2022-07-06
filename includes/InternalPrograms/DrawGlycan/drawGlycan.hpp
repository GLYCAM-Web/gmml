#ifndef GMML_INCLUDES_INTERNALPROGRAMS_DRAWGLYCAN_DRAWGLYCAN_HPP
#define GMML_INCLUDES_INTERNALPROGRAMS_DRAWGLYCAN_DRAWGLYCAN_HPP

// OG Feb 2022
// This class only exists so I can wrap it into gems without swig having to know about SequenceManipulator and the template Node class
// I don't know how to wrap that class, so I'm trying to keep it "on the other side of the line" when wrapping.

#include <string>
#include "includes/Abstract/builder.hpp"
#include "includes/InputSet/CondensedSequence/graphVizDotConfig.hpp"

namespace CondensedSequence
{
class DrawGlycan : public Abstract::Builder
{
public:
    //////////////////////////////////////////////////////////
    //                       CONSTRUCTORS                   //
    //////////////////////////////////////////////////////////
    DrawGlycan(GraphVizDotConfig configs, std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH");
    //////////////////////////////////////////////////////////
    //                       ACCESSORS                      //
    //////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    //                       FUNCTIONS                      //
    //////////////////////////////////////////////////////////
private:
    //////////////////////////////////////////////////////////
    //                       ATTRIBUTES                     //
    //////////////////////////////////////////////////////////
};
} // namespace
#endif // GMML_INCLUDES_INTERNALPROGRAMS_DRAWGLYCAN_DRAWGLYCAN_HPP
