#ifndef GMML_INCLUDES_INTERNALPROGRAMS_SEQUENCE_SEQUENCE_HPP_
#define GMML_INCLUDES_INTERNALPROGRAMS_SEQUENCE_SEQUENCE_HPP_

// OG Feb 2022
// This class only exists so I can wrap it into gems without swig having to know about SequenceManipulator and the
// template Node class I don't know how to wrap that class, so I'm trying to keep it "on the other side of the line"
// when wrapping.

#include <string>
#include "includes/Abstract/absBuilder.hpp"

namespace CondensedSequence
{
    class Sequence
    {
      public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTORS                   //
        //////////////////////////////////////////////////////////
        Sequence(std::string condensedSequence = "DManp[2S,3Me]a1-6DManpa1-6[DGlcpNAcb1-4][DNeup5Aca2-6DGalpb1-"
                                                 "4DGlcpNAc[3S]b1-2DManpa1-3]DManpb1-4DGlcpNAc[6Me]b1-4DGlcpNAcb1-OH");

        //////////////////////////////////////////////////////////
        //                       ACCESSORS                      //
        //////////////////////////////////////////////////////////
        inline std::string getInputSequence()
        {
            return inputSequence_;
        }

        inline std::string getInterpretedSequence()
        {
            return interpretedSequence_;
        }

        inline std::string getIndexOrdered()
        {
            return indexOrdered_;
        }

        inline std::string getIndexLabeled()
        {
            return indexLabeled_;
        }

        //////////////////////////////////////////////////////////
        //                       MUTATORS                       //
        //////////////////////////////////////////////////////////
        inline void setInputSequence(const std::string s)
        {
            inputSequence_ = s;
        }

        inline void setInterpretedSequence(const std::string s)
        {
            interpretedSequence_ = s;
        }

        inline void setIndexOrdered(const std::string s)
        {
            indexOrdered_ = s;
        }

        inline void setIndexLabeled(const std::string s)
        {
            indexLabeled_ = s;
        }

      private:
        //////////////////////////////////////////////////////////
        //                       ATTRIBUTES                     //
        //////////////////////////////////////////////////////////
        std::string inputSequence_       = "";
        std::string interpretedSequence_ = "";
        std::string indexOrdered_        = "";
        std::string indexLabeled_        = "";
    };
} // namespace CondensedSequence
#endif /* INCLUDES_INTERNALPROGRAMS_SEQUENCE_SEQUENCE_HPP_ */
