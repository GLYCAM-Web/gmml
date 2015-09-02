#ifndef CONDENSEDSEQUENCE_HPP
#define CONDENSEDSEQUENCE_HPP

#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "../common.hpp"

namespace CondensedSequenceSpace
{
    class CondensedSequenceResidue;
    class CondensedSequence
    {
        public:
            //////////////////////////////////////////////////////////
            //                    TYPE DEFINITION                   //
            //////////////////////////////////////////////////////////
            /*! \typedef
              * List of condensed sequence residues
              */
            typedef std::vector<CondensedSequenceResidue*> CondensedSequenceResidueVector;
            /*! \typedef
              * List of condensed sequence residues
              */
            typedef std::vector<gmml::CondensedSequenceTokenType> CondensedSequenceTokenTypeVector;

            //////////////////////////////////////////////////////////
            //                       CONSTRUCTOR                    //
            //////////////////////////////////////////////////////////
            /*! \fn
              * Default constructor
              */
            CondensedSequence();

            //////////////////////////////////////////////////////////
            //                       ACCESSOR                       //
            //////////////////////////////////////////////////////////
            CondensedSequenceResidueVector GetResidues();
            CondensedSequenceTokenTypeVector GetTokens();

            //////////////////////////////////////////////////////////
            //                       MUTATOR                        //
            //////////////////////////////////////////////////////////
            void SetResidues(CondensedSequenceResidueVector residues);
            void AddResidue(CondensedSequenceResidue* residue);
            void SetTokens(CondensedSequenceTokenTypeVector tokens);
            void AddToken(gmml::CondensedSequenceTokenType token);

            //////////////////////////////////////////////////////////
            //                        FUNCTIONS                     //
            //////////////////////////////////////////////////////////

            //////////////////////////////////////////////////////////
            //                       DISPLAY FUNCTION               //
            //////////////////////////////////////////////////////////
            /*! \fn
              * A function to print out the condensed sequence contents in a structural format
              * Print out the information in a defined structure
              * @param out An output stream, the print result will be written in the given output stream
              */
            void Print(std::ostream& out = std::cout);

        private:
            //////////////////////////////////////////////////////////
            //                       ATTRIBUTES                     //
            //////////////////////////////////////////////////////////
            CondensedSequenceResidueVector residues_;
            CondensedSequenceTokenTypeVector tokens_;

    };
}

#endif // CONDENSEDSEQUENCE_HPP
