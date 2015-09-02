#include "../../../includes/FileSet/CondensedSequenceSpace/condensedsequence.hpp"
#include "../../../includes/common.hpp"

using namespace std;
using namespace gmml;
using namespace CondensedSequenceSpace;

//////////////////////////////////////////////////////////
//                       CONSTRUCTOR                    //
//////////////////////////////////////////////////////////
CondensedSequence::CondensedSequence()
{
}

//////////////////////////////////////////////////////////
//                         ACCESSOR                     //
//////////////////////////////////////////////////////////
CondensedSequence::CondensedSequenceResidueVector CondensedSequence::GetResidues()
{
    return residues_;
}
CondensedSequence::CondensedSequenceTokenTypeVector CondensedSequence::GetTokens()
{
    return tokens_;
}

//////////////////////////////////////////////////////////
//                          MUTATOR                     //
//////////////////////////////////////////////////////////
void CondensedSequence::SetResidues(CondensedSequenceResidueVector residues)
{
    residues_.clear();
    for(CondensedSequenceResidueVector::iterator it = residues.begin(); it != residues.end(); it++)
    {
        residues_.push_back(*it);
    }
}
void CondensedSequence::AddResidue(CondensedSequenceResidue* residue)
{
    residues_.push_back(residue);
}
void CondensedSequence::SetTokens(CondensedSequenceTokenTypeVector tokens)
{
    tokens_.clear();
    for(CondensedSequenceTokenTypeVector::iterator it = tokens.begin(); it != tokens.end(); it++)
    {
        tokens_.push_back(*it);
    }
}
void CondensedSequence::AddToken(gmml::CondensedSequenceTokenType token)
{
    tokens_.push_back(token);
}

//////////////////////////////////////////////////////////
//                        FUNCTIONS                     //
//////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
//                      DISPLAY FUNCTION                //
//////////////////////////////////////////////////////////
void CondensedSequence::Print(ostream &out)
{}





