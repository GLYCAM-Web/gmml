#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_ASSEMBLY_BUILDER_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_ASSEMBLY_BUILDER_HPP

#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"
#include "includes/MolecularModeling/residue.hpp"

namespace CondensedSequence
{
    class AssemblyBuilder : public SequenceManipulator
    {
    public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        AssemblyBuilder(std::string inputSequence, std::string prepFilePath);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        std::vector<MolecularModeling::Residue> GenerateResidues(std::string prepFilePath);
    private:
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        inline std::map<std::string, PrepFileSpace::PrepFileResidue*>* GetPrepResidueMap() {return &prepResidueMap_;}
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////

        inline void SetPrepResidueMap(std::map<std::string, PrepFileSpace::PrepFileResidue*> prepMap) { prepResidueMap_ = prepMap;}
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void RecurveGenerateResidues(ParsedResidue *currentChild, MolecularModeling::Residue &parent, std::vector<MolecularModeling::Residue> &createdResidues); 
        void BondResiduesDeduceAtoms(MolecularModeling::Residue& parentResidue, MolecularModeling::Residue& childResidue, std::string linkageLabel);
        std::string GetGlycamResidueName(ParsedResidue &residue);
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResidueMap_;
        std::vector<MolecularModeling::Residue> createdResidues_;
    };
}
#endif