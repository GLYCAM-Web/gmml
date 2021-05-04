#ifndef SEQUENCE_ASSEMBLY_HPP
#define SEQUENCE_ASSEMBLY_HPP

#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"
#include "includes/MolecularModeling/residue.hpp"

namespace CondensedSequence
{
    class SequenceAssembly : public SequenceManipulator
    {
    public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
        SequenceAssembly(std::string inputSequence, std::string prepFilePath);
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

        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResidueMap_;
        std::vector<MolecularModeling::Residue> createdResidues_;
    };
}
#endif