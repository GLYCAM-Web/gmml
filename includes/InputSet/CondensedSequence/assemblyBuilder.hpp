#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_ASSEMBLYBUILDER_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_ASSEMBLYBUILDER_HPP

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
        AssemblyBuilder(std::string inputSequence, std::string prepFilePath, MolecularModeling::Assembly *inputAssembly);
        //////////////////////////////////////////////////////////
        //                       ACCESSOR                       //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       MUTATOR                        //
        //////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
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
        void EnsureIntegralCharge(double charge);
        void GenerateResidues(MolecularModeling::Assembly *inputAssembly);
        void RecurveGenerateResidues(ParsedResidue *currentChild, MolecularModeling::Residue &parent, MolecularModeling::Assembly* assembly); 
        void BondResiduesDeduceAtoms(MolecularModeling::Residue& parentResidue, MolecularModeling::Residue& childResidue, std::string linkageLabel);
        std::string GetGlycamResidueName(ParsedResidue &residue);
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResidueMap_;
        std::vector<MolecularModeling::Residue> createdResidues_;
    };
}
#endif // GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_ASSEMBLYBUILDER_HPP
