#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"
namespace CondensedSequence
{
    class Carbohydrate : public SequenceManipulator
    {
    public:
        //////////////////////////////////////////////////////////
        //                       CONSTRUCTOR                    //
        //////////////////////////////////////////////////////////
    	Carbohydrate(std::string inputSequence, std::string prepFilePath);
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
//        inline std::map<std::string, PrepFileSpace::PrepFileResidue*>* GetPrepResidueMap() {return &prepResidueMap_;}
//        //////////////////////////////////////////////////////////
//        //                       MUTATOR                        //
//        //////////////////////////////////////////////////////////
//        inline void SetPrepResidueMap(std::map<std::string, PrepFileSpace::PrepFileResidue*> prepMap) { prepResidueMap_ = prepMap;}
        //////////////////////////////////////////////////////////
        //                       FUNCTIONS                      //
        //////////////////////////////////////////////////////////
        void EnsureIntegralCharge(double charge);
//        void RecurveGenerateResidues(ParsedResidue *currentChild, MolecularModeling::Residue &parent, MolecularModeling::Assembly* assembly);
        void BondResiduesDeduceAtoms(cds::Residue* parentResidue, cds::Residue* childResidue);
        std::vector<std::string> GetGlycamNamesOfResidues() const;
        std::string GetGlycamResidueName(ParsedResidue *residue) const;
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        //std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResidueMap_;
    };
}
#endif
