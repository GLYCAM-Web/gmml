#ifndef GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP
#define GMML_INCLUDES_INPUTSET_CONDENSEDSEQUENCE_CARBOHYDRATE_HPP

#include "includes/InputSet/CondensedSequence/sequenceManipulator.hpp"
#include "includes/ParameterSet/PrepFile/prepFile.hpp"
#include "includes/CentralDataStructure/Shapers/residueLinkage.hpp"
#include "includes/Abstract/absBuilder.hpp"
#include <vector>
namespace CondensedSequence
{
    class Carbohydrate : public SequenceManipulator, public Abstract::absBuilder
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
    	void Generate3DStructureFiles(std::string fileOutputDirectory = "unspecified", std::string outputFileNaming = "structure");
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
    	void ApplyDeoxy(ParsedResidue* deoxyResidue);
    	void MoveAtomsFromPrepResidueToParsedResidue(prep::PrepFile& prepResidues, ParsedResidue* parsedResidue);
    	void DerivativeChargeAdjustment(ParsedResidue* parsedResidue);
        void EnsureIntegralCharge(double charge);
//        void RecurveGenerateResidues(ParsedResidue *currentChild, MolecularModeling::Residue &parent, MolecularModeling::Assembly* assembly);
        void ConnectAndSetGeometry(cds::Residue* parentResidue, cds::Residue* childResidue);
        std::vector<std::string> GetGlycamNamesOfResidues() const;
        std::string GetGlycamResidueName(ParsedResidue *residue) const;
        void FigureOutResidueLinkages(cds::Residue* from_this_residue1, cds::Residue* to_this_residue2);
        void SetDefaultShapeUsingMetadata();
        void ResolveOverlaps();
        //////////////////////////////////////////////////////////
        //                 PRIVATE MEMBERS                      //
        //////////////////////////////////////////////////////////
        //std::map<std::string, PrepFileSpace::PrepFileResidue*> prepResidueMap_;
        std::string inputSequenceString_;
        std::vector<cds::ResidueLinkage> glycosidicLinkages_;
    };
}
#endif
